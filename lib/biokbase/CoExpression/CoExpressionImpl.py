#BEGIN_HEADER
import os
import sys
import traceback
import argparse
import json
import logging
import time
from  pprint import pprint
import string
import subprocess
from os import environ
from ConfigParser import ConfigParser
import re

# 3rd party imports
import requests

# KBase imports
import biokbase.workspace.client 
import biokbase.Transform.script_utils as script_utils 

def empty_results(err_msg, expr, workspace_service_url, param, logger, ws):
    if 'description' not in expr: 
        expr['description'] = "Filtered Expression Matrix"
    expr['description'] += " : Empty Expression Matrix by '{0}' method; {1}".format(param['method'], err_msg)

    expr['feature_mapping'] = {}
    expr['data'] = {'row_ids' : [], 'col_ids' : [], 'values' : []}

    ws.save_objects({'workspace' : param['workspace_name'], 'objects' : [{'type' : 'KBaseFeatureValues.ExpressionMatrix',
                                                                          'data' : expr,
                                                                          'name' : (param['out_expr_object_name'])}]})

def empty_cluster_results(err_msg, expr, workspace_service_url, param, logger, ws):

    #clrst = {'feature_clusters' : [{'id_to_pos' : {} }], 
    clrst = {'feature_clusters' : [], 
             'report' : { 
                 'checkTypeDetected' : '',
                 'checkUsed' : '',
                 'checkDescriptions' : [],
                 'checkResults' : [],
                 'messages' : [],
                 'warnings' : [],
                 'errors' : [err_msg]
                        }
            }

    ws.save_objects({'workspace' : param['workspace_name'], 'objects' : [{'type' : 'KBaseFeatureValues.FeatureClusters',
                                                                          'data' : clrst,

                                                                          'name' : (param['out_object_name'])}]})

def clean_up_expr_matrix(fn, logger):
    

#END_HEADER


class CoExpression:
    '''
    Module Name:
    CoExpression

    Module Description:
    Co-Expression Service APIs 

 This module provides services in support of the coexpression network. 
 The modules supports retrieval of the following information:
 1. Identify differentially expressed genes
 2. WGCNA clustering
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    FVE_2_TSV = 'trns_transform_KBaseFeatureValues_ExpressionMatrix_to_TSV'
    TSV_2_FVE = 'trns_transform_TSV_Exspression_to_KBaseFeatureValues_ExpressionMatrix'
    RAWEXPR_DIR = 'raw_dir'
    FLTRD_DIR = 'fltr_dir'
    CLSTR_DIR = 'clstr_dir'
    FINAL_DIR = 'final_dir'
    EXPRESS_FN = 'expression.tsv'
    SAMPLE_FN = 'sample.tsv'
    COEX_FILTER = 'coex_filter'
    COEX_CLUSTER = 'coex_cluster2'
    FLTRD_FN = 'filtered.tsv'
    CLSTR_FN = 'clusters.tsv'
    FINAL_FN = 'filtered.json'
    PVFDT_FN = 'pv_distribution.json'
    GENELST_FN = 'selected.tsv'
    __WS_URL = 'https://ci.kbase.us/services/ws'
    __HS_URL = 'https://ci.kbase.us/services/handle_service'
    __SHOCK_URL = 'https://ci.kbase.us/services/shock-api'
    __PUBLIC_SHOCK_NODE = 'true'
    logger = None
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        pprint(config)
        if 'ws_url' in config:
              self.__WS_URL = config['ws_url']
        if 'shock_url' in config:
              self.__SHOCK_URL = config['shock_url']
	if 'hs_url' in config:
	      self.__HS_URL = config['hs_url']
        if 'cluster_dir' in config:
              self.__CLSTR_DIR = config['cluster_dir']
        if 'final_dir' in config:
              self.__FINAL_DIR = config['final_dir']
        if 'coex_filter' in config:
              self.__COEX_FILTER = config['coex_filter']
        if 'coex_cluster' in config:
              self.__COEX_CLUSTER = config['coex_cluster']
	if 'force_shock_node_2b_public' in config: # expect 'true' or 'false' string
	      self.__PUBLIC_SHOCK_NODE = config['force_shock_node_2b_public']
	
        # logging
        self.logger = logging.getLogger('CoExpression')
        if 'log_level' in config:
              self.logger.setLevel(config['log_level'])
        else:
              self.logger.setLevel(logging.INFO)
        streamHandler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter("%(asctime)s - %(filename)s - %(lineno)d - %(levelname)s - %(message)s")
        formatter.converter = time.gmtime
        streamHandler.setFormatter(formatter)
        self.logger.addHandler(streamHandler)
        self.logger.info("Logger was set")
        #END_CONSTRUCTOR
        pass

    def diff_p_distribution(self, ctx, args):
        # ctx is the context object
        # return variables are: result
        #BEGIN diff_p_distribution
        try:
            os.makedirs(self.RAWEXPR_DIR)
        except:
            pass
        try:
            os.makedirs(self.FLTRD_DIR)
        except:
            pass
        try:
            os.makedirs(self.FINAL_DIR)
        except:
            pass
 
        if self.logger is None:
            self.logger = script_utils.stderrlogger(__file__)
        
        result = {}
        self.logger.info("Starting conversion of KBaseFeatureValues.ExpressionMatrix to TSV")
        token = ctx['token']
 
        eenv = os.environ.copy()
        eenv['KB_AUTH_TOKEN'] = token

        param = args
 
 
        from biokbase.workspace.client import Workspace
        ws = Workspace(url=self.__WS_URL, token=token)
        expr = ws.get_objects([{'workspace': param['workspace_name'], 'name' : param['object_name']}])[0]['data']
 
 
        cmd_dowload_cvt_tsv = [self.FVE_2_TSV, '--workspace_service_url', self.__WS_URL, 
                                          '--workspace_name', param['workspace_name'],
                                          '--object_name', param['object_name'],
                                          '--working_directory', self.RAWEXPR_DIR,
                                          '--output_file_name', self.EXPRESS_FN
                              ]
 
        # need shell in this case because the java code is depending on finding the KBase token in the environment
        #  -- copied from FVE_2_TSV
        tool_process = subprocess.Popen(" ".join(cmd_dowload_cvt_tsv), stderr=subprocess.PIPE, shell=True, env=eenv)
        stdout, stderr = tool_process.communicate()
        
        if stdout is not None and len(stdout) > 0:
            self.logger.info(stdout)
 
        if stderr is not None and len(stderr) > 0:
            self.logger.info(stderr)
 
        self.logger.info("Identifying differentially expressed genes")
 
        ## Prepare sample file
        # detect num of columns
        with open("{0}/{1}".format(self.RAWEXPR_DIR, self.EXPRESS_FN), 'r') as f:
          fl = f.readline()
        ncol = len(fl.split('\t'))
        
        # force to use ANOVA if the number of sample is two
        if(ncol == 3): param['method'] = 'anova'
 
        with open("{0}/{1}".format(self.RAWEXPR_DIR, self.SAMPLE_FN), 'wt') as s:
          s.write("0")
          for j in range(1,ncol-1):
            s.write("\t{0}".format(j))
          s.write("\n")
 
 
        ## Run coex_filter
        cmd_coex_filter = [self.COEX_FILTER, '-i', "{0}/{1}".format(self.RAWEXPR_DIR, self.EXPRESS_FN), '-o', "{0}/{1}".format(self.FLTRD_DIR, self.FLTRD_FN),
                           '-m', param['method'], '-s', "{0}/{1}".format(self.RAWEXPR_DIR, self.SAMPLE_FN),
                           '-x', "{0}/{1}".format(self.RAWEXPR_DIR, self.GENELST_FN), '-t', 'y', '-j', self.PVFDT_FN]
        if 'num_features' in param:
          cmd_coex_filter.append("-n")
          cmd_coex_filter.append(str(param['num_features']))
 
        if 'p_value' in param:
          cmd_coex_filter.append("-p")
          cmd_coex_filter.append(str(param['p_value']))
 
 
        tool_process = subprocess.Popen(cmd_coex_filter, stderr=subprocess.PIPE)
        stdout, stderr = tool_process.communicate()
        
        if stdout is not None and len(stdout) > 0:
            self.logger.info(stdout)
 
        if stderr is not None and len(stderr) > 0:
            self.logger.info(stderr)
 
        ## loading pvalue distribution FDT
        pvfdt = {};
        with open(self.PVFDT_FN, 'r') as myfile:
           pvfdt = json.load(myfile)
 
 
        fig_properties = {"xlabel" : "-log2(p-value)", "ylabel" : "Number of features", "xlog_mode" : "-log2", "ylog_mode" : "none", "title" : "Histogram of P-values", "plot_type" : "histogram"}
        #"data_ref" : "4997/1/1"
        sstatus = ws.save_objects({'workspace' : param['workspace_name'], 'objects' : [{'type' : 'MAK.FloatDataTable',
                                                                              'data' : pvfdt,
                                                                              'name' : (param['out_data_object_name'])}]})

        data_ref = "{0}/{1}/{2}".format(sstatus[0][6], sstatus[0][0], sstatus[0][4])
        fig_properties['data_ref'] = data_ref

        sstatus = ws.save_objects({'workspace' : param['workspace_name'], 'objects' : [{'type' : 'MAK.FloatDataTable',
                                                                              'data' : fig_properties,
                                                                              'name' : (param['out_figure_object_name'])}]})
 
        result = fig_properties
        #END diff_p_distribution

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method diff_p_distribution return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def filter_genes(self, ctx, args):
        # ctx is the context object
        # return variables are: result
        #BEGIN filter_genes
        try:
            os.makedirs(self.RAWEXPR_DIR)
        except:
            pass
        try:
            os.makedirs(self.FLTRD_DIR)
        except:
            pass
        try:
            os.makedirs(self.FINAL_DIR)
        except:
            pass
 
        if self.logger is None:
            self.logger = script_utils.stderrlogger(__file__)
        
        result = {}
        self.logger.info("Starting conversion of KBaseFeatureValues.ExpressionMatrix to TSV")
        token = ctx['token']
 
        eenv = os.environ.copy()
        eenv['KB_AUTH_TOKEN'] = token

        param = args
 
 
        from biokbase.workspace.client import Workspace
        ws = Workspace(url=self.__WS_URL, token=token)
        expr = ws.get_objects([{'workspace': param['workspace_name'], 'name' : param['object_name']}])[0]['data']
 
 
        cmd_dowload_cvt_tsv = [self.FVE_2_TSV, '--workspace_service_url', self.__WS_URL, 
                                          '--workspace_name', param['workspace_name'],
                                          '--object_name', param['object_name'],
                                          '--working_directory', self.RAWEXPR_DIR,
                                          '--output_file_name', self.EXPRESS_FN
                              ]
 
        # need shell in this case because the java code is depending on finding the KBase token in the environment
        #  -- copied from FVE_2_TSV
        tool_process = subprocess.Popen(" ".join(cmd_dowload_cvt_tsv), stderr=subprocess.PIPE, shell=True, env=eenv)
        stdout, stderr = tool_process.communicate()
        
        if stdout is not None and len(stdout) > 0:
            self.logger.info(stdout)
 
        if stderr is not None and len(stderr) > 0:
            self.logger.info(stderr)
 
        self.logger.info("Identifying differentially expressed genes")
 
        ## Prepare sample file
        # detect num of columns
        with open("{0}/{1}".format(self.RAWEXPR_DIR, self.EXPRESS_FN), 'r') as f:
          fl = f.readline()
        ncol = len(fl.split('\t'))
        
        # force to use ANOVA if the number of sample is two
        if(ncol == 3): param['method'] = 'anova'
 
        with open("{0}/{1}".format(self.RAWEXPR_DIR, self.SAMPLE_FN), 'wt') as s:
          s.write("0")
          for j in range(1,ncol-1):
            s.write("\t{0}".format(j))
          s.write("\n")
 
 
        ## Run coex_filter
        cmd_coex_filter = [self.COEX_FILTER, '-i', "{0}/{1}".format(self.RAWEXPR_DIR, self.EXPRESS_FN), '-o', "{0}/{1}".format(self.FLTRD_DIR, self.FLTRD_FN),
                           '-m', param['method'], '-s', "{0}/{1}".format(self.RAWEXPR_DIR, self.SAMPLE_FN),
                           '-x', "{0}/{1}".format(self.RAWEXPR_DIR, self.GENELST_FN), '-t', 'y']
        if 'num_features' in param:
          cmd_coex_filter.append("-n")
          cmd_coex_filter.append(str(param['num_features']))
 
        if 'p_value' in param:
          cmd_coex_filter.append("-p")
          cmd_coex_filter.append(str(param['p_value']))
 
        if 'p_value' not in param and 'num_features' not in param:
          self.logger.error("One of p_value or num_features must be defined");
          return empty_results("One of p_value or num_features must be defined", expr,self.__WS_URL, param, self.logger, ws)
          #sys.exit(2) #TODO: No error handling in narrative so we do graceful termination
 
        #if 'p_value' in param and 'num_features' in param:
        #  self.logger.error("Both of p_value and num_features cannot be defined together");
        #  sys.exit(3)
 
        tool_process = subprocess.Popen(cmd_coex_filter, stderr=subprocess.PIPE)
        stdout, stderr = tool_process.communicate()
        
        if stdout is not None and len(stdout) > 0:
            self.logger.info(stdout)
 
        if stderr is not None and len(stderr) > 0:
            self.logger.info(stderr)
 
        ## Header correction
        try:
            with open("{0}/{1}".format(self.FLTRD_DIR, self.FLTRD_FN), 'r') as ff:
                fe = ff.readlines()
            with open("{0}/{1}".format(self.FLTRD_DIR, self.FLTRD_FN), 'w') as ff:
                ff.write(fl) # use original first line that has correct header information
                fe.pop(0)
                ff.writelines(fe)
        except:
            self.logger.error("Output was not found");
            return empty_results("Increase p_value or specify num_features", expr,self.__WS_URL, param, self.logger, ws)
            
        
        ## checking genelist
        with open("{0}/{1}".format(self.RAWEXPR_DIR, self.GENELST_FN),'r') as glh:
          gl = glh.readlines()
        gl = [x.strip('\n') for x in gl]
 
        if(len(gl) < 1) :
          self.logger.error("No genes are selected")
          return empty_results("Increase p_value or specify num_features", expr,self.__WS_URL, param, self.logger, ws)
          #sys.exit(4)
 
        ## Upload FVE
        # change workspace to be the referenced object's workspace_name because it may not be in the same working ws due to referencing
        # Updates: change missing genome handling strategy by copying reference to working workspace
        cmd_upload_expr = [self.TSV_2_FVE, '--workspace_service_url', self.__WS_URL, 
                                          '--object_name', param['out_expr_object_name'],
                                          '--working_directory', self.FINAL_DIR,
                                          '--input_directory', self.FLTRD_DIR,
                                          '--output_file_name', self.FINAL_FN
                              ]
        tmp_ws = param['workspace_name']
        if 'genome_ref' in expr:
            obj_infos = ws.get_object_info_new({"objects": [{'ref':expr['genome_ref']}]})[0]
 
            if len(obj_infos) < 1:
                self.logger.error("Couldn't find {0} from the workspace".format(expr['genome_ref']))
                raise Exception("Couldn't find {0} from the workspace".format(expr['genome_ref']))
 
            #tmp_ws = "{0}".format(obj_infos[7])
            self.logger.info("{0} => {1} / {2}".format(expr['genome_ref'], obj_infos[7], obj_infos[1]))
            if obj_infos[7] != param['workspace_name']:
                #we need to copy it from the other workspace
                try:
                  self.logger.info("trying to copy the referenced genome object : {0}".format(expr['genome_ref']))
                  ws.copy_object({'from' : {'ref' : expr['genome_ref']},'to' : {'workspace': param['workspace_name'], 'name' : obj_infos[1]}})
                  # add genome_object_name only after successful copy
                  cmd_upload_expr.append('--genome_object_name')
                  cmd_upload_expr.append(obj_infos[1])
                except:
                  # no permission or any issues... then, give up providing genome reference
                  self.logger.info("".join(traceback.format_exc()))
                  pass
            else:
                # it is local... we can simply add reference without copying genome
                cmd_upload_expr.append('--genome_object_name')
                cmd_upload_expr.append(obj_infos[1])
 
        # updated ws name
        cmd_upload_expr.append('--workspace_name')
        cmd_upload_expr.append(tmp_ws)
 
        self.logger.info(" ".join(cmd_upload_expr))
 
        tool_process = subprocess.Popen(" ".join(cmd_upload_expr), stderr=subprocess.PIPE, shell=True, env=eenv)
        stdout, stderr = tool_process.communicate()
        
        if stdout is not None and len(stdout) > 0:
            self.logger.info(stdout)
 
        if stderr is not None and len(stderr) > 0:
            self.logger.info(stderr)
 
        
        with open("{0}/{1}".format(self.FINAL_DIR,self.FINAL_FN),'r') as et:
          eo = json.load(et)
 
        if 'description' not in expr: 
            expr['description'] = "Filtered Expression Matrix"
        expr['description'] += " : Filtered by '{1}' method ".format(expr['description'], param['method'])
 
        if 'feature_mapping' in expr and 'feature_mapping' in eo:
            expr['feature_mapping'] = eo['feature_mapping']
        expr['data'] = eo['data']
 
        ws.save_objects({'workspace' : param['workspace_name'], 'objects' : [{'type' : 'KBaseFeatureValues.ExpressionMatrix',
                                                                              'data' : expr,
                                                                              'name' : (param['out_expr_object_name'])}]})
 
        ## Upload FeatureSet
        fs ={'elements': {}}
        fs['description'] = "FeatureSet identified by filtering method '{0}' ".format(param['method'])
 
        fs['description'] += "from {0}/{1}".format(param['workspace_name'], param['object_name'])
 
        for g in gl:
          if 'genome_ref' in expr:
            fs['elements'][g] = [expr['genome_ref']]
          else:
            fs['elements'][g] = []
 
        ws.save_objects({'workspace' : param['workspace_name'], 'objects' : [{'type' : 'KBaseCollections.FeatureSet',
                                                                              'data' : fs,
                                                                              'name' : (param['out_fs_object_name'])}]})
        result = {'workspace_name' : param['workspace_name'], 'out_expr_object_name' : param['out_expr_object_name'], 'out_fs_object_name' : param['out_fs_object_name']}
        #END filter_genes

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method filter_genes return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def const_coex_net_clust(self, ctx, args):
        # ctx is the context object
        # return variables are: result
        #BEGIN const_coex_net_clust
        try:
            os.makedirs(self.RAWEXPR_DIR)
        except:
            pass
        try:
            os.makedirs(self.CLSTR_DIR)
        except:
            pass
        try:
            os.makedirs(self.FINAL_DIR)
        except:
            pass
 
        if self.logger is None:
            self.logger = script_utils.stderrlogger(__file__)
        
        result = {}
        self.logger.info("Starting conversion of KBaseFeatureValues.ExpressionMatrix to TSV")
        token = ctx['token']

        param = args
 
        from biokbase.workspace.client import Workspace
        ws = Workspace(url=self.__WS_URL, token=token)
        expr = ws.get_objects([{'workspace': param['workspace_name'], 'name' : param['object_name']}])[0]['data']
 
 
        eenv = os.environ.copy()
        eenv['KB_AUTH_TOKEN'] = token
        cmd_dowload_cvt_tsv = [self.FVE_2_TSV, '--workspace_service_url', self.__WS_URL, 
                                          '--workspace_name', param['workspace_name'],
                                          '--object_name', param['object_name'],
                                          '--working_directory', self.RAWEXPR_DIR,
                                          '--output_file_name', self.EXPRESS_FN
                              ]
 
        # need shell in this case because the java code is depending on finding the KBase token in the environment
        #  -- copied from FVE_2_TSV
        tool_process = subprocess.Popen(" ".join(cmd_dowload_cvt_tsv), stderr=subprocess.PIPE, shell=True, env=eenv)
        stdout, stderr = tool_process.communicate()
        
        if stdout is not None and len(stdout) > 0:
            self.logger.info(stdout)
 
        if stderr is not None and len(stderr) > 0:
            self.logger.info(stderr)
            #raise Exception(stderr)
 
        self.logger.info("Coexpression clustering analysis")
 
        ## Prepare sample file
        # detect num of columns
        with open("{0}/{1}".format(self.RAWEXPR_DIR, self.EXPRESS_FN), 'r') as f:
          fl = f.readline()
        ncol = len(fl.split('\t'))
        
        with open("{0}/{1}".format(self.RAWEXPR_DIR, self.SAMPLE_FN), 'wt') as s:
          s.write("0")
          for j in range(1,ncol-1):
            s.write("\t{0}".format(j))
          s.write("\n")
 
 
        ## Run coex_cluster
        cmd_coex_cluster = [self.COEX_CLUSTER, '-t', 'y',
                           '-i', "{0}/{1}".format(self.RAWEXPR_DIR, self.EXPRESS_FN), 
                           '-o', "{0}/{1}".format(self.CLSTR_DIR, self.CLSTR_FN)]
 
        for p in ['net_method', 'minRsq', 'maxmediank', 'maxpower', 'clust_method', 'minModuleSize', 'detectCutHeight']:
           if p in param:
             cmd_coex_cluster.append("--{0}".format(p))
             cmd_coex_cluster.append(str(param[p]))
  
 
        #sys.exit(2) #TODO: No error handling in narrative so we do graceful termination
 
        #if 'p_value' in param and 'num_features' in param:
        #  self.logger.error("Both of p_value and num_features cannot be defined together");
        #  sys.exit(3)
 
        tool_process = subprocess.Popen(cmd_coex_cluster, stderr=subprocess.PIPE)
        stdout, stderr = tool_process.communicate()
        
        if stdout is not None and len(stdout) > 0:
            self.logger.info(stdout)
 
        if stderr is not None and len(stderr) > 0:
            if re.search(r'^There were \d+ warnings \(use warnings\(\) to see them\)', stderr):
              self.logger.info(stderr)
            else:
              self.logger.error(stderr)
              raise Exception(stderr)
 
        
        # build index for gene list
        pos_index ={expr['data']['row_ids'][i]: i for i in range(0, len(expr['data']['row_ids']))}
 
 
        # parse clustering results
        cid2genelist = {}
        with open("{0}/{1}".format(self.CLSTR_DIR, self.CLSTR_FN),'r') as glh:
            glh.readline() # skip header
            for line in glh:
                gene, cluster = line.replace('"','').split("\t")
                if cluster not in cid2genelist:
                    cid2genelist[cluster] = []
                cid2genelist[cluster].append(gene)
 
        if(len(cid2genelist) < 1) :
          self.logger.error("Clustering failed")
          return empty_results("Error: No cluster output", expr,self.__WS_URL, param, self.logger, ws)
          #sys.exit(4)
 
        self.logger.info("Uploading the results onto WS")
        feature_clusters = []
        for cluster in cid2genelist:
            feature_clusters.append( { "id_to_pos" : { gene : pos_index[gene] for gene in cid2genelist[cluster]}})
                
 
        ## Upload Clusters
        feature_clusters ={"original_data": "{0}/{1}".format(param['workspace_name'],param['object_name']),
                           "feature_clusters": feature_clusters}
 
        ws.save_objects({'workspace' : param['workspace_name'], 'objects' : [{'type' : 'KBaseFeatureValues.FeatureClusters',
                                                                          'data' : feature_clusters,
                                                                          'name' : (param['out_object_name'])}]})
        result = {'workspace_name' : param['workspace_name'], 'out_object_name' : param['out_object_name']}
        #END const_coex_net_clust

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method const_coex_net_clust return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]
