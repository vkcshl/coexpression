#!/usr/bin/env python

# standard library imports
import os
import sys
import traceback
import argparse
import json
import logging
import time
import pprint 
import string
import subprocess
from os import environ
from ConfigParser import ConfigParser

# 3rd party imports
import requests

# KBase imports
import biokbase.workspace.client 
import biokbase.Transform.script_utils as script_utils 


FVE_2_TSV = 'trns_transform_KBaseFeatureValues_ExpressionMatrix_to_TSV'
TSV_2_FVE = 'trns_transform_TSV_Exspression_to_KBaseFeatureValues_ExpressionMatrix'
RAWEXPR_DIR = 'raw_dir'
FLTRD_DIR = 'fltr_dir'
FINAL_DIR = 'final_dir'
EXPRESS_FN = 'expression.tsv'
SAMPLE_FN = 'sample.tsv'
COEX_FILTER = 'coex_filter'
FLTRD_FN = 'filtered.tsv'
FINAL_FN = 'filtered.json'
GENELST_FN = 'selected.tsv'


def empty_results(err_msg, expr, workspace_service_url, param, logger, ws):
    if 'description' not in expr: 
        expr['description'] = "Filtered Expression Matrix"
    expr['description'] += " : Empty Expression Matrix by '{0}' method; {1}".format(param['method'], err_msg)

    expr['feature_mapping'] = {}
    expr['data'] = {'row_ids' : [], 'col_ids' : [], 'values' : []}

    ws.save_objects({'workspace' : param['workspace_name'], 'objects' : [{'type' : 'KBaseFeatureValues.ExpressionMatrix',
                                                                          'data' : expr,
                                                                          'name' : (param['out_expr_object_name'])}]})

    ## Upload FeatureSet
    # Empty element feature set was not properly handled by widget, so skip it
    #fs ={'elements': {}}
    #fs['description'] = "Empty FeatureSet by '{0}' method; {1}".format(param['method'], err_msg)

    #ws.save_objects({'workspace' : param['workspace_name'], 'objects' : [{'type' : 'KBaseCollections.FeatureSet',
    #                                                                      'data' : fs,
    #                                                                      'name' : (param['out_fs_object_name'])}]})

def run_filter_genes(workspace_service_url=None, param_file = None, level=logging.INFO, logger = None):
    """
    Narrative Job Wrapper script to execute coex_filter
    
    Args:
        workspace_service_url:  A url for the KBase Workspace service 
        param_file: parameter file
        object_name: Name of the object in the workspace 
        level: Logging level, defaults to logging.INFO.
    
    Returns:
        Output is written back in WS
    
    Authors:
        Shinjae Yoo
    
    """ 

    try:
        os.makedirs(RAWEXPR_DIR)
    except:
        pass
    try:
        os.makedirs(FLTRD_DIR)
    except:
        pass
    try:
        os.makedirs(FINAL_DIR)
    except:
        pass

    if logger is None:
        logger = script_utils.stderrlogger(__file__)
    
    logger.info("Starting conversion of KBaseFeatureValues.ExpressionMatrix to TSV")
    token = os.environ.get("KB_AUTH_TOKEN")

    with open(param_file) as paramh:
      param = json.load(paramh)


    from biokbase.workspace.client import Workspace
    ws = Workspace(url=workspace_service_url, token=os.environ['KB_AUTH_TOKEN'])
    expr = ws.get_objects([{'workspace': param['workspace_name'], 'name' : param['object_name']}])[0]['data']


    cmd_dowload_cvt_tsv = [FVE_2_TSV, '--workspace_service_url', workspace_service_url, 
                                      '--workspace_name', param['workspace_name'],
                                      '--object_name', param['object_name'],
                                      '--working_directory', RAWEXPR_DIR,
                                      '--output_file_name', EXPRESS_FN
                          ]

    # need shell in this case because the java code is depending on finding the KBase token in the environment
    #  -- copied from FVE_2_TSV
    tool_process = subprocess.Popen(" ".join(cmd_dowload_cvt_tsv), stderr=subprocess.PIPE, shell=True)
    stdout, stderr = tool_process.communicate()
    
    if stdout is not None and len(stdout) > 0:
        logger.info(stdout)

    if stderr is not None and len(stderr) > 0:
        logger.info(stderr)

    logger.info("Identifying differentially expressed genes")

    ## Prepare sample file
    # detect num of columns
    with open("{0}/{1}".format(RAWEXPR_DIR, EXPRESS_FN), 'r') as f:
      fl = f.readline()
    ncol = len(fl.split('\t'))
    
    # force to use ANOVA if the number of sample is two
    if(ncol == 3): param['method'] = 'anova'

    with open("{0}/{1}".format(RAWEXPR_DIR, SAMPLE_FN), 'wt') as s:
      s.write("0")
      for j in range(1,ncol-1):
        s.write("\t{0}".format(j))
      s.write("\n")


    ## Run coex_filter
    cmd_coex_filter = [COEX_FILTER, '-i', "{0}/{1}".format(RAWEXPR_DIR, EXPRESS_FN), '-o', "{0}/{1}".format(FLTRD_DIR, FLTRD_FN),
                       '-m', param['method'], '-s', "{0}/{1}".format(RAWEXPR_DIR, SAMPLE_FN),
                       '-x', "{0}/{1}".format(RAWEXPR_DIR, GENELST_FN), '-t', 'y']
    if 'num_features' in param:
      cmd_coex_filter.append("-n")
      cmd_coex_filter.append(param['num_features'])

    if 'p_value' in param:
      cmd_coex_filter.append("-p")
      cmd_coex_filter.append(param['p_value'])

    if 'p_value' not in param and 'num_features' not in param:
      logger.error("One of p_value or num_features must be defined");
      return empty_results("One of p_value or num_features must be defined", expr,workspace_service_url, param, logger, ws)
      #sys.exit(2) #TODO: No error handling in narrative so we do graceful termination

    #if 'p_value' in param and 'num_features' in param:
    #  logger.error("Both of p_value and num_features cannot be defined together");
    #  sys.exit(3)

    tool_process = subprocess.Popen(cmd_coex_filter, stderr=subprocess.PIPE)
    stdout, stderr = tool_process.communicate()
    
    if stdout is not None and len(stdout) > 0:
        logger.info(stdout)

    if stderr is not None and len(stderr) > 0:
        logger.info(stderr)

    ## Header correction
    try:
        with open("{0}/{1}".format(FLTRD_DIR, FLTRD_FN), 'r') as ff:
            fe = ff.readlines()
        with open("{0}/{1}".format(FLTRD_DIR, FLTRD_FN), 'w') as ff:
            ff.write(fl) # use original first line that has correct header information
            fe.pop(0)
            ff.writelines(fe)
    except:
        logger.error("Output was not found");
        return empty_results("Increase p_value or specify num_features", expr,workspace_service_url, param, logger, ws)
        
    
    ## checking genelist
    with open("{0}/{1}".format(RAWEXPR_DIR, GENELST_FN),'r') as glh:
      gl = glh.readlines()
    gl = [x.strip('\n') for x in gl]

    if(len(gl) < 1) :
      logger.error("No genes are selected")
      return empty_results("Increase p_value or specify num_features", expr,workspace_service_url, param, logger, ws)
      #sys.exit(4)

    ## Upload FVE
    # change workspace to be the referenced object's workspace_name because it may not be in the same working ws due to referencing
    # Updates: change missing genome handling strategy by copying reference to working workspace
    cmd_upload_expr = [TSV_2_FVE, '--workspace_service_url', workspace_service_url, 
                                      '--object_name', param['out_expr_object_name'],
                                      '--working_directory', FINAL_DIR,
                                      '--input_directory', FLTRD_DIR,
                                      '--output_file_name', FINAL_FN
                          ]
    tmp_ws = param['workspace_name']
    if 'genome_ref' in expr:
        obj_infos = ws.get_object_info_new({"objects": [{'ref':expr['genome_ref']}]})[0]

        if len(obj_infos) < 1:
            logger.error("Couldn't find {0} from the workspace".format(expr['genome_ref']))
            raise Exception("Couldn't find {0} from the workspace".format(expr['genome_ref']))

        #tmp_ws = "{0}".format(obj_infos[7])
        logger.info("{0} => {1} / {2}".format(expr['genome_ref'], obj_infos[7], obj_infos[1]))
        if obj_infos[7] != param['workspace_name']:
            #we need to copy it from the other workspace
            try:
              logger.info("trying to copy the referenced genome object : {0}".format(expr['genome_ref']))
              ws.copy_object({'from' : {'ref' : expr['genome_ref']},'to' : {'workspace': param['workspace_name'], 'name' : obj_infos[1]}})
              # add genome_object_name only after successful copy
              cmd_upload_expr.append('--genome_object_name')
              cmd_upload_expr.append(obj_infos[1])
            except:
              # no permission or any issues... then, give up providing genome reference
              logger.info("".join(traceback.format_exc()))
              pass
        else:
            # it is local... we can simply add reference without copying genome
            cmd_upload_expr.append('--genome_object_name')
            cmd_upload_expr.append(obj_infos[1])

    # updated ws name
    cmd_upload_expr.append('--workspace_name')
    cmd_upload_expr.append(tmp_ws)

    logger.info(" ".join(cmd_upload_expr))

    tool_process = subprocess.Popen(" ".join(cmd_upload_expr), stderr=subprocess.PIPE, shell=True)
    stdout, stderr = tool_process.communicate()
    
    if stdout is not None and len(stdout) > 0:
        logger.info(stdout)

    if stderr is not None and len(stderr) > 0:
        logger.info(stderr)

    
    with open("{0}/{1}".format(FINAL_DIR,FINAL_FN),'r') as et:
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


# called only if script is run from command line
if __name__ == "__main__":	
    import sys

    parser = argparse.ArgumentParser(prog='njs-run-coex.py', 
                                     description='NJS Service Wrapper Script',
                                     epilog='Authors: Shinjae Yoo')
    parser.add_argument('-s', '--service_url', help='Service url', action='store', type=str, default='impl', nargs='?', required=False)
    parser.add_argument('-w', '--ws_url', help='Workspace url', action='store', type=str, default='https://kbase.us/services/ws/', nargs='?', required=False)
    parser.add_argument('-c','--command', help ='Command name', action='store', type=str, nargs='?', required=True)
    parser.add_argument('-p','--param_file', help ='Input parameter file name', action='store', type=str, nargs='?', required=True)
    parser.add_argument('-t','--token', help ='token', action='store', type=str, nargs='?', default=None, required=False)

    args = parser.parse_args()

    logger = script_utils.stderrlogger(__file__)
    try:
        ret_json = run_filter_genes(args.ws_url, args.param_file, logger=logger)
        
        #logger.info("Writing out JSON.")
        #with open(args.output_filename, "w") as outFile:
        #    outFile.write(json.dumps(ret_json,sort_keys = True, indent = 4))
        
   	logger.info("Execution completed.")
    except:
        logger.exception("".join(traceback.format_exc()))
        sys.exit(1)
    
    sys.exit(0)
