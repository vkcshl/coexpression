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

    if 'num_features' not in param and 'p_value' in param:
      cmd_coex_filter.append("-p")
      cmd_coex_filter.append(param['p_value'])

    if 'p_value' not in param and 'num_features' not in param:
      logger.error("One of p_value or num_features must be defined");
      sys.exit(2)

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
    with open("{0}/{1}".format(FLTRD_DIR, FLTRD_FN), 'r') as ff:
        fe = ff.readlines()
    with open("{0}/{1}".format(FLTRD_DIR, FLTRD_FN), 'w') as ff:
        ff.write(fl) # use original first line that has correct header information
        fe.pop(0)
        ff.writelines(fe)
    

    ## Upload FVE
    from biokbase.workspace.client import Workspace
    ws = Workspace(url=workspace_service_url, token=os.environ['KB_AUTH_TOKEN'])
    expr = ws.get_objects([{'workspace': param['workspace_name'], 'name' : param['object_name']}])[0]['data']
    
    cmd_upload_expr = [TSV_2_FVE, '--workspace_service_url', workspace_service_url, 
                                      '--workspace_name', param['workspace_name'],
                                      '--object_name', "{0}_filtered".format(param['object_name']),
                                      '--working_directory', FINAL_DIR,
                                      '--input_directory', FLTRD_DIR,
                                      '--output_file_name', FINAL_FN
                          ]
    if 'genome_ref' in expr:
        cmd_upload_expr.append('--genome_object_name')
        cmd_upload_expr.append(expr['genome_ref'])

    tool_process = subprocess.Popen(" ".join(cmd_upload_expr), stderr=subprocess.PIPE, shell=True)
    stdout, stderr = tool_process.communicate()
    
    if stdout is not None and len(stdout) > 0:
        logger.info(stdout)

    if stderr is not None and len(stderr) > 0:
        logger.info(stderr)

    
    with open("{0}/{1}".format(FINAL_DIR,FINAL_FN),'r') as et:
      eo = json.load(et)

    if 'description' in expr: expr['description'] = "{0}, coex_filter by {1}".format(expr['description'], " ".join(cmd_coex_filter))
    if 'feature_mapping' in expr:
        expr['feature_mapping'] = eo['feature_mapping']
    expr['data'] = eo['data']

    ws.save_objects({'workspace' : param['workspace_name'], 'objects' : [{'type' : 'KBaseFeatureValues.ExpressionMatrix',
                                                                          'data' : expr,
                                                                          'name' : (param['out_expr_object_name'])}]})

    ## Upload FeatureSet
    fs ={'description':'Differentially expressed genes generated by {0}'.format(" ".join(cmd_coex_filter)),
         'elements': {}}
    
    with open("{0}/{1}".format(RAWEXPR_DIR, GENELST_FN),'r') as glh:
      gl = glh.readlines()
    gl = [x.strip('\n') for x in gl]

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
