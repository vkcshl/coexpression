#!/usr/bin/python
import argparse
import sys
import os
import time
import traceback
import sys
import ctypes
import subprocess
from subprocess import Popen, PIPE
import os
from optparse import OptionParser
from biokbase.workspace.client import Workspace

desc1 = '''
NAME
      coex-filter-genes -- select differentially expressed genes

SYNOPSIS      
      
'''

desc2 = '''
DESCRIPTION
  coex-filter-genes provides the function to identify differentially expressed genes given an expression series/experiment. An expression series/experiment contains a list of expression samples. A expression sample is the measurement of mRNA abundance in a biological sample. The design of expression profiling usually includes replicates. The replicates allows us to differ the non-relevent expression variation and the relevent expression variation.
  The replicate information is manully extracted by KBase developers. Only a part of samples has been assigned to a replicate group. For those samples without an assignment, the variation of its expression abundance is used directly.
  filter_genes now has two methods to identify differentially expressed genes: ANOVA and lor(from limma r package). The output of this function is a list of genes
  All the data is feteched from KBase workspace. The output is stored back into KBase workspace.
'''

desc3 = '''
EXAMPLES
      Filter genes with ANOVA
      > coex-filter-genes --ws_url 'https://kbase.us/services/ws' --ws_id KBasePublicExpression  --in_id 'my_series' --out_id 'filtered_series' --filter_method anova --p_value 0.01 
      > coex-filter-genes -u 'https://kbase.us/services/ws' -w KBasePublicExpression  -i 'my_series' -o 'filtered_series' -m anova -p 0.01 
      
      Filter genes with LOR
      > coex-filter-genes --ws_url 'https://kbase.us/services/ws' --ws_id KBasePublicExpression  --in_id 'my_series' -out_id 'filtered_series' --filter_method lor --p_value 0.01 
      > coex-filter-genes -u 'https://kbase.us/services/ws' -w KBasePublicExpression  -i 'my_series' -o 'filtered_series' -m lor -p 0.01


SEE ALSO
      coex_filter

AUTHORS
Shinjae Yoo, Gang Fang, Fei He, Daifeng Wang.
'''


def filter_expression (args) :
    ###
    # download ws object and convert them to csv
    wsd = Workspace(url=args.ws_url, token=os.environ.get('KB_AUTH_TOKEN'))
    lseries = wsd.get_object({'id' : args.inobj_id,
                  'type' : 'KBaseExpression.ExpressionSeries', 
                  'workspace' : args.ws_id})['data']

    if lseries is None:
        raise COEXException("Object {} not found in workspace {}".format(args.inobj_id, args.ws_id))

    samples, sids, genome_id = {}, [], ""
    # assume only one genome id
    for gid in sorted(lseries['genome_expression_sample_ids_map'].keys()):
        genome_id = gid
        for samid in lseries['genome_expression_sample_ids_map'][gid]:
            sids.append({'ref': samid})
        samples = wsd.get_objects(sids)
        break

    cif = open(args.exp_fn, 'w')
    header = ",".join([s['data']['source_id'] for s in samples])
    cif.write(header + "\n")

    # find common gene list
    gids = set(samples[0]['data']['expression_levels'].keys())  # each sample has same gids
    for s in samples:
        gids = gids.intersection(set(s['data']['expression_levels'].keys()))
    for gid in sorted(gids):
        line = gid + ","
        line += ",".join([str(s['data']['expression_levels'][gid]) for s in samples])
        cif.write(line + "\n")
    cif.close()

    sif = open(args.rp_smp_fn, 'w')
    sample = ",".join(map(str, range(len(samples))))
    sif.write(sample + "\n")
    sif.close()

    ###
    # execute filtering
    flt_cmd_lst = ['coex_filter', "-i", args.exp_fn]
    if (args.method     is not None): 
        flt_cmd_lst.append('-m')
        flt_cmd_lst.append(args.method)
    if (args.p_value    is not None): 
        flt_cmd_lst.append('-p')
        flt_cmd_lst.append(args.p_value)
    if (args.num_genes  is not None): 
        flt_cmd_lst.append('-n')
        flt_cmd_lst.append(args.num_genes)
    if (args.flt_out_fn is not None): 
        flt_cmd_lst.append('-o')
        flt_cmd_lst.append(args.flt_out_fn)
    if (args.rp_smp_fn  is not None): 
        flt_cmd_lst.append('-s')
        flt_cmd_lst.append(args.rp_smp_fn)

    p1 = Popen(flt_cmd_lst, stdout=PIPE)
    out_str = p1.communicate()
    # print output message for error tracking
    if out_str[0] is not None : print out_str[0]
    if out_str[1] is not None : print >> sys.stderr, out_str[1]
    flt_cmd = " ".join(flt_cmd_lst)
   
    ###
    # put it back to workspace
    elm = {};
    fif = open(args.flt_out_fn, 'r')
    fif.readline(); # skip header
    
    nsamples = len(samples)
    for i in range(nsamples): elm[i] = {}
    
    for line in fif :
        line.strip();
        values = line.split(',')
        gene_id = values[0].replace("\"", "")
        for i in range(nsamples): elm[i][gene_id] = float(values[i + 1])
 
    data_list = [];
    sid_list =[];
    for i in range(nsamples) :
        samples[i]['data']['expression_levels'] = elm[i]
        if samples[i]['data']['title'] is None: samples[i]['data']['title'] = " Filtered by coex-filter-genes" 
        else : samples[i]['data']['title'] += " filtered by coex-filter-genes"
        if samples[i]['data']['description'] is None : samples[i]['data']['description'] = "Generated by " + flt_cmd
        else : samples[i]['data']['description'] += " Generated by " + flt_cmd
        samples[i]['data']['id']+=".filtered";
        samples[i]['data']['source_id']+=".filtered";
        data_list.append({'type' : 'KBaseExpression.ExpressionSample', 'data' : samples[i]['data'], 'name' : samples[i]['data']['id']})
    sv_rst = wsd.save_objects({'workspace' : args.ws_id, 'objects' : data_list})
    for i in range(nsamples):sid_list.append(str(sv_rst[i][6]) + "/" + str(sv_rst[i][0]) + "/" + str(sv_rst[i][4]))
 
    data_list = [];
    # assume only one genome id
    lseries['genome_expression_sample_ids_map'][genome_id] = sid_list
    lseries['title'] += " filtered by coex_filter for " + genome_id
    lseries['source_id'] += ".filtered"
    lseries['id'] = args.outobj_id
    data_list.append({'type' : 'KBaseExpression.ExpressionSeries', 'data' : lseries, 'name' : lseries['id'], 'meta' : {'org.series' : args.inobj_id}})
    wsd.save_objects({'workspace' : args.ws_id, 'objects' : data_list})

    if(args.del_tmps is "true") :
        os.remove(args.exp_fn)
        os.remove(args.rp_smp_fn)
        os.remove(args.flt_out_fn)
 

if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='coex-filter-genes', epilog=desc3)
    parser.add_argument('-u', '--ws_url', help='Workspace url', action='store', dest='ws_url', default='https://kbase.us/services/ws')
    parser.add_argument('-w', '--ws_id', help='Workspace id', action='store', dest='ws_id', default=None, required=True)
    parser.add_argument('-i', '--in_id', help='Input Series object id', action='store', dest='inobj_id', default=None, required=True)
    parser.add_argument('-o', '--out_id', help='Output Series object id', action='store', dest='outobj_id', default=None, required=True)
    parser.add_argument('-m', '--filter_method', help='Filtering methods (\'anova\' for ANOVA or \'lor\' for log-odd ratio', action='store', dest='method', default='anova')
    parser.add_argument('-n', '--num_genes', help='The number of genes to be selected', action='store', dest='num_genes', default=None)
    parser.add_argument('-p', '--p_value', help='The p-value cut-off', action='store', dest='p_value', default=None)
    parser.add_argument('-e', '--expression_fn', help='Expression file name (temporary file)', action='store', dest='exp_fn', default='expression.csv')
    parser.add_argument('-r', '--replicate_sample_id_fn', help='Replicate sample id file name (temporary file)', action='store', dest='rp_smp_fn', default='sample.csv')
    parser.add_argument('-f', '--filter_out_fn', help='Filtering output file name (temporary file)', action='store', dest='flt_out_fn', default='filtered.csv')
    parser.add_argument('-d', '--del_tmp_files', help='Delete temporary files', action='store', dest='del_tmps', default='true')
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    if(args.p_value is None and args.num_genes is None) :
        print "Either p_value or num_genes has to be specified";
        exit(1);

    # main loop
    filter_expression(args)
    exit(0);
