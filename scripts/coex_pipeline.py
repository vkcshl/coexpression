#!/usr/bin/env python

import subprocess
import sys
import argparse
import os


# GENERATE THE HELP DOCUMENTATION #
desc = \
"""
NAME
    coex_pipeline -- pipeline to identify highly differentially expressed genes, construct gene co-expression network, and cluster genes into co-expression modules

SYNOPSIS
    coex_pipeline [-ifqnusmcvrkptelzwoxy]

DESCRIPTION
    Pipeline for computing gene coexpression network and clustering genes into co-expression modules using:
    coex_filter.r => coex_network.r => (coex_cluster.r | coex_cluster2.r)
    
    Please read the help description for the above functions for further help.
    
    NOTE: It is highly recommended that you use the long options to avoid confusion, since there are so many of them.
    
    Inputs: none, unless you do not specify '--sample_index', and then you will have to answer (y|n) to the prompt(s).  It is highly recommended that you input the '--sample_index' option.
    
    Outputs: command line arguments that were run, the output of each step in the pipeline, and writes specified files to disk.
    
    Options:
    --skip_filter             skips coex_filter.r
    --skip_network            skips coex_net.r
    --skip_cluster            skips coex_cluster.r
    --skip_cluster2_hclust    skips coex_cluster2.r using hclust method
    --skip_cluster2_wgcna     skips coex_cluster2.r using wgcna method
    
    Python-generated argument list:
"""

epil = \
"""
EXAMPLES
    $ ./coex_pipeline.py -d arabidopsis/cultured-protoplast-cell -i data.csv --filter_p 0.0001 --filter_method a --sample_index sample_id.csv --filter_output datafilter.csv --corr_threshold 0.8 --network_output coexpression_network_edge_list.csv --min_module_size 50 --max_median_k 2000 --hclust_output coex_module_h.csv --wgcna_output coex_module_w.csv
    $ ./coex_pipeline.py -d arabidopsis/fruit -i data.csv --filter_p 0.05 --sample_index sample_id.csv --filter_output datafilter.csv --corr_threshold 0.8 --network_output coexpression_network_edge_list.csv --min_module_size 50 --max_power 100 --minRsq 0.75 --hclust_output coex_module_h.csv --skip_cluster2_wgcna
    $ ./coex_pipeline.py -d poplar/xylem -i data.csv

SEE ALSO
    coex_filter.r
    coex_net.r
    coex_cluster.r
    coex_cluster2.r

AUTHORS
    Eric Pan, Daifeng Wang, Gang Fang, Mark Gerstein, Yale University
"""

# OBTAIN AND PARSE THE ARGUMENTS #

parser = argparse.ArgumentParser(description = desc, epilog = epil, formatter_class = argparse.RawTextHelpFormatter, usage = argparse.SUPPRESS)

# Set the working directory for this script
parser.add_argument('-d', '--directory', dest = 'directory', required = False, help = 'Directory containing the input file and where you want the files saved.  Leaving this out means that you must specify locations in the input and output files, or all intermediate files will be saved in your current directory, replacing the existing files.')

# Arguments for coex_filter.r
parser.add_argument('-i', '--input', dest = 'inputFile', required = False, default = "data.csv", help = "Input file name that stores original gene expression data matrix. Each row corresponds to a gene, and each column corresponds to a replicate. The column names are replicate names. The row names are gene names.")
parser.add_argument('-f', '--filter_method', dest = 'filterMethod', required = False, help = "Method to identify highly differential expressed genes. Choices include anova and log-odd ratio (lor).")
parser.add_argument('-q', '--filter_p', dest = 'filterP', required = False, help = "Maximum p-value up to which genes significantly highly differentially express")
parser.add_argument('-n', '--top_number', dest = 'topNumber', required = False, help = "Number of top genes that most highly differentially express")
parser.add_argument('-u', '--filter_output', dest = 'filterOutput', required = False, default = 'datafilter.csv', help = "Output file name that stores gene expression data matrix including selected genes only.  NOTE: This is the same as the input to coex_net.r")
parser.add_argument('-s', '--sample_index', dest = 'sampleIndex', required = False, help = "RECOMMENDED unless no replicates in data: The file that stores a numeric vector to indicate the sample indices that replicates correspond. For example, if sample_index is equal to c(1,1,2,2,3,3), the first two columns of data correspond to Sample 1, the third and fourth columns correspond to Sample 2, and the fifth and sixth columns correspond to Sample 3. Default is that every two columns correspond to a different sample.")

# Arguments for coex_net.r
parser.add_argument('-m', '--network_method', dest = 'networkMethod', required = False, help = "Method to identify highly differential expressed genes. Choices include anova and log-odd ratio (lor).")
parser.add_argument('-a', '--gene_list_1', dest = 'geneList1', required = False, help = "The first set of the genes of interest.  Leaving this out means all of the genes will be used.")
parser.add_argument('-b', '--gene_list_2', dest = 'geneList2', required = False, help = "The second set of the genes of interest.  Leaving this out means all of the genes will be used.")
parser.add_argument('-c', '--corr_threshold', dest = 'corrThreshold', required = False, help = "Maximum threshold to set elements of Pearson correlation matrix to be zero when method='simple'.")
parser.add_argument('-v', '--p_threshold', dest = 'pThreshold', required = False, help = "Maximum p-value up to which genes significantly highly differentially express")
parser.add_argument('-r', '--minRsq', dest = 'minRsq', required = False, help = "Minimum threshold for R2 that measures the fitness of gene co-expression network to scale-free topology in WGCNA. See pickSoftThreshold() of WGCNA for details.")
parser.add_argument('-k', '--max_median_k', dest = 'maxMedianK', required = False, help = "Maximum median connections for genes in network. See pickSoftThreshold() of WGCNA for details.")
parser.add_argument('-p', '--max_power', dest = 'maxPower', required = False, help = "Maximum power to decide the soft threshold. See pickSoftThreshold() of WGCNA for details.")
parser.add_argument('-t', '--output_type', dest = 'outputType', required = False, help = "Type of output. When method = 'edge' or 'e', the function outputs edge list of gene co-expression network consisting of gene pairs with Pearson correlation greater than corr_threshold. When method = 'adjmat' or 'a', the function outputs the adjacency matrix of constructed gene co-expression network.")
parser.add_argument('-e', '--network_output', dest = 'networkOutput', required = False, default = 'coexpression_network_edge_list.csv', help = "Output file name.  NOTE: This is the same as the input to coex_cluster.r or coex_cluster2.r, depending on what you selected for '--output_type'")

# Arguments for coex_cluster.r
parser.add_argument('-l', '--cluster_method', dest = 'clusterMethod', required = False, help = "Method to cluster genes into co-expression modules. When method = 'hclust' or 'h',the function uses hierarchical clustering. When method = 'WGCNA' or 'W', the function uses WGCNA. If users want to use 'hclust', the original adjacent matrix without any cutcoff (e.g., -c -1.2) is recommended.")
parser.add_argument('-z', '--min_module_size', dest = 'minModuleSize', required = False, help = "Minimum size of modules when method is WGCNA.")
parser.add_argument('-w', '--detect_cut_height',  dest = 'detectCutHeight', required = False, help = "Maximum heights to join modules in clustering.")
parser.add_argument('-o', '--output', dest = 'outputFile', required = False, default = 'coexpression_modules.csv', help = "Output file name that stores the clustering results. The first column includes gene names, and the second column includes their modules.")

# Arguments for coex_cluster2.r
parser.add_argument('-x', '--hclust_output', dest = 'hclustOutput', required = False, default = 'coex_module_h.csv', help = "Output file name that stores the clustering results from the hierarchial clustering method.")
parser.add_argument('-y', '--wgcna_output', dest = 'wgcnaOutput', required = False, default = 'coex_module_w.csv', help = "Output file name that stores the clustering results from the WGCNA method.")

[args, others] = parser.parse_known_args()
args = vars(args)

# Set up directories 
directory = args['directory']
if (directory is not None):
  lastChar = directory[-1:]
  if (lastChar != '/' and lastChar != '\\'):
    directory += '/'
    args['directory'] = directory    # this line is most likely unnecessary
  if (args['inputFile']   is not None): args['inputFile'] = directory + args['inputFile']
  if (args['sampleIndex'] is not None): args['sampleIndex'] = directory + args['sampleIndex']
  if (args['geneList1']   is not None): args['geneList1'] = directory + args['geneList1']
  if (args['geneList2']   is not None): args['geneList2'] = directory + args['geneList2']
  args['filterOutput']  = directory + args['filterOutput']
  args['networkOutput'] = directory + args['networkOutput']
  args['outputFile']    = directory + args['outputFile']
  args['hclustOutput']  = directory + args['hclustOutput']
  args['wgcnaOutput']   = directory + args['wgcnaOutput']



# GENERATE THE COMMANDS FOR THE DESIRED R FUNCTIONS #

if (not('--skip_filter' in others)):
  # Create the command line for coex_filter.r
  filterCommand = 'coex_filter -i ' + args['inputFile']
  if (args['filterMethod'] is not None): filterCommand += (' -m ' + args['filterMethod'])
  if (args['filterP']      is not None): filterCommand += (' -p ' + args['filterP'])
  if (args['topNumber']    is not None): filterCommand += (' -n ' + args['topNumber'])
  if (args['filterOutput'] is not None): filterCommand += (' -o ' + args['filterOutput'])       # somewhat superfluous, but looks nice
  if (args['sampleIndex']  is not None): filterCommand += (' -s ' + args['sampleIndex'])
  if ('--human_input'   in others):      filterCommand += (' -u y')
  if ('--no_replicates' in others):      filterCommand += (' -r n')
  if ('--not_default'   in others):      filterCommand += (' -d n')

if (not('--skip_network' in others)):
  # Create the command line for coex_net.r
  netCommand = 'coex_net -i ' + args['filterOutput']
  if (args['networkMethod'] is not None): netCommand += (' -m ' + args['networkMethod'])
  if (args['geneList1']     is not None): netCommand += (' -a ' + args['geneList1'])
  if (args['geneList2']     is not None): netCommand += (' -b ' + args['geneList2'])
  if (args['corrThreshold'] is not None): netCommand += (' -c ' + args['corrThreshold'])
  if (args['pThreshold']    is not None): netCommand += (' -q ' + args['pThreshold'])
  if (args['minRsq']        is not None): netCommand += (' -r ' + args['minRsq'])
  if (args['maxMedianK']    is not None): netCommand += (' -k ' + args['maxMedianK'])
  if (args['maxPower']      is not None): netCommand += (' -p ' + args['maxPower'])
  if (args['outputType']    is not None): netCommand += (' -t ' + args['outputType'])
  if (args['networkOutput'] is not None): netCommand += (' -o ' + args['networkOutput'])        # somewhat superfluous, but looks nice

if (not('--skip_cluster' in others)):
  # Create the command line for coex_cluster.r
  clusterCommand = 'coex_cluster -i ' + args['networkOutput']
  if (args['clusterMethod']   is not None): clusterCommand += (' -m ' + args['clusterMethod'])
  if (args['minModuleSize']   is not None): clusterCommand += (' -s ' + args['minModuleSize'])
  if (args['detectCutHeight'] is not None): clusterCommand += (' -d ' + args['detectCutHeight'])
  if (args['outputFile']      is not None): clusterCommand += (' -o ' + args['outputFile'])     # somewhat superfluous, but looks nice

if (not('--skip_hclust' in others)):
  # Create the command line for coex_cluster2.r using the hierarchial clustering method
  cluster2CommandH = 'coex_cluster2 -i ' + args['filterOutput'] + ' -n s -c hclust'
  if (args['minRsq']          is not None): cluster2CommandH += (' -r ' + args['minRsq'])
  if (args['maxMedianK']      is not None): cluster2CommandH += (' -k ' + args['maxMedianK'])
  if (args['maxPower']        is not None): cluster2CommandH += (' -p ' + args['maxPower'])
  if (args['minModuleSize']   is not None): cluster2CommandH += (' -s ' + args['minModuleSize'])
  if (args['detectCutHeight'] is not None): cluster2CommandH += (' -d ' + args['detectCutHeight'])
  if (args['hclustOutput']    is not None): cluster2CommandH += (' -o ' + args['hclustOutput'])

if (not('--skip_wgcna' in others)):
  # Create the command line for coex_cluster2.r using the WGCNA method
  cluster2CommandW = 'coex_cluster2 -i ' + args['filterOutput'] + ' -n w -c w'
  if (args['minRsq']          is not None): cluster2CommandW += (' -r ' + args['minRsq'])
  if (args['maxMedianK']      is not None): cluster2CommandW += (' -k ' + args['maxMedianK'])
  if (args['maxPower']        is not None): cluster2CommandW += (' -p ' + args['maxPower'])
  if (args['minModuleSize']   is not None): cluster2CommandW += (' -s ' + args['minModuleSize'])
  if (args['detectCutHeight'] is not None): cluster2CommandW += (' -d ' + args['detectCutHeight'])
  if (args['wgcnaOutput']     is not None): cluster2CommandW += (' -o ' + args['wgcnaOutput'])



# RUN THE DESIRED R FUNCTIONS FROM THE SHELL #

separator = '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
print(separator)

# coex_filter.r
if (not('--skip_filter' in others)):
  print(filterCommand)
  print(separator)
  subprocess.call(filterCommand, shell = True)
  print(separator)

# coex_net.r
if (not('--skip_network' in others)):
  print(netCommand)
  print(separator)
  subprocess.call(netCommand, shell = True)
  print(separator)

# coex_cluster.r OR coex_cluster2.r
# determine whether to run coex_cluster.r or coex_cluster2.r based on if the '--network_output' was 'adjmat' or 'a'
if (args['networkOutput'] == 'adjmat' or args['networkOutput'] == 'a'):
  if (not('--skip_cluster' in others)):
    # run coex_cluster.r
    print(clusterCommand)
    print(separator)
    subprocess.call(clusterCommand, shell = True)
    print(separator)
else:
  if (not('--skip_cluster2_hclust' in others)):
    # run coex_cluster2.r for the hclust method
    print(cluster2CommandH)
    print(separator)
    subprocess.call(cluster2CommandH, shell = True)
    print(separator)
  if (not('--skip_cluster2_wgcna' in others)):
    print(cluster2CommandW)
    print(separator)
    subprocess.call(cluster2CommandW, shell = True)
    print(separator)

print("Done!")
