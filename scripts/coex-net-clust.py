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
      coex-net-clust -- Construct coexpression network and identify clusters

SYNOPSIS

'''

desc2 = '''
DESCRIPTION

  This command provides the function to build coexpression network and identify the functional modules among it.
  A functional module is a network cluster with enrichment of certain biological function. const_coex_net_clust first construct coexpression network. Then, it identifys the clusters among the network. Finally, it identifys the GeneOntology enrichment for the genes in each cluster.


'''

desc3 = '''
EXAMPLES
      Using PCC to build coexpression network. Using hierarchical clustering to identify the clusters among the network. Enrichemnt test of gene ontology terms for each cluster will be automatically calculated.
      > coex-net-clust --ws_url 'https://kbase.us/services/ws' --ws_id 'plane83:Fei_coex_demo'  --in_id 'filtered_series' --out_id 'my_net'  --cut_off 0.8 --net_method simple --clust_method hclust  --num_module 15 
      > coex-net-clust -u 'https://kbase.us/services/ws' -w 'plane83:Fei_coex_demo'  -i 'filtered_series' --out_id 'my_net'  -c 0.8 -n simple -m hclust  -k 15 
     
      Using WGCNA to build coexpression network and to identify the clusters among the network. Enrichemnt test of gene ontology terms for each cluster will be automatically calculated.
      > coex-net-clust -u 'https://kbase.us/services/ws' -w 'plane83:Fei_coex_demo'  -i 'filtered_series' --out_id 'my_net'  -c 0.8 -n WGCNA -m WGCNA  -k 15 

SEE ALSO
      coex_filter

AUTHORS
Shinjae Yoo, Gang Fang, Fei He, Daifeng Wang.
      
'''
class Node:
    nodes = []
    edges = []
    ugids = {}
    igids = {}
    gid2nt = {}
    clst2genes = {}

    def __init__(self, unodes = [], uedges=[]):
      self._register_nodes(unodes)
      self._register_edges(uedges)
  
    def get_node_id(self, node, nt = "GENE"):
      if not node in self.ugids.keys() :
          #print node + ":" + nt
          self.ugids[node] = len(self.ugids)
          self.nodes.append( {
            'entity_id' : node,
            'name' : node,
            'user_annotations' : {},
            'type' : nt,
            'id' : 'kb|netnode.' + `self.ugids[node]`,
            'properties' : {}
          } )
          self.igids['kb|netnode.' + `self.ugids[node]`] = node
          self.gid2nt[node] = nt
      return "kb|netnode." + `self.ugids[node]`

    def add_edge(self, strength, ds_id, node1, nt1, node2, nt2, confidence):
      #print node1 + "<->" + node2
      self.edges.append( {
          'name' : 'interacting gene pair',
          'properties' : {},
          'strength' : float(strength),
          'dataset_id' : ds_id,
          'directed' : 'false',
          'user_annotations' : {},
          'id' : 'kb|netedge.'+`len(self.edges)`,
          'node_id1' : self.get_node_id(node1, nt1),
          'node_id2' : self.get_node_id(node2, nt2),
          'confidence' : float(confidence)
      })
      if(nt1 == 'CLUSTER'):
        if not node1 in self.clstr2genes.keys() : self.clst2genes[node1] = {}
        if(nt2 == 'GENE'):
          self.clst2gene[node1][node2] = 1
      else:
        if(nt2 == 'CLUSTER'):
          if not node2 in self.clst2genes.keys() : self.clst2genes[node2] = {}
          self.clst2genes[node2][node1] = 1
   
    def _register_nodes(self, unodes):
      self.nodes = unodes
      self.ugids = {}
      for node in self.nodes:
        nnid = node['id']
        nnid = nnid.replace("kb|netnode.","");
        self.ugids[node['entity_id']] = nnid
        self.igids[node['id']] = node['entity_id']
        self.gid2nt[node['entity_id']] = node['type']

    def _register_edges(self, uedges):
      self.edges = uedges
      for edge in self.edges:
        node1 = self.igids[edge['node_id1']];
        nt1  = self.gid2nt[node1];
        node2 = self.igids[edge['node_id2']];
        nt2  = self.gid2nt[node2];
        if(nt1 == 'CLUSTER'):
          if not node1 in self.clstr2genes.keys() : self.clst2genes[node1] = {}
          if(nt2 == 'GENE'):
            self.clst2genes[node1][node2] = 1
        else:
          if(nt2 == 'CLUSTER'):
            if not node2 in self.clst2genes.keys() : self.clst2genes[node2] = {}
            self.clst2genes[node2][node1] = 1
        

    def get_gene_list(self, cnode):
      if(cnode in self.clst2genes.keys()) : return self.clst2genes[cnode].keys()
      return []
     


def net_clust (args) :
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
    gids = samples[0]['data']['expression_levels'].keys()  # each sample has same gids
    for gid in sorted(gids):
        line = gid + ","
        line += ",".join([str(s['data']['expression_levels'][gid]) for s in samples])
        cif.write(line + "\n")
    cif.close()


    ###
    # generate network and cluster
    net_cmd_lst = ['coex_net', '-i', args.exp_fn]
    if (args.nmethod    is not None): 
        net_cmd_lst.append("-m")
        net_cmd_lst.append(args.nmethod)
    if (args.cut_off    is not None): 
        net_cmd_lst.append("-c")
        net_cmd_lst.append(args.cut_off)
    if (args.net_fn     is not None):
        net_cmd_lst.append("-o")
        net_cmd_lst.append(args.net_fn)
    p1 = Popen(net_cmd_lst, stdout=PIPE)
    out_str = p1.communicate()
    if out_str[0] is not None : print out_str[0]
    if out_str[1] is not None : print >> sys.stderr, out_str[1]
    net_cmd = " ".join(net_cmd_lst)
   
   
    clust_cmd_lst = ['coex_cluster2', '-i', args.exp_fn]
    if (args.cmethod    is not None):
        clust_cmd_lst.append("-c")
        clust_cmd_lst.append(args.cmethod)
    if (args.nmethod    is not None):
        clust_cmd_lst.append("-n")
        clust_cmd_lst.append(args.nmethod)
    if (args.k          is not None):
        clust_cmd_lst.append("-s")
        clust_cmd_lst.append(args.k)
    if (args.clust_fn   is not None):
        clust_cmd_lst.append("-o")
        clust_cmd_lst.append(args.clust_fn)
    p1 = Popen(clust_cmd_lst, stdout=PIPE)
    out_str = p1.communicate()
    if out_str[0] is not None : print out_str[0]
    if out_str[1] is not None : print >> sys.stderr, out_str[1]
    clust_cmd = " ".join(clust_cmd_lst)

   
    ###
    # Create network object
    #generate Networks datasets
    net_ds_id = args.inobj_id + ".net"
    clt_ds_id = args.inobj_id + ".clt"
 
    datasets = [
      {
        'network_type' : 'FUNCTIONAL_ASSOCIATION',
        'taxons' : [ genome_id ],
        'source_ref' : 'WORKSPACE',
        'name' : net_ds_id,
        'id' : clt_ds_id,
        'description' : "Coexpression network object of " + args.inobj_id,
        'properties' : {
          'original_data_type' : 'workspace',
          'original_ws_id' : args.ws_id,
          'original_obj_id' : args.inobj_id,
          'coex_net_cmd' : net_cmd
        }
      },
      {
        'network_type' : 'FUNCTIONAL_ASSOCIATION',
        'taxons' : [ genome_id ],
        'source_ref' : 'WORKSPACE',
        'name' : clt_ds_id,
        'id' : clt_ds_id,
        'description' : "Coexpression cluster object of " + args.inobj_id,
        'properties' : {
          'original_data_type' : 'workspace',
          'original_ws_id' : args.ws_id,
          'original_obj_id' : args.inobj_id,
          'coex_clust_cmd' : clust_cmd
        }
      }
    ]
 
 
    # process coex network file
    nc = Node()
 
    cnf = open(args.net_fn,'r');
    cnf.readline(); # skip header
    for line in cnf :
        line.strip();
        line = line.replace('"','')
        values = line.split(',')
        if values[0] != values[1] : nc.add_edge(float(values[2]), net_ds_id, values[0], 'GENE', values[1], 'GENE', 0.0) #we add edges meaningful
 
 
    # process coex cluster file
    cnf = open(args.clust_fn,'r')
    cnf.readline(); # skip header
    for line in cnf :
        line = line.strip();
        line = line.replace('"','')
        values = line.split(',')
        nc.add_edge(1.0, clt_ds_id, values[0], 'GENE', "cluster." + values[1], 'CLUSTER', 0.0)
 
    # generate Networks object
    net_object = {
      'datasets' : datasets,
      'nodes' : nc.nodes,
      'edges' : nc.edges,
      'user_annotations' : {},
      'name' : 'Coexpression Network',
      'id' : args.outobj_id,
      'properties' : {
        'graphType' : 'edu.uci.ics.jung.graph.SparseMultigraph'
      }
    }
 
    # Store results object into workspace
    wsd.save_objects({'workspace' : args.ws_id, 'objects' : [{'type' : 'KBaseNetworks.Network', 'data' : net_object, 'name' : args.outobj_id, 'meta' : {'org_obj_id' : args.inobj_id, 'org_ws_id' : args.ws_id}}]})
 
    if(args.del_tmps is "true") :
        os.remove(args.exp_fn)
        os.remove(args.net_fn)
        os.remove(args.clust_fn)
 

if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='coex-net-clust', epilog=desc3)
    parser.add_argument('-u', '--ws_url', help='Workspace url', action='store', dest='ws_url', default='https://kbase.us/services/ws')
    parser.add_argument('-w', '--ws_id', help='Workspace id', action='store', dest='ws_id', default=None, required=True)
    parser.add_argument('-i', '--in_id', help='Input Series object id', action='store', dest='inobj_id', default=None, required=True)
    parser.add_argument('-o', '--out_id', help='Output network object id', action='store', dest='outobj_id', default=None, required=True)
    parser.add_argument('-m', '--clust_method', help='Clustering method (\'hclust\' for hierachical clustering or \'WGCNA\' for WGCNA method', action='store', dest='cmethod', default='hclust')
    parser.add_argument('-n', '--net_method', help='Network construction method (\'simple\' for PCC or \'WGCNA\' for WGCNA method', action='store', dest='nmethod', default='simple')
    parser.add_argument('-c', '--cut_off', help='The edge cut-off value', action='store', dest='cut_off', default=None)
    parser.add_argument('-k', '--num_module', help='The number of module to be generated', action='store', dest='k', default=None, required=True)
    parser.add_argument('-e', '--expression_fn', help='Expression file name (temporary file)', action='store', dest='exp_fn', default='expression.csv')
    parser.add_argument('-t', '--net_out_fn', help='Network output file name (temporary file)', action='store', dest='net_fn', default='net.csv')
    parser.add_argument('-l', '--clust_out_fn', help='Cluster output file name (temporary file)', action='store', dest='clust_fn', default='clust.csv')
    parser.add_argument('-d', '--del_tmp_files', help='Delete temporary files', action='store', dest='del_tmps', default='true')
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # main loop
    net_clust(args)
    exit(0);
