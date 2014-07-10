/*
Co-Expression Service APIs 

 This module provides services for plant expression data in support of the coexpression
 network and ontology driven data needs of the plant sciences community. This version of
 the modules supports retrieval of the following information:
 1. Retrieval of GEO sample ID list for given EO (environmental ontology) and/or PO (plant ontology -plant tissues/organs of interest).
 2. Retrieval of the expression values for given GEO sample ID list.  
 3. For given expression values tables, it computes co-expression clusters or network (CLI only).

It will serve queries for tissue and condition specific gene expression co-expression network for biologically interesting genes/samples. Users can search differentially expressed genes in different tissues or in numerous experimental conditions or treatments (e.g various biotic or abiotic stresses). Currently the metadata annotation is provided for a subset of gene expression experiments from the NCBI GEO microarray experiments for Arabidopsis and Poplar. The samples of these experiments are manually annotated using plant ontology (PO) [http://www.plantontology.org/] and environment ontology (EO) [http://obo.cvs.sourceforge.net/viewvc/obo/obo/ontology/phenotype/environment/environment_ontology.obo]

*/

module CoExpression 
{
  authentication required;

  typedef structure {
    string ws_id;
    string inobj_id; /* series object id */
    string outobj_id;
    string p_value;
    string method;
    string num_genes;
  } FilterGenesParams;

  funcdef filter_genes(FilterGenesParams args) returns (list<string> job_id);

  typedef structure {
    string ws_id;
    string inobj_id; /* series object id */
    string outobj_id;
    string cut_off;
    string net_method;
    string clust_method;
    string num_modules;
  } ConstCoexNetClustParams;

  funcdef const_coex_net_clust(ConstCoexNetClustParams args) returns (list<string> job_id);
};
