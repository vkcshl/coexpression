/*
Co-Expression Service APIs 

 This module provides services in support of the coexpression network. 
 The modules supports retrieval of the following information:
 1. Identify differentially expressed genes
 2. WGCNA clustering

*/

module CoExpression 
{
  /*You need to use KBase auth service to get the authentication*/
  authentication required;

  typedef string obj_ref;

  /*
      Designed to support for p-value distribution plot to give idea on what would be good p-value cut off but it can be used for other purpose.
      @optional xlabel ylabel xlog_mode ylog_mode title description plot_type properties 
  */
  typedef structure {
    string xlabel; /* labels for xaxis */
    string ylabel; /* labels for xaxis */
    string xlog_mode; /* x axis log mode : "none", "log2", "log10", "-log2", "-log10", "no^-log2" for switching between no log and -log2, "no^log10" for switching between "none" and "log10", etc */
    string ylog_mode; /* y axis log mode : "none", "log2", "log10", "-log2", "-log10", "no^-log2" for switching between no log and -log2, "no^log10" for switching between "none" and "log10", etc */
    string title; /* title of plot */
    string description; /* description of plot */
    string plot_type; /* histogram, scatter, bar, treemap, pie, auto, user_select etc */
    mapping<string, string> properties; /* additional properties */
    obj_ref data_ref; /* data object such as MAK.FloatDataTable */
  } FigureProperties;

  /* @id ws KBaseFeatureValues.ExpressionMatrix */
  typedef string ws_expression_matrix_id;

  typedef structure {	  
    string workspace_name; /* workspace name*/
    string inobj_id; /*inobj_id is expression series object id */
    string outobj_id; /*outobj_id is the output object id*/
    string p_value;/*p_value is the p-value of the statistical significance of differential expression*/
    string method;/*method is the method used for identification of differentially expressed genes*/
    string num_genes;/*num_gene is user for specify how many differentially expressed genes are needed*/
  } FilterGenesParams;

  typedef structure {	  
    string workspace_name; /* workspace name*/
    string out_expr_object_name; /*output expression object name*/
    string out_fs_object_name; /*output feature set object name*/
  } FilterGenesResult;

  /* Description of diff_p_distribution:
    diff_p_distribution provides the function to generate p-value distribution as png or MAK.FloatDataTable object.
  */   
  async funcdef diff_p_distribution(FilterGenesParams args) returns (FigureProperties result);
  
  /* Description of filter_genes: 
  filter_genes provides the function to identify differentially expressed genes given an expression series/experiment. An expression series/experiment contains a list of expression samples. A expression sample is the measurement of mRNA abundance in a biological sample. The design of expression profiling usually includes replicates. The replicates allows us to differ the non-relevent expression variation and the relevent expression variation.
  The replicate information is manully extracted by KBase developers. Only a part of samples has been assigned to a replicate group. For those samples without an assignment, the variation of its expression abundance is used directly.
  filter_genes now has two methods to identify differentially expressed genes: ANOVA and lor(from limma r package). The output of this function is a reduced FV expression matrix and a FV FeatureSet*/
  
  async funcdef filter_genes(FilterGenesParams args) returns (FilterGenesResult result);

  typedef structure {
    string workspace_name; /* workspace name*/
    ws_expression_matrix_id inobj_id; /* series object id */
    string outobj_id; /*outobj_id is the output object id*/
    string cut_off; /*cut_off is the statistical threshold to define a coexpression relationship*/
    string net_method; /*net_method is the method to construct coexpression network. Currently, two methods have been implemented: WGCNA(Weighted Gene Co-Expression Network) and simple PCC-based approach*/
    string clust_method; /*clust_method is the method to identify network clusters. Currently, two methods have been implemented:  WGCNA(Weighted Gene Co-Expression Network) and hclust(Hierarchical clustering)*/
    string num_modules; /*num_modules is used to define the number of modules need to be identified from the network*/
  } ConstCoexNetClustParams;

  typedef structure {
    string workspace_name; /* workspace name*/
    string out_object_name; /*output object name*/
  } ConstCoexNetClustResult;



  /*Description of const_coex_net_clust
  const_coex_net_clust provides the function to build coexpression network and identify the functional modules among it.
  A functional module is a network cluster with enrichment of certain biological function. const_coex_net_clust first construct coexpression network. Then, it identifys the clusters among the network. Finally, it identifys the GeneOntology enrichment for the genes in each cluster.
  */
  async funcdef const_coex_net_clust(ConstCoexNetClustParams args) returns (ConstCoexNetClustResult result);

  
};



