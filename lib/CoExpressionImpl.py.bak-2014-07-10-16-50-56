#BEGIN_HEADER
#END_HEADER


class CoExpression:
    '''
    Module Name:
    CoExpression

    Module Description:
    Co-Expression Service APIs 

 This module provides services for plant expression data in support of the coexpression
 network and ontology driven data needs of the plant sciences community. This version of
 the modules supports retrieval of the following information:
 1. Retrieval of GEO sample ID list for given EO (environmental ontology) and/or PO (plant ontology -plant tissues/organs of interest).
 2. Retrieval of the expression values for given GEO sample ID list.  
 3. For given expression values tables, it computes co-expression clusters or network (CLI only).

It will serve queries for tissue and condition specific gene expression co-expression network for biologically interesting genes/samples. Users can search differentially expressed genes in different tissues or in numerous experimental conditions or treatments (e.g various biotic or abiotic stresses). Currently the metadata annotation is provided for a subset of gene expression experiments from the NCBI GEO microarray experiments for Arabidopsis and Poplar. The samples of these experiments are manually annotated using plant ontology (PO) [http://www.plantontology.org/] and environment ontology (EO) [http://obo.cvs.sourceforge.net/viewvc/obo/obo/ontology/phenotype/environment/environment_ontology.obo]
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        #END_CONSTRUCTOR
        pass

    def filter_genes(self, args):
        # self.ctx is set by the wsgi application class
        # return variables are: job_id
        #BEGIN filter_genes
        #END filter_genes

        #At some point might do deeper type checking...
        if not isinstance(job_id, list):
            raise ValueError('Method filter_genes return value ' +
                             'job_id is not type list as required.')
        # return the results
        return [job_id]

    def const_coex_net_clust(self, args):
        # self.ctx is set by the wsgi application class
        # return variables are: job_id
        #BEGIN const_coex_net_clust
        #END const_coex_net_clust

        #At some point might do deeper type checking...
        if not isinstance(job_id, list):
            raise ValueError('Method const_coex_net_clust return value ' +
                             'job_id is not type list as required.')
        # return the results
        return [job_id]
