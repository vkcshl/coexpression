#!/usr/bin/env Rscript

invR = function(p, df) {
  t = qt(p/2, df, lower.tail = FALSE)
  r = sqrt(t^2 / (df + t^2))
  return (r)
}

coex_net = function(data, geneList1 = NULL, geneList2 = NULL, method = 'simple', output_type = 'edge', outFileName = "", corr_thld = NA, p_thld = NULL, minRsq = 0.8, maxmediank = 40, maxpower = 50) {
  if (is.na(method)) { method = 'simple' }
  if (is.na(output_type)) { method = 'edge' }
  if (is.na(corr_thld) && is.null(p_thld)) { corr_thld = 0.75 }
  if (!(is.null(p_thld) || is.na(p_thld))) {
    df = dim(data)[2] - 2
    cmin = invR(p_thld, df)
    if (!is.na(corr_thld)) {
      corr_thld = max(corr_thld, cmin)
    } else {
      corr_thld = cmin
    }
  }
  if (is.na(minRsq)) { minRsq = 0.8 }
  if (is.na(maxmediank)) { maxmediank = 40 }
  if (is.na(maxpower)) { maxpower = 50 }
  if (method == 'simple' || method == 's') {
    adjmat=cor(t(data))
  } else if (method == 'WGCNA' || method == 'w') {
    datExpr = t(data)
    suppressPackageStartupMessages(library('WGCNA', quiet = TRUE)); options(stringsAsFactors=FALSE)
    allowWGCNAThreads()
    powers = c(1:maxpower)
    sft = pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)  
    sft_table = sft$fitIndices
    select_cond = sft_table$SFT.R.sq > minRsq & sft_table$median.k <= maxmediank
    if (sum(select_cond) == 0) { stop("No satisfied power found. Please decrease minRsq or increase maxpower.") }
    softPower = min(sft_table[select_cond, 'Power'])
    adjmat = adjacency(datExpr, power = softPower, type = 'signed', corFnc = "cor", corOptions = "use = 'p',method = 'pearson'")
  } else {
    stop("Please indicate a correct method. See help")
  }
  
  if (output_type == 'edge' || output_type == 'e') {
    if (is.na(corr_thld)) {
      adjmat_cutoff = as.matrix(adjmat * ((abs(adjmat) > corr_thld) * 1))
    } else {
      adjmat_cutoff = as.matrix(adjmat * ((adjmat > corr_thld) * 1))
    }
    suppressPackageStartupMessages(library('igraph', quiet = TRUE)); options(stringsAsFactors=FALSE)
    curr_graph = graph.adjacency(adjmat_cutoff, mode = 'undirected', weighted = T, diag = F)
    curr_edgelist = get.edgelist(curr_graph, names = TRUE)
    
    # (post-filter) METHOD 2 to filter genes of interest is removing edge pairs that are not part of the genes of interest
    if (!is.null(geneList1) || !is.null(geneList2)) {
      if (is.null(geneList1)) { geneList1 = rownames(data) }
      if (is.null(geneList2)) { geneList2 = rownames(data) }
      row1 = curr_edgelist[,1] %in% geneList1
      row2 = curr_edgelist[,2] %in% geneList2
      goodRows = row1 & row2
      curr_edgelist = curr_edgelist[goodRows,]
    }
    
    if(is.null(outFileName)) {
      outFileName = paste('co-expression_network_edge_list_method=',method,'.csv',sep="")
    }
    #write.csv(curr_edgelist,file=outFileName,row.names=F)
    return (curr_edgelist)
    
  }else if(output_type=='adjmat'||output_type=='a'){
  
    # (pre-filter) METHOD 1 to filter genes of interest is shaving down the adjacency matrix, which doesn't work with the 'igraph' methods
    if (!is.null(geneList1) || !is.null(geneList2))
    {
      if (is.null(geneList1)) { geneList1 = rownames(data) }
      if (is.null(geneList2)) { geneList2 = rownames(data) }
      goodRows = rownames(adjmat) %in% geneList1
      goodCols = colnames(adjmat) %in% geneList2
      adjmat = adjmat[goodRows, goodCols]
    }
    
    if (is.null(outFileName)) {
      outFileName = paste('co-expression_network_adjacency_matrix_method=', method, '.csv', sep = "")
    }
    #write.csv(adjmat,file=outFileName,row.names=T)
    return (adjmat)
  } else {
    stop("Please indicate a correct type of output. See help")
  }
}

filter_data = function(data, genes1, genes2) {
  genes = union(genes1, genes2)
  gNum1 = length(genes)
  rnames = rownames(data)
  bools = genes %in% rnames
  genes = genes[bools]
  gNum2 = length(genes);
  deleted = gNum1 - gNum2
  print(paste('NOTE: ', deleted, 'genes were not found in the data matrix and deleted.'))
  data = data[genes,]
  return (data)
}

coex_cluster = function(adjmat, method = 'WGCNA', outFileName = "", minModuleSize = 30, detectCutHeight = 0.99,data) {
  if (is.na(method)) { method = 'WGCNA' }
  if (is.na(minModuleSize)) { minModuleSize = 30 }
  if (is.na(detectCutHeight)) { detectCutHeight = 0.99 }
  adjmat = as.matrix(adjmat)
  if (dim(adjmat)[1] != dim(adjmat)[2]) { stop("Adjacent matrix should be symmetric") }
  if (method == 'hclust' || method == 'h') {
    suppressPackageStartupMessages(library('WGCNA', quiet = TRUE))
    suppressPackageStartupMessages(library('ape', quiet = TRUE))
    suppressPackageStartupMessages(library('reshape', quiet = TRUE))
    #allowWGCNAThreads()
    dissimilarity = 1 - adjmat
    distance = as.dist(dissimilarity)
    hclusters = hclust(distance)
    modulenames = cutree(hclusters, k = minModuleSize) #modulenames=cutree(hclusters,h=detectCutHeight)
    phy=as.phylo(hclusters)
    write.tree(phy,file=paste('hclust_tree_method=', method, '.txt', sep=""))
  } else if (method == 'WGCNA' || method == 'w') {
    suppressPackageStartupMessages(library('WGCNA', quiet = TRUE))
    suppressPackageStartupMessages(library('reshape', quiet = TRUE))
    #allowWGCNAThreads()
    TOM = TOMsimilarity(adjmat, verbose = 0)
    dissTOM = 1 - TOM
    geneTree = flashClust(as.dist(dissTOM), method = "average")
    dynamicMods = cutreeDynamic(dendro = geneTree, cutHeight = detectCutHeight, deepSplit = TRUE, minClusterSize = minModuleSize, method = 'tree')
    modulenames = labels2colors(dynamicMods)
    #par(pty="m");plotDendroAndColors(geneTree, dynamicColors,rowText=dynamicColors, dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
#     nSamples=nrow(datExpr)
     
#     modNames=substring(names(MEs$eigengenes),3)
#     geneModuleMembership=as.data.frame(cor(datExpr,MEs$eigengenes,use="p"))
#     MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
#     # Define column names for geneModuleMembership.
#     names(geneModuleMembership) = paste("MM", modNames, sep="")
#     # Define column names for MMPvalue.
#     names(MMPvalue) = paste("p.MM", modNames, sep="")
  } else {
    stop("Please indicate a correct method. See help")
  }
  print(table(modulenames))
  MEs=moduleEigengenes(t(data),modulenames)
  eigen_adjmat=cor(MEs$eigengenes,use='p')
  x = eigen_adjmat; x[upper.tri(eigen_adjmat)] = 0; mod_edgelist = melt(x); mod_edgelist = mod_edgelist[mod_edgelist$value!=0&mod_edgelist[,1]!=mod_edgelist[,2],]
  colnames(mod_edgelist)=c('node 1','node 2','weight')
  write.csv(mod_edgelist, file = paste('module_network_edgelist_method=', method, '.csv', sep=""), row.names = F) 
  clusterInfo = data.frame(Gene = colnames(adjmat), module = modulenames)
  if (is.na(outFileName)) {
    outFileName = paste('coex_cluster_method=', method, '.csv', sep="")
  }
  write.csv(clusterInfo, file = outFileName, row.names = FALSE)  
}

# argument parsing
suppressPackageStartupMessages(library("optparse"))
option_list=list( 
  make_option(c("-i", "--input"),dest="file_name",type="character", 
              help="REQUIRED: Input file is the gene expression data matrix in .csv format. Each row corresponds to a gene, and each column corresponds to a sample. The first column are gene names. The first row are gene names."), 
  make_option(c("-a", "--genes1"), type="character", default=NULL,
              help="The first set of the genes of interest.  Leaving this out means all of the genes will be used."), 
  make_option(c("-b", "--genes2"), type="character", default=NULL,
              help="The second set of the genes of interest.  Leaving this out means all of the genes will be used."), 
  make_option(c("-n", "--net_method"),type="character",default='simple', 
              help="Method to construct gene co-expression network. When net_method = ‘simple’ or ‘s’, the function constructs gene co-expression network using Pearson correlation matrix. When net_method = “WGCNA” or ‘w’, the function constructs gene co-expression network using Pearson correlations coefficient between genes, and “signed” network by WGCNA. [default \"%default\"]"), 
  make_option(c("-r", "--minRsq"), type="double", default=0.8,
              help="Minimum threshold for R2 that measures the fitness of gene co-expression network to scale-free topology in WGCNA. See pickSoftThreshold() of WGCNA for details. [default %default]"), 
  make_option(c("-k", "--maxmediank"), type="double", default=40,
              help="Maximum median connections for genes in network. See pickSoftThreshold() of WGCNA for details. [default %default]"), 
  make_option(c("-p", "--maxpower"), type="double", default=50,
              help="Maximum power to decide the soft threshold. See pickSoftThreshold() of WGCNA for details. [default %default]"), 
  make_option(c("-c", "--clust_method"),type="character",default='WGCNA', 
              help="Method to cluster genes into co-expression modules. When clust_method = ‘hclust’ or ‘h’,the function uses hierarchical clustering. When clust_method = “WGCNA” or ‘w’, the function uses WGCNA. [default \"%default\"]"), 
  make_option(c("-s", "--minModuleSize"), type="double", default=50,
              help="Minimum size of modules when the method is WGCNA. [default %default]"), 
  make_option(c("-d", "--detectCutHeight"), type="double", default=0.99,
              help="Maximum heights to join modules in clustering. [default %default]"), 
  make_option(c("-o", "--output"),dest="outFileName",type="character",default="coex_modules.csv",
              help="Output file that stores the clustering results. The first column is gene, and the second column is the corresponding module. [default \"%default\"]") 
) 

# pass arguments
description_text = "\nDESCRIPTION
\n\tGiven the gene expression table, this function clusters genes into co-expression modules 
using original hierarchical clustering or weighted gene co-expression network analysis (WGCNA, http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/Rpackages/WGCNA). 
\nInput: a table of gene expression data matrix in .csv file
\nOutput: a table of genes and their co-expression modules in .csv file"
usage_text = "\nNAME
\n\tcoex_cluster2 -- To cluster genes into co-expression modules from gene expression data matrix
\nSYNOPSIS
\n\t%prog [-inrkpcsdo]"
epilogue_text = "\nEXAMPLES 
\n\t$coex_cluster2  -i data.csv -c WGCNA -n WGCNA –o coex_modules.csv
\n\t$coex_cluster2  -i data.csv -c hclust -n simple –o coex_modules.csv
\nSEE ALSO 
\n\tcoex_filter 
\n\tcoex_net 
\n\tcoex_cluster
\nAUTHORS 
\n\tDaifeng Wang, Gang Fang, Mark Gerstein, Eric Pan, Yale University\n"
opt_obj = OptionParser(usage=usage_text,option_list=option_list,description=description_text,epilogue=epilogue_text)
opt = parse_args(opt_obj, print_help_and_exit = FALSE)
file_name = opt$file_name
net_method = opt$net_method
outFileName = opt$outFileName
minModuleSize = opt$minModuleSize
detectCutHeight = opt$detectCutHeight
clust_method = opt$clust_method
outFileName = opt$outFileName
minRsq = opt$minRsq
maxmediank = opt$maxmediank
maxpower = opt$maxpower
if (opt$help) {
  print_help(opt_obj)
  quit(status = 0)
}

# prepare data
options(stringsAsFactors=FALSE)
if (is.null(file_name)) {stop("please give your input file name $./coex_net.r -i [your input file name]") }
data = as.data.frame(read.csv(file_name, header = T, row.names = 1, stringsAsFactors = FALSE))
g1file = opt$genes1
g2file = opt$genes2
if (!is.null(g1file) && !is.null(g2file))
{
  g1 = read.csv(g1file)
  g2 = read.csv(g2file)
  genes = unique(rbind(g1, g2)[,1])
  gNum1 = length(genes)
  rnames = rownames(data)
  bools = genes %in% rnames
  genes = genes[bools]
  gNum2 = length(genes);
  deleted = gNum1 - gNum2;
  print(paste('NOTE: ', deleted, 'genes were not found in the data matrix and deleted.'))
  data = data[genes,]
  g1 = g1[,1]
  g2 = g2[,1]
} else {
  g1 = NULL
  g2 = NULL
}

# Running the function

#coex_cluster2(data=data,outFileName=outFileName,clust_method=clust_method,net_method=net_method,minRsq=minRsq,maxmediank=maxmediank,maxpower=maxpower,minModuleSize=minModuleSize,detectCutHeight=detectCutHeight)

adjmat = coex_net(data, geneList1 = g1, geneList2 = g2, output_type = 'adjmat', method = net_method, minRsq = minRsq, maxmediank = maxmediank, maxpower = maxpower)
coex_cluster(adjmat, method = clust_method, outFileName = outFileName, minModuleSize = minModuleSize, detectCutHeight = detectCutHeight,data=data)
