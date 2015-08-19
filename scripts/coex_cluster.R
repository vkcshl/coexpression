#!/usr/bin/env Rscript

coex_cluster = function(adjmat, method = 'WGCNA', outFileName = "", minModuleSize = 30, detectCutHeight = 0.99) {
  if (is.na(method)) { method = 'WGCNA' }
  if (is.na(minModuleSize)) { minModuleSize = 30 }
  if (is.na(detectCutHeight)) { detectCutHeight = 0.99 }
  adjmat = as.matrix(adjmat)
  if (dim(adjmat)[1] != dim(adjmat)[2]) { stop("Adjacent matrix should be symmetric") }
  if (method == 'hclust' || method == 'h') {
    dissimilarity = 1 - adjmat
    distance = as.dist(dissimilarity)
    hclusters = hclust(distance)
    modulenames = cutree(hclusters, k = minModuleSize) #modulenames=cutree(hclusters,h=detectCutHeight)
  } else if (method == 'WGCNA' || method == 'w') {
    suppressPackageStartupMessages(library('WGCNA', quiet = TRUE))
    allowWGCNAThreads()
    TOM = TOMsimilarity(adjmat, verbose = 0)
    dissTOM = 1 - TOM
    geneTree = flashClust(as.dist(dissTOM), method = "average")
    dynamicMods = cutreeDynamic(dendro = geneTree, cutHeight = detectCutHeight, deepSplit = TRUE, minClusterSize = minModuleSize, method = 'tree')
    modulenames = labels2colors(dynamicMods)
    #par(pty="m");plotDendroAndColors(geneTree, dynamicColors,rowText=dynamicColors, dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
#     nSamples=nrow(datExpr)
#     MEs=moduleEigengenes(datExpr,dynamicColors)
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
  clusterInfo = data.frame(Gene = colnames(adjmat), module = modulenames)
  if (is.na(outFileName)) {
    outFileName = paste('coex_cluster_method=', method, '.csv', sep="")
  }
  write.csv(clusterInfo, file = outFileName, row.names = FALSE)  
}

# argument parsing
suppressPackageStartupMessages(library("optparse"))
option_list = list( 
  make_option(c("-i", "--input"), dest = "file_name", type = "character", 
              help = "REQUIRED: Input file name that stores the adjacency matrix of constructed gene co-expression network for clustering. The row and column names of adjacency matrix correspond gene names."), 
  make_option(c("-m", "--method"), type = "character", default = 'WGCNA', 
              help = "Method to cluster genes into co-expression modules. When method = ‘hclust’ or ‘h’,the function uses hierarchical clustering. When method = “WGCNA” or ‘W’, the function uses WGCNA. [default \"%default\"] If users want to use “hclust”, the original adjacent matrix without any cutcoff (e.g., -c -1.2) is recommended."), 
  make_option(c("-s", "--minModuleSize"), type = "double", default = 50,
              help = "Minimum size of modules when method is WGCNA. [default %default]"), 
  make_option(c("-d", "--detectCutHeight"), type = "double", default = 0.99,
              help = "Maximum heights to join modules in clustering. [default %default]"), 
  make_option(c("-o", "--output"), dest = "outFileName", type = "character", default = "coex_modules.csv",
              help = "Output file name that stores the clustering results. The first column includes gene names, and the second column includes their modules. [default \"%default\"]") 
) 

# pass arguments
description_text = "\nDESCRIPTION
\n\tGiven the adjacency matrix of gene co-expression network, this function clusters genes into co-expression modules 
using original hierarchical clustering or weighted gene co-expression network analysis (WGCNA, http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/Rpackages/WGCNA). 
\nInput: adjacency matrix of gene co-expression network in .csv file
\nOutput: a table of genes and their co-expression modules in .csv file"
usage_text="\nNAME
\n\tcoex_cluster -- cluster genes inot co-expression modules from adjacency matrix of gene co-expressio network
\nSYNOPSIS
\n\t%prog [-imsdo]"
epilogue_text = "\nEXAMPLES 
\n\t$coex_cluster  -i adjmat.csv -m WGCNA -s 30 -d 0.999 –o coex_modules.csv 
\nSEE ALSO 
\n\tcoex_filter 
\n\tcoex_net 
\n\tcoex_cluster2
\nAUTHORS 
\n\tDaifeng Wang, Gang Fang, Mark Gerstein, Eric Pan, Yale University\n"
opt_obj = OptionParser(usage=usage_text,option_list=option_list,description=description_text,epilogue=epilogue_text)
opt=parse_args(opt_obj, print_help_and_exit = FALSE)
file_name = opt$file_name
method = opt$method
outFileName = opt$outFileName
minModuleSize = opt$minModuleSize
detectCutHeight = opt$detectCutHeight
if (opt$help) {
  print_help(opt_obj)
  quit(status = 0)
}

# prepare adjmat;
options(stringsAsFactors = FALSE)
if (is.null(file_name)) { stop("please give your input file name $./coex_cluster.r -i [your input file name]") }
adjmat = as.data.frame(read.csv(file_name, header = TRUE, row.names = 1, stringsAsFactors = FALSE))

# Running the function
coex_cluster(adjmat, method = method, outFileName = outFileName, minModuleSize = minModuleSize, detectCutHeight = detectCutHeight)
