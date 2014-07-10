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
    suppressPackageStartupMessages(library('igraph', quiet = TRUE)); 
    suppressPackageStartupMessages(library('reshape', quiet = TRUE));
    options(stringsAsFactors=FALSE)
    curr_graph = graph.adjacency(adjmat_cutoff, mode = 'undirected', weighted = T, diag = F)
    curr_edgelist = get.edgelist(curr_graph, names = TRUE)
    x = adjmat_cutoff; x[upper.tri(x)] = 0; curr_edgelist = melt(x); curr_edgelist = curr_edgelist[curr_edgelist$value > 0,]
    colnames(curr_edgelist)=c('Gene1','Gene2','Weight')
    
    
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

# argument parsing
suppressPackageStartupMessages(library("optparse"))
option_list=list( 
  make_option(c("-i", "--input"), dest = "file_name", type = "character", 
              help = "REQUIRED: Input file that stores the gene expression data matrix that users want to construct gene co-expression network and do clustering. Each row corresponds to a gene, and each column corresponds to a replicate. The column names are replicate names. The row names are gene names."), 
  make_option(c("-a", "--genes1"), type = "character", default = NULL,
              help = "The first set of the genes of interest.  Leaving this out means all of the genes will be used."), 
  make_option(c("-b", "--genes2"), type = "character", default = NULL,
              help = "The second set of the genes of interest.  Leaving this out means all of the genes will be used."), 
  make_option(c("-m", "--method"), type = "character", default = 'simple', 
              help = "Method to construct gene co-expression network. When method = ‘simple’ or ‘s’, the function constructs gene co-expression network using Pearson correlation matrix whose elements less than corr_thld are set to be zero. When method = “WGCNA” or ‘W’, the function constructs gene co-expression network using Pearson correlations between genes, and “signed” network by WGCNA. [default \"%default\"]"), 
  make_option(c("-c", "--corr_thld"), type = "double", default = 0.75,
              help = "Maximum threshold to set elements of Pearson correlation matrix when method=’simple’. [default %default]"), 
  make_option(c("-q", "--p_thld"), type = "double", default = NULL,
              help = "Maximum threshold for p-value from Pearson correlation when method=’simple’. [default %default]"), 
  make_option(c("-r", "--minRsq"), type = "double", default = 0.8,
              help = "Minimum threshold for R2 that measures the fitness of gene co-expression network to scale-free topology in WGCNA. See pickSoftThreshold() of WGCNA for details. [default %default]"), 
  make_option(c("-k", "--maxmediank"), type = "double", default = 40,
              help = "Maximum median connections for genes in network. See pickSoftThreshold() of WGCNA for details. [default %default]"), 
  make_option(c("-p", "--maxpower"), type = "double", default = 50,
              help = "Maximum power to decide the soft threshold. See pickSoftThreshold() of WGCNA for details. [default %default]"), 
  make_option(c("-t", "--type"), dest = "output_type", type = "character", default = 'edge', 
              help = "Type of output. When method = ‘edge’ or ‘e’, the function outputs edge list of gene co-expression network consisting of gene pairs with Pearson correlation greater than corr_thld. When method = “adjmat” or ‘a’, the function outputs the adjacency matrix of the constructed gene co-expression network. [default \"%default\"]"), 
  make_option(c("-o", "--output"), dest="outFileName", type = "character", default="coexpression_network.csv",
              help = "Output file. [default \"%default\"]") 
) 

# pass arguments
description_text = "\nDESCRIPTION
\n\tFor a data matrix that represents gene expression across samples, this function constructs the gene co-expression network
using correlation matrix or weighted gene co-expression network analysis (WGCNA, http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/Rpackages/WGCNA). 
\nInput: a table of gene expression data matrix in .csv file
\nOutput: edge list or adjacency matrix of gene co-expression network in .csv file"
usage_text = "\nNAME
\n\tcoex_net -- To construct gene co-expression network
\nSYNOPSIS
\n\t%prog [-imcrkpto]"
epilogue_text = "\nEXAMPLES 
\n\t$coex_net  -i data.csv -c 0.8 –o edge_list.csv 
\n\t$coex_net  -i data.csv -m WGCNA -p 100 -r 0.7 -t adjmat –o adjmat.csv
\nSEE ALSO 
\n\tcoex_filter 
\n\tcoex_cluster 
\n\tcoex_cluster2
\nAUTHORS 
\n\tDaifeng Wang, Gang Fang, Mark Gerstein, Eric Pan, Yale University\n"
opt_obj = OptionParser(usage=usage_text,option_list=option_list,description=description_text,epilogue=epilogue_text)
opt = parse_args(opt_obj, print_help_and_exit = FALSE)
file_name = opt$file_name
method = opt$method
outFileName = opt$outFileName
output_type = opt$output_type
corr_thld = opt$corr_thld
p_thld = opt$p_thld
minRsq = opt$minRsq
maxmediank = opt$maxmediank
maxpower = opt$maxpower
if (opt$help) {
  print_help(opt_obj)
  quit(status = 0)
}

# prepare data
options(stringsAsFactors = FALSE)
if (is.null(file_name)) { stop("please give your input file name $./coex_net.r -i [your input file name]") }
data = as.data.frame(read.csv(file_name,header=T,row.names=1,stringsAsFactors=FALSE))
g1file = opt$genes1
g2file = opt$genes2
if (!is.null(g1file) && !is.null(g2file)) {
  g1 = read.csv(g1file)
  g2 = read.csv(g2file)
  genes1 = unlist(read.csv(g1file))
  genes2 = unlist(read.csv(g2file))
  data = filter_data(data, genes1, genes2)
} else if (!is.null(g1file)) {
  genes1 = unlist(read.csv(g1file))
} else if (!is.null(g2file)) {
  genes2 = unlist(read.csv(g2file))
} else {
  genes1 = NULL
  genes2 = NULL
}

# Running the function
res = coex_net(data, geneList1 = genes1, geneList2 = genes2, outFileName = outFileName, output_type = output_type, method = method, corr_thld = corr_thld, p_thld = p_thld, minRsq = minRsq, maxmediank = maxmediank, maxpower = maxpower)

# Save the results accordingly
if (output_type == 'e' || output_type == 'edge')
{
  write.csv(res, file = outFileName, row.names = FALSE)
} else if (output_type == 'adjmat' || output_type == 'a') {
  write.csv(res, file = outFileName, row.names = TRUE)
}
