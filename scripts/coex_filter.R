#!/usr/bin/env Rscript

coex_filter = function(data,method="lor",p_threshold=1.0,topnumber=0,outFileName="",sample_index,resp='y',tsv='y', genelistFileName = 'selected.txt', plotfn = 'pv_distribution.png', jsonfn = 'pv_distribution.json'){
  if (is.na(sample_index[1])) { sample_index = rep(c(1:(dim(data)[2]/2)), each = 2) }
  if (is.na(p_threshold)) { p_threshold = 0.05 }
  if (is.na(topnumber)) { topnumber = 0 }
  if (is.na(method)) { method = 'lor' }
  if (is.na(resp)) { resp = 'y' }
  if (dim(data)[2] != length(sample_index)) { stop("the length of sample index should be equal to the number of data columns") }
  gene_p = rep(0, times = nrow(data))
  if (method == 'lor' || method == 'l') {
    suppressPackageStartupMessages(library('limma', quiet = T))
    DES = as.factor(unlist(c(sapply(unique(sample_index),function(x) { rep(paste('S', x, sep = ''), sum(sample_index == x)) }))))
    design = model.matrix(~DES)
    #y = voom(data, design)
    #linear_fit = lmFit(y, design)
    linear_fit = lmFit(data, design)
    fit = eBayes(linear_fit)
    topgenes = topTable(fit, coef = 2, n = nrow(data))
    if (resp == 'n') { gene_p = topgenes$P.Val} else { gene_p = topgenes$adj.P.Val }
    #genelist = topgenes$ID ## for R v2
    genelist = row.names(topgenes) ## for R v3
  } else if (method == 'anova' || method == 'a') {
    for (i in 1:length(gene_p)) {
      gene_anv = aov(as.numeric(data[i,]) ~ sample_index)
      gene_p[i] = summary(gene_anv)[[1]]$"Pr(>F)"[1]#gene_anv$"Pr(>F)"[1]
      if(is.na(gene_p[i]) | is.nan(gene_p[i]) | is.null(gene_p[i])) {
        gene_p[i] = 1.1
      }
      genelist=rownames(data)
    }
  }
  if (topnumber > length(gene_p)) {
    stop('Requested top genes are more than total genes.')
  } else if (sum(gene_p < p_threshold, na.rm = TRUE) == 0) {
    stop('No highly differentially expressed gene was found. Please increase p-value threshold.')
  }
  if( p_threshold < 1) { 
    selected_genes = genelist[gene_p < p_threshold]
    gene_p_selected = gene_p[gene_p < p_threshold]
    if( topnumber > 0 & topnumber <= length(selected_genes)) {
      selected_genes = selected_genes[which(order(gene_p_selected, decreasing = FALSE) <= topnumber)]
    }
  }
  else if (topnumber > 0 & topnumber <= length(gene_p)) {
    selected_genes = genelist[which(order(gene_p, decreasing = FALSE) <= topnumber)] 
  }
  else {
    stop('Please specify either p-value or number of genes.')
  }
  
  if (length(unique(sample_index)) != dim(data)[2]) {
    data_preproc = sapply(unique(sample_index),function(x) {
      if (sum(sample_index == x) > 1) {
        rowMeans(data[,which(sample_index==x)])
      } else {
        data[,x]
      }
    })
    data_filter = data_preproc[selected_genes,]
    if(length(selected_genes) == 1) {
      data_filter = t(as.data.frame((data_preproc[selected_genes,])))
      rownames(data_filter) = selected_genes
    }
  } else {
    data_filter = data[selected_genes,]
    if(length(selected_genes) == 1) {
      data_filter = t(as.data.frame((data[selected_genes,])))
      rownames(data_filter) = selected_genes
    }
  }
  
  if (is.na(outFileName)) {
    outFileName = "Gene_expression_data_preprocessed_for_clustering.csv"
  }  


#generate a plot to show the distribution of p-values. Added by Fei, Jan 21, 2016
  iw=as.matrix(gene_p)
  iw=-log10(iw)
  png(filename=plotfn)
  hist(iw,50,main='Histogram of P-values',xlab='-log10(P-value)',col='green')
  dev.off()
#end
  rownames(iw) = genelist
  colnames(iw) = c('p-value')
  #iwjs = toJSON(as.data.frame((iw)))
  datas = toJSON((t(iw)))
  column_ids = toJSON(genelist)
  js = paste('{"name" : "Histogram of P-values", "row_labels":["p-value"], "column_labels":', column_ids, ',"data":', datas, "}")
  
  write(js, jsonfn)



  if(tsv == 'y') {
    write.table(data_filter, file = outFileName, sep="\t", row.names = TRUE)
  } else {
    write.csv(data_filter, file = outFileName, row.names = TRUE)
  }
  write(selected_genes, file = genelistFileName, sep="\t")
  if (topnumber > 0) {
    cat(paste(nrow(data_filter), ' highly differentially expressed genes are selected.\n', sep = ''))
  } else {
    cat(paste(nrow(data_filter),' highly differentially expressed genes are selected with p < ', p_threshold, '.\n', sep = ''))
  }
}

# argument parsing
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("jsonlite"))
option_list = list( 
  make_option(c("-i", "--input"), dest = "file_name", type = "character", 
              help = "REQUIRED: Input file that stores original gene expression data matrix. Each row corresponds to a gene, and each column corresponds to a sample. The column names are sample names. The row names are gene names."), 
  make_option(c("-m", "--method"), type = "character", default = 'lor', 
              help = "Method to identify highly differentially expressed genes. Choices include anova and log-odd ratio (lor). [default \"%default\"]"), 
  make_option(c("-p", "--p_threshold"), type = "double", default = 1.01,
              help = "Maximum p-value up to which genes are significantly differentially expressed [default %default]"), 
  make_option(c("-n", "--topnumber"), type = "integer", default = 0,
              help = "Number of top genes that are most highly differentially expressed [default %default]"),
  make_option(c("-o", "--output"), dest = "outFileName", type = "character",
              default = "Gene_expression_data_filtered_for_clustering.csv", 
              help = "Output file that stores gene expression data matrix including selected genes only [default \"%default\"]"), 
  make_option(c("-s", "--sample_index"), type = "character", default = NA,
              help = "RECOMMENDED unless no replicates in data: The file that stores a numeric vector to indicate the sample indices that replicates correspond. For example, if sample_index is equal to c(1,1,2,2,3,3), the first two columns of data correspond to Sample 1, the third and fourth columns correspond to Sample 2, and the fifth and sixth columns correspond to Sample 3."),
  make_option(c("-u", "--human_input"), type = "character", default = 'n',
              help = "Accept human input regarding replicates? 'y' to allow, 'n' to automatically use defaults. [y/n]  (Default is 'n')"),
  make_option(c("-r", "--has_replicates"), type = "character", default = 'y',
              help = "Does your data have replicates (this helps avoid user input during program runtime)?  [y/n]  (Default is 'y')"),
  make_option(c("-d", "--default_action"),type="character",default='y',
              help = "Use the default action when dealing with replicates?  [y/n]  (Default is 'y')"),
  make_option(c("-x", "--genelistFileName"), type = "character", default = 'selected.txt', 
              help = "Output selected gene list file name. [default \"%default\"]"), 
  make_option(c("-l", "--plotpv"), type = "character", default = 'pv_distribution.png', 
              help = "Output p-value distribution file name. [default \"%default\"]"), 
  make_option(c("-j", "--jsonpv"), type = "character", default = 'pv_distribution.json', 
              help = "Output p-value list json file name. [default \"%default\"]"), 
  make_option(c("-t", "--tab_delim"),type="character",default='n',
              help = "Use tab as a deliminator?  [y/n]  (Default is 'n')")
)

# pass arguments
description_text = "\nDESCRIPTION
\n\tThis function identifies genes that are highly differentially expressed across samples using ANOVA (http://www.maths.bath.ac.uk/~jjf23/book/) 
or Log-odd ratios (http://www.bioconductor.org/packages/2.12/bioc/html/limma.html).
Users can identify the highly differentially expressed genes across samples either setting up a p-value threshold or choosing the top number.
\nInput: a table of gene expression data matrix in .csv file
\nOutput: a table of highly differentially expressed genes along with expression value, i.e., mean of replicates per sample in .csv file"
usage_text = "\nNAME
\n\tcoex_filter -- identify highly differeitially expressed genes
\nSYNOPSIS
\n\t%prog [-impnos]"
epilogue_text = "\nEXAMPLES 
\n\t$coex_filter -i data.csv –m lor -p 0.01 -s sample_id.csv -o datafilter.csv 
\n\t$coex_filter -i data.csv –m anova -n 200 -s sample_id.csv -o datafilter.csv
\nSEE ALSO 
\n\tcoex_net 
\n\tcoex_cluster 
\n\tcoex_cluster2
\nAUTHORS 
\n\tDaifeng Wang, Gang Fang, Mark Gerstein, Eric Pan, Yale University\n"
opt_obj = OptionParser(usage=usage_text,option_list=option_list,description=description_text,epilogue=epilogue_text)
opt = parse_args(opt_obj, print_help_and_exit = FALSE)
file_name = opt$file_name
method = opt$method
p_threshold = opt$p_threshold
topnumber = opt$topnumber
outFileName = opt$outFileName
sample_index = opt$sample_index
humanInput = opt$human_input
hasReplicates = opt$has_replicates
defaultAction = opt$default_action
tab_delim = opt$tab_delim
glfn = opt$genelistFileName
plotfn = opt$plotpv
jsonfn = opt$jsonpv
if (opt$help) {
  print_help(opt_obj)
  quit(status = 0)
}

# prepare data
options(stringsAsFactors = FALSE)
if (is.null(file_name)) { stop("please give your input file name $./coex_filter.r -i [your input file name]") }
if(tab_delim == 'y'){ # inefficient but minimal testing
  data = as.data.frame(read.csv(file_name, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE))
} else {
  data = as.data.frame(read.csv(file_name, header = TRUE, row.names = 1, stringsAsFactors = FALSE))
}
if (is.na(sample_index)) {
  if (humanInput == 'y') {
    cat("Does your data have replicates? [y/n]: ") 
    resp=readLines("stdin", 1, warn = FALSE)
  }
  else
  {
    resp = hasReplicates
  }
  if (resp == 'y') {
    if (humanInput == 'y') {
      cat("Would you like to use the default setting for replicates? [y/n]: ") 
      resp2 = readLines("stdin", 1, warn = FALSE)
    }
    else
    {
      resp2 = defaultAction
    }
    if (resp2 == 'y') {
      cat('Function is running by default...')
    } else {
      stop("please indicate columns of each sample. See SAMPLE_INDEX using $./coex_filter_optparse.r -h")
    }
  } else if (resp == 'n') {
    cat('Function treats each column to a different sample. Running now...\n')
    sample_index = c(c(1:(dim(data)[2])), c(1:(dim(data)[2])))
    data = cbind(data, data)
  } else {
    stop('please answer does your data have replicates? [y/n]:')
  }
} else {
  resp = 'y'
  if(tab_delim == 'y'){ # inefficient but minimal testing
    sample_index = as.numeric(read.csv(sample_index, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
  } else {
    sample_index = as.numeric(read.csv(sample_index, header = FALSE, stringsAsFactors = FALSE))
  }
  if (length(unique(sample_index)) == ncol(data)) {
    sample_index = c(c(1:(dim(data)[2])), c(1:(dim(data)[2])))
    data = cbind(data, data)
  }
}
# Running the function
coex_filter(data, method = method, p_threshold = p_threshold, topnumber = topnumber, sample_index = sample_index, outFileName = outFileName, resp = resp, tsv=tab_delim, genelistFileName = glfn, plotfn = plotfn, jsonfn = jsonfn)
