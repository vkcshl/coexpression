source("http://bioconductor.org/biocLite.R") 
biocLite(c("RSQLite"), lib = "[% rlib %]" )
biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"), lib = "[% rlib %]") 
install.packages ("WGCNA", repos = "http://cran.wustl.edu", lib = "[% rlib %]")
