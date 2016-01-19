# install new version
source("http://bioconductor.org/biocLite.R")
## bioc v 3.0 or later: upgrade packages
biocLite(lib = "[% rlib %]")
#biocLite(ask=FALSE)
#try(biocLite("BiocUpgrade", ask=FALSE, lib = "[% rlib %]" ))

## Now ready to install WGCNA
biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore", "limma"), lib = "[% rlib %]")
install.packages ("WGCNA", repos = "http://cran.wustl.edu", lib = "[% rlib %]")
