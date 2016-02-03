#try(update.packages(checkBuilt=TRUE, ask=FALSE, repos = "http://cran.wustl.edu", lib = "[% rlib %]"))
# upgrade any existing CRAN package to the latest
try(update.packages(checkBuilt=TRUE, ask=FALSE, repos = "http://cran.wustl.edu", lib = "[% rlib %]"))
# remove bio conductor, so that it would have less issues
try(remove.packages("BiocInstaller"))
install.packages("optparse", repos = "http://cran.wustl.edu", lib = "[% rlib %]")
install.packages("reshape", repos = "http://cran.wustl.edu", lib = "[% rlib %]")
install.packages("flashClust", repos = "http://cran.wustl.edu", lib = "[% rlib %]")
install.packages("igraph", repos = "http://cran.wustl.edu", lib = "[% rlib %]")
install.packages("jsonlite", repos = "http://cran.wustl.edu", lib = "[% rlib %]")
