#try(update.packages(checkBuilt=TRUE, ask=FALSE, repos = "http://cran.wustl.edu", lib = "[% rlib %]"))
# upgrade any existing CRAN package to the latest
try(update.packages(checkBuilt=TRUE, ask=FALSE, repos = "http://cran.wustl.edu"))
# remove bio conductor, so that it would have less issues
try(remove.packages("BiocInstaller", ask=FALSE))
