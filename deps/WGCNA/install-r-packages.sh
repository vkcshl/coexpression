mkdir -p $(TARGET)/lib/R/library
tpage --define rlib=$(TARGET)/lib ./deps/r-packages.R  | $(KB_RUNTIME)/bin/R --vanilla --slave

