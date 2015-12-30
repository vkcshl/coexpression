#!/bin/bash
dest=${TARGET-/usr/}
runt=${KB_RUNTIME-/usr}/bin
srcd="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "using $dest as installation directory";
echo "using $runt as kb_runtime directory";
echo "using $srcd as package directory";

###
# Warning: the following installation package assumes the latest version installation is stable on WGCNA, 
# which might not be the best if any new piece is broken.
# TODO: Fix all the version of the dependencies

c=$(grep "deb http://cran.rstudio.com/bin/linux/ubuntu trusty" /etc/apt/sources.list | wc -l);
if [ "$c" == 0 ]; then
    echo "Upgrading R version to 3.2 or latest"
    sudo sh -c 'echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list'
    gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
    gpg -a --export E084DAB9 | sudo apt-key add -
    sudo apt-get update
    sudo apt-get -y remove r-base-core r-base-dev
    # ensure the installed packages gone for clean install
    # by default, R would not update the installed packages with the latest version,
    # which causes version compatibility issues.
    sudo rm -rf /usr/lib/R
    sudo rm -rf /usr/local/lib/R
    sudo apt-get -y install r-base
fi

mkdir -p $dest/lib/R/library # for sanity and it actually does not use the created folder
# the kb runtime execution environment has the following variable was set
export R_LIBS=$dest/lib
if [ -e $runt/R ]; then 
	# /kb/runtime case
	tpage --define rlib=$dest/lib "$srcd/r-packages.R"        | $runt/R --vanilla --slave
	tpage --define rlib=$dest/lib "$srcd/r-wgcna-packages.R"  | $runt/R --vanilla --slave
else # docker does not have R on /kb/runtime
	# system default case
	tpage --define rlib=$dest/lib "$srcd/r-packages.R"        | /usr/bin/R --vanilla --slave
	tpage --define rlib=$dest/lib "$srcd/r-wgcna-packages.R"  | /usr/bin/R --vanilla --slave
fi
