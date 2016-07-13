#pbbamr [![CircleCI](https://circleci.com/gh/PacificBiosciences/pbbamr.svg?style=svg)](https://circleci.com/gh/PacificBiosciences/pbbamr)

An R library to view and examine PacBio BAM format files.  Cannot be used on 
windows computers.


##Installation

```
install.packages("devtools")
library(devtools)
install_github("PacificBiosciences/pbbamr", build_vignettes = FALSE)
```
On some machines, a configuration conflict may lead to the
message `Problem with the SSL CA cert (path? access rights?)`.  If this appears,
simply run the following command before those given above:

    httr::set_config( httr::config( ssl_verifypeer = 0L ) )


## Documentation

[The online vignette](http://htmlpreview.github.io/?http://github.com/PacificBiosciences/pbbamr/blob/master/vignettes/pbbamr.html)

[Version 0.3.0 changes](http://rpubs.com/evolvedmicrobe/195680)


## The origin of this package

This package is designed to be a stand alone library that can access PacBio BAM
files.  In order to be stand alone, it has hard forked copies of the following
libraries:

* Boost (1.58)
* HTSLIB (P4, 12/10/2015)
* PBBAM (415cb1a1bcad853978ed37eda3fd4a23aed8e490)

Directly copying the files was done to enable a cleaner build, but none of the
files were edited and can likely be recopied using the set of
commands described below to refill the files for the package.

The package relies on Rcpp, and outside of just copying the files, generating
the Makefile and Makevars file was the only real challenge in adding the
dependencies.  Both htslib and pbbam are intended to be linked in as static
libraries, avoiding any dynamic loading issues, and the new make files are as
such.

###How to update pbbam

The following commands followed by a git commit
can bring you up to date.

```
# Move pbbam over
mkdir ~/git/pbbamr/src/pbbam/
cp -r ~/git/pbbam/src/ ~/git/pbbamr/src/pbbam/
cp -r ~/git/pbbam/include/pbbam ~/git/pbbamr/inst/include/
# ditch the swig stuff.
rm -r ~/git/pbbamr/src/pbbam/swig
```

### Move htslib package over 
```
cd /Users/nigel/p4/ndelaney_Nigels-MacBook-Pro_1663/depot/software/smrtanalysis/bioinformatics/staging/PostPrimary/
mkdir ~/git/pbbamr/src/htslib
cp -r htslib/ ~/git/pbbamr/src/htslib/
mkdir ~/git/pbbamr/inst/include/htslib
mv ~/git/pbbamr/src/htslib/* ~/git/pbbamr/inst/include/htslib/
# Delete the tests
rm -rf ~/git/pbbamr/src/htslib/test
```


### Subset the boost headers we need
```
# bcp is a utility to do this and you may need to download (brew install boost-bcp)
# Run this command to detect the dependencies we need and add them to 
# the include directory
bcp --scan --boost=/Users/nigel/git/bh/local/boost_1_58_0/ src/pbbam/*.cpp inst/include/
bcp --scan --boost=/Users/nigel/git/bh/local/boost_1_58_0/ src/pbbam/*.h inst/include/
bcp --scan --boost=/Users/nigel/git/bh/local/boost_1_58_0/ inst/include/pbbam/*.h inst/include/
bcp --scan --boost=/Users/nigel/git/bh/local/boost_1_58_0/ inst/include/pbbam/internal/*.h inst/include/
# For some reason this didn't copy
cp  boost_1_58_0/boost/utility/string_ref.hpp ~/git/pbbamr/inst/include/boost/utility/
```
