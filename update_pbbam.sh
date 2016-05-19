#!/bin/sh
#mkdir ~/git/pbbamr/src/pbbam/
cp -r ~/git/pbbam/src/ ~/git/pbbamr/src/pbbam/
cp -r ~/git/pbbam/include/pbbam ~/git/pbbamr/inst/include/
# ditch the swig stuff.
rm -r ~/git/pbbamr/src/pbbam/swig