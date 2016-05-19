#!/bin/sh

# Remove the old stuff
#mkdir ~/git/pbbamr/src/pbbam/
find ~/git/pbbamr/src/pbbam -name "*.cpp" -exec rm -rf {} \;
find ~/git/pbbamr/src/pbbam -name "*.h" -exec rm -rf {} \;
rm -rf ~/git/pbbamr/inst/include/pbbam
# Copy the new stuff
cp -r ~/git/pbbam/src/ ~/git/pbbamr/src/pbbam/
cp -r ~/git/pbbam/include/pbbam ~/git/pbbamr/inst/include/
# ditch the swig stuff.
rm -r ~/git/pbbamr/src/pbbam/swig