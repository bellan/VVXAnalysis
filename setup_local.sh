#!/bin/sh

########################################################
# Creates a local setup for VVXAnalysis                #
#                                                      #
# Author: A. Mecca (alberto.mecca@cern.ch)             #
########################################################

set -e
set -u

branch_ZZ=Run2UL_22
branch_VVX=Run2UltraLegacy

# Clone ZZAnalysis
git clone -b ${branch_ZZ} https://github.com/CJLST/ZZAnalysis.git

# Clone VVXAnalysis (try with ssh key, then with plain https)
git clone -b ${branch_VVX} git@github.com:bellan/VVXAnalysis.git || git clone -b ${branch_VVX} https://github.com:bellan/VVXAnalysis.git

# cmake and compile VVXAnalysis
cd VVXAnalysis/TreeAnalysis
cmake CMakeLists.txt
make -j

# Ensure that cmsstyle is installed
python3 -c "import cmsstyle" || python3 -m pip install --user cmsstyle
