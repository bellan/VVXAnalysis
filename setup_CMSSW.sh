#!/bin/sh

########################################################
# Creates a CMSSW environment and sets up VVXAnalysis  #
#                                                      #
# Author: A. Mecca (alberto.mecca@cern.ch)             #
########################################################

set -e
set -u

CMSSW_VERSION=10_6_26
branch_ZZ=Run2UL_22
branch_VVX=Run2UltraLegacy
checkoutscript=checkout_10X.csh

# Create the CMSSW area
cmsrel CMSSW_${CMSSW_VERSION}

# Fetch the ZZAnalysis setup script
wget https://raw.githubusercontent.com/CJLST/ZZAnalysis/${branch_ZZ}/${checkoutscript}
chmod u+x ${checkoutscript}

# Move to ${CMSSW_BASE}/src and cmsenv
cd CMSSW_${CMSSW_VERSION}/src
cmsenv

# Execute the ZZAnalysis setup script
../../${checkoutscript}

# Clone VVXAnalysis (try with ssh key, then with plain https)
git clone -b ${branch_VVX} git@github.com:bellan/VVXAnalysis.git || git clone -b ${branch_VVX} https://github.com:bellan/VVXAnalysis.git

# Patch CommonLHETools/LHEHandler so that it does not crash with unknown weights
patch CommonLHETools/LHEHandler/src/LHEHandler.cc <<EOF
914d913
<         throw cms::Exception("LHEWeights") << "Don't know what to do with alternate weight id = " << wgtid << " (weightstype == " << weightstype << ")";
EOF

# Compile with SCRAM
scram b -j

# Ensure that a generic python exists; if not, make it a symlink to the python3 used by this CMSSW release
command -v python >/dev/null || ln -s $(which python3) ${CMSSW_BASE}/bin/${SCRAM_ARCH}/python

# Ensure that cmsstyle is installed
python3 -c "import cmsstyle" || python3 -m pip install --user cmsstyle
