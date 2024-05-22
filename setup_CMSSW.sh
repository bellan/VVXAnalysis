#!/bin/sh

########################################################
# Creates a CMSSW environment and sets up VVXAnalysis  #
#                                                      #
# Author: A. Mecca (alberto.mecca@cern.ch)             #
########################################################

set -e
set -u

CMSSW_VERSION=13_0_16
branch_ZZ=Run3
branch_VVX=el8-py3

# Create the CMSSW area
cmsrel CMSSW_${CMSSW_VERSION}

# Fetch the ZZAnalysis setup script
wget https://raw.githubusercontent.com/CJLST/ZZAnalysis/${branch_ZZ}/checkout_13X.csh
chmod u+x checkout_13X.csh

# Move to ${CMSSW_BASE}/src and cmsenv
cd CMSSW_${CMSSW_VERSION}/src
cmsenv

# Execute the ZZAnalysis setup script
../../checkout_13X.csh

# Clone VVXAnalysis (try with ssh key, then with plain https)
git clone -b ${branch_VVX} git@github.com:bellan/VVXAnalysis.git || git clone -b ${branch_VVX} https://github.com:bellan/VVXAnalysis.git

# Patch CommonLHETools/LHEHandler so that it does not crash with unknown weights
patch CommonLHETools/LHEHandler/src/LHEHandler.cc <<EOF
914d913
<         throw cms::Exception("LHEWeights") << "Don't know what to do with alternate weight id = " << wgtid << " (weightstype == " << weightstype << ")";
EOF
(eval $(MelaAnalytics/setup.sh env); cd CommonLHETools; scram b -f -j)

# Compile
scram b -f -j

# python is python3, specifically the one used by this CMSSW release
command -v python >/dev/null || ln -s $(which python3) ${CMSSW_BASE}/bin/${SCRAM_ARCH}/python

# Ensure that cmsstyle is installed
python3 -c "import cmsstyle" || python3 -m pip install --user cmsstyle
