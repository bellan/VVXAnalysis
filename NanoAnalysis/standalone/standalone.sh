#!/bin/bash

CWD=$PWD
MyProject=$1
echo "checkout of PhysicsTools/NanoAODTools"
wget https://raw.githubusercontent.com/cms-sw/cmssw/master/PhysicsTools/NanoAODTools/standalone/checkoutStandalone.sh
bash checkoutStandalone.sh -d $MyProject
cd $MyProject
echo "checkout of ZZAnalysis"
git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis
(cd ZZAnalysis; git checkout Run3)
cd ..
echo "checkout of VVXAnalysis"
git clone https://github.com/bellan/VVXAnalysis.git VVXAnalysis
(cd VVXAnalysis; git checkout Run3NanoAOD)
echo "Configuring the environment"
source NanoAnalysis/standalone/env_standalone.sh build
