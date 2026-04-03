#!/bin/bash

CWD=$PWD
MyProject=$1
if [ ! -d $PWD/$MyProject ]; then
    mkdir $MyProject
fi 
echo "checkout of PhysicsTools/NanoAODTools"
wget https://raw.githubusercontent.com/cms-sw/cmssw/master/PhysicsTools/NanoAODTools/standalone/checkoutStandalone.sh
bash checkoutStandalone.sh -d $MyProject
cd $MyProject
echo "checkout of ZZAnalysis"
git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis
(cd ZZAnalysis; git checkout Run3)
echo "checkout of VVXAnalysis"
git clone https://github.com/bellan/VVXAnalysis.git VVXAnalysis
(cd VVXAnalysis; git checkout Run3NanoAOD)
echo "Configuring the environment"
source VVXAnalysis/NanoAnalysis/standalone/env_standalone.sh
init ZZAnalysis NanoAnalysis build
init VVXAnalysis NanoAnalysis build
init VVXAnalysis NanoAnalysis
init PhysicsTools NanoAODTools

