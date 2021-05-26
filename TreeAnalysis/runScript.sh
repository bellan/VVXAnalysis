#!/bin/bash

# Calls the Analyzer (via run.py) over all the samples in "samples/"
# Last Author: A. Mecca (alberto.mecca@edu.unito.it)
# Usage: without arguments, it uses default values. Otherwise it takes the
# first as the Analyzer and passes the others unchanged as parameters to run.py

ANALYZER=VZZAnalyzer #default values
EXTRA_ARGS='-rMC -l137100'  # '-rMC -n10000'

if [ $# -gt 0 ]; then
	ANALYZER=$1
fi
if [ $# -gt 1 ]; then
	EXTRA_ARGS="${@:2}" #All the aguments except the first
fi

for SAMPLE in samples/2016/* ; do
	SAMPLE=${SAMPLE##*/}
	SAMPLE=${SAMPLE%.root}
	./python/run.py ${ANALYZER} ${SAMPLE} ${EXTRA_ARGS}
	#echo "${ANALYZER} ${SAMPLE} ${EXTRA_ARGS}"
done

#./python/run.py ZZjAnalyzer  ZZTo4e 
#./python/run.py ZZjAnalyzer  ZZTo4mu 
#./python/run.py ZZjAnalyzer  ZZTo2e2mu
#./python/run.py ZZjAnalyzer  ZZTo4eJJ_SMHContinInterf_H125.6
#./python/run.py ZZjAnalyzer  ZZTo2e2muJJ_SMHContinInterf_H125.6
#./python/run.py ZZjAnalyzer  ZZTo4muJJ_SMHContinInterf_H125.6
#./python/run.py ZZjAnalyzer    ggTo4e_SMHContinInterf-MCFM67_H125.6
#./python/run.py ZZjAnalyzer    ggTo4mu_SMHContinInterf-MCFM67_H125.6
#./python/run.py ZZjAnalyzer    ggTo2e2mu_SMHContinInterf-MCFM67_H125.6
#./python/run.py ZZjAnalyzer    ZZZJets
#./python/run.py ZZjAnalyzer    WZZJets
#./python/run.py ZZjAnalyzer    ttH126
#./python/run.py ZZjAnalyzer    ZH
#./python/run.py ZZjAnalyzer    WH
#./python/run.py ZZjAnalyzer    VBFH126
#./python/run.py ZZjAnalyzer    WZ
#./python/run.py ZZjAnalyzer    WWWJets
#./python/run.py ZZjAnalyzer    TTTo2L2Nu2B
#./python/run.py ZZjAnalyzer    TTWJets
#./python/run.py ZZjAnalyzer    TTGJets
#./python/run.py ZZjAnalyzer    TTWWJets
#./python/run.py ZZjAnalyzer    WWZJets
#./python/run.py ZZjAnalyzer    TTZJets
#./python/run.py ZZjAnalyzer    WWJets
#./python/run.py ZZjAnalyzer    WGToLNuG
#./python/run.py ZZjAnalyzer    WWGJets 

#./python/run.py ZZjAnalyzer  MuEGA
#./python/run.py ZZjAnalyzer  MuEGB
#./python/run.py ZZjAnalyzer  MuEGC
#./python/run.py ZZjAnalyzer  MuEGD
#./python/run.py ZZjAnalyzer DoubleMuA
#./python/run.py ZZjAnalyzer DoubleMuB
#./python/run.py ZZjAnalyzer DoubleMuC
#./python/run.py ZZjAnalyzer DoubleMuD
#./python/run.py ZZjAnalyzer DoubleEleA
#./python/run.py ZZjAnalyzer DoubleEleB
#./python/run.py ZZjAnalyzer DoubleEleC
#./python/run.py ZZjAnalyzer DoubleEleD

#./python/run.py ZZjAnalyzer  ZZTo4e -r CR2P2F 
#./python/run.py ZZjAnalyzer  ZZTo4mu -r CR2P2F
#./python/run.py ZZjAnalyzer  ZZTo2e2mu -r CR2P2F
#./python/run.py ZZjAnalyzer  ZZTo4eJJ_SMHContinInterf_H125.6  -r CR2P2F
#./python/run.py ZZjAnalyzer  ZZTo2e2muJJ_SMHContinInterf_H125.6 -r CR2P2F
#./python/run.py ZZjAnalyzer  ZZTo4muJJ_SMHContinInterf_H125.6 -r CR2P2F
#./python/run.py ZZjAnalyzer    ggTo4e_SMHContinInterf-MCFM67_H125.6 -r CR2P2F
#./python/run.py ZZjAnalyzer    ggTo4mu_SMHContinInterf-MCFM67_H125.6 -r CR2P2F
#./python/run.py ZZjAnalyzer    ggTo2e2mu_SMHContinInterf-MCFM67_H125.6 -r CR2P2F
#./python/run.py ZZjAnalyzer    ZZZJets -r CR2P2F
#./python/run.py ZZjAnalyzer    WZZJets -r CR2P2F
#./python/run.py ZZjAnalyzer    ttH126 -r CR2P2F
#./python/run.py ZZjAnalyzer    ZH -r CR2P2F
#./python/run.py ZZjAnalyzer    WH -r CR2P2F
#./python/run.py ZZjAnalyzer    VBFH126 -r CR2P2F
#./python/run.py ZZjAnalyzer    WZ -r CR2P2F
#./python/run.py ZZjAnalyzer    WWWJets -r CR2P2F
#./python/run.py ZZjAnalyzer    TTTo2L2Nu2B -r CR2P2F
#./python/run.py ZZjAnalyzer    TTWJets -r CR2P2F
#./python/run.py ZZjAnalyzer    TTGJets -r CR2P2F
#./python/run.py ZZjAnalyzer    TTWWJets -r CR2P2F
#./python/run.py ZZjAnalyzer    WWZJets -r CR2P2F
#./python/run.py ZZjAnalyzer    TTZJets -r CR2P2F
#./python/run.py ZZjAnalyzer    WWJets -r CR2P2F
#./python/run.py ZZjAnalyzer    WGToLNuG -r CR2P2F
#./python/run.py ZZjAnalyzer    WWGJets  -r CR2P2F

#./python/run.py ZZjAnalyzer  MuEGA -r CR2P2F
#./python/run.py ZZjAnalyzer  MuEGB -r CR2P2F
#./python/run.py ZZjAnalyzer  MuEGC -r CR2P2F
#./python/run.py ZZjAnalyzer  MuEGD -r CR2P2F
#./python/run.py ZZjAnalyzer DoubleMuA -r CR2P2F
#./python/run.py ZZjAnalyzer DoubleMuB -r CR2P2F
#./python/run.py ZZjAnalyzer DoubleMuC -r CR2P2F
#./python/run.py ZZjAnalyzer DoubleMuD -r CR2P2F
#./python/run.py ZZjAnalyzer DoubleEleA -r CR2P2F
#./python/run.py ZZjAnalyzer DoubleEleB -r CR2P2F
#./python/run.py ZZjAnalyzer DoubleEleC -r CR2P2F
#./python/run.py ZZjAnalyzer DoubleEleD -r CR2P2F
 
#./python/run.py ZZjAnalyzer  ZZTo4e -r CR3P1F
#./python/run.py ZZjAnalyzer  ZZTo4mu -r CR3P1F
#./python/run.py ZZjAnalyzer  ZZTo2e2mu -r CR3P1F
#./python/run.py ZZjAnalyzer  ZZTo4eJJ_SMHContinInterf_H125.6 -r CR3P1F
#./python/run.py ZZjAnalyzer  ZZTo2e2muJJ_SMHContinInterf_H125.6 -r CR3P1F
#./python/run.py ZZjAnalyzer  ZZTo4muJJ_SMHContinInterf_H125.6 -r CR3P1F
#./python/run.py ZZjAnalyzer    ggTo4e_SMHContinInterf-MCFM67_H125.6 -r CR3P1F
#./python/run.py ZZjAnalyzer    ggTo4mu_SMHContinInterf-MCFM67_H125.6 -r CR3P1F
#./python/run.py ZZjAnalyzer    ggTo2e2mu_SMHContinInterf-MCFM67_H125.6 -r CR3P1F
#./python/run.py ZZjAnalyzer    ZZZJets -r CR3P1F
#./python/run.py ZZjAnalyzer    WZZJets -r CR3P1F
#./python/run.py ZZjAnalyzer    ttH126 -r CR3P1F
#./python/run.py ZZjAnalyzer    ZH -r CR3P1F
#./python/run.py ZZjAnalyzer    WH -r CR3P1F
#./python/run.py ZZjAnalyzer    VBFH126 -r CR3P1F
#./python/run.py ZZjAnalyzer    WZ -r CR3P1F
#./python/run.py ZZjAnalyzer    WWWJets -r CR3P1F
#./python/run.py ZZjAnalyzer    TTTo2L2Nu2B -r CR3P1F
#./python/run.py ZZjAnalyzer    TTWJets -r CR3P1F
#./python/run.py ZZjAnalyzer    TTGJets -r CR3P1F
#./python/run.py ZZjAnalyzer    TTWWJets -r CR3P1F
#./python/run.py ZZjAnalyzer    WWZJets -r CR3P1F
#./python/run.py ZZjAnalyzer    TTZJets -r CR3P1F
#./python/run.py ZZjAnalyzer    WWJets -r CR3P1F
#./python/run.py ZZjAnalyzer    WGToLNuG -r CR3P1F
#./python/run.py ZZjAnalyzer    WWGJets -r CR3P1F 

#./python/run.py ZZjAnalyzer  MuEGA -r CR3P1F
#./python/run.py ZZjAnalyzer  MuEGB -r CR3P1F
#./python/run.py ZZjAnalyzer  MuEGC -r CR3P1F
#./python/run.py ZZjAnalyzer  MuEGD -r CR3P1F
#./python/run.py ZZjAnalyzer DoubleMuA -r CR3P1F
#./python/run.py ZZjAnalyzer DoubleMuB -r CR3P1F
#./python/run.py ZZjAnalyzer DoubleMuC -r CR3P1F
#./python/run.py ZZjAnalyzer DoubleMuD -r CR3P1F
#./python/run.py ZZjAnalyzer DoubleEleA -r CR3P1F
#./python/run.py ZZjAnalyzer DoubleEleB -r CR3P1F
#./python/run.py ZZjAnalyzer DoubleEleC -r CR3P1F
#./python/run.py ZZjAnalyzer DoubleEleD -r CR3P1F


