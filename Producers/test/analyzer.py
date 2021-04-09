### ----------------------------------------------------------------------
### Based on ZZ->4l strategy.
###----------------------------------------------------------------------


### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------

### ----------------------------------------------------------------------

### ----------------------------------------------------------------------
### Execute the analysis
### ----------------------------------------------------------------------
import os


LEPTON_SETUP = 2018


PyFilePath = os.environ['CMSSW_BASE'] + "/src/VVXAnalysis/Producers/python/"
execfile(PyFilePath + "analyzer_VVjj.py")
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"
execfile(PyFilePath + "prod/pyFragments/RecoProbabilities.py")


### ----------------------------------------------------------------------

### ----------------------------------------------------------------------
### Replace module prameters
### ----------------------------------------------------------------------

process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(''))

 
#2016
if LEPTON_SETUP == 2016:
    process.source.fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV709_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/20000/2637D760-16D4-E711-8612-0026B9277A4C.root')

#2017
if LEPTON_SETUP == 2017:
    process.source.fileNames = cms.untracked.vstring(#'/store/mc/RunIIFall17MiniAODv2/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/54B0F730-A142-E811-BB8F-901B0E6459CA.root'
        '/store/mc/RunIIFall17MiniAODv2/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/C4D6693C-1542-E811-9B5A-1866DA87967B.root'
)

#2018
if LEPTON_SETUP == 2018:
    process.source.fileNames = cms.untracked.vstring('/store/mc/RunIIAutumn18MiniAOD/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v2/100000/6DB7F3BB-8FB7-7F43-9288-9B9A0700F467.root'
#'/store/mc/RunIIAutumn18MiniAOD/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/270000/E5E2F122-AA57-5248-8177-594EC87DD494.root'
    )


process.maxEvents.input = 10


