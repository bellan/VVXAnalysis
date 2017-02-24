

### ----------------------------------------------------------------------
### Based on ZZ->4l strategy.
###----------------------------------------------------------------------

from ZZAnalysis.AnalysisStep.defaults import declareDefault

### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------

### ----------------------------------------------------------------------


### ----------------------------------------------------------------------
### Execute the analysis
### ----------------------------------------------------------------------
import os

#declareDefault("APPLYTRIG", False, globals())
#declareDefault("SKIM_REQUIRED",False,globals()) # To look into events outside SR and CRs. e.g. Look into MC signal at gen level.
#declareDefault("IsMC",False,globals())  # To run data. Remember that you sholu need the json file for run only the right events.
#PD= "DoubleEle" # Choose the right one. Look the data you are running. "DoubleEle","DoubleMu","MuEG","SingleElectron","SingleMuon" 

PyFilePath = os.environ['CMSSW_BASE'] + "/src/VVXAnalysis/Producers/python/"
execfile(PyFilePath + "analyzer_ZZjj.py")


process.calibratedPatElectrons.isSynchronization = cms.bool(True)
process.calibratedMuons.isSynchronization = cms.bool(True)

### ----------------------------------------------------------------------
### Replace module prameters
### ----------------------------------------------------------------------

process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(


#2017
#data
#'/store/data/Run2016H/DoubleMuon/MINIAOD/PromptReco-v1/000/281/085/00000/9292AFD0-897F-E611-B2F1-02163E014284.root'
#'/store/data/Run2016E/DoubleMuon/MINIAOD/23Sep2016-v1/100000/0049B6D2-278C-E611-AEB4-0025902D944E.root'
#'/store/data/Run2016C/DoubleEG/MINIAOD/23Sep2016-v1/50000/00B319ED-4A86-E611-8BC7-0CC47A4D7634.root'
#'/store/data/Run2016B/DoubleEG/MINIAOD/23Sep2016-v3/00000/0060C751-C097-E611-9FE6-FA163EFD4308.root'
#MC
# '/store/mc/RunIISummer16MiniAODv2/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/22F32262-3FC5-E611-B373-D4AE526DEDB7.root'
'/store/mc/RunIISummer16MiniAODv2/ZZTo4L_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/3C331910-FED5-E611-AE8B-0025905A60CA.root'
#'/store/mc/RunIISummer16MiniAODv2/ZZTo4L_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/5CE2841E-0AD6-E611-BF31-0025905B85BC.root',
#'/store/mc/RunIISummer16MiniAODv2/ZZTo4L_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/62F0599D-8DD6-E611-9B36-0CC47AA992C8.root'
#' /store/mc/RunIISummer16MiniAODv2/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/040B9D16-B5BD-E611-A70F-0CC47A78A3F8.root'
#'/store/mc/RunIISummer16MiniAODv2/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/221CC46F-2FC6-E611-8FFC-0CC47A1E0488.root',
#'/store/mc/RunIISummer16MiniAODv2/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/4E41610B-14C6-E611-97B1-001E67580704.root'


)) 
#   ),eventsToProcess = cms.untracked.VEventRange('1:49904:14039512')) #To look in singles events
#   ),eventsToProcess = cms.untracked.VEventRange('1:54597:15360049')) #To look in singles events



#process.filltrees = cms.EndPath(process.treePlanter *process.dumpUserData *process.printTree)
#process.analysis = cms.EndPath(process.dumpUserData)
#process.filltrees = cms.EndPath(process.treePlanter *process.printTree)
