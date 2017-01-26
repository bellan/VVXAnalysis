### ----------------------------------------------------------------------
### Based on ZZ->4l strategy.
###----------------------------------------------------------------------

from ZZAnalysis.AnalysisStep.defaults import declareDefault

### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------

ELEREGRESSION = "None"

### ----------------------------------------------------------------------


### ----------------------------------------------------------------------
### Execute the analysis
### ----------------------------------------------------------------------
import os

declareDefault("APPLYTRIG", False, globals())
declareDefault("SKIM_REQUIRED",False,globals()) # To look into events outside SR and CRs. e.g. Look into MC signal at gen level.
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
#'/store/data/Run2016B/DoubleEG/MINIAOD/23Sep2016-v3/00000/0060C751-C097-E611-9FE6-FA163EFD4308.root'
#MC
' /store/mc/RunIISummer16MiniAODv2/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/040B9D16-B5BD-E611-A70F-0CC47A78A3F8.root'
#'/store/mc/RunIISummer16MiniAODv2/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/221CC46F-2FC6-E611-8FFC-0CC47A1E0488.root'
 
#2016
#Data
#'/store/data/Run2016B/DoubleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/FC69F6B0-F019-E611-90A4-02163E0142B5.root'
#'/store/data/Run2016B/DoubleEG/MINIAOD/23Sep2016-v3/00000/0060C751-C097-E611-9FE6-FA163EFD4308.root'
#MC
#'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/40000/00200284-F15C-E611-AA9B-002590574776.root'
#'/store/mc/RunIISpring16MiniAODv2/ZZJJTo4L_EWK_13TeV-madgraph-pythia8/MINIAODSIM/premix_withHLT_80X_mcRun2_asymptotic_v14-v1/80000/403C4E51-2A72-E611-B589-002590DE7230.root'
#'/store/mc/RunIISpring16MiniAODv2/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/0AA68890-7F1A-E611-B8B0-0CC47A7FC6E6.root'
#'/store/mc/RunIISpring16MiniAODv2/VBFToHiggs0PMContinToZZTo2e2muJJ_M125_GaSM_13TeV_phantom_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/20000/021115D1-983F-E611-A4B2-0CC47A4D768E.root'
#'/store/mc/RunIISpring16MiniAODv2/GluGluToHiggs0MContinToZZTo2e2mu_M125_10GaSM_13TeV_MCFM701_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/50000/02D23032-EE5B-E611-8E99-A0369F301924.root'
#'/store/mc/RunIISpring16MiniAODv1/GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v4/10000/0419BFA2-8C3E-E611-8C5E-00259021A0E2.root'
#'/store/mc/RunIISpring16MiniAODv2/ZZTo4L_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/00000/38182EE2-2352-E611-B861-00259073E3F2.root'
#'/store/mc/RunIISpring16MiniAODv2/GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/30000/9067B609-3C1A-E611-9E4D-0025904C7FBE.root'
#'/store/mc/RunIISpring16MiniAODv1/ttZJets_13TeV_madgraphMLM/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/00F31E7E-6301-E611-A1A8-0017A4770C2C.root'
#'/store/mc/RunIISpring16MiniAODv1/ttZJets_13TeV_madgraphMLM/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/00E65BA2-9502-E611-8EC3-0025907B4F8A.root'
# '/store/mc/RunIISpring16MiniAODv2/VBFToContinZZTo4eJJ_13TeV_phantom_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/60000/0000B6D4-962D-E611-877D-00259074AEB2.root'
# '/store/mc/RunIISpring16MiniAODv1/VBFToContinZZTo4eJJ_13TeV_phantom_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/10000/1431B7BB-6D33-E611-B3BC-0017A4770C00.root'
# '/store/mc/RunIISpring16MiniAODv1/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/046E1689-FD0D-E611-BE90-002590791D60.root'
# '/store/mc/RunIISpring16MiniAODv2/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/0CD58AA0-D12A-E611-850E-00266CFCCD94.root'
))
#    ),eventsToProcess = cms.untracked.VEventRange('1:31907:6167667','1:31907:6167734')) #To look in singles events
#),eventsToProcess = cms.untracked.VEventRange('1:4103554')) #To look in singles events
#    ),eventsToProcess = cms.untracked.VEventRange('273158:37:57587558')) #To look in singles events


#process.filltrees = cms.EndPath(process.treePlanter *process.dumpUserData *process.printTree)
#process.analysis = cms.EndPath(process.dumpUserData)
#process.filltrees = cms.EndPath(process.treePlanter *process.printTree)
