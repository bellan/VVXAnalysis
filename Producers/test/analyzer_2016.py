### ----------------------------------------------------------------------
### Based on ZZ->4l strategy.
###----------------------------------------------------------------------


### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------

ELEREGRESSION = "None"

### ----------------------------------------------------------------------


### ----------------------------------------------------------------------
### Execute the analysis
### ----------------------------------------------------------------------
import os

PyFilePath = os.environ['CMSSW_BASE'] + "/src/VVXAnalysis/Producers/python/"
execfile(PyFilePath + "analyzer_ZZjj.py")

### ----------------------------------------------------------------------


### ----------------------------------------------------------------------
### Replace module prameters
### ----------------------------------------------------------------------

process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(

        #2016
'/store/mc/RunIISpring16MiniAODv2/GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/30000/9067B609-3C1A-E611-9E4D-0025904C7FBE.root'
#'/store/mc/RunIISpring16MiniAODv1/ttZJets_13TeV_madgraphMLM/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/00F31E7E-6301-E611-A1A8-0017A4770C2C.root'
#'/store/mc/RunIISpring16MiniAODv1/ttZJets_13TeV_madgraphMLM/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/00E65BA2-9502-E611-8EC3-0025907B4F8A.root'
# '/store/mc/RunIISpring16MiniAODv2/VBFToContinZZTo4eJJ_13TeV_phantom_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/60000/0000B6D4-962D-E611-877D-00259074AEB2.root'
# '/store/mc/RunIISpring16MiniAODv1/VBFToContinZZTo4eJJ_13TeV_phantom_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/10000/1431B7BB-6D33-E611-B3BC-0017A4770C00.root'
# '/store/mc/RunIISpring16MiniAODv1/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/046E1689-FD0D-E611-BE90-002590791D60.root'
# '/store/mc/RunIISpring16MiniAODv2/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/0CD58AA0-D12A-E611-850E-00266CFCCD94.root'
))
    #),eventsToProcess = cms.untracked.VEventRange('1:31907:6167667','1:31907:6167734')) #To look in singles events

process.maxEvents.input = -1
