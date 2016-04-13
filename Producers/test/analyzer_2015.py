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

print "SELSETUP",SELSETUP
print "BESTCAND_AMONG",BESTCAND_AMONG,"\n"
### ----------------------------------------------------------------------


### ----------------------------------------------------------------------
### Replace module prameters
### ----------------------------------------------------------------------

process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(

        '/store/mc/RunIIFall15MiniAODv2/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/00AE5257-01BA-E511-AA8F-002590596484.root'
        #'/store/mc/RunIIFall15MiniAODv1/ZZTo4L_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/064FDFEF-3DA7-E511-8340-001E6739850C.root'

        ))
    #),eventsToProcess = cms.untracked.VEventRange('1:31907:6167667','1:31907:6167734')) #To look in singles events

process.maxEvents.input = 100
