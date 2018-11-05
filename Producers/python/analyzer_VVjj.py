from ZZAnalysis.AnalysisStep.defaults import declareDefault

SIGNALDEFINITION = int('101',2)  # -1 means get everything, 1 means the request of having a HZZ pair with the  mass in the chosen windows. 11 means the request of having a ZZ pair with the  mass in the chosen windows. For other topology see the README under VVXAnalysis/Commons.

declareDefault("PD","",globals())
declareDefault("XSEC",-1.,globals())
declareDefault("MCFILTER","",globals())
declareDefault("SKIM_REQUIRED",True,globals())
declareDefault("KINREFIT", False, globals())
declareDefault("BESTCANDCOMPARATOR", "byBestZ1bestZ2", globals())
declareDefault("APPLYTRIG", True, globals()) 
declareDefault("VVMODE", 1, globals())
declareDefault("VVDECAYMODE", 0, globals())
declareDefault("ADDLHEKINEMATICS", False, globals())

declareDefault("APPLY_QCD_GGF_UNCERT", False, globals() )

# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

execfile(PyFilePath + "MasterPy/ZZ4lAnalysis.py")         # 2012 reference analysis

### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------


### ----------------------------------------------------------------------
### Output root file
### ----------------------------------------------------------------------

process.TFileService=cms.Service('TFileService', fileName=cms.string('VVXAnalysis.root'))



### ------------------------------- Analyses -----------------------------

### ZZjj paths
VVjj_search_path = os.environ['CMSSW_BASE'] + "/src/VVXAnalysis/Producers/python/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------



#process.goodIsoMuons = cms.EDProducer("CandSelector",
#                                      src = cms.InputTag("appendPhotons:muons"),
#                                      cut = cms.string("userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"))

#process.goodIsoElectrons = cms.EDProducer("CandSelector",
#                                      src = cms.InputTag("appendPhotons:electrons"),
#                                      cut = cms.string("userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"))

process.VVjjEventTagger = cms.EDFilter("EventTagger",
                                       Topology = cms.int32(-1), 
                                       src = cms.InputTag("softLeptons"),
                                       TightSelection = cms.string("userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')")
                                       )

process.path(process.VVjjEventTagger)


## ZZ->4l
execfile(VVjj_search_path + "analyzer_ZZjj.py")

## WZ->3lnu
#execfile(VVjj_search_path + "analyzer_WZjj.py") 

## VZ->2l2j
#execfile(VVjj_search_path + "analyzer_VZjj.py")

## WW->2l2nu
#execfile(VVjj_search_path + "analyzer_WWjj.py")

### ----------------------------------------------------------------------








### ------------------------------------------------------------------------- ###
### Fill the tree for the analysis
### ------------------------------------------------------------------------- ###

 
process.treePlanter = cms.EDAnalyzer("TreePlanter",
                                     sampleName   = cms.string(SAMPLENAME),
                                     setup        = cms.int32(LEPTON_SETUP),
                                     sampleType   = cms.int32(SAMPLE_TYPE),
                                     PD           = cms.string(PD),
                                     skimPaths    = cms.vstring(SkimPaths),
                                     SkimRequired = cms.untracked.bool(SKIM_REQUIRED),
                                     MCFilterPath = cms.string(MCFILTER),
                                     isMC         = cms.untracked.bool(IsMC),
                                     VVMode = cms.int32(int(VVMODE)),
                                     VVDecayMode = cms.int32(int(VVDECAYMODE)),
                                     signalDefinition = cms.int32(SIGNALDEFINITION),
                                     JetAlgo        = cms.string("AK4PFchs"),        
                                     jetResFile_pt  = cms.FileInPath('JRDatabase/textFiles/Fall15_25nsV2_MC/Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt'),
                                     jetResFile_phi = cms.FileInPath('JRDatabase/textFiles/Fall15_25nsV2_MC/Fall15_25nsV2_MC_PhiResolution_AK4PFchs.txt'),
                                     jetResFile_SF  = cms.FileInPath('JRDatabase/textFiles/Fall15_25nsV2_MC/Fall15_25nsV2_MC_SF_AK4PFchs.txt'),
                                     muons        = cms.InputTag("postCleaningMuons"),     # all good isolated muons BUT the ones coming from ZZ decay
                                     electrons    = cms.InputTag("postCleaningElectrons"), # all good isolated electrons BUT the ones coming from ZZ decay
                                     jets         = cms.InputTag("disambiguatedJets"),     # jets which do not contains leptons from ZZ or other good isolated leptons are removed
                                     Vhad         = cms.InputTag("VhadCand"),
                                     ZZ           = cms.InputTag("ZZFiltered"),            # only the best ZZ->4l candidate that pass the FULL selection
                                     ZL           = cms.InputTag("ZlCand"),
                                     MET          = cms.InputTag("slimmedMETs"),
                                     Vertices     = cms.InputTag("goodPrimaryVertices"),                                    
                                     XSection     = cms.untracked.double(XSEC)
     )


### ------------------------------------------------------------------------- ###
### Run the TreePlanter
### ------------------------------------------------------------------------- ###

process.filltrees = cms.EndPath(cms.ignore(process.zzTrigger) + process.srCounter + process.cr2P2FCounter + process.cr3P1FCounter + process.treePlanter)

########################################################################################################################################################################

########################################################
### Tools for further debuging                       ###
########################################################

process.printTree = cms.EDAnalyzer("ParticleListDrawer",
                                   maxEventsToPrint = cms.untracked.int32(-1),
                                   printVertex = cms.untracked.bool(False),
                                   src = cms.InputTag("prunedGenParticles")
                                   )


process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
     dumpTrigger = cms.untracked.bool(True),

     muonSrcs =  cms.PSet(
        slimmedMuons = cms.InputTag("slimmedMuons"),
        calibratedMuons   = cms.InputTag("calibratedMuons"),
        muons        = cms.InputTag("appendPhotons:muons"),
        postCleaningMuons  = cms.InputTag("postCleaningMuons")
     ),

     electronSrcs = cms.PSet(
        slimmedElectron        = cms.InputTag("slimmedElectrons"),
        calibratedPatElectrons = cms.InputTag("calibratedPatElectrons"),
        electrons              = cms.InputTag("appendPhotons:electrons"),
        postCleaningElectrons  = cms.InputTag("postCleaningElectrons")
     ),
     candidateSrcs = cms.PSet(
        Z   = cms.InputTag("ZCand"),
        ZZ  = cms.InputTag("ZZCand"),
        ZLL = cms.InputTag("ZLLCand"),  
        ZL  = cms.InputTag("ZlCand") 
        ),
      jetSrc = cms.PSet(
        cleanJets          = cms.InputTag("cleanJets"),
        JetsWithLeptonsRemover = cms.InputTag("JetsWithLeptonsRemover"),
        disambiguatedJets = cms.InputTag("disambiguatedJets")
        )
)


#process.filltrees = cms.Path(process.preselection * process.genCategory * process.treePlanter * process.printTree)
#process.filltrees = cms.EndPath(process.treePlanter *process.dumpUserData)
#process.filltrees = cms.EndPath(process.treePlanter *process.printTree)

########################################################
