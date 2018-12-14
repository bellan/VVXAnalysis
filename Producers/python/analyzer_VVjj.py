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

execfile(PyFilePath + "MasterPy/ZZ4lAnalysis.py")         # 2016 reference analysis

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



#process.VVjjEventTagger = cms.EDFilter("EventTagger",
#                                       Topology = cms.int32(-1), 
#                                       src = cms.InputTag("softLeptons"),
#                                       TightSelection = cms.string("userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')")
#                                       )

#process.eventTagger = cms.Path(process.VVjjEventTagger)


process.genParticlesFromHardProcess = cms.EDFilter("GenParticleSelector",
                                           filter = cms.bool(False),
                                           src = cms.InputTag("prunedGenParticles"),
                                           #acceptance cut on leptons?
                                           cut = cms.string('status == 1 && (isPromptFinalState && fromHardProcessFinalState && abs(pdgId) >= 11 && abs(pdgId) <= 16) ||  abs(pdgId) == 22'),        
                                           stableOnly = cms.bool(True)
                                           )

# FIXME! They need to be disambiguated from leptons!!
process.selectedGenJets = cms.EDFilter("GenJetSelector",
                                       filter = cms.bool(False),
                                       src = cms.InputTag("slimmedGenJets"),
                                       cut = cms.string('pt > 20 && abs(eta) < 4.7'),
                                       )



process.genPath = cms.Path(process.genParticlesFromHardProcess + process.selectedGenJets)



## ZZ->4l
execfile(VVjj_search_path + "analyzer_ZZjj.py")

## WZ->3lnu
#execfile(VVjj_search_path + "analyzer_WZjj.py") 

## VZ->2l2j
#execfile(VVjj_search_path + "analyzer_VZjj.py")

## WW->2l2nu
#execfile(VVjj_search_path + "analyzer_WWjj.py")

### ----------------------------------------------------------------------




### ......................................................................... ###
### Build collections of muons and electrons that pass a quality criteria (isGood + isolation) and that are NOT selected to form the ZZ best candidate that pass the full selection
### ......................................................................... ###


# Muons cleaning. First, create a muon collection from the best ZZ candidate grand daughters
process.muonsFromZZ = cms.EDProducer("PATMuonsFromCompositeCandidates", src =  cms.InputTag("ZZFiltered"), SplitLevel = cms.int32(1))

# Electrons cleaning. First, create a electron collection from the best ZZ candidate grand daughters
process.electronsFromZZ = cms.EDProducer("PATElectronsFromCompositeCandidates", src =  cms.InputTag("ZZFiltered"), SplitLevel = cms.int32(1))


process.postCleaningMuons = cms.EDFilter("PATMuonSelector", src = cms.InputTag("appendPhotons:muons"),
                                         cut = cms.string("pt > 10 && userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"))


process.postCleaningElectrons = cms.EDFilter("PATElectronSelector", src = cms.InputTag("appendPhotons:electrons"),
                                             cut = cms.string("pt > 10 && userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"))


process.muonsToBeRemovedFromJets = cms.EDProducer("PATMuonMerger",
                                                  src = cms.VInputTag(cms.InputTag("muonsFromZZ"), cms.InputTag("postCleaningMuons")))

process.electronsToBeRemovedFromJets = cms.EDProducer("PATElectronMerger",
                                                      src = cms.VInputTag(cms.InputTag("electronsFromZZ"), cms.InputTag("postCleaningElectrons")))




### ......................................................................... ###
# Remove from the event the jets that have leptons from the ZZ best candidate.
# Jets are also checked against other good isolated leptons not coming from the ZZ best candidate. 
# The jets, to be stored in the event, must pass the preselction (specified below by the user) AND the looseID + PU veto that is implemented in the code (same algo as for H->ZZ VBF selection)
### ......................................................................... ###

## FIXME: Logic need to be recheck as of new FSR strategy has been implemented
process.disambiguatedJets = cms.EDProducer("JetsWithLeptonsRemover",
                                           JetPreselection      = cms.string("pt > 20"),
                                           DiBosonPreselection  = cms.string(""),
                                           MuonPreselection     = cms.string(""),
                                           ElectronPreselection = cms.string(""),
                                           MatchingType         = cms.string("byDeltaR"), 
                                           Jets      = cms.InputTag("dressedJets"),
                                           Muons     = cms.InputTag("muonsToBeRemovedFromJets"),
                                           Electrons = cms.InputTag("electronsToBeRemovedFromJets"),
                                           Diboson   = cms.InputTag(""),
                                           cleanFSRFromLeptons = cms.bool(True)
                                           )


# Number of disambiguated jets
process.jetCounterFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("disambiguatedJets"), minNumber = cms.uint32(0))

process.jetCleaning = cms.Path(process.muonsFromZZ * process.postCleaningMuons * process.muonsToBeRemovedFromJets 
                               + process.electronsFromZZ * process.postCleaningElectrons * process.electronsToBeRemovedFromJets                               
                               + process.disambiguatedJets
                               + process.jetCounterFilter)





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
                                     muons        = cms.InputTag("postCleaningMuons"),     # all good isolated muons BUT the ones coming from ZZ decay
                                     electrons    = cms.InputTag("postCleaningElectrons"), # all good isolated electrons BUT the ones coming from ZZ decay
                                     jets         = cms.InputTag("disambiguatedJets"),     # jets which do not contains leptons from ZZ or other good isolated leptons are removed
                                     Vhad         = cms.InputTag(""),
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
