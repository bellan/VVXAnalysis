from ZZAnalysis.AnalysisStep.defaults import declareDefault

SIGNALDEFINITION = int('1',2)  # -1 means get everything, 1 means the request of having a ZZ pair with the  mass in the chosen windows. For other topology see the README under VVXAnalysis/Commons.

declareDefault("PD","",globals())
declareDefault("MCFILTER","",globals())
declareDefault("XSEC",-1,globals()) # was -1. FIXME
declareDefault("GENXSEC", -1, globals()) # FIXME
declareDefault("GENBR", -1, globals()) # FIXME
declareDefault("SKIM_REQUIRED",True,globals())
declareDefault("KINREFIT", False, globals())
declareDefault("BESTCANDCOMPARATOR", "byBestZ1bestZ2", globals())
declareDefault("APPLYTRIG", True, globals()) 
declareDefault("VVMODE", 1, globals())
declareDefault("VVDECAYMODE", 0, globals())
declareDefault("ADDLHEKINEMATICS", False, globals())
#declareDefault("PROCESS_CR", False, globals())
declareDefault("ADDZTREE", False, globals())


declareDefault("APPLY_QCD_GGF_UNCERT", False, globals() )

# K factors
declareDefault("APPLY_K_NNLOQCD_ZZGG", 0, globals()) # 0: Do not; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
declareDefault("APPLY_K_NNLOQCD_ZZQQB", False, globals())
declareDefault("APPLY_K_NLOEW_ZZQQB", False, globals())

#failed events
declareDefault("SKIP_EMPTY_EVENTS", True, globals())
declareDefault("FAILED_TREE_LEVEL", 0, globals())

if FAILED_TREE_LEVEL and not SKIP_EMPTY_EVENTS:
    raise ValueError(
                     "Inconsistent options: FAILED_TREE_LEVEL={}, SKIP_EMPTY_EVENTS={}\n"
                     "If you want to write a failed tree, set SKIP_EMPTY_EVENTS=True"
                     .format(FAILED_TREE_LEVEL, SKIP_EMPTY_EVENTS)
                    )


# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Import ZZAnalysis
### ----------------------------------------------------------------------

execfile(PyFilePath + "MasterPy/ZZ4lAnalysis.py")      

### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------


### ----------------------------------------------------------------------
### Output root file
### ----------------------------------------------------------------------

#process.TFileService=cms.Service('TFileService', fileName=cms.string('VVXAnalysis.root'))
process.TFileService=cms.Service('TFileService', fileName=cms.string('ZZ4lAnalysis.root'))



### ------------------------------- Analyses -----------------------------

### ZZjj paths
VVjj_search_path = os.environ['CMSSW_BASE'] + "/src/VVXAnalysis/Producers/python/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

### ---------------------------------- MC --------------------------------
if IsMC:
    
    genparticles_cut = '(status == 1 && isPromptFinalState && fromHardProcessFinalState && abs(pdgId) <= 22)'
    # genquarks  = '(isPromptFinalState && fromHardProcessFinalState && abs(pdgId) <= 6)'  # && (mother.pdgId == 23) || abs(mother.pdgId) == 24)
    
    process.genParticlesFromHardProcess = cms.EDFilter("GenParticleSelector", filter = cms.bool(False), stableOnly = cms.bool(False),
                                                       src = cms.InputTag("prunedGenParticles"),
                                                       cut = cms.string(genparticles_cut))
    
    # process.genPhotons = cms.EDFilter("GenParticleSelector", filter = cms.bool(False), stableOnly = cms.bool(False),
    #                                   src = cms.InputTag("prunedGenParticles"),
    #                                   cut = cms.string('status == 1 && pdgId == 22'))
    
    process.genTaus = cms.EDFilter("GenParticleSelector", filter = cms.bool(False), stableOnly = cms.bool(False),
                                   src = cms.InputTag("prunedGenParticles"),
                                   cut = cms.string('isPromptDecayed && abs(pdgId) == 15'))


    # FIXME! They need to be disambiguated from leptons!! RB: done in the signal definition!
    process.selectedGenJets = cms.EDFilter("GenJetSelector", filter = cms.bool(False),
                                           src = cms.InputTag("slimmedGenJets"),
                                           cut = cms.string('pt > 20 && abs(eta) < 4.7'))
                                    

    process.selectedGenJetsAK8 = cms.EDFilter("GenJetSelector", filter = cms.bool(False),
                                              src = cms.InputTag("slimmedGenJetsAK8"),
                                              cut = cms.string('pt > 20 && abs(eta) < 4.7'))
                                          
    process.genPath = cms.Path(process.genParticlesFromHardProcess
                               # + process.genPhotons
                               + process.selectedGenJets
                               + process.selectedGenJetsAK8
                               + process.genTaus)


### ---------- If it is MC, run also the signal definition path -----------

    # Empty sequence to attach the signal filter (if specified in the CSV file)
    process.mcSelectionCounter = cms.EDProducer("EventCountProducer") # not really needeed... it is mainly an hack to get the path executed
    process.signalFilters = cms.Sequence(process.mcSelectionCounter) 
    process.mcSelection   = cms.Path(process.signalFilters)
    MCFILTER = "mcSelection"

    process.genCategory =  cms.EDFilter("VVXGenFilterCategory",
                                Topology       = cms.int32(SIGNALDEFINITION), 
                                src            = cms.InputTag("genParticlesFromHardProcess"),
                                GenJets        = cms.InputTag("selectedGenJets"),
                                GenJetsAK8     = cms.InputTag("selectedGenJetsAK8"),
                                )


    process.kFactor = cms.EDProducer('kfactorProducer',
                                     isMC  = cms.untracked.bool(IsMC),
                                     src   = cms.InputTag("prunedGenParticles")) # RB: switch to genParticlesFromHardProcess ??
    
 
    process.signalCounter    = cms.EDProducer("EventCountProducer")
    process.signalDefinition = cms.Path(process.genCategory * process.kFactor * process.signalCounter)

### ---------------------------------------------------------------------



### ------------------------------- Photons -----------------------------
process.load("RecoEgamma/PhotonIdentification/photonIDValueMapProducer_cff")  # this creates a process.photonIDValueMapProducer

process.filledPhotons = cms.EDProducer("PhillerVVX",
                                       # photons = cms.InputTag("filteredPhotons"),
                                       photons = cms.InputTag("slimmedPhotons"),
                                       chIsoValueMap = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation")
)

process.filteredPhotons = cms.EDFilter("PATPhotonSelector",
                                       src = cms.InputTag("slimmedPhotons"),
                                       # src = cms.InputTag("filledPhotons"),
                                       cut = cms.string("pt > 15 && abs(eta) < 3.0")
                                      )

process.photonSelection = cms.Path(
    # process.photonIDValueMapProducer *
    # process.filledPhotons *
    process.filteredPhotons
)
### ---------------------------------------------------------------------





### ------------------------------- AK8 jets -----------------------------
AK8_JEC_tag = None
if IsMC:
    if   (SAMPLE_TYPE == 2016):
        #AK8_JEC_tag    = 'JetCorrectorParametersCollection_Summer16_07Aug2017_V11_MC_AK8PFPuppi', 
        #AK8_JEC_tag    = 'JetCorrectorParametersCollection_Summer19UL16APV_V7_MC_AK8PFPuppi', # APV
        AK8_JEC_tag    = 'JetCorrectorParametersCollection_Summer19UL16_V7_MC_AK8PFPuppi' # NON APV
    elif (SAMPLE_TYPE == 2017):
        #AK8_JEC_tag    = 'JetCorrectorParametersCollection_Fall17_17Nov2017_V32_94X_MC_AK8PFPuppi', #FIXME: need to be tested
        AK8_JEC_tag    = 'JetCorrectorParametersCollection_Summer19UL17_V5_MC_AK8PFPuppi'
    elif (SAMPLE_TYPE == 2018):
        #AK8_JEC_tag    = 'JetCorrectorParametersCollection_Autumn18_V19_MC_AK8PFPuppi',
        AK8_JEC_tag    = 'JetCorrectorParametersCollection_Summer19UL18_V5_MC_AK8PFPuppi'

else:
    if   (SAMPLE_TYPE == 2016):
        #AK8_JEC_tag    = 'JetCorrectorParametersCollection_Summer16_07Aug2017All_V11_DATA_AK8PFPuppi', #for 80X/Moriond17
        AK8_JEC_tag    = 'JetCorrectorParametersCollection_Summer19UL16_RunBCDEFGH_Combined_V7_DATA_AK8PFPuppi'
    elif (SAMPLE_TYPE == 2017):
        #AK8_JEC_tag    = 'JetCorrectorParametersCollection_Fall17_17Nov2017_V32_94X_DATA_AK8PFPuppi', #FIXME: need to be tested
        AK8_JEC_tag    = 'JetCorrectorParametersCollection_Summer19UL17_RunBCDEF_V5_DATA_AK8PFPuppi'
    elif (SAMPLE_TYPE == 2018):
        #AK8_JEC_tag    = 'JetCorrectorParametersCollection_Autumn18_RunABCD_V19_DATA_AK8PFPuppi',
        AK8_JEC_tag    = 'JetCorrectorParametersCollection_Summer19UL18_V5_DATA_AK8PFPuppi'


if AK8_JEC_tag is not None:
    process.jec.toGet.append(cms.PSet( record = cms.string('JetCorrectionsRecord'),
                                       tag    = cms.string(AK8_JEC_tag),
                                       label  = cms.untracked.string('AK8PFPuppi')
                                   ))
else:
    print "UNKNOWN YEAR", SAMPLE_TYPE    



process.patJetCorrFactorsReapplyJECAK8 = updatedPatJetCorrFactors.clone(
    src     = cms.InputTag("slimmedJetsAK8"),
    levels  = ['L1FastJet','L2Relative','L3Absolute'],
    payload = 'AK8PFPuppi'
)
if not IsMC:
    process.patJetCorrFactorsReapplyJECAK8.levels = ['L1FastJet','L2Relative','L3Absolute', 'L2L3Residual']

process.patJetsReapplyJECAK8 = updatedPatJets.clone(
    jetSource = cms.InputTag("slimmedJetsAK8"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJECAK8") )
)


from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.goodJetsAK8 = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                   filterParams = pfJetIDSelector.clone(),
                                   src = cms.InputTag("patJetsReapplyJECAK8"),
                                   #src = cms.InputTag("selectedUpdatedPatJetsAK8WithDeepTags"),
                                   filter = cms.bool(False) )


process.correctedJetsAK8 = cms.EDProducer("CorrJetsProducer",
                                          year    = cms.int32  (LEPTON_SETUP),
                                          jets    = cms.InputTag( "goodJetsAK8" ), # FIXME check with Roberto, it was cleanJetsFat/AK8
                                          vertex  = cms.InputTag( "goodPrimaryVertices" ), 
                                          rho     = cms.InputTag( "fixedGridRhoFastjetAll"   ),
                                          payload = cms.string  ( "AK8PFPuppi" ),
                                          isData  = cms.bool    (  not IsMC )) # FIXME check with Roberto


from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from RecoBTag.ONNXRuntime.pfDeepBoostedJet_cff import _pfDeepBoostedJetTagsAll, _pfDeepBoostedJetTagsProbs, _pfDeepBoostedJetTagsMetaDiscrs, _pfMassDecorrelatedDeepBoostedJetTagsProbs, _pfMassDecorrelatedDeepBoostedJetTagsMetaDiscrs
from RecoBTag.ONNXRuntime.pfParticleNet_cff import _pfParticleNetJetTagsAll

updateJetCollection(
    process,
    jetSource = cms.InputTag('correctedJetsAK8'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    rParam = 0.8,
    jetCorrections = ('AK8PFPuppi', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'), #, 'L2Relative', 'L3Absolute', 'L2L3Residual'
    btagDiscriminators = _pfDeepBoostedJetTagsAll + _pfParticleNetJetTagsAll + ['pfDeepFlavourJetTags:probb','pfDeepFlavourJetTags:probc','pfDeepFlavourJetTags:probg','pfDeepFlavourJetTags:problepb','pfDeepFlavourJetTags:probbb','pfDeepFlavourJetTags:probuds'],
    postfix='AK8WithDeepTags',
    printWarning = True
   )

patAlgosToolsTask = getPatAlgosToolsTask(process)
process.outpathPAT = cms.EndPath(patAlgosToolsTask)


#process.fatJets = cms.Sequence(process.patJetCorrFactorsReapplyJECAK8 + process.patJetsReapplyJECAK8 + process.goodJetsAK8 + process.correctedJetsAK8 + process.disambiguatedJetsAK8)
process.fatJets = cms.Path(process.patJetCorrFactorsReapplyJECAK8 + process.patJetsReapplyJECAK8 + process.goodJetsAK8 + process.correctedJetsAK8)
### ---------------------------------------------------------------------

## targetting ZZ->4l
execfile(VVjj_search_path + "4_leptons_regions.py") 

## targetting WZ->3lnu
execfile(VVjj_search_path + "3_leptons_regions.py")

## VZ->2l2j
execfile(VVjj_search_path + "2_leptons_regions.py")

## WW->2l2nu
#execfile(VVjj_search_path + "analyzer_WWjj.py")

### ----------------------------------------------------------------------

process.counters = cms.Task(process.SR4PCounter , process.CR3P1FCounter   , process.CR2P2FCounter, process.SR4P1LCounter,
                            process.SRHZZCounter, process.CR3P1FHZZCounter, process.CR2P2FHZZCounter,
                            process.SR3PCounter, process.CR110Counter, process.CR101Counter, process.CR011Counter, process.CR100Counter, process.CR001Counter, process.CR010Counter, process.CR000Counter, process.SR3P1LCounter,
                            process.CRLFRCounter,
                            process.SR2PCounter, process.SR2P1LCounter)



### ......................................................................... ###
### Build collections of muons and electrons that pass a quality criteria (isGood + isolation) and that are NOT selected to form the ZZ best candidate that pass the full selection
### ......................................................................... ###


process.pogMuons     = cms.EDFilter("PATMuonSelector", 
                                    src = cms.InputTag("appendPhotons:muons"),
                                    #cut = cms.string("pt > 10 && userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"))
### ID as PKS
                                    cut = cms.string("pt > 10 && abs(eta) < 2.4 && passed('CutBasedIdTight') && (pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt < 0.15"))

process.pogElectrons = cms.EDFilter("PATElectronSelector", 
                                    src = cms.InputTag("appendPhotons:electrons"),
                                    #cut = cms.string("pt > 10 && userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"))
                                    cut = cms.string("pt > 10 && abs(eta) < 2.5 && userInt('cutBasedElectronID-Fall17-94X-V2-medium') == 1023 && ((abs(eta) < 1.479 && userFloat('dxy')<0.05 && userFloat('dz')<0.1) || (abs(eta)>1.479 && userFloat('dxy')<0.1 && userFloat('dz')<0.2))"))  # userInt() == 1023 <--> electronID('cutBasedElectronID_Fall17_94X_V2_medium'), since it is a bitset containing 9 cuts



process.pogIdLeptons = cms.Path(process.pogMuons + process.pogElectrons)


process.muonsToBeRemovedFromJets = cms.EDProducer("PATMuonCollectionMerger",
                                                  src = cms.VInputTag(cms.InputTag("muonsFromZZ"), cms.InputTag("muonsFromZW"), cms.InputTag("muonsFromZV")))

process.electronsToBeRemovedFromJets = cms.EDProducer("PATElectronCollectionMerger",
                                                      src = cms.VInputTag(cms.InputTag("electronsFromZZ"), cms.InputTag("electronsFromZW"), cms.InputTag("electronsFromZV")))




### ......................................................................... ###
# Remove from the event the jets that have leptons from the ZZ best candidate.
# Jets are also checked against other good isolated leptons not coming from the ZZ best candidate. 
# The jets, to be stored in the event, must pass the preselection (specified below by the user) AND the looseID + PU veto that is implemented in the code (same algo as for H->ZZ VBF selection)
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


process.disambiguatedJetsAK8 = cms.EDProducer("JetsWithLeptonsRemover",
                                              JetPreselection      = cms.string("pt > 20"),
                                              DiBosonPreselection  = cms.string(""),
                                              MuonPreselection     = cms.string(""),
                                              ElectronPreselection = cms.string(""),
                                              MatchingType         = cms.string("byDeltaR"), 
                                              DeltaRCut = cms.untracked.double(0.8),
                                              Jets      = cms.InputTag("correctedJetsAK8:corrJets"), # need to create AK8 dressed jets???
                                              Muons     = cms.InputTag("muonsToBeRemovedFromJets"),
                                              Electrons = cms.InputTag("electronsToBeRemovedFromJets"),
                                              Diboson   = cms.InputTag(""),
                                              cleanFSRFromLeptons = cms.bool(True)
                                              )







# Number of disambiguated jets
process.jetCounterFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("disambiguatedJets"), minNumber = cms.uint32(0))

process.jetCleaning = cms.Path(  process.muonsToBeRemovedFromJets 
                               + process.electronsToBeRemovedFromJets
                               + process.disambiguatedJets + process.disambiguatedJetsAK8
                               + process.jetCounterFilter)


### ------------------------------------------------------------------------- ###
### Fill the tree for the analysis
### ------------------------------------------------------------------------- ###

 
process.treePlanter = cms.EDAnalyzer("TreePlanter",
                                     sampleName   = cms.string(SAMPLENAME),
                                     setup        = cms.int32(LEPTON_SETUP),
                                     sampleType   = cms.int32(SAMPLE_TYPE),
                                     PD           = cms.string(PD),
                                     dataTag      = cms.string(DATA_TAG),                  #added for recognizing UL16 pre/post VFP
                                     skimPaths    = cms.vstring(SkimPaths),
                                     SkimRequired = cms.untracked.bool(SKIM_REQUIRED),
                                     MCFilterPath = cms.string(MCFILTER),
                                     isMC         = cms.untracked.bool(IsMC),
                                     VVMode = cms.int32(int(VVMODE)),
                                     VVDecayMode = cms.int32(int(VVDECAYMODE)),
                                     signalDefinition = cms.int32(SIGNALDEFINITION),
                                     AddLHEKinematics = cms.bool(ADDLHEKINEMATICS),
                                     muons        = cms.InputTag("pogMuons"),     # For comparison
                                     electrons    = cms.InputTag("pogElectrons"), # For comparison
                                     photons      = cms.InputTag("filteredPhotons"),       # all photons that pass pt cut
                                     jets         = cms.InputTag("disambiguatedJets"),     # jets which do not contains leptons from ZZ or other good isolated leptons are removed
                                     jetsAK8      = cms.InputTag("disambiguatedJetsAK8"),  # jets which do not contains leptons from ZZ or other good isolated leptons are removed
                                     Z            = cms.InputTag("selectedZCand"),
                                     Vhad         = cms.InputTag(""),
                                     ZZ           = cms.InputTag("ZZFiltered"),            # only the best ZZ->4l candidate that pass the FULL selection
                                     ZL           = cms.InputTag("ZlSelected"),
                                     ZW           = cms.InputTag("bareZWCand"),
                                     MET          = cms.InputTag("slimmedMETs"),
                                     Vertices     = cms.InputTag("goodPrimaryVertices"),                                    
                                     XSection     = cms.untracked.double(XSEC)
     )

### ------------------------------------------------------------------------- ###
### Triggers
### ------------------------------------------------------------------------- ###

triggerFilterForVVX = cms.EDFilter("VVXTriggerFilter",
                                   channelType  = cms.string("UNDEF"),
                                   isMC         = cms.untracked.bool(IsMC),
                                   setup        = cms.int32(LEPTON_SETUP),
                                   sampleType   = cms.int32(SAMPLE_TYPE),
                                   PD           = cms.string(PD),
                                   skimPaths    = cms.vstring(SkimPaths),
                                   MCFilterPath = cms.string(MCFILTER)
                               )

process.triggerForZZ = triggerFilterForVVX.clone()
process.triggerForZZ.channelType = cms.string("ZZ")

process.triggerForZW = triggerFilterForVVX.clone()
process.triggerForZW.channelType = cms.string("ZW")

process.triggerForZV = triggerFilterForVVX.clone()
process.triggerForZV.channelType = cms.string("ZV")

process.triggerForZL = triggerFilterForVVX.clone()
process.triggerForZL.channelType = cms.string("ZL")


########################################################
### Tools for further debuging                       ###
########################################################

process.printTree = cms.EDAnalyzer("ParticleListDrawer",
                                   maxEventsToPrint = cms.untracked.int32(-1),
                                   printVertex = cms.untracked.bool(False),
                                   src = cms.InputTag("prunedGenParticles")
                                   )


process.dumpUserData =  cms.EDAnalyzer("dumpData", #"dumpUserData",
    dumpTrigger = cms.untracked.bool(True),
    options = cms.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound')),
    muonSrcs =  cms.PSet(
        # slimmedMuons      = cms.InputTag("slimmedMuons"),
        # calibratedMuons   = cms.InputTag("calibratedMuons"),
        # muons             = cms.InputTag("appendPhotons:muons"),
        # postCleaningMuons = cms.InputTag("postCleaningMuons")
      ),
    electronSrcs = cms.PSet(
        # slimmedElectron        = cms.InputTag("slimmedElectrons"),
        # calibratedPatElectrons = cms.InputTag("calibratedPatElectrons"),
        # electrons              = cms.InputTag("appendPhotons:electrons"),
        # postCleaningElectrons  = cms.InputTag("postCleaningElectrons")
    ),
    candidateSrcs = cms.PSet(
        # Z   = cms.InputTag("ZCand"),
        # ZZ  = cms.InputTag("ZZCand"),
        # ZLL = cms.InputTag("ZLLCand"),
        # ZL  = cms.InputTag("ZlCand")
    ),
    photonSrc = cms.InputTag("filteredPhotons"),
    # full5x5SigmaIEtaIEtaMap   = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"),
    # phoChargedIsolation       = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
    # phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
    # phoPhotonIsolation        = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
    jetSrc = cms.PSet(
    #     cleanJets         = cms.InputTag("cleanJets"),
    #     JetsWithLeptonsRemover = cms.InputTag("JetsWithLeptonsRemover"),
        disambiguatedJets = cms.InputTag("disambiguatedJets"),
    #     correctedJetsAK8  = cms.InputTag("correctedJetsAK8"),
        disambiguatedJetsAK8 = cms.InputTag("disambiguatedJetsAK8")
    )
    # jetSrc = cms.InputTag("disambiguatedJetsAK8")
)


### ------------------------------------------------------------------------- ###
### Run the TreePlanter
### ------------------------------------------------------------------------- ###

process.filltrees = cms.EndPath(cms.ignore(process.triggerForZZ) + cms.ignore(process.triggerForZW) + cms.ignore(process.triggerForZV) + cms.ignore(process.triggerForZL) +
                                # process.dumpUserData +
                                process.treePlanter,
                                process.counters)


#process.filltrees = cms.EndPath(cms.ignore(process.zzTrigger) + process.srCounter + process.cr2P2FCounter + process.cr3P1FCounter + process.treePlanter + process.dumpUserData)

#process.filltrees = cms.Path(process.preselection * process.genCategory * process.treePlanter * process.printTree)
#process.filltrees = cms.EndPath(process.treePlanter *process.dumpUserData)

########################################################

