### ----------------------------------------------------------------------
### Based on ZZ->4l strategy.
###----------------------------------------------------------------------


LEPTON_SETUP = 2012
JET_SETUP    = 2012
#PD = ""
#MCFILTER = ""
ELECORRTYPE   = "Paper" # "None", "Moriond", or "Paper"
ELEREGRESSION = "Paper" # "None", "Moriond", "PaperNoComb", or "Paper" 
APPLYMUCORR = True

SIGNALDEFINITION = int('1',2)  # -1 means get everything, 1 means the request of having a ZZ pair with the  mass in the choosedn windows. For other topology see the README under VVXAnalysis/Commons.

try:
    IsMC
except NameError:
    IsMC = True

try:
    LEPTON_SETUP
except NameError:
    LEPTON_SETUP = 2012 # define the set of effective areas, rho corrections, etc.

try:
    JET_SETUP
except NameError:
    JET_SETUP = 2012 # define the MVA for the jet PU Id

try:
    PD
except NameError:
    PD = ""             # "" for MC, "DoubleEle", "DoubleMu", or "MuEG" for data 

try:
    MCFILTER
except NameError:
    MCFILTER = ""

try:
    XSEC
except NameError:
    XSEC = -1

try:
    FULL_NOCUTS
except NameError:
    FULL_NOCUTS = False



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
SkimPaths.append("preselection")

### ----------------------------------------------------------------------
  ### Replace parameters
### ----------------------------------------------------------------------
process.source.fileNames = cms.untracked.vstring(
    # '/store/cmst3/group/cmgtools/CMG/WZZ_8TeV-aMCatNLO-herwig/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_13_1_YZF.root'    
    # 'root://lxcms00//data3/2013/HZZ_cmgTuple/BE539_H1258TeV.root' #533 V5_15_0 version
    #'/store/cmst3/group/cmgtools/CMG/WZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_10_1_nLP.root'
    #'/store/cmst3/group/cmgtools/CMG/WZZ_8TeV-aMCatNLO-herwig/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_100_1_GEb.root'
    #'/store/cmst3/group/cmgtools/CMG/ZZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_10_1_UV1.root'
    #'/store/cmst3/user/cmgtools/CMG/ZZTo2e2mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_100_1_irQ.root'
    #'/store/cmst3/user/cmgtools/CMG/ZZTo4mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_100_1_UR6.root'
    #'/store/cmst3/user/cmgtools/CMG/ZZTo2e2tau_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_10_1_E2X.root'
    #'/store/cmst3/user/cmgtools/CMG/ZZTo2e2muJJ_SMHContinInterf_M-125p6_8TeV-phantom-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_100_1_tZF.root'
    #'/store/cmst3/user/cmgtools/CMG/ZZTo2e2tau_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_50_1_Fka.root',
    #'/store/cmst3/user/cmgtools/CMG/ZZTo2e2tau_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_6_1_EfC.root'
    '/store/cmst3/user/cmgtools/CMG//ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_6_1_zVE.root'
    #'file:/tmp/bellan/cmgTuple_6_1_zVE.root'

    # '/store/cmst3/user/cmgtools/CMG/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_1006_0_kDH.root',
    # '/store/cmst3/user/cmgtools/CMG/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_1015_0_zrA.root',
    # '/store/cmst3/user/cmgtools/CMG/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_1021_0_2Im.root',
    # '/store/cmst3/user/cmgtools/CMG/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_1023_0_nwv.root',
    # '/store/cmst3/user/cmgtools/CMG/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_1025_0_h1F.root',
    # '/store/cmst3/user/cmgtools/CMG/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_1026_0_J5E.root'

    )

process.maxEvents.input = 1000

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

### ----------------------------------------------------------------------
### Output root file
### ----------------------------------------------------------------------
#process.source.eventsToProcess = cms.untracked.VEventRange("1:5666496")

process.TFileService=cms.Service('TFileService', fileName=cms.string('ZZjjAnalysis.root'))


#####################################################################################################################################################


### ------------------------------------------------------------------------- ###
### Activate some skimming and cleaning of the collections
### ------------------------------------------------------------------------- ###





### ......................................................................... ###
### Build collections of muons and electrons that pass a quality criteria (isGood + isolation) and that are NOT selected to form the ZZ best candidate that pass the full selection
### ......................................................................... ###



# Muons cleaning. First, create a muon collection from the best ZZ candidate grand daughters
process.muonsFromZZ = cms.EDProducer("PATMuonsFromCompositeCandidates", src =  cms.InputTag("ZZFiltered"), SplitLevel = cms.int32(1))

# Muons cleaning. Second, remove from the muon collection the muons that come from the best ZZ candidate that pass the full selection (previous collection).
# The newly produced collection is also filtered in muon quality and isolation
process.postCleaningMuons = cms.EDProducer("PATMuonCleaner",
                                           # pat electron input source
                                           src = cms.InputTag("appendPhotons:muons"),
                                           # preselection (any string-based cut for pat::Muons)
                                           preselection = cms.string("pt > 10 && userFloat('isGood') && userFloat('CombRelIsoPF') < 0.4"),
                                           # overlap checking configurables
                                           checkOverlaps = cms.PSet(
        muons = cms.PSet(
            src       = cms.InputTag("muonsFromZZ"), # Start from loose lepton def
            algorithm = cms.string("byDeltaR"),
            preselection        = cms.string(""), 
            deltaR              = cms.double(0.05),  
            checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
            pairCut             = cms.string(""),
            requireNoOverlaps   = cms.bool(True), # overlaps don't cause the electron to be discared
            )
        ),
                                           # finalCut (any string-based cut for pat::Muons)
                                           finalCut = cms.string(''),
                                           )





# Electrons cleaning. First, create a electron collection from the best ZZ candidate grand daughters
process.electronsFromZZ = cms.EDProducer("PATElectronsFromCompositeCandidates", src =  cms.InputTag("ZZFiltered"), SplitLevel = cms.int32(1))

# Electrons cleaning. Second, remove from the electron collection the electrons that come from the best ZZ candidate that pass the full selection (previous collection).
# The newly produced collection is also filtered in electron quality and isolation
process.postCleaningElectrons = cms.EDProducer("PATElectronCleaner",
                                               # pat electron input source
                                               src = cms.InputTag("appendPhotons:electrons"),
                                               # preselection (any string-based cut for pat::Electron)
                                               preselection = cms.string("pt > 10 && userFloat('isGood') && userFloat('CombRelIsoPF') < 0.4"),
                                               # overlap checking configurables
                                               checkOverlaps = cms.PSet(
        electrons = cms.PSet(
            src       = cms.InputTag("electronsFromZZ"), # Start from loose lepton def
            algorithm = cms.string("byDeltaR"),
            preselection        = cms.string(""), #
            deltaR              = cms.double(0.05),  
            checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
            pairCut             = cms.string(""),
            requireNoOverlaps   = cms.bool(True), # overlaps don't cause the electron to be discared
            )
        ),
                                               # finalCut (any string-based cut for pat::Electron)
                                               finalCut = cms.string(''),
                                               )




### ......................................................................... ###
# Remove from the event the jets that have leptons from the ZZ best candidate.
# Jets are also checked against other good isolated leptons not coming from the ZZ best candidate. 
# The jets, to be stored in the event, must pass the preselction (specified below by the user) AND the looseID + PU veto that is implemented in the code (same algo as for H->ZZ VBF selection)
### ......................................................................... ###

process.disambiguatedJets = cms.EDProducer("JetsWithLeptonsRemover",
                                           Setup               = cms.int32(JET_SETUP),
                                           JetPreselection     = cms.string("pt > 10"),
                                           DiBosonPreselection = cms.string(""),
                                           MatchingType        = cms.string("byDeltaR"), 
                                           Jets      = cms.InputTag("cmgPFJetSel"),
                                           Muons     = cms.InputTag("postCleaningMuons"),
                                           Electrons = cms.InputTag("postCleaningElectrons"),
                                           Diboson   = cms.InputTag("ZZFiltered"),
                                           EnergyFractionAllowed = cms.double(0), # maximum energy fraction carried by the lepton in the jet, to accept a jet as non from lepton                             
                                           DebugPlots= cms.untracked.bool(False)
                                           )




### ......................................................................... ###
### Build the W->jj candidate out of the previously disambiguated jet collection, restricted to the central regions
### ......................................................................... ###

process.centralJets = cms.EDFilter("EtaPtMinCMGPFJetSelector", 
                                   src = cms.InputTag("disambiguatedJets"),
                                   ptMin   = cms.double(10),
                                   etaMin = cms.double(-2.4),
                                   etaMax = cms.double(2.4)
                                   )

process.bareVhadCand = cms.EDProducer("CandViewShallowCloneCombiner",
                                      decay = cms.string('centralJets centralJets'),
                                      cut = cms.string('mass > 20'), # protect against ghosts
                                      checkCharge = cms.bool(False))

process.VhadCand = cms.EDProducer("WCandidateFiller",
                                  src = cms.InputTag("bareVhadCand"))

process.VhadSequence = cms.Sequence(process.centralJets * process.bareVhadCand * process.VhadCand)





### ......................................................................... ###
### Clean the 4 lepton candidates to select only the best possible candidate, requiring that passes the FULL selection
### ......................................................................... ###

Z1MASS            = "daughter('Z1').mass>60 && daughter('Z1').mass<120"
Z2MASS            = "daughter('Z2').mass>60 && daughter('Z2').mass<120"
FULLSELTIGHT      = (FULLSEL + "&&" + Z1MASS + "&&" + Z2MASS)

process.MMMMCand.flags.FullSelTight = cms.string(FULLSELTIGHT)
process.EEEECand.flags.FullSelTight = cms.string(FULLSELTIGHT)
process.EEMMCand.flags.FullSelTight = cms.string(FULLSELTIGHT)


process.MMMMFiltered = cms.EDFilter("PATCompositeCandidateSelector",
                                    src = cms.InputTag("MMMMCand"),
                                    cut = cms.string("userFloat('isBestCand') && userFloat('FullSelTight')")
                                    )

process.EEEEFiltered = cms.EDFilter("PATCompositeCandidateSelector",
                                    src = cms.InputTag("EEEECand"),
                                    cut = cms.string("userFloat('isBestCand') && userFloat('FullSelTight')")
                                    )

process.EEMMFiltered = cms.EDFilter("PATCompositeCandidateSelector",
                                    src = cms.InputTag("EEMMCand"),
                                    cut = cms.string("userFloat('isBestCand') && userFloat('FullSelTight')")
                                    )

### ......................................................................... ###
### Clean the Z+2 lepton candidates to select only the best possible candidate for that region, requiring that at least one lepton fails the FULL selection
### ......................................................................... ###


# Z1MASS_LARGE is used to define the best Z1 candidate, to be aligned with what done for Higgs studies.
# Z1MASS is then used for the actual cut.

CR_Z1MASS_LARGE  = "daughter(0).mass>40 && daughter(0).mass<120"
CR_Z1MASS  = "daughter(0).mass>60 && daughter(0).mass<120"
CR_Z2MASS  = "daughter(1).mass>60 && daughter(1).mass<120"
# Value here below are the ones used for the H->ZZ analysis and here for cross-check for dedicated studies. 
#Z1MASS  = "daughter(0).mass>40 && daughter(0).mass<120"
#Z2MASS  = "daughter(1).mass>12 && daughter(1).mass<120"

Z2LL_OS = "(daughter(1).daughter(0).pdgId + daughter(1).daughter(1).pdgId) == 0" #Z2 = l+l-

PASSD0 = "(userFloat('d1.d0.isGood') && userFloat('d1.d0.combRelIsoPFFSRCorr') < 0.4)"
PASSD1 = "(userFloat('d1.d1.isGood') && userFloat('d1.d1.combRelIsoPFFSRCorr') < 0.4)"
FAILD0 = "!" + PASSD0
FAILD1 = "!" + PASSD1
BOTHPASS = PASSD0 + "&&" + PASSD1
BOTHFAIL = FAILD0 + "&&" + FAILD1
PASSD0_XOR_PASSD1 = "((" + PASSD0 + "&&" + FAILD1 + ") || (" + PASSD1 + "&&" + FAILD0 + "))"
PASSD0_OR_PASSD1  = "(" + PASSD0 + "||" + PASSD1 + ")"

# CR Base
CR_BESTZLLos = (CR_BESTCANDBASE       + "&&" +  
                CR_Z1MASS_LARGE       + "&&" +  
                "userFloat('d0.isBestZ') &&" + 
                Z2LL_OS               + "&&" +  
                Z2SIP)


# CR 3P1F 
CR_BESTZLLos_3P1F = (CR_BESTZLLos    + "&&" +  
                     PASSD0_OR_PASSD1)

CR_ZLLosSEL_3P1F  = (CR_BASESEL + "&&" + PASSD0_XOR_PASSD1 + "&&" + CR_Z1MASS + "&&" + CR_Z2MASS)



# CR 2P2F
CR_BESTZLLos_2P2F = (CR_BESTZLLos)

CR_ZLLosSEL_2P2F  = (CR_BASESEL + "&&" + BOTHFAIL + "&&" + CR_Z1MASS + "&&" + CR_Z2MASS)


process.ZLLCand.bestCandAmong.isBestCRZLLos_3P1F = cms.string(CR_BESTZLLos_3P1F)
process.ZLLCand.flags.SelZLL_3P1F = cms.string(CR_ZLLosSEL_3P1F)

process.ZLLCand.bestCandAmong.isBestCRZLLos_2P2F = cms.string(CR_BESTZLLos_2P2F)
process.ZLLCand.flags.SelZLL_2P2F = cms.string(CR_ZLLosSEL_2P2F)

process.ZLLFiltered3P1F = cms.EDFilter("PATCompositeCandidateSelector",
                                       src = cms.InputTag("ZLLCand"),
                                       cut = cms.string("userFloat('isBestCRZLLos_3P1F') && userFloat('SelZLL_3P1F')")
                                       )

process.ZLLFiltered2P2F = cms.EDFilter("PATCompositeCandidateSelector",
                                       src = cms.InputTag("ZLLCand"),
                                       cut = cms.string("userFloat('isBestCRZLLos_2P2F') && userFloat('SelZLL_2P2F')")
                                       )


# Merger of all ZZ final states.
process.ZZFiltered = cms.EDProducer("PATCompositeCandidateMergerWithPriority",
                                    src = cms.VInputTag(cms.InputTag("MMMMFiltered"), cms.InputTag("EEEEFiltered"), cms.InputTag("EEMMFiltered"),
                                                        cms.InputTag("ZLLFiltered2P2F"), cms.InputTag("ZLLFiltered3P1F")),
                                    priority = cms.vint32(1,1,1,0,0)
                                    )

### ------------------------------------------------------------------------- ###
### Define the post reconstruction cleaning sequence
### ------------------------------------------------------------------------- ###

process.postRecoCleaning = cms.Sequence( process.MMMMFiltered
                                         + process.EEEEFiltered
                                         + process.EEMMFiltered
                                         + process.ZLLFiltered2P2F
                                         + process.ZLLFiltered3P1F
                                         + process.ZZFiltered
                                         + process.muonsFromZZ*process.postCleaningMuons 
                                         + process.electronsFromZZ*process.postCleaningElectrons
                                         + process.disambiguatedJets
                                         + process.VhadSequence
                                         )



### ------------------------------------------------------------------------- ###
### Run the preselection
### ------------------------------------------------------------------------- ###

# Skim counters
process.prePreselectionCounter       = cms.EDProducer("EventCountProducer")
#process.preSkimSignalCounter         = cms.EDProducer("EventCountProducer")
process.postPreselectionCounter      = cms.EDProducer("EventCountProducer")


### Some filters

# Number of disambiguated jets
process.jetCounterFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("disambiguatedJets"), minNumber = cms.uint32(0))

# Select only events with one such candidate
process.zzCounterFilter  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("ZZFiltered"), minNumber = cms.uint32(1))

### Path that pre-select the higher level objects that will input the TreePlanter
process.preselection = cms.Path( process.prePreselectionCounter
                                 * process.CR
                                 * process.postRecoCleaning 
                                 * process.zzCounterFilter * process.jetCounterFilter
                                 * process.postPreselectionCounter)


# Some counters and paths functional to the fake lepton background estimation and signal region check

process.cand2P2F       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("ZZFiltered"), cut = cms.string(BOTHFAIL))
process.cand2P2FFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand2P2F"), minNumber = cms.uint32(1))
process.cr2P2F         = cms.Path(process.cand2P2F * process.cand2P2FFilter)
process.cr2P2FCounter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("cr2P2F","preselection"))

process.cand3P1F       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("ZZFiltered"), cut = cms.string(PASSD0_XOR_PASSD1))
process.cand3P1FFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand3P1F"), minNumber = cms.uint32(1))
process.cr3P1F         = cms.Path(process.cand3P1F * process.cand3P1FFilter)
process.cr3P1FCounter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("cr3P1F","preselection"))

process.candSR       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("ZZFiltered"), cut = cms.string(BOTHPASS))
process.candSRFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("candSR"), minNumber = cms.uint32(1))
process.sr           = cms.Path(process.candSR * process.candSRFilter)
process.srCounter    = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("sr","preselection"))


### If it is MC, run also the signal definition path
if IsMC:
    # Empty sequence to attach the signal filter (if specified in the CSV file)
    process.mcSelectionCounter = cms.EDProducer("EventCountProducer") # not really needeed... it is mainly an hack to get the path executed
    process.signalFilters = cms.Sequence(process.mcSelectionCounter) 
    process.mcSelection   = cms.Path(process.signalFilters)
    MCFILTER = "mcSelection"

    process.genCategory =  cms.EDFilter("ZZGenFilterCategory",
                                        Topology       = cms.int32(SIGNALDEFINITION), 
                                        ParticleStatus = cms.int32(1), 
                                        src            = cms.InputTag("genParticlesPruned"),
                                        GenJets        = cms.InputTag("genJetSel"),
                                        )
    process.signalCounter    = cms.EDProducer("EventCountProducer")
    process.signalDefinition = cms.Path(process.genCategory * process.signalCounter)



### ------------------------------------------------------------------------- ###
### Fill the tree for the analysis
### ------------------------------------------------------------------------- ###

 
process.treePlanter = cms.EDAnalyzer("TreePlanter",
                                     sampleName   = cms.string(SAMPLENAME),
                                     setup        = cms.int32(LEPTON_SETUP),
                                     sampleType   = cms.int32(SAMPLE_TYPE),
                                     PD           = cms.string(PD),
                                     skimPaths    = cms.vstring(SkimPaths),
                                     SkimRequired = cms.untracked.bool(FULL_NOCUTS),
                                     MCFilterPath = cms.string(MCFILTER),
                                     isMC         = cms.untracked.bool(IsMC),
                                     signalDefinition = cms.int32(SIGNALDEFINITION),
                                     muons        = cms.InputTag("postCleaningMuons"),     # all good isolated muons BUT the ones coming from ZZ decay
                                     electrons    = cms.InputTag("postCleaningElectrons"), # all good isolated electrons BUT the ones coming from ZZ decay
                                     jets         = cms.InputTag("disambiguatedJets"),     # jets which contains leptons from ZZ or other good isolated leptons are removed
                                     Vhad         = cms.InputTag("VhadCand"),
                                     ZZ           = cms.InputTag("ZZFiltered"),            # only the best ZZ candidate that pass the FULL selection
                                     ZL           = cms.InputTag("ZlCand"),    
                                     MET          = cms.InputTag("cmgPFMET"),
                                     Vertices     = cms.InputTag("goodPrimaryVertices"),                                    
                                     XSection     = cms.untracked.double(XSEC)
                                     )



### ------------------------------------------------------------------------- ###
### Run the TreePlanter
### ------------------------------------------------------------------------- ###

process.filltrees = cms.EndPath(process.srCounter + process.cr2P2FCounter + process.cr3P1FCounter + process.treePlanter)


########################################################################################################################################################################




########################################################
### Tools for further debuging                       ###
########################################################

process.printTree = cms.EDAnalyzer("ParticleListDrawer",
                                   maxEventsToPrint = cms.untracked.int32(-1),
                                   printVertex = cms.untracked.bool(False),
                                   src = cms.InputTag("genParticlesPruned")
                                   )

process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
                                       dumpTrigger = cms.untracked.bool(True),
                                       muonSrc = cms.InputTag("appendPhotons:muons"), 
                                       electronSrc = cms.InputTag("appendPhotons:electrons"),
                                       candidateSrcs = cms.PSet( Zmm   = cms.InputTag("MMCand"),
                                                                 Zee   = cms.InputTag("EECand"),
                                                                 #MMMM  = cms.InputTag("MMMMCand"),
                                                                 #EEEE  = cms.InputTag("EEEECand"),
                                                                 ZLL  = cms.InputTag("ZLLCand"),
                                                                 #ZLLCandPreFiltered = cms.InputTag("ZLLCandPreFiltered"),
                                                                 EEMM  = cms.InputTag("EEMMCand"),
                                                                 )
                                       )

#process.filltrees = cms.Path(process.preselection * process.genCategory * process.treePlanter * process.printTree)
#process.filltrees = cms.EndPath(process.srCounter + process.cr2P2FCounter + process.cr3P1FCounter +process.treePlanter *process.dumpUserData)
#process.filltrees = cms.EndPath(process.treePlanter *process.printTree)

########################################################
