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

#MCFILTER  = "signaldefinition"

SIGNALDEFINITION = int('1',2)  # -1 means get everything, 1 means the request of having a ZZ pair with the  mass in the choosedn windows. For other topology see the README under VVXAnalysis/Commons.

CONTROLREGION = '2P2F'

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
    #'/store/cmst3/user/cmgtools/CMG/ZZTo2e2mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_100_1_irQ.root'
    #'/store/cmst3/user/cmgtools/CMG//DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_997_0_MQw.root',
    #'/store/cmst3/user/cmgtools/CMG//DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_998_0_V7X.root',
    #'/store/cmst3/user/cmgtools/CMG//DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_999_0_dXJ.root',
    #'/store/cmst3/user/cmgtools/CMG//DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_99_1_z8J.root',
    #'/store/cmst3/user/cmgtools/CMG//DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_9_1_blv.root'
    
    '/store/cmst3/user/cmgtools/CMG/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_1006_0_kDH.root',
    '/store/cmst3/user/cmgtools/CMG/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_1015_0_zrA.root',
    '/store/cmst3/user/cmgtools/CMG/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_1021_0_2Im.root',
    '/store/cmst3/user/cmgtools/CMG/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_1023_0_nwv.root',
    '/store/cmst3/user/cmgtools/CMG/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_1025_0_h1F.root',
    '/store/cmst3/user/cmgtools/CMG/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_1026_0_J5E.root'
    )

#process.source.skipEvents = cms.untracked.uint32(13328)
#process.source.skipEvents = cms.untracked.uint32(12788)

process.maxEvents.input = -1

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

### ----------------------------------------------------------------------
### Output root file
### ----------------------------------------------------------------------


process.TFileService=cms.Service('TFileService', fileName=cms.string('ZZjjAnalysisCR{0:s}.root'.format(CONTROLREGION)))


#####################################################################################################################################################


### ------------------------------------------------------------------------- ###
### Fill the tree for the analysis
### ------------------------------------------------------------------------- ###

 
process.treePlanter = cms.EDAnalyzer("TreePlanter",
                                     sampleName   = cms.string(SAMPLENAME),
                                     setup        = cms.int32(LEPTON_SETUP),
                                     sampleType   = cms.int32(SAMPLE_TYPE),
                                     PD           = cms.string(PD),
                                     skimPaths    = cms.vstring(SkimPaths),
                                     MCFilterPath = cms.string(MCFILTER),
                                     isMC         = cms.untracked.bool(IsMC),
                                     signalDefinition = cms.int32(SIGNALDEFINITION),
                                     muons        = cms.InputTag("postCleaningMuons"),     # all good isolated muons BUT the ones coming from ZZ decay
                                     electrons    = cms.InputTag("postCleaningElectrons"), # all good isolated electrons BUT the ones coming from ZZ decay
                                     jets         = cms.InputTag("disambiguatedJets"),     # jets which contains leptons from ZZ or other good isolated leptons are removed
                                     Zmm          = cms.InputTag("MMCand"),
                                     Zee          = cms.InputTag("EECand"),
                                     Wjj          = cms.InputTag("WCand"),
                                     ZZ4m         = cms.InputTag(""),          # only the best ZZ->4mu candidate that pass the FULL selection
                                     ZZ4e         = cms.InputTag(""),          # only the best ZZ->4e candidate that pass the FULL selection
                                     ZZ2e2m       = cms.InputTag(""),          # only the best ZZ->2e2mu candidate that pass the FULL selection
                                     Zll          = cms.InputTag("ZLLFiltered"), 
                                     MET          = cms.InputTag("cmgPFMET"),
                                     Vertices     = cms.InputTag("goodPrimaryVertices"),                                    
                                     XSection     = cms.untracked.double(XSEC)
                                     )


### ------------------------------------------------------------------------- ###
### Activate some skimming and cleaning of the collections
### ------------------------------------------------------------------------- ###





### ......................................................................... ###
### Build collections of muons and electrons that pass a quality criteria (isGood + isolation) and that are NOT selected to form the ZZ best candidate that pass the full selection
### ......................................................................... ###



# Muons cleaning. First, create a muon collection from the best ZZ candidate grand daughters
process.muonsFromZZ = cms.EDProducer("PATMuonsFromCompositeCandidates", src =  cms.InputTag("ZLLFiltered"), SplitLevel = cms.int32(1))

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
process.electronsFromZZ = cms.EDProducer("PATElectronsFromCompositeCandidates", src =  cms.InputTag("ZLLFiltered"), SplitLevel = cms.int32(1))

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
                                           JetPreselection     = cms.string("pt > 20"),
                                           DiBosonPreselection = cms.string(""),
                                           MatchingType        = cms.string("byDeltaR"), 
                                           Jets      = cms.InputTag("cmgPFJetSel"),
                                           Muons     = cms.InputTag("postCleaningMuons"),
                                           Electrons = cms.InputTag("postCleaningElectrons"),
                                           Diboson   = cms.InputTag("ZLLFiltered"),
                                           EnergyFractionAllowed = cms.double(0), # maximum energy fraction carried by the lepton in the jet, to accept a jet as non from lepton                             
                                           DebugPlots= cms.untracked.bool(False)
                                           )




### ......................................................................... ###
### Build the W->jj candidate out of the previously disambiguated jet collection, restricted to the central regions
### ......................................................................... ###

process.centralJets = cms.EDFilter("EtaPtMinCMGPFJetSelector", 
                                   src = cms.InputTag("disambiguatedJets"),
                                   ptMin   = cms.double(30),
                                   etaMin = cms.double(-2.4),
                                   etaMax = cms.double(2.4)
                                   )

process.bareWCand = cms.EDProducer("CandViewShallowCloneCombiner",
                                   decay = cms.string('centralJets centralJets'),
                                   cut = cms.string('mass > 20'), # protect against ghosts
                                   checkCharge = cms.bool(False))

process.WCand = cms.EDProducer("WCandidateFiller",
                               src = cms.InputTag("bareWCand"))

process.WjjSequence = cms.Sequence(process.centralJets * process.bareWCand * process.WCand)





### ......................................................................... ###
### Clean the 4 lepton candidates to select only the best possible candidate, requiring that passes the FULL selection
### ......................................................................... ###

# Z1MASS_LARGE is used to define the best Z1 candidate, to be aligned with what done for Higgs studies.
# Z1MASS is then used for the actual cut.

Z1MASS_LARGE  = "daughter(0).mass>40 && daughter(0).mass<120"
Z1MASS  = "daughter(0).mass>60 && daughter(0).mass<120"
Z2MASS  = "daughter(1).mass>60 && daughter(1).mass<120"
# Value here below are the ones used for the H->ZZ analysis and here for cross-check for dedicated studies. 
#Z1MASS  = "daughter(0).mass>40 && daughter(0).mass<120"
#Z2MASS  = "daughter(1).mass>12 && daughter(1).mass<120"

LL_OS = "(daughter(1).daughter(0).pdgId + daughter(1).daughter(1).pdgId) == 0" #Z2 = l+l-

PASSD0 = "(userFloat('d1.d0.isGood') && userFloat('d1.d0.combRelIsoPFFSRCorr') < 0.4)"
PASSD1 = "(userFloat('d1.d1.isGood') && userFloat('d1.d1.combRelIsoPFFSRCorr') < 0.4)"
FAILD0 = "!" + PASSD0
FAILD1 = "!" + PASSD1
BOTHFAIL = FAILD0 + "&&" + FAILD1
PASSD0_XOR_PASSD1 = "((" + PASSD0 + "&&" + FAILD1 + ") || (" + PASSD1 + "&&" + FAILD0 + "))"
PASSD0_OR_PASSD1  = "(" + PASSD0 + "||" + PASSD1 + ")"

# CR 3P1F
BESTZLL_3P1F   = (CR_BESTCANDBASE       + "&&" +  
                  Z1MASS_LARGE          + "&&" +  
                  "userFloat('d0.isBestZ') &&" + 
                  LL_OS                 + "&&" +  
                  PASSD0_OR_PASSD1      + "&&" +  
                  Z2SIP)

BESTZ2mLL_3P1F = (CR_BESTCANDBASE         + "&&" +  
                  Z1MASS_LARGE            + "&&" +  
                  "userFloat('d0.isBestZmm') &&" + 
                  LL_OS                   + "&&" +  
                  PASSD0_OR_PASSD1        + "&&" +  
                  Z2SIP)

BESTZ2eLL_3P1F = (CR_BESTCANDBASE         + "&&" +  
                  Z1MASS_LARGE            + "&&" +  
                  "userFloat('d0.isBestZee') &&" + 
                  LL_OS                   + "&&" +  
                  PASSD0_OR_PASSD1        + "&&" +  
                  Z2SIP)


ZLLSEL_3P1F  = (CR_BASESEL + "&&" + PASSD0_XOR_PASSD1 + "&&" + Z1MASS + "&&" + Z2MASS)



# CR 2P2F
BESTZLL_2P2F   = (CR_BESTCANDBASE       + "&&" +  
                  Z1MASS_LARGE          + "&&" +  
                  "userFloat('d0.isBestZ') &&" + 
                  LL_OS                 + "&&" +  
                  Z2SIP)

BESTZ2mLL_2P2F = (CR_BESTCANDBASE         + "&&" +  
                  Z1MASS_LARGE            + "&&" +  
                  "userFloat('d0.isBestZmm') &&" + 
                  LL_OS                   + "&&" +  
                  Z2SIP)

BESTZ2eLL_2P2F = (CR_BESTCANDBASE         + "&&" +  
                  Z1MASS_LARGE            + "&&" +  
                  "userFloat('d0.isBestZee') &&" + 
                  LL_OS                   + "&&" +  
                  Z2SIP)


ZLLSEL_2P2F  = (CR_BASESEL + "&&" + BOTHFAIL + "&&" + Z1MASS + "&&" + Z2MASS)



BESTZ2mLL = BESTZ2mLL_3P1F
BESTZ2eLL = BESTZ2eLL_3P1F
BESTZLL   = BESTZLL_3P1F
ZLLSEL    = ZLLSEL_3P1F


if CONTROLREGION == '2P2F':
    BESTZ2mLL = BESTZ2mLL_2P2F
    BESTZ2eLL = BESTZ2eLL_2P2F
    BESTZLL = BESTZLL_2P2F
    ZLLSEL  = ZLLSEL_2P2F
elif not CONTROLREGION == '2P2F' and not CONTROLREGION == '3P1F' :
    print "Do not know what tho do with {0:s} control region. Collapsing into 3P1F one.".format(CONTROLREGION)
    

process.ZLLCand.bestCandAmong.isBestCandCRZ2mLL = cms.string(BESTZ2mLL)
process.ZLLCand.bestCandAmong.isBestCandCRZ2eLL = cms.string(BESTZ2eLL)
process.ZLLCand.bestCandAmong.isBestCandZLL = cms.string(BESTZLL)
process.ZLLCand.flags.SelZLL = cms.string(ZLLSEL)



### ......................................................................... ###
### Clean the Z+2 lepton candidates to select only the best possible candidate for that region, requiring that at least one lepton fails the FULL selection
### ......................................................................... ###

process.Z2mLLFiltered = cms.EDFilter("PATCompositeCandidateSelector",
                                     src = cms.InputTag("ZLLCand"),
                                     cut = cms.string("userFloat('isBestCandCRZ2mLL') && userFloat('SelZLL')")
                                     )

process.Z2eLLFiltered = cms.EDFilter("PATCompositeCandidateSelector",
                                     src = cms.InputTag("ZLLCand"),
                                     cut = cms.string("userFloat('isBestCandCRZ2eLL') && userFloat('SelZLL')")
                                     )

#process.ZLLFiltered = cms.EDFilter("PATCompositeCandidateSelector",
#                                   src = cms.InputTag("ZLLCand"),
#                                   cut = cms.string("userFloat('isBestCandZLL') && userFloat('SelZLL')")
#                                   )

# Merger of all ZZ final states.
process.ZLLFiltered = cms.EDProducer("PATCompositeCandidateMerger",
                                     src = cms.VInputTag(cms.InputTag("Z2mLLFiltered"), cms.InputTag("Z2eLLFiltered"))
                                     )



### ------------------------------------------------------------------------- ###
### Define the post reconstruction cleaning sequence
### ------------------------------------------------------------------------- ###

process.postRecoCleaning = cms.Sequence( process.Z2mLLFiltered
                                         + process.Z2eLLFiltered
                                         + process.ZLLFiltered
                                         + process.muonsFromZZ*process.postCleaningMuons 
                                         + process.electronsFromZZ*process.postCleaningElectrons
                                         + process.disambiguatedJets
                                         + process.WjjSequence
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
process.zllCounterFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("ZLLFiltered"), minNumber = cms.uint32(1))

# Empty sequence to attach the signal filter (if specified in the CSV file)
process.signalFilters    = cms.Sequence() 


# Some counters and paths functional to the fake lepton background estimation

process.cand2P2F =  cms.EDFilter("PATCompositeCandidateSelector",
                                 src = cms.InputTag("ZLLFiltered"),
                                 cut = cms.string(BOTHFAIL)
                                 )

process.cand3P1F =  cms.EDFilter("PATCompositeCandidateSelector",
                                 src = cms.InputTag("ZLLFiltered"),
                                 cut = cms.string(PASSD0_XOR_PASSD1)
                                 )

process.cand2P2FFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand2P2F"), minNumber = cms.uint32(1))
process.cand3P1FFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand3P1F"), minNumber = cms.uint32(1))


#process.cr2P2FCounter      = cms.EDProducer("SelectedEventCountProducer", name = cms.vstring("cr2P2F","preselection"))
#process.cr3P1FCounter      = cms.EDProducer("SelectedEventCountProducer", name = cms.vstring("cr3P1F","preselection"))
process.cr2P2FCounter      = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("cr2P2F","preselection"))
process.cr3P1FCounter      = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("cr3P1F","preselection"))


#process.crDefinitions = cms.Sequence(cms.ignore(process.cand2P2F) * cms.ignore(process.cand2P2FFilter) + 
#                                     cms.ignore(process.cand3P1F) * cms.ignore(process.cand3P1FFilter)) 


### Path that pre-select the higher level objects that will input the TreePlanter
process.preselection = cms.Path( process.prePreselectionCounter
                                 * process.CR
                                 * process.signalFilters
                                 * process.postRecoCleaning 
                                 * process.zllCounterFilter * process.jetCounterFilter
#                                 * process.crDefinitions
                                 * process.postPreselectionCounter)


process.cr2P2F = cms.Path(process.cand2P2F * process.cand2P2FFilter)
process.cr3P1F = cms.Path(process.cand3P1F * process.cand3P1FFilter)



### If it is MC, run also the signal definition path
if IsMC:
    process.genCategory =  cms.EDFilter("ZZGenFilterCategory",
                                        Topology = cms.int32(SIGNALDEFINITION),
                                        ParticleStatus = cms.int32(1), 
                                        src            = cms.InputTag("genParticlesPruned"),
                                        GenJets        = cms.InputTag("genJetSel"),
                                        )
    process.signalCounter    = cms.EDProducer("EventCountProducer")
    process.signalDefinition = cms.Path(process.genCategory * process.signalCounter)


### ------------------------------------------------------------------------- ###
### Run the TreePlanter
### ------------------------------------------------------------------------- ###

process.filltrees = cms.EndPath(process.cr2P2FCounter + process.cr3P1FCounter + process.treePlanter)

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
                                                                 MMMM  = cms.InputTag("MMMMCand"),
                                                                 EEEE  = cms.InputTag("EEEECand"),
                                                                 #EEMM  = cms.InputTag("EEMMCand"),
                                                                 EEMM  = cms.InputTag("ZLLFiltered"),
                                                                 )
                                       )

#process.filltrees = cms.Path(process.preselection * process.genCategory * process.treePlanter * process.printTree)
#process.filltrees = cms.EndPath(process.cr2P2FCounter + process.cr3P1FCounter + process.treePlanter + process.dumpUserData)

########################################################
