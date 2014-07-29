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
    # '/store/cmst3/group/cmgtools/CMG/WZZ_8TeV-aMCatNLO-herwig/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_13_1_YZF.root'    
    # 'root://lxcms00//data3/2013/HZZ_cmgTuple/BE539_H1258TeV.root' #533 V5_15_0 version
    #'/store/cmst3/group/cmgtools/CMG/WZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_10_1_nLP.root'
    #'/store/cmst3/group/cmgtools/CMG/WZZ_8TeV-aMCatNLO-herwig/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_100_1_GEb.root'
    # '/store/cmst3/group/cmgtools/CMG/ZZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_10_1_UV1.root'
    '/store/cmst3/user/cmgtools/CMG/ZZTo2e2mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_100_1_irQ.root'
    #'/store/cmst3/user/cmgtools/CMG//ZZTo4mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_100_1_UR6.root'
    #'/store/cmst3/user/cmgtools/CMG/VBF_phantom_8TeV/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_42.root'
    )


process.maxEvents.input = -1

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

### ----------------------------------------------------------------------
### Output root file
### ----------------------------------------------------------------------


process.TFileService=cms.Service('TFileService', fileName=cms.string('ZZjjAnalysis.root'))


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
                                     muons        = cms.InputTag("postCleaningMuons"),     # all good isolated muons BUT the ones coming from ZZ decay
                                     electrons    = cms.InputTag("postCleaningElectrons"), # all good isolated electrons BUT the ones coming from ZZ decay
                                     jets         = cms.InputTag("disambiguatedJets"),     # jets which contains leptons from ZZ or other good isolated leptons are removed
                                     Zmm          = cms.InputTag("MMCand"),
                                     Zee          = cms.InputTag("EECand"),
                                     Wjj          = cms.InputTag("WCand"),
                                     ZZ4m         = cms.InputTag("MMMMFiltered"),          # only the best ZZ->4mu candidate that pass the FULL selection
                                     ZZ4e         = cms.InputTag("EEEEFiltered"),          # only the best ZZ->4e candidate that pass the FULL selection
                                     ZZ2e2m       = cms.InputTag("EEMMFiltered"),          # only the best ZZ->2e2mu candidate that pass the FULL selection
                                     Zll          = cms.InputTag("ZLLCand"), 
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
                                           JetPreselection     = cms.string("pt > 20"),
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

process.ZZFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
                                  src = cms.InputTag("ZZCand"),
                                  cut = cms.string("userFloat('isBestCand') && userFloat('FullSel')")
                                  )

process.MMMMFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
                                    src = cms.InputTag("MMMMCand"),
                                    cut = cms.string("userFloat('isBestCand') && userFloat('FullSel')")
                                    )

process.EEEEFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
                                    src = cms.InputTag("EEEECand"),
                                    cut = cms.string("userFloat('isBestCand') && userFloat('FullSel')")
                                    )

process.EEMMFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
                                    src = cms.InputTag("EEMMCand"),
                                    cut = cms.string("userFloat('isBestCand') && userFloat('FullSel')")
                                    )

### ------------------------------------------------------------------------- ###





### ------------------------------------------------------------------------- ###
### Define the post reconstruction cleaning sequence
### ------------------------------------------------------------------------- ###

process.postRecoCleaning = cms.Sequence( process.ZZFiltered
                                         + process.MMMMFiltered
                                         + process.EEEEFiltered
                                         + process.EEMMFiltered
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
process.zzCounterFilter  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("ZZFiltered"), minNumber = cms.uint32(1))

# Empty sequence to attach the signal filter (if specified in the CSV file)
process.signalFilters    = cms.Sequence() 

### Path that pre-select the higher level objects that will input the TreePlanter
process.preselection = cms.Path( process.prePreselectionCounter
                                 * process.signalFilters
                                 * process.postRecoCleaning 
                                 * process.zzCounterFilter * process.jetCounterFilter
                                 * process.postPreselectionCounter)



### If it is MC, run also the signal definition path
if IsMC:
    process.genCategory =  cms.EDFilter("GenFilterCategory",
                                        src = cms.InputTag("genParticlesPruned"),
                                        Category = cms.int32(-1),
                                        SignalDefinition = cms.int32(3))
    
    process.signalDefinition = cms.Path(process.genCategory)



### Run the CR sequence too ###
process.CRPath = cms.Path(process.CR)



### ------------------------------------------------------------------------- ###
### Run the TreePlanter
### ------------------------------------------------------------------------- ###

process.filltrees = cms.EndPath(process.treePlanter)

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
                                                                 EEMM  = cms.InputTag("EEMMCand"),
                                                                 )
                                       )

#process.filltrees = cms.Path(process.preselection * process.genCategory * process.treePlanter * process.printTree)
#process.filltrees = cms.EndPath(process.treePlanter *process.dumpUserData)

########################################################
