### ----------------------------------------------------------------------
### Based on ZZ->4l strategy.
###----------------------------------------------------------------------

LEPTON_SETUP = 2015
#PD = ""
#MCFILTER = ""
ELECORRTYPE   = "None" # "None", "Moriond", or "Paper"
ELEREGRESSION = "None" # "None", "Moriond", "PaperNoComb", or "Paper" 
APPLYMUCORR = False # ??? FIXME

SIGNALDEFINITION = int('1',2)  # -1 means get everything, 1 means the request of having a ZZ pair with the  mass in the choosedn windows. For other topology see the README under VVXAnalysis/Commons.

try:
    IsMC
except NameError:
    IsMC = True

try:
    LEPTON_SETUP
except NameError:
    LEPTON_SETUP = 2015

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
    '/store/mc/Phys14DR/ZZTo4L_Tune4C_13TeV-powheg-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/04CD96C9-E269-E411-9D64-00266CF9ADA0.root'
    #'/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/3295EF7C-2070-E411-89C4-7845C4FC35DB.root'
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
                                           JetPreselection     = cms.string("pt > 30"),
                                           DiBosonPreselection = cms.string(""),
                                           MatchingType        = cms.string("byDeltaR"), 
                                           Jets      = cms.InputTag("slimmedJets"),
                                           Muons     = cms.InputTag("postCleaningMuons"),
                                           Electrons = cms.InputTag("postCleaningElectrons"),
                                           Diboson   = cms.InputTag("ZZFiltered"),
                                           EnergyFractionAllowed = cms.double(0), # maximum energy fraction carried by the lepton in the jet, to accept a jet as non from lepton                             
                                           DebugPlots= cms.untracked.bool(False)
                                           )




### ......................................................................... ###
### Build the W->jj candidate out of the previously disambiguated jet collection, restricted to the central regions
### ......................................................................... ###

process.centralJets = cms.EDFilter("EtaPtMinCandViewSelector", 
                                   src = cms.InputTag("disambiguatedJets"),
                                   ptMin   = cms.double(30),
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

process.ZZCand.flags.FullSelTight = cms.string(FULLSELTIGHT)

process.ZZFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
                                  src = cms.InputTag("ZZCand"),
                                  cut = cms.string("userFloat('isBestCand') && userFloat('FullSelTight')")
                                  )

### ......................................................................... ###
### Clean the Z+2 lepton candidates to select only the best possible candidate for that region, requiring that at least one lepton fails the FULL selection
### ......................................................................... ###

#process.ZLLFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
#                                   src = cms.InputTag("ZLLCand"),
#                                   cut = cms.string("userFloat('isBestCand') && userFloat('FullSel')")
#                                   )



### ------------------------------------------------------------------------- ###
### Define the post reconstruction cleaning sequence
### ------------------------------------------------------------------------- ###

process.postRecoCleaning = cms.Sequence( process.ZZFiltered
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
                                 * process.postRecoCleaning 
                                 * process.zzCounterFilter * process.jetCounterFilter
                                 * process.postPreselectionCounter)


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
                                        src            = cms.InputTag("prunedGenParticles"),
                                        GenJets        = cms.InputTag("slimmedGenJets"),
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
                                     #SkimRequired = cms.untracked.bool(False),
                                     MCFilterPath = cms.string(MCFILTER),
                                     isMC         = cms.untracked.bool(IsMC),
                                     signalDefinition = cms.int32(SIGNALDEFINITION),
                                     JECFileName  = cms.string("VVXAnalysis/Producers/test/Winter14_V5_DATA_Uncertainty_AK5PF.txt"),
                                     muons        = cms.InputTag("postCleaningMuons"),     # all good isolated muons BUT the ones coming from ZZ decay
                                     electrons    = cms.InputTag("postCleaningElectrons"), # all good isolated electrons BUT the ones coming from ZZ decay
                                     jets         = cms.InputTag("disambiguatedJets"),     # jets which do not contains leptons from ZZ or other good isolated leptons
                                     Z            = cms.InputTag("ZCand"),
                                     Vhad         = cms.InputTag("VhadCand"),
                                     ZZ           = cms.InputTag("ZZFiltered"),            # only the best ZZ->4l candidate that pass the FULL selection
                                     Zll          = cms.InputTag(""), 
                                     MET          = cms.InputTag("slimmedMETs"),
                                     Vertices     = cms.InputTag("goodPrimaryVertices"),                                    
                                     XSection     = cms.untracked.double(XSEC)
                                     )



### ------------------------------------------------------------------------- ###
### Run the TreePlanter
### ------------------------------------------------------------------------- ###


process.thisIsTheEnd = cms.EndPath(process.treePlanter)

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
#process.filltrees = cms.EndPath(process.treePlanter *process.printTree)

########################################################
