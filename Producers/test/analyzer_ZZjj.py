### ----------------------------------------------------------------------
### Based on ZZ->4l strategy.
###----------------------------------------------------------------------


LEPTON_SETUP = 2012
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
    '/store/cmst3/group/cmgtools/CMG/WZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_10_1_nLP.root'
    #'/store/cmst3/group/cmgtools/CMG/WZZ_8TeV-aMCatNLO-herwig/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_100_1_GEb.root'
    # '/store/cmst3/group/cmgtools/CMG/ZZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_10_1_UV1.root'
    # '/store/cmst3/user/cmgtools/CMG/ZZTo2e2mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_100_1_irQ.root'
    #'/store/cmst3/user/cmgtools/CMG//ZZTo4mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_100_1_UR6.root'
    #'/store/cmst3/user/cmgtools/CMG/VBF_phantom_8TeV/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_42.root'
    )


process.maxEvents.input = 1000

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000


### ----------------------------------------------------------------------
### Analyzer for Plots
### ----------------------------------------------------------------------


# jet-jet pairs
process.disambiguatedJets = cms.EDProducer("JetsWithLeptonsRemover",
                                           Preselection = cms.string("pt > 20"),
                                           Jets  = cms.InputTag("cmgPFJetSel"),
                                           Muons = cms.InputTag("appendPhotons:muons"),
                                           Electrons = cms.InputTag("appendPhotons:electrons"),
                                           EnergyFractionAllowed = cms.double(0)) # maximum energy fraction carried by the lepton in the jet, to accept a jet as non from lepton


process.centralJets = cms.EDFilter("EtaPtMinCMGPFJetSelector", 
                                   src = cms.InputTag("disambiguatedJets"),
                                   ptMin   = cms.double(20),
                                   etaMin = cms.double(-2.5),
                                   etaMax = cms.double(2.5)
                                   )


process.bareWCand = cms.EDProducer("CandViewShallowCloneCombiner",
                                   decay = cms.string('centralJets centralJets'),
                                   cut = cms.string('mass > 0'), # protect against ghosts
                                   checkCharge = cms.bool(False))

process.WCand = cms.EDProducer("WCandidateFiller",
                               src = cms.InputTag("bareWCand"))


### Simple plot for test

process.wPt= cms.EDAnalyzer("CandViewHistoAnalyzer",
                            src = cms.InputTag("WCand"),
                            histograms = cms.VPSet(cms.PSet(min = cms.untracked.double(0.0),
                                                            max = cms.untracked.double(1000.0),
                                                            nbins = cms.untracked.int32(100),   
                                                            name = cms.untracked.string("W_pT"),
                                                            description = cms.untracked.string("W pT [GeV/c]"),
                                                            plotquantity = cms.untracked.string("pt"))))

# Build W->jj candidate
process.WjjSequence = cms.Sequence(process.disambiguatedJets * process.centralJets * process.bareWCand * process.WCand)#  * process.wPt)

# Use the same sequence as in H->ZZ
process.Candidates += process.WjjSequence 



# Fill the tree for the analysis
process.treePlanter = cms.EDAnalyzer("TreePlanter",
                                     channel = cms.untracked.string('aChannel'),
                                     setup        = cms.int32(LEPTON_SETUP),
                                     sampleType   = cms.int32(SAMPLE_TYPE),
                                     PD           = cms.string(PD),
                                     skimPaths    = cms.vstring(SkimPaths),
                                     MCFilterPath = cms.string(MCFILTER),
                                     isMC         = cms.untracked.bool(IsMC),
                                     muons        = cms.InputTag("appendPhotons:muons"),
                                     electrons    = cms.InputTag("appendPhotons:electrons"),
                                     jets         = cms.InputTag("disambiguatedJets"),
                                     Zmm          = cms.InputTag("MMCand"),
                                     Zee          = cms.InputTag("EECand"),
                                     Wjj          = cms.InputTag("WCand"),
                                     ZZmmmm       = cms.InputTag("MMMMCand"),
                                     ZZeeee       = cms.InputTag("EEEECand"),
                                     ZZeemm       = cms.InputTag("EEMMCand"),
                                     MET          = cms.InputTag("cmgPFMET"),
                                     Vertices     = cms.InputTag("goodPrimaryVertices"),                                    
                                     XSection     = cms.untracked.double(XSEC)
                                     )


process.genCategory =  cms.EDFilter("GenFilterCategory",
                                    src = cms.InputTag("genParticlesPruned"),
                                    Category = cms.int32(-1),
                                    SignalDefinition = cms.int32(3))

### Activate some skimming ###
process.preSkimCounter = cms.EDProducer("EventCountProducer")

process.ZZFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
                                  src = cms.InputTag("ZZCand"),
                                  cut = cms.string("userFloat('isBestCand')")
                                  )

### Select only events with one such candidate
process.zzCounterFilter= cms.EDFilter("CandViewCountFilter",
                                      src = cms.InputTag("ZZFiltered"),
                                      minNumber = cms.uint32(1)
                                      )


process.jetPtFilter = cms.EDFilter("PtMinCandViewSelector",
                                   src = cms.InputTag("disambiguatedJets"),
                                   ptMin = cms.double(20)
                                   )
process.jetCounterFilter = cms.EDFilter("CandViewCountFilter",
                                        src = cms.InputTag("jetPtFilter"),
                                        minNumber = cms.uint32(1),
                                        )
process.postSkimCounter = cms.EDProducer("EventCountProducer")


process.preselection = cms.Path(process.preSkimCounter*
                                (process.ZZFiltered*process.zzCounterFilter+
                                process.jetPtFilter*process.jetCounterFilter)*
                                process.postSkimCounter)
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
     candidateSrcs = cms.PSet(
        Zmm   = cms.InputTag("MMCand"),
        Zee   = cms.InputTag("EECand"),
#        Z     = cms.InputTag("ZCand"),
        MMMM  = cms.InputTag("MMMMCand"),
        EEEE  = cms.InputTag("EEEECand"),
        EEMM  = cms.InputTag("EEMMCand"),
     )
)





process.signalDefinition = cms.Path(process.genCategory)
#process.filltrees = cms.Path(process.preselection * process.genCategory * process.treePlanter * process.printTree)

process.filltrees = cms.EndPath(process.treePlanter *process.dumpUserData)

### ----------------------------------------------------------------------
### Output root file (monitoring histograms)
### ----------------------------------------------------------------------
process.TFileService=cms.Service('TFileService',
                                fileName=cms.string('ZZWAnalysis.root')
                                )

