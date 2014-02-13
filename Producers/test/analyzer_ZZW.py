### ----------------------------------------------------------------------
### Based on ZZ->4l strategy.
###----------------------------------------------------------------------


LEPTON_SETUP = 2012
PD = ""
MCFILTER = ""
ELECORRTYPE   = "Paper" # "None", "Moriond", or "Paper"
ELEREGRESSION = "Paper" # "None", "Moriond", "PaperNoComb", or "Paper" 
APPLYMUCORR = True

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





#samples = [('WZZJets','cmgtools_group','/WZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "")]
#samples = [
#    ('WZZJets', 'cmgtools_group', '/WZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, '', 0.01968)
    # ('ZZ4mu',         'cmgtools', '/ZZTo4mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',         'cmgTuple.*root', 4,  ""),
    # ('ZZ4e',          'cmgtools', '/ZZTo4e_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                     'cmgTuple.*root', 4,  ""),
    # ('ZZ2mu2tau',     'cmgtools', '/ZZTo2mu2tau_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                'cmgTuple.*root', 4, ""),
    # ('ZZ2e2tau',      'cmgtools', '/ZZTo2e2tau_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                 'cmgTuple.*root', 4, ""),
    # ('ZZ2e2mu',       'cmgtools', '/ZZTo2e2mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                  'cmgTuple.*root', 4, ""),
    # ('ZZ4tau',        'cmgtools', '/ZZTo4tau_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                   'cmgTuple.*root', 4, ""),           
    # ('ZZ4mu_ext',     'cmgtools', '/ZZTo4mu_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',                'cmgTuple.*root', 4,  ""),
    # ('ZZ4e_ext',      'cmgtools', '/ZZTo4e_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',                 'cmgTuple.*root', 6,  ""),
    # ('ZZ2mu2tau_ext', 'cmgtools', '/ZZTo2mu2tau_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',            'cmgTuple.*root', 20, ""),
    # ('ZZ2e2tau_ext',  'cmgtools', '/ZZTo2e2tau_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',             'cmgTuple.*root', 20, ""),
    # ('ZZ2e2mu_ext',   'cmgtools', '/ZZTo2e2mu_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',              'cmgTuple.*root', 12, ""),
    # ('ZZ4tau_ext',    'cmgtools', '/ZZTo4tau_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',               'cmgTuple.*root', 40, ""),
    # ('VBFH126','cmgtools','/VBF_HToZZTo4L_M-126_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',  'cmgTuple.*root', 15, ""),
    # ('ttH126','cmgtools_group','/TTbarH_HToZZTo4L_M-126_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 1, ""),
    # ('WH126','cmgtools_group','/WH_HToZZTo4L_M-126_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',      'cmgTuple.*root', 2, ""),
    # ('ZH126','cmgtools_group','/ZH_HToZZTo4L_M-126_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',      'cmgTuple.*root', 2, ""),
    # ('DYJetsToLLTuneZ2M50-NoB','cmgtools','/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 4, "DYJets_NoB"),
    # ('DYJetsToLLTuneZ2M10-NoB','cmgtools','/DYJetsToLL_M-10To50filter_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 4, "DYJets_NoB"),
    # ('DYJetsToLLTuneZ2M50-B','cmgtools','/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 4, "DYJets_B"),
    # ('DYJetsToLLTuneZ2M10-B','cmgtools','/DYJetsToLL_M-10To50filter_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 4, "DYJets_B"),
    # ('TTTo2L2Nu2B','cmgtools','/TTTo2L2Nu2B_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 4, ""), 
    # ('ZZJetsTo4L','cmgtools','/ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 4, ""),
    # ('ZZZJets','cmgtools_group','/ZZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    # ('WZ','cmgtools_group','/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 4, ""),
    # ('TTZJets','cmgtools_group','/TTZJets_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, "")
    # ('powheg15jhuGenV3H126','cmgtools','/SMHiggsToZZTo4L_M-126_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v2/AODSIM/PAT_CMG_V5_15_0/','cmgTuple.*root', 2, ""),
    # ('minloH126', 'cmgtools', '/GluGluToHToZZTo4L_M-126_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, "MCAllEvents,MH126"),
    # ('WZZ_aMCatNLO','cmgtools_group','/WZZ_8TeV-aMCatNLO-herwig/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    # ('WWJets','cmgtools_group','/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0','cmgTuple.*root', 8, "")
#                    ]


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
  ### Replace parameters
### ----------------------------------------------------------------------
# process.source.fileNames = cms.untracked.vstring(
#     #'root://lxcms00//data3/2013/HZZ_cmgTuple/BE539_H1258TeV.root' #533 V5_15_0 version
#     '/store/cmst3/group/cmgtools/CMG/WZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_10_1_nLP.root'
#     #'/store/cmst3/group/cmgtools/CMG/WZZ_8TeV-aMCatNLO-herwig/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_100_1_GEb.root'
#     #'/store/cmst3/group/cmgtools/CMG/ZZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_10_1_UV1.root'
#     )


process.maxEvents.input = -1

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


### ----------------------------------------------------------------------
### Analyzer for Plots
### ----------------------------------------------------------------------


# jet-jet pairs
process.disambiguatedJets = cms.EDProducer("JetsWithLeptonsRemover",
                                           Preselection = cms.string("pt > 25 && abs(eta) < 2.5"),
                                           Jets  = cms.InputTag("cmgPFJetSel"),
                                           Muons = cms.InputTag("appendPhotons:muons"),
                                           Electrons = cms.InputTag("appendPhotons:electrons"),
                                           EnergyFractionAllowed = cms.double(0)) # maximum energy fraction carried by the lepton in the jet, to accept a jet as non from lepton

process.bareWCand = cms.EDProducer("CandViewShallowCloneCombiner",
                                   decay = cms.string('disambiguatedJets disambiguatedJets'),
                                   cut = cms.string('mass > 0'), # protect against ghosts
                                   checkCharge = cms.bool(True))

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
process.WjjSequence = cms.Sequence(process.disambiguatedJets * process.bareWCand * process.WCand)#  * process.wPt)


process.Candidates = cms.Path(process.muons             +
                              process.electrons         + process.cleanSoftElectrons +
                              process.appendPhotons     +
                              process.softLeptons       +
                              # Build Boson candidates
                              process.bareEECand        + process.EECand +  
                              process.bareMMCand        + process.MMCand +
                              process.WjjSequence
                              )

# Fill the tree for the analysis
process.treePlanter = cms.EDAnalyzer("TreePlanter",
                                     setup      = cms.int32(LEPTON_SETUP),
                                     sampleType = cms.int32(SAMPLE_TYPE),
                                     muons      = cms.InputTag("appendPhotons:muons"),
                                     electrons  = cms.InputTag("appendPhotons:electrons"),
                                     jets       = cms.InputTag("disambiguatedJets"),
                                     Zmm        = cms.InputTag("MMCand"),
                                     Zee        = cms.InputTag("EECand"),
                                     Wjj        = cms.InputTag("WCand"),
                                     MET        = cms.InputTag("cmgPFMET"),
                                     Vertices   = cms.InputTag("goodPrimaryVertices"),                                    
                                     isMC       = cms.untracked.bool(True),
                                     XSection   = cms.untracked.double(XSEC)
                                     )


process.genCategory =  cms.EDFilter("GenFilterCategory",
                                    src = cms.InputTag("genParticlesPruned"),
                                    Category = cms.int32(-1),
                                    SignalDefinition = cms.int32(3))

#process.filltrees = cms.Path(process.printTree + process.genCategory * process.treePlanter)
process.filltrees = cms.Path(process.genCategory * process.treePlanter)


### ----------------------------------------------------------------------
### Output root file (monitoring histograms)
### ----------------------------------------------------------------------
process.TFileService=cms.Service('TFileService',
                                fileName=cms.string('ZZWAnalysis.root')
                                )
