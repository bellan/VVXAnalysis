1### NOTE: This is prepared to run on the newest PDFs with LHAPDF >=3.8.4
### so it requires local installation of LHAPDF libraries in order to run 
### out of the box. Otherwise, substitute the PDF sets by older sets

import FWCore.ParameterSet.Config as cms


### ----------------------------------------------------------------------
### Based on ZZ->4l strategy.
###----------------------------------------------------------------------
LEPTON_SETUP = 2012
JET_SETUP = 2012
#PD = ""
#MCFILTER = ""
ELECORRTYPE = "Paper" # "None", "Moriond", or "Paper"
ELEREGRESSION = "Paper" # "None", "Moriond", "PaperNoComb", or "Paper"
APPLYMUCORR = True
#MCFILTER = "signaldefinition"
SIGNALDEFINITION = int('1',2) # -1 means get everything, 1 means the request of having a ZZ pair with the mass in the choosedn windows. For other topology see the README under VVXAnalysis/Commons.
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
    PD = "" # "" for MC, "DoubleEle", "DoubleMu", or "MuEG" for data

try:
    MCFILTER
except NameError:
    MCFILTER = ""

try:
    XSEC
except NameError:
    XSEC = -1


import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/VVXAnalysis/Producers/python/"
execfile(PyFilePath + "analyzer_ZZjj.py")
process.filltrees = cms.EndPath()

# Max events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

infile1='/store/cmst3/user/cmgtools/CMG/ZZTo4mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_100_1_UR6.root'
infile2= '/store/cmst3/user/cmgtools/CMG/ZZTo2e2mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_100_1_irQ.root'
infile3='/store/cmst3/group/cmgtools/CMG/ZZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_10_1_UV1.root'
# Input files (on disk)

weightfile='file:PdfWeight.root'

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(weightfile),
                            secondaryFileNames = cms.untracked.vstring(infile1)
                            )

process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.TFileService=cms.Service('TFileService', fileName=cms.string('PDFStudies.root'))


process.signalFilters += process.genCategory

# Collect uncertainties for rate and acceptance
process.pdfSystematics = cms.EDAnalyzer("PdfSystematicsAnalyzerZZ",
                                        isMC         = cms.untracked.bool(IsMC),
                                        setup        = cms.int32(LEPTON_SETUP),
                                        sampleType   = cms.int32(SAMPLE_TYPE),
                                        PD           = cms.string(PD),
                                        skimPaths    = cms.vstring(SkimPaths),
                                        MCFilterPath = cms.string(MCFILTER),
                                        FilterNames = cms.vstring('sr','preselection','zzTrigger'),
                                        PdfWeightTags = cms.untracked.VInputTag("pdfWeights:CT10"
                                                                                , "pdfWeights:MSTW2008nlo68cl"
                                                                                , "pdfWeights:NNPDF20"
                                                                                )
                                        )

process.genEventCounter  = cms.EDProducer("EventCountProducer")

process.theEnd = cms.EndPath(process.genEventCounter*cms.ignore(process.zzTrigger)*process.pdfSystematics)


process.MessageLogger = cms.Service("MessageLogger",
                                    cout = cms.untracked.PSet(
        default = cms.untracked.PSet(limit = cms.untracked.int32(100)),
        threshold = cms.untracked.string('INFO')
        ),
                                    destinations = cms.untracked.vstring('cout')
                                    )
