1### NOTE: This is prepared to run on the newest PDFs with LHAPDF >=3.8.4
### so it requires local installation of LHAPDF libraries in order to run 
### out of the box. Otherwise, substitute the PDF sets by older sets

import FWCore.ParameterSet.Config as cms


### ----------------------------------------------------------------------
### Based on ZZ->4l strategy.
###----------------------------------------------------------------------

LEPTON_SETUP = 2015
#PD = ""
#MCFILTER = ""
ELECORRTYPE   = "None" # "None", "Moriond", or "Paper"
ELEREGRESSION = "None" # "None", "Moriond", "PaperNoComb", or "Paper" 
APPLYMUCORR = False # ??? FIXME

SIGNALDEFINITION = int('1',2)  # -1 means get everything, 1 means the request of having a ZZ pair with the  mass in the chosen windows. For other topology see the README under VVXAnalysis/Commons.

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

try:
    SKIM_REQUIRED
except NameError:
    SKIM_REQUIRED = True



import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/VVXAnalysis/Producers/python/"
execfile(PyFilePath + "analyzer_ZZjj.py")
process.filltrees = cms.EndPath()

# Max events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(),
                            secondaryFileNames = cms.untracked.vstring()
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




# Printouts
process.MessageLogger = cms.Service("MessageLogger",
      cout = cms.untracked.PSet(
            default = cms.untracked.PSet(limit = cms.untracked.int32(100)),
            threshold = cms.untracked.string('INFO')
      ),
      destinations = cms.untracked.vstring('cout')
)
