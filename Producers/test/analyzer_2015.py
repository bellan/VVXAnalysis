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



# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/VVXAnalysis/Producers/python/"
execfile(PyFilePath + "analyzer_ZZjj.py")


### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------
SkimPaths.append("preselection")

### ----------------------------------------------------------------------
  ### Replace parameters
### ----------------------------------------------------------------------
process.source.fileNames = cms.untracked.vstring(
#    '/store/mc/RunIISpring15DR74/GluGluToZZTo4mu_BackgroundOnly_13TeV_MCFM/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/1CAF3EBF-9118-E511-8143-02163E0136C2.root'
    '/store/mc/RunIISpring15DR74/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/02B9AF23-1B1A-E511-937A-0025905B8576.root'

    #'/store/mc/Phys14DR/ZZTo4L_Tune4C_13TeV-powheg-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/04CD96C9-E269-E411-9D64-00266CF9ADA0.root'
    #'/store/mc/Phys14DR/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/3295EF7C-2070-E411-89C4-7845C4FC35DB.root'
    #'/store/mc/Spring14miniaod/WH_ZH_HToZZ_4LFilter_M-125_13TeV_pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/0ABFFCEF-AA09-E411-8022-001E673970FD.root'
    )

process.maxEvents.input = 100
