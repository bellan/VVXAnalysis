### ----------------------------------------------------------------------
### Based on ZZ->4l strategy.
###----------------------------------------------------------------------

LEPTON_SETUP = 2015
#PD = ""
#MCFILTER = ""
ELECORRTYPE   = "None" # "None", "Moriond", or "Paper"
ELEREGRESSION = "None" # "None", "Moriond", "PaperNoComb", or "Paper" 
APPLYMUCORR = False # ??? FIXME

#FSRMODE = "RunII" #CHECK

SIGNALDEFINITION = int('1',2)  # -1 means get everything, 1 means the request of having a ZZ pair with the  mass in the chosen windows. For other topology see the README under VVXAnalysis/Commons.

try:
    IsMC
except NameError:
    IsMC = True
#    IsMC = False

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


#SELSETUP="allCutsAtOnceButMZ2"

# Get absolute path
import os

SELSETUP = "allCutsAtOnceZZFullRange" #Setup for ZZFull

PyFilePath = os.environ['CMSSW_BASE'] + "/src/VVXAnalysis/Producers/python/"
execfile(PyFilePath + "analyzer_ZZjj.py")

print "SELSETUP",SELSETUP
print "BESTCAND_AMONG",BESTCAND_AMONG,"\n"


#print "CR_BESTZLLos_3P1F",CR_BESTZLLos_3P1F
#print "CR_BESTZLLos_2P2F",CR_BESTZLLos_2P2F

# CR 3P1F
process.ZLLCand.bestCandAmong.isBestCRZLLos_3P1F = cms.string(CR_BESTZLLos_3P1F.replace("daughter(1).mass>12","daughter(1).mass>4"))

# CR 2P2F
process.ZLLCand.bestCandAmong.isBestCRZLLos_2P2F = cms.string(CR_BESTZLLos_2P2F.replace("daughter(1).mass>12","daughter(1).mass>4"))


#print "\nNew CR_BESTZLLos_3P1F",CR_BESTZLLos_3P1F.replace("daughter(1).mass>12","daughter(1).mass>4")


### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------
SkimPaths.append("preselection")

### ----------------------------------------------------------------------
  ### Replace parameters
### ----------------------------------------------------------------------
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(

#'/store/mc/RunIISpring15DR74/GluGluToZZTo4mu_BackgroundOnly_13TeV_MCFM/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/1CAF3EBF-9118-E511-8143-02163E0136C2.root'
#'/store/mc/RunIISpring15MiniAODv2/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/624BC18A-CA6D-E511-A73B-003048F3511E.root'
#'/store/mc/RunIISpring15MiniAODv2/GluGluToZZTo4e_BackgroundOnly_13TeV_MCFM/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/186EDBC4-DE72-E511-8720-0025905A60A6.root',
#'/store/mc/RunIIFall15MiniAODv1/ZZTo4L_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/064FDFEF-3DA7-E511-8340-001E6739850C.root'
#'/store/mc/RunIISpring15MiniAODv2/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v2/40000/00768931-1C76-E511-81DA-00266CFAEA48.root'

'/store/mc/RunIISpring15MiniAODv2/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v2/40000/76E4A530-1C76-E511-A951-00266CFAE764.root'


#'/store/mc/RunIISpring15MiniAODv2/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v3/60000/00181849-176A-E511-8B11-848F69FD4C94.root'

#'/store/data/Run2015D/DoubleMuon/MINIAOD/PromptReco-v4/000/258/159/00000/0C6D4AB0-6F6C-E511-8A64-02163E0133CD.root'

#'/store/mc/RunIISpring15MiniAODv2/GluGluToZZTo4e_BackgroundOnly_13TeV_MCFM/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/16D6F7C9-DE72-E511-B3BE-0026189438E2.root'

    #news2016
 #   '/store/mc/RunIIFall15MiniAODv2/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/024F30C4-96B9-E511-91A0-24BE05C616C1.root'

#'/store/mc/RunIIFall15MiniAODv1/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0C0945C1-86B2-E511-B553-0CC47A78A440.root'

    ))
    #),eventsToProcess = cms.untracked.VEventRange('1:31907:6167667','1:31907:6167734'))
process.maxEvents.input = -1