import ROOT
from ROOT import gSystem

# global MCSet
# MCSet = "Mad"

# def SetGlobalVariable(SetIn):

#     global MCSet
#     MCSet = SetIn

DataSamples = [{"sample":'data',"name":'Data'}]

BkgSamples  = [{"sample":'WWZJets',"name":'Irreducible background'},{"sample":'TTWWJets',"name":'Irreducible background'},{"sample":'TTZJets',"name":'Irreducible background'}]

SignalZZ_qq_Pow = [{"sample":'ZZTo2e2mu',"name":'qq/qg #rightarrow ZZ(+jets)'},{"sample":'ZZTo4e',"name":'qq/qg #rightarrow ZZ(+jets)'},{"sample":'ZZTo4mu',"name":'qq/qg #rightarrow ZZ(+jets)'}]

SignalZZ_VBS = [{"sample":'ZZTo2e2muJJ_SMHContinInterf_H125.6',"name":'qq/qg #rightarrow ZZ(+jets)'},{"sample":'ZZTo4eJJ_SMHContinInterf_H125.6',"name":'qq/qg #rightarrow ZZ(+jets)'},{"sample":'ZZTo4muJJ_SMHContinInterf_H125.6',"name":'qq/qg #rightarrow ZZ(+jets)'}]

SignalZZ_gg=[{"sample":'ggTo2e2mu_SMHContinInterf-MCFM67_H125.6',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)'},{"sample":'ggTo4e_SMHContinInterf-MCFM67_H125.6',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)'},{"sample":'ggTo4mu_SMHContinInterf-MCFM67_H125.6',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)'}]

#Signal_Others = [{"sample":'ZZZJets',"color":ROOT.kAzure-6,"name":'Other ZZ processes'},{"sample":'WZZJets',"color":ROOT.kAzure-6,"name":'Other ZZ processes'},{"sample":'WH126',"color":ROOT.kAzure-6,"name":'Other ZZ processes'},{"sample":'ZH126',"color":ROOT.kAzure-6,"name":'Other ZZ processes'},{"sample":'ttH126',"color":ROOT.kAzure-6,"name":'Other ZZ processes'}] 

Signal_Others = [{"sample":'ZZZJets',"color":ROOT.kAzure-6,"name":'Other ZZ processes'},{"sample":'WZZJets',"color":ROOT.kAzure-6,"name":'Other ZZ processes'},{"sample":'ttH126',"color":ROOT.kAzure-6,"name":'Other ZZ processes'}] 

SignalZZ_qq_Mad =  [{"sample":'ZZJetsTo4L',"name":'qq/qg #rightarrow ZZ(+jets)'}]
#SignalZZ_qq_Mad =  [{"sample":'ZZJetsTo4L_Normal',"name":'qq/qg #rightarrow ZZ(+jets)'}]

#Powheg set
SignalSamples_Mad =  SignalZZ_qq_Mad + SignalZZ_gg + SignalZZ_VBS + Signal_Others
SignalSamples_Pow =  SignalZZ_qq_Pow + SignalZZ_gg + SignalZZ_VBS + Signal_Others

DataSamples =  [{"sample":'data',"name":'Data'}]

inputdir_RMC="./results/ZZRecoAnalyzer_SR/"
inputdir_MC = "./results/ZZMCAnalyzer_MC/"
inputdir_CR="./results/ZZRecoAnalyzer_CR/" #DEL

# CenterData = "./FinalResults_"+MCSet+"/Data.root"
# UpData = "./FinalResults_"+MCSet+"/DataUp.root"
# DownData = "./FinalResults_"+MCSet+"/DataDown.root"

# CenterDataUF = "./FinalResults_"+MCSet+"/DataUnfold.root"
# UpDataUF = "./FinalResults_"+MCSet+"/DataUnfoldUp.root"
# DownDataUF = "./FinalResults_"+MCSet+"/DataUnfoldDown.root"

# MCSample = "./FinalResults_"+MCSet+"/MC.root"
# MCRecoSample = "./FinalResults_"+MCSet+"/MCReco.root"

BRele = 0.03363
BRmu = 0.03366
Lumi=19712 #pb-1

GlobSistList = [{"name":"Trig","value":0.015},{"name":"IsoId","value":0.015},{"name":"Lumi","value":0.026},{"name":"Acc","value":0.05},{"name":"Eff","value":0.015}]
#GlobSistList = [{"name":"Trig","value":0.015},{"name":"Lumi","value":0.026},{"name":"Acc","value":0.05}]

#DiffSistListUnfold = ("RedBkg","IrrBkg","qqgg","MCgen","UnfDataOverGenMC")
DiffSistListUnfold = ("RedBkg","IrrBkg","qqgg","MCgen")

DiffSistListJetsUnfold = ("JES_ModData","JER")

DiffSistList = ("Red","Irr") 
#DiffSistList = ("Red","Irr","sFactor") 
