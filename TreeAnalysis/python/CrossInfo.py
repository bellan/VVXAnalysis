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

Signal_Others = [{"sample":'ZZZJets',"color":ROOT.kAzure-6,"name":'Other ZZ processes'},{"sample":'WZZJets',"color":ROOT.kAzure-6,"name":'Other ZZ processes'},{"sample":'ttH126',"color":ROOT.kAzure-6,"name":'Other ZZ processes'}] 

SignalZZ_qq_Mad =  [{"sample":'ZZJetsTo4L',"name":'qq/qg #rightarrow ZZ(+jets)'}]

#MadGraph and Powheg sets of samples
SignalSamples_Mad =  SignalZZ_qq_Mad + SignalZZ_gg + SignalZZ_VBS + Signal_Others
SignalSamples_Pow =  SignalZZ_qq_Pow + SignalZZ_gg + SignalZZ_VBS + Signal_Others

DataSamples =  [{"sample":'data',"name":'Data'}]

inputdir_RMC="./results/ZZRecoAnalyzer_SR/"
inputdir_MC = "./results/ZZMCAnalyzer_MC/"
inputdir_CR="./results/ZZRecoAnalyzer_CR/" 

BRele = 0.03363
BRmu = 0.03366
Lumi=19712 #pb-1

#List of global systematic uncertainties in the normal (wide) fiducial region:
#GlobSistList = [{"name":"Trig","value":0.015},{"name":"Lumi","value":0.026},{"name":"Acc","value":0.05}] 

#List of global systematic uncertainties in the tight fiducial region (the only difference is in the acceptance value):
GlobSistList = [{"name":"Trig","value":0.015},{"name":"Lumi","value":0.026},{"name":"Acc","value":0.01}] 

#List used if Unfold == False 
DiffSistList = ("Red","Irr","sFactor")

#List of systematic uncertainties propagated through the unfolding procedure
DiffSistListUnfold = ("RedBkg","IrrBkg","qqgg","MCgen","UnfDataOverGenMC","SFSq") 
DiffSistListJetsUnfold = ("JES_ModData","JES_ModMat","JER")#not to include for the estimate of the inclusive cross-section as the integral of the jet multiplicity.

#The following numbers are used to estimate the inclusive cross-section using the likelihood method
xs_2e2m = 10.15 #taken from 1507.06257v1 at NNLO (for tight region)
xs_4m = 4.90 #taken from 1507.06257v1 at NNLO (for tight region)
xs_4e =5.16 #taken from 1507.06257v1 at NNLO (for tight region)

xs_OS_2e2m =  1.#FIXME
xs_OS_4m   =  1.#FIXME
xs_OS_4e   =  1.#FIXME
