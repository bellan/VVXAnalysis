import ROOT
from ROOT import gSystem


DataSamples = [{"sample":'data',"name":'Data'}]

BkgSamples  = [{"sample":'WWZ',"name":'Irreducible background'},{"sample":'TTZToLL',"name":'Irreducible background'}] # ,{"sample":'TTWWJets',"name":'Irreducible background'}  To be added

SignalZZ_qq_Pow = [{"sample":'ZZTo4l',"name":'qq/qg #rightarrow ZZ(+jets)',"kfac":1.1}]
SignalZZ_qq_Mad = [{"sample":'ZZTo4lamcatnlo',"name":'qq/qg #rightarrow ZZ(+jets)',"kfac":1.1}]

SignalZZ_gg=[{"sample":'ggZZ2e2mu',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)',"kfac":1.7},{"sample":'ggZZ4e',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)',"kfac":1.7},{"sample":'ggZZ4mu',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)',"kfac":1.7}] #FIX use  gg sample with higgs

SignalZZ_VBS = [{"sample":'ZZTo2e2muJJ',"name":'qq/qg #rightarrow ZZ(+jets)'},{"sample":'ZZTo4eJJ',"name":'qq/qg #rightarrow ZZ(+jets)'},{"sample":'ZZTo4muJJ',"name":'qq/qg #rightarrow ZZ(+jets)'}] 

#Signal_Others = [{"sample":'ZZZJets',"color":ROOT.kAzure-6,"name":'Other ZZ processes'},{"sample":'WZZJets',"color":ROOT.kAzure-6,"name":'Other ZZ processes'},{"sample":'ttH126',"color":ROOT.kAzure-6,"name":'Other ZZ processes'}] 


#Powheg Set
SignalSamples_Pow =  SignalZZ_qq_Pow + SignalZZ_gg + SignalZZ_VBS 

#Amcatnlo Set
SignalSamples_Mad =  SignalZZ_qq_Mad  + SignalZZ_gg + SignalZZ_VBS

DataSamples =  [{"sample":'data',"name":'Data'}]

inputdir_RMC ="./results/ZZRecoAnalyzer_SR/"
inputdir_MC  = "./results/ZZMCAnalyzer_MC/"
inputdir_CR  ="./results/ZZRecoAnalyzer_CR/"

BRele = 0.03363
BRmu  = 0.03366
Lumi  = 12892
#Lumi   = 15941

xs_wide = 16.5

# Theoretic xs for the tight fiducial region MCFM

xs_tight = {"2e2m":16.37,"4m":8.19,"4e":8.19,"4l":32.75} #To be corrected with MCFM values

GlobSystList           = [{"name":"Trig","value":0.02},{"name":"Lumi","value": 0.062}] 
DiffSystList           = [{"name":"Red","longname":"Reducible background"},{"name":"Irr","longname":"Irreducible background"},{"name":"sFactor","longname":"Scale Factor"},{"name":"MCgen","longname":"Monte Carlo choice"}]

DiffSystList           = [{"name":"Red","longname":"Reducible background"},{"name":"Irr","longname":"Irreducible background"},{"name":"sFactor","longname":"Scale Factor"},{"name":"MCgen","longname":"Monte Carlo choice"}]

DiffSystListUnfold           = [{"name":"RedBkg","longname":"Reducible background"},{"name":"IrrBkg","longname":"Irreducible background"},{"name":"SFSq","longname":"Scale Factor"},{"name":"MCgen","longname":"Monte Carlo choice"}]
#DiffSystListUnfold     = ("RedBkg","IrrBkg","qqgg","MCgen","SFSq")
DiffSystListJetsUnfold = [{"name":"JES_ModData","longname":"JES data correction"},{"name":"JES_ModMat","longname":"JES Mat correction"},{"name":"JER","longname":"Jet energy resolution"}]
#DiffSystListUnfold     = ("RedBkg","IrrBkg","qqgg","MCgen","UnfDataOverGenMC","SFSq")
