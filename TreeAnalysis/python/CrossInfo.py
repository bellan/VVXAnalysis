import ROOT
from ROOT import gSystem

DataSamples = [{"sample":'data',"name":'Data'}]
BkgSamples  = [{"sample":'WWZ',"name":'Irreducible background'},{"sample":'TTZToLL',"name":'Irreducible background'}] 
SignalZZ_qq_Pow = [{"sample":'ZZTo4l',"name":'qq/qg #rightarrow ZZ(+jets)',"kfac":1.1}]
SignalZZ_qq_Mad = [{"sample":'ZZTo4lamcatnlo',"name":'qq/qg #rightarrow ZZ(+jets)',"kfac":1.1}]

SignalZZ_gg=[{"sample":'ggTo2e2mu_Contin_MCFM701',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)',"kfac":1.7},{"sample":'ggTo4e_Contin_MCFM701',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)',"kfac":1.7},{"sample":'ggTo4mu_Contin_MCFM701',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)',"kfac":1.7}]

SignalZZ_VBS = [{"sample":'ZZTo2e2muJJ',"name":'qq/qg #rightarrow ZZ(+jets)'},{"sample":'ZZTo4eJJ',"name":'qq/qg #rightarrow ZZ(+jets)'},{"sample":'ZZTo4muJJ',"name":'qq/qg #rightarrow ZZ(+jets)'}] 

Signal_Others = [{"sample":'ZZZ',"color":ROOT.kAzure-6,"name":'Other ZZ processes'},{"sample":'WZZ',"color":ROOT.kAzure-6,"name":'Other ZZ processes'}] 

SignalZZ_gg=[{"sample":'gg_4l',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)',"kfac":1.7}]
SignalZZ_VBS = [{"sample":'qq_4l2j',"name":'qq/qg #rightarrow ZZ(+jets)'}]

#Powheg Set
SignalSamples_Pow =  SignalZZ_qq_Pow + SignalZZ_gg + SignalZZ_VBS #+ Signal_Others 

#Amcatnlo Set
SignalSamples_Mad =  SignalZZ_qq_Mad  + SignalZZ_gg + SignalZZ_VBS #+ Signal_Others 

SignalSamples_Compact = [{"sample":'Signal',"color":ROOT.kAzure-5,"name":'all'}]
BkgSamples_Compact = [{"sample":'Irr',"color":ROOT.kAzure-5,"name":'Irr'}]

DataSamples =  [{"sample":'data',"name":'Data'}]
VarList = ["Mass","nJets","nJets_Central","Mjj","Mjj_Central","Deta","Deta_Central","PtJet1","PtJet2","EtaJet1","EtaJet2","dRZZ" ,"PtZZ"]

inputdir_RMC ="./results/ZZRecoAnalyzer_SR/"
inputdir_MC  = "./results/ZZMCAnalyzer_MC/"
inputdir_CR  ="./results/ZZRecoAnalyzer_CR/"

BRele  = 0.03363
BRmu   = 0.03366
Lumi   = 35900

xs_wide = 16.5

# Theoretic xs for the tight fiducial region MCFM

xs_tight_exp  = {"2e2m":21.0,"4m":10.8,"4e":9.9,"4l":42.2}  # ZZ inclusive official.  Powheg 
xs_tight      = {"2e2m":20.1630,"4m":10.0270,"4e":10.0270,"4l":40.2170}  #Born level. ZZ Inclusive official. Pow*1.1 + MCFM*1.7 


Xs_OS_2e2m =  32.64
xs_OS_4m   =  16.32
xs_OS_4e   =  16.32

GlobSystList           = [{"name":"Trig","value":0.02},{"name":"Lumi","value": 0.026}] 

PdfSyst_fid=[{"name":"Pdf","value":0.00}] 
PdfSyst    =[{"name":"Pdf","value":0.00}] 

DiffSystList           = [{"name":"Red","longname":"Reducible background"},{"name":"Irr","longname":"Irreducible background"},{"name":"sFactor","longname":"Scale Factor"},{"name":"MCgen","longname":"Monte Carlo choice"},{"name":"Pu","longname":"Pile up"}] #add pdf

DiffSystListJets = [{"name":"JES_ModMat","longname":"JES Mat correction"},{"name":"JER","longname":"Jet energy resolution"}]

DiffSystListUnfold     = [{"name":"RedBkg","longname":"Reducible background"},{"name":"IrrBkg","longname":"Irreducible background"},{"name":"EleSFSq","longname":"Electron Scale Factor"},{"name":"MuSFSq","longname":"Muon Scale Factor"},{"name":"MCgen","longname":"Monte Carlo choice"},{"name":"qqgg","longname":"qq gg cross section"},{"name":"Pu","longname":"Pileup"},{"name":"PDF","longname":"PDF"},{"name":"As","longname":"\Alpha_S"}]

#DiffSystListUnfold     = [{"name":"RedBkg","longname":"Reducible background"},{"name":"IrrBkg","longname":"Irreducible background"},{"name":"SFSq","longname":"Scale Factor"},{"name":"MCgen","longname":"Monte Carlo choice"},{"name":"qqgg","longname":"qq gg cross section"}]

DiffSystListJetsUnfold = [{"name":"JES_ModMat","longname":"JES Mat correction"},{"name":"JER","longname":"Jet energy resolution"}]


