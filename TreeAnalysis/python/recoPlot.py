#! /usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import CrossInfo
from CrossInfo import* 
from ROOT import TH1F,TCanvas, TLegend
from recoPlotUtils import*
import sys,ast
import math
import operator
import CMS_lumi, tdrstyle
from PersonalInfo import*

regions = ['SR4P','CR4L','CR2P2F','CR3P1F','SR_HZZ','CR3P1F_HZZ','CR2P2F_HZZ']

parser = OptionParser(usage="usage: %prog <final state> [options]")

parser.add_option("-r", "--region", dest="region",
                  default="SR4P",
                  help="Region type are {0:s}. Default is SR.".format(', '.join(regions)))

parser.add_option("-f", "--finalstate", dest="FinalState",
                  default="4l",
                  help="Final state are 4l, 4m, 2e2m and 4e. Default is 4l")

parser.add_option("-s", "--state", dest="state",
                  default="Is",
                  help="State type are Fs for finalstate and Is for initial state. Default is Is.")

parser.add_option("-c", "--category", dest="category",
                  default="All",
                  help="Category type are All, Sig, IrrBkg, RedBkg, or combination like IsSig+IrrBkg")

parser.add_option("-d", "--data", dest="DoData",
                  action="store_true",
                  default=True,
                  help="Data is True or False to plot or not the data. Default is True")

parser.add_option("-t", "--type", dest="Type",
                  default="Mass",
                  help= "type type to choose the  plot you want. Mass, Jets, DeltaEta, mjj")

parser.add_option("-F", "--Fiducial", dest="DoFiducial",
                  default="False",
                  help= "type fiducial region to plot the fiducial region plot")

parser.add_option("-S", "--Save", dest="SavePlot",
                  action="store_true",
                  default=False,
                  help="Save plot option, default is False")

parser.add_option("-p", "--plot", dest="showPlot",
                  action="store_true",
                  default=False,
                  help="Show plot option, default is False")

parser.add_option("-m", "--mcset", dest="mcSet",
                  default="mad",
                  help= "Monte Carlo Set, pow for Powheg, mad for amcatnlo")

parser.add_option("-l", "--lumiProj", dest="LumiProj",
                  default="",
                  help="Lumi projection")

parser.add_option("-D", "--Dir", dest="Dir",
                  default="test",
                  help="Directory where save plots")

parser.add_option("-A", "--Analysis", dest="Analysis",
                  default="ZZ",
                  help="Analysis. Default is ZZ. Othe oprtion is VBS")

parser.add_option("-y", "--year", dest="year",
                  default="2016",
                  help= "valid inputs are 2016, 2017, 2018")


#REMEMBER ADD DEFINTION PLOT


(options, args) = parser.parse_args()

DoData     = options.DoData
region     = options.region
state      = options.state
category   = options.category
Type       = options.Type
Save       = options.SavePlot
mcSet      = options.mcSet
LumiProj   = options.LumiProj
Dir        = options.Dir
Analysis   = options.Analysis
showPlot   = options.showPlot
year       = options.year

tdrstyle.setTDRStyle()

if not showPlot: ROOT.gROOT.SetBatch(True)

if category=="IrrBkg" or category=="RedBkg":      c1 = TCanvas( 'c1', Type+"_"+mcSet, 200, 10, 800, 800 )
else:                       c1 = TCanvas( 'c1', Type+"_"+mcSet , 200, 10, 900, 1200 )

#leg = ROOT.TLegend(.78,.10,.8,.22);

FinState= options.FinalState

Addfake=True

InfoType_zz = {"Mass":["m_{4l} [GeV]","m_{4\ell}",10],"Mjj":["m_{jj} [GeV]","m_{JJ}",20],"Z1Mass":["Z1 Mass","m_{2\ell}",10,],"Z2Mass":["Z2 Mass","m_{2\ell}",10,],"Z1lep0_sip":["Z1 lep 0 Sip","Sip",4],"Z1lep0_iso":["Z1 lep 0 Iso","Iso",4],"Z0lep0_pt":["Z1 lep 0 pT","p_{T}",4],"nJets":["N_{jets} (|#eta^{jet}| < 4.7)","N_{jets} (|#eta^{jet}| < 4.7)",1],"nJets_central":["N_{jets} (|#eta^{jet}| < 4.7)","N_{jets} (|#eta^{jet}| < 4.7)",1],"z":["z1","z1",1],"PtJet1":["p_{T}^{jet1} [GeV]","p_{T}^{jet}",1],"EtaJet1":["#eta^{jet1}","#eta^{jet}",9],"PtJet2":["p_{T}^{jet2} [GeV]","p_{T}^{jet}",1],"EtaJet2":["#eta^{jet2}","#eta^{jet}",10],"Z1pt":["Z1 p_{T}","p_{T}",20],"Z2pt":["Z2 p_{T}","p_{T}",10],"Z1z":["Z1 z","z_{Z_{1}}",7],"Z2z":["Z2 z","z_{Z_{2}}",7],"ptJRatio":["","#Sigma p_{T}/# Sum  ",2],"ptRatio":["","#Sum p_{T}",2],"PtZZ":["p_{T}^{4\\ell}","Sum p_{T}",20],"deltaEtaJJ":["|#eta_{jj}|","|#eta_{jj}|",2],"Dphi":["#Delta #phi_{jj}","#Delta #phi_{jj}",10],"Deta":["|#Delta#eta_{jj}|","#Delta #eta_{jj}",5],"Mjj_Central":["m_{jj}","m_{jj}",20],"Deta_Central":["#Delta #eta_{jj}","#Delta #eta_{jj}",5],"Deta2Jet":["#Delta #eta_{jj}, 2 jet","#Delta #eta_{jj} =2 jet",5],"Deta_noCentral":["#Delta #eta_{jj}, >2 jet","#Delta #eta_{jj} > 2 jet",5],"Deta_1noCentral":["#Delta #eta_{jj}, >2 jet","#Delta #eta_{jj} > 2 jet",5],"PtJet1_noCentral":["#eta Jet","#eta^{jet}",9],"EtaJet1_noCentral":["#eta Jet","#eta^{jet}",10]}

InfoType_vbs = {"Mass":["m_{4\ell}","m_{4\ell}",40],"Mjj":["m_{jj}","m_{JJ}",20],"Z1Mass":["Z1 Mass","m_{2\ell}",10,],"Z2Mass":["Z2 Mass","m_{2\ell}",10,],"Z1lep0_sip":["Z1 lep 0 Sip","Sip",4],"Z1lep0_iso":["Z1 lep 0 Iso","Iso",4],"Z0lep0_pt":["Z1 lep 0 pT","p_{T}",4],"nJets":["# jets","# jets",1],"nJets_central":["# jets","# jets",1],"z":["z1","z1",1],"PtJet1":["pT Jet","p_{T}^{jet}",10],"EtaJet1":["#eta Jet","#eta^{jet}",10],"PtJet2":["pT Jet","p_{T}^{jet}",10],"EtaJet2":["#eta Jet","#eta^{jet}",10],"Z1pt":["Z1 p_{T}","p_{T}",20],"Z2pt":["Z2 p_{T}","p_{T}",10],"Z1z":["Z1 z","z_{Z_{1}}",7],"Z2z":["Z2 z","z_{Z_{2}}",7],"ptJRatio":["","#Sigma p_{T}/# Sum  ",2],"ptRatio":["","#Sum p_{T}",2],"PtZZ":["p_{T}^{4\\ell}","Sum p_{T}",60],"deltaEtaJJ":["|#eta_{jj}|","|#eta_{jj}|",2],"Dphi":["#Delta #phi_{jj}","#Delta #phi_{jj}",10],"Deta":["#Delta #eta_{jj}","#Delta #eta_{jj}",5],"Mjj_Central":["m_{jj}","m_{jj}",20],"Deta_Central":["#Delta #eta_{jj}","#Delta #eta_{jj}",5]}


if Analysis =="ZZ": InfoType=InfoType_zz
else: InfoType=InfoType_vbs

##font[3]{l}
#set the tdr style

#change the CMS_lumi variables  (see CMS_lumi.py)                                                                                  
#CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"                               
#CMS_lumi.lumi_8TeV = "19.7 fb^{-1}"                                                                                                                                                                                                          

lumi = round(Lumi/1000.,1)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "{0}".format(lumi)+" fb^{-1} (13 TeV)\n"

iPos = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 0

try:
    os.stat("./Plot/RecoPlots/")
except:
    os.mkdir("./Plot/RecoPlots/")
    
Var = Type
if Analysis =="ZZ":
    InputDir = year+"/ZZjAnalyzer_"
elif Analysis =="VBS":
    InputDir = "VBSAnalyzer_"
else: sys.exit("Wrong analysis. Chose between ZZ and VBS")

if LumiProj!="":  InputDir+=LumiProj+"fbm1_"

#if region=="SR" and  category =="All" and state=="Is":  category =  "AllIs"


if category=="Sig":
    Addfake=False

if category=="RedBkg":
    hData = GetFakeRate("results/"+InputDir,"ZZTo"+FinState+"_"+Var,"data",InfoType[Type][2])            
    (hMC,leg)=GetMCPlot("results/"+InputDir+region.replace("CR4L","SR4P")+"/",category,"ZZTo"+FinState+"_"+Var,False,mcSet,InfoType[Type][2])          
elif category=="IrrBkg":      
    (hMC,leg)=GetMCPlot("results/"+InputDir+region+"/",category,"ZZTo"+FinState+"_"+Var,False,mcSet,InfoType[Type][2])    
    DoData=False
elif category=="Bkg":      
    (hMC,leg)=GetMCPlot("results/"+InputDir+region+"/",category,"ZZTo"+FinState+"_"+Var,True,mcSet,InfoType[Type][2])    

elif category=="CR4L":      
    (hMC,leg)=GetMCPlot("results/"+InputDir+region+"/",category,"ZZTo"+FinState+"_"+Var,False,mcSet,InfoType[Type][2])    
    (hData,histodata)=GetDataPlot("results/"+InputDir+region+"/","ZZTo"+FinState+"_"+Var,region,InfoType[Type][2])
else:
    if state=="Fs":
        (hMC,leg)=GetMCPlot_fstate("results/"+InputDir+region+"/",category ,"ZZTo"+FinState+"_"+Var,Addfake,mcSet,InfoType[Type][2]) 
    else:
        (hMC,leg)=GetMCPlot("results/"+InputDir+region+"/",category,"ZZTo"+FinState+"_"+Var,Addfake,mcSet,InfoType[Type][2])    
    (hData,histodata)=GetDataPlot("results/"+InputDir+region+"/","ZZTo"+FinState+"_"+Var,region,InfoType[Type][2])
    print("int",histodata.Integral())
    
YMaxMC = YMax=hMC.GetMaximum()

if category=="IrrBkg" or category == "RedBkg": YMaxData = 0. # YMaxData = YMax=hData.GetMaximum()
else:
    YMaxData = ROOT.TMath.MaxElement(hData.GetN(),hData.GetEYhigh()) + ROOT.TMath.MaxElement(hData.GetN(),hData.GetY())
        
#print hData.GetN(),ROOT.TMath.MaxElement(hData.GetN(),hData.GetEYhigh()),ROOT.TMath.MaxElement(hData.GetN(),hData.GetY()),hData.GetEYhigh(),YMaxData,YMaxMC,YMaxData>YMaxMC,DoData

if category != "Sig":
    hMCErr =  copy.deepcopy(hMC.GetStack().Last())
    YMaxMC = hMCErr.GetBinContent(hMCErr.GetMaximumBin())+hMCErr.GetBinError(hMCErr.GetMaximumBin())

if YMaxData>YMaxMC and DoData:
    YMax = YMaxData
else: YMax = YMaxMC

YMax*=1.37

# if len(InfoType[Type])>3: 
#     hMC.SetMaximum(InfoType[Type][3])
# else: hMC.SetMaximum(YMax+YMax*0.2)



if category !="IrrBkg" and category !="RedBkg":
    HistoData = copy.deepcopy(histodata)
    c1.cd()
    pad1 = ROOT.TPad ('hist', '', 0., 0.22, 1.0, 1.0)#0.35
    pad1.SetTopMargin    (0.10)
    pad1.SetRightMargin  (0.06)#0.10
    pad1.SetLeftMargin   (0.1)
    pad1.SetBottomMargin (1.5) 
    pad1.Draw()
    
    c1.cd()

    pad2 = ROOT.TPad ('rat', 'Data/MC ratio', 0., 0.0,  1., 0.28)#0.15
    pad2.SetTopMargin (0.01)
    pad2.SetRightMargin (0.06)#0.10
    pad2.SetLeftMargin (0.1)
    pad2.SetBottomMargin(0.3);
    pad2.Draw()
    
    pad1.cd()

    histodata.Divide(hMC.GetStack().Last())    
else:
    c1.SetLeftMargin(0.15)
    c1.SetRightMargin(0.08)
hMC.SetMaximum(YMax)
hMC.Draw("hist")

if category !="IrrBkg" and category !="RedBkg": 
    hMC.GetHistogram().GetXaxis().SetLabelSize(0)
else: 
    hMC.SetName(InfoType[Type][1])
    hMC.GetHistogram().GetXaxis().SetTitle(InfoType[Type][0])
    hMC.GetHistogram().GetXaxis().SetTitleSize(0.05)
    hMC.GetHistogram().GetYaxis().SetTitleSize(0.05)
    hMC.GetHistogram().GetXaxis().SetTitleOffset(1)
    hMC.GetHistogram().GetYaxis().SetTitle("Events")

if category != "Sig":
    hMCErr.SetFillStyle(3005)
    hMCErr.SetMarkerStyle(1)
    hMCErr.SetFillColor(ROOT.kBlack)
    hMCErr.Draw("sameE2")
    leg.AddEntry(hMCErr, "MC Err","f")
if DoData:
    hData.SetMarkerStyle(20)
    hData.SetMarkerSize(.9)
    hData.Draw("samep")
    leg.AddEntry(hData, "Data","lpe")
    #HistoData.Draw("same text")
    if category !="IrrBkg" and category !="RedBkg":    HistoData.Draw("same")

if Type=="nJets":
    hMC.GetHistogram().GetXaxis().SetTitle("N_{jets} (|#eta^{jet}| < 4.7)")
    hMC.GetHistogram().GetXaxis().SetBinLabel(1,"0 ")
    hMC.GetHistogram().GetXaxis().SetBinLabel(2,"1 ")
    hMC.GetHistogram().GetXaxis().SetBinLabel(3,"2 ")
    hMC.GetHistogram().GetXaxis().SetBinLabel(4,"3 ")
    hMC.GetHistogram().GetXaxis().SetBinLabel(5,">3 ")
    if category != "IrrBkg" and category !="RedBkg":
        histodata.GetXaxis().SetTitle("N_{jets} (|#eta^{jet}| < 4.7)")
        histodata.GetXaxis().SetBinLabel(1,"0 ")
        histodata.GetXaxis().SetBinLabel(2,"1 ")
        histodata.GetXaxis().SetBinLabel(3,"2 ")
        histodata.GetXaxis().SetBinLabel(4,"3 ")
        histodata.GetXaxis().SetBinLabel(5,">3 ")


leg.Draw("same")    

CMS_lumi.CMS_lumi(c1, iPeriod, iPos)

#print Ratio.GetNbinsX(),hMC.GetStack().Last().GetNbinsX()

if category !="IrrBkg" and category !="RedBkg":

    pad2.cd()
    
    Line = ROOT.TLine(hMC.GetXaxis().GetXmin(),1,hMC.GetXaxis().GetXmax(),1) 
    Line.SetLineWidth(2)
    histodata.GetXaxis().SetTitle(InfoType[Type][0])
    histodata.GetXaxis().SetLabelSize(0.08)
    histodata.GetXaxis().SetTitleSize(0.08)
    #histodata.Draw("AXIS")
    histodata.SetMarkerStyle(20)
    histodata.Draw("E1")
    Line.Draw("same")

if category == 'All':
    Title= "All"
elif category=="AllIs":
    Title="Initial_state"
elif category == 'Sig':
    Title = "Signal"
elif category == 'IrrBkg':
    Title = "Irreducible_Bkg"
elif category == "data":
    Title = "data"
elif category == "RedBkg":
    Title = "Reducible_Bkg"
elif category == "All13TeV":
    Title = "All13TeV"
if state == "Fs":
    Title = "Final_State"
elif category == "CR4L":
    Title = "CR4L"
elif category == "Sig13TeV":
    Title = "Sig13TeV"


Title=Var+"_"+Title+"_"+mcSet+"_"+region+"_"+FinState

ROOT.gStyle.SetOptStat(0);   
ROOT.gStyle.SetOptTitle(0)
c1.Update()

c1.SetTitle(Title)
c1.SaveAs("Plot/RecoPlots/"+Title+".root")        
c1.SaveAs("Plot/RecoPlots/"+Title+".png")        
c1.SaveAs("Plot/RecoPlots/"+Title+".eps")        
c1.SaveAs("Plot/RecoPlots/"+Title+".pdf")  


if Save:

    Dir+="/"

    try:
        os.stat(PersonalFolder+Dir)
    except:
        os.mkdir(PersonalFolder+Dir)
        os.system("cp "+PersonalFolder+"index.php "+PersonalFolder+Dir )
    try:
        os.stat(PersonalFolder+Dir+"/Reco/")
    except:
        os.mkdir(PersonalFolder+Dir+"/Reco/")
        os.system("cp "+PersonalFolder+"index.php "+PersonalFolder+Dir+"/Reco/" )

    if(LumiProj!=""): 
        c1.SaveAs("~/www/PlotsVV/13TeV_"+LumiProj+"fb/"+Title+".png")        
        c1.SaveAs("~/www/PlotsVV/13TeV_"+LumiProj+"fb/"+Title+".pdf")        
    else:
        c1.SaveAs(PersonalFolder+Dir+"Reco/"+Title+".png")        
        c1.SaveAs(PersonalFolder+Dir+"Reco/"+Title+".root")        
        c1.SaveAs(PersonalFolder+Dir+"Reco/"+Title+".pdf")        
        c1.SaveAs("~/../../../../../eos/user/g/gpinnaan/"+Dir+"Reco/"+Title+".pdf")        


