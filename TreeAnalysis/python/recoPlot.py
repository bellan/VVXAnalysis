#! /usr/bin/env python
from optparse import OptionParser
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import TH1F,TCanvas, TLegend
#from ZZj4lReco_Mods_New import*
from recoPlotUtils import*
import sys,ast
import math
import operator


regions = ['SR','CR','CR2P2F','CR3P1F','SR_HZZ','CR3P1F_HZZ','CR2P2F_HZZ']

parser = OptionParser(usage="usage: %prog <final state> [options]")

parser.add_option("-r", "--region", dest="region",
                  default="SR",
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

parser.add_option("-m", "--mcset", dest="mcSet",
                  default="pow",
                  help= "Monte Carlo Set, pow for Powheg, mad for amcatnlo")

parser.add_option("-v", "--visual", dest="doVisulize",
                  default="True",
                  help= "Monte Carlo Set, pow for Powheg, mad for amcatnlo")

parser.add_option("-l", "--lumiProj", dest="LumiProj",
                  default="",
                  help="Lumi projection")


#REMEMBER ADD DEFINTION PLOT


(options, args) = parser.parse_args()

DoData=options.DoData
region =options.region
state = options.state
category = options.category
Type = options.Type
Save = options.SavePlot
mcSet = options.mcSet
doVisulize = ast.literal_eval(options.doVisulize)
LumiProj = options.LumiProj

print mcSet

c1 = TCanvas( 'c1', 'Reco Plots', 200, 10, 900, 700 )
leg = ROOT.TLegend(.78,.10,.92,.22);

FinState= options.FinalState

Addfake=True

InfoType = {"Mass":["M_{4l}","M_{4l}",10],"mJJ":["mjj","mJJ",4],"Z1Mass":["Z1 Mass","M_{2l}",10,],"Z2Mass":["Z2 Mass","M_{2l}",10,],"Z1lep0_sip":["Z1 lep 0 Sip","Sip",3],"Z1lep0_iso":["Z1 lep 0 Iso","Iso",3],"Z0lep0_pt":["Z1 lep 0 pT","p_{T}",4],"nJets":["# jets","# jets",1],"nJetsOrig":["# jets","# jets",1],"z":["z1","z1",1],"ptJet":["pT Jet","p_{T}^{jet}",4],"etaJet":["#eta Jet","#eta^{jet}",9],"Z1pt":["Z1 p_{T}","p_{T}",10],"Z2pt":["Z2 p_{T}","p_{T}",10],"Z1z":["Z1 z","z_{Z_{1}}",4],"Z2z":["Z2 z","z_{Z_{2}}",4],"ptJRatio":["","#Sigma p_{T}/# Sum  ",2],"ptRatio":["","#Sum p_{T}",2]} 


try:
    os.stat("./Plot/RecoPlots/")
except:
    os.mkdir("./Plot/RecoPlots/")
    
Var = Type
InputDir = "ZZjAnalyzer_"
    
if LumiProj!="":  InputDir+=LumiProj+"fbm1_"

#if region=="SR" and  category =="All" and state=="Is":  category =  "AllIs"

if category=="Sig":
    Addfake=False

if category=="RedBkg":
    hData = GetFakeRate("results/"+InputDir,"ZZTo"+FinState+"_"+Var,"data",InfoType[Type][2])            
    (hMC,leg)=GetMCPlot("results/"+InputDir+region.replace("CR","SR")+"/",category,"ZZTo"+FinState+"_"+Var,False,mcSet,InfoType[Type][2])          
elif category=="IrrBkg":      
    (hMC,leg)=GetMCPlot("results/"+InputDir+region+"/",category,"ZZTo"+FinState+"_"+Var,False,mcSet,InfoType[Type][2])    
    DoData=False
elif category=="Bkg":      
    (hMC,leg)=GetMCPlot("results/"+InputDir+region+"/",category,"ZZTo"+FinState+"_"+Var,True,mcSet,InfoType[Type][2])    

elif category=="CR":      
    (hMC,leg)=GetMCPlot("results/"+InputDir+region+"/",category,"ZZTo"+FinState+"_"+Var,False,mcSet,InfoType[Type][2])    
    hData=GetDataPlot("results/"+InputDir+region+"/","ZZTo"+FinState+"_"+Var,region,InfoType[Type][2])
else:
    if state=="Fs":
        (hMC,leg)=GetMCPlot_fstate("results/"+InputDir+region+"/",category ,"ZZTo"+FinState+"_"+Var,Addfake,mcSet,InfoType[Type][2]) 
    else:
        (hMC,leg)=GetMCPlot("results/"+InputDir+region+"/",category,"ZZTo"+FinState+"_"+Var,Addfake,mcSet,InfoType[Type][2])    
      
    hData=GetDataPlot("results/"+InputDir+region+"/","ZZTo"+FinState+"_"+Var,region,InfoType[Type][2])
    
YMaxMC = YMax=hMC.GetMaximum()

if category=="IrrBkg" or category == "RedBkg": YMaxData = 0. # YMaxData = YMax=hData.GetMaximum()
else:  YMaxData = ROOT.TMath.MaxElement(hData.GetN(),hData.GetEYhigh()) + ROOT.TMath.MaxElement(hData.GetN(),hData.GetY())


#print hData.GetN(),ROOT.TMath.MaxElement(hData.GetN(),hData.GetEYhigh()),ROOT.TMath.MaxElement(hData.GetN(),hData.GetY()),hData.GetEYhigh(),YMaxData,YMaxMC,YMaxData>YMaxMC,DoData

if YMaxData>YMaxMC and DoData:
    YMax = YMaxData
else: YMax = YMaxMC

print Type,len(InfoType[Type])
if len(InfoType[Type])>3: 
    hMC.SetMaximum(InfoType[Type][3])
else: hMC.SetMaximum(YMax+YMax*0.2)

print "MAX",YMax

hMC.Draw("hist")

hMC.GetHistogram().GetXaxis().SetTitle(InfoType[Type][0])
hMC.SetName(InfoType[Type][1])

if DoData:
    hData.SetMarkerStyle(20)
    hData.SetMarkerSize(.9)
    hData.Draw("samep")

if Type=="nJets":
    hMC.GetHistogram().GetXaxis().SetTitle("#Jets")
    hMC.GetHistogram().GetXaxis().SetBinLabel(1,"0 Jets")
    hMC.GetHistogram().GetXaxis().SetBinLabel(2,"1 Jets")
    hMC.GetHistogram().GetXaxis().SetBinLabel(3,"2 Jets")
    hMC.GetHistogram().GetXaxis().SetBinLabel(4,"3 Jets")
    hMC.GetHistogram().GetXaxis().SetBinLabel(4,">3 Jets")


leg.Draw("same")    

print category

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
elif category == "CR":
    Title = "CR"
elif category == "Sig13TeV":
    Title = "Sig13TeV"


Title=Var+"_"+Title+"_"+mcSet+"_"+region+"_"+FinState

ROOT.gStyle.SetOptStat(0);   
c1.Update()

c1.SetTitle(Title)

c1.SaveAs("Plot/RecoPlots/"+Title+".root")        
c1.SaveAs("Plot/RecoPlots/"+Title+".png")        
c1.SaveAs("Plot/RecoPlots/"+Title+".eps")        
c1.SaveAs("Plot/RecoPlots/"+Title+".pdf")  

if Save:
    if(LumiProj!=""): 
        c1.SaveAs("~/www/PlotsVV/13TeV_"+LumiProj+"fb/"+Title+".png")        
        c1.SaveAs("~/www/PlotsVV/13TeV_"+LumiProj+"fb/"+Title+".pdf")        
    else:
        c1.SaveAs("~/www/PlotsVV/13TeV/"+Title+".png")        
        c1.SaveAs("~/www/PlotsVV/13TeV/"+Title+".pdf")        


fout  = ROOT.TFile("results/data/"+Title+".root","update")
hwrite  = copy.deepcopy(hMC.GetStack().Last())
hwrite.SetName(Type+"_"+LumiProj)
hwrite.Write()

if doVisulize: CloseVar=input("digit anything you want to end the script \n")        
