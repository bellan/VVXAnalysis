#! /usr/bin/env python
from optparse import OptionParser
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import TH1F,TCanvas, TLegend
from ZZj4lReco_Mods import*
import sys,ast
import math
import operator


regions = ['SR','CR','CR2P2F','CR3P1F','CR3P1F_HZZ','CR2P2F_HZZ']

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
                  default="True",
                  help="Data is True or False to plot or not the data. Default is True")

parser.add_option("-t", "--type", dest="Type",
                  default="Mass",
                  help= "type type to choose the  plot you want. Mass, Jets, DeltaEta, mjj")

parser.add_option("-F", "--Fiducial", dest="DoFiducial",
                  default="False",
                  help= "type fiducial region to plot the fiducial region plot")


#REMEMBER ADD DEFINTION PLOT


(options, args) = parser.parse_args()

DoData=ast.literal_eval(options.DoData)
region =options.region
state = options.state
category = options.category
Type = options.Type

c2 = TCanvas( 'c2', 'c2', 200, 10, 900, 700 )
leg = ROOT.TLegend(.65,.12,.85,.26);

FinState= options.FinalState

Addfake=True


try:
    os.stat("./Plot/RecoPlots/")
except:
    os.mkdir("./Plot/RecoPlots/")


if Type=="Mass": 
    Var = "Mass"
    InputDir = "ZZjAnalyzer_"
    
elif Type=="Jets":
    #    Var = "Jets_01"
    Var = "Jets_JERCentralSmear_01"
    InputDir = "ZZRecoAnalyzer_"

if region=="SR" and  category =="All" and state=="Is":  category =  "AllIs"

if category=="Sig" or category=="All13TeV":
    Addfake=False

if category=="RedBkg":
    hData = GetFakeRate("results/"+InputDir,"ZZTo"+FinState+"_"+Var,"data")            
    (hMC,leg)=GetMCPlot("results/"+InputDir+region+"/",category,"ZZTo"+FinState+"_"+Var,False)          
elif category=="IrrBkg":      
    (hMC,leg)=GetMCPlot("results/"+InputDir+region+"/",category,"ZZTo"+FinState+"_"+Var,False)    
    DoData=False
elif category=="Bkg":      
    (hMC,leg)=GetMCPlot("results/"+InputDir+region+"/",category,"ZZTo"+FinState+"_"+Var,True)    

else:
    if state=="Fs":
        (hMC,leg)=GetMCPlot_fstate("results/"+InputDir+region+"/",category ,"ZZTo"+FinState+"_"+Var,Addfake) 
    else:
        (hMC,leg)=GetMCPlot("results/"+InputDir+region+"/",category,"ZZTo"+FinState+"_"+Var,Addfake)    
      
    hData=GetDataPlot("results/"+InputDir+region+"/","ZZTo"+FinState+"_"+Var,region)

if category=="RedBkg" or category=="IrrBkg" or category=="All13TeV": YMax=hMC.GetMaximum()
else: YMax = ROOT.TMath.MaxElement(hData.GetN(),hData.GetY()); 


if YMax>hMC.GetMaximum() and DoData: 
    YErrMax = ROOT.TMath.MaxElement(hData.GetN(),hData.GetEYhigh()); 
    hMC.SetMaximum(YMax+YErrMax*1.1)

hMC.Draw("hist")
print hMC.GetNhists()

if Type=="Mass":
    hMC.GetHistogram().GetXaxis().SetTitle("M_{4l}")
    hMC.SetName("M_{4l}")

elif Type=="Jets":
    hMC.GetHistogram().GetXaxis().SetTitle("#Jets")
    hMC.GetHistogram().GetXaxis().SetBinLabel(1,"0 Jets")
    hMC.GetHistogram().GetXaxis().SetBinLabel(2,"1 Jets")
    hMC.GetHistogram().GetXaxis().SetBinLabel(3,"2 Jets")
    hMC.GetHistogram().GetXaxis().SetBinLabel(4,">2 Jets")
    
    #hMC.GetYaxis().SetTitle("Events/10GeV")

    hMC.SetTitle("")
    
    if DoData is True:
        leg.AddEntry(hData,"Data","lpe")
        if category=="RedBkg":  hData.Draw("same")
        else: hData.Draw("p")


elif Type=="Def":
    (hDef,leg2)=GetSignalDefPlot("results/"+InputDir+region+"/",category)
    hDef.Draw("hist")
    leg2.Draw("same")
    hDef.GetHistogram().GetXaxis().SetTitle("M_{4l}")
    hDef.GetYaxis().SetTitle("Events/10GeV")
    hDef.SetTitle("")
   
else:
    print "Choose Type between Norm and Def"

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

Title=Var+"_"+Title+"_Pow"

ROOT.gStyle.SetOptStat(0);   
c2.Update()

c2.SetTitle(Title)

c2.SaveAs("Plot/RecoPlots/"+Title+".root")        
c2.SaveAs("Plot/RecoPlots/"+Title+".png")        
c2.SaveAs("Plot/RecoPlots/"+Title+".eps")        
c2.SaveAs("Plot/RecoPlots/"+Title+".pdf")  

CloseVar=input("digit anything you want to end the script \n")        
