#! /usr/bin/env python

##################################
## G. Pinna (UNITO) - Jun 2015 ###
##################################

from optparse import OptionParser
import ROOT,copy
ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import gSystem, TCanvas, TH1,  TPad, gStyle, TLegend, TGraphAsymmErrors
from collections import OrderedDict
from CrossSection import*
import sys,ast,os
import math
import operator

parser = OptionParser(usage="usage: %prog [options]")


parser.add_option("-s", "--save", dest="SavePlot",
                  default="False",
                  help="Save plot option, default is False")

parser.add_option("-u", "--unfold", dest="UseUnfold",
                  default="True",
                  help="Use unfolded data or not, default is True")

parser.add_option("-r", "--recoMC", dest="UseMCReco",
                  default="False",
                  help="Use reconstructed MC samples or not, default is False")

parser.add_option("-t", "--Type", dest="Type",
                  default="Mass",
                  help="Type of differential Cross secion. ['Mass','Jets']. Default is Mass")

parser.add_option("-i", "--Inclusive", dest="DoInclusive",
                  default="False",
                  help="Compute inclusive Cross secion. Default is False ")

parser.add_option("-m", "--MCSet", dest="MCSet",
                  default="Pow",
                  help="Choose MC Set between Pow (Powheg) and Mad (MadGraph). Default is Pow")


(options, args) = parser.parse_args()


SavePlot  =ast.literal_eval(options.SavePlot)
UseUnfold =ast.literal_eval(options.UseUnfold)
UseMCReco =ast.literal_eval(options.UseMCReco)
DoInclusive =ast.literal_eval(options.DoInclusive)
MCSet = options.MCSet

Type = options.Type


try:
    os.stat("./Plot/CrossSection/")
except:
    os.mkdir("./Plot/CrossSection/")


if DoInclusive:
    print Red("\n\n##############################################################\n################# ZZ4l Inclusive Cross Section ###############\n##############################################################\n")
else:  
    print Red("\n\n#################################################################################\n############## ZZ4l Differential Cross Section as function of {0} ##############\n#################################################################################\n".format(Type))


SetGlobalVariable(MCSet)

hMCList = getCrossPlot_MC(Type)

if DoInclusive:
    TotalCross()
    sys.exit(0)

if UseMCReco:
    (hDataList,hDataListUp,hDataListDown) = getCrossPlot_Data(UseUnfold,Type,0,True) 
   
else:
    (hDataList,hDataListUp,hDataListDown)  = getCrossPlot_Data(UseUnfold,Type,0,False) 

for i in range(0,4):
    #print "{0} Tot Cross {1} {2:.2f} \n".format(hDataList[i]["name"],(15-len(hDataList[i]["name"]))*" ", hDataList[i]["state"].Integral(1,-1)) # Check total cross section without normalization
    hDataList[i]["state"].SetMarkerSize(1.2)
    hDataList[i]["state"].SetMarkerColor(1)
    hDataList[i]["state"].SetMarkerStyle(20)
    hDataList[i]["state"].SetLineColor(1)
    
    grSist = getSistGraph(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"])
    grSist.SetFillStyle(3001)
    grSist.SetFillColor(ROOT.kRed)
 
    Err=ROOT.Double(0.)
    
    LastBin = hMCList[i]["state"].GetXaxis().GetLast()    

    hRatio= ROOT.TH1F()
    hRatio = hDataList[i]["state"].Clone("hRatioNew")
    hRatio.SetStats(0)
    hRatio.Divide(hMCList[i]["state"])
    hRatio.SetMinimum(0);  
    hRatio.SetMaximum(3); 
    if Type=="Mass":
        xName = "M_{4l} [GeV]"
        hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d M_{ZZ} [pb/GeV]")

    elif Type=="Jets":
        xName = "#Jets"
        hRatio.GetXaxis().SetTitle("#Jets")
        hRatio.GetXaxis().SetBinLabel(1,"0 Jets")
        hRatio.GetXaxis().SetBinLabel(2,"1 Jets")
        hRatio.GetXaxis().SetBinLabel(3,"2 Jets")
        hRatio.GetXaxis().SetBinLabel(4,">2 Jets")

        hMCList[i]["state"].GetXaxis().SetTitle("#Jets")
        hMCList[i]["state"].GetXaxis().SetBinLabel(1,"0 Jets")
        hMCList[i]["state"].GetXaxis().SetBinLabel(2,"1 Jets")
        hMCList[i]["state"].GetXaxis().SetBinLabel(3,"2 Jets")
        hMCList[i]["state"].GetXaxis().SetBinLabel(4,">2 Jets")

        hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d #jets [pb]")        

    hRatio.GetXaxis().SetTitle(xName)
    
    hRatio.GetYaxis().SetTitle("Data/MC")
    
    Max=hDataList[i]["state"].GetMaximum()

    if Type == "Mass": Max=0.08 #Max=Max+Max/5.
    else: Max=Max+Max/9.
    
    hMCList[i]["state"].SetMaximum(Max) 
    hMCList[i]["state"].SetFillColor(ROOT.kCyan-3)
    hMCList[i]["state"].SetLineColor(1)
    hMCList[i]["state"].GetXaxis().SetTitle(xName)
    hMCList[i]["state"].SetTitle(hMCList[i]["name"]+" Cross Section")
    
    c1 = TCanvas( 'c1', hMCList[i]["name"] , 200, 10, 900, 700 )
    
    leg = ROOT.TLegend(.65,.52,.85,.66);


    leg.AddEntry(grSist,"Syst","f")
    leg.AddEntry(hDataList[i]["state"],"Data + Stat","lep")
    leg.AddEntry(hMCList[i]["state"],"MC","lpf")


    fInt = ROOT.TF1("constant","1",hMCList[i]["state"].GetXaxis().GetXmin(),    hMCList[i]["state"].GetXaxis().GetXmax());
        
    c1.cd()
    pad1 = ROOT.TPad ('hist', '', 0., 0.30, 1.0, 1.0)
    pad1.SetTopMargin (0.10)
    pad1.SetRightMargin (0.10)
    pad1.SetLeftMargin (0.10)
    pad1.Draw()
        
    c1.cd()
    pad2 = ROOT.TPad ('rat', 'Data/MC ratio', 0., 0.0,  1., 0.3)
    pad2.SetTopMargin (0.10)
    pad2.SetRightMargin (0.10)
    pad2.SetLeftMargin (0.10)
    pad2.Draw()
    
    pad1.cd()
    hMCList[i]["state"].Draw("hist")
    
    if not UseMCReco:
        grSist.Draw("same2")
    hDataList[i]["state"].Draw("sameE1")
    leg.Draw("same")
    
    pad2.cd()
    hRatio.Draw("E1")
    fInt.Draw("same")
    
    ROOT.gStyle.SetOptStat(0);   
    c1.Update()

    if UseUnfold: PlotType="_Unfolded" 
    else: PlotType=""

    if UseMCReco: PlotType = PlotType+"_RecoMC"

    MCSetStr = ""
    if MCSet == "Pow":     MCSetStr = "Powheg"
    elif MCSet == "Mad":     MCSetStr = "MadGraph"

    if SavePlot:
        c1.SaveAs("Plot/CrossSection/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".png")        
        c1.SaveAs("Plot/CrossSection/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".eps")        
        c1.SaveAs("Plot/CrossSection/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".pdf")        
        c1.SaveAs("Plot/CrossSection/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".root")        
        
   
     



