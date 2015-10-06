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
                  help="Type of differential Cross secion. ['Mass','Jets','Mjj','Deta',CentralJets,CentralMjj,CentralDeta]. Default is Mass")

parser.add_option("-i", "--Inclusive", dest="DoInclusive",
                  default="False",
                  help="Compute inclusive Cross secion. Default is False ")

parser.add_option("-m", "--MCSet", dest="MCSetIn",
                  default="Pow",
                  help="Choose MC Set between Pow (Powheg) and Mad (MadGraph). Default is Pow")

parser.add_option("-n", "--Normalized", dest="DoNormalized",
                  default="True",
                  help="Compute normilized cross section. Default is True")

(options, args) = parser.parse_args()


SavePlot  =ast.literal_eval(options.SavePlot)
UseUnfold =ast.literal_eval(options.UseUnfold)
UseMCReco =ast.literal_eval(options.UseMCReco)
DoInclusive =ast.literal_eval(options.DoInclusive)
MCSetIn = options.MCSetIn
Type = options.Type
DoNormalized =ast.literal_eval(options.DoNormalized)

try:
    os.stat("./Plot")
except:
    os.mkdir("./Plot")
try:
    os.stat("./Plot/CrossSection")
except:
    os.mkdir("./Plot/CrossSection")

if DoInclusive:
    print Red("\n\n##############################################################\n################# ZZ4l Inclusive Cross Section ###############\n##############################################################\n")
else:  
    print Red("\n\n#################################################################################\n############## ZZ4l Differential Cross Section as function of {0} ##############\n#################################################################################\n".format(Type))


#SetGlobalVariable(MCSet) 
#CrossInfo.MCSet=MCSetIn


#import CrossInfo
#from CrossInfo import*

hMCList = getCrossPlot_MC(MCSetIn,Type,DoNormalized)

if DoInclusive:
    TotalCross(MCSetIn,Type)
    sys.exit(0)

if UseMCReco:
    (hDataList,hDataListUp,hDataListDown) = getCrossPlot_Data(MCSetIn,UseUnfold,Type,0,True,DoNormalized) 

else:
    (hDataList,hDataListUp,hDataListDown)  = getCrossPlot_Data(MCSetIn,UseUnfold,Type,0,False,DoNormalized) 


for i in range(0,4):
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
   
    if Type=="Mass":
        xName = "M_{4l} [GeV]"
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("frac{1}{#sigma_{fid}}#frac{d#sigma_{fid}}{d M_{ZZ}} [1/GeV]")
        else: hMCList[i]["state"].GetYaxis().SetTitle("#frac{d#sigma}{d M_{ZZ}} [pb/GeV]")

    elif "Jets" in Type:
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
        
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}}#frac{d#sigma_{fid}}{d #jets}")        
        else: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d #jets [pb]") 

    elif "Deta" in Type:
        xName = "#Delta #eta_{jj}"
        hRatio.GetXaxis().SetTitle(xName)
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}}#frac{d#sigma_{fid}}{d #Delta #eta_{jj}} [1/GeV]")
        else: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d #Delta #eta_{jj}[pb]")  

    elif "Mjj" in Type:
        xName = "m_{jj} [Gev|"
        hRatio.GetXaxis().SetTitle("#Jets")
        hMCList[i]["state"].GetXaxis().SetTitle("m_{jj}")
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}}#frac{d#sigma_{fid}}{d  m_{jj}} [1/GeV]")      
        else: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d m_{jj} [pb/GeV]")  

    elif "PtJet" in Type:
        xName = "p_{T} [Gev|"
        hRatio.GetXaxis().SetTitle("#Jets")
        hMCList[i]["state"].GetXaxis().SetTitle("p_{T}")
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{d p_{T}^{jet}} [1/GeV]")             
        else:  hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d p_{T}^{jet} [pb/GeV]") 

    elif "EtaJet" in Type:
        xName = "#eta [Gev|"
        hRatio.GetXaxis().SetTitle("#Jets")
        hMCList[i]["state"].GetXaxis().SetTitle("#eta")
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}}#frac{d#sigma_{fid}}{d #eta^{jet}} [1/GeV]")             
        else: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d #eta^{jet} [pb/GeV]")

    hRatio.GetXaxis().SetTitle(xName)
    
    hRatio.GetYaxis().SetTitle("Data/MC")
    
    Max= hDataList[i]["state"].GetMaximum()+hDataList[i]["state"].GetBinError(hDataList[i]["state"].GetMaximumBin())
    MaxMC=hMCList[i]["state"].GetMaximum()


    if Type == "Mass": 
        if DoNormalized:
            if "fr" in MCSetIn: Max=0.011 #Max=Max+Max/5.
            else: Max = 0.011
        else: 
            if "fr" in MCSetIn: Max=0.0002 #Max=Max+Max/5.
            else: Max = 0.1
    elif "Mjj" in Type: 
        if DoNormalized:
            Max=0.0060 
    elif "Deta" in Type: 
        if DoNormalized:
            Max=1.
        else: Max=0.4
        #hRatio.SetMinimum(hRatio.GetMinimum()-hRatio.GetMinimum()/8.)
        #hRatio.SetMinimum(0.)
        #hRatio.SetMaximum(hRatio.GetMaximum()+hRatio.GetMaximum()/8.)
        #hRatio.SetMaximum(1.5)
        
    elif Type=="PtJet1":
        if DoNormalized: Max=0.04
        else: Max=0.1
        
    elif Type=="PtJet2":
        if DoNormalized: Max=0.020
        else: Max=0.015

    elif Type=="EtaJet1":
        Max=2.

    elif Type=="EtaJet2":
        if DoNormalized:Max = 1.
        else: Max=0.50

    else: Max=Max+Max/9.

    
    hMCList[i]["state"].SetMaximum(Max)
    hMCList[i]["state"].SetMinimum(0)
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
#        hDataList[i]["state"].Draw("E1")
    leg.Draw("same")
    
    pad2.cd()
    hRatio.Draw("E1")
    fInt.Draw("same")
    
    ROOT.gStyle.SetOptStat(0);   
    c1.SetLogy()    
    #c1.Update()

    if UseUnfold: PlotType="_Unfolded" 
    else: PlotType=""

    if UseMCReco: PlotType = PlotType+"_RecoMC"

    MCSetStr = ""
    if MCSetIn == "Pow":     MCSetStr = "Powheg"
    elif MCSetIn == "Mad":     MCSetStr = "MadGraph"
    elif MCSetIn == "fr_Mad":     MCSetStr = "fr_MadGraph"
    elif MCSetIn == "fr_Pow":     MCSetStr = "fr_Powheg"

    if SavePlot:
        c1.SaveAs("Plot/CrossSection/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".png")        
        c1.SaveAs("Plot/CrossSection/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".eps")        
        c1.SaveAs("Plot/CrossSection/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".pdf")        
        c1.SaveAs("Plot/CrossSection/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".root")        

        if DoNormalized:
            c1.SaveAs("~/www/VBS/CrossSections/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+"_norm.png")        
            c1.SaveAs("~/www/VBS/CrossSections/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+"_norm.eps")        
            c1.SaveAs("~/www/VBS/CrossSections/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+"_norm.pdf")        
            c1.SaveAs("~/www/VBS/CrossSections/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+"_norm.root")       
        else:
            c1.SaveAs("~/www/VBS/CrossSections/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".png")        
            c1.SaveAs("~/www/VBS/CrossSections/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".eps")        
            c1.SaveAs("~/www/VBS/CrossSections/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".pdf")        
            c1.SaveAs("~/www/VBS/CrossSections/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".root")  
   
     
