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
from  PersonalInfo import*

parser = OptionParser(usage="usage: %prog [options]")


parser.add_option("-a", "--analysis", dest="analysis",
                  default="ZZ",
                  help="Analysis option, default is ZZ. Others are ZZFull and HZZ")

parser.add_option("-s", "--save", dest="SavePlot",
                  action= "store_true",
                  default=False,
                  help="Save plot option, default is False")

parser.add_option("-u", "--unfold", dest="UseUnfold",
                  action= "store_true",
                  default=False,
                  help="Use unfolded data or not, default is True")

parser.add_option("-r", "--recoMC", dest="UseMCReco",
                  action= "store_true",
                  default=False,
                  help="Use reconstructed MC samples or not, default is False")

parser.add_option("-t", "--Type", dest="Type",
                  default="Mass",
                  help="Type of differential Cross secion. ['Mass','Jets','Mjj','Deta',CentralJets,CentralMjj,CentralDeta]. Default is Mass")

parser.add_option("-i", "--Inclusive", dest="DoInclusive",
                  action= "store_true",
                  default=False,
                  help="Compute inclusive Cross secion. Default is False ")

parser.add_option("-S", "--Set", dest="MCSetIn",
                  default="Mad",
                  help="Choose MC Set between Pow (Powheg) and Mad (MadGraph). Default is Mad")

parser.add_option("-n", "--Normalized", dest="DoNormalized",
                  action= "store_true",
                  default=False,
                  help="Compute normalized cross section. Default is True")

parser.add_option("-f", "--Fiducial", dest="doFiducial",
                  action  = "store_true",
                  default = False,
                  help="Compute fiducial cross section. Default is False")

parser.add_option("-d", "--Dir", dest="Dir",
                  default = "test",
                  help="Choose a specific directory name for save plots. Default is ''")

(options, args) = parser.parse_args()


SavePlot     =  options.SavePlot
UseUnfold    =  options.UseUnfold
UseMCReco    =  options.UseMCReco
DoInclusive  =  options.DoInclusive
MCSetIn      =  options.MCSetIn
Type         =  options.Type
DoNormalized =  options.DoNormalized
DoFiducial   =  options.doFiducial
analysis     =  options.analysis
Dir          =  options.Dir

Dir+="/"

try:
    os.stat("./Plot")
except:
    os.mkdir("./Plot")
try:
    os.stat("./Plot/CrossSection")
except:
    os.mkdir("./Plot/CrossSection")

if SavePlot:
    try:
        os.stat(PersonalFolder+Dir)
    except:
        os.mkdir(PersonalFolder+Dir)
        os.system("cp "+PersonalFolder+"index.php "+PersonalFolder+Dir )

    try:
        os.stat(PersonalFolder+Dir+"/CrossSection/")
    except:
        os.mkdir(PersonalFolder+Dir+"/CrossSection/")
        os.system("cp "+PersonalFolder+"index.php "+PersonalFolder+Dir+"/CrossSection/" )


if DoInclusive:
    print Red("\n\n##############################################################\n################# ZZ4l Inclusive Cross Section ###############\n##############################################################\n")
else:  
    print Red("\n\n#################################################################################\n############## ZZ4l Differential Cross Section as function of {0} ##############\n#################################################################################\n".format(Type))




if DoInclusive:
    hMCList = getCrossPlot_MC(MCSetIn,"Total",analysis,DoNormalized,DoFiducial)
    TotalCross(MCSetIn,"Tot",analysis,DoFiducial,UseMCReco)
    sys.exit(0)

hMCList = getCrossPlot_MC(MCSetIn,Type,analysis,DoNormalized,DoFiducial)

if UseMCReco:
    (hDataList,hDataListUp,hDataListDown) = getCrossPlot_Data(MCSetIn,UseUnfold,Type,analysis,0,True,DoNormalized,DoFiducial) 

else:
    (hDataList,hDataListUp,hDataListDown)  = getCrossPlot_Data(MCSetIn,UseUnfold,Type,analysis,0,False,DoNormalized,DoFiducial) 



for i in range(0,4):
    hDataList[i]["state"].SetMarkerSize(1.2)
    hDataList[i]["state"].SetMarkerColor(1)
    hDataList[i]["state"].SetMarkerStyle(20)
    hDataList[i]["state"].SetLineColor(1)
    
    grSyst = getSystGraph(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],DoFiducial,hDataListDown[i]["name"])
    grSyst.SetFillStyle(3001)
    grSyst.SetFillColor(ROOT.kRed)
 
    Err=ROOT.Double(0.)
    
    LastBin = hMCList[i]["state"].GetXaxis().GetLast()    

    hRatio= ROOT.TH1F()
    hRatio = hDataList[i]["state"].Clone("hRatioNew")
    hRatio.SetStats(0)
    hRatio.Divide(hMCList[i]["state"])
   
    if Type=="Mass":
        xName = "M_{4l} [GeV]"
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("frac{1}{#sigma_{fid}}#frac{d#sigma_{fid}}{d M_{ZZ}} [1/GeV]")
        elif DoFiducial: hMCList[i]["state"].GetYaxis().SetTitle("#frac{d#sigma}{d M_{ZZ}} [fb/GeV]")
        else:  hMCList[i]["state"].GetYaxis().SetTitle("#frac{d#sigma}{d M_{ZZ}} [pb/GeV]")

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
        elif DoFiducial: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d #jets [fb]") 
        else: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d #jets [pb]") 
        WriteJetsNorm(hMCList,Type,DoFiducial) #Write partial xsec for variable with jet>0

    elif "Deta" in Type:
        xName = "#Delta #eta_{jj}"
        hRatio.GetXaxis().SetTitle(xName)
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}}#frac{d#sigma_{fid}}{d #Delta #eta_{jj}} [1/GeV]")
        elif DoFiducial: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d #Delta #eta_{jj}[fb]")  
        else:            hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d #Delta #eta_{jj}[pb]")  

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
        elif DoFiducial: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d p_{T}^{jet} [fb/GeV]") 
        else:  hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d p_{T}^{jet} [pb/GeV]") 

    elif "EtaJet" in Type:
        xName = "#eta"
        hRatio.GetXaxis().SetTitle("#eta")
        hMCList[i]["state"].GetXaxis().SetTitle("#eta")
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}}#frac{d#sigma_{fid}}{d #eta^{jet}} [1/GeV]")             
        elif DoFiducial: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d #eta^{jet} [fb/GeV]")
        else:            hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d #eta^{jet} [pb/GeV]")

    elif "dRZZ" in Type:
        xName = "#Delta R_{ZZ}"
        hRatio.GetXaxis().SetTitle("#Delta R_{ZZ}")
        hMCList[i]["state"].GetXaxis().SetTitle("#Delta R_{ZZ}")
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}}#frac{d#sigma_{fid}}{d #Delta r^{jets}} [1/GeV]")             
        elif DoFiducial: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d #Delta R_{ZZ} [fb/GeV]")
        else:            hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d #Delta R_{ZZ} [pb/GeV]")


    hRatio.GetXaxis().SetTitle(xName)
    
    hRatio.GetYaxis().SetTitle("Data/MC")
    

    # "Type":{fr:{norm,notnorm},notfr:{norm,notnorm}}
    MaxList = {"Mass":{True:{True:0.009,False:0.35},False:{True:0.009,False:0.15}},"Jets":{True:{True:0.8,False:46.4},False:{True:0.85,False:13}},"Jets_Central":{True:{True:.82,False:50},False:{True:0.85,False:13}},"Mjj":{True:{True:0.004,False:0.02},False:{True:0.004, False:0.0086}},"Mjj_Central":{True:{True:0.005,False:0.018},False:{True:0.005,False:0.006}},"Deta":{True:{True:0.6,False:2.2},False:{True:0.45,False:1.2}},"Deta_Central":{True:{True:0.6,False:1.98},False:{True:0.7,False:.9}} ,"PtJet1":{True:{True:0.03,False:0.48},False:{True:0.06,False:0.16}},"PtJet2":{True:{True:0.016,False:0.08},False:{True:0.016,False:0.025}},"EtaJet1":{True:{True:0.5,False:7.4},False:{True:0.6,False:2.8}},"EtaJet2":{True:{True:0.5,False:2.6},False:{True:0.58,False:0.9}},"dRZZ":{True:{True:0.65,False:23.6},False:{True:0.7,False:12.}}}

    # Max   = hDataList[i]["state"].GetMaximum()+hDataList[i]["state"].GetBinError(hDataList[i]["state"].GetMaximumBin()) #Automatic max
    Max = MaxList[Type][DoFiducial][DoNormalized]

    hMCList[i]["state"].SetMaximum(Max)
    hMCList[i]["state"].SetMinimum(0)
    hMCList[i]["state"].SetFillColor(ROOT.kCyan-3)
    hMCList[i]["state"].SetLineColor(1)
    hMCList[i]["state"].GetXaxis().SetTitle(xName)
    hMCList[i]["state"].SetTitle(hMCList[i]["name"]+" Cross Section")
    
    c1 = TCanvas( 'c1', hMCList[i]["name"] , 200, 10, 900, 700 )
  

    leg = ROOT.TLegend(.65,.52,.85,.66);

    leg.AddEntry(grSyst,"Syst","f")
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
        grSyst.Draw("same2")
    hDataList[i]["state"].Draw("sameE1")
    leg.Draw("same")
    
    pad2.cd()
    hRatio.Draw("E1")
    fInt.Draw("same")
    
    ROOT.gStyle.SetOptStat(0);   
    c1.SetLogy()    

    if UseUnfold: PlotType="_Unfolded" 
    else: PlotType=""

    if UseMCReco: PlotType = PlotType+"_RecoMC"

    MCSetStr = ""
    if   MCSetIn == "Pow":        MCSetStr = "Powheg"
    elif MCSetIn == "Mad":        MCSetStr = "MadGraph"
    elif MCSetIn == "fr_Mad":     MCSetStr = "fr_MadGraph"
    elif MCSetIn == "fr_Pow":     MCSetStr = "fr_Powheg"


    Kind = ""
    if DoFiducial:   Kind += "_tight"
    if DoNormalized: Kind += "_norm"

    if SavePlot:
        c1.SaveAs("Plot/CrossSection/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+Kind+".png")
        c1.SaveAs("Plot/CrossSection/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+Kind+".root")
        c1.SaveAs(PersonalFolder+Dir+"CrossSection/diffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+Kind+".png")        
        c1.SaveAs(PersonalFolder+Dir+"CrossSection/diffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+Kind+".eps")        
        c1.SaveAs(PersonalFolder+Dir+"CrossSection/diffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+Kind+".pdf")        
        c1.SaveAs(PersonalFolder+Dir+"CrossSection/diffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+Kind+".root")       
   
     
