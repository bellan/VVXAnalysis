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
from LL import*
import sys,ast,os
import math
import operator
from  PersonalInfo import*
import CMS_lumi, tdrstyle
import LatexUtils
from LatexUtils import*
import PlotStyle
from PlotStyle import*

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
                  help="Type of differential Cross secion. ['Mass','nJets','Mjj','Deta',CentralnJets,CentralMjj,CentralDeta]. Default is Mass")

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

parser.add_option("-L", "--printLatex", dest="printLatex",
                  action="store_true",
                  default=False,
                  help="Print tabular LateX code")

parser.add_option("-l", "--dolog", dest="doLog",
                  action="store_true",
                  default=False,
                  help="plot in log scale")

parser.add_option("-p", "--plot", dest="showPlot",
                  action="store_true",
                  default=False,
                  help="Show plot option, default is False")

parser.add_option("-R", "--Ratio8TeV", dest="do8TeVRatio",
                  action="store_true",
                  default=False,
                  help="8 TeV ratio ange")


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
printLatex   =  options.printLatex
doLog        =  options.doLog
showPlot     =  options.showPlot
do8TeVRatio  =  options.do8TeVRatio
Dir+="/"

if not showPlot: ROOT.gROOT.SetBatch(True)

if printLatex or DoNormalized: 
    LumiErr = GlobSystList[1]["value"] 
    GlobSystList[1]["value"]=0.

if DoFiducial: GlobSystList+=PdfSyst_fid
else: GlobSystList+=PdfSyst


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
        os.stat(PersonalFolder+Dir+"/"+Type)
    except:
        os.mkdir(PersonalFolder+Dir+"/"+Type)
        os.system("cp "+PersonalFolder+"index.php "+PersonalFolder+Dir+"/"+Type )
    try:
        os.stat(PersonalFolder+Dir+"/"+Type+"/CrossSection/")
    except:
        os.mkdir(PersonalFolder+Dir+"/"+Type+"/CrossSection/")
        os.system("cp "+PersonalFolder+"index.php "+PersonalFolder+Dir+"/"+Type+"/CrossSection/" )


if DoInclusive:
    print Red("\n\n##############################################################\n################# ZZ4l Inclusive Cross Section ###############\n##############################################################\n")
else:  
    print Red("\n\n#################################################################################\n############## ZZ4l Differential Cross Section as function of {0} ##############\n#################################################################################\n".format(Type))

#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetHatchesSpacing(hatchSpace)
ROOT.gStyle.SetHatchesLineWidth(1)                                                                                  
       
#ROOT.gStyle.SetErrorX(0)                                  

lumi = round(Lumi/1000.,1)
CMS_lumi.writeExtraText = True
#CMS_lumi.extraText = "Preliminary"
extratext1 = "anti-k_{T} (R=0.4) Jets"
if "Central"  not in Type:  extratext2 = "p_{T}^{j}>30 GeV, |#eta_{j}|<4.7"
else: extratext2 = "p_{T}^{j}>30 GeV, |#eta_{j}|<2.4"

#extratext3 = 
CMS_lumi.extraText = "#splitline{"+extratext1+"}{"+extratext2+"}"
#CMS_lumi.extraText = "#splitline{a}{b}"
CMS_lumi.lumi_sqrtS = "{0}".format(lumi)+" fb^{-1} (13 TeV)\n"
CMS_lumi.cmsTextSize = 1.
CMS_lumi.relPosY    = 0.78

if "Mjj" in Type:  
   # CMS_lumi.relPosX = 0.63
   # CMS_lumi.relPosY = 0.71
    CMS_lumi.relPosX = 0.55   
    CMS_lumi.relPosY = 0.55
else:
    CMS_lumi.relPosY    = 0.78

iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 0


if DoInclusive:
    (hMCList,a,b) = getCrossPlot_MC(MCSetIn,"Total",analysis,DoNormalized,DoFiducial)
    TotalCross(MCSetIn,"Tot",analysis,DoFiducial,UseMCReco)
    sys.exit(0)

if MCSetIn   == "Pow": mcop = "Mad"
elif MCSetIn == "Mad" : mcop = "Pow"

(hMCList,hMCListUp,hMCListDown)          = getCrossPlot_MC(MCSetIn,Type,analysis,DoNormalized,DoFiducial)
(hMCList_op,hMCListUp_op,hMCListDown_op) = getCrossPlot_MC(mcop,Type,analysis,DoNormalized,DoFiducial)


if UseMCReco:
    (hDataList,hDataListUp,hDataListDown) = getCrossPlot_Data(MCSetIn,UseUnfold,Type,analysis,0,True,DoNormalized,DoFiducial) 

else:
    (hDataList,hDataListUp,hDataListDown)  = getCrossPlot_Data(MCSetIn,UseUnfold,Type,analysis,0,False,DoNormalized,DoFiducial) 

if Type=="nJets" or Type=="nIncJets":
    WriteJetsNorm(hDataList,Type,DoFiducial) #Write partial xsec for variable with jet>0    
    CrossDic = collections.OrderedDict()

    print "Cross Section per jet multiplicity"
    for h,hup,hdown,hmc,hmcup,hmcdown in zip(hDataList,hDataListUp,hDataListDown,hMCList,hMCListUp,hMCListDown):
        print Red(h["name"])
        for bin in range(1,h["state"].GetNbinsX()+1):
            print "{0} jets:  {1:.2f} +- {2:.2f}(stat) + {3:.2f} - {4:.2f} (syst)\n".format(bin,h["state"].GetBinContent(bin),h["state"].GetBinError(bin),hup["state"].GetBinContent(bin)-h["state"].GetBinContent(bin),h["state"].GetBinContent(bin)-hdown["state"].GetBinContent(bin))
            if h["name"]=="4l" and printLatex:
                CrossDic[bin]={}  
                CrossDic[bin]["cross"]="{0:.3f}".format(h["state"].GetBinContent(bin))
                CrossDic[bin]["stat"]="{0:.3f}".format(h["state"].GetBinError(bin))
                CrossDic[bin]["up"]="{0:.3f}".format(hup["state"].GetBinContent(bin)-h["state"].GetBinContent(bin))
                CrossDic[bin]["down"]="{0:.3f}".format(h["state"].GetBinContent(bin)-hdown["state"].GetBinContent(bin))
                CrossDic[bin]["lumi"]="{0:.3f}".format(h["state"].GetBinContent(bin)*LumiErr)
                CrossDic[bin]["theor"]="{0:.3f}".format(hmc["state"].GetBinContent(bin))
                CrossDic[bin]["theorUp"]="{0:.3f}".format(hmcup["state"].GetBinContent(bin)-hmc["state"].GetBinContent(bin))
                CrossDic[bin]["theorDn"]="{0:.3f}".format(hmc["state"].GetBinContent(bin)-hmcdown["state"].GetBinContent(bin))

for i in range(0,4):
    hDataList[i]["state"].SetMarkerSize(1.3)
    hDataList[i]["state"].SetMarkerColor(1)
    hDataList[i]["state"].SetMarkerStyle(20)
    hDataList[i]["state"].SetLineColor(1)

    if not DoFiducial and hDataList[i]["name"]=="4l":
        grSyst = getUncGraph(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],MCSetIn,hDataList[i]["name"],True,"syst",Type)
    else:
        print "Hei ",MCSetIn,hDataList[i]["name"]
        grSyst = getUncGraph(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],MCSetIn,hDataList[i]["name"],True,"statsyst",Type)

    grMCSyst    = getUncGraph(hMCList[i]["state"],hMCListUp[i]["state"],hMCListDown[i]["state"],MCSetIn,hMCList[i]["name"],False,"syst",Type)    
    grMCSyst_op = getUncGraph(hMCList_op[i]["state"],hMCListUp_op[i]["state"],hMCListDown_op[i]["state"],MCSetIn,hMCList_op[i]["name"],False,"syst",Type)

    grSyst.SetFillStyle(DataFillStyle)
    grSyst.SetFillColorAlpha(DataFillColor,DataAlpha)
    grSyst.SetLineColorAlpha(DataFillColor,DataAlpha-corrCont)

    grMCSyst.SetFillColorAlpha(MG5FillColor,MG5Alpha) 
    grMCSyst.SetFillStyle(MG5FillStyle)
    grMCSyst.SetLineColor(MG5Color)
    grMCSyst.SetMarkerColor(MG5Color)
    grMCSyst.SetLineStyle(MG5LineStyle)
    grMCSyst.SetLineWidth(MG5LineWidth)

    grMCSyst_op.SetFillColorAlpha(PowFillColor,PowAlpha)
    grMCSyst_op.SetLineColor(PowColor)
    grMCSyst_op.SetFillStyle(PowFillStyle)
    grMCSyst_op.SetLineStyle(PowLineStyle)
    grMCSyst_op.SetMarkerColor(PowColor)
    grMCSyst_op.SetLineWidth(PowLineStyle)
    grMCSyst_op.SetLineWidth(PowLineWidth)

    Err=ROOT.Double(0.)    
    LastBin = hMCList[i]["state"].GetXaxis().GetLast()    

    if "Pow" in MCSetIn: 
        mcref_entry = "Powheg+MCFM+Pythia8"
#        mcop_entry  = "MadGraph5_aMC@NLO+MCFM+Pythia8"
        mcop_entry  = "MG5_aMC@NLO+MCFM+Pythia8"
    else:
        mcref_entry = "MG5_aMC@NLO+MCFM+Pythia8"
#        mcref_entry = "MG5_aMC@NLO+MCFM+Pythia8"
        mcop_entry  = "Powheg+MCFM+Pythia8"
    
###########################################################
   
    #Ratio and uncertainties on data points for the reference set of samples 
          
    #statistical uncertainty

    (hRatio,hRatio_up,hRatio_down) = LL(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],hMCList[i]["state"],Type,hDataList[i]["name"],"stat",True)     
    grStat_ratio = getUncGraph(hRatio,hRatio_up,hRatio_down,MCSetIn,hDataList[i]["name"],True,"ratio",Type)  

    #systematic uncertainty
    (hRatio_ss,hRatio_up_ss,hRatio_down_ss) = LL(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],hMCList[i]["state"],Type,hDataList[i]["name"],"statsyst",True) 

    print "ratio ss up",hRatio_up_ss.GetBinContent(2) 
    grSyst_ratio = getUncGraph(hRatio_ss,hRatio_up_ss,hRatio_down_ss,MCSetIn,hDataList[i]["name"],True,"syst",Type)
    grSyst_ratio.SetFillStyle(DataFillStyle)
    grSyst_ratio.SetFillColorAlpha(DataFillColor,DataAlpha-corrMad)

    (hMCRatio_ss,hMCRatio_up_ss,hMCRatio_down_ss) = LL(hMCList[i]["state"],hMCListUp[i]["state"],hMCListDown[i]["state"],hMCList[i]["state"],Type,hMCList[i]["name"],"syst",True) 
    grMCSyst_ratio = getUncGraph(hMCRatio_ss,hMCRatio_up_ss,hMCRatio_down_ss,MCSetIn,hMCList[i]["name"],False,"syst",Type)  
    grMCSyst_ratio.SetFillColorAlpha(MG5Color,MG5Alpha)
    grMCSyst_ratio.SetFillStyle(MG5FillStyle)

    (hMCRatio_op_ss,hMCRatio_op_up_ss,hMCRatio_op_down_ss) = LL(hMCList_op[i]["state"],hMCListUp_op[i]["state"],hMCListDown_op[i]["state"],hMCList_op[i]["state"],Type,hMCList_op[i]["name"],"syst",True) 
    grMCSyst_ratio_op = getUncGraph(hMCRatio_op_ss,hMCRatio_op_up_ss,hMCRatio_op_down_ss,MCSetIn,hMCList_op[i]["name"],False,"syst",Type)  

    grMCSyst_ratio_op.SetFillColorAlpha(PowFillColor,PowAlpha)
    grMCSyst_ratio_op.SetFillStyle(PowFillStyle)

    maxtot = hRatio_up_ss.GetMaximum()   


###################################################################

    #Ratio and uncertainties on data for the opposite set of samples 
 
    #statistical uncertainty
    (hRatio_op,hRatio_op_up,hRatio_op_down) = LL(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],hMCList_op[i]["state"],Type,hDataList[i]["name"],"stat",True)     
    grStat_ratio_op = getUncGraph(hRatio_op,hRatio_op_up,hRatio_op_down,MCSetIn,hDataList[i]["name"],True,"ratio",Type)   
  
    #systematic uncertainty
    (hRatio_op_ss,hRatio_op_up_ss,hRatio_op_down_ss) = LL(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],hMCList_op[i]["state"],Type,hDataList[i]["name"],"statsyst",True)  
    grSyst_ratio_op = getUncGraph(hRatio_op_ss,hRatio_op_up_ss,hRatio_op_down_ss,MCSetIn,hDataList[i]["name"],True,"syst",Type)  

    grSyst_ratio_op.SetFillStyle(DataFillStyle)
    grSyst_ratio_op.SetFillColorAlpha(DataFillColor,DataAlpha-corrPow)

    maxtot_op =  hRatio_op_up_ss.GetMaximum() 
    maxratios_op = max(maxtot,maxtot_op)


    FidStr = "" #"_{fid}"
    FidStrJets = "" # "fid,"

    if DoFiducial: Unit = "fb"
    else:          Unit = "pb"

    if "Central" in Type: eta_cut = " (|#eta_{j}| < 2.4)"
    else: eta_cut = " (|#eta_{j}| < 4.7)"

    if Type=="Mass":
        xName = "m_{4 l} [GeV]"
        hRatio.GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma"+FidStr+"} #frac{d#sigma"+FidStr+"}{dm_{4l}} [1/GeV]")
        else: hMCList[i]["state"].GetYaxis().SetTitle("#frac{d#sigma}{dm_{4l}} ["+Unit+"/GeV]")

    elif Type=="PtZZ":
        xName = "p_{T}^{ZZ} [GeV]"
        hRatio.GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma"+FidStr+"} #frac{d#sigma"+FidStr+"}{dp_{T}^{4l}} [1/GeV]")
        else: hMCList[i]["state"].GetYaxis().SetTitle("#frac{d#sigma}{dp_{T}^{4l}} ["+Unit+"/GeV]")

    elif "nJets" in Type:
        hRatio.GetXaxis().SetBinLabel(1,"0")
        hRatio.GetXaxis().SetBinLabel(2,"1")
        hRatio.GetXaxis().SetBinLabel(3,"2")
        hRatio.GetXaxis().SetBinLabel(4,"#geq 3")

        hRatio_op.GetXaxis().SetBinLabel(1,"0")
        hRatio_op.GetXaxis().SetBinLabel(2,"1")
        hRatio_op.GetXaxis().SetBinLabel(3,"2")
        hRatio_op.GetXaxis().SetBinLabel(4,"#geq 3")

        hMCList[i]["state"].GetXaxis().SetBinLabel(1,"0")
        hMCList[i]["state"].GetXaxis().SetBinLabel(2,"1")
        hMCList[i]["state"].GetXaxis().SetBinLabel(3,"2")
        hMCList[i]["state"].GetXaxis().SetBinLabel(4,"#geq 3")
        
        xName = "N_{jets}"#+eta_cut #fixme
       
        hRatio.GetXaxis().SetTitle(xName)
        hMCList[i]["state"].GetXaxis().SetTitle(xName)

        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma"+FidStr+"} #frac{d#sigma"+FidStr+"}{dN_{jets}}")        
        else: hMCList[i]["state"].GetYaxis().SetTitle("#frac{d#sigma"+FidStr+"}{dN_{jets}} ["+Unit+"]")  

    elif "nIncJets" in Type:
        hRatio_op.GetXaxis().SetBinLabel(1,"#geq0")
        hRatio_op.GetXaxis().SetBinLabel(2,"#geq1")
        hRatio_op.GetXaxis().SetBinLabel(3,"#geq2")
        hRatio_op.GetXaxis().SetBinLabel(4,"#geq3")

        hMCList[i]["state"].GetXaxis().SetBinLabel(1,"#geq0")
        hMCList[i]["state"].GetXaxis().SetBinLabel(2,"#geq1")
        hMCList[i]["state"].GetXaxis().SetBinLabel(3,"#geq2")
        hMCList[i]["state"].GetXaxis().SetBinLabel(4,"#geq3")
        
        xName = "N_{jets}"#+eta_cut
       
        hRatio.GetXaxis().SetTitle(xName)
        hMCList[i]["state"].GetXaxis().SetTitle(xName)

        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma"+FidStr+"} #frac{d#sigma"+FidStr+"}{dN_{jets}}")        
        else: hMCList[i]["state"].GetYaxis().SetTitle("#frac{d#sigma"+FidStr+"}{dN_{jets}}["+Unit+"]")  

    elif "Deta" in Type:
        xName = "#Delta#eta_{jj}" #+ eta_cut
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{"+FidStrJets+"jets#geq2}} #frac{d#sigma"+FidStr+"}{d#Delta#eta_{jj}}")
        else: hMCList[i]["state"].GetYaxis().SetTitle("#frac{d#sigma"+FidStr+"}{d#Delta#eta_{jj}}["+Unit+"]")  

    elif "Mjj" in Type:
        xName = "m_{jj}" + " [GeV]" #+ eta_cut
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized:
            hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{"+FidStrJets+"jets#geq2}} #frac{d#sigma"+FidStr+"}{dm_{jj}} [1/GeV]")      
            #hMCList[i]["state"].GetYaxis().SetMaxDigits(1);
        else: hMCList[i]["state"].GetYaxis().SetTitle("#frac{d#sigma"+FidStr+"}{dm_{jj}} ["+Unit+"/GeV]")  

    elif "PtJet1" in Type:
        xName = "p_{T}^{j1} [GeV]"
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{"+FidStrJets+"jets#geq1}} #frac{d#sigma"+FidStr+"}{dp_{T}^{j1}} [1/GeV]")             
        else:  hMCList[i]["state"].GetYaxis().SetTitle("#frac{d#sigma"+FidStr+"}{dp_{T}^{j1}} ["+Unit+"/GeV]") 
    elif "PtJet2" in Type:
        xName = "p_{T}^{j2} [GeV]"
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{"+FidStrJets+"jets#geq2}} #frac{d#sigma"+FidStr+"}{dp_{T}^{j2}} [1/GeV]")             
        else:  hMCList[i]["state"].GetYaxis().SetTitle("#frac{d#sigma"+FidStr+"}{dp_{T}^{j2}} ["+Unit+"/GeV]") 
    elif "EtaJet1" in Type:
        xName = "|#eta_{j1}|"
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{"+FidStrJets+"jets#geq1}} #frac{d#sigma"+FidStr+"}{d|#eta_{j1}|}")             
        else: hMCList[i]["state"].GetYaxis().SetTitle("#frac{d#sigma"+FidStr+"}{d|#eta_{j1}|} ["+Unit+"]") 
    elif "EtaJet2" in Type:
        xName = "|#eta_{j2}|"
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{"+FidStrJets+"jets#geq2}} #frac{d#sigma"+FidStr+"}{d|#eta_{j2}|}")             
        else: hMCList[i]["state"].GetYaxis().SetTitle("#frac{d#sigma"+FidStr+"}{d|#eta_{j2}|} ["+Unit+"]")
    elif "dRZZ" in Type:
        xName = "#Delta R_{ZZ}"
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma"+FidStr+"} #frac{d#sigma"+FidStr+"}{d #Delta R_{ZZ}}")             
        else: hMCList[i]["state"].GetYaxis().SetTitle("#frac{d#sigma"+FidStr+"}{d #Delta R_{ZZ}} ["+Unit+"]")


    hRatio.GetXaxis().SetTitle(xName) 
    hRatio_op.GetXaxis().SetTitle(xName)
    ratioName = "Data/MC"
    hRatio.GetYaxis().SetTitle(ratioName) 
    hRatio_op.GetYaxis().SetTitle(ratioName)
 
    # "Type":{fr:{norm,notnorm},notfr:{norm,notnorm}}
    MaxList = {"Mass":{True:{True:0.011,False:0.42},False:{True:0.011,False:0.17}},"nJets":{True:{True:0.97,False:40.},False:{True:1.,False:18.}},"nIncJets":{True:{True:0.97,False:60.},False:{True:1.,False:18.}},"nJets_Central":{True:{True:.98,False:55},False:{True:0.99,False:16}},"Mjj":{True:{True:0.0036,False:0.028},False:{True:0.006, False:0.0086}},"Mjj_Central":{True:{True:0.005,False:0.018},False:{True:0.005,False:0.006}},"Deta":{True:{True:0.6,False:3.2},False:{True:0.72,False:1.25}},"Deta_Central":{True:{True:0.8,False:1.98},False:{True:0.7,False:.9}} ,"PtJet1":{True:{True:0.028,False:0.4},False:{True:0.04,False:0.16}},"PtJet2":{True:{True:0.038,False:0.21},False:{True:0.036,False:0.145}},"EtaJet1":{True:{True:0.64,False:9.2},False:{True:0.65,False:3.2}},"EtaJet2":{True:{True:0.55,False:2.9},False:{True:0.58,False:1.2}},"dRZZ":{True:{True:0.65,False:29.6},False:{True:0.7,False:12.}},"PtZZ":{True:{True:0.03,False:1.2},False:{True:0.032,False:.55}}}

    MaxRatioList_8TeV = {"Mass":{True:{True:1.9,False:0.2},False:{True:1.9,False:0.01}},"nJets":{True:{True:2.3,False:-0.6},False:{True:2.5,False:-0.55}},"nJets_Central":{True:{True:2.42,False:-0.41},False:{True:0.92,False:16}},"Mjj":{True:{True:2.2,False:0.03},False:{True:3.2, False:0.00}},"Mjj_Central":{True:{True:1.5,False:0.02},False:{True:0.005,False:0.006}},"Deta":{True:{True:2.05,False:0.05},False:{True:2.9,False:0}},"Deta_Central":{True:{True:1.8,False:0.01},False:{True:2.2,False:1e-3}} ,"PtJet1":{True:{True:9.5,False:-2.},False:{True:2.5,False:1e-4}},"PtJet2":{True:{True:3.4,False:0.03},False:{True:5.,False:0.1}},"EtaJet1":{True:{True:2.4,False:-0.3},False:{True:1.6,False:3e-2}},"EtaJet2":{True:{True:1.78,False:.2},False:{True:1.7,False:0.2}},"dRZZ":{True:{True:0.65,False:29.6},False:{True:0.7,False:12.}},"PtZZ":{True:{True:1.5,False:0.3},False:{True:1.9,False:0.2}}}

    # "Type":{norm:{up,notup},notnorm:{up,notup}}
    MaxRatioList = {"Mass":{True:{True:1.9,False:0.2},False:{True:1.9,False:0.01}},"nJets":{True:{True:1.6,False:.4},False:{True:1.65,False:0.3}},"nIncJets":{True:{True:1.5,False:.45},False:{True:1.65,False:0.3}},"nJets_Central":{True:{True:1.6,False:0.6},False:{True:1.7,False:5.2e-1}},"Mjj":{True:{True:2.3,False:0.03},False:{True:3.2, False:0.00}},"Mjj_Central":{True:{True:1.5,False:0.02},False:{True:0.005,False:0.006}},"Deta":{True:{True:2.2,False:-0.2},False:{True:2.9,False:-2}},"Deta_Central":{True:{True:1.8,False:0.01},False:{True:2.2,False:1e-3}} ,"PtJet1":{True:{True:2.8,False:0.19},False:{True:3.2,False:1e-4}},"PtJet2":{True:{True:3.6,False:0.32},False:{True:6.,False:0.1}},"EtaJet1":{True:{True:1.6,False:.1},False:{True:1.9,False:5e-5}},"EtaJet2":{True:{True:1.78,False:.2},False:{True:1.7,False:0.2}},"dRZZ":{True:{True:0.65,False:29.6},False:{True:0.7,False:12.}},"PtZZ":{True:{True:1.5,False:0.3},False:{True:1.9,False:0.2}}}

    # "Type":{not norm:{max,min},norm:{max,min}}
    LogList = {"Mass":{True:{True:0.1,False:2e-5},False:{True:2.5,False:0.0009}},"nJets":{True:{True:2.1,False:0.02},False:{True:.9e2,False:0.8}},"nIncJets":{True:{True:3.2,False:0.02},False:{True:4e2,False:0.5}},"nJets_Central":{True:{True:3.2,False:1e-2},False:{True:110.9,False:.5}},"Mjj":{True:{True:0.02,False:5e-5},False:{True:0.004, False:0.0086}},"Mjj_Central":{True:{True:0.005,False:0.018},False:{True:0.005,False:0.006}},"Deta":{True:{True:12.,False:5.5e-3},False:{True:45.45,False:3e-2}},"Deta_Central":{True:{True:20.5,False:0.004},False:{True:0.7,False:.9}} ,"PtJet1":{True:{True:.05,False:4e-5},False:{True:1.1,False:1e-3}},"PtJet2":{True:{True:0.09,False:0.0001},False:{True:.96,False:2.5e-5}},"EtaJet1":{True:{True:31.,False:3e-3},False:{True:2e2,False:7e-2}},"EtaJet2":{True:{True:0.65,False:2.},False:{True:0.58,False:0.9}},"dRZZ":{True:{True:0.65,False:29.6},False:{True:0.7,False:12.}},"PtZZ":{True:{True:0.03,False:1.2},False:{True:0.02,False:.30}}}

    # "Type":{Norm:{max,min},Not norm:{max,min}}
    LogList_wide = {"Mass":{True:{True:0.1,False:2e-5},False:{True:2.5,False:0.0009}},"nJets":{True:{True:5.,False:1.2e-2},False:{True:135.,False:1e-1}},"nIncJets":{True:{True:5.,False:1.2e-2},False:{True:135.,False:1e-1}},"nJets_Central":{True:{True:9.2,False:1e-2},False:{True:8e+1,False:.1}},"Mjj":{True:{True:0.005,False:0.028},False:{True:0.004, False:0.0086}},"Mjj_Central":{True:{True:0.005,False:0.018},False:{True:0.005,False:0.006}},"Deta":{True:{True:12.,False:5.5e-3},False:{True:45.45,False:3e-2}},"Deta_Central":{True:{True:20.5,False:0.004},False:{True:0.7,False:.9}} ,"PtJet1":{True:{True:1.8,False:4e-5},False:{True:1.5,False:1e-4}},"PtJet2":{True:{True:1.7,False:0.00001},False:{True:2.1,False:1e-4}},"EtaJet1":{True:{True:31.,False:3e-3},False:{True:2e2,False:7e-2}},"EtaJet2":{True:{True:0.65,False:2.},False:{True:0.58,False:0.9}},"dRZZ":{True:{True:0.65,False:29.6},False:{True:0.7,False:12}},"PtZZ":{True:{True:0.03,False:1.2},False:{True:0.02,False:.30}}}

    # Maxratios = max(maxratios_op,maxtot_MGatNLO)+0.80*max(maxratios_op,maxtot_MGatNLO)
    # Max   = hDataList[i]["state"].GetMaximum()+hDataList[i]["state"].GetBinError(hDataList[i]["state"].GetMaximumBin()) #Automatic max
    Max = MaxList[Type][DoFiducial][DoNormalized]
    if doLog:
        if DoFiducial:  Max = LogList[Type][DoNormalized][True]
        else:   Max = LogList_wide[Type][DoNormalized][True]
    hMCList[i]["state"].SetMaximum(Max)
    hMCList[i]["state"].SetMinimum(0) 
    hMCList[i]["state"].SetFillColorAlpha(MadFillColor,MadAlpha)
    hMCList[i]["state"].SetLineStyle(MadLineStyle)
    hMCList[i]["state"].SetLineStyle(MG5LineStyle)
    hMCList[i]["state"].SetLineWidth(MG5LineWidth) 
    hMCList[i]["state"].GetXaxis().SetNoExponent()    
    hMCList[i]["state"].GetXaxis().SetTitle("")
    hMCList[i]["state"].GetYaxis().SetTitleOffset(0.4) 
    hMCList[i]["state"].GetXaxis().SetLabelOffset(0.5) 
    hMCList[i]["state"].SetTitle("") 
    hMCList_op[i]["state"].SetLineColor(PowColor);
    hMCList_op[i]["state"].SetLineWidth(PowLineWidth);
    hMCList_op[i]["state"].SetLineStyle(PowLineStyle);
     
    c1 = TCanvas( 'c1', hMCList[i]["name"] , 200, 10, 900, 1200 )
  
#    leg = ROOT.TLegend(.57,.45,.94,.91);
#    leg = ROOT.TLegend(.47,.35,.94,.91);

    leg = ROOT.TLegend(.5,.56,.94,.92);
    leg.SetBorderSize(0);
    #leg.SetTextSize(0.025);
    leg.SetTextSize(0.03);
    leg.AddEntry(hDataList[i]["state"],"Unfolded data + stat. uncertainty","lep")
    leg.AddEntry(grSyst,"Stat. #oplus syst. uncertainty","f") 

    leg.AddEntry(grMCSyst,mcref_entry,"flp") 
    leg.AddEntry(grMCSyst_op,mcop_entry,"flp")

    fInt = ROOT.TF1("constant","1",hMCList[i]["state"].GetXaxis().GetXmin(),    hMCList[i]["state"].GetXaxis().GetXmax());    
    line    =  ROOT.TLine(hMCList[i]["state"].GetXaxis().GetXmin(),1,hMCList[i]["state"].GetXaxis().GetXmax(),1) 
    line_op =  ROOT.TLine(hMCList[i]["state"].GetXaxis().GetXmin(),1,hMCList[i]["state"].GetXaxis().GetXmax(),1) 

    line.SetLineColor(MG5Color)
    line.SetLineStyle(MG5LineStyle)
    line.SetLineWidth(MG5LineWidth)

    line_op.SetLineColor(PowColor)
    line_op.SetLineStyle(PowLineStyle)
    line_op.SetLineWidth(PowLineWidth)

    c1.cd()
    pad1 = ROOT.TPad ('hist', '', 0., 0.45, 1.0, 1.0)
    if doLog: pad1.SetLogy()
#    pad1.SetRightMargin (0.04)
    pad1.SetRightMargin (0.045)
    pad1.SetLeftMargin (0.17)
    pad1.SetBottomMargin (0.015) 
    pad1.Draw()

    c1.cd()
    pad2 = ROOT.TPad ('rat', 'Data/MC ratio', 0., 0.28,  1., 0.45)
    pad2.SetTopMargin (0.0)
    pad2.SetRightMargin (0.045)
    pad2.SetLeftMargin (0.17)
    pad2.SetBottomMargin(0.35);
    pad2.Draw()

    c1.cd()
    pad3 = ROOT.TPad ('rat2', 'Data/MC ratio', 0., 0.17,  1., 0.34)
    pad3.SetTopMargin (0.0)
    pad3.SetRightMargin (0.045)
    pad3.SetLeftMargin (0.17)
    pad3.SetBottomMargin(0.35);
    pad3.Draw()

    #needed because of strange effects on pdf plots
    c1.cd()
    pad5 = ROOT.TPad ('space', 'space', 0.5, 0.1,  1., 0.17)
    pad5.SetTopMargin (0.0)
    pad5.SetRightMargin (0.045)
    pad5.SetLeftMargin (0.17)
    pad5.SetBottomMargin(0.60)
    pad5.Draw()
    pad1.cd() 

    if doLog:
        if DoFiducial:    Min = LogList[Type][DoNormalized][False]
        else:  Min = LogList_wide[Type][DoNormalized][False]
        hMCList[i]["state"].SetMinimum(Min)
    
    hMCList[i]["state"].GetXaxis().SetLabelSize(0)
    hMCList_op[i]["state"].GetXaxis().SetNoExponent()

    hMCList[i]["state"].Draw(MadDrawStyle)
    hMCList_op[i]["state"].Draw("same][")

    hMCList[i]["state"].GetXaxis().SetTitle("");

    grMCSyst.Draw("same2")   
    grMCSyst_op.Draw("same2")

    if "Mjj" or "Pt" in Type:hMCList[i]["state"].GetYaxis().SetTitleOffset(1.2)
    else: hMCList[i]["state"].GetYaxis().SetTitleOffset(1)
    hMCList[i]["state"].GetYaxis().SetTitleSize(0.06)
   
    if not UseMCReco:
        grSyst.Draw("same2")

    hDataList[i]["state"].Draw("sameE1")

    print "DataHei",hDataList[i]["name"]
    for bi in range(1, hDataList[i]["state"].GetNbinsX()+1): print hDataList[i]["state"].GetBinContent(bi),hDataList[i]["state"].GetBinError(bi)

    print "#### Statistic test ###\n\nOnly Statistic uncertainty"
    hDataList[i]["state"].Chi2Test(hMCList[i]["state"],"WW P")


    hDataForChi2  = GetHistoForChi2(hDataList[i]["state"],grSyst)
    hMCforChi2    = GetHistoForChi2(hMCList[i]["state"],grMCSyst)
    hMCforChi2_op = GetHistoForChi2(hMCList_op[i]["state"],grMCSyst_op)

    print "\nAll uncertainties\n"

    hDataForChi2.Chi2Test(hMCforChi2,"WW P")
    print "\nKolomogorov  ",hDataList[i]["state"].KolmogorovTest(hMCforChi2)
    PValue    =hDataForChi2.Chi2Test(hMCforChi2,"WW")
    PValue_op =hDataForChi2.Chi2Test(hMCforChi2_op,"WW")

    
    leg.Draw("same")
    
    pad2.cd()
    hRatio.GetYaxis().SetLabelSize(0.10);  
    hRatio.GetXaxis().SetNoExponent()    

    if "nJets" in Type:     hRatio.GetYaxis().SetNdivisions(604,1); 
    else:    hRatio.GetYaxis().SetNdivisions(604,1); 
    hRatio.GetXaxis().SetLabelSize(0); 
    hRatio.GetXaxis().SetTitle(""); 
    hRatio.GetXaxis().SetTitle(""); 
    hRatio.GetYaxis().SetTitleOffset(0.3);
    hRatio.GetYaxis().SetTitleSize(0.13);
    hRatio.GetXaxis().SetTitleSize(0.13); 
      
    hRatio.SetTitleSize(0.13);  

    if Type == "Mass": hRatio.SetMinimum(0.1); 
    elif Type == "Central_Mjj" or "Pt" in Type: hRatio.SetMinimum(-1.);
    elif "nJets_Central" in Type: 
        hRatio.SetMaximum(4.2);
        hRatio.SetMinimum(-0.5);
    elif "nJets" in Type:
        hRatio.SetMaximum(3.5);
        hRatio.SetMinimum(-0.7);
    else: hRatio.SetMinimum(-0.5);
    
    if do8TeVRatio:
        hRatio.SetMaximum(MaxRatioList_8TeV[Type][DoFiducial][True])
        hRatio.SetMinimum(MaxRatioList_8TeV[Type][DoFiducial][False])
    else:
        print "Min",MaxRatioList[Type][DoFiducial][False]
        hRatio.SetMaximum(MaxRatioList[Type][DoNormalized][True])
        hRatio.SetMinimum(MaxRatioList[Type][DoNormalized][False])

    hRatio.Draw("AXIS") 
    line.Draw("same")

    grMCSyst_ratio.Draw("same2")
   
    ratioTitle = ROOT.TText();
    grStat_ratio.Draw("P same")
    grSyst_ratio.Draw("same2")
    grSyst_ratio.Draw("same2")
    ratioTitle.SetNDC();
    #ratioTitle.SetTextFont(1);
    ratioTitle.SetTextColor(1);
    ratioTitle.SetTextSize(0.10);
    ratioTitle.SetTextAlign(22);
    ratioTitle.SetTextAngle(0);
#    if "Pow" in MCSetIn: ratioTitle.DrawText(0.4, 0.90,"{0}. p-value = {1:.3f}".format(mcref_entry,PValue)); 
 #   else: ratioTitle.DrawText(0.48, 0.90,"{0}. p-value = {1:.3f}".format(mcref_entry,PValue));  
#    ratioTitle.DrawText(0.48, 0.90,"{0}.".format(mcref_entry));  

    ratioTitle.DrawText(0.365, 0.90,"{0}.".format(mcref_entry));   #fixme

    pad3.cd()
    hRatio_op.GetYaxis().SetLabelSize(0.10);  
    hRatio_op.GetYaxis().SetNdivisions(604,1); 
    if "Mjj" in Type:  hRatio_op.GetXaxis().SetNdivisions(-505) #fixme add options in PlotStyle.py
    if "Jets" in Type: hRatio_op.GetXaxis().SetLabelSize(0.22)
    else: hRatio_op.GetXaxis().SetLabelSize(.18)
    hRatio_op.GetXaxis().SetLabelOffset(0.05);
    hRatio_op.GetXaxis().SetTitleOffset(1.5);
    hRatio_op.GetYaxis().SetTitleOffset(0.3);
    hRatio_op.GetYaxis().SetTitleSize(0.13);
    hRatio_op.GetXaxis().SetTitleSize(0.18); 
    hRatio_op.GetXaxis().SetTitle("");
    hRatio_op.GetXaxis().SetNoExponent()    

    if do8TeVRatio:
        hRatio_op.SetMaximum(MaxRatioList_8TeV[Type][DoNormalized][True])
        hRatio_op.SetMinimum(MaxRatioList_8TeV[Type][DoNormalized][False])
    else:
        hRatio_op.SetMaximum(MaxRatioList[Type][DoNormalized][True])
        hRatio_op.SetMinimum(MaxRatioList[Type][DoNormalized][False])

    hRatio_op.Draw("AXIS")
    line_op.Draw("same")
    grMCSyst_ratio_op.Draw("same2")
    grSyst_ratio_op.Draw("same2") 
    grStat_ratio_op.Draw("P same")

    ratioTitle_op = ROOT.TText();
    ratioTitle_op.SetNDC();
    ratioTitle_op.SetTextColor(1);
    ratioTitle_op.SetTextSize(0.1);
    ratioTitle_op.SetTextAlign(22);
    ratioTitle_op.SetTextAngle(0);

    ratioTitle.DrawText(0.32, 0.90,"{0}.".format(mcop_entry));  
    #if "Pow" in MCSetIn: ratioTitle_op.DrawText(0.48, 0.90,"{0}. p-value = {1:.3f}".format(mcop_entry,PValue_op));
    #else: ratioTitle_op.DrawText(0.4, 0.90,"{0}. p-value = {1:.3f}".format(mcop_entry,PValue_op)); 

    grSyst_ratio_op.Draw("same2") 
 
    pad5.cd()
    ratioTitle = ROOT.TLatex(0.94,0.48,xName);
    ratioTitle.SetNDC();
    ratioTitle.SetTextColor(1);
    ratioTitle.SetTextSize(0.51); 
    ratioTitle.SetTextFont(42);
    ratioTitle.SetTextAlign(32);
    ratioTitle.SetTextAngle(0);
    ratioTitle.Draw();

    #write CMS Preliminary and luminosity
    CMS_lumi.CMS_lumi(pad1, iPeriod, iPos)

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

    if do8TeVRatio:
        Kind+="_8TeV"

    if SavePlot:
        if doLog:
            c1.SaveAs(PersonalFolder+Dir+"/"+Type+"/CrossSection/diffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+Kind+"_log"+".png")        
            c1.SaveAs(PersonalFolder+Dir+"/"+Type+"/CrossSection/diffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+Kind+"_log"+".eps")        
            c1.SaveAs(PersonalFolder+Dir+"/"+Type+"/CrossSection/diffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+Kind+"_log"+".pdf")        
            c1.SaveAs(PersonalFolder+Dir+"/"+Type+"/CrossSection/diffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+Kind+"_log"+".root")       
            c1.SaveAs(PersonalFolder+Dir+"/"+Type+"/CrossSection/diffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+Kind+"_log"+".C")       
        else:
            c1.SaveAs("Plot/CrossSection/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+Kind+".png")
            c1.SaveAs("Plot/CrossSection/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+Kind+".root")
            c1.SaveAs(PersonalFolder+Dir+"/"+Type+"/CrossSection/diffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+Kind+".png")      
            c1.SaveAs(PersonalFolder+Dir+"/"+Type+"/CrossSection/diffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+Kind+".eps")      
            c1.SaveAs(PersonalFolder+Dir+"/"+Type+"/CrossSection/diffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+Kind+".pdf")      
            c1.SaveAs(PersonalFolder+Dir+"/"+Type+"/CrossSection/diffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+Kind+".root")       

    print "Pvalue ",PValue
     
if printLatex and (Type=="nJets" or Type=="nIncJets"): 
    print "\nJet multuplicity cross section tabular:\n"
    Text =("Number of jets ($|\eta_{j}| < 4.7$)","Total cross-section [fb]","stat","syst","lumi")
    CrossLatex(Text,CrossDic,Type)
