#! /usr/bin/env python

##################################
## G. Pinna (UNITO) - Jun 2015 ###
##################################

from optparse import OptionParser
import ROOT,copy
ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import gSystem, TCanvas, TH1,  TPad, gStyle, TLegend, TGraphAsymmErrors, Math#,QuantFuncMathCore 
from collections import OrderedDict
from CrossSection_copy import*
from LL import*
import sys,ast,os
import math
import operator
import CMS_lumi, tdrstyle

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

parser.add_option("-v", "--GenMCVar", dest="MCvar",
                  default="central",
                  help="Up, central or down MC generated distributions")

parser.add_option("-c", "--DataOverMCratio", dest="DoDataOverMC",
                  default="True",
                  help="Compute normilized cross section. Default is True")

parser.add_option("-p", "--InclusiveXSFromIntegral", dest="DoInclInt",
                  default="False",
                  help="Compute the inclusive cross section as distribution integral. Defaiult is False")

parser.add_option("-q", "--DifferentialXSFromIntegral", dest="DoDiffInt",
                  default="False",
                  help="Compute the differentiaal (bin per bin) cross section. Defaiult is False")

(options, args) = parser.parse_args()


SavePlot  =ast.literal_eval(options.SavePlot)
UseUnfold =ast.literal_eval(options.UseUnfold)
UseMCReco =ast.literal_eval(options.UseMCReco)
DoInclusive =ast.literal_eval(options.DoInclusive)
MCSetIn = options.MCSetIn
Type = options.Type
DoNormalized =ast.literal_eval(options.DoNormalized)
MCvar = options.MCvar
DoDataOverMC =ast.literal_eval(options.DoDataOverMC)
DoInclInt = ast.literal_eval(options.DoInclInt)
DoDiffInt = ast.literal_eval(options.DoDiffInt)

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

#set the tdr style
tdrstyle.setTDRStyle()

#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
CMS_lumi.lumi_8TeV = "19.7 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "19.7 fb^{-1} (8 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 0

ROOT.gStyle.SetHatchesSpacing(0.3)
ROOT.gStyle.SetHatchesLineWidth(1)
if DoInclusive:
    TotalCross(MCSetIn,Type)
    sys.exit(0)


if MCSetIn == "Pow": mcop = "Mad"
elif MCSetIn == "Mad" : mcop = "Pow"
elif MCSetIn == "fr_Pow" : mcop = "fr_Mad"
elif MCSetIn == "fr_Mad" : mcop = "fr_Pow"
(hMCList,hMCListUp,hMCListDown) = getCrossPlot_MC(MCSetIn,Type,DoNormalized)
(hMCList_op,hMCListUp_op,hMCListDown_op) = getCrossPlot_MC(mcop,Type,DoNormalized)

if DoInclInt == True or DoDiffInt == True:  optInt = True
else: optInt = False

if "fr" in MCSetIn: 
    (hMCList_MGatNLO,hMCListUp_MGatNLO,hMCListDown_MGatNLO) = getCrossPlot_MC("fr_MGatNLO",Type,DoNormalized)
    fb =1000
    um = "fb"
else: (hMCList_MGatNLO,hMCListUp_MGatNLO,hMCListDown_MGatNLO) = getCrossPlot_MC("MGatNLO",Type,DoNormalized)
 
if UseMCReco:
    (hDataList,hDataListUp,hDataListDown) = getCrossPlot_Data(MCSetIn,UseUnfold,Type,0,True,DoNormalized,optInt) 
    fb = 1
    um = "pb"
else:
    (hDataList,hDataListUp,hDataListDown) = getCrossPlot_Data(MCSetIn,UseUnfold,Type,0,False,DoNormalized,optInt) 
    

#Inclusive cross-section obtained from the integral of the data distribution (remember to run SetHistos.py with DoInclInt == True, not to include JER and JES systematics)
if DoInclInt:
 
    for i in range(0,4): 
        grSist = getUncGraph(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],MCSetIn,hDataList[i]["name"],True,"syst",False,optInt)  
        print Red("Final Results for the")+(" "*(9-len("Final Results for the"))), 
        print Red(hDataList[i]["name"]+ " Final State")
        err_syst_up = 0
        err_syst_down = 0
        err_syst_up_sq = 0
        err_syst_down_sq = 0 
        staterr = 0
        lumierr =  0
        binning = hDataList[i]["state"].GetNbinsX()

        for b in range(1,binning+1):
            err_syst_up_sq += (hDataListUp[i]["state"].GetBinContent(b)- hDataList[i]["state"].GetBinContent(b))*(hDataListUp[i]["state"].GetBinContent(b)- hDataList[i]["state"].GetBinContent(b))  
            err_syst_down_sq += (hDataList[i]["state"].GetBinContent(b)- hDataListDown[i]["state"].GetBinContent(b))*(hDataList[i]["state"].GetBinContent(b)- hDataListDown[i]["state"].GetBinContent(b))

        err_syst_up = math.sqrt(err_syst_up_sq)
        err_syst_down = math.sqrt(err_syst_down_sq)
        staterr = ROOT.Double(0)
        xs_incl =  hDataList[i]["state"].IntegralAndError(1,-1,staterr,"width")*fb
        lumierr = xs_incl*lumiErr_perc

    
        print  " xs inclusive = {0:.2f} +-{1:.2f} (stat) +{2:.2f} -{3:.2f} (syst) +-{4:.2f} (lumi.) [{5}]\n".format(xs_incl,staterr*fb, err_syst_up*fb, err_syst_down*fb,lumierr,um) 

    sys.exit(0)



#Cross-section obtained as a function of the number of the bin (remember to run SetHistos.py with DoInclInt == False, to include JER and JES systematics)
if DoDiffInt:
    for i in range(0,4): 
        grSist = getUncGraph(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],MCSetIn,hDataList[i]["name"],True,"syst",False,optInt)  
        print Red("Final Results for the")+(" "*(9-len("Final Results for the"))), 
        print Red(hDataList[i]["name"]+ " Final State")
        
        binning = hDataList[i]["state"].GetNbinsX()
        for b in range(1,binning+1):
       
            if "fr" in MCSetIn:
                lumierr = 0
                lumierr =  hDataList[i]["state"].GetBinContent(b)*fb*lumiErr_perc
                print Blue("bin "+ str(b)), 
                print " {0:.3f} +-{1:.3f} (+{2:.3f} -{3:.3f}) (stat) +-{4:.3f} (+{5:.3f} -{6:.3f}) (syst) +- {7:.3f}% (+{8:.3f}% -{9:.3f}%)(syst %) +-{10:.3f} (lumi.)[fb]\n".format( hDataList[i]["state"].GetBinContent(b)*fb, hDataList[i]["state"].GetBinError(b)*fb, hDataList[i]["state"].GetBinErrorUp(b)*fb, hDataList[i]["state"].GetBinErrorLow(b)*fb,  grSist.GetErrorY(b-1)*fb, (hDataListUp[i]["state"].GetBinContent(b)- hDataList[i]["state"].GetBinContent(b))*fb,(hDataList[i]["state"].GetBinContent(b)- hDataListDown[i]["state"].GetBinContent(b))*fb, grSist.GetErrorY(b-1)*100/hDataList[i]["state"].GetBinContent(b), (hDataListUp[i]["state"].GetBinContent(b)- hDataList[i]["state"].GetBinContent(b))*100/hDataList[i]["state"].GetBinContent(b),(hDataList[i]["state"].GetBinContent(b)- hDataListDown[i]["state"].GetBinContent(b))*100/hDataList[i]["state"].GetBinContent(b),lumierr)

            else: 
                print Blue("bin "+ str(b)), 
                print " xs =", hDataList[i]["state"].GetBinContent(b), " +- ",  hDataList[i]["state"].GetBinError(b), " +- ",  grSist.GetErrorY(b-1), "pb"  
    sys.exit(0) 
       

#To obtain the final plots:
 #loop on final states
for i in range(0,4):
    
    hDataList[i]["state"].SetMarkerSize(1.3)
    hDataList[i]["state"].SetMarkerColor(1)
    hDataList[i]["state"].SetMarkerStyle(20)
    hDataList[i]["state"].SetLineColor(1)
    
    grSist = getUncGraph(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],MCSetIn,hDataList[i]["name"],True,"statsyst",False,optInt)  
    if "Mad" in MCSetIn: grSistVar = getUncGraph(hMCList_op[i]["state"],hMCListUp_op[i]["state"],hMCListDown_op[i]["state"],mcop,hMCList_op[i]["name"],False,"syst",False,optInt)  
    else: grSistVar = getUncGraph(hMCList[i]["state"],hMCListUp[i]["state"],hMCListDown[i]["state"],mcop,hMCList[i]["name"],False,"syst",False,optInt)  

    grSistVar_MGatNLO = getUncGraph(hMCList_MGatNLO[i]["state"],hMCListUp_MGatNLO[i]["state"],hMCListDown_MGatNLO[i]["state"],"MGatNLO",hMCList_MGatNLO[i]["name"],False,"syst",False,optInt) 


    grSist.SetFillStyle(3005)
    grSist.SetFillColor(ROOT.kBlack)
    grSistVar.SetFillStyle(3001)
    grSistVar.SetFillColorAlpha(ROOT.kRed-4,0.70)
    grSistVar.SetFillStyle(3350) 
    #grSistVar.SetFillStyle(3001) 
    grSistVar.SetLineStyle(1) 
    grSistVar.SetLineWidth(1)
    grSistVar.SetLineColorAlpha(ROOT.kRed-9,0.10) 
    grSistVar.SetLineColor(ROOT.kRed) 
    grSistVar.SetLineWidth(2)  
    grSistVar.SetLineStyle(1)  
    grSistVar_MGatNLO.SetFillStyle(3001)
    grSistVar_MGatNLO.SetFillColorAlpha(ROOT.kBlue,0.70) 
    #grSistVar_MGatNLO.SetFillColor(ROOT.kBlue) 
    # grSistVar_MGatNLO.SetLineStyle(1) 
    # grSistVar_MGatNLO.SetLineWidth(1)
    grSistVar_MGatNLO.SetLineColorAlpha(ROOT.kBlue,0.10)   
    grSistVar_MGatNLO.SetLineColor(ROOT.kBlue) 
    grSistVar_MGatNLO.SetLineWidth(2)  
    grSistVar_MGatNLO.SetLineStyle(1) 
    grSistVar_MGatNLO.SetFillStyle(3359)  
    #grSistVar_MGatNLO.SetFillStyle(3001)  
    #grSistVar_MGatNLO.SetLineWidth(2)
    Err=ROOT.Double(0.)
    
    LastBin = hMCList[i]["state"].GetXaxis().GetLast() 
   
    if "Pow" in MCSetIn: 
        mcref_entry = "Powheg+MCFM+Phantom"
        mcop_entry = "MadGraph+MCFM+Phantom"
    else:
        mcref_entry = "MadGraph+MCFM+Phantom"
        mcop_entry = "Powheg+MCFM+Phantom"
    
###########################################################
   
    #Ratio and uncertainties on data points for the reference set of samples 
          
    #statistical uncertainty
    (hRatio,hRatio_up,hRatio_down) = LL(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],hMCList[i]["state"],Type,hDataList[i]["name"],"stat",DoDataOverMC)     
    grStat_ratio = getUncGraph(hRatio,hRatio_up,hRatio_down,MCSetIn,hDataList[i]["name"],True,"syst",False,optInt)  
    
    #systematic uncertainty
    (hRatio_ss,hRatio_up_ss,hRatio_down_ss) = LL(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],hMCList[i]["state"],Type,hDataList[i]["name"],"statsyst",DoDataOverMC) 
    grSist_ratio = getUncGraph(hRatio_ss,hRatio_up_ss,hRatio_down_ss,MCSetIn,hDataList[i]["name"],True,"syst",False,optInt)  
    grSist_ratio.SetFillStyle(3005)
    grSist_ratio.SetFillColor(ROOT.kBlack)

    maxtot = hRatio_up_ss.GetMaximum()   


###################################################################

    #Ratio and uncertainties on data for the opposite set of samples 
 
    #statistical uncertainty
    (hRatio_op,hRatio_op_up,hRatio_op_down) = LL(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],hMCList_op[i]["state"],Type,hDataList[i]["name"],"stat",DoDataOverMC)     
    grStat_ratio_op = getUncGraph(hRatio_op,hRatio_op_up,hRatio_op_down,MCSetIn,hDataList[i]["name"],True,"syst",False,optInt)   
  
    #systematic uncertainty
    (hRatio_op_ss,hRatio_op_up_ss,hRatio_op_down_ss) = LL(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],hMCList_op[i]["state"],Type,hDataList[i]["name"],"statsyst",DoDataOverMC)  
    grSist_ratio_op = getUncGraph(hRatio_op_ss,hRatio_op_up_ss,hRatio_op_down_ss,MCSetIn,hDataList[i]["name"],True,"syst",False,optInt)  
    grSist_ratio_op.SetFillStyle(3005)
    grSist_ratio_op.SetFillColor(ROOT.kBlack)

    maxtot_op =  hRatio_op_up_ss.GetMaximum() 
    maxratios_op = max(maxtot,maxtot_op)


    #histograms for the theoretical uncertainty
          
    if "Mad" in MCSetIn: grSistMC_ratio_op = getUncGraphMC(hMCList_op[i]["state"],hMCListUp_op[i]["state"],hMCListDown_op[i]["state"],MCSetIn,hDataList[i]["name"])  
    else: grSistMC_ratio_op = getUncGraphMC(hMCList[i]["state"],hMCListUp[i]["state"],hMCListDown[i]["state"],MCSetIn,hDataList[i]["name"])

    grSistMC_ratio_op.SetFillStyle(3001)
    grSistMC_ratio_op.SetFillColorAlpha(ROOT.kRed -4,0.70)
    grSistMC_ratio_op.SetFillStyle(3350)
    #grSistMC_ratio_op.SetFillStyle(3001)

 ##############################################

    #Ratio and uncertainties on data for the MGatNLO set of samples 
  
    #statistical uncertainty
    (hRatio_MGatNLO,hRatio_MGatNLO_up,hRatio_MGatNLO_down) = LL(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],hMCList_MGatNLO[i]["state"],Type,hDataList[i]["name"],"stat",DoDataOverMC)     
 
    grStat_ratio_MGatNLO = getUncGraph(hRatio_MGatNLO,hRatio_MGatNLO_up,hRatio_MGatNLO_down,MCSetIn,hDataList[i]["name"],True,"syst",False,optInt)
 

    #systematic uncertainty
    (hRatio_MGatNLO_ss,hRatio_MGatNLO_up_ss,hRatio_MGatNLO_down_ss) = LL(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],hMCList_MGatNLO[i]["state"],Type,hDataList[i]["name"],"statsyst",DoDataOverMC)  
    grSist_ratio_MGatNLO = getUncGraph(hRatio_MGatNLO_ss,hRatio_MGatNLO_up_ss,hRatio_MGatNLO_down_ss,MCSetIn,hDataList[i]["name"],True,"syst",False,optInt)  

    # for l in range(1,binning+1):
    #     print "MGATNLO"
    #     print "data",hDataList[i]["state"].GetBinContent(l), "data error",hDataList[i]["state"].GetBinError(l)/hDataList[i]["state"].GetBinContent(l),"data up",hDataListUp[i]["state"].GetBinContent(l),"data down",hDataListDown[i]["state"].GetBinContent(l)
    #     print "mc",hMCList_MGatNLO[i]["state"].GetBinContent(l)
    #     print "ratio", hRatio_MGatNLO.GetBinContent(l),"ratio up",hRatio_MGatNLO_up.GetBinContent(l),"ratio down", hRatio_MGatNLO_down.GetBinContent(l)
    
    grSist_ratio_MGatNLO.SetFillStyle(3005)
    grSist_ratio_MGatNLO.SetFillColor(ROOT.kBlack)

    maxtot_MGatNLO =  hRatio_MGatNLO_up.GetMaximum()
  
    if "EtaJet2" in Type: maxratios = max(maxratios_op,maxtot_MGatNLO)+0.35*max(maxratios_op,maxtot_MGatNLO)
    elif "PtJet1" in Type: maxratios = max(maxratios_op,maxtot_MGatNLO)+0.30*max(maxratios_op,maxtot_MGatNLO)
    elif "PtJet2" in Type: maxratios = max(maxratios_op,maxtot_MGatNLO)+0.35*max(maxratios_op,maxtot_MGatNLO)
    elif "Mjj" in Type: maxratios = max(maxratios_op,maxtot_MGatNLO)+0.35*max(maxratios_op,maxtot_MGatNLO) 
    elif "Deta" in Type: maxratios = max(maxratios_op,maxtot_MGatNLO)+0.50*max(maxratios_op,maxtot_MGatNLO)
    else: maxratios = max(maxratios_op,maxtot_MGatNLO)+0.80*max(maxratios_op,maxtot_MGatNLO)
 
    #histograms for the theoretical uncertainty
    grSistMC_ratio_MGatNLO = getUncGraphMC(hMCList_MGatNLO[i]["state"], hMCListUp_MGatNLO[i]["state"],hMCListDown_MGatNLO[i]["state"],MCSetIn,hDataList[i]["name"])  

    grSistMC_ratio_MGatNLO.SetFillStyle(3001)
    grSistMC_ratio_MGatNLO.SetFillColorAlpha(ROOT.kBlue,0.70)
    grSistMC_ratio_MGatNLO.SetFillStyle(3359)
    #grSistMC_ratio_MGatNLO.SetFillStyle(3001)

#############################  

    if "Central" in Type: eta_cut = " (|#eta^{jet}| < 2.4)"
    else: eta_cut = " (|#eta^{jet}| < 4.7)"

    if Type=="Mass":
        xName = "m_{4l} [GeV]"
        hRatio_MGatNLO.GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{dm_{4l}} [1/GeV]")
        else: hMCList[i]["state"].GetYaxis().SetTitle("#frac{d#sigma}{dm_{4l}} [pb/GeV]")

    elif "Jets" in Type:
        hRatio_MGatNLO.GetXaxis().SetBinLabel(1,"0")
        hRatio_MGatNLO.GetXaxis().SetBinLabel(2,"1")
        hRatio_MGatNLO.GetXaxis().SetBinLabel(3,"2")
        hRatio_MGatNLO.GetXaxis().SetBinLabel(4,"#geq 3")

        hMCList[i]["state"].GetXaxis().SetBinLabel(1,"0")
        hMCList[i]["state"].GetXaxis().SetBinLabel(2,"1")
        hMCList[i]["state"].GetXaxis().SetBinLabel(3,"2")
        hMCList[i]["state"].GetXaxis().SetBinLabel(4,"#geq 3")
        
        xName = "N_{jets}"+eta_cut
       
        hRatio_MGatNLO.GetXaxis().SetTitle(xName)
        hMCList[i]["state"].GetXaxis().SetTitle(xName)

        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{dN_{jets}}")        
        else: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/dN_{jets} [pb]") 

    elif "Deta" in Type:
        xName = "#Delta#eta_{jj}" + eta_cut
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{d#Delta#eta_{jj}}")
        else: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d#Delta#eta_{jj}[pb]")  

    elif "Mjj" in Type:
        xName = "m_{jj}" + eta_cut+ " [GeV]"
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{dm_{jj}} [1/GeV]")      
        else: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/dm_{jj} [pb/GeV]")  

    elif "PtJet1" in Type:
        xName = "p_{T}^{jet1} [GeV]"
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{dp_{T}^{jet1}} [1/GeV]")             
        else:  hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/dp_{T}^{jet1} [pb/GeV]") 
    elif "PtJet2" in Type:
        xName = "p_{T}^{jet2} [GeV]"
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{dp_{T}^{jet2}} [1/GeV]")             
        else:  hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/dp_{T}^{jet2} [pb/GeV]") 
    elif "EtaJet1" in Type:
        xName = "|#eta^{jet1}|"
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{d|#eta^{jet1}|}")             
        else: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d|#eta^{jet1}| [pb]") 
    elif "EtaJet2" in Type:
        xName = "|#eta^{jet2}|"
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{d|#eta^{jet2}|}")             
        else: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d|#eta^{jet2}| [pb]")

    hRatio_MGatNLO.GetXaxis().SetTitle(xName)
    #hRatio_op.GetXaxis().SetTitle(xName)
    if DoDataOverMC == True: ratioName = "Data/MC"
    else: ratioName = "MC/Data"
    hRatio.GetYaxis().SetTitle(ratioName) 
    hRatio_op.GetYaxis().SetTitle(ratioName)
    hRatio_MGatNLO.GetYaxis().SetTitle(ratioName)
  

    Max= hDataList[i]["state"].GetMaximum()+hDataList[i]["state"].GetBinError(hDataList[i]["state"].GetMaximumBin())
    MaxMC=hMCList[i]["state"].GetMaximum()


    if Type == "Mass": 
        if DoNormalized:
            if "fr" in MCSetIn: Max=0.013 #Max=Max+Max/5.
            else: Max = 0.013
        else: 
            if "fr" in MCSetIn: Max=0.0003 #Max=Max+Max/5.
            else: Max = 0.1
    elif "Mjj" in Type: 
        if DoNormalized:
            Max=0.0060
        else: Max = 0.006  
    elif "Deta" in Type: 
        if DoNormalized:
            Max=1.
        else: Max=0.4
              
    elif Type=="PtJet1":
        if DoNormalized: Max=0.05
        else: Max=0.1
        
    elif Type=="PtJet2":
        if DoNormalized: Max=0.023
        else: Max=0.015

    elif Type=="EtaJet1":
       if DoNormalized: Max=1.
       else:  Max=1.

    elif Type=="EtaJet2":
        if DoNormalized:Max = 1.1
        else: Max=0.70
   
    elif "Jets" in Type:
        if DoNormalized: Max = 1.2
        else:
            if "fr" in MCSetIn: Max = 0.023
            else: Max=10.
   
    else: Max=Max+Max/9.

    #if previous ZZ cross-sections
    #Max_data= hDataList[i]["state"].GetMaximum()+hDataList[i]["state"].GetBinError(hDataList[i]["state"].GetMaximumBin())
    #MaxMC=hMCList[i]["state"].GetMaximum()
    #Max = max(Max_data,MaxMC)*2

    hMCList[i]["state"].SetMaximum(Max)
    hMCList[i]["state"].SetMinimum(0)
    #hMCList[i]["state"].SetFillColor(ROOT.kCyan-3)
    hMCList[i]["state"].SetFillColor(851)
    hMCList[i]["state"].SetLineColor(852)
    hMCList[i]["state"].SetLineStyle(1)
    hMCList[i]["state"].SetLineWidth(2)
    #hMCList[i]["state"].GetXaxis().SetTitle(xName)
    hMCList[i]["state"].GetXaxis().SetTitle("")
    hMCList[i]["state"].GetYaxis().SetTitleOffset(1.4) 
    hMCList[i]["state"].GetXaxis().SetLabelOffset(0.5)
    #hMCList[i]["state"].SetTitle(hMCList[i]["name"]+" Cross Section")
    hMCList[i]["state"].SetTitle("")
    hMCList_op[i]["state"].SetLineColor(ROOT.kRed); 
   # hMCList_op[i]["state"].SetFillStyle(307);
    hMCList_op[i]["state"].SetLineWidth(2);
    hMCList_op[i]["state"].SetLineStyle(1);
    hMCList_MGatNLO[i]["state"].SetLineColor(ROOT.kBlue);
    hMCList_MGatNLO[i]["state"].SetLineWidth(2);
    hMCList_MGatNLO[i]["state"].SetLineStyle(1);
    #hMCList_MGatNLO[i]["state"].SetFillStyle(308);

    #c1 = TCanvas( 'c1', hMCList[i]["name"] , 200, 10, 900, 700 )
    c1 = TCanvas( 'c1', hMCList[i]["name"] , 200, 10, 900, 1200 )
   
    leg = ROOT.TLegend(.51,.50,.800,.85);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.032);
    leg.AddEntry(hDataList[i]["state"],"Unfolded data + stat. uncertainty","lep")
    leg.AddEntry(grSist,"Total uncertainty","f") 
    leg.AddEntry(hMCList[i]["state"],mcref_entry,"f") 
    leg.AddEntry(grSistVar,mcop_entry,"fl")
    leg.AddEntry(grSistVar_MGatNLO, "MadGraph5_aMCatNLO+MCFM+Phantom","fl")
   
    if(Type == "Mass"): 
        line = ROOT.TLine(100,1,800,1)
        line_op = ROOT.TLine(100,1,800,1)
    elif(Type == "Jets" or Type == "CentralJets"): 
        line =  ROOT.TLine(0,1,4,1)
        line_op =  ROOT.TLine(0,1,4,1)
    elif(Type == "Mjj" or Type == "CentralMjj"): 
        line =  ROOT.TLine(0,1,800,1)
        line_op =  ROOT.TLine(0,1,800,1)
    elif(Type == "Deta"  or Type == "CentralDeta"): 
        line =  ROOT.TLine(0,1,4.7,1)
        line_op =  ROOT.TLine(0,1,4.7,1) 
    elif(Type =="PtJet1"):  
        line =  ROOT.TLine(30,1,300,1) 
        line_op =  ROOT.TLine(30,1,300,1) 
    elif(Type =="PtJet2" ):  
        line =  ROOT.TLine(30,1,200,1) 
        line_op =  ROOT.TLine(30,1,200,1) 
    elif(Type =="EtaJet1"or Type =="EtaJet2" ):  
        line =  ROOT.TLine(0,1,4.7,1)
        line_op =  ROOT.TLine(0,1,4.7,1) 
    #line.SetLineColor(852);
    line.SetLineColor(ROOT.kBlack);
    line.SetLineWidth(1);
    #line_op.SetLineColor(ROOT.kRed);
    line_op.SetLineColor(ROOT.kBlack);
    line_op.SetLineWidth(1);

    c1.cd()
    pad1 = ROOT.TPad ('hist', '', 0., 0.51, 1.0, 1.0)#0.35
    pad1.SetTopMargin (0.10)
    pad1.SetRightMargin (0.03)#0.10
    pad1.SetLeftMargin (0.17)
    pad1.SetBottomMargin (0.03) 
    pad1.Draw()
   
    c1.cd()
    pad2 = ROOT.TPad ('rat', 'Data/MC ratio', 0., 0.31,  1., 0.51)#0.15
    pad2.SetTopMargin (0.0)
    pad2.SetRightMargin (0.03)#0.10
    pad2.SetLeftMargin (0.17)
    pad2.SetBottomMargin(0.35);
    pad2.Draw()
     
    c1.cd()
    pad3 = ROOT.TPad ('rat2', 'Data/MC ratio', 0., 0.18,  1., 0.38)#0.02-0.22
    pad3.SetTopMargin (0.0)
    pad3.SetRightMargin (0.03)#0.10
    pad3.SetLeftMargin (0.17)
    pad3.SetBottomMargin(0.35);
    pad3.Draw()
 
    c1.cd()
    pad4 = ROOT.TPad ('rat3', 'Data/MC ratio', 0., 0.05,  1., 0.25)#0.02-0.22
    pad4.SetTopMargin (0.0)
    pad4.SetRightMargin (0.03)#0.10
    pad4.SetLeftMargin (0.17)
    pad4.SetBottomMargin(0.35)#0.35
    pad4.Draw()
    
    #needed because of strange effects on pdf plots
    c1.cd()
    pad5 = ROOT.TPad ('space', 'space', 0.5, 0.0,  1., 0.07)#0.02-0.22
    pad5.SetTopMargin (0.0)
    pad5.SetRightMargin (0.03)#0.10
    pad5.SetLeftMargin (0.17)
    pad5.SetBottomMargin(0.60)#0.35
    pad5.Draw()
   
    pad1.cd() 
    hMCList[i]["state"].Draw("hist")
    hMCList_op[i]["state"].Draw("hist, same")
    hMCList_MGatNLO[i]["state"].Draw("hist, same")
    if "Mjj" or "Pt" in Type:hMCList[i]["state"].GetYaxis().SetTitleOffset(1.2)
    else: hMCList[i]["state"].GetYaxis().SetTitleOffset(1)
    hMCList[i]["state"].GetYaxis().SetTitleSize(0.06)
   
    if not UseMCReco:
        grSistVar_MGatNLO.Draw("same2")
        grSistVar.Draw("same2")
        grSist.Draw("same2")
       
    hDataList[i]["state"].Draw("same E1")
    leg.Draw("same")
      
    pad2.cd()
    hRatio.GetYaxis().SetLabelSize(0.10);  
    #hRatio.GetXaxis().CenterLabels();
    hRatio.GetYaxis().SetNdivisions(604,1); 
    hRatio.GetXaxis().SetLabelSize(0.13); 
    hRatio.GetXaxis().SetLabelOffset(0.5);
    hRatio.GetXaxis().SetTitle(""); 
    hRatio.GetYaxis().SetTitleOffset(0.4);
    hRatio.GetYaxis().SetTitleSize(0.13);#0.13
    hRatio.GetXaxis().SetTitleSize(0.13); #0.13  

    hRatio.SetTitleSize(0.13);  
    if maxratios > 1.9: hRatio.SetMaximum(maxratios); 
    else: hRatio.SetMaximum(1.9); 
    if Type == "Mass": hRatio.SetMinimum(0.1); 
    elif Type == "CentralMjj" or "Pt" in Type: hRatio.SetMinimum(-1.);
    else: hRatio.SetMinimum(-0.5);


    hRatio.Draw() 
    line.Draw("same")
    grSist_ratio.Draw("same2")
    hRatio.Draw("same")
    grSist_ratio.Draw("same2")
    ratioTitle = ROOT.TText();
    ratioTitle.SetNDC();
    #ratioTitle.SetTextFont(1);
    ratioTitle.SetTextColor(1);
    ratioTitle.SetTextSize(0.10);
    ratioTitle.SetTextAlign(22);
    ratioTitle.SetTextAngle(0);
    if "Pow" in MCSetIn: ratioTitle.DrawText(0.345, 0.90,mcref_entry);
    else: ratioTitle.DrawText(0.36, 0.90,mcref_entry);  
    hRatio.Draw("same")
    grStat_ratio.Draw("P same")
    
    pad3.cd()
    hRatio_op.GetYaxis().SetLabelSize(0.10);  
    hRatio_op.GetYaxis().SetNdivisions(604,1); 
    if "Jets" in Type: hRatio_op.GetXaxis().SetLabelSize(0.22)#0.15; 
    else: hRatio_op.GetXaxis().SetLabelSize(0.15)#0.15; 
    hRatio_op.GetXaxis().SetLabelOffset(0.05);
    hRatio_op.GetXaxis().SetTitleOffset(1.5);#1.2; 
    hRatio_op.GetYaxis().SetTitleOffset(0.4);
    hRatio_op.GetYaxis().SetTitleSize(0.13);#0.13
    hRatio_op.GetXaxis().SetTitleSize(0.18); 
    
    if maxratios > 1.9: hRatio_op.SetMaximum(maxratios); 
    else: hRatio_op.SetMaximum(1.9); 
    if Type == "Mass": hRatio_op.SetMinimum(0.1); 
    elif Type == "CentralMjj" or "Pt" in Type: hRatio_op.SetMinimum(-1.);
    else: hRatio_op.SetMinimum(-0.5);
    
    hRatio_op.Draw()
    line_op.Draw("same")
    grSist_ratio_op.Draw("same2") 
    grStat_ratio_op.Draw("P same")
    hRatio_op.Draw("same") 
    if "Mad" in MCSetIn: grSistMC_ratio_op.Draw("same2")  #uncomment
    grSist_ratio_op.Draw("same2")#uncomment
   
    ratioTitle_op = ROOT.TText();
    ratioTitle_op.SetNDC();
    ratioTitle_op.SetTextColor(1);
    ratioTitle_op.SetTextSize(0.10);
    ratioTitle_op.SetTextAlign(22);
    ratioTitle_op.SetTextAngle(0);
    if "Pow" in MCSetIn: ratioTitle_op.DrawText(0.36, 0.90,mcop_entry);
    else: ratioTitle_op.DrawText(0.345, 0.90,mcop_entry); 
    hRatio_op.Draw("same")
    grSist_ratio_op.Draw("same2") 
    grStat_ratio_op.Draw("P same")
    ROOT.gStyle.SetOptStat(0);   
  
   
    pad4.cd()
    hRatio_MGatNLO.GetYaxis().SetLabelSize(0.10);  
    hRatio_MGatNLO.GetYaxis().SetNdivisions(604,1); 
    if "Jets" in Type: hRatio_MGatNLO.GetXaxis().SetLabelSize(0.22)#0.15; 
    else: hRatio_MGatNLO.GetXaxis().SetLabelSize(0.15)#0.15; 
    hRatio_MGatNLO.GetXaxis().SetLabelOffset(0.05);
    hRatio_MGatNLO.GetXaxis().SetTitleOffset(1.5);#1.2; 
    hRatio_MGatNLO.GetYaxis().SetTitleOffset(0.4);
    hRatio_MGatNLO.GetYaxis().SetTitleSize(0.13);#0.13
    hRatio_MGatNLO.GetXaxis().SetTitleSize(0.18); 
    if maxratios > 1.9: hRatio_MGatNLO.SetMaximum(maxratios); 
    else: hRatio_MGatNLO.SetMaximum(1.9); 
    if Type == "Mass": hRatio_MGatNLO.SetMinimum(0.1); 
    elif Type == "CentralMjj" or "Pt" in Type: hRatio_MGatNLO.SetMinimum(-1.);
    else: hRatio_MGatNLO.SetMinimum(-0.5);
    
    hRatio_MGatNLO.Draw()
    line_op.Draw("same")
    grStat_ratio_MGatNLO.Draw("P same") 
    grSist_ratio_MGatNLO.Draw("same2")
    hRatio_MGatNLO.Draw("same")  
    grSistMC_ratio_MGatNLO.Draw("same2")
    hRatio_MGatNLO.Draw("same")
    grStat_ratio_MGatNLO.Draw("P same") 

    ratioTitle_MGatNLO = ROOT.TText();
    ratioTitle_MGatNLO.SetNDC();
    ratioTitle_MGatNLO.SetTextColor(1);
    ratioTitle_MGatNLO.SetTextSize(0.10);
    ratioTitle_MGatNLO.SetTextAlign(22);
    ratioTitle_MGatNLO.SetTextAngle(0);
    ratioTitle_MGatNLO.DrawText(0.435, 0.90,"MadGraph5_aMCatNLO+MCFM+Phantom");
    
    ROOT.gStyle.SetOptStat(0); 

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
    if MCSetIn == "Pow":     MCSetStr = "Powheg"
    elif MCSetIn == "Mad":     MCSetStr = "MadGraph"
    elif MCSetIn == "fr_Mad":     MCSetStr = "fr_MadGraph"
    elif MCSetIn == "fr_Pow":     MCSetStr = "fr_Powheg"

    if SavePlot:
        c1.SaveAs("Plot/CrossSection/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".pdf")        
        c1.SaveAs("Plot/CrossSection/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".eps")        
        c1.SaveAs("Plot/CrossSection/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".png")        
        c1.SaveAs("Plot/CrossSection/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".root")        

        if DoNormalized:
            c1.SaveAs("~/www/VBS/test/CrossSections/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+"_norm.pdf")        
            c1.SaveAs("~/www/VBS/test/CrossSections/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+"_norm.eps")        
            c1.SaveAs("~/www/VBS/test/CrossSections/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+"_norm.png")        
            c1.SaveAs("~/www/VBS/test/CrossSections/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+"_norm.root")       
        else:
            c1.SaveAs("~/www/VBS/test/CrossSections/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".pdf")        
            c1.SaveAs("~/www/VBS/test/CrossSections/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".eps")        
            c1.SaveAs("~/www/VBS/test/CrossSections/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".png")        
            c1.SaveAs("~/www/VBS/test/CrossSections/DiffCrossSecZZTo"+hMCList[i]["name"]+Type+PlotType+"_"+MCSetStr+".root")  
   
     
        
