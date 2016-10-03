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

#set the tdr style
tdrstyle.setTDRStyle()

#change the CMS_lumi variables (see CMS_lumi.py)
#CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
#CMS_lumi.lumi_8TeV = "19.7 fb^{-1}"
lumi = round(Lumi/1000.,1)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "{0}".format(lumi)+" fb^{-1} (13 TeV)\n"

iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 0


if DoInclusive:
    hMCList = getCrossPlot_MC(MCSetIn,"Total",analysis,DoNormalized,DoFiducial)
    TotalCross(MCSetIn,"Tot",analysis,DoFiducial,UseMCReco)
    sys.exit(0)


if MCSetIn == "Pow": mcop = "Mad"
elif MCSetIn == "Mad" : mcop = "Pow"
elif MCSetIn == "fr_Pow" : mcop = "fr_Mad"
elif MCSetIn == "fr_Mad" : mcop = "fr_Pow"

hMCList = getCrossPlot_MC(MCSetIn,Type,analysis,DoNormalized,DoFiducial)
hMCList_op = getCrossPlot_MC(mcop,Type,analysis,DoNormalized,DoFiducial)


if UseMCReco:
    (hDataList,hDataListUp,hDataListDown) = getCrossPlot_Data(MCSetIn,UseUnfold,Type,analysis,0,True,DoNormalized,DoFiducial) 

else:
    (hDataList,hDataListUp,hDataListDown)  = getCrossPlot_Data(MCSetIn,UseUnfold,Type,analysis,0,False,DoNormalized,DoFiducial) 



for i in range(0,4):
    hDataList[i]["state"].SetMarkerSize(1.3)
    hDataList[i]["state"].SetMarkerColor(1)
    hDataList[i]["state"].SetMarkerStyle(20)
    hDataList[i]["state"].SetLineColor(1)
    
    grSyst = getUncGraph(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],MCSetIn,hDataList[i]["name"],True,"statsyst")
    grSyst.SetFillStyle(3005)
    grSyst.SetFillColor(ROOT.kBlack)
 
    Err=ROOT.Double(0.)
    
    LastBin = hMCList[i]["state"].GetXaxis().GetLast()    

    # hRatio= ROOT.TH1F()
    # hRatio = hDataList[i]["state"].Clone("hRatioNew")
    # hRatio.SetStats(0)
    # hRatio.Divide(hMCList[i]["state"])
   
    if "Pow" in MCSetIn: 
        mcref_entry = "Powheg+MCFM+Phantom"
        mcop_entry = "MadGraph5_aMCatNLO+MCFM+Phantom"
    else:
        mcref_entry = "MadGraph5_aMCatNLO+MCFM+Phantom"
        mcop_entry = "Powheg+MCFM+Phantom"
    
###########################################################
   
    #Ratio and uncertainties on data points for the reference set of samples 
          
    #statistical uncertainty
    (hRatio,hRatio_up,hRatio_down) = LL(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],hMCList[i]["state"],Type,hDataList[i]["name"],"stat",True)     
    grStat_ratio = getUncGraph(hRatio,hRatio_up,hRatio_down,MCSetIn,hDataList[i]["name"],True,"syst")  
    
    #systematic uncertainty
    (hRatio_ss,hRatio_up_ss,hRatio_down_ss) = LL(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],hMCList[i]["state"],Type,hDataList[i]["name"],"statsyst",True) 
    grSist_ratio = getUncGraph(hRatio_ss,hRatio_up_ss,hRatio_down_ss,MCSetIn,hDataList[i]["name"],True,"syst")  
    grSist_ratio.SetFillStyle(3005)
    grSist_ratio.SetFillColor(ROOT.kBlack)

    maxtot = hRatio_up_ss.GetMaximum()   

###################################################################

    #Ratio and uncertainties on data for the opposite set of samples 
 
    #statistical uncertainty
    (hRatio_op,hRatio_op_up,hRatio_op_down) = LL(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],hMCList_op[i]["state"],Type,hDataList[i]["name"],"stat",True)     
    grStat_ratio_op = getUncGraph(hRatio_op,hRatio_op_up,hRatio_op_down,MCSetIn,hDataList[i]["name"],True,"syst")   
  
    #systematic uncertainty
    (hRatio_op_ss,hRatio_op_up_ss,hRatio_op_down_ss) = LL(hDataList[i]["state"],hDataListUp[i]["state"],hDataListDown[i]["state"],hMCList_op[i]["state"],Type,hDataList[i]["name"],"statsyst",True)  
    grSist_ratio_op = getUncGraph(hRatio_op_ss,hRatio_op_up_ss,hRatio_op_down_ss,MCSetIn,hDataList[i]["name"],True,"syst")  
    grSist_ratio_op.SetFillStyle(3005)
    grSist_ratio_op.SetFillColor(ROOT.kBlack)

    maxtot_op =  hRatio_op_up_ss.GetMaximum() 
    maxratios_op = max(maxtot,maxtot_op)


    if DoFiducial: Unit = "fb"
    else:          Unit = "pb"

    if "Central" in Type: eta_cut = " (|#eta^{jet}| < 2.4)"
    else: eta_cut = " (|#eta^{jet}| < 4.7)"

    if Type=="Mass":
        xName = "m_{4l} [GeV]"
        hRatio.GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{dm_{4l}} [1/GeV]")
        else: hMCList[i]["state"].GetYaxis().SetTitle("#frac{d#sigma}{dm_{4l}} ["+Unit+"/GeV]")

    elif "nJets" in Type:
        hRatio_op.GetXaxis().SetBinLabel(1,"0")
        hRatio_op.GetXaxis().SetBinLabel(2,"1")
        hRatio_op.GetXaxis().SetBinLabel(3,"2")
        hRatio_op.GetXaxis().SetBinLabel(3,"3")
        hRatio_op.GetXaxis().SetBinLabel(4,"#geq 4")

        hMCList[i]["state"].GetXaxis().SetBinLabel(1,"0")
        hMCList[i]["state"].GetXaxis().SetBinLabel(2,"1")
        hMCList[i]["state"].GetXaxis().SetBinLabel(3,"2")
        hMCList[i]["state"].GetXaxis().SetBinLabel(3,"3")
        hMCList[i]["state"].GetXaxis().SetBinLabel(4,"#geq 4")
        
        xName = "N_{jets}"+eta_cut
       
        hRatio.GetXaxis().SetTitle(xName)
        hMCList[i]["state"].GetXaxis().SetTitle(xName)

        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{dN_{jets}}")        
        else: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/dN_{jets} ["+Unit+"]") 

    elif "Deta" in Type:
        xName = "#Delta#eta_{jj}" + eta_cut
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{d#Delta#eta_{jj}}")
        else: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d#Delta#eta_{jj}["+Unit+"]")  

    elif "Mjj" in Type:
        xName = "m_{jj}" + eta_cut+ " [GeV]"
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{dm_{jj}} [1/GeV]")      
        else: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/dm_{jj} ["+Unit+"/GeV]")  

    elif "PtJet1" in Type:
        xName = "p_{T}^{jet1} [GeV]"
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{dp_{T}^{jet1}} [1/GeV]")             
        else:  hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/dp_{T}^{jet1} ["+Unit+"/GeV]") 
    elif "PtJet2" in Type:
        xName = "p_{T}^{jet2} [GeV]"
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{dp_{T}^{jet2}} [1/GeV]")             
        else:  hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/dp_{T}^{jet2} ["+Unit+"/GeV]") 
    elif "EtaJet1" in Type:
        xName = "|#eta^{jet1}|"
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{d|#eta^{jet1}|}")             
        else: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d|#eta^{jet1}| ["+Unit+"]") 
    elif "EtaJet2" in Type:
        xName = "|#eta^{jet2}|"
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{d|#eta^{jet2}|}")             
        else: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d|#eta^{jet2}| ["+Unit+"]")
    elif "dRZZ" in Type:
        xName = "#Delta R_{ZZ}"
        hMCList[i]["state"].GetXaxis().SetTitle(xName)
        if DoNormalized: hMCList[i]["state"].GetYaxis().SetTitle("#frac{1}{#sigma_{fid}} #frac{d#sigma_{fid}}{d #Delta R_{ZZ}}")             
        else: hMCList[i]["state"].GetYaxis().SetTitle("d#sigma/d #Delta R_{ZZ} ["+Unit+"]")



    hRatio.GetXaxis().SetTitle(xName) 
    hRatio_op.GetXaxis().SetTitle(xName)
    ratioName = "Data/MC"
    hRatio.GetYaxis().SetTitle(ratioName) 
    hRatio_op.GetYaxis().SetTitle(ratioName)
    
    maxtot_MGatNLO = 0
    if "EtaJet2" in Type: maxratios = max(maxratios_op,maxtot_MGatNLO)+0.35*max(maxratios_op,maxtot_MGatNLO)
    elif "PtJet1" in Type: maxratios = max(maxratios_op,maxtot_MGatNLO)+0.30*max(maxratios_op,maxtot_MGatNLO)
    elif "PtJet2" in Type: maxratios = max(maxratios_op,maxtot_MGatNLO)+0.35*max(maxratios_op,maxtot_MGatNLO)
    elif "Mjj" in Type: maxratios = max(maxratios_op,maxtot_MGatNLO)+0.35*max(maxratios_op,maxtot_MGatNLO) 
    elif "Deta" in Type: maxratios = max(maxratios_op,maxtot_MGatNLO)+0.50*max(maxratios_op,maxtot_MGatNLO)
    else: maxratios = max(maxratios_op,maxtot_MGatNLO)+0.80*max(maxratios_op,maxtot_MGatNLO)
 

    # "Type":{fr:{norm,notnorm},notfr:{norm,notnorm}}
    MaxList = {"Mass":{True:{True:0.009,False:0.35},False:{True:0.009,False:0.15}},"nJets":{True:{True:0.8,False:28.5},False:{True:0.85,False:13}},"nJets_Central":{True:{True:.82,False:50},False:{True:0.85,False:13}},"Mjj":{True:{True:0.004,False:0.02},False:{True:0.004, False:0.0086}},"Mjj_Central":{True:{True:0.005,False:0.018},False:{True:0.005,False:0.006}},"Deta":{True:{True:0.6,False:2.2},False:{True:0.45,False:1.2}},"Deta_Central":{True:{True:0.6,False:1.98},False:{True:0.7,False:.9}} ,"PtJet1":{True:{True:0.03,False:0.48},False:{True:0.06,False:0.16}},"PtJet2":{True:{True:0.016,False:0.08},False:{True:0.016,False:0.025}},"EtaJet1":{True:{True:0.5,False:7.4},False:{True:0.6,False:2.8}},"EtaJet2":{True:{True:0.5,False:2.6},False:{True:0.58,False:0.9}},"dRZZ":{True:{True:0.65,False:23.6},False:{True:0.7,False:12.}}}

    # Max   = hDataList[i]["state"].GetMaximum()+hDataList[i]["state"].GetBinError(hDataList[i]["state"].GetMaximumBin()) #Automatic max
    Max = MaxList[Type][DoFiducial][DoNormalized]*1.5

    hMCList[i]["state"].SetMaximum(Max)
    hMCList[i]["state"].SetMinimum(0) 
    hMCList[i]["state"].SetFillColor(851)
    hMCList[i]["state"].SetLineColor(852)
    hMCList[i]["state"].SetLineStyle(1)
    hMCList[i]["state"].SetLineWidth(2) 
    hMCList[i]["state"].GetXaxis().SetTitle("")
    hMCList[i]["state"].GetYaxis().SetTitleOffset(1.4) 
    hMCList[i]["state"].GetXaxis().SetLabelOffset(0.5) 
    hMCList[i]["state"].SetTitle("") 
    hMCList_op[i]["state"].SetLineColor(ROOT.kRed); 
   # hMCList_op[i]["state"].SetFillStyle(307);
    hMCList_op[i]["state"].SetLineWidth(2);
    hMCList_op[i]["state"].SetLineStyle(1);
  #  hMCList[i]["state"].SetFillColor(ROOT.kCyan-3)
   # hMCList[i]["state"].SetLineColor(1)
    #hMCList[i]["state"].GetXaxis().SetTitle(xName)
    #hMCList[i]["state"].SetTitle(hMCList[i]["name"]+" Cross Section")
     
       
    c1 = TCanvas( 'c1', hMCList[i]["name"] , 200, 10, 900, 1200 )
  

    #leg = ROOT.TLegend(.65,.52,.85,.66);

    #leg.AddEntry(grSyst,"Syst","f")
    #leg.AddEntry(hDataList[i]["state"],"Data + Stat","lep")
    #leg.AddEntry(hMCList[i]["state"],"MC","lpf")
    leg = ROOT.TLegend(.51,.50,.80,.85);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.032);
    leg.AddEntry(hDataList[i]["state"],"Unfolded data + stat. uncertainty","lep")
    leg.AddEntry(grSyst,"Total uncertainty","f") 
    leg.AddEntry(hMCList[i]["state"],mcref_entry,"f") 
    leg.AddEntry(hMCList_op[i]["state"],mcop_entry,"l")

    fInt = ROOT.TF1("constant","1",hMCList[i]["state"].GetXaxis().GetXmin(),    hMCList[i]["state"].GetXaxis().GetXmax());
        

    if(Type == "Mass"): 
        line = ROOT.TLine(100,1,800,1)
        line_op = ROOT.TLine(100,1,800,1)
    elif(Type == "nJets" or Type == "nJets_Central"): 
        line =  ROOT.TLine(0,1,4,1)
        line_op =  ROOT.TLine(0,1,4,1)
    elif(Type == "Mjj" or Type == "Mjj_Central"): 
        line =  ROOT.TLine(0,1,800,1)
        line_op =  ROOT.TLine(0,1,800,1)
    elif(Type == "Deta"  or Type == "Deta_Central"): 
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
    elif(Type =="dRZZ" ):  
        line =  ROOT.TLine(30,1,200,1)  #FIXME
        line_op =  ROOT.TLine(30,1,200,1) 
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

    #needed because of strange effects on pdf plots
    c1.cd()
    pad5 = ROOT.TPad ('space', 'space', 0.5, 0.13,  1., 0.20)#0.02-0.22
    pad5.SetTopMargin (0.0)
    pad5.SetRightMargin (0.03)#0.10
    pad5.SetLeftMargin (0.17)
    pad5.SetBottomMargin(0.60)#0.35
    pad5.Draw()
  
    pad1.cd() 
    hMCList[i]["state"].Draw("hist")
    hMCList_op[i]["state"].Draw("hist, same")

    if "Mjj" or "Pt" in Type:hMCList[i]["state"].GetYaxis().SetTitleOffset(1.2)
    else: hMCList[i]["state"].GetYaxis().SetTitleOffset(1)
    hMCList[i]["state"].GetYaxis().SetTitleSize(0.06)
   
    if not UseMCReco:
        grSyst.Draw("same2")
    hDataList[i]["state"].Draw("sameE1")
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
    elif Type == "Central_Mjj" or "Pt" in Type: hRatio.SetMinimum(-1.);
    else: hRatio.SetMinimum(-0.5);
#fInt.Draw("same")
   
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
    else: ratioTitle.DrawText(0.435, 0.90,mcref_entry);  
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
    elif Type == "Central_Mjj" or "Pt" in Type: hRatio_op.SetMinimum(-1.);
    else: hRatio_op.SetMinimum(-0.5);
    
    hRatio_op.Draw()
    line_op.Draw("same")
    grSist_ratio_op.Draw("same2") 
    grStat_ratio_op.Draw("P same")
    hRatio_op.Draw("same") 
    #if "Mad" in MCSetIn: grSistMC_ratio_op.Draw("same2")  #uncomment
    grSist_ratio_op.Draw("same2")#uncomment
   
    ratioTitle_op = ROOT.TText();
    ratioTitle_op.SetNDC();
    ratioTitle_op.SetTextColor(1);
    ratioTitle_op.SetTextSize(0.10);
    ratioTitle_op.SetTextAlign(22);
    ratioTitle_op.SetTextAngle(0);
    if "Pow" in MCSetIn: ratioTitle_op.DrawText(0.435, 0.90,mcop_entry);
    else: ratioTitle_op.DrawText(0.345, 0.90,mcop_entry); 
    hRatio_op.Draw("same")
    grSist_ratio_op.Draw("same2") 
    grStat_ratio_op.Draw("P same")
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
    
    #ROOT.gStyle.SetOptStat(0);   
    #c1.SetLogy()    

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
   
     
