#! /usr/bin/env python
from optparse import OptionParser
import ROOT,copy
ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import TH1F,TCanvas, TLegend, TParameter
import CrossInfo
from CrossInfo import*
from array import*
import sys,ast,os
import math
import operator
import collections
from Colours import *


parser = OptionParser(usage="usage: %prog <final state> [options]")


parser.add_option("-t", "--type", dest="Type",
                  default="Mass",
                  help="Mass,Total or Jets.")

parser.add_option("-s", "--save", dest="SavePlot",
                  default="False",
                  help="Save plot option, default is False")

parser.add_option("-S", "--Set", dest="Set",
                  default="Pow",
                  help="MC samples Set, default is Pow (Powheg) the other one is Amc (Amcatnlo)")

(options, args) = parser.parse_args()

Type = options.Type
Set = options.Set
SavePlot  = ast.literal_eval(options.SavePlot)

if "Pow" in Set: SignalSamples = SignalSamples_Pow
#elif "Mad" in Set: SignalSamples = SignalSamples_Mad
elif "Amc" in Set: SignalSamples = SignalSamples_Amc
else: sys.exit(Set,"is a wrong MC set, chose between Pow and Amc")

try:
    os.stat("./Plot/Acceptance/")
except:
    os.mkdir("./Plot/Acceptance/")
    
def SetAcceptance(inputdir,SampleType, SavePlot,SistErr):

    if Set=="Pow": print "\nPowheg Monte Carlo set"
    elif  Set=="Amc": print "\nAmcatNLO Monte Carlo set"  

    if SistErr == "0":     FileOut =  ROOT.TFile("Acceptance_"+Set+".root","UPDATE") 
    elif  SistErr == "1":     FileOut =  ROOT.TFile("AcceptanceSFactorPlus_"+Set+".root","UPDATE") 
    elif  SistErr == "-1":     FileOut =  ROOT.TFile("AcceptanceSFactorMinus_"+Set+".root","UPDATE") 
    elif  SistErr == "+Sq":     FileOut =  ROOT.TFile("AcceptanceSFactorSqPlus_"+Set+".root","UPDATE") 
    elif  SistErr == "-Sq":     FileOut =  ROOT.TFile("AcceptanceSFactorSqMinus_"+Set+".root","UPDATE") 
    
 
    if SampleType=="Mass" or SampleType=="Total":
        Xbins = sorted([100,200,250,300,350,400,500,600,800])
    elif "Jets" in SampleType:
        Xbins = sorted([0,1,2,3,4]) 
    elif SampleType=="PtJet1":
        Xbins = sorted([30,50,100,200,300,500])
    elif SampleType=="PtJet2":
        Xbins = sorted([30,100,200,300,500]) 
    elif "EtaJet" in SampleType:
        Xbins = sorted([0,2,4,6])
    elif "Deta" in SampleType:
        Xbins = sorted([0,2.4,4.7])
    elif "Mjj" in SampleType:
        Xbins = sorted([0.,200,800])
    else: 
        print "Wrong Type"

    runArrayX = array('d',Xbins)

    c1 = TCanvas( 'c1', 'c1', 200, 10, 900, 700 )
    ROOT.gStyle.SetOptStat(0) 

    if SampleType=="Mass":
        leg = ROOT.TLegend(.15,.15,.30,.29); 
    else:
        if "fr" in Set:
            leg = ROOT.TLegend(.15,.15,.30,.29);    
        else:        
            leg = ROOT.TLegend(.15,.70,.30,.84);    

    Acc1=0. 
    Acc2=0.
    Acc3=0.


    FinStateAcc ={"2e2m":Acc1,"4m":Acc2,"4e":Acc3}
    for Fin in ("2e2m","4m","4e"):
        print Red("\n#### FINAL STATE {0:s} Scale factor sist {1} ####\n".format(Fin,SistErr))

        hMergReco = ROOT.TH1F("HAcc_"+Fin+"_"+SampleType,"Acceptance for "+Fin+" vs "+SampleType, len(runArrayX)-1, runArrayX)
        hMergGen = ROOT.TH1F("HGen_"+Fin+"_"+SampleType,"Gen events for "+Fin+" vs "+SampleType, len(runArrayX)-1, runArrayX)
        print Blue("Sample                                         Acceptance \n")
        for sample in SignalSamples:
            #print "sample ",sample
            fileIn = ROOT.TFile(inputdir+sample["sample"]+".root")
            
            if "fr" in Set:
                hGen_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"Gen_01_fr") #if tighter fiducial region
                if SistErr=="0":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenReco_01_fr")
                elif SistErr=="1":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFErrPlus_01_fr")
                elif SistErr=="-1":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFErrMinus_01_fr")
                elif SistErr=="+Sq":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFErrSqPlus_01_fr")
                elif SistErr=="-Sq":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFErrSqMinus_01_fr") 
                else:          
                    sys.exit("Wrong value for scale factor sistematic")
            
            else:
                hGen_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"Gen_01")
                if SistErr=="0":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenReco_01")
                elif SistErr=="1":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFErrPlus_01")
                elif SistErr=="-1":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFErrMinus_01")
                elif SistErr=="+Sq":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFErrSqPlus_01")
                elif SistErr=="-Sq":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFErrSqMinus_01") 
                else:          
                    sys.exit("Wrong value for scale factor sistematic")

            hGen  = copy.deepcopy(hGen_)          
            hReco = copy.deepcopy(hRec_)
            if hGen==None: 
                # print "in sample",sample["sample"],"ZZTo"+Fin+"_"+SampleType+"Gen_01"+"   Not found"
                continue
            if hReco==None:
                print sample["sample"],"No events in reco sample but events in gen sample"
                continue
            n = len(sample["sample"])
            spacestring = (45-n)*" "
            print "{0} {1} {2:.2f}".format(sample["sample"],spacestring,hReco.Integral(1,-1)/hGen.Integral(1,-1))

            hMergReco.Add(hReco)
            hMergGen.Add(hGen)             

        Nbins= hMergGen.GetNbinsX()
        
        NTotEv  = hMergGen.Integral(1,-1)
        NRecoEv = hMergReco.Integral(1,-1)
        Acc=NRecoEv/NTotEv
        FinStateAcc[Fin]=Acc

        
        print Blue("Total Acceptance {0} {1:.2f} \n".format(29*" ",Acc))
       
        hMergReco.Divide(hMergGen)
        FileOut.cd()
        if SampleType=="Mass": 
            TotAcc = TParameter(float)("TotAcc"+Fin,Acc)
            TotAcc.Write("",TotAcc.kOverwrite)
        hMergReco.Write("",hMergReco.kOverwrite)

        if SavePlot:
            c2 = TCanvas( 'c2', 'c2', 200, 10, 900, 700 )
            if SampleType=="Mass":
                hMergReco.SetXTitle("gen M_{ZZ} [Gev]")
                                             
            elif SampleType=="Jets": 
                hMergReco.SetXTitle("# Jets")
                hMergReco.GetXaxis().SetTitle("#Jets")
                hMergReco.GetXaxis().SetBinLabel(1,"0 Jets")
                hMergReco.GetXaxis().SetBinLabel(2,"1 Jets")
                hMergReco.GetXaxis().SetBinLabel(3,"2 Jets")
                hMergReco.GetXaxis().SetBinLabel(4,">2 Jets")
               
            elif SampleType=="CentralJets": 
                hMergReco.SetXTitle("# Central Jets")
                hMergReco.GetXaxis().SetTitle("#CentralJets")
                hMergReco.GetXaxis().SetBinLabel(1,"0 Jets")
                hMergReco.GetXaxis().SetBinLabel(2,"1 Jets")
                hMergReco.GetXaxis().SetBinLabel(3,"2 Jets")
                hMergReco.GetXaxis().SetBinLabel(4,">2 Jets")
               
            elif SampleType == "Deta":
                hMergReco.SetXTitle("gen #Delta#eta_{jj}")
             
            elif SampleType == "CentralDeta":
                hMergReco.SetXTitle("gen #Delta#eta_{jj} (|#eta^{jet}|<2.4)")
              
            elif SampleType == "Mjj":
                hMergReco.SetXTitle("gen M_{jj}")
              
            elif SampleType == "CentralMjj":
                hMergReco.SetXTitle("gen M_{jj} (|#eta^{jet}|<2.4)")
                          
            elif SampleType == "PtJet1":
                hMergReco.SetXTitle("gen p_{T}^{jet1}") 
                                      
            elif SampleType == "PtJet2":
                hMergReco.SetXTitle("gen p_{T}^{jet2}")
                                              
            elif SampleType == "EtaJet1":
                hMergReco.SetXTitle("gen #eta^{jet1}")
             
            elif SampleType == "EtaJet2":
                hMergReco.SetXTitle("gen #eta^{jet2}") 
             
            
            if "fr" in Set: 
                hMergReco.SetMaximum(1.);
                hMergReco.SetMinimum(0.29);   
                hMergReco.SetYTitle("#epsilon") 
            
            else:
                if SampleType == "Mass":
                    hMergReco.SetMaximum(0.8); 
                    hMergReco.SetMinimum(0.20);
                else:
                    hMergReco.SetMaximum(0.8);
                    hMergReco.SetMinimum(0.29);  
                hMergReco.SetYTitle("A #upoint #epsilon") 
            
            hMergReco.SetMarkerStyle(21)
            hMergReco.SetMarkerSize(1.2)
            hMergReco.SetLineColor(1)
            hMergReco.SetFillColor(ROOT.kAzure-4)
            hMergReco.Draw("hist");
         
            if SistErr=="0": 
                if "fr" in Set:  
                    c2.SaveAs("Plot/Acceptance/Acceptance_for_"+Fin+"_"+SampleType+"_fr.png")
                else:
                    c2.SaveAs("Plot/Acceptance/Acceptance_for_"+Fin+"_"+SampleType+".png")
                    c2.SaveAs("~/www/PlotsVV/13TeV/Acceptance/DiffAcceptance_"+Fin+"_"+SampleType+"_"+Set+".png")
            
            c1.cd()
            
            if Fin=="2e2m":
                hcopy1 = copy.deepcopy(hMergReco)
                hcopy1.SetMarkerColor(2)
                hcopy1.SetNameTitle("TotAcc","Differential Acceptance "+SampleType)
                leg.AddEntry(hcopy1,"ZZ #rightarrow 2e2#mu","lep")
                hcopy1.Draw()
            elif Fin=="4m":
                hcopy2 = copy.deepcopy(hMergReco)
                hcopy2.SetMarkerColor(3)
                leg.AddEntry(hcopy2,"ZZ #rightarrow 4#mu","lep")
                hcopy2.Draw("same");
            elif Fin=="4e":
                hcopy3 = copy.deepcopy(hMergReco)
                hcopy3.SetMarkerColor(4)
                leg.AddEntry(hcopy3,"ZZ #rightarrow 4e","lep")
                hcopy3.Draw("same");
                leg.Draw("same")
        
    if SistErr=="0":
        if "fr" in Set: 
            c1.SaveAs("./Plot/Acceptance/DiffAcceptance_"+SampleType+"_"+Set+"_fr.png") 
        else:        
            c1.SaveAs("./Plot/Acceptance/DiffAcceptance_"+SampleType+"_"+Set+".png") 
            c1.SaveAs("~/www/PlotsVV/13TeV/Acceptance/DiffAcceptance_"+SampleType+"_"+Set+".png") 

        
    return FinStateAcc

#PlusAcc    = SetAcceptance("results/ZZMCAnalyzer_MC/",Type,SavePlot,"1")  # Correlated Errors //To be added again
#MinusAcc   = SetAcceptance("results/ZZMCAnalyzer_MC/",Type,SavePlot,"-1") # Correlated Errors

PlusSqAcc   = SetAcceptance("results/ZZMCAnalyzer_MC/",Type,SavePlot,"+Sq")
MinusSqAcc  = SetAcceptance("results/ZZMCAnalyzer_MC/",Type,SavePlot,"-Sq")
CentralAcc  = SetAcceptance("results/ZZMCAnalyzer_MC/",Type,SavePlot,"0")


print Red(("Final total results for MC set {0} \n").format(Set))
for fin in ("2e2m","4m","4e"):
    print "acceptance for {0} {1} {2:.2f} + {3:.3f} %  - {4:.3f} % ".format(fin, (4-len(fin))*" ", CentralAcc[fin],(-1+PlusSqAcc[fin]/CentralAcc[fin])*100,(1-MinusSqAcc[fin]/CentralAcc[fin])*100)
    # print "Correlated errors {0} + {1:.5f} % - {2:.5f} %".format(10*" ",(-1+PlusAcc[fin]/CentralAcc[fin])*100,(1-MinusAcc[fin]/CentralAcc[fin])*100)
