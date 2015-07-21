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
                  help="Mass or Jets.")

parser.add_option("-s", "--save", dest="SavePlot",
                  default="False",
                  help="Save plot option, default is False")

parser.add_option("-S", "--Set", dest="Set",
                  default="Pow",
                  help="MC samples Set, default is Pow (Powheg) the other one is Mad (MadGraph)")

(options, args) = parser.parse_args()

Type = options.Type
Set = options.Set
SavePlot  = ast.literal_eval(options.SavePlot)

if Set=="Pow": SignalSamples = SignalSamples_Pow
elif Set=="Mad": SignalSamples = SignalSamples_Mad
else: 
    print Set,"is a wrong MC set, chose between Pow and Mad"
    sys.exit()

try:
    os.stat("./Plot/Acceptance/")
except:
    os.mkdir("./Plot/Acceptance/")


def SetAcceptance(inputdir,SampleType, SavePlot,SistErr):
    print "MC set",Set
    
    if SistErr == "0":     FileOut =  ROOT.TFile("Acceptance_"+Set+".root","UPDATE") 
    elif  SistErr == "1":     FileOut =  ROOT.TFile("AcceptanceSFactorPlus_"+Set+".root","UPDATE") 
    elif  SistErr == "-1":     FileOut =  ROOT.TFile("AcceptanceSFactorMinus_"+Set+".root","UPDATE") 
    elif  SistErr == "+Sq":     FileOut =  ROOT.TFile("AcceptanceSFactorSqPlus_"+Set+".root","UPDATE") 
    elif  SistErr == "-Sq":     FileOut =  ROOT.TFile("AcceptanceSFactorSqMinus_"+Set+".root","UPDATE") 

    if SampleType=="Mass":
        Xbins = sorted([100,200,250,300,350,400,500,600,800])
    elif SampleType =="Jets":
        Xbins = sorted([0,1,2,3,4])
    else: 
        print "Wrong Type"

    runArrayX = array('d',Xbins)

    c1 = TCanvas( 'c1', 'c1', 200, 10, 900, 700 )
    leg = ROOT.TLegend(.75,.62,.95,.76);    

    ROOT.gStyle.SetOptStat(0) 

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
            fileIn = ROOT.TFile(inputdir+sample["sample"]+".root")
            hGen_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"Gen_01")
          

            if SistErr=="0":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenReco_01")
            elif SistErr=="1":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFErrPlus_01")
            elif SistErr=="-1":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFErrMinus_01")
            elif SistErr=="+Sq":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFErrSqPlus_01")
            elif SistErr=="-Sq":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFErrSqMinus_01")
            else:          
                sys.exit("Wrong vaue for scale factor sistematic")

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
        TotAcc = TParameter(float)("TotAcc"+Fin,Acc)
        
        print Blue("Total Acceptance {0} {1:.2f} \n".format(29*" ",Acc))
       
        hMergReco.Divide(hMergGen)
        FileOut.cd()
        TotAcc.Write("",TotAcc.kOverwrite)
        hMergReco.Write("",hMergReco.kOverwrite)

        if SavePlot:
            c2 = TCanvas( 'c2', 'c2', 200, 10, 900, 700 )
            if SampleType=="Mass":
                hMergReco.SetXTitle("M_{ZZ} gen [Gev]")
                hMergReco.SetMaximum(0.515); 
                
            elif SampleType=="Jets": 
                hMergReco.SetXTitle("# Jets")
                hMergReco.GetXaxis().SetTitle("#Jets")
                hMergReco.GetXaxis().SetBinLabel(1,"0 Jets")
                hMergReco.GetXaxis().SetBinLabel(2,"1 Jets")
                hMergReco.GetXaxis().SetBinLabel(3,"2 Jets")
                hMergReco.GetXaxis().SetBinLabel(4,">2 Jets")
                hMergReco.SetMaximum(0.55); 
                hMergReco.SetMinimum(0.29); 
                
                
            hMergReco.SetYTitle("A #upoint #epsilon")
            hMergReco.SetMarkerStyle(21)
            hMergReco.SetMarkerSize(1.2)
            hMergReco.SetLineColor(1)
            hMergReco.SetFillColor(ROOT.kAzure-4)
            hMergReco.Draw("hist");
         
            if SistErr=="0":  c2.SaveAs("Plot/Acceptance/Acceptance_for_"+Fin+"_"+SampleType+".png")
            
            c1.cd()
            
            if Fin=="2e2m":
                hcopy1 = copy.deepcopy(hMergReco)
                hcopy1.SetMarkerColor(2)
                hcopy1.SetNameTitle("TotAcc","Differential Acceptance "+SampleType)
                leg.AddEntry(hcopy1,Fin,"lep")
                hcopy1.Draw()
            elif Fin=="4m":
                hcopy2 = copy.deepcopy(hMergReco)
                hcopy2.SetMarkerColor(3)
                leg.AddEntry(hcopy2,Fin,"lep")
                hcopy2.Draw("same");
            elif Fin=="4e":
                hcopy3 = copy.deepcopy(hMergReco)
                hcopy3.SetMarkerColor(4)
                leg.AddEntry(hcopy3,Fin,"lep")
                hcopy3.Draw("same");
                leg.Draw("same")
        
    if SistErr=="0":
        c1.SaveAs("./Plot/Acceptance/DiffAcceptance_"+SampleType+"_"+Set+".png")    
        #c1.SaveAs("./Plot/Acceptance/DiffAcceptance_"+SampleType+"_"+Set+".pdf")   
        #c1.SaveAs("./Plot/Acceptance/DiffAcceptance_"+SampleType+"_"+Set+".root") 
        #c1.SaveAs("./Plot/Acceptance/DiffAcceptance_"+SampleType+"_"+Set+".eps") 

    return FinStateAcc

#PlusAcc    = SetAcceptance("results/ZZMCAnalyzer_MC/",Type,SavePlot,"1")  # Correlated Errors
#MinusAcc   = SetAcceptance("results/ZZMCAnalyzer_MC/",Type,SavePlot,"-1") # Correlated Errors


PlusSqAcc   = SetAcceptance("results/ZZMCAnalyzer_MC/",Type,SavePlot,"+Sq")
MinusSqAcc  = SetAcceptance("results/ZZMCAnalyzer_MC/",Type,SavePlot,"-Sq")
CentralAcc  = SetAcceptance("results/ZZMCAnalyzer_MC/",Type,SavePlot,"0")


print Red(("Final total results for MC set {0} \n").format(Set))
for fin in ("2e2m","4m","4e"):
    print "acceptance for {0} {1} {2:.2f} + {3:.3f} %  - {4:.3f} % ".format(fin, (4-len(fin))*" ", CentralAcc[fin],(-1+PlusSqAcc[fin]/CentralAcc[fin])*100,(1-MinusSqAcc[fin]/CentralAcc[fin])*100)
    #print "Correlated errors {0} + {1:.5f} % - {2:.5f} %".format(10*" ",(-1+PlusAcc[fin]/CentralAcc[fin])*100,(1-MinusAcc[fin]/CentralAcc[fin])*100)
