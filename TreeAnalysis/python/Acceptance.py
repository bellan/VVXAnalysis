#! /usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
import ROOT,copy
ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import TH1F,TCanvas, TLegend, TParameter,TList
import CrossInfo
from CrossInfo import*
from PersonalInfo import*
from array import*
import sys,ast,os
import math
import operator
import collections
from Colours import *
from os.path import expanduser


parser = OptionParser(usage="usage: %prog <final state> [options]")

parser.add_option("-t", "--type", dest="Type",
                  default="Mass",
                  help="Mass, or nJets.")

parser.add_option("-s", "--save", dest="SavePlot",
                  action="store_true",
                  default=False,
                  help="Save plot option, default is False")

parser.add_option("-S", "--Set", dest="Set",
                  default="Mad",
                  help="MC samples Set, default is Mad (Amcatnlo) the other one is Pow (Powheg)")

parser.add_option("-A", "--Analysis", dest="Analysis",
                  default="ZZ",
                  help="Analysis, default is  ZZ, others are ZZFull and HZZ")

parser.add_option("-d", "--Dir", dest="Dir",
                  default = "test",
                  help="Choose a specific directory name for save plots. Default is ''")

parser.add_option("-p", "--plot", dest="showPlot",
                  action="store_true",
                  default=False,
                  help="Show plot option, default is False")


(options, args) = parser.parse_args()

Type      = options.Type
Set       = options.Set
SavePlot  = options.SavePlot
Analysis  = options.Analysis
Dir       = options.Dir
showPlot   = options.showPlot

if not showPlot: ROOT.gROOT.SetBatch(True)

Dir+="/"

if Analysis=="ZZ": Analysis=""
else: Analysis="_"+Analysis


if "Pow" in Set: SignalSamples = SignalSamples_Pow
elif "Mad" in Set: SignalSamples = SignalSamples_Mad
else: sys.exit(Set,"is a wrong MC set, chose between Pow and Mad")

#if Type != "Mass" and Type != "nJets": sys.exit("ERROR \nWrong Type, choose between Mass or nJets.")

try:
    os.stat("./Plot")
except:
    os.mkdir("./Plot")

try:
    os.stat("./Plot/Acceptance/")
except:
    os.mkdir("./Plot/Acceptance/")

try:
    os.stat("./Acceptance/")
except:
    os.mkdir("./Acceptance/")

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
        os.system("cp "+PersonalFolder+"index.php "+PersonalFolder+Dir+"/"+Type) 
    try:
        os.stat(PersonalFolder+Dir+"/"+Type+"/Acceptance/")
    except:
        os.mkdir(PersonalFolder+Dir+"/"+Type+"/Acceptance/")
        os.system("cp "+PersonalFolder+"index.php "+PersonalFolder+Dir+"/"+Type+"/Acceptance/" )


    
FileOut_0   =  ROOT.TFile("./Acceptance/Acceptance_"+Set+Analysis+"_"+Type+".root","RECREATE") 
FileOut_1   =  ROOT.TFile("./Acceptance/AcceptanceSFactorUp_"+Set+Analysis+"_"+Type+".root","RECREATE") 
FileOut_m1  =  ROOT.TFile("./Acceptance/AcceptanceSFactorDn_"+Set+Analysis+"_"+Type+".root","RECREATE") 
FileOut_Sq  =  ROOT.TFile("./Acceptance/AcceptanceSFactorSqUp_"+Set+Analysis+"_"+Type+".root","RECREATE") 
FileOut_mSq =  ROOT.TFile("./Acceptance/AcceptanceSFactorSqDn_"+Set+Analysis+"_"+Type+".root","RECREATE") 
    

FileOut_Acc1   =  ROOT.TFile("./Acceptance/AcceptanceTheorUp_"+Set+Analysis+"_"+Type+".root","RECREATE") 
FileOut_Accm1  =  ROOT.TFile("./Acceptance/AcceptanceTheorDn_"+Set+Analysis+"_"+Type+".root","RECREATE") 

def SetAcceptance(inputdir,SampleType, SavePlot,SistErr,Fr):

    if Fr =="Acc":    
        if    SistErr == "0":     FileOut = FileOut_0 
        elif  SistErr == "1":     FileOut = FileOut_Acc1 
        elif  SistErr == "-1":    FileOut = FileOut_Accm1 
    else:
        if    SistErr == "0":     FileOut = FileOut_0 
        elif  SistErr == "+Sq":   FileOut = FileOut_Sq 
        elif  SistErr == "-Sq":   FileOut = FileOut_mSq 
        elif  SistErr == "1":     FileOut = FileOut_1 
        elif  SistErr == "-1":    FileOut = FileOut_m1 


    if    Set=="Pow": print("\nPowheg Monte Carlo set")
    elif  Set=="Mad": print("\nMadgraph Monte Carlo set")  
 

    c1 = TCanvas( 'c1', 'c1', 200, 10, 900, 700 )
    ROOT.gStyle.SetOptStat(0) 

    if SampleType=="Mass":
        leg = ROOT.TLegend(.15,.15,.30,.29); 
    else:
        if Fr=="Eff_Tight_noOut":
            leg = ROOT.TLegend(.15,.15,.30,.29);    
        else:        
            leg = ROOT.TLegend(.15,.70,.30,.84);    

    Acc1=0. 
    Acc2=0.
    Acc3=0.

    Tot4lfr = 0;
    Tot4lReco = 0;
    Tot4l     = 0;

    FinStateAcc ={"2e2m":Acc1,"4m":Acc2,"4e":Acc3}
    for Fin in ("2e2m","4m","4e"):
        print(Red("\n#### FINAL STATE {0:s} Scale factor sist {1} ####\n".format(Fin,SistErr)))

        hMergReco = ROOT.TH1F()
        hMergGen  = ROOT.TH1F()
        
        print(Blue("Sample                                         Acceptance \n"))
        isFirst=True            
        for sample in SignalSamples:
            fileIn = ROOT.TFile(inputdir+sample["sample"]+".root")
            fileRecoIn = ROOT.TFile( ( (inputdir.replace("_MC","_SR")).replace("MC","Reco"))+sample["sample"]+".root")
            if Fr=="Tot":
                hGen_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"Gen_01")
                if SistErr=="0":      hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenReco_01")
                elif SistErr=="1":    hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFUp_01")
                elif SistErr=="-1":   hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFDn_01")
                elif SistErr=="+Sq":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFSqUp_01")
                elif SistErr=="-Sq":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFSqDn_01") 
                else:          
                    sys.exit("Wrong value for scale factor sistematic")

            
            elif Fr=="Eff_Tight":
                hGen_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"Gen_01_fr")
                if SistErr=="0":      hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenReco_01")
                elif SistErr=="1":    hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFUp_01")
                elif SistErr=="-1":   hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFDn_01")
                elif SistErr=="+Sq":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFSqUp_01")
                elif SistErr=="-Sq":  hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"GenRecoSFSqDn_01") 
                else:          
                    sys.exit("Wrong value for scale factor sistematic")


            elif Fr=="EffAndFake_Tight":
                hGen_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"Gen_01_fr")
                if SistErr=="0":      hRec_ = fileRecoIn.Get("ZZTo"+Fin+"_"+SampleType+"_01")
                elif SistErr=="1":    hRec_ = fileRecoIn.Get("ZZTo"+Fin+"_"+SampleType+"_SFUp_01")
                elif SistErr=="-1":   hRec_ = fileRecoIn.Get("ZZTo"+Fin+"_"+SampleType+"_SFDn_01")
                elif SistErr=="+Sq":  hRec_ = fileRecoIn.Get("ZZTo"+Fin+"_"+SampleType+"_SFSqUp_01")
                elif SistErr=="-Sq":  hRec_ = fileRecoIn.Get("ZZTo"+Fin+"_"+SampleType+"_SFSqDn_01") 
                else:          
                    sys.exit("Wrong value for scale factor sistematic")


            elif Fr=="Acc":                    
                if SistErr=="0": 
                    hGen_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"Gen_01")
                    hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"Gen_01_fr")
                if SistErr=="1":    
                    hGen_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"Gen_01_scaleUp")
                    hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"Gen_01_scaleUp_fr")
                if SistErr=="-1":
                    hGen_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"Gen_01_scaleDn")
                    hRec_ = fileIn.Get("ZZTo"+Fin+"_"+SampleType+"Gen_01_scaleDn_fr")


            hGen  = copy.deepcopy(hGen_)          
            hReco = copy.deepcopy(hRec_)
            if hGen==None: 
                print("in sample",sample["sample"],"ZZTo"+Fin+"_"+SampleType+"Gen_01"+"   Not found")
                continue
            if hReco==None:
                print(sample["sample"],"No events in reco sample but events in gen sample")
                continue
  
            n = len(sample["sample"])
            spacestring = (45-n)*" "
            print("{0} {1} {2:.4f}".format(sample["sample"],spacestring,hReco.Integral(1,-1)/hGen.Integral(1,-1)))
            #print "{0} {1} {2:.4f}, reco {3}, total {4}".format(sample["sample"],spacestring,hReco.Integral(1,-1)/hGen.Integral(1,-1),hReco.Integral(1,-1),hGen.Integral(1,-1))
            
            if isFirst:
                hMergReco = hReco
                hMergGen  = hGen
                isFirst=False            
            else:
                hMergReco.Add(hReco)
                hMergGen.Add(hGen)             
           
        Nbins   = hMergGen.GetNbinsX()
        NTotEv  = hMergGen.Integral(0,-1)
        NRecoEv = hMergReco.Integral(0,-1)


        Tot4lReco+=NRecoEv
        Tot4l    +=NTotEv

        Acc= NRecoEv/NTotEv
 #       print "Tot",NTotEv,"Reco",NRecoEv

        FinStateAcc[Fin]=Acc
        
        print(Blue("Total Acceptance {0} {1:.2f} \n".format(29*" ",Acc)))
        hMergReco.Divide(hMergGen)
        FileOut.cd()
        if SampleType=="Mass": 
            if Fr=="Eff_Tight":          TotAcc = TParameter(float)("TotAcc"+Fin+"_fr",Acc)
            elif Fr=="Tot":              TotAcc = TParameter(float)("TotAcc"+Fin+"_Tot",Acc)
            elif Fr=="Acc":              TotAcc = TParameter(float)("TotAcc"+Fin+"_Acc",Acc)
            elif Fr=="EffAndFake_Tight":    TotAcc = TParameter(float)("TotAcc"+Fin+"_frAndFake",Acc)
            else: sys.exit("Error: Wrong region ")
           
            TotAcc.Write("",TotAcc.kOverwrite)

        hMergReco.SetNameTitle("HAcc_"+Fin+"_"+SampleType,"Acceptance for "+Fin+" vs "+SampleType)
        hMergReco.SetName("H"+Fr+"_"+Fin+"_"+SampleType)
        hMergReco.Write("",hMergReco.kOverwrite)

        if SavePlot:
            c2 = TCanvas( 'c2', 'c2', 200, 10, 900, 700 )
            if SampleType=="Mass":
                hMergReco.SetXTitle("gen M_{ZZ} [Gev]")
                                             
            elif SampleType=="nJets": 
                hMergReco.SetXTitle("# nJets")
                hMergReco.GetXaxis().SetTitle("#nJets")
                hMergReco.GetXaxis().SetBinLabel(1,"0 Jets")
                hMergReco.GetXaxis().SetBinLabel(2,"1 Jets")
                hMergReco.GetXaxis().SetBinLabel(3,"2 Jets")
                hMergReco.GetXaxis().SetBinLabel(4,">2 Jets")

            elif SampleType=="nJets_Central": 
                hMergReco.SetXTitle("# Central nJets")
                hMergReco.GetXaxis().SetTitle("#CentralnJets")
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
              
            elif SampleType == "Mjj_Central":
                hMergReco.SetXTitle("gen M_{jj} (|#eta^{jet}|<2.4)")
                          
            elif SampleType == "PtJet1":
                hMergReco.SetXTitle("gen p_{T}^{jet1}") 
                                      
            elif SampleType == "PtJet2":
                hMergReco.SetXTitle("gen p_{T}^{jet2}")
                                              
            elif SampleType == "EtaJet1":
                hMergReco.SetXTitle("gen #eta^{jet1}")
             
            elif SampleType == "EtaJet2":
                hMergReco.SetXTitle("gen #eta^{jet2}") 
             
            
            if "Eff_Tight" in Fr: 
                hMergReco.SetMaximum(1.);
                hMergReco.SetMinimum(0.29);   
                hMergReco.SetYTitle("#epsilon") 
            
            elif Fr=="Acc":
                hMergReco.SetMaximum(0.9);
                hMergReco.SetMinimum(0.29);  
                hMergReco.SetYTitle("ACC.") 
            else:
                if SampleType == "Mass":
                    hMergReco.SetMaximum(0.9); 
                    hMergReco.SetMinimum(0.20);
                else:
                    hMergReco.SetMaximum(0.9);
                    hMergReco.SetMinimum(0.29);  
                hMergReco.SetYTitle("A #upoint #epsilon") 
            
            hMergReco.SetMarkerStyle(21)
            hMergReco.SetMarkerSize(1.2)
            hMergReco.SetLineColor(1)
            hMergReco.SetFillColor(ROOT.kAzure-4)
            hMergReco.Draw("hist");
         
#            if SistErr=="0": 
            c2.SaveAs("Plot/Acceptance/DiffAcceptance_"+Fin+"_"+SampleType+"_"+Set+Analysis+"_"+Fr+"_"+SistErr+".png")
            
            c1.cd()
            
            if Fin=="2e2m":
                hcopy1 = copy.deepcopy(hMergReco)
                hcopy1.SetMarkerColor(2)
                hcopy1.SetNameTitle("TotAcc","Differential Acceptance "+SampleType)
                leg.AddEntry(hcopy1,"ZZ #rightarrow 2e2#mu","lep")
                hcopy1.Draw("E1")
            elif Fin=="4m":
                hcopy2 = copy.deepcopy(hMergReco)
                hcopy2.SetMarkerColor(3)
                leg.AddEntry(hcopy2,"ZZ #rightarrow 4#mu","lep")
                hcopy2.Draw("sameE1");
            elif Fin=="4e":
                hcopy3 = copy.deepcopy(hMergReco)
                hcopy3.SetMarkerColor(4)
                leg.AddEntry(hcopy3,"ZZ #rightarrow 4e","lep")
                hcopy3.Draw("sameE1");
                leg.Draw("same")

    FileOut.cd
    Acc4l = Tot4lReco/Tot4l
    print("Acc 4l",SistErr,Acc4l)
    TotAcc = TParameter(float)("TotAcc4l_"+Fr,Acc4l)
    TotAcc.Write("",TotAcc.kOverwrite)
    c1.SaveAs("./Plot/Acceptance/DiffAcceptance_"+SampleType+"_"+Set+Analysis+"_"+Fr+".png") 
    if SavePlot:
        print(PersonalFolder+Dir+"/"+Type+"/Acceptance/DiffAcceptance_"+SampleType+"_"+Set+Analysis+"_"+Fr+"_"+SistErr+".png") 
        c1.SaveAs(PersonalFolder+Dir+"/"+Type+"/Acceptance/DiffAcceptance_"+SampleType+"_"+Set+Analysis+"_"+Fr+"_"+SistErr+".png") 
        c1.SaveAs(PersonalFolder+Dir+"/"+Type+"/Acceptance/DiffAcceptance_"+SampleType+"_"+Set+Analysis+"_"+Fr+"_"+SistErr+".pdf") 
    return FinStateAcc


#PlusAcc   = SetAcceptance("results/ZZMCAnalyzer_MC/",Type,SavePlot,"1")  # Correlated Errors //To be added again
#MinuAcc   = SetAcceptance("results/ZZMCAnalyzer_MC/",Type,SavePlot,"-1") # Correlated Errors


print(Yellow("Total"))
PlusSqAcc   = SetAcceptance("results/ZZMCAnalyzer_MC"+Analysis+"/",Type,SavePlot,"+Sq","Tot")            
MinusSqAcc  = SetAcceptance("results/ZZMCAnalyzer_MC"+Analysis+"/",Type,SavePlot,"-Sq","Tot")
CentralAcc  = SetAcceptance("results/ZZMCAnalyzer_MC"+Analysis+"/",Type,SavePlot,"0","Tot")

print(Yellow("\n Fiducial"))
PlusSqAcc_fr   = SetAcceptance("results/ZZMCAnalyzer_MC"+Analysis+"/",Type,SavePlot,"+Sq","Eff_Tight")
MinusSqAcc_fr  = SetAcceptance("results/ZZMCAnalyzer_MC"+Analysis+"/",Type,SavePlot,"-Sq","Eff_Tight")
CentralAcc_fr  = SetAcceptance("results/ZZMCAnalyzer_MC"+Analysis+"/",Type,SavePlot,"0","Eff_Tight")

print(Yellow("\n Fiducial and fake"))
PlusSqAcc_frAndFake   = SetAcceptance("results/ZZMCAnalyzer_MC"+Analysis+"/",Type,SavePlot,"+Sq","EffAndFake_Tight")
MinusSqAcc_frAndFake  = SetAcceptance("results/ZZMCAnalyzer_MC"+Analysis+"/",Type,SavePlot,"-Sq","EffAndFake_Tight")
CentralAcc_frAndFake  = SetAcceptance("results/ZZMCAnalyzer_MC"+Analysis+"/",Type,SavePlot,"0","EffAndFake_Tight")

print(Yellow("\n Acceptance Tr->Fr"))
CentralAcc_Acc  = SetAcceptance("results/ZZMCAnalyzer_MC"+Analysis+"/",Type,SavePlot,"0","Acc")
CentralAcc_Acc_up  = SetAcceptance("results/ZZMCAnalyzer_MC"+Analysis+"/",Type,SavePlot,"1","Acc")
CentralAcc_Acc_dn  = SetAcceptance("results/ZZMCAnalyzer_MC"+Analysis+"/",Type,SavePlot,"-1","Acc")

print(Red(("Final total results for MC set {0} \n").format(Set)))
print(Red("Total Region \nAcceptance x efficency"))
for fin in ("2e2m","4m","4e"):
    print("{0} {1} {2:.3f} + {3:.3f} %  - {4:.3f} % ".format(fin, (6-len(fin))*" ", CentralAcc[fin],(-1+PlusSqAcc[fin]/CentralAcc[fin])*100,(1-MinusSqAcc[fin]/CentralAcc[fin])*100))
    # print "Correlated errors {0} + {1:.5f} % - {2:.5f} %".format(10*" ",(-1+PlusAcc[fin]/CentralAcc[fin])*100,(1-MinusAcc[fin]/CentralAcc[fin])*100)

print(Red("Fiducial Region \nEfficency "))
for fin in ("2e2m","4m","4e"):
    print("{0} {1} {2:.4f} + {3:.3f}  - {4:.3f}".format(fin, (6-len(fin))*" ", CentralAcc_fr[fin],(PlusSqAcc_fr[fin]-CentralAcc_fr[fin]),(-MinusSqAcc_fr[fin]+CentralAcc_fr[fin])))

print(Red("Fiducial Region \nEfficency and fake correction "))
for fin in ("2e2m","4m","4e"):
   print("{0} {1} {2:.4f} + {3:.3f}  - {4:.3f}".format(fin, (6-len(fin))*" ", CentralAcc_frAndFake[fin],(PlusSqAcc_frAndFake[fin]-CentralAcc_frAndFake[fin]),(-MinusSqAcc_frAndFake[fin]+CentralAcc_frAndFake[fin])))


print(Red("From Tight to Total \nAcceptance "))
for fin in ("2e2m","4m","4e"):
    print("{0} {1} {2:.4f}".format(fin, (6-len(fin))*" ", CentralAcc_Acc[fin]))

print("")


print(Red("Fiducial Region \nEfficency "))
print(" \\begin{tabular}{ccc}")
print("  \\hline")
print("Final state & Efficiency  & acceptance \\\\")
print("  \\hline")
for fin in ("2e2m","4m","4e"):
    print("${0} $ & {1} {2:.1f} + {3:.1f}  - {4:.1f} \\% & {5:.1f} + {6:.1f} - {7:.1f} \\% \\\\".format(fin, (6-len(fin))*" ", 100*CentralAcc_fr[fin],100*(PlusSqAcc_fr[fin]-CentralAcc_fr[fin]),100*(-MinusSqAcc_fr[fin]+CentralAcc_fr[fin]), 100*CentralAcc_Acc[fin],  100*(CentralAcc_Acc_dn[fin]-CentralAcc_Acc[fin]), 100*(CentralAcc_Acc[fin]-CentralAcc_Acc_up[fin] )))
print("  \\hline")
print("\\end{tabular} ")
