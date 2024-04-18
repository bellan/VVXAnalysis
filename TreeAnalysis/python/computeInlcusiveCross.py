#! /usr/bin/env python

##################################
## G. Pinna (UNITO) - Jun 2015 ##
##################################

from __future__ import print_function
import ROOT,copy
from ROOT import gSystem, TCanvas, TH1,  TPad, gStyle, TLegend,TGraphAsymmErrors,Math,TArrayD
import collections 
import CrossInfo
from CrossInfo import*
from Colours import *
from optparse import OptionParser
import sys,ast,os
import math
import operator
import textwrap
import collections                                                                                    

import LikelihoodXS
from LikelihoodXS import*


parser = OptionParser(usage="usage: %prog <final state> [options]")

parser.add_option("-S", "--Set", dest="Set",
                  default="Pow",
                  help="MC samples Set, default is Pow (Powheg) the other one is Amc (Amcatnlo)")

parser.add_option("-A", "--Analysis", dest="Analysis",
                  default="ZZ",
                  help="Analysis, default is  ZZ, others are ZZFull and HZZ")


(options, args) = parser.parse_args()

Set = options.Set
Analysis  = options.Analysis

if Analysis!="ZZ":
    inputdir     = "./results/ZZRecoAnalyzer_SR_"+Analysis+"/"
    inputdir_CR  = "./results/ZZRecoAnalyzer_CR_"+Analysis+"/"
else:
    inputdir     = "./results/ZZRecoAnalyzer_SR/"
    inputdir_CR  = "./results/ZZRecoAnalyzer_CR/"


# def getRedBkg(FinState,Sign):

#         hFakeRate.Add(setErrorsEntries(fileFake.Get("ZZTo"+FinState+"_Mass_FRVar")),-1) #CHECK
#     if  Sign==-1:
#         hFakeRate = fileFake.Get("ZZTo"+FinState+"_Mass_01")
#         hFakeRate.Add(setErrorsEntries(fileFake.Get("ZZTo"+FinState+"_Mass_FRVar")),1)
#     else:
#         hFakeRate = fileFake.Get("ZZTo"+FinState+"_Mass_01")
#     Err=ROOT.Double(0.)
#     Integr= hFakeRate.IntegralAndError(0,-1,Err)
#     if hFakeRate.Integral(0,-1)<=0: 
#         hFakeRate.Scale(0)    

#     print "{0} contribution {1}> {2:.2f}\n".format(FinState,(33-len(FinState))*"-",hFakeRate.Integral(0,-1))
#     return  copy.deepcopy(hFakeRate)

#################################################################################################################

def setErrorsEntries(hVar):
    Nbins=hVar.GetNbinsX()
    for k in range(1,Nbins+1):	
        hVar.SetBinContent(k,math.sqrt(hVar.GetBinContent(k)))
#        print "integrale",hVar.Integral(0,-1)
    return hVar


##################################################################################################################

####################  Get, sum and set histograms from data or from MC reco as they were data ####################

##################################################################################################################


DataFile = ROOT.TFile(inputdir+"data.root")

hsum2e2mu = ROOT.TH1F()
hsum4e    = ROOT.TH1F()
hsum4mu   = ROOT.TH1F()

hSum = [{"hist":hsum2e2mu,"name":'2e2m'},{"hist":hsum4e,"name":'4e'},{"hist":hsum4mu,"name":'4m'}]
print(Red("\nData\n"))
for h in hSum:
    h1 = DataFile.Get("ZZTo"+h["name"]+"_Mass_01")
    h["hist"] = copy.deepcopy(h1)
    print("{0} contribution {1}> {2:.2f}\n".format(h["name"],(33-len(h["name"]))*"-",h1.Integral(0,-1)))

hIrredBkg2e2mu= ROOT.TH1F() 
hIrredBkg4mu = ROOT.TH1F() 
hIrredBkg4e = ROOT.TH1F() 

hIrredSum = [{"hist":hIrredBkg2e2mu,"name":'2e2m'},{"hist":hIrredBkg4e,"name":'4e'},{"hist":hIrredBkg4mu,"name":'4m'}]         

filesbkg ={}

for b in BkgSamples:
    filesbkg[b["sample"]] = ROOT.TFile(inputdir+b["sample"]+".root") 

print(Red("Irreducible Background\n"))

for h in hIrredSum:
    isFirst=1
    for b in BkgSamples:
        h1 = filesbkg[b["sample"]].Get("ZZTo"+h["name"]+"_Mass_01")
        if h1==None:
            print("For sample ", b["sample"], "h"+h["name"],"has no enetries or is a zombie")      
            continueb
        if isFirst:
            h["hist"]=copy.deepcopy(h1) 
            isFirst=0
            continue
        print("{0} contribution {1}> {2:.2f}\n".format(h["name"],(33-len(h["name"]))*"-",h["hist"].Integral(0,-1)))
        h["hist"].Add(h1)     

print(Red("Reducible Background\n"))


fileFake = ROOT.TFile(inputdir_CR+"data.root")

print(Blue("Central value\n" ))
hRed2e2mu = copy.deepcopy(fileFake.Get("ZZTo2e2m_Mass_01"))
hRed4e    = copy.deepcopy(fileFake.Get("ZZTo4e_Mass_01"))
hRed4mu   = copy.deepcopy(fileFake.Get("ZZTo4m_Mass_01"))

hRedSum = [{"hist":hRed2e2mu,"name":'2e2m'},{"hist":hRed4e,"name":'4e'},{"hist":hRed4mu,"name":'4m'}]

for h in hRedSum:
    Err=ROOT.Double(0.)
    print("{0} contribution {1}> {2:.2f} +- {3:.2f}\n".format(h["name"],(33-len(h["name"]))*"-",h["hist"].IntegralAndError(0,-1,Err),Err))

AccFile     = ROOT.TFile("./Acceptance/Acceptance_"+Set+"_Mass.root")
AccFile_Hi  = ROOT.TFile("./Acceptance/AcceptanceSFactorSqPlus_"+Set+"_Mass.root")
AccFile_Low = ROOT.TFile("./Acceptance/AcceptanceSFactorSqMinus_"+Set+"_Mass.root")

wspace = ROOT.RooWorkspace("w")

for hData,hRed,hIrr in zip(hSum,hRedSum,hIrredSum):

    hData["hist"].Rebin(8)
    hRed["hist"].Rebin(8)

    print("\nnobs_" +hData["name"],hData["hist"].Integral(0,-1))
    wspace.factory("nobs_"+hData["name"]+"["+str(hData["hist"].Integral(0,-1))+",0,300]")

    Err=ROOT.Double(0.)
    print("red_" +hData["name"],hRed["hist"].IntegralAndError(0,-1,Err),"+-",Err)
    wspace.factory("bkg_"+ hData["name"]+"["+str(hRed["hist"].Integral(0,-1))+"]")
    wspace.factory("bkg0_"+ hData["name"]+"["+str(hRed["hist"].Integral(0,-1))+",0,50]")
    wspace.factory("sigmaBkg_"+hData["name"]+"["+str(Err)+"]") 
 

    print("Irr_" +hData["name"],hIrr["hist"].IntegralAndError(0,-1,Err),"+-",Err)
    wspace.factory("irrBkg_"+ hData["name"]+"["+str(hIrr["hist"].Integral(0,-1))+"]")
    wspace.factory("irrBkg0_"+ hData["name"]+"["+str(hIrr["hist"].Integral(0,-1))+",0,50]")
    wspace.factory("sigmaIrrBkg_"+hData["name"]+"["+str(Err)+"]") 


    Acc     = AccFile.Get("TotAcc"+hData["name"]+"_fr").GetVal()
    Acc_Hi  = AccFile_Hi.Get("TotAcc"+hData["name"]+"_fr").GetVal()
    Acc_Low = AccFile_Low.Get("TotAcc"+hData["name"]+"_fr").GetVal()

    print("Acc",Acc,"Sigma Acc",Acc_Hi-Acc)
    wspace.factory("eff_"+hData["name"]+"["+str(Acc)+"]")
    wspace.factory("eff0_"+hData["name"]+"["+str(Acc)+",0,1.]")
    wspace.factory("sigmaEff_"+hData["name"]+"["+str(Acc_Hi-Acc)+"]")

    # Use diffent funcion to for nuisance parameter?
    # wspace.factory("sigmaEff"+hData["name"]+"["+str(Acc-Acc_Low))+"]")
    
#xs_tight = {"2e2m":16.37,"4m":8.19,"4e":8.19,"4l":32.75} #To be corrected with MCFM values

wspace.factory("xs_4m["+str(xs_tight["4m"])+",0,30]")
wspace.factory("xs_4e["+str(xs_tight["4e"])+",0,30]")
wspace.factory("xs_2e2m["+str(xs_tight["2e2m"])+",0,30]")

Lumi = Lumi/1000.
wspace.factory("Lumi["+str(Lumi)+",0,100]")

cross  = LikeCross(wspace)
#print "CROSS",cross
    


