#! /usr/bin/env python

##################################
## G. Pinna (UNITO) - Jun 2015 ##
##################################

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

Lumi = 2.568


if Analysis!="ZZ":
    inputdir     = "./results/ZZRecoAnalyzer_SR_"+Analysis+"/"
    inputdir_CR  = "./results/ZZRecoAnalyzer_CR_"+Analysis+"/"
else:
    inputdir     = "./results/ZZRecoAnalyzer_SR/"
    inputdir_CR  = "./results/ZZRecoAnalyzer_CR/"


def getRedBkg(FinState,Sign):

    fileFake = ROOT.TFile(inputdir_CR+"data.root")
    hFakeRate=ROOT.TH1F()
    if  Sign==1:
        hFakeRate = fileFake.Get("ZZTo"+FinState+"_Mass_01")
        hFakeRate.Add(setErrorsEntries(fileFake.Get("ZZTo"+FinState+"_Mass_FRVarLow")),-1) 
    if  Sign==-1:
        hFakeRate = fileFake.Get("ZZTo"+FinState+"_Mass_01")
        hFakeRate.Add(setErrorsEntries(fileFake.Get("ZZTo"+FinState+"_Mass_FRVarHigh")),1)
    else:
        hFakeRate = fileFake.Get("ZZTo"+FinState+"_Mass_01")
    Err=ROOT.Double(0.)
    Integr= hFakeRate.IntegralAndError(0,-1,Err)
    if hFakeRate.Integral(0,-1)<=0: 
        hFakeRate.Scale(0)    

    print Sign,"Total integral {0} contribution {1}> {2:.2f}\n\n".format(FinState,(33-len(FinState))*"-",hFakeRate.Integral(0,-1))
    return  copy.deepcopy(hFakeRate)

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

for h in hSum:
    h1 = DataFile.Get("ZZTo"+h["name"]+"_Mass_01")
    h["hist"] = copy.deepcopy(h1)
    print "\n",h["name"],"Total integral contribution ----------> ",h1.Integral(0,-1)        
    

hIrredBkg2e2mu= ROOT.TH1F() 
hIrredBkg4mu = ROOT.TH1F() 
hIrredBkg4e = ROOT.TH1F() 

hIrredSum = [{"hist":hIrredBkg2e2mu,"name":'2e2m'},{"hist":hIrredBkg4e,"name":'4e'},{"hist":hIrredBkg4mu,"name":'4m'}]         


for b in BkgSamples:
    filesbkg[b["sample"]] = ROOT.TFile(inputdir+b["sample"]+".root") 
    
    for h in hIrredSum:
        isFirst=1
        for b in BkgSamples:
            h1 = filesbkg[b["sample"]].Get("ZZTo"+h["name"]+"_"+TypeString+"_01")
            if h1==None:
                print "For sample ", b["sample"], "h"+h["name"],"has no enetries or is a zombie"      
                continue
            if isFirst:
                h["hist"]=copy.deepcopy(h1) 
                isFirst=0
                continue
            print "\nTotal integral {0} contribution {1}> {2:.2f}\n\n".format(h["name"],(33-len(h["name"]))*"-",h["hist"].Integral(0,-1))
            h["hist"].Add(h1)     


hRed2e2mu=getRedBkg("2e2m",0)
hRed4e=getRedBkg("4e",0)
hRed4mu=getRedBkg("4m",0)

hRedSum = [{"hist":hRed2e2mu,"name":'2e2m'},{"hist":hRed4e,"name":'4e'},{"hist":hRed4mu,"name":'4m'}]


hRed2e2mu_Hi = getRedBkg("2e2m",1)
hRed4e_Hi    = getRedBkg("4e",1)
hRed4mu_Hi   = getRedBkg("4m",1)

hRedSum_Hi   = [{"hist":hRed2e2mu_Hi,"name":'2e2m'},{"hist":hRed4e_Hi,"name":'4e'},{"hist":hRed4mu_Hi,"name":'4m'}]


hRed2e2mu_Low = getRedBkg("2e2m",-1)
hRed4e_Low    = getRedBkg("4e",-1)
hRed4mu_Low   = getRedBkg("4m",-1)

hRedSum_Low   = [{"hist":hRed2e2mu_Low,"name":'2e2m'},{"hist":hRed4e_Low,"name":'4e'},{"hist":hRed4mu_Low,"name":'4m'}]


AccFile = ROOT.TFile("./Acceptance/Acceptance_"+Set+".root")
AccFile_Hi = ROOT.TFile("./Acceptance/AcceptanceSFactorSqPlus_"+Set+".root")
AccFile_Low = ROOT.TFile("./Acceptance/AcceptanceSFactorSqMinus_"+Set+".root")


#Acc = AccFile.Get("TotAcc"+i["name"]+"_Tot").GetVal() #FIXME

wspace = ROOT.RooWorkspace("w")


for hData,hRed,hRed_Hi,hRed_Low,hIrr in zip(hSum,hRedSum,hRedSum_Hi,hRedSum_Low,hIrredSum):

    ErrIrr=ROOT.Double(0.)

    print "nobs_" +hData["name"],hData["hist"].Integral(0,-1)
    wspace.factory("nobs_"+hData["name"]+"["+str(hData["hist"].Integral(0,-1))+",0,50]")

    wspace.factory("bkg_red_"+ hData["name"]+"["+str(hRed["hist"].Integral(0,-1))+"]")
    wspace.factory("bkg0_red_"+ hData["name"]+"["+str(hRed["hist"].Integral(0,-1))+",0,10]")
    wspace.factory("sigmaBkg_Red_"+hData["name"]+"["+str(-hRed["hist"].Integral(0,-1)+hRed_Low["hist"].Integral(0,-1))+"]")

    # Use diffent funcion to for nuisance parameter?
    # wspace.factory("sigmaBkg_Red_"+hData["name"]+"["+str(hRed_Hi["hist"].Integral(0,-1)-hRed["hist"].Integral(0,-1))+"]") 

    wspace.factory("bkg_irr_"+ hData["name"]+"["+str(hIrr["hist"].Integral(0,-1))+"]")
    wspace.factory("bkg0_irr_"+ hData["name"]+"["+str(-hIrr["hist"].IntegralAndError(0,-1,ErrIrr))+",0,10]")
    wspace.factory("sigmaBkg_Irr_"+hData["name"]+"["+str(ErrIrr)+"]")
    
    print "red_" +hData["name"],hRed["hist"].Integral(0,-1)
    print "redHi_" +hData["name"],hRed_Hi["hist"].Integral(0,-1)
    print "redLow_" +hData["name"],hRed_Low["hist"].Integral(0,-1)
 
    print "Sigma redHi_" +hData["name"],-hRed_Hi["hist"].Integral(0,-1)+hRed["hist"].Integral(0,-1)
    print "Sigma redLow_" +hData["name"],hRed_Low["hist"].Integral(0,-1)-hRed["hist"].Integral(0,-1)

    Acc = AccFile.Get("TotAcc"+hData["name"]+"_fr").GetVal()
    Acc_Hi = AccFile_Hi.Get("TotAcc"+hData["name"]+"_fr").GetVal()
    Acc_Low = AccFile_Low.Get("TotAcc"+hData["name"]+"_fr").GetVal()

    wspace.factory("eff_"+hData["name"]+"["+str(Acc)+"]")
    wspace.factory("eff0_"+hData["name"]+"["+str(Acc)+",0,1.]")
    wspace.factory("sigmaEff_"+hData["name"]+"["+str(Acc_Hi-Acc)+"]")
    # Use diffent funcion to for nuisance parameter?
    # wspace.factory("sigmaEff"+hData["name"]+"["+str(Acc-Acc_Low))+"]")
    

wspace.factory("xs_4m["+str(xs_4m)+",0,30]")
wspace.factory("xs_4e["+str(xs_4e)+",0,30]")
wspace.factory("xs_2e2m["+str(xs_2e2m)+",0,30]")

wspace.factory("Lumi["+str(Lumi)+",0,100]")

cross  = LikeCross(wspace)
print "CROSS",cross
    


