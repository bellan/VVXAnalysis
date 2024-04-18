#! /usr/bin/env python

##################################
## G. Pinna (UNITO) - Jun 2015 ###
##################################

from __future__ import print_function
import ROOT,copy
from ROOT import gSystem, TCanvas, TH1,  TPad, gStyle, TLegend,TGraphAsymmErrors,Math,TArrayD
import CrossInfo
from CrossInfo import* 
import collections
from optparse import OptionParser
import sys
import math
import json
import operator
from Colours import *
import os.path
from  PersonalInfo import*
#from operator import add #del

##################################################################################################################

########################################### Add inclusive systematic #############################################

##################################################################################################################


def addGlobSyst2(hSum,hSumUp,hSumDown,doFiducial):
  
    SystList = GlobSystList

    sigmaTrig = {"4e":[],"4m":[],"2e2m":[],"4l":[]}
    for h,hUp,hDown in zip(hSum,hSumUp,hSumDown):                                                     
        Nbins = h["state"].GetNbinsX()
        print(h["name"])
        for i in range(1,Nbins+1):
            valUp = (hUp["state"].GetBinContent(i)-h["state"].GetBinContent(i))**2
            valDown = (h["state"].GetBinContent(i)-hDown["state"].GetBinContent(i))**2
            for syst in SystList:
#                print doFiducial, h["name"],syst["name"]"Trig": 
                if doFiducial and h["name"] == '4l' and syst["corrEll"]: 
            #for the tight region the trigger uncertainty is added when the 4l cross-section is computed
                    valUp+= 0.
                    valDown+= 0.
                else:
                    print(syst["name"],h["state"].GetBinContent(i)*syst["value"])
                    valUp+=(h["state"].GetBinContent(i)*syst["value"])**2
                    valDown+=(h["state"].GetBinContent(i)*syst["value"])**2
                    if syst["name"]=="Trig": sigmaTrig[h["name"]].append(h["state"].GetBinContent(i)*syst["value"])

            hUp["state"].SetBinContent(i,h["state"].GetBinContent(i)+math.sqrt(valUp))
            hDown["state"].SetBinContent(i,h["state"].GetBinContent(i)-math.sqrt(valDown))

    for sig4e,sig4m,sig2e2m in zip(sigmaTrig["4e"],sigmaTrig["4m"],sigmaTrig["2e2m"]):
        sigmaTrig["4l"].append(math.sqrt(sig4e*sig4e +sig4m*sig4m) +sig2e2m)
        #sigmaTrig["4l"].append(sig4e+sig4m+sig2e2m)


    for h,hup,hdn in zip(hSum,hSumUp,hSumDown):
        Nbins = h["state"].GetNbinsX()
        for i in range(1,Nbins+1):
            if hup["name"] == "4l":
                hup["state"].SetBinContent(i,h["state"].GetBinContent(i)+math.sqrt( pow(hup["state"].GetBinContent(i)-h["state"].GetBinContent(i),2)+pow(sigmaTrig["4l"][i-1],2) ) )
                hdn["state"].SetBinContent(i,h["state"].GetBinContent(i)-math.sqrt( pow(h["state"].GetBinContent(i)-hdn["state"].GetBinContent(i),2)+pow(sigmaTrig["4l"][i-1],2) ) )




def addGlobSyst(hCent,hUp,hDown,FinState,doFiducial):
  
    SystList = GlobSystList
    Nbins = hCent.GetNbinsX()
    for i in range(1,Nbins+1):
        valUp = (hUp.GetBinContent(i)-hCent.GetBinContent(i))**2
        valDown = (hCent.GetBinContent(i)-hDown.GetBinContent(i))**2
        for syst in SystList:
#for the tight region the trigger uncertainty is added when the 4l cross-section is computed
             if doFiducial and FinState == '4l' and syst["name"] == "Trig": 
                 valUp+= 0.
                 valDown+= 0.
             else:
                 print(syst["name"],hCent.GetBinContent(i)*syst["value"])
                 valUp+=(hCent.GetBinContent(i)*syst["value"])**2
                 valDown+=(hCent.GetBinContent(i)*syst["value"])**2
                 
        hUp.SetBinContent(i,hCent.GetBinContent(i)+math.sqrt(valUp))
        hDown.SetBinContent(i,hCent.GetBinContent(i)-math.sqrt(valDown))


##################################################################################################################

##################### Create TGraphAsymmetricError from up and down systematic distributions #####################

##################################################################################################################


def getUncGraph(HCent,HUp,HDown,MCSet,FinState,Data,Err,Type):
    
    NBins = HCent.GetNbinsX()
    
    #Add global systematic uncertainties: acceptance, luminosity and (if not tight region and not 4l) trigger  

    grErr = ROOT.TGraphAsymmErrors(HCent)
    NBins = HCent.GetNbinsX()

    for i in range(1,NBins+1):
            
        statUnc = 0     
        sistUnc_up = 0 
        sistUnc_down = 0 
        totUnc_up = 0 
        totUnc_down = 0
        if Data == True:
            statUnc_up   = HCent.GetBinErrorUp(i)
            statUnc_down = HCent.GetBinErrorLow(i)
        else:
            statUnc_up   = HCent.GetBinError(i)
            statUnc_down = HCent.GetBinError(i)

        sistUnc_up   = HUp.GetBinContent(i)   - HCent.GetBinContent(i)
        sistUnc_down = HCent.GetBinContent(i) - HDown.GetBinContent(i)

        if Err == "syst" or Err == "ratio":
            totUnc_up = sistUnc_up
            totUnc_down = sistUnc_down 
        elif Err == "statsyst":
            totUnc_up   = math.sqrt(statUnc_up*statUnc_up + sistUnc_up*sistUnc_up)
            totUnc_down = math.sqrt(statUnc_down*statUnc_down + sistUnc_down*sistUnc_down)
        
        # print i-1,"cent",HCent.GetBinContent(i)
        # print "stat up",statUnc_up,"syst up",sistUnc_up,"tot up",totUnc_up
        # print "stat down",statUnc_down,"syst down",sistUnc_down,"tot down",totUnc_down

        grErr.SetPointEYhigh(i-1,totUnc_up) 
        grErr.SetPointEYlow(i-1,totUnc_down)

        if "nJets del"in Type and Data and Err=="ratio":
            print("set 0")
            grErr.SetPointEXhigh(i-1,0.) 
            grErr.SetPointEXlow(i-1,0.)
       
    return copy.deepcopy(grErr)

##################################################################################################################

##################### Create TGraphAsymmetricError from up and down systematic distributions #####################

##################################################################################################################


def getCrossPlot_MC(MCSet,Type,analysis,doNormalize,DoFiducial):

    if DoFiducial:    fInMC  = ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/MC_"+Type+"_fr.root")  
    else:             fInMC  = ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/MC_"+Type+".root")  
    
    if Type == "Total": Type = "Mass"

    print(Red("######################### Monte Carlo #######################\n"))

    list2e2mu = []
    list4mu   = []
    list4e    = []
    list4l    = []

    listUpDn = [] # for up and down list. Just to put something. 

    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()
    hsum4l    = ROOT.TH1F()

    hSum = [{"state":hsum2e2mu,"name":'2e2m',"samples":list2e2mu},{"state":hsum4e,"name":'4e',"samples":list4e},{"state":hsum4mu,"name":'4m',"samples":list4mu}]
    if DoFiducial:
        hSum.append({"state":hsum4l,"name":'4l',"samples":list4l})
    
    hsum2e2muUp = ROOT.TH1F()
    hsum4eUp    = ROOT.TH1F()
    hsum4muUp   = ROOT.TH1F()

    hsum4lUp    = ROOT.TH1F()

    hSumUp = [{"state":hsum2e2muUp,"name":'2e2m',"samples":listUpDn},{"state":hsum4e,"name":'4e',"samples":listUpDn},{"state":hsum4muUp,"name":'4m',"samples":listUpDn}]    

    if DoFiducial:
        hSumUp.append({"state":hsum4lUp,"name":'4l',"samples":listUpDn})

    hsum2e2muDn = ROOT.TH1F()
    hsum4eDn    = ROOT.TH1F()
    hsum4muDn   = ROOT.TH1F()
    hsum4lDn    = ROOT.TH1F()

    hSumDn = [{"state":hsum2e2muDn,"name":'2e2m',"samples":listUpDn},{"state":hsum4eDn,"name":'4e',"samples":listUpDn},{"state":hsum4muDn,"name":'4m',"samples":listUpDn}]    

    if DoFiducial:
        hSumDn.append({"state":hsum4lDn,"name":'4l',"samples":listUpDn})

    normStr = ""
    if doNormalize: normStr = "_norm"
    for h,hup,hdn in  zip(hSum,hSumUp,hSumDn):

        if DoFiducial:
            
            h["state"]    = copy.deepcopy(fInMC.Get("ZZTo"+h["name"]+"_"+Type+"Gen_01"+"_fr"+normStr))
            hup["state"]  = copy.deepcopy(fInMC.Get("ZZTo"+h["name"]+"_"+Type+"Gen_01"+"_fr_Up"+normStr))
            hdn["state"]  = copy.deepcopy(fInMC.Get("ZZTo"+h["name"]+"_"+Type+"Gen_01"+"_fr_Dn"+normStr))

            if h["state"]  == None:  sys.exit("ERROR: Check MC histograms. ZZTo"+h["name"]+"_"+Type+"Gen_01"+"_fr"+normStr+" does not exist in "+fInMC.GetName())
            if hup["state"]== None:  sys.exit("ERROR: Check MC histograms. ZZTo"+h["name"]+"_"+Type+"Gen_01_scaleUp"+"_fr"+normStr+" does not exist  in file ./FinalResults_"+MCSet+"_"+analysis+"/MC_"+Type+"_fr.root")
            if hdn["state"]== None: sys.exit("ERROR: Check MC histograms. ZZTo"+h["name"]+"_"+Type+"Gen_01_scaleDn"+"_fr"+normStr+"does not exist") 
        else:  
            h["state"]    = copy.deepcopy(fInMC.Get("ZZTo"+h["name"]+"_"+Type+"Gen_01"+normStr))
            hup["state"]  = copy.deepcopy(fInMC.Get("ZZTo"+h["name"]+"_"+Type+"Gen_01_Up"+normStr))
            hdn["state"]  = copy.deepcopy(fInMC.Get("ZZTo"+h["name"]+"_"+Type+"Gen_01_Dn"+normStr))
            
            if  h["state"] is None: sys.exit("ZZTo"+h["name"]+"_"+Type+"Gen_01"+normStr,"Doesn't exist")
            elif  hup["state"] == None: sys.exit("ZZTo"+h["name"]+"_"+Type+"Gen_01_Up"+normStr+" Doesn't exist")
            
        if MCSet == "Mad":SignalSamples = SignalSamples_Mad
        else: SignalSamples = SignalSamples_Pow 

        for i in SignalSamples:
            h_s   = copy.deepcopy(fInMC.Get(i["sample"]+"_"+h["name"]))
            h_sup = copy.deepcopy(fInMC.Get(i["sample"]+"_"+h["name"]))
            h_sdn = copy.deepcopy(fInMC.Get(i["sample"]+"_"+h["name"]))
    
            if h_s==None: continue

            h["samples"].append(h_s)
            hup["samples"].append(h_sup)
            hdn["samples"].append(h_sdn)

        if h==None:
            print("ERROR no data for",h["name"])
            break

    if not DoFiducial:

        hTOTCross   = copy.deepcopy(hSum[0]["state"])
        hTOTCrossUp = copy.deepcopy(hSumUp[0]["state"])
        hTOTCrossDn = copy.deepcopy(hSumDn[0]["state"])

        hTOTElem   = {"state":hTOTCross,  "color":ROOT.kAzure-6,"name":'4l'}
        hTOTElemUp = {"state":hTOTCrossUp,"color":ROOT.kAzure-6,"name":'4l'}
        hTOTElemDn = {"state":hTOTCrossDn,"color":ROOT.kAzure-6,"name":'4l'}

        hSum.append(hTOTElem)
        hSumUp.append(hTOTElemUp)
        hSumDn.append(hTOTElemDn)

    for h,hup,hdn in zip(hSum,hSumUp,hSumDn):
        print(Red(h["name"]))
        if h["name"]!="4l":
            for (h_s,h_sup,h_sdn) in zip(h["samples"],hup["samples"],hdn["samples"]):
                setCrossSectionMC(h_s,h_sup,h_sdn,h["name"],Type,doNormalize,MCSet,DoFiducial)
        setCrossSectionMC(h["state"],hup["state"],hdn["state"],h["name"],Type,doNormalize,MCSet,DoFiducial) 

    return (hSum,hSumUp,hSumDn)


##############################################################################################################
def getCrossPlot_Data(MCSet,UseUnfold,Type,analysis,Sign,UseMCReco,doNormalize,doFiducial):
    
    Fr = ""
    if doFiducial: Fr="_fr"
    
    if UseMCReco:  print(Red("########################### MC RECO ########################\n"))
    else: print(Red("############################ DATA  #########################\n".format(Sign)))
    
    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()
    hsum4l   = ROOT.TH1F()
    
    if not doFiducial:  hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]    
    else:               hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'},{"state":hsum4l,"name":'4l'}]    
    
    hsum2e2muUp = ROOT.TH1F()
    hsum4eUp    = ROOT.TH1F()
    hsum4muUp   = ROOT.TH1F()
    hsum4lUp   = ROOT.TH1F()

    if not doFiducial: hSumUp = [{"state":hsum2e2muUp,"name":'2e2m'},{"state":hsum4eUp,"name":'4e'},{"state":hsum4muUp,"name":'4m'}]    
    else: hSumUp = [{"state":hsum2e2muUp,"name":'2e2m'},{"state":hsum4eUp,"name":'4e'},{"state":hsum4muUp,"name":'4m'},{"state":hsum4lUp,"name":'4l'}]    

    hsum2e2muDown = ROOT.TH1F()
    hsum4eDown    = ROOT.TH1F()
    hsum4muDown   = ROOT.TH1F()
    hsum4lDown   = ROOT.TH1F()

    if not doFiducial:    hSumDown = [{"state":hsum2e2muDown,"name":'2e2m'},{"state":hsum4eDown,"name":'4e'},{"state":hsum4muDown,"name":'4m'}]
    else:   hSumDown = [{"state":hsum2e2muDown,"name":'2e2m'},{"state":hsum4eDown,"name":'4e'},{"state":hsum4muDown,"name":'4m'},{"state":hsum4lDown,"name":'4l'}]    
  
    if UseMCReco:
        FInCenter = ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/MCReco.root")
        FInUp     = ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/DataUp.root") 
        FInDown   = ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/DataDown.root")  
    else:
        if UseUnfold:
            FInCenter =  ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/DataUnfold"+Fr+".root")
            FInUp     =  ROOT.TFile( "./FinalResults_"+MCSet+"_"+analysis+"/DataUnfoldUp"+Fr+".root")
            FInDown   =  ROOT.TFile( "./FinalResults_"+MCSet+"_"+analysis+"/DataUnfoldDown"+Fr+".root")  
        else: 
            FInCenter = ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/Data"+Fr+".root") 
            FInUp     = ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/DataUp"+Fr+".root")
            FInDown   = ROOT.TFile( "./FinalResults_"+MCSet+"_"+analysis+"/DataDown"+Fr+".root")  

    if doNormalize: normStr = "_norm"
    else: normStr = ""

    for h,hup,hdown in zip(hSum,hSumUp,hSumDown):
        if not doFiducial and  h["name"]=="4l": continue

        h["state"]     = copy.deepcopy(FInCenter.Get("ZZTo"+h["name"]+"_"+Type+"_01"))
        hup["state"]   = copy.deepcopy(FInUp.Get("ZZTo"+h["name"]+"_"+Type+"_01"+normStr))
        hdown["state"] = copy.deepcopy(FInDown.Get("ZZTo"+h["name"]+"_"+Type+"_01"+normStr))

        if h["state"]==None:   sys.exit( "ERROR no data for "+h["name"])

        setCrossSectionData(h["state"],h["name"],MCSet,doFiducial)
        setCrossSectionData(hup["state"],hup["name"],MCSet,doFiducial)
        setCrossSectionData(hdown["state"],hdown["name"],MCSet,doFiducial)


    if not doFiducial:
        hTOTCross= ROOT.TH1F()
        (hTOTCross,hTOTCrossUp,hTOTCrossDown)=combineCrossCorr(hSum,hSumUp,hSumDown) 

        hSum.append({"state":hTOTCross,"name":'4l'})
        hSumUp.append({"state":hTOTCrossUp,"name":'4l'})
        hSumDown.append({"state":hTOTCrossDown,"name":'4l'})

    addGlobSyst2(hSum,hSumUp,hSumDown,doFiducial)
    for h,hup,hdown in zip(hSum,hSumUp,hSumDown):

        print("{0} Tot Cross     {1} {2:.2f} \n".format(h["name"],(15-len(h["name"]))*" ", h["state"].Integral(0,-1)))     
        print("{0} Tot Cross up  {1} {2:.2f} \n".format(hup["name"],(15-len(hup["name"]))*" ", hup["state"].Integral(0,-1)))     
        print("{0} Tot Cross down{1} {2:.2f} \n".format(hdown["name"],(15-len(hdown["name"]))*" ", hdown["state"].Integral(0,-1)))     

        if h["state"]==None:
            print(h["state"]," has no enetries")

            #      NormUp   =   hup["state"].Integral(1,-1)/h["state"].Integral(1,-1) 
            #      NormDown = hdown["state"].Integral(1,-1)/h["state"].Integral(1,-1) 

        DoInclNormalized= False
        DoNormal = True

        if doNormalize:
            hup["state"].Scale(1./h["state"].Integral(),"width")
            hdown["state"].Scale(1./h["state"].Integral(),"width")
            h["state"].Scale(1./h["state"].Integral(),"width")

        elif DoNormal:
            # hup["state"].Scale((1000.)/(Lumi),"width")
            # hdown["state"].Scale((1000.)/(Lumi),"width")
            # h["state"].Scale((1000.)/(Lumi),"width")
            hup["state"].Scale(1.,"width")
            hdown["state"].Scale(1.,"width")
            h["state"].Scale(1.,"width")
                                                
        elif DoInclNormalized:
            norm = GetNorm(h["name"],Type,doFiducial)
            hup["state"].Scale(1./norm,"width")
            hdown["state"].Scale(1./norm,"width")
            h["state"].Scale(1./norm,"width")

        else:
            norm = GetNorm(h["name"],Type,doFiducial)
            hup["state"].Scale((norm)/h["state"].Integral(0,-1),"width")
            hdown["state"].Scale((norm)/h["state"].Integral(0,-1),"width")
            h["state"].Scale(norm/h["state"].Integral(0,-1),"width")
            print("NORM",norm)


    return hSum,hSumUp,hSumDown


##########################################################
################# Set Cross Section ######################

def setCrossSectionData(h1,FinState,MCSet,doFiducial):

    if   FinState=='4e':   BR=BRele*BRele
    elif FinState=='4m':   BR=BRmu*BRmu
    elif FinState=='2e2m': BR=2*BRmu*BRele

    if doFiducial: 
        # norm = GetNorm(FinState,Type,doFiducial) 
         #h1.Scale(norm/Integral,"width")  
         h1.Scale(1000./(Lumi)) #if you don't want to use official xs

    else:   h1.Scale(1./(Lumi*BR))
    
##########################################################
def setCrossSectionMC(h1,h1up,h1dn,FinState,Type,doNormalize,MCSet,doFiducial):
    
    if doFiducial: BR = 1
    else:
        if   FinState=='4e':   BR=BRele*BRele
        elif FinState=='4m':   BR=BRmu*BRmu
        elif FinState=='2e2m': BR=2*BRmu*BRele
        else: BR=2*BRmu*BRele

    name = (h1.GetName()).replace("_"+FinState,"")
    if "Gen" in name: name = "Total"
    if doFiducial:
        print("{0} {1} {2:.2f} + {3:.2f} -{4:.2f} [fb]\n".format(name,((25-len(name))*" "),1000*(h1.Integral(0,-1))/(Lumi),1000*((h1up.Integral(0,-1) - h1.Integral(0,-1))/(Lumi)),1000*((h1.Integral(0,-1) - h1dn.Integral(0,-1) )/(Lumi))))
        if FinState=="4l" and Type=="Total":
            AccFile = ROOT.TFile("./Acceptance/Acceptance_"+MCSet+"_"+Type+".root")
            Acc = AccFile.Get("TotAcc2e2m_Acc").GetVal()
            print("Cross section Wide region",(h1.Integral(0,-1))/(Lumi*Acc*(BRele*BRele+2*BRele*BRmu+BRmu*BRmu)),"[pb]")              
    else:
        print("{0} {1} {2:.2f} [pb]\n".format(name,((25-len(name))*" "),(h1.Integral(1,-1))/(Lumi*BR))) # Check total cross section witho
    
    Integral = h1.Integral() 

    if doNormalize:  
        h1.Scale(1.,"width") 
        h1up.Scale(1.,"width") 
        h1dn.Scale(1.,"width") 
    elif doFiducial:  
        h1.Scale(1000./(Lumi),"width")  
        h1up.Scale(1000./(Lumi),"width")  
        h1dn.Scale(1000./(Lumi),"width")  
    else:
        h1.Scale(1./(Lumi*BR),"width")  
        h1up.Scale(1./(Lumi*BR),"width")  
        h1dn.Scale(1./(Lumi*BR),"width")  
        
##################################################################################################################

##################### Combine different final state cross section distribution for final 4l ######################

####################### using weighted mean formula with both statistic ad sytematic errors ######################

##################################################################################################################


def combineCross(HList,HListUp,HListDown):

    HCross     = ROOT.TH1F()
    HCrossUp   = ROOT.TH1F()
    HCrossDown = ROOT.TH1F()

    #just equal to one
    HCross=copy.deepcopy(HList[1]["state"])
    HCrossUp=copy.deepcopy(HListUp[1]["state"])
    HCrossDown=copy.deepcopy(HListDown[1]["state"])

    Nbins= HList[1]["state"].GetNbinsX()

    for i in range(1,Nbins+1):
        Hlist  = list(zip(HList, HListUp,HListDown))
        #Sort List by entries magnitude, from higher to lower to skip 0 entries bins
        SortedHlist = sorted(Hlist,key=lambda value: value[0]["state"].GetBinContent(i),reverse = True)        

        Cross       = 0.
        ErrSystUp   = 0.
        ErrSystDown = 0.
        ErrStat     = 0.

        WeightStat     = 0.
        WeightTot      = 0.
        WeightSystUp   = 0.
        WeightSystDown = 0.


        for elem in SortedHlist:

            Entries = elem[0]["state"].GetBinContent(i)
            print(elem[0]["name"],elem[0]["state"].GetBinContent(i))
            if Entries == 0.:   break   # Because of sorting also others final state will be 0 so use break and no continue

            errStatSq     = (elem[0]["state"].GetBinError(i)/Entries)**2
            errSystUpSq   = ((elem[1]["state"].GetBinContent(i)-Entries)/Entries)**2
            errSystDownSq = ((elem[2]["state"].GetBinContent(i)-Entries)/Entries)**2
            errSystSq     = ((elem[1]["state"].GetBinContent(i)-elem[2]["state"].GetBinContent(i) )/(Entries*2.))**2 #Use the average of systematic up and down. is it a +, right? //DELhot

            weightStat =  1./errStatSq

            if errSystUpSq==0:    weightSystUp   = 0.
            else:                 weightSystUp   =  1./errSystUpSq

            if errSystDownSq==0:  weightSystDown = 0.
            else:                 weightSystDown =  1./errSystDownSq

            weightTot       = 1./(errStatSq+errSystSq)
            WeightStat     += weightStat
            WeightTot      += weightTot
            WeightSystUp   += weightSystUp
            WeightSystDown += weightSystDown 

            Cross          += weightTot*Entries

        Cross = Cross/WeightTot

        ErrStat     = math.sqrt(1./WeightStat)        
        ErrSystUp   = math.sqrt(1./WeightSystUp)
        if WeightSystDown ==0:  
            ErrSystDown = 0 
            print("ErrSystDown for bin",i,"is 0. Check it")
        else: ErrSystDown = math.sqrt(1./WeightSystDown)
        
        HCrossUp.SetBinContent(i,Cross+ErrSystUp*Cross)
        HCrossDown.SetBinContent(i,Cross-ErrSystDown*Cross)
        HCross.SetBinContent(i,Cross)
        HCross.SetBinError(i,ErrStat*Cross)      
 
    return HCross,HCrossUp,HCrossDown

def combineCrossCorr(HList,HListUp,HListDown):

    HCross     = ROOT.TH1F()
    HCrossUp   = ROOT.TH1F()
    HCrossDown = ROOT.TH1F()

    #just equal to one
    HCross=copy.deepcopy(HList[1]["state"])
    HCrossUp=copy.deepcopy(HListUp[1]["state"])
    HCrossDown=copy.deepcopy(HListDown[1]["state"])

    Nbins= HList[1]["state"].GetNbinsX()

    for i in range(1,Nbins+1):
        Hlist  = list(zip(HList, HListUp,HListDown))
        #Sort List by entries magnitude, from higher to lower to skip 0 entries bins
        SortedHlist = sorted(Hlist,key=lambda value: value[0]["state"].GetBinContent(i),reverse = True)        

        Cross       = 0.
        ErrSystUp   = 0.
        ErrSystDown = 0.
        ErrStat     = 0.

        WeightStat     = 0.
        WeightTot      = 0.
        WeightSystUp   = 0.
        WeightSystDown = 0.


        errStatSq = {"2e2m":0,"4m":0,"4e":0}
        errSystSq = {"2e2m":0,"4m":0,"4e":0}
        errSyst   = {"2e2m":0,"4m":0,"4e":0}
        errStat   = {"2e2m":0,"4m":0,"4e":0}
        Entries   = {"2e2m":0,"4m":0,"4e":0}

        for elem in SortedHlist:

            Entries[elem[0]["name"]] = elem[0]["state"].GetBinContent(i)
            #print  elem[0]["name"],elem[0]["state"].GetBinContent(i),(elem[0]["state"].GetBinError(i))
            if Entries[elem[0]["name"]] == 0.:   break   # Because of sorting also others final state will be 0 so use break and no continue

#            errStatSq[ elem[0]["name"]]     = (elem[0]["state"].GetBinError(i)/Entries[elem[0]["name"]])**2
            errStat[elem[0]["name"]]     = (elem[0]["state"].GetBinError(i))
            errStatSq[elem[0]["name"]]   = (elem[0]["state"].GetBinError(i))**2
            errSyst[elem[0]["name"]]     = ((elem[1]["state"].GetBinContent(i)-elem[2]["state"].GetBinContent(i) )/(2.)) #Use the average of systematic up and do\

        Cross,ErrTot,ErrStat = ToyMCAverage(Entries["2e2m"],errStat["2e2m"],errSyst["2e2m"],Entries["4m"],errStat["4m"],errSyst["4m"],Entries["4e"],errStat["4e"],errSyst["4e"],i)
        print("2e2m",Entries["2e2m"], errStat["2e2m"], errStat["2e2m"]/Entries["2e2m"],errSyst["2e2m"]/Entries["2e2m"])
        print("4m"  ,Entries["4m"],   errStat["4m"],   errStat["4m"]/Entries["4m"],    errSyst["4m"]/Entries["4m"])
        print("4e",  Entries["4e"],   errStat["4e"],   errStat["4e"]/Entries["4e"],    errSyst["4e"]/Entries["4e"])
        print(Cross,ErrTot)
        
        HCrossUp.SetBinContent(i,Cross+ErrTot)
        HCrossDown.SetBinContent(i,Cross-ErrTot)
        HCross.SetBinContent(i,Cross)
        HCross.SetBinError(i,ErrStat)      
 
    return HCross,HCrossUp,HCrossDown

def ToyMCAverage(val1,stat1,syst1,val2,stat2,syst2,val3,stat3,syst3,bin):

    ctoy  = ROOT.TCanvas("ctoy","ctoy", 200, 10,1200, 900) 
    r =   ROOT.TRandom()

    hToy     = ROOT.TH1F("hToy","",100,val1-stat1*4,val1+stat1*4)
    hToyStat = ROOT.TH1F("hToystat","",100,val1-stat1*4,val1+stat1*4)

    if val1<=0.01 or val2<=0.01 or val3<=0.01:

        if val1<=0.01: 
            val1  = val3
            stat1 = stat3
            syst1 = syst3
        elif val2<=0.01:
            val2  = val3
            stat2 = stat3
            syst2 = syst3

        for i in range(0,10000):
            
            corr = r.Gaus(0,1)                
            hToy.Fill(
                ( 
                    ((val1+r.Gaus(0,1)*stat1+corr*syst1)/(stat1*stat1+syst1*syst1))
                    +((val2+r.Gaus(0,1)*stat2+corr*syst2)/(stat2*stat2+syst2*syst2)))
                     /( 1/(stat1*stat1+syst1*syst1) +1/(stat2*stat2+syst2*syst2)) )
            
            hToyStat.Fill(
                (
                    ((val1+r.Gaus(0,1)*stat1)/(stat1*stat1 ))
                    +((val2+r.Gaus(0,1)*stat2)/(stat2*stat2))
                    )/( 1/(stat1*stat1) +1/(stat2*stat2)) )

    else:
        for i in range(0,100000):
            
            corr = r.Gaus(0,1)      
            hToy.Fill(
                ( 
                    ((val1+r.Gaus(0,1)*stat1+corr*syst1)/((stat1*stat1+syst1*syst1)/(val1*val1)))
                    +((val2+r.Gaus(0,1)*stat2+corr*syst2)/((stat2*stat2+syst2*syst2)/(val2*val2)))
                    +((val3+r.Gaus(0,1)*stat3+corr*syst3)/((stat3*stat3+syst3*syst3)/(val3*val3))))/( (val1*val1)/(stat1*stat1+syst1*syst1) +(val2*val2)/(stat2*stat2+syst2*syst2) +(val3*val3)/(stat3*stat3+syst3*syst3)  ) )
            
            hToyStat.Fill(
                (
                    (((val1+r.Gaus(0,1)*stat1)*val1*val1)/(stat1*stat1 ))
                    +(((val2+r.Gaus(0,1)*stat2)*val2*val2)/(stat2*stat2))
                    +(((val3+r.Gaus(0,1)*stat3)*val3*val3)/(stat3*stat3)) )/( (val1*val1)/(stat1*stat1) +(val2*val2)/(stat2*stat2) +(val3*val3)/(stat3*stat3)  ) )
       
    ctoy.Divide(2,1)

    ctoy.cd(1)
    hToy.SetLineColor(4)
    hToy.Draw("hist")
    f = hToy.Fit("gaus","S")

    ctoy.cd(2)
    hToyStat.Draw("hist")
    fStat = hToyStat.Fit("gaus","S")
    Type ="Mass"

    ctoy.SaveAs(PersonalFolder+"/test/ToyMC_"+Type+"bin_"+str(bin)+".png")

    return f.GetParams()[1],f.GetParams()[2],fStat.GetParams()[2]


def combineCrossFiducial(HList,HListUp,HListDown):

    HCross    = ROOT.TH1F()
    HCrossUp  = ROOT.TH1F()
    HCrossDown= ROOT.TH1F()

    #just copy one finalstate
    HCross     = copy.deepcopy(HList[1]["state"])
    HCrossUp   = copy.deepcopy(HListUp[1]["state"])
    HCrossDown = copy.deepcopy(HListDown[1]["state"])

    HCross.SetName(HCross.GetName().replace("4e", "4l"))
    HCrossUp.SetName(HCrossUp.GetName().replace("4e", "4l"))
    HCrossDown.SetName(HCrossDown.GetName().replace("4e", "4l"))
   
    Nbins= HList[1]["state"].GetNbinsX()

    for i in range(1,Nbins+1):

        Hlist  = list(zip(HList, HListUp,HListDown))
        #Sort List by entries magnitude, from higher to lower  to skip 0 entries bins
        SortedHlist = sorted(Hlist,key=lambda value: value[0]["state"].GetBinContent(i),reverse = True)        

        Cross       = 0
        ErrSystUp   = 0
        ErrSystDown = 0
        ErrStat     = 0

        for elem in SortedHlist:
            Entries = elem[0]["state"].GetBinContent(i)
            if Entries == 0.:  break   # Because of sorting also others final state will be 0 so use break and no continue
            
            #print elem[0]["state"].GetBinContent(i),elem[1]["state"].GetBinContent(i),elem[2]["state"].GetBinContent(i)

            errStatSq   = (elem[0]["state"].GetBinError(i))**2
            errSystUp   = (elem[1]["state"].GetBinContent(i)-Entries)**2
            errSystDown = (elem[2]["state"].GetBinContent(i)-Entries)**2

            Cross       += Entries 
            ErrStat     += errStatSq 
            ErrSystUp   += errSystUp
            ErrSystDown += errSystDown
       
        sigma_4e = Hlist[0][0]["state"].GetBinContent(i)
        sigma_4m = Hlist[1][0]["state"].GetBinContent(i)
        sigma_2e2m = Hlist[2][0]["state"].GetBinContent(i)

        #4e and 4m are not correlated while 2e2m is correlated with both.  
        sigma_4l_corrunc = math.sqrt(sigma_4e*sigma_4e +sigma_4m*sigma_4m) +sigma_2e2m  
        corrunc_4l = (GlobSystList[0]["value"]*sigma_4l_corrunc)**2  #trigger 

        ErrStatSq     = math.sqrt( ErrStat )  
        ErrSystUpSq   = math.sqrt( ErrSystUp   + corrunc_4l ) 
        ErrSystDownSq = math.sqrt( ErrSystDown + corrunc_4l )

        HCrossUp.SetBinContent(i,Cross+ErrSystUpSq)
        HCrossDown.SetBinContent(i,Cross-ErrSystDownSq)
        HCross.SetBinContent(i,Cross)
        HCross.SetBinError(i,ErrStatSq)      
 
    return HCross,HCrossUp,HCrossDown


def getHisto(Sign,MCSet,analysis,Type,doFiducial,UseMCReco):
    
    fIn = ROOT.TFile() 

    doFid = ""
    if doFiducial: doFid = "_fr" 
    
    if Sign==0: 
        if UseMCReco:
            fIn = ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/MCReco.root")
        else:
            fIn = ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/Data"+doFid+".root") 

    elif Sign==1: fIn= ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/DataUp"+doFid+".root") 
    elif Sign==-1: fIn= ROOT.TFile( "./FinalResults_"+MCSet+"_"+analysis+"/DataDown"+doFid+".root")
  
    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()
   
    hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]    

    for h in hSum:
        if (Type == "Mass") or (Type == "Jets"):
            h1 = fIn.Get("ZZTo"+h["name"]+"_"+Type+"_01")
            h["state"] = copy.deepcopy(h1)
        elif Type=="Tot":
            h1 = fIn.Get("ZZTo"+h["name"]+"_"+Type+doFid)
            h["state"] = copy.deepcopy(h1)
            ErrStat=ROOT.Double(0.)
        else:
            print(Red("\n WRONG TYPE: Mass or Jets only\n")) 
      
    return hSum


def TotalCross(MCSet,Type,analysis,doFiducial,UseMCReco):

    hData     = getHisto(0,MCSet,analysis,Type,doFiducial,UseMCReco)
    hDataUp   = getHisto(1,MCSet,analysis,Type,doFiducial,UseMCReco)
    hDataDown = getHisto(-1,MCSet,analysis,Type,doFiducial,UseMCReco)
   
    list2e2m = [0,0,0,0,0]
    list4e   = [0,0,0,0,0]
    list4m   = [0,0,0,0,0]
 
    CrossDic  =  collections.OrderedDict()
    CrossDic["2e2m"] = list2e2m
    CrossDic["4e"]   = list4e
    CrossDic["4m"]   = list4m
    
    for i,j,k in zip(hData,hDataUp,hDataDown):

        if doFiducial: BR = 1
        else:
            if   i["name"]=='4e':   BR=BRele*BRele
            elif i["name"]=='4m':   BR=BRmu*BRmu
            elif i["name"]=='2e2m': BR=2*BRmu*BRele
            if   i["name"]=='4l':   break
            
        ErrStat=ROOT.Double(0.)

        CrossDic[i["name"]][0] = i["state"].IntegralAndError(0,-1,ErrStat)/(BR*Lumi) 
        CrossDic[i["name"]][1] = ErrStat/(BR*Lumi)
        CrossDic[i["name"]][2] = j["state"].Integral(0,-1)/(BR*Lumi) - CrossDic[i["name"]][0]     
        CrossDic[i["name"]][3] = CrossDic[i["name"]][0] - k["state"].Integral(0,-1)/(BR*Lumi)
        
        MidSyst =  (abs(CrossDic[i["name"]][2]) + abs(CrossDic[i["name"]][3]))/2.
        CrossDic[i["name"]][4] = math.sqrt(CrossDic[i["name"]][1]**2+MidSyst**2) 

    print(Red("\n##################### RESULTS FOR INCLUSIVE CROSS SECTION ####################\n")) 

    if doFiducial: list4l = combineInclusiveCrossFiducial(CrossDic)
    else:          list4l = combineInclusiveCross(CrossDic)
    dic4l = {"4l":list4l}  
    CrossDic.update(dic4l)
 
    for fin,value in CrossDic.items():
        value[2]=value[2]**2
        value[3]=value[3]**2 
#        value[4]=value[4]**2 
        if doFiducial and fin == "4l": #only lumi is added #HOT
            value[2]= value[2]+(value[0]*0.026)**2
            value[3]= value[3]+(value[0]*0.026)**2
        else:            
            for syst in GlobSystList: #Add inclusive systematic errors
                value[2]+=(value[0]*syst["value"])**2
                value[3]+=(value[0]*syst["value"])**2 
                #value[4]+=(value[0]*syst["value"])**2 #CHECk
        value[4]= math.sqrt(value[2]+value[3])#CHECK
        #value[4]= math.sqrt(value[4])
        value[2]= math.sqrt(value[2])
        value[3]= math.sqrt(value[3])

        print(Blue("{0}".format(fin)))
        if doFiducial:
            print(" {0:.2f} +- {1:.2f} (stat) + {2:.2f} (syst) - {3:.2f} (syst) +- {4:.2f} (Total) [fb]\n".format(value[0]*1000,value[1]*1000,value[2]*1000,value[3]*1000,value[4]*1000))
            if fin=="4l":
                AccFile = ROOT.TFile("./Acceptance/Acceptance_"+MCSet+"_Mass.root")
                Acc = (AccFile.Get("TotAcc2e2m_Acc").GetVal()+AccFile.Get("TotAcc4e_Acc").GetVal()+AccFile.Get("TotAcc4m_Acc").GetVal())/3. 
                print("Cross Section in wide Region",value[0]/(Acc*(BRele*BRele+2*BRele*BRmu+BRmu*BRmu)),"\nAcceptance",Acc,"\n") 

        else: print(" {0:.2f} +- {1:.2f} (stat) + {2:.2f} (syst) - {3:.2f} (syst) +- {4:.2f} (Total) [pb] \n".format(value[0],value[1],value[2],value[3],value[4]))


def combineInclusiveCross(Dic):

    TotStat     = math.sqrt(1./(1./(Dic["2e2m"][1]*Dic["2e2m"][1])+ 1./(Dic["4e"][1]*Dic["4e"][1])+ 1./(Dic["4m"][1]*Dic["4m"][1])))

    TotSystUp   = math.sqrt(1./(1./(Dic["2e2m"][2]*Dic["2e2m"][2])+ 1./(Dic["4e"][2]*Dic["4e"][2])+ 1./(Dic["4m"][2]*Dic["4m"][2])))

    TotSystDown = math.sqrt(1./(1./(Dic["2e2m"][3]*Dic["2e2m"][3])+ 1./(Dic["4e"][3]*Dic["4e"][3])+ 1./(Dic["4m"][3]*Dic["4m"][3])))
         
    TotErr      = math.sqrt(1./(1./(Dic["2e2m"][4]*Dic["2e2m"][4])+ 1./(Dic["4e"][4]*Dic["4e"][4])+ 1./(Dic["4m"][4]*Dic["4m"][4])))

    WhTot       = 1./(1./(Dic["2e2m"][4]*Dic["2e2m"][4])+1./(Dic["4e"][4]*Dic["4e"][4])+1./(Dic["4m"][4]*Dic["4m"][4]))

    TotCross    = (Dic["2e2m"][0]/(Dic["2e2m"][4]*Dic["2e2m"][4])+Dic["4e"][0]/(Dic["4e"][4]*Dic["4e"][4])+Dic["4m"][0]/(Dic["4m"][4]*Dic["4m"][4]))*WhTot

    return [TotCross,TotStat,TotSystUp,TotSystDown,TotErr]


def combineInclusiveCrossFiducial(Dic):


    TotStat     = math.sqrt((Dic["2e2m"][1]*Dic["2e2m"][1])+(Dic["4e"][1]*Dic["4e"][1])+(Dic["4m"][1]*Dic["4m"][1]))

    TriggerUnc  = (0.015*(math.sqrt(Dic["4e"][0]*Dic["4e"][0]+Dic["4m"][0]*+Dic["4m"][0]) + Dic["2e2m"][0]))**2

    TotSystUp   = math.sqrt((Dic["2e2m"][2])**2+(Dic["4e"][2])**2+(Dic["4m"][2])**2 +TriggerUnc)

    TotSystDown = math.sqrt((Dic["2e2m"][3])**2+(Dic["4e"][3])**2+(Dic["4m"][3])**2 +TriggerUnc)

    TotSyst     = (TotSystUp + TotSystDown)/2.          

    TotErr      = math.sqrt(TotStat*TotStat+TotSyst*TotSyst)

    TotCross    = Dic["2e2m"][0]+Dic["4e"][0]+Dic["4m"][0]

    return [TotCross,TotStat,TotSystUp,TotSystDown,TotErr]

        
def WriteJetsNorm(hMCList,Type,DoFid):

    if "Central" in Type:
        if DoFid:  out_file = open("JetsNorm_Central_tight.json","w")    
        else:      out_file = open("JetsNorm_Central.json","w")    
    else:
        if DoFid:  out_file = open("JetsNorm_tight.json","w")    
        else:      out_file = open("JetsNorm.json","w")    
    out_file.write("{\n")
    for i in range(0,4):
        out_file.write('"'+hMCList[i]["name"]+'":')
        json.dump([hMCList[i]["state"].Integral(2,4),hMCList[i]["state"].Integral(3,4)],out_file)
        if i == 3: out_file.write("\n")
        else:     out_file.write(",\n")
    out_file.write("}")
    out_file.close()

def GetNorm(finState,Type,doFid):

    if "Jets" in Type or "Mass" in Type or "PtZZ" in Type or "dRZZ" in Type:
        if doFid:
            return xs_tight_exp[finState]
        else: return xs_wide 
    else:
        FileName = "JetsNorm"
        CommandName = "./python/ComputeCross.py -t Jets"
        if doFid:
            if "Central" in Type: 
                FileName += "_Central_tight.json"
                CommandName += "_Central -f"
            else:                 
                FileName += "_tight.json"
                CommandName += " -f"
        else:
            if "Central" in Type: 
                FileName += "_Central.json"
                CommandName += "_Central"
            else:                 
                FileName += ".json"

        if not os.path.exists(FileName): 
                sys.exit("ERROR: "+FileName+" doesn't exist yet. You need to run '"+CommandName+"' first to put normalization values to file")

        File = open(FileName,"r") 
        normjson = json.load(File)
        if "Jet1" in Type: return normjson[finState][0]
        elif "Jet2" in Type or "Deta" in Type or "Mjj" in Type: return normjson[finState][1]



def GetHistoForChi2(h,gr):
    
    NBins = h.GetNbinsX()
    hnew = copy.deepcopy(h)        
    for i in range(1,NBins+1):
        #   print gr.GetErrorYhigh(i-1),gr.GetErrorYhigh(i-1),abs(gr.GetErrorYhigh(i-1)+gr.GetErrorYhigh(i-1))/2
        hnew.SetBinError(i,abs(gr.GetErrorYhigh(i-1)+gr.GetErrorYhigh(i-1))/2)
    return hnew
