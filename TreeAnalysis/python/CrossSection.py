#! /usr/bin/env python

##################################
## G. Pinna (UNITO) - Jun 2015 ###
##################################

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

##################################################################################################################

########################################### Add inclusive systematic #############################################

##################################################################################################################

def addGlobSist(hCent,hUp,hDown):
      
    Nbins = hCent.GetNbinsX()
    for i in range(1,Nbins+1):
        valUp = (hUp.GetBinContent(i)-hCent.GetBinContent(i))**2
        valDown = (hCent.GetBinContent(i)-hDown.GetBinContent(i))**2
        #print "i",i, math.sqrt(valUp)      
        for sist in GlobSistList:
            valUp+=(hCent.GetBinContent(i)*sist["value"])**2
            valDown+=(hCent.GetBinContent(i)*sist["value"])**2

        #print i, math.sqrt(valUp)
        hUp.SetBinContent(i,hCent.GetBinContent(i)+math.sqrt(valUp))
        #print hCent.GetBinContent(i)+math.sqrt(valUp),"\n"
        hDown.SetBinContent(i,hCent.GetBinContent(i)-math.sqrt(valDown))

#######################for the Fiducial region
def addGlobSist4lFiducial(hCent,hUp,hDown):
      
    Nbins = hCent.GetNbinsX()
    for i in range(1,Nbins+1):
        valGlobal = (hCent.GetBinContent(i)*(0.026))**2 #lumi
        valUp = (hUp.GetBinContent(i)-hCent.GetBinContent(i))**2 + valGlobal
        valDown = (hCent.GetBinContent(i)-hDown.GetBinContent(i))**2 + valGlobal
        #print "i",i, math.sqrt(valUp)      
      
        #print i, math.sqrt(valUp)
        hUp.SetBinContent(i,hCent.GetBinContent(i)+math.sqrt(valUp))
        #print hCent.GetBinContent(i)+math.sqrt(valUp),"\n"
        hDown.SetBinContent(i,hCent.GetBinContent(i)-math.sqrt(valDown))

##################################################################################################################

##################### Create TGraphAsymmetricError from up and down systematic distributions #####################

##################################################################################################################


def getSistGraph(HCent,HUp,HDown,DoFiducial,FinState):

    if DoFiducial and FinState == '4l': addGlobSist4lFiducial(HCent,HUp,HDown)
    else: addGlobSist(HCent,HUp,HDown)
    grSist = ROOT.TGraphAsymmErrors(HCent)
    NBins = HCent.GetNbinsX()
    for i in range(1,NBins+1):
#        print i,HCent.GetBinContent(i),"+-", HCent.GetBinError(i),"+",HUp.GetBinContent(i)-HCent.GetBinContent(i),"-",HCent.GetBinContent(i)- HDown.GetBinContent(i)
        grSist.SetPointEYhigh(i-1,HUp.GetBinContent(i)-HCent.GetBinContent(i))
        grSist.SetPointEYlow(i-1,HCent.GetBinContent(i)- HDown.GetBinContent(i))

    return copy.deepcopy(grSist)

##################################################################################################################

def getCrossPlot_MC(MCSet,Type,analysis,DoNormalized,DoFiducial):

    fInMC  = ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/MC.root")  
   
    print Red("######################### Monte Carlo #######################\n")

    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()

    hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]    

    for h in hSum:
        if DoFiducial:
            h["state"]    = copy.deepcopy(fInMC.Get("ZZTo"+h["name"]+"_"+Type+"Gen_01_fr"))
        else:  h["state"] = copy.deepcopy(fInMC.Get("ZZTo"+h["name"]+"_"+Type+"Gen_01"))

        if h==None:
            print "ERROR no data for",h["name"]
            break

    hTOTCross = copy.deepcopy(hSum[0]["state"])
    if DoFiducial:
        hTOTCross.Add(hSum[1]["state"])
        hTOTCross.Add(hSum[2]["state"])

    hTOTElem = {"state":hTOTCross,"color":ROOT.kAzure-6,"name":'4l'}
    hSum.append(hTOTElem)
    
    for h in hSum:
        setCrossSectionMC(h["state"],h["name"],Type,DoNormalized,MCSet,DoFiducial)
    return hSum

##############################################################################################################
def getCrossPlot_Data(MCSet,UseUnfold,Type,analysis,Sign,UseMCReco,DoNormalized,doFiducial):

    if UseMCReco:  print Red("########################### MC RECO ########################\n")
    else: print Red("############################ DATA  #########################\n".format(Sign))
    
    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()

    hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]    

    hsum2e2muUp = ROOT.TH1F()
    hsum4eUp    = ROOT.TH1F()
    hsum4muUp   = ROOT.TH1F()

    hSumUp = [{"state":hsum2e2muUp,"name":'2e2m'},{"state":hsum4eUp,"name":'4e'},{"state":hsum4muUp,"name":'4m'}]    

    hsum2e2muDown = ROOT.TH1F()
    hsum4eDown    = ROOT.TH1F()
    hsum4muDown   = ROOT.TH1F()

    hSumDown = [{"state":hsum2e2muDown,"name":'2e2m'},{"state":hsum4eDown,"name":'4e'},{"state":hsum4muDown,"name":'4m'}]    
  
    if UseMCReco:
        FInCenter = ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/MCReco.root")
        FInUp     = ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/DataUp.root") 
        FInDown   = ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/DataDown.root")  
    else:
        if UseUnfold:
            FInCenter =  ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/DataUnfold.root")
            FInUp     =  ROOT.TFile( "./FinalResults_"+MCSet+"_"+analysis+"/DataUnfoldUp.root")
            FInDown   =  ROOT.TFile( "./FinalResults_"+MCSet+"_"+analysis+"/DataUnfoldDown.root")  
        else: 
            FInCenter = ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/Data.root") 
            FInUp     = ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/DataUp.root")
            FInDown   = ROOT.TFile( "./FinalResults_"+MCSet+"_"+analysis+"/DataDown.root")  

    for h,hup,hdown in zip(hSum,hSumUp,hSumDown):

        h["state"] = copy.deepcopy(FInCenter.Get("ZZTo"+h["name"]+"_"+Type+"_01"))
        hup["state"] = copy.deepcopy(FInUp.Get("ZZTo"+h["name"]+"_"+Type+"_01"))
        hdown["state"] = copy.deepcopy(FInDown.Get("ZZTo"+h["name"]+"_"+Type+"_01"))

        if h==None:
            print "ERROR no data for",h["name"]
            break

        NormUp   =   hup["state"].Integral(1,-1)/h["state"].Integral(1,-1)
        NormDown = hdown["state"].Integral(1,-1)/h["state"].Integral(1,-1)

        print Blue("Central  "),
        setCrossSectionData(h["state"],h["name"],MCSet,False)
        print Blue("Sist Up  "),
        setCrossSectionData(hup["state"],hup["name"],MCSet,False)
        print Blue("Sist Down"),
        setCrossSectionData(hdown["state"],hdown["name"],MCSet,False)
        
    hTOTCross= ROOT.TH1F()
    
    #if (DoNormalized == False) and ("fr" in MCSet): (hTOTCross,hTOTCrossUp,hTOTCrossDown)=combineCrossFiducial(hSum,hSumUp,hSumDown)   
    if doFiducial: (hTOTCross,hTOTCrossUp,hTOTCrossDown)=combineCrossFiducial(hSum,hSumUp,hSumDown)   
    else:          (hTOTCross,hTOTCrossUp,hTOTCrossDown)=combineCross(hSum,hSumUp,hSumDown) 

    for hTot,h4lTot,sistSt in zip([hSum,hSumUp,hSumDown],[hTOTCross,hTOTCrossUp,hTOTCrossDown],["central","Sist Up","Sist Down"]):
        hTOTElem = {"state":h4lTot,"name":'4l'}
        print Blue(sistSt)+(" "*(9-len(sistSt))),
        print "{0} Tot Cross {1} {2:.2f} \n".format("4l",(15-len("4l"))*" ", hTOTElem["state"].Integral(1,-1)) 
        hTot.append(hTOTElem)
        for i in hTot:
            if i["state"]==None:
                print i["state"]," has no enetries" 
                
            if  "Jets" not in Type:  
                if DoNormalized: 
                    i["state"].Scale(1./i["state"].Integral(0,-1),"width")
                else:  i["state"].Scale(1.,"width") 
            else:
                if DoNormalized: 
                    i["state"].Scale(1./i["state"].Integral(0,-1),"width")#normalization to the total integral of data distribution          
    return hSum,hSumUp,hSumDown


##########################################################
################# Set Cross Section ######################

def setCrossSectionData(h1,FinState,MCSet,doFiducial):

    if FinState=='4e':     BR=BRele*BRele
    elif FinState=='4m':   BR=BRmu*BRmu
    elif FinState=='2e2m': BR=2*BRmu*BRele
    
    if doFiducial: h1.Scale(1./(Lumi)) 
    else:          h1.Scale(1./(Lumi*BR))
    
    print "{0} Tot Cross {1} {2:.2f} \n".format(FinState,(15-len(FinState))*" ", h1.Integral(0,-1))
  
##########################################################
def setCrossSectionMC(h1,FinState,Type,DoNormalized,MCSet,doFiducial):
    
    if doFiducial: BR = 1
    else:
        if FinState=='4e':     BR=BRele*BRele
        elif FinState=='4m':   BR=BRmu*BRmu
        elif FinState=='2e2m': BR=2*BRmu*BRele
        else: BR=2*BRmu*BRele

    if doFiducial:
        print "{0} Tot Cross {1} {2:.6f} [fb]\n".format(FinState,(25-len(FinState))*" ",1000*(h1.Integral(1,-1))/(Lumi)) # Check total cross section withou
        if FinState=="4l" and Type=="Total":
            AccFile = ROOT.TFile("./Acceptance/Acceptance_"+MCSet+"_"+Type+".root")
            Acc = AccFile.Get("TotAcc2e2m_Acc").GetVal() #FIXME Use a weightd avarage?
            print " Wide region",(h1.Integral(1,-1))/(Lumi*Acc*(BRele*BRele+2*BRele*BRmu+BRmu*BRmu)),"[pb]"              
    else:
        print "{0} Tot Cross {1} {2:.6f} [pb]\n".format(FinState,(25-len(FinState))*" ", (h1.Integral(1,-1))/(Lumi*BR)) # Check total cross section witho
    
    
    Integral = h1.Integral(0,-1) #Use integral with overflows entries to scale with theoretical value which include also the overflow entries.
    
    #h1.Scale(1/(Lumi*BR),"width") # If you don't want to normalize at the official cross section value.
    if DoNormalized: 
        h1.Scale(1./Integral,"width") 
    else:
        if doFiducial and ("Mass" in Type or "Jets" in Type):
            h1.Scale(1000./(Lumi),"width")  #FIXME. Use h1.Scale(norm/Integral,"width") as far as the MCFM values are available    
        else:
            norm = GetNorm(FinState,Type,doFiducial)
            h1.Scale(norm/Integral,"width")        
      
        
##################################################################################################################

##################### Combine different final state cross section distribution for final 4l ######################

####################### using weighted mean formula with both statistic ad sytematic errors ######################

##################################################################################################################


def combineCross(HList,HListUp,HListDown):

    HCross=ROOT.TH1F()
    HCrossUp=ROOT.TH1F()
    HCrossDown=ROOT.TH1F()

    #just equal to one
    HCross=copy.deepcopy(HList[1]["state"])
    HCrossUp=copy.deepcopy(HListUp[1]["state"])
    HCrossDown=copy.deepcopy(HListDown[1]["state"])

    Nbins= HList[1]["state"].GetNbinsX()

    for i in range(1,Nbins+1):

        Hlist  = zip(HList, HListUp,HListDown)
        #Sort List by entries magnitude, from higher to lower to skip 0 entries bins
        SortedHlist = sorted(Hlist,key=lambda value: value[0]["state"].GetBinContent(i),reverse = True)        

        Cross = 0
        ErrSistUp = 0
        ErrSistDown = 0
        ErrStat = 0

        WeightStat = 0.
        WeightTot = 0.
        WeightSistUp = 0.
        WeightSistDown = 0.


        for elem in SortedHlist:
            Entries = elem[0]["state"].GetBinContent(i)

            if Entries == 0.:  
                WeightTot = 1.
                WeightStat = 1.
                WeightSistUp = 1.
                WeightSistDown = 1.

                break   # Because of sorting also others final state will be 0 so use break and no continue
            
            errStatSq     = (elem[0]["state"].GetBinError(i))**2
            errSistUpSq   = (elem[1]["state"].GetBinContent(i)-Entries)**2
            errSistDownSq = (elem[2]["state"].GetBinContent(i)-Entries)**2
            errSistSq     = ((elem[1]["state"].GetBinContent(i)+elem[2]["state"].GetBinContent(i))/2.)**2 #Use the average of sistematic up and down. is it a +, right?

            weightStat =  1./errStatSq

            if errSistUpSq==0:  weightSistUp = 0.
            else:  weightSistUp =  1./errSistUpSq

            if errSistDownSq==0:  weightSistDown = 0.
            else:  weightSistDown =  1./errSistDownSq

            weightTot = 1./(errStatSq+errSistSq)
            WeightStat += weightStat
            WeightTot += weightTot
            WeightSistUp += weightSistUp
            WeightSistDown += weightSistDown 
            Cross +=  weightTot*Entries

        Cross = Cross/WeightTot

        ErrStat = math.sqrt(1./WeightStat)        
        ErrSistUp =math.sqrt(1./WeightSistUp)
        ErrSistDown =math.sqrt(1./WeightSistDown)

        HCrossUp.SetBinContent(i,Cross+ErrSistUp)
        HCrossDown.SetBinContent(i,Cross-ErrSistDown)
        HCross.SetBinContent(i,Cross)
        HCross.SetBinError(i,ErrStat)      
 
    return HCross,HCrossUp,HCrossDown

def combineCrossFiducial(HList,HListUp,HListDown):

    HCross=ROOT.TH1F()
    HCrossUp=ROOT.TH1F()
    HCrossDown=ROOT.TH1F()

#just equal to one
    HCross=copy.deepcopy(HList[1]["state"])
    HCrossUp=copy.deepcopy(HListUp[1]["state"])
    HCrossDown=copy.deepcopy(HListDown[1]["state"])

    Nbins= HList[1]["state"].GetNbinsX()

    for i in range(1,Nbins+1):

        Hlist  = zip(HList, HListUp,HListDown)
        #Sort List by entries magnitude, from higher to lower  to skip 0 entries bins
        SortedHlist = sorted(Hlist,key=lambda value: value[0]["state"].GetBinContent(i),reverse = True)        

        Cross = 0
        ErrSistUp = 0
        ErrSistDown = 0
        ErrStat = 0

        for elem in SortedHlist:
            Entries = elem[0]["state"].GetBinContent(i)
            if Entries == 0.:  break   # Because of sorting also others final state will be 0 so use break and no continue
            
            #print elem[0]["state"].GetBinContent(i),elem[1]["state"].GetBinContent(i),elem[2]["state"].GetBinContent(i)

            errStatSq = (elem[0]["state"].GetBinError(i))**2
            errSistUp = (elem[1]["state"].GetBinContent(i)-Entries)**2
            errSistDown = (elem[2]["state"].GetBinContent(i)-Entries)**2
            errSistSq =  ((elem[1]["state"].GetBinContent(i)+elem[2]["state"].GetBinContent(i))/2.)**2 #Use the average of sistematic up and down

            Cross +=  Entries #Giusto?? FIXME
            ErrStat += errStatSq 
            ErrSistUp += errSistUp
            ErrSistDown += errSistDown
       
        sigma_4e = Hlist[0][1]["state"].GetBinContent(i)
        sigma_4m = Hlist[0][2]["state"].GetBinContent(i)
        sigma_2e2m = Hlist[0][0]["state"].GetBinContent(i)
        
        sigma_4l_corrunc = math.sqrt(sigma_4e*sigma_4e +sigma_4m*sigma_4m) +sigma_2e2m 
        corrunc_4l = (0.015*sigma_4l_corrunc)**2 + (0.03*sigma_4l_corrunc)**2 #trigger unc plus scale factor unc
        
        #print sigma_4e,sigma_4m, sigma_2e2m,sigma_4l_corrunc,corrunc_4l #DEL HOT
        
        ErrStatSq = math.sqrt(ErrStat)  
        ErrSistUpSq =math.sqrt(ErrSistUp + corrunc_4l) 
        ErrSistDownSq =math.sqrt(ErrSistDown + corrunc_4l)

        #ErrSistUpSq =math.sqrt(ErrSistUp+((0.078+0.029+0.044)*Cross)**2)#uncertainties coming from the unfolding summed in quadrature + lumi 7.8% (2.6*3), Triggher 2,9% (sqrt(2*(0.015)^2)+ 0.015), SF 4.4% (sqrt(2*(0.03)^2)+ 0.03)
        #ErrSistDownSq =math.sqrt(ErrSistDown+((0.078+0.029+0.044)*Cross)**2)

        HCrossUp.SetBinContent(i,Cross+ErrSistUpSq)
        HCrossDown.SetBinContent(i,Cross-ErrSistDownSq)
        HCross.SetBinContent(i,Cross)
        HCross.SetBinError(i,ErrStatSq)      
 
    return HCross,HCrossUp,HCrossDown


def getHisto(Sign,MCSet,analysis,Type,doFiducial,UseMCReco):
    

    fIn = ROOT.TFile() 
  

    if Sign==0: 
        if UseMCReco:
            fIn = ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/MCReco.root")
        else:
            fIn = ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/Data.root") 

    elif Sign==1: fIn= ROOT.TFile("./FinalResults_"+MCSet+"_"+analysis+"/DataUp.root") 
    elif Sign==-1: fIn= ROOT.TFile( "./FinalResults_"+MCSet+"_"+analysis+"/DataDown.root")
  
    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()
   
    hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]    


    doFid = ""
    if doFiducial: doFid = "_fr" 

    for h in hSum:
        if (Type == "Mass") or (Type == "Jets"):
            h1 = fIn.Get("ZZTo"+h["name"]+"_"+Type+"_01")
            h["state"] = copy.deepcopy(h1)
        elif Type=="Tot":
            h1 = fIn.Get("ZZTo"+h["name"]+"_"+Type+doFid)
            h["state"] = copy.deepcopy(h1)
            ErrStat=ROOT.Double(0.)
            #print h["state"].IntegralAndError(0,-1,ErrStat) #HOT
            #print "err stat",ErrStat
        else:
            print Red("\n WRONG TYPE: Mass or Jets only\n") 
      
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
        
        MidSist =  (abs(CrossDic[i["name"]][2])+ abs(CrossDic[i["name"]][3]))/2.
        CrossDic[i["name"]][4] = math.sqrt(CrossDic[i["name"]][1]**2+MidSist**2) #HOT

    print Red("\n##################### RESULTS FOR INCLUSIVE CROSS SECTION ####################\n") 

    if doFiducial: list4l = combineInclusiveCrossFiducial(CrossDic)
    else:          list4l = combineInclusiveCross(CrossDic)
    dic4l = {"4l":list4l}
    CrossDic.update(dic4l)
 
    for fin,value in CrossDic.iteritems():
        value[2]=value[2]**2
        value[3]=value[3]**2 
#        value[4]=value[4]**2 
        if DoFiducial and fin == "4l": #only lumi is added #HOT
            value[2]= value[2]+(value[0]*0.026)**2
            value[3]= value[3]+(value[0]*0.026)**2
        else:            
            for sist in GlobSistList: #Add inclusive systematic errors
                value[2]+=(value[0]*sist["value"])**2
                value[3]+=(value[0]*sist["value"])**2 #CHECK
                # value[4]+=(value[0]*sist["value"])**2 #CHECK #HOT #FIX
        #value[4]= math.sqrt(value[2]+value[3])#CHECK
        value[4]= math.sqrt(value[4])
        value[2]= math.sqrt(value[2])
        value[3]= math.sqrt(value[3])

        print Blue("{0}".format(fin))
        if doFiducial:
            print  " {0:.2f} +- {1:.2f} (stat) + {2:.2f} (sist) - {3:.2f} (sist) +- {4:.2f} (Total) [fb]\n".format(value[0]*1000,value[1]*1000,value[2]*1000,value[3]*1000,value[4]*1000)
            if fin=="4l":
                AccFile = ROOT.TFile("./Acceptance/Acceptance_"+MCSet+"_"+Type+".root")
                Acc = (AccFile.Get("TotAcc2e2m_Acc").GetVal()+AccFile.Get("TotAcc4e_Acc").GetVal()+AccFile.Get("TotAcc4m_Acc").GetVal())/3. 
                print "Wide Region",value[0]/(Acc*(BRele*BRele+2*BRele*BRmu+BRmu*BRmu)),"Acc",Acc 

        else: print " {0:.2f} +- {1:.2f} (stat) + {2:.2f} (sist) - {3:.2f} (sist) +- {4:.2f} (Total) [pb] \n".format(value[0],value[1],value[2],value[3],value[4])


def combineInclusiveCross(Dic):

    TotStat     = math.sqrt(1./(1./(Dic["2e2m"][1]*Dic["2e2m"][1])+ 1./(Dic["4e"][1]*Dic["4e"][1])+ 1./(Dic["4m"][1]*Dic["4m"][1])))

    TotSistUp   = math.sqrt(1./(1./(Dic["2e2m"][2]*Dic["2e2m"][2])+ 1./(Dic["4e"][2]*Dic["4e"][2])+ 1./(Dic["4m"][2]*Dic["4m"][2])))

    TotSistDown = math.sqrt(1./(1./(Dic["2e2m"][3]*Dic["2e2m"][3])+ 1./(Dic["4e"][3]*Dic["4e"][3])+ 1./(Dic["4m"][3]*Dic["4m"][3])))
         
    TotErr      = math.sqrt(1./(1./(Dic["2e2m"][4]*Dic["2e2m"][4])+ 1./(Dic["4e"][4]*Dic["4e"][4])+ 1./(Dic["4m"][4]*Dic["4m"][4])))

    WhTot       = 1./(1./(Dic["2e2m"][4]*Dic["2e2m"][4])+1./(Dic["4e"][4]*Dic["4e"][4])+1./(Dic["4m"][4]*Dic["4m"][4]))

    TotCross    = (Dic["2e2m"][0]/(Dic["2e2m"][4]*Dic["2e2m"][4])+Dic["4e"][0]/(Dic["4e"][4]*Dic["4e"][4])+Dic["4m"][0]/(Dic["4m"][4]*Dic["4m"][4]))*WhTot

    return [TotCross,TotStat,TotSistUp,TotSistDown,TotErr]


def combineInclusiveCrossFiducial(Dic):


    TotStat= math.sqrt((Dic["2e2m"][1]*Dic["2e2m"][1])+(Dic["4e"][1]*Dic["4e"][1])+(Dic["4m"][1]*Dic["4m"][1]))

    TriggerUnc = (0.015*(math.sqrt(Dic["4e"][0]*Dic["4e"][0]+Dic["4m"][0]*+Dic["4m"][0]) + Dic["2e2m"][0]))**2

    TotSistUp  = math.sqrt((Dic["2e2m"][2])**2+(Dic["4e"][2])**2+(Dic["4m"][2])**2 +TriggerUnc)

    TotSistDown= math.sqrt((Dic["2e2m"][3])**2+(Dic["4e"][3])**2+(Dic["4m"][3])**2 +TriggerUnc)

    TotSist = (TotSistUp + TotSistDown)/2          

    TotErr= math.sqrt(TotStat*TotStat+TotSist*TotSist)

    TotCross= Dic["2e2m"][0]+Dic["4e"][0]+Dic["4m"][0]

    return [TotCross,TotStat,TotSistUp,TotSistDown,TotErr]

    
    
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

    if "Jets" in Type or "Mass" in Type:
        if doFid:
            return xs_tight[finState]
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

