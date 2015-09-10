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
import operator
from Colours import *

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

##################################################################################################################

##################### Create TGraphAsymmetricError from up and down systematic distributions #####################

##################################################################################################################


def getSistGraph(HCent,HUp,HDown):

    addGlobSist(HCent,HUp,HDown)
    grSist = ROOT.TGraphAsymmErrors(HCent)
    NBins = HCent.GetNbinsX()
    for i in range(1,NBins+1):
        grSist.SetPointEYhigh(i-1,HUp.GetBinContent(i)-HCent.GetBinContent(i))
        grSist.SetPointEYlow(i-1,HCent.GetBinContent(i)- HDown.GetBinContent(i))

    return copy.deepcopy(grSist)


##############################################################################################################

def getCrossPlot_MC(MCSet,Type):

   
    #MCRecoSample = "./FinalResults_"+MCSet+"/MCReco.root"

    fInMC  = ROOT.TFile("./FinalResults_"+MCSet+"/MC.root")  
   
    print Red("######################### Monte Carlo #######################\n")

#    Type=Type+"Gen"
    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()

    hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]    

    for h in hSum:
        h["state"] = copy.deepcopy(fInMC.Get("ZZTo"+h["name"]+"_"+Type+"Gen_01"))

        if h==None:
            print "ERROR no data for",h["name"]
            break
    
        setCrossSectionMC(h["state"],h["name"],Type)

    hTOTCross = copy.deepcopy(hSum[1]["state"])
    hTOTElem = {"state":hTOTCross,"color":ROOT.kAzure-6,"name":'4l'}
    hSum.append(hTOTElem)

    return hSum

##############################################################################################################
def getCrossPlot_Data(MCSet,UseUnfold,Type,Sign,UseMCReco):

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



# CenterData = "./FinalResults_"+MCSet+"/Data.root"
# UpData = "./FinalResults_"+MCSet+"/DataUp.root"
# DownData =  "./FinalResults_"+MCSet+"/DataDown.root"

# CenterDataUF = "./FinalResults_"+MCSet+"/DataUnfold.root"
# UpDataUF = "./FinalResults_"+MCSet+"/DataUnfoldUp.root"
# DownDataUF = "./FinalResults_"+MCSet+"/DataUnfoldDown.root"

# fInCenter  = ROOT.TFile(CenterData)
# fInUp  =  ROOT.TFile(UpData)
# fInDown  = ROOT.TFile(DownData)  

# fInCenterUF  = ROOT.TFile(CenterDataUF)
# fInUpUF  =  ROOT.TFile(UpDataUF)
# fInDownUF  = ROOT.TFile(DownDataUF)  

# fInMC  = ROOT.TFile(MCSample)  
# fInMCReco  = ROOT.TFile(MCRecoSample)

# fInCenter  = ROOT.TFile("./FinalResults_"+MCSet+"/Data.root")
# fInUp  =  ROOT.TFile("./FinalResults_"+MCSet+"/DataUp.root")
# fInDown  = ROOT.TFile( "./FinalResults_"+MCSet+"/DataDown.root")  

# fInCenterUF  = ROOT.TFile("./FinalResults_"+MCSet+"/DataUnfold.root")
# fInUpUF  =  ROOT.TFile( "./FinalResults_"+MCSet+"/DataUnfoldUp.root")
# fInDownUF  = ROOT.TFile( "./FinalResults_"+MCSet+"/DataUnfoldDown.root")  

# fInMC  = ROOT.TFile(MCSample)  
# fInMCReco  = ROOT.TFile(MCRecoSample)

#     if UseMCReco:
#         FInCenter = fInMCReco
#         FInUp=fInUp
#         FInDown=fInDown
#     else:
#         if UseUnfold:
#             FInCenter=fInCenterUF                    
#             FInUp=fInUpUF
#             FInDown=fInDownUF
#         else: 
#             FInCenter=fInCenter
#             FInUp=fInUp
#             FInDown=fInDown

  
    if UseMCReco:
        FInCenter = ROOT.TFile("./FinalResults_"+MCSet+"/MCReco.root")
        FInUp     = ROOT.TFile("./FinalResults_"+MCSet+"/DataUp.root") 
        FInDown   = ROOT.TFile("./FinalResults_"+MCSet+"/DataDown.root")  
    else:
        if UseUnfold:
            FInCenter =  ROOT.TFile("./FinalResults_"+MCSet+"/DataUnfold.root")
            FInUp     =  ROOT.TFile( "./FinalResults_"+MCSet+"/DataUnfoldUp.root")
            FInDown   =  ROOT.TFile( "./FinalResults_"+MCSet+"/DataUnfoldDown.root")  
        else: 
            FInCenter = ROOT.TFile("./FinalResults_"+MCSet+"/Data.root") 
            FInUp     = ROOT.TFile("./FinalResults_"+MCSet+"/DataUp.root")
            FInDown   = ROOT.TFile( "./FinalResults_"+MCSet+"/DataDown.root")  

    for h,hup,hdown in zip(hSum,hSumUp,hSumDown):

        h["state"] = copy.deepcopy(FInCenter.Get("ZZTo"+h["name"]+"_"+Type+"_01"))
        hup["state"] = copy.deepcopy(FInUp.Get("ZZTo"+h["name"]+"_"+Type+"_01"))
        hdown["state"] = copy.deepcopy(FInDown.Get("ZZTo"+h["name"]+"_"+Type+"_01"))

        if h==None:
            print "ERROR no data for",h["name"]
            break
        print Blue("Central  "),
        setCrossSectionData(h["state"],h["name"])
        print Blue("Sist Up  "),
        setCrossSectionData(hup["state"],hup["name"])
        print Blue("Sist Down"),
        setCrossSectionData(hdown["state"],hdown["name"])
        
    hTOTCross= ROOT.TH1F()
    
    (hTOTCross,hTOTCrossUp,hTOTCrossDown)=combineCross(hSum,hSumUp,hSumDown)   
    

    for hTot,h4lTot,sistSt in zip([hSum,hSumUp,hSumDown],[hTOTCross,hTOTCrossUp,hTOTCrossDown],["central","Sist Up","Sist Down"]):
        hTOTElem = {"state":h4lTot,"name":'4l'}
        print Blue(sistSt)+(" "*(9-len(sistSt))),
        print "{0} Tot Cross {1} {2:.2f} \n".format("4l",(15-len("4l"))*" ", hTOTElem["state"].Integral(1,-1)) 
        hTot.append(hTOTElem)
    
        for i in hTot:
            if i["state"]==None:
                print i["state"]," has no enetries" 
            if "Jets" not in Type:
                i["state"].Scale(1.,"width")
    return hSum,hSumUp,hSumDown


#########################################################
################# Set Cros Section ######################

def setCrossSectionData(h1,FinState):

    if FinState=='4e':     BR=BRele*BRele
      
    elif FinState=='4m':   BR=BRmu*BRmu

    elif FinState=='2e2m': BR=2*BRmu*BRele

    h1.Scale(1/(Lumi*BR))

    print "{0} Tot Cross {1} {2:.2f} \n".format(FinState,(15-len(FinState))*" ", h1.Integral(1,-1))


#####################################################
def setCrossSectionMC(h1,FinState,Type):
    
    if FinState=='4e':      BR=BRele*BRele

    elif FinState=='4m':   BR=BRmu*BRmu

    elif FinState=='2e2m': BR=2*BRmu*BRele

    print "{0} Tot Cross {1} {2:.2f} \n".format(FinState,(25-len(FinState))*" ", (h1.Integral(1,-1))/(Lumi*BR)) # Check total cross section without normalization
    
    #h1.Scale(1/(Lumi*BR),"width") # If you don't want to normalize at the official cross section value.

    #Use integral with overflows entries to scale with theoretical value which include also the overflow entries.

    Integral = h1.Integral(0,-1) 
   # print "Intagral before",Integral
    if Type == "Mass":   
        h1.Scale(7.5/Integral,"width") 
    elif "Jets" in Type: 
        h1.Scale(7.5/Integral)     
       # print "Int34",h1.Integral(3,4)
    elif "PtJet1" in Type:
        h1.Scale(2.97096294176353171/Integral,"width") # Mjj and Deta don't contain all the events of theoretical cross section 7,5. 
    else:
        #h1.Scale(1/(BR*Lumi),"width") # Mjj and Deta don't contain all the events of theoretical cross section 7,5.
        #print "Int before",h1.Integral()
        h1.Scale(6.20670780539512634e-01/Integral,"width") # Mjj and Deta don't contain all the events of theoretical cross section 7,5. 
        #h1.Scale(6.20670780539512634e-01/Integral) # Mjj and Deta don't contain all the events of theoretical cross section 7,5. 
        print "Int",h1.Integral(0,-1)


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
        #Sort List by entries magnitude, from higher to lower  to skip 0 entries bins
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
            
          
            if Entries == 0.:  break   # Because of sorting also others final state will be 0 so use break and no continue

            errStatSq = (elem[0]["state"].GetBinError(i))**2
            errSistUpSq = (elem[1]["state"].GetBinContent(i)-Entries)**2
            errSistDownSq = (elem[2]["state"].GetBinContent(i)-Entries)**2
            errSistSq =  ((elem[1]["state"].GetBinContent(i)-elem[2]["state"].GetBinContent(i))/2.)**2 #Use the average of sistematic up and down

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


def getHisto(Sign):
    
    fIn = ROOT.TFile()
    if Sign==0:    fIn=fInCenter
    elif Sign==1:  fIn=fInUp
    elif Sign==-1: fIn=fInDown

    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()
    
    hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]    
    for h in hSum:
        h1 = fIn.Get("ZZTo"+h["name"]+"_Mass_01")
        h["state"] = copy.deepcopy(h1)
       
    return hSum


def TotalCross():

    hData = getHisto(0)
    hDataUp = getHisto(1)
    hDataDown = getHisto(-1)
    
    list2e2m = [0,0,0,0,0]
    list4e = [0,0,0,0,0]
    list4m = [0,0,0,0,0]
    
    CrossDic = collections.OrderedDict()
    CrossDic["2e2m"]=list2e2m
    CrossDic["4e"]=list4e
    CrossDic["4m"]=list4m
    

    for i,j,k in zip(hData,hDataUp,hDataDown):

        if i["name"]=='4e':     BR=BRele*BRele
        elif i["name"]=='4m':   BR=BRmu*BRmu
        elif i["name"]=='2e2m': BR=2*BRmu*BRele
        if i["name"]=='4l': break

        ErrStat=ROOT.Double(0.)
        CrossDic[i["name"]][0] = i["state"].IntegralAndError(0,-1,ErrStat)/(BR*Lumi) 
        CrossDic[i["name"]][1] = ErrStat/(BR*Lumi)
       
        CrossDic[i["name"]][2] = j["state"].Integral(0,-1)/(BR*Lumi) - CrossDic[i["name"]][0]     
        CrossDic[i["name"]][3] = CrossDic[i["name"]][0] - k["state"].Integral(0,-1)/(BR*Lumi)
        MidSist =  (CrossDic[i["name"]][2]+ CrossDic[i["name"]][3])/2.
        
        CrossDic[i["name"]][4] = math.sqrt(CrossDic[i["name"]][1]**2+MidSist**2)
    
    list4l = combineInclusiveCross(CrossDic)

    dic4l = {"4l":list4l}
    CrossDic.update(dic4l)
 
    print Red("\n##################### RESULTS FOR INCLUSIVE CROSS SECTION ####################\n") 
    for fin,value in CrossDic.iteritems():
        value[2]=value[2]**2
        value[3]=value[3]**2
        for sist in GlobSistList: #Add inclusive systematic errors
            value[2]+=(value[0]*sist["value"])**2
            value[3]+=(value[0]*sist["value"])**2 #CHEC
        value[4]= math.sqrt(value[2]+value[3])
        value[2]=math.sqrt(value[2])
        value[3]=math.sqrt(value[3])


        print Blue("{0}".format(fin)),
        print " {0:.2f} +- {1:.2f} (stat) + {2:.2f} (sist) - {3:.2f} (sist) +- {4:.2f} (Total)\n".format(value[0],value[1],value[2],value[3],value[4])
    

def combineInclusiveCross(Dic):

    TotStat= math.sqrt(1./(1./(Dic["2e2m"][1]*Dic["2e2m"][1])+ 1./(Dic["4e"][1]*Dic["4e"][1])+ 1./(Dic["4m"][1]*Dic["4m"][1])))

    TotSistUp= math.sqrt(1/(1/(Dic["2e2m"][2]*Dic["2e2m"][2])+ 1/(Dic["4e"][2]*Dic["4e"][2])+ 1/(Dic["4m"][2]*Dic["4m"][2])))

    TotSistDown= math.sqrt(1/(1/(Dic["2e2m"][3]*Dic["2e2m"][3])+ 1/(Dic["4e"][3]*Dic["4e"][3])+ 1/(Dic["4m"][3]*Dic["4m"][3])))
         
    TotErr= math.sqrt(1/(1/(Dic["2e2m"][4]*Dic["2e2m"][4])+ 1/(Dic["4e"][4]*Dic["4e"][4])+ 1/(Dic["4m"][4]*Dic["4m"][4])))

    WhTot= 1/(1/(Dic["2e2m"][4]*Dic["2e2m"][4])+1/(Dic["4e"][4]*Dic["4e"][4])+1/(Dic["4m"][2]*Dic["4m"][4]))

    TotCross = (Dic["2e2m"][0]/(Dic["2e2m"][4]*Dic["2e2m"][4])+Dic["4e"][0]/(Dic["4e"][4]*Dic["4e"][4])+Dic["4m"][0]/(Dic["4m"][2]*Dic["4m"][4]))*WhTot

    return [TotCross,TotStat,TotSistUp,TotSistDown,TotErr]



   
