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

##################################################################################################################

####################  Get, sum and set histograms from data or from MC reco as they were data ####################

##################################################################################################################

def getHisto(Type,isData,Sign,sist):
    if isData:
        Samples=DataSamples
    else:
        if not isData: print Red("\n######################## MC RECO ########################\n")
        if Set=="Pow":
            Samples=SignalSamples_Pow
        elif Set=="Mad":
            Samples=SignalSamples_Mad

    files={}
    filesbkg ={}
    

    inputdir = inputdir_RMC

    for s in Samples:
        files[s["sample"]] = ROOT.TFile(inputdir+s["sample"]+".root")

    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()

    hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]

    for h in hSum:
        if not isData: print Blue(h["name"])
        isFirst=1
        for s in Samples:
            h1 = files[s["sample"]].Get("ZZTo"+h["name"]+"_"+Type+"_01")
            if h1==None:
               #print "For sample ", s["sample"], "h"+h["name"],"has no enetries or is a zombie"       
                continue
            if isFirst:
                h["state"]=copy.deepcopy(h1) 
                isFirst=0
                if not isData: print "{0} {1}> {2:.2f}".format(s["sample"],(61-len(s["sample"]))*"-",h1.Integral(0,-1))
                continue
            if not isData: print "{0} {1}> {2:.2f}".format(s["sample"],(61-len(s["sample"]))*"-",h1.Integral(0,-1))
            h["state"].Add(h1) 
        if not isData:print "\nTotal integral {0} contribution {1}> {2:.2f}\n\n".format(h["name"],(33-len(h["name"]))*"-",h["state"].Integral(0,-1))

            
    hIrredBkg2e2mu= ROOT.TH1F() 
    hIrredBkg4mu = ROOT.TH1F() 
    hIrredBkg4e = ROOT.TH1F() 
    
    hIrredSum = [{"state":hIrredBkg2e2mu,"name":'2e2m'},{"state":hIrredBkg4e,"name":'4e'},{"state":hIrredBkg4mu,"name":'4m'}]         

    if Sign==0: print  Red("\n####################### Contribution to reducible background #######################\n")
    hfake2e2mu=getRedBkg("2e2m",Sign)
    hfake4e=getRedBkg("4e",Sign)
    hfake4mu=getRedBkg("4m",Sign)
    
    hFakeSum = [{"state":hfake2e2mu,"name":'2e2m'},{"state":hfake4e,"name":'4e'},{"state":hfake4mu,"name":'4m'}]
                

    if isData:

        if Type=="Jets":         TypeString=Type+"_JERCentralSmear"
        else: TypeString=Type

        if Sign==0: print Red("\n######### Contribution to Irreducible background  #########\n")

        for b in BkgSamples:
            filesbkg[b["sample"]] = ROOT.TFile(inputdir+b["sample"]+".root") 
    
        for h in hIrredSum:
            if Sign==0: print Blue(h["name"])
            isFirst=1
            for b in BkgSamples:
                h1 = filesbkg[b["sample"]].Get("ZZTo"+h["name"]+"_"+TypeString+"_01")
                if h1==None:
                    print "For sample ", b["sample"], "h"+h["name"],"has no enetries or is a zombie"       
                    continue
                #print "\n",h["name"],"Total integral contribution ----------> ",h1.Integral(0,-1)
                if isFirst:
                    h["state"]=copy.deepcopy(h1) 
                    isFirst=0
                    continue
                if Sign==0: print "{0} {1}> {2:.2f}".format(b["sample"],(61-len(b["sample"]))*"-",h1.Integral(0,-1))

                h["state"].Add(h1)     
                if Sign==0: print "\nTotal integral {0} contribution {1}> {2:.2f}\n\n".format(h["name"],(33-len(h["name"]))*"-",h["state"].Integral(0,-1))
                
    if sist=="sFactor":
        if Sign==1: AccFile = ROOT.TFile("AcceptanceSFactorSqPlus_"+Set+".root")
        elif Sign==-1:  AccFile = ROOT.TFile("AcceptanceSFactorSqMinus_"+Set+".root")
        
    else: AccFile = ROOT.TFile("Acceptance_"+Set+".root")
    for i,j,k in zip(hSum,hIrredSum,hFakeSum):
        
        Nbins = i["state"].GetNbinsX()      
        if i["state"]==None:
                print i["state"]," has no enetries" 
                continue
        if isData:
            if sist=="Irr":
                j["state"]=addSist(j["state"],Sign)
            elif sist=="Red":
                k["state"]=addSist(k["state"],Sign)

            i["state"].Add(j["state"],-1)
            i["state"].Add(k["state"],-1)
            
            Nbins = i["state"].GetNbinsX()      
          
            
            for l in range(1,Nbins+2):
                if i["state"].GetBinContent(l)<0:  i["state"].SetBinContent(l,0.)  

            ErrStat=ROOT.Double(0.)                

            #            print "Before",i["state"].IntegralAndError(1,-1,ErrStat)
            #            print ErrStat
          
            
        hAcc = AccFile.Get("HAcc_"+i["name"]+"_"+Type)
        i["state"].Divide(hAcc)

#            print "After",i["state"].IntegralAndError(1,Nbins+1,ErrStat)
#            print ErrStat

            #print i["state"].Integral(),i["state"].Integral(Nbins,Nbins+1)
            #print i["state"].Integral(),i["state"].Integral(0,Nbins+1)
    return hSum

#################################################################################################################

def getRedBkg(FinState,Sign):

    fileFake = ROOT.TFile(inputdir_CR+"data.root")
    hFakeRate=ROOT.TH1F()
    hFakeRate = fileFake.Get("ZZTo"+FinState+"_"+Type+"_01")
    Err=ROOT.Double(0.)
    Integr= hFakeRate.IntegralAndError(0,-1,Err)
    if Sign==0: print "\nTotal integral {0} contribution {1}> {2:.2f}\n\n".format(FinState,(33-len(FinState))*"-",hFakeRate.Integral(0,-1))
    return  copy.deepcopy(hFakeRate)

##################################################################################################################

################################# Add background error for systematic distribution ###############################

##################################################################################################################

def addSist(h,Sign):
    NBin = h.GetNbinsX()
    
    for i in range(1,NBin+1):
        h.AddBinContent(i,-Sign*h.GetBinError(i))
        if h.GetBinContent(i)<0: h.SetBinContent(i,0)
      
    return h 

##################################################################################################################

#########################################  Get unfolded histograms ###############################################

##################################################################################################################

def getHistoUnfold(Sign,sist):

    if sist != "": sist = "_"+sist
    fileUnfold = ROOT.TFile("./macros/UnfoldingMacros/UnfoldFolder_"+Set+"/UnfoldData_"+Type+sist+".root")
    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()

    hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]
    
    if Sign == -1: Signs = "_m"
    elif Sign == +1: Signs = "_p"
    else: Signs = ""

    for h in hSum:        
        h["state"] = copy.deepcopy(fileUnfold.Get("ZZTo"+h["name"]+"_"+Type+Signs))

        # take same histogram for not plus/minus systemtic.
        if h["state"]==None:
            h["state"] = copy.deepcopy(fileUnfold.Get("ZZTo"+h["name"]+"_"+Type))
        h["state"].SetName("ZZTo"+h["name"]+"_"+Type+"_01")

    return hSum


##################################################################################################################

########################### Get histograms plus or minus differential sistematicss ###############################

############################ nb: inclusive systematic are added in CrossSection.py ###############################

##################################################################################################################

def getSist(Sign,isUnfold,HData): 
   
    hFinSist = copy.deepcopy(HData)
    hSistList = {}
    print Red(("Systematics added on {0} hisogram:").format(Sign))
    if isUnfold:
        for sist in SistList:
            hSistList[sist] =  getHistoUnfold(Sign,sist)
            print sist
    else:
        for sist in SistList:
            hSistList[sist] =  getHisto(Type,True,Sign,sist)


    for i in range(0,3):
        NBin = hFinSist[i]["state"].GetNbinsX()

        for b in range(1,NBin+1):
            Content = 0

            # loop over histograms with systematics
            for sist in SistList:    
                
                #print sist,hSistList[sist][i]["state"].GetBinContent(b),"diff",hSistList[sist][i]["state"].GetBinContent(b)-hFinSist[i]["state"].GetBinContent(b)

                #square sum of systematics  
                Content += (hSistList[sist][i]["state"].GetBinContent(b)-hFinSist[i]["state"].GetBinContent(b))**2.

            hFinSist[i]["state"].AddBinContent(b, Sign*math.sqrt(Content))

            # print "Variation",  Sign*math.sqrt(Content),"New bin content", hFinSist[i]["state"].GetBinContent(b) 

    return hFinSist
##################################################################################################################

##################################### Get, sum and set histograms from Gen MC ####################################

##################################################################################################################


def getPlot_MC(Type):

    print Red("\n############################ MC SIGNAL ############################\n")

    files={}
    
    inputdir = inputdir_MC
    CrossType = Type+"Gen"

    if Set=="Pow":
        Samples=SignalSamples_Pow
    elif Set=="Mad":
        Samples=SignalSamples_Mad

    for s in Samples:
        files[s["sample"]] = ROOT.TFile(inputdir+s["sample"]+".root")

    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()

  
    hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]

    for h in hSum:
        print Blue(h["name"])
        isFirst=1
        for s in Samples:
           h1 = files[s["sample"]].Get("ZZTo"+h["name"]+"_"+CrossType+"_01")  
           if h1==None:
               #print "For sample ", s["sample"], "h"+h["name"],"has no enetries or is a zombie"       
               continue
           if isFirst:
               h["state"]=copy.deepcopy(h1)           
               print "{0} {1}-> {2:.2f} ".format(s["sample"],(61-len(s["sample"]))*"-",h1.Integral(0,-1))
               isFirst=0
               continue
           h["state"].Add(h1) 
           print "{0} {1}-> {2:.2f} ".format(s["sample"],(61-len(s["sample"]))*"-",h1.Integral(0,-1))
        print "\nTotal integral {0} contribution {1}> {2:.2f}\n\n".format(h["name"],(3-len(h["name"]))*"-",h["state"].Integral(0,-1))


    return copy.deepcopy(hSum) 


Type= sys.argv[1]
isUnfold = sys.argv[2]
isUnfold = ast.literal_eval(isUnfold)

Set = sys.argv[3]

# Sistematic lists defined in CrossInfo.py
if isUnfold:
    SistList = DiffSistListUnfold
    if Type=="Jets": SistList = SistList+DiffSistListJetsUnfold
else: SistList = DiffSistList



hMC =  getPlot_MC(Type)

FileOutMC =  ROOT.TFile("./FinalResults_"+Set+"/MC.root","update") 

for i in hMC:
    i["state"].Write("",i["state"].kOverwrite)


hMCReco =  getHisto(Type,False,0,"")

try:
    os.stat("./FinalResults_"+Set)
except:
    os.mkdir("./FinalResults_"+Set)       

FileOutMCReco =  ROOT.TFile("./FinalResults_"+Set+"/MCReco.root","update") 

for i in hMCReco:
    i["state"].Write("",i["state"].kOverwrite)
    
    
if isUnfold:

    hData =  getHistoUnfold(0,"")

    hDataUp = getSist(1,isUnfold,hData)
    
    hDataDown = getSist(-1,isUnfold,hData) 
        

    print Red("\n#############################  UNFOLD DATA  ############################")  
    
    FileOutData =  ROOT.TFile("./FinalResults_"+Set+"/DataUnfold.root","update") 
    for i in hData:
        i["state"].Write("",i["state"].kOverwrite)
        print "\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))
        
    print Red("\n#########################  DATA UNFOLD SIST UP #########################")  
    
    FileOutDataUp =  ROOT.TFile("./FinalResults_"+Set+"/DataUnfoldUp.root","update")      
    for i in hDataUp:
        i["state"].Write("",i["state"].kOverwrite)
        print "\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1)) 

    print Red("\n########################  DATA UNFOLD SIST DOWN ########################")              

    FileOutDataDown =  ROOT.TFile("./FinalResults_"+Set+"/DataUnfoldDown.root","update") 
    for i in hDataDown:
        i["state"].Write("",i["state"].kOverwrite)
        print "\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))

    FileOutData.Close()

else:

    hData = getHisto(Type,True,0,"")
    hDataUp = getSist(1,isUnfold,hData)
    hDataDown = getSist(-1,isUnfold,hData) 

    print Red("\n###############################  DATA  ###############################")  
    FileOutData =  ROOT.TFile("./FinalResults_"+Set+"/Data.root","update") 
    for i in hData:
        i["state"].Write("",i["state"].kOverwrite)
        print "\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))

    print Red("\n###########################  DATA SIST UP ############################")  
    FileOutDataUp =  ROOT.TFile("./FinalResults_"+Set+"/DataUp.root","update") 
    for i in hDataUp:
        i["state"].Write("",i["state"].kOverwrite)
        print "\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))
        
    print Red("\n##########################  DATA SIST DOWN ###########################")              
    FileOutDataDown =  ROOT.TFile("./FinalResults_"+Set+"/DataDown.root","update") 
    for i in hDataDown:
        i["state"].Write("",i["state"].kOverwrite)
        print "\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))

    print Red("\n##############################  MC RECO  #############################")               
         
    FileOutMCREco =  ROOT.TFile("./FinalResults_"+Set+"/DataMCReco.root","update") 
    for i in hMCReco:
        i["state"].Write("",i["state"].kOverwrite)
        print "\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))



    
