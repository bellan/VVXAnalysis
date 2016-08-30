#! /usr/bin/env python
##################################
## G. Pinna (UNITO) - Jun 2015 ##
##################################

from optparse import OptionParser
import ROOT,copy
from ROOT import gSystem, TCanvas, TH1,  TPad, gStyle, TLegend,TGraphAsymmErrors,Math,TArrayD
import collections 
import CrossInfo
from CrossInfo import*
from Colours import *
import sys,ast,os
import math
import operator
import textwrap


parser = OptionParser(usage="usage: %prog [options]")

parser.add_option("-t", "--type", dest="Type",
                  default="Mass",
                  help="Mass, Total or Jets.")

parser.add_option("-f", "--fiducial", dest="isFiducial",
                  action="store_true",
                  default=False,
                  help="is Fiducial distribution")

parser.add_option("-u", "--unfold", dest="isUnfold",
                  action="store_true",
                  default=False,
                  help="is Unfold distribution")


parser.add_option("-S", "--Set", dest="Set",
                  default="Mad",
                  help="MC samples Set, default is Mad (MadGraph) the other one is Pow (Powheg)")

parser.add_option("-A", "--Analysis", dest="Analysis",
                  default="ZZ",
                  help="Analysis, default is  ZZ, other is HZZ")


(options, args) = parser.parse_args()

Type = options.Type
Set = options.Set
Analysis  = options.Analysis
isFiducial = options.isFiducial
isUnfold = options.isUnfold


if Analysis!="ZZ":
    inputdir_RMC = "./results/ZZRecoAnalyzer_SR_"+Analysis+"/"
    inputdir_MC  = "./results/ZZMCAnalyzer_MC_"+Analysis+"/"
    inputdir_CR  = "./results/ZZRecoAnalyzer_CR_"+Analysis+"/"


##################################################################################################################

####################  Get, sum and set histograms from data or from MC reco as they were data ####################

##################################################################################################################

def getHisto(Type,isData,Sign,sist,isTot,isFiducial):
    if isData:
        Samples=DataSamples
    else:
        if not isData: print Red("\n######################## MC RECO ########################\n")
        if "Pow" in Set:
            Samples=SignalSamples_Pow
        elif "Mad" in Set:
            Samples=SignalSamples_Mad
        elif "Mad" in Set:
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
            if not isData: 
                print "{0} {1}> {2:.2f}".format(s["sample"],(61-len(s["sample"]))*"-",h1.Integral(0,-1))
            h["state"].Add(h1) 
        if not isData:print "\nTotal integral {0} contribution {1}> {2:.2f}\n\n".format(h["name"],(33-len(h["name"]))*"-",h["state"].Integral(0,-1))
        if isTot: h["state"].SetName((h["state"].GetName()).replace("Mass", "Total"))
            
    hIrredBkg2e2mu = ROOT.TH1F() 
    hIrredBkg4mu   = ROOT.TH1F() 
    hIrredBkg4e    = ROOT.TH1F() 
    
    hIrredSum = [{"state":hIrredBkg2e2mu,"name":'2e2m'},{"state":hIrredBkg4e,"name":'4e'},{"state":hIrredBkg4mu,"name":'4m'}]         

    hfake2e2mu =  ROOT.TH1F() 
    hfake4e    =  ROOT.TH1F() 
    hfake4mu   =  ROOT.TH1F() 

    hFakeSum = [{"state":hfake2e2mu,"name":'2e2m'},{"state":hfake4e,"name":'4e'},{"state":hfake4mu,"name":'4m'}]

    if isData:                

        if Sign==0: print  Red("\n####################### Contribution to reducible background #######################\n")

        for hfake in hFakeSum:
            hfake["state"] = getRedBkg(hfake["name"],Sign,sist)
        # hfake4e    = getRedBkg("4e",Sign,sist)
        # hfake4mu   = getRedBkg("4m",Sign,sist)

#        if "Jet" in Type or "Deta" in Type or "Mjj" in Type : TypeString=Type+"_JERSmear"
#        else: 
        TypeString=Type       
       
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
                    if Sign==0: print "{0} {1}> {2:.2f}".format(b["sample"],(61-len(b["sample"]))*"-",h1.Integral(0,-1))
                    continue
                if Sign==0: print "{0} {1}> {2:.2f}".format(b["sample"],(61-len(b["sample"]))*"-",h1.Integral(0,-1))
                h["state"].Add(h1)     
            if Sign==0: print "\nTotal integral {0} contribution {1}> {2:.2f}\n\n".format(h["name"],(33-len(h["name"]))*"-",h["state"].Integral(0,-1))

        if isTot:
            for hRed,hIrr,hData in zip(hFakeSum,hSum,hIrredSum):
                hRed["state"].Rebin(4)
                hData["state"].Rebin(4)
                hIrr["state"].Rebin(4)

    if sist=="sFactor":
        if   Sign==-1:  AccFile = ROOT.TFile("./Acceptance/AcceptanceSFactorSqPlus_"+Set+"_"+Type+".root")#SF errors non-correlated #Check
        elif Sign==+1:  AccFile = ROOT.TFile("./Acceptance/AcceptanceSFactorSqMinus_"+Set+"_"+Type+".root")#SF errors non-correlated
        #if Sign==1: AccFile = ROOT.TFile("AcceptanceSFactorPlus_"+Set+".root")
        #elif Sign==-1:  AccFile = ROOT.TFile("AcceptanceSFactorMinus_"+Set+".root")
    else: AccFile = ROOT.TFile("./Acceptance/Acceptance_"+Set+"_"+Type+".root")
 
    if sist=="" and isData:  print Red("\n###############################  DATA  ###############################")  
    for i,j,k in zip(hSum,hIrredSum,hFakeSum):
        
        Nbins = i["state"].GetNbinsX()      
        if i["state"]==None:
                print i["state"]," has no enetries" 
                continue
        ErrStat=ROOT.Double(0.)

        if sist=="":
            print Blue(i["name"])
            print "\nIntegral {0}> {1:.2f} +- {2:.2f}\n ".format((55-len(i["name"]))*"-",i["state"].IntegralAndError(0,-1,ErrStat),ErrStat)
        if isData:
            if sist=="Irr":
                j["state"]=addSist(j["state"],Sign)

            if sist=="Red":
                k["state"]=addSist(k["state"],Sign)

            i["state"].Add(j["state"],-1)
 
            for l in range(0,Nbins+2):
                  
                if k["state"].GetBinContent(l)<0 or math.isnan(k["state"].GetBinContent(l)): 
                    #print "negative or nan content in bin ",l,k["state"].GetBinContent(l)   
                    k["state"].SetBinContent(l,0.)  #To avoid negative bin.
                    k["state"].SetBinError(l,0.)  #To avoid negative bin.
                    
            i["state"].Add(k["state"],-1)

            Nbins = i["state"].GetNbinsX()        
            
            #for l in range(1,Nbins+2):
             #   if i["state"].GetBinContent(l)<0:  i["state"].SetBinContent(l,0.)  

            ErrStat=ROOT.Double(0.)                
  #          if sist=="": print "Integral after background subtraction {0}> {1:.2f} +- {2:.2f}\n ".format((26-len(i["name"]))*"-",i["state"].IntegralAndError(0,-1,ErrStat),ErrStat)
            if sist=="": print "Integral after background subtraction {0}> {1:.2f} +- {2:.2f}\n ".format(23*"-",i["state"].IntegralAndError(0,-1,ErrStat),ErrStat)
        if isTot:
            if isFiducial:   Acc = AccFile.Get("TotAcc"+i["name"]+"_fr").GetVal()
            else:            Acc = AccFile.Get("TotAcc"+i["name"]+"_Tot").GetVal() #FIXME
#            print Acc
            i["state"].Scale(1./Acc)

        else:
            if isFiducial:         hAcc = AccFile.Get("HEff_Tight"+i["name"]+"_"+Type)
            else:                  hAcc = AccFile.Get("HTot_"+i["name"]+"_"+Type)

            if hAcc== None: sys.exit("HAcc_"+i["name"]+"_"+Type+" is Null or doesn't exist")
            i["state"].Divide(hAcc)
        if sist=="": print "Integral after acceptance correction {0}> {1:.2f} +- {2:.2f}\n ".format((27-len(i["name"]))*"-",i["state"].IntegralAndError(0,-1,ErrStat),ErrStat)
    return hSum

#################################################################################################################

def setErrorsEntries(hVar):
    Nbins=hVar.GetNbinsX()
    NegEntries = False
    for k in range(1,Nbins+1):
        if hVar.GetBinContent(k)<0:
            hVar.SetBinContent(k,0)
            NegEntries =True
        else:  hVar.SetBinContent(k,math.sqrt(hVar.GetBinContent(k)))
    if  NegEntries: print "Negative entries in Reducible background variance histogram"
    return hVar
    
#################################################################################################################

def getRedBkg(FinState,Sign,sist):

    fileFake = ROOT.TFile(inputdir_CR+"data.root")
    hFakeRate=ROOT.TH1F()
    if sist=="Red" and Sign==1:
        hFakeRate = fileFake.Get("ZZTo"+FinState+"_"+Type+"_01")
        hFakeRate.Add(setErrorsEntries(fileFake.Get("ZZTo"+FinState+"_"+Type+"_FRVarLow")),-1) 
    if sist=="Red" and Sign==-1:
        hFakeRate = fileFake.Get("ZZTo"+FinState+"_"+Type+"_01")
        hFakeRate.Add(setErrorsEntries(fileFake.Get("ZZTo"+FinState+"_"+Type+"_FRVarHigh")),1)
    else:
        hFakeRate = fileFake.Get("ZZTo"+FinState+"_"+Type+"_01")
    Err=ROOT.Double(0.)
    Integr= hFakeRate.IntegralAndError(0,-1,Err)
    
    if Sign==0: print "Total integral {0} contribution {1}> {2:.2f}\n\n".format(FinState,(33-len(FinState))*"-",hFakeRate.Integral(0,-1))
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

def getHistoUnfold(Sign,sist,isFiducial):
    if sist != "": sist = "_"+sist
    if isFiducial: fileUnfold = ROOT.TFile("./UnfoldFolder_fr_"+Set+"/UnfoldData_"+Type+sist+".root") 
    else:          fileUnfold = ROOT.TFile("./UnfoldFolder_"+Set+"/UnfoldData_"+Type+sist+".root") 

    
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

def getSist(Sign,isUnfold,HData,isTot,isFiducial): 
   
    hFinSist = copy.deepcopy(HData)
    hSistList = {}
    print Red(("Systematics added on {0} hisogram:").format(Sign))
    if isUnfold:
        for sist in SistList:
            hSistList[sist] =  getHistoUnfold(Sign,sist,isFiducial)
    else:
        for sist in SistList:
            hSistList[sist] =  getHisto(Type,True,Sign,sist,isTot,isFiducial)

#            hSistList[sist].SetName((hSistList[sist].GetName()).replace("Mass", "Total"))


    for i in range(0,3):
        print  hSistList[sist][i]["name"]
        for sist in SistList:    
            if Sign==-1: 
                print "Error from {0} {1}-> {2:.3f} %".format(sist,(15-len(sist))*"-",(1-hSistList[sist][i]["state"].Integral(0,-1)/hFinSist[i]["state"].Integral(0,-1))*100)
            else: 
                print "Error from {0} {1}-> {2:.3f} %".format(sist,(15-len(sist))*"-",(1-hFinSist[i]["state"].Integral(0,-1)/hSistList[sist][i]["state"].Integral(0,-1))*100)

        NBin = hFinSist[i]["state"].GetNbinsX()

        for b in range(1,NBin+1):
            Content = 0

            # loop over histograms bins with systematics
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


def getPlot_MC(Type,isFiducial):

    print Red("\n############################ MC SIGNAL ############################\n")

    files={}

    # fileUnfold = ROOT.TFile("./UnfoldFolder_"+Set+"/UnfoldData_"+Type+".root")     
    # if fileUnfold.IsZombie():
    # sys.exit("Errors! File dosn't exist")
    
    inputdir = inputdir_MC
    CrossType = Type+"Gen"

    if "Pow" in Set:
        Samples=SignalSamples_Pow
    elif "Mad" in Set:
        Samples=SignalSamples_Mad

    for s in Samples:
        files[s["sample"]] = ROOT.TFile(inputdir+s["sample"]+".root")

    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()

  
    hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]

    isFromUnfold = False

    for h in hSum:
        print Blue(h["name"])

        if isFromUnfold:
            h["state"] = copy.deepcopy(fileUnfold.Get("ZZTo"+h["name"]+"_"+Type+"_GEN"))
            if isFiducial:
                h["state"].SetName( "ZZTo"+h["name"]+"_"+Type+"Gen_01")
            else:
                h["state"].SetName( "ZZTo"+h["name"]+"_"+Type+"Gen_01_fr") 
        else:
            isFirst=1
            for s in Samples:
                if isFiducial: h1 = files[s["sample"]].Get("ZZTo"+h["name"]+"_"+CrossType+"_01_fr")  
                else: h1 = files[s["sample"]].Get("ZZTo"+h["name"]+"_"+CrossType+"_01")  
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
        
       # if isFiducial: h["state"].SetName( h["state"].GetName()+"_fr" )
           
        print "\nTotal integral {0} contribution {1}> {2:.2f}\n\n".format(h["name"],(34-len(h["name"]))*"-",h["state"].Integral(0,-1))


    return copy.deepcopy(hSum) 

isTot = False

if Type=="Total":
    Type="Mass"
    isTot=True

try:
    os.stat("./FinalResults_"+Set+"_"+Analysis)
except:
    os.mkdir("./FinalResults_"+Set+"_"+Analysis)       


# Set sistematic lists defined in CrossInfo.py
if isUnfold:
    SistList = DiffSistListUnfold
    if "Jet" in Type or "Deta" in Type or "Mjj" in Type: SistList = SistList+DiffSistListJetsUnfold #Add Jet systematic
else: SistList = DiffSistList

hMC =  getPlot_MC(Type,isFiducial)
FileOutMC =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/MC.root","update") 

for i in hMC:
    i["state"].Write("",i["state"].kOverwrite)
    
    
if isUnfold:

    hData     = getHistoUnfold(0,"",isFiducial)
    hDataUp   = getSist(1,isUnfold,hData,False,isFiducial)   
    hDataDown = getSist(-1,isUnfold,hData,False,isFiducial) 
        
 
    print Red("\n#############################  UNFOLD DATA  ############################")  
    

    FileOutData =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/DataUnfold.root","update") 
    for i in hData:
        i["state"].Write("",i["state"].kOverwrite)
        print "\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))
        
    print Red("\n#########################  DATA UNFOLD SIST UP #########################")  
    
    FileOutDataUp =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/DataUnfoldUp.root","update")      
    for i in hDataUp:
        i["state"].Write("",i["state"].kOverwrite)
        print "\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1)) 

    print Red("\n########################  DATA UNFOLD SIST DOWN ########################")              

    FileOutDataDown =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/DataUnfoldDown.root","update") 
    for i in hDataDown:
        i["state"].Write("",i["state"].kOverwrite)
        print "\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))

    FileOutData.Close()
    
else:
    if isTot:
        if isFiducial:
            hDataTot     = getHisto("Mass",True,0,"",True,True)
            hDataUpTot   = getSist(1,isUnfold,hDataTot,True,True)
            hDataDownTot = getSist(-1,isUnfold,hDataTot,True,True)
            hMCRecoTot   = getHisto("Mass",False,0,"",True,True)
        else:
            hDataTot     = getHisto("Mass",True,0,"",True,False)
            hDataUpTot   = getSist(1,isUnfold,hDataTot,True,False)
            hDataDownTot = getSist(-1,isUnfold,hDataTot,True,False)
            hMCRecoTot   = getHisto("Mass",False,0,"",True,False)
    else:
        hData      = getHisto(Type,True,0,"",False,False)
        hDataUp    = getSist(1,isUnfold,hData,False,False)
        hDataDown  = getSist(-1,isUnfold,hData,False,False)
        hMCReco    = getHisto(Type,False,0,"",False,False)
 
    print Red("\n##############################  SUMMARY  #############################")  

    print Red("\n###############################  MC RECO  ###############################")  
    FileOutMCReco =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/MCReco.root","update") 

    if isTot:
        for i in hMCRecoTot:
            if isFiducial:  i["state"].Write("ZZTo"+i["name"]+"_Tot_fr",i["state"].kOverwrite)
            else: i["state"].Write("ZZTo"+i["name"]+"_Tot",i["state"].kOverwrite)
            print "\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))
    else:
        for i in hMCReco:
            i["state"].Write("",i["state"].kOverwrite)
            print "\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))

    print Red("\n###############################  DATA  ###############################")  


    FileOutData =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/Data.root","update") 
    if isTot:
        for i in hDataTot:
            if isFiducial:  i["state"].Write("ZZTo"+i["name"]+"_Tot_fr",i["state"].kOverwrite)
            else: i["state"].Write("ZZTo"+i["name"]+"_Tot",i["state"].kOverwrite)
            print "\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))
    else:
        for i in hData:
            i["state"].Write("",i["state"].kOverwrite)
            print "\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))
 

    print Red("\n###########################  DATA SIST UP ############################")  
    FileOutDataUp =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/DataUp.root","update") 
    
    if isTot: 
        for i in hDataUpTot:
            if isFiducial:  i["state"].Write("ZZTo"+i["name"]+"_Tot_fr",i["state"].kOverwrite)
            else:           i["state"].Write("ZZTo"+i["name"]+"_Tot",i["state"].kOverwrite)
            print "\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))
    else:
        for i in hDataUp:
            i["state"].Write("",i["state"].kOverwrite)
            print "\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))

    print Red("\n##########################  DATA SIST DOWN ###########################")              
    FileOutDataDown =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/DataDown.root","update") 
 
    if isTot:
        for i in hDataDownTot:
            if isFiducial:  i["state"].Write("ZZTo"+i["name"]+"_Tot_fr",i["state"].kOverwrite)
            else:            i["state"].Write("ZZTo"+i["name"]+"_Tot",i["state"].kOverwrite)

            print "\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))
    else:
        for i in hDataDown:
            i["state"].Write("",i["state"].kOverwrite)
            print "\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))
