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
import LatexUtils
from LatexUtils import*

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

parser.add_option("-L", "--printLatex", dest="printLatex",
                  action="store_true",
                  default=False,
                  help="Print LateX code of tabulars")



(options, args) = parser.parse_args()

Type       = options.Type
Set        = options.Set
Analysis   = options.Analysis
isFiducial = options.isFiducial
isUnfold   = options.isUnfold
printLatex = options.printLatex

SampleDic = collections.OrderedDict()
RedDic    = collections.OrderedDict()
IrrDic    = collections.OrderedDict()
DataDic   = collections.OrderedDict()
TotDic    = collections.OrderedDict()
SystDicUp = collections.OrderedDict()
SystDicDown = collections.OrderedDict()

if Analysis!="ZZ":
    inputdir_RMC = "./results/ZZRecoAnalyzer_SR_"+Analysis+"/"
    inputdir_MC  = "./results/ZZMCAnalyzer_MC_"+Analysis+"/"
    inputdir_CR  = "./results/ZZRecoAnalyzer_CR_"+Analysis+"/"


##################################################################################################################

####################  Get, sum and set histograms from data or from MC reco as they were data ####################

##################################################################################################################

def getHisto(Type,isData,Sign,syst,isTot,isFiducial):
    if isData:
        Samples=DataSamples
    else:
        if not isData: print Red("\n######################## MC RECO ########################\n")
        if "Pow" in Set:
            Samples=SignalSamples_Pow
        elif "Mad" in Set:
            Samples=SignalSamples_Mad
        else: sys.exit("Wrong Mont Carlo Set in GetHisto") 
            
    files={}
    filesbkg ={}
    
    inputdir = inputdir_RMC


    for s in Samples:
        files[s["sample"]] = ROOT.TFile(inputdir+s["sample"]+".root")
        if isData:  DataDic[s["sample"]] = {}
        else:       SampleDic[s["sample"]] = {}
    if not isData:    SampleDic["Total Signal"]={}
    
    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()

    hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]

    for h in hSum:
        if not isData: print Blue(h["name"])
        isFirst=1
        for s in Samples:
            sDic={}
            h1 = files[s["sample"]].Get("ZZTo"+h["name"]+"_"+Type+"_01")
            if h1==None:
                sDic["yield"]="-"
                sDic["Err"]="-"
                if isData: DataDic[s["sample"]][h["name"]]=sDic
                else:      SampleDic[s["sample"]][h["name"]]=sDic
               #print "For sample ", s["sample"], "h"+h["name"],"has no enetries or is a zombie"       
                continue
            if isFirst:
                h["state"]=copy.deepcopy(h1) 
                isFirst=0
                sDic["yield"]="{0:.2f}".format(h1.Integral(0,-1))
                sDic["Err"]="{0:.2f}".format(math.sqrt(h1.Integral(0,-1)))
                if isData: DataDic[s["sample"]][h["name"]]=sDic
                else:      SampleDic[s["sample"]][h["name"]]=sDic
                if not isData: print "{0} {1}> {2:.2f}".format(s["sample"],(61-len(s["sample"]))*"-",h1.Integral(0,-1))
                continue
            if not isData:
                sDic["yield"]="{0:.2f}".format(h1.Integral(0,-1))
                sDic["Err"]="{0:.2f}".format(math.sqrt(h1.Integral(0,-1)))
                SampleDic[s["sample"]][h["name"]]=sDic
                print "{0} {1}> {2:.2f}".format(s["sample"],(61-len(s["sample"]))*"-",h1.Integral(0,-1))

            h["state"].Add(h1)

        if not isData:print "\nTotal integral {0} contribution {1}> {2:.2f}\n\n".format(h["name"],(33-len(h["name"]))*"-",h["state"].Integral(0,-1))
        if isTot: h["state"].SetName((h["state"].GetName()).replace("Mass", "Total"))
        sDic={}
        sDic["yield"]="{0:.2f}".format(h["state"].Integral(0,-1))
        sDic["Err"]="{0:.2f}".format(math.sqrt(h["state"].Integral(0,-1)))
        if not isData:  SampleDic["Total Signal"][h["name"]]=sDic
    
#    print SampleDic
    hIrredBkg2e2mu = ROOT.TH1F() 
    hIrredBkg4mu   = ROOT.TH1F() 
    hIrredBkg4e    = ROOT.TH1F() 
    
    hIrredSum = [{"state":hIrredBkg2e2mu,"name":'2e2m'},{"state":hIrredBkg4e,"name":'4e'},{"state":hIrredBkg4mu,"name":'4m'}]         

    hfake2e2mu =  ROOT.TH1F() 
    hfake4e    =  ROOT.TH1F() 
    hfake4mu   =  ROOT.TH1F() 

    hFakeSum = [{"state":hfake2e2mu,"name":'2e2m'},{"state":hfake4e,"name":'4e'},{"state":hfake4mu,"name":'4m'}]

    if True: #fIXME                

        if Sign==0: print  Red("\n####################### Contribution to reducible background #######################\n")

        RedDic["Reducible background"]={}
        for hfake in hFakeSum:
            sDic = {}
            ErrFake=ROOT.Double(0.)
            hfake["state"] = getRedBkg(hfake["name"],Sign,syst)
            sDic["yield"]="{0:.2f}".format(hfake["state"].IntegralAndError(0,-1,ErrFake))
            sDic["Err"]="{0:.2f}".format(ErrFake)
            RedDic["Reducible background"][hfake["name"]]=sDic
#        print RedDic
#        if "Jet" in Type or "Deta" in Type or "Mjj" in Type : TypeString=Type+"_JERSmear" #Check
#        else: 
        TypeString=Type       
       
        if Sign==0: print Red("\n######### Contribution to Irreducible background  #########\n")
            
        for b in BkgSamples:
            filesbkg[b["sample"]] = ROOT.TFile(inputdir+b["sample"]+".root") 
            IrrDic[b["sample"]] = {}                
        IrrDic["Total Irreducible"]={}    
        for h in hIrredSum:
            if Sign==0: print Blue(h["name"])
            isFirst=1
            for b in BkgSamples:
                sDic = {}                    
                h1 = filesbkg[b["sample"]].Get("ZZTo"+h["name"]+"_"+TypeString+"_01")
                if h1==None:
                    sDic["yield"]="-"
                    sDic["Err"]="-"
                    IrrDic[b["sample"]][h["name"]]=sDic
                    print "For sample ", b["sample"], "h"+h["name"],"has no enetries or is a zombie"       
                    continue
                #print "\n",h["name"],"Total integral contribution ----------> ",h1.Integral(0,-1)
                if isFirst:
                    sDic["yield"]="{0:.2f}".format(h1.Integral(0,-1))
                    sDic["Err"]="{0:.2f}".format(math.sqrt(h1.Integral(0,-1)))
                    IrrDic[b["sample"]][h["name"]]=sDic
                    h["state"]=copy.deepcopy(h1) 
                    isFirst=0
                    if Sign==0: print "{0} {1}> {2:.2f}".format(b["sample"],(61-len(b["sample"]))*"-",h1.Integral(0,-1))
                    continue
                if Sign==0: print "{0} {1}> {2:.2f}".format(b["sample"],(61-len(b["sample"]))*"-",h1.Integral(0,-1))
                sDic["yield"]="{0:.2f}".format(h1.Integral(0,-1))
                sDic["Err"]="{0:.2f}".format(math.sqrt(h1.Integral(0,-1)))
                IrrDic[b["sample"]][h["name"]]=sDic
                h["state"].Add(h1)     
            sDic={}
            sDic["yield"]="{0:.2f}".format(h["state"].Integral(0,-1))
            sDic["Err"]="{0:.2f}".format(math.sqrt(h["state"].Integral(0,-1)))
            IrrDic["Total Irreducible"][h["name"]]=sDic
            if Sign==0:
                print "\nTotal integral {0} contribution {1}> {2:.2f}\n\n".format(h["name"],(33-len(h["name"]))*"-",h["state"].Integral(0,-1))
                #print IrrDic
        if isTot:
            for hRed,hIrr,hData in zip(hFakeSum,hSum,hIrredSum):
                hRed["state"].Rebin(4)
                hData["state"].Rebin(4)
                hIrr["state"].Rebin(4)

    if not isData:
        hTot2e2mu = ROOT.TH1F() 
        hTot4e    = ROOT.TH1F() 
        hTot4mu   = ROOT.TH1F() 
        
        hTotSum = [{"state":hTot2e2mu,"name":'2e2m'},{"state":hTot4e,"name":'4e'},{"state":hTot4mu,"name":'4m'}]         
        TotDic["Total expected"] = {}
        for h,hmc,hirr,hred in zip(hTotSum,hSum,hIrredSum,hFakeSum):
            h["state"]=hmc["state"]
            h["state"].Add(hirr["state"])
            h["state"].Add(hred["state"])
            sDic={}
            sDic["yield"]="{0:.2f}".format(h["state"].Integral(0,-1))
            sDic["Err"]="{0:.2f}".format(math.sqrt(h["state"].Integral(0,-1)))
            TotDic["Total expected"][h["name"]]=sDic

    if syst=="MCgen":
        if Set=="Pow":   AccFile = ROOT.TFile("./Acceptance/Acceptance_Mad_"+Type+".root")
        else:  AccFile = ROOT.TFile("./Acceptance/Acceptance_Pow_"+Type+".root")
    elif syst=="sFactor":
        if   Sign==-1:  AccFile = ROOT.TFile("./Acceptance/AcceptanceSFactorSqPlus_"+Set+"_"+Type+".root")#SF errors non-correlated #Check
        elif Sign==+1:  AccFile = ROOT.TFile("./Acceptance/AcceptanceSFactorSqMinus_"+Set+"_"+Type+".root")#SF errors non-correlated

        #if Sign==1: AccFile = ROOT.TFile("AcceptanceSFactorPlus_"+Set+".root")
        #elif Sign==-1:  AccFile = ROOT.TFile("AcceptanceSFactorMinus_"+Set+".root")
    else: AccFile = ROOT.TFile("./Acceptance/Acceptance_"+Set+"_"+Type+".root")
 
    if syst=="" and isData:  print Red("\n###############################  DATA  ###############################")  
    for i,j,k in zip(hSum,hIrredSum,hFakeSum):
        
        Nbins = i["state"].GetNbinsX()  

        # if syst=="":
        #     for l in range(0,Nbins+2):
        #         print l,i["state"].GetBinContent(l)
    
        if i["state"]==None:
                print i["state"]," has no enetries" 
                continue
        ErrStat=ROOT.Double(0.)

        if syst=="":
            print Blue(i["name"])
            print "\nIntegral {0}> {1:.2f} +- {2:.2f}\n ".format((55-len(i["name"]))*"-",i["state"].IntegralAndError(0,-1,ErrStat),ErrStat)
        if isData:
            if syst=="Irr":
                j["state"]=addSyst(j["state"],Sign)
             
            if syst=="Red":
                k["state"]=addSyst(k["state"],Sign)
#                print "Integral",k["state"].Integral(0,-1)
            i["state"].Add(j["state"],-1)                    
            i["state"].Add(k["state"],-1)

            for l in range(0,Nbins+2):
                if k["state"].GetBinContent(l)<0 or math.isnan(k["state"].GetBinContent(l)): 
#                    print k["name"],"Reducible bkg negative or nan content in bin ",l,k["state"].GetBinContent(l)   
                    k["state"].SetBinContent(l,0.)  #To avoid negative bin.
                    k["state"].SetBinError(l,0.)    #To avoid negative bin.
               
                if i["state"].GetBinContent(l)<0 or math.isnan(k["state"].GetBinContent(l)): 
 #                   print i["name"],"Data negative or nan content in bin ",l,i["state"].GetBinContent(l)   
                    i["state"].SetBinContent(l,0.)  #To avoid negative bin.                  
                    #i["state"].SetBinError(l,0.)    #To avoid negative bin.

            Nbins = i["state"].GetNbinsX()        
        

            ErrStat=ROOT.Double(0.)                
            if syst=="": print "Integral after background subtraction {0}> {1:.2f} +- {2:.2f}\n ".format(23*"-",i["state"].IntegralAndError(0,-1,ErrStat),ErrStat)
        if isTot:
            if isFiducial:   Acc = AccFile.Get("TotAcc"+i["name"]+"_fr").GetVal()
            else:            Acc = AccFile.Get("TotAcc"+i["name"]+"_Tot").GetVal() #FIXME
            i["state"].Scale(1./Acc)

        else:
            if isFiducial:  
                hAcc = AccFile.Get("HEff_Tight_"+i["name"]+"_"+Type)
                if hAcc== None: sys.exit("HEff_Tight_"+i["name"]+"_"+Type+" is Null or doesn't exist")
            else:                 
                hAcc = AccFile.Get("HTot_"+i["name"]+"_"+Type)
                if hAcc== None: sys.exit("HTot_"+i["name"]+"_"+Type+" is Null or doesn't exist")
            i["state"].Divide(hAcc)
        if syst=="": print "Integral after acceptance correction {0}> {1:.2f} +- {2:.2f}\n ".format((27-len(i["name"]))*"-",i["state"].IntegralAndError(0,-1,ErrStat),ErrStat)
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

def getRedBkg(FinState,Sign,syst): #FIX

    fileFake = ROOT.TFile(inputdir_CR+"data.root")
    hFakeRate=ROOT.TH1F()
    if syst=="Red" and Sign==1:
        hFakeRate = fileFake.Get("ZZTo"+FinState+"_"+Type+"_01")
        #hFakeRate.Add(setErrorsEntries(fileFake.Get("ZZTo"+FinState+"_"+Type+"_FRVarLow")),-1) 
    if syst=="Red" and Sign==-1:
        hFakeRate = fileFake.Get("ZZTo"+FinState+"_"+Type+"_01")
        ##hFakeRate.Add(setErrorsEntries(fileFake.Get("ZZTo"+FinState+"_"+Type+"_FRVarHigh")),1)
    else:
        hFakeRate = fileFake.Get("ZZTo"+FinState+"_"+Type+"_01")
    Err=ROOT.Double(0.)
    Integr= hFakeRate.IntegralAndError(0,-1,Err)
    if Sign==0: print "Total integral {0} contribution {1}> {2:.2f}\n\n".format(FinState,(33-len(FinState))*"-",hFakeRate.Integral(0,-1))
    #    print "Total integral {0} contribution {1}> {2:.2f}\n\n".format(FinState,(33-len(FinState))*"-",hFakeRate.Integral(0,-1))
    return  copy.deepcopy(hFakeRate)

##################################################################################################################

################################# Add background error for systematic distribution ###############################

##################################################################################################################

def addSyst(h,Sign):
    NBin = h.GetNbinsX()    
    for i in range(1,NBin+1):
        #print i,h.GetBinContent(i),"+-",h.GetBinError(i)
        h.AddBinContent(i,-Sign*h.GetBinError(i))
        if h.GetBinContent(i)<0: h.SetBinContent(i,0)
      
    return h 

##################################################################################################################

#########################################  Get unfolded histograms ###############################################

##################################################################################################################

def getHistoUnfold(Sign,syst,isFiducial):
    if syst != "": syst = "_"+syst
    if isFiducial: fileUnfold = ROOT.TFile("./UnfoldFolder_fr_"+Set+"/UnfoldData_"+Type+syst+".root") 
    else:          fileUnfold = ROOT.TFile("./UnfoldFolder_"+Set+"/UnfoldData_"+Type+syst+".root") 

    
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

########################### Get histograms plus or minus differential systematicss ###############################

############################ nb: inclusive systematic are added in CrossSection.py ###############################

##################################################################################################################

def getSyst(Sign,isUnfold,HData,isTot,isFiducial): 
   
    hFinSyst = copy.deepcopy(HData)
    hSystList = {}
    if Sign==+1:    print Red("Up Systematics")
    else:           print Red("Down Systematics")
    if isUnfold:
        for syst in SystList:
            hSystList[syst["name"]] =  getHistoUnfold(Sign,syst["name"],isFiducial)
    else:
        for syst in SystList:
            if Sign==+1:    SystDicUp[syst["longname"]] = {}   
            else:           SystDicDown[syst["longname"]] = {}   
            hSystList[syst["name"]] =  getHisto(Type,True,Sign,syst["name"],isTot,isFiducial)

#            hSystList[syst].SetName((hSystList[syst].GetName()).replace("Mass", "Total"))

#HHOT
#            sDic["yield"]="{0:.2f}".format(h1.Integral(0,-1))
# SampleDic[s["sample"]] = {} 

    for i in range(0,3):
        print  Blue(hSystList[syst["name"]][i]["name"])
        for syst in SystList:    
            sDic = {}
            if Sign==-1: 
                if (hSystList[syst["name"]][i]["state"].Integral(0,-1)/hFinSyst[i]["state"].Integral(0,-1)) > 1: 
#                    print syst,"is negative.If is MCgen is ok"
                    sDic["yield"]="-"
                    SystDicDown[syst["longname"]][hSystList[syst["name"]][i]["name"]]=sDic
                    continue
                sDic["yield"]="{0:.1f}".format((1-hSystList[syst["name"]][i]["state"].Integral(0,-1)/hFinSyst[i]["state"].Integral(0,-1))*100)
                print "{0} {1}-> {2:.3f} %".format(syst["longname"],(30-len(syst["longname"]))*"-",(1-hSystList[syst["name"]][i]["state"].Integral(0,-1)/hFinSyst[i]["state"].Integral(0,-1))*100)
                SystDicDown[syst["longname"]][hSystList[syst["name"]][i]["name"]]=sDic
                print hSystList[syst["name"]][i]["name"]
            else: 
                if (hFinSyst[i]["state"].Integral(0,-1)/hSystList[syst["name"]][i]["state"].Integral(0,-1)) > 1:
                    #                   print syst,"is negative. If si MCgen is ok"
                    sDic["yield"]="-"
                    SystDicUp[syst["longname"]][hSystList[syst["name"]][i]["name"]]=sDic
                    continue
                sDic["yield"]="{0:.1f}".format((1-hFinSyst[i]["state"].Integral(0,-1)/hSystList[syst["name"]][i]["state"].Integral(0,-1))*100)
                print "sDic",sDic
                print "{0} {1}-> {2:.3f} %".format(syst["longname"],(30-len(syst["longname"]))*"-",(1-hFinSyst[i]["state"].Integral(0,-1)/hSystList[syst["name"]][i]["state"].Integral(0,-1))*100)
                SystDicUp[syst["longname"]][hSystList[syst["name"]][i]["name"]]=sDic
                print SystDicUp[syst["longname"]]


        NBin = hFinSyst[i]["state"].GetNbinsX()
        print "SystDicUp",SystDicUp
      #  print "SystDicDown",SystDicDown

        for b in range(1,NBin+1):
            Content = 0
 
            # loop over histograms bins with systematics
            for syst in SystList:    

                if   Sign==-1 and hFinSyst[i]["state"].GetBinContent(b) != 0 and  (hSystList[syst["name"]][i]["state"].GetBinContent(b)/hFinSyst[i]["state"].GetBinContent(b)) > 1: continue
                elif Sign==+1 and hSystList[syst["name"]][i]["state"].GetBinContent(b) !=0 and (hFinSyst[i]["state"].GetBinContent(b)/hSystList[syst["name"]][i]["state"].GetBinContent(b)) > 1: continue

                #square sum of systematics  
                Content += (hSystList[syst["name"]][i]["state"].GetBinContent(b)-hFinSyst[i]["state"].GetBinContent(b))**2.

            hFinSyst[i]["state"].AddBinContent(b, Sign*math.sqrt(Content))

            # print "Variation",  Sign*math.sqrt(Content),"New bin content", hFinSyst[i]["state"].GetBinContent(b) 

    return hFinSyst
##################################################################################################################

##################################### Get, sum and set histograms from Gen MC ####################################

##################################################################################################################


def getPlot_MC(Type,isFiducial):

    print Red("\n############################ MC SIGNAL ############################\n")

    files={}

    # fileUnfold = ROOT.TFile("./UnfoldFolder_"+Set+"/UnfoldData_"+Type+".root")     
    # if fileUnfold.IsZombie():
    #     sys.exit("Errors! File dosn't exist")
    
    inputdir = inputdir_MC
    CrossType = Type+"Gen"

    if "Pow" in Set:
        Sample=SignalSamples_Pow
    elif "Mad" in Set:
        Samples=SignalSamples_Mad

    for s in Samples:
        files[s["sample"]] = ROOT.TFile(inputdir+s["sample"]+".root")

    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()

    list2e2mu = []
    list4mu   = []
    list4e    = []
  
    hSum = [{"state":hsum2e2mu,"name":'2e2m',"samples":list2e2mu},{"state":hsum4e,"name":'4e',"samples":list4e},{"state":hsum4mu,"name":'4m',"samples":list4mu}]

    isFromUnfold = False

    for h in hSum:
        print Blue(h["name"])
        if isFromUnfold:
            h["state"] = copy.deepcopy(fileUnfold.Get("ZZTo"+h["name"]+"_"+Type+"_GEN"))
            if isFiducial:
                h["state"].SetName( "ZZTo"+h["name"]+"_"+Type+"Gen_01_fr")
            else:
                h["state"].SetName( "ZZTo"+h["name"]+"_"+Type+"Gen_01") 
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
                    h1.SetName(s["sample"]+"_"+h["name"])
                    h["samples"].append(h1)
                    print "{0} {1}-> {2:.2f} ".format(s["sample"],(61-len(s["sample"]))*"-",h1.Integral(0,-1))
                    isFirst=0
                    continue
                h["state"].Add(h1) 
                h1.SetName(s["sample"]+"_"+h["name"])
                h["samples"].append(h1)
  
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


# Set systematic lists defined in CrossInfo.py
if isUnfold:
    SystList = DiffSystListUnfold
    if "Jet" in Type or "Deta" in Type or "Mjj" in Type: SystList = SystList+DiffSystListJetsUnfold #Add Jet systematic
else: SystList = DiffSystList

hMC =  getPlot_MC(Type,isFiducial)
#FileOutMC =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/MC.root","update") 


Fr = ""
if isFiducial: Fr="_fr"

if isTot: FileOutMC =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/MC_Total"+Fr+".root","RECREATE") 
else:     FileOutMC =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/MC_"+Type+Fr+".root","RECREATE") 

for i in hMC:
    i["state"].Write("",i["state"].kOverwrite)
    for h in i["samples"]:
        h.Write("",i["state"].kOverwrite)
    
    
if isUnfold:

    hData     = getHistoUnfold(0,"",isFiducial)
    hDataUp   = getSyst(1,isUnfold,hData,False,isFiducial)   
    hDataDown = getSyst(-1,isUnfold,hData,False,isFiducial) 
        
 
    print Red("\n#############################  UNFOLD DATA  ############################")  
    

    FileOutData =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/DataUnfold"+Fr+".root","update") 
    for i in hData:
        i["state"].Write("",i["state"].kOverwrite)
        print "\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))
        
    print Red("\n#########################  DATA UNFOLD SYST UP #########################")  
    
    FileOutDataUp =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/DataUnfoldUp"+Fr+".root","update")      
    for i in hDataUp:
        i["state"].Write("",i["state"].kOverwrite)
        print "\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1)) 

    print Red("\n########################  DATA UNFOLD SYST DOWN ########################")              

    FileOutDataDown =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/DataUnfoldDown"+Fr+".root","update") 
    for i in hDataDown:
        i["state"].Write("",i["state"].kOverwrite)
        print "\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))

    FileOutData.Close()
    
else:
    if isTot:
        hDataTot     = getHisto("Mass",True,0,"",True,isFiducial)
        hDataUpTot   = getSyst(1,isUnfold,hDataTot,True,isFiducial)
        hDataDownTot = getSyst(-1,isUnfold,hDataTot,True,isFiducial)
        hMCRecoTot   = getHisto("Mass",False,0,"",True,isFiducial)
    else:
        hData      = getHisto(Type,True,0,"",False,isFiducial) 
        hDataUp    = getSyst(1,isUnfold,hData,False,isFiducial)
        hDataDown  = getSyst(-1,isUnfold,hData,False,isFiducial)
        hMCReco    = getHisto(Type,False,0,"",False,isFiducial)

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


    FileOutData =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/Data"+Fr+".root","update") 
    if isTot:
        for i in hDataTot:
            if isFiducial:  i["state"].Write("ZZTo"+i["name"]+"_Tot_fr",i["state"].kOverwrite)
            else: i["state"].Write("ZZTo"+i["name"]+"_Tot",i["state"].kOverwrite)
            print "\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))
    else:
        for i in hData:
            i["state"].Write("",i["state"].kOverwrite)
            print "\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))
 

    print Red("\n###########################  DATA SYST UP ############################")  
    FileOutDataUp =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/DataUp"+Fr+".root","update") 
    
    if isTot: 
        for i in hDataUpTot:
            if isFiducial:  i["state"].Write("ZZTo"+i["name"]+"_Tot_fr",i["state"].kOverwrite)
            else:           i["state"].Write("ZZTo"+i["name"]+"_Tot",i["state"].kOverwrite)
            print "\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))
    else:
        for i in hDataUp:
            i["state"].Write("",i["state"].kOverwrite)
            print "\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))

    print Red("\n##########################  DATA SYST DOWN ###########################")              
    FileOutDataDown =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/DataDown"+Fr+".root","update") 
 
    if isTot:
        for i in hDataDownTot:
            if isFiducial:  i["state"].Write("ZZTo"+i["name"]+"_Tot_fr",i["state"].kOverwrite)
            else:            i["state"].Write("ZZTo"+i["name"]+"_Tot",i["state"].kOverwrite)

            print "\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))
    else:
        for i in hDataDown:
            i["state"].Write("",i["state"].kOverwrite)
            print "\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))

if printLatex: 
    print "\nYields tabular:\n"
    Text =("Sample","2e2mu" ,"4mu" ,"4e")
    YieldLatex(Text,SampleDic,RedDic,IrrDic,TotDic,DataDic)
    print "\n\nSyst tabular:\n"
    Text =("Systematic","2e2mu" ,"4mu" ,"4e")
    SystLatex(Text,SystDicUp,SystDicDown)
