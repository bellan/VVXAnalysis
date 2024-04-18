#! /usr/bin/env python
##################################
## G. Pinna (UNITO) - Jun 2015 ##
##################################

from __future__ import print_function
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

parser.add_option("-n", "--normalize", dest="doNormalized",
                  action="store_true",
                  default=False,
                  help="set histograms for normalized distribution")

parser.add_option("-S", "--Set", dest="Set",
                  default="Mad",
                  help="MC samples Set, default is Mad (MadGraph) the other one is Pow (Powheg)")

parser.add_option("-A", "--Analysis", dest="Analysis",
                  default="ZZ",
                  help="Analysis, default is  ZZ, other is HZZ and VBS")

parser.add_option("-L", "--printLatex", dest="printLatex",
                  action="store_true",
                  default=False,
                  help="Print LateX code of tabulars")

parser.add_option("-v", "--verbose", dest="verbose",
                  action="store_true",
                  default=False,
                  help="Print LateX code of tabulars")


(options, args) = parser.parse_args()

Type         = options.Type
Set          = options.Set
Analysis     = options.Analysis
isFiducial   = options.isFiducial
isUnfold     = options.isUnfold
printLatex   = options.printLatex
verbose      = options.verbose
doNormalized = options.doNormalized

SampleDic = collections.OrderedDict()
RedDic    = collections.OrderedDict()
IrrDic    = collections.OrderedDict()
DataDic   = collections.OrderedDict()
TotDic    = collections.OrderedDict()

SystDicUp   = collections.OrderedDict()
SystDicDown = collections.OrderedDict()

SampleDic_bins = collections.OrderedDict()
RedDic_bins    = collections.OrderedDict()
IrrDic_bins    = collections.OrderedDict()
DataDic_bins   = collections.OrderedDict()
TotDic_bins    = collections.OrderedDict()

bin1 = 0;
bin2 = -1;


#ResultsFolder = "OfficialUnfoldResults"
ResultsFolder = ""

if Analysis=="VBS":
    inputdir_RMC = "./results/VBSRecoAnalyzer_SR/"
    inputdir_MC  = "./results/VBSMCAnalyzer_MC/"
    inputdir_CR  = "./results/VBSRecoAnalyzer_CR/"

    SignalSamples_Mad  = SignalZZ_VBS
    BkgSamples = [{"sample":'VBSbkg',"name":'Irreducible Background'}]
    #  BkgSamples += SignalZZ_gg 
    #  BkgSamples += SignalZZ_qq_Mad 
    DiffSystList.remove({"name":"Red","longname":"Reducible background"})

##################################################################################################################

####################  Get, sum and set histograms from data or from MC reco as they were data ####################

##################################################################################################################

def getHisto(Type,isData,Sign,syst,isTot,isFiducial):
    if isData:
        Samples=DataSamples
        lbin.setSampleTypeStatus("Data")
    else:
        if not isData: print(Red("\n######################## MC RECO ########################\n"))
        lbin.setSampleTypeStatus("Signal")
        if "Pow" in Set:
            Samples=SignalSamples_Pow
        elif "Mad" in Set:
            Samples=SignalSamples_Mad
        else: sys.exit("Wrong Mont Carlo Set in GetHisto") 

    DataSyst = "";
    if syst =="JES_Data": 
        if Sign==1: DataSyst="_JESDataUpSmear"
        else: DataSyst="_JESDataDownSmear"

    files={}
    filesbkg ={}
    
    inputdir = inputdir_RMC

    for s in Samples:
        files[s["sample"]] = ROOT.TFile(inputdir+s["sample"]+".root")
        if isData:  
            DataDic[s["sample"]] = {}
        else:      
            SampleDic[s["sample"]] = {}

    if not isData:    SampleDic["Total Signal"]={}
    
    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()

    hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]

    for h in hSum:
        if not isData: print(Blue(h["name"]))
        isFirst=True
        for s in Samples:
            sDic={}
            h1 = files[s["sample"]].Get("ZZTo"+h["name"]+"_"+Type+DataSyst+"_01")
            if h1==None:
                sDic["yield"]="-"
                sDic["Err"]="-"
                if isData:
                    lbin.setSampleTypeStatus("Data")
                    DataDic[s["sample"]][h["name"]]=sDic
                else: 
                    lbin.setSampleTypeStatus(s["sample"])
                    SampleDic[s["sample"]][h["name"]]=sDic
                for bin in range(1,h["state"].GetNbinsX()+1):  lbin.setBin(bin,h["name"],"yield",h["state"].GetBinContent())

                continue
            if isFirst:
                h["state"]=copy.deepcopy(h1) 
                isFirst=False
                sDic["yield"]="{0:.2f}".format(h1.Integral(bin1,bin2))
                sDic["Err"]="{0:.2f}".format(math.sqrt(h1.Integral(bin1,bin2)))
                if isData: 
                    DataDic[s["sample"]][h["name"]]=sDic ##del
                    lbin.setSampleStatus("Data")
                else:      
                    SampleDic[s["sample"]][h["name"]]=sDic #ddel
                    lbin.setSampleStatus(s["sample"])
                if not isData: print("{0} {1}> {2:.2f}".format(s["sample"],(61-len(s["sample"]))*"-",h1.Integral(0,-1)))
                for bin in range(1,h1.GetNbinsX()+1):  lbin.setBin(bin,h["name"],"yield",h1.GetBinContent(bin))
                continue
            if not isData:
                sDic["yield"]="{0:.2f}".format(h1.Integral(bin1,bin2)) #ddel
                sDic["Err"]="{0:.2f}".format(math.sqrt(h1.Integral(bin1,bin2))) #ddel
                SampleDic[s["sample"]][h["name"]]=sDic
                lbin.setSampleStatus(s["sample"])
                print("{0} {1}> {2:.2f}".format(s["sample"],(61-len(s["sample"]))*"-",h1.Integral(0,-1)))
                for bin in range(1,h1.GetNbinsX()+1):  lbin.setBin(bin,h["name"],"yield",h1.GetBinContent(bin))
            h["state"].Add(h1)

        if not isData:print("\nTotal integral {0} contribution {1}> {2:.2f}\n\n".format(h["name"],(33-len(h["name"]))*"-",h["state"].Integral(0,-1)))

        if isTot: h["state"].SetName((h["state"].GetName()).replace("Mass", "Total"))
        sDic={}
        sDic["yield"]="{0:.2f}".format(h["state"].Integral(bin1,bin2))
        sDic["Err"]="{0:.2f}".format(math.sqrt(h["state"].Integral(bin1,bin2)))
        if not isData:  SampleDic["Total Signal"][h["name"]]=sDic

    hIrredBkg2e2mu = ROOT.TH1F() 
    hIrredBkg4mu   = ROOT.TH1F() 
    hIrredBkg4e    = ROOT.TH1F() 
    
    hIrredSum = [{"state":hIrredBkg2e2mu,"name":'2e2m'},{"state":hIrredBkg4e,"name":'4e'},{"state":hIrredBkg4mu,"name":'4m'}]         

    hfake2e2mu =  ROOT.TH1F() 
    hfake4e    =  ROOT.TH1F() 
    hfake4mu   =  ROOT.TH1F() 

    hFakeSum = [{"state":hfake2e2mu,"name":'2e2m'},{"state":hfake4e,"name":'4e'},{"state":hfake4mu,"name":'4m'}]

    if Sign==0: print(Red("\n####################### Contribution to reducible background #######################\n"))

    RedDic["Reducible background"]={}
    lbin.setSampleTypeStatus("Background")
    lbin.setSampleStatus("Reducible")
    for hfake in hFakeSum:
        sDic = {}
        ErrFake=ROOT.Double(0.)
        if Analysis!="VBS":  hfake["state"] = getRedBkgNoZZ(hfake["name"],Sign,syst)
        else: 
            hfake["state"]=hSum[0]["state"]
            hfake["state"].Reset()
        sDic["yield"]="{0:.2f}".format(hfake["state"].IntegralAndError(bin1,bin2,ErrFake))
        sDic["Err"]="{0:.2f}".format(ErrFake)
        RedDic["Reducible background"][hfake["name"]]=sDic
        for bin in range(1,hfake["state"].GetNbinsX()+1):  lbin.setBin(bin,hfake["name"],"yield",hfake["state"].GetBinContent(bin))
    TypeString=Type       
    
    if Sign==0: print(Red("\n######### Contribution to Irreducible background  #########\n"))
    lbin.setSampleStatus("Irreducible")
    for b in BkgSamples:
        filesbkg[b["sample"]] = ROOT.TFile(inputdir+b["sample"]+".root") 
        IrrDic[b["sample"]]   = {}                
    IrrDic["Total Irreducible"]={}        
    for h in hIrredSum:
        if Sign==0: print(Blue(h["name"]))
        isFirst= True
        for b in BkgSamples:
            sDic = {}           
            lbin.setSampleStatus(b["sample"])
            h1 = filesbkg[b["sample"]].Get("ZZTo"+h["name"]+"_"+TypeString+"_01")
            ErrIrr=ROOT.Double(0.)
            if h1==None:
                sDic["yield"]="-"
                sDic["Err"]="-"
                IrrDic[b["sample"]][h["name"]]=sDic
                print("For sample ", b["sample"], "h"+h["name"],"has no enetries or is a zombie")       
                continue
              
            if isFirst:
                sDic["yield"]="{0:.2f}".format(h1.IntegralAndError(bin1,bin2,ErrIrr))
                sDic["Err"]="{0:.2f}".format(ErrIrr)
                IrrDic[b["sample"]][h["name"]]=sDic
                h["state"]=copy.deepcopy(h1) 
                isFirst=False
                for bin in range(1,h1.GetNbinsX()+1):  lbin.setBin(bin,h["name"],"yield",h1.GetBinContent(bin))
                if Sign==0: print("{0} {1}> {2:.2f} +- {3:2f}".format(b["sample"],(61-len(b["sample"]))*"-",h1.IntegralAndError(0,-1,ErrIrr),ErrIrr))
                continue
            if Sign==0: print("{0} {1}> {2:.2f} +- {3:2f}".format(b["sample"],(61-len(b["sample"]))*"-",h1.IntegralAndError(0,-1,ErrIrr),ErrIrr))
            sDic["yield"]="{0:.2f}".format(h1.IntegralAndError(bin1,bin2,ErrIrr))
            sDic["Err"]="{0:.2f}".format(ErrIrr)
            IrrDic[b["sample"]][h["name"]]=sDic
            for bin in range(1,h1.GetNbinsX()+1):  lbin.setBin(bin,h["name"],"yield",h1.GetBinContent(bin))
            h["state"].Add(h1)     
        sDic={}
        sDic["yield"]="{0:.2f}".format(h["state"].Integral(bin1,bin2))
        sDic["Err"]="{0:.2f}".format(math.sqrt(h["state"].Integral(bin1,bin2)))
        IrrDic["Total Irreducible"][h["name"]]=sDic
        if Sign==0:
            print("\nTotal integral {0} contribution {1}> {2:.2f}+- {3:2f}\n\n".format(h["name"],(33-len(h["name"]))*"-",h["state"].IntegralAndError(0,-1,ErrIrr),ErrIrr))
        
    if isTot:
        Rebin = 8
        for hRed,hData,hIrr in zip(hFakeSum,hSum,hIrredSum):
            hRed["state"].Rebin(Rebin)
            hData["state"].Rebin(Rebin)
            hIrr["state"].Rebin(Rebin)

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
            sDic["yield"]="{0:.2f}".format(h["state"].Integral(bin1,bin2))
            sDic["Err"]="{0:.2f}".format(math.sqrt(h["state"].Integral(bin1,bin2)))
            TotDic["Total expected"][h["name"]]=sDic

    if syst=="MCgen":
        if Set=="Pow":   AccFile = ROOT.TFile("./Acceptance/Acceptance_Mad_"+Type+".root")
        else:  AccFile = ROOT.TFile("./Acceptance/Acceptance_Pow_"+Type+".root")
    elif syst=="sFactor":
        if   Sign==-1:  AccFile = ROOT.TFile("./Acceptance/AcceptanceSFactorSqUp_"+Set+"_"+Type+".root")#SF errors non-correlated 
        elif Sign==+1:  AccFile = ROOT.TFile("./Acceptance/AcceptanceSFactorSqDn_"+Set+"_"+Type+".root")#SF errors non-correlated

        #if Sign==1: AccFile = ROOT.TFile("AcceptanceSFactorPlus_"+Set+".root")
        #elif Sign==-1:  AccFile = ROOT.TFile("AcceptanceSFactorMinus_"+Set+".root")
    else: AccFile = ROOT.TFile("./Acceptance/Acceptance_"+Set+"_"+Type+".root")
 
    if syst=="" and isData:  print(Red("\n###############################  DATA  ###############################"))  
    for i,j,k in zip(hSum,hIrredSum,hFakeSum):
        
        Nbins = i["state"].GetNbinsX()  
    
        if i["state"]==None:
                print(i["state"]," has no enetries") 
                continue
        ErrStat=ROOT.Double(0.)

        if syst=="":
            print(Blue(i["name"]))
            print("\nIntegral {0}> {1:.2f} +- {2:.2f}\n ".format((55-len(i["name"]))*"-",i["state"].IntegralAndError(0,-1,ErrStat),ErrStat))
        if isData:

            if syst=="Irr":
                j["state"]=addSyst(j["state"],Sign)
             
            if syst=="Red":
                k["state"].Scale(1-Sign*.3)
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


            Nbins = i["state"].GetNbinsX()        
        

            ErrStat=ROOT.Double(0.)                
            if syst=="": print("Integral after background subtraction {0}> {1:.2f} +- {2:.2f}\n ".format(23*"-",i["state"].IntegralAndError(0,-1,ErrStat),ErrStat))
        if isTot:
#            if isFiducial:   Acc = AccFile.Get("TotAcc"+i["name"]+"_fr").GetVal()
            if isFiducial:   Acc = AccFile.Get("TotAcc"+i["name"]+"_frAndFake").GetVal()
            else:            Acc = AccFile.Get("TotAcc"+i["name"]+"_Tot").GetVal() 
            i["state"].Scale(1./Acc)

        else:
            if isFiducial:  
                hAcc = AccFile.Get("HEff_Tight_"+i["name"]+"_"+Type)
                if hAcc== None: sys.exit("HEff_Tight_"+i["name"]+"_"+Type+" is Null or doesn't exist")
            else:                 
                hAcc = AccFile.Get("HTot_"+i["name"]+"_"+Type)
                if hAcc== None: sys.exit("HTot_"+i["name"]+"_"+Type+" is Null or doesn't exist")
            i["state"].Divide(hAcc)
        if syst=="": print("Integral after acceptance correction {0}> {1:.2f} +- {2:.2f}\n ".format((27-len(i["name"]))*"-",i["state"].IntegralAndError(0,-1,ErrStat),ErrStat))
    return hSum

#################################################################################################################

def setErrorsEntries(hVar): # to correct negative entries
    Nbins=hVar.GetNbinsX()
    NegEntries = False
    for k in range(1,Nbins+1):
        if hVar.GetBinContent(k)<0:
            hVar.SetBinContent(k,0)
            NegEntries =True
        else:  hVar.SetBinContent(k,math.sqrt(hVar.GetBinContent(k)))
    if  NegEntries: print("Negative entries in Reducible background variance histogram")
    return hVar
    
#################################################################################################################

def getRedBkgNoZZ(FinState,Sign,syst):

    fileFake2P2F = ROOT.TFile(inputdir_CR.replace("R/","R")+"2P2F/data.root")
    fileFake3P1F = ROOT.TFile(inputdir_CR.replace("R/","R")+"3P1F/data.root")
    fileFakeqqZZ3P1F = ROOT.TFile(inputdir_CR.replace("R/","R")+"3P1F/ZZTo4lamcatnlo.root")
    fileFakeggZZ3P1F = ROOT.TFile(inputdir_CR.replace("R/","R")+"3P1F/gg_4l.root")

    hFakeRate=ROOT.TH1F()
    if syst=="Red" and Sign==1:      TypeRed=Type+"_RedDn"
    elif syst=="Red" and Sign==-1:   TypeRed=Type+"_RedUp"
    else:                            TypeRed=Type

    hFakeRate       = fileFake3P1F.Get("ZZTo"+FinState+"_"+TypeRed+"_01")
    hFakeRate2P2F   = fileFake2P2F.Get("ZZTo"+FinState+"_"+TypeRed+"_01")
    hFakeRateggZZ3P1F = fileFakeggZZ3P1F.Get("ZZTo"+FinState+"_"+TypeRed+"_01")
    hFakeRateqqZZ3P1F = fileFakeqqZZ3P1F.Get("ZZTo"+FinState+"_"+TypeRed+"_01")

    hFakeRate.Add(hFakeRateqqZZ3P1F,-1)
    hFakeRate.Add(hFakeRateggZZ3P1F,-1)
    hFakeRate.Add(hFakeRate2P2F)

    Err=ROOT.Double(0.)
    Integr= hFakeRate.IntegralAndError(0,-1,Err)
    if Sign==0: print("Total integral {0} contribution {1}> {2:.2f}\n".format(FinState,(33-len(FinState))*"-",hFakeRate.Integral(0,-1)))
    if syst=="Red": print(Sign,"Total integral {0} contribution {1}> {2:.2f}\n\n".format(FinState,(33-len(FinState))*"-",hFakeRate.Integral(0,-1)))
    return  copy.deepcopy(hFakeRate)

#################################################################################################################                                                                    

def getRedBkg(FinState,Sign,syst):

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
    if Sign==0: print("Total integral {0} contribution {1}> {2:.2f}\n".format(FinState,(33-len(FinState))*"-",hFakeRate.Integral(0,-1)))
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
    if isFiducial: 
        fileUnfold    =  ROOT.TFile("."+ResultsFolder+"/UnfoldFolder_fr_"+Set+"/UnfoldData_"+Type+syst+".root") 
        fileUnfoldCen =  ROOT.TFile("."+ResultsFolder+"/UnfoldFolder_fr_"+Set+"/UnfoldData_"+Type+".root") 
    else:  
        fileUnfold    =  ROOT.TFile("."+ResultsFolder+"/UnfoldFolder_"+Set+"/UnfoldData_"+Type+syst+".root") 
        fileUnfoldCen =  ROOT.TFile("."+ResultsFolder+"/UnfoldFolder_"+Set+"/UnfoldData_"+Type+".root") 
    
    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()

    hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]
    
    for h in hSum:        
        if Sign ==0:
            h["state"] = copy.deepcopy(fileUnfold.Get("ZZTo"+h["name"]+"_"+Type))
        else:
            hcen  = copy.deepcopy(fileUnfoldCen.Get("ZZTo"+h["name"]+"_"+Type))
            hup   = copy.deepcopy(fileUnfold.Get("ZZTo"+h["name"]+"_"+Type+"_p"))
            hdown = copy.deepcopy(fileUnfold.Get("ZZTo"+h["name"]+"_"+Type+"_m"))

            # take same histogram for not plus/minus systemtic.
            if hup==None:
                hup   = copy.deepcopy(fileUnfold.Get("ZZTo"+h["name"]+"_"+Type))
                hdown = hup
            if doNormalized:  
                hup.Scale( hcen.Integral(0,-1)/hup.Integral(0,-1) )
                hdown.Scale( hcen.Integral(0,-1)/hdown.Integral(0,-1) )
            h["state"] = GetRightSyst(hcen,hup,hdown,Sign)

        h["state"].SetName("ZZTo"+h["name"]+"_"+Type+"_01")

    return hSum

def GetRightSyst(h,hup,hdown,sign):
    newH = copy.deepcopy(h)
    nBins = h.GetNbinsX()
    for bin in range(1,nBins+1):
        if sign ==1:
            #if(hup.GetBinContent(bin)-h.GetBinContent(bin) >0 and hdown.GetBinContent(bin)-h.GetBinContent(bin) >0): print "warning"
            newCont = max(hup.GetBinContent(bin)-h.GetBinContent(bin),hdown.GetBinContent(bin)-h.GetBinContent(bin))
            if newCont>0: newH.SetBinContent(bin,newCont+h.GetBinContent(bin))
            else:         newH.SetBinContent(bin,h.GetBinContent(bin))
        if sign ==-1:
            newCont = max(h.GetBinContent(bin)-hdown.GetBinContent(bin),h.GetBinContent(bin)-hup.GetBinContent(bin))
            if newCont>0: newH.SetBinContent(bin,h.GetBinContent(bin)-newCont)
            else:         newH.SetBinContent(bin,h.GetBinContent(bin))

    return newH



##################################################################################################################

########################### Get histograms plus or minus differential systematicss ###############################

############################ nb: inclusive systematic are added in CrossSection.py ###############################

##################################################################################################################

def getSyst(Sign,isUnfold,HData,isTot,isFiducial): 

    hFinSyst = copy.deepcopy(HData)

    # 4l combination
    h4l = copy.deepcopy(hFinSyst[0]) # copy of one channel
    h4l["name"]="4l"
    h4l["state"].SetName((hFinSyst[0]["state"].GetName()).replace("2e2m","4l"))
    # set to 0 entries and errors
    for b in range(0,h4l["state"].GetNbinsX()+2): h4l["state"].SetBinContent(b,0)
    for b in range(0,h4l["state"].GetNbinsX()+2): h4l["state"].SetBinError(b,0)

    h4lFinSyst = copy.deepcopy(h4l)

    # sum of all the channels
    for i in range(0,3):
        h4lFinSyst["state"].Add(HData[i]["state"])

    if Sign==+1:    print(Red("Up Systematics"))
    else:           print(Red("Down Systematics"))

    if isUnfold:
        for syst in SystList:
            if Sign==+1:    SystDicUp[syst["longname"]] = {}
            else:           SystDicDown[syst["longname"]] = {}
            syst["histo"] =  getHistoUnfold(Sign,syst["name"],isFiducial)
    else:
        for syst in SystList:
            if Sign==+1:    SystDicUp[syst["longname"]] = {}   
            else:           SystDicDown[syst["longname"]] = {}   

            syst["histo"] = getHisto(Type,True,Sign,syst["name"],isTot,isFiducial)

    # create 4l systemtic histogram.
    for systType in SystList:

        systType["histo"].append(copy.deepcopy(h4l))
        if systType["corr"] != 1:  
            systType["histo"][3]["state"].Add(h4lFinSyst["state"])
            hUnCorrSyst = copy.deepcopy(h4l)
        for i in range(0,3):
            if systType["corr"]==1:  
                systType["histo"][3]["state"].Add(systType["histo"][i]["state"])
            else:
                hUnc = copy.deepcopy(systType["histo"][i]["state"])
                hUnc.Add(HData[i]["state"],-1)
                hUnc.Multiply(hUnc)             
                hUnCorrSyst["state"].Add(hUnc)

        if systType["corr"] != 1: 

            for b in range(0,h4l["state"].GetNbinsX()+2): 
                systType["histo"][3]["state"].SetBinContent( b, systType["histo"][3]["state"].GetBinContent(b)+Sign*math.sqrt(hUnCorrSyst["state"].GetBinContent(b)) )

    hFinSyst.append(h4lFinSyst)

    for i in range(0,4):
        print(Blue(hFinSyst[i]["name"]))
        for syst in SystList:
            lSyst.setSampleStatus(syst["longname"])
            sDic = {}
            if Sign==-1: 
                if (syst["histo"][i]["state"].Integral(0,-1)/hFinSyst[i]["state"].Integral(0,-1) ) > 1 and 0: 
                    print(syst,"is negative.If is MCgen is ok")
                    sDic["yield"]="-"
                    SystDicDown[syst["longname"]][syst["histo"][i]["name"]]=sDic
                    continue
                sDic["yield"]="{0:.3f}".format((1-syst["histo"][i]["state"].Integral(0,-1)/hFinSyst[i]["state"].Integral(0,-1))*100)
                print("{0} {1}-> {2:.3f} %".format(syst["longname"],(30-len(syst["longname"]))*"-",(1-syst["histo"][i]["state"].Integral(0,-1)/hFinSyst[i]["state"].Integral(0,-1))*100))

                for bin in range(1,syst["histo"][i]["state"].GetNbinsX()+1):
                    if(syst["histo"][i]["state"].GetBinContent(bin) != 0): 
                        lSyst.setSystBin(bin,hFinSyst[i]["name"],Sign,((Sign*syst["histo"][i]["state"].GetBinContent(bin)-1*Sign*hFinSyst[i]["state"].GetBinContent(bin))/hFinSyst[i]["state"].GetBinContent(bin)))
                        if verbose:     print("{0} {1}-> {2:.3f} %".format(bin,(30-len(syst["longname"]))*"-",(1-syst["histo"][i]["state"].GetBinContent(bin)/hFinSyst[i]["state"].GetBinContent(bin))*100))
                SystDicDown[syst["longname"]][syst["histo"][i]["name"]]=sDic
            else: 
                if (hFinSyst[i]["state"].Integral(0,-1)/syst["histo"][i]["state"].Integral(0,-1)) > 1 and 0:
                    sDic["yield"]="-"
                    SystDicUp[syst["longname"]][syst["histo"][i]["name"]]=sDic
                    continue
                sDic["yield"]="{0:.1f}".format((1-hFinSyst[i]["state"].Integral(0,-1)/syst["histo"][i]["state"].Integral(0,-1))*100)

                print("{0} {1}-> {2:.3f} %".format(syst["longname"],(30-len(syst["longname"]))*"-",(1-hFinSyst[i]["state"].Integral(0,-1)/syst["histo"][i]["state"].Integral(0,-1))*100))
     
                for bin in range(1,syst["histo"][i]["state"].GetNbinsX()+1): 
                    #print bin,syst["histo"][i]["state"].GetBinContent(bin),hFinSyst[i]["state"].GetBinContent(bin)
                    if(syst["histo"][i]["state"].GetBinContent(bin) != 0):
                        lSyst.setSystBin(bin,hFinSyst[i]["name"],Sign,((Sign*syst["histo"][i]["state"].GetBinContent(bin)-1*Sign*hFinSyst[i]["state"].GetBinContent(bin))/hFinSyst[i]["state"].GetBinContent(bin))) 
                        if verbose:  print("{0} {1}-> {2:.3f} %".format(bin,(30-len(syst["longname"]))*"-",( ( -1+(syst["histo"][i]["state"].GetBinContent(bin))/hFinSyst[i]["state"].GetBinContent(bin))*100) ))
                SystDicUp[syst["longname"]][syst["histo"][i]["name"]]=sDic

        NBin = hFinSyst[i]["state"].GetNbinsX()
        
        if verbose:
            if hFinSyst[i]["name"]=="4l": print(Red("Overall systematics"))
        for b in range(1,NBin+1):

            Content = 0

            # loop over histograms bins with systematics
            for syst in SystList:    
                if Sign==-1 and hFinSyst[i]["state"].GetBinContent(b) != 0 and  (syst["histo"][i]["state"].GetBinContent(b)/hFinSyst[i]["state"].GetBinContent(b)) > 1: continue
                elif Sign==+1 and syst["histo"][i]["state"].GetBinContent(b) !=0 and (hFinSyst[i]["state"].GetBinContent(b)/syst["histo"][i]["state"].GetBinContent(b)) > 1: continue
                #square sum of systematics  
                Content += (syst["histo"][i]["state"].GetBinContent(b)-hFinSyst[i]["state"].GetBinContent(b))**2.

                hFinSyst[i]["state"].AddBinContent(b, Sign*math.sqrt(Content))
                if hFinSyst[i]["name"]=="4l": print(b,(math.sqrt(Content)/hFinSyst[i]["state"].GetBinContent(b))*100,"\%")

    return hFinSyst
##################################################################################################################

##################################### Get, sum and set histograms from Gen MC ####################################

##################################################################################################################


def getPlot_MC(Type,isFiducial,sign):

    files={}
    
    inputdir  = inputdir_MC
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
    hsum4l    = ROOT.TH1F()

    list2e2mu = []
    list4mu   = []
    list4e    = []
    list4l    = []

    hSum = [{"state":hsum2e2mu,"name":'2e2m',"samples":list2e2mu},{"state":hsum4e,"name":'4e',"samples":list4e},{"state":hsum4mu,"name":'4m',"samples":list4mu}]
    hSum_4l = {"state":hsum4l,"name":'4l',"samples":list4l}

    ScalVar = [{"state":{},"name":"_mf1mr1"},{"state":{},"name":"_mf1mr2"},{"state":{},"name":"_mf1mr0p5"},{"state":{},"name":"_mf2mr1"},{"state":{},"name":"_mf2mr2"},{"state":{},"name":"_mf0p5mr1"},{"state":{},"name":"_mf0p5mr0p5"}]
    SystVarUp = [{"state":{},"name":"_pdfUp"},{"state":{},"name":"_asMZUp"}]
    SystVarDn = [{"state":{},"name":"_pdfDn"},{"state":{},"name":"_asMZDn"}]

    kFacUp = {"name":"_kFacUp","state":{}}
    kFacDn = {"name":"_kFacDn","state":{}}


    for h in hSum:
        print(Blue(h["name"]))

        if sign !=0:
            isFirst = True
            
            for s in Samples:
                if isFiducial: frStr = "_fr"
                else: frStr =""
                
                if "ZZTo4l" in s["sample"]:
                    for sVar in ScalVar:
                        sVar["state"][h["name"]] = files[s["sample"]].Get("ZZTo"+h["name"]+"_"+Type+"Gen"+"_01"+frStr+sVar["name"])
                    for sVarUp,sVarDn in zip(SystVarUp,SystVarDn):
                        sVarUp["state"][h["name"]] =  files[s["sample"]].Get("ZZTo"+h["name"]+"_"+Type+"Gen"+"_01"+frStr+sVarUp["name"])
                        sVarDn["state"][h["name"]] =  files[s["sample"]].Get("ZZTo"+h["name"]+"_"+Type+"Gen"+"_01"+frStr+sVarDn["name"])

                    hCent = copy.deepcopy(files[s["sample"]].Get("ZZTo"+h["name"]+"_"+CrossType+"_01"+frStr))
                    print("{0} {1}-> {2:.2f} ".format(s["sample"],(61-len(s["sample"]))*"-",hCent.Integral(0,-1)))           
                    continue
                h1 = files[s["sample"]].Get("ZZTo"+h["name"]+"_"+CrossType+"_01"+frStr)                      
                if h1==None:
                    print("For sample ", s["sample"], "h"+h["name"],"has no enetries or is a zombie")       
                    continue


                if "gg_" in s["sample"]:
                    print("Ecco",h["name"])
                    kFacUp["state"][h["name"]] = copy.deepcopy(h1)
                    kFacDn["state"][h["name"]] = copy.deepcopy(h1)
                    kFacUp["state"][h["name"]].Scale(0.1)
                    kFacDn["state"][h["name"]].Scale(-0.1)

                if isFirst:
                    h["state"]=copy.deepcopy(h1)           
                    h1.SetName(s["sample"]+"_"+h["name"])
                    h["samples"].append(h1)
                    print("{0} {1}-> {2:.2f} ".format(s["sample"],(61-len(s["sample"]))*"-",h1.Integral(0,-1)))
                    isFirst=False
                    continue

                h["state"].Add(h1) 
                h1.SetName(s["sample"]+"_"+h["name"])
                h["samples"].append(h1)
                print("{0} {1}-> {2:.2f} ".format(s["sample"],(61-len(s["sample"]))*"-",h1.Integral(0,-1)))           

            for sVarUp,sVarDn in zip(SystVarUp,SystVarDn):
                sVarUp["state"][h["name"]].Add(h["state"])
                sVarDn["state"][h["name"]].Add(h["state"])
            for sVar in ScalVar:
                sVar["state"][h["name"]].Add(h["state"])
                #if doNormalized: 
                 #   sVarDn["state"].Scale(1/sVarDn["state"].Integral())
                 #   sVarUp["state"].Scale(1/sVarUp["state"].Integral())
            h["state"].Add(hCent) 

            kFacUp["state"][h["name"]].Add(h["state"])
            kFacDn["state"][h["name"]].Add(h["state"])

        else:
            isFirst = True
            for s in Samples:
                if isFiducial: frStr = "_fr"
                else: frStr =""
                h1 = files[s["sample"]].Get("ZZTo"+h["name"]+"_"+CrossType+"_01"+frStr)                      
                if h1==None:
                    print("For sample ", s["sample"], "h"+h["name"],"has no enetries or is a zombie")       
                    continue
                if isFirst:
                    h["state"]=copy.deepcopy(h1)           
                    h1.SetName(s["sample"]+"_"+h["name"])
                    h["samples"].append(h1)
                    print("{0} {1}-> {2:.2f} ".format(s["sample"],(61-len(s["sample"]))*"-",h1.Integral(0,-1)))
                    isFirst=False
                    continue
                h["state"].Add(h1) 
                h1.SetName(s["sample"]+"_"+h["name"])
                h["samples"].append(h1)
                print("{0} {1}-> {2:.2f} ".format(s["sample"],(61-len(s["sample"]))*"-",h1.Integral(0,-1)))           
            if doNormalized: h["state"].Scale(1/h["state"].Integral(0,-1))
            print("\nTotal integral {0} contribution {1}> {2:.2f}\n\n".format(h["name"],(34-len(h["name"]))*"-",h["state"].Integral(0,-1)))

    hSum_4l["state"]  = copy.deepcopy(hSum[0]["state"])
    hSum_4l["state"].SetName(hSum_4l["state"].GetName().replace("2e2m","4l"))
    hSum_4l["state"].Add(hSum[1]["state"])
    hSum_4l["state"].Add(hSum[2]["state"])
    
    hSum.append(hSum_4l)

    if doNormalized:
        for h in hSum:
            h["state"].Scale(1/h["state"].Integral(0,-1))
    if sign != 0:
        isFirst = True
        for sVar in ScalVar:
            ScalVar4l  = copy.deepcopy(sVar["state"]["2e2m"])
            ScalVar4l.SetName(ScalVar4l.GetName().replace("2e2m","4l"))
            ScalVar4l.Add(sVar["state"]["4e"])
            ScalVar4l.Add(sVar["state"]["4m"])
            sVar["state"]["4l"] = ScalVar4l
            
        for sVarUp,sVarDn in zip(SystVarUp,SystVarDn):
            SystVar4lUp  = copy.deepcopy(sVarUp["state"]["2e2m"])
            SystVar4lUp.SetName(SystVar4lUp.GetName().replace("2e2m","4l"))
            SystVar4lUp.Add(sVarUp["state"]["4e"])
            SystVar4lUp.Add(sVarUp["state"]["4m"])
            
            sVarUp["state"]["4l"] = SystVar4lUp
            
            SystVar4lDn  = copy.deepcopy(sVarDn["state"]["2e2m"])
            SystVar4lDn.SetName(SystVar4lDn.GetName().replace("2e2m","4l"))
            SystVar4lDn.Add(sVarDn["state"]["4e"])
            SystVar4lDn.Add(sVarDn["state"]["4m"])
            
            sVarDn["state"]["4l"] = SystVar4lDn


        kFacUp["state"]["4l"] = copy.deepcopy(kFacUp["state"]["2e2m"])
        kFacUp["state"]["4l"].SetName(kFacUp["state"]["4l"].GetName().replace("2e2m","4l"))
        kFacUp["state"]["4l"].Add(kFacUp["state"]["4e"])
        kFacUp["state"]["4l"].Add(kFacUp["state"]["4m"])

        kFacDn["state"]["4l"] = copy.deepcopy(kFacDn["state"]["2e2m"])
        kFacDn["state"]["4l"].SetName(kFacDn["state"]["4l"].GetName().replace("2e2m","4l"))
        kFacDn["state"]["4l"].Add(kFacDn["state"]["4e"])
        kFacDn["state"]["4l"].Add(kFacDn["state"]["4m"])

        if doNormalized:
            for sVarUp,sVarDn in zip(SystVarUp,SystVarDn):            

                for fin,hsyst in sVarUp["state"].items():
                    hsyst.Scale(1/hsyst.Integral(0,-1))
                for fin,hsyst in sVarDn["state"].items():
                    hsyst.Scale(1/hsyst.Integral(0,-1))

            for fin,hsyst in kFacUp["state"].items():
                hsyst.Scale(1/hsyst.Integral(0,-1))

            for fin,hsyst in kFacDn["state"].items():
                hsyst.Scale(1/hsyst.Integral(0,-1))


        SystVarDn.append(kFacDn)
        SystVarUp.append(kFacUp)


        ScaleUp = {"name":"_scaleUp","state":{}}
        ScaleDn = {"name":"_scaleDn","state":{}}

        
        for h in hSum:

            if doNormalized: h["state"].Scale(1/h["state"].Integral(0,-1))

            isFirst = True
            for sVar in ScalVar:
                if doNormalized: sVar["state"][h["name"]].Scale(1/sVar["state"][h["name"]].Integral(0,-1))
                if isFirst: 
                    ScaleUp["state"][h["name"]] = copy.deepcopy(sVar["state"][h["name"]] )
                    ScaleDn["state"][h["name"]] = copy.deepcopy(sVar["state"][h["name"]] )
                    isFirst = False
                    continue
                for bin in range(1,ScaleUp["state"][h["name"]].GetNbinsX()+1):
                    maxVal = max(ScaleUp["state"][h["name"]].GetBinContent(bin),sVar["state"][h["name"]].GetBinContent(bin) )
                    minVal = min(ScaleDn["state"][h["name"]].GetBinContent(bin),sVar["state"][h["name"]].GetBinContent(bin) )

                    ScaleUp["state"][h["name"]].SetBinContent(bin,maxVal)
                    ScaleDn["state"][h["name"]].SetBinContent(bin,minVal)

        SystVarDn.append(ScaleDn)
        SystVarUp.append(ScaleUp)

        for h in hSum:
            for hSystUp,hSystDn in zip(SystVarUp,SystVarDn): 
                for bin in range(1,hSystUp["state"][h["name"]].GetNbinsX()+1):

                    maxVal = max( hSystDn["state"][h["name"]].GetBinContent(bin),hSystUp["state"][h["name"]].GetBinContent(bin) ) 
                    minVal = min( hSystDn["state"][h["name"]].GetBinContent(bin),hSystUp["state"][h["name"]].GetBinContent(bin) ) 
                    hSystUp["state"][h["name"]].SetBinContent(bin,maxVal)
                    hSystDn["state"][h["name"]].SetBinContent(bin,minVal)


            for bin in range(1,h["state"].GetNbinsX()+1):
                binEntryUp = 0
                binEntryDn = 0
                for hSystUp,hSystDn in zip(SystVarUp,SystVarDn):

                    binEntryUp += (-h["state"].GetBinContent(bin)+ hSystUp["state"][h["name"]].GetBinContent(bin))*(-h["state"].GetBinContent(bin)+ hSystUp["state"][h["name"]].GetBinContent(bin))
                    binEntryDn += (h["state"].GetBinContent(bin) - hSystDn["state"][h["name"]].GetBinContent(bin))*(h["state"].GetBinContent(bin) - hSystDn["state"][h["name"]].GetBinContent(bin))

                if sign == +1: h["state"].SetBinContent(bin,h["state"].GetBinContent(bin)+math.sqrt(binEntryUp))
                else:  h["state"].SetBinContent(bin,h["state"].GetBinContent(bin)-math.sqrt(binEntryDn))

            if sign   == +1:  h["state"].SetName(h["state"].GetName()+"_Up") 
            elif sign == -1:  h["state"].SetName(h["state"].GetName()+"_Dn") 
            else: sys.exit("Wrong sign")

            print("\nSyst Total integral {0} contribution {1}> {2:.4f}\n\n".format(h["name"],(34-len(h["name"]))*"-",h["state"].Integral(0,-1)))                

    return copy.deepcopy(hSum) 


def GetMCsyst(finState,fileIn,sign):

    ScalVar = ["_mf1mr1","_mf1mr2","_mf1mr0p5","_mf2mr1","_mf2mr2","_mf0p5mr1","_mf0p5mr0p5"]
    if isFiducial: frStr = "_fr"
    else: frStr =""
    hOut = ROOT.TH1F()
    isFirst = True
    for Var in ScalVar:
        h1 = fileIn.Get("ZZTo"+h["name"]+"_"+Type+"Gen"+"_01"+Var+frStr)  
        if h1==None: sys.exit("ZZTo"+finState+"_"+Type+"Gen"+"_01"+Var+frStr+" doesn't exixt")
        if isFirst:
            hOut = h1
        else:
            for bin in range(0,hOut.GetNbinsX()+1):
                if sign==+1:  hOut.SetBinContent(bin,max(hOut.GetBinContent(bin),h1.GetBinContent(bin)))
                elif sign==-1:  hOut.SetBinContent(bin,min(hOut.GetBinContent(bin),h1.GetBinContent(bin)))    
    return hOut

isTot = False

if Type=="Total":
    Type  = "Mass"
    isTot = True

try:
    os.stat("./FinalResults_"+Set+"_"+Analysis)
except:
    os.mkdir("./FinalResults_"+Set+"_"+Analysis)       

# Set systematic lists defined in CrossInfo.py
if isUnfold:
    SystList = DiffSystListUnfold
    if "Jet" in Type or "Deta" in Type or "Mjj" in Type: SystList = SystList+DiffSystListJetsUnfold #Add Jet systematic
else: 
    SystList = DiffSystList
    if "Jet" in Type or "Deta" in Type or "Mjj" in Type: SystList = SystList+DiffSystListJets #Add Jet systematic

print(Red("\n############################ MC SIGNAL ############################\n"))

Variables = ["Mass","nJets","nIncJets","nJets_Central","Mjj","Mjj_Central","Deta","Deta_Central","PtJet1","PtJet2","EtaJet1","EtaJet2","dRZZ","PtZZ"]
if Type not in Variables:
    print("Wrong variable choose between", end=' ')
    for var in Variables: print(var, end=' ')
    sys.exit()

if doNormalized:
    hMC    =  getPlot_MC(Type,isFiducial,0)
    print(Red("### Syst Up ###\n "))
    hMC_Up =  getPlot_MC(Type,isFiducial,1)
    print(Red("### Syst Down ###\n"))
    hMC_Dn =  getPlot_MC(Type,isFiducial,-1)
else:
    hMC    =  getPlot_MC(Type,isFiducial,0)
    print(Red("### Syst Up ###\n "))
    hMC_Up =  getPlot_MC(Type,isFiducial,1)
    print(Red("### Syst Down ###\n"))
    hMC_Dn =  getPlot_MC(Type,isFiducial,-1)

lbin  = Latex("yield","perbin")
lSyst = Latex("Syst" ,"Syst")

lSyst.setSampleTypeStatus("Syst")

Fr = ""
if isFiducial: Fr="_fr"


if isTot: FileOutMC =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/MC_Total"+Fr+".root","update") 
else:     FileOutMC =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/MC_"+Type+Fr+".root","update") 


for i,j,k in zip(hMC,hMC_Up,hMC_Dn):
    if doNormalized:
        i["state"].SetName(i["state"].GetName()+"_norm")
        j["state"].SetName(j["state"].GetName()+"_norm")
        k["state"].SetName(k["state"].GetName()+"_norm")
    i["state"].Write("",i["state"].kOverwrite)
    j["state"].Write("",j["state"].kOverwrite)
    k["state"].Write("",k["state"].kOverwrite)
    for h,hup,hdn in zip(i["samples"],j["samples"],k["samples"]):
        h.Write("",h.kOverwrite)
        
if isUnfold:
    
    hData     = getHistoUnfold(0,"",isFiducial)
    hDataUp   = getSyst(1,isUnfold,hData,False,isFiducial)   
    hDataDown = getSyst(-1,isUnfold,hData,False,isFiducial) 

    hFinSyst = copy.deepcopy(hData)

    # 4l combination fo central value
    h4l = copy.deepcopy(hFinSyst[0]) # copy of one channel
    h4l["name"]="4l"
    h4l["state"].SetName((hFinSyst[0]["state"].GetName()).replace("2e2m","4l"))
    for b in range(0,h4l["state"].GetNbinsX()+2): h4l["state"].SetBinContent(b,0)
    for b in range(0,h4l["state"].GetNbinsX()+2): h4l["state"].SetBinError(b,0)

    h4lFin = copy.deepcopy(h4l)

    for i in range(0,3):
        h4lFin["state"].Add(hData[i]["state"])

    hData.append(h4lFin)
    

    print("central",hData[3]["name"],hData[3]["state"].Integral(),hDataUp[3]["state"].Integral(),hDataDown[3]["state"].Integral())

    print(Red("\n#############################  UNFOLD DATA  ############################"))  
    
    FileOutData =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/DataUnfold"+Fr+".root","update") 
    for i in hData:
        i["state"].Write("",i["state"].kOverwrite)
        if i["name"]=="4l" and verbose:
            print(Red("Stat uncertainties"))
            for bin in range(1,i["state"].GetNbinsX()+1):
                print(bin,(i["state"].GetBinError(bin)/i["state"].GetBinContent(bin))*100,"%")
        print("\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1)))


    print(Red("\n#########################  DATA UNFOLD SYST UP #########################"))  
    
    FileOutDataUp =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/DataUnfoldUp"+Fr+".root","update")      
    for i in hDataUp:
        if doNormalized:  i["state"].SetName(i["state"].GetName()+"_norm")
        i["state"].Write("",i["state"].kOverwrite)
        print("\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1))) 

    print(Red("\n########################  DATA UNFOLD SYST DOWN ########################"))              

    FileOutDataDown =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/DataUnfoldDown"+Fr+".root","update") 
    for i in hDataDown:
        if doNormalized:  i["state"].SetName(i["state"].GetName()+"_norm")
        i["state"].Write("",i["state"].kOverwrite)
        print("\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1)))

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


    print(Red("\n##############################  SUMMARY  #############################"))  

    print(Red("\n###############################  MC RECO  ###############################"))  
    FileOutMCReco =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/MCReco.root","update") 

    if isTot:
        for i in hMCRecoTot:
            if isFiducial:  i["state"].Write("ZZTo"+i["name"]+"_Tot_fr",i["state"].kOverwrite)
            else: i["state"].Write("ZZTo"+i["name"]+"_Tot",i["state"].kOverwrite)
            print("\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1)))
    else:
        for i in hMCReco:
            i["state"].Write("",i["state"].kOverwrite)
            print("\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1)))

    print(Red("\n###############################  DATA  ###############################"))  


    FileOutData =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/Data"+Fr+".root","update") 
    DataTot = 0
    if isTot:
        for i in hDataTot:
            if isFiducial:  i["state"].Write("ZZTo"+i["name"]+"_Tot_fr",i["state"].kOverwrite)
            else: i["state"].Write("ZZTo"+i["name"]+"_Tot",i["state"].kOverwrite)
            print("\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1)))
    else:
        for i in hData:
            i["state"].Write("",i["state"].kOverwrite)
            print("\n Total integral {0} contribution {1}> {2:.2f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1)))


 

    print(Red("\n###########################  DATA SYST UP ############################"))  
    FileOutDataUp =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/DataUp"+Fr+".root","update") 
    
    if isTot: 
        for i in hDataUpTot:
            if isFiducial:  i["state"].Write("ZZTo"+i["name"]+"_Tot_fr",i["state"].kOverwrite)
            else:           i["state"].Write("ZZTo"+i["name"]+"_Tot",i["state"].kOverwrite)
            print("\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1)))
    else:
        for i in hDataUp:
            i["state"].Write("",i["state"].kOverwrite)
            print("\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1)))

    print(Red("\n##########################  DATA SYST DOWN ###########################"))              
    FileOutDataDown =  ROOT.TFile("./FinalResults_"+Set+"_"+Analysis+"/DataDown"+Fr+".root","update") 
 
    if isTot:
        for i in hDataDownTot:
            if isFiducial:  i["state"].Write("ZZTo"+i["name"]+"_Tot_fr",i["state"].kOverwrite)
            else:            i["state"].Write("ZZTo"+i["name"]+"_Tot",i["state"].kOverwrite)

            print("\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1)))
    else:
        for i in hDataDown:
            i["state"].Write("",i["state"].kOverwrite)
            print("\n Total integral {0} contribution {1}> {2:.3f}".format(i["name"],(32-len(i["name"]))*"-",i["state"].Integral(0,-1)))

if printLatex: 
    print("\nYields tabular:\n")
    Text =("Sample","2e2mu" ,"4mu" ,"4e")
    YieldLatex(Text,SampleDic,RedDic,IrrDic,TotDic,DataDic)

    print("\n\nSyst tabular:\n")
    Text =("Systematic","2e2mu" ,"4mu" ,"4e")
    SystLatex(Text,SystDicUp,SystDicDown,isTot)
    if not isUnfold:  lbin.printPerBinCode()
    else: lSyst.printSystRanges()
