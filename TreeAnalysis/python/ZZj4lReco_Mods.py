#! /usr/bin/env python

import ROOT, copy, sys
from ROOT import gSystem, TCanvas, TH1,  TPad, gStyle, TLegend, THStack, TGraphAsymmErrors,Math
from readSampleInfo import *
from collections import OrderedDict
from Colours import *
import re

def GetTypeofsamples(category,Set):
    
    signal_qq_Pow = [{"sample":'ZZTo2e2mu',"color":ROOT.kAzure-4,"name":'qq/qg #rightarrow ZZ(+jets)'},{"sample":'ZZTo4e',"color":ROOT.kAzure-4,"name":'qq/qg #rightarrow ZZ(+jets)'},{"sample":'ZZTo4mu',"color":ROOT.kAzure-4,"name":'qq/qg #rightarrow ZZ(+jets)'}]
    
    signal_qq_Mad = [{"sample":'ZZJetsTo4L',"color":ROOT.kAzure-4,"name":'qq/qg/gg #rightarrow ZZ(+jets)'}]
    
    signal_gg=[{"sample":'ggTo2e2mu_SMHContinInterf-MCFM67_H125.6',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)'},{"sample":'ggTo4e_SMHContinInterf-MCFM67_H125.6',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)'},{"sample":'ggTo4mu_SMHContinInterf-MCFM67_H125.6',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)'}]
    
    signal_VBS = [{"sample":'ZZTo2e2muJJ_SMHContinInterf_H125.6',"color":ROOT.kAzure-6,"name":'other ZZ processes'},{"sample":'ZZTo4eJJ_SMHContinInterf_H125.6',"color":ROOT.kAzure-6,"name":'other ZZ processes'},{"sample":'ZZTo4muJJ_SMHContinInterf_H125.6',"color":ROOT.kAzure-6,"name":'other ZZ processes'}]
    
    signal_other = [{"sample":'ZZZJets',"color":ROOT.kAzure-6,"name":'other ZZ processes'},{"sample":'WZZJets',"color":ROOT.kAzure-6,"name":'other ZZ processes'},{"sample":'WH126',"color":ROOT.kAzure-6,"name":'ZZ processes'},{"sample":'ZH126',"color":ROOT.kAzure-6,"name":'other ZZ processes'},{"sample":'ttH126',"color":ROOT.kAzure-6,"name":'other ZZ processes'}] # check if there is everything
    
    if Set=="pow":    signal_tot = signal_qq_Pow+signal_gg+signal_VBS+signal_other
    elif Set=="mad":  signal_tot = signal_qq_Mad+signal_gg+signal_VBS+signal_other
    else: sys.exit("Check MC set name. Possibilty are pow and mad")    
    
    irred_background_tot = [{"sample":'WWZJets',"color":ROOT.kOrange-3,"name":'Irreducible background'},{"sample":'TTWWJets',"color":ROOT.kOrange-3,"name":'Irreducible background'},{"sample":'TTZJets',"color":ROOT.kOrange-3,"name":'Irreducible background'}]    


    background_red = [{"sample":'WWWJets',"color":ROOT.kBlue-9,"name":'Triboson'},{"sample":'WZ',"color":ROOT.kOrange-4,"name":'WZ'},{"sample":'TTTo2L2Nu2B',"color":ROOT.kRed-2,"name":'tt'},{"sample":'TTWJets',"color":ROOT.kRed-2,"name":'tt'},{"sample":'TTGJets',"color":ROOT.kRed-2,"name":'tt'},{"sample":'DYJetsToLLTuneZ2M50',"color":ROOT.kGreen-5,"name":'DY'},{"sample":'DYJetsToLLTuneZ2M10',"color":ROOT.kGreen-5,"name":'DY'}]

    data = [{"sample":'data',"color":ROOT.kBlack,"name":'Data'}]

    #13 TeV

    signal_qq_13 =  [{"sample":'ZZTo4l',"color":ROOT.kAzure-4,"name":'qq/qg #rightarrow ZZ(+jets)'}]
 
    signal_gg_13 = [{"sample":'ggZZ4e',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)'},{"sample":'ggZZ2e2mu',"color":ROOT.kAzure-5,"name":'ZZ processes'},{"sample":'ggZZ4mu',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)'}]

    signal_others_13 =  [{"sample":'ZH125',"color":ROOT.kAzure-6,"name":'other ZZ processes'},{"sample":'VBFH125',"color":ROOT.kAzure-6,"name":'other ZZ processes'},{"sample":'gH125',"color":ROOT.kAzure-6,"name":'other ZZ processes'},{"sample":'WplusH125',"color":ROOT.kAzure-6,"name":'other ZZ processes'},{"sample":'WminusH125',"color":ROOT.kAzure-6,"name":'other ZZ processes'}]

    bkg_13 = [{"sample":'TTTo2L2Nu2B',"color":ROOT.kOrange-3,"name":'Reducible background'},{"sample":'DYJetsToLLTuneCM50',"color":ROOT.kOrange-3,"name":'Reducible background'}]

    signal_tot_13 = signal_qq_13 +  signal_gg_13


    typeofsamples = []

    if category == 'Sig':
        typeofsamples = signal_tot 
    elif category == 'All':
        typeofsamples = signal_tot
    elif category == "RedBkg":
       typeofsamples =  background_red
    elif category == 'IrrBkg' or category == 'Bkg' :
        typeofsamples = irred_background_tot
    elif category== "AllIs":
         typeofsamples = irred_background_tot+signal_tot
    elif category == "data":
       typeofsamples = data

    #13TeV

    elif category == 'Sig13TeV':
        typeofsamples = signal_tot_13+signal_others_13
    elif category == 'All13TeV':
        typeofsamples = bkg_13+signal_tot_13+signal_others_13
    elif category == 'Bkg13TeV':
        typeofsamples = bkg_13

    else: sys.exit("ERROR, check Category") 

    return typeofsamples
    

####### Extract MC plot #########
def GetMCPlot(inputdir, category, plot,Addfake):
    print Red("\n#########################################\n############## Monte Carlo ##############\n#########################################\n")
    print category
    leg = TLegend(0.61,0.56,0.85,0.81)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.025)
 
    typeofsamples = GetTypeofsamples(category,"pow")
    files = {}
    stack = ROOT.THStack("stack",plot+"_stack")
   
    if Addfake:
        hfake = GetFakeRate(inputdir.replace("SR/",""),plot,"data") 
        stack.Add(hfake)
        leg.AddEntry(hfake,"Reducible background","f")
     
    LastColor = ROOT.kBlack

    for sample in typeofsamples:
        files[sample["sample"]] = ROOT.TFile(inputdir+sample["sample"]+".root")
        totalMC = 0

    print Red("\n######### Contribution of every sample #########\n")

    for sample in typeofsamples:

        h = files[sample["sample"]].Get(plot)
        print plot
        if not h:
            print sample["sample"],"has no enetries or is a zombie"
            continue
        print sample["sample"],"..........................", h.Integral()

        totalMC += h.Integral()
        if "Mass" in plot:        h.Rebin(2)
        h.SetLineColor(sample["color"])
        h.SetFillColor(sample["color"])
        h.SetMarkerStyle(21)
        h.SetMarkerColor(sample["color"])
        stack.Add(h)
        if LastColor!=sample["color"]:
            leg.AddEntry(h,sample["name"],"f")
        LastColor=sample["color"]
    print "\n Total MC .......................... {0:.2f}".format(totalMC)
    print "____________________________________ "       
    return (copy.deepcopy(stack),copy.deepcopy(leg))


#################################################

def GetDataPlot(inputdir, plot, Region):
    print "\n",""
    print Red("\n############################################\n################### DATA ###################\n############################################\n")
    files = {}
    typeofsamples = GetTypeofsamples("data","pow")
    hdata=ROOT.TH1F()
    plot = plot.replace("_JERCentralSmear","")
    for sample in typeofsamples:
        files[sample["sample"]] = ROOT.TFile(inputdir+sample["sample"]+".root")
    
    isFirst=1
    for sample in typeofsamples:
        h = files[sample["sample"]].Get(plot)
        if not h:
            print sample['sample'],'has no entries or is a zombie'
            continue

        print sample["sample"], "..........................",h.Integral()       
        if isFirst:
            hdata=copy.deepcopy(h)
            isFirst=0
            continue
        hdata.Add(h)   
    if 1+inputdir.find("CR2P2F"): 
        print "Dir",inputdir,inputdir.find("CR2P2F")                
        hdata.Add(hdata,-2)     
    #fdata = ROOT.TFile(inputdir+"data.root")
    #hdata = fdata.Get(plot)
    hdata.SetMarkerColor(ROOT.kBlack)
    hdata.SetLineColor(ROOT.kBlack)
    hdata.SetMarkerStyle(20)
    hdata.SetMarkerSize(.8)
    if "Mass" in plot:    hdata.Rebin(2)
    print "Total data ..........................",hdata.Integral()
    print " __________________ "   
    DataGraph=SetError(hdata,Region)
    DataGraph.SetMarkerStyle(20)
    DataGraph.SetMarkerSize(.9)
    return copy.deepcopy(DataGraph)

###############################################################

def GetMCPlot_fstate(inputdir, category, plot,Addfake):
    print Red("\n#########################################\n############## Monte Carlo ##############\n#########################################\n")
   
    print Red("\n######### Contribution to Signal #########\n")    
    leg = TLegend(0.61,0.56,0.85,0.81)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.025)
    files={}
    filesbkg ={}

    if "13TeV" in category:
        bkgsamples = GetTypeofsamples("Bkg13TeV","pow")    
        typeofsamples = GetTypeofsamples("Sig13TeV","pow") 
    else: 
        bkgsamples = GetTypeofsamples("IrrBkg","pow")    
        typeofsamples = GetTypeofsamples("Sig","pow") 

    for s in typeofsamples:
        files[s["sample"]] = ROOT.TFile(inputdir+s["sample"]+".root")

    for b in bkgsamples:
        filesbkg[b["sample"]] = ROOT.TFile(inputdir+b["sample"]+".root")

    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()
    hsum4l    = ROOT.TH1F()

    Var  = plot.replace("ZZTo4l","")

    hsum = [{"state":hsum2e2mu,"color":ROOT.kAzure-4,"name":'2e2m'},{"state":hsum4e,"color":ROOT.kAzure-5,"name":'4e'},{"state":hsum4mu,"color":ROOT.kAzure-6,"name":'4m'},{"state":hsum4l,"color":ROOT.kAzure-6,"name":'4l'}]

    for h in hsum:
        print Blue("### "+h["name"]+" ###")
        NoSamples = "For "+h["name"]+" there are no events in "
        isFirst=1

        for s in typeofsamples:
            print "ZZTo"+h["name"]+Var
            hsamp = files[s["sample"]].Get("ZZTo"+h["name"]+Var)  
            if hsamp==None:
                NoSamples+=s["sample"]+" "
                continue
            if isFirst:
                h["state"]=copy.deepcopy(hsamp)            
                isFirst=0
                continue 

            ErrStat=ROOT.Double(0.)
            print "{0} contribution {1} {2:.3f} +- {3: .3f} \n".format(s["sample"],(40-len(s["sample"]))*" ",hsamp.IntegralAndError(1,-1,ErrStat),ErrStat)
            h["state"].Add(hsamp)        
        print NoSamples,"\n\n" 
        
    print Blue("### Signal ###")  
    for h in hsum:
        print ("Total contribution "+h["name"]+" {0} {1:.3f} +- {2: .3f} \n").format((32-len(h["name"]))*" ",h["state"].IntegralAndError(1,-1,ErrStat),ErrStat)

    print Red("\n######### Contribution to Irreducible Background#########\n")    


    bsum2e2mu = ROOT.TH1F()
    bsum4e    = ROOT.TH1F()
    bsum4mu   = ROOT.TH1F()
    bsum4l   = ROOT.TH1F()
   
    bsum = [{"state":bsum2e2mu,"color":ROOT.kAzure-4,"name":'2e2m'},{"state":bsum4e,"color":ROOT.kAzure-5,"name":'4e'},{"state":bsum4mu,"color":ROOT.kAzure-6,"name":'4m'},{"state":bsum4l,"color":ROOT.kAzure-6,"name":'4l'}]

    for hbkg in bsum:
        print Blue("### "+hbkg["name"]+" ###")
        NoSamples = "For "+hbkg["name"]+" there are no events in "
        isFirst=1
        for b in bkgsamples:
            hb = filesbkg[b["sample"]].Get("ZZTo"+hbkg["name"]+Var)  
            if hb==None:
                print "For sample ", b["sample"], "has no enetries or is a zombie"       
                NoSamples+=b["sample"]+" "
                continue
            if isFirst:
                hbkg["state"]=copy.deepcopy(hb)            
                print "{0} contribution {1} {2:.3f} +- {3: .3f} \n".format(b["sample"],(40-len(b["sample"]))*" ",hb.IntegralAndError(1,-1,ErrStat),ErrStat)
                isFirst=0
                continue 

            ErrStat=ROOT.Double(0.)
            print "{0} contribution {1} {2:.3f} +- {3: .3f} \n".format(b["sample"],(40-len(b["sample"]))*" ",hb.IntegralAndError(1,-1,ErrStat),ErrStat)
            hbkg["state"].Add(hb)        
        print NoSamples,"\n\n" 

    for hbkg in bsum:
        print ("Total contribution "+hbkg["name"]+" {0} {1:.3f} +- {2: .3f} \n").format((32-len(hbkg["name"]))*" ",hbkg["state"].IntegralAndError(1,-1,ErrStat),ErrStat)
        if "Mass" in plot:    hbkg["state"].Rebin(2)
       
    stack = ROOT.THStack("stack",plot+"_stack")   
    stack.Add(bsum[3]["state"])
    leg.AddEntry(bsum[3]["state"], "Irreducible background","f")
    bsum[3]["state"].SetLineColor(b["color"])
    bsum[3]["state"].SetFillColor(b["color"])           
  
    if Addfake:
        print Red("\n######### Contribution to Reducible Background#########\n")    
        for i in ["2e2m","4e","4m","4l"]:
            print Blue("### "+i+" ###")
            hfake = GetFakeRate(inputdir.replace("SR/",""),"ZZTo"+i+Var,"data") 
            if i=="4l":
                stack.Add(hfake)
                leg.AddEntry(hfake,"Reducible background","f")
        
    # if category!="FsSig" and category!="Sig":
    #     stack.Add(hIrredBkg)
    #     leg.AddEntry(hIrredBkg, "Irreducible background","f")

    print Red("\n######### Signal samples for every final state #########\n")

    LastColor = ROOT.kBlack
    for i in hsum:
        if i["name"]=="4l": continue
        if i["state"]==None:
            print i["state"]," has no enetries" 
            continue
  
        i["state"].SetLineColor(i["color"])
        i["state"].SetFillColor(i["color"])
        if "Mass" in plot:        i["state"].Rebin(2) # FIX MEEEEE
       
        stack.Add(i["state"])
       
        if LastColor!=i["color"]:
            leg.AddEntry(i["state"],i["name"],"f")
        LastColor=i["color"]
    
        
    return (copy.deepcopy(stack),copy.deepcopy(leg))


###############################################################


def GetFakeRate(inputdir,plot, method):
#    hData = GetFakeRate("results/ZZjAnalyzer_","ZZTo"+FinState+Var,"data")            

    if method=="MC":
        print " Non MC method yet"
        #fileFake = ROOT.TFile(inputdir+"data.root")
        return 0.
    else:
        fileFake = ROOT.TFile(inputdir+"CR/data.root")
    hFakeRate=ROOT.TH1F()
    plot = plot.replace("_JERCentralSmear","")

    hFakeRate = fileFake.Get(plot)
    Err=ROOT.Double(0.)
    Integr= hFakeRate.IntegralAndError(0,-1,Err)
    if "Mass" in plot:    hFakeRate.Rebin(2)
    hFakeRate.SetFillColor(ROOT.kGray)
    hFakeRate.SetLineColor(ROOT.kBlack)
    hFakeRate.SetMarkerStyle(21)
    hFakeRate.SetMarkerSize(.5)


    #    print "\n","Total integral",FinState,"contribution ----------> ",Integr," +- ",Err,"\n"
    print "\n","Total integral contribution ----------> ",Integr," +- ",Err,"\n"
    return  copy.deepcopy(hFakeRate) 

########################################################

def GetSignalDefPlot(inputdir,category):

    leg = TLegend(0.61,0.56,0.85,0.81)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.025)

    typeofsamples=GetTypeofsamples(category)
    files = {}

    for s in typeofsamples:
        files[s["sample"]] = ROOT.TFile(inputdir+s["sample"]+".root")
    
    stack = ROOT.THStack("stack","SR: signal definition for "+category)
    hSig = ROOT.TH1F()
    hNoSig = ROOT.TH1F()

    isFirst =1
    for s in typeofsamples:
        hs = files[s["sample"]].Get("PassDef")
        if hs==None:
            print "For sample ", s["sample"],"PassDef has no entries"
            continue
        if isFirst:
            hSig = copy.deepcopy(hs)
            isFirst=0
            continue
        
        hSig.Add(hs)
    print "Total events passing signal defition ", hSig.Integral()

    isFirst =1
    for s in typeofsamples:
        hn = files[s["sample"]].Get("NoPassDef")
        if hn==None:
            print "For sample ", s["sample"],"NoPassDef has no entries"
            continue
        if isFirst:
            hNoSig = copy.deepcopy(hn)
            isFirst=0
            continue
        
        hNoSig.Add(hn)
    print "Total events not passing signal defition ", hNoSig.Integral()
       
    hSig.SetFillColor(ROOT.kAzure-4)
    hSig.SetLineColor(ROOT.kAzure-4)
    hNoSig.SetFillColor(ROOT.kAzure-6)
    hNoSig.SetLineColor(ROOT.kAzure-6)      
    if "Mass" in plot:    hSig.Rebin(2)
    if "Mass" in plot:    hNoSig.Rebin(2)
    stack.Add(hSig)
    stack.Add(hNoSig)
    leg.AddEntry(hSig,"Pass S def","f")
    leg.AddEntry(hNoSig,"No Pass S def","f")

    return (copy.deepcopy(stack),copy.deepcopy(leg))

#######################################################

def SetError(Histo,Region):
    q=(1-0.6827)/2.
    Graph=ROOT.TGraphAsymmErrors(Histo)
    Nbins= Histo.GetNbinsX()
    for i in range(1,Nbins):
       # if Region=="CR3P1F" or Region=="CR2P2F":
        
        N=Histo.GetBinContent(i)
        if N==0 : statMin=0
        else: statMin = (N-ROOT.Math.chisquared_quantile_c(1-q,2*N)/2.)
        statPlus = ROOT.Math.chisquared_quantile_c(q,2*(N+1))/2.-N
        Graph.SetPointError(i-1,0,0,statMin,statPlus)
    return Graph
