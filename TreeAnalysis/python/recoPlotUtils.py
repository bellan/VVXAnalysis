#! /usr/bin/env python

import ROOT, copy, sys
from ROOT import gSystem, TCanvas, TH1,  TPad, gStyle, TLegend, THStack, TGraphAsymmErrors,Math
from readSampleInfo import *
from collections import OrderedDict
from Colours import *
import ctypes
import re

def GetTypeofsamples(category,Set):
    
    signal_qq_pow =  [{"sample":'ZZTo4l',"color":ROOT.kAzure-4,"name":'qq/qg #rightarrow ZZ(+jets)'}]
    signal_qq_mad =  [{"sample":'ZZTo4lamcatnlo',"color":ROOT.kAzure-4,"name":'qq/qg #rightarrow ZZ(+jets)'}] 

    if Set=="pow":   signal_qq = signal_qq_pow
    elif Set=="mad": signal_qq = signal_qq_mad
    else: sys.exit("Wrong Set, choose pow or mad")

    signal_gg = [{"sample":'ggZZ',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)'}]
    signal_VBS = [{"sample":'ZZ4lJJ',"color":ROOT.kAzure-6,"name":'VBS'}]
    signal_others =  [{"sample":'HZZ',"color":ROOT.kAzure-7,"name":'higgs'}]
    
    #signal_gg = [{"sample":'ggZZ4e',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)'},{"sample":'ggZZ2e2mu',"color":ROOT.kAzure-5,"name":'ZZ processes'},{"sample":'ggZZ4mu',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)'}]
#    signal_gg = [{"sample":'ggTo4e_Contin_MCFM701',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)'},{"sample":'ggTo2e2mu_Contin_MCFM701',"color":ROOT.kAzure-5,"name":'ZZ processes'},{"sample":'ggTo4mu_Contin_MCFM701',"color":ROOT.kAzure-5,"name":'gg #rightarrow ZZ(+jets)'}]
     #signal_VBS = [{"sample":'ZZTo4eJJ',"color":ROOT.kAzure-6,"name":'VBS'},{"sample":'ZZTo2e2muJJ',"color":ROOT.kAzure-6,"name":'VBS'},{"sample":'ZZTo4muJJ',"color":ROOT.kAzure-6,"name":'VBS'}]
#    signal_VBS = [{"sample":'ZZJJTo4L_ewk',"color":ROOT.kAzure-6,"name":'VBS'}]

#    signal_others =  [{"sample":'ZZZ',"color":ROOT.kAzure-7,"name":'other ZZ processes'},{"sample":'WZZ',"color":ROOT.kAzure-7,"name":'other ZZ processes'}]#,{"sample":'ggH125',"color":ROOT.kAzure-6,"name":'other ZZ processes'},{"sample":'VBFH125',"color":ROOT.kAzure-6,"name":'other ZZ processes'}]
    
    #signal_others =  [{"sample":'ggH',"color":ROOT.kAzure-7,"name":'higgs'}]

    bkg_red = [{"sample":'WZ',"color":ROOT.kRed+2,"name":'WZ'},{"sample":'TTJets',"color":ROOT.kRed-4,"name":'t#bar{t}'},{"sample":'TTWJets',"color":ROOT.kRed-4,"name":'TT'},{"sample":'TTGJets',"color":ROOT.kRed-4,"name":'TT'},{"sample":'WWW',"color":ROOT.kGreen-1,"name":'others'},{"sample":'DYJetsToLL_M50',"color":ROOT.kGreen-5,"name":'DY'}]

    bkg_irr = [{"sample":'WWZ',"color":ROOT.kOrange,"name":'WWZ'},{"sample":'TTZJets_M10_MLM',"color":ROOT.kOrange,"name":'TTZ'}]
    bkg_irr_divided = [{"sample":'WWZ',"color":ROOT.kOrange,"name":'WWZ'},{"sample":'TTZJets_M10_MLM',"color":ROOT.kOrange-5,"name":'TTZ'}]
    data = [{"sample":'data',"color":ROOT.kBlack,"name":'Data'}]

    signal_tot = signal_qq + signal_gg  + signal_VBS + signal_others

    typeofsamples = []

    if category == 'Sig':
        typeofsamples = signal_tot 
    elif category == 'All':
        typeofsamples = signal_tot 
    elif category == "RedBkg":
       typeofsamples =  bkg_red
    elif category == 'IrrBkg' or category == 'Bkg' :
        typeofsamples = bkg_irr_divided
    elif category== "AllIs":
         typeofsamples = signal_tot
    elif category == "data":
       typeofsamples = data
    elif category == "CR4P":
       typeofsamples = bkg_red+signal_tot
    else: sys.exit("ERROR, check Category") 

    return typeofsamples
    

####### Extract MC plot #########
def GetMCPlot(inputdir, category, plot,Addfake,MCSet,rebin):
    print Red("\n#########################################\n############## Monte Carlo ##############\n#########################################\n")
    print "Category",category
    print "plot",plot
    leg = TLegend(0.6,0.52,0.79,0.87)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.025)

    bkgsamples = GetTypeofsamples("IrrBkg",MCSet)  
    typeofsamples = GetTypeofsamples(category,MCSet)

    files = {}
    filesbkg = {}
    stack = ROOT.THStack("stack",plot+"_stack")
    ErrStat = ctypes.c_double(0.) #ROOT.Double(0.)

    if category != "RedBkg" and category != "IrrBkg":
        print Red("\n######### Contribution to Irreducible Background#########\n")    
        
        for b in bkgsamples:
            filesbkg[b["sample"]] = ROOT.TFile(inputdir+b["sample"]+".root")
            Var     = plot.replace("ZZTo4l","")
            bsum    = ROOT.TH1F()    
            print Blue("### "+plot+" ###")
            NoSamples = "For "+plot+" there are no events in "
            isFirst=1
        for b in bkgsamples:
            hb = filesbkg[b["sample"]].Get(plot)  
            if hb==None:
                print "For sample ", b["sample"], "has no enetries or is a zombie"       
                NoSamples+=b["sample"]+" "
                continue
            if isFirst:
                bsum=copy.deepcopy(hb)            
                print "{0} contribution {1} {2:.3f} +- {3: .3f} \n".format(b["sample"],(40-len(b["sample"]))*" ",hb.IntegralAndError(1,-1,ErrStat),ErrStat.value)
                isFirst=0
                continue 
            
            print "{0} contribution {1} {2:.3f} +- {3: .3f} \n".format(b["sample"],(40-len(b["sample"]))*" ",hb.IntegralAndError(1,-1,ErrStat),ErrStat.value)
            bsum.Add(hb)        
                
        if rebin != 1:  bsum.Rebin(rebin)
        #bsum.Rebin(rebin)
        bsum.SetName("hirr")
        stack.Add(bsum)
            
        leg.AddEntry(bsum, "Irreducible background","f")
        bsum.SetLineColor(b["color"])
        bsum.SetFillColor(b["color"])           
            

    if Addfake:
        print Red("\n######### Contribution to Reducible Background#########\n")    
        hfake = GetFakeRate(inputdir.replace("SR4P/",""),plot,"data",rebin) 
        stack.Add(hfake)
        leg.AddEntry(hfake,"Reducible background","f")
     
    LastColor = ROOT.kBlack

    for sample in typeofsamples:

        #if category=="RedBkg":   files[sample["sample"]] = ROOT.TFile(inputdir+"reducible_background_from_"+sample["sample"]+".root") /To take mc reducible bkg with the data drive method
        if category=="RedBkg": 
            inputdir = inputdir.replace("CR4P","SR4P")
            files[sample["sample"]] = ROOT.TFile(inputdir+sample["sample"]+".root")
        else:                    files[sample["sample"]] = ROOT.TFile(inputdir+sample["sample"]+".root")

        totalMC = 0

    print Red("\n######### Contribution to signal  #########\n")

    for sample in typeofsamples:

        h = files[sample["sample"]].Get(plot)
        if not h:
            print sample["sample"],"has no enetries or is a zombie"
            continue

        print "{0} contribution {1} {2:.3f} +- {3: .3f} \n".format(sample["sample"],(40-len(sample["sample"]))*" ",h.IntegralAndError(0,-1,ErrStat),ErrStat.value)

        totalMC += h.Integral()

        if rebin!=1: h.Rebin(rebin) 

        if "IrrBkg" in category:    h.SetLineColor(1)
        else:   h.SetLineColor(sample["color"])

        h.SetFillColor(sample["color"])
        h.SetMarkerStyle(21)
        h.SetMarkerColor(sample["color"])

        if "CR2P2F" in inputdir:
            h.Scale(-1)
        stack.Add(h)
        if LastColor!=sample["color"]:
            leg.AddEntry(h,sample["name"],"f")
        LastColor=sample["color"]

    print "\n Total MC .......................... {0:.2f}".format(totalMC)
    print "____________________________________ "       
    return (copy.deepcopy(stack),copy.deepcopy(leg))


#################################################

def GetDataPlot(inputdir, plot, Region,rebin):
    print "\n",""
    print Red("\n############################################\n################### DATA ###################\n############################################\n")
    files = {}
    typeofsamples = GetTypeofsamples("data","pow")
    hdata=ROOT.TH1F()
#    plot = plot.replace("_JERCentralSmear","") #FIXME
    for sample in typeofsamples:
        files[sample["sample"]] = ROOT.TFile(inputdir+sample["sample"]+".root")
    
    isFirst=1
    for sample in typeofsamples:
        h = files[sample["sample"]].Get(plot)
        if not h:
            print sample['sample'],'has no entries or is a zombie'
            continue

        print sample["sample"], "..........................",h.Integral(0,-1)       
        if isFirst:
            hdata=copy.deepcopy(h)
            isFirst=0
            continue
        hdata.Add(h)   
    if 1+inputdir.find("CR2P2F"): 
        print "Dir",inputdir,inputdir.find("CR2P2F")                
        hdata.Add(hdata,-2)     
    #hdata = fdata.Get(plot)
    hdata.SetMarkerColor(ROOT.kBlack)
    hdata.SetLineColor(ROOT.kBlack)
    hdata.SetMarkerStyle(20)
    hdata.SetMarkerSize(.8)

    if rebin!=1: hdata.Rebin(rebin) 

    print "Total data ..........................",hdata.Integral(0,-1)
    print "_________________________ "   
    DataGraph=SetError(hdata,Region,False)
    DataGraph.SetMarkerStyle(20)
    DataGraph.SetMarkerSize(.9)
    return copy.deepcopy(DataGraph),hdata

###############################################################

def GetMCPlot_fstate(inputdir, category, plot,Addfake,MCSet,rebin):
    print Red("\n#########################################\n############## Monte Carlo ##############\n#########################################\n")
   
    print Red("\n######### Contribution to Signal #########\n")    
    leg = TLegend(0.51,0.56,0.85,0.81)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.025)
    files={}
    filesbkg ={}

    bkgsamples = GetTypeofsamples("IrrBkg",MCSet)  
    typeofsamples = GetTypeofsamples("Sig",MCSet) 

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
        ErrStat=ctypes.c_double(0.)
        print "ZZTo"+h["name"]+Var
        for s in typeofsamples:
            hsamp = files[s["sample"]].Get("ZZTo"+h["name"]+Var)  
            if hsamp==None:
                NoSamples+=s["sample"]+" "
                continue
            if isFirst:
                h["state"]=copy.deepcopy(hsamp)            
                isFirst=0
                print "{0} {1} {2:.3f} +- {3: .3f}".format(s["sample"],(40-len(s["sample"]))*" ",hsamp.IntegralAndError(1,-1,ErrStat),ErrStat.value)
                continue 

            print "{0} {1} {2:.3f} +- {3: .3f}".format(s["sample"],(40-len(s["sample"]))*" ",hsamp.IntegralAndError(1,-1,ErrStat),ErrStat.value)
            h["state"].Add(hsamp)        
        print "\n",NoSamples,"\n\n" 
        
    print Blue("### Signal ###")  
    for h in hsum:
        print ("Total contribution "+h["name"]+" {0} {1:.3f} +- {2: .3f} \n").format((32-len(h["name"]))*" ",h["state"].IntegralAndError(1,-1,ErrStat),ErrStat.value)

    stack = ROOT.THStack("stack",plot+"_stack")   


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
            
            ErrStat=ctypes.c_double(0.)
            print "{0} contribution {1} {2:.3f} +- {3: .3f} \n".format(b["sample"],(40-len(b["sample"]))*" ",hb.IntegralAndError(1,-1,ErrStat),ErrStat.value)
            hbkg["state"].Add(hb)        
            print NoSamples,"\n\n" 
            
            for hbkg in bsum:
                print ("Total contribution "+hbkg["name"]+" {0} {1:.3f} +- {2: .3f} \n").format((32-len(hbkg["name"]))*" ",hbkg["state"].IntegralAndError(1,-1,ErrStat),ErrStat.value)
                
                if rebin != 1:  hbkg["state"].Rebin(rebin)
        

    stack.Add(bsum[3]["state"])
    leg.AddEntry(bsum[3]["state"], "Irreducible background","f")
    bsum[3]["state"].SetLineColor(b["color"])
    bsum[3]["state"].SetFillColor(b["color"])           
  
    if Addfake:
        print Red("\n######### Contribution to Reducible Background#########\n")    
        for i in ["2e2m","4e","4m","4l"]:
            print Blue("### "+i+" ###")
            hfake = GetFakeRate(inputdir.replace("SR4P/",""),"ZZTo"+i+Var,"data",rebin) 
            if i=="4l":
                stack.Add(hfake)
                leg.AddEntry(hfake,"Reducible background","f")
    

    print Red("\n######### Signal samples for every final state #########\n")

    LastColor = ROOT.kBlack
    for i in hsum:
        if i["name"]=="4l": continue
        if i["state"]==None:
            print i["state"]," has no enetries" 
            continue
  
        i["state"].SetLineColor(i["color"])
        i["state"].SetFillColor(i["color"])
        if rebin !=1:   i["state"].Rebin(rebin) # FIX MEEEEE
       
        stack.Add(i["state"])
       
        if LastColor!=i["color"]:
            leg.AddEntry(i["state"],i["name"],"f")
        LastColor=i["color"]
        
    return (copy.deepcopy(stack),copy.deepcopy(leg))


###############################################################


def GetFakeRate(inputdir,plot, method,rebin):


#    print Red("\n######### Reducible Background #########\n")


    if method=="MC":
        print " Non MC method yet"
        #fileFake = ROOT.TFile(inputdir+"data.root")
        return 0.
    else:
        fileFake = ROOT.TFile(inputdir+"CR4P/data.root")
    hFakeRate=ROOT.TH1F()
    #if "Jet" in plot:
     #   plot = plot.replace("_JERCentralSmear","")

    hFakeRate = fileFake.Get(plot)
    Err=ctypes.c_double(0.)
    Integr= hFakeRate.IntegralAndError(0,-1,Err)
    if rebin != 1: hFakeRate.Rebin(rebin)
    hFakeRate.SetFillColor(ROOT.kGray)
    hFakeRate.SetLineColor(ROOT.kGray)
    hFakeRate.SetMarkerStyle(21)
    hFakeRate.SetMarkerSize(.5)

    ErrStat = ctypes.c_double(0.)
    #    print "\n","Total integral",FinState,"contribution ----------> ",Integr," +- ",Err,"\n"
#    print "{0} contribution {1} {2:.3f} +- {3: .3f} \n".format(b["sample"],(40-len(b["sample"]))*" ",hb.IntegralAndError(1,-1,ErrStat),ErrStat)

    print "Total contribution {0} {1:.3f} +- {2: .3f} \n".format(38*" ",hFakeRate.IntegralAndError(1,-1,ErrStat),ErrStat.value)

    return  copy.deepcopy(hFakeRate) 

########################################################

def GetSignalDefPlot(inputdir,category):

    leg = TLegend(0.61,0.56,0.82,0.81)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.025)

    typeofsamples=GetTypeofsamples(category)
    files = {}

    for s in typeofsamples:
        files[s["sample"]] = ROOT.TFile(inputdir+s["sample"]+".root")
    
    stack = ROOT.THStack("stack","SR4P: signal definition for "+category)
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
    if rebin!=1:
        hSig.Rebin(rebin)
        hNoSig.Rebin(rebin)
    stack.Add(hSig)
    stack.Add(hNoSig)
    leg.AddEntry(hSig,"Pass S def","f")
    leg.AddEntry(hNoSig,"No Pass S def","f")

    return (copy.deepcopy(stack),copy.deepcopy(leg))

#######################################################

def SetError(Histo,Region,Set0Error):
    q=(1-0.6827)/2.
    Graph=ROOT.TGraphAsymmErrors(Histo)
    Nbins= Histo.GetNbinsX()
    for i in range(1,Nbins):
       # if Region=="CR3P1F" or Region=="CR2P2F":
        if False:
            N=Histo.GetBinContent(i)
            if N==0 : 
                statMin=0
                if Set0Error: statPlus = ROOT.Math.chisquared_quantile_c(q,2*(N+1))/2.-N
                else: statPlus = 0
            else: 
                statMin = (N-ROOT.Math.chisquared_quantile_c(1-q,2*N)/2.)
                statPlus = ROOT.Math.chisquared_quantile_c(q,2*(N+1))/2.-N
                Graph.SetPointError(i-1,0,0,statMin,statPlus)
    return Graph


def setCMSStyle(self, canvas, author='N. Woods', textRight=True, dataType='Preliminary Simulation', energy=13, intLumi=19710.):
    '''
    Set plotting defaults to something appropriate for CMS Analysis Notes
    intLumi is given in pb^-1 and converted to fb^-1, unless it is less than 1 fb^-1
    If intLumi is nonpositive, it is not printed
    '''
    # Make sure that if there's an exponent on the X axis, it's visible but not on top of the axis title
    self.fixXExponent(canvas)

    # Make sure that temperature plot scales don't run off the right side
    self.fixZScale(canvas)

    # Put "Preliminary" or similar on the plots
    if dataType:
        CMS_lumi.relPosX = 0.12
        CMS_lumi.extraText = dataType
    else:
        CMS_lumi.writeExtraText = False

    # Put sqrt(s) on plots
    try:
        energy = [int(energy)]
    except TypeError:
        assert isinstance(energy,list) and all(isinstance(e, int) for e in energy), \
            "Energy must be an integer or list of integers"

    try:
        intLumi = [float(intLumi)]
    except TypeError:
        assert isinstance(intLumi, list) and all(isinstance(e, float) for il in intLumi), \
            "Integrated Luminosity must be a float  or list of floats"
    assert len(intLumi) == len(energy), "Must have exactly one integrated luminosity per energy"

    iPeriod = 0
    for i, e in enumerate(energy):
        iL = intLumi[i]
        if iL > 0.:
            if iL >= 1000.:
                iL /= 1000. # convert to fb^-1
                unit = "fb^{-1}"
            else:
                unit = "pb^{-1}"
            iLStr = makeNumberPretty(iL, 1)
        else:
            iLStr = ""
            unit = ""

        if e == 13:
            iPeriod += 4
            CMS_lumi.lumi_13TeV = CMS_lumi.lumi_13TeV.replace("20.1","%s"%iLStr).replace("fb^{-1}", unit)
        elif energy == 8:
            iPeriod += 2
            CMS_lumi.lumi_8TeV = CMS_lumi.lumi_8TeV.replace("19.7","%.1f"%iLStr).replace("fb^{-1}", unit)
        if energy == 7:
            iPeriod += 1
            CMS_lumi.lumi_7TeV = CMS_lumi.lumi_7TeV.replace("5.1","%.1f"%iLStr).replace("fb^{-1}", unit)

    # Put "CMS preliminary simulation" or whatever above the left side of the plot
    iPos = 0

    # Draw all that stuff
    CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)

    # Put author name and "Preliminary Exam" and a box in the top right corner of the frame
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(kBlack)
    latex.SetTextFont(61)
    latex.SetTextSize(0.03)
    latex.SetTextAlign(12)
    latex.DrawLatex(0.01, 0.05, author)
#        latex.DrawLatex(0.01, 0.02, "U. Wisconsin Preliminary Exam")

#         # Make frame and tick marks thicker
#         gStyle.SetFrameLineWidth(3)
#         gStyle.SetLineWidth(3)
