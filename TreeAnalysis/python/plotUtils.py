#!/usr/bin/env python2

import ROOT, copy, sys, os
from copy import deepcopy
import math
from readSampleInfo import *
from collections import OrderedDict
from Colours import *
import ctypes
from plotUtils23 import TFileContext, addIfExisting, PlotNotFoundError, InputDir, InputFile, set_overflow_range
import samplesByRegion # getSamplesByRegion, data_obs, ZZG, WZG, ...

##############################################
# Utilities, function definitions, and stuff #
##############################################

def getPlot_impl(filename, var):
    try:
        with TFileContext(filename) as tf:
            retrieved = deepcopy(tf.Get(var))
    except (IOError, OSError) as e:
        # ROOT prints an error message on its own
        retrieved = None
    else:
        if(not retrieved):
            retrieved = None
        if(retrieved is None):
            print('WARN: plot "{}" not present in file "{}"'.format(var, filename))
    return retrieved

def getPlot(plot, sample, region, inputdir='results', year='2016', analyzer='VVGammaAnalyzer'):
    theInputDir  = InputDir(inputdir, year=year, region=region, analyzer=analyzer)
    theInputFile = InputFile(theInputDir, '{:s}.root'.format(sample))
    if year == 'Run2':
        plots = []
        for y in ['2016preVFP', '2016postVFP', '2017', '2018']:
            filename = theInputFile.path(year=y)
            h = getPlot_impl(filename, plot)
            plots.append( h )
        return addIfExisting(*plots)
    else:
        filename = theInputFile.path()
        return getPlot_impl(filename, plot)

def iterate_bins(h, **kwargs):
    loX = 0 if kwargs.get("underX") else 1
    loY = 0 if kwargs.get("underY") else 1
    loZ = 0 if kwargs.get("underZ") else 1
    upX = 0 if kwargs.get("overfX") else 1
    upY = 0 if kwargs.get("overfY") else 1
    upZ = 0 if kwargs.get("overfZ") else 1
    for i in range(loX, h.GetNbinsX() + 1 + upX):
        for j in range(loY, h.GetNbinsY() + 1 + upY):
            for k in range(loZ, h.GetNbinsZ() + 1 + upZ):
                yield h.GetBin(i, j, k)

def is4Lregion(region):
    return region in ['SR4P', 'CR3P1F' , 'CR2P2F' , 'SR4P_1L', 'SR4P_1P', 'CR4P_1F', 'CR4L']
def is3Lregion(region):
    return region in ['SR3P', 'CR110'  , 'CR101'  , 'CR011'  , 'CR100'  , 'CR001'  , 'CR010', 'CR000', 'SR3P_1L', 'SR3P_1P', 'CR3P_1F', 'CRLFR', 'CR3L']
def is2Lregion(region):
    return region in ['SR2P', 'SR2P_1L', 'CR2P_1F']

###################################################
# Functions used to retrieve and manipulate plots #
###################################################

def GetTypeofsamples(category,Set):
    
    signal_qq_pow =  [{"sample":'ZZTo4l'        ,"color":ROOT.kAzure-4,"name":'qq/qg #rightarrow ZZ(+jets)'}]
    signal_qq_mad =  [{"sample":'ZZTo4lamcatnlo',"color":ROOT.kAzure-4,"name":'qq/qg #rightarrow ZZ(+jets)'}] 

    if   Set=="pow": signal_qq = signal_qq_pow
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
    elif category == "CR4L":
       typeofsamples = bkg_red+signal_tot
    else: sys.exit("ERROR, check Category") 

    return typeofsamples


def getPlotFromSample(inputdir, sample, plot, verbosity, forcePositive, note=None):
    _nameFormat = "{:24.24s}"
    errStat = ctypes.c_double(0.)
    totalIntegral = totalError = 0
    h = None

    if(inputdir.year == 'Run2'): years = ('2016preVFP', '2016postVFP', '2017', '2018')
    else:                        years = (inputdir.year,)
    multiyear = len(years) > 1
    isReversed = forcePositive and inputdir.region in ['CR2P2F','CR100','CR010','CR001']

    for fname in sample['files']:
        integralFile = errorFile = 0.
        for year in years:
            rootfilename = os.path.join(inputdir.path(year=year), fname+".root")
            fname_year = fname if not multiyear else fname+' '+year
            if(not os.path.exists(rootfilename)):
                if(verbosity >= 2):
                    print _nameFormat.format(fname_year), "No file" + ("" if(verbosity < 3) else " (%s)"%(rootfilename))
                continue

            with TFileContext(rootfilename) as fhandle:
                h_current = fhandle.Get(plot)

                if(not h_current):
                    if(verbosity >= 2):
                        print _nameFormat.format(fname_year), "No histo" + ("" if(verbosity < 3) else " (%s)"%(plot)) + " in file" + ("" if(verbosity < 4) else " (%s)"%(rootfilename))
                    continue

                if isReversed:
                    h_current.Scale(-1)

                integral = h_current.IntegralAndError(0, -1, errStat)  # Get overflow events too
                integralFile += integral
                errorFile = math.sqrt(errorFile**2 + errStat.value**2)

                if(h is None): h = copy.deepcopy(h_current)
                else         : h.Add(h_current)

        if(verbosity >= 2 and h is not None):
            if(note is not None): fname_print = fname + ' ' + note
            else:                 fname_print = fname
            print (_nameFormat+" {: 10.2f} +- {: 10.2f}").format(fname_print, integralFile, errorFile)
        totalIntegral += integralFile
        totalError    += errorFile

    return h, (totalIntegral, totalError)

###############################################################

def GetMCPlot_fstate(inputdir, category, plot,Addfake,MCSet,rebin):
    print Red("\n#########################################\n############## Monte Carlo ##############\n#########################################\n")
   
    print Red("\n######### Contribution to Signal #########\n")    
    leg = ROOT.TLegend(0.51,0.56,0.85,0.81)
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
                print "For sample ", b["sample"], "has no entries or is a zombie"       
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
            hfake = GetFakeRate(inputdir.replace("SR4P/",""), {'name':"ZZTo"+i+Var, 'rebin':rebin}, "data")
            if i=="4l":
                stack.Add(hfake)
                leg.AddEntry(hfake,"Reducible background","f")
    

    print Red("\n######### Signal samples for every final state #########\n")

    LastColor = ROOT.kBlack
    for i in hsum:
        if i["name"]=="4l": continue
        if i["state"]==None:
            print i["state"]," has no entries" 
            continue
  
        i["state"].SetLineColor(i["color"])
        i["state"].SetFillColor(i["color"])
        if rebin !=1:   i["state"].Rebin(rebin) # FIX MEEEEE
       
        stack.Add(i["state"])
       
        if LastColor!=i["color"]:
            leg.AddEntry(i["state"],i["name"],"f")
        LastColor=i["color"]
        
    return stack, leg

########################################################

def GetSignalDefPlot(inputdir,category):

    leg = ROOT.TLegend(0.61,0.56,0.82,0.81)
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
    # See also https://twiki.cern.ch/twiki/bin/viewauth/CMS/PoissonErrorBars
    h_copy = Histo.Clone(Histo.GetName()+'_copy')
    h_copy.SetBinErrorOption(ROOT.TH1.kPoisson)
    return ROOT.TGraphAsymmErrors(h_copy)


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
#         ROOT.gStyle.SetFrameLineWidth(3)
#         ROOT.gStyle.SetLineWidth(3)
