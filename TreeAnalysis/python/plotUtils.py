#! /usr/bin/env python

import ROOT, copy, sys, os
from copy import deepcopy
from errno import EEXIST
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

# Emulate os.makekdirs(..., exists_ok=True) for python2
def makedirs_ok(*args):
    try: os.makedirs(*args)
    except OSError as e:
        if(e.errno != EEXIST): raise e  # Catch only "File esists"

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


####### Extract predictions plot #########
def GetPredictionsPlot(inputdir, plotInfo, predType, MCSet, forcePositive=False, verbosity=1):
    region = inputdir.region
    plot = plotInfo['name']
    rebin = plotInfo.get('rebin', 1)
    overflow  = plotInfo.get('draw_overflow' , False)
    underflow = plotInfo.get('draw_underflow', False)

    controlRegions = []
    if region == 'SR4P':
        controlRegions = ['CR2P2F', 'CR3P1F']
    elif region == 'SR3P':
        controlRegions = ['CR000', 'CR001', 'CR010', 'CR011', 'CR100', 'CR101', 'CR110']
    elif region in ['CR001', 'CR010', 'CR100']:
        controlRegions = ['CR000']
    elif region == 'CR011': controlRegions = ['CR000', 'CR001', 'CR010']
    elif region == 'CR101': controlRegions = ['CR000', 'CR001', 'CR100']
    elif region == 'CR110': controlRegions = ['CR000', 'CR010', 'CR100']
    else:
        if(predType in ['fromCR', 'fakeMC']):
            print 'WARN: no rule for fake-lepton bkg for region "{}"'.format(region)  # "You should know what you are doing"

    useFakeLeptonsFromData = predType in ('fromCR', 'lepCR', 'fullCR')
    useFakePhotonsFromData = predType in ('fullCR', 'phoCR') and plotInfo.get('fake_photons') is not None
    
    if(verbosity == 1):
        print Red("\n############## "+    plot     +" ##############")
    elif(verbosity >= 2):
        print Red("\n###############"+'#'*len(plot)+"###############"
                  "\n############## "+    plot     +" ##############"
                  "\n###############"+'#'*len(plot)+"###############")

    leg = ROOT.TLegend(0.5,0.52,0.79,0.87)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.025)
    leg.SetFillStyle(0)

    samples = samplesByRegion.getSamplesByRegion(region, MCSet, predType)

    stack = ROOT.THStack("stack",plot+"_stack")
    ErrStat = ctypes.c_double(0.)

    if useFakeLeptonsFromData:
        if(verbosity >= 1):
            print Green("\nNon-prompt leptons background")
        hfakes = []
        for CR in controlRegions:
            newdir = copy.deepcopy(inputdir)
            newdir.region = CR
            hfakes.append( GetFakeRate(newdir, plotInfo, "data", MCSet, verbosity=verbosity) )
        hfake = addIfExisting(*hfakes)
        if(hfake is None):
            raise PlotNotFoundError('Fake lepton plot not found for ' + plotInfo['name'])

        hfake.SetLineColor(ROOT.kBlack)
        set_overflow_range(hfake, underflow=underflow, overflow=overflow)
        stack.Add(hfake)
        leg.AddEntry(hfake,"Non-prompt leptons","f")

    elif predType == 'fakeMC':  # Hack: use MCs in CRs as if they were data
        if(verbosity >= 1):
            print Green('\nNon-prompt leptons from MC in control regions')
        hfake = None
        newdir = copy.deepcopy(inputdir)
        for controlRegion in controlRegions:
            newdir.region = controlRegion
            hfakeTmp, _ = GetPredictionsPlot(newdir, plotInfo, 'fullMC', MCSet, forcePositive=forcePositive)
            if hfakeTmp is None: continue
            if not hfakeTmp.GetStack():
                print("WARN: fakeMC stack is null!")
                continue
            if hfakeTmp.GetStack().GetEntries() == 0:
                print("WARN: got 0 predictions from fakeMC")
                continue
            if hfake is None:
                hfake = copy.deepcopy(hfakeTmp.GetStack().Last())
            else:
                hfake.Add(hfakeTmp.GetStack().Last())
        hfake.SetLineColor(ROOT.kBlack)
        stack.Add(hfake)
        leg.AddEntry(hfake,"Non-prompt lept (MC)","f")

    if useFakePhotonsFromData:
        if(verbosity >= 1):
            print Green("\nNon-prompt photons background")
        fakeName = plotInfo['fake_photons']
        hfakePho, (integral, _) = getPlotFromSample(inputdir, samplesByRegion.data_obs, fakeName, verbosity, forcePositive)
        if(not hfakePho):
            raise PlotNotFoundError('Missing non-prompt photon plot {} (in {})'.format(fakeName, inputdir.path()))
        hfakePho.SetLineColor(ROOT.kBlack)
        hfakePho.SetFillColor(ROOT.kGreen-8)
        set_overflow_range(hfakePho, underflow=underflow, overflow=overflow)
        stack.Add(hfakePho)
        leg.AddEntry(hfakePho, "Non-prompt photons", "f")
    
    totalMC = 0
    
    if(verbosity >= 1):
        print Red("\n######### Contribution to {0:s}  #########\n".format(region))
    
    for sample in samples:
        h = None
        do_prompt_ph    = not sample.get('skip_prompt_ph'   , False)
        do_nonprompt_ph = not sample.get('skip_nonprompt_ph', False)
        splitPromptPh = (sample.get('split_prompt_ph') or not do_prompt_ph or not do_nonprompt_ph) and plotInfo.get('split_prompt_ph')

        if(splitPromptPh):
            split_pattern = plotInfo.get('split_prompt_ph_pattern', plot+'_%s')

            if(do_prompt_ph):
                h_prompt, (integralPrompt, _) = getPlotFromSample(inputdir, sample, split_pattern % ('prompt'), verbosity, forcePositive, note='prompt')
            else:
                h_prompt, integralPrompt = None, 0

            if(do_nonprompt_ph and not useFakePhotonsFromData):
                h_nonpro, (integralNonpro, _) = getPlotFromSample(inputdir, sample, split_pattern % ('nonpro'), verbosity, forcePositive, note='nonpro')
            else:
                h_nonpro, integralNonpro = None, 0
            totalMC += integralPrompt + integralNonpro

            for h in [h_prompt, h_nonpro]:
                if(h is None):
                    continue
                h.Scale(sample.get("kfactor", 1.))
                if rebin!=1: h.Rebin(rebin)
                set_overflow_range(h, underflow=underflow, overflow=overflow)

                h.SetLineColor(ROOT.kBlack)
                h.SetFillColor(sample["color"])
                h.SetMarkerStyle(21)

            if(h_nonpro):
                if(do_prompt_ph):  # Change color only if both prompt and nonprompt are present
                    h_nonpro.SetFillStyle(3002)
                leg.AddEntry(h_nonpro, sample["name"]+' nonprompt', "f")
                stack.Add(h_nonpro)
            if(h_prompt):
                leg.AddEntry(h_prompt, sample["name"]+' prompt', "f")
                stack.Add(h_prompt)

        else:
            h, (integral, _) = getPlotFromSample(inputdir, sample, plot, verbosity, forcePositive)
            totalMC += integral

            if(h is None):
                continue

            h.Scale(sample.get("kfactor", 1.))
            if rebin!=1: h.Rebin(rebin)
            set_overflow_range(h, underflow=underflow, overflow=overflow)

            h.SetLineColor(ROOT.kBlack)
            leg.AddEntry(h,sample["name"],"f")

            h.SetFillColor(sample["color"])
            h.SetMarkerStyle(21)
            h.SetMarkerColor(sample["color"])

            stack.Add(h)
    
    if(verbosity >= 1):
        print "\n Total MC .......................... {0:.2f}".format(totalMC)
        print "____________________________________ "       
    return stack, leg


#################################################
def GetClosureStack(region, inputDir, plotInfo, forcePositive=False, verbosity=1):
    plot  = plotInfo['name']

    if  (verbosity >= 1):
        print Red("\n############## "+    plot     +" ##############")
    elif(verbosity >= 2):
        print Red("\n###############"+'#'*len(plot)+"###############"
                  "\n############## "+    plot     +" ##############"
                  "\n###############"+'#'*len(plot)+"###############")
    leg = ROOT.TLegend(0.6,0.52,0.79,0.87)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.025)
    leg.SetFillStyle(0)

    stack = ROOT.THStack("stack",plot+"_stack")
    
    if   region == 'SR4P':
        samples_prompt = samplesByRegion.ZZG
    elif region in ['SR3P', 'CR3P1F']:
        samples_prompt = samplesByRegion.ZZG + samplesByRegion.WZG
    elif region in ['SR2P', 'CR2P2F', 'CR110', 'CR101', 'CR011']:
        samples_prompt = samplesByRegion.ZG
    
    isReversed = forcePositive and region in ['CR2P2F','CR100','CR010','CR001']

    for sample in samples_prompt:
        sample.update({'title': sample['name' ]   })  # TEMP, must change convention also in samplesByRegion
        sample.update({'name' : sample['files'][0]})
    # plot ~ PhFRClosure_KtoVL_reweighted_mZZG
    plot_reweight = plot.replace('PASS', 'reweighted')
    
    totalMC = 0
    ErrStat = ctypes.c_double(0.)

    with TFileContext(os.path.join(inputDir, 'data.root'), 'READ') as tf:
        hFakeData   = tf.Get(plot_reweight)
        hFakeData.SetDirectory(0)  # Prevent ROOT from deleting stuff under my nose
        if(isReversed): hFakeData.Scale(-1)
    if  (verbosity >= 2):
        print Red("\n######### Nonprompt photon background for {0:s}  #########\n".format(region))
        integral = hFakeData.IntegralAndError(0,-1,ErrStat)  # Get overflow events too
        print "{0:16.16} {1:.3f} +- {2: .3f}".format('data', integral, ErrStat.value)

    for sample_prompt in samples_prompt:
        with TFileContext(os.path.join(inputDir, sample_prompt['files'][0]+'.root'), 'READ') as tf:
            hFakePrompt = tf.Get(plot_reweight)
            hPrompt     = tf.Get(plot)
            hFakePrompt.SetDirectory(0)
            hPrompt    .SetDirectory(0)
            if(isReversed):
                hPrompt.Scale(-1)
                hFakePrompt.Scale(-1)
                pass
        if  (verbosity >= 2):
            integral = hFakePrompt.IntegralAndError(0,-1,ErrStat)  # Get overflow events too
            print "{0:16.16} {1:.3f} +- {2: .3f}".format(sample_prompt['files'][0], integral, ErrStat.value)

        hFakeData.Add(hFakePrompt, -1)  # subtract prompt contribution from "fail" region; it is already weighted by the FR
        sample_prompt.update({'hist': hPrompt})

    samples = samples_prompt + [{'name':'fake-photons', 'color':ROOT.kGray, 'title':'nonprompt #gamma', 'hist':hFakeData}]

    if  (verbosity >= 1):
        print Red("\n######### Contribution to {0:s}  #########\n".format(region))

    for sample in samples:
        h = sample['hist']
        if not h:
            if  (verbosity >= 2):
                print "{0:16.16s}".format(sample['name']), "No entries or is a zombie"
            continue
        
        #h.Scale(sample.get("kfactor", 1.))

        integral = h.IntegralAndError(0,-1,ErrStat)  # Get overflow events too
        if  (verbosity >= 2):
            print "{0:16.16} {1:.3f} +- {2: .3f}".format(sample['name'], integral, ErrStat.value)
        totalMC += integral

        h.Rebin(plotInfo.get('rebin', 1))

        h.SetLineColor  (ROOT.kBlack)  # h.SetLineColor(sample["color"])
        h.SetFillColor  (sample["color"])
        h.SetMarkerColor(sample["color"])
        h.SetMarkerStyle(21)
        stack.Add(h)
        leg.AddEntry(h, sample['title'], "f")

    # print '>>> legend', [(p.GetLabel(), p.GetObject().GetEntries()) for p in leg.GetListOfPrimitives()]
    # print '>>> stack ', [(p.GetName(), p.GetTitle(), p.GetEntries()) for p in stack.GetHists()]
    if  (verbosity >= 1):
        print "\n Total background .......................... {0:.2f}".format(totalMC)
        print "____________________________________ "       
    return stack, leg


#################################################

def GetDataPlot(inputdir, plotInfo, forcePositive=False, verbosity=1):
    plot = plotInfo['name']
    overflow  = plotInfo.get('draw_overflow' , False)
    underflow = plotInfo.get('draw_underflow', False)

    if  (verbosity >= 1):
        print Red("\n###################    DATA    ###################\n")
    sample = samplesByRegion.data_obs
    hdata = None
    
    isFirst=1
    for fname in sample['files']:
        h = getPlot(plot, fname, inputdir.region, inputdir.basedir, year=inputdir.year, analyzer=inputdir.analyzer)
        if(not h):
            continue

        if forcePositive and any(cr in inputdir.region for cr in ['CR2P2F','CR100','CR010','CR001']):
            h.Scale(-1)
        
        if not h:
            if  (verbosity >= 1):
                print fname, 'has no entries or is a zombie'
            continue
        
        if  (verbosity >= 1):
            print fname, "in", inputdir.region, "..........................",h.Integral(0,-1)
        if hdata is None:
            hdata = copy.deepcopy(h)
        else:
            hdata.Add(h)

    if(hdata is None):
        raise PlotNotFoundError('no data plot "{}" in {}'.format(plot, inputdir))
    hdata.SetMarkerColor(ROOT.kBlack)
    hdata.SetLineColor(ROOT.kBlack)
    hdata.SetMarkerStyle(20)
    hdata.SetMarkerSize(.8)

    hdata.Rebin(plotInfo.get('rebin', 1))
    set_overflow_range(hdata, underflow=underflow, overflow=overflow)

    if  (verbosity >= 1):
        print "Total data in {0:s} region .......................... {1:.2f}".format(inputdir.region, hdata.Integral(0,-1))
        print "_________________________ "   
    DataGraph=SetError(hdata, inputdir.region, False)
    DataGraph.SetMarkerStyle(20)
    DataGraph.SetMarkerSize(.9)
    return DataGraph, hdata

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


###############################################################


def GetFakeRate(inputdir, plotInfo, method, MCSet='mad', verbosity=1):
    plot = plotInfo['name']
    region = inputdir.region

    hFakeRate, (integral, error) = getPlotFromSample(inputdir, samplesByRegion.data_obs, plot, verbosity=verbosity, forcePositive=False, note=None)
    if(hFakeRate is None):
        return None

    hFakeRate.Rebin(plotInfo.get('rebin', 1))

    hFakeRate.SetFillColor(ROOT.kGray)
    hFakeRate.SetLineColor(ROOT.kGray)
    hFakeRate.SetMarkerStyle(21)
    hFakeRate.SetMarkerSize(.5)
    

    if method=="MC":  # MC subtraction of prompt processes from CRs
        assert region is not None, "Must provide the region to subtract prompt background"
        samples = samplesByRegion.getSamplesByRegion(region, MCSet, predType)

        listPromptMC = []
        for sample in samples:
            h, _ = getPlotFromSample(inputdir, sample, plot, verbosity=1, forcePositive=False, note=None)
            listPromptMC.append(h)
        hPromptMC = addIfExisting(*listPromptMC)

        if(hPromptMC is not None):
            hFakeRate.Add(hPromptMC, -1)

        Err2    = ctypes.c_double(0.)
        Integr2 = hFakeRate.IntegralAndError(0,-1,Err2)
        if(Integr2 * integral < 0):
            print "WARN: data-driven background changed sign after prompt MC subtraction!"
            print "      samples used:", [sample["name"] for sample in samples]
        print "data-promptMC ({:6.6s})\t {:.3f} +- {: .3f}".format(region, Integr2 , Err2.value)
    else:
        print "data ({:6.6s}) \t {:.3f} +- {: .3f}".format(region, integral, error)
    return hFakeRate

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
    q=(1-0.6827)/2.
    Graph=ROOT.TGraphAsymmErrors(Histo)

    # Add manually under-/overflow points if needed
    if(Histo.GetXaxis().GetLast() > Histo.GetXaxis().GetNbins()):
        n = Graph.GetN()
        xaxis = Histo.GetXaxis()
        x = xaxis.GetBinCenter (xaxis.GetLast())
        y = Histo.GetBinContent(xaxis.GetLast())
        ey = math.sqrt(y)
        Graph.AddPoint(x, y)
        Graph.SetPointEYhigh(n, ey)
        Graph.SetPointEYlow (n, ey)
    if(Histo.GetXaxis().GetFirst() < 1):
        n = Graph.GetN()
        xaxis = Histo.GetXaxis()
        x = xaxis.GetBinCenter (0)
        y = Histo.GetBinContent(0)
        ey = math.sqrt(y)
        Graph.AddPoint(x, y)
        Graph.SetPointEYhigh(n, ey)
        Graph.SetPointEYlow (n, ey)

    for i in range(Histo.GetXaxis().GetFirst(), Histo.GetXaxis().GetLast()+1):
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
#         ROOT.gStyle.SetFrameLineWidth(3)
#         ROOT.gStyle.SetLineWidth(3)
