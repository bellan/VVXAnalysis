from array import array
from math import sqrt
import sys
import os
import copy
import logging
import ctypes
import numpy as np
from Colours import Red, Green
import ROOT
import cmsstyle

import samplesByRegion
from utils23 import NONPROMPT_CAP

if(sys.version_info.major < 3):
    import errno
    class FileNotFoundError(OSError):
        pass

class TFileContext(object):
    def __init__(self, *args):
        self.tfile = ROOT.TFile(*args)
        if(not (self.tfile and self.tfile.IsOpen())):
            raise FileNotFoundError(args[0] if len(args) > 0 else '')

    def __enter__(self):
        return self.tfile

    def __exit__(self, exc_type, exc_value, traceback):
        self.tfile.Close()


class PlotNotFoundError(Exception):
    pass


class InputDir:
    def __init__(self, basedir, year, region, analyzer):
        self.basedir  = basedir
        self.year     = year
        self.region   = region
        self.analyzer = analyzer

    def path(self, **kwargs):
        '''
        Returns a string corresponding to the path.
        Accepts keyword args for quick substitutions
        '''
        return os.path.join(
            kwargs.get('basedir', self.basedir),
            kwargs.get('year'   , self.year   ),
            '{}_{}'.format(
                kwargs.get('analyzer', self.analyzer),
                kwargs.get('region'  , self.region  )
            )
        )

    def __str__(self):
        return self.path()

    def __repr__(self):
        return str(vars(self))


class InputFile:
    def __init__(self, inputDir, fname):
        self.inputDir = inputDir
        self.fname = fname

    def path(self, **kwargs):
        return os.path.join(
            self.inputDir.path(**kwargs),
            kwargs.get('fname', self.fname)
        )

    def __str__(self):
        return self.path()

    def __repr__(self):
        return str(vars(self))


def addIfExisting(*args):
    result = None
    for a in [ a for a in args if a is not None ]:
        if result is None:
            result = a
        else:
            result.Add(a)
    return result


def debug_hist(hist):
    string = ''
    for bx in range(1, hist.GetNbinsX()+1):
        for by in range(1, hist.GetNbinsY()+1):
            b = hist.GetBin(bx, by)
            val = hist.GetBinContent(b)
            string += ' '.join(('\tx: %d,  y: %d' %(bx, by)
                                , '  bin: %2d' %(b)
                                , ' [%.0f, %3.0f]' %(hist     .GetXaxis().GetBinLowEdge(bx), hist     .GetXaxis().GetBinUpEdge(bx))
                                , ' [%.2f, %.2f]'  %(hist     .GetYaxis().GetBinLowEdge(by), hist     .GetYaxis().GetBinUpEdge(by))
                                , ' value: %+6.3g' %(val)
                                ))+'\n'
    return string


def rebin2D(hist_orig, x_bins=None, y_bins=None, verbose=False):
    x_bins_orig = array('d', hist_orig.GetXaxis().GetXbins())
    y_bins_orig = array('d', hist_orig.GetYaxis().GetXbins())
    if(x_bins is None): x_bins = x_bins_orig
    if(y_bins is None): y_bins = y_bins_orig
    assert all(edge in x_bins_orig for edge in x_bins), 'x_bins must not contain edges which are not in the original histogram: '+str(x_bins_orig)
    assert all(edge in y_bins_orig for edge in y_bins), 'y_bins must not contain edges which are not in the original histogram: '+str(y_bins_orig)

    hist = ROOT.TH2F(hist_orig.GetName()+'_rebin', hist_orig.GetTitle(),
                     len(x_bins) - 1, x_bins,
                     len(y_bins) - 1, y_bins)

    if(verbose): print('*****', hist_orig.GetName(), '*****')
    for bx in range(1, hist_orig.GetNbinsX()+1):
        x = hist_orig.GetXaxis().GetBinCenter(bx)
        for by in range(1, hist_orig.GetNbinsY()+1):
            y = hist_orig.GetYaxis().GetBinCenter(by)
            b = hist_orig.GetBin(bx, by)
            val = hist_orig.GetBinContent(b)
            err = hist_orig.GetBinError(b)
            b_new = hist.FindFixBin(x, y)
            previous_val_new = hist.GetBinContent(b_new)
            val_new = hist.GetBinContent(b_new) + val
            err_new = sqrt( hist.GetBinError(b_new)**2 + err**2 )  # Manually track error because Fill() assumes that each call is a single entry --> err^2 = summ(w_i^2)
            hist.SetBinContent(b_new, val_new)
            hist.SetBinError  (b_new, err_new)
            if(verbose):
                print('\tx:', bx, ' y:', by
                      , '  bin: %2d' %(b)
                      , ' [%.0f, %3.0f]' %(hist_orig.GetXaxis().GetBinLowEdge(bx), hist_orig.GetXaxis().GetBinUpEdge(bx))
                      , ' [%.2f, %.2f]'  %(hist_orig.GetYaxis().GetBinLowEdge(by), hist_orig.GetYaxis().GetBinUpEdge(by))
                      , ' value: %+6.3g' %(val)
                      , ' - new'
                      , ' bin: %2d' %(b_new)
                      , ' value: %+6.3g + %+6.3g = %+6.3g' %(previous_val_new, val, hist.GetBinContent(b_new))
                      # , '+- %5.3g' %(err),
                      # , ' -  err_new: %5.3g' %(err_new)
                      )

    if(verbose):
        print('    * rebinned *')
        print(debug_hist(hist))
    return hist


def retrieve_bin_edges(axis):
    '''
    Return an array.array of nbins+1 bin edges both in case of fixed and variable bins
    '''
    if( len(axis.GetXbins()) > 0 ):
        # variable bin size
        edges = array('d', axis.GetXbins())
    else:
        # fixed bin size
        edges = array('d', np.linspace(axis.GetXmin(), axis.GetXmax(), axis.GetNbins()+1))
    return edges


def get_plots_singleyear(inputdir, sample, plots):
    fname = InputFile(inputdir, sample+".root").path()
    retrieved = []
    with TFileContext(fname) as rFile:
        for plot in plots:
            h = rFile.Get(plot)
            if(not h):
                logging.warning('Could not get "%s" from "%s"' % (plot, fname))
                retrieved.append(None)
            else:
                h.SetDirectory(0)
                retrieved.append(h)
                if(logging.getLogger().isEnabledFor(logging.DEBUG)):
                    error  = ctypes.c_double(0)
                    ndim = h.GetDimension()
                    if  (ndim == 1): integr = h.IntegralAndError(0, -1, error)
                    elif(ndim == 2): integr = h.IntegralAndError(0, -1, 0, -1, error)
                    elif(ndim == 3): integr = h.IntegralAndError(0, -1, 0, -1, 0, -1, error)
                    else:            integr = float('inf')  # py2 portability
                    logging.debug('{:s} - integral: {:7.1f} +- {:5.1f}'.format(plot, integr, error.value))
    return retrieved


def get_plots(inputdir, sample, plots):
    if inputdir.year == 'Run2':
        inputdir_copy = copy.deepcopy(inputdir)  # Avoid overwriting the original
        plot_matrix = []  # dimension 0 = year, dimension 1 = requested plot names
        for year in ('2016preVFP', '2016postVFP', '2017', '2018'):
            inputdir_copy.year = year
            plot_matrix.append( get_plots_singleyear(inputdir_copy, sample, plots) )

        # Sum along dimension 0
        plot_list = [addIfExisting(*l) for l in zip(*plot_matrix)]
        return plot_list
    else:
        return get_plots_singleyear(inputdir, sample, plots)


def set_overflow_range(h, underflow=False, overflow=True):
    bx_min = (0 if underflow else 1)
    bx_max = (1 if overflow  else 0) + h.GetXaxis().GetNbins()
    h.GetXaxis().SetRange(bx_min, bx_max)


def deduce_region_label(plotname, region):
    split = plotname.split('_')
    if  (any(part.startswith('loose', 'failReweight') for part in split)):
        return region+'_1P'
    elif(any(part.startswith('fail') for part in split)):
        return region+'_1F'
    elif(any(part.startswith('veryLoose') for part in split)):
        return region+'_1L'
    else:
        return region


# Stuff ported from plotUtils (py2-only) and made py3-ready
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
                    print(_nameFormat.format(fname_year) + " No file" + ("" if(verbosity < 3) else " (%s)"%(rootfilename)))
                continue

            with TFileContext(rootfilename) as fhandle:
                h_current = fhandle.Get(plot)

                if(not h_current):
                    if(verbosity >= 2):
                        print(_nameFormat.format(fname_year) + " No histo" + ("" if(verbosity < 3) else " (%s)"%(plot)) + " in file" + ("" if(verbosity < 4) else " (%s)"%(rootfilename)))
                    continue

                if isReversed:
                    h_current.Scale(-1)

                integral = h_current.IntegralAndError(0, -1, errStat)  # Get overflow events too
                integralFile += integral
                errorFile = sqrt(errorFile**2 + errStat.value**2)

                if(h is None): h = copy.deepcopy(h_current)
                else         : h.Add(h_current)

        if(verbosity >= 2 and h is not None):
            if(note is not None): fname_print = fname + ' ' + note
            else:                 fname_print = fname
            print ((_nameFormat+" {: 10.2f} +- {: 10.2f}").format(fname_print, integralFile, errorFile))
        totalIntegral += integralFile
        totalError    += errorFile

    return h, (totalIntegral, totalError)


def GetFakeRate(inputdir, plotInfo, method, MCSet='mad', verbosity=1):
    plot = plotInfo['name']
    region = inputdir.region

    hFakeRate, (integral, error) = getPlotFromSample(inputdir, samplesByRegion.data_obs, plot, verbosity=verbosity, forcePositive=False, note=None)
    if(hFakeRate is None):
        return None

    hFakeRate.Rebin(plotInfo.get('rebin', 1))

    hFakeRate.SetFillColor(samplesByRegion.fake_leptons['color'])
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
            print("WARN: data-driven background changed sign after prompt MC subtraction!")
            print("      samples used: {}".format([sample["name"] for sample in samples]))
        print("data-promptMC ({:6.6s})\t {:.3f} +- {: .3f}".format(region, Integr2 , Err2.value))
    else:
        print("data ({:6.6s}) \t {:.3f} +- {: .3f}".format(region, integral, error))
    return hFakeRate


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
            print('WARN: no rule for fake-lepton bkg for region "{}"'.format(region))  # "You should know what you are doing"

    useFakeLeptonsFromData = predType in ('fromCR', 'lepCR', 'fullCR')
    useFakePhotonsFromData = predType in ('fullCR', 'phoCR') and plotInfo.get('fake_photons') is not None
    
    if(verbosity == 1):
        print(Red("\n############## "+    plot     +" ##############"))
    elif(verbosity >= 2):
        print(Red("\n###############"+'#'*len(plot)+"###############"
                  "\n############## "+    plot     +" ##############"
                  "\n###############"+'#'*len(plot)+"###############"))

    leg = ROOT.TLegend(0.32,0.5,0.8,0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    samples = samplesByRegion.getSamplesByRegion(region, MCSet, predType, special=plotInfo.get('special', False))

    stack = ROOT.THStack("stack",plot+"_stack")
    ErrStat = ctypes.c_double(0.)

    if useFakeLeptonsFromData:
        if(verbosity >= 1):
            print(Green("\nNon-prompt leptons background"))
        hfakes = []
        for CR in controlRegions:
            newdir = copy.deepcopy(inputdir)
            newdir.region = CR

            # A "temporary" hack to get the correct fake_leptons plot for SYS_mZZGloose(-nonpro)_central
            # for the special "compare" plot (data-driven vs MC predictions in SR4P_1P)
            # In the MC the plot is SYS_mZZGloose-nonpro_central, in fake_leptons it's SYS_mZZGloose_central
            info_copy = copy.deepcopy(plotInfo)
            if(plotInfo.get('special')):
               info_copy['name'] = plotInfo['name'].replace('-nonpro', '')

            hfakes.append( GetFakeRate(newdir, info_copy, "data", MCSet, verbosity=verbosity) )
        hfake = addIfExisting(*hfakes)
        if(hfake is None):
            raise PlotNotFoundError('Fake lepton plot not found for ' + plotInfo['name'])

        hfake.SetLineColor(ROOT.kBlack)
        set_overflow_range(hfake, underflow=underflow, overflow=overflow)
        stack.Add(hfake)
        leg.AddEntry(hfake, NONPROMPT_CAP+" l", "f")

    elif predType == 'fakeMC':  # Hack: use MCs in CRs as if they were data
        if(verbosity >= 1):
            print(Green('\nNon-prompt leptons from MC in control regions'))
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
        leg.AddEntry(hfake, NONPROMPT_CAP+" l (MC)", "f")

    if useFakePhotonsFromData:
        if(verbosity >= 1):
            print(Green("\nNon-prompt photons background"))
        fakeName = plotInfo['fake_photons']
        hfakePho, (integral, _) = getPlotFromSample(inputdir, samplesByRegion.data_obs, fakeName, verbosity, forcePositive)
        if(not hfakePho):
            raise PlotNotFoundError('Missing non-prompt photon plot {} (in {})'.format(fakeName, inputdir.path()))
        hfakePho.SetLineColor(ROOT.kBlack)
        hfakePho.SetFillColor(samplesByRegion.fake_photons['color'])
        set_overflow_range(hfakePho, underflow=underflow, overflow=overflow)
        stack.Add(hfakePho)
        leg.AddEntry(hfakePho, NONPROMPT_CAP+" #gamma", "f")

    totalMC = 0
    totalMCerr = 0
    
    if(verbosity >= 1):
        print(Red("\n######### Contribution to {0:s}  #########\n".format(region)))
    
    for sample in samples:
        h = None
        do_prompt_ph    = not sample.get('skip_prompt_ph'   , False)
        do_nonprompt_ph = not sample.get('skip_nonprompt_ph', False)
        splitPromptPh = (sample.get('split_prompt_ph') or not do_prompt_ph or not do_nonprompt_ph) and plotInfo.get('split_prompt_ph')

        if(splitPromptPh):
            split_pattern = plotInfo.get('split_prompt_ph_pattern', plot+'_%s')

            if(do_prompt_ph):
                h_prompt, (integralPrompt, errorPrompt) = getPlotFromSample(inputdir, sample, split_pattern % ('prompt'), verbosity, forcePositive, note='prompt')
            else:
                h_prompt, integralPrompt, errorPrompt = None, 0, 0

            if(do_nonprompt_ph and not useFakePhotonsFromData):
                h_nonpro, (integralNonpro, errorNonpro) = getPlotFromSample(inputdir, sample, split_pattern % ('nonpro'), verbosity, forcePositive, note='nonpro')
            else:
                h_nonpro, integralNonpro, errorNonpro = None, 0, 0
            totalMC += integralPrompt + integralNonpro
            totalMCerr = sqrt(totalMCerr**2 + errorPrompt**2 + errorNonpro**2)

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
                leg.AddEntry(h_nonpro, sample["name"]+' non-prompt', "f")
                stack.Add(h_nonpro)
            if(h_prompt):
                leg.AddEntry(h_prompt, sample["name"]+' prompt', "f")
                stack.Add(h_prompt)

        else:
            h, (integral, error) = getPlotFromSample(inputdir, sample, plot, verbosity, forcePositive)
            totalMC += integral
            totalMCerr = sqrt(totalMCerr**2 + error**2)

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
        print("\n Total MC .......................... {0:.2f} +- {1:.2f}".format(totalMC, totalMCerr))
        print("____________________________________")
    return stack, leg


def GetClosureStack(region, inputDir, plotInfo, forcePositive=False, verbosity=1):
    '''
    Compare simulation (stack) and data-driven (points) fake photon predictions
    In order to use a cmsstyle.cmsLeg, this should be a python3-only function (no cmsstyle for python2)
    '''
    plot  = plotInfo['name']

    if  (verbosity >= 1):
        print(Red("\n############## "+    plot     +" ##############"))
    elif(verbosity >= 2):
        print(Red("\n###############"+'#'*len(plot)+"###############"
                  "\n############## "+    plot     +" ##############"
                  "\n###############"+'#'*len(plot)+"###############"))
    leg = ROOT.TLegend(0.6,0.52,0.79,0.87, "", "brNDC")
    leg.SetTextSize(0.03)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)

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
        print(Red("\n######### Nonprompt photon background for {0:s}  #########".format(region)))
        integral = hFakeData.IntegralAndError(0,-1,ErrStat)  # Get overflow events too
        print("{0:16.16} {1:.3f} +- {2: .3f}".format('data', integral, ErrStat.value))

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
            print("{0:16.16} {1:.3f} +- {2: .3f}".format(sample_prompt['files'][0], integral, ErrStat.value))

        hFakeData.Add(hFakePrompt, -1)  # subtract prompt contribution from "fail" region; it is already weighted by the FR
        sample_prompt.update({'hist': hPrompt})

    samples = samples_prompt + [{'name':'fake-photons', 'color':ROOT.kGray, 'title':'non-prompt #gamma', 'hist':hFakeData}]

    if  (verbosity >= 1):
        print(Red("\n######### Contribution to {0:s}  #########".format(region)))

    for sample in samples:
        h = sample['hist']
        if not h:
            if  (verbosity >= 2):
                print("{0:16.16s} No entries or is a zombie".format(sample['name']))
            continue
        
        h.Scale(sample.get("kfactor", 1.))

        integral = h.IntegralAndError(0,-1,ErrStat)  # Get overflow events too
        if  (verbosity >= 2):
            print("{0:16.16} {1:.3f} +- {2: .3f}".format(sample['name'], integral, ErrStat.value))
        totalMC += integral

        h.Rebin(plotInfo.get('rebin', 1))

        h.SetLineColor  (ROOT.kBlack)  # h.SetLineColor(sample["color"])
        h.SetFillColor  (sample["color"])
        h.SetMarkerColor(sample["color"])
        h.SetMarkerStyle(21)
        stack.Add(h)
        leg.AddEntry(h, sample['title'], "f")

    if  (verbosity >= 1):
        print("\n Total background .......................... {0:.2f}".format(totalMC))
        print("____________________________________")
    return stack, leg


def SetError(Histo,Region,Set0Error):
    '''Creates a TGraphAsymmErrors from a TH1, with the appropriate error bars'''
    # See also https://twiki.cern.ch/twiki/bin/viewauth/CMS/PoissonErrorBars
    h_copy = Histo.Clone(Histo.GetName()+'_copy')
    h_copy.SetBinErrorOption(ROOT.TH1.kPoisson)
    tga = ROOT.TGraphAsymmErrors(h_copy)

    return tga


def GetDataPlot(inputdir, plotInfo, forcePositive=False, verbosity=1):
    '''
    Retrieve the data histogram (plotInfo) from the results in inputdir.
    Returns a TGraphAsymmErrors (to be drawn with the MC stack) and the
    original TH1 (to be used e.g. to calculate a data/MC ratio).
    '''
    plot = plotInfo['name']
    overflow  = plotInfo.get('draw_overflow' , False)
    underflow = plotInfo.get('draw_underflow', False)

    if  (verbosity >= 1):
        print(Red("\n###################    DATA    ###################\n"))
    sample = samplesByRegion.data_obs
    hdata = None
    
    isFirst=1
    for fname in sample['files']:
        h = get_plots(inputdir, fname, [plot])[0]
        if(not h):
            continue

        if forcePositive and any(cr in inputdir.region for cr in ['CR2P2F','CR100','CR010','CR001']):
            h.Scale(-1)
        
        if not h:
            if  (verbosity >= 1):
                print('{} has no entries or is a zombie'.format(fname))
            continue
        
        if  (verbosity >= 1):
            print("{} in {} .......................... {}". format(fname, inputdir.region, h.Integral(0,-1)))
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
        c_err = ctypes.c_double(0)
        print("Total data in {0:s} region .......................... {1:.2f} +- {2:.2f}".format(inputdir.region, hdata.IntegralAndError(0,-1, c_err), c_err.value))
        print("_________________________")

    return hdata


def graph_and_ratio(histodata, hStackSum, xedges, bx_min, bx_max, unblind=True, remove_zeros=False):
    '''
    Returns two TGraphAsymmErrors: one the upper pad, and the ratio of histodata and the MC sum
    '''
    tgaData = ROOT.TGraphAsymmErrors()
    tgaData.SetName('ratio')
    if(not unblind):
        return ROOT.TGraphAsymmErrors(histodata), tgaData

    # Create new data and MC histograms with possibly the under-/overflow bins
    # This is the only way to include them in the TGraph resulting from Divide()
    tmpdata = ROOT.TH1F(histodata.GetName()+'_tmpdata', histodata.GetTitle(), len(xedges)-1, xedges)
    tmpMC   = ROOT.TH1F(histodata.GetName()+'_tmpMC'  , hStackSum.GetTitle(), len(xedges)-1, xedges)
    b_new = 0  # bin index in the new histograms
    for b_old in range(bx_min, bx_max+1):
        b_new += 1
        tmpdata.SetBinContent(b_new, histodata.GetBinContent(b_old))
        tmpdata.SetBinError  (b_new, histodata.GetBinError  (b_old))
        tmpMC  .SetBinContent(b_new, hStackSum.GetBinContent(b_old))
        tmpMC  .SetBinError  (b_new, hStackSum.GetBinError  (b_old))

    graphData = SetError(tmpdata, '', False)

    # Remove points where data = 0
    if(remove_zeros):
        remove_zeros_tg(graphData)

    tgaData.Divide(tmpdata, tmpMC, 'pois')

    del tmpdata, tmpMC

    return graphData, tgaData

def remove_zeros_tg(tg):
    to_remove = []
    for i in range(tg.GetN()):
        y = tg.GetPointY(i)
        if(y <= 1e-6):
            to_remove.append(i)
    for i in reversed(to_remove):
        # loop backwards to avoid changing the point indices while deleting
        tg.RemovePoint(i)


def integral_and_error(h, binx1=0, binx2=-1, option=""):
    '''
    Wrapper around TH1::IntegralAndError to avoid dealing with a raw double*
    '''
    err = ctypes.c_double(0.)
    integral = h.IntegralAndError(binx1, binx2, err, option)
    return integral, err.value


def cmsDiCanvas_fromTH1(name, h, r, y_scale=1, range_include_err=False, **kwargs):
    cmsargs = dict()
    x_min, x_max = getTAxisLimits(h.GetXaxis())
    cmsargs['x_min'] = kwargs.get('x_min', x_min)
    cmsargs['x_max'] = kwargs.get('x_max', x_max)

    cmsargs['y_min'] = kwargs.get('y_min', h.GetMinimum())
    cmsargs['y_max'] = kwargs.get('y_max', h.GetMaximum()*y_scale)

    if(not ('r_min' in kwargs and 'r_max' in kwargs)):
        for argname in ('y_scale', 'min_lo', 'max_lo', 'min_hi', 'max_hi'):
            # massage arg names for clamp_expnd_r()
            if(argname+'_r' in kwargs): kwargs[argname] = kwargs.pop(argname+'_r')
        if(r.GetN() > 0):
            r_min, r_max = get_range_tga(r, include_err=range_include_err)
        else:
            r_min, r_max = 0, 8
        r_min, r_max = clamp_expnd_r(r_min, r_max, **kwargs)
    r_min = kwargs.get('r_min', r_min)
    r_max = kwargs.get('r_max', r_max)

    cmsargs['nameXaxis'] = kwargs.get('nameXaxis', h.GetXaxis().GetTitle())
    cmsargs['nameYaxis'] = kwargs.get('nameYaxis', h.GetYaxis().GetTitle())
    cmsargs['nameRatio'] = kwargs.get('nameRatio', r.GetYaxis().GetTitle())
    if('iPos' in kwargs): cmsargs['iPos'] = kwargs['iPos']

    c = cmsstyle.cmsDiCanvas('canvas_%s' %(name), r_min=r_min, r_max=r_max, **cmsargs)

    return c


def getTAxisLimits(axis):
    return \
        axis.GetBinLowEdge(1), \
        axis.GetBinLowEdge(axis.GetNbins()+1)


def get_range_tga(g, include_err=False):
    '''
    Get the y range needed to draw a TGraphAsymmErrors
    '''
    name = g.GetName()

    if(include_err):
        np = g.GetN()
        y_max = max( p for p in (g.GetPointY(i)+g.GetErrorYhigh(i) for i in range(np)) )
        y_min = min( p for p in (g.GetPointY(i)-g.GetErrorYlow (i) for i in range(np)) )
    else:
        buf = array('d', g.GetY())
        y_max = max( buf )
        y_min = min( buf )
    logging.debug('%s range (raw): [%.3g, %.3g]', name, y_min, y_max)

    return y_min, y_max


def clamp_expnd_r(lo, hi, y_scale=0.1, min_lo=0., max_lo=0.9, min_hi=1.1, max_hi=100, name='[range]', **kwargs):
    '''
    Massage the range [lo, hi]: enlarge it by (hi-lo)*y_scale, 
    and clamp both (lo|hi) between [min_(lo|hi), max_(lo|hi)]
    '''

    def clamp(v, min_y, max_y):
        return min(max(v, min_y), max_y)

    delta = hi - lo
    lo = clamp(lo - y_scale * delta, min_lo, max_lo)
    hi = clamp(hi + y_scale * delta, min_hi, max_hi)
    logging.debug('%s range (fix): [%.3g, %.3g]', name, lo, hi)

    return lo, hi


def debugTGA(g):
    for i in range(g.GetN()):
        print('%5.3g [%.3g, %.3g]' %(g.GetPointX(i), g.GetErrorXlow(i), g.GetErrorXhigh(i)))
