from array import array
from math import sqrt
import os
import ROOT

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
