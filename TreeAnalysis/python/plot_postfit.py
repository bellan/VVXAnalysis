#!/usr/bin/env python
from argparse import ArgumentParser
import logging
import os
from math import sqrt
from copy import deepcopy
import ROOT
import cmsstyle
from array import array
from subprocess import run
from ctypes import c_double

from utils23 import config_logging, lumi_dict
from plotUtils23 import TFileContext, addIfExisting, cmsDiCanvas_fromTH1, getTAxisLimits
from PersonalInfo import personalFolder
import samplesByRegion


_varinfo = {
    'mZZG': {'bins': array('d', range(0, 1100, 100)), 'xtitle': 'm_{4l#gamma} [GeV]'},
    'pt'  : {'bins': array('d', [20., 25., 35., 50., 80., 120.]), 'xtitle': 'p_{T}^{#gamma} [GeV]'},
}

_samplesinfo = {
    'signal'      :{'color': samplesByRegion.ZZG[0]['color']},
    'qqZZ'        :{'color': samplesByRegion.qqZZ_pow[0]['color']},
    'ggZZ'        :{'color': samplesByRegion.ggZZ[0]['color']},
    'fake_photons':{'color': samplesByRegion.fake_photons['color']},
    'fake_leptons':{'color': samplesByRegion.fake_leptons['color']},
    'rare_bkg'    :{'color': samplesByRegion.rare_4l[0]['color']}
}

_SHAPE_LABELS = {'fit_s': 'Post-fit', 'fit_b': 'Bkg. only fit', 'prefit': 'Pre-fit'}


def main(args):
    logging.debug('args = %s', args)
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    with TFileContext(args.workspace) as tf:
        h_years_map = get_hists(tf, shapes=args.shapes)
    logging.debug('retrieved keys = %s', h_years_map.keys())

    # Deduce if it is the result of a triboson card by the folder name
    if(args.isTriboson is None):
        args.isTriboson = os.path.split(args.workspace)[0].endswith('_triboson')
        logging.debug('Setting triboson label: %s', args.isTriboson)

    outname = os.path.join(args.out,
                           os.path.splitext(os.path.basename(args.workspace))[0]\
                           .replace('fitDiagnostics_','')\
                           .replace('Run2_', args.year+'_')
                           +'_'+args.shapes
                           )

    if(args.year == 'Run2'):
        if(args.cut_n_count):
            h_map = join_hists_year(h_years_map)
        else:
            h_map = sum_hists(h_years_map)
    else:
        h_map = h_years_map[args.year]
        h_map['data'] = tga2hist(h_map['data'])

    # Set x bin edges, since Combine discards this information when creating the workspace
    varinfo = _varinfo['mZZG']
    if('_ptloose' in args.workspace or '_ptwp' in args.workspace):
        varinfo = _varinfo['pt']

    if(not args.cut_n_count):
        h_map = fix_binning(h_map, varinfo['bins'])

    hdata = h_map.pop('data')
    hdata.GetXaxis().SetTitle(varinfo['xtitle'])
    h_map_grouped = group_hists(h_map, isTriboson=args.isTriboson)

    info_list = sort_h_map(h_map_grouped)
    logging.debug('Plotting these MCs = %s', ['%s ("%s")' %(i['name'], i['title']) for i in info_list])

    # Customize style
    cmsstyle.SetExtraText("")
    cmsstyle.setCMSStyle()
    cmsstyle.SetLumi(138)
    ROOT.gStyle.SetLabelSize(0.045, "X")

    err = plot(hdata, info_list, outname=outname, **vars(args))

    return err


def parse_args():
    parser = ArgumentParser(epilog='For the PAS, --yscale was: 1.8 (inclusive), 2.3 (triboson)')
    parser.add_argument('workspace', metavar='rootfile', help='FitDiagnostics output (ROOT file)')
    parser.add_argument('-y', '--year', dest='year',
                        default='Run2',
                        help= 'valid inputs are 2016preVFP, 2016postVFP, 2017, 2018, Run2')
    parser.add_argument(      '--triboson', action='store_true', dest='isTriboson', default=None,
                              help='Set the legend entry for ZZG (default:%(default)s)')
    parser.add_argument(      '--no-triboson', action='store_false', dest='isTriboson')
    parser.add_argument(      '--region-label', action='store_true', default=False, help='Draw the region label (default: %(default)s)')
    parser.add_argument(      '--no-region-label', action='store_false', dest='region_label')
    parser.add_argument(      '--postfit-label', action='store_true', help='Add '+'/'.join([v for k,v in _SHAPE_LABELS.items()])+' to the region label')
    parser.add_argument(      '--yscale', type=float, default=1.8,
                              help='Factor that scales y_max in the upper plot (default: %(default)s)')
    parser.add_argument(      '--y_max', type=float, default=None, help='Set y_max in the upper plot (default: %(default)s)')
    parser.add_argument(      '--r_max', type=float, default=None, help='Set r_max in the lower plot (default: %(default)s)')
    parser.add_argument(      '--shapes', choices=_SHAPE_LABELS.keys(), default='fit_s',
                              help='Name of the folder in the FitDiagnostics file that contains the histograms (default: %(default)s)')
    parser.add_argument(      '--cut-n-count', action='store_true',
                              help='In case there is only one bin per year')
    parser.add_argument('-o', '--out', default=personalFolder, help='Output directory for plots (default:%(default)s)')
    parser.add_argument(      '--ext', default=['png'], nargs='+', help='Format(s) for the images produced (default: %(default)s)')
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')

    return parser.parse_args()


def plot(hdata, info_list, isTriboson=False, outname='postfit', ext=['png'], yscale=1.8, **kwargs):
    for b in range(0, hdata.GetNbinsX()+2):
        logging.debug('>>> %d: %f', b, hdata.GetBinContent(b))

    stack = mk_stack(info_list)

    ratio = ROOT.TGraphAsymmErrors()
    ratio.SetName('ratio')
    logging.debug('data: %s', hdata)
    logging.debug('MC  : %s', stack.GetStack().Last())
    ratio.Divide(hdata, stack.GetStack().Last(), 'pois')

    # Create the canvas
    dicanvas_kwargs = dict(y_min=0, y_scale=yscale, min_hi_r=2., max_lo_r=0., range_include_err=True,
                           nameYaxis='Events', nameRatio='Data/Pred.', iPos=0)
    if(args.y_max is not None): dicanvas_kwargs['y_max'] = args.y_max
    if(args.r_max is not None): dicanvas_kwargs['r_max'] = args.r_max
    canvas = cmsDiCanvas_fromTH1(args.shapes, hdata, ratio,
                                 **dicanvas_kwargs)
    if(hdata.GetXaxis().IsAlphanumeric()):
        logging.info('alphanumeric axis')
        canv_hist = cmsstyle.GetcmsCanvasHist(canvas.cd(2))
        xaxis_canv = canv_hist.GetXaxis()
        xaxis_canv.SetAlphanumeric()
        xaxis_d = hdata.GetXaxis()
        xaxis_canv.Set(xaxis_d.GetNbins(), xaxis_d.GetXmin(), xaxis_d.GetXmax())
        for b in range(1, xaxis_d.GetNbins()+1):
            xaxis_canv.SetBinLabel(b, xaxis_d.GetBinLabel(b))
        canvas.RedrawAxis()
        cmsstyle.UpdatePad(canvas)
    canvas.cd()

    # The legend needs to be created after the canvas, otherwise it won't be drawn
    legend = mk_legend(info_list)

    ### Upper pad ###
    canvas.cd(1)

    # Region and postfit label
    if(args.region_label or args.postfit_label):
        pad = ROOT.gPad
        y_top = 1 - pad.GetTopMargin() - 0.025
        y_bot = y_top - 0.05 - 0.06*(int(args.region_label) + int(args.postfit_label))

        region_text = ROOT.TPaveText(
            pad.GetLeftMargin()+0.05, y_top,
            0.5, y_bot, "NB NDC")
        region_text.SetTextAlign(ROOT.ETextAlign.kHAlignLeft + ROOT.ETextAlign.kVAlignTop)
        region_text.SetTextSize(.05)
        region_text.SetFillColor(ROOT.kWhite)

        if(args.region_label):
            region_text.AddText( 'Triboson ZZ#gamma' if isTriboson else 'Inclusive pp #rightarrow 4l#gamma' )

        if(args.postfit_label):
            region_text.AddText( _SHAPE_LABELS[kwargs['shapes']])

        region_text.Draw('same')

    # Error band in the upper canvas
    hMCErr = deepcopy(stack.GetStack().Last())

    hMCErr.SetFillStyle(3005)
    hMCErr.SetMarkerStyle(1)
    hMCErr.SetFillColor(ROOT.kBlack)
    legend.AddEntry(hMCErr, "Stat. only", "f")

    # Style data
    hdata.SetLineColor(ROOT.kBlack)
    hdata.SetMarkerStyle(20)
    hdata.SetMarkerSize(.8)
    hdata.SetBinErrorOption(ROOT.TH1.kPoisson)
    legend.AddEntry(hdata, 'data', 'lpe')

    # Draw
    stack.Draw('SAMEHIST')
    hMCErr.Draw("SAMEE2")
    hdata.Draw('SAMEPE0X0')

    ### Lower pad ###
    canvas.cd(2)

    # Line y=1 in the ratio plot
    x_min, x_max = getTAxisLimits(hdata.GetXaxis())
    logging.debug('x_min=%.3g, x_max=%.3g', x_min, x_max)
    ref_line = ROOT.TLine(x_min, 1, x_max, 1)
    cmsstyle.cmsDrawLine(ref_line, lcolor=ROOT.kBlack, lstyle=ROOT.kDotted)

    # Ratio
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(.8)

    # Draw
    ratio.Draw('PE')

    for e in ext:
        if e == 'root': continue
        outfname = '.'.join([outname, e])
        canvas.SaveAs(outfname)
        if('workspace' in kwargs):
            campaign = os.path.split( os.path.dirname(kwargs['workspace']) )[-1]
            cmd = ['exiftool', '-overwrite_original', '-Keywords=%s'%(campaign), outfname]
            logging.info('running: %s', ' '.join(cmd))
            run(cmd) # subprocess; willingly ignore errors

    return 0


def get_hists(tf, shapes='fit_s'):
    '''
    Retrieve histograms (and a TGraphAsymmErrors for data) from the output of
    Combine's FitDiagnostics

    Return schema: {year: {process: <TH1F>, "data": TGraphAsymmErrors}}
    '''
    # In each year there are:
    # - "data": <TGraphAsymmErrors>
    # - several MC: <TH1F>
    # - total(_signal|_background): <TH1F>
    # - total_covar: <TH2F>

    h_map = dict()

    shapes_dir = tf.Get('shapes_'+shapes)
    for k_year in shapes_dir.GetListOfKeys():
        # Each bin corresponds to a year
        year = k_year.GetName().lstrip('y')
        logging.debug('year = %s', year)
        assert k_year.IsFolder(), 'k_year is of type %s' %(k_year.GetClassName())

        for k in k_year.ReadObj().GetListOfKeys():
            # logging.debug('    k = %s (%s)', k.GetName(), k.GetClassName())
            name = k.GetName()
            if(name.startswith('total')):
               # logging.debug('        skipped')
               continue

            obj = k.ReadObj()
            if(isinstance(obj, ROOT.TH1)): obj.SetDirectory(0) # disable ROOT's broken garbage collector

            # Schema: sample_name -> year -> hist
            h_map.setdefault(year, {})[name] = obj

    return h_map


def fix_binning(h_map_in, bin_edges):
    '''
    Create a new dictionary of TH1F, with the x axis set as per `bin_edges`,
    for each histogram that is in the first argument.
    '''
    nb = len(bin_edges) - 1
    buf = array('d', bin_edges)
    h_map_out = dict()

    for name, h_old in h_map_in.items():
        # Check that the supplied bin edges are ok
        assert h_old.GetNbinsX() == nb, 'Wrong number of bins: %d (expected %d)' %(h_old.GetNbinsX(), nb)

        h_new = ROOT.TH1F(h_old.GetName(), h_old.GetTitle(), len(buf)-1, buf)
        for bx in range(nb+1):
            h_new.SetBinContent(bx, h_old.GetBinContent(bx))
            h_new.SetBinError  (bx, h_old.GetBinError  (bx))

        h_map_out[name] = h_new

    return h_map_out


def sum_hists(in_map):
    '''
    Get the total yield for each group of processes for the whole Run2.
    The sum of the TGraphAsymmErrors for "data" is summed in a TH1F as well.

    Return schema: {sample_group: <TH1F>}
    '''
    out_map = dict()
    data_x = None
    data_y = None # manual sum of TGraphs

    for _, processes in in_map.items():
        for proc, hist in processes.items():
            if(proc == 'data'):
                buf = array('d', hist.GetY())
                if(data_y is None):
                    data_y = buf
                else:
                    for i in range(len(data_y)):
                        data_y[i] += buf[i]
                continue

            if(proc in out_map):
                out_map[proc].Add(hist)
            else:
                out_map[proc] = hist
                if(data_x is None):
                    nb = hist.GetNbinsX()
                    data_x = array('d', [0.]*nb)
                    hist.GetXaxis().GetLowEdge(data_x)
                    data_x.append(hist.GetBinLowEdge(nb+1))

    # Handle data specially
    logging.debug('data_x = %s', data_x)
    logging.debug('data_y = %s', data_y)

    data = ROOT.TH1F('data', '', len(data_x)-1, data_x)
    for b in range(len(data_x)-1):
        data.SetBinContent(b+1, data_y[b])
        # data.SetBinError  (b+1, sqrt(data_y[b]))

    # data = ROOT.TGraphAsymmErrors(len(data_y), data_x, data_y)
    # data.SetName('data')
    data.SetBinErrorOption(ROOT.TH1.kPoisson)
    out_map['data'] = data

    return out_map


def join_hists_year(in_map):
    '''
    Get the yield for each group of processes for each year and put them in
    a single histogram. The TGraphAsymmErrors for "data" is converted to a TH1F as well.

    Return schema: {sample: <TH1F>}
    '''
    out_map = dict()
    years_sorted = sorted(in_map.keys(),
                          key=lambda y: (
                              int(y[:4]), # if only Python's atoi() behaved like C
                              1 if 'post' in y else -1 if 'pre' in y else 0
                          ))

    for year, processes in in_map.items():
        for proc, hist in processes.items():
            if(not proc in out_map):
                hnew = ROOT.TH1F(hist.GetName(), hist.GetTitle(), len(years_sorted),0,len(years_sorted))
                for b,y in enumerate(years_sorted):
                    hnew.GetXaxis().SetBinLabel(b+1, y)
                out_map[proc] = hnew

            if(proc == 'data'):
                tot = sum( array('d', hist.GetY()) )
                b = hnew.GetXaxis().FindFixBin(year)
                out_map[proc].SetBinContent(b, tot)
            else:
                err = c_double()
                tot = hist.IntegralAndError(0, -1, err)
                hnew = out_map[proc]
                b = hnew.GetXaxis().FindFixBin(year)
                hnew.SetBinContent(b, tot)
                hnew.SetBinError  (b, err)

    out_map['data'].SetBinErrorOption(ROOT.TH1.kPoisson)

    return out_map


def mk_stack(info_list):
    stack = ROOT.THStack("stack", "stack")

    for info in info_list:
        logging.debug('info: %s', info)
        title= info['title']
        hist = info['h']
        hist.SetFillColor(info['color'])
        hist.SetLineColor(ROOT.kBlack)
        logging.debug('%s -> %s', hist, title)
        stack.Add(hist)

    return stack


def mk_legend(info_list):
    ymax = .92
    ymin = ymax - 0.05*(len(info_list)+2)  # +2: MC stat, data
    # logging.debug('nhist = %d+2 - y = [%.2f, %.2f]', len(info_list), ymin, ymax)
    legend = cmsstyle.cmsLeg(.55, ymin, .90, ymax, textSize=.03)
    for info in info_list:
        title= info['title']
        hist = info['h']
        logging.debug('%s -> %s', hist, title)
        legend.AddEntry(hist, title, 'f')

    return legend


def sort_h_map(h_map):
    info_list = []
    for proc, data in h_map.items():
        data['name'] = proc
        info_list.append(data)
    info_list.sort(key=lambda x: x.get('key', 99), reverse=True)

    return info_list


def group_hists(h_map_ungrouped, isTriboson=False):
    '''
    Sums hist of the same group and assigns titles for the legend
    and assign them a color
    '''
    h_map = {}
    for sample, hist in h_map_ungrouped.items():
        if('-' in sample):
            base, extra = sample.split('-')
            nonpro = (extra == 'nonpro')
            extra_t = ' non-prompt' if nonpro else ''
            extra_k = 1 if nonpro else 0
        else:
            base = sample
            extra = ''
            extra_t = ''
            extra_k = 0

        if  (base == 'ZZGTo4LG' or base == 'signal'):
            title = 'ZZ#gamma' if isTriboson else '4l #gamma'
            # if(extra == 'nonpro'): title += ' OSD'
            h_map.setdefault(base, dict(
                title=title, color=_samplesinfo['signal']['color'], key=extra_k, hlist=[]
            ))['hlist'].append(hist)
        elif(base == 'ZZTo4l'):
            h_map[sample] = dict(h=hist, title='qq #rightarrow ZZ'+extra_t, color=_samplesinfo['qqZZ']['color'], key=4+extra_k)
        elif(base.startswith('ggTo')):
            h_map.setdefault('ggTo4l'+extra, dict(
                             title='gg #rightarrow ZZ'+extra_t, color=_samplesinfo['ggZZ']['color'], key=6+extra_k, hlist=[]
                             ))['hlist'].append(hist)
        elif(base == 'fake_photons'):
            h_map[sample] = dict(h=hist, title='Non-prompt #gamma', color=_samplesinfo['fake_photons']['color'], key=8)
        elif(base == 'fake_leptons'):
            h_map[sample] = dict(h=hist, title='Non-prompt l', color=_samplesinfo['fake_leptons']['color'], key=9)
        else:
            h_map.setdefault('rare_bkg', dict(
                             title='Rare backgrounds', color=_samplesinfo['rare_bkg']['color'], key=2, hlist=[]
                             ))['hlist'].append(hist)

    for _, data in h_map.items():
        if('hlist' in data):
            logging.debug('grouping "%s"', _)
            data['h'] = addIfExisting(*data.pop('hlist'))
    return h_map


def tga2hist(tga):
    data_y = array('d', tga.GetY())
    # Combine discards x bin information, xaxis goes from 0 to n in n steps.
    # This is fixed in fix_binning()

    n = len(data_y)
    h = ROOT.TH1F(tga.GetName(), '', n, 0, n)
    for b in range(n):
        h.SetBinContent(b+1, data_y[b])

    return h


if __name__ == '__main__':
    args = parse_args()
    config_logging(args.loglevel)

    exit(main(args))
