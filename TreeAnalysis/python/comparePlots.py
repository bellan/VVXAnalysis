#!/usr/bin/env python2

import os
from argparse import ArgumentParser
import ctypes
import logging
import re
import ROOT

from plotUtils23 import TFileContext, InputDir, InputFile, get_plots
# from plotUtils import GetPredictionsPlot
from variablesInfo import getVariablesInfo
from utils23 import lumi_dict
import CMS_lumi, tdrstyle

def parse_args():
    regions = ['SR4P', 'CR3P1F' , 'CR2P2F' , 'SR4P_1L', 'SR4P_1P', 'CR4P_1F', 'CR4L',    
               'SR3P', 'CR110'  , 'CR101'  , 'CR011'  , 'CR100'  , 'CR001'  , 'CR010', 'CR000', 'SR3P_1L', 'SR3P_1P', 'CR3P_1F', 'CRLFR', 'CR3L',
               'SR2P', 'SR2P_1L', 'SR2P_1P', 'CR2P_1F', 
               'SR_HZZ', 'CR2P2F_HZZ', 'CR3P1F_HZZ', 'CR_HZZ']
    
    parser = ArgumentParser()
    
    parser.add_argument('-i', '--inputdir'      , default='results', help='Base input directory, where the results of the analyzers are')
    parser.add_argument('-r', '--region'        , default='SR4P', help='Default: %(default)s')
    parser.add_argument('-y', '--year'          , default='2018', choices=lumi_dict.keys(), help='Default: %(default)s')
    parser.add_argument('-A', '--analysis'      , default='VVGammaAnalyzer', dest='analysis', help='Default: %(default)s')
    parser.add_argument('-t', '--variables'     , default=None, dest='var_include', type=re.compile, help='Only plot names that match this regex will be used. Defaults to everything')
    parser.add_argument('-s', '--skip'          , default=None, dest='var_skip'   , type=re.compile, help='Plots names that match this regex will be skipped')
    parser.add_argument(      '--rebin'         , default=1, type=int)
    parser.add_argument(      '--no-split'      , dest='split_prompt' , action='store_false', help='Never split prompt/nonprompt')
    parser.add_argument(      '--no-title'      , dest='do_title'     , action='store_false', help='Do not paint the title on the canvas (default: False)')
    parser.add_argument(      '--no-norm'       , dest='do_norm'      , action='store_false', help='Do not normalize shapes')
    parser.add_argument(      '--ext'           , default=['png'], dest='extensions', nargs='+', help='Format(s) for output images (default: %(default)s)')
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')

    return parser.parse_args()


def cutStudy(var, plotInfo, inputdir, do_title=True, do_norm=True, extensions=['png'], rebin=1, **kwargs):
    canvas = ROOT.TCanvas(var+'_canvas', var, 1600, 1200)

    canvas.cd()
    pad1 = ROOT.TPad('hist', '', 0., 0.28, 1.0, 1.0)
    pad1.SetTopMargin   (0.10)
    pad1.SetRightMargin (0.04)
    pad1.SetLeftMargin  (0.13)
    pad1.SetBottomMargin(0.02)
    pad1.Draw()
    
    canvas.cd()
    pad2 = ROOT.TPad('ratio', 'S/B ratio', 0., 0.0,  1., 0.28)
    pad2.SetTopMargin   (0.02)
    pad2.SetRightMargin (0.04)
    pad2.SetLeftMargin  (0.13)
    pad2.SetBottomMargin(0.3);
    pad2.Draw()

    pad1.cd()

    splitPromptPh = plotInfo.get('split_prompt_ph')
    if(splitPromptPh):
        split_pattern = plotInfo.get('split_prompt_ph_pattern', var+'_%s')
        logging.info('Splitting prompt/nonprompt photons')
        to_get = [split_pattern %('prompt'), split_pattern %('nonpro')]
    else:
        to_get = [var]
    logging.info('to_get: %s', to_get)

    error = ctypes.c_double(0.)
    sigstack  = None
    bkgstacks = []
    legend_x0 = 0.95
    legend_y0 = 0.89
    legend = ROOT.TLegend(legend_x0 - 0.20,
                          legend_y0 - 0.06*len(hardcoded_config['groups'])*(2 if splitPromptPh else 1)
                          , legend_x0, legend_y0)

    nbinsx = None
    for group, gdata in hardcoded_config['groups'].items():
        logging.debug('Starting group: %s', group)
        groupstack = ROOT.THStack(group+'_stack', var if do_title else '')
        # legend.AddEntry(groupstack, group, "LPF")

        hist_counter = 0
        for fname in gdata['files']:
            hists = get_plots(inputdir, fname, to_get)
            for hname, hist in zip(to_get, hists):
                if(not hist):
                    continue
                if(hist.GetDimension() > 1):
                    logging.error('histograms with dimension %d > 1 are not supported', hist.GetDimension())
                    return 1

                hist_counter += 1
                hist.Rebin(rebin)
                if  (nbinsx is None):
                    nbinsx = hist.GetNbinsX()
                elif(nbinsx != hist.GetNbinsX()):
                    raise RuntimeError('Histograms have number of bins (previous: %d, now: %d)' %(nbinsx, hist.GetNbinsX()))

                hist.SetLineColor(gdata['color'])
                hist.SetLineStyle(hist_counter)
                hist.SetLineWidth(2)
                if(splitPromptPh and 'fill_style_prompt' in gdata and 'fill_style_nonpro' in gdata):
                    hist.SetFillColor(gdata['color'])
                    hist.SetFillStyle(gdata['fill_style_prompt'] if 'prompt' in hist.GetName() else gdata['fill_style_nonpro'])
                elif('fill_style' in gdata):
                    hist.SetFillColor(gdata['color'])
                    hist.SetFillStyle(gdata['fill_style'])

                groupstack.Add(hist)
                legend.AddEntry(hist, fname if not splitPromptPh else fname+' '+('nonpro' if 'nonpro' in hname else 'prompt'), "L")
                del hist

        if(sigstack is None): sigstack = groupstack
        else: bkgstacks.append(groupstack)
        logging.debug('Finished group "%s"', group)

    logging.debug('Finished retrieving plots')

    ymax = 0
    for stack in [sigstack]+bkgstacks:
        if(do_norm):
            logging.debug('Normalizing stack %s', stack)
            last = stack.GetStack().Last()
            integral = last.IntegralAndError(1, nbinsx, error)
            for h in stack.GetHists():
                h.Scale(1./integral)
            stack.Modified()

        ymax = max(ymax, 1.1*stack.GetStack().Last().GetBinContent(stack.GetStack().Last().GetMaximumBin()))

    sigstack.SetMaximum(ymax)
    sigstack.Draw('hist')
    sigstack.GetHistogram().GetXaxis().SetLabelSize(0)
    if(do_norm):
        sigstack.GetHistogram().GetYaxis().SetTitle('Normalized shapes')
    sigstack.GetHistogram().GetYaxis().SetTitleOffset(1.)
    logging.debug('Upper title offset %f', sigstack.GetHistogram().GetYaxis().GetTitleOffset())
    for stack in bkgstacks:
        stack.Draw('hist same')
    legend.Draw('same')

    hsig = sigstack.GetStack().Last()
    hbkg = bkgstacks[0].GetStack().Last()
    for stack in bkgstacks[1:]:
        hbkg.Add(stack.GetStack().Last())

    CMS_lumi.CMS_lumi(pad1, iPeriod=0, iPosX=0)

    pad2.cd()
    ratio = ROOT.TGraphAsymmErrors()
    ratio.Divide(hsig, hbkg, 'pois')
    ymax_ratio = max(y+ey for y,ey in zip(ratio.GetY(), ratio.GetEYhigh()))
    ymax_ratio = min(ymax_ratio, 4)
    # ymin_ratio = min(y-ey for y,ey in zip(ratio.GetY(), ratio.GetEYlow()))
    
    frame = hsig.Clone('frame_ratio')
    frame.GetYaxis().SetRangeUser(0, ymax_ratio*1.1)
    frame.GetXaxis().SetLabelSize(0.12)
    frame.GetXaxis().SetTitleSize(0.12)
    frame.GetYaxis().SetLabelSize(0.08)
    frame.GetYaxis().SetTitleSize(0.12)
    frame.GetYaxis().SetTitleOffset(0.2)
    frame.GetYaxis().SetTitle('S/B')
    frame.Draw('AXIS')
    
    ROOT.gStyle.SetOptStat(0)
    
    line_1 = ROOT.TLine(hsig.GetXaxis().GetBinLowEdge(hsig.GetXaxis().GetFirst()), 1, hsig.GetXaxis().GetXmax(), 1)
    line_1.SetLineWidth(2)
    line_1.SetLineStyle(7)
    line_1.Draw('SAME')

    ratio.Draw('PE')

    for ext in extensions:
        canvas.SaveAs('%s.%s' %(var, ext))

    return 0


hardcoded_config = {
    'groups': {
        "ZZG" : {'files':['ZZGTo4LG'      ], 'title':'ZZ#gamma' , 'color':ROOT.kRed    , 'fill_style':3345, 'fill_style_prompt': 3315, 'fill_style_nonpro':3351},
        "qqZZ": {'files':['ZZTo4l'        ], 'title':'qqZZ'     , 'color':ROOT.kBlue-4 , 'fill_style':3354, 'fill_style_prompt': 3365, 'fill_style_nonpro':3356}# , 'kfactor': 1.325/1.256}
    }
}


def search_variable(var_list, re_include, re_skip, analysis='unknown', region='unknown'):
    if re_include is None:
        variables = var_list
    else:
        variables = [ var for var in var_list if re_include.search(var) ]  # Allow for regexp to be specified from command line
        if len(variables) == 0:
            logging.warning('no variables matching regex "{}" for {} in {}'.format(re_include.pattern, analysis, region))

    if len(variables) > 0 and re_skip is not None:
        variables = [ var for var in variables if not re_skip.search(var) ]
        if len(variables) == 0:
            logging.warning('using regex {} all variables are skipped for {} in {}'.format(re_skip.pattern, analysis, region))

    return variables


def main():
    args = parse_args()

    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    ROOT.gROOT.SetBatch(True)
    tdrstyle.setTDRStyle()
    CMS_lumi.writeExtraText = True
    CMS_lumi.extraText = " Private work"
    CMS_lumi.lumi_sqrtS = "{0:.3g} fb^{{-1}} (13 TeV)\n".format( lumi_dict[args.year]['value']/1000. )
    CMS_lumi.relPosX = 0.12

    inputdir = InputDir(basedir=args.inputdir, year=args.year, region=args.region, analyzer=args.analysis)

    varInfo = getVariablesInfo(args.analysis, args.region)
    variables = search_variable(varInfo.keys(), args.var_include, args.var_skip, analysis=args.analysis, region=args.region)

    if(len(variables) == 0):
        # Try reading the TKeys from the signal file
        fname = hardcoded_config['groups']['ZZG']['files'][0]
        rootfile = InputFile(inputdir, fname+'.root')
        if(rootfile.inputDir.year == 'Run2'): rootfile.inputDir.year = '2018'
        logging.info('Search in getVariablesInfo unsuccesful. Trying with the TKeys in the file "%s"', rootfile)

        with TFileContext(rootfile.path()) as tf:
            keys = [k.GetName() for k in tf.GetListOfKeys()]
        variables = search_variable(keys, args.var_include, args.var_skip, analysis=args.analysis, region=args.region)

        for variable in variables:
            varInfo[variable] = {}  # Assign default config

    if(args.split_prompt == False):
        for _, info in varInfo.items(): info['split_prompt_ph'] = False

    logging.info('variables: %s', variables)

    argsdict = {k:v for k,v in vars(args).items() if not k in ('inputdir', 'year', 'region', 'analysis', 'loglevel', 'var_include', 'var_skip')}

    status = 0
    for variable in variables:
        status += cutStudy(variable, plotInfo=varInfo[variable], inputdir=inputdir, **argsdict)
    return status


if __name__ == '__main__':
    exit(main())
