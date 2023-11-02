#!/usr/bin/env python3

from __future__ import print_function
import sys, os
from argparse import ArgumentParser
import logging
from ctypes import c_double
import ROOT
from utils23 import lumi_dict
from plotUtils23 import TFileContext, InputDir

def fix_negative_bins(hnum, hden):
    for b in range(0, hden.GetNbinsX()+1):
        num = hnum.GetBinContent(b)
        den = hden.GetBinContent(b)
        if(num > den):
            hnum.SetBinContent(b, den)


def efficiency(tfile, label, var, outdir, canvas, c_debug):
    effNameFormat = "Eff_{label:s}_{ND:s}_{var:s}"  # to be formatted
    effNameFormat = effNameFormat.format(label=label, var=var, ND='{ND:s}')
    effTitle = effNameFormat.replace('_{ND:s}', '')
    numerator   = effNameFormat.format(ND="NUM")
    denominator = effNameFormat.format(ND="DEN")

    logging.info("### %s / %s ###", numerator, denominator)
    hnum = tfile.Get(numerator)
    hden = tfile.Get(denominator)
    if(not hnum or not hden):
        if(not hnum): logging.error('could not get "{}" from "{}"'.format(numerator  , tfile.GetName()) )
        if(not hden): logging.error('could not get "{}" from "{}"'.format(denominator, tfile.GetName()) )
        return

    err = c_double(0)
    int_num = hnum.IntegralAndError(0, -1, err)
    err_num = err.value
    int_den = hden.IntegralAndError(0, -1, err)
    err_den = err.value

    logging.debug("bins   :  num: %d - den: %d"    , hnum.GetNbinsX() , hden.GetNbinsX() )
    logging.debug("entries:  num: %.0f - den: %.0f", hnum.GetEntries(), hden.GetEntries())
    logging.debug("integr :  num: %.2f - den: %.2f", int_num          , int_den          )
    logging.info(" <eff>  :  %.3f += %.3f", int_num/int_den, ((err_num/int_den)**2 + (err_den*int_num/int_den**2)**2)**0.5)

    fix_negative_bins(hnum, hden)

    # c_debug.cd()
    # hden.SetTitle( effTitle )
    # hden.Draw("hist")
    # hnum.SetLineColor(ROOT.kRed)
    # hnum.Draw("same")
    # c_debug.SaveAs(outdir+'debug_'+hden.GetTitle()+'.png')

    canvas.cd()
    geff = ROOT.TGraphAsymmErrors(hnum, hden, "")
    geff.SetTitle( effTitle )
    geff.GetYaxis().SetTitle('Efficiency')
    geff.SetLineColor(ROOT.kBlack)
    geff.SetMarkerStyle(ROOT.kFullTriangleUp)
    geff.SetMarkerColor(geff.GetLineColor())
    geff.Draw("APLE")
    geff.GetXaxis().SetTitle(hden.GetXaxis().GetTitle())
    canvas.SaveAs(outdir+geff.GetTitle()+'.png')


# def resolution(tfile, label, var, c, c_debug):
#     resNameFormat = "Res_%s_%s".format(label, var)
def resolution2D(tfile, label, var, canvas):
    resNameFormat = "Res_{label}_{var}".format(label=label, var=var)
    print("###", resNameFormat, "###")
    hRes = tfile.Get(resNameFormat)
    if(not hRes):
        print('ERROR: could not get "{}" from {}'.format(resNameFormat, tfile.GetName()) )
        return

    assert hRes.Class().InheritsFrom("TH2")
    canvas.cd()
    profile = hRes.ProfileX()
    profile.SetTitle(resNameFormat)
    profile.Draw()
    canvas.SaveAs(outdir+profile.GetTitle()+'.png')


def main():
    parser = ArgumentParser()
    parser.add_argument('-i', '--inputdir', default='results', help='Top directory containing input (default: %(default)s)')
    parser.add_argument('-y', '--year'    , default='2018', choices=lumi_dict.keys(), help='Default: %(default)s')
    parser.add_argument('-r', '--region'  , default='SR4P', help='Default: %(default)s')
    parser.add_argument('-A', '--analyzer', default='VVGammaAnalyzer', help='Default: %(default)s')
    parser.add_argument('-s', '--sample'  , default='ZZGTo4LG', help='Default: %(default)s')
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')

    args = parser.parse_args()
    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    ROOT.gROOT.SetBatch(True)

    inputDir = InputDir(args.inputdir, year=args.year, region=args.region, analyzer=args.analyzer)
    filename = os.path.join(inputDir.path(), args.sample+'.root')

    # output
    outdir = 'Plot/Efficiencies/'
    os.makedirs(outdir, exist_ok=True)

    # (numerator, denominator)
    labels = [
        # "AK4_quarks",
        # "AK4_genAK4",
        # "genAK4_quarks",
        "AK4_genJets",
        "goodPhotons_gen",
        "kinPhotons_gen",
        "goodPhotons_genKin",
        "goodPhotons_genDRl",
        "goodPhotons_genPrompt",
        "goodPhotons_genKinPrompt",
        "goodPhotons_genKinDRl"
    ]
    effVars = {'pt', 'E', 'eta', 'N'}
    # resVars = {'dR', 'E', 'pt'}
    # resVars2D = {'EvsE'}

    c = ROOT.TCanvas("c", "The canvas", 0, 0, 1000, 1000)
    c_debug = ROOT.TCanvas("c_debug", "Num and den", 0, 0, 1000, 1000)

    with TFileContext(filename) as tfile:
        for label in labels:
            for var in effVars:
               efficiency(tfile, label, var, outdir, c, c_debug)

            # for var in resVars:
            #     resolution(tfile, label, var, c)

            # for var in resVars2D:
            #     resolution2D(tfile, label, var, c)

    return 0

if __name__ == '__main__':
    exit(main())
