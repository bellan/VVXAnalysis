#!/usr/bin/env python

################################################################
# Small helper script to produce images with a nice formatting #
# from a TH2 stored in a rootfile (e.g. for scale factors)     #
#                                                              #
# Authtor: A. Mecca (alberto.mecca@cern.ch)                    #
################################################################

import os
import sys
from argparse import ArgumentParser
import logging
import ROOT

# Find TreeAnalysis/python
def get_VVXdir():
    if('CMSSW_BASE' in os.environ):
        return os.path.join(environ['CMSSW_BASE'], 'src', 'VVXAnalysis')
    else:
        wd = os.getcwd()
        if('VVXAnalysis' in wd):
            return os.path.join(wd.split('VVXAnalysis')[0], 'VVXAnalysis')
    # Could not find VVXAnalysis
    return None


def print_SF_hist(tf, histname, logx=False, logy=False, style='', do_title=True, extensions=(), text_size=1., **kwargs):
    canvas = ROOT.TCanvas('canvas', 'canvas', 1600, 900)
    canvas.SetLogx(logx)
    canvas.SetLogy(logy)
    canvas.SetRightMargin(0.11)

    h = tf.Get(histname)
    if(not h):
        raise PlotNotFoundError(histname)
    else:
        logging.info('retrieved %s', histname)

    h.GetXaxis().SetNoExponent(logx)
    h.GetYaxis().SetNoExponent(logy)
    h.GetXaxis().SetMoreLogLabels(logx)
    h.GetYaxis().SetMoreLogLabels(logy)

    h.SetMarkerSize(text_size)  # For some reason ROOT developers decided that painted text's size should have the same multiplier as the marker size

    if(not do_title):
        h.SetTitle('')
    h.Draw(style)

    # Extra stuff
    h.GetXaxis().SetTitleSize(.05)
    h.GetYaxis().SetTitleOffset(0.8)

    for ext in extensions:
        canvas.SaveAs('{:s}_{:s}.{:s}'.format(os.path.splitext(tf.GetName())[0], histname, ext))
    del canvas


def main():
    parser = ArgumentParser()
    parser.add_argument('fnames', nargs='+', metavar='FILE')
    parser.add_argument('-n', '--histname', default='EGamma_SF2D', help='Name of the TH2 inside the file (default: %(default)s)')
    parser.add_argument(      '--logx'    , dest='logx', action="store_true")
    parser.add_argument(      '--logy'    , dest='logy', action="store_true")
    parser.add_argument(      '--format'  , default='.2f', help='printf-style format string to be passed to SetPaintTextFormat (default: %(default)s)')
    parser.add_argument(      '--text-size', type=float, default=1., metavar='SIZE', help='Factor that scales the text printed on each bin (default: %(default)f)')
    parser.add_argument(      '--style'   , default='colz texte', help='To be passed to TH1::Draw() (default: %(default)s)')
    parser.add_argument(      '--no-title', dest='do_title',action='store_false', help='Do not paint the title on the canvas (default: False)')
    parser.add_argument(      '--extensions', nargs='+', choices=('pdf', 'ps', 'eps', 'svg', 'tex', 'png', 'jpg', 'tiff'), default=('png', 'pdf',))
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')

    args = parser.parse_args()
    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    # Add TreeAnalysis/python to the path to import some utilities
    VVXdir = get_VVXdir()
    if(VVXdir is None):
        logging.fatal('Could not find VVXAnalysis')
    else:
        sys.path.append(os.path.join(VVXdir, 'TreeAnalysis', 'python'))

    from plotUtils23 import TFileContext

    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetPaintTextFormat(args.format.replace('%', ''))

    for fname in args.fnames:
        with TFileContext(fname) as tf:
            logging.info('opened %s', fname)
            print_SF_hist(tf, **vars(args))

if __name__ == '__main__':
    exit(main())
