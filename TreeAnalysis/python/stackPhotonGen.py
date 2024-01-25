#!/usr/bin/env python

from argparse import ArgumentParser
import logging
import ROOT
from plotUtils23 import TFileContext

def stackGenPhoton(tf, hname):
    outname = (hname %('')).strip('_')
    canvas = ROOT.TCanvas('canvas', outname, 1600, 900)
    canvas.SetRightMargin(0.05)
    
    h_promptm = tf.Get(hname %('promptm'))
    h_matched = tf.Get(hname %('matched'))
    h_nomatch = tf.Get(hname %('nomatch'))

    style={
        'promptm': {'color': ROOT.kGreen },
        'matched': {'color': ROOT.kYellow},
        'nomatch': {'color': ROOT.kOrange}
    }

    h_promptm.SetFillColor(style['promptm']['color'])
    h_matched.SetFillColor(style['matched']['color'])
    h_nomatch.SetFillColor(style['nomatch']['color'])
    for h in [h_promptm, h_matched, h_nomatch]:
        h.SetLineColor(ROOT.kBlack)

    stack  = ROOT.THStack('stack', ';%s;%s' %(h_promptm.GetXaxis().GetTitle(), h_promptm.GetYaxis().GetTitle()))
    legend = ROOT.TLegend(.78, .73, .93, .88)
    stack.Add(h_promptm)
    stack.Add(h_matched)
    stack.Add(h_nomatch)
    legend.AddEntry(h_promptm, 'promptm')
    legend.AddEntry(h_matched, 'matched')
    legend.AddEntry(h_nomatch, 'nomatch')
    stack.Draw('hist')
    legend.Draw('same')

    for ext in ('png', 'pdf'):
        canvas.SaveAs(outname+'.'+ext)

def main():
    parser = ArgumentParser()

    parser.add_argument('fname', metavar='FILE')
    parser.add_argument(      '--wp' , choices=('kin', 'veryLoose', 'loose'), default='kin')
    parser.add_argument('-t', '--var', choices=('DRLep', 'DRJet', 'FPJet', 'DPJet'), default='FPJet')
    parser.add_argument(      '--list', action='store_true', dest='do_list')
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')    
    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s',
                        level=args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel))

    hname = "PhGenStudy_%s_%s_%s" %(args.var, args.wp, '%s')
    logging.info('hname: %s', hname)
    with TFileContext(args.fname) as tf:
        if(args.do_list):
            print('\n'.join(k.GetName() for k in tf.GetListOfKeys() if k.GetName().startswith('PhGenStudy')))
            return 0
        stackGenPhoton(tf, hname=hname)
    return 0

if __name__ == '__main__':
    exit(main())
