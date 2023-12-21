#!/usr/bin/env python

###############################################################
# Read a 2D histogram from a JSON, create a TH2F from it,     #
# and save it to a rootfile.                                  #
#                                                             #
# Author: A. Mecca (alberto.mecca@cern.ch)                    #
# Date: 21/12/2023                                            #
###############################################################


import os
import json
from array import array
from argparse import ArgumentParser
import logging

import ROOT

from plotUtils23 import TFileContext

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("inputfile"        , metavar="FILE")
    parser.add_argument("-o", "--output"   , metavar="FILE", default=None          , help="Specify output file. Overrides --outputdir")
    parser.add_argument("-O", "--outputdir", metavar="DIR" , default=os.path.curdir, help="Output directory")
    parser.add_argument(      "--mode"     , default="RECREATE", help="Option passed to TFile's constructor. Choices: [CREATE, RECREATE, UPDATE]. Default: %(default)s")
    parser.add_argument(      "--histname" , default="FR", help="Histogram name in the output file")
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')
    args = parser.parse_args()

    if(args.output is None):
        basename = os.path.basename(args.inputfile)
        basename_noext = os.path.splitext(basename)[0]
        args.output = os.path.join(args.outputdir, basename_noext+'.root')

    return args

def main():
    args = parse_args()

    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    ROOT.gROOT.SetBatch(True)

    with open(args.inputfile) as f:
        json_data = json.load(f)

    xbins = array('d', json_data["xbins"])
    ybins = array('d', json_data["ybins"])
    logging.debug('xbins: %s', xbins)
    logging.debug('ybins: %s', ybins)
    nbinsx = len(xbins)-1
    nbinsy = len(ybins)-1
    h = ROOT.TH2F(args.histname, ";%s;%s"%(json_data["xvar"], json_data["yvar"]),
                  len(xbins)-1, xbins,
                  len(ybins)-1, ybins
                  )

    # Note: the source JSONs were written with the row order inverted
    for j,row in enumerate(json_data["values"]):
        for i,cell in enumerate(row):
            # logging.debug("%d-%d: %.2f +- %.2f", i, j, *cell)
            bx = i+1
            by = nbinsx-j+1
            b = h.GetBin(bx, by)

            # The source has only 2 significant digits
            if(cell[1] <= 0):
                cell[1] = 0.005
                logging.warning("Bin error <= 0 in bin %d (%d, %d) -> %s: [%g, %g]  %s: [%g, %g]",
                                b, bx, by,
                                json_data["xvar"], h.GetXaxis().GetBinLowEdge(bx), h.GetXaxis().GetBinUpEdge(bx),
                                json_data["yvar"], h.GetYaxis().GetBinLowEdge(by), h.GetYaxis().GetBinUpEdge(by)
                                )

            h.SetBinContent(b, cell[0])
            h.SetBinError  (b, cell[1])

    with TFileContext(args.output, args.mode) as tf:
        h.Write()

    logging.info('Output in: %s', args.output)


if __name__ == "__main__":
    exit(main())
