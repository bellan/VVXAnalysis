#!/usr/bin/env python3
from argparse import ArgumentParser
import logging
import ROOT
import cmsstyle

class TFileContext(object):
    def __init__(self, *args):
        self.tfile = ROOT.TFile(*args)
        if(not (self.tfile and self.tfile.IsOpen())):
            raise FileNotFoundError(args[0] if len(args) > 0 else '')

    def __enter__(self):
        return self.tfile

    def __exit__(self, exc_type, exc_value, traceback):
        self.tfile.Close()


def main(args):
    logging.debug('args = %s', args)
    plot_name = 'covariance_fit_%s' %(args.fit_name)
    logging.debug('plot_name = %s', plot_name)

    # Get correlation matrix
    with TFileContext(args.filename) as fin:
        # Get the TH2 with the correlation matrix
        h = fin.Get(plot_name)
        h.SetDirectory(0)
    logging.debug('h = %s', h)

    # Style options and canvas
    ROOT.gROOT.SetBatch(True)
    cmsstyle.SetExtraText(args.extra_text)
    cmsstyle.SetLumi(args.lumi)
    logging.debug('gStyle = %s', ROOT.gStyle.GetName())

    canv = cmsstyle.cmsCanvas(
        'canv',
        **get_canvas_params(h),
        square=False,
        extraSpace=.1,
        iPos=0,
        with_z_axis = True
    )
    # NOTE: the canvas draws numeric axes, so the TH2 (alphanumeric) must be overwrite them (e.g. NOT use "same").
    # Then the lumi and CMS logo must be re-drawn on top of it.

    # Draw
    canv.cd()
    ROOT.gStyle.SetPaintTextFormat('.2f')
    ROOT.gStyle.SetPalette(ROOT.kViridis)
    h.SetMarkerSize(.7)
    label_size = get_max_label_size(h)
    logging.debug('label_size: %f', label_size)
    h.GetXaxis().SetLabelSize(label_size)
    h.GetYaxis().SetLabelSize(label_size)
    h.GetZaxis().SetLabelSize(.03)
    logging.debug('x label size: %f', h.GetXaxis().GetLabelSize())
    logging.debug('y label size: %f', h.GetYaxis().GetLabelSize())
    h.Draw('colz text')

    # Re-draw the lumi
    cmsstyle.CMS_lumi(canv, iPosX=0)

    # Move the palette (z axis) to a reasonable location...
    palette = h.GetListOfFunctions().FindObject('palette')
    x_palette = .9
    y_palette = .2
    palette.SetX1NDC(x_palette)
    palette.SetX2NDC(x_palette + 0.0345)
    palette.SetY1NDC(y_palette)
    palette.SetY2NDC(y_palette + 0.7)
    ROOT.gPad.Modified()
    ROOT.gPad.Update()
    logging.debug('palette = %s', palette)
    logging.debug('palette x1: %f  x2: %f (%f) - y1: %f  y2: %f (%f)',
                  palette.GetX1NDC(),
                  palette.GetX2NDC(),
                  palette.GetX2NDC() - palette.GetX1NDC(),
                  palette.GetY1NDC(),
                  palette.GetY2NDC(),
                  palette.GetY2NDC() - palette.GetY1NDC()
                  )

    # Save canvas
    for ext in ('pdf', 'png'):
        canv.SaveAs(plot_name+'.'+ext)


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('filename', metavar='FILE', help='ROOT file that contains the fit (as TH2)')
    parser.add_argument('-l', '--lumi', default='138')
    parser.add_argument('-t', '--extra-text', default='Preliminary')
    parser.add_argument(      '--b-only', action='store_const', dest='fit_name', const='b', default='s', help='Plot the background-only fit')
    parser.add_argument(      '--sig'   , action='store_const', dest='fit_name', const='s', help='Plot the fit with the signal (default)')
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')
    return parser.parse_args()


def get_canvas_params(h):
    '''
    Set the axis limits according to the plot limits
    '''
    xaxis = h.GetXaxis()
    yaxis = h.GetYaxis()
    return {
        'x_min': xaxis.GetBinLowEdge(1),
        'x_max': xaxis.GetBinLowEdge(xaxis.GetNbins()),
        'y_min': yaxis.GetBinLowEdge(1),
        'y_max': yaxis.GetBinLowEdge(yaxis.GetNbins()),
        'nameXaxis': xaxis.GetTitle(),
        'nameYaxis': yaxis.GetTitle(),
    }


def get_max_label_size(h):
    space = 0.69
    max_len_x = max(len(l) for l in h.GetXaxis().GetLabels())

    return space / max_len_x


if __name__ == '__main__':
    args = parse_args()
    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    exit(main(args))
