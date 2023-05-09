#!/usr/bin/env python3

from os import path
import copy
from ctypes import c_double, c_int
import ROOT

class TFileContext(object):
    def __init__(self, *args):
        self.tfile = ROOT.TFile(*args)
    def __enter__(self):
        return self.tfile
    def __exit__(self, type, value, traceback):
        self.tfile.Close()


def getPlots(inputdir, sample, plots, verbose=0):
    fname = path.join(inputdir, sample+".root")
    if(not path.exists(fname)):
        print('WARN: file "{}" does not exist'.format(fname))
        return [None for plot in plots]

    retrieved = []
    with TFileContext(fname) as rFile:
        for plot in plots:
            h = rFile.Get(plot)
            if(not h):
                if(verbose > 0): print('WARN: Could not get "%s" from "%s"' % (plot, fname))
                retrieved.append(None)
            else:
                retrieved.append( copy.deepcopy(h) )
                if(verbose > 1):
                    ignore = c_double(0)
                    n = h.IntegralAndError(0, -1, 0, -1, ignore)
                    print('\t\t{:s} - entries: {:6.0f}'.format(plot, n))
    return retrieved


def getSlices(h2, name=None, direction='X'):
    if(not h2.Class().InheritsFrom("TH2")):
        print('ERROR: "{:s}" does not inherit from TH2'.format(h2.GetName()))
        return
    if(name is None):
        name = h2.GetName()

    if  (direction=='Y'):
        axis_2D_proj = h2.GetYaxis()
        # axis_2D_draw = h2.GetXaxis()
        nbins_proj = h2.GetNbinsY()
        # nbins_draw = h2.GetNbinsX()
        proj = h2.ProjectionX
    elif(direction=='X'):
        axis_2D_proj = h2.GetXaxis()
        # axis_2D_draw = h2.GetYaxis()
        nbins_proj = h2.GetNbinsX()
        # nbins_draw = h2.GetNbinsY()
        proj = h2.ProjectionY

    h1s = []

    for j in range(0, nbins_proj+2):
        if(j == 0):
            binTitle = '%s < {}'.format(axis_2D_proj.GetBinUpEdge(j))
        elif(j == nbins_proj+1):
            binTitle = '{} > %s'.format(axis_2D_proj.GetBinLowEdge(j))
        else:
            binTitle = '{} < %s < {}'.format(axis_2D_proj.GetBinLowEdge(j), axis_2D_proj.GetBinUpEdge(j))
        binTitle = binTitle % (axis_2D_proj.GetTitle())
        # print(f'\t{j=} - {binTitle=}')
        h1 = proj(h2.GetName()+'_proj{}_{:d}'.format(direction, j), j, j, 'e')  # TODO projectionY
        h1.SetTitle(' - '.join([h2.GetTitle(), binTitle]))
        h1s.append(h1)

    return h1s


def computeRange(pad):
    assert pad.Class().InheritsFrom('TVirtualPad')

    xmin, ymin, xmax, ymax = c_double(0), c_double(0), c_double(0), c_double(0)
    pad.GetRangeAxis(xmin, ymin, xmax, ymax)

    return xmin.value, ymin.value, xmax.value, ymax.value
