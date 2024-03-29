#!/usr/bin/env python

######################################################################
# Get lepton fake rates from ZZAnalysis (as TGraphs)                 #
# and store them in a TH2D (x:|eta|, y:pt)                           #
# One file per year, which contains two TH2Ds (eletrons and muons)   #
#                                                                    #
# Author: A. Mecca (alberto.mecca@cern.ch)                           #
######################################################################

from __future__ import print_function
import os
from array import array
from ctypes import c_double
from argparse import ArgumentParser

import ROOT

class TFileContext(object):
    def __init__(self, *args):
        self.tfile = ROOT.TFile(*args)
    def __enter__(self):
        return self.tfile
    def __exit__(self, type, value, traceback):
        self.tfile.Close()

_allowed_years = (2016, 2017, 2018)

parser = ArgumentParser()
parser.add_argument('-y', '--years', nargs='+', choices=_allowed_years, default=_allowed_years, type=int)
parser.add_argument('-p', '--print-graphs', dest='do_print'    , action='store_true')
parser.add_argument(      '--no-root'     , dest='do_rootfiles', action='store_false', default=True)
parser.add_argument(      '--no-title'    , dest='do_title'    ,action='store_false', help='Do not paint the title')
args = parser.parse_args()

if(args.do_print):
    ROOT.gROOT.SetBatch(True)

# Locate the directory in which VVXAnalysis and ZZAnalysis are
if(os.environ.get('CMSSW_BASE', False)):
    basepath = os.path.join(environ['CMSSW_BASE'], 'src')
else:
    try:
        cwd_split = os.getcwd().split(os.sep)
        idx = cwd_split.index('VVXAnalysis')
        basepath = '/'.join(cwd_split[:idx])
    except ValueError:
        print('ERROR: could not locate top directory of VVXAnalysis')
        exit(2)

# Check that the input dir exists
ZZpath = os.path.join(basepath, 'ZZAnalysis')
assert os.path.isdir(ZZpath), 'Could not determine the locaton of ZZAnalysis'
fr_dir = os.path.join(ZZpath, 'AnalysisStep/data/FakeRates')
assert os.path.isdir(fr_dir), 'FakeRates folder not found in ZZAnalysis'

# Check that the output dir exists
out_dir = os.path.join(basepath, 'VVXAnalysis', 'Commons', 'data')
assert os.path.isdir(out_dir), 'Cannot find VVXAnalysis/Commons/data to store output rootfiles'


def print_graphs(fin, year):
    ele_EB = fin.Get("FR_OS_electron_EB")
    ele_EE = fin.Get("FR_OS_electron_EE")
    muo_EB = fin.Get("FR_OS_muon_EB")
    muo_EE = fin.Get("FR_OS_muon_EE")
    if(not (ele_EB and ele_EE and muo_EB and muo_EE)):
        print('Could not get all 4 FakeRate histograms from "{}"!'.format(fin.GetName()))
        return

    c = ROOT.TCanvas('c', 'canvas', 1600, 1200)
    c.SetRightMargin(0.05)
    xmin = ymin = 100
    xmax = ymax = -100
    for g in (ele_EE, ele_EB, muo_EB, muo_EE):
        xmin_ = c_double(0.)
        ymin_ = c_double(0.)
        xmax_ = c_double(0.)
        ymax_ = c_double(0.)
        g.ComputeRange(xmin_, ymin_, xmax_, ymax_)
        xmin = min(xmin, xmin_.value)
        ymin = min(ymin, ymin_.value)
        xmax = max(xmax, xmax_.value)
        ymax = max(ymax, ymax_.value)
    yrange = ymax - ymin
    assert yrange > 0
    ymin_draw = ymin - 0.1*yrange
    ymax_draw = ymax + 0.1*yrange

    ele_EB.SetMinimum(0)  # (ymin_draw)
    ele_EB.SetMaximum(.16)# (ymax_draw)

    legend = ROOT.TLegend(*list(map(lambda x,y: x+y, [0.3, 0.7]*2, [0, 0, .22, .18])))

    ele_EB.SetMarkerStyle(ROOT.kFullSquare)
    ele_EE.SetMarkerStyle(ROOT.kOpenSquare)
    muo_EB.SetMarkerStyle(ROOT.kFullCircle)
    muo_EE.SetMarkerStyle(ROOT.kOpenCircle)
    ele_EB.SetMarkerSize(1.2)
    ele_EE.SetMarkerSize(1.2)
    muo_EB.SetMarkerSize(1.2)
    muo_EE.SetMarkerSize(1.2)
    ele_EB.SetLineStyle(1)
    ele_EE.SetLineStyle(1)
    muo_EB.SetLineStyle(1)
    muo_EE.SetLineStyle(1)
    ele_EB.SetLineColor(ROOT.kRed)
    ele_EE.SetLineColor(ROOT.kRed)
    muo_EB.SetLineColor(ROOT.kBlue)
    muo_EE.SetLineColor(ROOT.kBlue)
    ele_EB.SetMarkerColor(ele_EB.GetLineColor())
    ele_EE.SetMarkerColor(ele_EE.GetLineColor())
    muo_EB.SetMarkerColor(muo_EB.GetLineColor())
    muo_EE.SetMarkerColor(muo_EE.GetLineColor())

    ele_EB.SetTitle('Lepton fake rates - {}'.format(year) if args.do_title else '')
    ele_EB.GetXaxis().SetTitle('p_{T} [GeV/c]')
    ele_EB.GetYaxis().SetTitle('Fake Rate')
    legend.AddEntry(ele_EB, 'electron barrel', 'lp')
    legend.AddEntry(ele_EE, 'electron endcap', 'lp')
    legend.AddEntry(muo_EB, 'muon barrel'    , 'lp')
    legend.AddEntry(muo_EE, 'muon endcap'    , 'lp')
    ele_EB.Draw('APE')
    ele_EE.Draw('PE')
    muo_EB.Draw('PE')
    muo_EE.Draw('PE')

    legend.Draw('same')
    c.SaveAs('leptonFakeRate_{}.eps'.format(year))

def write_new_file(fout, fin):
    ele_EB = fin.Get("FR_OS_electron_EB")
    ele_EE = fin.Get("FR_OS_electron_EE")
    muo_EB = fin.Get("FR_OS_muon_EB")
    muo_EE = fin.Get("FR_OS_muon_EE")
    if(not (ele_EB and ele_EE and muo_EB and muo_EE)):
        print('Could not get all 4 FakeRate histograms from "{}"!'.format(fin.GetName()))
        return

    print(fin.GetName(), '-->', fout.GetName())

    # Create the TH2s with the bin edges corresponding to points +- errors from the TGraphs
    def get_xedges(tgraph):
        nx = tgraph.GetN()
        points_x  = array('d', tgraph.GetX() )
        points_ex = array('d', tgraph.GetEX())

        xedges = array('d')

        for x, ex in zip(points_x, points_ex):
            print(x, '+-', ex)
            xedges.append(x - ex)
        xedges.append(points_x[-1] + points_ex[-1])
        return nx, xedges

    ny_e, yedges_e = get_xedges(ele_EB)
    ny_m, yedges_m = get_xedges(muo_EB)
    nx_e, xedges_e = 2, array('d', [0., 1.459, 2.5])
    nx_m, xedges_m = 2, array('d', [0., 1.2, 2.4])

    fout.cd()
    th2_e = ROOT.TH2D("fakeRate_e", "Electron fake rate scale factor;|#eta|;p_{T} [GeV/c]", nx_e, xedges_e, ny_e, yedges_e)
    th2_m = ROOT.TH2D("fakeRate_m",     "Muon fake rate scale factor;|#eta|;p_{T} [GeV/c]", nx_m, xedges_m, ny_m, yedges_m)

    # Fill the TH2s
    def fill_th2(th2, tgr_eb, tgr_ee):
        values_EB = array('d', tgr_eb.GetY() )
        values_EE = array('d', tgr_ee.GetY() )
        errors_EB = array('d', tgr_eb.GetEY())
        errors_EE = array('d', tgr_ee.GetEY())

        for by in range(1, th2.GetNbinsY()+1):
            # Bin #1 (in the th2) is the FIRST --> index #0 of the array (from the TGraph)
            print('>>> by:', by)
            b_EB = th2.GetBin(1, by)
            th2.SetBinContent(b_EB, values_EB[by-1])
            th2.SetBinError  (b_EB, errors_EB[by-1])
            print( '\tb_EB: {:2d}  y: {:.3f}  ey: {:.3f}'.format(b_EB, values_EB[by-1], errors_EE[by-1]) )

            b_EE = th2.GetBin(2, by)
            th2.SetBinContent(b_EE, values_EE[by-1])
            th2.SetBinError  (b_EE, errors_EE[by-1])
            print( '\tb_EE: {:2d}  y: {:.3f}  ey: {:.3f}'.format(b_EE, values_EE[by-1], errors_EE[by-1]) )

    fill_th2(th2_e, ele_EB, ele_EE)
    fill_th2(th2_m, muo_EB, muo_EE)

    th2_e.Write()
    th2_m.Write()


# Loop over the years
for year in args.years:
    fr_fname = 'newData_FakeRates_OS_{}.root'.format(year)
    fr_fpath = os.path.join(fr_dir, fr_fname)
    if(not os.path.isfile(fr_fpath)):
        print('ERROR: Cannot not find "{}" in "{}"'.format(fr_fname, fr_dir))
        continue

    out_fname = 'leptonFakeRates_{}.root'.format(year)
    out_fpath = os.path.join(out_dir, out_fname)

    if(args.do_print):
        with TFileContext(fr_fpath, 'READ') as fin:
            print_graphs(fin, year)

    if(args.do_rootfiles):
        with TFileContext(fr_fpath, 'READ') as fin, TFileContext(out_fpath, 'RECREATE') as fout:
            write_new_file(fout, fin)
