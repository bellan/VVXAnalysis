#!/usr/bin/env python3

from os import path
from argparse import ArgumentParser
from array import array

import ROOT
import plotUtils3  # getPlots, computeRange


if __name__ == '__main__':
    ROOT.gROOT.SetBatch(True)
    parser = ArgumentParser()

    parser.add_argument('-i', '--inputdir', default='results')
    parser.add_argument('-y', '--year'    , default=2018)
    parser.add_argument('-r', '--region'  , default='SR4P')
    parser.add_argument('-A', '--analyzer', default='VVGammaAnalyzer')
    parser.add_argument('-s', '--sample'  , default='ZZGTo4LG')
    parser.add_argument(      '--logy'    , action='store_true')
    parser.add_argument('-v', '--verbose' , dest='verbosity', action='count', default=1)
    parser.add_argument('-q', '--quiet'   , dest='verbosity', action='store_const', const=0)

    args = parser.parse_args()

    inpath = path.join(args.inputdir, str(args.year), '_'.join([args.analyzer, args.region]))
    outpath = path.join('Plot', 'slices', '_'.join(args.inputdir.split('_')[1:]))
    print(f'{inpath=}')
    print(f'{outpath=}')

    cprof = ROOT.TCanvas('cprof', 'Profile', 2000, 1500)
    c_all = ROOT.TCanvas('c_all', 'All profiles', 2000, 1500)

    plots = plotUtils3.getPlots(inpath, args.sample, ('lead_loose_dRl_pt',), args.verbosity)

    _style = [
        {"color":ROOT.kRed     , "marker":ROOT.kFullCircle},
        {"color":ROOT.kBlue    , "marker":ROOT.kFullTriangleUp},
        {"color":ROOT.kGreen   , "marker":ROOT.kFullTriangleDown},
        {"color":ROOT.kBlack   , "marker":ROOT.kFullSquare},
        {"color":ROOT.kMagenta , "marker":ROOT.kFullCrossX},
        {"color":ROOT.kCyan    , "marker":ROOT.kFullStar},
        {"color":ROOT.kYellow+3, "marker":ROOT.kFullCross},
        {"color":ROOT.kOrange+4, "marker":ROOT.kOpenTriangleUp},
        {"color":ROOT.kViolet-5, "marker":ROOT.kOpenTriangleDown},
        {"color":ROOT.kTeal+3  , "marker":ROOT.kOpenSquare},
        {"color":ROOT.kSpring+3, "marker":ROOT.kOpenCircle},
        {"color":ROOT.kMagenta+4,"marker":ROOT.kOpenCross}
    ]

    for h2 in plots:
        # Individual bins
        cprof.cd()
        min_X, max_X, min_Y, max_Y = 1, 0, 1, 0

        ROOT.gStyle.SetOptStat('euio')

        slices_X = plotUtils3.getSlices(h2, direction='X')
        for h1 in slices_X:
            # cprof.Clear()
            # h1.Draw('hist')
            # cprof.SaveAs( path.join(outpath, h1.GetName()+'.png') )
            minval = h1.GetMinimum()  # h1.GetBinContent(h1.GetMinimumBin())
            if(minval != 0):
                min_X = min(min_X, minval)
            max_X = max(max_X, h1.GetMaximum())

        slices_Y = plotUtils3.getSlices(h2, direction='Y')
        for h1 in slices_Y:
            # cprof.Clear()
            # h1.Draw('hist')
            # cprof.SaveAs( path.join(outpath, h1.GetName()+'.png') )
            minval = h1.GetMinimum()  # h1.GetBinContent(h1.GetMinimumBin())
            if(minval != 0):
                min_Y = min(min_X, minval)
            max_Y = max(max_Y, h1.GetMaximum())


        # Canvas with all the plots
        print(f'>>> {min_X=}, {max_X=}, {min_Y=}, {max_Y=}')
        c_all.cd()
        if(args.logy):
            c_all.SetLogy()
            scale_X = 10**(0.0012*(max_X/min_X))
            scale_Y = 10**(0.0012*(max_Y/min_Y))
            print(f'>>> {scale_X=} - {scale_Y=}')
        else:
            scale_X = scale_Y = 1.6

        ROOT.gStyle.SetOptStat("")

        legend_w = 0.2
        legend_dh = 0.025

        c_all.Clear()
        legend_X = ROOT.TLegend(legend_w, legend_dh*len(slices_X))
        for i, s in enumerate(slices_X[:len(_style)]):
            s.SetMaximum(max_X * scale_X)
            s.SetMinimum(min_X / scale_X)
            style = _style[i]
            s.SetMarkerStyle(style['marker'])
            s.SetMarkerColor(style['color' ])
            s.SetLineColor  (style['color' ])
            s.GetYaxis().SetMoreLogLabels()
            s.GetYaxis().SetNoExponent()
            legend_X.AddEntry(s, s.GetTitle().split(' - ')[-1])
            s.SetTitle(h2.GetTitle())            
            s.Draw('L same' if i != 0 else 'L')
            s.Draw('E0 same')

        # xmin, ymin, xmax, ymax = plotUtils3.computeRange(c_all)
        # print(f'>>> {xmin=}, {ymin=}, {xmax=}, {ymax=}')

        legend_X.Draw('same')

        c_all.SaveAs(path.join(outpath, h2.GetName()+'_allX.png'))


        legend_Y = ROOT.TLegend(legend_w, min(legend_dh*len(slices_Y), 0.25))
        for i, s in enumerate(slices_Y[:len(_style)]):
            s.SetMaximum(max_Y * scale_Y)
            s.SetMinimum(min_Y / scale_Y)
            style = _style[i]
            s.SetMarkerStyle(style['marker'])
            s.SetMarkerColor(style['color' ])
            s.SetLineColor  (style['color' ])
            s.GetYaxis().SetMoreLogLabels()
            s.GetYaxis().SetNoExponent()
            legend_Y.AddEntry(s, s.GetTitle().split(' - ')[-1])
            s.SetTitle(h2.GetTitle())
            s.Draw('E0 same'  if i != 0 else 'E0')

        legend_Y.Draw('same')

        c_all.SaveAs(path.join(outpath, h2.GetName()+'_allY.png'))

