#! /usr/bin/env python

import ROOT, copy

from ROOT import TH1F


def write(particle,region,outname,fout):
    f = ROOT.TFile("results_1/VVXAnalyzer_MC/data.root")
    hn = f.Get("FakeRate_num_"+particle+"_"+region+"_pt")
    hd = f.Get("FakeRate_denom_"+particle+"_"+region+"_pt")
    hn.Divide(hd)
    hn.SetTitle(outname)
    fout.cd()
    hn.Write(outname)


particles = ['muons','electrons']
regions   = ['barrel','endcap']

fout = ROOT.TFile("fakeRates.root", "RECREATE")

for particle in particles:
    for region in regions:
        p = ''
        r = ''
        if particle == 'muons': p = 'mu'
        if particle == 'electrons': p = 'el'
        if region   == 'barrel': r = 'B'
        if region   == 'endcap': r = 'E'

        write(particle,region,'h1D_FR'+p+'_E'+r,fout)

fout.Close()
