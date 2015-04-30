#! /usr/bin/env python

import ROOT, copy

from ROOT import TH1F


def write(particle,region,outname,fout):
    f = ROOT.TFile("results/VVXAnalyzer_MC/data.root")
    hn = f.Get("FakeRate_num_"+particle+"_"+region+"_pt")
    hd = f.Get("FakeRate_denom_"+particle+"_"+region+"_pt")

    
    fWZ = ROOT.TFile("results_1/VVXAnalyzer_MC/WZ.root")
    hnWZ = fWZ.Get("FakeRate_num_"+particle+"_"+region+"_pt")
    hdWZ = fWZ.Get("FakeRate_denom_"+particle+"_"+region+"_pt")

    hFake = hn.Clone("FakeRate_"+outname)
    hFake.Divide(hd)
    hFake.SetTitle("FakeRate_"+outname)
    fout.cd()
    hFake.Write("FakeRate_"+outname)

    hFake_NoWZ = hn.Clone("FakeRate_NoWZ_"+outname)
    hFake_NoWZ.Add(hnWZ,-1)

    hDen_NoWZ = hd.Clone("Den_NoWZ_"+outname)
    hDen_NoWZ.Add(hdWZ,-1)
    
    hFake_NoWZ.Divide(hDen_NoWZ)
    hFake_NoWZ.SetTitle("FakeRate_NoWZ_"+outname)
    hFake_NoWZ.Write("FakeRate_NoWZ_"+outname)


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
