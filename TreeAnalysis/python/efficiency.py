#!/usr/bin/env python3

from __future__ import print_function
import sys, os
import ROOT

ROOT.gROOT.SetBatch(True)

class TFileContext(object):
    def __init__(self, *args):
        # print('>>>Opening file with args:', args)
        self.tfile = ROOT.TFile(*args)
    def __enter__(self):
        return self.tfile
    def __exit__(self, type, value, traceback):
        # print('<<<Closing TFile "%s"' % (self.tfile.GetName()))
        self.tfile.Close()

# input
config = {
    'results_folder': 'results',  # 'rsync_results/Updated/EXT'
    'analyzer': 'VVGammaAnalyzer',
    'region': 'SR2P',
    'year': 2016,
    'sample': 'ZZTo2Q2L'
}

filename = '{results_folder}/{year}/{analyzer}_{region}/{sample}.root'.format(**config)
assert os.path.exists(filename), 'Could not find "%s"' %(filename)

# output
outdir = 'Plot/Efficiencies/'
os.makedirs(outdir, exist_ok=True)

# (numerator, denominator)
effName = "Eff_{label:s}_{ND:s}_{var:s}"  # to be formatted
labels = [
    "AK4_quarks",
    "AK4_genAK4",
    "genAK4_quarks"
    # ("photons", "kin")
]
# effVars = {'pt':{}, 'E', 'eta', 'N'}
# resVars = {'dR', 'E', 'pt'}
# resVars2D = {'EvsE'}

c = ROOT.TCanvas("c", "The canvas", 0, 0, 1000, 1000)
c_debug = ROOT.TCanvas("c_debug", "Num and den", 0, 0, 1000, 1000)

def fix_negative_bins(hnum, hden):
    for b in range(0, hden.GetNbinsX()+1):
        num = hnum.GetBinContent(b)
        den = hden.GetBinContent(b)
        if(num > den):
            hnum.SetBinContent(b, den)


def efficiency(tfile, label, var, canvas, canvas_debug):
    effNameFormat = effName.format(label=label, var=var, ND='{ND:s}')
    effTitle = effNameFormat.replace('_{ND:s}', '')
    numerator   = effNameFormat.format(ND="NUM")
    denominator = effNameFormat.format(ND="DEN")
    
    print("###", numerator, "/", denominator, "###")
    hnum = tfile.Get(numerator)
    hden = tfile.Get(denominator)
    if(not hnum or not hden):
        if(not hnum): print('ERROR: could not get "{}" from "{}"'.format(numerator  , tfile.GetName()) )
        if(not hden): print('ERROR: could not get "{}" from "{}"'.format(denominator, tfile.GetName()) )
        return
           
    print(">>>bins   :  num:", hnum.GetNbinsX() , "- den:", hden.GetNbinsX())
    print(">>>entries:  num:", hnum.GetEntries(), "- den:", hden.GetEntries())
    
    fix_negative_bins(hnum, hden)
    
    c_debug.cd()
    hden.SetTitle( effTitle )
    hden.Draw("hist")
    hnum.SetLineColor(ROOT.kRed)
    hnum.Draw("same")
    c_debug.SaveAs(outdir+'debug_'+hden.GetTitle()+'.png')
    
    c.cd()
    geff = ROOT.TGraphAsymmErrors(hnum, hden, "")
    geff.SetTitle( effTitle )
    geff.GetYaxis().SetTitle('Efficiency')
    geff.SetLineColor(ROOT.kBlack)
    geff.SetMarkerStyle(ROOT.kFullTriangleUp)
    geff.SetMarkerColor(geff.GetLineColor())
    geff.Draw("APLE")
    geff.GetXaxis().SetTitle(hden.GetXaxis().GetTitle())
    
    c.SaveAs(outdir+geff.GetTitle()+'.png')
    

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


with TFileContext(filename) as tfile:
    for label in labels:
        for var in effVars:
           efficiency(tfile, label, var, c, c_debug)
        
        # for var in resVars:
        #     resolution(tfile, label, var, c)
        
        for var in resVars2D:
            resolution2D(tfile, label, var, c)
