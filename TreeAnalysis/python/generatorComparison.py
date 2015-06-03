#! /usr/bin/env python

import ROOT, copy, sys, ast

from optparse import OptionParser

from tdrstyle import *

inputdir = './'

def getPlot(filename,plotName, color, style):
    f = ROOT.TFile(inputdir+filename)
    h = f.Get(plotName)
    h.SetLineColor(color)
    h.SetMarkerColor(color)
    h.SetMarkerStyle(style)
    h.SetMarkerSize(0.5)
    return copy.deepcopy(h)




def plot(plotName, log = False):
    hMGNLO  = getPlot("mgnlo.root",plotName,ROOT.kRed, 20)
    hPowheg = getPlot("powheg.root",plotName,ROOT.kBlue,24)
    hMGLO   = getPlot("mglo.root",plotName,8,21)

    legend = ROOT.TLegend(0.54,0.78,0.79,0.91)
    if 'Z0mass' in plotName or 'Z1mass' in plotName: legend = ROOT.TLegend(0.20,0.78,0.45,0.91)

    legend.AddEntry(hMGNLO,"MadGraph (NLO)","lp")
    legend.AddEntry(hMGLO,"MadGraph (LO)","lp")
    legend.AddEntry(hPowheg,"Powheg (NLO)","lp")

    legend.SetTextSize(0.05)
    legend.SetTextFont(42)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetLineColor(ROOT.kWhite)
    legend.SetShadowColor(ROOT.kWhite)
    
    

    c1 = ROOT.TCanvas(plotName,plotName)
    c1.Divide(2,2)
    c1.cd(1)
    hMGNLO.Draw()
    hMGLO.Draw("same")
    hPowheg.Draw("same")
    legend.Draw("same")


    line = ROOT.TLine(hMGNLO.GetXaxis().GetXmin(),1,hMGNLO.GetXaxis().GetXmax(),1)
    line.SetLineColor(ROOT.kRed)

    l = ROOT.TLatex()
    l.SetNDC()
    l.SetTextSize(0.05)
    l.SetTextFont(42)
    

    # MGNLO/Powheg
    c1.cd(3)
    hRatio_MGNLO_Powheg = hMGNLO.Clone("ratio_MGNLO_Powheg")
    hRatio_MGNLO_Powheg.Divide(hPowheg)
    hRatio_MGNLO_Powheg.SetLineColor(ROOT.kBlack) 
    hRatio_MGNLO_Powheg.SetMarkerColor(ROOT.kBlack) 
    hRatio_MGNLO_Powheg.SetMarkerStyle(21) 
    hRatio_MGNLO_Powheg.GetYaxis().SetRangeUser(0,5)
    hRatio_MGNLO_Powheg.Draw()
    line.Draw("same")
    l.DrawLatex(0.25,0.75,"MadGraph (NLO) / Powheg")

    # MGNLO/MGLO
    c1.cd(2)
    hRatio_MGNLO_MGLO = hMGNLO.Clone("ratio_MGNLO_MGLO")
    hRatio_MGNLO_MGLO.Divide(hMGLO)
    hRatio_MGNLO_MGLO.SetLineColor(ROOT.kBlack) 
    hRatio_MGNLO_MGLO.SetMarkerColor(ROOT.kBlack) 
    hRatio_MGNLO_MGLO.SetMarkerStyle(21) 
    hRatio_MGNLO_MGLO.GetYaxis().SetRangeUser(0,5)
    hRatio_MGNLO_MGLO.Draw()
    line.Draw("same")
    l.DrawLatex(0.25,0.75,"MadGraph (NLO) / MadGrap (LO)")

    # Powheg/MGLO
    c1.cd(4)
    hRatio_Powheg_MGLO = hPowheg.Clone("ratio_Powheg_MGLO")
    hRatio_Powheg_MGLO.Divide(hMGLO)
    hRatio_Powheg_MGLO.SetLineColor(ROOT.kBlack) 
    hRatio_Powheg_MGLO.SetMarkerColor(ROOT.kBlack) 
    hRatio_Powheg_MGLO.SetMarkerStyle(21) 
    hRatio_Powheg_MGLO.GetYaxis().SetRangeUser(0,5)
    hRatio_Powheg_MGLO.Draw()
    line.Draw("same")
    l.DrawLatex(0.25,0.75,"Powheg / MadGraph (LO)")


    c1.cd(1)
    maxYtmp = hMGLO.GetMaximum() if hMGLO.GetMaximum() > hPowheg.GetMaximum() else hPowheg.GetMaximum()
    maxY = hMGNLO.GetMaximum() if hMGNLO.GetMaximum() > maxYtmp else maxYtmp
    hMGNLO.GetYaxis().SetRangeUser(0,maxY*1.10)
    if log: 
        minYtmp = hMGLO.GetMinimum() if hMGLO.GetMinimum() > hPowheg.GetMinimum() else hPowheg.GetMinimum()
        minY = hMGNLO.GetMinimum() if hMGNLO.GetMinimum() > minYtmp else minYtmp
        hMGNLO.GetYaxis().SetRangeUser(minY if not minY == 0 else 0.0001, maxY*1.20)
        c1.GetPad(1).SetLogy()

    c1.Update()
    c1.SaveAs('{0:s}{1:s}.png'.format(plotName, '_log' if log else ''))
    #input()

if __name__ == '__main__':
    setTDRStyle()

    parser = OptionParser(usage="usage: %prog <plot name> [options]")
    parser.add_option("-l", "--log", dest="log",
                  action="store_true",
                  help="Use this option if you want the comparison plot in log scale")

    (options, args) = parser.parse_args()
    plotName = args[0] 
    log = False if options.log is None else options.log

    if plotName == "all": 
        f = ROOT.TFile(inputdir+'mgnlo.root')
        
        next = ROOT.TIter(f.GetListOfKeys())
        key = next()
        while key:
            plot(key.GetName(), log)
            key = next()

    else:
        plot(plotName, log)

