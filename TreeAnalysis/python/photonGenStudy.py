#!/usr/bin/env python3

from os import path, mkdir
from copy import deepcopy
import ROOT
from ctypes import c_double


class TFileContext(object):
    def __init__(self, *args):
        self.tfile = ROOT.TFile(*args)
    def __enter__(self):
        return self.tfile
    def __exit__(self, type, value, traceback):
        self.tfile.Close()


def getPlots(inputdir, sample, plots):
    fname = path.join(inputdir, sample+".root")
    retrieved = []
    with TFileContext(fname) as rFile:
        for plot in plots:
            h = rFile.Get(plot)
            if(not h):
                print('Warning: Could not get "%s" from "%s"' % (plot, fname))
                retrieved.append(None)
            else:
                retrieved.append( deepcopy(h) )
    return retrieved


def integralTH1(h, direction=+1):
    if( (not h.Class().InheritsFrom('TH1')) or h.Class().InheritsFrom('TH2') or h.Class().InheritsFrom('TH3')):
        raise TypeError('Expected TH1')
    
    r = h.Clone(h.GetName() + "_integral")
    err = c_double(0.)
    # print('>>> Integral of', h.GetName())
    for b in range(0, h.GetNbinsX() + 2):
        if(direction >= 0):
            integral = h.IntegralAndError(b, -1, err)
        else:
            integral = h.IntegralAndError(0, b, err)
        # print('\t', b, integral, '+-', err.value)
        r.SetBinContent(b, integral)
        r.SetBinError(b, err.value)
    # print()
    return r
    

def transformTH1(h, f, ferr, label='transformed'):
    if( (not h.Class().InheritsFrom('TH1')) or h.Class().InheritsFrom('TH2') or h.Class().InheritsFrom('TH3')):
        raise TypeError('Expected TH1')

    r = h.Clone(h.GetName() + '_'+label)
    # print('>>> Transformation of', h.GetName())
    for b in range(0, h.GetNbinsX() + 2):
        r.SetBinContent( b, f(h.GetBinContent(b)) )
        r.SetBinError( b, ferr(h.GetBinContent(b)) * h.GetBinError(b) )
        # print('\t', b, h.GetBinContent(b), '+-', h.GetBinError(b), '-->', f(h.GetBinContent(b)), '+-', ferr(h.GetBinContent(b)) * h.GetBinError(b))
        
    return r


_rebinning = {
    'DRJet': 4
}
_xlabels = {
    'DRPhLep': r'\Delta\!R(\upgamma, \ell_{closest})'
}
_line_styles = (1, 2, 7)


def stackPlots(ZZpromptmO, ZZmatchedO, ZZnomatchO, ZZGpromptmO, ZZGmatchedO, ZZGnomatchO, variable, phID, region):
    # copy plots to avoid accidentaly modifying them
    ZZpromptm  = ZZpromptmO.Clone( ZZpromptmO.GetName() +'_copy')
    ZZmatched  = ZZmatchedO.Clone( ZZmatchedO.GetName() +'_copy')
    ZZnomatch  = ZZnomatchO.Clone( ZZnomatchO.GetName() +'_copy')
    ZZGpromptm = ZZGpromptmO.Clone(ZZGpromptmO.GetName()+'_copy')
    ZZGmatched = ZZGmatchedO.Clone(ZZGmatchedO.GetName()+'_copy')
    ZZGnomatch = ZZGnomatchO.Clone(ZZGnomatchO.GetName()+'_copy')
    
    def normalize_stack(stack):
        print('>>> nomalizing',stack.GetName())
        htot = stack.GetStack().Last().Clone('htot')
        for b in range(0, htot.GetNcells()): htot.SetBinError(b, 0)
        for h in stack.GetStack(): h.Divide(htot)
        stack.SetMaximum(1.)

    def fill_stack(h_promptm, h_matched, h_nomatch, sample, colors, do_normalize=True):
        stack = ROOT.THStack( '_'.join(('stack', sample, variable, phID.capitalize())), ' '.join((sample, variable, phID.capitalize())) )
        legend = ROOT.TLegend(0.7,0.7,0.87,0.87, '' ,'brNDC')
        _labels = ("match prompt", "match noprompt", "no match")
        for i,h in enumerate([h_promptm, h_matched, h_nomatch]):
            h.SetFillColor(colors[i])
            h.SetLineColor(ROOT.kBlack)
            h.SetLineWidth(2)
            stack.Add(h, 'hist')
            legend.AddEntry(h, _labels[i], 'lf')
            
        # for h in stack.GetStack():
        #     print('>>> name:', h.GetName(), '- title:', h.GetTitle(), '- x:', h.GetXaxis().GetTitle(), '- y:', h.GetYaxis().GetTitle())
        if(do_normalize):
            normalize_stack(stack)
        stack.Draw()
        legend.Draw("same")
        stack.GetXaxis().SetTitle(h_nomatch.GetXaxis().GetTitle())
        stack.GetYaxis().SetTitle('Fraction of events')
        return stack, legend  # Prevent garbage collection from destroying them

    do_normalize = variable != 'm4lG'
    
    canvas_stackZZ  = ROOT.TCanvas("canvas_stackZZ", "stackZZ", 1600, 900)
    canvas_stackZZ.cd()
    stackZZ, legendZZ   = fill_stack(ZZpromptm , ZZmatched , ZZnomatch , 'ZZ' , colors=(ROOT.kGreen -3, ROOT.kYellow+1, ROOT.kOrange+7), do_normalize=do_normalize)
    canvas_stackZZ.SaveAs(  path.join('Plot', 'photonGenStudy',  '{}_{}_{}_stackZZ.png'.format(variable, phID, region)) )
    
    canvas_stackZZG = ROOT.TCanvas("canvas_stackZZG", "stackZZG", 1600, 900)
    canvas_stackZZG.cd()
    stackZZG, legendZZG = fill_stack(ZZGpromptm, ZZGmatched, ZZGnomatch, 'ZZG', colors=(ROOT.kGreen -3, ROOT.kYellow+1, ROOT.kOrange+7), do_normalize=do_normalize)
    canvas_stackZZG.SaveAs( path.join('Plot', 'photonGenStudy', '{}_{}_{}_stackZZG.png'.format(variable, phID, region)) )


def myPlots(variable, phID, region):
    # basedir = 'rsync_results/Updated/EXT/2016/VVGammaAnalyzer_{}'.format(region)
    basedir = 'results/2016/VVGammaAnalyzer_{}'.format(region)
    name = 'PhGenStudy_{variable}_{phID}_{matched}'.format(variable=variable, phID=phID, matched='{}')
    ZZpromptm , ZZmatched , ZZnomatch  = getPlots(basedir, 'ZZTo4l'  , [name.format('matchedPrompt'), name.format('matched'), name.format('nomatch')] )
    ZZGpromptm, ZZGmatched, ZZGnomatch = getPlots(basedir, 'ZZGTo4LG', [name.format('matchedPrompt'), name.format('matched'), name.format('nomatch')] )
    
    if(ZZnomatch is None or ZZGnomatch is None):
        print('ERROR: one of "nomatch" plots is missing. Skipping', region, variable, phID)
        return
    if(ZZpromptm  is None):
        # PyRoot is not made to manipulate arrays: TAxis::GetXbins().GetArray() is simply unusable here!
        # ZZpromptm  = ROOT.TH1F('PhGenStudy_DRLep_loose_matchedPrompt',phID.capitalize(), ZZnomatch. GetNbinsX(), ZZnomatch. GetXaxis().GetXbins().GetArray())
        ZZpromptm = ZZnomatch.Clone('PhGenStudy_DRLep_loose_matchedPrompt')
        for b in range(0, ZZnomatch.GetNcells()+2):
            ZZpromptm.SetBinContent(b, 0)
            ZZpromptm.SetBinError(b, 0)
    if(ZZmatched  is None):
        ZZmatched = ZZnomatch.Clone('PhGenStudy_DRLep_loose_matched')
        for b in range(0, ZZnomatch.GetNcells()+2):
            ZZmatched.SetBinContent(b, 0)
            ZZmatched.SetBinError(b, 0)
    if(ZZGpromptm is None):
        ZZGpromptm = ZZGnomatch.Clone('PhGenStudy_DRLep_loose_matchedPrompt')
        for b in range(0, ZZGnomatch.GetNcells()+2):
            ZZGpromptm.SetBinContent(b, 0)
            ZZGpromptm.SetBinError(b, 0)
    if(ZZGmatched is None):
        ZZGmatched = ZZGnomatch.Clone('PhGenStudy_DRLep_loose_matched')
        for b in range(0, ZZGnomatch.GetNcells()+2):
            ZZGmatched.SetBinContent(b, 0)
            ZZGmatched.SetBinError(b, 0)
    
    for h in (ZZpromptm, ZZmatched, ZZnomatch, ZZGpromptm, ZZGmatched, ZZGnomatch):
        # print(h)
        h.SetTitle('{} photons (in {})'.format(phID.capitalize(), region))
        axis = h.GetYaxis()
        axis.SetMoreLogLabels()
        axis.SetNoExponent()
        h.Rebin(_rebinning.get(variable, 1))
        print('>>> rebinning', variable, _rebinning.get(variable, 1))

    ################################################################################
    # Stack Plots        
    stackPlots(ZZpromptm, ZZmatched, ZZnomatch, ZZGpromptm, ZZGmatched, ZZGnomatch, variable, phID, region)

    ################################################################################
    # Old Plots
    yMax  = max([max(h.GetBinContent(h.GetMaximumBin()), h.GetBinContent(h.GetNbinsX() + 1))
                 for h in (ZZpromptm, ZZmatched, ZZnomatch, ZZGpromptm, ZZGmatched, ZZGnomatch)])
    yMax *= 1.1
    yMin = 0 #10e-5 #max(ZZGnomatch.GetMinimum() * 0.1, 10e-4)
        
    canvas = ROOT.TCanvas('canvas', 'canvas', 1600, 900)
    canvas.cd()
    
    # Set max/min and draw axis
    frame = deepcopy(ZZnomatch)
    xlabel = _xlabels.get(variable, None)
    if(xlabel is not None):
        frame.GetXaxis().SetTitle(xlabel)
    frame.GetXaxis().SetRange(1, frame.GetXaxis().GetNbins()+1)
    frame.GetYaxis().SetRangeUser(yMin, yMax)
    frame.Draw("axis")
    
    # For some reason the title is not drawn
    title = '{} photons (in {}): {}'.format(phID.capitalize(), region, variable)
    titlePave = ROOT.TPaveText(.25, .92, .75, .98, 'NB NDC')
    titlePave.AddText(title)
    titlePave.SetFillColor(ROOT.kWhite)
    titlePave.Draw()

    def style_plots(h_promptm, h_matched, h_nomatch, legend, sample, color):
        h_promptm.SetMarkerStyle(ROOT.kFullSquare      )
        h_matched.SetMarkerStyle(ROOT.kFullTriangleUp  )
        h_nomatch.SetMarkerStyle(ROOT.kOpenTriangleDown)
        _labels = ("match prompt", "match noprompt", "no match")
        for i,h in enumerate([h_promptm, h_matched, h_nomatch]):
            h.SetLineColor(color)
            h.SetMarkerColor(color)
            h.SetLineStyle(_line_styles[i])
            legend.AddEntry(h, sample+' '+_labels[i], 'pl')
            h.Draw('CP hist same')
        
    legend = ROOT.TLegend(0.70,0.60,0.88,0.87, '' ,'brNDC')
    style_plots(ZZpromptm , ZZmatched , ZZnomatch , legend, sample='ZZ' , color=ROOT.kRed )
    style_plots(ZZGpromptm, ZZGmatched, ZZGnomatch, legend, sample='ZZG', color=ROOT.kBlue)
    legend.Draw('same')

    canvas.SaveAs( path.join('Plot', 'photonGenStudy', '{}_{}_{}_lin.png'.format(variable, phID, region)) )
    
    canvas.SetLogy()
    frame.GetYaxis().SetRangeUser(10e-5, yMax *4)
    canvas.Update()
    canvas.SaveAs( path.join('Plot', 'photonGenStudy', '{}_{}_{}_log.png'.format(variable, phID, region)) )
    
    
    ################################################################################
    ############################## Significance study ##############################
    ################################################################################
    
    ##### Sums #####
    canvas_sum = ROOT.TCanvas('canvas_sum', 'sum', 1600, 900)
    canvas_sum.cd()
    legend_sum = ROOT.TLegend(0.7,0.77,0.87,0.87, '' ,'brNDC')
    hZZ  = ZZpromptm  + ZZmatched  + ZZnomatch
    hZZG = ZZGpromptm + ZZGmatched + ZZGnomatch

    hZZG.SetLineColor(ROOT.kRed )
    hZZ .SetLineColor(ROOT.kBlue)
    hZZG.SetMarkerStyle(ROOT.kFullTriangleUp)
    hZZ .SetMarkerStyle(ROOT.kFullTriangleDown)
    hZZG.SetMarkerColor(hZZG.GetLineColor())
    hZZ .SetMarkerColor(hZZ .GetLineColor())
    legend_sum.AddEntry(hZZG, 'ZZ#gamma inclusive', 'lp')
    legend_sum.AddEntry(hZZ , 'ZZ   inclusive'    , 'lp')
    
    hZZ.Draw('pe')
    hZZG.Draw('pe same')
    legend_sum.Draw('same')

    canvas_sum.SetGrid()
    canvas_sum.SaveAs( path.join('Plot', 'photonGenStudy', '{}_{}_{}_sum.png'.format(variable, phID, region)) )

    ##### Integral #####
    # canvas_int = ROOT.TCanvas('canvas_int', 'integral', 1600, 900)
    # canvas_int.cd()
    # legend_int = ROOT.TLegend(0.7,0.77,0.87,0.87, '' ,'brNDC')

    hBplusS = integralTH1(hZZ) + integralTH1(hZZG)
    hSqrt = transformTH1(hBplusS, lambda x: x**0.5, lambda x: 1/(2*x**0.5), 'sqrt')
    # hSqrt.SetMarkerStyle(ROOT.kFullTriangleUp)
    # hBplusS.SetLineColor(ROOT.kMagenta+2)
    # hSqrt  .SetLineColor(ROOT.kGreen+2  )
    # hBplusS.SetMarkerColor(hBplusS.GetLineColor())
    # hSqrt  .SetMarkerColor(hSqrt  .GetLineColor())
    # legend_int.AddEntry(hBplusS, 'B+S'      , 'lp')
    # legend_int.AddEntry(hSqrt  , 'sqrt(B+S)', 'lp')

    # hBplusS.SetTitle('Integrals')
    # hBplusS.SetMinimum(0)
    # hBplusS.Draw('pe')
    # hSqrt.Draw('pe same')
    # legend_int.Draw('same')

    # canvas_int.SetGrid()
    # canvas_int.SaveAs( path.join('Plot', 'photonGenStudy', '{}_{}_{}_integral.png'.format(variable, phID, region)) )

    ##### Significance #####
    canvas_sig = ROOT.TCanvas('canvas_int', 'integral', 1600, 900)
    canvas_sig.cd()

    hSignif = integralTH1(hZZG) / hSqrt

    # print('>>> hSignif')
    # for b in range(0, hSignif.GetNbinsX() + 2):
    #     print(b, hSignif.GetBinContent(b), '+-', hSignif.GetBinError(b))

    hSignif.SetTitle(title + ' - S/#sqrt{S+B}')
    hSignif.SetLineColor(ROOT.kBlack)
    hSignif.SetMarkerStyle(ROOT.kFullTriangleUp)
    hSignif.SetMarkerColor(hSignif.GetLineColor())
    hSignif.SetMarkerSize(1.5)
    hSignif.SetMaximum(hSignif.GetBinContent(hSignif.GetMaximumBin()) * 1.5)
    hSignif.SetMinimum(0.)
    hSignif.Draw('PE hist C')

    canvas_sig.SetGrid()
    canvas_sig.SaveAs( path.join('Plot', 'photonGenStudy', '{}_{}_{}_significance.png'.format(variable, phID, region)) )


if __name__ == '__main__':
    try:
        mkdir('Plot/photonGenStudy')
    except FileExistsError:
        pass

    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch()
    
    # myPlots('DRLep', 'loose', 'SR4P')
    # myPlots('DRJet', 'loose', 'SR4P')
    for region in ['SR4P']:
        for variable in ['DRJet', 'DRPhLep', 'm4lG', 'nJets', 'ptPh']:
            for phID in ['kin', 'veryLoose', 'fail', 'loose']:
                myPlots(variable, phID, 'SR4P')