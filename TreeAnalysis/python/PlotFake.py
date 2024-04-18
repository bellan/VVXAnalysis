#! /usr/bin/env python
from __future__ import print_function
import ROOT, copy
import recoPlotUtils
import sys
from ROOT import gSystem, TCanvas, TH1,  TPad, gStyle, TLegend, THStack, TGraphAsymmErrors,Math, TF1
from Colours import *

def plot(plot,obs):

    c1 = TCanvas( 'c1', 'c1', 200, 10, 900, 700 )
    
    f = ROOT.TFile("results/FakeRateAnalyzer_MC/data.root")
    h = f.Get(plot)
    h.Rebin(4)
    Max = h.GetMaximum()+ h.GetMaximum()/3.

    stack = ROOT.THStack("stack_num","Stack_"+plot)   
    hRatio = copy.deepcopy(h)
    leg = TLegend(0.81,0.77,0.95,0.96)

#    mclist = [{"sample":'WZTo3LNu',"color":ROOT.kOrange-4,"name":'WZ'},{"sample":'TTTo2L2Nu2B',"color":ROOT.kRed-2,"name":'tt'},{"sample":'TTJets',"color":ROOT.kRed-2,"name":'tt'},{"sample":'DYJetsToLLTuneCM50',"color":ROOT.kGreen-5,"name":'DY'}]

    mclist = [{"sample":'WZ',"color":ROOT.kOrange-4,"name":'WZ'},{"sample":'TTTo2L2Nu',"color":ROOT.kRed-2,"name":'tt'},{"sample":'DYJetsToLL_M50',"color":ROOT.kGreen-5,"name":'DY'}]

    LastColor=""
    for sample in mclist:
        print(sample["sample"])
        fMC = ROOT.TFile("results/FakeRateAnalyzer_MC/"+sample["sample"]+".root")
        h_mc_or = fMC.Get(plot)
        h_mc = copy.deepcopy(h_mc_or)
        h_mc.Rebin(4)
#       h_mc.SetNameTitle(sample["name"])
        h_mc.SetFillColor(sample["color"])           
        h_mc.SetLineColor(sample["color"])
        if LastColor!=sample["color"]:
            leg.AddEntry(h_mc,sample["name"],"f")
        LastColor=sample["color"]
        stack.Add(h_mc)

    leg.AddEntry(h,"data","lpe")
    stack.SetMaximum(Max)

    fInt = ROOT.TF1("constant","1",h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax());

    hRatio.SetStats(0)
    hRatio.Divide(stack.GetStack().Last())

    pad1 = ROOT.TPad ('hist', '', 0., 0.30, 1.0, 1.0)
    pad1.SetTopMargin (0.10)
    pad1.SetRightMargin (0.10)
    pad1.SetLeftMargin (0.10)
    pad1.Draw()

    c1.cd()
    pad2 = ROOT.TPad ('rat', 'Data/MC ratio', 0., 0.0,  1., 0.3)
    pad2.SetTopMargin (0.10)
    pad2.SetRightMargin (0.10)
    pad2.SetLeftMargin (0.10)
    pad2.Draw()
    
    pad1.cd()

    stack.Draw("hist")
    stack.GetHistogram().GetXaxis().SetTitle(obs) 
    h.SetMarkerStyle(20)
    h.Draw("sameE1")
    leg.Draw("same")

    pad2.cd()
    hRatio.SetMarkerStyle(20)
    hRatio.Draw("E1")
    fInt.Draw("same")


    c1.SaveAs("~/www/PlotsVV/13TeV/FakeRate/"+plot+".png")        
    pad1.cd()
#    c1.cd()
    pad1.SetLogy()    
    c1.SaveAs("~/www/PlotsVV/13TeV/FakeRate/"+plot+"_Log.png")        

    return (h,stack)

def createHistos(plot):
    f = ROOT.TFile("fakeRates.root")
    f_mc = ROOT.TFile("fakeRates_MC.root")
    h_ = f.Get(plot)
    h = copy.deepcopy(h_)
   # Max = h.GetMaximum()+ h.GetMaximum()/5.

    stack = ROOT.THStack("stack_num","Stack_"+plot)   
    leg = TLegend(0.81,0.77,0.95,0.96)

    h_mc_or = f_mc.Get(plot)
    h_mc = copy.deepcopy(h_mc_or)
    leg.AddEntry(h,"data","lpe")
    leg.AddEntry(h_mc,"MC","lpe")
    h.SetMaximum(0.5)

    return(h,h_mc,leg)

def Plot(List,obs):

    #h1= ROOT.TH1F()
    h1=copy.deepcopy(List[0])
    h2=List[1]
    leg=List[2]

    c1 = TCanvas( 'c1', 'c1', 200, 10, 900, 700 )
    fInt = ROOT.TF1("constant","1",h1.GetXaxis().GetXmin(),h1.GetXaxis().GetXmax());

    hRatio = copy.deepcopy(h2)   
   
    if "Stack" in type(h1).__name__:
        hRatio.Divide(h1.GetStack().Last())
        hRatio.SetStats(0)  
    elif "TH1" in type(h1).__name__:
        hRatio.Divide(h1)
        hRatio.SetStats(0)  

    if "TGraph" not in type(h1).__name__: 
        pad1 = ROOT.TPad ('hist', '', 0., 0.30, 1.0, 1.0)
        pad1.SetTopMargin (0.10)
        pad1.SetRightMargin (0.10)
        pad1.SetLeftMargin (0.10)
        
        c1.cd()
        pad2 = ROOT.TPad ('rat', 'Data/MC ratio', 0., 0.0,  1., 0.3)
        pad2.SetTopMargin (0.10)
        pad2.SetRightMargin (0.10)
        pad2.SetLeftMargin (0.10)
    
    if "TH1" in type(h1).__name__:
        h1.Draw("hist")
        h1.GetXaxis().SetTitle(obs) 
        pad1.Draw()
        pad2.Draw() 
        h2.Draw("sameE1")

    elif "TGraph" in type(h1).__name__: 
        h1.SetMarkerSize(1)
        h1.SetMarkerStyle(20)
        h1.Draw("AP")
        h1.GetHistogram().GetYaxis().SetTitle(obs) 
        h2.SetMarkerSize(1)
        h2.SetMarkerStyle(21)
        h2.SetMarkerColor(2)
        h1.GetHistogram().GetXaxis().SetTitle("p_{T}") 
#        h2.GetXaxis().SetTitle("p_{T}")
        h2.Draw("sameP")

    else: sys.exit(type(h1).__name__+" is not contemplate")
    
    if "TGraph" not in type(h1).__name__: 
        pad1.cd()
        h2.SetMarkerStyle(20)

    leg.Draw("same")

    if "TGraph" not in type(h1).__name__: 
        pad2.cd()
        hRatio.SetMarkerStyle(20)
        hRatio.Draw("E1")
        fInt.Draw("same")

    c1.SaveAs("~/www/PlotsVV/13TeV/FakeRate/"+h1.GetName()+".png")        

#    c1.cd()
    if "TGraph" not in type(h1).__name__: 
        pad1.cd()        
        pad1.SetLogy()    
        c1.SaveAs("~/www/PlotsVV/13TeV/FakeRate/"+h1.GetName()+"_Log.png")        

particles = ['muons','electrons']
regions   = ['barrel','endcap']


for particle in particles:
    p = ''
    if particle == 'muons': p = 'mu'
    elif particle == 'electrons': p = 'el'
#    plot("MT_"+p,"m_{T} GeV")    
 #   plot("MTNoHF_"+p,"m_{T} GeV")    

    for region in regions:
        print(particle,region)
        r = ''
        if region   == 'barrel': r = 'B'
        if region   == 'endcap': r = 'E'


        #Plot(createHistos("NoWZ_h1D_FR"+p+"_E"+r),"fr")
        Plot(createHistos("grFakeRate_NoWZ_h1D_FR"+p+"_E"+r),"fr")

        #plot("FakeRate_num_"+particle+"_"+region+"_pt","p_{T} GeV")
        #plot("FakeRate_denom_"+particle+"_"+region+"_pt","p_{T} GeV")

        #plot("ZL_"+particle+"_"+region+"_Iso","Iso")
        #plot("ZL_"+particle+"_"+region+"_Iso","sip")


