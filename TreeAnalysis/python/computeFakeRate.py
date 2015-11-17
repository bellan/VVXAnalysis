#! /usr/bin/env python

import ROOT, copy

from ROOT import TH1F

def GetError(n, k, isUp, level=0.68540158589942957):

    alpha = (1.0 - level)/2
    if k<=0: xlow = 0.
    else: xlow = ROOT.Math.beta_quantile(alpha,k,n-k+1.)
    if k==n: xhigh = 1.
    else:  xhigh = ROOT.Math.beta_quantile(1. - alpha,k+1.,n-k)
    if isUp==1:
        #print "HighErr",xhigh
        return xhigh
    elif isUp==0:
        #print "LowErr",xlow
        return xlow

def GetGrEff(hpass, htotal):

    Nbin= hpass.GetNbinsX()
    for bin in  range(1,Nbin+2):

        if hpass.GetBinContent(bin)<0: hpass.SetBinContent(bin,0.)
        if htotal.GetBinContent(bin)<0: htotal.SetBinContent(bin,0.)
        #print  "bin",bin,"pass",hpass.GetBinContent(bin),"total", htotal.GetBinContent(bin)
    hpass_= copy.deepcopy(hpass)
    hpass_.Divide(htotal)
    for bin in range(1,Nbin+1):  
        if hpass_.GetBinContent(bin)==0:  hpass_.SetBinContent(bin,0)
    grEff = ROOT.TGraphAsymmErrors(hpass,htotal,"e0")

    for bin in range(1,Nbin+1):   
        #print "bin",bin,"pass",hpass.GetBinContent(bin),"total",htotal.GetBinContent(bin),"\nVal",hpass_.GetBinContent(bin) 

        XVal=ROOT.Double(0.)
        YVal=ROOT.Double(0.)
        grEff.GetPoint(bin-1,XVal,YVal)

        if bin==1 and YVal==0:
            grEff.SetPointEYhigh(bin-1, 0)
            grEff.SetPointEYlow(bin-1, 0)
        else:
            grEff.SetPointEYhigh(bin-1, GetError(htotal.GetBinContent(bin), hpass.GetBinContent(bin),1))
            grEff.SetPointEYlow(bin-1, GetError(htotal.GetBinContent(bin), hpass.GetBinContent(bin),0))

        grEff.GetPoint(bin-1,XVal,YVal)
        #print "Gr",XVal,YVal,grEff.GetErrorY(bin-1),grEff.GetErrorYhigh(bin-1)
        #print "\n"

        grEff.SetMarkerStyle(20)
        grEff.SetMarkerSize(.3)    
    return grEff



def write(particle,region,outname,fout):
    print "\n",particle,region
    f = ROOT.TFile("results/FakeRateAnalyzer_MC/data.root")
    hn = f.Get("FakeRate_num_"+particle+"_"+region+"_pt")
    hd = f.Get("FakeRate_denom_"+particle+"_"+region+"_pt")

    
    fWZ = ROOT.TFile("results/FakeRateAnalyzer_MC/WZ.root")
    hnWZ = fWZ.Get("FakeRate_num_"+particle+"_"+region+"_pt")
    hdWZ = fWZ.Get("FakeRate_denom_"+particle+"_"+region+"_pt")

    hFake =  copy.deepcopy(hn)
    hFake.Divide(hd)

    hFake.SetName("hFakeRate_"+outname)   
    hFake.SetTitle("hFakeRate_"+outname)
    
    grFake = GetGrEff(hn,hd)
    grFake.SetName("grFakeRate_"+outname)   
    grFake.SetTitle("grFakeRate_"+outname)
    

    hn.Add(hnWZ,-1)
    hd.Add(hdWZ,-1)

    hFake_NoWZ = copy.deepcopy(hn)
    hFake_NoWZ.Divide(hd)

    hFake_NoWZ.SetName("FakeRate_NoWZ_"+outname)
    hFake_NoWZ.SetTitle("FakeRate_NoWZ_"+outname)

       
    grFake_NoWZ = GetGrEff(hn,hd)
    grFake_NoWZ.SetName("grFakeRate_NoWZ_"+outname)
    grFake_NoWZ.SetTitle("grFakeRate_NoWZ_"+outname)

    fout.cd()

    hFake.Write(outname)
    grFake.Write("grFakeRate_"+outname)

    hFake_NoWZ.Write("NoWZ_"+outname)
    grFake_NoWZ.Write("grFakeRate_NoWZ_"+outname)

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

#write("electrons",'endcap','h1D_FRel_EE',fout)
#write("muons",'endcap','h1D_FRmu_EE',fout)
fout.Close()

