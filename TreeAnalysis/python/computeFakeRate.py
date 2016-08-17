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
    
    hpass_= copy.deepcopy(hpass)
    hpass_.Divide(htotal)
    for bin in range(1,Nbin+1):  
        if hpass_.GetBinContent(bin)==0:  hpass_.SetBinContent(bin,0)
    grEff = ROOT.TGraphAsymmErrors(hpass,htotal,"e0")

    for bin in range(1,Nbin+1):   
        #print "bin ",bin," passed events ",hpass.GetBinContent(bin),"total events ",htotal.GetBinContent(bin),"\nRatio ",hpass_.GetBinContent(bin) 

        XVal=ROOT.Double(0.)
        YVal=ROOT.Double(0.)
        grEff.GetPoint(bin-1,XVal,YVal)

        if bin==1 and YVal==0:
            grEff.SetPointEYhigh(bin-1, 0)
            grEff.SetPointEYlow(bin-1, 0)
        else:
            grEff.SetPointEYhigh(bin-1,-YVal+GetError(htotal.GetBinContent(bin), hpass.GetBinContent(bin),1))
            grEff.SetPointEYlow(bin-1, YVal-GetError(htotal.GetBinContent(bin), hpass.GetBinContent(bin),0))

        grEff.GetPoint(bin-1,XVal,YVal)
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
    
    Nbin= hn.GetNbinsX()

    print "\nData before WZ subtraction\n"
    for bin in range(1,Nbin+1):
        print "bin",bin
        print "pass",hn.GetBinContent(bin)
        print "total",hd.GetBinContent(bin)

    print "\nWZ\n"
    for bin in range(1,Nbin+1):
        print "bin",bin
        print "pass",hnWZ.GetBinContent(bin)
        print "total",hdWZ.GetBinContent(bin)
        

    # To avoid problems such passed bin content > total bin content
    hn.SetBinContent(Nbin+1,0)
    hd.SetBinContent(Nbin+1,0) 

    hn.Add(hnWZ,-1)
    hd.Add(hdWZ,-1)

    print hn.Integral(0,-1),hd.Integral(0,-1)

    print "\ndata after WZ subtraction\n"
    for bin in range(1,Nbin+1):
        print "bin",bin
        print "pass",hn.GetBinContent(bin)
        print "total",hd.GetBinContent(bin)


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


def write_MC(particle,region,outname,fout):
    print "\n",particle,region

    mclist = [{"sample":'TTJets',"color":ROOT.kRed-2,"name":'tt'},{"sample":'DYJetsToLL_M50',"color":ROOT.kGreen-5,"name":'DY'}]

    stack_n = ROOT.THStack("stack_num","Stack_"+"FakeRate_num_"+particle+"_"+region+"_pt")   
    stack_d = ROOT.THStack("stack_den","Stack_"+"FakeRate_denom_"+particle+"_"+region+"_pt")   
    
    for sample in mclist:
        print sample["sample"]
        fMC = ROOT.TFile("results/FakeRateAnalyzer_MC/"+sample["sample"]+".root")
        h_mc_n = fMC.Get("FakeRate_num_"+particle+"_"+region+"_pt")
        h_mc_d = fMC.Get("FakeRate_denom_"+particle+"_"+region+"_pt")
        #h_mc = copy.deepcopy(h_mc_or)
        stack_n.Add(h_mc_n)
        stack_d.Add(h_mc_d)

    hFake =  copy.deepcopy(stack_n.GetStack().Last())
    hFake_d = copy.deepcopy(stack_d.GetStack().Last())

    grFake = GetGrEff(hFake,hFake_d)

    hFake.Divide(hFake_d)

    Nbin= hFake_d.GetNbinsX()

    # print "\ndata first\n"
    # for bin in range(1,Nbin+2):
    #     print "bin",bin
    #     print "pass",hn.GetBinContent(bin)
    #     print "total",hd.GetBinContent(bin)


    # print "\nWZ\n"
    # for bin in range(1,Nbin+2):
    #     print "bin",bin
    #     print "pass",hnWZ.GetBinContent(bin)
    #     print "total",hdWZ.GetBinContent(bin)        

    hFake.SetName("FakeRate_NoWZ_"+outname)
    hFake.SetTitle("FakeRate_NoWZ_"+outname)
       
    grFake.SetName("grFakeRate_NoWZ_"+outname)
    grFake.SetTitle("grFakeRate_NoWZ_"+outname)

    fout.cd()

    hFake.Write("NoWZ_"+outname)
    grFake.Write("grFakeRate_NoWZ_"+outname)

particles = ['muons','electrons']
regions   = ['barrel','endcap']

fout = ROOT.TFile("fakeRates.root", "RECREATE")
fout_MC = ROOT.TFile("fakeRates_MC.root", "RECREATE")


for particle in particles:
    for region in regions:
        p = ''
        r = ''
        if particle == 'muons': p = 'mu'
        if particle == 'electrons': p = 'el'
        if region   == 'barrel': r = 'B'
        if region   == 'endcap': r = 'E'
        write(particle,region,'h1D_FR'+p+'_E'+r,fout)
        write_MC(particle,region,'h1D_FR'+p+'_E'+r,fout_MC)

#write("electrons",'endcap','h1D_FRel_EE',fout)
#write("muons",'endcap','h1D_FRmu_EE',fout)
fout.Close()

