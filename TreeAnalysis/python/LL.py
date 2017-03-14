#! /usr/bin/env python

import ROOT,copy
import CrossInfo
import math
from CrossInfo import* 
from ROOT import RooStats
#pragma link C++ class RooTreeData ;
#pragma link C++ class RooTreeData::PlotOpt ;
#pragma link C++ class RooTruthModel ;


def LL(hData,hDataUp,hDataDown,hMC,Type,FinState,UncType,DataOverMC):
 
    hRatio= ROOT.TH1F()
    hRatio_up= ROOT.TH1F() 
    hRatio_down= ROOT.TH1F()
    hRatio = copy.deepcopy(hData)
    hRatio_up =copy.deepcopy(hData) 
    hRatio_down = copy.deepcopy(hData) 

    NBins = hData.GetNbinsX() 
    for j in range(1,NBins+1):
        d = 0 
        sd = 0
        sd_stat = 0
        sd_syst = 0
        mc = 0 
        smc = 0
        errlo = 0
        errhi =0
        errstd = 0
        d = hData.GetBinContent(j)
        sd_stat = hData.GetBinError(j)
        sd_syst = (hDataUp.GetBinContent(j)-hDataDown.GetBinContent(j))/2#mean of up and down systematic uncertainties 
        if UncType == "stat": sd = sd_stat
        elif UncType == "syst": sd = sd_syst
        elif UncType == "statsyst": sd = math.sqrt(sd_stat*sd_stat+sd_syst*sd_syst)
        mc =  hMC.GetBinContent(j)
        smc = hMC.GetBinError(j)

        if DataOverMC == 1:
            r0 = d/mc 
            r = ROOT.RooRealVar("r","r",r0,0.,15.)
            x = ROOT.RooRealVar("x","x",mc*0.9,mc*1.1);
            #x = ROOT.RooRealVar("x","x",mc*0.5,mc*1.5);#FIXME!!!!!
            x0 = ROOT.RooRealVar("x0","x0",mc)
            sx = ROOT.RooRealVar("sx","sx",smc)
            y0 = ROOT.RooRealVar("y0","y0",d) 
            sy = ROOT.RooRealVar("sy","sy",sd) 
        
        else:#FIXME
            if d != 0: r0 = mc/d 
            else: r0 = 0
            r = ROOT.RooRealVar("r","r",r0,0.,15.)
            x = ROOT.RooRealVar("x","x",d*0.9,d*1.1);
            x0 = ROOT.RooRealVar("x0","x0",d)
            sx = ROOT.RooRealVar("sx","sx",sd)
            y0 = ROOT.RooRealVar("y0","y0",mc) 
            sy = ROOT.RooRealVar("sy","sy",smc) 

        rx = ROOT.RooProduct("rx","rx",ROOT.RooArgList(r,x))
        
        g1 = ROOT.RooGaussian("g1","g1",x,x0,sx)
        g2 = ROOT.RooGaussian("g2","g2",rx,y0,sy)
        #g2 = ROOT.RooPoisson("g2","g2",rx,y0) 
        
        LL = ROOT.RooProdPdf("LL","LL",g1,g2)
        
        obs = ROOT.RooArgSet(x0,y0) #observables
        poi = ROOT.RooArgSet(r) #parameters of interest
        data = ROOT.RooDataSet("data", "data", obs)
        data.add(obs) #actually add the data
        
        res = ROOT.RooFitResult(LL.fitTo(data,ROOT.RooFit.Minos(poi),ROOT.RooFit.Save(),ROOT.RooFit.Hesse(0)))
        
        #print "res status = "#, res.status()
        if res.status()==0:
            #r.Print()
            #x.Print()
            errlo =r.getErrorLo()
            errhi = r.getErrorHi()
        else: 
            print "Likelihood maximization failed"
            
        nll = LL.createNLL(data)
        frame = r.frame()
        pll = nll.createProfile(poi)
        #pll.plotOn(frame)#,RooFit::LineColor(ROOT::kRed))
        #frame.Draw()
 
        print j, r.getVal(), r0
        #hRatio.SetBinContent(j,r0)
 
        rval = r.getVal()
        r.setVal(0.)
        if errlo == 0 and pll.getVal()<0.5: errlo = -rval
        #else: break
        #print "sd",sd,"smc",smc,"mc",mc  
        stderr = r0*math.sqrt(((sd/d)*(sd/d)+(smc/mc)*(smc/mc)))#to compare with
        print "profile likelihood ratio results",j, errhi, errlo, stderr, stderr+errlo
      
        #hRatio_Errup.SetBinContent(j,math.fabs(stderr - errhi))
        #hRatio_Errdown.SetBinContent(j,math.fabs(stderr+errlo))
        hRatio.SetBinContent(j,r0)
        hRatio_up.SetBinContent(j,r0+errhi)
        hRatio_down.SetBinContent(j,r0+errlo)
        
    # hRatio.SetName(Type+"_"+FinState+"_"+"hRatio") 
    # hRatio_Errup.SetName(Type+"_"+FinState+"_"+"hRatio_Errup") 
    # hRatio_Errdown.SetName(Type+"_"+FinState+"_"+"hRatio_Errdown") 
    # hRatio.Write()
    # hRatio_Errup.Write()
    # hRatio_Errdown.Write()
    return hRatio, hRatio_up, hRatio_down

