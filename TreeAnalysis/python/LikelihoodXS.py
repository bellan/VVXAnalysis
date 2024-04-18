#! /usr/bin/env python

from __future__ import print_function
import ROOT
import CrossInfo
from CrossInfo import* 
from ROOT import RooStats
import copy

## to time the macro

def LikeCross(wspace):
           
    #  make Poisson model * Gaussian constraint
    wspace.factory("prod:s_4m(mu[1.,0.,1000],xs_4m,eff_4m,Lumi)")  
    wspace.factory("prod:s_2e2m(mu[1.,0.,1000],xs_2e2m,eff_2e2m,Lumi)")  
    wspace.factory("prod:s_4e(mu[1.,0.,1000],xs_4e,eff_4e,Lumi)")  
        
    wspace.factory("sum:nexp_4m(s_4m,bkg_4m,irrBkg_4m)")
    wspace.factory("sum:nexp_2e2m(s_2e2m,bkg_2e2m, irrBkg_2e2m)")
    wspace.factory("sum:nexp_4e(s_4e,bkg_4e,irrBkg_4e)")
 
    wspace.Print()
    
    #Poisson of (n | s+b)
    wspace.factory("Poisson:pdf_4m(nobs_4m,nexp_4m)")
    wspace.factory("Poisson:pdf_2e2m(nobs_2e2m,nexp_2e2m)")
    wspace.factory("Poisson:pdf_4e(nobs_4e,nexp_4e)")
    
    wspace.factory("Gaussian:constraint_4m(bkg0_4m,bkg_4m,sigmaBkg_4m)")
    wspace.factory("Gaussian:constraint_2e2m(bkg0_2e2m,bkg_2e2m,sigmaBkg_2e2m)")
    wspace.factory("Gaussian:constraint_4e(bkg0_4e,bkg_4e,sigmaBkg_4e)")

    wspace.factory("Gaussian:constraintIrrBkg_4m(irrBkg0_4m,irrBkg_4m,sigmaIrrBkg_4m)")
    wspace.factory("Gaussian:constraintIrrBkg_2e2m(irrBkg0_2e2m,irrBkg_2e2m,sigmaIrrBkg_2e2m)")
    wspace.factory("Gaussian:constraintIrrBkg_4e(irrBkg0_4e,irrBkg_4e,sigmaIrrBkg_4e)")
    
    wspace.factory("Gaussian:constraintEff_4m(eff0_4m,eff_4m,sigmaEff_4m)")
    wspace.factory("Gaussian:constraintEff_2e2m(eff0_2e2m,eff_2e2m,sigmaEff_2e2m)")
    wspace.factory("Gaussian:constraintEff_4e(eff0_4e,eff_4e,sigmaEff_4e)")
    
    wspace.factory("PROD:model_4m(pdf_4m,constraint_4m,constraintIrrBkg_4m,constraintEff_4m)")
    wspace.factory("PROD:model_2e2m(pdf_2e2m,constraint_2e2m,constraintIrrBkg_2e2m,constraintEff_2e2m)")
    wspace.factory("PROD:model_4e(pdf_4e,constraint_4e,constraintIrrBkg_4e,constraintEff_4e)")
    
    wspace.factory("PROD:model(model_4m,model_2e2m,model_4e)")
            
    BR_4m   = BRmu*BRmu
    BR_2e2m = 2*BRmu*BRele
    BR_4e   = BRele*BRele    
    BR      = BR_4m+BR_2e2m+BR_4e
    
    print("BR sum",BR," BRmean ",(BR_4m + BR_4e + BR_2e2m)/3.," BR_4m ",BR_4m," BR_2e2m ",BR_2e2m," BR_4e ",BR_4e)
    
    # wspace.var("BR_4m").setVal(BR_4m)
    # wspace.var("BR_2e2m").setVal(BR_2e2m)
    # wspace.var("BR_4e").setVal(BR_4e)
    
    wspace.var("xs_4m").setConstant(True) # needed for being treated as global observables
    wspace.var("xs_2e2m").setConstant(True) # needed for being treated as global observables
    wspace.var("xs_4e").setConstant(True) # needed for being treated as global observables
    
    wspace.var("Lumi").setConstant(True) 
           
    wspace.defineSet("obs","nobs_4m,nobs_4e,nobs_2e2m")  #observables
    wspace.defineSet("poi","mu"); #parameters of interest
    wspace.defineSet("np","bkg_4m,bkg_4e,bkg_2e2m,irrBkg_4m,irrBkg_4e,irrBkg_2e2m,eff_4m,eff_4e,eff_2e2m") #nuisance_lumi
        
    wspace.Print()
    
    
  # make data set with the namber of observed events
    print("---Make data Set---")       
    data = ROOT.RooDataSet("data","", ROOT.RooArgSet(wspace.var("nobs_4m"),wspace.var("nobs_2e2m"),wspace.var("nobs_4e")));
     
    data.add(ROOT.RooArgSet(wspace.var("nobs_4m") ))
    data.add(ROOT.RooArgSet(wspace.var("nobs_2e2m") ))
    data.add(ROOT.RooArgSet(wspace.var("nobs_4e") ))

    print("---Make Frame---")       
#    frame = wspace.var("mu").frame(0.95,1.1);
    frame = wspace.var("mu").frame(1.1,1.25);
    
    print("---Make nll---")       
    nll = wspace.pdf("model").createNLL(data)
    print("Initial NLL value",nll.getVal())
    nll.plotOn(frame,ROOT.RooFit.ShiftToZero()); #the ShiftToZero option puts the minimum at 0 on the y-axis

    #nllcopy = copy.deepcopy(nll)

    #nllcopy.plotOn(frame,ROOT.RooFit.ShiftToZero()); #the ShiftToZero option puts the minimum at 0 on the y-axis
    #   nllcopy.plotOn(frame)
    mini = ROOT.RooMinuit(nll);
    mini.minos(wspace.set("poi"));

    print("NLL value",nll.getVal())       
    print("---Create profile---")   

    pll = nll.createProfile(wspace.set("poi"))
    pll.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kRed))
    
    frame.SetMaximum(1.) 
    c = ROOT.TCanvas( 'c', 'c', 200, 10, 900, 700 )
    frame.Draw()
    c.SaveAs("Plot/likelihood.png")
    c.SaveAs("~/www/PlotsVV/13TeV/likelihood.png")
    mu = wspace.var("mu").getValV()
    mu_hi = wspace.var("mu").getErrorHi()
    mu_lo = wspace.var("mu").getErrorLo()

    #   print "mu",mu,"cross",mu*34.34,"+",(mu+mu_hi)*34.34,"-",(mu+mu_lo)*34.34
    #   print "Total Cross ",(mu*34.34)/(0.5338*BR)

    TotThCross = xs_tight["4l"] #xs_2e2m +xs_4m + xs_4e
#    Acc = xs_4e/xs_OS_4e #FIXME 
    Acc = xs_tight["4e"]/xs_OS_4e #FIXME 
    print("Acc",Acc)
    #print "mu",mu,"cross",mu*TotThCross,"+",(mu+mu_hi)*TotThCross,"-",(mu+mu_lo)*TotThCross

    print("mu",mu,"cross",mu*TotThCross,"+",(mu_hi)*TotThCross,(mu_lo)*TotThCross)
    print("Total Cross ",(mu*TotThCross)/(Acc*BR))

    return mu*TotThCross




