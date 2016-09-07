#include <iostream>

#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooCBShape.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooGaussian.h>
#include <RooProduct.h>
#include <RooProdPdf.h>
#include <RooGenericPdf.h>
#include <RooLandau.h>

void LL(){

  //y0 = 0.000135096401209 sigma_y0 = 0.000103896581837 x0 = 0.000446013873443 sigma_x0 =1.81384394011e-06
  //0.014108652249 0.0168368471049 0.0219755396247 0.000120423865262 1.5575931164 1.55759310722 3.41637854038
  //0.072569437325 0.084063541977 0.0376693978906 0.000284216132439 0.51908074913 0.519080758095 1.12037749267
 // double d = 0.014108652249;
 //  double sd = 0.0168368471049;
 //  double mc = 0.0219755396247;
 //  double smc = 0.000120423865262;
 //  double r0 = d/mc;

  double d = 0.072569437325;
  double sd =  0.084063541977;
  double mc =  0.0376693978906;
  double smc =  0.00028421613243;
  double r0 = d/mc;

  RooRealVar x("x","x",mc*0.9,mc*1.1);
  RooRealVar x0("x0","x0",mc);
  RooRealVar sx("sx","sx",smc);

  RooRealVar r("r","r",r0,0.,5.);
  RooRealVar y0("y0","y0",d); 
  RooRealVar sy("sy","sy",sd); 
  
  RooProduct rx("rx","rx",RooArgList(r,x));

  RooGaussian g1("g1","g1",x,x0,sx);
  RooGaussian g2("g2","g2",rx,y0,sy);

  RooProdPdf LL("LL","LL",g1,g2);

  RooArgSet obs(x0,y0); //observables
  RooArgSet poi(r); //parameters of interest
  RooDataSet data("data", "data", obs);
  data.add(obs); //actually add the data


  RooFitResult* res = LL.fitTo(data,RooFit::Minos(poi),RooFit::Save(),RooFit::Hesse(false));
  if(res->status()==0) {
    r.Print();
    x.Print();
    cout << r.getErrorLo() << " " << r.getErrorHi() << endl;
  } else {
    cout << "Likelihood maximization failed" << endl;
  }
  
  RooAbsReal* nll = LL.createNLL(data); 
  RooPlot* frame = r.frame();
  RooAbsReal* pll = nll->createProfile(poi);
  pll->plotOn(frame);//,RooFit::LineColor(ROOT::kRed));
  frame->Draw();

  r.setVal(0.);
  cout << pll->getVal() << endl; 

  return;
    
    


}
