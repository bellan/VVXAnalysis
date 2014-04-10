#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraph.h"

#include <iostream>
#include <vector>

using namespace std;


void MacroCompare() {
  
  gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginZZWGenAnalysisPlugins.so");
  
  // on LXPLUS:
  gSystem->AddIncludePath(" -I$CMSSW_BASE/src/ -I$CMSSW_RELEASE_BASE/src");
  
  // on laptop
  //gSystem->AddIncludePath(" -I../../../../ -I$CMSSW_RELEASE_BASE/src");
  
  std::string path = std::string(gROOT->GetMacroPath())+"$CMSSW_BASE/src:$CMSSW_BASE/src/VVXAnalysis/Producers/test/GenAnalysis";
  gROOT->SetMacroPath(path.c_str());
  
  
  gROOT->LoadMacro("H6f.h+");
  gROOT->LoadMacro("Hbos.h+");
  gROOT->LoadMacro("Hjets.h+");

  
  TFile* fQED = new TFile("QED6_SIGNAL.root");
  TFile* fWZZ = new TFile("WZZ_SIGNAL.root");
  
    
  //  gDirectory->cd("Rint:/");

  TCanvas* c1 =  new TCanvas("Signal(2f->6f) - Signal(WZZ on-shell) Bosons Mass ","Signal(2f->6f) - Signal(WZZ on-shell) Bosons Mass",900,700);
  // c1->Divide(2,1);
  
  TLegend* L1 = new TLegend(0.67, 0.61,0.88,0.76);
  
  
  float L = 300.;
  
  float xsec_QED6 = 52.76;
  float xsec_Wp = 0.02985;
  float xsec_Wm = 0.01525;
  
  int Ngen_QED6 = 9998*5;
  int Ngen_Wp = 99989;
  int Ngen_Wm = 9997;

  float w_QED6 = (L*xsec_QED6)/Ngen_QED6;
  float w_WZZ = (L*(xsec_Wp+xsec_Wm))/(Ngen_Wp+Ngen_Wm);
 
  cout << "\n%%%%%%%%%%%%\nWEIGHTS\nw_QED6= " << w_QED6 << "\nw_WZZ= " << w_WZZ << endl; 


  //=================== 6f =========================

  TH1F* m6f_QED = new TH1F("6fMass_QED", "6fMass_QED", 3000, 0, 3000);
  m6f_QED = (TH1F*)fQED->Get("myAnalyzer/all6fMass_0");
  TH1F* m6f_WZZ = new TH1F("6fMass_WZZ", "6fMass_WZZ", 3000, 0, 3000);
  m6f_WZZ = (TH1F*)fWZZ->Get("myAnalyzer/all6fMass");

  m6f_QED->Scale(w_QED6);
  m6f_WZZ->Scale(w_WZZ);
 
  m6f_WZZ->SetLineColor(kRed);

  float errQED = sqrt(m6f_QED->GetEntries())*w_QED6;
  float errWZZ = sqrt(m6f_WZZ->GetEntries())*w_WZZ;
  
  float Integral_WZZ = m6f_WZZ->Integral(0, m6f_WZZ->GetNbinsX()+1);
  float Integral_QED = m6f_QED->Integral(0, m6f_QED->GetNbinsX()+1);
  float scaling = Integral_WZZ/Integral_QED;

  cout << "\nINTEGRAL_WZZ= "  << Integral_WZZ  << "\tError_WZZ= "    << errWZZ << endl;
  cout << "\nINTEGRAL_QED6= " << Integral_QED << "\t\tError_QED6= " << errQED << endl;
  cout << "scaling= " << scaling << endl;

  m6f_QED->Scale(scaling);
 
  cout << "\nINTEGRAL_WZZ= "  << Integral_WZZ  << endl;
  cout << "\nINTEGRAL_QED6= " << m6f_QED->Integral(0, m6f_QED->GetNbinsX()+1) << endl;

  c1->cd(0);
  m6f_QED->SetTitle("6fMass");
  m6f_QED->Rebin(2);
  m6f_QED->DrawClone();
  m6f_WZZ->DrawClone("same");


  L1->AddEntry(m6f_QED, "2f -> 6f", "l");
  L1->AddEntry(m6f_WZZ, "WZZ", "l");
  L1->Draw();

  TStyle::gStyle->SetOptStat(0);


//   //=================== 4l =========================


//   TH1F* Mass4l_0 = new TH1F("4lMass0", "4lMass0", 3000, 0, 3000);
//   Mass4l_0 = (TH1F*)f0->Get("myAnalyzer/all4lMass");
//   TH1F* Mass4l_1 = new TH1F("4lMass1", "4lMass1", 3000, 0, 3000);
//   Mass4l_1 = (TH1F*)f1->Get("myAnalyzer/all4lMass");
//   TH1F* Mass4l_2 = new TH1F("4lMass2", "4lMass2", 3000, 0, 3000);
//   Mass4l_2 = (TH1F*)f2->Get("myAnalyzer/all4lMass");
//   TH1F* Mass4l_3 = new TH1F("4lMass3", "4lMass3", 3000, 0, 3000);
//   Mass4l_3 = (TH1F*)f3->Get("myAnalyzer/all4lMass");
//   TH1F* Mass4l_4 = new TH1F("4lMass4", "4lMass4", 3000, 0, 3000);
//   Mass4l_4 = (TH1F*)f4->Get("myAnalyzer/all4lMass");
  
//   TH1F* Mass4l_Wp = new TH1F("4lMassWp", "4lMassWp", 3000, 0, 3000);
//   Mass4l_Wp = (TH1F*)fWp->Get("myAnalyzer/all4lMass");
//   TH1F* Mass4l_Wm = new TH1F("4lMassWm", "4lMassWm", 3000, 0, 3000);
//   Mass4l_Wm = (TH1F*)fWm->Get("myAnalyzer/all4lMass");


//   Mass4l_0->Scale(w_QED6);
//   Mass4l_1->Scale(w_QED6);
//   Mass4l_2->Scale(w_QED6);
//   Mass4l_3->Scale(w_QED6);
//   Mass4l_4->Scale(w_QED6);
//   Mass4l_Wp->Scale(w_QED6);
//   Mass4l_Wm->Scale(w_QED6);  
//   Mass4l_Wp->SetLineColor(kRed);
//   Mass4l_Wm->SetLineColor(kRed);
 
//   Mass4l_0->Scale(scaling);
//   Mass4l_1->Scale(scaling);
//   Mass4l_2->Scale(scaling);
//   Mass4l_3->Scale(scaling);
//   Mass4l_4->Scale(scaling);

//   TH1F QED6_4l = (*(Mass4l_0)+*(Mass4l_1)+*(Mass4l_2)+*(Mass4l_3)+*(Mass4l_4));
//   TH1F WZZ_4l  = (*(Mass4l_Wp)+*(Mass4l_Wm));
  
//   c1->cd(2);
//   QED6_4l.SetTitle("4lMass");
//   QED6_4l.Rebin(2);
//   QED6_4l.DrawClone();
//   WZZ_4l.DrawClone("same");



//   //=================== Z0_Pt =========================


//   H6f* Bosons_ZZW_WpZZ = new H6f("ZZW6f",fWp);
//   H6f* ZZW_WmZZ = new H6f("ZZW6f",fWm);

//   Hbos* Bosons_QED6_0 = new Hbos("myAnalyzer/Bosons",f0);
//   Hbos* Bosons_QED6_1 = new Hbos("myAnalyzer/Bosons",f1);
//   Hbos* Bosons_QED6_2 = new Hbos("myAnalyzer/Bosons",f2);
//   Hbos* Bosons_QED6_3 = new Hbos("myAnalyzer/Bosons",f3);
//   Hbos* Bosons_QED6_4 = new Hbos("myAnalyzer/Bosons",f4);





  /////////////////////////////////////////////////////////
  //                                                     //
  //              Signal-Background QED6                 //
  //                                                     //
  /////////////////////////////////////////////////////////

 

}
