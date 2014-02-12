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

TH1F* Integral(TH1F* h);
TH1F* Divide(TH1F* h_num, TH1F* h_den);

void macro_ZZW_gen() {
 
  gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginVVVAnalysisAnalysisStepPlugins.so");

  gROOT->LoadMacro("../interface/H6f.h+");
  gROOT->LoadMacro("../interface/Hbos.h+");
  gROOT->LoadMacro("../interface/Hjets.h+");
  gROOT->LoadMacro("../interface/Boson.h+");

    TFile* f0_0 = new TFile("QED6_0_3_cat0.root");  
    TFile* f0_1 = new TFile("QED6_1_3_cat0.root");
    TFile* f0_2 = new TFile("QED6_2_3_cat0.root");
    TFile* f0_3 = new TFile("QED6_3_3_cat0.root");
    TFile* f0_4 = new TFile("QED6_4_3_cat0.root"); 

    TFile* f1_0 = new TFile("QED6_0_3_cat1.root");  
    TFile* f1_1 = new TFile("QED6_1_3_cat1.root");
    TFile* f1_2 = new TFile("QED6_2_3_cat1.root");
    TFile* f1_3 = new TFile("QED6_3_3_cat1.root");
    TFile* f1_4 = new TFile("QED6_4_3_cat1.root"); 
    
    TFile* f2_0 = new TFile("QED6_0_3_cat2.root");  
    TFile* f2_1 = new TFile("QED6_1_3_cat2.root");
    TFile* f2_2 = new TFile("QED6_2_3_cat2.root");
    TFile* f2_3 = new TFile("QED6_3_3_cat2.root");
    TFile* f2_4 = new TFile("QED6_4_3_cat2.root"); 

  gDirectory->cd("Rint:/");
   
  Hbos* Bosons_QED6_0_cat0 = new Hbos("Bosons",f0_0);
  Hbos* Bosons_QED6_1_cat0 = new Hbos("Bosons",f0_1);
  Hbos* Bosons_QED6_2_cat0 = new Hbos("Bosons",f0_2);
  Hbos* Bosons_QED6_3_cat0 = new Hbos("Bosons",f0_3);
  Hbos* Bosons_QED6_4_cat0 = new Hbos("Bosons",f0_4);
 
  Hbos* Bosons_QED6_0_cat1 = new Hbos("Bosons",f1_0);
  Hbos* Bosons_QED6_1_cat1 = new Hbos("Bosons",f1_1);
  Hbos* Bosons_QED6_2_cat1 = new Hbos("Bosons",f1_2);
  Hbos* Bosons_QED6_3_cat1 = new Hbos("Bosons",f1_3);
  Hbos* Bosons_QED6_4_cat1 = new Hbos("Bosons",f1_4);

  Hbos* Bosons_QED6_0_cat2 = new Hbos("Bosons",f2_0);
  Hbos* Bosons_QED6_1_cat2 = new Hbos("Bosons",f2_1);
  Hbos* Bosons_QED6_2_cat2 = new Hbos("Bosons",f2_2);
  Hbos* Bosons_QED6_3_cat2 = new Hbos("Bosons",f2_3);
  Hbos* Bosons_QED6_4_cat2 = new Hbos("Bosons",f2_4);

  Hbos* Bosons_QED6_0_cat0_Cut = new Hbos("BosonsCut",f0_0);
  Hbos* Bosons_QED6_1_cat0_Cut = new Hbos("BosonsCut",f0_1);
  Hbos* Bosons_QED6_2_cat0_Cut = new Hbos("BosonsCut",f0_2);
  Hbos* Bosons_QED6_3_cat0_Cut = new Hbos("BosonsCut",f0_3);
  Hbos* Bosons_QED6_4_cat0_Cut = new Hbos("BosonsCut",f0_4);

  Hbos* Bosons_QED6_0_cat1_Cut = new Hbos("BosonsCut",f1_0);
  Hbos* Bosons_QED6_1_cat1_Cut = new Hbos("BosonsCut",f1_1);
  Hbos* Bosons_QED6_2_cat1_Cut = new Hbos("BosonsCut",f1_2);
  Hbos* Bosons_QED6_3_cat1_Cut = new Hbos("BosonsCut",f1_3);
  Hbos* Bosons_QED6_4_cat1_Cut = new Hbos("BosonsCut",f1_4);

  Hbos* Bosons_QED6_0_cat2_Cut = new Hbos("BosonsCut",f2_0);
  Hbos* Bosons_QED6_1_cat2_Cut = new Hbos("BosonsCut",f2_1);
  Hbos* Bosons_QED6_2_cat2_Cut = new Hbos("BosonsCut",f2_2);
  Hbos* Bosons_QED6_3_cat2_Cut = new Hbos("BosonsCut",f2_3);
  Hbos* Bosons_QED6_4_cat2_Cut = new Hbos("BosonsCut",f2_4);
  
  TH1F* lostEv_0 = new TH1F("lostEv_0","lostEv_0",3, 0., 3.);
  lostEv_0 = (TH1F*)f0_0->Get("MyAnalyzer/lostEv");
  TH1F* lostEv_1 = new TH1F("lostEv_1","lostEv_1",3, 0., 3.);
  lostEv_1 = (TH1F*)f0_1->Get("MyAnalyzer/lostEv");
  TH1F* lostEv_2 = new TH1F("lostEv_2","lostEv_2",3, 0., 3.);
  lostEv_2 = (TH1F*)f0_2->Get("MyAnalyzer/lostEv");
  TH1F* lostEv_3 = new TH1F("lostEv_3","lostEv_3",3, 0., 3.);
  lostEv_3 = (TH1F*)f0_3->Get("MyAnalyzer/lostEv");
  TH1F* lostEv_4 = new TH1F("lostEv_4","lostEv_4",3, 0., 3.);
  lostEv_4 = (TH1F*)f0_4->Get("MyAnalyzer/lostEv");

  float L = 300.;

  float xsec_QED6 = 52.76;

  int Ngen_QED6 = 9998*5;

  float w_QED6 = (L*xsec_QED6)/Ngen_QED6;
 
  cout << "\n%%%%%%%%%%%%\nWEIGHTS\nw_QED6= " << w_QED6  <<endl; 
 

  Bosons_QED6_0_cat0->Scale(w_QED6);
  Bosons_QED6_1_cat0->Scale(w_QED6);
  Bosons_QED6_2_cat0->Scale(w_QED6);
  Bosons_QED6_3_cat0->Scale(w_QED6);
  Bosons_QED6_4_cat0->Scale(w_QED6);

  Bosons_QED6_0_cat1->Scale(w_QED6);
  Bosons_QED6_1_cat1->Scale(w_QED6);
  Bosons_QED6_2_cat1->Scale(w_QED6);
  Bosons_QED6_3_cat1->Scale(w_QED6);
  Bosons_QED6_4_cat1->Scale(w_QED6);

  Bosons_QED6_0_cat2->Scale(w_QED6);
  Bosons_QED6_1_cat2->Scale(w_QED6);
  Bosons_QED6_2_cat2->Scale(w_QED6);
  Bosons_QED6_3_cat2->Scale(w_QED6);
  Bosons_QED6_4_cat2->Scale(w_QED6);
 
  
  TCanvas* myc0 = new TCanvas("Z0 Mass", "Z0 Mass", 900,700);
  myc0->Divide(2,1);

  TLegend* leg0 = new TLegend(0.67, 0.61,0.88,0.76);
 
  
  //Canvas 0----------------SIGNAL QED6 - ZZW (produced on shell) - ZZZ---------------------//
  myc->cd(1);
  TH1F Z0Mass = *(Bosons_QED6_0_cat0->hZ0Mass) + *(Bosons_QED6_1_cat0->hZ0Mass) + *(Bosons_QED6_2_cat0->hZ0Mass) + *(Bosons_QED6_3_cat0->hZ0Mass) + *(Bosons_QED6_4_cat0->hZ0Mass);  
  Z0Mass.SetTitle("Z0Mass");
  Z0Mass.DrawClone();
  TH1F Z0MassCut = *(Bosons_QED6_0_cat0_Cut->hZ0Mass) + *(Bosons_QED6_1_cat0_Cut->hZ0Mass) + *(Bosons_QED6_2_cat0_Cut->hZ0Mass) + *(Bosons_QED6_3_cat0_Cut->hZ0Mass) + *(Bosons_QED6_4_cat0_Cut->hZ0Mass);  
  Z0MassCut.DrawClone("same");
  
  leg0->AddEntry(Bosons_QED6_0_cat0->hZ0Mass, "NoCut", "l");
  leg0->AddEntry(Bosons_QED6_0_cat0_Cut->hZ0Mass, "Cut", "l");
  leg0->Draw();
 

 TStyle::gStyle->SetOptStat(0);
  
 }


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// INTEGRAL //
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TH1F* Integral(TH1F* h) {
  
  TH1F* h1=h->Clone("h1");	      
  h1->Reset();
  
  int lastbin = h->GetNbinsX();
  
  for (int bin=1;bin<=lastbin; bin++) {
    float int_value = h->Integral(bin,lastbin);
    h1->Fill(bin,int_value);  
  }
  
  return h1;
  
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// DIVIDE //
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TH1F* Divide(TH1F* h_num, TH1F* h_den) {
  
  TH1F* h_div= h_den->Clone();
  h_div->Reset();
  
  int lastbin_den = h_den->GetNbinsX();
  
  for (int bin_i=1;bin_i<=lastbin_den;bin_i++) {
    float num = h_num->GetBinContent(bin_i);
    float den = h_den->GetBinContent(bin_i);
    if (den!=0) {
      float div = num/sqrt(den);
      h_div->Fill(bin_i,div);
    }    
  }
  
  return h_div;
  
}
