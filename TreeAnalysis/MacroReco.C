#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TAxis.h"

#include <iostream>
#include <vector>

using namespace std;

void MacroReco() {

  TFile* ZZ             = new TFile("/afs/cern.ch/work/MSSW_5_3_9/src/VVXAnalysis/TreeAnalysis/results_test/ZZWAnalysis_baseline/ZZ.root");
  TFile* H              = new TFile("/afs/cern.ch/work/MSSW_5_3_9/src/VVXAnalysis/TreeAnalysis/results_test/ZZWAnalysis_baseline/H.root");
  TFile* Z              = new TFile("/afs/cern.ch/work/MSSW_5_3_9/src/VVXAnalysis/TreeAnalysis/results_test/ZZWAnalysis_baseline/Z.root");
  TFile* tt             = new TFile("/afs/cern.ch/work/MSSW_5_3_9/src/VVXAnalysis/TreeAnalysis/results_test/ZZWAnalysis_baseline/tt.root");
  TFile* triboson       = new TFile("/afs/cern.ch/work/MSSW_5_3_9/src/VVXAnalysis/TreeAnalysis/results_test/ZZWAnalysis_baseline/triboson.root");
  TFile* WZ             = new TFile("/afs/cern.ch/work/MSSW_5_3_9/src/VVXAnalysis/TreeAnalysis/results_test/ZZWAnalysis_baseline/WZ.root");
  TFile* WZZSignal      = new TFile("/afs/cern.ch/work/MSSW_5_3_9/src/VVXAnalysis/TreeAnalysis/results_test/ZZWAnalysis_baseline/WZZJetsSignal.root");
  TFile* WZZNoSignal    = new TFile("/afs/cern.ch/work/MSSW_5_3_9/src/VVXAnalysis/TreeAnalysis/results_test/ZZWAnalysis_baseline/WZZJetsNoSignal.root");
  TFile* ZZJetsSignal   = new TFile("/afs/cern.ch/work/MSSW_5_3_9/src/VVXAnalysis/TreeAnalysis/results_test/ZZWAnalysis_baseline/ZZJetsTo4L_Signal.root");
  TFile* ZZJetsNoSignal = new TFile("/afs/cern.ch/work/MSSW_5_3_9/src/VVXAnalysis/TreeAnalysis/results_test/ZZWAnalysis_baseline/ZZJetsTo4L_NoSignal.root");

  TCanvas *c = new TCanvas("m6f - BosonsPt", "m6f - BosonsPt", 700, 900);
  c->Divide(2,2); 
  TLegend *l = new TLegend(0.67, 0.61,0.88,0.76);

  TCanvas *c1 = new TCanvas("Delta Phi", "Delta Phi", 700, 900);
  c1->Divide(2,2); 
  TLegend *l1 = new TLegend(0.67, 0.61,0.88,0.76);
  
  TH1F *m6f_ZZ             = (TH1F*)ZZ->Get("6f_Mass");
  TH1F *m6f_H              = (TH1F*)H->Get("6f_Mass");
  TH1F *m6f_Z              = (TH1F*)Z->Get("6f_Mass");
  TH1F *m6f_tt             = (TH1F*)tt->Get("6f_Mass;1");
  TH1F *m6f_triboson       = (TH1F*)triboson->Get("6f_Mass");
  TH1F *m6f_WZ             = (TH1F*)WZ->Get("6f_Mass");
  TH1F *m6f_WZZSignal      = (TH1F*)WZZSignal->Get("6f_Mass");
  TH1F *m6f_WZZNoSignal    = (TH1F*)WZZNoSignal->Get("6f_Mass");
  TH1F *m6f_ZZJetsSignal   = (TH1F*)ZZJetsSignal->Get("6f_Mass");
  TH1F *m6f_ZZJetsNoSignal = (TH1F*)ZZJetsNoSignal->Get("6f_Mass");
  
  TH1F *WPt_ZZ             = (TH1F*)ZZ->Get("W_Pt");
  TH1F *WPt_H              = (TH1F*)H->Get("W_Pt");
  TH1F *WPt_Z              = (TH1F*)Z->Get("W_Pt");
  TH1F *WPt_tt             = (TH1F*)tt->Get("W_Pt");
  TH1F *WPt_triboson       = (TH1F*)triboson->Get("W_Pt");
  TH1F *WPt_WZ             = (TH1F*)WZ->Get("W_Pt");
  TH1F *WPt_WZZSignal      = (TH1F*)WZZSignal->Get("W_Pt");
  TH1F *WPt_WZZNoSignal    = (TH1F*)WZZNoSignal->Get("W_Pt");
  TH1F *WPt_ZZJetsSignal   = (TH1F*)ZZJetsSignal->Get("W_Pt");
  TH1F *WPt_ZZJetsNoSignal = (TH1F*)ZZJetsNoSignal->Get("W_Pt");

  TH1F *Z0Pt_ZZ             = (TH1F*)ZZ->Get("Z0_Pt");
  TH1F *Z0Pt_H              = (TH1F*)H->Get("Z0_Pt");
  TH1F *Z0Pt_Z              = (TH1F*)Z->Get("Z0_Pt");
  TH1F *Z0Pt_tt             = (TH1F*)tt->Get("Z0_Pt");
  TH1F *Z0Pt_triboson       = (TH1F*)triboson->Get("Z0_Pt");
  TH1F *Z0Pt_WZ             = (TH1F*)WZ->Get("Z0_Pt");
  TH1F *Z0Pt_WZZSignal      = (TH1F*)WZZSignal->Get("Z0_Pt");
  TH1F *Z0Pt_WZZNoSignal    = (TH1F*)WZZNoSignal->Get("Z0_Pt");
  TH1F *Z0Pt_ZZJetsSignal   = (TH1F*)ZZJetsSignal->Get("Z0_Pt");
  TH1F *Z0Pt_ZZJetsNoSignal = (TH1F*)ZZJetsNoSignal->Get("Z0_Pt");
 
  TH1F *Z1Pt_ZZ             = (TH1F*)ZZ->Get("Z1_Pt");
  TH1F *Z1Pt_H              = (TH1F*)H->Get("Z1_Pt");
  TH1F *Z1Pt_Z              = (TH1F*)Z->Get("Z1_Pt");
  TH1F *Z1Pt_tt             = (TH1F*)tt->Get("Z1_Pt");
  TH1F *Z1Pt_triboson       = (TH1F*)triboson->Get("Z1_Pt");
  TH1F *Z1Pt_WZ             = (TH1F*)WZ->Get("Z1_Pt");
  TH1F *Z1Pt_WZZSignal      = (TH1F*)WZZSignal->Get("Z1_Pt");
  TH1F *Z1Pt_WZZNoSignal    = (TH1F*)WZZNoSignal->Get("Z1_Pt");
  TH1F *Z1Pt_ZZJetsSignal   = (TH1F*)ZZJetsSignal->Get("Z1_Pt");
  TH1F *Z1Pt_ZZJetsNoSignal = (TH1F*)ZZJetsNoSignal->Get("Z1_Pt");
 
  ///////////////////// m6f ////////////////////////////

  c->cd(1);
  gPad->SetLogy();

  m6f_ZZ->Rebin(8);
  m6f_H->Rebin(8);
  m6f_Z->Rebin(8);
  m6f_tt->Rebin(8);
  m6f_triboson->Rebin(8);
  m6f_WZ->Rebin(8); 
  m6f_WZZSignal->Rebin(8);      
  m6f_WZZNoSignal->Rebin(8);    
  m6f_ZZJetsSignal->Rebin(8);   
  m6f_ZZJetsNoSignal->Rebin(8);
  
  m6f_ZZ->SetLineColor(kGreen);
  m6f_H->SetLineColor(kBlue);
  m6f_Z->SetLineColor(kYellow);
  m6f_tt->SetLineColor(kViolet);
  m6f_triboson_noSkim->SetLineColor(kRed);
  m6f_WZ->SetLineColor(kRed); 
  m6f_WZZSignal->SetLineColor(kRed);
  m6f_WZZNoSignal-> SetLineColor(kRed);   
  m6f_ZZJetsSignal->SetLineColor(kRed);
  m6f_ZZJetsNoSignal->SetLineColor(kRed);
 
  
  m6f_ZZ->GetYaxis()->SetRangeUser(0.00000001,0.01);

  m6f_ZZ->Draw();
  m6f_H->Draw("same");
  m6f_Z->Draw("same");
  m6f_tt->Draw("same");
  m6f_triboson->Draw("same");
  m6f_WZ->Draw("same");
  m6f_WZZSignal->Draw("same");   
  m6f_WZZNoSignal->Draw("same");    
  m6f_ZZJetsSignal->Draw("same");
  m6f_ZZJetsNoSignal->Draw("same");

  
  l->AddEntry(m6f_WZZSignal, "WZZSignal", "l");
  l->AddEntry(m6f_ZZJetsSignal, "ZZJetsSignal", "l");
  l->AddEntry(m6f_ZZ, "ZZ", "l");
  l->AddEntry(m6f_H, "H", "l");
  l->AddEntry(m6f_Z, "Z", "l");
  l->AddEntry(m6f_tt, "tt", "l");
  l->AddEntry(m6f_triboson, "triboson", "l");
  l->AddEntry(m6f_WZZNoSignal, "WZZNoSignal", "l");
  l->AddEntry(m6f_WZ, "WZ", "l");
  l->AddEntry(m6f_ZZJetsNoSignal, "ZZJetsNoSignal", "l");
 

  
  l->Draw();

  ///////////////////// WPt ////////////////////////////  

  c->cd(2);
  gPad->SetLogy();

  WPt_ZZ_noSkim->Rebin(4); 
  WPt_H_noSkim->Rebin(4);
  WPt_Z_noSkim->Rebin(4);
  WPt_tt_noSkim->Rebin(4);
  WPt_triboson_noSkim->Rebin(4); 
  WPt_WZZ_noSkim->Rebin(4);
  
  WPt_ZZ_noSkim->SetLineColor(kGreen);
  WPt_H_noSkim->SetLineColor(kBlue);
  WPt_Z_noSkim->SetLineColor(kYellow);
  WPt_tt_noSkim->SetLineColor(kViolet);
  WPt_triboson_noSkim->SetLineColor(kRed);
  WPt_WZZ_noSkim->SetLineColor(kBlack);

  WPt_ZZ_noSkim->GetYaxis()->SetRangeUser(0.00000001,0.01);

  WPt_ZZ_noSkim->Draw();
  WPt_H_noSkim->Draw("same");
  WPt_Z_noSkim->Draw("same");
  WPt_tt_noSkim->Draw("same");
  WPt_triboson_noSkim->Draw("same");
  WPt_WZZ_noSkim->Draw("same");
  
  




  gStyle->SetOptStat(0);

}
