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

  TFile* ZZ = new TFile("/afs/cern.ch/user/b/bellan/public/for_Giulia/results/zzwAnalysis_baseline/ZZ.root");
  TFile* H = new TFile("/afs/cern.ch/user/b/bellan/public/for_Giulia/results/zzwAnalysis_baseline/H.root");
  TFile* Z = new TFile("/afs/cern.ch/user/b/bellan/public/for_Giulia/results/zzwAnalysis_baseline/Z.root");
  TFile* tt = new TFile("/afs/cern.ch/user/b/bellan/public/for_Giulia/results/zzwAnalysis_baseline/tt.root");
  TFile* triboson = new TFile("/afs/cern.ch/user/b/bellan/public/for_Giulia/results/zzwAnalysis_baseline/triboson.root");
  TFile* WZZ = new TFile("/afs/cern.ch/user/b/bellan/public/for_Giulia/results/zzwAnalysis_baseline/WZZJets.root");

  TCanvas *c = new TCanvas("Ricostruiti", "Ricostruiti", 700,900);
  c->Divide(2,2);

  TLegend* l = new TLegend(0.67, 0.61,0.88,0.76);

  TH1F *m6f_ZZ = (TH1F*)ZZ->Get("6f_Mass");
  TH1F *m6f_H = (TH1F*)H->Get("6f_Mass");
  TH1F *m6f_Z = (TH1F*)Z->Get("6f_Mass");
  TH1F *m6f_tt = (TH1F*)tt->Get("6f_Mass;1");
  TH1F *m6f_triboson = (TH1F*)triboson->Get("6f_Mass");
  TH1F *m6f_WZZ = (TH1F*)WZZ->Get("6f_Mass");

  TH1F *WPt_ZZ = (TH1F*)ZZ->Get("W_Pt");
  TH1F *WPt_H = (TH1F*)H->Get("W_Pt");
  TH1F *WPt_Z  = (TH1F*)Z->Get("W_Pt");
  TH1F *WPt_tt = (TH1F*)tt->Get("W_Pt");
  TH1F *WPt_triboson = (TH1F*)triboson->Get("W_Pt");
  TH1F *WPt_WZZ = (TH1F*)WZZ->Get("W_Pt");

  TH1F *Z0Pt_ZZ = (TH1F*)ZZ->Get("Z0_Pt");
  TH1F *Z0Pt_H = (TH1F*)H->Get("Z0_Pt");
  TH1F *Z0Pt_Z  = (TH1F*)Z->Get("Z0_Pt");
  TH1F *Z0Pt_tt = (TH1F*)tt->Get("Z0_Pt");
  TH1F *Z0Pt_triboson = (TH1F*)triboson->Get("Z0_Pt");
  TH1F *Z0Pt_WZZ = (TH1F*)WZZ->Get("Z0_Pt");

  TH1F *Z1Pt_ZZ = (TH1F*)ZZ->Get("Z1_Pt");
  TH1F *Z1Pt_H = (TH1F*)H->Get("Z1_Pt");
  TH1F *Z1Pt_Z  = (TH1F*)Z->Get("Z1_Pt");
  TH1F *Z1Pt_tt = (TH1F*)tt->Get("Z1_Pt");
  TH1F *Z1Pt_triboson = (TH1F*)triboson->Get("Z1_Pt");
  TH1F *Z1Pt_WZZ = (TH1F*)WZZ->Get("Z1_Pt");

  ///////////////////// m6f ////////////////////////////
  c->cd(1);
  gPad->SetLogy();

  m6f_ZZ->Rebin(8);
  m6f_H ->Rebin(8);
  m6f_Z->Rebin(8);
  m6f_tt->Rebin(8);
  m6f_triboson->Rebin(8);
  m6f_WZZ->Rebin(8);
  
  m6f_ZZ->SetLineColor(kGreen);
  m6f_H->SetLineColor(kBlue);
  m6f_Z->SetLineColor(kYellow);
  m6f_tt->SetLineColor(kViolet);
  m6f_triboson->SetLineColor(kRed);
  m6f_WZZ->SetLineColor(kBlack);
  
  m6f_ZZ->GetYaxis()->SetRangeUser(0.00000001,0.01);

  m6f_ZZ->Draw();
  m6f_H->Draw("same");
  m6f_Z->Draw("same");
  m6f_tt->Draw("same");
  m6f_triboson->Draw("same");
  m6f_WZZ->Draw("same");
  
  l->AddEntry(m6f_WZZ, "ZZW", "l");
  l->AddEntry(m6f_ZZ, "ZZ", "l");
  l->AddEntry(m6f_H, "H", "l");
  l->AddEntry(m6f_Z, "Z", "l");
  l->AddEntry(m6f_tt, "tt", "l");
  l->AddEntry(m6f_triboson, "triboson", "l");
  
  
  l->Draw();

  ///////////////////// WPt ////////////////////////////  
  c->cd(2);
  gPad->SetLogy();

  WPt_ZZ->Rebin(4); 
  WPt_H->Rebin(4);
  WPt_Z->Rebin(4);
  WPt_tt->Rebin(4);
  WPt_triboson->Rebin(4); 
  WPt_WZZ->Rebin(4);
  
  WPt_ZZ->SetLineColor(kGreen);
  WPt_H->SetLineColor(kBlue);
  WPt_Z->SetLineColor(kYellow);
  WPt_tt->SetLineColor(kViolet);
  WPt_triboson->SetLineColor(kRed);
  WPt_WZZ->SetLineColor(kBlack);

  WPt_ZZ->GetYaxis()->SetRangeUser(0.00000001,0.01);

  WPt_ZZ->Draw();
  WPt_H->Draw("same");
  WPt_Z->Draw("same");
  WPt_tt->Draw("same");
  WPt_triboson->Draw("same");
  WPt_WZZ->Draw("same");
  

  ///////////////////// Z0Pt ////////////////////////////
  c->cd(3);
  gPad->SetLogy();

  Z0Pt_ZZ->Rebin(4); 
  Z0Pt_H->Rebin(4);
  Z0Pt_Z->Rebin(4);
  Z0Pt_tt->Rebin(4);
  Z0Pt_triboson->Rebin(4); 
  Z0Pt_WZZ->Rebin(4);

  Z0Pt_ZZ->SetLineColor(kGreen);
  Z0Pt_H->SetLineColor(kBlue);
  Z0Pt_Z->SetLineColor(kYellow);
  Z0Pt_tt->SetLineColor(kViolet);
  Z0Pt_triboson->SetLineColor(kRed);
  Z0Pt_WZZ->SetLineColor(kBlack);
  
  Z0Pt_ZZ->GetYaxis()->SetRangeUser(0.00000001,0.01);

  Z0Pt_ZZ->Draw();
  Z0Pt_H->Draw("same");
  Z0Pt_Z->Draw("same");
  Z0Pt_tt->Draw("same");
  Z0Pt_triboson->Draw("same");
  Z0Pt_WZZ->Draw("same");

  ///////////////////// Z1Pt ////////////////////////////
  c->cd(4);
  gPad->SetLogy();

  Z1Pt_ZZ->Rebin(4); 
  Z1Pt_H->Rebin(4);
  Z1Pt_Z->Rebin(4);
  Z1Pt_tt->Rebin(4);
  Z1Pt_triboson->Rebin(4); 
  Z1Pt_WZZ->Rebin(4);

  Z1Pt_ZZ->SetLineColor(kGreen);
  Z1Pt_H->SetLineColor(kBlue);
  Z1Pt_Z->SetLineColor(kYellow);
  Z1Pt_tt->SetLineColor(kViolet);
  Z1Pt_triboson->SetLineColor(kRed);
  Z1Pt_WZZ->SetLineColor(kBlack);
  
  Z1Pt_ZZ->GetYaxis()->SetRangeUser(0.00000001,0.01);

  Z1Pt_ZZ->Draw();
  Z1Pt_H->Draw("same");
  Z1Pt_Z->Draw("same");
  Z1Pt_tt->Draw("same");
  Z1Pt_triboson->Draw("same");
  Z1Pt_WZZ->Draw("same");

  TStyle::gStyle->SetOptStat(0);

}
