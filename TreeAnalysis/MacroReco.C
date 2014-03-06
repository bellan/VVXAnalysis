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

  TFile* ZZ_noSkim       = new TFile("/afs/cern.ch/user/b/bellan/public/for_Giulia/results/zzwAnalysis_baseline/ZZ.root");
  TFile* H_noSkim        = new TFile("/afs/cern.ch/user/b/bellan/public/for_Giulia/results/zzwAnalysis_baseline/H.root");
  TFile* Z_noSkim        = new TFile("/afs/cern.ch/user/b/bellan/public/for_Giulia/results/zzwAnalysis_baseline/Z.root");
  TFile* tt_noSkim       = new TFile("/afs/cern.ch/user/b/bellan/public/for_Giulia/results/zzwAnalysis_baseline/tt.root");
  TFile* triboson_noSkim = new TFile("/afs/cern.ch/user/b/bellan/public/for_Giulia/results/zzwAnalysis_baseline/triboson.root");
  TFile* WZZ_noSkim      = new TFile("/afs/cern.ch/user/b/bellan/public/for_Giulia/results/zzwAnalysis_baseline/WZZJets.root");


  TCanvas *c = new TCanvas("Ricostruiti", "Ricostruiti", 700, 900);
  c->Divide(2,2);
 
  TLegend* l = new TLegend(0.67, 0.61,0.88,0.76);

  TH1F *m6f_ZZ_noSkim       = (TH1F*)ZZ_noSkim->Get("6f_Mass");
  TH1F *m6f_H_noSkim        = (TH1F*)H_noSkim->Get("6f_Mass");
  TH1F *m6f_Z_noSkim        = (TH1F*)Z_noSkim->Get("6f_Mass");
  TH1F *m6f_tt_noSkim       = (TH1F*)tt_noSkim->Get("6f_Mass;1");
  TH1F *m6f_triboson_noSkim = (TH1F*)triboson_noSkim->Get("6f_Mass");
  TH1F *m6f_WZZ_noSkim      = (TH1F*)WZZ_noSkim->Get("6f_Mass");

  TH1F *WPt_ZZ_noSkim       = (TH1F*)ZZ_noSkim->Get("W_Pt");
  TH1F *WPt_H_noSkim        = (TH1F*)H_noSkim->Get("W_Pt");
  TH1F *WPt_Z_noSkim        = (TH1F*)Z_noSkim->Get("W_Pt");
  TH1F *WPt_tt_noSkim       = (TH1F*)tt_noSkim->Get("W_Pt");
  TH1F *WPt_triboson_noSkim = (TH1F*)triboson_noSkim->Get("W_Pt");
  TH1F *WPt_WZZ_noSkim      = (TH1F*)WZZ_noSkim->Get("W_Pt");

  TH1F *Z0Pt_ZZ_noSkim       = (TH1F*)ZZ_noSkim->Get("Z0_Pt");
  TH1F *Z0Pt_H_noSkim        = (TH1F*)H_noSkim->Get("Z0_Pt");
  TH1F *Z0Pt_Z_noSkim        = (TH1F*)Z_noSkim->Get("Z0_Pt");
  TH1F *Z0Pt_tt_noSkim       = (TH1F*)tt_noSkim->Get("Z0_Pt");
  TH1F *Z0Pt_triboson_noSkim = (TH1F*)triboson_noSkim->Get("Z0_Pt");
  TH1F *Z0Pt_WZZ_noSkim      = (TH1F*)WZZ_noSkim->Get("Z0_Pt");

  TH1F *Z1Pt_ZZ_noSkim       = (TH1F*)ZZ_noSkim->Get("Z1_Pt");
  TH1F *Z1Pt_H_noSkim        = (TH1F*)H_noSkim->Get("Z1_Pt");
  TH1F *Z1Pt_Z_noSkim        = (TH1F*)Z_noSkim->Get("Z1_Pt");
  TH1F *Z1Pt_tt_noSkim       = (TH1F*)tt_noSkim->Get("Z1_Pt");
  TH1F *Z1Pt_triboson_noSkim = (TH1F*)triboson_noSkim->Get("Z1_Pt");
  TH1F *Z1Pt_WZZ_noSkim      = (TH1F*)WZZ_noSkim->Get("Z1_Pt");


  ///////////////////// m6f ////////////////////////////

  c->cd(1);
  gPad->SetLogy();

  m6f_ZZ_noSkim->Rebin(8);
  m6f_H_noSkim->Rebin(8);
  m6f_Z_noSkim->Rebin(8);
  m6f_tt_noSkim->Rebin(8);
  m6f_triboson_noSkim->Rebin(8);
  m6f_WZZ_noSkim->Rebin(8);
  
  m6f_ZZ_noSkim->SetLineColor(kGreen);
  m6f_H_noSkim->SetLineColor(kBlue);
  m6f_Z_noSkim->SetLineColor(kYellow);
  m6f_tt_noSkim->SetLineColor(kViolet);
  m6f_triboson_noSkim->SetLineColor(kRed);
  m6f_WZZ_noSkim->SetLineColor(kBlack);
  
  m6f_ZZ_noSkim->GetYaxis()->SetRangeUser(0.00000001,0.01);

  m6f_ZZ_noSkim->Draw();
  m6f_H_noSkim->Draw("same");
  m6f_Z_noSkim->Draw("same");
  m6f_tt_noSkim->Draw("same");
  m6f_triboson_noSkim->Draw("same");
  m6f_WZZ_noSkim->Draw("same");
  
  l->AddEntry(m6f_WZZ_noSkim, "ZZW", "l");
  l->AddEntry(m6f_ZZ_noSkim, "ZZ", "l");
  l->AddEntry(m6f_H_noSkim, "H", "l");
  l->AddEntry(m6f_Z_noSkim, "Z", "l");
  l->AddEntry(m6f_tt_noSkim, "tt", "l");
  l->AddEntry(m6f_triboson_noSkim, "triboson", "l");
  
  
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
  

  ///////////////////// Z0Pt ////////////////////////////
  c->cd(3);
  gPad->SetLogy();

  Z0Pt_ZZ_noSkim->Rebin(4); 
  Z0Pt_H_noSkim->Rebin(4);
  Z0Pt_Z_noSkim->Rebin(4);
  Z0Pt_tt_noSkim->Rebin(4);
  Z0Pt_triboson_noSkim->Rebin(4); 
  Z0Pt_WZZ_noSkim->Rebin(4);

  Z0Pt_ZZ_noSkim->SetLineColor(kGreen);
  Z0Pt_H_noSkim->SetLineColor(kBlue);
  Z0Pt_Z_noSkim->SetLineColor(kYellow);
  Z0Pt_tt_noSkim->SetLineColor(kViolet);
  Z0Pt_triboson_noSkim->SetLineColor(kRed);
  Z0Pt_WZZ_noSkim->SetLineColor(kBlack);
  
  Z0Pt_ZZ_noSkim->GetYaxis()->SetRangeUser(0.00000001,0.01);

  Z0Pt_ZZ_noSkim->Draw();
  Z0Pt_H_noSkim->Draw("same");
  Z0Pt_Z_noSkim->Draw("same");
  Z0Pt_tt_noSkim->Draw("same");
  Z0Pt_triboson_noSkim->Draw("same");
  Z0Pt_WZZ_noSkim->Draw("same");

  ///////////////////// Z1Pt ////////////////////////////
  c->cd(4);
  gPad->SetLogy();

  Z1Pt_ZZ_noSkim->Rebin(4); 
  Z1Pt_H_noSkim->Rebin(4);
  Z1Pt_Z_noSkim->Rebin(4);
  Z1Pt_tt_noSkim->Rebin(4);
  Z1Pt_triboson_noSkim->Rebin(4); 
  Z1Pt_WZZ_noSkim->Rebin(4);

  Z1Pt_ZZ_noSkim->SetLineColor(kGreen);
  Z1Pt_H_noSkim->SetLineColor(kBlue);
  Z1Pt_Z_noSkim->SetLineColor(kYellow);
  Z1Pt_tt_noSkim->SetLineColor(kViolet);
  Z1Pt_triboson_noSkim->SetLineColor(kRed);
  Z1Pt_WZZ_noSkim->SetLineColor(kBlack);
  
  Z1Pt_ZZ_noSkim->GetYaxis()->SetRangeUser(0.00000001,0.01);

  Z1Pt_ZZ_noSkim->Draw();
  Z1Pt_H_noSkim->Draw("same");
  Z1Pt_Z_noSkim->Draw("same");
  Z1Pt_tt_noSkim->Draw("same");
  Z1Pt_triboson_noSkim->Draw("same");
  Z1Pt_WZZ_noSkim->Draw("same");

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
  //                                                                        //
  //                     Skimmed samples                                    // 
  //                                                                        //
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

  TFile* WZZgen = new TFile("/afs/cern.ch/work/g/ggonella/CMSSW_5_3_9/src/VVXAnalysis/TreeAnalysis/results/ZZWAnalyzer_baseline/WZZJets.root");
  TFile* WZZ    = new TFile("/afs/cern.ch/work/g/ggonella/CMSSW_5_3_9/src/VVXAnalysis/TreeAnalysis/results/ZZWAnalyzer_baseline/WZZJets.root");
  
  TH1F *Z0genPt_WZZ = (TH1F*)WZZgen->Get("Z0Gen_Pt");
  TH1F *Z1genPt_WZZ = (TH1F*)WZZgen->Get("Z1Gen_Pt");
  TH1F *WgenPt_WZZ  = (TH1F*)WZZgen->Get("WGen_Pt" );
  
  TH1F *Z0Pt_WZZ = (TH1F*)WZZ->Get("Z0_Pt");
  TH1F *Z1Pt_WZZ = (TH1F*)WZZ->Get("Z1_Pt");
  TH1F *WPt_WZZ  = (TH1F*)WZZ->Get("W_Pt" );
  
  TCanvas * c1 = new TCanvas ("Pt genPart - recoPart", "Pt genPart - recoPart", 1100, 600);
  c1->Divide(3,1);
  
  TLegend* l1 = new TLegend(0.89, 0.89, 0.70, 0.78);
  
  ///////////////////// comparison: Pt genPart - recoPart  ////////////////////////////

  c1->cd(1);

  Z0genPt_WZZ->SetLineColor(kRed);
  Z0Pt_WZZ->SetLineColor(kBlue);
  Z0genPt_WZZ->SetTitle("Z0 Pt");

  Z0genPt_WZZ->Draw();
  Z0Pt_WZZ->Draw("same");
    
  l1->AddEntry(Z0genPt_WZZ, "Z0gen", "l");
  l1->AddEntry(Z0Pt_WZZ, "Z0reco", "l");

  l1->Draw();

  c1->cd(2);
  
  Z1genPt_WZZ->SetLineColor(kRed);
  Z1Pt_WZZ->SetLineColor(kBlue);
  Z1genPt_WZZ->SetTitle("Z1 Pt");

  Z1genPt_WZZ->Draw(); 
  Z1Pt_WZZ->Draw("same");
    
  c1->cd(3);
  
  WgenPt_WZZ->SetLineColor(kRed);
  WPt_WZZ->SetLineColor(kBlue);
  WgenPt_WZZ->SetTitle("W Pt");

  WgenPt_WZZ->Draw(); 
  WPt_WZZ->Draw("same");
  

  gStyle->SetOptStat(0);

}
