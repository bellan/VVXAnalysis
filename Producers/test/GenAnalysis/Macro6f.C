#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TObject.h"

#include <iostream>
#include <vector>

using namespace std;

void Macro6f() {

  gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginZZWGenAnalysisPlugins.so");
  
  gSystem->AddIncludePath(" -I$CMSSW_BASE/src/ -I$CMSSW_RELEASE_BASE/src");
  
  
  std::string path = std::string(gROOT->GetMacroPath())+"$CMSSW_BASE/src:$CMSSW_BASE/src/VVXAnalysis/Producers/test/GenAnalysis";
  gROOT->SetMacroPath(path.c_str());
  
 
  TFile* QED6_S = new TFile("QED6_SIGNAL.root");
  TFile* QED6_B = new TFile("QED6_BACKGR.root");
  
  TCanvas* c =  new TCanvas("Mass 6 fermions","Mass 6 fermions", 900, 700);
  c->Divide(2,1);
  
  TLegend* leg = new TLegend(0.67, 0.61,0.88,0.76);
  

  float L = 300.;
  
  float xsec_QED6 = 52.76;
  int Ngen_QED6   = 9998*5;
  float w_QED6    = (L*xsec_QED6)/Ngen_QED6;
  
  TH1F* m6f_sign = new TH1F("m6f_sign", "m6f_sign", 3000, 0, 3000);
  m6f_sign = (TH1F*)QED6_S->Get("myAnalyzer/all6fMass_0");
  TH1F* m6f_1 = new TH1F("m6f_1", "m6f_1", 3000, 0, 3000);
  m6f_1 = (TH1F*)QED6_B->Get("myAnalyzer/all6fMass_1");
  TH1F* m6f_2 = new TH1F("m6f_2", "m6f_2", 3000, 0, 3000);
  m6f_2 = (TH1F*)QED6_B->Get("myAnalyzer/all6fMass_2");
  TH1F* m6f_3 = new TH1F("m6f_3", "m6f_3", 3000, 0, 3000);
  m6f_3 = (TH1F*)QED6_B->Get("myAnalyzer/all6fMass_3");
  TH1F* m6f_4 = new TH1F("m6f_4", "m6f_4", 3000, 0, 3000);
  m6f_4 = (TH1F*)QED6_B->Get("myAnalyzer/all6fMass_4");
  TH1F* m6f_5 = new TH1F("m6f_5", "m6f_5", 3000, 0, 3000);
  m6f_5 = (TH1F*)QED6_B->Get("myAnalyzer/all6fMass_5");
  TH1F* m6f_6 = new TH1F("m6f_6", "m6f_6", 3000, 0, 3000);
  m6f_6 = (TH1F*)QED6_B->Get("myAnalyzer/all6fMass_6");
  TH1F* m6f_7 = new TH1F("m6f_7", "m6f_7", 3000, 0, 3000);
  m6f_7 = (TH1F*)QED6_B->Get("myAnalyzer/all6fMass_7");
  TH1F* m6f_8 = new TH1F("m6f_8", "m6f_8", 3000, 0, 3000);
  m6f_8 = (TH1F*)QED6_B->Get("myAnalyzer/all6fMass_8");
  TH1F* m6f_9 = new TH1F("m6f_9", "m6f_9", 3000, 0, 3000);
  m6f_9 = (TH1F*)QED6_B->Get("myAnalyzer/all6fMass_9");

  TH1F* my1 = new TH1F("my1", "my1", 3000, 0, 3000);
  TH1F* my2 = new TH1F("my2", "my2", 3000, 0, 3000);
  TH1F* my3 = new TH1F("my3", "my3", 3000, 0, 3000);
  TH1F* my4 = new TH1F("my4", "my4", 3000, 0, 3000);
  TH1F* my5 = new TH1F("my5", "my5", 3000, 0, 3000);
  TH1F* my6 = new TH1F("my6", "my6", 3000, 0, 3000);
  TH1F* my7 = new TH1F("my7", "my7", 3000, 0, 3000);
  TH1F* my8 = new TH1F("my8", "my8", 3000, 0, 3000);
  TH1F* my9 = new TH1F("my9", "my9", 3000, 0, 3000);  
  
  m6f_sign->Scale(w_QED6);
  m6f_1->Scale(w_QED6);
  m6f_2->Scale(w_QED6);
  m6f_3->Scale(w_QED6);
  m6f_4->Scale(w_QED6);
  m6f_5->Scale(w_QED6);
  m6f_6->Scale(w_QED6);
  m6f_7->Scale(w_QED6);
  m6f_8->Scale(w_QED6);
  m6f_9->Scale(w_QED6);
  

  my1 = (TH1F*)m6f_3->Clone();
  my2 = (TH1F*)my1->Clone();
  my2->Add(m6f_2,1);
  my3 = (TH1F*)my2->Clone();
  my3->Add(m6f_1,1);
  my4 = (TH1F*)my3->Clone();
  my4->Add(m6f_8,1);
  my5 = (TH1F*)my4->Clone();
  my5->Add(m6f_5,1);
  my6 = (TH1F*)my5->Clone();
  my6->Add(m6f_6,1);
  my7 = (TH1F*)my6->Clone();
  my7->Add(m6f_7,1);
  my8 = (TH1F*)my7->Clone();
  my8->Add(m6f_4,1);
  my9 = (TH1F*)my8->Clone();
  my9->Add(m6f_9,1);
  
  c->cd(0);
  c->SetLogy();
  my1->SetFillColor(kGreen);
  my2->SetFillColor(kRed);  
  my3->SetFillColor(kBlue);
  my4->SetFillColor(kPink);
  my5->SetFillColor(kBlack)
  my6->SetFillColor(kYellow);
  my7->SetFillColor(kAzure);
  my8->SetFillColor(kViolet);
  my9->SetFillColor(kOrange);
  m6f_sign->SetLineColor(kBlack);

  my9->SetTitle("Mass 6 fermions");
  my9->SetMaximum(1000.);
  my9->Rebin(2);
  my9->Draw();
  my8->Draw("same");
  my7->Draw("same");
  my6->Draw("same");
  my5->Draw("same");
  my4->Draw("same");
  my3->Draw("same");
  my2->Draw("same");
  my1->Draw("same");
  m6f_sign->Draw("same");

 
  leg->AddEntry(my1, "1 backgr", "f");
  leg->AddEntry(my2, "2 backgr", "f");
  leg->AddEntry(my3, "3 backgr", "f");
  leg->AddEntry(my4, "4 backgr", "f");
  leg->AddEntry(my5, "5 backgr", "f");
  leg->AddEntry(my6, "6 backgr", "f");
  leg->AddEntry(my7, "7 backgr", "f");
  leg->AddEntry(my8, "8 backgr", "f");
  leg->AddEntry(my9, "9 backgr", "f");
  leg->AddEntry(m6f_sign, "signal", "l");
  leg->Draw();

  TStyle::gStyle->SetOptStat(0);
}
