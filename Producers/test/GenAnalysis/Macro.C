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

void Macro() {

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

  int numb;

  cout << "Type a number to analyze data from \n1: MC history \n2: Real signal, MadGraph pairing \n3: Real signal, real pairing" << endl;
  cin >> numb;

  if (numb == 1 ) {
    
    TFile* f0 = new TFile("QED6_0_1_GenAn.root");  
    TFile* f1 = new TFile("QED6_1_1_GenAn.root");
    TFile* f2 = new TFile("QED6_2_1_GenAn.root");
    TFile* f3 = new TFile("QED6_3_1_GenAn.root");
    TFile* f4 = new TFile("QED6_4_1_GenAn.root");  
    TFile* fWp = new TFile("W+ZZ_1_GenAn.root");
    TFile* fWm = new TFile("W-ZZ_1_GenAn.root");

  }
  
  if ( numb == 2 ) {
    
    TFile* f0 = new TFile("QED6_0_2_GenAn.root"); 
    TFile* f1 = new TFile("QED6_1_2_GenAn.root");
    TFile* f2 = new TFile("QED6_2_2_GenAn.root");
    TFile* f3 = new TFile("QED6_3_2_GenAn.root");
    TFile* f4 = new TFile("QED6_4_2_GenAn.root");  
    TFile* fWp = new TFile("W+ZZ_2_GenAn.root");
    TFile* fWm = new TFile("W-ZZ_2_GenAn.root");
    
  }
  
  if ( numb == 3 ) {
    
    TFile* f0 = new TFile("QED6_0_3_GenAn.root");
    TFile* f1 = new TFile("QED6_1_3_GenAn.root");
    TFile* f2 = new TFile("QED6_2_3_GenAn.root");
    TFile* f3 = new TFile("QED6_3_3_GenAn.root");
    TFile* f4 = new TFile("QED6_4_3_GenAn.root");  
    TFile* fWp = new TFile("W+ZZ_3_GenAn.root");
    TFile* fWm = new TFile("W-ZZ_3_GenAn.root");
    
  }
  
  gDirectory->cd("Rint:/");
  
  H6f* ZZW_QED6_0 = new H6f("ZZW6f",f0);
  H6f* ZZW_QED6_1 = new H6f("ZZW6f",f1);
  H6f* ZZW_QED6_2 = new H6f("ZZW6f",f2);
  H6f* ZZW_QED6_3 = new H6f("ZZW6f",f3);
  H6f* ZZW_QED6_4 = new H6f("ZZW6f",f4);
  H6f* Backgr_QED6_0 = new H6f("Backgr",f0);
  H6f* Backgr_QED6_1 = new H6f("Backgr",f1);
  H6f* Backgr_QED6_2 = new H6f("Backgr",f2);
  H6f* Backgr_QED6_3 = new H6f("Backgr",f3);
  H6f* Backgr_QED6_4 = new H6f("Backgr",f4);
  H6f* ZZZ_QED6_0 = new H6f("ZZZ6f",f0);
  H6f* ZZZ_QED6_1 = new H6f("ZZZ6f",f1);
  H6f* ZZZ_QED6_2 = new H6f("ZZZ6f",f2);
  H6f* ZZZ_QED6_3 = new H6f("ZZZ6f",f3);
  H6f* ZZZ_QED6_4 = new H6f("ZZZ6f",f4);
  H6f* Backgr_QED6_Cut_0 = new H6f("BackgrCut",f0);
  H6f* Backgr_QED6_Cut_1 = new H6f("BackgrCut",f1);
  H6f* Backgr_QED6_Cut_2 = new H6f("BackgrCut",f2);
  H6f* Backgr_QED6_Cut_3 = new H6f("BackgrCut",f3);
  H6f* Backgr_QED6_Cut_4 = new H6f("BackgrCut",f4);
  H6f* AllCut_0 = new H6f("AllCut",f0);
  H6f* AllCut_1 = new H6f("AllCut",f1);
  H6f* AllCut_2 = new H6f("AllCut",f2);
  H6f* AllCut_3 = new H6f("AllCut",f3);
  H6f* AllCut_4 = new H6f("AllCut",f4);
  H6f* ZZW_WpZZ = new H6f("ZZW6f",fWp);
  H6f* ZZW_WmZZ = new H6f("ZZW6f",fWm);
  H6f* ZZW_WpZZ_Cut = new H6f("ZZWCut",fWp); 
  H6f* ZZW_WmZZ_Cut = new H6f("ZZWCut",fWm); 
  Hbos* Bosons_QED6_0 = new Hbos("Bosons",f0);
  Hbos* Bosons_QED6_1 = new Hbos("Bosons",f1);
  Hbos* Bosons_QED6_2 = new Hbos("Bosons",f2);
  Hbos* Bosons_QED6_3 = new Hbos("Bosons",f3);
  Hbos* Bosons_QED6_4 = new Hbos("Bosons",f4);
  Hbos* Bosons_WpZZ = new Hbos("Bosons",fWp);
  Hbos* Bosons_WmZZ = new Hbos("Bosons",fWm);
  Hjets* Backgr_Jets_0 = new Hjets("BackgrJets",f0);
  Hjets* Backgr_Jets_1 = new Hjets("BackgrJets",f1);
  Hjets* Backgr_Jets_2 = new Hjets("BackgrJets",f2);
  Hjets* Backgr_Jets_3 = new Hjets("BackgrJets",f3);
  Hjets* Backgr_Jets_4 = new Hjets("BackgrJets",f4);
  Hjets* Backgr_Jets_Cut_0 = new Hjets("BackgrCutJets",f0);
  Hjets* Backgr_Jets_Cut_1 = new Hjets("BackgrCutJets",f1);
  Hjets* Backgr_Jets_Cut_2 = new Hjets("BackgrCutJets",f2);
  Hjets* Backgr_Jets_Cut_3 = new Hjets("BackgrCutJets",f3);
  Hjets* Backgr_Jets_Cut_4 = new Hjets("BackgrCutJets",f4);  
  Hjets* WpZZ_Jets = new Hjets("AllJets",fWp);
  Hjets* WpZZ_Jets_Cut = new Hjets("CutJets",fWp);
  Hjets* WmZZ_Jets = new Hjets("AllJets",fWm);
  Hjets* WmZZ_Jets_Cut = new Hjets("CutJets",fWm);
  TH1F* lostEv_0 = new TH1F("lostEv_0","lostEv_0",3, 0., 3.);
  lostEv_0 = (TH1F*)f0->Get("genAnalyzer/lostEvEtaRange");
  TH1F* lostEv_1 = new TH1F("lostEv_1","lostEv_1",3, 0., 3.);
  lostEv_1 = (TH1F*)f1->Get("genAnalyzer/lostEvEtaRange");
  TH1F* lostEv_2 = new TH1F("lostEv_2","lostEv_2",3, 0., 3.);
  lostEv_2 = (TH1F*)f2->Get("genAnalyzer/lostEvEtaRange");
  TH1F* lostEv_3 = new TH1F("lostEv_3","lostEv_3",3, 0., 3.);
  lostEv_3 = (TH1F*)f3->Get("genAnalyzer/lostEvEtaRange");
  TH1F* lostEv_4 = new TH1F("lostEv_4","lostEv_4",3, 0., 3.);
  lostEv_4 = (TH1F*)f4->Get("genAnalyzer/lostEvEtaRange");

  float L = 300.;

  float xsec_QED6 = 52.76;
  float xsec_Wp = 0.02985;
  float xsec_Wm = 0.01525;

  int Ngen_QED6 = 9998*5;
  int Ngen_Wp = 99989;
  int Ngen_Wm = 9997;

  float w_QED6 = (L*xsec_QED6)/Ngen_QED6;
  float w_WpZZ = (L*xsec_Wp)/Ngen_Wp;
  float w_WmZZ = (L*xsec_Wm)/Ngen_Wm;
 
  cout << "\n%%%%%%%%%%%%\nWEIGHTS\nw_QED6= " << w_QED6 << "\nw_WpZZ= " << w_WpZZ << "\nw_WmZZ=" << w_WmZZ <<endl; 

  ZZW_QED6_0->Scale(w_QED6);
  ZZW_QED6_1->Scale(w_QED6);
  ZZW_QED6_2->Scale(w_QED6);
  ZZW_QED6_3->Scale(w_QED6);
  ZZW_QED6_4->Scale(w_QED6);
  Backgr_QED6_0->Scale(w_QED6);
  Backgr_QED6_1->Scale(w_QED6);
  Backgr_QED6_2->Scale(w_QED6);
  Backgr_QED6_3->Scale(w_QED6);
  Backgr_QED6_4->Scale(w_QED6);
  ZZZ_QED6_0->SetLineColor(kGreen);
  ZZZ_QED6_1->SetLineColor(kGreen);
  ZZZ_QED6_2->SetLineColor(kGreen);
  ZZZ_QED6_3->SetLineColor(kGreen);
  ZZZ_QED6_4->SetLineColor(kGreen);
  ZZZ_QED6_0->Scale(w_QED6);
  ZZZ_QED6_1->Scale(w_QED6);
  ZZZ_QED6_2->Scale(w_QED6);
  ZZZ_QED6_3->Scale(w_QED6);
  ZZZ_QED6_4->Scale(w_QED6);
  Backgr_QED6_Cut_0->Scale(w_QED6);
  Backgr_QED6_Cut_1->Scale(w_QED6);
  Backgr_QED6_Cut_2->Scale(w_QED6);
  Backgr_QED6_Cut_3->Scale(w_QED6);
  Backgr_QED6_Cut_4->Scale(w_QED6);
  ZZW_WpZZ->SetLineColor(kRed);
   ZZW_WpZZ->Scale(w_WpZZ);
  ZZW_WmZZ->SetLineColor(kRed);
  ZZW_WmZZ->Scale(w_WmZZ);
  ZZW_WpZZ_Cut->SetLineColor(kRed); 
  ZZW_WpZZ_Cut->Scale(w_WpZZ);
  ZZW_WmZZ_Cut->SetLineColor(kRed);
  ZZW_WmZZ_Cut->Scale(w_WmZZ);
  Bosons_QED6_0->Scale(w_QED6);
  Bosons_QED6_1->Scale(w_QED6);
  Bosons_QED6_2->Scale(w_QED6);
  Bosons_QED6_3->Scale(w_QED6);
  Bosons_QED6_4->Scale(w_QED6);
  Bosons_WpZZ->SetLineColor(kRed);  
  Bosons_WpZZ->Scale(w_WpZZ);
  Bosons_WmZZ->SetLineColor(kRed);
  Bosons_WmZZ->Scale(w_WmZZ);  
  Backgr_Jets_0->Scale(w_QED6);
  Backgr_Jets_1->Scale(w_QED6);
  Backgr_Jets_2->Scale(w_QED6);
  Backgr_Jets_3->Scale(w_QED6);
  Backgr_Jets_4->Scale(w_QED6);
  Backgr_Jets_Cut_0->Scale(w_QED6);
  Backgr_Jets_Cut_1->Scale(w_QED6);
  Backgr_Jets_Cut_2->Scale(w_QED6);
  Backgr_Jets_Cut_3->Scale(w_QED6);
  Backgr_Jets_Cut_4->Scale(w_QED6);
  WpZZ_Jets->SetLineColor(kGreen);
  WpZZ_Jets->Scale(w_WpZZ);
  WmZZ_Jets->SetLineColor(kGreen);
  WmZZ_Jets->Scale(w_WmZZ);
  WpZZ_Jets_Cut->SetLineColor(kGreen);
  WpZZ_Jets_Cut->Scale(w_WpZZ);
  WmZZ_Jets_Cut->SetLineColor(kGreen);
  WmZZ_Jets_Cut->Scale(w_WmZZ);

  TCanvas* c1 =  new TCanvas("Signal(2f->6f) - Signal(ZZW on-shell) Bosons Mass ","Signal(2f->6f) - Signal(ZZW on-shell) Bosons Mass",900,700);
  c1->Divide(3,1);
  
  TCanvas* c2 = new TCanvas("Signal(2f->6f) - Signal(ZZW on-shell) Mass6f & Pts","Signal(2f->6f) - Signal(ZZW on-shell) Mass6f & Pts", 900, 700);
  c2->Divide(2,2); 

  TCanvas* c3 = new TCanvas("Signal - Background","Signal - Background",900,700);
  c3->Divide(3,1); 

  TCanvas* c4 = new TCanvas("SignalQED6 - ZZW - ZZZ","SignalQED6 - ZZW - ZZZ",900,700);
  c4->Divide(2,2);
  
  TCanvas* c5 = new TCanvas("BackgroundQED6 - ZZW CUT", "BackgroundQED6 - ZZW CUT",900,700);
  c5->Divide(2,2); 

  TCanvas* c6 = new TCanvas("Leptons Pt (lep4)", "Leptons Pt (lep4)",900,700);
  c6->Divide(2,2);
  
  TCanvas* c7 = new TCanvas("Leptons Pt (lep3)", "Leptons Pt (lep3)",900,700);
  c7->Divide(2,2);

  TCanvas* c8 = new TCanvas("Leptons Pt (lep2)", "Leptons Pt (lep2)",900,700);
  c8->Divide(2,2);
  
  TCanvas* c9 = new TCanvas("Leptons Pt (lep1)", "Leptons Pt (lep1)",900,700);
  c9->Divide(2,2);
  
  TCanvas* c10 = new TCanvas("Jets Pt (j2)", "Jets Pt (j2)",900,700);
  c10->Divide(2,2);
  
  TCanvas* c11 = new TCanvas("Jets Pt (j1)", "Jets Pt (j1)",900,700);
  c11->Divide(2,2);
  
  TCanvas* c12 = new TCanvas("Jets Kinematics: Deta", "Jets Kinematics: Deta",900,600);
  c12->Divide(2,1);

  TCanvas* c13 = new TCanvas("Jets Kinematics: Dphi", "Jets Kinematics: Dphi",900,600);
  c13->Divide(2,1);
  
  TCanvas* c14 = new TCanvas("Jets Kinematics: DR", "Jets Kinematics: DR",900,600);
  c14->Divide(2,1);  

  TCanvas* c15 = new TCanvas("Z Pts", "Z Pts",900,600);
  c15->Divide(2,1);

  TLegend* L1 = new TLegend(0.67, 0.61,0.88,0.76);
  TLegend* L2 = new TLegend(0.67, 0.61,0.88,0.76);
  TLegend* L4 = new TLegend(0.67, 0.61,0.88,0.76);
  TLegend* L5 = new TLegend(0.67, 0.61,0.88,0.76);
  TLegend* L6 = new TLegend(0.70, 0.66,0.88,0.76);
  TLegend* L7 = new TLegend(0.70, 0.66,0.88,0.76);
  TLegend* L8 = new TLegend(0.70, 0.66,0.88,0.76);
  TLegend* L9 = new TLegend(0.70, 0.66,0.88,0.76);
  TLegend* L10 = new TLegend(0.70, 0.66,0.88,0.76);
  TLegend* L11 = new TLegend(0.70, 0.66,0.88,0.76);
  TLegend* L12 = new TLegend(0.70, 0.66,0.88,0.76);
  TLegend* L13 = new TLegend(0.70, 0.66,0.88,0.76);
  TLegend* L14 = new TLegend(0.70, 0.66,0.88,0.76);
  TLegend* L15 = new TLegend(0.70, 0.66,0.88,0.76);
 
  ///////////////////////////////////////// Comaprisons 2f->6f - ZZW on-shell //////////////////////////////////

 //Canvas 1: -------------------------Signal(2f->6f) - Signal(ZZW on-shell) Mass ---------------------//

  //----2l0 Mass-----
  c1->cd(1);
  TH1F ZZW_QED6sum2l0 = (*(ZZW_QED6_0->h2l0Mass)+*(ZZW_QED6_1->h2l0Mass)+*(ZZW_QED6_2->h2l0Mass)+*(ZZW_QED6_3->h2l0Mass)+*(ZZW_QED6_4->h2l0Mass));  
  ZZW_QED6sum2l0.SetTitle("2l0Mass");  
  ZZW_QED6sum2l0.DrawClone();
  TH1F Wsum2l0 = (*(ZZW_WpZZ->h2l0Mass)+*(ZZW_WmZZ->h2l0Mass));
  Wsum2l0.DrawClone("same");

  L1->AddEntry(ZZW_QED6_0->h6fMass, "2f->6f", "l");
  L1->AddEntry(ZZW_WpZZ->h6fMass, "ZZW on-shell", "l");
  L1->Draw();
  
  //----2l1 Mass------
  c1->cd(2);
  TH1F ZZW_QED6sum2l1 = (*(ZZW_QED6_0->h2l1Mass)+*(ZZW_QED6_1->h2l1Mass)+*(ZZW_QED6_2->h2l1Mass)+*(ZZW_QED6_3->h2l1Mass)+*(ZZW_QED6_4->h2l1Mass));  
  ZZW_QED6sum2l1.SetTitle("2l1Mass");  
  ZZW_QED6sum2l1.DrawClone();
  TH1F Wsum2l1 = (*(ZZW_WpZZ->h2l1Mass)+*(ZZW_WmZZ->h2l1Mass));
  Wsum2l1.DrawClone("same");
  
  //----jj Mass------
  c1->cd(3);
  TH1F ZZW_QED6sumjj = (*(ZZW_QED6_0->hjjMass)+*(ZZW_QED6_1->hjjMass)+*(ZZW_QED6_2->hjjMass)+*(ZZW_QED6_3->hjjMass)+*(ZZW_QED6_4->hjjMass));  
  ZZW_QED6sumjj.SetTitle("jjMass");  
  ZZW_QED6sumjj.DrawClone();
  TH1F Wsumjj = (*(ZZW_WpZZ->hjjMass)+*(ZZW_WmZZ->hjjMass));
  Wsumjj.DrawClone("same");


  //Canvas 2: -------------------------Signal(2f->6f) - Signal(ZZW on-shell) Mass 6f & Pts ---------------------//
  
  TH1F ZZW_QED6sum6f_in = (*(ZZW_QED6_0->h6fMass)+*(ZZW_QED6_1->h6fMass)+*(ZZW_QED6_2->h6fMass)+*(ZZW_QED6_3->h6fMass)+*(ZZW_QED6_4->h6fMass));  
  float IntegralSignalQED6_in =  ZZW_QED6sum6f_in.Integral(0, ZZW_QED6sum6f_in.GetNbinsX()+1);
  TH1F Wsum6f = (*(ZZW_WpZZ->h6fMass)+*(ZZW_WmZZ->h6fMass));
  float IntegralSignalZZW = Wsum6f.Integral(0, Wsum6f.GetNbinsX()+1);
  
  float scalefactor = IntegralSignalZZW/IntegralSignalQED6_in;

  cout <<"%%%%%%%%%%%%%"<< endl;
  cout << "Integral QED6= " << IntegralSignalQED6_in << endl;
  cout << "Integral ZZW= " << IntegralSignalZZW  << endl;
  cout << "Scale factor= " << scalefactor << endl;

  ZZW_QED6_0->Scale(scalefactor);
  ZZW_QED6_1->Scale(scalefactor);
  ZZW_QED6_2->Scale(scalefactor);
  ZZW_QED6_3->Scale(scalefactor);
  ZZW_QED6_4->Scale(scalefactor);
 
  TH1F ZZW_QED6sum6f = (*(ZZW_QED6_0->h6fMass)+*(ZZW_QED6_1->h6fMass)+*(ZZW_QED6_2->h6fMass)+*(ZZW_QED6_3->h6fMass)+*(ZZW_QED6_4->h6fMass));  
  float IntegralSignalQED6 =  ZZW_QED6sum6f.Integral(0, ZZW_QED6sum6f.GetNbinsX()+1);
  //  TH1F Wsum6f = (*(ZZW_WpZZ->h6fMass)+*(ZZW_WmZZ->h6fMass));
  //float IntegralSignalZZW = Wsum6f.Integral(0, Wsum6f.GetNbinsX()+1);

  cout <<"------ After  scaling ------"<< endl;
  cout << "Integral QED6= " << IntegralSignalQED6 << endl;
  cout << "Integral ZZW= " << IntegralSignalZZW  << endl;
  cout <<"%%%%%%%%%%%%%"<< endl;

  c2->cd(1);
  ZZW_QED6sum6f.SetTitle("6fMass");
  ZZW_QED6sum6f.Rebin(4);
  Wsum6f.Rebin(4);
  ZZW_QED6sum6f.DrawClone();
  Wsum6f.DrawClone("same");

  L2->AddEntry(ZZW_QED6_0->h2l0Pt, "2f->6f", "l");
  L2->AddEntry(ZZW_WpZZ->h2l0Pt, "ZZW on-shell", "l");
  L2->Draw();
  
  c2->cd(2);
  TH1F ZZW_QED6sum2l0Pt = (*(ZZW_QED6_0->h2l0Pt)+*(ZZW_QED6_1->h2l0Pt)+*(ZZW_QED6_2->h2l0Pt)+*(ZZW_QED6_3->h2l0Pt)+*(ZZW_QED6_4->h2l0Pt));  
  TH1F Wsum2l0Pt = (*(ZZW_WpZZ->h2l0Pt)+*(ZZW_WmZZ->h2l0Pt));
  ZZW_QED6sum2l0Pt.SetTitle("2l0Pt"); 
  ZZW_QED6sum2l0Pt.Rebin(4);
  Wsum2l0Pt.Rebin(4);
  ZZW_QED6sum2l0Pt.DrawClone();
  Wsum2l0Pt.DrawClone("same");

  c2->cd(3);
  TH1F ZZW_QED6sum2l1Pt = (*(ZZW_QED6_0->h2l1Pt)+*(ZZW_QED6_1->h2l1Pt)+*(ZZW_QED6_2->h2l1Pt)+*(ZZW_QED6_3->h2l1Pt)+*(ZZW_QED6_4->h2l1Pt));  
  TH1F Wsum2l1Pt = (*(ZZW_WpZZ->h2l1Pt)+*(ZZW_WmZZ->h2l1Pt));
  ZZW_QED6sum2l1Pt.SetTitle("2l1Pt");  
  ZZW_QED6sum2l1Pt.Rebin(4);
  Wsum2l1Pt.Rebin(4);
  ZZW_QED6sum2l1Pt.DrawClone();
  Wsum2l1Pt.DrawClone("same");

  c2->cd(4);
  TH1F ZZW_QED6sumjjPt = (*(ZZW_QED6_0->hjjPt)+*(ZZW_QED6_1->hjjPt)+*(ZZW_QED6_2->hjjPt)+*(ZZW_QED6_3->hjjPt)+*(ZZW_QED6_4->hjjPt));  
  TH1F WsumjjPt = (*(ZZW_WpZZ->hjjPt)+*(ZZW_WmZZ->hjjPt));
  ZZW_QED6sumjjPt.SetTitle("jjPt");  
  ZZW_QED6sumjjPt.Rebin(6);
  WsumjjPt.Rebin(6);
  ZZW_QED6sumjjPt.DrawClone();
  WsumjjPt.DrawClone("same");



  //Canvas 3: ------------------------ Signal - Background 2l0Mass, 2l1Mass, Mass6f ---------------------//

  ZZW_QED6_0->Scale(1/scalefactor);
  ZZW_QED6_1->Scale(1/scalefactor);
  ZZW_QED6_2->Scale(1/scalefactor);
  ZZW_QED6_3->Scale(1/scalefactor);
  ZZW_QED6_4->Scale(1/scalefactor);

  c3->cd(1);
  gPad->SetLogy();
  TH1F Backgr_QED6sum2l0 = (*(Backgr_QED6_0->h2l0Mass)+*(Backgr_QED6_1->h2l0Mass)+*(Backgr_QED6_2->h2l0Mass)+*(Backgr_QED6_3->h2l0Mass)+*(Backgr_QED6_4->h2l0Mass));  
  Backgr_QED6sum2l0.SetTitle("2l0Mass");  
  Backgr_QED6sum2l0.DrawClone();
  // TH1F ZZW_QED6sum2l0 = (*(ZZW_QED6_0->h2l0Mass)+*(ZZW_QED6_1->h2l0Mass)+*(ZZW_QED6_2->h2l0Mass)+*(ZZW_QED6_3->h2l0Mass)+*(ZZW_QED6_4->h2l0Mass));
  ZZW_QED6sum2l0.SetLineColor(kRed);
  ZZW_QED6sum2l0.DrawClone("same");


  c3->cd(2);
  gPad->SetLogy();
  TH1F Backgr_QED6sum2l1 = (*(Backgr_QED6_0->h2l1Mass)+*(Backgr_QED6_1->h2l1Mass)+*(Backgr_QED6_2->h2l1Mass)+*(Backgr_QED6_3->h2l1Mass)+*(Backgr_QED6_4->h2l1Mass));  
  Backgr_QED6sum2l1.SetTitle("2l1Mass");  
  Backgr_QED6sum2l1.DrawClone();
  // TH1F ZZW_QED6sum2l1 = (*(ZZW_QED6_0->h2l1Mass)+*(ZZW_QED6_1->h2l1Mass)+*(ZZW_QED6_2->h2l1Mass)+*(ZZW_QED6_3->h2l1Mass)+*(ZZW_QED6_4->h2l1Mass));  
  ZZW_QED6sum2l1.SetLineColor(kRed);
  ZZW_QED6sum2l1.DrawClone("same");

  c3->cd(3); 
  gPad->SetLogy();
  TH1F Backgr_QED6sum6f = (*(Backgr_QED6_0->h6fMass)+*(Backgr_QED6_1->h6fMass)+*(Backgr_QED6_2->h6fMass)+*(Backgr_QED6_3->h6fMass)+*(Backgr_QED6_4->h6fMass));  
  Backgr_QED6sum6f.SetTitle("6fMass");
  Backgr_QED6sum6f.DrawClone();
  //  TH1F ZZW_QED6sum6f = (*(ZZW_QED6_0->h6fMass)+*(ZZW_QED6_1->h6fMass)+*(ZZW_QED6_2->h6fMass)+*(ZZW_QED6_3->h6fMass)+*(ZZW_QED6_4->h6fMass));
  ZZW_QED6sum6f.SetLineColor(kRed);
  ZZW_QED6sum6f.DrawClone("same");

 

  //Canvas 4: ----------------Signal QED6 - ZZW - ZZZ---------------------//
  c4->cd(1);
  //TH1F ZZW_QED6sum6f = (*(ZZW_QED6_0->h6fMass)+*(ZZW_QED6_1->h6fMass)+*(ZZW_QED6_2->h6fMass)+*(ZZW_QED6_3->h6fMass)+*(ZZW_QED6_4->h6fMass));  
  ZZW_QED6sum6f.SetTitle("6fMass");
  ZZW_QED6sum6f.DrawClone();
  //TH1F Wsum6f = (*(ZZW_WpZZ->h6fMass)+*(ZZW_WmZZ->h6fMass));
  TH1F ZZZ_QED6sum6f = (*(ZZZ_QED6_0->h6fMass)+*(ZZZ_QED6_1->h6fMass)+*(ZZZ_QED6_2->h6fMass)+*(ZZZ_QED6_3->h6fMass)+*(ZZZ_QED6_4->h6fMass));  
  ZZZ_QED6sum6f.DrawClone("same");
  Wsum6f.DrawClone("same");
  
  L4->AddEntry(ZZW_QED6_0->h6fMass, "QED6", "l");
  L4->AddEntry(ZZW_WpZZ->h6fMass, "ZZW", "l");
  L4->AddEntry(ZZZ_QED6_0->h6fMass, "ZZZ", "l");
  L4->Draw();

  c4->cd(2);
  //  TH1F ZZW_QED6sum2l0 = (*(ZZW_QED6_0->h2l0Mass)+*(ZZW_QED6_1->h2l0Mass)+*(ZZW_QED6_2->h2l0Mass)+*(ZZW_QED6_3->h2l0Mass)+*(ZZW_QED6_4->h2l0Mass));  
  ZZW_QED6sum2l0.SetTitle("2l0Mass");  
  ZZW_QED6sum2l0.DrawClone();
  //  TH1F Wsum2l0 = (*(ZZW_WpZZ->h2l0Mass)+*(ZZW_WmZZ->h2l0Mass));
  TH1F ZZZ_QED6sum2l0 = (*(ZZZ_QED6_0->h2l0Mass)+*(ZZZ_QED6_1->h2l0Mass)+*(ZZZ_QED6_2->h2l0Mass)+*(ZZZ_QED6_3->h2l0Mass)+*(ZZZ_QED6_4->h2l0Mass));  
  ZZZ_QED6sum2l0.DrawClone("same");
  Wsum2l0.DrawClone("same");
  
  c4->cd(3);
  // TH1F ZZW_QED6sum2l1 = (*(ZZW_QED6_0->h2l1Mass)+*(ZZW_QED6_1->h2l1Mass)+*(ZZW_QED6_2->h2l1Mass)+*(ZZW_QED6_3->h2l1Mass)+*(ZZW_QED6_4->h2l1Mass));  
  ZZW_QED6sum2l1.SetTitle("2l1Mass");  
  ZZW_QED6sum2l1.DrawClone();
  //  TH1F Wsum2l1 = (*(ZZW_WpZZ->h2l1Mass)+*(ZZW_WmZZ->h2l1Mass));
  TH1F ZZZ_QED6sum2l1 = (*(ZZZ_QED6_0->h2l1Mass)+*(ZZZ_QED6_1->h2l1Mass)+*(ZZZ_QED6_2->h2l1Mass)+*(ZZZ_QED6_3->h2l1Mass)+*(ZZZ_QED6_4->h2l1Mass));  
  ZZZ_QED6sum2l1.DrawClone("same");
  Wsum2l1.DrawClone("same");
  
  c4->cd(4);
  //  TH1F ZZW_QED6sumjj = (*(ZZW_QED6_0->hjjMass)+*(ZZW_QED6_1->hjjMass)+*(ZZW_QED6_2->hjjMass)+*(ZZW_QED6_3->hjjMass)+*(ZZW_QED6_4->hjjMass));  
  ZZW_QED6sumjj.SetTitle("jjMass");  
  ZZW_QED6sumjj.DrawClone();
  //  TH1F Wsumjj = (*(ZZW_WpZZ->hjjMass)+*(ZZW_WmZZ->hjjMass));
  TH1F ZZZ_QED6sumjj = (*(ZZZ_QED6_0->hjjMass)+*(ZZZ_QED6_1->hjjMass)+*(ZZZ_QED6_2->hjjMass)+*(ZZZ_QED6_3->hjjMass)+*(ZZZ_QED6_4->hjjMass));  
  ZZZ_QED6sumjj.DrawClone("same");
  Wsumjj.DrawClone("same");
  
  
//Canvas 5: --------------Background QED6 - ZZW (produced on shell)  CUT----------- m_6f>300 && Pt_lep >7 && Pt_j1>40 && Pt_j2>25-----------//

  c5->cd(1);
  TH1F Cut_QED6sum6f = (*(Backgr_QED6_Cut_0->h6fMass)+*(Backgr_QED6_Cut_1->h6fMass)+*(Backgr_QED6_Cut_2->h6fMass)+*(Backgr_QED6_Cut_3->h6fMass)+*(Backgr_QED6_Cut_4->h6fMass));  
  Cut_QED6sum6f.SetTitle("6fMass");  
  Cut_QED6sum6f.DrawClone();
  TH1F Wsum6fCut = (*(ZZW_WpZZ_Cut->h6fMass)+*(ZZW_WmZZ_Cut->h6fMass));
  Wsum6fCut.DrawClone("same");

  L5->AddEntry(Backgr_QED6_Cut_0->h6fMass, "B", "l");
  L5->AddEntry(ZZW_WpZZ_Cut->h6fMass, "S", "l");
  L5->Draw();
  
  c5->cd(2);
  TH1F Cut_QED6sum2l0 = (*(Backgr_QED6_Cut_0->h2l0Mass)+*(Backgr_QED6_Cut_1->h2l0Mass)+*(Backgr_QED6_Cut_2->h2l0Mass)+*(Backgr_QED6_Cut_3->h2l0Mass)+*(Backgr_QED6_Cut_4->h2l0Mass));  
  Cut_QED6sum2l0.SetTitle("2l0Mass");  
  Cut_QED6sum2l0.DrawClone();
  TH1F Wsum2l0Cut = (*(ZZW_WpZZ_Cut->h2l0Mass)+*(ZZW_WmZZ_Cut->h2l0Mass));
  Wsum2l0Cut.DrawClone("same");
  
  c5->cd(3);
  TH1F Cut_QED6sum2l1= (*(Backgr_QED6_Cut_0->h2l1Mass)+*(Backgr_QED6_Cut_1->h2l1Mass)+*(Backgr_QED6_Cut_2->h2l1Mass)+*(Backgr_QED6_Cut_3->h2l1Mass)+*(Backgr_QED6_Cut_4->h2l1Mass));  
  Cut_QED6sum2l1.SetTitle("2l1Mass");  
  Cut_QED6sum2l1.DrawClone();
  TH1F Wsum2l1Cut = (*(ZZW_WpZZ_Cut->h2l1Mass)+*(ZZW_WmZZ_Cut->h2l1Mass));
  Wsum2l1Cut.DrawClone("same");
  
  c5->cd(4);
  TH1F Cut_QED6sumjj= (*(Backgr_QED6_Cut_0->hjjMass)+*(Backgr_QED6_Cut_1->hjjMass)+*(Backgr_QED6_Cut_2->hjjMass)+*(Backgr_QED6_Cut_3->hjjMass)+*(Backgr_QED6_Cut_4->hjjMass));  
  Cut_QED6sumjj.SetTitle("jjMass");  
  Cut_QED6sumjj.DrawClone();
  TH1F WsumjjCut = (*(ZZW_WpZZ_Cut->hjjMass)+*(ZZW_WmZZ_Cut->hjjMass));
  WsumjjCut.DrawClone("same");
  

  float errZZW_QED6_0  = sqrt(ZZW_QED6_0->h6fMass->GetEntries())*w_QED6;
  float errZZW_QED6_1  = sqrt(ZZW_QED6_1->h6fMass->GetEntries())*w_QED6;
  float errZZW_QED6_2  = sqrt(ZZW_QED6_2->h6fMass->GetEntries())*w_QED6;
  float errZZW_QED6_3  = sqrt(ZZW_QED6_3->h6fMass->GetEntries())*w_QED6;
  float errZZW_QED6_4  = sqrt(ZZW_QED6_4->h6fMass->GetEntries())*w_QED6;
  float errBkgr_QED6_0 = sqrt(Backgr_QED6_0->h6fMass->GetEntries())*w_QED6;
  float errBkgr_QED6_1 = sqrt(Backgr_QED6_1->h6fMass->GetEntries())*w_QED6;
  float errBkgr_QED6_2 = sqrt(Backgr_QED6_2->h6fMass->GetEntries())*w_QED6;
  float errBkgr_QED6_3 = sqrt(Backgr_QED6_3->h6fMass->GetEntries())*w_QED6;
  float errBkgr_QED6_4 = sqrt(Backgr_QED6_4->h6fMass->GetEntries())*w_QED6;
  float errCut_QED6_0  = sqrt(Backgr_QED6_Cut_0->h6fMass->GetEntries())*w_QED6;
  float errCut_QED6_1  = sqrt(Backgr_QED6_Cut_1->h6fMass->GetEntries())*w_QED6;
  float errCut_QED6_2  = sqrt(Backgr_QED6_Cut_2->h6fMass->GetEntries())*w_QED6;
  float errCut_QED6_3  = sqrt(Backgr_QED6_Cut_3->h6fMass->GetEntries())*w_QED6;
  float errCut_QED6_4  = sqrt(Backgr_QED6_Cut_4->h6fMass->GetEntries())*w_QED6;
  float errWp = sqrt(ZZW_WpZZ->h6fMass->GetEntries())*w_WpZZ;
  float errWm = sqrt(ZZW_WmZZ->h6fMass->GetEntries())*w_WmZZ;
  float errWp_Cut = sqrt(ZZW_WpZZ_Cut->h6fMass->GetEntries())*w_WpZZ;
  float errWm_Cut = sqrt(ZZW_WmZZ_Cut->h6fMass->GetEntries())*w_WmZZ;
  
  
  cout << "\n%%%%%%%%%%%%\nAFTER SCALING with WEIGHTS" <<endl;

  cout << "\nSIGNAL_2l0= " << ZZW_QED6sum2l0.Integral(0, ZZW_QED6sum2l0.GetNbinsX()+1) 
       << "\nSIGNAL_2l1= " << ZZW_QED6sum2l1.Integral(0, ZZW_QED6sum2l1.GetNbinsX()+1)
       << "\t\tErrorSignal= " << sqrt(errZZW_QED6_0*errZZW_QED6_0 + errZZW_QED6_1*errZZW_QED6_1 + errZZW_QED6_2*errZZW_QED6_2 + errZZW_QED6_3*errZZW_QED6_3 + errZZW_QED6_4*errZZW_QED6_4)
       << "\nBackground = " << Backgr_QED6sum6f.Integral(0,Backgr_QED6sum6f.GetNbinsX()+1) 
       << "\t\tError = " << sqrt(errBkgr_QED6_0*errBkgr_QED6_0 + errBkgr_QED6_1*errBkgr_QED6_1 + errBkgr_QED6_2*errBkgr_QED6_2 + errBkgr_QED6_3*errBkgr_QED6_3 + errBkgr_QED6_4*errBkgr_QED6_4) 
       << "\nBackground_Cut= " << Cut_QED6sum6f.Integral(0,Cut_QED6sum6f.GetNbinsX()+1) 
       << "\t\tError= " << sqrt(errCut_QED6_0*errCut_QED6_0 + errCut_QED6_1*errCut_QED6_1 + errCut_QED6_2*errCut_QED6_2 + errCut_QED6_3*errCut_QED6_3 + errCut_QED6_4*errCut_QED6_4)    
       <<"\nSIGNAL_onshell= " << Wsum6f.Integral(0, Wsum6f.GetNbinsX()+1) 
       << "\t\t\tErrorSignalOnShell= " << sqrt(errWp*errWp + errWm*errWm)    
       << "\nSignal_Cut= " << Wsum6fCut.Integral(0, Wsum6fCut.GetNbinsX()+1)
       << "\t\tError= " << sqrt(errWp_Cut*errWp_Cut + errWm_Cut*errWm_Cut)  << endl;
  


  //////////////////////////////////////////////////////////
///////////===========Pt cut analysis=============////////////


  TH1F* bl4 = (TH1F*)(*(Backgr_QED6_Cut_0->h4lPt_4)+*(Backgr_QED6_Cut_1->h4lPt_4)+*(Backgr_QED6_Cut_2->h4lPt_4)+*(Backgr_QED6_Cut_3->h4lPt_4)+*(Backgr_QED6_Cut_4->h4lPt_4));
  TH1F* sl4 = (TH1F*)(*(ZZW_WpZZ_Cut->h4lPt_4)+*(ZZW_WmZZ_Cut->h4lPt_4));

  TH1F* bl3 = (TH1F*)(*(Backgr_QED6_Cut_0->h4lPt_3)+*(Backgr_QED6_Cut_1->h4lPt_3)+*(Backgr_QED6_Cut_2->h4lPt_3)+*(Backgr_QED6_Cut_3->h4lPt_3)+*(Backgr_QED6_Cut_4->h4lPt_3));
  TH1F* sl3 = (TH1F*)(*(ZZW_WpZZ_Cut->h4lPt_3)+*(ZZW_WmZZ_Cut->h4lPt_3));
  
  TH1F* bl2 = (TH1F*)(*(Backgr_QED6_Cut_0->h4lPt_2)+*(Backgr_QED6_Cut_1->h4lPt_2)+*(Backgr_QED6_Cut_2->h4lPt_2)+*(Backgr_QED6_Cut_3->h4lPt_2)+*(Backgr_QED6_Cut_4->h4lPt_2));
  TH1F* sl2 = (TH1F*)(*(ZZW_WpZZ_Cut->h4lPt_2)+*(ZZW_WmZZ_Cut->h4lPt_2));
  
  TH1F* bl1 = (TH1F*)(*(Backgr_QED6_Cut_0->h4lPt_1)+*(Backgr_QED6_Cut_1->h4lPt_1)+*(Backgr_QED6_Cut_2->h4lPt_1)+*(Backgr_QED6_Cut_3->h4lPt_1)+*(Backgr_QED6_Cut_4->h4lPt_1));
  TH1F* sl1 = (TH1F*)(*(ZZW_WpZZ_Cut->h4lPt_1)+*(ZZW_WmZZ_Cut->h4lPt_1));
  
  TH1F* bj2 = (TH1F*)(*(Backgr_QED6_Cut_0->hjPt_2)+*(Backgr_QED6_Cut_1->hjPt_2)+*(Backgr_QED6_Cut_2->hjPt_2)+*(Backgr_QED6_Cut_3->hjPt_2)+*(Backgr_QED6_Cut_4->hjPt_2));
  TH1F* sj2 = (TH1F*)(*(ZZW_WpZZ_Cut->hjPt_2)+*(ZZW_WmZZ_Cut->hjPt_2));
  
  TH1F* bj1 = (TH1F*)(*(Backgr_QED6_Cut_0->hjPt_1)+*(Backgr_QED6_Cut_1->hjPt_1)+*(Backgr_QED6_Cut_2->hjPt_1)+*(Backgr_QED6_Cut_3->hjPt_1)+*(Backgr_QED6_Cut_4->hjPt_1));
  TH1F* sj1 = (TH1F*)(*(ZZW_WpZZ_Cut->hjPt_1)+*(ZZW_WmZZ_Cut->hjPt_1));


//%%%%%%%% Canvas 6: ------QED6_4lPt (lep4) - ZZW_4lPt------%%%%%%%%
  c6->cd(1);
  bl4->SetTitle("B (QED6) - S (ZZW)");
  bl4->SetName("B (QED6) - S (ZZW)");
  bl4->DrawClone();  
  sl4->DrawClone("same");
  

  L6->AddEntry(Backgr_QED6_Cut_0->h4lPt_4, "B (QED6)", "l");
  L6->AddEntry(ZZW_WpZZ_Cut->h4lPt_4, "S (ZZW)", "l");
  L6->Draw();

//-----------Integral functions---------
  c6->cd(2);
  TH1F* int_bl4 = Integral(bl4);
  TH1F* int_sl4 = Integral(sl4);  
  int_bl4->SetTitle("Int(B) - Int(S)");
  int_bl4->SetName("Int(B) - Int(S)");
  int_bl4->Draw();
  int_sl4->Draw("same");

//-----------S/sqrt(B)------------------ 
  c6->cd(3);  
  TH1F* cut0 = Divide(int_sl4,int_bl4);  
  cut0->SetTitle("S/sqrt(B)");
  cut0->SetName("S/sqrt(B)");
  cut0->Draw();

//----------S/sqrt(S+B)-----------------
  c6->cd(4);  
  TH1F* int_sbl4sum = (TH1F*)(*(int_sl4)+*(int_bl4));
  TH1F* cut1 = Divide(int_sl4,int_sbl4sum);  
  cut1->SetTitle("S/sqrt(S+B)");
  cut1->SetName("S/sqrt(S+B)");
  cut1->Draw();


//%%%%%%%% Canvas 7: ------QED6_4lPt (lep3) - ZZW_4lPt------%%%%%%%%

  c7->cd(1);
  bl3->SetTitle("B (QED6) - S (ZZW)");
  bl3->SetName("B (QED6) - S (ZZW)");
  bl3->DrawClone();  
  sl3->DrawClone("same");
  

  L7->AddEntry(Backgr_QED6_Cut_0->h4lPt_3, "B (QED6)", "l");
  L7->AddEntry(ZZW_WpZZ_Cut->h4lPt_3, "S (ZZW)", "l");
  L7->Draw();

//-----------Integral functions---------
  c7->cd(2);
  TH1F* int_bl3 = Integral(bl3);
  TH1F* int_sl3 = Integral(sl3);  
  int_bl3->SetTitle("Int(B) - Int(S)");
  int_bl3->SetName("Int(B) - Int(S)");
  int_bl3->Draw();
  int_sl3->Draw("same");

//-----------S/sqrt(B)------------------ 
  c7->cd(3);  
  TH1F* cut2 = Divide(int_sl3,int_bl3);  
  cut2->SetTitle("S/sqrt(B)");
  cut2->SetName("S/sqrt(B)");
  cut2->Draw();

//----------S/sqrt(S+B)-----------------
  c7->cd(4);  
  TH1F* int_sbl3sum = (TH1F*)(*(int_sl3)+*(int_bl3));
  TH1F* cut3 = Divide(int_sl3,int_sbl3sum);  
  cut3->SetTitle("S/sqrt(S+B)");
  cut3->SetName("S/sqrt(S+B)");
  cut3->Draw();



//%%%%%%%% Canvas 8: ------QED6_4lPt (lep2) - ZZW_4lPt------%%%%%%%%

  c8->cd(1);
  bl2->SetTitle("B (QED6) - S (ZZW)");
  bl2->SetName("B (QED6) - S (ZZW)");
  bl2->DrawClone();  
  sl2->DrawClone("same");
  

  L8->AddEntry(Backgr_QED6_Cut_0->h4lPt_2, "B (QED6)", "l");
  L8->AddEntry(ZZW_WpZZ_Cut->h4lPt_2, "S (ZZW)", "l");
  L8->Draw();

//-----------Integral functions---------
  c8->cd(2);
  TH1F* int_bl2 = Integral(bl2);
  TH1F* int_sl2 = Integral(sl2);  
  int_bl2->SetTitle("Int(B) - Int(S)");
  int_bl2->SetName("Int(B) - Int(S)");
  int_bl2->Draw();
  int_sl2->Draw("same");

//-----------S/sqrt(B)------------------ 
  c8->cd(3);  
  TH1F* cut4 = Divide(int_sl2,int_bl2);  
  cut4->SetTitle("S/sqrt(B)");
  cut4->SetName("S/sqrt(B)");
  cut4->Draw();

//----------S/sqrt(S+B)-----------------
  c8->cd(4);  
  TH1F* int_sbl2sum = (TH1F*)(*(int_sl2)+*(int_bl2));
  TH1F* cut5 = Divide(int_sl2,int_sbl2sum);  
  cut5->SetTitle("S/sqrt(S+B)");
  cut5->SetName("S/sqrt(S+B)");
  cut5->Draw();


//%%%%%%%% Canvas 9: ------QED6_4lPt (lep1) - ZZW_4lPt------%%%%%%%%

  c9->cd(1);
  bl1->SetTitle("B (QED6) - S (ZZW)");
  bl1->SetName("B (QED6) - S (ZZW)");
  bl1->DrawClone();  
  sl1->DrawClone("same");
  

  L9->AddEntry(Backgr_QED6_Cut_0->h4lPt_1, "B (QED6)", "l");
  L9->AddEntry(ZZW_WpZZ_Cut->h4lPt_1, "S (ZZW)", "l");
  L9->Draw();

//-----------Integral functions---------
  c9->cd(2);
  TH1F* int_bl1 = Integral(bl1);
  TH1F* int_sl1 = Integral(sl1);  
  int_bl1->SetTitle("Int(B) - Int(S)");
  int_bl1->SetName("Int(B) - Int(S)");
  int_bl1->Draw();
  int_sl1->Draw("same");

//-----------S/sqrt(B)------------------ 
  c9->cd(3);  
  TH1F* cut6 = Divide(int_sl1,int_bl1);  
  cut6->SetTitle("S/sqrt(B)");
  cut6->SetName("S/sqrt(B)");
  cut6->Draw();

//----------S/sqrt(S+B)-----------------
  c9->cd(4);  
  TH1F* int_sbl1sum = (TH1F*)(*(int_sl1)+*(int_bl1));
  TH1F* cut7 = Divide(int_sl1,int_sbl1sum);  
  cut7->SetTitle("S/sqrt(S+B)");
  cut7->SetName("S/sqrt(S+B)");
  cut7->Draw();


//%%%%%%%% Canvas 10: ------QED6_4lPt (j2) - ZZW_4lPt------%%%%%%%%

  c10->cd(1);
  bj2->SetTitle("B (QED6) - S (ZZW)");
  bj2->SetName("B (QED6) - S (ZZW)");
  bj2->DrawClone();  
  sj2->DrawClone("same");
  

  L10->AddEntry(Backgr_QED6_Cut_0->hjPt_2, "B (QED6)", "l");
  L10->AddEntry(ZZW_WpZZ_Cut->hjPt_2, "S (ZZW)", "l");
  L10->Draw();

//-----------Integral functions---------
  c10->cd(2);
  TH1F* int_bj2 = Integral(bj2);
  TH1F* int_sj2 = Integral(sj2);  
  int_bj2->SetTitle("Int(B) - Int(S)");
  int_bj2->SetName("Int(B) - Int(S)");
  int_bj2->Draw();
  int_sj2->Draw("same");

//-----------S/sqrt(B)------------------ 
  c10->cd(3);  
  TH1F* cut4 = Divide(int_sj2,int_bj2);  
  cut4->SetTitle("S/sqrt(B)");
  cut4->SetName("S/sqrt(B)");
  cut4->Draw();

//----------S/sqrt(S+B)-----------------
  c10->cd(4);  
  TH1F* int_sbj2sum = (TH1F*)(*(int_sj2)+*(int_bj2));
  TH1F* cut5 = Divide(int_sj2,int_sbj2sum);  
  cut5->SetTitle("S/sqrt(S+B)");
  cut5->SetName("S/sqrt(S+B)");
  cut5->Draw();



//%%%%%%%% Canvas 11: ------QED6_4lPt (j1) - ZZW_4lPt------%%%%%%%%

  c11->cd(1);
  bj1->SetTitle("B (QED6) - S (ZZW)");
  bj1->SetName("B (QED6) - S (ZZW)");
  bj1->DrawClone();  
  sj1->DrawClone("same");
  

  L11->AddEntry(Backgr_QED6_Cut_0->hjPt_1, "B (QED6)", "l");
  L11->AddEntry(ZZW_WpZZ_Cut->hjPt_1, "S (ZZW)", "l");
  L11->Draw();

//-----------Integral functions---------
  c11->cd(2);
  TH1F* int_bj1 = Integral(bj1);
  TH1F* int_sj1 = Integral(sj1);  
  int_bj1->SetTitle("Int(B) - Int(S)");
  int_bj1->SetName("Int(B) - Int(S)");
  int_bj1->Draw();
  int_sj1->Draw("same");

//-----------S/sqrt(B)------------------ 
  c11->cd(3);  
  TH1F* cut6 = Divide(int_sj1,int_bj1);  
  cut6->SetTitle("S/sqrt(B)");
  cut6->SetName("S/sqrt(B)");
  cut6->Draw();

//----------S/sqrt(S+B)-----------------
  c11->cd(4);  
  TH1F* int_sbj1sum = (TH1F*)(*(int_sj1)+*(int_bj1));
  TH1F* cut7 = Divide(int_sj1,int_sbj1sum);  
  cut7->SetTitle("S/sqrt(S+B)");
  cut7->SetName("S/sqrt(S+B)");
  cut7->Draw();



  ///////////////////////////////////////////////////////
///////////===========Jets angular analysis=============////////////

TH1F* bDeta =(TH1F*)(*(Backgr_Jets_0->hjjDeta)+*(Backgr_Jets_1->hjjDeta)+*(Backgr_Jets_2->hjjDeta)+*(Backgr_Jets_3->hjjDeta)+*(Backgr_Jets_4->hjjDeta));
TH1F* sDeta =(TH1F*)(*(WpZZ_Jets->hjjDeta)+*(WmZZ_Jets->hjjDeta));

TH1F* bDphi =(TH1F*)(*(Backgr_Jets_0->hjjDphi)+*(Backgr_Jets_1->hjjDphi)+*(Backgr_Jets_2->hjjDphi)+*(Backgr_Jets_3->hjjDphi)+*(Backgr_Jets_4->hjjDphi));
TH1F* sDphi =(TH1F*)(*(WpZZ_Jets->hjjDphi)+*(WmZZ_Jets->hjjDphi));

TH1F* bDR =(TH1F*)(*(Backgr_Jets_0->hjjDR)+*(Backgr_Jets_1->hjjDR)+*(Backgr_Jets_2->hjjDR)+*(Backgr_Jets_3->hjjDR)+*(Backgr_Jets_4->hjjDR));
TH1F* sDR =(TH1F*)(*(WpZZ_Jets->hjjDR)+*(WmZZ_Jets->hjjDR));

TH1F* bDeta_Cut =(TH1F*)(*(Backgr_Jets_Cut_0->hjjDeta)+*(Backgr_Jets_Cut_1->hjjDeta)+*(Backgr_Jets_Cut_2->hjjDeta)+*(Backgr_Jets_Cut_3->hjjDeta)+*(Backgr_Jets_Cut_4->hjjDeta));
TH1F* sDeta_Cut =(TH1F*)(*(WpZZ_Jets_Cut->hjjDeta)+*(WmZZ_Jets_Cut->hjjDeta));

TH1F* bDphi_Cut =(TH1F*)(*(Backgr_Jets_Cut_0->hjjDphi)+*(Backgr_Jets_Cut_1->hjjDphi)+*(Backgr_Jets_Cut_2->hjjDphi)+*(Backgr_Jets_Cut_3->hjjDphi)+*(Backgr_Jets_Cut_4->hjjDphi));
TH1F* sDphi_Cut =(TH1F*)(*(WpZZ_Jets_Cut->hjjDphi)+*(WmZZ_Jets_Cut->hjjDphi));

TH1F* bDR_Cut =(TH1F*)(*(Backgr_Jets_Cut_0->hjjDR)+*(Backgr_Jets_Cut_1->hjjDR)+*(Backgr_Jets_Cut_2->hjjDR)+*(Backgr_Jets_Cut_3->hjjDR)+*(Backgr_Jets_Cut_4->hjjDR));
TH1F* sDR_Cut =(TH1F*)(*(WpZZ_Jets_Cut->hjjDR)+*(WmZZ_Jets_Cut->hjjDR));


//Canvas 12: ----------------Deta-----------------

 c12->cd(1);
 bDeta->SetTitle("Deta");
 bDeta->DrawClone();
 sDeta->DrawClone("same");


//------Cut--------
 c12->cd(2);
 bDeta_Cut->SetTitle("Deta_Cut");
 bDeta_Cut->DrawClone();
 sDeta_Cut->DrawClone("same");

 L12->AddEntry(Backgr_Jets_0->hjjDeta, "B", "l");
 L12->AddEntry(WpZZ_Jets->hjjDeta, "S", "l");
 L12->Draw();


//Canvas 13: ----------------Dphi----------------

 c13->cd(1);
 bDphi->SetTitle("Dphi");
 bDphi->DrawClone();
 sDphi->DrawClone("same");

//------Cut--------
 c13->cd(2);
 bDphi_Cut->SetTitle("Dphi_Cut");
 bDphi_Cut->DrawClone();
 sDphi_Cut->DrawClone("same");

 L13->AddEntry(Backgr_Jets_0->hjjDphi, "B", "l");
 L13->AddEntry(WpZZ_Jets->hjjDphi, "S", "l");
 L13->Draw();

//Canvas 14: ----------------DR------------------

 c14->cd(1);
 bDR->SetTitle("DR");
 bDR->DrawClone();
 sDR->DrawClone("same");

//------Cut--------
 c14->cd(2);
 bDR_Cut->SetTitle("DR_Cut");
 bDR_Cut->DrawClone();
 sDR_Cut->DrawClone("same");

 L14->AddEntry(Backgr_Jets_0->hjjDR, "B", "l");
 L14->AddEntry(WpZZ_Jets->hjjDR, "S", "l");
 L14->Draw();


 //////===========Bosons=============////////////

 //Canvas 15: ----------------ZZPts------------------

 TH1F* QED6ZPt_1 =(TH1F*)(*(Bosons_QED6_0->hZPt_1)+*(Bosons_QED6_1->hZPt_1)+*(Bosons_QED6_2->hZPt_1)+*(Bosons_QED6_3->hZPt_1)+*(Bosons_QED6_4->hZPt_1));
 TH1F* QED6ZPt_2 =(TH1F*)(*(Bosons_QED6_0->hZPt_2)+*(Bosons_QED6_1->hZPt_2)+*(Bosons_QED6_2->hZPt_2)+*(Bosons_QED6_3->hZPt_2)+*(Bosons_QED6_4->hZPt_2));
 TH1F* WZZ_ZPt_1 =(TH1F*)(*(Bosons_WpZZ->hZPt_1)+*(Bosons_WmZZ->hZPt_1));
 TH1F* WZZ_ZPt_2 =(TH1F*)(*(Bosons_WpZZ->hZPt_2)+*(Bosons_WmZZ->hZPt_2));
 
 c15->cd(1);
 QED6ZPt_1->SetTitle("ZPt_1");
 QED6ZPt_1->DrawClone();
 WZZ_ZPt_1->DrawClone("same"); 

 L15->AddEntry(Bosons_QED6_0->hZPt_1, "QED6", "l");
 L15->AddEntry(Bosons_WpZZ->hZPt_1, "ZZW", "l");
 L15->Draw();
 
 c15->cd(2);
 QED6ZPt_2->SetTitle("ZPt_2");
 QED6ZPt_2->DrawClone();
 WZZ_ZPt_2->DrawClone("same");

 ///////////===========ETA CUT:LOST EVENTS and SIGNAL=============////////////

 int lostev = (lostEv_0->GetEntries())+(lostEv_1->GetEntries())+(lostEv_2->GetEntries())+(lostEv_3->GetEntries())+(lostEv_4->GetEntries());

 cout << "\n%%%%%%%%%%%% \nETA CUT for leptons: -2.5 < eta < 2.5" << endl;
 cout << "Generated events: " << Ngen_QED6 << "\nLost events: " << lostev << endl;
 

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
