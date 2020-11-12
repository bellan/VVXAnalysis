/*
Macro to analyze efficiency and resolution for a VZZ sample analyzed by VZZaQGCAnalyzer. Path will probably need changing.
Author: Marozzo Giovanni Battista
Date: 2020/11/25
*/

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define TEST_MODE

using namespace std;

void efficiency(string histname1,string histname2,string axistitle, string title, string efftitle, double ymax, TFile *result, TCanvas *c){
  TVirtualPad *pd1=c->cd(1);
  TH1F *hist1= (TH1F*)result->Get(histname1.c_str());
  TH1F *hist2= (TH1F*)result->Get(histname2.c_str());
  TH1F *hist3= (TH1F*)hist2->Clone("efficiency");
  
  hist1->SetTitle(title.c_str());
  hist1->GetYaxis()->SetRangeUser(0,ymax);
  hist1->GetXaxis()->SetTitle(axistitle.c_str());
  hist1->SetLineColor(1);
  hist2->SetLineColor(2);
  hist1->SetLineWidth(2);
  hist2->SetLineWidth(2);
  hist1->Draw();
  hist2->Draw("same");
  TLegend *legend= new TLegend(0.7,0.77,0.98,0.94,"");
  legend->AddEntry(hist1,"generated bosons","l");
  legend->AddEntry(hist2,"reconstructed bosons","l");
  legend->Draw();

  TVirtualPad *pd2=c->cd(2);
  hist3->Divide(hist1);
  hist3->SetTitle(efftitle.c_str());
  hist3->GetYaxis()->SetRangeUser(0,1);
  hist3->GetXaxis()->SetTitle(axistitle.c_str());
  hist3->GetYaxis()->SetTitle("efficiency");
  hist3->SetLineColor(1);
  hist3->SetLineWidth(2);
  hist3->Draw();}


void VZZaQGCEfficiency(){

  gErrorIgnoreLevel = kFatal; //to ignore errors due to nonexistent tree branches

  TFile *result = TFile::Open("./results/2018/VZZaQGCAnalyzer_MC/WZZ.root");
  TCanvas *c1 = new TCanvas("c1","canvas",0,0,1000,500);
  c1->Divide(2,1);
  TCanvas *c2 = new TCanvas("c2","canvas",0,0,1000,500);
  c2->Divide(2,1);
  TCanvas *c3 = new TCanvas("c3","canvas",0,0,1000,500);
  c3->Divide(2,1);
  TCanvas *c4 = new TCanvas("c4","canvas",0,0,1000,500);
  c4->Divide(2,1);
  TCanvas *c5 = new TCanvas("c5","canvas",0,0,1000,500);
  c5->Divide(2,1);
  TCanvas *c6 = new TCanvas("c6","canvas",0,0,1000,500);
  c6->Divide(2,1);
  TCanvas *c7 = new TCanvas("c7","canvas",0,0,1000,500);
  c7->Divide(2,1);
  TCanvas *c8 = new TCanvas("c8","canvas",0,0,1000,500);
  c8->Divide(2,1);
  TCanvas *c9 = new TCanvas("c9","canvas",0,0,1000,500);
  c9->Divide(2,1);
  TCanvas *c10 = new TCanvas("c10","canvas",0,0,1000,500);
  c10->Divide(2,1);
  TCanvas *c11 = new TCanvas("c11","canvas",0,0,1000,500);
  c11->Divide(2,1);
  TCanvas *c12 = new TCanvas("c12","canvas",0,0,1000,500);
  c12->Divide(2,1);

  efficiency("mass of well generated Z1","mass of well reconstructed Z1","mass (Gev/c^2)","Generated and reconstructed Z1 mass","Efficiency as a function of Z1 mass",1000,result,c1);
  
  efficiency("mass of well generated Z2","mass of well reconstructed Z2","mass (Gev/c^2)","Generated and reconstructed Z2 mass","efficiency as a function of Z2 mass",1000,result,c2);
  
  efficiency("pt of well generated Z1","pt of well reconstructed Z1","pt (Gev/c)","Generated and reconstructed Z1 pt","efficiency as a function of Z1 pt",600,result,c3);
  
  efficiency("pt of well generated Z2","pt of well reconstructed Z2","pt (Gev/c)","Generated and reconstructes Z2 pt","efficiency as a function of Z2 pt",600,result,c4);
  
  efficiency("energy of well generated Z1","energy of well reconstructed Z1","E (Gev)","Generated and reconstructed Z1 energy","efficiency as a function of Z1 energy",1200,result,c5);
  
  efficiency("energy of well generated Z2","energy of well reconstructed Z2","E (Gev)","Generated and reconstructed Z2 energy","efficiency as a function of Z2 energy",1200,result,c6);
  
  efficiency("eta of well generated Z1","eta of well reconstructed Z1","eta","Generated and reconstructed Z1 eta","efficiency as a function of Z1 eta",200,result,c7);
  
  efficiency("eta of well generated Z2","eta of well reconstructed Z2","eta","Generated and reconstructed Z2 eta","efficiency as a function of Z2 eta",200,result,c8);
  
  efficiency("energy of good Z1 major lepton","energy of well reconstructed Z1 major lepton","E (Gev)","Generated and reconstructed Z1 major lepton energy","efficiency as a function of Z1 major lepton energy",400,result,c9);
  
  efficiency("energy of good Z2 major lepton","energy of well reconstructed Z2 major lepton","E (Gev)","Generated and reconstructed Z2 major lepton energy","efficiency as a function of Z2 major lepton energy",400,result,c10);
  
  efficiency("energy of good Z1 minor lepton","energy of well reconstructed Z1 minor lepton","E (Gev)","Generated and reconstructed Z1 minor lepton energy","efficiency as a function of Z1 minor lepton energy",400,result,c11);
  
  efficiency("energy of good Z2 minor lepton","energy of well reconstructed Z2 minor lepton","E (Gev)","Generated and reconstructed Z2 minor lepton energy","efficiency as a function of Z2 minor lepton energy",400,result,c12);
}
