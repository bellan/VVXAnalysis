#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define TEST_MODE

using namespace std;

void VZZaQGCEfficiency(){

  TFile *result = TFile::Open("/home/giovanni/Desktop/tesi/VVXAnalysis/TreeAnalysis/results/2018/VZZaQGCAnalyzer_MC/WZZ.root");
  TCanvas *c1 = new TCanvas("c1","canvas",0,0,1000,1000);
  TCanvas *c2 = new TCanvas("c2","canvas",0,0,1000,1000);
  TCanvas *c3 = new TCanvas("c3","canvas",0,0,1000,1000);
  TCanvas *c4 = new TCanvas("c4","canvas",0,0,1000,1000);
  TCanvas *c5 = new TCanvas("c5","canvas",0,0,1000,1000);
  TCanvas *c6 = new TCanvas("c6","canvas",0,0,1000,1000);
  TCanvas *c7 = new TCanvas("c7","canvas",0,0,1000,1000);
  TCanvas *c8 = new TCanvas("c8","canvas",0,0,1000,1000);
  
  c1->cd(); 
  TH1F *massgen = (TH1F*)result->Get("massa bosoni generati bene");
  TH1F *massric = (TH1F*)result->Get("massa bosoni ricostruiti bene");
  massric->SetTitle("efficienza in funzione della massa");
  TH1F *massefficiency = (TH1F*)massric->Clone("massefficiency");
  massefficiency->Divide(massgen);
  massefficiency->GetYaxis()->SetRangeUser(0,1);
  massefficiency->GetXaxis()->SetTitle("massa (GeV/c^2)");
  massefficiency->GetYaxis()->SetTitle("efficienza");
  massefficiency->SetLineColor(1);
  massefficiency->SetLineWidth(2);
  massefficiency->Draw();
  c5->cd();
  massgen->SetTitle("massa bosoni generati e ricostruiti");
  massgen->GetYaxis()->SetRangeUser(0,1000);
  massgen->GetXaxis()->SetTitle("massa (GeV/c^2)");
  massgen->SetLineColor(1);
  massric->SetLineColor(2);
  massgen->SetLineWidth(2);
  massric->SetLineWidth(2);
  massgen->Draw();
  massric->Draw("same");

  c2->cd();
  TH1F *ptgen = (TH1F*)result->Get("pt bosoni generati bene");
  TH1F *ptric = (TH1F*)result->Get("pt bosoni ricostruiti bene");
  ptric->SetTitle("efficienza in funzione della pt");
  TH1F *ptefficiency = (TH1F*)ptric->Clone("ptefficiency");
  ptefficiency->Divide(ptgen);
  ptefficiency->GetYaxis()->SetRangeUser(0,1);
  ptefficiency->GetXaxis()->SetTitle("pt (GeV/c)");
  ptefficiency->GetYaxis()->SetTitle("efficienza");
  ptefficiency->SetLineColor(1);
  ptefficiency->SetLineWidth(2);
  ptefficiency->Draw();
  c6->cd();
  ptgen->SetTitle("pt bosoni generati e ricostruiti");
  ptgen->GetYaxis()->SetRangeUser(0,600);
  ptgen->GetXaxis()->SetTitle("pt (GeV/c)");
  ptgen->SetLineColor(1);
  ptric->SetLineColor(2);
  ptgen->SetLineWidth(2);
  ptric->SetLineWidth(2);
  ptgen->Draw();
  ptric->Draw("same");

  c3->cd();
  TH1F *energygen = (TH1F*)result->Get("energia bosoni generati bene");
  TH1F *energyric = (TH1F*)result->Get("energia bosoni ricostruiti bene");
  energyric->SetTitle("efficienza in funzione dell'energia");
  TH1F *energyefficiency = (TH1F*)energyric->Clone("energyefficiency");
  energyefficiency->Divide(energygen);
  energyefficiency->GetYaxis()->SetRangeUser(0,1);
  energyefficiency->GetXaxis()->SetTitle("E (GeV)");
  energyefficiency->GetYaxis()->SetTitle("efficienza");
  energyefficiency->SetLineColor(1);
  energyefficiency->SetLineWidth(2);
  energyefficiency->Draw();
  c7->cd();
  energygen->SetTitle("energia bosoni generati e ricostruiti");
  energygen->GetYaxis()->SetRangeUser(0,1200);
  energygen->GetXaxis()->SetTitle("energia(GeV)");
  energygen->SetLineColor(1);
  energyric->SetLineColor(2);
  energygen->SetLineWidth(2);
  energyric->SetLineWidth(2);
  energygen->Draw();
  energyric->Draw("same");

  c4->cd(); 
  TH1F *etagen = (TH1F*)result->Get("eta bosoni generati bene");
  TH1F *etaric = (TH1F*)result->Get("eta bosoni ricostruiti bene");
  etaric->SetTitle("efficienza in funzione della eta");
  TH1F *etaefficiency = (TH1F*)etaric->Clone("etaefficiency");
  etaefficiency->Divide(etagen);
  etaefficiency->GetYaxis()->SetRangeUser(0,1);
  etaefficiency->GetXaxis()->SetTitle("eta");
  etaefficiency->GetYaxis()->SetTitle("efficienza");
  etaefficiency->SetLineColor(1);
  etaefficiency->SetLineWidth(2);
  etaefficiency->Draw();
  c8->cd();
  etagen->SetTitle("eta bosoni generati e ricostruiti");
  etagen->GetYaxis()->SetRangeUser(0,200);
  etagen->GetXaxis()->SetTitle("eta");
  etagen->SetLineColor(1);
  etaric->SetLineColor(2);
  etagen->SetLineWidth(2);
  etaric->SetLineWidth(2);
  etagen->Draw();
  etaric->Draw("same");
}
