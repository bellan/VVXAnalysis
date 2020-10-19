#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define TEST_MODE

using namespace std;

void VZZDistributionComparison(){

  TCanvas *c1 = new TCanvas("c1","canvas",0,0,1000,1000);
  TCanvas *c2 = new TCanvas("c2","canvas",0,0,1000,1000);
  TCanvas *c3 = new TCanvas("c3","canvas",0,0,1000,1000);
  TCanvas *c4 = new TCanvas("c4","canvas",0,0,1000,1000);
  TCanvas *c5 = new TCanvas("c5","canvas",0,0,1000,1000);
  TCanvas *c6 = new TCanvas("c6","canvas",0,0,1000,1000);

  TFile *result1 = TFile::Open("./results/2018/VVXnocutsAnalyzer_MC/WZZ.root");
  TH1F *masstot1 = (TH1F*)result1->Get("massa tribosoni");
  TH1F *massacoppie1 = (TH1F*)result1->Get("massa bosoni a coppie");
  TH1F *energiatot1 = (TH1F*)result1->Get("energia totale bosoni");
  TH1F *energiamax1 = (TH1F*)result1->Get("energia bosone piu' energetico");
  TH1F *ptmax1 = (TH1F*)result1->Get("pt bosone maggiore");
  TH1F *angolo1 = (TH1F*)result1->Get("angolo relativo bosoni");

  TFile *result2 = TFile::Open("./results/2016/VVXnocutsAnalyzer_MC/test.root");
  TH1F *masstot2 = (TH1F*)result2->Get("massa tribosoni");
  TH1F *massacoppie2 = (TH1F*)result2->Get("massa bosoni a coppie");
  TH1F *energiatot2 = (TH1F*)result2->Get("energia totale bosoni");
  TH1F *energiamax2 = (TH1F*)result2->Get("energia bosone piu' energetico");
  TH1F *ptmax2 = (TH1F*)result2->Get("pt bosone maggiore");
  TH1F *angolo2 = (TH1F*)result2->Get("angolo relativo bosoni");
  
  c1->cd();
  masstot1->SetTitle("massa tribosoni");
  masstot1->GetYaxis()->SetRangeUser(0,700);
  masstot1->GetXaxis()->SetTitle("massa (GeV/c^2)");
  masstot1->SetLineColor(1);
  masstot2->SetLineColor(2);
  masstot1->SetLineWidth(2);
  masstot2->SetLineWidth(2);
  cout<<masstot1->Integral()<<endl;
  cout<<masstot2->Integral()<<endl;
  masstot2->Scale(8);
  masstot1->Draw();
  masstot2->Draw("same");

  c2->cd();
  massacoppie1->SetTitle("massa bosoni a coppie");
  massacoppie1->GetYaxis()->SetRangeUser(0,2000);
  massacoppie1->GetXaxis()->SetTitle("massa (GeV/c^2)");
  massacoppie1->SetLineColor(1);
  massacoppie2->SetLineColor(2);
  massacoppie1->SetLineWidth(2);
  massacoppie2->SetLineWidth(2);
  massacoppie2->Scale(8);
  massacoppie1->Draw();
  massacoppie2->Draw("same");

  c3->cd();
  energiatot1->SetTitle("energia totale bosoni");
  energiatot1->GetYaxis()->SetRangeUser(0,700);
  energiatot1->GetXaxis()->SetTitle("energia (GeV)");
  energiatot1->SetLineColor(1);
  energiatot2->SetLineColor(2);
  energiatot1->SetLineWidth(2);
  energiatot2->SetLineWidth(2);
  energiatot2->Scale(8);
  energiatot1->Draw();
  energiatot2->Draw("same");

  c4->cd();
  energiamax1->SetTitle("energia bosone piu' energetico");
  energiamax1->GetYaxis()->SetRangeUser(0,400);
  energiamax1->GetXaxis()->SetTitle("energia (GeV)");
  energiamax1->SetLineColor(1);
  energiamax2->SetLineColor(2);
  energiamax1->SetLineWidth(2);
  energiamax2->SetLineWidth(2);
  energiamax2->Scale(8);
  energiamax1->Draw();
  energiamax2->Draw("same");

  c5->cd();
  ptmax1->SetTitle("pt massima bosoni");
  ptmax1->GetYaxis()->SetRangeUser(0,600);
  ptmax1->GetXaxis()->SetTitle("pt (GeV/c)");
  ptmax1->SetLineColor(1);
  ptmax2->SetLineColor(2);
  ptmax1->SetLineWidth(2);
  ptmax2->SetLineWidth(2);
  ptmax2->Scale(8);
  ptmax1->Draw();
  ptmax2->Draw("same");

  c6->cd();
  angolo1->SetTitle("angolo relativo bosoni");
  angolo1->GetYaxis()->SetRangeUser(0,900);
  angolo1->GetXaxis()->SetTitle("angolo (rad)");
  angolo1->SetLineColor(1);
  angolo2->SetLineColor(2);
  angolo1->SetLineWidth(2);
  angolo2->SetLineWidth(2);
  angolo2->Scale(8);
  angolo1->Draw();
  angolo2->Draw("same"); 
}
