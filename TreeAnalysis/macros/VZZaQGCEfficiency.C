#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define TEST_MODE

using namespace std;

void efficiency(string nomegrafico1,string nomegrafico2,string titoloasse, string titolo, string titoloeff, double ymax, TFile *result, TCanvas *c){
  TVirtualPad *pd1=c->cd(1);
  TH1F *grafico1= (TH1F*)result->Get(nomegrafico1.c_str());
  TH1F *grafico2= (TH1F*)result->Get(nomegrafico2.c_str());
  TH1F *grafico3= (TH1F*)grafico2->Clone("efficienza");
  
  grafico1->SetTitle(titolo.c_str());
  grafico1->GetYaxis()->SetRangeUser(0,ymax);
  grafico1->GetXaxis()->SetTitle(titoloasse.c_str());
  grafico1->SetLineColor(1);
  grafico2->SetLineColor(2);
  grafico1->SetLineWidth(2);
  grafico2->SetLineWidth(2);
  grafico1->Draw();
  grafico2->Draw("same");
  TLegend *legend= new TLegend(0.7,0.77,0.98,0.94,"");
  legend->AddEntry(grafico1,"bosoni generati","l");
  legend->AddEntry(grafico2,"bosoni ricostruiti","l");
  legend->Draw();

  TVirtualPad *pd2=c->cd(2);
  grafico3->Divide(grafico1);
  grafico3->SetTitle(titoloeff.c_str());
  grafico3->GetYaxis()->SetRangeUser(0,1);
  grafico3->GetXaxis()->SetTitle(titoloasse.c_str());
  grafico3->GetYaxis()->SetTitle("efficienza");
  grafico3->SetLineColor(1);
  grafico3->SetLineWidth(2);
  grafico3->Draw();}


void VZZaQGCEfficiency(){

  gErrorIgnoreLevel = kFatal;

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

  efficiency("massa Z1 generati bene","massa Z1 ricostruiti bene","massa (Gev/c^2)","massa Z1 generati e ricostruiti","efficienza in funzione della massa di Z1",1000,result,c1);
  
  efficiency("massa Z2 generati bene","massa Z2 ricostruiti bene","massa (Gev/c^2)","massa Z2 generati e ricostruiti","efficienza in funzione della massa di Z2",1000,result,c2);
  
  efficiency("pt Z1 generati bene","pt Z1 ricostruiti bene","pt (Gev/c)","pt Z1 generati e ricostruiti","efficienza in funzione della pt di Z1",600,result,c3);
  
  efficiency("pt Z2 generati bene","pt Z2 ricostruiti bene","pt (Gev/c)","pt Z2 generati e ricostruiti","efficienza in funzione della pt di Z2",600,result,c4);
  
  efficiency("energia Z1 generati bene","energia Z1 ricostruiti bene","E (Gev)","energia Z1 generati e ricostruiti","efficienza in funzione dell'energia di Z1",1200,result,c5);
  
  efficiency("energia Z2 generati bene","energia Z2 ricostruiti bene","E (Gev)","energia Z2 generati e ricostruiti","efficienza in funzione dell'energia di Z2",1200,result,c6);
  
  efficiency("eta Z1 generati bene","eta Z1 ricostruiti bene","eta","eta Z1 generati e ricostruiti","efficienza in funzione dell'eta di Z1",200,result,c7);
  
  efficiency("eta Z2 generati bene","eta Z2 ricostruiti bene","eta","eta Z2 generati e ricostruiti","efficienza in funzione dell'eta di Z2",200,result,c8);
  
  efficiency("E leptone maggiore Z1 buono","E leptone maggiore Z1 buono ricostruito","E (Gev)","energia leptoni maggiori di Z1 generati e ricostruiti","efficienza in funzione dell'energia leptone maggiore di Z1",400,result,c9);
  
  efficiency("E leptone maggiore Z2 buono","E leptone maggiore Z2 buono ricostruito","E (Gev)","energia leptoni maggiori di Z2 generati e ricostruiti","efficienza in funzione dell'energia leptone maggiore di Z2",400,result,c10);
  
  efficiency("E leptone minore Z1 buono","E leptone minore Z1 buono ricostruito","E (Gev)","energia leptoni minori di Z1 generati e ricostruiti","efficienza in funzione dell'energia leptone minore di Z1",400,result,c11);
  
  efficiency("E leptone minore Z2 buono","E leptone minore Z2 buono ricostruito","E (Gev)","energia leptoni minori di Z2 generati e ricostruiti","efficienza in funzione dell'energia leptone minore di Z2",400,result,c12);
}
