#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <sstream>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace std;

void DistributionComparison(string nomegrafico, string titolo,string nomeasse, double ymax, TFile *result1, TFile *result2, double scale){
  TH1F *grafico1= (TH1F*)result1->Get(nomegrafico.c_str());
  TH1F *grafico2= (TH1F*)result2->Get(nomegrafico.c_str());
  grafico1->SetTitle(titolo.c_str());
  grafico1->GetYaxis()->SetRangeUser(0,ymax);
  grafico1->GetXaxis()->SetTitle(nomeasse.c_str());
  grafico1->SetLineColor(1);
  grafico2->SetLineColor(2);
  grafico1->SetLineWidth(2);
  grafico2->SetLineWidth(2);
  grafico2->Scale(scale);
  grafico1->Draw();
  grafico2->Draw("same");
  TLegend *legend= new TLegend(0.75,0.75,0.98,0.95);
  legend->AddEntry(grafico1,"dati originali","l");
  legend->AddEntry(grafico2,"dati MadGraph","l");
  ostringstream stream;
  stream << grafico1->Chi2Test(grafico2,"WW");
  string insert= "p value: "+stream.str();
  legend->AddEntry((TObject*)0,insert.c_str(),"");
  legend->Draw();
}
  

void VZZDistributionComparison(){

  gErrorIgnoreLevel = kWarning; //per ignorare messaggi "Info"

  TCanvas *c1 = new TCanvas("c1","canvas",0,0,1000,1000);
  TCanvas *c2 = new TCanvas("c2","canvas",0,0,1000,1000);
  TCanvas *c3 = new TCanvas("c3","canvas",0,0,1000,1000);
  TCanvas *c4 = new TCanvas("c4","canvas",0,0,1000,1000);
  TCanvas *c5 = new TCanvas("c5","canvas",0,0,1000,1000);
  TCanvas *c6 = new TCanvas("c6","canvas",0,0,1000,1000);
  TCanvas *c7 = new TCanvas("c7","canvas",0,0,1000,1000);
  TCanvas *c8 = new TCanvas("c8","canvas",0,0,1000,1000);
  TCanvas *c9 = new TCanvas("c9","canvas",0,0,1000,1000);
  TCanvas *c10 = new TCanvas("c10","canvas",0,0,1000,1000);
  TCanvas *c11 = new TCanvas("c11","canvas",0,0,1000,1000);

  TFile *result1 = TFile::Open("./results/2018/VVXnocutsAnalyzer_MC/WZZ.root");
  TFile *result2 = TFile::Open("./results/2016/VVXnocutsAnalyzer_MC/test.root");

  double entries1=6162.0;
  double entries2=7780.0;
  double scale=entries1/entries2;

  c1->cd();
  DistributionComparison("massa tribosoni","massa tribosoni","massa (GeV/c^2)",700,result1,result2,scale);
  c2->cd();
  DistributionComparison("massa bosoni a coppie","massa bosoni a coppie","massa (GeV/c^2)",2000,result1,result2,scale);
  c3->cd();
  DistributionComparison("energia totale bosoni","energia totale","energia (GeV)",700,result1,result2,scale);
  c4->cd();
  DistributionComparison("energia bosone piu' energetico","energia bosone piu' energetico","energia (GeV)",400,result1,result2,scale);
  c5->cd();
  DistributionComparison("pt bosone maggiore","pt massima","pt (GeV/c)",600,result1,result2,scale);
  c6->cd();
  DistributionComparison("angolo relativo bosoni","angolo relativo bosoni","angolo (rad)",900,result1,result2,scale);

  double entries3=8575.0;
  double entries4=1350.0;
  double scale2=entries3/entries4;
  
  c7->cd();
  DistributionComparison("E leptone maggiore","energia leptone piu' energetico","energia (GeV)",2000,result1,result2,scale2);
  c8->cd();
  DistributionComparison("E leptone minore","energia leptone meno energetico","energia (GeV)",2800,result1,result2,scale2);
  c9->cd();
  DistributionComparison("massa bosoni leptonici","massa bosoni leptonici","massa (GeV/c^2)",3300,result1,result2,scale2);
  c10->cd();
  DistributionComparison("energia bosoni leptonici","energia bosoni leptonici","energia (GeV)",1700,result1,result2,scale2);
  c11->cd();
  DistributionComparison("pt bosoni leptonici","pt bosoni leptonici","pt (GeV/c)",1500,result1,result2,scale2);
}
