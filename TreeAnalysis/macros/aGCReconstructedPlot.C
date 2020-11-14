/*
Macro to plot the expected results for the presence of anomalous couplings in reconstructed VZZ variables. Path will probably need changing
Author:Marozzo Giovanni Battista
Date: 2020/11/25
*/

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLatex.h>

using namespace std;

const double sigma=0.05565;
const double sampleevents=58308;

void singlevariablePlot(string histname,string title, string axistitle, TFile *result){
  
  TH1F *histSM= (TH1F*)result->Get(histname.c_str());
  string histname2="weighted "+histname;
  TH1F *histaGC= (TH1F*)result->Get(histname2.c_str());
  
  histSM->SetTitle(title.c_str());
  histSM->GetYaxis()->SetRangeUser(0,10);
  histSM->GetXaxis()->SetTitle(axistitle.c_str());
  histSM->SetLineColor(1);
  histaGC->SetLineColor(2);
  histSM->SetLineWidth(2);
  histaGC->SetLineWidth(2);
  histSM->Draw();
  histaGC->Draw("same");
  
  TLegend *legend= new TLegend(0.7,0.77,0.98,0.94,"");
  ostringstream streamSM,streamaGC;
  streamSM << 1000*sigma*histSM->Integral()/sampleevents;
  string legendentrySM="SM data; #sigma= "+streamSM.str()+" fb";
  legend->AddEntry(histSM,legendentrySM.c_str(),"l");
  streamaGC << 1000*sigma*histaGC->Integral()/sampleevents;
  string legendentryaGC="aGC data; #sigma= "+streamaGC.str()+" fb";
  legend->AddEntry(histaGC,legendentryaGC.c_str(),"l");
  legend->Draw();

  ostringstream stream1,stream2,stream3;
  stream1 << histSM->Chi2Test(histaGC,"WWCHI2");
  string insert1= "chi2: "+stream1.str();
  legend->AddEntry((TObject*)0,insert1.c_str(),"");
  stream2 << histSM->Chi2Test(histaGC,"WWCHI2/NDF");
  string insert2= "chi2/NDF: "+stream2.str();
  legend->AddEntry((TObject*)0,insert2.c_str(),"");
  legend->Draw();
  stream3 << histSM->Chi2Test(histaGC,"WW");
  string insert3= "p value: "+stream3.str();
  legend->AddEntry((TObject*)0,insert3.c_str(),"");
  legend->Draw();
}

void aGCReconstructedPlot(string sample){

  gErrorIgnoreLevel = kError; //to ignore "Info" and "Warning" messages
  
  TCanvas *c1 = new TCanvas("c1","canvas",0,0,1000,500);
  /*
  TCanvas *c2 = new TCanvas("c2","canvas",0,0,1000,500);
  TCanvas *c3 = new TCanvas("c3","canvas",0,0,1000,500);
  */
  
  string filename="./results/2018/VZZaQGCAnalyzer_MC/"+sample+".root";
  TFile *result=TFile::Open(filename.c_str());
  
  c1->cd();
  singlevariablePlot("energy of all bosons 2","Energy of all bosons","E (GeV)",result);
  /*
  c2->cd();
  singlevariablePlot("mass of tribosons","Mass of tribosons","m (GeV/c^2)",result);
  c3->cd();
  singlevariablePlot("energy of major Z1 leptons","Energy of major Z1 leptons","E (GeV)",result);
  */
}
  
