/**
 *  Macro to get efficiency from results of WZZAnalyzer 
 *		
 *  $Date: 2019/09/12 9:47:08 
 *  $Revision: 0.0 $
 *
 *  \author C. Tarricone cristiano.tarrico@edu.unito.it
 */

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define TEST_MODE

using namespace std;

void WZZEfficiencyAnalysis(){
  //cout<<"Opening \""<<sampleName<<".root\"\n";
  TFile* result = TFile::Open("~/VVXAnalysis/TreeAnalysis/results/WZZAnalyzer_MC/WZZ.root");
  TH1F* hDen = (TH1F*)result->Get("genMuonEta_den");
  TH1F* hNum = (TH1F*)result->Get("genMuonEta_num");
  //TGraphAsymmErrors *hEff = new TGraphAsymmErrors(hNum, hDen, "cp")
  TGraphAsymmErrors* hEff = new TGraphAsymmErrors(hNum, hDen, "cp");
  hEff->SetTitle("Electrons_efficiency_vs_eta");
  hEff->GetYaxis()->SetRangeUser(0.,1.01);
  TCanvas *cDrawing = new TCanvas("Muons_efficiency_vs_eta","Muons_efficiency_vs_eta", 10,0,1280,1024);
  cDrawing->cd();
  hEff->Draw("AP");
  result->Close("R");

}
