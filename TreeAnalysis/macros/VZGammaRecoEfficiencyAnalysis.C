/**
 *  Macro to get efficiency of Bosons reconstruction from results of VZGAnalyzer 
 *
 *
 *	Usage: 	root [-l] [-b] [-q] 'VZGammaRecoEfficiencyAnalysis.C("<sample name>")'		-l do not display banner 	-b run in background		-q close after finishing
 *	e.g.: 	root -l 'VZGammaRecoEfficiencyAnalysis.C("WZGTo2L2jG")'
 *	It may become necessary to change the path and/or to expand the samples list	
 *			
 *  $Date: 2023/06/26 10:17:08 
 *  $Revision: 2.0 $
 *
 *  \author C. Tarricone cristiano.tarricone@cern.ch
 */

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TPaveText.h>
#include <map>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define TEST_MODE

using namespace std;

void VZGammaRecoEfficiencyAnalysis(TString requestedSample) {
  
  TString path = "/eos/home-c/ctarrico/Frameworks/CMSSW_10_6_26/src/VVXAnalysis/TreeAnalysis/results/2017/VZGAnalyzer_MC/";

  vector<TString> samples = {"WZGTo2L2jG", "ZZGTo2L2jG"};
  vector<TString> parNames = {"W","Z"};
  vector<TString> typeNames = {"Pt"};
  vector<TString> algNames = {"m6W4Z", "m7W3Z", "m8W2Z", "m9W1Z"};//{"mW", "mZ", "mWZ", "m8W2Z"};

  map<TString, TString> algConv;
  algConv[TString("mW")] = TString("V cand w/mass closest to m_{W}");
  algConv[TString("mZ")] = TString("V cand w/mass closest to m_{Z}");
  algConv[TString("mWZ")] = TString("V cand w/mass closest to the closest m_{V}");
  algConv[TString("m8W2Z")] = TString("V cand w/mass closer 0.8*m_{W}+0.2*m_{Z}");
  algConv[TString("m6W4Z")] = TString("V cand w/mass closer 0.6*m_{W}+0.4*m_{Z}");
  algConv[TString("m7W3Z")] = TString("V cand w/mass closer 0.7*m_{W}+0.3*m_{Z}");
  algConv[TString("m9W1Z")] = TString("V cand w/mass closer 0.9*m_{W}+0.1*m_{Z}");


  TString sampleName, name;
  if(requestedSample != "") {
    for(int i = 0; i < samples.size(); i++) {
      if(requestedSample == samples.at(i)) {
        sampleName = requestedSample;
        name=parNames.at(i);
        cout<<name<<"\n";
      }
    }
  }

  if(sampleName == "") {
    cout<<"Unknown sample \""<<sampleName<<"\"\n";
    return;
  }

  cout<<"Opening \""<<sampleName<<".root\"\n";
  TFile* result = TFile::Open(path + sampleName + ".root");

  TCanvas *cDrawing = new TCanvas(name + "ReconstructionEfficiency_vs_" + typeNames.at(0), name + "ReconstructionEfficiency_vs_" + typeNames.at(0), 10, 0, 1280, 1024);

  TLegend* legend = new TLegend(0.65, 0.75, 0.9, 0.9);
  legend->SetBorderSize(1);

  int nGraphs = 0;

  // Array di colori per i grafici dei rapporti
  Color_t colors[] = {kRed, kBlue, kGreen, kOrange, kMagenta};

  foreach(TString& alg, algNames) {
    foreach(TString& type, typeNames) {
      TH1F* hNum = (TH1F*)result->Get(type+"_"+alg+"_num");
      if(hNum == nullptr) {
        cout<<"Could not open gen"<<type<<"_"<<alg<<"_num""\"\n";
        continue;
      }

      TH1F* hDen = (TH1F*)result->Get(type+"_den");
      if(hDen == nullptr) {
        cout<<"Could not open gen"<<name<<type<<"_"<<sampleName<<"_"<<alg<<"_den""\"\n";
        continue;
      }

      TGraphAsymmErrors* hEff = new TGraphAsymmErrors(hNum, hDen, "cp");
      hEff->SetTitle("VB Reconstruction Efficiency");
      hEff->GetYaxis()->SetRangeUser(0.4, 1.01);

      hEff->SetLineColor(colors[nGraphs]); // Imposta il colore del grafico del rapporto

      hEff->Draw(nGraphs == 0 ? "AP" : "P SAME");

      hEff->GetXaxis()->SetTitle("p_{T}^{VCand} (GeV/c)");
      hEff->GetXaxis()->SetTitleSize(0.035);
      hEff->GetYaxis()->SetTitle("VB Reconstruction Efficiency");
      hEff->GetYaxis()->SetTitleSize(0.035);

      TString legendEntry = algConv[alg];
      legend->AddEntry(hEff, legendEntry, "lp");
      nGraphs++;
    }
  }

  legend->Draw();

  //cDrawing->SetLogy();

  cDrawing->Modified();
  cDrawing->Update();

  TPaveText *paveTLeft = new TPaveText(0.065, 0.87, 0.9, 0.95, "NDCNDC");
  paveTLeft->SetFillColor(0);
  paveTLeft->SetFillStyle(0);
  paveTLeft->SetBorderSize(0);
  paveTLeft->SetTextAlign(11);
  paveTLeft->SetTextFont(62);
  paveTLeft->SetTextSize(0.02);
  paveTLeft->AddText("CMS Preliminary");

  TPaveText *paveTRight = new TPaveText(0.80, 0.87, 0.9, 0.95, "NDCNDC");
  paveTRight->SetFillColor(0);
  paveTRight->SetFillStyle(0);
  paveTRight->SetBorderSize(0);
  paveTRight->SetTextAlign(11);
  paveTRight->SetTextFont(62);
  paveTRight->SetTextSize(0.02);
  paveTRight->AddText("Simulation");

  cDrawing->cd();
  paveTLeft->Draw();
  paveTLeft->SetBit(kCanDelete);
  paveTRight->Draw();
  paveTRight->SetBit(kCanDelete);



  TString fname = "Refinement_" + name + "RecoEff_" + typeNames.at(0) + ".png";
  cDrawing->SaveAs("efficiencyPlots/" + fname);

  result->Close("R");
}
