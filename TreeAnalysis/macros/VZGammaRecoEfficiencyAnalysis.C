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
  
  TString path = "/eos/home-c/ctarrico/Frameworks/CMSSW_10_6_26/src/VVXAnalysis/TreeAnalysis/results/2016preVFP/VZGAnalyzer_SR2P/";

  vector<TString> samples = {"WZGTo2L2jG", "ZZGTo2L2jG"};
  vector<TString> parNames = {"W","Z"};
  vector<TString> typeNames = {"Pt"}; 
  vector<TString> algNames = {"mjjBased", "mixSmth01", "mixSmth05", "mixSmth10", "mixSmth20", "mixSmth100", "qglBased"};//{"mW", "mZ", "mWZ", "m8W2Z"};

  map<TString, TString> algConv;
  algConv[TString("mjjBased")]   = TString("DiJet mass based");
  algConv[TString("qglBased")]   = TString("DiJet QGL based");
  algConv[TString("mixSmth01")]  = TString("Weight transition par. = 1/10");
  algConv[TString("mixSmth05")]  = TString("Weight transition par. = 1/2");
  algConv[TString("mixSmth10")]  = TString("Weight transition par. = 1");
  algConv[TString("mixSmth20")]  = TString("Weight transition par. = 2");
  algConv[TString("mixSmth100")] = TString("Weight transition par. = 10");


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

  TCanvas *cDrawing = new TCanvas(name + "ReconstructionEfficiency_vs_" + typeNames.at(0), name + "ReconstructionEfficiency_vs_" + typeNames.at(0), 0, 0, 800, 800);

  TLegend* legend = new TLegend(0.15, 0.15, 0.45, 0.45);
  legend->SetBorderSize(1);

  int nGraphs = 0;

  // Array di colori per i grafici dei rapporti
  Color_t colors[] = {kRed, kOrange+2, kOrange-3, kAzure, kSpring, kBlue, kViolet  };
  Int_t markerStyles[] = {20,22,23,24,26,27,21  };

  foreach(TString& alg, algNames) {
    foreach(TString& type, typeNames) {
      TH1F* hNum = (TH1F*)result->Get(type+"_"+alg+"_num");
      //      hNum->Rebin(2);
      if(hNum == nullptr) {
        cout<<"Could not open gen"<<type<<"_"<<alg<<"_num""\"\n";
        continue;
      }

      TH1F* hDen = (TH1F*)result->Get(type+"_den");
      //      hDen->Rebin(2);
      if(hDen == nullptr) {
        cout<<"Could not open gen"<<name<<type<<"_"<<sampleName<<"_"<<alg<<"_den""\"\n";
        continue;
      }

      TGraphAsymmErrors* hEff = new TGraphAsymmErrors(hNum, hDen, "n");
      hEff->SetTitle(name+" Reconstruction Efficiency");
      hEff->GetYaxis()->SetRangeUser(0., 1.01);

      hEff->SetLineColor(colors[nGraphs]);
      hEff->SetMarkerStyle(markerStyles[nGraphs]);
      hEff->SetMarkerColor(colors[nGraphs]);

      hEff->Draw(nGraphs == 0 ? "AP" : "P SAME");

      hEff->GetXaxis()->SetTitle("p_{T} V had. Cand.");
      hEff->GetXaxis()->SetTitleSize(0.035);
      hEff->GetYaxis()->SetTitle(name+" Reconstruction Efficiency");
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
