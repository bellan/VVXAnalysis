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

TString resDir="resultsForReweighting";
//TString year="2018";
const std::vector<TString> years = {"2016preVFP","2016postVFP","2017","2018"};


using namespace std;

void extractWeights() {

  gROOT->SetBatch(kTRUE);
  TCanvas *cDrawing = new TCanvas("", "", 0, 0, 800, 800);

  TLegend* legend = new TLegend(0.15, 0.15, 0.45, 0.45);
  legend->SetBorderSize(1);


  // Array di colori per i grafici dei rapporti
  Color_t colors[] = {kRed, kSpring, kBlue, kViolet  };
  Int_t markerStyles[] = {20,22,23,24};

  for(int i=0; i<4; i++){
    TString year=years[i];
    cout<<"---------------"<<year<<"---------------"<<endl;
    TString path = "../"+resDir+"/"+year+"/VZGAnalyzer_SR2P/";

    TFile *DYFile = new TFile(path+"DYJetsToLL_M50.root", "READ");
    TFile *DataFile = new TFile(path+"data_obs.root", "READ");

    if (!DYFile || DYFile->IsZombie()) {
      cout << "Error opening " << path+"DYJetsToLL_M50.root" << endl;
      return;
    }
    if (!DataFile || DataFile->IsZombie()) {
      cout << "Error opening " << path+"data_obs.root" << endl;
      return;
    }

    TH1F* hNum = (TH1F*)DataFile->Get("REWGT_numDATA");
    TH1F* hDen = (TH1F*)DYFile->Get("REWGT_DY-GEN");
    if (!hNum || !hDen) {
      cout << "Error: Histogram not found!" << endl;
      return;
    }
    /*
    for (int j = 1; j <= hDen->GetNbinsX(); j++) {
      cout << "Bin " << j << "= " <<hNum->GetBinContent(j)<<" / "<<hDen->GetBinContent(j)<<" = "<<hNum->GetBinContent(j)/hDen->GetBinContent(j)<< endl;

      if (hDen->GetBinContent(j) <= 0) {
        cout << "Warning: Bin " << j << " in hDen is zero!" << endl;
        hDen->SetBinContent(j, 1e-6); // Evita divisioni per zero
      }
    }
    */
    cout <<"{ ";
    for (int j = 1; j <= hDen->GetNbinsX(); j++) {
      cout << hNum->GetBinContent(j)/hDen->GetBinContent(j)<<", ";

      if (hDen->GetBinContent(j) <= 0) {
        cout << "Warning: Bin " << j << " in hDen is zero!" << endl;
        hDen->SetBinContent(j, 1e-6); // Evita divisioni per zero
      }
    }
    cout <<"} "<<endl;
    cout<< "denominator n Bins = "<<hDen->GetNbinsX()<<endl;
    cout<< "numerator   n Bins = "<<hNum->GetNbinsX()<<endl;

    hNum->Sumw2();
    hDen->Sumw2();

    TH1F *hEff = (TH1F*)hNum->Clone("hEff");
    hEff->Divide(hDen);

    hEff->SetTitle("");
    hEff->GetYaxis()->SetRangeUser(0., 3);

    hEff->SetLineColor(colors[i]);
    hEff->SetMarkerStyle(markerStyles[i]);
    hEff->SetMarkerColor(colors[i]);

    hEff->Draw(i == 0 ? "AP" : "P SAME");

    hEff->GetXaxis()->SetTitle("gen. Z p_{T} [GeV]");
    hEff->GetXaxis()->SetTitleSize(0.035);
    hEff->GetYaxis()->SetTitle("weights");
    hEff->GetYaxis()->SetTitleSize(0.035);

    legend->AddEntry(hEff, year+"weights", "lp");
    
    DYFile->Close();
    DataFile->Close();

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
  paveTLeft->AddText("CMS Private Work");

  TPaveText *paveTRight = new TPaveText(0.80, 0.87, 0.9, 0.95, "NDCNDC");
  paveTRight->SetFillColor(0);
  paveTRight->SetFillStyle(0);
  paveTRight->SetBorderSize(0);
  paveTRight->SetTextAlign(11);
  paveTRight->SetTextFont(62);
  paveTRight->SetTextSize(0.02);
  paveTRight->AddText("Simulation (13 TeV)");

  cDrawing->cd();
  paveTLeft->Draw();
  paveTLeft->SetBit(kCanDelete);
  paveTRight->Draw();
  paveTRight->SetBit(kCanDelete);



  TString fname = "DY_Zpt_reweight";
  cDrawing->SaveAs("efficiencyPlots/" + fname +".png");
  cDrawing->SaveAs("efficiencyPlots/" + fname +".pdf");

}
