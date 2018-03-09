#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMath.h>

using namespace std;

// --------------------- 1D graphs --------------------- 
void Draw1D(TString name, TString back_path = "QCD", TString sig_path = "EW"){
  // analysis results files
  TString path = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/WZ";
  TFile *EW = TFile::Open(path + sig_path + ".root");
  TFile *QCD = TFile::Open(path + back_path + ".root");

  // canvas
  TCanvas *cDrawing = new TCanvas("cDrawing", name, 10, 0, 1280, 1024);
  gStyle->SetOptStat(0);
  cDrawing->SetGrid();

  // signal histogram
  TH1F *hEW = (TH1F*)EW->Get(name);
  hEW->SetMarkerColor(kGreen+2);
  hEW->SetMarkerSize(0.6);
  hEW->SetMarkerStyle(21);
  hEW->SetLineColor(kGreen+2);

  // backgroung histogram
  TH1F *hQCD = (TH1F*)QCD->Get(name);
  hQCD->SetMarkerColor(kRed);
  hQCD->SetMarkerSize(0.6);
  hQCD->SetMarkerStyle(21);
  hQCD->SetLineColor(kRed);

  // legend
  TLegend *legend = new TLegend(0.75, 0.80, 0.95, 0.95);
  legend->AddEntry(hEW, "signal");
  legend->AddEntry(hQCD, "background");

  cout << "\nEW's " << name << " integral is  " << hEW->Integral(0,-1) << endl;
  cout << "QCD's " << name << " integral is " << hQCD->Integral(0,-1) << endl << endl;

  // drawing
  cDrawing->cd();
  hEW-> SetMaximum(1.1 * max(hEW->GetMaximum(), hQCD->GetMaximum()));
  hEW->Draw();
  hQCD->Draw("same");
  legend->Draw();
}

void Draw1DNorm(TString name, TString back_path = "QCD", TString sig_path = "EW"){
  // analysis results files
  TString path = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/WZ";
  TFile *EW = TFile::Open(path + sig_path + ".root");
  TFile *QCD = TFile::Open(path + back_path + ".root");

  // canvas
  TCanvas *cDrawing = new TCanvas("cDrawing", name, 10, 0, 1280, 1024);
  gStyle->SetOptStat(0);
  cDrawing->SetGrid();

  // signal histogram
  TH1F *hEW = (TH1F*)EW->Get(name);
  hEW->Scale(1/hEW->Integral(0, -1));
  hEW->SetMarkerSize(0.6);
  hEW->SetMarkerStyle(21);
  hEW->SetMarkerColor(kGreen+2);
  hEW->SetLineColor(kGreen+2);

  // backgroung histogram
  TH1F *hQCD = (TH1F*)QCD->Get(name);
  hQCD->Scale(1/hQCD->Integral(0, -1));
  hQCD->SetMarkerSize(0.6);
  hQCD->SetMarkerStyle(21);
  hQCD->SetMarkerColor(kRed);
  hQCD->SetLineColor(kRed);

  // legend
  TLegend *legend = new TLegend(0.75, 0.80, 0.95, 0.95);
  legend->AddEntry(hEW, "signal");
  legend->AddEntry(hQCD, "background");

  // drawing
  cDrawing->cd();
  hEW-> SetMaximum(1.1 * max(hEW->GetMaximum(), hQCD->GetMaximum()));
  hEW->Draw();
  hQCD->Draw("same");
  legend->Draw();
}

// --------------------- 2D graphs ---------------------

void Draw2D(TString name, TString back_path = "QCD", TString sig_path = "EW"){
  // analysis results files
  TString path = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/WZ";
  TFile *EW = TFile::Open(path + sig_path + ".root");
  TFile *QCD = TFile::Open(path + back_path + ".root");
  
  // canvas
  TCanvas *cDrawing2 = new TCanvas("cDrawing2", name, 10, 0, 2048, 820);
  gStyle->SetOptStat(0);
  cDrawing2->Divide(2,1);

  // signal histogram
  TH2F *hEW2 = (TH2F*)EW->Get(name);
  hEW2->SetTitle("Signal");

  // backgroung histogram
  TH2F *hQCD2 = (TH2F*)QCD->Get(name);
  hQCD2->SetTitle("Background");

  // drawing
  cDrawing2->cd(1);
  cDrawing2->cd(1)->SetGridx();
  cDrawing2->cd(1)->SetGridy();
  hEW2->Draw("colz");

  cDrawing2->cd(2);
  cDrawing2->cd(2)->SetGridx();
  cDrawing2->cd(2)->SetGridy();
  hQCD2->Draw("colz");
}

void SSB(TString name, TString back_path = "QCD", TString sig_path = "EW"){
  // analysis results files
  TString path = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/WZ";
  TFile *EW = TFile::Open(path + sig_path + ".root");
  TFile *QCD = TFile::Open(path + back_path + ".root");

  // canvas
  TCanvas *cDrawing = new TCanvas("cDrawingSSB", name, 10, 0, 1280, 1024);
  gStyle->SetOptStat(0);
  cDrawing->SetGrid();

  // signal histogram
  TH1F *hEW = (TH1F*)EW->Get(name);

  // backgroung histogram
  TH1F *hQCD = (TH1F*)QCD->Get(name);

  // signal/sqrt(background) histogram
  TH1F *hSSB = new TH1F("hSSB", name + " #rightarrow #frac{S}{#sqrt{S + B}}", hEW->GetNbinsX(), hEW->GetBinLowEdge(1), hEW->GetBinLowEdge(hEW->GetNbinsX()) + hEW->GetBinWidth(1));

  float signal;
  float background;
  float weight;
  
  for(int i = 1; i <= hEW->GetNbinsX(); i++){
    signal = hEW->Integral(i, -1);
    background = hQCD->Integral(i, -1);
    weight = signal/sqrt(signal + background);

    hSSB->Fill(hEW->GetBinCenter(i), weight);
  }

  // drawing
  cDrawing->cd();
  hSSB->Draw("hist");

  cout << "\nThe maximum is in x = " << hSSB->GetBinLowEdge(hSSB->GetMaximumBin()) << endl << endl;
}
