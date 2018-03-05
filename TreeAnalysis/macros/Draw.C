#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMath.h>

using namespace std;

// --------------------- 1D graphs --------------------- 
void Draw1D(TString nome, TString perfondo = "QCD", TString persegnale = "EW"){
  // analysis results files
  TString per = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/WZ";
  TFile *EW = TFile::Open(per + persegnale + ".root");
  TFile *QCD = TFile::Open(per + perfondo + ".root");

  // canvas
  TCanvas *cDrawing = new TCanvas("cDrawing", nome, 10, 0, 1280, 1024);
  gStyle->SetOptStat(0);
  cDrawing->SetGrid();

  // signal histogram
  TH1F *hEW = (TH1F*)EW->Get(nome);
  hEW->SetMarkerColor(kGreen+2);
  hEW->SetMarkerSize(0.6);
  hEW->SetMarkerStyle(21);
  hEW->SetLineColor(kGreen+2);

  // backgroung histogram
  TH1F *hQCD = (TH1F*)QCD->Get(nome);
  hQCD->SetMarkerColor(kRed);
  hQCD->SetMarkerSize(0.6);
  hQCD->SetMarkerStyle(21);
  hQCD->SetLineColor(kRed);

  // legend
  TLegend *legend = new TLegend(0.75, 0.80, 0.95, 0.95);
  legend->AddEntry(hEW, "signal");
  legend->AddEntry(hQCD, "background");

  cout << "EW's " << nome << " integral is  " << hEW->Integral(0,-1) << endl;
  cout << "QCD's " << nome << " integral is " << hQCD->Integral(0,-1) << endl;

  // drawing
  cDrawing->cd();
  hEW-> SetMaximum(1.1 * max(hEW->GetMaximum(), hQCD->GetMaximum()));
  hEW->Draw();
  hQCD->Draw("same");
  legend->Draw();
}

void Draw1DNorm(TString nome, TString perfondo = "QCD", TString persegnale = "EW"){
  // analysis results files
  TString per = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/WZ";
  TFile *EW = TFile::Open(per + persegnale + ".root");
  TFile *QCD = TFile::Open(per + perfondo + ".root");

  // canvas
  TCanvas *cDrawing = new TCanvas("cDrawing", nome, 10, 0, 1280, 1024);
  gStyle->SetOptStat(0);
  cDrawing->SetGrid();

  // signal histogram
  TH1F *hEW = (TH1F*)EW->Get(nome);
  hEW->Scale(1/hEW->Integral(0, -1));
  hEW->SetMarkerSize(0.6);
  hEW->SetMarkerStyle(21);
  hEW->SetMarkerColor(kGreen+2);
  hEW->SetLineColor(kGreen+2);

  // backgroung histogram
  TH1F *hQCD = (TH1F*)QCD->Get(nome);
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

void Draw2D(TString nome, TString perfondo = "QCD", TString persegnale = "EW"){
  // analysis results files
  TString per = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/WZ";
  TFile *EW = TFile::Open(per + persegnale + ".root");
  TFile *QCD = TFile::Open(per + perfondo + ".root");

  /*
  // canvas 1
  TCanvas *cDrawing = new TCanvas("cDrawing", nome, 10, 0, 1280, 1024);
  gStyle->SetOptStat(0);
  cDrawing->SetGrid();

  // signal histogram
  TH2F *hEW = (TH2F*)EW->Get(nome);
  hEW->SetMarkerColor(kGreen+2);
  hEW->SetLineColor(kGreen+2);

  // backgroung histogram
  TH2F *hQCD = (TH2F*)QCD->Get(nome);
  hQCD->SetMarkerColor(kRed);
  hQCD->SetLineColor(kRed);

  // legend
  TLegend *legend = new TLegend(0.75, 0.80, 0.95, 0.95);
  legend->AddEntry(hEW, "signal");
  legend->AddEntry(hQCD, "background");

  // drawing
  cDrawing->cd();
  hEW->Draw("surf");
  hQCD->Draw("surfsame");
  legend->Draw();  
  */
  
  // canvas 2
  TCanvas *cDrawing2 = new TCanvas("cDrawing2", nome, 10, 0, 2048, 820);
  gStyle->SetOptStat(0);
  cDrawing2->Divide(2,1);

  // signal histogram
  TH2F *hEW2 = (TH2F*)EW->Get(nome);
  hEW2->SetTitle("Signal");

  // backgroung histogram
  TH2F *hQCD2 = (TH2F*)QCD->Get(nome);
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

void Draw1Dbc(TString nome, TString perfondo = "QCD", TString persegnale = "EW"){
  // analysis results files
  TString per = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/WZ";
  TFile *EW = TFile::Open(per + persegnale + "bc.root");
  TFile *QCD = TFile::Open(per + perfondo + "bc.root");

  // canvas
  TCanvas *cDrawingBC = new TCanvas("cDrawingBC", nome, 10, 0, 1280, 1024);
  gStyle->SetOptStat(0);
  cDrawingBC->SetGrid();

  // signal histogram
  TH1F *hEW = (TH1F*)EW->Get(nome);
  hEW->SetMarkerColor(kGreen+2);
  hEW->SetMarkerSize(0.6);
  hEW->SetMarkerStyle(21);
  hEW->SetLineColor(kGreen+2);

  // backgroung histogram
  TH1F *hQCD = (TH1F*)QCD->Get(nome);
  hQCD->SetMarkerColor(kRed);
  hQCD->SetMarkerSize(0.6);
  hQCD->SetMarkerStyle(21);
  hQCD->SetLineColor(kRed);

  // legend
  TLegend *legend = new TLegend(0.75, 0.80, 0.95, 0.95);
  legend->AddEntry(hEW, "signal");
  legend->AddEntry(hQCD, "background");

  cout << "EW's " << nome << " integral is  " << hEW->Integral(0,-1) << endl;
  cout << "QCD's " << nome << " integral is " << hQCD->Integral(0,-1) << endl;

  // drawing
  cDrawingBC->cd();
  hEW-> SetMaximum(1.1 * max(hEW->GetMaximum(), hQCD->GetMaximum()));
  hEW->Draw();
  hQCD->Draw("same");
  legend->Draw();
}

void Draw1DNormbc(TString nome, TString perfondo = "QCD", TString persegnale = "EW"){
  // analysis results files
  TString per = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/WZ";
  TFile *EW = TFile::Open(per + persegnale + "bc.root");
  TFile *QCD = TFile::Open(per + perfondo + "bc.root");

  // canvas
  TCanvas *cDrawingBC = new TCanvas("cDrawingBC", nome, 10, 0, 1280, 1024);
  gStyle->SetOptStat(0);
  cDrawingBC->SetGrid();

  // signal histogram
  TH1F *hEW = (TH1F*)EW->Get(nome);
  hEW->Scale(1/hEW->Integral(0, -1));
  hEW->SetMarkerSize(0.6);
  hEW->SetMarkerStyle(21);
  hEW->SetMarkerColor(kGreen+2);
  hEW->SetLineColor(kGreen+2);

  // backgroung histogram
  TH1F *hQCD = (TH1F*)QCD->Get(nome);
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
  cDrawingBC->cd();
  hEW-> SetMaximum(1.1 * max(hEW->GetMaximum(), hQCD->GetMaximum()));
  hEW->Draw();
  hQCD->Draw("same");
  legend->Draw();
}

// --------------------- 2D graphs ---------------------

void Draw2Dbc(TString nome, TString perfondo = "QCD", TString persegnale = "EW"){
  // analysis results files
  TString per = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/WZ";
  TFile *EW = TFile::Open(per + persegnale + "bc.root");
  TFile *QCD = TFile::Open(per + perfondo + "bc.root");
  
  // canvas
  TCanvas *cDrawingBC = new TCanvas("cDrawing2BC", nome, 10, 0, 2048, 820);
  gStyle->SetOptStat(0);
  cDrawingBC->Divide(2,1);

  // signal histogram
  TH2F *hEW2 = (TH2F*)EW->Get(nome);
  hEW2->SetTitle("Signal");

  // backgroung histogram
  TH2F *hQCD2 = (TH2F*)QCD->Get(nome);
  hQCD2->SetTitle("Background");

  // drawing
  cDrawingBC->cd(1);
  cDrawingBC->cd(1)->SetGridx();
  cDrawingBC->cd(1)->SetGridy();
  hEW2->Draw("colz");

  cDrawingBC->cd(2);
  cDrawingBC->cd(2)->SetGridx();
  cDrawingBC->cd(2)->SetGridy();
  hQCD2->Draw("colz");
}

void SSB(TString nome, TString perfondo = "QCD", TString persegnale = "EW"){
  // analysis results files
  TString per = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/WZ";
  TFile *EW = TFile::Open(per + persegnale + ".root");
  TFile *QCD = TFile::Open(per + perfondo + ".root");

  // canvas
  TCanvas *cDrawing = new TCanvas("cDrawingSSB", nome, 10, 0, 1280, 1024);
  gStyle->SetOptStat(0);
  cDrawing->SetGrid();

  // signal histogram
  TH1F *hEW = (TH1F*)EW->Get(nome);

  // backgroung histogram
  TH1F *hQCD = (TH1F*)QCD->Get(nome);

  // signal/sqrt(background) histogram
  TH1F *hSSB = new TH1F("hSSB", nome + " #rightarrow #frac{S}{#sqrt{S + B}}", hQCD->GetNbinsX(), hEW->GetBinLowEdge(1), hEW->GetBinLowEdge(hEW->GetNbinsX()) + hEW->GetBinWidth(1));

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

  cout << "The maximum is in x = " << hSSB->GetBinLowEdge(hSSB->GetMaximumBin()) << endl;
}
