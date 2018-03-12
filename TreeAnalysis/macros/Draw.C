#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMath.h>

using namespace std;

// --------------------- 1D graphs --------------------- 
void Draw1D(TString name, int cuts = 0, TString back_path = "QCD", TString sig_path = "EW"){
  
  // analysis results files
  TString path = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/WZ";
  TFile *EW = TFile::Open(path + sig_path + ".root");
  TFile *QCD = TFile::Open(path + back_path + ".root");
    
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
  
  switch(cuts){
  case 0:{    
    // canvas
    TCanvas *cDrawing = new TCanvas("cDrawing", name, 10, 0, 1280, 1024);
    gStyle->SetOptStat(0);
    cDrawing->SetGrid();
    
    cout << "\nEW's " << name << " integral is  " << hEW->Integral(0,-1) << endl;
    cout << "QCD's " << name << " integral is " << hQCD->Integral(0,-1) << endl << endl;
    
    // drawing
    cDrawing->cd();
    hEW->SetMaximum(1.1 * max(hEW->GetMaximum(), hQCD->GetMaximum()));
    hEW->Draw();
    hQCD->Draw("same");
    legend->Draw();
  }
    break;
    
  case 1:{
    // analysis results files
    TFile *EWbc = TFile::Open(path + sig_path + "bc.root");
    TFile *QCDbc = TFile::Open(path + back_path + "bc.root");
    
    // canvas
    TCanvas *cDrawing = new TCanvas("cDrawing", name, 10, 0, 1920, 769);
    cDrawing->Divide(2,1);
    gStyle->SetOptStat(0);
    
    // signal histograms
    hEW->SetTitle("After cuts");
    
    TH1F *hEWbc = (TH1F*)EWbc->Get(name);
    hEWbc->SetMarkerColor(kGreen+2);
    hEWbc->SetMarkerSize(0.6);
    hEWbc->SetMarkerStyle(21);
    hEWbc->SetLineColor(kGreen+2);
    hEWbc->SetTitle("Before cuts");
    
    // backgroung histograms    
    TH1F *hQCDbc = (TH1F*)QCDbc->Get(name);
    hQCDbc->SetMarkerColor(kRed);
    hQCDbc->SetMarkerSize(0.6);
    hQCDbc->SetMarkerStyle(21);
    hQCDbc->SetLineColor(kRed);
    
    // legend    
    TLegend *legendbc = new TLegend(0.75, 0.80, 0.95, 0.95);
    legendbc->AddEntry(hEWbc, "signal");
    legendbc->AddEntry(hQCDbc, "background");

    cout << "\nIntegral before cuts:" << endl; 
    cout << "Signal's " << name << " integral is      " << hEWbc->Integral(0,-1) << endl;
    cout << "Background's " << name << " integral is " << hQCDbc->Integral(0,-1) << endl << endl;

    cout << "\nIntegral after cuts:" << endl; 
    cout << "Signal's " << name << " integral is      " << hEW->Integral(0,-1) << endl;
    cout << "Background's " << name << " integral is " << hQCD->Integral(0,-1) << endl << endl;
    
    // drawing
    cDrawing->cd(1);
    cDrawing->cd(1)->SetGrid();
    hEWbc-> SetMaximum(1.1 * max(hEWbc->GetMaximum(), hQCDbc->GetMaximum()));
    hEWbc->Draw();
    hQCDbc->Draw("same");
    legendbc->Draw();
    
    cDrawing->cd(2);
    cDrawing->cd(2)->SetGrid();
    hEW-> SetMaximum(1.1 * max(hEW->GetMaximum(), hQCD->GetMaximum()));
    hEW->Draw();
    hQCD->Draw("same");
    legend->Draw();
  }
    break;
    
  default:{
    cout << "Wrong initialization: Draw1D(graph_to_compare, cuts, background, signal)" << endl;
    cout << "Set cuts = 1 to compare graphs before and after all cuts." << endl;
    
    return;
  }
  }
}

void Draw1DNorm(TString name, int cuts = 0, TString back_path = "QCD", TString sig_path = "EW"){
  // analysis results files
  TString path = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/WZ";
  TFile *EW = TFile::Open(path + sig_path + ".root");
  TFile *QCD = TFile::Open(path + back_path + ".root");

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

  switch(cuts){
  case 0:{
    // canvas
    TCanvas *cDrawing = new TCanvas("cDrawing", name, 10, 0, 1280, 1024);
    gStyle->SetOptStat(0);
    cDrawing->SetGrid();
    
    // drawing
    cDrawing->cd();
    hEW-> SetMaximum(1.1 * max(hEW->GetMaximum(), hQCD->GetMaximum()));
    hEW->Draw();
    hQCD->Draw("same");
    legend->Draw();
  }
    break;
    
  case 1:{
    // analysis results files
    TFile *EWbc = TFile::Open(path + sig_path + "bc.root");
    TFile *QCDbc = TFile::Open(path + back_path + "bc.root");
    
    // canvas
    TCanvas *cDrawing = new TCanvas("cDrawing", name, 10, 0, 1920, 769);
    cDrawing->Divide(2,1);
    gStyle->SetOptStat(0);
    
    // signal histograms
    hEW->SetTitle("After cuts");
    
    TH1F *hEWbc = (TH1F*)EWbc->Get(name);
    hEWbc->Scale(1/hEWbc->Integral(0, -1));
    hEWbc->SetMarkerColor(kGreen+2);
    hEWbc->SetMarkerSize(0.6);
    hEWbc->SetMarkerStyle(21);
    hEWbc->SetLineColor(kGreen+2);
    hEWbc->SetTitle("Before cuts");
    
    // backgroung histograms    
    TH1F *hQCDbc = (TH1F*)QCDbc->Get(name);
    hQCDbc->Scale(1/hQCDbc->Integral(0, -1));
    hQCDbc->SetMarkerColor(kRed);
    hQCDbc->SetMarkerSize(0.6);
    hQCDbc->SetMarkerStyle(21);
    hQCDbc->SetLineColor(kRed);
    
    // legend    
    TLegend *legendbc = new TLegend(0.75, 0.80, 0.95, 0.95);
    legendbc->AddEntry(hEWbc, "signal");
    legendbc->AddEntry(hQCDbc, "background");

    cout << "\nIntegral before cuts:" << endl; 
    cout << "Signal's " << name << " integral is      " << hEWbc->Integral(0,-1) << endl;
    cout << "Background's " << name << " integral is " << hQCDbc->Integral(0,-1) << endl << endl;

    cout << "\nIntegral after cuts:" << endl; 
    cout << "Signal's " << name << " integral is      " << hEW->Integral(0,-1) << endl;
    cout << "Background's " << name << " integral is " << hQCD->Integral(0,-1) << endl << endl;
    
    // drawing
    cDrawing->cd(1);
    cDrawing->cd(1)->SetGrid();
    hEWbc-> SetMaximum(1.1 * max(hEWbc->GetMaximum(), hQCDbc->GetMaximum()));
    hEWbc->Draw();
    hQCDbc->Draw("same");
    legendbc->Draw();
    
    cDrawing->cd(2);
    cDrawing->cd(2)->SetGrid();
    hEW-> SetMaximum(1.1 * max(hEW->GetMaximum(), hQCD->GetMaximum()));
    hEW->Draw();
    hQCD->Draw("same");
    legend->Draw();    
  }
    break;
    
  default:{
    cout << "Wrong initialization: Draw1DNorm(graph_to_compare, cuts, background, signal)" << endl;
    cout << "Set cuts = 1 to compare graphs before and after all cuts." << endl;
    
    return;
  }    
  }
}

// --------------------- 2D graphs ---------------------

void Draw2D(TString name, int cuts = 0, TString back_path = "QCD", TString sig_path = "EW"){
  // analysis results files
  TString path = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/WZ";
  TFile *EW = TFile::Open(path + sig_path + ".root");
  TFile *QCD = TFile::Open(path + back_path + ".root");

  switch(cuts){
  case 0:{
    // canvas
    TCanvas *cDrawing2 = new TCanvas("cDrawing2", name, 10, 0, 1920, 769);
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
    break;

  case 1:{
    // analysis results files
    TFile *EWbc = TFile::Open(path + sig_path + "bc.root");
    TFile *QCDbc = TFile::Open(path + back_path + "bc.root");
    
    // canvas
    TCanvas *cDrawing2 = new TCanvas("cDrawing2", name, 10, 0, 1186, 950);
    gStyle->SetOptStat(0);
    cDrawing2->Divide(2,2);
    
    // signal histogram
    TH2F *hEW2 = (TH2F*)EW->Get(name);
    hEW2->SetTitle("Signal after cuts");
    TH2F *hEW2bc = (TH2F*)EWbc->Get(name);
    hEW2->SetTitle("Signal before cuts");
    
    // backgroung histogram
    TH2F *hQCD2 = (TH2F*)QCD->Get(name);
    hQCD2->SetTitle("Background after cuts");
    TH2F *hQCD2bc = (TH2F*)QCDbc->Get(name);
    hQCD2->SetTitle("Background before cuts");
    
    // drawing
    cDrawing2->cd(1);
    cDrawing2->cd(1)->SetGridx();
    cDrawing2->cd(1)->SetGridy();
    hEW2bc->Draw("colz");
    
    cDrawing2->cd(2);
    cDrawing2->cd(2)->SetGridx();
    cDrawing2->cd(2)->SetGridy();
    hQCD2bc->Draw("colz");
    
    cDrawing2->cd(3);
    cDrawing2->cd(3)->SetGridx();
    cDrawing2->cd(3)->SetGridy();
    hEW2->Draw("colz");
    
    cDrawing2->cd(4);
    cDrawing2->cd(4)->SetGridx();
    cDrawing2->cd(4)->SetGridy();
    hQCD2->Draw("colz");
  }
    break;
    
  default:{
    cout << "Wrong initialization: Draw2D(graph_to_compare, cuts, background, signal)" << endl;
    cout << "Set cuts = 1 to compare graphs before and after all cuts." << endl;
    
    return;
  }    
  }
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
