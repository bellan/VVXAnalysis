#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMath.h>
#include <THStack.h>

using namespace std;

struct test{
  template<class T>
  bool operator()(T const &a, T const &b) const{
    return a > b;
  }
};


// ************************************ Draw1D ************************************

void Draw1D(TString name, int cuts = 0, TString back_path = "WZQCD", TString sig_path = "WZEW"){
  
  // analysis results files
  TString path = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/";
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
  legend->AddEntry(hEW, sig_path);
  legend->AddEntry(hQCD, back_path);
  
  switch(cuts){
  case 0:{    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ case 0 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // canvas
    TCanvas *cDrawing = new TCanvas("cDrawing", name, 10, 0, 1280, 1024);
    gStyle->SetOptStat(0);
    cDrawing->SetGrid();
    
    cout << "\nEW's " << name << " integral is  " << hEW->Integral(0,-1) << endl;
    cout << "QCD's " << name << " integral is " << hQCD->Integral(0,-1) << endl << endl;
    cout << "signal over sqrt(background) = " << hEW->Integral(0,-1)/sqrt(hQCD->Integral(0,-1)) << endl;
    cout << "signal over sqrt(signal + background) = " << hEW->Integral(0,-1)/sqrt(hEW->Integral(0,-1) + hQCD->Integral(0,-1)) << endl << endl;
    
    // drawing
    cDrawing->cd();
    hEW->SetMaximum(1.1 * max(hEW->GetMaximum(), hQCD->GetMaximum()));
    hEW->Draw();
    hQCD->Draw("same");
    legend->Draw();
  }
    break;
    
  case 1:{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ case 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    legendbc->AddEntry(hEWbc, sig_path);
    legendbc->AddEntry(hQCDbc, back_path);

    cout << "\nIntegral before cuts:" << endl; 
    cout << "Signal's " << name << " integral is      " << hEWbc->Integral(0,-1) << endl;
    cout << "Background's " << name << " integral is " << hQCDbc->Integral(0,-1) << endl;
    cout << "signal over sqrt(background) = " << hEWbc->Integral(0,-1)/sqrt(hQCDbc->Integral(0,-1)) << endl;
    cout << "signal over sqrt(signal + background) = " << hEWbc->Integral(0,-1)/sqrt(hEWbc->Integral(0,-1) + hQCDbc->Integral(0,-1)) << endl << endl;

    cout << "\nIntegral after cuts:" << endl; 
    cout << "Signal's " << name << " integral is      " << hEW->Integral(0,-1) << endl;
    cout << "Background's " << name << " integral is " << hQCD->Integral(0,-1) << endl;
    cout << "signal over sqrt(background) = " << hEW->Integral(0,-1)/sqrt(hQCD->Integral(0,-1)) << endl;
    cout << "signal over sqrt(signal + background) = " << hEW->Integral(0,-1)/sqrt(hEW->Integral(0,-1) + hQCD->Integral(0,-1)) << endl << endl;
    
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

  case 2:{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ case 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    delete hEW;
    delete hQCD;
    delete EW;
    delete QCD;
    delete legend;
    
    TFile *filetemp;
    vector<TString> samples = {"DYJetsToLL_M50", "TTJets", "WZ", "TTZToLL", "ZZTo4l", "TTWJets", "WZEW"};
    vector<Color_t> colors = {kRed, kOrange-3, kOrange-2, kGreen+1, kAzure+10, kBlue, kViolet+1, kMagenta};
    vector<TFile *> files;
    vector<TH1F *> histos;
    vector<float> maxs;
    
    // canvas
    TCanvas *cDrawing = new TCanvas("cDrawing", name, 10, 0, 1280, 1024);
    gStyle->SetOptStat(0);

    // legend
    TLegend *legenda2 = new TLegend(0.75, 0.80, 0.95, 0.95);

    // analysis results file + histograms + legend
    for(int i = 0; i < (int)samples.size(); i++){
      filetemp = TFile::Open(path + samples[i] + ".root");
      files.push_back(filetemp);
      
      histos.push_back((TH1F*)files[i]->Get(name));
      histos[i]->SetMarkerSize(0.6);
      histos[i]->SetMarkerStyle(21);
      histos[i]->SetMarkerColor(colors[i]);
      histos[i]->SetLineColor(colors[i]);
      legenda2->AddEntry(histos[i], samples[i]);
	
      maxs.push_back(histos[i]->GetMaximum());
    }

    sort(maxs.begin(), maxs.end(), test());

    histos[0]->SetMaximum(1.1 * maxs[0]);
    
    // drawing
    cDrawing->cd();
    histos[0]->Draw();
    legenda2->Draw();

    for(int i = 1; i < (int)samples.size(); i++){
      histos[i]->Draw("same");
    }
    
  }
    break;
    
  case 3:{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ case 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    delete hEW;
    delete hQCD;
    delete EW;
    delete QCD;
    delete legend;
    
    TFile *filetemp;
    vector<TString> samples = {"DYJetsToLL_M50", "TTJets", "WZ", "TTZToLL", "ZZTo4l", "TTWJets", "WZEW"};
    vector<Color_t> colors = {kRed, kOrange-3, kOrange-2, kGreen+1, kAzure+10, kBlue, kViolet+1};
    vector<TFile *> files;
    vector<TH1F *> histos;
    
    // canvas
    TCanvas *cStack = new TCanvas("cStack", name, 10, 0, 1280, 1024);
    gStyle->SetOptStat(0);
    THStack *hstack = new THStack("hstack", name);

    // legend
    TLegend *legenda2 = new TLegend(0.75, 0.80, 0.95, 0.95);

    // analysis results file + histograms + legend
    for(int i = 0; i < (int)samples.size(); i++){
      filetemp = TFile::Open(path + samples[i] + ".root");
      files.push_back(filetemp);
      
      histos.push_back((TH1F*)files[i]->Get(name));
      histos[i]->SetLineColor(colors[i]);
      histos[i]->SetFillColor(colors[i]);
      legenda2->AddEntry(histos[i], samples[i]);
      hstack->Add(histos[i]);
    }

    // drawing
    cStack->cd();
    hstack->Draw("hist");
    legenda2->Draw();
  }
    break;
    
  default:{
    cout << "Wrong initialization: Draw1D(graph_to_compare, case, background, signal)" << endl;
    cout << "Set case = 0 to compare signal and one background" << endl;
    cout << "Set case = 1 to compare signal and one background before and after all cuts." << endl;
    cout << "Set case = 2 to compare all samples" << endl;
    cout << "Set case = 3 to get stack graphs" << endl;
    
    return;
  }
  }
}


// ********************************** Draw1DNorm **********************************

void Draw1DNorm(TString name, int cuts = 0, TString back_path = "WZQCD", TString sig_path = "WZEW"){
  // analysis results files
  TString path = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/";
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
  legend->AddEntry(hEW, sig_path);
  legend->AddEntry(hQCD, back_path);

  switch(cuts){
  case 0:{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ case 0 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ case 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    legendbc->AddEntry(hEWbc, sig_path);
    legendbc->AddEntry(hQCDbc, back_path);

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

  case 2:{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ case 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    delete hEW;
    delete hQCD;
    delete EW;
    delete QCD;
    delete legend;
    
    TFile *filetemp;
    vector<TString> samples = {"DYJetsToLL_M50", "TTJets", "WZ", "TTZToLL", "ZZTo4l", "TTWJets", "WZEW"};
    vector<Color_t> colors = {kRed, kOrange-3, kOrange-2, kGreen+1, kAzure+10, kBlue, kViolet+1, kMagenta};
    vector<TFile *> files;
    vector<TH1F *> histos;
    vector<float> maxs;
    
    // canvas
    TCanvas *cDrawing = new TCanvas("cDrawing", name, 10, 0, 1280, 1024);
    gStyle->SetOptStat(0);

    // legend
    TLegend *legenda2 = new TLegend(0.75, 0.80, 0.95, 0.95);

    // analysis results file + histograms + legend
    for(int i = 0; i < (int)samples.size(); i++){
      filetemp = TFile::Open(path + samples[i] + ".root");
      files.push_back(filetemp);
      
      histos.push_back((TH1F*)files[i]->Get(name));
      histos[i]->Scale(1/histos[i]->Integral(0, -1));
      histos[i]->SetMarkerSize(0.6);
      histos[i]->SetMarkerStyle(21);
      histos[i]->SetMarkerColor(colors[i]);
      histos[i]->SetLineColor(colors[i]);
      legenda2->AddEntry(histos[i], samples[i]);
	
      maxs.push_back(histos[i]->GetMaximum());
    }

    sort(maxs.begin(), maxs.end(), test());

    histos[0]->SetMaximum(1.1 * maxs[0]);
    // drawing
    cDrawing->cd();
    histos[0]->Draw();
    legenda2->Draw();

    for(int i = 1; i < (int)samples.size(); i++){
      histos[i]->Draw("same");
    }
    
  }
    break;
    
  case 3:{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ case 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    delete hEW;
    delete hQCD;
    delete EW;
    delete QCD;
    delete legend;
    
    TFile *filetemp;
    vector<TString> samples = {"DYJetsToLL_M50", "TTJets", "WZ", "TTZToLL", "ZZTo4l", "TTWJets", "WZEW"};
    vector<Color_t> colors = {kRed, kOrange-3, kOrange-2, kGreen+1, kAzure+10, kBlue, kViolet+1};
    vector<TFile *> files;
    vector<TH1F *> histos;
    
    // canvas
    TCanvas *cStack = new TCanvas("cStack", name, 10, 0, 1280, 1024);
    gStyle->SetOptStat(0);
    THStack *hstack = new THStack("hstack", name);

    // legend
    TLegend *legenda2 = new TLegend(0.75, 0.80, 0.95, 0.95);

    // analysis results file + histograms + legend
    for(int i = 0; i < (int)samples.size(); i++){
      filetemp = TFile::Open(path + samples[i] + ".root");
      files.push_back(filetemp);
      
      histos.push_back((TH1F*)files[i]->Get(name));
      histos[i]->Scale(1/histos[i]->Integral(0, -1));
      histos[i]->SetLineColor(colors[i]);
      histos[i]->SetFillColor(colors[i]);
      legenda2->AddEntry(histos[i], samples[i]);
      hstack->Add(histos[i]);
    }

    // drawing
    cStack->cd();
    hstack->Draw("hist");
    legenda2->Draw();
  }
    break;
    
  default:{
    cout << "Wrong initialization: Draw1D(graph_to_compare, case, background, signal)" << endl;
    cout << "Set case = 0 to compare signal and one background" << endl;
    cout << "Set case = 1 to compare graphs before and after all cuts." << endl;
    cout << "Set case = 2 to compare all samples" << endl;
    cout << "Set case = 3 to get stack graphs" << endl;
    
    return;
  }
  }
}


// ************************************ Draw2D ************************************

void Draw2D(TString name, int cuts = 0, TString back_path = "WZQCD", TString sig_path = "WZEW"){
  // analysis results files
  TString path = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/";
  TFile *EW = TFile::Open(path + sig_path + ".root");
  TFile *QCD = TFile::Open(path + back_path + ".root");

  switch(cuts){
  case 0:{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ case 0 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // canvas
    TCanvas *cDrawing2 = new TCanvas("cDrawing2", name, 10, 0, 1920, 769);
    gStyle->SetOptStat(0);
    cDrawing2->Divide(2,1);
    
    // signal histogram
    TH2F *hEW2 = (TH2F*)EW->Get(name);
    hEW2->SetTitle(sig_path);
    
    // backgroung histogram
    TH2F *hQCD2 = (TH2F*)QCD->Get(name);
    hQCD2->SetTitle(back_path);
    
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
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ case 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    cout << "Wrong initialization: Draw2D(graph_to_compare, case, background, signal)" << endl;
    cout << "Set case = 0 to compare signal and one background" << endl;
    cout << "Set case = 1 to compare graphs before and after all cuts." << endl;
    
    return;
  }    
  }
}


// ************************************* SSB **************************************

void SSB(TString name, int cut = 1, TString back_path = "WZQCD", TString sig_path = "WZEW"){
  
  // analysis results files
  TString path = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/";
    
  // canvas
  TCanvas *cDrawing = new TCanvas("cDrawingSSB", name, 10, 0, 1280, 1024);
  gStyle->SetOptStat(0);
  cDrawing->SetGrid();
    
  float signal;
  float background;
  float weight;

  TH1F *hSSB;
  
  switch(cut){
  case 0:{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ case 0 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // analysis results files
    TFile *EW = TFile::Open(path + sig_path + ".root");
    TFile *QCD = TFile::Open(path + back_path + ".root");
    
    // signal histogram
    TH1F *hEW = (TH1F*)EW->Get(name);
    
    // backgroung histogram
    TH1F *hQCD = (TH1F*)QCD->Get(name);
    
    // signal/sqrt(background) histogram
    TH1F *hSSB = new TH1F("hSSB", name + " #rightarrow #frac{S}{#sqrt{B}}", hEW->GetNbinsX(), hEW->GetBinLowEdge(1), hEW->GetBinLowEdge(hEW->GetNbinsX()) + hEW->GetBinWidth(hEW->GetNbinsX()));
    
    for(int i = 1; i <= hEW->GetNbinsX(); i++){
      signal = hEW->Integral(i, -1);
      background = hQCD->Integral(i, -1);
      weight = signal/sqrt(background);
      
      hSSB->Fill(hEW->GetBinCenter(i), weight);
    }
    
    // drawing
    cDrawing->cd();
    hSSB->Draw("hist");
    
    cout << "\nThe maximum is in x = " << hSSB->GetBinLowEdge(hSSB->GetMaximumBin()) << endl << endl;
  }
    break;

  case 1:{ 
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ case 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    TFile *filetemp;
    vector<TString> samples = {"DYJetsToLL_M50", "TTJets", "WZ", "TTZToLL", "ZZTo4l", "TTWJets", "WZEW"};
    vector<TFile *> files;
    vector<TH1F *> histos;
    
    // analysis results files + original histograms
    for(int i = 0; i < (int)samples.size(); i++){
      filetemp = TFile::Open(path + samples[i] + ".root");
      files.push_back(filetemp);
      
      histos.push_back((TH1F*)files[i]->Get(name));
    }
    
    // signal/sqrt(background) histogram
    TH1F *hSSB = new TH1F("hSSB", name + " #rightarrow #frac{S}{#sqrt{#sum B}}", histos[0]->GetNbinsX(), histos[0]->GetBinLowEdge(1), histos[0]->GetBinLowEdge(histos[0]->GetNbinsX()) + histos[0]->GetBinWidth(histos[0]->GetNbinsX()));
    
    for(int i = 1; i <= histos[0]->GetNbinsX(); i++){

      for(int j = 0; j <= samples.size() - 2; j++){
	background = background + histos[j]->Integral(i, -1);
      }
      
      signal = histos[samples.size() - 1]->Integral(i, -1);
      weight = signal/sqrt(background);
      hSSB->Fill(histos[0]->GetBinCenter(i), weight);
    }
    
    // drawing
    cDrawing->cd();
    hSSB->Draw("hist");
    
    cout << "\nThe maximum is in x = " << hSSB->GetBinLowEdge(hSSB->GetMaximumBin()) << endl << endl;
  }
    break;

  default:{
    cout << "Wrong initialization: SSB(graph, case, background, signal)" << endl;
    cout << "Set case = 0 to get the significance plot for one background" << endl;
    cout << "Set case = 1 to get the significance plot for all background" << endl;
    
    return;
  }
  }
}


// ************************************* Cuts *************************************

void Cuts(int cuts = 0, TString back_path = "WZQCD", TString sig_path = "WZEW"){
  TString path = "~/Tesi/VVXAnalysis/TreeAnalysis/results/WZAnalyzer_MC/";
  float weight;
  vector<TString> namebin = {"All events", "ZLCand's size > 0", "3^{rd} lep pass full sel", "60 < m_{Z} < 120 #frac{GeV}{c^{2}}", "jets' size > 2", "p_{T, MET} > 30 #frac{GeV}{c}", "30 < m_{T, W^{#pm}} < 500 #frac{GeV}{c^{2}}", "m_{J_{0}J_{1}} > 258,5 #frac{GeV}{c^{2}}", "|#Delta#eta| > 2,1", "Z_{llJ_{0}} > 0,6"};
  
  // canvas1
  TCanvas *cCuts = new TCanvas("cCuts", "RecoCuts1", 10, 0, 1280, 1024);
  gStyle->SetOptStat(0);
  cCuts->SetGrid();
  
  // canvas 2
  TCanvas *cCuts2 = new TCanvas("cCuts2", "RecoCuts2", 10, 0, 1820, 1024);
  gStyle->SetOptStat(0);
  cCuts2->SetGrid();
  
  // canvas 3
  TCanvas *cCuts3 = new TCanvas("cCuts3", "RecoCuts3", 10, 0, 1820, 1024);
  gStyle->SetOptStat(0);
  cCuts3->SetGrid();
  
  // legend
  TLegend *legend = new TLegend(0.75, 0.80, 0.95, 0.95);
  TLegend *legend2 = new TLegend(0.75, 0.80, 0.95, 0.95);

  switch(cuts){
  case 0:{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ case 0 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TFile *EW = TFile::Open(path + sig_path + ".root");
    TFile *QCD = TFile::Open(path + back_path + ".root");

    // -------------------------------- Cuts 1 --------------------------------
    // signal histogram
    TH1F *hEW = (TH1F*)EW->Get("RecoCuts");
    hEW->SetMarkerColor(kGreen+2);
    hEW->SetMarkerSize(0.6);
    hEW->SetMarkerStyle(21);
    hEW->SetLineColor(kGreen+2);
    
    // backgroung histogram
    TH1F *hQCD = (TH1F*)QCD->Get("RecoCuts");
    hQCD->SetMarkerColor(kRed);
    hQCD->SetMarkerSize(0.6);
    hQCD->SetMarkerStyle(21);
    hQCD->SetLineColor(kRed);

    for(int i = 1; i <= (int)namebin.size(); i++){
      hEW->GetXaxis()->SetBinLabel(i, namebin[i - 1]);   
    }
    
    // drawing
    cCuts->cd();
    hEW->SetMaximum(1.1 * max(hEW->GetMaximum(), hQCD->GetMaximum()));
    hEW->Draw();
    hQCD->Draw("same");

    // -------------------------------- Cuts 2 --------------------------------
    // signal histogram
    TH1F *hEW2 = new TH1F("hEW2", "Recocuts", hEW->GetNbinsX(), hEW->GetBinLowEdge(1), hEW->GetBinLowEdge(hEW->GetNbinsX()) + hEW->GetBinWidth(hEW->GetNbinsX()));
    hEW2->SetLineColor(kGreen+2);
    hEW2->SetFillColor(kGreen+2);
    hEW2->SetFillStyle(3002);
    
    // backgroung histogram
    TH1F *hQCD2 = new TH1F("hQCD2", "Recocuts", hQCD->GetNbinsX(), hQCD->GetBinLowEdge(1), hQCD->GetBinLowEdge(hQCD->GetNbinsX()) + hQCD->GetBinWidth(hQCD->GetNbinsX()));
    hQCD2->SetLineColor(kRed);
    hQCD2->SetFillColor(kRed);
    hQCD2->SetFillStyle(3004);
    
    weight = hEW->GetBinContent(1)/hEW->GetBinContent(1);
    hEW2->Fill(0., weight);
    
    weight = hQCD->GetBinContent(1)/hQCD->GetBinContent(1);
    hQCD2->Fill(0., weight);
    
    for(int i = 2; i < hEW->GetNbinsX() + 1; i++){
      weight =  hEW->GetBinContent(i)/hEW->GetBinContent(1);
      hEW2->Fill(i-1, weight);
      
      weight = hQCD->GetBinContent(i)/hQCD->GetBinContent(1);
      hQCD2->Fill(i-1, weight);
    }

    for(int i = 1; i <= (int)namebin.size(); i++){
      hEW->GetXaxis()->SetBinLabel(i, namebin[i - 1]);     
    }
    
    // drawing
    cCuts2->cd();
    hEW2->Draw("histo");
    hQCD2->Draw("histosame");

    // -------------------------------- Cuts 3 --------------------------------
    // signal histogram
    TH1F *hEW3 = new TH1F("hEW3", "Recocuts", hEW->GetNbinsX(), hEW->GetBinLowEdge(1), hEW->GetBinLowEdge(hEW->GetNbinsX()) + hEW->GetBinWidth(hEW->GetNbinsX()));
    hEW3->SetLineColor(kGreen+2);
    hEW3->SetFillColor(kGreen+2);
    hEW3->SetFillStyle(3002);
    
    // backgroung histogram
    TH1F *hQCD3 = new TH1F("hQCD3", "Recocuts", hQCD->GetNbinsX(), hQCD->GetBinLowEdge(1), hQCD->GetBinLowEdge(hQCD->GetNbinsX()) + hQCD->GetBinWidth(hQCD->GetNbinsX()));
    hQCD3->SetLineColor(kRed);
    hQCD3->SetFillColor(kRed);
    hQCD3->SetFillStyle(3004);
    
    weight = hEW->GetBinContent(1)/hEW->GetBinContent(1);
    hEW3->Fill(0., weight);
    
    weight = hQCD->GetBinContent(1)/hQCD->GetBinContent(1);
    hQCD3->Fill(0., weight);
    
    for(int i = 2; i < 13; i++){
      if(hQCD->GetBinContent(i-1) > 0){
	weight = hQCD->GetBinContent(i)/hQCD->GetBinContent(i-1);
	hQCD3->Fill(i-1, weight);
      }
      else{
	hQCD3->Fill(i-1, 0);
      }
      
      if(hEW->GetBinContent(i-1) > 0){
	weight = hEW->GetBinContent(i)/hEW->GetBinContent(i-1);
	hEW3->Fill(i-1, weight);
      }
      else{
	hEW3->Fill(i-1, 0);
      }
    }

    for(int i = 1; i <= (int)namebin.size(); i++){
      hEW->GetXaxis()->SetBinLabel(i, namebin[i - 1]);   
    }
    
    // drawing
    cCuts3->cd();
    hEW3->Draw("histo");
    hQCD3->Draw("histosame");

    
    // legend
    legend->AddEntry(hEW, sig_path);
    legend->AddEntry(hQCD, back_path);

    legend2->AddEntry(hEW, sig_path);
    legend2->AddEntry(hQCD, back_path);
  }
    break;
    
  case 1:{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ case 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    TFile *filetemp;    
    vector<TString> samples = {"100/DYJetsToLL_M50", "100/TTJets", "100/WZ", "100/TTZToLL", "100/ZZTo4l", "100/TTWJets", "100/WZEW"};
    vector<Color_t> colors = {kRed, kOrange-3, kOrange-2, kGreen+1, kAzure+10, kBlue, kViolet+1, kMagenta};
    vector<Style_t> style = {3004, 3005, 3006, 3007, 3016, 3020, 3002};
    vector<TFile *> files;
    vector<TH1F *> histos;
    vector<TH1F *> histos2;
    vector<TH1F *> histos3;
    vector<float> maxs;

    // -------------------------------- Cuts 1 --------------------------------
    // analysis results file + histograms + legend
    for(int i = 0; i < (int)samples.size(); i++){
      filetemp = TFile::Open(path + samples[i] + ".root");
      files.push_back(filetemp);
      
      histos.push_back((TH1F*)files[i]->Get("RecoCuts"));
      histos[i]->SetMarkerSize(0.6);
      histos[i]->SetMarkerStyle(21);
      histos[i]->SetMarkerColor(colors[i]);
      histos[i]->SetLineColor(colors[i]);
      legend->AddEntry(histos[i], samples[i]);
	
      maxs.push_back(histos[i]->GetMaximum());
    }

    sort(maxs.begin(), maxs.end(), test());

    // -------------------------------- Cuts 2 --------------------------------
    THStack *hsCuts2 = new THStack("hsCuts2", "RecoCuts");
    
    for(int i = 0; i < (int)samples.size(); i++){
      histos2.push_back(new TH1F("histos2", "RecoCuts", histos[i]->GetNbinsX(), histos[i]->GetBinLowEdge(1), histos[i]->GetBinLowEdge(histos[i]->GetNbinsX()) + histos[i]->GetBinWidth(histos[i]->GetNbinsX())));
      histos2[i]->SetLineColor(colors[i]);
      histos2[i]->SetFillColor(colors[i]);

      weight = histos[i]->GetBinContent(1)/histos[i]->GetBinContent(1);
      histos2[i]->Fill(0., weight);

      for(int j = 2; j < histos[i]->GetNbinsX() + 1; j++){
	weight = histos[i]->GetBinContent(j)/histos[i]->GetBinContent(1);
	histos2[i]->Fill(j-1, weight);
      }

      hsCuts2->Add(histos2[i]);
      legend2->AddEntry(histos2[i], samples[i]);
    }

    // -------------------------------- Cuts 3 --------------------------------
    THStack *hsCuts3 = new THStack("hsCuts3", "RecoCuts");
    
    for(int i = 0; i < (int)samples.size(); i++){
      histos3.push_back(new TH1F("histos3", "RecoCuts", histos[i]->GetNbinsX(), histos[i]->GetBinLowEdge(1), histos[i]->GetBinLowEdge(histos[i]->GetNbinsX()) + histos[i]->GetBinWidth(histos[i]->GetNbinsX())));
      histos3[i]->SetLineColor(colors[i]);
      histos3[i]->SetFillColor(colors[i]);

      weight = histos[i]->GetBinContent(1)/histos[i]->GetBinContent(1);
      histos3[i]->Fill(0., weight);

      for(int j = 2; j < histos[i]->GetNbinsX() + 1; j++){
	if(histos[i]->GetBinContent(j) > 0){
	  weight = histos[i]->GetBinContent(j)/histos[i]->GetBinContent(j-1);
	  histos3[i]->Fill(j-1, weight);
	}
	else{
	  histos3[i]->Fill(j-1, 0);
	}
      }

      hsCuts3->Add(histos3[i]);
    }

    
    // drawing    
    for(int i = 1; i <= (int)namebin.size(); i++){
      histos[0]->GetXaxis()->SetBinLabel(i, namebin[i - 1]);
      histos2[0]->GetXaxis()->SetBinLabel(i, namebin[i - 1]);
      histos3[0]->GetXaxis()->SetBinLabel(i, namebin[i - 1]);
    }
    
    cCuts->cd();
    histos[0]->SetMaximum(1.1 * maxs[0]);
    histos[0]->Draw();
    for(int i = 1; i < (int)samples.size(); i++){
      histos[i]->Draw("same");
    }
    
    cCuts2->cd();
    hsCuts2->Draw("histonostackb");
    
    cCuts3->cd();
    hsCuts3->Draw("histonostackb");
    
  }
    break;
    
  default:{
    cout << "Wrong initialization: Cuts(case, background, signal)" << endl;
    cout << "Set case = 0 to compare the cuts for two samples" << endl;
    cout << "Set case = 1 to compare the cuts for all samples" << endl;
  }
  }
  
  // drawing the legend
  cCuts->cd();
  legend->Draw();
  
  cCuts2->cd();
  legend2->Draw();
  
  cCuts3->cd();
  legend2->Draw();
}
