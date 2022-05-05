#include "TH2F.h"
#include "TLine.h"
#include "TLatex.h"

#include "../../Commons/interface/Constants.h"

void XiPlotter(TH2F *hist){
  TLine *line1x = new TLine(0.04,0.04,0.04,0.2);
  TLine *line2x = new TLine(0.2,0.04,0.2,0.2);
  TLine *line1y = new TLine(0.04,0.04,0.2,0.04);
  TLine *line2y = new TLine(0.04,0.2,0.2,0.2);

  hist->Draw();
			   
  line1x->SetLineColor(kBlack);
  line1x->SetLineWidth(3);
  line1x->Draw("same");
  line2x->SetLineColor(kBlack);
  line2x->SetLineWidth(3);
  line2x->Draw("same");
  line1y->SetLineColor(kBlack);
  line1y->SetLineWidth(3);
  line1y->Draw("same");
  line2y->SetLineColor(kBlack);
  line2y->SetLineWidth(3);
  line2y->Draw("same");
}


void MatchingPlotter(TH2F *histmatch){
  histmatch->SetTitle("2D matching distribution (with pileup)");
  histmatch->GetXaxis()->SetTitle("#bf{1-m_{ZZ}/m_{pp}}");
  histmatch->GetYaxis()->SetTitle("#bf{y_{pp}- y_{ZZ}}");
  gStyle->SetOptStat(0);
  histmatch->DrawCopy("colz");
  
  TLegend *l = new TLegend(0.6,0.7,0.95,0.8);
  l->AddEntry(histmatch,"a^{0}_{Z}/#Lambda^{2} = 5*10^{-5} GeV^{-2}","f");
  l->Draw("same");

  TLine *line1 = new TLine(-0.5,0.4,0,-0.1);
  line1->SetLineColor(kYellow);
  line1->SetLineWidth(4);
  line1->Draw("same");
  TLine *line2 = new TLine(-1.5,0.9,-0.5,0.4);
  line2->SetLineColor(kYellow);
  line2->SetLineWidth(4);
  line2->Draw("same");
  TLine *line3 = new TLine(-2,1.1,-1.5,0.9);
  line3->SetLineColor(kYellow);
  line3->SetLineWidth(4);
  line3->Draw("same");
  TLine *line4 = new TLine(-0.5,0.6,0,0.1);
  line4->SetLineColor(kYellow);
  line4->SetLineWidth(4);
  line4->Draw("same");
  TLine *line5 = new TLine(-1.5,1.1,-0.5,0.6);
  line5->SetLineColor(kYellow);
  line5->SetLineWidth(4);
  line5->Draw("same");
  TLine *line6 = new TLine(-2,1.3,-1.5,1.1);
  line6->SetLineColor(kYellow);
  line6->SetLineWidth(4);
  line6->Draw("same");
  TLine *line7 = new TLine(-0.5,-0.4,0,0.1);
  line7->SetLineColor(kYellow);
  line7->SetLineWidth(4);
  line7->Draw("same");
  TLine *line8 = new TLine(-1.5,-0.9,-0.5,-0.4);
  line8->SetLineColor(kYellow);
  line8->SetLineWidth(4);
  line8->Draw("same");
  TLine *line9 = new TLine(-2,-1.1,-1.5,-0.9);
  line9->SetLineColor(kYellow);
  line9->SetLineWidth(4);
  line9->Draw("same");
  TLine *line10 = new TLine(-0.5,-0.6,0,-0.1);
  line10->SetLineColor(kYellow);
  line10->SetLineWidth(4);
  line10->Draw("same");
  TLine *line11 = new TLine(-1.5,-1.1,-0.5,-0.6);
  line11->SetLineColor(kYellow);
  line11->SetLineWidth(4);
  line11->Draw("same");
  TLine *line12 = new TLine(-2,-1.3,-1.5,-1.1);
  line12->SetLineColor(kYellow);
  line12->SetLineWidth(4);
  line12->Draw("same");
  TLine *line13 = new TLine(-2,1.1,-2,1.3);
  line13->SetLineColor(kYellow);
  line13->SetLineWidth(4);
  line13->Draw("same");
  TLine *line14 = new TLine(-2,-1.3,-2,-1.1);
  line14->SetLineColor(kYellow);
  line14->SetLineWidth(4);
  line14->Draw("same");
  
  TLine *line15 = new TLine(-0.1,0,0,0.1);
  line15->SetLineColor(kRed);
  line15->SetLineWidth(4);
  line15->Draw("same");
  TLine *line16 = new TLine(-0.1,0,0,-0.1);
  line16->SetLineColor(kRed);
  line16->SetLineWidth(4);
  line16->Draw("same"); 
  TLine *line17 = new TLine(0,0.1,0.05,0.1);
  line17->SetLineColor(kRed);
  line17->SetLineWidth(4);
  line17->Draw("same");
  TLine *line18 = new TLine(0,-0.1,0.05,-0.1);
  line18->SetLineColor(kRed);
  line18->SetLineWidth(4);
  line18->Draw("same"); 
  TLine *line19 = new TLine(0.05,0.1,0.15,0);
  line19->SetLineColor(kRed);
  line19->SetLineWidth(4);
  line19->Draw("same");
  TLine *line20 = new TLine(0.05,-0.1,0.15,0);
  line20->SetLineColor(kRed);
  line20->SetLineWidth(4);
  line20->Draw("same");
  
}

//(TObject*)0
void mppPlotter(TH1F *mpp){
  //mpp->SetTitle("Mass of ZZ pairs with n>=2 reconstructed protons");
  mpp->GetXaxis()->SetTitle("ZZ mass (GeV)");
  mpp->GetYaxis()->SetTitle("# of events");
  mpp->Draw();
}

void yppPlotter(TH1F *ypp){
  ypp->SetTitle("Rapidity of matched ZZ pairs");
  ypp->GetXaxis()->SetTitle("expected ZZ rapidity");
  ypp->GetYaxis()->SetTitle("# of events");
  ypp->Draw();
}

void GeneralPlotter(TH1F *th){
  //th->SetTitle("backward proton multiplicity without pileup");
  th->GetXaxis()->SetTitle("# of possible matches");
  th->GetYaxis()->SetTitle("# of events");
  th->DrawCopy();
}

void PPZZPlotter(){
  TFile *filenoPUP = new TFile("results/2018/PPZZAnalyzer_SR4P/ZZTo4lnoPUP2018.root","READ");
  TFile *filewithPUP = new TFile("results/2018/PPZZAnalyzer_SR4P/ZZTo4lwithPUP2018.root","READ");
  TFile *file = new TFile("results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root","READ");
  TH2F *histxi= (TH2F*)file->Get("th2xi");
  //XiPlotter(histxi);
  TH2F *histmatch = (TH2F*)file->Get("th2");
  //MatchingPlotter(histmatch);
  //TCanvas *c = new TCanvas();
  TH2F *histmatch2 = (TH2F*)file->Get("th2good");
  //MatchingPlotter(histmatch2);
  TH1F *mpphist = (TH1F*)filenoPUP->Get("mZZtwoprotons");
  //mppPlotter(mpphist);
  TH1F *ypphist = (TH1F*)file->Get("ypp");
  //yppPlotter(ypphist);
  TH1F *mZZ = (TH1F*)file->Get("goodmZZ");
  mppPlotter(mZZ);
  TH1F *yZZ = (TH1F*)file->Get("goodyZZ");
  //yppPlotter(yZZ);

  TH1F *n = (TH1F*)filewithPUP->Get("nmatches");
  //GeneralPlotter(n);

  //TH2F *histmatchnoPUP = (TH2F*)filenoPUP->Get("th2");
  //TH2F *histmatchwithPUP = (TH2F*)filewithPUP->Get("th2");
  //histmatchnoPUP->Add(histmatchwithPUP,-1);
  //MatchingPlotter(histmatchnoPUP);
  
}
  
