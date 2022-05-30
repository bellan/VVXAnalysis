#include "TH2F.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPaveLabel.h"

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

double f1(double *x, double *par){
  double xx=x[0];
  return par[0]*log(par[1]*(1-xx));
}

void MatchingPlotter(TH2F *histmatch, double kmin, double kmax){
  histmatch->SetTitle("2D pp-ZZ matching distribution");
  histmatch->GetXaxis()->SetTitle("#bf{m_{match}=1-m_{ZZ}/m_{pp}}");
  histmatch->GetYaxis()->SetTitle("#bf{y_{match}= y_{pp}- y_{ZZ}}");
  gStyle->SetOptStat(0);
  histmatch->DrawCopy("colz");
  
  //TPaveLabel *label = new TPaveLabel(-0.35,1.145,0.5,1.495,"a^{0}_{Z}/#Lambda^{2} = 0.9*10^{-5} GeV^{-2}");
  TPaveLabel *label = new TPaveLabel(-0.35,1.145,0.5,1.495,"qq->ZZ->4l");
  label->SetFillColor(kSpring-4);
  label->SetTextColor(kOrange+10);
  label->SetTextSize(0.5);
  label->Draw("same");
  
  TLine *line1 = new TLine(-0.05,-0.05,-0.05,0.05);
  line1->SetLineColor(kRed);
  line1->SetLineWidth(6);
  line1->Draw("same");
  TLine *line2 = new TLine(-0.05,0.05,0,0.1);
  line2->SetLineColor(kRed);
  line2->SetLineWidth(6);
  line2->Draw("same");
  TLine *line3 = new TLine(-0.05,-0.05,0,-0.1);
  line3->SetLineColor(kRed);
  line3->SetLineWidth(6);
  line3->Draw("same");
  TLine *line4 = new TLine(0,0.1,0.05,0.1);
  line4->SetLineColor(kRed);
  line4->SetLineWidth(6);
  line4->Draw("same");
  TLine *line5 = new TLine(0,-0.1,0.05,-0.1);
  line5->SetLineColor(kRed);
  line5->SetLineWidth(6);
  line5->Draw("same");
  TLine *line6 = new TLine(0.05,0.1,0.1,0.05);
  line6->SetLineColor(kRed);
  line6->SetLineWidth(6);
  line6->Draw("same");
  TLine *line7 = new TLine(0.05,-0.1,0.1,-0.05);
  line7->SetLineColor(kRed);
  line7->SetLineWidth(6);
  line7->Draw("same");
  TLine *line8 = new TLine(0.1,0.05,0.1,-0.05);
  line8->SetLineColor(kRed);
  line8->SetLineWidth(6);
  line8->Draw("same");
  

  TF1 *f = new TF1("myfunc",f1,-2.01,-0.04,2);
  f->SetLineWidth(6);
  f->SetLineColor(kOrange-2);

  double y0,y1,y2,y3,y4,y5,y6,y7;
  
  f->SetParameter(0,1);
  f->SetParameter(1,kmax);
  y0=f->Eval(-0.05);
  y4=f->Eval(-2);
  f->DrawCopy("same");
  f->SetParameter(1,kmin);
  y1=f->Eval(-0.05);
  y5=f->Eval(-2);
  f->DrawCopy("same");
  f->SetParameter(0,-1);
  y2=f->Eval(-0.05);
  y6=f->Eval(-2);
  f->DrawCopy("same");
  f->SetParameter(1,kmax);
  y3=f->Eval(-0.05);
  y7=f->Eval(-2);
  f->DrawCopy("same");

  
  TLine *line9 = new TLine(-2,y4,-2,y5);
  line9->SetLineColor(kOrange-2);
  line9->SetLineWidth(6);
  line9->Draw("same");
  TLine *line10 = new TLine(-2,y6,-2,y7);
  line10->SetLineColor(kOrange-2);
  line10->SetLineWidth(6);
  line10->Draw("same");
  TLine *line11 = new TLine(-0.05,y0,-0.05,y1);
  line11->SetLineColor(kOrange-2);
  line11->SetLineWidth(6);
  line11->Draw("same");
  TLine *line12 = new TLine(-0.05,y2,-0.05,y3);
  line12->SetLineColor(kOrange-2);
  line12->SetLineWidth(6);
  line12->Draw("same");

  line1->Draw("same");
  
}

//(TObject*)0
void mppPlotter(TH1F *mpp){
  mpp->SetTitle("ZZ invariant mass");
  gStyle->SetOptStat(0);
  mpp->GetXaxis()->SetTitle("mass (GeV)");
  mpp->GetYaxis()->SetTitle("# of events");
  TPaveLabel *label = new TPaveLabel(-0.35,1.145,0.5,1.495,"qq->ZZ->4l");
  label->SetFillColor(kSpring-4);
  label->SetTextColor(kOrange+10);
  label->SetTextSize(0.5);
  mpp->DrawCopy("same");
  label->Draw("same");
}

void yppPlotter(TH1F *ypp){
  gStyle->SetOptStat(0);
  ypp->SetTitle("Matched ZZ rapidity");
  ypp->GetXaxis()->SetTitle("y");
  ypp->GetYaxis()->SetTitle("# of events");
  ypp->Draw();
  TLegend *l = new TLegend(0.6,0.7,0.95,0.8);
  l->AddEntry((TObject*)0,"a^{c}_{Z}/#Lambda^{2} = 5*10^{-5} GeV^{-2}","f");
  l->Draw("same");
}

void GeneralPlotter(TH1F *th){
  th->SetTitle("backwards proton multiplicity without pileup");
  gStyle->SetOptStat(0);
  th->GetXaxis()->SetTitle("multiplicity");
  th->GetYaxis()->SetTitle("# of events");
  th->DrawCopy();
  TLegend *l = new TLegend(0.6,0.7,0.95,0.8);
  l->AddEntry((TObject*)0,"a^{0}_{Z}/#Lambda^{2} = 5*10^{-5} GeV^{-2}","f");
  l->Draw("same");
}

void PPZZPlotter(){
  
  TFile *filenoPUP = new TFile("results/2018/PPZZAnalyzer_SR4P/ZZTo4lnoPUP2018.root","READ");
  TFile *filewithPUP = new TFile("results/2018/PPZZAnalyzer_SR4P/ZZTo4lwithPUP2018.root","READ");
  TFile *file = new TFile("results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root","READ");
  //TFile *file = new TFile("results/2018/PPZZAnalyzer_SR4P/ZZTo4la0zD2.root","READ");

  double theWeight= ((TH1F*)file->Get("weight_full"))->GetMean();
  
  TH2F *histxi= (TH2F*)file->Get("th2xi");
  //XiPlotter(histxi);
  TH2F *histmatch = (TH2F*)file->Get("th2goodC");
  MatchingPlotter(histmatch,0.95,1.15);
  //TCanvas *c = new TCanvas();
  TH2F *histmatch2 = (TH2F*)file->Get("th2good");
  //MatchingPlotter(histmatch2);
  TH1F *mpphist = (TH1F*)file->Get("mZZ");
  //mppPlotter(mpphist);
  TH1F *ypphist = (TH1F*)file->Get("ypp");
  //yppPlotter(ypphist);
  TH1F *mZZ = (TH1F*)file->Get("goodmZZ");
  //mppPlotter(mZZ);
  //TCanvas *c = new TCanvas();
  TH1F *yZZ = (TH1F*)file->Get("goodyZZ");
  //yppPlotter(yZZ);

  TH1F *n = (TH1F*)file->Get("backprotonmult");
  //GeneralPlotter(n);
  /*
  TH2F *histmatchonlyPUP = (TH2F*)file->Get("th2good");
  TH2F *histmatchwithPUP = (TH2F*)file->Get("th2goodC");
  histmatchwithPUP->Add(histmatchonlyPUP,-1);
  for(int i=1;i<=3600;i++){
    if(histmatchwithPUP->GetBinContent(i)<0) histmatchwithPUP->SetBinContent(i,0);}
  MatchingPlotter(histmatchwithPUP,0.95,1.15);
  */
  
  TH1 *counteromicron = (TH2F*)file->Get("counteromicron2");
  TH1 *counterdelta = (TH2F*)file->Get("counterdelta2");
  cout<<endl<<"Algorithm A (check if delta region, if not min distance):"<<endl;
  if(counterdelta){cout<<"# of events in the delta signal region: "<<counterdelta->GetEntries()*theWeight<<endl;
    cout<<"# of events in the omicron signal region: "<<(counteromicron->GetEntries()-counterdelta->GetEntries())*theWeight<<endl<<endl;}
  else{cout<<"# of events in the delta signal region: 0"<<endl;
    cout<<"# of events in the omicron signal region: "<<(counteromicron->GetEntries())*theWeight<<endl<<endl;}
  
  TH1 *counteromicronxi = (TH2F*)file->Get("counteromicronxi");
  TH1 *counterdeltaxi = (TH2F*)file->Get("counterdeltaxi");
  cout<<"Algorithm B (max xi):"<<endl;
  cout<<"# of events in the delta signal region: "<<counterdeltaxi->GetEntries()*theWeight<<endl;
  cout<<"# of events in the omicron signal region: "<<(counteromicronxi->GetEntries()-counterdeltaxi->GetEntries())*theWeight<<endl<<endl;
   
  TH1 *counteromicronC = (TH2F*)file->Get("counteromicronC");
  cout<<"Algorithm C (check if delta region, if not max xi):"<<endl;
  cout<<"# of events in the delta signal region: "<<counterdelta->GetEntries()*theWeight<<endl;
  cout<<"# of events in the omicron signal region: "<<counteromicronC->GetEntries()*theWeight<<endl<<endl;
  
  TH1 *counteromicron4e = (TH2F*)file->Get("counteromicron4e");
  TH1 *counterdelta4e = (TH2F*)file->Get("counterdelta4e");
  TH1 *counteromicron4mu = (TH2F*)file->Get("counteromicron4mu");
  TH1 *counterdelta4mu = (TH2F*)file->Get("counterdelta4mu");
  TH1 *counteromicron2e2mu = (TH2F*)file->Get("counteromicron2e2mu");
  TH1 *counterdelta2e2mu = (TH2F*)file->Get("counterdelta2e2mu");
  cout<<"4mu:"<<endl;
  cout<<"# of events in the delta signal region: "<<counterdelta4mu->GetEntries()*theWeight<<endl;
  cout<<"# of events in the omicron signal region: "<<counteromicron4mu->GetEntries()*theWeight<<endl<<endl;
  cout<<"4e:"<<endl;
  cout<<"# of events in the delta signal region: "<<counterdelta4e->GetEntries()*theWeight<<endl;
  cout<<"# of events in the omicron signal region: "<<counteromicron4e->GetEntries()*theWeight<<endl<<endl;
  cout<<"2e2mu:"<<endl;
  cout<<"# of events in the delta signal region: "<<counterdelta2e2mu->GetEntries()*theWeight<<endl;
  cout<<"# of events in the omicron signal region: "<<counteromicron2e2mu->GetEntries()*theWeight<<endl<<endl;
  
}
  
