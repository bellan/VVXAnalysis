#include "TH2F.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPaveLabel.h"
#include "TString.h"
#include "../../Commons/interface/Constants.h"

double f1(double *x, double *par);
void MatchingPlotter(TH2F *histmatch, double kmin, double kmax);
void mppPlotter(TH1F *mpp);
void yppPlotter(TH1F *ypp);
void GeneralPlotter(TH1F *th);


void PPZZPlotter2(TString sample){
  
  TFile *filenoPUP = new TFile("results/2018/PPZZAnalyzer_SR4P/ZZTo4lnoPUP2018.root","READ");
  TFile *filewithPUP = new TFile("results/2018/PPZZAnalyzer_SR4P/ZZTo4lwithPUP2018.root","READ");
  TFile *file = new TFile("results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root","READ");
  //TFile *file = new TFile("results/2018/PPZZAnalyzer_SR4P/ZZTo4la0zD2.root","READ");

  double theWeight;
  
  if(sample=="a0z") theWeight=0.003783742; //a0z signal
  if(sample=="qqZZ") theWeight=0.000570889; //qqZZ
  if(sample=="ggZZ") theWeight=0.000242160969322; //ggZZ

  
  TH2F *histxi= (TH2F*)file->Get("th2xi");
  //XiPlotter(histxi);
  TH2F *histmatch = (TH2F*)file->Get("th2goodC");
  MatchingPlotter(histmatch,0.95,1.15);
  /*
  TH1F *mpphist = (TH1F*)file->Get("mpp");
  TH1F *mZZ = (TH1F*)file->Get("mZZ");
  mpphist->SetTitle("m_{pp} vs m_{ZZ}");
  mpphist->GetXaxis()->SetTitle("m (GeV)");
  TCanvas *c1=new TCanvas();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  mpphist->Draw("hist");
  mZZ->SetLineColor(kRed);
  mZZ->Draw("samehist");
  TLegend *l = new TLegend(0.7,0.8,0.9,0.9);
  l->AddEntry(mpphist,"m_{pp}","l");
  l->AddEntry(mZZ,"m_{ZZ}","l");
  l->Draw("same");
  */
  TH1F *ypphist = (TH1F*)file->Get("ypp");
  //yppPlotter(ypphist);
  //TH1F *mZZ = (TH1F*)file->Get("goodmZZ");
  //mppPlotter(mZZ);
  TH1F *yZZ = (TH1F*)file->Get("goodyZZ");
  //yppPlotter(yZZ);
  
  TH1F *aco = (TH1F*)file->Get("acoplanarity");
  //GeneralPlotter(aco);
  
  TH1 *counteromicron = (TH2F*)file->Get("counteromicron2");
  TH1 *counterdelta = (TH2F*)file->Get("counterdelta2");
  cout<<endl<<"Algorithm A (check if delta region, if not min distance):"<<endl;
  if(counterdelta){cout<<"# of events in the delta signal region: "<<counterdelta->GetEntries()*theWeight<<endl;
    cout<<"# of events in the omicron signal region: "<<(counteromicron->GetEntries()-counterdelta->GetEntries())*theWeight<<endl<<endl;}
  else if(counteromicron) {cout<<"# of events in the delta signal region: 0"<<endl;
    cout<<"# of events in the omicron signal region: "<<(counteromicron->GetEntries())*theWeight<<endl<<endl;}
  else {cout<<"# of events in the delta signal region: 0"<<endl;
    cout<<"# of events in the omicron signal region: 0"<<endl<<endl;}
  
  TH1 *counteromicronxi = (TH2F*)file->Get("counteromicronxi");
  TH1 *counterdeltaxi = (TH2F*)file->Get("counterdeltaxi");
  cout<<"Algorithm B (max xi):"<<endl;
  if(counterdeltaxi) cout<<"# of events in the delta signal region: "<<counterdeltaxi->GetEntries()*theWeight<<endl;
  else cout<<"# of events in the delta signal region: 0"<<endl;
  if(counteromicronxi) cout<<"# of events in the omicron signal region: "<<(counteromicronxi->GetEntries()-counterdeltaxi->GetEntries())*theWeight<<endl<<endl;
  else cout<<"# of events in the omicron signal region: 0"<<endl<<endl;
  
  TH1 *counteromicronC = (TH2F*)file->Get("counteromicronC");
  cout<<"Algorithm C (check if delta region, if not max xi):"<<endl;
  if(counterdelta) cout<<"# of events in the delta signal region: "<<counterdelta->GetEntries()*theWeight<<endl;
  else cout<<"# of events in the delta signal region: 0"<<endl;
  if(counteromicronC) cout<<"# of events in the omicron signal region: "<<counteromicronC->GetEntries()*theWeight<<endl<<endl;
  else cout<<"# of events in the omicron signal region: 0"<<endl<<endl;
  /*
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
  */
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
  
  TPaveLabel *label = new TPaveLabel(-0.35,1.145,0.5,1.495,"a^{0}_{Z}/#Lambda^{2} = 0.9*10^{-5} GeV^{-2}");
  //TPaveLabel *label = new TPaveLabel(-0.35,1.145,0.5,1.495,"gg->ZZ->4l");
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
  th->SetTitle("ZZ acoplanarity");
  gStyle->SetOptStat(0);
  th->GetXaxis()->SetTitle("acoplanarity");
  th->GetYaxis()->SetTitle("# of events");
  th->DrawCopy();
  TPaveLabel *label = new TPaveLabel(0.7,4,1,4.45,"gg->ZZ->4l");
  label->SetFillColor(kSpring-4);
  label->SetTextColor(kOrange+10);
  label->SetTextSize(0.5);
  label->Draw("same");
}




void sbvsnprotons(){
  
  float sb[]={242.9907,241.1025,242.7663,227.4062};
  float sberror[]={16.2353706845161,16.1092100895071,16.5180270078873,15.9215044414918};
  float nprotons[]={5.,4.,3.,2.};
  float errnprotons[]={0.,0.,0.,0.};

  TGraphErrors *tge = new TGraphErrors(4,nprotons,sb,errnprotons,sberror);
  tge->SetMarkerSize(5);
  tge->SetMarkerColor(kRed);
  tge->Draw("AP");

}

void sbvskmin(){
  
  float sb[]={42.031,43.099,46.658,47.670,49.327,49.722,48.832,46.649,41.206};
  float sberror[]={0.951,0.994,1.093,1.128,4.666,4.674,4.520,4.438,1.132};
  float kmin[]={0.80,0.85,0.90,0.93,0.94,0.95,0.96,0.97,1.00};
  float kminerror[]={0.,0.,0.,0.,0.,0.,0.,0.,0.};

  TGraphErrors *tge = new TGraphErrors(9,kmin,sb,kminerror,sberror);
  tge->SetTitle("S/#sqrt{B} vs k_{min}");
  tge->GetXaxis()->SetTitle("k_{min}");
  tge->GetYaxis()->SetTitle("S/#sqrt{B}");
  tge->SetMarkerSize(0.8);
  tge->SetMarkerColor(kRed);
  tge->SetMarkerStyle(21);
  tge->SetLineColor(kBlack);
  tge->SetLineWidth(1);
  tge->Draw("AP");

  TLine *line = new TLine(0.95,39,0.95,44);
  line->SetLineColor(kBlue);
  line->SetLineWidth(3);
  line->Draw("same");
  TLine *line2 = new TLine(0.95,54.7,0.95,55.5);
  line2->SetLineColor(kBlue);
  line2->SetLineWidth(3);
  line2->Draw("same");

}

void sbvskmax(){
  
  float sb[]={53.603,54.328,53.241,53.305,51.748,50.489,49.722,45.073,43.157,42.289};
  float sberror[]={1.362,5.403,5.057,5.097,4.755,4.499,4.399,1.064,1.009,0.981};
  float kmax[]={1.09,1.1,1.11,1.12,1.13,1.14,1.15,1.20,1.25,1.30};
  float kmaxerror[]={0.,0.,0.,0.,0.,0.,0.,0.,0.};

  TGraphErrors *tge = new TGraphErrors(10,kmax,sb,kmaxerror,sberror);
  tge->SetTitle("S/#sqrt{B} vs k_{max}");
  tge->GetXaxis()->SetTitle("k_{max}");
  tge->GetYaxis()->SetTitle("S/#sqrt{B}");
  tge->SetMarkerSize(0.8);
  tge->SetMarkerColor(kRed);
  tge->SetMarkerStyle(21);
  tge->SetLineColor(kBlack);
  tge->SetLineWidth(1);
  tge->Draw("AP");

  TLine *line = new TLine(1.1,40,1.1,48);
  line->SetLineColor(kBlue);
  line->SetLineWidth(3);
  line->Draw("same");
  TLine *line2 = new TLine(1.1,60,1.1,61.5);
  line2->SetLineColor(kBlue);
  line2->SetLineWidth(3);
  line2->Draw("same");
  
}
