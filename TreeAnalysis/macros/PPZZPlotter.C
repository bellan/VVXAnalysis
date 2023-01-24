#include "TH2F.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPaveLabel.h"
#include "TString.h"
#include "../../Commons/interface/Constants.h"
#include "UnfoldingMacros/CMS_lumi.h"
#include "UnfoldingMacros/CMS_lumi.C"
#include "TFile.h"

double f1(double *x, double *par);
void MatchingPlotter(TH2F *histmatch, double kmin, double kmax);
void mppPlotter(TH1F *mpp);
void yppPlotter(TH1F *ypp);
void GeneralPlotter(TH1F *th);


void PPZZPlotter(){
  
  TFile *filenoPUP = new TFile("results/2017/PPZZAnalyzer_SR4P/ZZTo4lnoPUP.root","READ");
  TFile *file = new TFile("results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root","READ");
  TFile *fileBG = new TFile("results/2018/PPZZAnalyzer_SR4P/BG.root","READ");
  TFile *fileqqZZ = new TFile("results/2018/PPZZAnalyzer_SR4P/qqZZ.root","READ");
  TFile *fileggZZ = new TFile("results/2018/PPZZAnalyzer_SR4P/ggZZ.root","READ");
  TFile *fileVBS = new TFile("results/2018/PPZZAnalyzer_SR4P/VBS.root","READ");
  TFile *filea0z = new TFile("results/2018/PPZZAnalyzer_SR4P/a0z.root","READ");
  TFile *fileaCz = new TFile("results/2018/PPZZAnalyzer_SR4P/aCz.root","READ");
  TFile *filemix = new TFile("results/2018/PPZZAnalyzer_SR4P/mix.root","READ");
  TFile *filedata = new TFile("results/2018/PPZZAnalyzer_SR4P/DATA.root","READ");
  
  TH2F *histxi= (TH2F*)file->Get("th2xi");
  //XiPlotter(histxi);
  TH2F *histmatch = (TH2F*)file->Get("th2goodC");
  TH2F *histmatchZeta = (TH2F*)filea0z->Get("th2goodCZeta");
  //histmatch->GetXaxis()->SetRangeUser(-0.5,0.5);
  //histmatch->GetYaxis()->SetRangeUser(-0.5,0.5);
  /*
  TFile *f= new TFile;
  f=TFile::Open("/home/giovanni/PPtoPPWWjets/PPtoPPWWjets/analysismacros/NewMixingDistributions2018.root","r");
  TH1D *hn45 = (TH1D *)f->Get("mtpl_multi_S45");
  TH1D *hn56 = (TH1D *)f->Get("mtpl_multi_S56");
  TH1D *hxi45 = (TH1D *)f->Get("xi_multi_S45");
  TH1D *hxi56 = (TH1D *)f->Get("xi_multi_S56");
  gStyle->SetOptStat(0);
  hxi45->SetTitle("");
  hxi45->GetXaxis()->SetTitle("#xi");
  hxi45->GetYaxis()->SetTitle("# of events");
  hxi56->GetXaxis()->SetTitle("#xi");
  hxi56->GetYaxis()->SetTitle("# of events");
  hxi56->SetTitle("");
  hxi45->DrawCopy();
  hxi56->SetLineColor(kRed);
  hxi56->DrawCopy("same");

  TLegend *legend = new TLegend(0.48,0.7,0.9,0.9);
  hn45->SetLineColor(kBlack);
  legend->AddEntry(hn45,"sector 45 (forward)","l");
  hn56->SetLineColor(kRed);
  TH1F *htest = new TH1F();
  htest->SetLineColor(kRed);
  legend->AddEntry(htest,"sector 56 (backward)","l");
  legend->Draw("same");

  f->Close();
  delete f;
  */
  
  TCanvas *c_p = new TCanvas ("c_p","c_p");
  string p_png = "images/matching2017noPUP.png";
  MatchingPlotter(histmatchZeta,0.8,1.15);
  CMS_lumi(c_p,18,5);
  c_p->Print(p_png.c_str());
  //cout<<"Numero di punti sul grafico: "<<histmatch->GetEntries()<<endl;
  //cout<<"Numero di punti sul grafico: "<<histmatchZeta->GetEntries()<<endl;
  
  /*
  TH2F *histmatchqqZZ = (TH2F*)fileqqZZ->Get("th2goodC");
  TH2F *histmatchggZZ = (TH2F*)fileggZZ->Get("th2goodC");
  TH2F *histmatchVBS = (TH2F*)fileVBS->Get("th2goodC");
  histmatchqqZZ->Add(histmatchggZZ,1);
  histmatchqqZZ->Add(histmatchVBS,1);
  TCanvas *c_p = new TCanvas ("c_p","c_p");
  string p_png = "prova.png";
  MatchingPlotter(histmatchqqZZ,0.95,1.10);
  CMS_lumi(c_p,18,5);
  c_p->Print(p_png.c_str());
  */
  /*
  TH1F *acosignal = (TH1F*) file->Get("acoplanarity");
  TH1F *acoqqZZ = (TH1F*)fileqqZZ->Get("acoplanarity");
  TH1F *acoggZZ = (TH1F*)fileggZZ->Get("acoplanarity");
  TH1F *acoVBS = (TH1F*)fileVBS->Get("acoplanarity");
  acoqqZZ->Add(acoggZZ,1);
  acoqqZZ->Add(acoVBS,1);
  acoqqZZ->SetTitle("ZZ acoplanarity");
  TCanvas *c_p = new TCanvas ("c_p","c_p");
  string p_png = "images/acoplanarityBG.png";
  GeneralPlotter(acoqqZZ);
  CMS_lumi(c_p,18,5);
  c_p->Print(p_png.c_str());
  */
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
  
  TCanvas *c_p = new TCanvas ("c_p","c_p");
  string p_png = "images/backmultcomparison.png";
  TH1F *forwnoPUP = (TH1F*)filenoPUP->Get("forwprotonmult");
  TH1F *forw = (TH1F*)file->Get("forwprotonmult");
  TH1F *backnoPUP = (TH1F*)filenoPUP->Get("backprotonmult");
  TH1F *back = (TH1F*)file->Get("backprotonmult");
  backnoPUP->SetLineColor(kRed);
  GeneralPlotter(backnoPUP);
  //forw->SetLineColor(kRed);
  back->Draw("same");
  CMS_lumi(c_p,18,5);
  c_p->Print(p_png.c_str());
*/
  /*
  TCanvas *c_p = new TCanvas ("c_p","c_p");
  string p_png = "images/xiback.png";
  TH1F *xiforwnoPUP = (TH1F*)filenoPUP->Get("xiforw");
  TH1F *xiforw = (TH1F*)file->Get("xiforw");
  TH1F *xibacknoPUP = (TH1F*)filenoPUP->Get("xiback");
  TH1F *xiback = (TH1F*)file->Get("xiback");

  xiforw->Add(xiforwnoPUP,-1);
  xiback->Add(xibacknoPUP,-1);
  
  xiforwnoPUP->Scale(1/xiforwnoPUP->Integral());
  xiforw->Scale(1/xiforw->Integral());
  xibacknoPUP->Scale(1/xibacknoPUP->Integral());
  xiback->Scale(1/xiback->Integral());
  xiforwnoPUP->SetLineColor(kRed);
  xiforw->SetLineColor(kBlack);
  xibacknoPUP->SetLineColor(kRed);
  xiback->SetLineColor(kBlack);
  xiforw->GetYaxis()->SetRangeUser(0,0.15);
  xiback->GetYaxis()->SetRangeUser(0,0.15);
  GeneralPlotter(xiback);
  //forw->SetLineColor(kRed);
  xibacknoPUP->Draw("same");
  CMS_lumi(c_p,18,5);
  c_p->Print(p_png.c_str());
  */
  /*
  double forwc=0,forwnoPUPc=0,backc=0,backnoPUPc=0;
  for(int i=1;i<7;i++){
    forwnoPUPc+=forwnoPUP->GetBinContent(i)*(i-1);
    forwc+=forw->GetBinContent(i)*(i-1);
    backnoPUPc+=backnoPUP->GetBinContent(i)*(i-1);
    backc+=back->GetBinContent(i)*(i-1);
  }
  
  cout<<"FORW: "<<forwnoPUPc<<" / "<<forwc<<endl;
  cout<<"BACK: "<<backnoPUPc<<" / "<<backc<<endl;
  
  //TH1F *ypphist = (TH1F*)file->Get("ypp");
  //yppPlotter(ypphist);
  TH1F *mZZ = (TH1F*)file->Get("mZZ");
  //mppPlotter(mZZ);
  TH1F *yZZ = (TH1F*)file->Get("goodyZZ");
  //yppPlotter(yZZ);
  */
  /*
  TCanvas *c_p = new TCanvas ("c_p","c_p");
  string p_png = "images/mZZcomparison2017.png";
  TH1F *mZZ1 = (TH1F*)file->Get("mZZ");
  TH1F *mZZ2 = (TH1F*)file->Get("mZZoneproton");
  TH1F *mZZ3 = (TH1F*)file->Get("mZZtwoprotons");
  mZZ1->SetLineColor(kBlack);
  mZZ2->SetLineColor(kBlue);
  mZZ3->SetLineColor(kRed);
  mZZ1->SetTitle("");
  mZZ1->GetYaxis()->SetRangeUser(0,13);
  mppPlotter(mZZ1);
  //gPad->SetLogy();
  mZZ2->SetTitle("");
  mZZ3->SetTitle("");
  mZZ2->Draw("same");
  mZZ3->Draw("same");
  CMS_lumi(c_p,18,5);
  c_p->Print(p_png.c_str());
  */
  
  TH1 *counterdelta = (TH2F*)file->Get("counterdelta2");
  /*
  TH1 *counteromicron = (TH2F*)file->Get("counteromicron2");
  cout<<endl<<"Algorithm A (check if delta region, if not min distance):"<<endl;
  if(counterdelta){cout<<"# of events in the delta signal region: "<<counterdelta->GetBinContent(1)<<endl;
    cout<<"# of events in the omicron signal region: "<<counteromicron->GetBinContent(1)-counterdelta->GetBinContent(1)<<endl<<endl;}
  else if(counteromicron) {cout<<"# of events in the delta signal region: 0"<<endl;
    cout<<"# of events in the omicron signal region: "<<counteromicron->GetBinContent(1)<<endl<<endl;}
  else {cout<<"# of events in the delta signal region: 0"<<endl;
    cout<<"# of events in the omicron signal region: 0"<<endl<<endl;}
  
  TH1 *counteromicronxi = (TH2F*)file->Get("counteromicronxi");
  TH1 *counterdeltaxi = (TH2F*)file->Get("counterdeltaxi");
  cout<<"Algorithm B (max xi):"<<endl;
  if(counterdeltaxi) cout<<"# of events in the delta signal region: "<<counterdeltaxi->GetBinContent(1)<<endl;
  else cout<<"# of events in the delta signal region: 0"<<endl;
  if(counteromicronxi) cout<<"# of events in the omicron signal region: "<<counteromicronxi->GetBinContent(1)-counterdeltaxi->GetBinContent(1)<<endl<<endl;
  else cout<<"# of events in the omicron signal region: 0"<<endl<<endl;
  */
  TH1 *counteromicronC = (TH2F*)file->Get("counteromicronC");
  TH1 *counterzeta = (TH2F*)file->Get("counterzeta");
  cout<<"Algorithm C (check if delta region, if not max xi):"<<endl;
  if(counterdelta) cout<<"# of events in the delta signal region: "<<counterdelta->GetBinContent(1)<<" +- "<<counterdelta->GetBinError(1)<<endl;
  else cout<<"# of events in the delta signal region: 0"<<endl;
  if(counteromicronC) cout<<"# of events in the omicron signal region: "<<counteromicronC->GetBinContent(1)<<" +- "<<counteromicronC->GetBinError(1)<<endl;
  else cout<<"# of events in the omicron signal region: 0"<<endl;
  if(counterzeta) cout<<"# of events in the zeta signal region: "<<counterzeta->GetBinContent(1)<<" +- "<<counterzeta->GetBinError(1)<<endl;
  else cout<<"# of events in the zeta signal region: 0"<<endl<<endl;
  /*
  TH1 *counteromicron4e = (TH2F*)file->Get("counteromicron4e");
  TH1 *counterdelta4e = (TH2F*)file->Get("counterdelta4e");
  TH1 *counteromicron4mu = (TH2F*)file->Get("counteromicron4mu");
  TH1 *counterdelta4mu = (TH2F*)file->Get("counterdelta4mu");
  TH1 *counteromicron2e2mu = (TH2F*)file->Get("counteromicron2e2mu");
  TH1 *counterdelta2e2mu = (TH2F*)file->Get("counterdelta2e2mu");
  cout<<"4mu:"<<endl;
  cout<<"# of events in the delta signal region: "<<counterdelta4mu->GetBinContent(1)<<endl;
  cout<<"# of events in the omicron signal region: "<<counteromicron4mu->GetBinContent(1)<<endl<<endl;
  cout<<"4e:"<<endl;
  cout<<"# of events in the delta signal region: "<<counterdelta4e->GetBinContent(1)<<endl;
  cout<<"# of events in the omicron signal region: "<<counteromicron4e->GetBinContent(1)<<endl<<endl;
  cout<<"2e2mu:"<<endl;
  cout<<"# of events in the delta signal region: "<<counterdelta2e2mu->GetBinContent(1)<<endl;
  cout<<"# of events in the omicron signal region: "<<counteromicron2e2mu->GetBinContent(1)<<endl<<endl;
  */
}


double f1(double *x, double *par){
  double xx=x[0];
  return par[0]*log(par[1]*(1-xx));
}

void MatchingPlotter(TH2F *histmatch, double kmin, double kmax){
  histmatch->SetTitle("");
  histmatch->GetXaxis()->SetTitle("#bf{m_{match}=1-m_{ZZ}/m_{pp}}");
  histmatch->GetYaxis()->SetTitle("#bf{y_{match}= y_{pp}- y_{ZZ}}");
  gStyle->SetOptStat(0);
  histmatch->DrawCopy("colz");
  
  TPaveLabel *label = new TPaveLabel(-0.35,1.145,1,1.495,"a^{0}_{Z}/#Lambda^{2} = 0.9*10^{-5} GeV^{-2}");
  //TPaveLabel *label = new TPaveLabel(0.1,0.378,0.5,0.498,"a^{0}_{Z}/#Lambda^{2} = 0.9*10^{-5} GeV^{-2}");
  //TPaveLabel *label = new TPaveLabel(-2.4975,1.145,-1.5,1.495,"background");
  label->SetFillColor(kSpring-4);
  label->SetTextColor(kOrange+10);
  label->SetTextSize(0.5);
  label->Draw("same");
  
  TLine *line1 = new TLine(-0.05,-0.05,-0.05,0.05);
  line1->SetLineColor(kRed); //kRed
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
  f->SetLineColor(kOrange-2); //kOrange-2

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
  mpp->SetTitle("");
  gStyle->SetOptStat(0);
  mpp->GetXaxis()->SetTitle("mass (GeV)");
  mpp->GetYaxis()->SetTitle("# of events");
  mpp->GetYaxis()->SetRangeUser(0,16.25);
  mpp->DrawCopy("");
  
  TPaveLabel *label = new TPaveLabel(4570,5.69,8000,6.49,"a^{0}_{Z}/#Lambda^{2} = 0.9*10^{-5} GeV^{-2}");
  label->SetFillColor(kSpring-4);
  label->SetTextColor(kOrange+10);
  label->SetTextSize(0.5);
  label->Draw("same");
  
  TLegend *legend = new TLegend(0.54,0.60,0.9,0.8);
  TH1F *th2 = new TH1F(*mpp);
  TH1F *th3 = new TH1F(*mpp);
  th2->SetLineColor(kBlue);
  th3->SetLineColor(kRed);
  //legend->AddEntry(mpp,"a^{0}_{Z}/#Lambda^{2} = 0.9*10^{-5} GeV^{-2}","l");
  //legend->AddEntry(th2,"a^{C}_{Z}/#Lambda^{2} = 4.0*10^{-5} GeV^{-2}","l");
  //legend->AddEntry(th3,"both","l");
  legend->AddEntry(mpp,"all ZZ pairs","l");
  legend->AddEntry(th2,">=1 proton","l");
  legend->AddEntry(th3,">=2 protons","l");
  legend->Draw("same");
}

void yppPlotter(TH1F *ypp){
  gStyle->SetOptStat(0);
  ypp->SetTitle("");
  ypp->GetXaxis()->SetTitle("y");
  ypp->GetYaxis()->SetTitle("# of events");
  ypp->Draw();
  
  TPaveLabel *label = new TPaveLabel(0.81,11.31,5,12.97,"a^{0}_{Z}/#Lambda^{2} = 0.9*10^{-5} GeV^{-2}");
  label->SetFillColor(kSpring-4);
  label->SetTextColor(kOrange+10);
  label->SetTextSize(0.5);
  label->Draw("same");
  
  TLegend *legend = new TLegend(0.565,0.65,0.9,0.8);
  TH1F *th2 = new TH1F(*ypp);
  TH1F *th3 = new TH1F(*ypp);
  th2->SetLineColor(kBlue);
  th3->SetLineColor(kRed);
  legend->AddEntry(ypp,"all ZZ pairs","l");
  legend->AddEntry(th2,">=1 proton","l");
  legend->AddEntry(th3,">=2 protons","l");
  legend->Draw("same");
}

void GeneralPlotter(TH1F *th){
  th->SetTitle("");
  gStyle->SetOptStat(0);
  th->GetXaxis()->SetTitle("acoplanarity");
  th->GetYaxis()->SetTitle("# of events");
  th->DrawCopy("");
  //TPaveLabel *label = new TPaveLabel(0.00275,36.5,0.005,42.8,"a^{0}_{Z}/#Lambda^{2} = 0.9*10^{-5} GeV^{-2}");
  TPaveLabel *label = new TPaveLabel(0.55,275,1,321,"background");
  label->SetFillColor(kSpring-4);
  label->SetTextColor(kOrange+10);
  label->SetTextSize(0.5);
  label->Draw("same");
  /*
  TLegend *legend = new TLegend(0.54,0.65,0.9,0.784);
  TH1F *th2 = new TH1F(*th);
  th2->SetLineColor(kRed);
  th->SetLineColor(kBlack);
  legend->AddEntry(th2,"no pileup","l");
  legend->AddEntry(th,"pileup","l");
  legend->Draw("same");
  */
}




void sbvsnprotons(){
  
  float sb[]={36.6,35.4,34.7,31.7};
  float sberror[]={2.3,2.2,2.2,2.1};
  float nprotons[]={5.,4.,3.,2.};
  float errnprotons[]={0.,0.,0.,0.};

  TGraphErrors *tge = new TGraphErrors(4,nprotons,sb,errnprotons,sberror);
  tge->SetMarkerSize(0.8);
  tge->SetMarkerColor(kRed);
  tge->SetMarkerStyle(21);
  tge->SetTitle("S/#sqrt{B} vs max protons per arm");
  tge->GetXaxis()->SetTitle("max protons per arm");
  tge->GetYaxis()->SetTitle("S/#sqrt{B}");
  tge->Draw("AP");

  TPaveLabel *label = new TPaveLabel(1.8,38,3,40,"a^{0}_{Z}/#Lambda^{2} = 0.9*10^{-5} GeV^{-2}");
  label->SetFillColor(kSpring-4);
  label->SetTextColor(kOrange+10);
  label->SetTextSize(0.5);
  label->Draw("same");

}

void sbvskmin(){

  TCanvas *c_p = new TCanvas ("c_p","c_p");
  string p_png = "images/kminselection.png";
  
  float sb1[]={35.3,36.1,39.3,40.5,41.9,42.1,41.4,39.8,35.2};
  float sberror1[]={2.2,2.3,2.8,3.0,3.3,3.4,3.3,3.2,3.1};

  float sb2[]={52.0,53.4,58.1,59.4,62.1,62.0,60.9,58.9,50.1};
  float sberror2[]={3.2,3.4,4.1,4.4,4.9,4.9,4.9,4.8,4.4};

  float sb3[]={172.6,176.5,191.9,197.5,205.6,206.1,199.2,192.2,167.3};
  float sberror3[]={10.5,11.4,13.7,14.7,16.1,16.4,16.0,15.5,14.7};

  
  float kmin[]={0.80,0.85,0.90,0.93,0.94,0.95,0.96,0.97,1.00};
  float kminerror[]={0.,0.,0.,0.,0.,0.,0.,0.,0.};

  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("#omicron S/#sqrt{B} vs k_{min}");

  TGraphErrors *tge = new TGraphErrors(10,kmin,sb1,kminerror,sberror1);
  tge->SetTitle("");
  tge->SetMarkerSize(0.8);
  tge->SetMarkerColor(kRed);
  tge->SetMarkerStyle(21);
  tge->SetLineColor(kBlack);
  tge->SetLineWidth(1);
  tge->GetYaxis()->SetRangeUser(20,50);
  tge->GetXaxis()->SetRangeUser(0.78,1.02);
  tge->GetXaxis()->SetTitle("k_{min}");
  tge->GetYaxis()->SetTitle("S/#sqrt{B}");
  tge->Draw("AP");
  
  TGraphErrors *tge2 = new TGraphErrors(10,kmin,sb2,kminerror,sberror2);
  tge2->SetTitle("S/#sqrt{B} vs k_{max}");
  tge2->SetMarkerSize(0.8);
  tge2->SetMarkerColor(kRed);
  tge2->SetMarkerStyle(21);
  tge2->SetLineColor(kBlack);
  tge2->SetLineWidth(1); 
  
  TGraphErrors *tge3 = new TGraphErrors(10,kmin,sb3,kminerror,sberror3);
  tge3->SetTitle("S/#sqrt{B} vs k_{max}");
  tge3->SetMarkerSize(0.8);
  tge3->SetMarkerColor(kRed);
  tge3->SetMarkerStyle(21);
  tge3->SetLineColor(kBlack);
  tge3->SetLineWidth(1);
  tge3->GetYaxis()->SetRangeUser(0,300);

  mg->Add(tge);
  mg->Add(tge2);
  mg->Add(tge3);

  mg->GetXaxis()->SetTitle("k_{max}");
  mg->GetYaxis()->SetTitle("S/#sqrt{B}");

  mg->Draw("AP");

  TPaveLabel *label = new TPaveLabel(0.7765,45.7,0.89,49.9,"a^{0}_{Z}/#Lambda^{2} = 0.9*10^{-5} GeV^{-2}");
  label->SetFillColor(kSpring-4);
  label->SetTextColor(kOrange+10);
  label->SetTextSize(0.5);
  label->Draw("same");
  
  
  TLine *line = new TLine(0.95,21,0.95,38);
  line->SetLineColor(kBlue);
  line->SetLineWidth(3);
  line->Draw("same");
  TLine *line2 = new TLine(0.95,46,0.95,50);
  line2->SetLineColor(kBlue);
  line2->SetLineWidth(3);
  line2->Draw("same");

  CMS_lumi(c_p,18,5);
  c_p->Print(p_png.c_str());

}

void sbvskmax(){

  TCanvas *c_p = new TCanvas ("c_p","c_p");
  string p_png = "images/kmaxselection.png";
  
  float sb1[]={45.4,45.6,46.1,45.2,45.4,43.9,42.7,42.1,38.4,37.0,36.2};
  float sberror1[]={4.6,4.3,4.3,4.0,4.0,3.7,3.5,3.4,2.7,2.5,2.4};

  float sb2[]={67.4,67.6,67.9,66.6,66.2,64.5,62.7,62.0,55.6,53.4,52.3};
  float sberror2[]={6.8,6.4,6.3,5.9,5.8,5.4,5.1,4.9,4.0,3.6,3.4};

  float sb3[]={228.6,229.0,228.1,221.8,220.2,214.8,209.3,206.1,184.9,178.7,175.4};
  float sberror3[]={22.8,21.7,21.0,19.7,19.2,18.1,17.1,16.4,13.2,12.1,11.4};

  
  float kmax[]={1.08,1.09,1.1,1.11,1.12,1.13,1.14,1.15,1.20,1.25,1.30};
  float kmaxerror[]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("");

  TGraphErrors *tge = new TGraphErrors(10,kmax,sb1,kmaxerror,sberror1);
  tge->SetTitle("");
  tge->SetMarkerSize(0.8);
  tge->SetMarkerColor(kRed);
  tge->SetMarkerStyle(21);
  tge->SetLineColor(kBlack);
  tge->SetLineWidth(1);
  tge->GetYaxis()->SetRangeUser(0,70);
  
  TGraphErrors *tge2 = new TGraphErrors(10,kmax,sb2,kmaxerror,sberror2);
  tge2->SetTitle("S/#sqrt{B} vs k_{max}");
  tge2->SetMarkerSize(0.8);
  tge2->SetMarkerColor(kRed);
  tge2->SetMarkerStyle(21);
  tge2->SetLineColor(kBlack);
  tge2->SetLineWidth(1); 
  
  TGraphErrors *tge3 = new TGraphErrors(10,kmax,sb3,kmaxerror,sberror3);
  tge3->SetTitle("S/#sqrt{B} vs k_{max}");
  tge3->SetMarkerSize(0.8);
  tge3->SetMarkerColor(kRed);
  tge3->SetMarkerStyle(21);
  tge3->SetLineColor(kBlack);
  tge3->SetLineWidth(1);
  tge3->GetYaxis()->SetRangeUser(0,300);

  mg->Add(tge);
  //mg->Add(tge2);
  //mg->Add(tge3);

  mg->GetXaxis()->SetTitle("k_{max}");
  mg->GetYaxis()->SetTitle("S/#sqrt{B}");

  mg->Draw("AP");

  TPaveLabel *label = new TPaveLabel(1.16,48.3,1.2585,51.15,"a^{0}_{Z}/#Lambda^{2} = 0.9*10^{-5} GeV^{-2}");
  label->SetFillColor(kSpring-4);
  label->SetTextColor(kOrange+10);
  label->SetTextSize(0.5);
  label->Draw("same");
  
  
  TLine *line = new TLine(1.1,34,1.1,41.3);
  line->SetLineColor(kBlue);
  line->SetLineWidth(3);
  line->Draw("same");
  TLine *line2 = new TLine(1.1,50.5,1.1,51);
  line2->SetLineColor(kBlue);
  line2->SetLineWidth(3);
  line2->Draw("same");
  
  CMS_lumi(c_p,18,5);
  c_p->Print(p_png.c_str());
  
}
