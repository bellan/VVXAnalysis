#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <sstream>
#include <string>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace std;

double SM_crosssec=0.000141309273192;
double scale(double crosssec){return crosssec/SM_crosssec;}

//cross sections for anomalous phenomena

double sec1=0;
double sec2=5.21541181034e-06;
double sec3=0;
double sec4=7.84718458183e-06;
double sec5=-5.4669356329e-07;
double sec6=6.6865180464e-05;

//scale factors for anomalous phenomenona

double scale1=scale(sec1);
double scale2=scale(sec2);
double scale3=scale(sec3);
double scale4=scale(sec4);
double scale5=scale(sec5);
double scale6=scale(sec6);

//generates the filename input for other functions

string filename(string samplename){return "./results/2016/VVXnocutsAnalyzer_MC/"+samplename+".root";}

//plots a single distribution

TH1F* DistributionPlot(TFile *result, TLegend *legend, string histname, string title, string axisname, string legendname, double ymax, double colour){
  TH1F *hist= (TH1F*)result->Get(histname.c_str());
  hist->SetTitle(title.c_str());
  hist->GetYaxis()->SetRangeUser(0,ymax);
  hist->GetXaxis()->SetTitle(axisname.c_str());
  hist->SetLineColor(colour);
  hist->SetLineWidth(2);
  legend->AddEntry(hist,legendname.c_str(),"l");
  hist->Draw();
  return hist;}

//plots the distribution that sums SM, linear and quadratic contributions

  TH1F *DistributionSummer(TFile *result1, TFile *result2, double scale1, TFile *result3, double scale2, TLegend *legend, string histname, string title, string axisname, string legendname, double ymax, double colour){
  TH1F *hist1= (TH1F*)result1->Get(histname.c_str());
  TH1F *hist2= (TH1F*)result2->Get(histname.c_str());
  TH1F *hist3= (TH1F*)result3->Get(histname.c_str());
  TH1F *histsum= new TH1F(*hist1);
  hist1->SetTitle(title.c_str());
  hist1->GetYaxis()->SetRangeUser(0,ymax);
  hist1->GetXaxis()->SetTitle(axisname.c_str());
  hist1->SetLineWidth(2);
  histsum->Add(hist1,hist2,1,scale1);
  histsum->Add(histsum,hist3,1,scale2*scale2);
  histsum->Scale(1/(1+scale1+scale2*scale2));
  histsum->SetLineColor(colour);
  legend->AddEntry(histsum,legendname.c_str());
  histsum->Draw("same");
  return histsum;}

//compares SM and anomalous couplings distributions for a single variable

  void DistributionComparison(string histname, string title,string axisname, double ymax, TFile *result1, TFile *result2, TFile *result3){
  TLegend *legend= new TLegend(0.75,0.75,0.98,0.95);
  TH1F *hist=DistributionPlot(result1,legend,histname,title,axisname,"SM data",ymax,kBlack);
  /*
  TH1F *histsum=DistributionSummer(result1, result2, scale1, result3, scale2, legend, histname, title, axisname, "cW=0.65", ymax, kMagenta);
  TH1F *histsum=DistributionSummer(result1, result2, scale3, result3, scale4, legend, histname, title, axisname, "cW=0.8", ymax, kGreen);
  */
  TH1F *histsum=DistributionSummer(result1, result2, scale5, result3, scale6, legend, histname, title, axisname, "cW=1", ymax, kRed);

  double chi=hist->Chi2Test(histsum,"WWCHI2");
  ostringstream stream1,stream2;
  stream1 << hist->Chi2Test(histsum,"WWCHI2/NDF");
  string insert1= "chi2/NDF: "+stream1.str();
  legend->AddEntry((TObject*)0,insert1.c_str(),"");
  stream2 << hist->Chi2Test(histsum,"WW");
  string insert2= "p value: "+stream2.str();
  legend->AddEntry((TObject*)0,insert2.c_str(),"");
  legend->Draw();
}  

//compares distributions for multiple variables (modifiable as needed)

  void VZZDistributionComparison(string samplename1, string samplename2, string samplename3){
    
  gErrorIgnoreLevel = kWarning; //to ignore "Info" messages

  string filename1= filename(samplename1);
  string filename2= filename(samplename2);
  string filename3= filename(samplename3);
  TFile *result1 = TFile::Open(filename1.c_str());
  TFile *result2 = TFile::Open(filename2.c_str());
  TFile *result3 = TFile::Open(filename3.c_str());

  TCanvas *c1 = new TCanvas("c1","canvas",0,0,1000,1000);
  TCanvas *c2 = new TCanvas("c2","canvas",0,0,1000,1000);
  TCanvas *c3 = new TCanvas("c3","canvas",0,0,1000,1000);
  TCanvas *c4 = new TCanvas("c4","canvas",0,0,1000,1000);
  TCanvas *c5 = new TCanvas("c5","canvas",0,0,1000,1000);
  TCanvas *c6 = new TCanvas("c6","canvas",0,0,1000,1000);
  TCanvas *c7 = new TCanvas("c7","canvas",0,0,1000,1000);
  TCanvas *c8 = new TCanvas("c8","canvas",0,0,1000,1000);
  
  c1->cd();
  DistributionComparison("mass of tribosons","Triboson mass","mass (GeV/c^2)",1600,result1,result2,result3);
  c2->cd();
  DistributionComparison("mass of coupled bosons","Coupled bosons mass","mass (GeV/c^2)",2500,result1,result2,result3);
  c3->cd();
  DistributionComparison("energy of all bosons","Total boson energy","energy (GeV)",1600,result1,result2,result3);
  c4->cd();
  DistributionComparison("energy of major bosons","Major boson energy","energy (GeV)",1000,result1,result2,result3);
  c5->cd();
  DistributionComparison("pt of major bosons","Maximum boson pt","pt (GeV/c)",1200,result1,result2,result3);
  c6->cd();
  DistributionComparison("boson relative angle","Boson relative angle","angle (rad)",900,result1,result2,result3);
  c7->cd();
  DistributionComparison("energy of leptonic bosons","Leptonic boson energy","energy (GeV)",2800,result1,result2,result3);
  c8->cd();
  DistributionComparison("pt of leptonic bosons","Leptonic boson pt","pt (GeV/c)",3600,result1,result2,result3);
}

//compares distribution of a variable (with chi squared) for two different samples

void TwoDistributionComparison(string histname, string title,string axisname, string legendname1, string legendname2, double ymax, string samplename1, string samplename2, string year1, string year2, double scale){

  string filename1= "./results/"+year1+"/VVXnocutsAnalyzer_MC/"+samplename1+".root";
  string filename2= "./results/"+year2+"/VVXnocutsAnalyzer_MC/"+samplename2+".root";
  TFile *result1 = TFile::Open(filename1.c_str());
  TFile *result2 = TFile::Open(filename2.c_str());
  
  TH1F *hist1= (TH1F*)result1->Get(histname.c_str());
  TH1F *hist2= (TH1F*)result2->Get(histname.c_str());
  hist1->SetTitle(title.c_str());
  hist1->GetYaxis()->SetRangeUser(0,ymax);
  hist1->GetXaxis()->SetTitle(axisname.c_str());
  hist1->SetLineColor(1);
  hist2->SetLineColor(2);
  hist1->SetLineWidth(2);
  hist2->SetLineWidth(2);
  hist2->Scale(scale);
  hist1->Draw();
  hist2->Draw("same");
  TLegend *legend= new TLegend(0.75,0.75,0.98,0.95);
  legend->AddEntry(hist1,legendname1.c_str(),"l");
  legend->AddEntry(hist2,legendname2.c_str(),"l");
  ostringstream stream;
  stream << hist1->Chi2Test(hist2,"WW");
  string insert= "p value: "+stream.str();
  legend->AddEntry((TObject*)0,insert.c_str(),"");
  legend->Draw();
}
