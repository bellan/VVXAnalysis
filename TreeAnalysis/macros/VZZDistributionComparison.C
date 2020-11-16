/*
Macro to analyze variable distributions for anomalous couplings, as analyzed by VVXnocutsAnalyzer.  Path will probably need changing.
Author: Marozzo Giovanni Battista
Date: 2020/11/25
*/

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <sstream>
#include <string>
#include <TLatex.h>

using namespace std;

double scale(double sig,double sig_SM){return sig/sig_SM;}

//cross sections for SM

const double sig_WZZ_SM=0.000141309273192;    //pb
const double sig_ZZZ_SM=0.000132818264676;    //pb

//cross sections for anomalous phenomena

//cW WZZ
const double sig_WZZ_cW_LI=-5.4669356329e-07; //pb
const double sig_WZZ_cW_QU=6.6865180464e-05;  //pb
//cW ZZZ
const double sig_ZZZ_cW_LI=2.5001339964e-08;  //pb
const double sig_ZZZ_cW_QU=3.64336224597e-06; //pb

//generates the filename input for other functions

string filename(string samplename){return "./results/2016/VVXnocutsAnalyzer_MC/"+samplename+".root";}

//plots a single distribution

TH1F* DistributionPlot(TFile *result, TLegend *legend, string histname, string title, string axisname, string legendname, double ymax, double colour, double sig_SM){
  TH1F *hist= (TH1F*)result->Get(histname.c_str());
  hist->SetTitle(title.c_str());
  hist->GetYaxis()->SetRangeUser(0,ymax);
  hist->GetXaxis()->SetTitle(axisname.c_str());
  hist->SetLineColor(colour);
  hist->SetLineWidth(2);
  ostringstream stream;
  stream << 1000*sig_SM;
  string legendentry=legendname+"; #sigma= "+stream.str()+" fb";
  legend->AddEntry(hist,legendentry.c_str(),"l");
  return hist;}

//plots the distribution that sums SM, linear and quadratic contributions

TH1F *DistributionSummer(double parameter, TFile *result1, TFile *result2, double scale1, TFile *result3, double scale2, TLegend *legend, string histname, string title, string axisname, double ymax, double colour,double sig_SM){
  TH1F *hist1= (TH1F*)result1->Get(histname.c_str());
  TH1F *hist2= (TH1F*)result2->Get(histname.c_str());
  TH1F *hist3= (TH1F*)result3->Get(histname.c_str());
  TH1F *histsum= new TH1F(*hist1);
  hist1->SetTitle(title.c_str());
  hist1->GetYaxis()->SetRangeUser(0,ymax);
  hist1->GetXaxis()->SetTitle(axisname.c_str());
  hist1->SetLineWidth(2);
  histsum->Add(hist1,hist2,1,scale1*parameter);
  histsum->Add(histsum,hist3,1,scale2*parameter*parameter);
  histsum->SetLineColor(colour);
  ostringstream stream1,stream2;
  stream1 << parameter;
  stream2 << 1000*sig_SM*(1+scale1*parameter+scale2*parameter*parameter);
  string insert= "cW= "+stream1.str()+"; #sigma= "+stream2.str()+" fb";
  legend->AddEntry(histsum,insert.c_str(),"l");
  return histsum;}

//compares SM and anomalous couplings distributions for a single variable for WZZ or ZZZ

TH1F *DistributionComparison(double parameter, string histname, string title,string axisname, double ymax, TFile *result1, TFile *result2, TFile *result3, double scale1, double scale2, double sig_SM){
  
  TLegend *legend= new TLegend(0.75,0.75,0.98,0.95);
  legend->AddEntry((TObject*)0,"#it{L}=150 fb^{-1}","");
  TH1F *hist=DistributionPlot(result1,legend,histname,title,axisname,"SM data",ymax,kBlack,sig_SM);
  TH1F *histsum=DistributionSummer(parameter,result1, result2, scale1, result3, scale2, legend, histname, title, axisname, ymax, kRed,sig_SM);
  hist->Scale(150*sig_SM*pow(10,-1));
  histsum->Scale(150*sig_SM*pow(10,-1));
  hist->Draw();
  histsum->Draw("same");

  ostringstream stream1,stream2,stream3;
  stream1 << hist->Chi2Test(histsum,"WWCHI2");
  string insert1= "chi2: "+stream1.str();
  legend->AddEntry((TObject*)0,insert1.c_str(),"");
  stream2 << hist->Chi2Test(histsum,"WWCHI2/NDF");
  string insert2= "chi2/NDF: "+stream2.str();
  legend->AddEntry((TObject*)0,insert2.c_str(),"");
  legend->Draw();
  stream3 << hist->Chi2Test(histsum,"WW");
  string insert3= "p value: "+stream3.str();
  legend->AddEntry((TObject*)0,insert3.c_str(),"");
  legend->Draw();
  TH1F *weighthist= new TH1F(*histsum);
  weighthist->Divide(hist);
  return weighthist;
}


//compares distributions for multiple variables (modifiable as needed)

void VZZDistributionComparison(string sample, double parameter){

  gErrorIgnoreLevel = kWarning; //to ignore "Info" messages

  string filename1= filename(sample+"_SM");
  string filename2= filename(sample+"_cW_LI");
  string filename3= filename(sample+"_cW_QU");
  TFile *result1 = TFile::Open(filename1.c_str());
  TFile *result2 = TFile::Open(filename2.c_str());
  TFile *result3 = TFile::Open(filename3.c_str());
  
  TCanvas *c1 = new TCanvas("c1","canvas",0,0,1000,1000);
  
  //picks the correct cross sections for the scale factors
  
  double sig_SM,sig_LI,sig_QU;
  if(sample=="WZZ"){
    sig_SM=sig_WZZ_SM;
    sig_LI=sig_WZZ_cW_LI;
    sig_QU=sig_WZZ_cW_QU;}
  if(sample=="ZZZ"){
    sig_SM=sig_ZZZ_SM;
    sig_LI=sig_ZZZ_cW_LI;
    sig_QU=sig_ZZZ_cW_QU;}
  
  double scale1=scale(sig_LI,sig_SM);
  double scale2=scale(sig_QU,sig_SM);

  c1->cd();
  TH1F *weighthist=DistributionComparison(parameter,"energy of all bosons","Total boson energy","energy (GeV)",1000,result1,result2,result3,scale1,scale2,sig_SM);
  TFile *file=new TFile("aGCweighthist.root","recreate");
  weighthist->Write();
}
