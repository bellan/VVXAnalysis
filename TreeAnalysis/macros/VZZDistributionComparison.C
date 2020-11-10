#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <sstream>
#include <string>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace std;

double scale(double sig,double sig_SM){return sig/sig_SM;}

//cross sections for SM

const double sig_WZZ_SM=0.000141309273192;
const double sig_ZZZ_SM=0;

//cross sections for anomalous phenomena

//cW WZZ
const double sig_WZZ_cW_LI=-5.4669356329e-07;
const double sig_WZZ_cW_QU=6.6865180464e-05;
//cW ZZZ
const double sig_ZZZ_cW_LI=0;
const double sig_ZZZ_cW_QU=3.64336224597e-06;
//cHW WZZ
const double sig_WZZ_cHW_LI=1.37528873887e-07;
const double sig_WZZ_cHW_QU=6.63411412681e-07;
//cHW ZZZ
const double sig_ZZZ_cHW_LI=0;
const double sig_ZZZ_cHW_QU=6.60299819742e-07;

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
  return hist;}

//plots the distribution that sums SM, linear and quadratic contributions

TH1F *DistributionSummer(string operatorname, double parameter, TFile *result1, TFile *result2, double scale1, TFile *result3, double scale2, TLegend *legend, string histname, string title, string axisname, double ymax, double colour){
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
  /*
  histsum->Scale(1/(1+scale1*parameter+scale2*parameter*parameter));
  */
  histsum->SetLineColor(colour);
  ostringstream stream;
  stream << parameter;
  string insert= operatorname+"= "+stream.str();
  legend->AddEntry(histsum,insert.c_str());
  return histsum;}

//compares SM and anomalous couplings distributions for a single variable for WZZ or ZZZ

void DistributionComparison(string operatorname, double parameter, string histname, string title,string axisname, double ymax, TFile *result1, TFile *result2, TFile *result3, double scale1, double scale2){
  
  TLegend *legend= new TLegend(0.75,0.75,0.98,0.95);
  TH1F *hist=DistributionPlot(result1,legend,histname,title,axisname,"SM data",ymax,kBlack);
  TH1F *histsum=DistributionSummer(operatorname,parameter,result1, result2, scale1, result3, scale2, legend, histname, title, axisname, ymax, kRed);
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
}

//compares SM and anomalous couplings distributions for a single variable for VZZ

void DistributionComparisonVZZ(string operatorname, double parameter, string histname, string title,string axisname, double ymax, TFile *result1, TFile *result2, TFile *result3, TFile *result4, TFile *result5, TFile *result6, double sig1, double sig2, double sig3, double sig4, double sig5, double sig6){

  double scale1=sig2/sig1;
  double scale2=sig3/sig1;
  double scale3=sig5/sig4;
  double scale4=sig6/sig4;
  
  TLegend *legend= new TLegend(0.75,0.75,0.98,0.95);
  TH1F *histWZZ=DistributionPlot(result1,legend,histname,title,axisname,"SM data",ymax,kBlack);
  TH1F *histsumWZZ=DistributionSummer(operatorname,parameter,result1, result2, scale1, result3, scale2, legend, histname, title, axisname, ymax, kRed);
  TH1F *histZZZ=DistributionPlot(result4,legend,histname,title,axisname,"SM data",ymax,kBlack);
  TH1F *histsumZZZ=DistributionSummer(operatorname,parameter,result4, result5, scale3, result6, scale4, legend, histname, title, axisname, ymax, kRed);

  TH1F *hist=new TH1F(*histWZZ);
  hist->Add(histWZZ,histZZZ,1,sig2/sig1);
  TH1F *histsum=new TH1F(*histsumWZZ);
  histsum->Add(histsumWZZ,histsumZZZ,1,(1+(sig5/sig4)*parameter+(sig6/sig4)*parameter*parameter)/(1+(sig2/sig1)*parameter+(sig3/sig1)*parameter*parameter));
  
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
}

//compares distributions for multiple variables (modifiable as needed)

void VZZDistributionComparison(string sample, string operatorname, double parameter){

  gErrorIgnoreLevel = kWarning; //to ignore "Info" messages

  TCanvas *c1 = new TCanvas("c1","canvas",0,0,1000,1000);
  TCanvas *c2 = new TCanvas("c2","canvas",0,0,1000,1000);

  //for VZZ
  
  if(sample=="VZZ"){
    
    string filename1= filename("WZZ_SM");
    string filename2= filename("WZZ_"+operatorname+"_LI");
    string filename3= filename("WZZ_"+operatorname+"_QU");
    string filename4= filename("ZZZ_SM");
    string filename5= filename("ZZZ_"+operatorname+"_LI");
    string filename6= filename("ZZZ_"+operatorname+"_QU");
    TFile *result1 = TFile::Open(filename1.c_str());
    TFile *result2 = TFile::Open(filename2.c_str());
    TFile *result3 = TFile::Open(filename3.c_str());
    TFile *result4 = TFile::Open(filename4.c_str());
    TFile *result5 = TFile::Open(filename5.c_str());
    TFile *result6 = TFile::Open(filename6.c_str());

    double sig1=sig_WZZ_SM;
    double sig4=sig_ZZZ_SM;
    double sig2,sig3,sig5,sig6;
    if(operatorname=="cW"){
      sig2=sig_WZZ_cW_LI;
      sig3=sig_WZZ_cW_QU;
      sig5=sig_ZZZ_cW_LI;
      sig6=sig_ZZZ_cW_QU;}
    if(operatorname=="cHW"){
      sig2=sig_WZZ_cHW_LI;
      sig3=sig_WZZ_cHW_QU;
      sig5=sig_ZZZ_cHW_LI;
      sig6=sig_ZZZ_cHW_QU;}
    
    c1->cd();
    DistributionComparisonVZZ(operatorname,parameter,"mass of coupled bosons","Coupled bosons mass","mass (GeV/c^2)",2700,result1,result2,result3,result4,result5,result6,sig1,sig2,sig3,sig4,sig5,sig6);
    c2->cd();
    DistributionComparisonVZZ(operatorname,parameter,"energy of all bosons","Total boson energy","energy (GeV)",1000,result1,result2,result3,result4,result5,result6,sig1,sig2,sig3,sig4,sig5,sig6);
    return;}

  string filename1= filename(sample+"_SM");
  string filename2= filename(sample+"_"+operatorname+"_LI");
  string filename3= filename(sample+"_"+operatorname+"_QU");
  TFile *result1 = TFile::Open(filename1.c_str());
  TFile *result2 = TFile::Open(filename2.c_str());
  TFile *result3 = TFile::Open(filename3.c_str());
  
  /*
  TCanvas *c3 = new TCanvas("c3","canvas",0,0,1000,1000);
  TCanvas *c4 = new TCanvas("c4","canvas",0,0,1000,1000);
  TCanvas *c5 = new TCanvas("c5","canvas",0,0,1000,1000);
  TCanvas *c6 = new TCanvas("c6","canvas",0,0,1000,1000);
  TCanvas *c7 = new TCanvas("c7","canvas",0,0,1000,1000);
  TCanvas *c8 = new TCanvas("c8","canvas",0,0,1000,1000);
  */

  //picks the correct cross sections for the scale factors
  
  double sig_SM,sig_LI,sig_QU;
  if(sample=="WZZ"){
    sig_SM=sig_WZZ_SM;
    if(operatorname=="cW"){
      sig_LI=sig_WZZ_cW_LI;
      sig_QU=sig_WZZ_cW_QU;}
    if(operatorname=="cHW"){
      sig_LI=sig_WZZ_cHW_LI;
      sig_QU=sig_WZZ_cHW_QU;}}
  if(sample=="ZZZ"){
    sig_SM=sig_ZZZ_SM;
    if(operatorname=="cW"){
      sig_LI=sig_ZZZ_cW_LI;
      sig_QU=sig_ZZZ_cW_QU;}
    if(operatorname=="cHW"){
      sig_LI=sig_ZZZ_cHW_LI;
      sig_QU=sig_ZZZ_cHW_QU;}}
  
  double scale1=scale(sig_LI,sig_SM);
  double scale2=scale(sig_QU,sig_SM);

  c1->cd();
  DistributionComparison(operatorname,parameter,"mass of coupled bosons","Coupled bosons mass","mass (GeV/c^2)",2700,result1,result2,result3,scale1,scale2);
  c2->cd();
  DistributionComparison(operatorname,parameter,"energy of all bosons","Total boson energy","energy (GeV)",1000,result1,result2,result3,scale1,scale2);
  /*
  c3->cd();
  DistributionComparison(operatorname,parameter,"mass of tribosons","Triboson mass","mass (GeV/c^2)",1600,result1,result2,result3,scale1,scale2);
  c4->cd();
  DistributionComparison(operatorname,parameter,"energy of major bosons","Major boson energy","energy (GeV)",1000,result1,result2,result3,scale1,scale2);
  c5->cd();
  DistributionComparison(operatorname,parameter,"pt of major bosons","Maximum boson pt","pt (GeV/c)",1200,result1,result2,result3,scale1,scale2);
  c6->cd();
  DistributionComparison(operatorname,parameter,"boson relative angle","Boson relative angle","angle (rad)",900,result1,result2,result3,scale1,scale2);
  c7->cd();
  DistributionComparison(operatorname,parameter,"energy of leptonic bosons","Leptonic boson energy","energy (GeV)",2800,result1,result2,result3,scale1,scale2);
  c8->cd();
  DistributionComparison(operatorname,parameter,"pt of leptonic bosons","Leptonic boson pt","pt (GeV/c)",3600,result1,result2,result3,scale1,scale2);
  */
}
