#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <string>
#include <sstream>
using std::cout;
using std::endl;

//#include "TRandom.h"
#include <TROOT.h>
#include <TFile.h>
#include <TStyle.h>
#include "TH1D.h"
#include "TMatrix.h"
//#include "TArray.h"
#include "TH2.h"
#include "TLegend.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
#include "TCanvas.h"
//#include "../RooUnfold-1.1.1/src/TSVDUnfold.h"
#include "TPad.h"
#endif

//This macro tests the unfolding procedure on MC samples

void Unfold_MCtest(string var = "Mass",string finalstate = "4e", bool bayes = 0, int it = 4, bool MadMatrix =1, bool MadDistribution = 1, bool FullSample = 0)
{
  
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  
  TH1 * h_true; 
  TH1 * h_measured; 
  TH1 * h_measured_unf; 
  TH1 * h_true_unf; 
  TH2 * h_Resmat;
  TH1 * hReco;
  // TMatrixD * cov;
  // TVectorD * Vcov_stat; 
  // TVectorD * Vcov_unf; 
  string  filePath = "/../../";
  string fileName = filePath+var+"_test/matrices.root";
  string fileName_Pow =filePath +var+"_test/matrices_Pow.root"; 

  TFile *file_mad = new TFile(fileName.c_str());
  TFile *file_pow = new TFile(fileName_Pow.c_str());
 
  //TFile *output = new TFile("testUnfoldData.root", "UPDATE");
 
  string matrixName; 
  string histoName; 
  string histoNamegen; 
  string histoNameGen_unf;
  string histoName_unf;

 if(FullSample == 0){
   matrixName = "ResMat_qqggJJ_"+var+"_ZZTo" + finalstate + "_st_0";
   histoName = var+"_qqggJJ_ZZTo" + finalstate + "_st_0";
   histoNamegen = var+"Gen_qqggJJ_ZZTo" + finalstate + "_st_0";
   histoName_unf =  var+"_qqggJJ_ZZTo" + finalstate + "_st_1"; //same as histoName, MC test n1
   histoNameGen_unf = var+"Gen_qqggJJ_ZZTo" + finalstate + "_st_1";
 }
 
 else{
   matrixName = "ResMat_qqggJJ_"+var+"_ZZTo" + finalstate + "_st_01";
   histoName = var+"_qqggJJ_ZZTo" + finalstate + "_st_01";
   histoNamegen = var+"Gen_qqggJJ_ZZTo" + finalstate + "_st_01";
   histoName_unf =  var+"_qqggJJ_ZZTo" + finalstate + "_st_01"; //same as histoName, MC test n1
   histoNameGen_unf = var+"Gen_qqggJJ_ZZTo" + finalstate + "_st_01";
 } 
  
 if(MadMatrix == 1){
   h_measured = (TH1*) file_mad->Get(histoName.c_str());  
   h_true = (TH1*) file_mad->Get(histoNamegen.c_str());
   h_Resmat = (TH2*)file_mad->Get(matrixName.c_str()); 
 }

 else{
   h_measured = (TH1*) file_pow->Get(histoName.c_str());  
   h_true = (TH1*) file_pow->Get(histoNamegen.c_str());
   h_Resmat = (TH2*)file_pow->Get(matrixName.c_str()); 
 }

 if(MadDistribution == 1){
   h_measured_unf = (TH1*) file_mad->Get(histoName_unf.c_str());  
   h_true_unf = (TH1*) file_mad->Get(histoNameGen_unf.c_str());
 }
 
 else{
   h_measured_unf = (TH1*) file_pow->Get(histoName_unf.c_str());  
   h_true_unf = (TH1*) file_pow->Get(histoNameGen_unf.c_str());
 }


  RooUnfoldResponse response(h_measured, h_true, h_Resmat, "response", "response"); 
  RooUnfoldBayes unfold_bayes(&response, h_measured_unf,it);
  RooUnfoldSvd   unfold_svd(&response, h_measured_unf, it);
  
  if (bayes ==1){
    hReco= (TH1*) unfold_bayes.Hreco(RooUnfold::kCovariance);
  }
  else { 
    hReco= (TH1*) unfold_svd.Hreco(RooUnfold::kCovariance);
  }

  string fs;
  string XaxisTitle;
  if(finalstate == "4m") fs = "4#mu";
  else if(finalstate == "2e2m") fs = "2e2#mu";
  else fs = finalstate;
  
  if(var == "Mass") XaxisTitle = "m_{" + fs + "}"; 
  else if(var == "Jets")  XaxisTitle = "N Jets (" + fs + " final state)";

  string YaxisTitle = "Events";
  string YaxisTitle2 = "Unfolded/True";
  TCanvas *c = new TCanvas ("c","c");
  TLegend *leg = new TLegend(0.65,0.65,0.45,0.85); 
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42); 

  c->cd();
  
  TPad  *pad1 = new TPad("pad1","", 0., 0.2, 1.0, 1.0);
  pad1->SetTopMargin (0.10);
  pad1->SetRightMargin (0.10);
  pad1->SetLeftMargin (0.10);
  pad1->Draw();
  
  c->cd();
  TPad  *pad2 = new TPad("pad2", "", 0., 0.0,  1.0, 0.2);
  pad2->SetTopMargin (0.10);
  pad2->SetRightMargin (0.10);
  pad2->SetLeftMargin (0.10);
  pad2->Draw(); 
    
  pad1->cd();
  //h_true_unf->SetTitle(title.c_str());
  h_true_unf->SetTitle("");
  h_true_unf->GetXaxis()->SetRange(0,25);
  h_true_unf->GetXaxis()->SetTitle(XaxisTitle.c_str()); 
  h_true_unf->GetYaxis()->SetTitle(YaxisTitle.c_str());
  // float max = h_true_unf->GetBinContent(2)+20;
  // h_true_unf->SetMaximum(max);
  h_true_unf->SetMinimum(0); 
  h_true_unf->SetMarkerColor(8);
  h_true_unf->SetMarkerStyle(8);
  h_true_unf->SetLineColor(8); 
  h_true_unf->SetLineWidth(1);
  h_true_unf->Draw("HIST E");
  
  hReco->SetLineColor(kBlue);
  hReco->SetLineWidth(1);
  hReco->SetMarkerColor(kBlue);
  hReco->SetMarkerStyle(8);
  hReco->Draw("E SAME"); 
  h_measured_unf->Draw("HIST E SAME");
  h_measured_unf->SetMarkerColor(2);
  h_measured_unf->SetMarkerStyle(8);
  h_measured_unf->SetLineColor(2);
  h_measured_unf->SetLineWidth(1);

  leg->AddEntry(hReco,"unfolded distribution","lep"); 
  leg->Draw(); 
  leg->AddEntry(h_true_unf,"true distribution","lep"); 
  leg->Draw("SAME");
  leg->AddEntry(h_measured_unf,"reco distribution","lep"); 
  leg->Draw("SAME");
  
  TH1 * hReco_r = (TH1*) hReco->Clone();
    
  for(int k =1;k<9;k++){
    float unf=0;
    float tr = 0;
    float ratio =0;
    float err_unf=0;
    float err_tr = 0;
    float err_ratio =0;
    unf = hReco->GetBinContent(k);
    tr = h_true_unf->GetBinContent(k); 
    err_unf = hReco->GetBinError(k);
    err_tr = h_true_unf->GetBinError(k);
    ratio = unf/tr;
    err_ratio = sqrt((err_tr/tr)*(err_tr/tr)+(err_unf/unf)*(err_unf/unf));
    hReco_r->SetBinContent(k,ratio); 
    hReco_r->SetBinError(k,err_ratio);
  }

  pad2->cd();
  hReco_r->SetTitle(""); 
  hReco_r->SetMarkerColor(1);
  hReco_r->SetLineColor(1); 
  hReco_r->GetXaxis()->SetTitle(XaxisTitle.c_str()); 
  hReco_r->GetYaxis()->SetTitle(YaxisTitle2.c_str());
  hReco_r-> SetMaximum(1.5); 
  hReco_r-> SetMinimum(0.5);
  hReco_r->Draw("E");

  string pdf;
  string png;
  string sample;
  string matrix;
  string distr;
  if(FullSample == 0) sample = "HalfSample";
  else  sample = "FullSample";
  
  if(MadMatrix == 1) matrix = "MadMatrix";
 else matrix = "PowMatrix";

 if(MadDistribution == 1)distr = "MadDistr";
 else distr = "PowDistr";

  pdf = var+"_ZZTo" + finalstate + "_" +var+"_"+ matrix + "_" + distr + "_" + sample + ".pdf";
  png = var+"_ZZTo" + finalstate + "_" +var+"_"+ matrix + "_" + distr + "_" + sample + ".png";
  // c->Print(pdf.c_str());
  // c->Print(png.c_str());
   // // // output->cd();   
  // // // hReco->Write(UnfHistoName.c_str());
  // // // h_measured_unf->Write(HistoName.c_str());
  // // // h_true_unf->Write(TrueHistoName.c_str());
  // // // output->Close();


}

void MakeAllTest(string var = "Mass",string finalstate = "4e"){
  
  Unfold_MCtest(var.c_str(),finalstate.c_str(), 0, 4,1, 1, 1);
  Unfold_MCtest(var.c_str(),finalstate.c_str(), 0, 4,1, 1, 0);
  Unfold_MCtest(var.c_str(),finalstate.c_str(), 0, 4,1, 0, 1);
  Unfold_MCtest(var.c_str(),finalstate.c_str(), 0, 4,0, 1, 1);
  Unfold_MCtest(var.c_str(),finalstate.c_str(), 0, 4,0, 0, 1);
  Unfold_MCtest(var.c_str(),finalstate.c_str(), 0, 4,0, 0, 0); 
}

void MakeAllFinalStates (string var = "Mass"){//,bool MadOnPow =0,bool PowOnMad =0){
  MakeAllTest(var.c_str(),"4e");
  MakeAllTest(var.c_str(),"4m");
  MakeAllTest(var.c_str(),"2e2m");
}

#ifndef __CINT__
int main () { Unfold_MCtest(); return 0; }  // Main program when run stand-alone
#endif
