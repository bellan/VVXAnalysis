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
#endif

void RooUnfold_Mass_data(string finalstate = "4e", string matrix = "4e", bool bayes = 1, int it = 2)
{
  
#ifdef __CINT__
  gSystem->Load("../RooUnfold-1.1.1/libRooUnfold");
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
  TMatrixD * cov;
  TVectorD * Vcov_stat; 
  TVectorD * Vcov_unf; 
  
  // TFile *file = new TFile("/afs/cern.ch/user/l/lfinco/work/VVScattering/CMSSW_5_3_11/src/VVXAnalysis/TreeAnalysis/macros/output.root");
  // TFile *data = new TFile("/afs/cern.ch/user/l/lfinco/work/VVScattering/CMSSW_5_3_11/src/VVXAnalysis/TreeAnalysis/macros/DataToUnfold.root");
  TFile *file = new TFile("../output.root");
  TFile *data = new TFile("../DataToUnfold.root");
  TFile *output = new TFile("testUnfoldData.root", "UPDATE");
 
  string matrixName = "ResMat_qqgg_ZZTo" + matrix + "_st_01";
  string histoName = "Mass_qqgg_ZZTo" + matrix + "_st_01";
  string histoNamegen = "MassGen_qqgg_ZZTo" + matrix + "_st_01";
  string histoName_unf = "DataminusBkg_Mass_ZZTo"+finalstate;
  string histoNameGen_unf = "MassGen_qqgg_ZZTo" + finalstate + "_st_01";

  h_measured = (TH1*) file->Get(histoName.c_str());  
  h_true = (TH1*) file->Get(histoNamegen.c_str());
  h_Resmat = (TH2*)file->Get(matrixName.c_str()); 
  h_measured_unf = (TH1*) data->Get(histoName_unf.c_str());  
  h_true_unf = (TH1*) file->Get(histoNameGen_unf.c_str());

  RooUnfoldResponse response(h_measured, h_true, h_Resmat, "response", "response"); 
  RooUnfoldBayes unfold_bayes(&response, h_measured_unf,it);
  RooUnfoldSvd   unfold_svd(&response, h_measured_unf, it);
  
  if (bayes ==1){
    //unfold_bayes.IncludeSystematics(RooUnfold::kCovariance);      
    hReco= (TH1*) unfold_bayes.Hreco(RooUnfold::kCovariance); 
    std::cout << "3" << std::endl;
    unfold_bayes.PrintTable (cout, h_true_unf,RooUnfold::kCovariance);
    cov = (TMatrixD*) unfold_bayes.Ereco(RooUnfold::kCovariance);
    Vcov_stat = (TVectorD*) unfold_bayes.ErecoV(RooUnfold::kNoError);
    Vcov_unf =(TVectorD*) unfold_bayes.ErecoV(RooUnfold::kCovariance);
  }
  else { 
    // unfold_svd.IncludeSystematics(RooUnfold::kCovariance);    
    hReco= (TH1*) unfold_svd.Hreco(RooUnfold::kCovariance);
    unfold_svd.PrintTable (cout, h_true_unf,RooUnfold::kCovariance);
    cov =  (TMatrixD*)unfold_svd.Ereco(RooUnfold::kCovariance);
    Vcov_stat = (TVectorD*)unfold_svd.ErecoV(RooUnfold::kNoError);
    Vcov_unf =(TVectorD*) unfold_svd.ErecoV(RooUnfold::kCovariance);
   
  } 
  cov->Print();
  Vcov_stat->Print();
  Vcov_unf->Print();
  //std::cout << hReco->GetBinContent(1) << " " << hReco->GetBinError(1) << std::endl;

  TCanvas *c = new TCanvas ("c","c");
  TLegend *leg = new TLegend(0.65,0.65,0.6,0.85); 
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
 
  string title = "Unfolding m_{" + finalstate + "} distribution using " + matrix + "-response matrix";
  string XaxisTitle = "m_{" + finalstate + "}";
  hReco->Draw("E TEXT"); 
  hReco->SetTitle(title.c_str());
  hReco->GetXaxis()->SetRange(0,25);
  hReco->GetXaxis()->SetTitle(XaxisTitle.c_str());
 
  h_measured_unf->Draw("E SAME");
  h_measured_unf->SetLineColor(2);
  h_true_unf->SetLineColor(8);
  h_true_unf->Draw("E SAME");
  
  leg->AddEntry(hReco,"unfolded distribution","l"); 
  leg->Draw(); 
  leg->AddEntry(h_true_unf,"true distribution","l"); 
  leg->Draw("SAME");
  leg->AddEntry(h_measured_unf,"reco distribution","l"); 
  leg->Draw("SAME");

  string UnfHistoName = "ZZTo"+ finalstate +"_Mass_01";
  string HistoName = "ZZTo"+ finalstate +"_Mass_RECO_MC_01";
  string TrueHistoName = "ZZTo"+ finalstate +"_Mass_GEN_01";
  
  output->cd();   
  hReco->Write(UnfHistoName.c_str());
  h_measured_unf->Write(HistoName.c_str());
  h_true_unf->Write(TrueHistoName.c_str());
  output->Close();
}

void MakeAllFinalStates (bool bayes = 1, int it = 2){
  RooUnfold_Mass_data("4e","4e",bayes,it);
  RooUnfold_Mass_data("4m","4m",bayes,it);
  RooUnfold_Mass_data("2e2m","2e2m",bayes,it);
}

#ifndef __CINT__
int main () { RooUnfold_Mass_data(); return 0; }  // Main program when run stand-alone
#endif
