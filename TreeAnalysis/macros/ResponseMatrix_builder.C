#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAttLine.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TTree.h>
#include <iostream>
#include <string>
#include <sstream> 

void ResponseMatrix_builder(string dataset = "01", string finalstate = "4l", int xs_qq = 0, int xs_gg =0)//bool eta = false)
{
  //This macro builds the response matrix (measured X truth), the "measured" histo and the "truth" histo, that are the projections of "response" onto the X-axis and Y-axis respectively, but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and in "truth" for unmeasured events (inefficiency). using  and the TH1 measured and truth 

  //dataset == 01 full dataset
  //dataset == 0 first half of the dataset
  //dataset == 1 second half of the dataset

  //xs_qq == 0 matrix and histos built with central values of qq cross sections
  //xs_qq == 1 matrix and histos built with sigma_qq + 4.44% (theoretical uncertainty on sigma_qq)
  //xs_qq == -1 matrix and histos built with sigma_qq - 4.44% (theoretical uncertainty on sigma_qq)

  //xs_gg == 0 matrix and histos built with central values of gg cross sections
  //xs_gg == 1 matrix and histos built with sigma_qq + 25.36% (theoretical uncertainty on sigma_gg)
  //xs_gg == -1 matrix and histos built with sigma_gg - 25.36% (theoretical uncertainty on sigma_gg)
 
  gROOT->Reset();  

  gROOT->SetStyle("Plain");   
  
  gStyle->SetOptStat(0);

  TFile *output = new TFile("matrices.root", "UPDATE");
 
  //Reco samples (response matrices and signal region distributions) //ZZRecoAnalyzer not yet done
  TFile *ggZZTo2e2mu_r = new TFile("../results/ZZRecoAnalyzer_SR/ggTo2e2mu_SMHContinInterf-MCFM67_H125.6.root");
  TFile *ggZZTo4e_r = new TFile("../results/ZZRecoAnalyzer_SR/ggTo4e_SMHContinInterf-MCFM67_H125.6.root");
  TFile *ggZZTo4mu_r = new TFile("../results/ZZRecoAnalyzer_SR/ggTo4mu_SMHContinInterf-MCFM67_H125.6.root");
  TFile *ZZJetsTo4l_r = new TFile("../results/ZZRecoAnalyzer_SR/ZZJetsTo4L.root");
  
  //Truth samples (signal definition distributions) 
  TFile *ggZZTo2e2mu_g = new TFile("../results/ZZMCAnalyzer_MC/ggTo2e2mu_SMHContinInterf-MCFM67_H125.6.root");
  TFile *ggZZTo4e_g = new TFile("../results/ZZMCAnalyzer_MC/ggTo4e_SMHContinInterf-MCFM67_H125.6.root");
  TFile *ggZZTo4mu_g = new TFile("../results/ZZMCAnalyzer_MC/ggTo4mu_SMHContinInterf-MCFM67_H125.6.root");
  TFile *ZZJetsTo4l_g = new TFile("../results/ZZMCAnalyzer_MC/ZZJetsTo4L.root");
  
  //rescue file
  TFile *rescue = new TFile("../results/ZZRecoAnalyzer_SR/tt.root");

  TH2 * h_Resmat;
  TH2 * h_Resmat_4l;
  TH2 * h_Resmat_gg4mu;
  TH2 * h_Resmat_gg4e;
  TH2 * h_Resmat_gg2e2mu;  
  TH2 * h_Resmat_4lTot;
  TH2 * h_Resmat_ggTot;
 
  TH1 * h_4l;
  TH1 * h_gg4mu;
  TH1 * h_gg4e;
  TH1 * h_gg2e2mu; 
  TH1 * h_4l_gen;
  TH1 * h_gg4mu_gen;
  TH1 * h_gg4e_gen;
  TH1 * h_gg2e2mu_gen; 
  TH1 * h_4l_c;
  TH1 * h_4l_gen_c;
  TH2 * h_Resmat_normTot;

  TH1 *h_safe; 
  TH1 *h_safe_tmp;
  TH2 *h_Resmat_safe; 
  TH2 *h_Resmat_safe_tmp;
  
  string safeMatrixName = "ResMat_ZZTo2e2m_Mass_01";
  string safeHistoName = "ZZTo2e2m_Mass_01"; 
 
  h_Resmat_safe_tmp = (TH2*) rescue->Get(safeMatrixName.c_str()); 
  h_safe_tmp = (TH1*) rescue->Get(safeHistoName.c_str()); 
  h_Resmat_safe = (TH2*)h_Resmat_safe_tmp->Clone(safeMatrixName.c_str()); 
  h_safe = (TH1*)h_safe_tmp->Clone(safeHistoName.c_str());

  for(int k=1; k<10; k++){
    for(int l=1; l<10; l++){
      h_safe->SetBinContent(l,0.); 
      h_Resmat_safe->SetBinContent(l,k,0.);
    }
  }
  
  float unc_qq = 0; 
  float unc_gg = 0; 
  float totalint = 0; 
  
  string matrixName = "ResMat_ZZTo" + finalstate + "_Mass_" + dataset;
  string histoName_reco = "ZZTo" + finalstate + "_Mass_" +dataset;
  string histoName_gen =  "ZZTo" + finalstate + "_MassGen_" +dataset;

  h_Resmat_gg4mu = (TH2*) ggZZTo4mu_r->Get(matrixName.c_str()); 
  h_Resmat_gg4e = (TH2*) ggZZTo4e_r->Get(matrixName.c_str()); 
  h_Resmat_gg2e2mu = (TH2*) ggZZTo2e2mu_r->Get(matrixName.c_str()); 
  h_Resmat_4l = (TH2*) ZZJetsTo4l_r->Get(matrixName.c_str());

  h_gg4mu = (TH1*) ggZZTo4mu_r->Get(histoName_reco.c_str()); 
  h_gg4e = (TH1*) ggZZTo4e_r->Get(histoName_reco.c_str()); 
  h_gg2e2mu = (TH1*) ggZZTo2e2mu_r->Get(histoName_reco.c_str()); 
  h_4l = (TH1*) ZZJetsTo4l_r->Get(histoName_reco.c_str()); 

  h_gg4mu_gen = (TH1*) ggZZTo4mu_g->Get(histoName_gen.c_str()); 
  h_gg4e_gen = (TH1*) ggZZTo4e_g->Get(histoName_gen.c_str()); 
  h_gg2e2mu_gen = (TH1*) ggZZTo2e2mu_g->Get(histoName_gen.c_str()); 
  h_4l_gen = (TH1*) ZZJetsTo4l_g->Get(histoName_gen.c_str()); 
  
  if(h_Resmat_gg4mu == NULL) h_Resmat_gg4mu = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_gg4e == NULL) h_Resmat_gg4e = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_gg2e2mu == NULL) h_Resmat_gg2e2mu = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_gg4mu == NULL) h_gg4mu = (TH2*) h_safe->Clone("h_safe");
  if(h_gg4e == NULL) h_gg4e = (TH2*) h_safe->Clone("h_safe");
  if(h_gg2e2mu == NULL) h_gg2e2mu = (TH2*) h_safe->Clone("h_safe");
  if(h_gg4mu_gen == NULL) h_gg4mu_gen = (TH2*) h_safe->Clone("h_safe");
  if(h_gg4e_gen == NULL) h_gg4e_gen = (TH2*) h_safe->Clone("h_safe");
  if(h_gg2e2mu_gen == NULL) h_gg2e2mu_gen = (TH2*) h_safe->Clone("h_safe");

  h_gg4mu->Add(h_gg4e);
  h_gg4mu->Add(h_gg2e2mu);

  h_gg4mu_gen->Add(h_gg4e_gen);
  h_gg4mu_gen->Add(h_gg2e2mu_gen);

  h_Resmat_gg4mu -> Add(h_Resmat_gg4e, 1);
  h_Resmat_gg4mu -> Add(h_Resmat_gg2e2mu, 1);
 
  h_Resmat_4lTot = (TH2*)h_Resmat_4l->Clone("h_Resmat_4l");
  h_Resmat_ggTot = (TH2*)h_Resmat_gg4mu->Clone("h_Resmat_gg"); 
 
  unc_qq = 0.0444; //pdf: 3.4%  scale: 2.85% 
  unc_gg = 0.2536;//pdf: 7.10%  scale: 24.35% 
 
  if(xs_qq == 0) {
    h_Resmat_4lTot->Scale(1);
    h_4l->Scale(1);
    h_4l_gen->Scale(1);
  } 
  else if(xs_qq == 1) {
    h_Resmat_4lTot->Scale(1+unc_qq); 
    h_4l->Scale(1+unc_qq); 
    h_4l_gen->Scale(1+unc_qq);
  }
  else if(xs_qq == -1){
    h_Resmat_4lTot->Scale(1-unc_qq); 
    h_4l->Scale(1-unc_qq); 
    h_4l_gen->Scale(1-unc_qq);
  }
  else std::cout << "Error: xs_qq must be -1, 0 or 1" << std::endl;

  if(xs_gg == 0) {
    h_Resmat_ggTot->Scale(1);
    h_gg4mu->Scale(1);
    h_gg4mu_gen->Scale(1);
  }
  else if(xs_gg == 1) {
    h_Resmat_ggTot->Scale(1+unc_gg); 
    h_gg4mu->Scale(1+unc_gg); 
    h_gg4mu_gen->Scale(1+unc_gg);
  }
  else if(xs_gg == -1) {
    h_Resmat_ggTot->Scale(1-unc_gg); 
    h_gg4mu->Scale(1-unc_gg); 
    h_gg4mu_gen->Scale(1-unc_gg);
  }
  else std::cout << "Error: xs_gg must be -1, 0 or 1" << std::endl;

  h_Resmat = (TH2*)h_Resmat_4lTot->Clone("h_Resmat_4lTot");
  h_Resmat->Add(h_Resmat_ggTot);  
  h_Resmat->Sumw2(); 

  h_4l_c = (TH1*) h_4l ->Clone("h_4l");
  h_4l_c->Add(h_gg4mu);
  h_4l_c->Sumw2();

  h_4l_gen_c = (TH1*) h_4l_gen ->Clone("h_4l_gen");
  h_4l_gen_c->Add(h_gg4mu_gen);
  h_4l_gen_c->Sumw2();

  totalint = h_Resmat->Integral();

  h_Resmat_normTot = (TH2*)h_Resmat->Clone("h_Resmat");
  h_Resmat_normTot->Sumw2();
  
  h_Resmat_normTot->Scale(1/totalint);

  std::cout << "total integral " << totalint << std::endl;
 
  string unc;
  if(xs_qq == 0){
    if(xs_gg == 0) unc = "_st_"; //standard, no variations  
    else if(xs_gg == 1) unc = "_ggp_";
    else if(xs_gg == -1) unc = "_ggm_";
  }
  else if(xs_qq == 1){
    if(xs_gg == 0) unc = "_qqp_"; //standard, no variations  
    else if(xs_gg == 1) unc = "_qqp_ggp_";
    else if(xs_gg == -1) unc = "_qqp_ggm_";
  }
  else if(xs_qq == -1){
    if(xs_gg == 0) unc = "_qqm_"; //standard, no variations  
    else if(xs_gg == 1) unc = "_qqm_ggp_";
    else if(xs_gg == -1) unc = "_qqm_ggm_";
  }

  //  string title = "Response Matrix m_{" + binning + "_" + finalstate + "} " + unc+dataset;   
 

  string matrixNameFile = "ResMat_qqgg_ZZTo" + finalstate + unc + dataset; 
  string matrixNormTotNameFile =  "ResMat_qqgg_normTot_ZZTo" + finalstate + unc + dataset; 
  string histoName_recoFile = "Mass_qqgg_ZZTo" + finalstate + unc + dataset; 
  string histoName_genFile =  "MassGen_qqgg_ZZTo" + finalstate + unc + dataset; 

  h_Resmat->SetTitle(matrixNameFile.c_str());
  h_Resmat_normTot->SetTitle(matrixNormTotNameFile.c_str());
  h_4l_c->SetTitle(histoName_recoFile.c_str());
  h_4l_gen_c->SetTitle(histoName_genFile.c_str());

  output->cd(); 
  h_Resmat->Write(matrixNameFile.c_str());
  h_Resmat_normTot->Write(matrixNormTotNameFile.c_str());
  h_4l_c->Write(histoName_recoFile.c_str());
  h_4l_gen_c->Write(histoName_genFile.c_str());
  
  output->Close();
  
}

void BuildAllMatrices(string dataset = "01"){
  
  for(int i=-1; i<2; i++){
    for(int j=-1; j<2; j++){
      //ResponseMatrix_builder(dataset.c_str(),"4l", i, j);
      ResponseMatrix_builder(dataset.c_str(),"4m", i, j);
      ResponseMatrix_builder(dataset.c_str(),"4e", i, j);
      ResponseMatrix_builder(dataset.c_str(),"2e2m", i, j);
    }
  }
}
