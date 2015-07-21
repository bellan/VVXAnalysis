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

void ResponseMatrix_builder(string var = "Mass", string dataset = "01", string finalstate = "4e", int xs_qq = 0, int xs_gg =0, bool weight = 0)//bool eta = false)
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

  //If the option weight is equal to 1, the macro builds the same distributions, but weighted for the (unfolded_dat/generated_MC ratio)
 
  gROOT->Reset();  
  gROOT->SetStyle("Plain");   
  gStyle->SetOptStat(0);
  string outputName = var+"_test/matrices.root";
  string outputName_W = var+"_test/weightedMatrices.root";
  TFile *output;
  system(("mkdir "+var+"_test").c_str());  

  if(weight ==0)  output = new TFile(outputName.c_str(), "UPDATE");
  else  output = new TFile(outputName_W.c_str(), "UPDATE");

  //Reco samples (response matrices and signal region distributions) //ZZRecoAnalyzer not yet done
  TFile *ggZZTo2e2mu_r = new TFile("../results/ZZRecoAnalyzer_SR/ggTo2e2mu_SMHContinInterf-MCFM67_H125.6.root");
  TFile *ggZZTo4e_r = new TFile("../results/ZZRecoAnalyzer_SR/ggTo4e_SMHContinInterf-MCFM67_H125.6.root");
  TFile *ggZZTo4mu_r = new TFile("../results/ZZRecoAnalyzer_SR/ggTo4mu_SMHContinInterf-MCFM67_H125.6.root");
  TFile *ZZJetsTo4l_r = new TFile("../results/ZZRecoAnalyzer_SR/ZZJetsTo4L.root");
  TFile *ZZTo2e2muJJ_r = new TFile("../results/ZZRecoAnalyzer_SR/ZZTo2e2muJJ_SMHContinInterf_H125.6.root");
  TFile *ZZTo4eJJ_r = new TFile("../results/ZZRecoAnalyzer_SR/ZZTo4eJJ_SMHContinInterf_H125.6.root");
  TFile *ZZTo4muJJ_r = new TFile("../results/ZZRecoAnalyzer_SR/ZZTo4muJJ_SMHContinInterf_H125.6.root");

  //Truth samples (signal definition distributions) 
  TFile *ggZZTo2e2mu_g = new TFile("../results/ZZMCAnalyzer_MC/ggTo2e2mu_SMHContinInterf-MCFM67_H125.6.root");
  TFile *ggZZTo4e_g = new TFile("../results/ZZMCAnalyzer_MC/ggTo4e_SMHContinInterf-MCFM67_H125.6.root");
  TFile *ggZZTo4mu_g = new TFile("../results/ZZMCAnalyzer_MC/ggTo4mu_SMHContinInterf-MCFM67_H125.6.root");
  TFile *ZZJetsTo4l_g = new TFile("../results/ZZMCAnalyzer_MC/ZZJetsTo4L.root");
  TFile *ZZTo2e2muJJ_g = new TFile("../results/ZZMCAnalyzer_MC/ZZTo2e2muJJ_SMHContinInterf_H125.6.root");
  TFile *ZZTo4eJJ_g = new TFile("../results/ZZMCAnalyzer_MC/ZZTo4eJJ_SMHContinInterf_H125.6.root");
  TFile *ZZTo4muJJ_g = new TFile("../results/ZZMCAnalyzer_MC/ZZTo4muJJ_SMHContinInterf_H125.6.root");

  //rescue file
  TFile *rescue = new TFile("../results/ZZRecoAnalyzer_SR/tt.root");
 
  TH2 * h_Resmat;
  TH2 * h_Resmat_4l;
  TH2 * h_Resmat_gg4mu;
  TH2 * h_Resmat_gg4e;
  TH2 * h_Resmat_gg2e2mu;  
  TH2 * h_Resmat_4muJJ;
  TH2 * h_Resmat_4eJJ;
  TH2 * h_Resmat_2e2muJJ; 
  TH2 * h_Resmat_4lTot;
  TH2 * h_Resmat_ggTot;
  TH2 * h_Resmat_JJTot;

  TH1 * h_4l;
  TH1 * h_gg4mu;
  TH1 * h_gg4e;
  TH1 * h_gg2e2mu;
  TH1 * h_4muJJ;
  TH1 * h_4eJJ;
  TH1 * h_2e2muJJ;
 
  TH1 * h_4l_gen;
  TH1 * h_gg4mu_gen;
  TH1 * h_gg4e_gen;
  TH1 * h_gg2e2mu_gen;
  TH1 * h_4muJJ_gen;
  TH1 * h_4eJJ_gen;
  TH1 * h_2e2muJJ_gen;  

  TH1 * h_4l_c;
  TH1 * h_4l_gen_c;
  TH2 * h_Resmat_normTot;

  TH1 *h_safe; 
  TH1 *h_safe_tmp;
  TH2 *h_Resmat_safe; 
  TH2 *h_Resmat_safe_tmp;
  
  string safeMatrixName;
  string safeHistoName;
  string matrixName;
  string histoName_reco;
  string histoName_gen;
  string variable;
  string W;
  int b = 0;  

  if(var == "Mass") {
    variable = var;
    b=9;
  }
  else if(var == "Jets"){
    variable = "Jets_JERCentralSmear";
    b=6;
  }

  if(weight == 0) W = "";  
  else W = "W_";
  
  safeMatrixName = "ResMat_ZZTo2e2m_" + variable +"_"+ W + "01";
  safeHistoName = "ZZTo2e2m_" + variable +"_"+ W +"01"; 
  matrixName = "ResMat_ZZTo" + finalstate + "_" + variable+"_"+ W + dataset;
  histoName_reco = "ZZTo" + finalstate + "_" + variable+"_"+ W +dataset;
  histoName_gen =  "ZZTo" + finalstate + "_" + var + "Gen_" + W  +dataset;
 
  h_Resmat_safe_tmp = (TH2*) rescue->Get(safeMatrixName.c_str()); 
  h_safe_tmp = (TH1*) rescue->Get(safeHistoName.c_str()); 
  h_Resmat_safe = (TH2*)h_Resmat_safe_tmp->Clone(safeMatrixName.c_str()); 
  h_safe = (TH1*)h_safe_tmp->Clone(safeHistoName.c_str());

  for(int k=1; k<b; k++){
    for(int l=1; l<b; l++){
      h_safe->SetBinContent(l,0.);
      h_safe->SetBinError(l,0.);  
      h_Resmat_safe->SetBinContent(l,k,0.);
      h_Resmat_safe->SetBinError(l,k,0.);
    }
  }
 
  float unc_qq = 0; 
  float unc_gg = 0; 
  float totalint = 0; 
   
  h_Resmat_gg4mu = (TH2*) ggZZTo4mu_r->Get(matrixName.c_str()); 
  h_Resmat_gg4e = (TH2*) ggZZTo4e_r->Get(matrixName.c_str()); 
  h_Resmat_gg2e2mu = (TH2*) ggZZTo2e2mu_r->Get(matrixName.c_str()); 
  h_Resmat_4muJJ = (TH2*) ZZTo4muJJ_r->Get(matrixName.c_str()); 
  h_Resmat_4eJJ = (TH2*) ZZTo4eJJ_r->Get(matrixName.c_str()); 
  h_Resmat_2e2muJJ = (TH2*) ZZTo2e2muJJ_r->Get(matrixName.c_str()); 
  h_Resmat_4l = (TH2*) ZZJetsTo4l_r->Get(matrixName.c_str());
 
  h_gg4mu = (TH1*) ggZZTo4mu_r->Get(histoName_reco.c_str()); 
  h_gg4e = (TH1*) ggZZTo4e_r->Get(histoName_reco.c_str()); 
  h_gg2e2mu = (TH1*) ggZZTo2e2mu_r->Get(histoName_reco.c_str()); 
  h_4muJJ = (TH1*) ZZTo4muJJ_r->Get(histoName_reco.c_str()); 
  h_4eJJ = (TH1*) ZZTo4eJJ_r->Get(histoName_reco.c_str()); 
  h_2e2muJJ = (TH1*) ZZTo2e2muJJ_r->Get(histoName_reco.c_str()); 
  h_4l = (TH1*) ZZJetsTo4l_r->Get(histoName_reco.c_str()); 
 
  h_gg4mu_gen = (TH1*) ggZZTo4mu_g->Get(histoName_gen.c_str()); 
  h_gg4e_gen = (TH1*) ggZZTo4e_g->Get(histoName_gen.c_str()); 
  h_gg2e2mu_gen = (TH1*) ggZZTo2e2mu_g->Get(histoName_gen.c_str()); 
  h_4muJJ_gen = (TH1*) ZZTo4muJJ_g->Get(histoName_gen.c_str()); 
  h_4eJJ_gen = (TH1*) ZZTo4eJJ_g->Get(histoName_gen.c_str()); 
  h_2e2muJJ_gen = (TH1*) ZZTo2e2muJJ_g->Get(histoName_gen.c_str()); 
  h_4l_gen = (TH1*) ZZJetsTo4l_g->Get(histoName_gen.c_str()); 
 
  if(h_Resmat_gg4mu == NULL) h_Resmat_gg4mu = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_gg4e == NULL) h_Resmat_gg4e = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_gg2e2mu == NULL) h_Resmat_gg2e2mu = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe"); 
  if(h_Resmat_4muJJ == NULL) h_Resmat_4muJJ = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_4eJJ == NULL) h_Resmat_4eJJ = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_2e2muJJ == NULL) h_Resmat_2e2muJJ = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_gg4mu == NULL) h_gg4mu = (TH1*) h_safe->Clone("h_safe");
  if(h_gg4e == NULL) h_gg4e = (TH1*) h_safe->Clone("h_safe");
  if(h_gg2e2mu == NULL) h_gg2e2mu = (TH1*) h_safe->Clone("h_safe"); 
  if(h_4muJJ == NULL) h_4muJJ = (TH1*) h_safe->Clone("h_safe");
  if(h_4eJJ == NULL) h_4eJJ = (TH1*) h_safe->Clone("h_safe");
  if(h_2e2muJJ == NULL) h_2e2muJJ = (TH1*) h_safe->Clone("h_safe");
  if(h_gg4mu_gen == NULL) h_gg4mu_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_gg4e_gen == NULL) h_gg4e_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_gg2e2mu_gen == NULL) h_gg2e2mu_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_4muJJ_gen == NULL) h_4muJJ_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_4eJJ_gen == NULL) h_4eJJ_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_2e2muJJ_gen == NULL) h_2e2muJJ_gen = (TH1*) h_safe->Clone("h_safe");

  
  h_gg4mu->Add(h_gg4e);
  h_gg4mu->Add(h_gg2e2mu);
  h_gg4mu_gen->Add(h_gg4e_gen);
  h_gg4mu_gen->Add(h_gg2e2mu_gen);
  h_Resmat_gg4mu -> Add(h_Resmat_gg4e, 1);
  h_Resmat_gg4mu -> Add(h_Resmat_gg2e2mu, 1);

  h_4muJJ->Add(h_4eJJ);
  h_4muJJ->Add(h_2e2muJJ);
  h_4muJJ_gen->Add(h_4eJJ_gen);
  h_4muJJ_gen->Add(h_2e2muJJ_gen);

  h_Resmat_4muJJ -> Add(h_Resmat_4eJJ, 1);
  h_Resmat_4muJJ -> Add(h_Resmat_2e2muJJ, 1);
   
  h_Resmat_4lTot = (TH2*)h_Resmat_4l->Clone("h_Resmat_4l");
  h_Resmat_ggTot = (TH2*)h_Resmat_gg4mu->Clone("h_Resmat_gg"); 
  h_Resmat_JJTot = (TH2*)h_Resmat_4muJJ->Clone("h_Resmat_JJ"); 

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
  h_Resmat->Add(h_Resmat_JJTot);
  h_Resmat->Sumw2(); 
  h_4l_c = (TH1*) h_4l ->Clone("h_4l"); 
  h_4l_c->Add(h_gg4mu); 
  h_4l_c->Add(h_4muJJ);
  h_4l_c->Sumw2();
  h_4l_gen_c = (TH1*) h_4l_gen ->Clone("h_4l_gen");
  h_4l_gen_c->Add(h_gg4mu_gen);
  h_4l_gen_c->Add(h_4muJJ_gen);
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
 
  string matrixNameFile = "ResMat_qqggJJ_"+var+"_ZZTo" + finalstate + unc + dataset; 
  string matrixNormTotNameFile =  "ResMat_qqggJJ_"+var+"_normTot_ZZTo" + finalstate + unc + dataset; 
  string histoName_recoFile = var+"_qqggJJ_ZZTo" + finalstate + unc + dataset; 
  string histoName_genFile =  var+"Gen_qqggJJ_ZZTo" + finalstate + unc + dataset; 
  string histoName_recoFile_err = var+"_statErr_qqggJJ_ZZTo" + finalstate + unc + dataset; 

  h_Resmat->SetTitle(matrixNameFile.c_str());
  h_Resmat_normTot->SetTitle(matrixNormTotNameFile.c_str());
  h_4l_c->SetTitle(histoName_recoFile.c_str());
  h_4l_gen_c->SetTitle(histoName_genFile.c_str());
 

  TH1 * h_4l_err = (TH1*) h_4l_c ->Clone("h_4l_c");
  // for(int j=1; j<10; j++){
  //   h_4l_err->SetBinContent(j,0.); 
  // } 

  //h_4l_err->Scale(2);
  float err_stat = 0;
  float bin = 0;
   for(int i =1; i <5; i++){
     err_stat = 0;
     bin = 0;
     bin = h_4l->GetBinContent(i);
     err_stat = sqrt(bin);
     h_4l_err->SetBinError(i,err_stat);
     cout << bin << " " <<err_stat << " " <<   h_4l_err->GetBinError(i)<< " " << h_4l_err->GetEntries() <<endl;
  }

   // h_4l_err->Draw();
   
  output->cd(); 
  h_Resmat->Write(matrixNameFile.c_str());
  h_Resmat_normTot->Write(matrixNormTotNameFile.c_str());
  h_4l_c->Write(histoName_recoFile.c_str());
  h_4l_gen_c->Write(histoName_genFile.c_str());
  h_4l_err->Write(histoName_recoFile_err.c_str());
  output->Close();
  
}

void BuildAllMatrices(string var = "Mass", string dataset = "01"){
  // int p =0;
  //int q = 0;
  for(int p=-1; p<2; p++){
    for(int q=-1; q<2; q++){
      //ResponseMatrix_builder(dataset.c_str(),"4l", p, q);
      ResponseMatrix_builder(var.c_str(),dataset.c_str(),"4m", p, q,0);
      ResponseMatrix_builder(var.c_str(),dataset.c_str(),"4e", p, q,0);
      ResponseMatrix_builder(var.c_str(),dataset.c_str(),"2e2m", p, q,0);
    }
  }
}
void BuildAllWeightedMatrices(string var = "Mass", string dataset = "01"){
  int p =0;
  int q = 0;
  //for(int p=-1; p<2; p++){
  // for(int q=-1; q<2; q++){
      //ResponseMatrix_builder(dataset.c_str(),"4l", p, q);
      ResponseMatrix_builder(var.c_str(),dataset.c_str(),"4m", p, q,1);
      ResponseMatrix_builder(var.c_str(),dataset.c_str(),"4e", p, q,1);
      ResponseMatrix_builder(var.c_str(),dataset.c_str(),"2e2m", p, q,1);
      // }
      //}
}
void PlotMatrix(string var = "Mass",string fs = "4e", string dataset = "01", string unc = "st",bool weight = 0){

  gROOT->Reset();  
  gROOT->SetStyle("Plain");   
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("5.3f");
 
  TFile *file;
  string W;
  string title;
  string xAxis;
  string yAxis;
  string finalstate;
  double max = 0;  
  string fileName = var+"_test/matrices.root";
  string fileName_W = var+"_test/weightedMatrices.root";

  if(fs == "4m") finalstate = "4#mu";
  else if(fs == "2e2m") finalstate = "2e2#mu";
  else finalstate = fs;
  
  if(weight == 0){
    W = ""; 
    file = new TFile(fileName.c_str());
  } 
  else {
    W = "W";
    file = new TFile(fileName_W.c_str());
  }
 
  string matrixName =  "ResMat_qqggJJ_"+var+"_ZZTo" + fs + "_" + unc+ "_" + dataset;
  TH2D *matrix = (TH2D*) file->Get(matrixName.c_str());
  TCanvas *c = new TCanvas("c","c"); 
  
  if(var =="Mass"){
    title = "Response Matrix m_{"+finalstate+"} (dataset: " + dataset+" unc: "+ unc+")";
    xAxis = "reco Njets";
    yAxis = "gen Njets";
    max = matrix->GetBinContent(2,2)/2+3;
  }
  else if(var =="Jets"){
    title = "Response Matrix N jets - "+finalstate+"final state (dataset: " + dataset+" unc: "+ unc+")";
    xAxis = "reco Njets";
    yAxis = "gen Njets"; 
    max = matrix->GetBinContent(1,1)/2;
  }
  
  matrix->SetMaximum(max);
  // matrix->SetTitle(title.c_str());
  matrix->SetTitle("");
  matrix->GetXaxis()->SetTitle(xAxis.c_str());
  matrix->GetYaxis()->SetTitle(yAxis.c_str());
  matrix->Draw("COLZ(TEXT)");
   
 
  string png ="ResMat_qqggJJ_"+var+"_ZZTo" + fs + "_" + unc+ "_" + dataset  +"_" + W + ".png";
  string pdf ="ResMat_qqggJJ_"+var+"_ZZTo" + fs + "_" + unc+ "_" + dataset  +"_" + W + ".pdf";
 
  c->Print(png.c_str());
  c->Print(pdf.c_str());
}
