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

void ResponseMatrix_builder(string dataset = "01", string finalstate = "4e", string unc = "JESDown")//bool eta = false)
{
  //This macro builds the response matrix (measured X truth), the "measured" histo and the "truth" histo, that are the projections of "response" onto the X-axis and Y-axis respectively, but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and in "truth" for unmeasured events (inefficiency). using  and the TH1 measured and truth 

  //dataset == 01 full dataset
  //dataset == 0 first half of the dataset
  //dataset == 1 second half of the dataset

   
  gROOT->Reset();  
  gROOT->SetStyle("Plain");   
  gStyle->SetOptStat(0);

  TFile *output = new TFile("Jets_test/matrices_JESJER.root", "UPDATE");

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
  
  string safeMatrixName = "ResMat_ZZTo2e2m_Jets_JERCentralSmear_01";
  string safeHistoName = "ZZTo2e2m_Jets_JERCentralSmear_01"; 
  
  h_Resmat_safe_tmp = (TH2*) rescue->Get(safeMatrixName.c_str()); 
  h_safe_tmp = (TH1*) rescue->Get(safeHistoName.c_str()); 
  h_Resmat_safe = (TH2*)h_Resmat_safe_tmp->Clone(safeMatrixName.c_str()); 
  h_safe = (TH1*)h_safe_tmp->Clone(safeHistoName.c_str());

  for(int k=1; k<6; k++){
    for(int l=1; l<6; l++){
      h_safe->SetBinContent(l,0.);
      //h_safe->SetBinError(l,0.);  
      h_Resmat_safe->SetBinContent(l,k,0.);
      // h_Resmat_safe->SetBinError(l,k,0.);
    }
  }
 
   float totalint = 0; 
  
  string matrixName = "ResMat_ZZTo" + finalstate + "_Jets_"+unc+"Smear_" + dataset;
  string histoName_reco = "ZZTo" + finalstate + "_Jets_"+unc+"Smear_" + dataset;
  string histoName_gen =  "ZZTo" + finalstate + "_JetsGen_"+ dataset;
 
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
   cout << "7" << endl;
  h_Resmat_4lTot = (TH2*)h_Resmat_4l->Clone("h_Resmat_4l");
  h_Resmat_ggTot = (TH2*)h_Resmat_gg4mu->Clone("h_Resmat_gg"); 
  h_Resmat_JJTot = (TH2*)h_Resmat_4muJJ->Clone("h_Resmat_JJ"); 

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

  //  string title = "Response Matrix m_{" + binning + "_" + finalstate + "} " + unc+dataset;   
 
  string matrixNameFile = "ResMat_qqggJJ_Jets_ZZTo" + finalstate + "_" + unc + "_" + dataset; 
  string matrixNormTotNameFile =  "ResMat_qqggJJ_Jets_normTot_ZZTo" + finalstate + "_" + unc + "_" + dataset; 
  string histoName_recoFile = "Jets_qqggJJ_ZZTo" + finalstate + "_" + unc + "_" + dataset; 
  string histoName_genFile =  "JetsGen_qqggJJ_ZZTo" + finalstate + "_" + unc + "_" + dataset; 
  string histoName_recoFile_err = "Jets_statErr_qqggJJ_ZZTo" + finalstate + "_" + unc + "_" + dataset; 

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

void BuildAllMatrices(string dataset = "01", string finalstate = "4e"){

  ResponseMatrix_builder(dataset.c_str(),finalstate.c_str(), "JESDown");
  ResponseMatrix_builder(dataset.c_str(),finalstate.c_str(), "JESUp");
  ResponseMatrix_builder(dataset.c_str(),finalstate.c_str(), "JERDown");
  ResponseMatrix_builder(dataset.c_str(),finalstate.c_str(), "JERUp");
  
}

void BuildAllFinalstate(string dataset = "01") {

  BuildAllMatrices(dataset.c_str(),"4e");
  BuildAllMatrices(dataset.c_str(),"4m");
  BuildAllMatrices(dataset.c_str(),"2e2m");
  
}

void PlotMatrix(string finalstate = "4e", string dataset = "01", string unc = "JESDown"){

  gROOT->Reset();  
  gROOT->SetStyle("Plain");   
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("5.3f");

  TFile *file = new TFile("Jets/matrices_JESJER.root"); 
  string matrixName =  "ResMat_qqggJJ_Jets_ZZTo" + finalstate + "_" + unc+ "_" + dataset;
  string title = "Response Matrix N jets - "+finalstate+"final state (dataset: " + dataset+" unc: "+ unc+")";
  string xAxis = "reco Njets";
  string yAxis = "gen Njets";


  TH2D * matrix = (TH2D*) file->Get(matrixName.c_str());
  TCanvas *c = new TCanvas("c","c"); 
  double max = matrix->GetBinContent(1,1)/2;
  matrix->SetMaximum(max);
  // matrix->SetTitle(title.c_str());
  matrix->SetTitle("");
  matrix->GetXaxis()->SetTitle(xAxis.c_str());
  matrix->GetYaxis()->SetTitle(yAxis.c_str());
  matrix->Draw("COLZ(TEXT)");
  
  string png ="ResMat_qqggJJ_Jets_ZZTo" + finalstate + "_" + unc+ "_" + dataset  +".png";
  string eps ="ResMat_qqggJJ_Jets_ZZTo" + finalstate + "_" + unc+ "_" + dataset  +".eps";
 
  c->Print(eps.c_str());

}
