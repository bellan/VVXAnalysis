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

void PurityAndStability(string var = "Mass", string dataset = "01", string finalstate = "4e", bool Madgraph =1)//bool eta = false)
{ 

  //dataset == 01 full dataset
  //dataset == 0 first half of the dataset
  //dataset == 1 second half of the dataset
  
  gROOT->Reset();  
  gROOT->SetStyle("Plain");   
  gStyle->SetOptStat(0); 
  TFile *output; 
  TFile *matrixFile; 
  string matrixName = "ResMat_qqggJJ_"+var+"_ZZTo"+finalstate+"_st_"+dataset;
  string histoName = var+"_qqggJJ_ZZTo"+finalstate+"_st_"+dataset;
  string fileName = var+"_test/MatrixPurityStability.root"; 
  string fileName_Pow =  var+"_test/MatrixPurityStability_Pow.root";
  string matrixFileName = var+"_test/matrices.root";
  string matrixFileName_Pow = var+"_test/matrices_Pow.root";

  if(Madgraph ==1){
    output = new TFile(fileName.c_str(), "UPDATE"); 
    matrixFile = new TFile(matrixFileName.c_str());
  }
  
  else{
    output = new TFile(fileName_Pow.c_str(), "UPDATE"); 
    matrixFile = new TFile(matrixFileName_Pow.c_str());
  }

  TH2 * h_Resmat; 
  TH1 * h_purity;
  TH1 * h_stability; 
  float binContent = 0; 
  float binError = 0;
  float all_reco = 0;
  float all_gen = 0;
  float p_err =0;  
  float s_err =0;  
  double reco_err =0;  
  double gen_err =0;  
  int b = 0;

  h_Resmat = (TH2*) matrixFile->Get(matrixName.c_str()); 
  h_purity = (TH1*) matrixFile->Get(histoName.c_str()); 
  h_stability = (TH1*) matrixFile->Get(histoName.c_str()); 
  
  if(var == "Mass") b=9;
 
  else if(var == "Jets") b=5;
   
  for(int i =1; i<b; i++){
    reco_err=0;
    gen_err=0;
    binContent=0;  
    binError=0;    
    all_reco=0;
    all_gen=0; 
    p_err=0;
    s_err =0;
    
    binError = h_Resmat->GetBinError(i,i);
    binContent = h_Resmat->GetBinContent(i,i);
    all_reco = h_Resmat->IntegralAndError(i,i,1,b,reco_err);
    all_gen = h_Resmat->IntegralAndError(1,b,i,i,gen_err);
    p_err = sqrt((binError/binContent)*(binError/binContent) + (reco_err/all_reco)*(reco_err/all_reco))*(binContent/all_reco);    
    s_err = sqrt((binError/binContent)*(binError/binContent) + (gen_err/all_gen)*(gen_err/all_gen))*(binContent/all_gen);    


    std::cout <<"bin " << i <<  " all_reco  = " << all_reco << " +- " << reco_err << " p= "<<  binContent/all_reco << " +- " << p_err << std::endl;
    std::cout <<"bin " << i <<" all_gen  = " << all_gen <<" +- " << gen_err << " s= " << binContent/all_gen <<" +- " << s_err  << std::endl;
    
    if(all_reco!=0){
      h_purity->SetBinContent(i,binContent/all_reco);
      h_purity->SetBinError(i,p_err);
    }
    else h_purity->SetBinContent(i,0.);
    if(all_gen!=0) {
      h_stability->SetBinContent(i,binContent/all_gen);
      h_stability->SetBinError(i,s_err);
    }
    else h_stability->SetBinContent(i,0.);
  }

  h_purity->SetMaximum(1.);
  h_purity->SetMinimum(0.); 
  h_stability->SetMaximum(1.);
  h_stability->SetMinimum(0.);

  string purityName = "hPurity_ZZTo" + finalstate + "_" + dataset;
  string stabilityName = "hStability_ZZTo" + finalstate + "_" + dataset;

 output->cd(); 
 h_purity->Write(purityName.c_str()); 
 h_stability->Write(stabilityName.c_str()); 
 output->Close();

}
void MakeAllFinalStates(string var = "Mass",string dataset = "01", bool Madgraph = 1){
 
  //PurityAndStability(dataset.c_str(),"4l",Madgraph);
  PurityAndStability(var.c_str(),dataset.c_str(),"4m",Madgraph);
  PurityAndStability(var.c_str(),dataset.c_str(),"4e",Madgraph);
  PurityAndStability(var.c_str(),dataset.c_str(),"2e2m",Madgraph);

}

void Plot(string var = "Mass", string finalstate = "4m"){
  
  gROOT->Reset();  
  gROOT->SetStyle("Plain");   
  gStyle->SetOptStat(0); 

  string fileNameMad = var+"_test/MatrixPurityStability.root";
  string fileNamePow =  var+"_test/MatrixPurityStability_Pow.root";
  
  TFile * madgraph = new TFile(fileNameMad.c_str()); 
  TFile * powheg = new TFile(fileNamePow.c_str()); 
  
  string pName = "hPurity_ZZTo"+finalstate+"_01";
  string sName = "hStability_ZZTo"+finalstate+"_01";
 
  TH1 * p_mad = (TH1*) madgraph->Get(pName.c_str());
  TH1 * s_mad = (TH1*) madgraph->Get(sName.c_str());
  TH1 * p_pow = (TH1*) powheg->Get(pName.c_str());
  TH1 * s_pow = (TH1*) powheg->Get(sName.c_str());

  TCanvas *c_p = new TCanvas ("c_p","c_p");
  TLegend *leg_p = new TLegend(0.65,0.35,0.6,0.43); 
  leg_p->SetFillColor(kWhite);
  leg_p->SetBorderSize(0);
  leg_p->SetTextSize(0.03);
 
  //  string pTitle = "Purity distribution of ZZ#rightarrow" + finalstate + " final state"; 
  //string sTitle = "Stability distribution of ZZ#rightarrow" + finalstate + " final state";
  string sTitle = "";  
  string pTitle = "";
  string XaxisTitle; 
  string fs;

  if(finalstate == "4m") fs = "4#mu";
  else if(finalstate == "2e2m") fs = "2e2#mu";
  else finalstate = fs;
 
  if(var =="Mass"){
    XaxisTitle = "m_{"+fs+"}";
  }
  else if(var =="Jets"){
    XaxisTitle = "Njets ("+fs+"-final state)";
  }
 
  p_mad->Draw("E");  
  p_mad->SetLineColor(1);
  p_mad->SetTitle(pTitle.c_str());
  //hReco->GetXaxis()->SetRange(0,25);
  p_mad->GetXaxis()->SetTitle(XaxisTitle.c_str());
 
  p_pow->Draw("E SAME");
  p_pow->SetLineColor(2);
  
  leg_p->AddEntry(p_mad,"Madgraph set","l"); 
  leg_p->Draw(); 
  leg_p->AddEntry(p_pow,"Powheg set","l"); 
  leg_p->Draw("SAME");
  
  TCanvas *c_s = new TCanvas ("c_s","c_s");
  TLegend *leg_s = new TLegend(0.65,0.35,0.6,0.43); 
  leg_s->SetFillColor(kWhite);
  leg_s->SetBorderSize(0);
  leg_s->SetTextSize(0.03);

  s_mad->Draw("E");  
  s_mad->SetLineColor(1);
  s_mad->SetTitle(sTitle.c_str());
  //hReco->GetXaxis()->SetRange(0,25);
  s_mad->GetXaxis()->SetTitle(XaxisTitle.c_str());
 
  s_pow->Draw("E SAME");
  s_pow->SetLineColor(2);
  
  leg_s->AddEntry(s_mad,"Madgraph set","l"); 
  leg_s->Draw(); 
  leg_s->AddEntry(s_pow,"Powheg set","l"); 
  leg_s->Draw("SAME");

  string p_png = "purity_MadPow_"+finalstate+".png";
  string s_png = "stability_MadPow_"+finalstate+".png";

  // c_p->Print(p_png.c_str());
  // c_s->Print(s_png.c_str());
  
}

void PlotPurityAndStability(string var = "Mass",string finalstate = "4m", bool madgraph =1){
  
  gROOT->Reset();  
  gROOT->SetStyle("Plain");   
  gStyle->SetOptStat(0); 

  string fileNameMad = var+"_test/MatrixPurityStability.root";
  string fileNamePow =  var+"_test/MatrixPurityStability_Pow.root";
  
  TFile *file;
  if(madgraph ==1) file = new TFile(fileNameMad.c_str()); 
  else file= new TFile(fileNamePow.c_str()); 
  //TFile * madgraph = new TFile("MatrixPurityStability.root"); 
  //TFile * powheg = new TFile("MatrixPurityStability_PowhegSet.root"); 

  string pName = "hPurity_ZZTo"+finalstate+"_01";
  string sName = "hStability_ZZTo"+finalstate+"_01";
 
  TH1 * p_mad = (TH1*) file->Get(pName.c_str());
  TH1 * s_mad = (TH1*) file->Get(sName.c_str());
  TH1 * p_pow = (TH1*) file->Get(pName.c_str());
  TH1 * s_pow = (TH1*) file->Get(sName.c_str());

  TCanvas *c_p = new TCanvas ("c_p","c_p");
  TLegend *leg_p = new TLegend(0.70,0.35,0.6,0.45); 
  leg_p->SetFillColor(kWhite);
  leg_p->SetBorderSize(0);
  leg_p->SetTextSize(0.04);
  leg_p->SetTextFont(42); 

  //string pTitle = "Purity distribution of ZZ#rightarrow" + finalstate + " final state";
  //string sTitle = "";  
  // string pTitle = "";
  string XaxisTitle; 
  string fs;

  if(finalstate == "4m") fs = "4#mu";
  else if(finalstate == "2e2m") fs = "2e2#mu";
  else finalstate = fs;
 
  if(var =="Mass"){
    XaxisTitle = "m_{"+fs+"}";
  }
  else if(var =="Jets"){
    XaxisTitle = "Njets ("+fs+"-final state)";
  }

  p_mad->Draw("E");  
  p_mad->SetLineColor(1);
  p_mad->SetMarkerStyle(8);
  //p_mad->SetTitle(pTitle.c_str());
  p_mad->SetTitle("");
  //hReco->GetXaxis()->SetRange(0,25);
  p_mad->GetXaxis()->SetTitle(XaxisTitle.c_str());
 
  s_mad->Draw("E SAME");
  s_mad->SetLineColor(2);
  s_mad->SetMarkerStyle(8);
  s_mad->SetMarkerColor(2);

  leg_p->AddEntry(p_mad,"Purity","p"); 
  leg_p->Draw(); 
  leg_p->AddEntry(s_mad,"Stability","p"); 
  leg_p->Draw("SAME");
  
  string p_pdf = "PurityStability_"+finalstate+"_Njets.pdf";

  // c_p->Print(p_pdf.c_str());
   
}
