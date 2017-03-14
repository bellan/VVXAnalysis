#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <string>
#include <sstream>
using std::cout;
using std::endl;

#include <sys/stat.h>
#include <dirent.h>
#include <TROOT.h>
#include <TFile.h>
#include <TStyle.h>
#include "TH1D.h"
#include "TMatrix.h"
#include "TH2.h"
#include "TLegend.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TPad.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#endif

//This macro tests the unfolding procedure on MC samples


std::vector<std::string> Variables   = {"Mass","nJets","nJets_Central","Mjj","Mjj_Central","Deta","Deta_Central","PtJet1","PtJet2","EtaJet1","EtaJet2","dRZZ","PtZZ"};

void Unfold_MCtest(string var = "Mass",string finalstate = "4e", bool bayes = 0, int it = 4, bool MadMatrix =1, bool MadDistribution = 1, bool FullSample = 0,bool tightregion = 0, string date = "test")
{
  
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  std::string SavePage = "~/www/PlotsVV/13TeV/";


  struct stat st;

  //  if(stat((""+SavePage+"+date+"/").c_str(),&st) != 0){
    system(("mkdir "+SavePage+date+"/").c_str()); 
    system(("cp "+SavePage+"index.php "+SavePage+date+"/").c_str()); 
    //}
    //  if(satat((""+SavePage+date+"/"+ var+"/").c_str(),&st) != 0){
    system(("mkdir "+SavePage+date+"/"+ var+"/").c_str()); 
    system(("cp "+SavePage+"index.php "+SavePage+date+"/"+ var+"/").c_str());
    //  }
    //  if(stat((""+SavePage+"+date+"/"+ var+"/MCTests").c_str(),&st) != 0){
    system(("mkdir "+SavePage+date+"/"+ var+"/MCTests").c_str()); 
    system(("cp "+SavePage+"index.php "+SavePage+date+"/"+ var+"/MCTests").c_str());
    //  }

  setTDRStyle();
  string fs;
  if(finalstate == "4m")        fs = "4#mu";
  else if(finalstate == "2e2m") fs = "2e2#mu";
  else fs = finalstate;
  
  string label = fs + " channel";
  int iPeriod = 4; 
  int iPos = 11; 
  writeExtraText = true;    
  extraText  = "Simulation";
  extraText2  = label.c_str();
  string tightfr;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
   
  TH1 * h_true; 
  TH1 * h_measured; 
  TH1 * h_measured_unf; 
  TH1 * h_true_unf; 
  TH2 * h_Resmat;
  TH1 * hReco;
 
  string filePath = "../";
  string fileName_Mad = filePath + var +"_test/matrices"+tightfr+"_Mad.root";
  string fileName_Pow = filePath + var +"_test/matrices"+tightfr+"_Pow.root"; 

  TFile *file_mad = new TFile(fileName_Mad.c_str());
  TFile *file_pow = new TFile(fileName_Pow.c_str());
 
  string matrixName; 
  string histoName; 
  string histoNamegen; 
  string histoNameGen_unf;
  string histoName_unf;

 if(FullSample == 0){
   matrixName       = "ResMat_qqggJJ_"+var+"_ZZTo" + finalstate + "_st_0";
   histoName        = var+"_qqggJJ_ZZTo" + finalstate + "_st_0";
   histoNamegen     = var+"Gen_qqggJJ_ZZTo" + finalstate + "_st_0";
   histoName_unf    = var+"_qqggJJ_ZZTo" + finalstate + "_st_1"; 
   histoNameGen_unf = var+"Gen_qqggJJ_ZZTo" + finalstate + "_st_1";
 }
 
 else{
   matrixName = "ResMat_qqggJJ_"+var+"_ZZTo" + finalstate + "_st_01";
   histoName = var+"_qqggJJ_ZZTo" + finalstate + "_st_01";
   histoNamegen = var+"Gen_qqggJJ_ZZTo" + finalstate + "_st_01";
   histoName_unf =  var+"_qqggJJ_ZZTo" + finalstate + "_st_01"; 
   histoNameGen_unf = var+"Gen_qqggJJ_ZZTo" + finalstate + "_st_01";
 } 
  
 if(MadMatrix == 1){
   h_measured = (TH1*) file_mad->Get(histoName.c_str());  
   h_true     = (TH1*) file_mad->Get(histoNamegen.c_str());
   h_Resmat   = (TH2*) file_mad->Get(matrixName.c_str()); 
 }

 else{
   h_measured = (TH1*) file_pow->Get(histoName.c_str());  
   h_true     = (TH1*) file_pow->Get(histoNamegen.c_str());
   h_Resmat   = (TH2*) file_pow->Get(matrixName.c_str()); 
 }

 if(MadDistribution == 1){
   h_measured_unf = (TH1*) file_mad->Get(histoName_unf.c_str());  
   h_true_unf     = (TH1*) file_mad->Get(histoNameGen_unf.c_str());
 }
 
 else{
   h_measured_unf = (TH1*) file_pow->Get(histoName_unf.c_str());  
   h_true_unf     = (TH1*) file_pow->Get(histoNameGen_unf.c_str());
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

  string XaxisTitle;
  if(var == "Mass") XaxisTitle = "m_{4l} [GeV]"; 
  else if(var == "nJets")  XaxisTitle = "N jets (|#eta^{jet}| < 4.7)";
  else if(var == "nJets_Central")  XaxisTitle = "N jets (|#eta^{jet}| < 2.4)"; 
  else if(var == "Mjj") XaxisTitle ="m_{jj} (|#eta^{jet}| < 4.7) [GeV]";
  else if(var == "Mjj_Central") XaxisTitle ="m_{jj} (|#eta^{jet}| < 2.4) [GeV]";
  else if(var == "Deta") XaxisTitle ="#Delta#eta_{jj} (|#eta^{jet}| < 4.7)";
  else if(var == "Deta_Central") XaxisTitle ="#Delta#eta_{jj} (|#eta^{jet}| < 2.4)";
  else if(var =="PtJet1")  XaxisTitle = "p_{T}^{jet1} [GeV]";
  else if(var =="PtJet2")  XaxisTitle = "p_{T}^{jet2} [GeV]";
  else if(var =="EtaJet1")  XaxisTitle = "|#eta^{jet1}|";
  else if(var =="EtaJet2")  XaxisTitle = "|#eta^{jet2}|";
  else if(var =="dRZZ")  XaxisTitle = "#DeltaR(Z_{1},Z_{2})";
 
  string YaxisTitle = "Events";
  string YaxisTitle2 = "Unfolded/True";
  TCanvas *c = new TCanvas ("c","c");
  TLegend *leg = new TLegend(0.80,0.65,0.60,0.85); 
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42); 

  c->cd();
  
  TPad  *pad1 = new TPad("pad1","", 0., 0.22, 1.0, 1.0);
  pad1->SetTopMargin (0.10);
  pad1->SetRightMargin (0.10);
  pad1->SetLeftMargin (0.10);
  pad1->Draw();
  
  c->cd(); 
  TPad  *pad2 = new TPad("pad2", "", 0., 0.0,  1.0, 0.28);
  pad2->SetTopMargin (0.10);
  pad2->SetRightMargin (0.10);
  pad2->SetLeftMargin (0.10);
  pad2->SetBottomMargin (0.35);
  pad2->Draw(); 
    
  pad1->cd();
  h_true_unf->SetTitle("");
  //  h_true_unf->GetXaxis()->SetRange(0,25);
  h_true_unf->GetXaxis()->SetTitle(XaxisTitle.c_str()); 
  h_true_unf->GetYaxis()->SetTitle(YaxisTitle.c_str());
  
  float max = 0;

  if(var == "Mjj" ||var == "Mjj_Central" ||var == "Deta" ||var == "Deta_Central"){ 
    //max = h_true_unf->GetMaximum()+2;
    //h_true_unf->SetMaximum(max);
  }
  else if(var == "PtJet1"||var == "PtJet2"){
    //max = hReco->GetMaximum()+5;
    //h_true_unf->SetMaximum(max);
  }
  else if(var == "EtaJet1"||var == "EtaJet2"){
    // max = hReco->GetMaximum()+5;
    //h_true_unf->SetMaximum(max);
  }
  h_true_unf->SetMinimum(0); 
  h_true_unf->SetMarkerColor(8);
  h_true_unf->SetMarkerStyle(8);
  h_true_unf->SetLineColor(8); 
  h_true_unf->SetLineWidth(1);
  h_true_unf->Draw("HIST");
  h_true_unf->GetXaxis()->SetLabelOffset(0.5); 
  max = h_true_unf->GetMaximum()*1.5;
  h_true_unf->SetMaximum(max);
   h_true_unf->GetYaxis()->SetTitleSize(0.04); 
  h_true_unf->Draw("HIST");
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
  hReco_r->GetXaxis()->SetLabelSize(0.13); 
  if(var == "nJets" || var == "nJets_Central"){
   hReco_r->GetXaxis()->SetBinLabel(1,"0");
   hReco_r->GetXaxis()->SetBinLabel(2,"1");
   hReco_r->GetXaxis()->SetBinLabel(3,"2");
   hReco_r->GetXaxis()->SetBinLabel(4,">2");
   hReco_r->GetXaxis()->SetLabelFont(42);
    //h_true->GetXaxis()->SetLabelOffset(0.02);
   hReco_r->GetXaxis()->SetLabelSize(0.17);
  }

  hReco_r->GetYaxis()->SetLabelSize(0.10);  
  hReco_r->GetYaxis()->SetNdivisions(306); 
  hReco_r->GetXaxis()->SetLabelOffset(0.05);
  hReco_r->GetXaxis()->SetTitleOffset(1.3); 
  hReco_r->GetYaxis()->SetTitleOffset(0.4);
  hReco_r->GetYaxis()->SetTitleSize(0.10);
  hReco_r->GetXaxis()->SetTitleSize(0.13);
  hReco_r->SetMarkerColor(1);
  hReco_r->SetLineColor(1); 
  hReco_r->GetXaxis()->SetTitle(XaxisTitle.c_str()); 
  hReco_r->GetYaxis()->SetTitle(YaxisTitle2.c_str());
  hReco_r-> SetMaximum(2.); 
  hReco_r-> SetMinimum(0.);
  hReco_r->Draw("E");

  lumiTextSize     = 0.4;
  cmsTextSize      = 0.48;
  extraOverCmsTextSize  = 0.80;//0.63;
  CMS_lumi(pad1, iPeriod, iPos );
  //pad1->Update();
  //pad1->RedrawAxis();
  //pad1->GetFrame()->Draw();
  
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

 string Algo = "";
 if(bayes) Algo = "_Bayes";
 else Algo = "_SVD";

 pdf = SavePage+date+"/"+ var+"/MCTests/"+var+"_ZZTo" + finalstate + "_"+ matrix + "_" + distr + "_" + sample + tightfr +Algo+"_"+std::to_string(it)+".pdf";
 png = SavePage+date+"/"+ var+"/MCTests/"+var+"_ZZTo" + finalstate + "_"+ matrix + "_" + distr + "_" + sample + tightfr +Algo+"_"+std::to_string(it)+".png";
  c->Print(pdf.c_str());
  c->Print(png.c_str());
  
  //output->cd();   
  //hReco->Write(UnfHistoName.c_str());
  //h_measured_unf->Write(HistoName.c_str());
  //h_true_unf->Write(TrueHistoName.c_str());
  file_mad->Close();
  file_pow->Close();



}

void MakeAllTest(string var = "Mass",string finalstate = "4e",bool bayes = 0, int it = 4, bool tightregion = 0, string date = "test"){
  
  Unfold_MCtest(var.c_str(),finalstate.c_str(), bayes, it,1, 1, 1,tightregion,date);
  Unfold_MCtest(var.c_str(),finalstate.c_str(), bayes, it,1, 1, 0,tightregion,date);
  Unfold_MCtest(var.c_str(),finalstate.c_str(), bayes, it,1, 0, 1,tightregion,date);
  Unfold_MCtest(var.c_str(),finalstate.c_str(), bayes, it,0, 1, 1,tightregion,date);
  Unfold_MCtest(var.c_str(),finalstate.c_str(), bayes, it,0, 0, 1,tightregion,date);
  Unfold_MCtest(var.c_str(),finalstate.c_str(), bayes, it,0, 0, 0,tightregion,date); 
 
}

void MakeAllFinalStates (string var = "Mass",bool bayes = 0, int it = 4, bool tightregion = 0, string date = "test"){
  MakeAllTest(var.c_str(),"4e",bayes,it,tightregion,date);
  MakeAllTest(var.c_str(),"4m",bayes,it,tightregion,date);
  MakeAllTest(var.c_str(),"2e2m",bayes,it,tightregion,date);
}


void MakeAllTests (bool bayes = 0, int it = 4, string date = "test"){
  foreach(const std::string &var, Variables){
    MakeAllFinalStates(var,bayes,it,0,date);
    MakeAllFinalStates(var,bayes,it,1,date);
  }
}

void MakeAllOfficialTests (string date = "test"){
  MakeAllFinalStates("Mass",1,4,0,date);
  MakeAllFinalStates("Mass",1,4,1,date);
  MakeAllFinalStates("nJets",1,4,0,date);
  MakeAllFinalStates("nJets",1,4,1,date);
  MakeAllFinalStates("nJets_Central",1,4,0,date);
  MakeAllFinalStates("nJets_Central",1,4,1,date);
  MakeAllFinalStates("Mjj_Central",1,2,0,date);
  MakeAllFinalStates("Mjj_Central",1,2,1,date);  
  MakeAllFinalStates("Deta_Central",1,2,0,date);
  MakeAllFinalStates("Deta_Central",1,2,1,date); 
  MakeAllFinalStates("Mjj",1,2,0,date);
  MakeAllFinalStates("Mjj",1,2,1,date);  
  MakeAllFinalStates("Deta",1,2,0,date);
  MakeAllFinalStates("Deta",1,2,1,date);
  MakeAllFinalStates("PtJet1",1,4,0,date);
  MakeAllFinalStates("PtJet1",1,4,1,date);
  MakeAllFinalStates("PtJet2",1,4,0,date);
  MakeAllFinalStates("PtJet2",1,4,1,date);
  MakeAllFinalStates("EtaJet1",1,3,0,date);
  MakeAllFinalStates("EtaJet1",1,3,1,date);
  MakeAllFinalStates("EtaJet2",1,3,0,date);
  MakeAllFinalStates("EtaJet2",1,3,1,date);
  MakeAllFinalStates("PtZZ",1,4,0,date);
  MakeAllFinalStates("PtZZ",1,4,1,date);
  MakeAllFinalStates("dRZZ",1,4,0,date);
  MakeAllFinalStates("dRZZ",1,4,1,date);

  MakeAllFinalStates("nJets",0,2,1,date);
  MakeAllFinalStates("Mjj",0,2,0,date);
  MakeAllFinalStates("PtJet1",0,2,0,date);
}


#ifndef __CINT__
int main () { Unfold_MCtest(); return 0; }  // Main program when run stand-alone
#endif
