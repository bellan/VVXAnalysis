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
#include "TLine.h"
#include "TH2.h"
#include "TLegend.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
#include "TCanvas.h"
//#include <boost/filesystem.hpp>
#endif
//==============================================================================
// Global definitions
//==============================================================================
TH1 * h_true; 
TH1 * h_measured; 
TH1 * h_measured_unf; 
TH2 * h_Resmat;
TH1 * h_unfolded;
TH2 * h_Resmat_p;
TH1 * h_true_p;
TH1 * h_measured_p;
TH1 * h_unfolded_p; 
TH2 * h_Resmat_m;
TH1 * h_true_m;
TH1 * h_measured_m;
TH1 * h_unfolded_m;
TH1 * h_data_p;
TH1 * h_data_m;

TMatrixD * cov;
TVectorD * Vcov_stat; 
TVectorD * Vcov_unf; 

TFile *data;
TFile *matrix;
TFile * output;

std::string filePath;
std::string matrixFileName;
std::string dataFileName;
std::string outputFileName;
std::string matrixName; 
std::string histoName; 
std::string histoNameGen;
std::string histoName_unf;
std::string unfHistoName;
std::string recoHistoName;
std::string trueHistoName;
std::string matrix_p;
std::string histoReco_p;
std::string histoGen_p;
std::string matrix_m;
std::string histoReco_m;
std::string histoGen_m; 
std::string UnfHistoName_p;
std::string UnfHistoName_m;
std::string data_p;
std::string data_m;
std::string XaxisTitle;
std::string YaxisTitle;

/*To unfold Mass(Jets) distributions, first download the RooUnfold package following these instructions:
 http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html
Then create the folder "UnfoldFolder" and in ../RooUnfold-1.1.1/RooUnfold (after "make"):
root
.L ../../UnfoldingMacros/Unfold.cpp+ 
Unfold_data_All("Mass") (to unfold of all data distributions)
DoAllRatios("Mass") (to build Unfdata/GenMC ratio => run again run.py for all samples in order to create weighted distributions for data/MC systematic uncertainty)
DoAllSystematics("Mass") (to create distributions for all the systematic sources)
(or, simply: AllYouNeed())*/

void Unfold_data(std::string var = "Mass", std::string fs = "4e"){
  
 
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  system("mkdir ../../UnfoldingMacros/UnfoldFolder/");  
  filePath = "../../";
  matrixFileName = filePath + var + "_test/matrices.root";
  dataFileName = filePath +var + "_test/DataToUnfold.root";
  string outputFolderName =  filePath+"UnfoldingMacros/UnfoldFolder";
  outputFileName = filePath+"UnfoldingMacros/UnfoldFolder/UnfoldData_"+ var + ".root";
   
  matrix = new TFile(matrixFileName.c_str());
  data = new TFile(dataFileName.c_str());
  output = new TFile(outputFileName.c_str(), "UPDATE");
  
  matrixName = "ResMat_qqggJJ_" + var + "_ZZTo" + fs + "_st_01";
  histoName = var +"_qqggJJ_ZZTo" + fs + "_st_01";
  histoNameGen = var + "Gen_qqggJJ_ZZTo" + fs + "_st_01";
  histoName_unf = "DataminusBkg_" + var + "_ZZTo"+fs;
 
  h_measured = (TH1*) matrix->Get(histoName.c_str());  
  h_true = (TH1*) matrix->Get(histoNameGen.c_str());
  h_Resmat = (TH2*)matrix->Get(matrixName.c_str()); 
  h_measured_unf = (TH1*) data->Get(histoName_unf.c_str());  

  // h_measured->Draw();
  // h_true ->Draw("same");
  // h_Resmat->Draw("same");
  RooUnfoldResponse response(h_measured, h_true, h_Resmat, "response", "response"); 
  RooUnfoldSvd   unfold_svd(&response, h_measured_unf, 4); 
  RooUnfoldBayes unfold_bayes(&response, h_measured_unf,2);  
    
  if(var == "Mass") h_unfolded= (TH1*) unfold_svd.Hreco(RooUnfold::kCovariance);
  else if(var == "Jets") h_unfolded= (TH1*) unfold_bayes.Hreco(RooUnfold::kCovariance);
 
  unfHistoName = "ZZTo"+ fs +"_" + var;
  recoHistoName = "ZZTo"+ fs +"_"+var+"_RECO_MC";
  trueHistoName = "ZZTo"+ fs +"_"+var+"_GEN";
  string finalstate;
  if(fs == "4m") finalstate = "4#mu";
  else if(fs == "2e2m") finalstate = "2e2#mu";
  else finalstate = fs;

  if(var == "Mass") XaxisTitle = "m_{" + finalstate + "}"; 
  else if(var == "Jets")  XaxisTitle = "N Jets (" + finalstate + " final state)";

  YaxisTitle = "Events";
  string YaxisTitle2 = "Unfolded/True";

  TCanvas *c = new TCanvas ("c","c");
  TLegend *leg = new TLegend(0.85,0.65,0.65,0.85); 
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
  h_true->SetTitle("");
  h_true->GetXaxis()->SetRange(0,25);
  h_true->GetXaxis()->SetTitle(XaxisTitle.c_str()); 
  h_true->GetYaxis()->SetTitle(YaxisTitle.c_str());

  float max = 0;
  if(var == "Mass") max =  h_unfolded->GetBinContent(2)+20; 
  else if(var == "Jets")  max =  h_unfolded->GetBinContent(1)+40; 

  h_true->SetMaximum(max);
  h_true->SetMinimum(0); 
  h_true->SetMarkerColor(8);
  h_true->SetMarkerStyle(8);
  h_true->SetLineColor(8); 
  h_true->SetLineWidth(1);
  h_true->Draw("HIST E");
  
  h_unfolded->SetLineColor(kBlue);
  h_unfolded->SetLineWidth(1);
  h_unfolded->SetMarkerColor(kBlue);
  h_unfolded->SetMarkerStyle(8);
  h_unfolded->Draw("E SAME"); 
  h_measured_unf->Draw("HIST E SAME");
  h_measured_unf->SetMarkerColor(2);
  h_measured_unf->SetMarkerStyle(8);
  h_measured_unf->SetLineColor(2);
  h_measured_unf->SetLineWidth(1);

  leg->AddEntry(h_unfolded,"unfolded data","lep"); 
  leg->Draw(); 
  leg->AddEntry(h_true,"MC truth","lep"); 
  leg->Draw("SAME");
  leg->AddEntry(h_measured_unf,"data","lep"); 
  leg->Draw("SAME");

  TH1 * h_unfolded_r = (TH1*) h_unfolded->Clone();
  
  for(int k =1;k<9;k++){
    float unf=0;
    float tr = 0;
    float ratio =0;
    float err_unf=0;
    float err_tr = 0;
    float err_ratio =0;
    unf = h_unfolded->GetBinContent(k);
    tr = h_true->GetBinContent(k); 
    err_unf = h_unfolded->GetBinError(k);
    err_tr = h_true->GetBinError(k);
    ratio = unf/tr;
    err_ratio = sqrt((err_tr/tr)*(err_tr/tr)+(err_unf/unf)*(err_unf/unf));
    h_unfolded_r->SetBinContent(k,ratio); 
    h_unfolded_r->SetBinError(k,err_ratio);
  }

  pad2->cd();  
 
  h_unfolded_r->SetTitle(""); 
  h_unfolded_r->SetMarkerColor(1);
  h_unfolded_r->SetLineColor(1); 
  // h_unfolded_r->GetXaxis()->SetTitle(XaxisTitle.c_str()); 
  //  h_unfolded_r->GetYaxis()->SetTitle(YaxisTitle2.c_str());
  h_unfolded_r->GetYaxis()->SetTitle("");
  h_unfolded_r-> SetMaximum(2.5); 
  h_unfolded_r-> SetMinimum(-0.5);
  h_unfolded_r->Draw("E");
  TLine *line;

  if(var == "Mass")line = new TLine(100,1,800,1);
  else if(var == "Jets") line =  new TLine(0,1,4,1);
  line->SetLineColor(kRed);
  line->Draw("SAME");
  h_unfolded_r->Draw("E SAME");

  string png = var + "_ZZTo" + fs + ".png"; 
  string pdf = var + "_ZZTo" + fs + ".pdf";
  //c->Print(pdf.c_str());
  //c->Print(png.c_str());

  output->cd();   
  h_unfolded->Write(unfHistoName.c_str());
  h_measured_unf->Write(recoHistoName.c_str());
  h_true->Write(trueHistoName.c_str());
  output->Close();
}

void Unfold_data_All(std::string var = "Mass"){
  Unfold_data(var.c_str(),"4e");
  Unfold_data(var.c_str(),"4m");
  Unfold_data(var.c_str(),"2e2m");
}

void DoUnfoldedDataOverGenMCRatio(string var = "Mass", string fs = "4e"){
  
  system("mkdir ../../UnfoldingMacros/UnfoldFolder/");  
  filePath = "../../"; 
  dataFileName = filePath+"UnfoldingMacros/UnfoldFolder/UnfoldData_"+ var + ".root";
  outputFileName =filePath+"/UnfoldingMacros/UnfoldFolder/Ratio_UnfoldedDataOverGenMC.root"; 
 
  data = new TFile(dataFileName.c_str());
  output = new TFile(outputFileName.c_str(), "UPDATE");
  unfHistoName = "ZZTo"+ fs +"_" + var;
  recoHistoName = "ZZTo"+ fs +"_"+var+"_RECO_MC";
  trueHistoName = "ZZTo"+ fs +"_"+var+"_GEN";
  
  h_true = (TH1*) data->Get(trueHistoName.c_str());
  h_unfolded = (TH1*) data->Get(unfHistoName.c_str());
  
  TH1 * h_ratio  = (TH1*)h_unfolded->Clone("h_unfolded");
  float unf = 0;
  float gen = 0;
  float unf_err = 0;
  float gen_err = 0;
  float r = 0;
  float err = 0;
  int b = 0; 
  string finalstate;

  if(fs == "4m") finalstate = "4#mu";
  else if(fs == "2e2m") finalstate = "2e2#mu";
  else finalstate = fs;
 
  if(var =="Mass"){
    b=9;
    XaxisTitle = "m_{"+finalstate+"}";
  }
  else if(var =="Jets"){
    b=5;
    XaxisTitle = "Njets ("+finalstate+"-final state)";
  }

  for(int i = 1; i<b; i++){
    unf = 0;
    gen = 0;
    unf_err = 0;
    gen_err = 0;
    r = 0;
    err = 0;
    unf = h_unfolded->GetBinContent(i);
    gen = h_true->GetBinContent(i);  
    unf_err = h_unfolded->GetBinError(i);
    gen_err = h_true->GetBinError(i);
    r = unf/gen;
    err = sqrt((unf_err/unf)*(unf_err/unf)+(gen_err/gen)*(gen_err/gen))*r;
    h_ratio->SetBinContent(i,r); 
    h_ratio->SetBinError(i,err);
   
    cout << r <<" +- " << err << " " <<endl;

  }
   h_ratio->SetMaximum(2);
   h_ratio->GetXaxis()->SetTitle(XaxisTitle.c_str());
   string HistoRatio = "ZZTo"+ fs +"_"+ var +"_Ratio";
   
   output->cd();
   h_ratio->Write(HistoRatio.c_str());
   output->Close();
}

void DoAllRatios(string var = "Mass"){
  DoUnfoldedDataOverGenMCRatio(var.c_str(),"4e");
  DoUnfoldedDataOverGenMCRatio(var.c_str(),"4m");
  DoUnfoldedDataOverGenMCRatio(var.c_str(),"2e2m");

 }
void DoMCGenSystematic(string var = "Mass", string fs = "4e"){

 system("mkdir ../../UnfoldingMacros/UnfoldFolder/");  
  filePath = "../../";
  matrixFileName = filePath + var + "_test/matrices_Pow.root";
  dataFileName = filePath +var + "_test/DataToUnfold.root";
  outputFileName = filePath+"UnfoldingMacros/UnfoldFolder/UnfoldData_"+ var + "_MCgen.root";

  matrix = new TFile(matrixFileName.c_str());
  data = new TFile(dataFileName.c_str());
  output = new TFile(outputFileName.c_str(), "UPDATE");
  
  matrixName = "ResMat_qqggJJ_" + var + "_ZZTo" + fs + "_st_01";
  histoName = var +"_qqggJJ_ZZTo" + fs + "_st_01";
  histoNameGen = var + "Gen_qqggJJ_ZZTo" + fs + "_st_01";
  histoName_unf = "DataminusBkg_" + var + "_ZZTo"+fs;
 
  h_measured = (TH1*) matrix->Get(histoName.c_str());  
  h_true = (TH1*) matrix->Get(histoNameGen.c_str());
  h_Resmat = (TH2*)matrix->Get(matrixName.c_str()); 
  h_measured_unf = (TH1*) data->Get(histoName_unf.c_str());  

  RooUnfoldResponse response(h_measured, h_true, h_Resmat, "response", "response"); 
  RooUnfoldSvd   unfold_svd(&response, h_measured_unf, 4); 
  RooUnfoldBayes unfold_bayes(&response, h_measured_unf,2);  
    
  if(var == "Mass") h_unfolded= (TH1*) unfold_svd.Hreco(RooUnfold::kCovariance);
  else if(var == "Jets") h_unfolded= (TH1*) unfold_bayes.Hreco(RooUnfold::kCovariance);
 
  unfHistoName = "ZZTo"+ fs +"_" + var;
  
  output->cd();   
  h_unfolded->Write(unfHistoName.c_str());
  output->Close();
}


void DoQqggSystematic(string var = "Mass", string fs = "4e"){
 
  system("mkdir ../../UnfoldingMacros/UnfoldFolder/");  
  filePath = "../../";
  matrixFileName = filePath + var + "_test/matrices.root";
  dataFileName = filePath +var + "_test/DataToUnfold.root";
  outputFileName = filePath+"UnfoldingMacros/UnfoldFolder/UnfoldData_"+ var + "_qqgg.root";

  matrix = new TFile(matrixFileName.c_str());
  data = new TFile(dataFileName.c_str());
  output = new TFile(outputFileName.c_str(), "UPDATE");
  
  matrix_p = "ResMat_qqggJJ_" + var + "_ZZTo" + fs + "_qqp_ggm_01";
  histoReco_p = var + "_qqggJJ_ZZTo" + fs +"_qqp_ggm_01";
  histoGen_p = var + "Gen_qqggJJ_ZZTo" + fs +"_qqp_ggm_01";
  matrix_m = "ResMat_qqggJJ_" + var + "_ZZTo" + fs + "_qqm_ggp_01";
  histoReco_m = var + "_qqggJJ_ZZTo" + fs +"_qqm_ggp_01";
  histoGen_m = var + "Gen_qqggJJ_ZZTo" + fs +"_qqm_ggp_01";
  histoName_unf = "DataminusBkg_" + var + "_ZZTo"+fs;

  h_measured_p = (TH1*) matrix->Get(histoReco_p.c_str());  
  h_true_p = (TH1*) matrix->Get(histoGen_p.c_str());
  h_Resmat_p = (TH2*)matrix->Get(matrix_p.c_str()); 
  h_measured_m = (TH1*) matrix->Get(histoReco_m.c_str());  
  h_true_m = (TH1*) matrix->Get(histoGen_m.c_str());
  h_Resmat_m = (TH2*)matrix->Get(matrix_m.c_str()); 
  h_measured_unf = (TH1*) data->Get(histoName_unf.c_str());

  RooUnfoldResponse response_p(h_measured_p, h_true_p, h_Resmat_p, "response_p", "response_p"); 
  RooUnfoldResponse response_m(h_measured_m, h_true_m, h_Resmat_m, "response_m", "response_m"); 
  
  RooUnfoldSvd unfold_svd_p(&response_p, h_measured_unf, 4);
  RooUnfoldSvd unfold_svd_m(&response_m, h_measured_unf, 4);
  
  RooUnfoldBayes unfold_bayes_p(&response_p, h_measured_unf, 2);
  RooUnfoldBayes unfold_bayes_m(&response_m, h_measured_unf, 2);
    
  if(var == "Mass") {
    h_unfolded_p = (TH1*) unfold_svd_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_svd_m.Hreco(RooUnfold::kCovariance);
  }
  else if(var == "Jets"){ 
    h_unfolded_p = (TH1*) unfold_bayes_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes_m.Hreco(RooUnfold::kCovariance);
  }
 
  UnfHistoName_p = "ZZTo"+ fs + "_" + var + "_p";
  UnfHistoName_m = "ZZTo"+ fs + "_" + var + "_m";

  output->cd();
  h_unfolded_p->Write(UnfHistoName_p.c_str());
  h_unfolded_m->Write(UnfHistoName_m.c_str());   
  output->Close();

 }


void DoIrrBkgSystematic(string var = "Mass", string fs = "4e"){

#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  system("mkdir ../../UnfoldingMacros/UnfoldFolder/");  
  filePath = "../../";
  matrixFileName = filePath + var + "_test/matrices.root";
  dataFileName = filePath +var + "_test/DataToUnfold_syst.root";
  outputFileName = filePath+"UnfoldingMacros/UnfoldFolder/UnfoldData_"+ var + "_IrrBkg.root";
 
  matrix = new TFile(matrixFileName.c_str());
  data = new TFile(dataFileName.c_str());
  output = new TFile(outputFileName.c_str(), "UPDATE");
 
  matrixName = "ResMat_qqggJJ_" + var + "_ZZTo" + fs + "_st_01";
  histoName = var +"_qqggJJ_ZZTo" + fs + "_st_01";
  histoNameGen = var + "Gen_qqggJJ_ZZTo" + fs + "_st_01";
  data_p =  "DataminusBkg_irrp_"+var+"_ZZTo"+fs;
  data_m =  "DataminusBkg_irrm_"+var+"_ZZTo"+fs;
  
  h_measured = (TH1*) matrix->Get(histoName.c_str());  
  h_true = (TH1*) matrix->Get(histoNameGen.c_str());
  h_Resmat = (TH2*)matrix->Get(matrixName.c_str()); 
  h_data_p = (TH1*) data->Get(data_p.c_str());  
  h_data_m = (TH1*) data->Get(data_m.c_str()); 

  RooUnfoldResponse response(h_measured, h_true, h_Resmat, "response", "response"); 
  RooUnfoldSvd unfold_svd_p(&response, h_data_p, 4);
  RooUnfoldSvd unfold_svd_m(&response, h_data_m, 4);
  RooUnfoldBayes unfold_bayes_p(&response, h_data_p, 2);
  RooUnfoldBayes unfold_bayes_m(&response, h_data_m, 2);
  
  if(var == "Mass") {
    h_unfolded_p = (TH1*) unfold_svd_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_svd_m.Hreco(RooUnfold::kCovariance);
  }
  else if(var == "Jets"){ 
    h_unfolded_p = (TH1*) unfold_bayes_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes_m.Hreco(RooUnfold::kCovariance);
  }
  
  UnfHistoName_p = "ZZTo"+ fs + "_" + var + "_p";
  UnfHistoName_m = "ZZTo"+ fs + "_" + var + "_m";

  output->cd();
  h_unfolded_p->Write(UnfHistoName_p.c_str());
  h_unfolded_m->Write(UnfHistoName_m.c_str());   
  output->Close();
}


void DoRedBkgSystematic(string var = "Mass", string fs = "4e"){

#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
 
  system("mkdir ../../UnfoldingMacros/UnfoldFolder/");  
  filePath = "../../";
  matrixFileName = filePath + var + "_test/matrices.root";
  dataFileName = filePath +var + "_test/DataToUnfold_syst.root";
  outputFileName = filePath+"UnfoldingMacros/UnfoldFolder/UnfoldData_"+ var + "_RedBkg.root";
 
  matrix = new TFile(matrixFileName.c_str());
  data = new TFile(dataFileName.c_str());
  output = new TFile(outputFileName.c_str(), "UPDATE");
 
  matrixName = "ResMat_qqggJJ_" + var + "_ZZTo" + fs + "_st_01";
  histoName = var +"_qqggJJ_ZZTo" + fs + "_st_01";
  histoNameGen = var + "Gen_qqggJJ_ZZTo" + fs + "_st_01";
  data_p =  "DataminusBkg_redp_"+var+"_ZZTo"+fs;
  data_m =  "DataminusBkg_redm_"+var+"_ZZTo"+fs;
  
  h_measured = (TH1*) matrix->Get(histoName.c_str());  
  h_true = (TH1*) matrix->Get(histoNameGen.c_str());
  h_Resmat = (TH2*)matrix->Get(matrixName.c_str()); 
  h_data_p = (TH1*) data->Get(data_p.c_str());  
  h_data_m = (TH1*) data->Get(data_m.c_str()); 

  RooUnfoldResponse response(h_measured, h_true, h_Resmat, "response", "response"); 
  RooUnfoldSvd unfold_svd_p(&response, h_data_p, 4);
  RooUnfoldSvd unfold_svd_m(&response, h_data_m, 4);
  RooUnfoldBayes unfold_bayes_p(&response, h_data_p, 2);
  RooUnfoldBayes unfold_bayes_m(&response, h_data_m, 2);
  
  if(var == "Mass") {
    h_unfolded_p = (TH1*) unfold_svd_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_svd_m.Hreco(RooUnfold::kCovariance);
  }
  else if(var == "Jets"){ 
    h_unfolded_p = (TH1*) unfold_bayes_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes_m.Hreco(RooUnfold::kCovariance);
  }
  
  UnfHistoName_p = "ZZTo"+ fs + "_" + var + "_p";
  UnfHistoName_m = "ZZTo"+ fs + "_" + var + "_m";

  output->cd();
  h_unfolded_p->Write(UnfHistoName_p.c_str());
  h_unfolded_m->Write(UnfHistoName_m.c_str());   
  output->Close();


}

void DoUnfOverGenSystematic(string var = "Mass", string fs = "4e"){
  
  system("mkdir ../../UnfoldingMacros/UnfoldFolder/");  
  filePath = "../../";
  matrixFileName = filePath + var + "/weightedMatrices.root";
  dataFileName = filePath +var + "_test/DataToUnfold.root";
  outputFileName = filePath+"UnfoldingMacros/UnfoldFolder/UnfoldData_"+ var + "_UnfDataOverGenMC.root";

  matrix = new TFile(matrixFileName.c_str());
  data = new TFile(dataFileName.c_str());
  output = new TFile(outputFileName.c_str(), "UPDATE");
  
  matrixName = "ResMat_qqggJJ_" + var + "_ZZTo" + fs + "_st_01";
  histoName = var +"_qqggJJ_ZZTo" + fs + "_st_01";
  histoNameGen = var + "Gen_qqggJJ_ZZTo" + fs + "_st_01";
  histoName_unf = "DataminusBkg_" + var + "_ZZTo"+fs;
 
  h_measured = (TH1*) matrix->Get(histoName.c_str());  
  h_true = (TH1*) matrix->Get(histoNameGen.c_str());
  h_Resmat = (TH2*)matrix->Get(matrixName.c_str()); 
  h_measured_unf = (TH1*) data->Get(histoName_unf.c_str());  

  RooUnfoldResponse response(h_measured, h_true, h_Resmat, "response", "response"); 
  RooUnfoldSvd   unfold_svd(&response, h_measured_unf, 4); 
  RooUnfoldBayes unfold_bayes(&response, h_measured_unf,2);  
    
  if(var == "Mass") h_unfolded= (TH1*) unfold_svd.Hreco(RooUnfold::kCovariance);
  else if(var == "Jets") h_unfolded= (TH1*) unfold_bayes.Hreco(RooUnfold::kCovariance);
 
  unfHistoName = "ZZTo"+ fs +"_" + var;
  
  output->cd();   
  h_unfolded->Write(unfHistoName.c_str());
  output->Close();

}

void DoJERSystematic(string var = "Jets", string fs = "4e"){

#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  system("mkdir ../../UnfoldingMacros/UnfoldFolder/");  
  filePath = "../../";
  matrixFileName = filePath + var + "_test/matrices_JESJER.root";
  dataFileName = filePath +var + "_test/DataToUnfold.root";
  outputFileName = filePath+"UnfoldingMacros/UnfoldFolder/UnfoldData_"+ var + "_JER.root";

  matrix = new TFile(matrixFileName.c_str());
  data = new TFile(dataFileName.c_str());
  output = new TFile(outputFileName.c_str(), "UPDATE");

  matrix_p = "ResMat_qqggJJ_" + var + "_ZZTo" + fs + "_JERUp_01";
  histoReco_p = var + "_qqggJJ_ZZTo" + fs +"_JERUp_01";
  histoGen_p = var + "Gen_qqggJJ_ZZTo" + fs +"_JERUp_01";
  matrix_m = "ResMat_qqggJJ_" + var + "_ZZTo" + fs +"_JERDown_01";  
  histoReco_m = var + "_qqggJJ_ZZTo" + fs +"_JERDown_01";
  histoGen_m = var + "Gen_qqggJJ_ZZTo" + fs +"_JERDown_01"; 
  histoName_unf = "DataminusBkg_" + var + "_ZZTo"+fs;

  h_measured_p = (TH1*) matrix->Get(histoReco_p.c_str());  
  h_true_p = (TH1*) matrix->Get(histoGen_p.c_str());
  h_Resmat_p = (TH2*)matrix->Get(matrix_p.c_str()); 
  h_measured_m = (TH1*) matrix->Get(histoReco_m.c_str());  
  h_true_m = (TH1*) matrix->Get(histoGen_m.c_str());
  h_Resmat_m = (TH2*)matrix->Get(matrix_m.c_str()); 
  h_measured_unf = (TH1*) data->Get(histoName_unf.c_str());

  RooUnfoldResponse response_p(h_measured_p, h_true_p, h_Resmat_p, "response_p", "response_p"); 
  RooUnfoldResponse response_m(h_measured_m, h_true_m, h_Resmat_m, "response_m", "response_m"); 
  
  RooUnfoldSvd unfold_svd_p(&response_p, h_measured_unf, 4);
  RooUnfoldSvd unfold_svd_m(&response_m, h_measured_unf, 4);
  
  RooUnfoldBayes unfold_bayes_p(&response_p, h_measured_unf, 2);
  RooUnfoldBayes unfold_bayes_m(&response_m, h_measured_unf, 2);
    
  if(var == "Mass") {
    h_unfolded_p = (TH1*) unfold_svd_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_svd_m.Hreco(RooUnfold::kCovariance);
  }
  else if(var == "Jets"){ 
    h_unfolded_p = (TH1*) unfold_bayes_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes_m.Hreco(RooUnfold::kCovariance);
  }
 
  UnfHistoName_p = "ZZTo"+ fs + "_" + var + "_p";
  UnfHistoName_m = "ZZTo"+ fs + "_" + var + "_m";

  output->cd();
  h_unfolded_p->Write(UnfHistoName_p.c_str());
  h_unfolded_m->Write(UnfHistoName_m.c_str());   
  output->Close();
  
}


void DoJESSystematic_ModMat(string var = "Jets", string fs = "4e"){

#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  system("mkdir ../../UnfoldingMacros/UnfoldFolder/");  
  filePath = "../../";
  matrixFileName = filePath + var + "_test/matrices_JESJER.root";
  dataFileName = filePath +var + "_test/DataToUnfold.root";
  outputFileName = filePath+"UnfoldingMacros/UnfoldFolder/UnfoldData_"+ var + "_JES_ModMat.root";

  matrix = new TFile(matrixFileName.c_str());
  data = new TFile(dataFileName.c_str());
  output = new TFile(outputFileName.c_str(), "UPDATE");

  matrix_p = "ResMat_qqggJJ_" + var + "_ZZTo" + fs + "_JESUp_01";
  histoReco_p = var + "_qqggJJ_ZZTo" + fs +"_JESUp_01";
  histoGen_p = var + "Gen_qqggJJ_ZZTo" + fs +"_JESUp_01";
  matrix_m = "ResMat_qqggJJ_" + var + "_ZZTo" + fs +"_JESDown_01";  
  histoReco_m = var + "_qqggJJ_ZZTo" + fs +"_JESDown_01";
  histoGen_m = var + "Gen_qqggJJ_ZZTo" + fs +"_JESDown_01"; 
  histoName_unf = "DataminusBkg_" + var + "_ZZTo"+fs;

  h_measured_p = (TH1*) matrix->Get(histoReco_p.c_str());  
  h_true_p = (TH1*) matrix->Get(histoGen_p.c_str());
  h_Resmat_p = (TH2*)matrix->Get(matrix_p.c_str()); 
  h_measured_m = (TH1*) matrix->Get(histoReco_m.c_str());  
  h_true_m = (TH1*) matrix->Get(histoGen_m.c_str());
  h_Resmat_m = (TH2*)matrix->Get(matrix_m.c_str()); 
  h_measured_unf = (TH1*) data->Get(histoName_unf.c_str());

  RooUnfoldResponse response_p(h_measured_p, h_true_p, h_Resmat_p, "response_p", "response_p"); 
  RooUnfoldResponse response_m(h_measured_m, h_true_m, h_Resmat_m, "response_m", "response_m"); 
  
  RooUnfoldSvd unfold_svd_p(&response_p, h_measured_unf, 4);
  RooUnfoldSvd unfold_svd_m(&response_m, h_measured_unf, 4);
  
  RooUnfoldBayes unfold_bayes_p(&response_p, h_measured_unf, 2);
  RooUnfoldBayes unfold_bayes_m(&response_m, h_measured_unf, 2);
    
  if(var == "Mass") {
    h_unfolded_p = (TH1*) unfold_svd_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_svd_m.Hreco(RooUnfold::kCovariance);
  }
  else if(var == "Jets"){ 
    h_unfolded_p = (TH1*) unfold_bayes_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes_m.Hreco(RooUnfold::kCovariance);
  }
 
  UnfHistoName_p = "ZZTo"+ fs + "_" + var + "_p";
  UnfHistoName_m = "ZZTo"+ fs + "_" + var + "_m";

  output->cd();
  h_unfolded_p->Write(UnfHistoName_p.c_str());
  h_unfolded_m->Write(UnfHistoName_m.c_str());   
  output->Close();

}

void DoJESSystematic_ModData(string var = "Jets", string fs = "4e"){

#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  system("mkdir ../../UnfoldingMacros/UnfoldFolder/");  
  filePath = "../../";
  matrixFileName = filePath + var + "_test/matrices.root";
  dataFileName = filePath +var + "_test/DataToUnfold_JES.root";
  outputFileName = filePath+"UnfoldingMacros/UnfoldFolder/UnfoldData_"+ var + "_JES_ModData.root";

  matrix = new TFile(matrixFileName.c_str());
  data = new TFile(dataFileName.c_str());
  output = new TFile(outputFileName.c_str(), "UPDATE");

  matrixName = "ResMat_qqggJJ_" + var + "_ZZTo" + fs + "_st_01";
  histoName = var +"_qqggJJ_ZZTo" + fs + "_st_01";
  histoNameGen = var + "Gen_qqggJJ_ZZTo" + fs + "_st_01";
  data_p =  "DataminusBkg_"+var+"_ZZTo"+fs+ "_JESUp";
  data_m =  "DataminusBkg_"+var+"_ZZTo"+fs+ "_JESDown";
  
  h_measured = (TH1*) matrix->Get(histoName.c_str());  
  h_true = (TH1*) matrix->Get(histoNameGen.c_str());
  h_Resmat = (TH2*)matrix->Get(matrixName.c_str()); 
  h_data_p = (TH1*) data->Get(data_p.c_str());  
  h_data_m = (TH1*) data->Get(data_m.c_str()); 

  RooUnfoldResponse response(h_measured, h_true, h_Resmat, "response", "response"); 
  RooUnfoldSvd unfold_svd_p(&response, h_data_p, 4);
  RooUnfoldSvd unfold_svd_m(&response, h_data_m, 4);
  RooUnfoldBayes unfold_bayes_p(&response, h_data_p, 2);
  RooUnfoldBayes unfold_bayes_m(&response, h_data_m, 2);
 
  if(var == "Mass") {
    h_unfolded_p = (TH1*) unfold_svd_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_svd_m.Hreco(RooUnfold::kCovariance);
  }
  else if(var == "Jets"){ 
    h_unfolded_p = (TH1*) unfold_bayes_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes_m.Hreco(RooUnfold::kCovariance);
  }
 
  UnfHistoName_p = "ZZTo"+ fs + "_" + var + "_p";
  UnfHistoName_m = "ZZTo"+ fs + "_" + var + "_m";

  output->cd();
  h_unfolded_p->Write(UnfHistoName_p.c_str());
  h_unfolded_m->Write(UnfHistoName_m.c_str());   
  output->Close();

}


void  DoAllSystematics(string var = "Mass"){
  Unfold_data(var.c_str(), "4e");
  Unfold_data(var.c_str(), "4m");
  Unfold_data(var.c_str(), "2e2m");

  DoMCGenSystematic(var.c_str(), "4e");
  DoMCGenSystematic(var.c_str(), "4m");
  DoMCGenSystematic(var.c_str(), "2e2m");
  
  DoQqggSystematic(var.c_str(), "4e"); 
  DoQqggSystematic(var.c_str(), "4m");
  DoQqggSystematic(var.c_str(), "2e2m");
  
  DoIrrBkgSystematic(var.c_str(), "4e"); 
  DoIrrBkgSystematic(var.c_str(), "4m"); 
  DoIrrBkgSystematic(var.c_str(), "2e2m");
  
  DoRedBkgSystematic(var.c_str(), "4e"); 
  DoRedBkgSystematic(var.c_str(), "4m"); 
  DoRedBkgSystematic(var.c_str(), "2e2m"); 
  
  DoUnfOverGenSystematic(var.c_str(), "4e");
  DoUnfOverGenSystematic(var.c_str(), "4m");
  DoUnfOverGenSystematic(var.c_str(), "2e2m");
  
  if(var == "Jets"){
    DoJERSystematic(var.c_str(), "4e"); 
    DoJERSystematic(var.c_str(), "4m");
    DoJERSystematic(var.c_str(), "2e2m"); 
    
    DoJESSystematic_ModMat(var.c_str(), "4e"); 
    DoJESSystematic_ModMat(var.c_str(), "4m");
    DoJESSystematic_ModMat(var.c_str(), "2e2m");  
    
    DoJESSystematic_ModData(var.c_str(), "4e"); 
    DoJESSystematic_ModData(var.c_str(), "4m");
    DoJESSystematic_ModData(var.c_str(), "2e2m"); 
  }
}

void PlotResults(string var = "Jets", string fs = "4e", string syst = "MCgen"){
  
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  TFile * unfSyst;
  string systFileName;
  string systHistoName_p;
  string systHistoName_m;
  TH1 *h_true_pow;
  
  system("mkdir ../../UnfoldingMacros/UnfoldFolder/");  
  filePath = "../../";
  matrixFileName = filePath + var + "_test/matrices_Pow.root";
  dataFileName = filePath+"UnfoldingMacros/UnfoldFolder/UnfoldData_"+ var + ".root";
  systFileName =  filePath+"UnfoldingMacros/UnfoldFolder/UnfoldData_"+ var + "_"+syst + ".root";

  data = new TFile(dataFileName.c_str());
  unfSyst = new TFile(systFileName.c_str());
  matrix = new TFile(matrixFileName.c_str());

  unfHistoName = "ZZTo"+ fs +"_" + var;
  recoHistoName = "ZZTo"+ fs +"_"+var+"_RECO_MC";
  trueHistoName = "ZZTo"+ fs +"_"+var+"_GEN";
  systHistoName_p = "ZZTo"+ fs +"_" + var + "_p";
  systHistoName_m =  "ZZTo"+ fs +"_" + var + "_m";
  histoNameGen = var + "Gen_qqggJJ_ZZTo" + fs + "_st_01";
  
  h_true_pow = (TH1*) matrix->Get(histoNameGen.c_str());
  h_true = (TH1*) data->Get(trueHistoName.c_str());
  h_measured_unf = (TH1*) data->Get(recoHistoName.c_str());  
  h_unfolded = (TH1*) data->Get(unfHistoName.c_str());
 
  float max;
  float min;
  TH1* h_max; 
  TH1* h_min; 
  float syst_percentage[9];

  std::cout << "===================================================================================" << std::endl;
  std::cout << "                                 " << syst.c_str() << std::endl;
  std::cout << "===================================================================================" << std::endl;
  if(syst == "MCgen" || syst == "UnfDataOverGenMC"){
    h_unfolded_p = (TH1*) unfSyst->Get(unfHistoName.c_str());
    h_max = (TH1*)h_unfolded_p->Clone("");;
    h_min = (TH1*)h_unfolded_p->Clone("");;
    for(int i = 1; i<9; i++){
      max = 0;
      min = 0;
      syst_percentage[i-1]=0;
      max = h_unfolded->GetBinContent(i) + fabs(h_unfolded_p->GetBinContent(i)- h_unfolded->GetBinContent(i));
      min =  h_unfolded->GetBinContent(i) - fabs(h_unfolded_p->GetBinContent(i)- h_unfolded->GetBinContent(i));
    
      h_max->SetBinContent(i,max); 
      h_min->SetBinContent(i,min);
      syst_percentage[i-1]=0.5*(max-min)*100/(h_unfolded->GetBinContent(i));

      if(max!=0.) std::cout  << i << "-bin: min = "<< min << " max = " << max << " " << syst_percentage[i-1]  << "% " << syst.c_str() << " systematic uncertainty" << std::endl;
 
    }
  }
  else if(syst=="qqgg"||syst=="IrrBkg"||syst=="RedBkg"||syst=="JER"||syst=="JES_ModMat"||syst=="JES_ModData"){
    h_unfolded_p = (TH1*) unfSyst->Get(systHistoName_p.c_str());
    h_unfolded_m = (TH1*) unfSyst->Get(systHistoName_m.c_str());
    h_max = (TH1*)h_unfolded_p->Clone("");;
    h_min = (TH1*)h_unfolded_p->Clone("");;
    for(int i = 1; i<9; i++){
      max = 0;
      min = 0;
      syst_percentage[i-1]=0;
      max = TMath::Max(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
      min = TMath::Min(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
      // std::cout << min << " " << max << std::endl;
      h_max->SetBinContent(i,max); 
      h_min->SetBinContent(i,min); 
      syst_percentage[i-1]=0.5*(max-min)*100/(h_unfolded->GetBinContent(i)); 
     
      if(max!=0.) std::cout  << i << "-bin: min = "<< min << " max = " << max << " " << syst_percentage[i-1]  << "% " << syst.c_str() << " systematic uncertainty" << std::endl;
    }
  } 
  else std::cout << " Wrong systematic uncertainty!!!!" << std::endl;

  float M = syst_percentage[0];
  float m = syst_percentage[0];
  for(int i = 1; i<8; i++) {
      if(syst_percentage[i] > M) M = syst_percentage[i];
      if(syst_percentage[i] < m) m = syst_percentage[i];
     }

  std::cout << "Max syst = " << M << "%   Min syst = " << m << "%" << endl;

  TCanvas *c = new TCanvas ("c","c");
  TLegend *leg = new TLegend(0.7,0.65,0.6,0.85); 
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);

  string XaxisTitle;
  float max_y =0;
  if(var == "Mass"){
    if(fs == "4e")  XaxisTitle = "m_{4e}";
    else if(fs == "4m")  XaxisTitle = "m_{4#mu}";
    else if(fs == "2e2m") XaxisTitle = "m_{2e2#mu}"; 
    max_y = h_max->GetBinContent(2)+20;//*2+100;
  }
  else if(var == "Jets"){
    if(fs == "4e")  XaxisTitle = "N jets (4e final state)";
    else if(fs == "4m")  XaxisTitle = "N jets (4#mu final state)";
    else if(fs == "2e2m") XaxisTitle = "N jets (2e2#mu final state)"; 
    max_y = h_max->GetBinContent(1)+50;//*2+100;
  }  
  // h_max->GetXaxis()->SetRange(70.,800.);
  //std::cout << h_max <<std::endl;
  h_max->Draw("HIST");
  //h_max->GetXaxis()->SetRange(70.,800.); 
  h_max->SetFillColor(2);
  h_max->SetFillStyle(3144);
  h_max->SetLineColor(10);
  h_max->SetTitle("");
  h_max->GetXaxis()->SetRange(0,25);
  h_max->GetXaxis()->SetTitle(XaxisTitle.c_str()); 
  h_max->SetMaximum(max_y);
  //h_max->SetMinimum(0); 
  h_min->Draw("HIST SAME");
  h_min->SetFillColor(10);
  h_min->SetLineColor(10);
  
  h_unfolded->Draw("E SAME"); 
  h_unfolded->SetLineColor(1); 
  //h_unfolded->SetLineStyle(2);
  h_unfolded->SetLineWidth(1);
  h_unfolded->SetMarkerStyle(7);

  h_true->SetLineColor(8);
  h_true->SetLineStyle(2);  
  h_true->SetLineWidth(2);
  h_true->Draw("E HIST SAME");
  h_true_pow->SetLineColor(9);
  h_true_pow->SetLineStyle(2);  
  h_true_pow->SetLineWidth(2);
  h_true_pow->Draw("E HIST SAME");
  h_unfolded->Draw("E SAME"); 
  gPad->RedrawAxis();
  

  string syst_uncertainty;
  if(syst == "MCgen") syst_uncertainty = "MC generator syst. unc.";
  else if(syst == "UnfDataOverGenMC") syst_uncertainty = "Data/MC syst. unc.";
  else if(syst=="qqgg") syst_uncertainty = "#sigma_{qq}/#sigma_{gg} syst. unc.";
  else if(syst=="IrrBkg") syst_uncertainty = " Irr. Bkg. syst. unc.";
  else if(syst=="RedBkg") syst_uncertainty = " Red. Bkg. syst. unc.";
  else if(syst=="JER") syst_uncertainty = "JER syst. unc.";
  else if(syst=="JES_ModMat") syst_uncertainty = "JES syst. unc. (modified matrix)";
  else if(syst=="JES_ModData") syst_uncertainty = "JES syst. unc. (modified data)";
  else std::cout << "wrong systematic uncerainty!!!!" << std::endl;

  leg->AddEntry(h_unfolded,"unfolded data","lp"); 
  leg->Draw(); 
  leg->AddEntry(h_true,"MC truth (Madgraph Set)","l"); 
  leg->Draw("SAME"); 
  leg->AddEntry(h_true_pow,"MC truth (Powheg Set)","l"); 
  leg->Draw("SAME");
  leg->AddEntry(h_max,syst_uncertainty.c_str(),"f"); 
  leg->Draw("SAME");

  string png = "ZZTo"+fs+"_"+var+"_"+syst+".png";
  string pdf ="ZZTo"+fs+"_"+var+"_"+syst+".pdf";
  // c->Print(png.c_str());
  // c->Print(pdf.c_str());


}

void AllPlots_fs(string var = "Jets", string fs = "4e"){
  std::cout << "===================================================================================" << std::endl;
  std::cout << "===================================================================================" << std::endl;
  std::cout << "                                 " << fs.c_str() << std::endl;
  std::cout << "===================================================================================" << std::endl;
  PlotResults(var.c_str(),fs.c_str(),"MCgen");
  PlotResults(var.c_str(),fs.c_str(),"qqgg");
  PlotResults(var.c_str(),fs.c_str(),"IrrBkg");
  PlotResults(var.c_str(),fs.c_str(),"RedBkg");
  PlotResults(var.c_str(),fs.c_str(),"UnfDataOverGenMC");

  if(var =="Jets"){
    PlotResults(var.c_str(),fs.c_str(),"JER"); 
    PlotResults(var.c_str(),fs.c_str(),"JES_ModMat"); 
    PlotResults(var.c_str(),fs.c_str(),"JES_ModData");
  }
}

void AllPlots(string var = "Jets"){
  AllPlots_fs(var.c_str(),"4e");
  AllPlots_fs(var.c_str(),"4m");
  AllPlots_fs(var.c_str(),"2e2m");
}

void AllYouNeed(){
  Unfold_data_All("Mass");
  DoAllRatios("Mass");
  DoAllSystematics("Mass");
  Unfold_data_All("Jets");
  DoAllRatios("Jets");
  DoAllSystematics("Jets");
}
#ifndef __CINT__
int main () { DoAllSystematics(); return 0; }  // Main program when run stand-alone
#endif
