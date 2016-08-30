#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <string>
#include <sstream>
using std::cout;
using std::endl;

//#include "TRandom.h"
#include <TROOT.h>
#include <TLatex.h>
#include <TFile.h>
#include <TStyle.h>
#include "TH1D.h"
#include "TMatrix.h"
#include "TStyle.h"
#include "TLine.h"
#include "TFrame.h"
#include "TH2.h"
#include "TLegend.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
#include "TCanvas.h"
//#include <boost/filesystem.hpp>
#include "tdrstyle.C"
#include "CMS_lumi.C"
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
std::string recoMCHistoName;
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
std::string MCgen;

/*To unfold distributions, first download the RooUnfold package following these instructions:
 http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html
Then, in ../RooUnfold-1.1.1/RooUnfold (after "make"):
root
.L ../UnfoldingMacros/Unfold.cpp+ 
AllYouNeed_Var2("Variable","Folder_for_Plots")
*/

void Unfold_data(string var = "Mass", string fs = "4e", bool mad =1,bool tightregion = 0, string date = "test"){
  

#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  
  string tightfr;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";
 
  filePath = "../../../";
  if(mad ==1){
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Mad/").c_str());  
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Mad.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + ".root";
    MCgen = "_Mad";
  }
  else{
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Pow/").c_str());  
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Pow.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + ".root";
    MCgen = "_Pow";
  }

  system(("mkdir ~/www/VBS/"+date).c_str()); 
  system(("cp ~/www/VBS/index.php ~/www/VBS/"+date).c_str());
  system(("mkdir ~/www/VBS/"+date+"/"+ var).c_str()); 
  system(("cp ~/www/VBS/index.php ~/www/VBS/"+date+"/"+ var).c_str());

  dataFileName = "../" +var + "_test/DataToUnfold.root";
  //dataFileName = filePath +var + "_test/DataToUnfoldFake.root";
  
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
  RooUnfoldBayes unfold_bayes(&response, h_measured_unf,4);  
  RooUnfoldBayes unfold_bayes4(&response, h_measured_unf,4); 
    
  if(var == "Mass") h_unfolded= (TH1*) unfold_svd.Hreco(RooUnfold::kCovariance);
  else if  (var == "dRZZ")  h_unfolded= (TH1*) unfold_bayes4.Hreco(RooUnfold::kCovariance);
  else h_unfolded= (TH1*) unfold_bayes.Hreco(RooUnfold::kCovariance);
 
  unfHistoName = "ZZTo"+ fs +"_" + var;
  recoHistoName = "ZZTo"+ fs +"_"+var+"_RECO"; 
  recoMCHistoName = "ZZTo"+ fs +"_"+var+"_RECO_MC";
  trueHistoName = "ZZTo"+ fs +"_"+var+"_GEN";
  string finalstate;
  if(fs == "4m") finalstate = "4#mu";
  else if(fs == "2e2m") finalstate = "2e2#mu";
  else finalstate = fs;

  if(var == "Mass") XaxisTitle = "m_{" + finalstate + "} [GeV]"; 
  else if(var == "Jets")  XaxisTitle = "N Jets (" + finalstate + " final state)";
  else if(var == "Jets_Central")  XaxisTitle = "N Central Jets (" + finalstate + " final state)";
  else if(var == "Mjj" ||var == "Mjj_Central")  XaxisTitle = "m_{jj} (" + finalstate + " final state) [GeV]";
  else if(var == "Deta"||var == "Deta_Central")  XaxisTitle = "#Delta#eta_{jj} (" + finalstate + " final state)";
  else if(var =="PtJet1")  XaxisTitle = "p_{T}^{jet1} ("+ finalstate+" final state) [GeV]";
  else if(var =="PtJet2")  XaxisTitle = "p_{T}^{jet2} ("+ finalstate+" final state) [GeV]";
  else if(var =="EtaJet1")  XaxisTitle = "|#eta^{jet1}| ("+ finalstate+" final state)";
  else if(var =="EtaJet2")  XaxisTitle = "|#eta^{jet2}| ("+ finalstate+" final state)";
  else if(var =="dRZZ")  XaxisTitle = "#DeltaR(Z_{1},Z_{2}) ("+ finalstate+" final state)";
  YaxisTitle = "Events";
  string YaxisTitle2 = "Unfolded/True";
 
  float max = 0;
  if(var == "Mass") max =  h_unfolded->GetBinContent(2)+h_unfolded->GetBinContent(2)*0.15+h_unfolded->GetBinError(2); 
  else if(var == "Jets"||var == "Jets_Central")  max =  h_unfolded->GetBinContent(1)+h_unfolded->GetBinContent(1)*0.15+h_unfolded->GetBinError(1); 
  else if(var == "Mjj"|| var == "Mjj_Central") max =  h_true->GetBinContent(1)+h_true->GetBinContent(1)*0.20+h_unfolded->GetBinError(1); 
  else if(var == "Deta"|| var == "Deta_Central" ) max =  h_true->GetBinContent(1)+h_true->GetBinContent(1)*0.15+h_unfolded->GetBinError(1); 
  else if(var == "PtJet1") max =  h_true->GetBinContent(1)+h_true->GetBinContent(1)*0.15+h_unfolded->GetBinError(1); 
  else if(var == "PtJet2") max =  h_true->GetBinContent(1)+h_true->GetBinContent(1)*0.15+h_unfolded->GetBinError(1); 
  else if(var == "EtaJet1") max =  h_true->GetBinContent(1)+h_true->GetBinContent(1)*0.15+h_unfolded->GetBinError(1);
  else if(var == "EtaJet2") max =  h_true->GetBinContent(1)+h_true->GetBinContent(1)*0.15+h_unfolded->GetBinError(1); 
  else if(var == "dRZZ") max =  h_true->GetBinContent(4)+h_true->GetBinContent(4)+h_unfolded->GetBinError(1);
  else  max =  h_true->GetBinContent(1)+h_true->GetBinContent(1)*0.15+h_unfolded->GetBinError(1);
 
  TCanvas *c = new TCanvas ("c","c");
  TLegend *leg = new TLegend(0.90,0.65,0.70,0.85); 
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
  pad2->SetTopMargin(0.10);
  pad2->SetRightMargin(0.10);
  pad2->SetLeftMargin(0.10); 
  pad2->SetBottomMargin(0.35);
  pad2->Draw(); 
    
  pad1->cd();
  h_true->SetTitle("");
  h_true->GetXaxis()->SetRange(0,25);
  h_true->GetXaxis()->SetTitle(XaxisTitle.c_str()); 
  h_true->GetYaxis()->SetTitle(YaxisTitle.c_str());
  h_true->GetXaxis()->SetLabelOffset(0.5); 
 
  h_true->SetLineColor(kBlue); 
  h_true->SetLineWidth(1);
  h_true->SetLineStyle(7);
  h_true->SetMaximum(max);
  h_true->SetMinimum(0); 
  h_true->Draw("HIST E");
  h_measured->SetLineColor(kRed); 
  h_measured->SetLineStyle(7);
  h_measured->SetLineWidth(1);
  h_measured->Draw("HIST E SAME");
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
  
 
  leg->AddEntry(h_true,"MC truth","l"); 
  leg->Draw("SAME");
  leg->AddEntry(h_measured,"MC reco","l"); 
  leg->Draw("SAME");
  leg->AddEntry(h_unfolded,"unfolded data","lep"); 
  leg->Draw(); 
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
  h_unfolded_r->GetYaxis()->SetLabelSize(0.10);  
  h_unfolded_r->GetYaxis()->SetNdivisions(306); 
  h_unfolded_r->GetXaxis()->SetLabelSize(0.13); 
  h_unfolded_r->GetXaxis()->SetLabelOffset(0.05);
  h_unfolded_r->GetXaxis()->SetTitle(XaxisTitle.c_str()); 
  h_unfolded_r->GetYaxis()->SetTitle(YaxisTitle2.c_str());
  h_unfolded_r->GetXaxis()->SetTitleOffset(1.2); 
  h_unfolded_r->GetYaxis()->SetTitleOffset(0.4);
  h_unfolded_r->GetYaxis()->SetTitleSize(0.10);
  h_unfolded_r->GetXaxis()->SetTitleSize(0.13);
    h_unfolded_r-> SetMaximum(2.5); 
  h_unfolded_r-> SetMinimum(-0.5);
  h_unfolded_r->Draw("E");
  TLine *line;

  if(var == "Mass")line = new TLine(100,1,800,1);
  else if(var == "Jets" ||var == "Jets_Central") line =  new TLine(0,1,4,1);
  else if(var == "Mjj" || var == "Mjj_Central") line =  new TLine(0,1,800,1);
  else if(var == "Deta"  || var == "Deta_Central") line =  new TLine(0,1,4.7,1);
  else if(var =="PtJet1"||var =="PtJet2" )  line =  new TLine(30,1,500,1); 
  else if(var =="EtaJet1"||var =="EtaJet2" )  line =  new TLine(0,1,4.7,1); 
  else if(var == "dRZZ") line = new TLine(0,1,6,1);
  line->SetLineColor(kRed);
  line->Draw("SAME");
  h_unfolded_r->Draw("E SAME");
  
  // string png = "~/www/VBS/"+date+"/"+ var+"/"+var + "_ZZTo" + fs + MCgen+tightfr +".png"; 
  // string pdf = "~/www/VBS/"+date+"/"+ var+"/"+var + "_ZZTo" + fs + MCgen+tightfr +".pdf";
  // c->SaveAs(pdf.c_str());
  // c->SaveAs(png.c_str());

  output->cd();   
  h_unfolded->Write(unfHistoName.c_str());
  h_measured_unf->Write(recoHistoName.c_str());
  h_measured->Write(recoMCHistoName.c_str());
  h_true->Write(trueHistoName.c_str());
  output->Close();
}

void Unfold_data_All(std::string var = "Mass", bool mad = 1,bool tightregion = 0, string date = "test"){
  Unfold_data(var.c_str(),"4e",mad, tightregion,date.c_str());
  Unfold_data(var.c_str(),"4m",mad, tightregion,date.c_str());
  Unfold_data(var.c_str(),"2e2m",mad, tightregion,date.c_str());
}

void DoUnfoldedDataOverGenMCRatio(string var = "Mass", string fs = "4e", bool mad =1,bool tightregion = 0){
  
  string tightfr;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";

  filePath = "../../../";  

  if(mad ==1){
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Mad/").c_str());  
    dataFileName = filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + ".root";
    outputFileName =filePath+"UnfoldFolder"+tightfr+"_Mad/Ratio_UnfoldedDataOverGenMC.root"; 
    MCgen = "_Mad";
  }
  else{
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Pow/").c_str());  
    dataFileName = filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + ".root";
    outputFileName =filePath+"UnfoldFolder"+tightfr+"_Pow/Ratio_UnfoldedDataOverGenMC.root"; 
    MCgen = "_Pow";
  }

  data = new TFile(dataFileName.c_str());
  output = new TFile(outputFileName.c_str(), "UPDATE");
  unfHistoName = "ZZTo"+ fs +"_" + var;
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
 else if(var =="Jets_Central"){
    b=5;
    XaxisTitle = "N central jets ("+finalstate+"-final state)";
 }
 else if(var =="Mjj" || var =="Mjj_Central"){
    b=3;
    XaxisTitle = "m_{jj} ("+finalstate+"-final state)";
  }
 else if(var =="Deta"  || var =="Deta_Central"){
    b=3;
    XaxisTitle = "#Delta#eta_{jj} ("+finalstate+"-final state)";
  }
 else if(var =="PtJet1")  {
   b=6;
   XaxisTitle = "p_{T}^{jet1} ("+ finalstate+" final state)";
 }
 else if(var =="PtJet2") {
   b=6;
   XaxisTitle = "p_{T}^{jet2} ("+ finalstate+" final state)";
 } 
 else if(var =="EtaJet1")  {
   b=6;
   XaxisTitle = "|#eta^{jet1}| ("+ finalstate+" final state)";
 } 
 else if(var =="EtaJet2")  {
   b=6;
   XaxisTitle = "|#eta^{jet2}| ("+ finalstate+" final state)";
 } 
 else if(var =="dRZZ")  {
   b=7;
   XaxisTitle = "#DeltaR(Z_{1},Z_{2}) ("+ finalstate+" final state)";
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
   string HistoRatio = "ZZTo"+ fs +"_"+ var+"_Ratio";
   
   output->cd();
   h_ratio->Write(HistoRatio.c_str());
   output->Close();
}

void DoAllRatios(string var = "Mass", bool mad =1, bool tightregion =0){
  DoUnfoldedDataOverGenMCRatio(var.c_str(),"4e",mad,tightregion);
  DoUnfoldedDataOverGenMCRatio(var.c_str(),"4m",mad,tightregion);
  DoUnfoldedDataOverGenMCRatio(var.c_str(),"2e2m",mad,tightregion);

 }
void DoMCGenSystematic(string var = "Mass", string fs = "4e",bool mad =1 ,bool tightregion = 0){
  
  string tightfr;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";

  filePath = "../../../";
  
  if(mad ==1){
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Mad/").c_str());  
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Pow.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + "_MCgen.root";
    MCgen = "_Mad";
  }
  else{
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Pow/").c_str());  
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Mad.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + "_MCgen.root";
    MCgen = "_Pow";
  }


  dataFileName  = "../" +var + "_test/DataToUnfold.root";
  matrix        = new TFile(matrixFileName.c_str());
  data          = new TFile(dataFileName.c_str());
  output        = new TFile(outputFileName.c_str(), "UPDATE");
  
  matrixName    = "ResMat_qqggJJ_" + var + "_ZZTo" + fs + "_st_01";
  histoName     =  var +"_qqggJJ_ZZTo" + fs + "_st_01";
  histoNameGen  =  var + "Gen_qqggJJ_ZZTo" + fs + "_st_01";
  histoName_unf = "DataminusBkg_" + var + "_ZZTo"+ fs;
 
  h_measured = (TH1*) matrix->Get(histoName.c_str());
  h_true = (TH1*) matrix->Get(histoNameGen.c_str());
  h_Resmat = (TH2*)matrix->Get(matrixName.c_str()); 
  h_measured_unf = (TH1*) data->Get(histoName_unf.c_str()); 

  RooUnfoldResponse response(h_measured, h_true, h_Resmat, "response", "response"); 
  RooUnfoldSvd   unfold_svd(&response, h_measured_unf, 4); 
  RooUnfoldBayes unfold_bayes(&response, h_measured_unf,4);  
  RooUnfoldBayes unfold_bayes4(&response, h_measured_unf,4); 

  if(var == "Mass") h_unfolded= (TH1*) unfold_svd.Hreco(RooUnfold::kCovariance);
  else if(var == "dRZZ") h_unfolded= (TH1*) unfold_bayes4.Hreco(RooUnfold::kCovariance);
  else h_unfolded= (TH1*) unfold_bayes.Hreco(RooUnfold::kCovariance);
   
  unfHistoName = "ZZTo"+ fs +"_" + var;
  output->cd();   
  h_unfolded->Write(unfHistoName.c_str());
  output->Close();
}


void DoQqggSystematic(string var = "Mass", string fs = "4e", bool mad =1 ,bool tightregion = 0){
  
  string tightfr;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";
 
  filePath = "../../../";
  if(mad ==1){
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Mad/").c_str());  
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Mad.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + "_qqgg.root";
  }
  else{
      system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Pow/").c_str());       
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Pow.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + "_qqgg.root";
  }
  dataFileName = "../" +var + "_test/DataToUnfold.root";
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
  
  RooUnfoldBayes unfold_bayes_p(&response_p, h_measured_unf, 4);
  RooUnfoldBayes unfold_bayes_m(&response_m, h_measured_unf, 4);
    
  RooUnfoldBayes unfold_bayes4_p(&response_p, h_measured_unf, 4);
  RooUnfoldBayes unfold_bayes4_m(&response_m, h_measured_unf, 4);

  if(var == "Mass") {
    h_unfolded_p = (TH1*) unfold_svd_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_svd_m.Hreco(RooUnfold::kCovariance);
  }
  else if(var == "dRZZ"){
    h_unfolded_p = (TH1*) unfold_bayes4_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes4_m.Hreco(RooUnfold::kCovariance);
  }
  else{ 
    h_unfolded_p = (TH1*) unfold_bayes_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes_m.Hreco(RooUnfold::kCovariance);
  }
 
  h_unfolded_p->Draw();
  h_unfolded_m->Draw("same");

  double max =0;
  double min = 0;
  TH1 * h_max = (TH1*) h_unfolded_p->Clone("h_max");
  TH1 * h_min = (TH1*) h_unfolded_p->Clone("h_max");

   for(int i = 1; i<10; i++){
      max = 0;
      min = 0;
     
      max = TMath::Max(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
      min = TMath::Min(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
      cout << "max " << max << " min " << min<< endl;
      
      h_max->SetBinContent(i,max); 
      h_min->SetBinContent(i,min);
     
    }


  UnfHistoName_p = "ZZTo"+ fs + "_" + var + "_p";
  UnfHistoName_m = "ZZTo"+ fs + "_" + var + "_m";

  output->cd();
  h_max->Write(UnfHistoName_p.c_str());
  h_min->Write(UnfHistoName_m.c_str());
  output->Close();

 }


void DoIrrBkgSystematic(string var = "Mass", string fs = "4e", bool mad = 1,bool tightregion = 0){
  
  string tightfr;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";

#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
 
  filePath = "../../../";
  if(mad ==1){
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Mad/").c_str());  
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Mad.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + "_IrrBkg.root";
  }
  else{
      system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Pow/").c_str());       
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Pow.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + "_IrrBkg.root";
  }
  dataFileName = "../" +var + "_test/DataToUnfold_syst.root";
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
  RooUnfoldBayes unfold_bayes_p(&response, h_data_p, 4);
  RooUnfoldBayes unfold_bayes_m(&response, h_data_m, 4); 
  RooUnfoldBayes unfold_bayes4_p(&response, h_data_p, 4);
  RooUnfoldBayes unfold_bayes4_m(&response, h_data_m, 4);
  
  if(var == "Mass") {
    h_unfolded_p = (TH1*) unfold_svd_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_svd_m.Hreco(RooUnfold::kCovariance);
  }
  else if(var == "dRZZ"){
    h_unfolded_p = (TH1*) unfold_bayes4_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes4_m.Hreco(RooUnfold::kCovariance);
  }
  else{ 
    h_unfolded_p = (TH1*) unfold_bayes_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes_m.Hreco(RooUnfold::kCovariance);
  }
  
  double max =0;
  double min = 0;
  TH1 * h_max = (TH1*) h_unfolded_p->Clone("h_max");
  TH1 * h_min = (TH1*) h_unfolded_p->Clone("h_max");

   for(int i = 1; i<10; i++){
      max = 0;
      min = 0;
     
      max = TMath::Max(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
      min = TMath::Min(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
            
      h_max->SetBinContent(i,max); 
      h_min->SetBinContent(i,min);
     
    }


  UnfHistoName_p = "ZZTo"+ fs + "_" + var + "_p";
  UnfHistoName_m = "ZZTo"+ fs + "_" + var + "_m";

  output->cd();
  h_max->Write(UnfHistoName_p.c_str());
  h_min->Write(UnfHistoName_m.c_str());
  output->Close();

}


void DoRedBkgSystematic(string var = "Mass", string fs = "4e", bool mad = 1,bool tightregion = 0){
  
  string tightfr;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";

#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
 
 filePath = "../../../";
  if(mad ==1){
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Mad/").c_str());  
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Mad.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + "_RedBkg.root";
  }
  else{
      system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Pow/").c_str());       
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Pow.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + "_RedBkg.root";
  }
  dataFileName = "../" +var + "_test/DataToUnfold_syst.root";
 
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
  RooUnfoldBayes unfold_bayes_p(&response, h_data_p, 4);
  RooUnfoldBayes unfold_bayes_m(&response, h_data_m, 4);
  RooUnfoldBayes unfold_bayes4_p(&response, h_data_p, 4);
  RooUnfoldBayes unfold_bayes4_m(&response, h_data_m, 4);
  
  if(var == "Mass") {
    h_unfolded_p = (TH1*) unfold_svd_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_svd_m.Hreco(RooUnfold::kCovariance);
  } 
  else if(var == "dRZZ") {
    h_unfolded_p = (TH1*) unfold_bayes4_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes4_m.Hreco(RooUnfold::kCovariance);
  }
  else{ 
    h_unfolded_p = (TH1*) unfold_bayes_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes_m.Hreco(RooUnfold::kCovariance);
  }

  double max =0;
  double min = 0;
  TH1 * h_max = (TH1*) h_unfolded_p->Clone("h_max");
  TH1 * h_min = (TH1*) h_unfolded_p->Clone("h_max");

   for(int i = 1; i<10; i++){
      max = 0;
      min = 0;
     
      max = TMath::Max(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
      min = TMath::Min(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
        
      h_max->SetBinContent(i,max); 
      h_min->SetBinContent(i,min);
     
    }


  UnfHistoName_p = "ZZTo"+ fs + "_" + var + "_p";
  UnfHistoName_m = "ZZTo"+ fs + "_" + var + "_m";

  output->cd();
  h_max->Write(UnfHistoName_p.c_str());
  h_min->Write(UnfHistoName_m.c_str());
  output->Close();

}

void DoUnfOverGenSystematic(string var = "Mass", string fs = "4e", bool mad = 1,bool tightregion = 0){
  
  string tightfr;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";
  
  filePath = "../../../";
  
  if(mad ==1){ 
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Mad/").c_str());  
    matrixFileName = "../" + var + "_test/weightedMatrices"+tightfr+"_Mad.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + "_UnfDataOverGenMC.root";
  }  
  else{
      system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Pow/").c_str());       
    matrixFileName = "../" + var + "_test/weightedMatrices"+tightfr+"_Pow.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + "_UnfDataOverGenMC.root";
  }
  dataFileName = "../" +var + "_test/DataToUnfold.root";
  matrix   = new TFile(matrixFileName.c_str());
  data     = new TFile(dataFileName.c_str());
  output   = new TFile(outputFileName.c_str(), "UPDATE");

 
  std::cout<<"matrixFileName "<<matrixFileName.c_str()<<std::endl;
  matrixName    = "ResMat_qqggJJ_" + var + "_ZZTo" + fs + "_st_01";
  histoName     = var +"_qqggJJ_ZZTo" + fs + "_st_01";
  histoNameGen  = var + "Gen_qqggJJ_ZZTo" + fs + "_st_01";
  histoName_unf = "DataminusBkg_" + var + "_ZZTo"+fs;
 

  h_measured = (TH1*) matrix->Get(histoName.c_str());  
  h_true = (TH1*) matrix->Get(histoNameGen.c_str());
  h_Resmat = (TH2*)matrix->Get(matrixName.c_str()); 
  h_measured_unf = (TH1*) data->Get(histoName_unf.c_str());  

  RooUnfoldResponse response(h_measured, h_true, h_Resmat, "response", "response"); 
  RooUnfoldSvd      unfold_svd(&response, h_measured_unf, 4); 
  RooUnfoldBayes    unfold_bayes(&response, h_measured_unf,4);   
  RooUnfoldBayes    unfold_bayes4(&response, h_measured_unf,4); 

  if(var == "Mass") h_unfolded= (TH1*) unfold_svd.Hreco(RooUnfold::kCovariance);
  else if(var == "dRZZ") h_unfolded= (TH1*) unfold_bayes4.Hreco(RooUnfold::kCovariance);
  else  h_unfolded= (TH1*) unfold_bayes.Hreco(RooUnfold::kCovariance);
 
  unfHistoName = "ZZTo"+ fs +"_" + var;
  
  output->cd();   
  h_unfolded->Write(unfHistoName.c_str());
  output->Close();

}

void DoJERSystematic(string var = "Jets", string fs = "4e", bool mad = 1,bool tightregion = 0){
  
  string tightfr;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";

#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  filePath = "../../../";
  if(mad ==1){
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Mad/").c_str());  
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_JESJER_Mad.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + "_JER.root";
  }
  else{ 
     system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Pow/").c_str());    
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_JESJER_Pow.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + "_JER.root";
  }
  dataFileName = "../" +var + "_test/DataToUnfold.root";
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
  
  RooUnfoldBayes unfold_bayes_p(&response_p, h_measured_unf, 4);
  RooUnfoldBayes unfold_bayes_m(&response_m, h_measured_unf, 4);

  RooUnfoldBayes unfold_bayes4_p(&response_p, h_measured_unf, 4);
  RooUnfoldBayes unfold_bayes4_m(&response_m, h_measured_unf, 4);
    
  if(var == "Mass") {
    h_unfolded_p = (TH1*) unfold_svd_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_svd_m.Hreco(RooUnfold::kCovariance);
  }
  else if(var == "dRZZ"){ 
    h_unfolded_p = (TH1*) unfold_bayes4_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes4_m.Hreco(RooUnfold::kCovariance);
  }
  else{ 
    h_unfolded_p = (TH1*) unfold_bayes_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes_m.Hreco(RooUnfold::kCovariance);
  }

  double max =0;
  double min = 0;
  TH1 * h_max = (TH1*) h_unfolded_p->Clone("h_max");
  TH1 * h_min = (TH1*) h_unfolded_p->Clone("h_max");

   for(int i = 1; i<10; i++){
      max = 0;
      min = 0;
     
      max = TMath::Max(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
      min = TMath::Min(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
            
      h_max->SetBinContent(i,max); 
      h_min->SetBinContent(i,min);
     
    }


  UnfHistoName_p = "ZZTo"+ fs + "_" + var + "_p";
  UnfHistoName_m = "ZZTo"+ fs + "_" + var + "_m";

  output->cd();
  h_max->Write(UnfHistoName_p.c_str());
  h_min->Write(UnfHistoName_m.c_str());
  output->Close();
   
}

void DoJESSystematic_ModMat(string var = "Jets", string fs = "4e", bool mad = 1,bool tightregion = 0){
  
  string tightfr;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";

#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  filePath = "../../../";
  if(mad ==1){
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Mad/").c_str());  
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_JESJER_Mad.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + "_JES_ModMat.root";
  }
  else{ 
     system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Pow/").c_str());     
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_JESJER_Pow.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + "_JES_ModMat.root";
  }
  dataFileName = "../" +var + "_test/DataToUnfold.root";
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
  
  RooUnfoldBayes unfold_bayes_p(&response_p, h_measured_unf, 4);
  RooUnfoldBayes unfold_bayes_m(&response_m, h_measured_unf, 4);

  RooUnfoldBayes unfold_bayes4_p(&response_p, h_measured_unf, 4);
  RooUnfoldBayes unfold_bayes4_m(&response_m, h_measured_unf, 4);
    
  if(var == "Mass") {
    h_unfolded_p = (TH1*) unfold_svd_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_svd_m.Hreco(RooUnfold::kCovariance);
  } 
  else if(var == "dRZZ"){ 
    h_unfolded_p = (TH1*) unfold_bayes4_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes4_m.Hreco(RooUnfold::kCovariance);
  }
  else { 
    h_unfolded_p = (TH1*) unfold_bayes_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes_m.Hreco(RooUnfold::kCovariance);
  }
  
  double max =0;
  double min = 0;
  TH1 * h_max = (TH1*) h_unfolded_p->Clone("h_max");
  TH1 * h_min = (TH1*) h_unfolded_p->Clone("h_max");

   for(int i = 1; i<10; i++){
      max = 0;
      min = 0;
     
      max = TMath::Max(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
      min = TMath::Min(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
           
      h_max->SetBinContent(i,max); 
      h_min->SetBinContent(i,min);
     
    }


  UnfHistoName_p = "ZZTo"+ fs + "_" + var + "_p";
  UnfHistoName_m = "ZZTo"+ fs + "_" + var + "_m";

  output->cd();
  h_max->Write(UnfHistoName_p.c_str());
  h_min->Write(UnfHistoName_m.c_str());
  output->Close();

}

void DoJESSystematic_ModData(string var = "Jets", string fs = "4e",bool mad = 1,bool tightregion = 0){
  
  string tightfr;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";

#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  
  filePath = "../../../";
  if(mad ==1){
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Mad/").c_str());  
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Mad.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + "_JES_ModData.root";
  }
  else{
     system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Pow/").c_str());     
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Pow.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + "_JES_ModData.root";
  }
  dataFileName = "../" +var + "_test/DataToUnfold_JES.root";
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

  RooUnfoldBayes unfold_bayes_p(&response, h_data_p, 4);
  RooUnfoldBayes unfold_bayes_m(&response, h_data_m, 4);

  RooUnfoldBayes unfold_bayes4_p(&response, h_data_p, 4);
  RooUnfoldBayes unfold_bayes4_m(&response, h_data_m, 4);
 
  if(var == "Mass") {
    h_unfolded_p = (TH1*) unfold_svd_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_svd_m.Hreco(RooUnfold::kCovariance);
  } 
  else if(var == "dRZZ"){ 
    h_unfolded_p = (TH1*) unfold_bayes4_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes4_m.Hreco(RooUnfold::kCovariance);
  }
  else{ 
    h_unfolded_p = (TH1*) unfold_bayes_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes_m.Hreco(RooUnfold::kCovariance);
  }

  double max =0;
  double min = 0;
  TH1 * h_max = (TH1*) h_unfolded_p->Clone("h_max");
  TH1 * h_min = (TH1*) h_unfolded_p->Clone("h_max");

   for(int i = 1; i<10; i++){
      max = 0;
      min = 0;
     
      max = TMath::Max(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
      min = TMath::Min(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
            
      h_max->SetBinContent(i,max); 
      h_min->SetBinContent(i,min);
     
    }


  UnfHistoName_p = "ZZTo"+ fs + "_" + var + "_p";
  UnfHistoName_m = "ZZTo"+ fs + "_" + var + "_m";

  output->cd();
  h_max->Write(UnfHistoName_p.c_str());
  h_min->Write(UnfHistoName_m.c_str());
  output->Close();

}

void DoSFSystematic(string var = "Jets", string fs = "4e", bool mad = 1,bool tightregion = 0, bool corr = 0){
  
  string tightfr;
  string errtype;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";
  if(corr == 0)  errtype = "Sq";
  else errtype = "";

#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  filePath = "../../../";
  if(mad ==1){
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Mad/").c_str());  
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_SF_Mad.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + "_SF"+errtype+".root";
  }
  else{ 
     system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Pow/").c_str());     
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_SF_Pow.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + "_SF"+errtype+".root";
  }
  dataFileName = "../" +var + "_test/DataToUnfold.root";
  matrix = new TFile(matrixFileName.c_str());
  data = new TFile(dataFileName.c_str());
  output = new TFile(outputFileName.c_str(), "UPDATE");

  if(corr == 0){

    matrix_p       = "ResMat_qqggJJ_" + var + "_ZZTo" + fs + "_SFErrSqPlus_01";
    histoReco_p    = var + "_qqggJJ_ZZTo" + fs +"_SFErrSqPlus_01";
    histoGen_p     = var + "Gen_qqggJJ_ZZTo" + fs +"_SFErrSqPlus_01";
    matrix_m       = "ResMat_qqggJJ_" + var + "_ZZTo" + fs +"_SFErrSqMinus_01";  
    histoReco_m    = var + "_qqggJJ_ZZTo" + fs +"_SFErrSqMinus_01";
    histoGen_m     = var + "Gen_qqggJJ_ZZTo" + fs +"_SFErrSqMinus_01"; 
    histoName_unf  = "DataminusBkg_" + var + "_ZZTo"+fs;
  
  }
  else{
    matrix_p = "ResMat_qqggJJ_" + var + "_ZZTo" + fs + "_SFErrPlus_01";
    histoReco_p = var + "_qqggJJ_ZZTo" + fs +"_SFErrPlus_01";
    histoGen_p = var + "Gen_qqggJJ_ZZTo" + fs +"_SFErrPlus_01";
    matrix_m = "ResMat_qqggJJ_" + var + "_ZZTo" + fs +"_SFErrMinus_01";  
    histoReco_m = var + "_qqggJJ_ZZTo" + fs +"_SFErrMinus_01";
    histoGen_m = var + "Gen_qqggJJ_ZZTo" + fs +"_SFErrMinus_01"; 
    histoName_unf = "DataminusBkg_" + var + "_ZZTo"+fs; 
    }
  
  std::cout<<"matrix name "<<matrixFileName.c_str()<<std::endl;
  h_measured_p = (TH1*) matrix->Get(histoReco_p.c_str());
  std::cout<<histoReco_p.c_str()<<std::endl;  
  h_true_p = (TH1*) matrix->Get(histoGen_p.c_str());
  std::cout<<histoGen_p.c_str()<<std::endl;  
  h_Resmat_p = (TH2*)matrix->Get(matrix_p.c_str()); 
  std::cout<<matrix_p.c_str()<<std::endl;  
  h_measured_m = (TH1*) matrix->Get(histoReco_m.c_str());  
  h_true_m = (TH1*) matrix->Get(histoGen_m.c_str());
  h_Resmat_m = (TH2*)matrix->Get(matrix_m.c_str()); 
  h_measured_unf = (TH1*) data->Get(histoName_unf.c_str());


  RooUnfoldResponse response_p(h_measured_p, h_true_p, h_Resmat_p, "response_p", "response_p"); 
  RooUnfoldResponse response_m(h_measured_m, h_true_m, h_Resmat_m, "response_m", "response_m"); 
  RooUnfoldSvd unfold_svd_p(&response_p, h_measured_unf, 4);
  RooUnfoldSvd unfold_svd_m(&response_m, h_measured_unf, 4);
  RooUnfoldBayes unfold_bayes_p(&response_p, h_measured_unf, 4);
  RooUnfoldBayes unfold_bayes_m(&response_m, h_measured_unf, 4);
  RooUnfoldBayes unfold_bayes4_p(&response_p, h_measured_unf, 4);
  RooUnfoldBayes unfold_bayes4_m(&response_m, h_measured_unf, 4);
    
  if(var == "Mass") {
    h_unfolded_p = (TH1*) unfold_svd_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_svd_m.Hreco(RooUnfold::kCovariance);
  } 
  else if(var == "dRZZ"){ 
    h_unfolded_p = (TH1*) unfold_bayes4_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes4_m.Hreco(RooUnfold::kCovariance);
  }
  else{ 
    h_unfolded_p = (TH1*) unfold_bayes_p.Hreco(RooUnfold::kCovariance);
    h_unfolded_m = (TH1*) unfold_bayes_m.Hreco(RooUnfold::kCovariance);
  }
 
  double max =0;
  double min = 0;
  TH1 * h_max = (TH1*) h_unfolded_p->Clone("h_max");
  TH1 * h_min = (TH1*) h_unfolded_p->Clone("h_max");

   for(int i = 1; i<10; i++){
      max = 0;
      min = 0;
     
      max = TMath::Max(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
      min = TMath::Min(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
            
      h_max->SetBinContent(i,max); 
      h_min->SetBinContent(i,min);
     
    }


  UnfHistoName_p = "ZZTo"+ fs + "_" + var + "_p";
  UnfHistoName_m = "ZZTo"+ fs + "_" + var + "_m";

  output->cd();
  h_max->Write(UnfHistoName_p.c_str());
  h_min->Write(UnfHistoName_m.c_str());
  output->Close();
}

void  DoAllSystematics(string var = "Mass", bool mad =1, bool tightregion =0){
  
  DoMCGenSystematic(var.c_str(), "4e",mad,tightregion);
  DoMCGenSystematic(var.c_str(), "4m",mad,tightregion);
  DoMCGenSystematic(var.c_str(), "2e2m",mad,tightregion);
  DoQqggSystematic(var.c_str(), "4e",mad,tightregion); 
  DoQqggSystematic(var.c_str(), "4m",mad,tightregion);
  DoQqggSystematic(var.c_str(), "2e2m",mad,tightregion);
  DoIrrBkgSystematic(var.c_str(), "4e",mad,tightregion); 
  DoIrrBkgSystematic(var.c_str(), "4m",mad,tightregion); 
  DoIrrBkgSystematic(var.c_str(), "2e2m",mad,tightregion);
  DoRedBkgSystematic(var.c_str(), "4e",mad,tightregion); 
  DoRedBkgSystematic(var.c_str(), "4m",mad,tightregion); 
  DoRedBkgSystematic(var.c_str(), "2e2m",mad,tightregion); 
  
  //  DoUnfOverGenSystematic(var.c_str(), "4e",mad,tightregion);
  //DoUnfOverGenSystematic(var.c_str(), "4m",mad,tightregion);
  //sDoUnfOverGenSystematic(var.c_str(), "2e2m",mad,tightregion);
  
  if(var == "Jets"||var == "Mjj"||var == "Deta" || var == "Jets_Central"||var == "Mjj_Central"||var == "Deta_Central"||var == "PtJet1"||var == "PtJet2"||var == "EtaJet1"||var == "EtaJet2"){
    DoJERSystematic(var.c_str(), "4e",mad,tightregion); 
    DoJERSystematic(var.c_str(), "4m",mad,tightregion);
    DoJERSystematic(var.c_str(), "2e2m",mad,tightregion);
  
    DoJESSystematic_ModMat(var.c_str(), "4e",mad,tightregion); 
    DoJESSystematic_ModMat(var.c_str(), "4m",mad,tightregion);
    DoJESSystematic_ModMat(var.c_str(), "2e2m",mad,tightregion);  
  
    DoJESSystematic_ModData(var.c_str(), "4e",mad,tightregion); 
    DoJESSystematic_ModData(var.c_str(), "4m",mad,tightregion);
    DoJESSystematic_ModData(var.c_str(), "2e2m",mad,tightregion); 
  }
  
  DoSFSystematic(var.c_str(),"4e",mad,tightregion,0); 
  DoSFSystematic(var.c_str(),"4m",mad,tightregion,0); 
  DoSFSystematic(var.c_str(),"2e2m",mad,tightregion,0);
  // DoSFSystematic(var.c_str(),"4e",mad,tightregion,1);  
  // DoSFSystematic(var.c_str(),"4m",mad,tightregion,1); 
  // DoSFSystematic(var.c_str(),"2e2m",mad,tightregion,1); 
}

void PlotResults(string var = "Jets", string fs = "4e", string syst = "MCgen", bool mad =1,bool tightregion =0, string date = "test"){
  
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  setTDRStyle();

  string finalstate;
  if(fs == "4m") finalstate = "4#mu";
  else if(fs == "2e2m") finalstate = "2e2#mu";
  else finalstate = fs;
  
  string label = finalstate + " channel";
  int iPeriod = 4; 
  int iPos = 11; 
  writeExtraText = true;    
  extraText  = "Preliminary";
  extraText2  = label.c_str();
 
  string tightfr;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";
  
  TFile * unfSyst;
  string systFileName;
  string systHistoName_p;
  string systHistoName_m;
  TH1 *h_true_othMC;
 
  system(("mkdir ~/www/VBS/"+date+"/"+ var).c_str());
  system(("cp ~/www/VBS/index.php ~/www/VBS/"+date+"/"+ var).c_str()); 
  system(("mkdir ~/www/VBS/"+date+"/"+ var+"/Systematics").c_str()); 
  system(("cp ~/www/VBS/index.php ~/www/VBS/"+date+"/"+ var+"/Systematics").c_str());


  filePath = "../../../";
  if(mad ==1){
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Mad/").c_str());  
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Pow.root";
    dataFileName = filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + ".root";
    systFileName =  filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + "_"+syst + ".root";
    MCgen = "_Mad";
  }
  else{
     system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Pow/").c_str());     
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Mad.root";
    dataFileName = filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + ".root";
    systFileName =  filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + "_"+syst + ".root";
    MCgen = "_Pow";
  }
  data = new TFile(dataFileName.c_str());
  unfSyst = new TFile(systFileName.c_str());
  matrix = new TFile(matrixFileName.c_str());
  std::cout<<"HHHE"<<std::endl;
  unfHistoName = "ZZTo"+ fs +"_" + var;
  recoHistoName = "ZZTo"+ fs +"_"+var+"_RECO";
  trueHistoName = "ZZTo"+ fs +"_"+var+"_GEN";
  systHistoName_p = "ZZTo"+ fs +"_" + var + "_p";
  systHistoName_m =  "ZZTo"+ fs +"_" + var + "_m";
  histoNameGen = var + "Gen_qqggJJ_ZZTo" + fs + "_st_01";
  
  h_true_othMC = (TH1*) matrix->Get(histoNameGen.c_str());
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
    h_max = (TH1*)h_unfolded_p->Clone("");
    h_min = (TH1*)h_unfolded_p->Clone("");
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
  else if(syst=="qqgg"||syst=="IrrBkg"||syst=="RedBkg"||syst=="JER"||syst=="JES_ModMat"||syst=="JES_ModData"||syst=="SF"||syst=="SFSq"){
    h_unfolded_p = (TH1*) unfSyst->Get(systHistoName_p.c_str());
    h_unfolded_m = (TH1*) unfSyst->Get(systHistoName_m.c_str());
    h_max = (TH1*)h_unfolded_p->Clone("");;
    h_min = (TH1*)h_unfolded_p->Clone("");;
    for(int i = 1; i<9; i++){
      
      syst_percentage[i-1]=0;
      
      Double_t vec[] ={h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i),h_unfolded->GetBinContent(i)};
      min =vec[0];
      max =0;
    
      for(int j=0;j<3;j++)
	{
	  if(vec[j]<min) min = vec[j];
	  if(vec[j]>max) max = vec[j];
	}
      
      h_max->SetBinContent(i,max); 
      h_min->SetBinContent(i,min); 
      syst_percentage[i-1]=0.5*(max-min)*100/(h_unfolded->GetBinContent(i)); 
     
      if(max!=0.) std::cout  << i << "-bin: min = "<< min << " max = " << max << " " << syst_percentage[i-1]  << "% " << syst.c_str() << " systematic uncertainty" << std::endl;
    }
  } 
  else std::cout << " Wrong systematic uncertainty!!!!" << std::endl;

  float M = syst_percentage[0];
  float m = syst_percentage[0];
  for(int i = 0; i<8; i++) {
      if(syst_percentage[i] > M) M = syst_percentage[i];
      if(syst_percentage[i] < m) m = syst_percentage[i];
     }

  std::cout << "Max syst = " << M << "%   Min syst = " << m << "%" << endl;

  double totUnc_min = fabs(h_unfolded->Integral()-h_min->Integral())/h_unfolded->Integral();
  double totUnc_max = fabs(h_unfolded->Integral()-h_max->Integral())/h_unfolded->Integral();
  cout << "tot_integral = " << h_unfolded->Integral() << " + " <<  totUnc_max*100 << "% - " << 
    totUnc_min*100 << "%" <<endl;

  TCanvas *c = new TCanvas ("c","c"); 
  TLegend *leg = new TLegend(0.55,0.65,0.45,0.85); 
  
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);

  string XaxisTitle;
  float max_y =0;

 if(var == "Mass"){
    XaxisTitle = "m_{4l} [GeV]";
    YaxisTitle = "Events/GeV";  
    if (tightregion == 1)max_y = 125;
    else max_y = 350;
  }
  else if(var == "Jets")  {
    XaxisTitle = "N jets (|#eta^{jet}| < 4.7)";
    YaxisTitle = "Events";
    if (tightregion == 1)  max_y = 300;
    else max_y = 700; 
  }
  else if(var == "Jets_Central")  {
    XaxisTitle = "N jets (|#eta^{jet}| < 2.4)"; 
    YaxisTitle = "Events"; 
    if (tightregion == 1)  max_y = 300;
    else max_y = 700; 
  }
  else if(var == "Mjj") {
    XaxisTitle ="m_{jj} (|#eta^{jet}| < 4.7) [GeV]";
    YaxisTitle = "Events/GeV";  
    max_y =  h_true->GetBinContent(1)+h_true->GetBinContent(1)*1;
  }
  else if(var == "Mjj_Central") {
    XaxisTitle ="m_{jj} (|#eta^{jet}| < 2.4) [GeV]"; 
    YaxisTitle = "Events/GeV"; 
    max_y =  h_true->GetBinContent(1)+h_true->GetBinContent(1)*1;
  }
  else if( var == "Deta"){
    XaxisTitle ="#Delta#eta_{jj} (|#eta^{jet}| < 4.7)";
    YaxisTitle = "Events/(bin width)"; 
    if(mad == 1) max_y =  h_true->GetBinContent(1)+h_true->GetBinContent(1)*0.70; 
    else max_y =  h_true_othMC->GetBinContent(1)+h_true_othMC->GetBinContent(1)*0.70;  
  }
  else if(var == "Deta_Central"){ 
    XaxisTitle ="#Delta#eta_{jj} (|#eta^{jet}| < 2.4)";
    YaxisTitle = "Events/(bin width)"; 
    if(mad == 1) max_y =  h_true->GetBinContent(1)+h_true->GetBinContent(1)*0.70; 
    else max_y =  h_true_othMC->GetBinContent(1)+h_true_othMC->GetBinContent(1)*0.70; 
  }
  else if(var =="PtJet1") {
    XaxisTitle = "p_{T}^{jet1} [GeV]";
    YaxisTitle = "Events/GeV"; 
    max_y =  h_true->GetBinContent(1)+h_true->GetBinContent(1)*0.70;
  }
  else if(var =="PtJet2") { 
    XaxisTitle = "p_{T}^{jet2} [GeV]"; 
    YaxisTitle = "Events/GeV"; 
    if(mad == 1) max_y =  h_true->GetBinContent(1)+h_true->GetBinContent(1)*0.70; 
    else max_y =  h_true_othMC->GetBinContent(1)+h_true_othMC->GetBinContent(1)*0.70; 
   
  }
  else if(var =="EtaJet1"){  
    XaxisTitle = "|#eta^{jet1}|";
    YaxisTitle = "Events/(bin width)"; 
    if(mad == 1) max_y =  h_true->GetBinContent(1)+h_true->GetBinContent(1)*0.70; 
    else max_y =  h_true_othMC->GetBinContent(1)+h_true_othMC->GetBinContent(1)*0.70;  
  }
  else if(var =="EtaJet2"){  
    XaxisTitle = "|#eta^{jet2}|";
    YaxisTitle = "Events/(bin width)";  
    if(mad == 1) max_y =  h_true->GetBinContent(1)+h_true->GetBinContent(1)*0.70; 
    else max_y =  h_true_othMC->GetBinContent(1)+h_true_othMC->GetBinContent(1)*0.70;  
  }
  else if(var =="dRZZ") { 
    XaxisTitle = "#DeltaR(Z_{1},Z_{2})"; 
    YaxisTitle = "Events/(bin width)"; 
    max_y =  h_true->GetBinContent(4)+h_true->GetBinContent(4);
  }
   
  h_max->Scale(1,"width");
  h_min->Scale(1,"width");
  h_unfolded->Scale(1,"width");
  h_true->Scale(1,"width");
  h_true_othMC->Scale(1,"width");
  
 
  int binmax = h_max->GetMaximumBin();
  double binwidth = h_max->GetBinWidth(binmax);  
  
  h_max->Draw("HIST");
  h_max->SetFillColor(2);
  h_max->SetFillStyle(3001);
  h_max->SetLineColor(10);
  h_max->SetTitle("");
  h_max->GetXaxis()->SetRange(0,25);
  h_max->GetXaxis()->SetTitle(XaxisTitle.c_str());  
  h_max->GetYaxis()->SetTitle(YaxisTitle.c_str()); 
  if(var =="Mjj"||var =="Mjj_Central") h_max->GetYaxis()->SetTitleOffset(1.8);
  else  h_max->GetYaxis()->SetTitleOffset(1.4);
  h_max->GetXaxis()->SetTitleOffset(1.2);
  if(var == "Mjj_Central" && fs == "2e2m") h_max->SetMaximum(0.1);  
  else h_max->SetMaximum(max_y/binwidth);
  h_max->SetMinimum(0.);
  h_min->Draw("HIST SAME");
  h_min->SetFillColor(10);
  h_min->SetLineColor(10);
  
  h_unfolded->Draw("E SAME"); 
  h_unfolded->SetLineColor(1); 
  //h_unfolded->SetLineStyle(2);
  h_unfolded->SetLineWidth(1);
  h_unfolded->SetMarkerStyle(20);
  h_unfolded->SetMarkerColor(1);
  h_unfolded->SetMarkerSize(1);

  h_true->SetLineColor(8);
  h_true->SetLineStyle(1);  
  h_true->SetLineWidth(1);
  h_true->Draw("E HIST SAME");
  h_true_othMC->SetLineColor(9);
  h_true_othMC->SetLineStyle(1);  
  h_true_othMC->SetLineWidth(1);
  h_true_othMC->Draw("E HIST SAME");
  h_unfolded->Draw("E SAME"); 
  gPad->RedrawAxis();
  
  lumiTextSize     = 0.7;
  cmsTextSize      = 0.7;
  extraOverCmsTextSize  = 0.80;
  CMS_lumi(c, iPeriod, iPos );

  string syst_uncertainty;
  if(syst == "MCgen") syst_uncertainty = "MC generator syst. unc.";
  else if(syst == "UnfDataOverGenMC") syst_uncertainty = "Data/MC syst. unc.";
  else if(syst=="qqgg") syst_uncertainty = "#sigma_{qq}/#sigma_{gg} syst. unc.";
  else if(syst=="IrrBkg") syst_uncertainty = " Irr. Bkg. syst. unc.";
  else if(syst=="RedBkg") syst_uncertainty = " Red. Bkg. syst. unc.";
  else if(syst=="JER") syst_uncertainty = "JER syst. unc.";
  else if(syst=="JES_ModMat") syst_uncertainty = "JES syst. unc. (modified matrix)";
  else if(syst=="JES_ModData") syst_uncertainty = "JES syst. unc. (modified data)";
  else if(syst=="SFSq") syst_uncertainty = "Scale Factor syst. unc.";
    else if(syst=="SF") syst_uncertainty = "Scale Factor syst. unc. (correlated)";
  else std::cout << "wrong systematic uncerainty!!!!" << std::endl;

  leg->AddEntry(h_unfolded,"Unfolded data","lp"); 
  leg->Draw(); 
  if(mad ==1){
    leg->AddEntry(h_true,"MC truth (MadGraph Set)","l"); 
    leg->Draw("SAME"); 
    leg->AddEntry(h_true_othMC,"MC truth (Powheg Set)","l"); 
    leg->Draw("SAME");
  }
  else{ 
    leg->AddEntry(h_true_othMC,"MC truth (MadGraph Set)","l"); 
    leg->Draw("SAME"); 
    leg->AddEntry(h_true,"MC truth (Powheg Set)","l"); 
    leg->Draw("SAME");
  }
  leg->AddEntry(h_max,syst_uncertainty.c_str(),"f"); 
  leg->Draw("SAME");

  string png = " ~/www/VBS/"+date+"/"+ var+"/Systematics/ZZTo"+fs+"_"+var+"_"+syst+MCgen+tightfr+".png";
  string pdf = " ~/www/VBS/"+date+"/"+ var+"/Systematics/ZZTo"+fs+"_"+var+"_"+syst+MCgen+tightfr+".pdf";
  c->SaveAs(png.c_str());
  c->SaveAs(pdf.c_str());


}

void AllPlots_fs(string var = "Jets", string fs = "4e", bool mad =1,bool tightregion =0,string date = "test"){
  std::cout << "===================================================================================" << std::endl;
  std::cout << "===================================================================================" << std::endl;
  std::cout << "                                 " << fs.c_str() << std::endl;
  std::cout << "===================================================================================" << std::endl;
  PlotResults(var.c_str(),fs.c_str(),"MCgen",mad, tightregion, date.c_str());
  PlotResults(var.c_str(),fs.c_str(),"qqgg",mad, tightregion, date.c_str());
  PlotResults(var.c_str(),fs.c_str(),"IrrBkg",mad, tightregion, date.c_str());
  PlotResults(var.c_str(),fs.c_str(),"RedBkg",mad, tightregion, date.c_str());
  //  PlotResults(var.c_str(),fs.c_str(),"UnfDataOverGenMC",mad,tightregion, date.c_str());

  if(var == "Jets"||var == "Mjj"||var == "Deta"||var == "Jets_Central"||var == "Mjj_Central"||var == "Deta_Central"||var == "PtJet1"||var == "PtJet2"||var == "EtaJet1"||var == "EtaJet2"){
    PlotResults(var.c_str(),fs.c_str(),"JER",mad, tightregion, date.c_str()); 
    PlotResults(var.c_str(),fs.c_str(),"JES_ModMat",mad, tightregion, date.c_str()); 
    PlotResults(var.c_str(),fs.c_str(),"JES_ModData",mad, tightregion, date.c_str()); 
  }
  PlotResults(var.c_str(),fs.c_str(),"SFSq",mad, tightregion, date.c_str()); 
  //  PlotResults(var.c_str(),fs.c_str(),"SF",mad, tightregion, date.c_str());
}

void AllPlots(string var = "Jets", bool mad =1,bool tightregion =0,string date = "test"){
  AllPlots_fs(var.c_str(),"4e",mad, tightregion, date.c_str());
  AllPlots_fs(var.c_str(),"4m",mad, tightregion, date.c_str());
  AllPlots_fs(var.c_str(),"2e2m",mad, tightregion, date.c_str());
}
void AllSystematicPlotsForAllVariables(bool mad =1,bool tightregion =0,string date = "test"){
  
  AllPlots("Jets",mad, tightregion, date);
  AllPlots("Jets_Central",mad, tightregion, date);
  AllPlots("Mjj",mad, tightregion, date);
  AllPlots("Mjj_Central",mad, tightregion, date);
  AllPlots("Deta",mad, tightregion, date);
  AllPlots("Deta_Central",mad, tightregion, date);
  AllPlots("PtJet1",mad, tightregion, date);
  AllPlots("PtJet2",mad, tightregion, date);
  AllPlots("EtaJet1",mad, tightregion, date);
  AllPlots("EtaJet2",mad, tightregion, date);
}
void AllYouNeed(bool mad =1, bool tightregion = 0){
  //Unfold_data_All("Mass",mad,tightregion);
  //DoAllRatios("Mass",mad,tightregion);
  //DoAllSystematics("Mass",mad,tightregion);
  Unfold_data_All("Jets",mad,tightregion);
  DoAllRatios("Jets",mad,tightregion);
  DoAllSystematics("Jets",mad,tightregion);
  Unfold_data_All("Mjj",mad,tightregion);
  DoAllRatios("Mjj",mad,tightregion);
  DoAllSystematics("Mjj",mad,tightregion);
  Unfold_data_All("Deta",mad,tightregion);
  DoAllRatios("Deta",mad,tightregion);
  DoAllSystematics("Deta",mad,tightregion);
  Unfold_data_All("Jets_Central",mad,tightregion);
  DoAllRatios("Jets_Central",mad,tightregion);
  DoAllSystematics("Jets_Central",mad,tightregion);
  Unfold_data_All("Mjj_Central",mad,tightregion);
  DoAllRatios("Mjj_Central",mad,tightregion);
  DoAllSystematics("Mjj_Central",mad,tightregion);
  Unfold_data_All("Deta_Central",mad,tightregion);
  DoAllRatios("Deta_Central",mad,tightregion);
  DoAllSystematics("Deta_Central",mad,tightregion);
}

void PlotUnfoldData4L(string var = "Jets", bool mad =1,bool tightregion =0, string date = "test"){
  
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  setTDRStyle();
  int iPeriod = 2; 
  int iPos = 11; 
  writeExtraText = true;    
  extraText  = "Preliminary";
  extraText2 = "";
  string tightfr;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";
  
  TFile * unfData;
  string unfDataFileName;  
  filePath = "../../../";
  system(("mkdir ~/www/VBS/"+date).c_str());
  system(("cp ~/www/VBS/index.php ~/www/VBS/"+date).c_str());
  system(("mkdir ~/www/VBS/"+date+"/"+ var).c_str()); 
  system(("cp ~/www/VBS/index.php ~/www/VBS/"+date+"/"+ var).c_str());
  
  if(mad ==1){
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Mad/").c_str());  
    unfDataFileName = filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + ".root";
    MCgen = "_Mad";
  }
  else{
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Pow/").c_str());        
    unfDataFileName = filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + ".root";
    MCgen = "_Pow";
  }
  dataFileName = "../" +var + "_test/DataToUnfold.root";
  
  unfData = new TFile(unfDataFileName.c_str());
  data = new TFile(dataFileName.c_str());


  TH1 * h_unfdata4e = (TH1*) unfData->Get(("ZZTo4e_" + var).c_str());
  TH1 * h_unfdata4m = (TH1*) unfData->Get(("ZZTo4m_" + var).c_str());
  TH1 * h_unfdata2e2m = (TH1*) unfData->Get(("ZZTo2e2m_" + var).c_str());
  TH1 * h_data4e = (TH1*) unfData->Get(("ZZTo4e_" + var+"_RECO").c_str());
  TH1 * h_data4m = (TH1*) unfData->Get(("ZZTo4m_" + var+"_RECO").c_str());
  TH1 * h_data2e2m = (TH1*) unfData->Get(("ZZTo2e2m_" + var+"_RECO").c_str());
  TH1 * h_reco4e = (TH1*) unfData->Get(("ZZTo4e_" + var+"_RECO_MC").c_str());
  TH1 * h_reco4m = (TH1*) unfData->Get(("ZZTo4m_" + var+"_RECO_MC").c_str());
  TH1 * h_reco2e2m = (TH1*) unfData->Get(("ZZTo2e2m_" + var+"_RECO_MC").c_str());
  TH1 * h_gen4e = (TH1*) unfData->Get(("ZZTo4e_" + var+"_GEN").c_str());
  TH1 * h_gen4m = (TH1*) unfData->Get(("ZZTo4m_" + var+"_GEN").c_str());;
  TH1 * h_gen2e2m = (TH1*) unfData->Get(("ZZTo2e2m_" + var+"_GEN").c_str());;
  
  
  h_unfdata4e->Add(h_unfdata4m);
  h_unfdata4e->Add(h_unfdata2e2m);
  h_data4e->Add(h_data4m);
  h_data4e->Add(h_data2e2m);
  h_reco4e->Add(h_reco4m);
  h_reco4e->Add(h_reco2e2m);
  h_gen4e->Add(h_gen4m);
  h_gen4e->Add(h_gen2e2m);
  
  if(var == "Mass"){
    XaxisTitle = "m_{4l} [GeV]";
    YaxisTitle = "Events/GeV"; 
  }
  else if(var == "Jets")  {
    XaxisTitle = "N jets (|#eta^{jet}| < 4.7)";
    YaxisTitle = "Events"; 
  }
  else if(var == "Jets_Central")  {
    XaxisTitle = "N jets (|#eta^{jet}| < 2.4)"; 
    YaxisTitle = "Events"; 
  }
  else if(var == "Mjj") {
    XaxisTitle ="m_{jj} (|#eta^{jet}| < 4.7) [GeV]";
    YaxisTitle = "Events/GeV"; 
  }
  else if(var == "Mjj_Central") {
    XaxisTitle ="m_{jj} (|#eta^{jet}| < 2.4) [GeV]"; 
    YaxisTitle = "Events/GeV"; 
  }
  else if( var == "Deta"){
    XaxisTitle ="#Delta#eta_{jj} (|#eta^{jet}| < 4.7)";
    YaxisTitle = "Events/(bin width)"; 
  }
  else if(var == "Deta_Central"){ 
    XaxisTitle ="#Delta#eta_{jj} (|#eta^{jet}| < 2.4)";
    YaxisTitle = "Events/(bin width)"; 
  }
  else if(var =="PtJet1") {
    XaxisTitle = "p_{T}^{jet1} [GeV]";
    YaxisTitle = "Events/GeV"; 
  }
  else if(var =="PtJet2") { 
    XaxisTitle = "p_{T}^{jet2} [GeV]"; 
    YaxisTitle = "Events/GeV"; 
  }
  else if(var =="EtaJet1"){  
    XaxisTitle = "|#eta^{jet1}|";
    YaxisTitle = "Events/(bin width)"; 
  }
  else if(var =="EtaJet2"){  
    XaxisTitle = "|#eta^{jet2}|";
    YaxisTitle = "Events/(bin width)"; 
  }
  else if(var =="dRZZ") { 
    XaxisTitle = "#DeltaR(Z_{1},Z_{2})"; 
    YaxisTitle = "Events/(bin width)"; 
  }
  
  string YaxisTitle2 = "Unfolded/True";
  
  TCanvas *c = new TCanvas ("c","c");
  TLegend *leg = new TLegend(0.80,0.65,0.60,0.85); 
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42); 
  
  c->cd();
  
  TPad  *pad1 = new TPad("pad1","", 0., 0.19, 1.0, 1.0);
  pad1->SetTopMargin (0.10);
  pad1->SetRightMargin (0.10);
  pad1->SetLeftMargin (0.10);
  pad1->Draw();
  
  c->cd();
  TPad  *pad2 = new TPad("pad2", "", 0., 0.0,  1.0, 0.28);
  pad2->SetTopMargin(0.10);
  pad2->SetRightMargin(0.10);
  pad2->SetLeftMargin(0.10); 
  pad2->SetBottomMargin(0.35);
  pad2->Draw(); 
  
  pad1->cd();
  h_gen4e->SetTitle("");
  h_gen4e->GetXaxis()->SetRange(0,25);
  h_gen4e->GetXaxis()->SetTitle(XaxisTitle.c_str()); 
  h_gen4e->GetYaxis()->SetTitle(YaxisTitle.c_str());
  h_gen4e->GetYaxis()->SetTitleOffset(1.4);
  h_gen4e->GetXaxis()->SetLabelOffset(0.5);
  
  h_gen4e->Scale(1,"width");
  h_reco4e->Scale(1,"width");
  h_unfdata4e->Scale(1,"width");
  h_data4e->Scale(1,"width");
  
  float max = 0;
  if(var == "Mass"){
    if (tightregion == 1)max = 4.5;
    else max = 9;
  }
  else if(var == "Jets" || var == "Jets_Central"){
    if (tightregion == 1)  max = 450;
    else max = 900;
  } 
  else if(var == "Mjj"|| var == "Mjj_Central") max =  h_gen4e->GetBinContent(1)+h_gen4e->GetBinContent(1)*0.20+h_unfdata4e->GetBinError(1); 
  else if(var == "Deta"|| var == "Deta_Central" ) max =  h_gen4e->GetBinContent(1)+h_gen4e->GetBinContent(1)*0.15+h_unfdata4e->GetBinError(1); 
  else if(var == "PtJet1") max =  h_gen4e->GetBinContent(1)+h_gen4e->GetBinContent(1)*0.15+h_unfdata4e->GetBinError(1); 
  else if(var == "PtJet2") max =  h_gen4e->GetBinContent(1)+h_gen4e->GetBinContent(1)*0.15+h_unfdata4e->GetBinError(1); 
  else if(var == "EtaJet1") max =  h_gen4e->GetBinContent(1)+h_gen4e->GetBinContent(1)*0.15+h_unfdata4e->GetBinError(1);
  else if(var == "EtaJet2") max =  h_gen4e->GetBinContent(1)+h_gen4e->GetBinContent(1)*0.15+h_unfdata4e->GetBinError(1);
  else if(var == "dRZZ") max =  h_gen4e->GetBinContent(4)+h_gen4e->GetBinContent(4)+h_unfdata4e->GetBinError(1);
  else  max =  h_gen4e->GetBinContent(1)+h_gen4e->GetBinContent(1)*0.15+h_unfdata4e->GetBinError(1);
  
  h_gen4e->SetMaximum(max);
  h_gen4e->SetMinimum(0); 
  h_gen4e->SetLineColor(kBlue); 
  h_gen4e->SetLineWidth(1);
  h_gen4e->SetLineStyle(7);
  h_gen4e->Draw("HIST E");
  h_reco4e->SetLineColor(kRed);
  h_reco4e->SetLineStyle(7);
  h_reco4e->SetLineWidth(1);
  h_reco4e->Draw("HIST E SAME");
  h_unfdata4e->SetLineColor(kBlue);
  h_unfdata4e->SetLineWidth(1);
  h_unfdata4e->SetMarkerColor(kBlue);
  h_unfdata4e->SetMarkerStyle(8);
  h_unfdata4e->Draw("E SAME"); 
  h_data4e->Draw("HIST E SAME");
  h_data4e->SetMarkerColor(2);
  h_data4e->SetMarkerStyle(8);
  h_data4e->SetLineColor(2);
  h_data4e->SetLineWidth(1);
  
  
  leg->AddEntry(h_gen4e,"MC truth","l"); 
  leg->Draw("SAME");
  leg->AddEntry(h_reco4e,"MC reco","l"); 
  leg->Draw("SAME");
  leg->AddEntry(h_unfdata4e,"unfolded data","lep"); 
  leg->Draw(); 
  leg->AddEntry(h_data4e,"data","lep"); 
  leg->Draw("SAME");
  
  TH1 * h_unfdata4e_r = (TH1*) h_unfdata4e->Clone();
  
  for(int k =1;k<9;k++){
    float unf=0;
    float tr = 0;
    float ratio =0;
    float err_unf=0;
    float err_tr = 0;
    float err_ratio =0;
    unf = h_unfdata4e->GetBinContent(k);
    tr = h_gen4e->GetBinContent(k); 
    err_unf = h_unfdata4e->GetBinError(k);
    err_tr = h_gen4e->GetBinError(k);
    ratio = unf/tr;
    err_ratio = sqrt((err_tr/tr)*(err_tr/tr)+(err_unf/unf)*(err_unf/unf));
    h_unfdata4e_r->SetBinContent(k,ratio); 
    h_unfdata4e_r->SetBinError(k,err_ratio);
  }
  
  
  pad2->cd();  
  h_unfdata4e_r->SetTitle("");
  h_unfdata4e_r->SetMarkerColor(1);
  h_unfdata4e_r->SetLineColor(1); 
  h_unfdata4e_r->GetYaxis()->SetLabelSize(0.10);  
  h_unfdata4e_r->GetYaxis()->SetNdivisions(306); 
  h_unfdata4e_r->GetXaxis()->SetLabelSize(0.13); 
  h_unfdata4e_r->GetXaxis()->SetLabelOffset(0.05);
  h_unfdata4e_r->GetXaxis()->SetTitle(XaxisTitle.c_str()); 
  h_unfdata4e_r->GetYaxis()->SetTitle(YaxisTitle2.c_str());
  h_unfdata4e_r->GetXaxis()->SetTitleOffset(1.2); 
  h_unfdata4e_r->GetYaxis()->SetTitleOffset(0.5);
  h_unfdata4e_r->GetYaxis()->SetTitleSize(0.10);
  h_unfdata4e_r->GetXaxis()->SetTitleSize(0.13);
  
  if(var == "Jets" || var == "Jets_Central"){
    h_unfdata4e_r->GetXaxis()->SetBinLabel(1,"0");
    h_unfdata4e_r->GetXaxis()->SetBinLabel(2,"1");
    h_unfdata4e_r->GetXaxis()->SetBinLabel(3,"2");
    h_unfdata4e_r->GetXaxis()->SetBinLabel(4,">2");  
    h_unfdata4e_r->GetXaxis()->SetLabelFont(42);
  }
  h_unfdata4e_r-> SetMaximum(3.); 
  h_unfdata4e_r-> SetMinimum(-1.);
  h_unfdata4e_r->Draw("E");
  TLine *line;
  
  if(var == "Mass")line = new TLine(100,1,800,1);
  else if(var == "Jets" ||var == "Jets_Central" ) line =  new TLine(0,1,4,1);
  else if(var == "Mjj"||var == "Mjj_Central") line =  new TLine(0,1,800,1);
  else if(var == "Deta"||var == "Deta_Central") line =  new TLine(0,1,4.7,1);
  else if(var =="PtJet1"||var =="PtJet2" )  line =  new TLine(30,1,500,1); 
  else if(var =="EtaJet1"||var =="EtaJet2" )  line =  new TLine(0,1,4.7,1);
  else if(var == "dRZZ")line = new TLine(0,1,6,1); 
  line->SetLineColor(kRed);
  line->Draw("SAME");
  h_unfdata4e_r->Draw("E SAME");

  lumiTextSize     = 0.4;
  cmsTextSize      = 0.48;
  extraOverCmsTextSize  = 0.80;
  CMS_lumi(pad1, iPeriod, iPos );
  pad1->Update();
  pad1->RedrawAxis();
  pad1->GetFrame()->Draw();

  string png = "~/www/VBS/"+date+"/"+var+"/"+var + "_ZZTo4l"+MCgen+tightfr+".png"; 
  string pdf ="~/www/VBS/"+date+"/"+var+"/"+var + "_ZZTo4l"+MCgen+tightfr+".pdf";
  c->SaveAs(pdf.c_str());
  c->SaveAs(png.c_str());
  }

void AllPlot4L(bool mad =1, bool tightregion =0,string date = "test")
{
  PlotUnfoldData4L("Jets",mad, tightregion, date);
  PlotUnfoldData4L("Jets_Central",mad, tightregion, date);
  PlotUnfoldData4L("Mjj",mad, tightregion, date);
  PlotUnfoldData4L("Mjj_Central",mad, tightregion, date);
  PlotUnfoldData4L("Deta",mad, tightregion, date);
  PlotUnfoldData4L("Deta_Central",mad, tightregion, date);
  PlotUnfoldData4L("PtJet1",mad, tightregion, date);
  PlotUnfoldData4L("PtJet2",mad, tightregion, date);
  PlotUnfoldData4L("EtaJet1",mad, tightregion, date);
  PlotUnfoldData4L("EtaJet2",mad, tightregion, date);
  }
void PlotUnfoldData(string var = "Jets", string fs = "4e", bool mad =1,bool tightregion =0, string date = "test"){
  
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  setTDRStyle();
  
  string finalstate;
  if(fs == "4m") finalstate = "4#mu";
  else if(fs == "2e2m") finalstate = "2e2#mu";
  else finalstate = fs;
  
  string label = finalstate + " channel";
  int iPeriod = 4; 
  int iPos = 11; 
  writeExtraText = true;    
  extraText  = "Preliminary";
  extraText2  = label.c_str();
  string tightfr;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";
  
  TLatex latex;
  latex.SetTextFont(extraTextFont);
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    
  latex.SetTextSize(0.04);    
  //latex.SetTextAlign(12); 

  TFile * unfData;
  string unfDataFileName;
  
  filePath = "../../../"; 
  system(("mkdir ~/www/VBS/"+date).c_str());
  system(("cp ~/www/VBS/index.php ~/www/VBS/"+date).c_str());
  system(("mkdir ~/www/VBS/"+date+"/"+ var).c_str()); 
  system(("cp ~/www/VBS/index.php ~/www/VBS/"+date+"/"+ var).c_str());
  
  if(mad ==1){
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Mad/").c_str());  
    unfDataFileName = filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + ".root";
    MCgen = "_Mad";
  }
  else{
      system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Pow/").c_str());        
    unfDataFileName = filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + ".root";
    MCgen = "_Pow";
  }
  dataFileName = "../" +var + "_test/DataToUnfold.root";
  
  unfData = new TFile(unfDataFileName.c_str());
  data = new TFile(dataFileName.c_str());
    
  TH1 * h_unfdata = (TH1*) unfData->Get(("ZZTo"+fs+"_" + var).c_str());
  TH1 * h_data = (TH1*) unfData->Get(("ZZTo"+fs+"_" + var+"_RECO").c_str());
  TH1 * h_reco = (TH1*) unfData->Get(("ZZTo"+fs+"_" + var+"_RECO_MC").c_str());
  TH1 * h_gen = (TH1*) unfData->Get(("ZZTo"+fs+"_" + var+"_GEN").c_str());
  
  if(var == "Mass"){
    XaxisTitle = "m_{4l} [GeV]";
    YaxisTitle = "Events/GeV"; 
  }
  else if(var == "Jets")  {
    XaxisTitle = "N jets (|#eta^{jet}| < 4.7)";
    YaxisTitle = "Events"; 
  }
  else if(var == "Jets_Central")  {
    XaxisTitle = "N jets (|#eta^{jet}| < 2.4)"; 
    YaxisTitle = "Events"; 
  }
  else if(var == "Mjj") {
    XaxisTitle ="m_{jj} (|#eta^{jet}| < 4.7) [GeV]";
    YaxisTitle = "Events/GeV"; 
  }
  else if(var == "Mjj_Central") {
    XaxisTitle ="m_{jj} (|#eta^{jet}| < 2.4) [GeV]"; 
    YaxisTitle = "Events/GeV"; 
  }
  else if( var == "Deta"){
    XaxisTitle ="#Delta#eta_{jj} (|#eta^{jet}| < 4.7)";
    YaxisTitle = "Events/(bin width)"; 
  }
  else if(var == "Deta_Central"){ 
    XaxisTitle ="#Delta#eta_{jj} (|#eta^{jet}| < 2.4)";
    YaxisTitle = "Events/(bin width)"; 
  }
  else if(var =="PtJet1") {
    XaxisTitle = "p_{T}^{jet1} [GeV]";
    YaxisTitle = "Events/GeV"; 
  }
  else if(var =="PtJet2") { 
    XaxisTitle = "p_{T}^{jet2} [GeV]"; 
    YaxisTitle = "Events/GeV"; 
  }
  else if(var =="EtaJet1"){  
    XaxisTitle = "|#eta^{jet1}|";
    YaxisTitle = "Events/(bin width)"; 
  }
  else if(var =="EtaJet2"){  
    XaxisTitle = "|#eta^{jet2}|";
    YaxisTitle = "Events/(bin width)"; 
  }
  else if(var =="dRZZ") { 
    XaxisTitle = "#DeltaR(Z_{1},Z_{2})"; 
    YaxisTitle = "Events/(bin width)"; 
  }

  string YaxisTitle2 = "Unfolded/True";

  TCanvas *c = new TCanvas ("c","c");  
  TLegend *leg = new TLegend(0.80,0.65,0.60,0.85); 
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42); 

  TPad  *pad1 = new TPad("pad1","", 0., 0.19, 1.0, 1.0);
  pad1->SetTopMargin (0.10);
  pad1->SetRightMargin (0.10);
  pad1->SetLeftMargin (0.10);
  pad1->Draw();
  
  c->cd();
  TPad  *pad2 = new TPad("pad2", "", 0., 0.0,  1.0, 0.28);
  pad2->SetTopMargin(0.10);
  pad2->SetRightMargin(0.10);
  pad2->SetLeftMargin(0.10); 
  pad2->SetBottomMargin(0.35);
  pad2->Draw(); 

  pad1->cd();
  h_gen->SetTitle("");
  h_gen->GetXaxis()->SetRange(0,25);
  h_gen->GetXaxis()->SetTitle(XaxisTitle.c_str()); 
  h_gen->GetYaxis()->SetTitle(YaxisTitle.c_str());
  h_gen->GetYaxis()->SetTitleOffset(1.4);
  h_gen->GetXaxis()->SetLabelOffset(0.5);

  h_gen->Scale(1,"width");
  h_reco->Scale(1,"width");
  h_unfdata->Scale(1,"width");
  h_data->Scale(1,"width");
  
  float max = 0;
  if(var == "Mass") max =  h_unfdata->GetBinContent(2)+h_unfdata->GetBinContent(2)*0.50+h_unfdata->GetBinError(2); 
  else if(var == "Jets"||var == "Jets_Central")  max =  h_unfdata->GetBinContent(1)+h_unfdata->GetBinContent(1)*0.50+h_unfdata->GetBinError(1); 
  else if(var == "Mjj"|| var == "Mjj_Central") max =  h_gen->GetBinContent(1)+h_gen->GetBinContent(1)*0.50+h_unfdata->GetBinError(1); 
  else if(var == "Deta"|| var == "Deta_Central" ) max =  h_gen->GetBinContent(1)+h_gen->GetBinContent(1)*0.50+h_unfdata->GetBinError(1); 
  else if(var == "PtJet1") max =  h_gen->GetBinContent(1)+h_gen->GetBinContent(1)*0.50+h_unfdata->GetBinError(1); 
  else if(var == "PtJet2") max =  h_gen->GetBinContent(1)+h_gen->GetBinContent(1)*0.50+h_unfdata->GetBinError(1); 
  else if(var == "EtaJet1") max =  h_gen->GetBinContent(1)+h_gen->GetBinContent(1)*0.80+h_unfdata->GetBinError(1);
  else if(var == "EtaJet2") max =  h_gen->GetBinContent(1)+h_gen->GetBinContent(1)*0.80+h_unfdata->GetBinError(1);
  else if(var == "dRZZ") max =  h_gen->GetBinContent(4)+h_gen->GetBinContent(4)+h_unfdata->GetBinError(1);
  else  max =  h_gen->GetBinContent(1)+h_gen->GetBinContent(1)*0.50+h_unfdata->GetBinError(1);

  h_gen->SetMaximum(max);
  h_gen->SetMinimum(0); 
  h_gen->SetLineColor(kBlue); 
  h_gen->SetLineWidth(1);
  h_gen->SetLineStyle(7);
  h_gen->Draw("HIST E");
  h_reco->SetLineColor(kRed);
  h_reco->SetLineStyle(7);
  h_reco->SetLineWidth(1);
  h_reco->Draw("HIST E SAME");
  h_unfdata->SetLineColor(kBlue);
  h_unfdata->SetLineWidth(1);
  h_unfdata->SetMarkerColor(kBlue);
  h_unfdata->SetMarkerStyle(8);
  h_unfdata->Draw("E SAME"); 
  h_data->Draw("HIST E SAME");
  h_data->SetMarkerColor(2);
  h_data->SetMarkerStyle(8);
  h_data->SetLineColor(2);
  h_data->SetLineWidth(1);
  
  leg->AddEntry(h_gen,"MC truth","l"); 
  leg->Draw("SAME");
  leg->AddEntry(h_reco,"MC reco","l"); 
  leg->Draw("SAME");
  leg->AddEntry(h_unfdata,"unfolded data","lep"); 
  leg->Draw(); 
  leg->AddEntry(h_data,"data","lep"); 
  leg->Draw("SAME");

  TH1 * h_unfdata_r = (TH1*) h_unfdata->Clone();
  
  for(int k =1;k<9;k++){
    float unf=0;
    float tr = 0;
    float ratio =0;
    float err_unf=0;
    float err_tr = 0;
    float err_ratio =0;
    unf = h_unfdata->GetBinContent(k);
    tr = h_gen->GetBinContent(k); 
    err_unf = h_unfdata->GetBinError(k);
    err_tr = h_gen->GetBinError(k);
    ratio = unf/tr;
    err_ratio = sqrt((err_tr/tr)*(err_tr/tr)+(err_unf/unf)*(err_unf/unf));
    h_unfdata_r->SetBinContent(k,ratio); 
    h_unfdata_r->SetBinError(k,err_ratio);
  }
  pad2->cd();  
  h_unfdata_r->SetTitle("");
  h_unfdata_r->SetMarkerColor(1);
  h_unfdata_r->SetLineColor(1); 
  h_unfdata_r->GetYaxis()->SetLabelSize(0.10);  
  h_unfdata_r->GetYaxis()->SetNdivisions(306); 
  h_unfdata_r->GetXaxis()->SetLabelSize(0.13); 
  h_unfdata_r->GetXaxis()->SetLabelOffset(0.05);
  h_unfdata_r->GetXaxis()->SetTitle(XaxisTitle.c_str()); 
  h_unfdata_r->GetYaxis()->SetTitle(YaxisTitle2.c_str());
  h_unfdata_r->GetXaxis()->SetTitleOffset(1.2); 
  h_unfdata_r->GetYaxis()->SetTitleOffset(0.4);
  h_unfdata_r->GetYaxis()->SetTitleSize(0.10);
  h_unfdata_r->GetXaxis()->SetTitleSize(0.13);
  h_unfdata_r-> SetMaximum(3.); 
  h_unfdata_r-> SetMinimum(-1.);
  h_unfdata_r->Draw("E");
  
  if(var == "Jets" || var == "Jets_Central"){
    h_unfdata_r->GetXaxis()->SetBinLabel(1,"0");
    h_unfdata_r->GetXaxis()->SetBinLabel(2,"1");
    h_unfdata_r->GetXaxis()->SetBinLabel(3,"2");
    h_unfdata_r->GetXaxis()->SetBinLabel(4,">2");  
    h_unfdata_r->GetXaxis()->SetLabelFont(42);
  }

  TLine *line;

  if(var == "Mass")line = new TLine(100,1,800,1);
  else if(var == "Jets" ||var == "Jets_Central" ) line =  new TLine(0,1,4,1);
  else if(var == "Mjj"||var == "Mjj_Central") line =  new TLine(0,1,800,1);
  else if(var == "Deta"||var == "Deta_Central") line =  new TLine(0,1,4.7,1);
  else if(var =="PtJet1"||var =="PtJet2" )  line =  new TLine(30,1,500,1); 
  else if(var =="EtaJet1"||var =="EtaJet2" )  line =  new TLine(0,1,4.7,1); 
  else if(var == "dRZZ")line = new TLine(0,1,6,1);

  line->SetLineColor(kRed);
  line->Draw("SAME");
  h_unfdata_r->Draw("E SAME");

  lumiTextSize     = 0.4;
  cmsTextSize      = 0.48;
  extraOverCmsTextSize  = 0.80;//0.63;
  CMS_lumi(pad1, iPeriod, iPos );
  pad1->Update();
  pad1->RedrawAxis();
  pad1->GetFrame()->Draw();
  
  string png = "~/www/VBS/"+date+"/"+var+"/"+var + "_ZZTo"+fs+MCgen+tightfr+"_binwidth.png"; 
  string pdf = "~/www/VBS/"+date+"/"+var+"/"+var + "_ZZTo"+fs+MCgen+tightfr+"_binwidth.pdf";
  c->SaveAs(pdf.c_str());
  c->SaveAs(png.c_str());
 }
void AllPlotBinWidth(bool mad =1, bool tightregion =0,string date = "test")
{
  PlotUnfoldData4L("Mass",0, tightregion, date);
  PlotUnfoldData4L("Jets",mad, tightregion, date);
  PlotUnfoldData4L("Jets_Central",mad, tightregion, date);
  PlotUnfoldData4L("Mjj",mad, tightregion, date);
  PlotUnfoldData4L("Mjj_Central",mad, tightregion, date);
  PlotUnfoldData4L("Deta",mad, tightregion, date);
  PlotUnfoldData4L("Deta_Central",mad, tightregion, date);
  PlotUnfoldData4L("PtJet1",mad, tightregion, date);
  PlotUnfoldData4L("PtJet2",mad, tightregion, date);
  PlotUnfoldData4L("EtaJet1",mad, tightregion, date);
  PlotUnfoldData4L("EtaJet2",mad, tightregion, date);
  //PlotUnfoldData4L("dRZZ",mad, tightregion, date); 
   
  PlotUnfoldData("Mass","4e",0, tightregion, date);
  PlotUnfoldData("Jets","4e",mad, tightregion, date);
  PlotUnfoldData("Jets_Central","4e",mad, tightregion, date);
  PlotUnfoldData("Mjj","4e",mad, tightregion, date);
  PlotUnfoldData("Mjj_Central","4e",mad, tightregion, date);
  PlotUnfoldData("Deta","4e",mad, tightregion, date);
  PlotUnfoldData("Deta_Central","4e",mad, tightregion, date);
  PlotUnfoldData("PtJet1","4e",mad, tightregion, date);
  PlotUnfoldData("PtJet2","4e",mad, tightregion, date);
  PlotUnfoldData("EtaJet1","4e",mad, tightregion, date);
  PlotUnfoldData("EtaJet2","4e",mad, tightregion, date);
  //PlotUnfoldData("dRZZ","4e",mad, tightregion, date);

  PlotUnfoldData("Mass","4m",0, tightregion, date);
  PlotUnfoldData("Jets","4m",mad, tightregion, date);
  PlotUnfoldData("Jets_Central","4m",mad, tightregion, date);
  PlotUnfoldData("Mjj","4m",mad, tightregion, date);
  PlotUnfoldData("Mjj_Central","4m",mad, tightregion, date);
  PlotUnfoldData("Deta","4m",mad, tightregion, date);
  PlotUnfoldData("Deta_Central","4m",mad, tightregion, date);
  PlotUnfoldData("PtJet1","4m",mad, tightregion, date);
  PlotUnfoldData("PtJet2","4m",mad, tightregion, date);
  PlotUnfoldData("EtaJet1","4m",mad, tightregion, date);
  PlotUnfoldData("EtaJet2","4m",mad, tightregion, date);
  //PlotUnfoldData("dRZZ","4m",mad, tightregion, date);

  PlotUnfoldData("Mass","2e2m",0, tightregion, date);
  PlotUnfoldData("Jets","2e2m",mad, tightregion, date);
  PlotUnfoldData("Jets_Central","2e2m",mad, tightregion, date);
  PlotUnfoldData("Mjj","2e2m",mad, tightregion, date);
  PlotUnfoldData("Mjj_Central","2e2m",mad, tightregion, date);
  PlotUnfoldData("Deta","2e2m",mad, tightregion, date);
  PlotUnfoldData("Deta_Central","2e2m",mad, tightregion, date);
  PlotUnfoldData("PtJet1","2e2m",mad, tightregion, date);
  PlotUnfoldData("PtJet2","2e2m",mad, tightregion, date);
  PlotUnfoldData("EtaJet1","2e2m",mad, tightregion, date);
  PlotUnfoldData("EtaJet2","2e2m",mad, tightregion, date); 
  //PlotUnfoldData("dRZZ","2e2m",mad, tightregion, date);
  }

void PlotUnfoldData_fs(string var = "Jets", bool mad =1,bool tightregion =0, string date = "test"){ 
  PlotUnfoldData(var,"4e",mad, tightregion, date);
  PlotUnfoldData(var,"4m",mad, tightregion, date);
  PlotUnfoldData(var,"2e2m",mad, tightregion, date);
}

void AllYouNeed_Var1(string var = "Jets", bool mad =1 ,bool tightregion =0, string date ="test")
{  
  cout << "********************** UNFOLDING DATA *************************" <<endl;
  Unfold_data_All(var,mad,tightregion,date);
  cout << "********************** PLOTTING DISTRIBUTIONS *************************" <<endl;
  PlotUnfoldData_fs(var,mad,tightregion,date);
  cout << "********************** PLOTTING 4L DISTRIBUTIONS *************************" << endl;
  PlotUnfoldData4L(var,mad,tightregion,date);
  cout << "********************** CALCULATING UNF/MC RATIO *************************" << endl;
  DoAllRatios(var,mad,tightregion);
  cout << "********************** CALCULATING SYSTEMATICS *************************" << endl;
  DoAllSystematics(var,mad,tightregion);
  cout << "********************** PLOTTING SYSTEMATICS *************************" << endl;
  AllPlots(var,mad,tightregion,date);
}

void AllYouNeed_Var2(string var = "Jets", string date ="test"){
  AllYouNeed_Var1(var,1,0,date);
  AllYouNeed_Var1(var,0,0,date);
  AllYouNeed_Var1(var,0,1,date);
  AllYouNeed_Var1(var,1,1,date);
  
}


void FakeUnfold_data(string var = "Mass", string fs = "4e", bool mad =1,bool tightregion = 0, string date = "test", bool bayes = 1, int kreg = 2){
  

#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0); 
  setTDRStyle();
  string finalstate;
  if(fs == "4m") finalstate = "4#mu";
  else if(fs == "2e2m") finalstate = "2e2#mu";
  else finalstate = fs;
  
  string label = finalstate + " channel";
  int iPeriod = 4; 
  int iPos = 11; 
  writeExtraText = true;    
  extraText  = "Simulation";
  extraText2  = label.c_str();

  string tightfr;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";
 
  filePath = "../../../";
   
  if(mad ==1){
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Mad/").c_str());  
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Mad.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + ".root";
    MCgen = "_Mad";
  }
  else{
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Pow/").c_str());  
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Pow.root";
    outputFileName = filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + ".root";
    MCgen = "_Pow";
  }
  
  system(("mkdir ~/www/VBS/"+date+"/"+ var).c_str()); 
  system(("cp ~/www/VBS/index.php ~/www/VBS/"+date+"/"+ var).c_str());

  dataFileName = "../" +var + "_test/DataToUnfoldFake.root";
  
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
  RooUnfoldSvd   unfold_svd(&response, h_measured_unf, kreg); 
  RooUnfoldBayes unfold_bayes(&response, h_measured_unf,kreg);  
 
  if(bayes ==0) {h_unfolded= (TH1*) unfold_svd.Hreco(RooUnfold::kCovariance);
  }
  else if(bayes == 1){
    h_unfolded= (TH1*) unfold_bayes.Hreco(RooUnfold::kCovariance); 
  }
  else cout << "WRONG ALGORITHM!!!" << endl;
  
  unfHistoName = "ZZTo"+ fs +"_" + var;
  recoHistoName = "ZZTo"+ fs +"_"+var+"_RECO"; 
  recoMCHistoName = "ZZTo"+ fs +"_"+var+"_RECO_MC";
  trueHistoName = "ZZTo"+ fs +"_"+var+"_GEN";
  
  if(var == "Mass") XaxisTitle = "m_{" + finalstate + "} [GeV]"; 
  else if(var == "Jets")  XaxisTitle = "N jets (|#eta^{jet}|<4.7)";
  else if(var == "Jets_Central")  XaxisTitle = "N Central Jets (|#eta^{jet}|<2.4)";
  else if(var == "Mjj")  XaxisTitle = "m_{jj} (|#eta^{jet}|<4.7) [GeV]";
  else if(var == "Deta")  XaxisTitle = "#Delta#eta_{jj}(|#eta^{jet}|<4.7) "; 
  else if(var == "Mjj_Central")  XaxisTitle = "m_{jj} (|#eta^{jet}|<2.4) [GeV]";
  else if(var == "Deta_Central")  XaxisTitle = "#Delta#eta_{jj}(|#eta^{jet}|<2.4) ";
  else if(var =="PtJet1")  XaxisTitle = "p_{T}^{jet1} [GeV]";
  else if(var =="PtJet2")  XaxisTitle = "p_{T}^{jet2} [GeV]";
  else if(var =="EtaJet1")  XaxisTitle = "|#eta^{jet1}| ";
  else if(var =="EtaJet2")  XaxisTitle = "|#eta^{jet2}| ";
  else if(var =="dRZZ")  XaxisTitle = "#DeltaR(Z_{1},Z_{2}) ";
  
  YaxisTitle = "Events";
  string YaxisTitle2 = "Unfolded/True";
 
  TCanvas *c = new TCanvas ("c","c");
  TLegend *leg = new TLegend(0.80,0.65,0.60,0.85); 
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42); 

  h_true->SetTitle("");
  h_true->GetXaxis()->SetRange(0,25);
  h_true->GetXaxis()->SetTitle(XaxisTitle.c_str()); 
  h_true->GetYaxis()->SetTitle(YaxisTitle.c_str());

  double max = h_true->GetMaximum()+ h_true->GetMaximum()*0.35;
  h_true->SetLineColor(kBlue); 
  h_true->SetLineWidth(1);
  h_true->SetLineStyle(7);
  h_true->SetMaximum(max);
  h_true->SetMinimum(0);   
  h_true->Draw("HIST E");
  if(var == "Jets" || var == "Jets_Central"){
    h_true->GetXaxis()->SetBinLabel(1,"0");
    h_true->GetXaxis()->SetBinLabel(2,"1");
    h_true->GetXaxis()->SetBinLabel(3,"2");
    h_true->GetXaxis()->SetBinLabel(4,">2");  
    h_true->GetXaxis()->SetLabelFont(42);
    h_true->GetXaxis()->SetLabelSize(0.055);
  }
  h_true->GetXaxis()->SetTitleOffset(1.3);
  h_true->GetYaxis()->SetTitleOffset(1.6);
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
  
 
  leg->AddEntry(h_true,"MC truth","l"); 
  leg->Draw("SAME");
  leg->AddEntry(h_unfolded,"unfolded data","lep"); 
  leg->Draw(); 
  leg->AddEntry(h_measured_unf,"data","lep"); 
  leg->Draw("SAME");
 
  lumiTextSize     = 0.7;
  cmsTextSize      = 0.7;
  extraOverCmsTextSize  = 0.80;
  CMS_lumi(c, iPeriod, iPos );
  
  string algo;
  if(bayes == 0){
    if(kreg ==2) algo = "_SVD_2";
    else if(kreg ==4) algo = "_SVD_4";
  }
  else if(bayes ==1 && kreg == 2) algo = "_bayes_2"; 
  else if(bayes ==1 && kreg == 4) algo = "_bayes_4";

  string png = "~/www/VBS/"+date+"/"+ var+"/"+var + "_ZZTo" + fs + MCgen+tightfr +algo+".png"; 
  string pdf = "~/www/VBS/"+date+"/"+ var+"/"+var + "_ZZTo" + fs + MCgen+tightfr +algo+ ".pdf";
  c->SaveAs(pdf.c_str());
  c->SaveAs(png.c_str());

  // output->cd();   
  // h_unfolded->Write(unfHistoName.c_str());
  // h_measured_unf->Write(recoHistoName.c_str());
  // h_measured->Write(recoMCHistoName.c_str());
  // h_true->Write(trueHistoName.c_str());
  // output->Close();
}

void FakePlot(string var = "Mass", bool mad =1, string date = "test", bool bayes = 1, int kreg = 2)
{
  FakeUnfold_data(var,"4e",mad,0,date,bayes,kreg);
  FakeUnfold_data(var,"4m",mad,0,date,bayes,kreg);
  FakeUnfold_data(var,"2e2m",mad,0,date,bayes,kreg);
  FakeUnfold_data(var,"4e",mad,1,date,bayes,kreg);
  FakeUnfold_data(var,"4m",mad,1,date,bayes,kreg);
  FakeUnfold_data(var,"2e2m",mad,1,date,bayes,kreg);

  }


void Systematics_value(string var = "Jets", string fs = "4e", string syst = "MCgen", bool mad =1,bool tightregion =0){
  
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
 
  string tightfr;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";
  
  TFile * unfSyst;
  string systFileName;
  string systHistoName_p;
  string systHistoName_m;
  TH1 *h_true_othMC;
 
  filePath = "../../../";
  if(mad ==1){
    //system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Mad/").c_str());  
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Pow.root";
    dataFileName = filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + ".root";
    systFileName =  filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + "_"+syst + ".root";
    MCgen = "_Mad";
  }
  else{
    //system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Pow/").c_str());     
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Mad.root";
    dataFileName = filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + ".root";
    systFileName =  filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + "_"+syst + ".root";
    MCgen = "_Pow";
  }
  data = new TFile(dataFileName.c_str());
  unfSyst = new TFile(systFileName.c_str());
  matrix = new TFile(matrixFileName.c_str());

  unfHistoName = "ZZTo"+ fs +"_" + var;
  recoHistoName = "ZZTo"+ fs +"_"+var+"_RECO";
  trueHistoName = "ZZTo"+ fs +"_"+var+"_GEN";
  systHistoName_p = "ZZTo"+ fs +"_" + var + "_p";
  systHistoName_m =  "ZZTo"+ fs +"_" + var + "_m";
  histoNameGen = var + "Gen_qqggJJ_ZZTo" + fs + "_st_01";
  
  h_true_othMC = (TH1*) matrix->Get(histoNameGen.c_str());
  h_true = (TH1*) data->Get(trueHistoName.c_str());
  h_measured_unf = (TH1*) data->Get(recoHistoName.c_str());  
  h_unfolded = (TH1*) data->Get(unfHistoName.c_str());
 
  float max;
  float min;
  TH1* h_max; 
  TH1* h_min; 
  float syst_percentage[9];
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

      // if(max!=0.) std::cout  << i << "-bin: min = "<< min << " max = " << max << " " << syst_percentage[i-1]  << "% " << syst.c_str() << " systematic uncertainty" << std::endl;
 
    }
  }
  else if(syst=="qqgg"||syst=="IrrBkg"||syst=="RedBkg"||syst=="JER"||syst=="JES_ModMat"||syst=="JES_ModData"||syst=="SF"||syst=="SFSq"){
    h_unfolded_p = (TH1*) unfSyst->Get(systHistoName_p.c_str());
    h_unfolded_m = (TH1*) unfSyst->Get(systHistoName_m.c_str());
    h_max = (TH1*)h_unfolded_p->Clone("");;
    h_min = (TH1*)h_unfolded_p->Clone("");;
    for(int i = 1; i<9; i++){
      
      syst_percentage[i-1]=0;
      
      Double_t vec[] ={h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i),h_unfolded->GetBinContent(i)};
      min =vec[0];
      max =0;
    
      for(int j=0;j<3;j++)
	{
	  if(vec[j]<min) min = vec[j];
	  if(vec[j]>max) max = vec[j];
	}
      
      // max = TMath::Max(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
      // min = TMath::Min(h_unfolded_p->GetBinContent(i),h_unfolded_m->GetBinContent(i));
      // std::cout << min << " " << max << std::endl;
      h_max->SetBinContent(i,max); 
      h_min->SetBinContent(i,min); 
      syst_percentage[i-1]=0.5*(max-min)*100/(h_unfolded->GetBinContent(i)); 
     
      // if(max!=0.) std::cout  << i << "-bin: min = "<< min << " max = " << max << " " << syst_percentage[i-1]  << "% " << syst.c_str() << " systematic uncertainty" << std::endl;
    }
  } 
  else std::cout << " Wrong systematic uncertainty!!!!" << std::endl;

  float M = syst_percentage[0];
  float m = syst_percentage[0];
  for(int i = 0; i<8; i++) {
      if(syst_percentage[i] > M) M = syst_percentage[i];
      if(syst_percentage[i] < m) m = syst_percentage[i];
     }

  double totUnc_min = fabs(h_unfolded->Integral()-h_min->Integral())/h_unfolded->Integral();
  double totUnc_max = fabs(h_unfolded->Integral()-h_max->Integral())/h_unfolded->Integral();
  cout << syst << ": + " <<  totUnc_max*100 << "% - " << totUnc_min*100 << "%" <<endl;
  std::cout << "Max syst = " << M << "%   Min syst = " << m << "%" << endl;

 }

void AllSystematicsValues(string var = "Jets", string fs = "4e", bool mad =1,bool tightregion =0){

  std::cout << "===========================================" << std::endl;
  std::cout << "                    " << fs.c_str() << std::endl;
  std::cout << "===========================================" << std::endl;
  Systematics_value(var.c_str(),fs.c_str(),"MCgen",mad, tightregion);
  Systematics_value(var.c_str(),fs.c_str(),"qqgg",mad, tightregion);
  Systematics_value(var.c_str(),fs.c_str(),"IrrBkg",mad, tightregion);
  Systematics_value(var.c_str(),fs.c_str(),"RedBkg",mad, tightregion);
  Systematics_value(var.c_str(),fs.c_str(),"UnfDataOverGenMC",mad,tightregion);

  if(var == "Jets"||var == "Mjj"||var == "Deta"||var == "Jets_Central"||var == "Mjj_Central"||var == "Deta_Central"||var == "PtJet1"||var == "PtJet2"||var == "EtaJet1"||var == "EtaJet2"){
    Systematics_value(var.c_str(),fs.c_str(),"JER",mad, tightregion); 
    Systematics_value(var.c_str(),fs.c_str(),"JES_ModMat",mad, tightregion); 
    Systematics_value(var.c_str(),fs.c_str(),"JES_ModData",mad, tightregion); 
  }
  Systematics_value(var.c_str(),fs.c_str(),"SFSq",mad, tightregion); 
  Systematics_value(var.c_str(),fs.c_str(),"SF",mad, tightregion);
}

void AllSystematicsValues_fs(string var = "Jets", bool mad =1,bool tightregion =0){
  std::cout << "===========================================" << std::endl;
  std::cout << "                    " << var.c_str() << std::endl;
  AllSystematicsValues(var,"4e",mad,tightregion);
  AllSystematicsValues(var,"4m",mad,tightregion);
  AllSystematicsValues(var,"2e2m",mad,tightregion);
}

void PlotResults_4l(string var = "Jets", string syst = "MCgen", bool mad =1,bool tightregion =0, string date = "test"){
  
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif 
  
  gROOT->SetStyle("Plain");  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  setTDRStyle();
  int iPeriod = 4; 
  int iPos = 11; 
  writeExtraText = true;    
  extraText  = "Preliminary";
  extraText2  = "";
  string tightfr;
  if(tightregion == 1) tightfr = "_fr";
  else tightfr = "";
  
  TFile * unfSyst;
  string systFileName;
  string systHistoName_p;
  string systHistoName_m;
  TH1 *h_true_othMC;
  system(("mkdir ~/www/VBS/"+date).c_str());
  system(("cp ~/www/VBS/index.php ~/www/VBS/"+date).c_str()); 
  system(("mkdir ~/www/VBS/"+date+"/"+ var).c_str());
  system(("cp ~/www/VBS/index.php ~/www/VBS/"+date+"/"+ var).c_str()); 
  system(("mkdir ~/www/VBS/"+date+"/"+ var+"/Systematics").c_str()); 
  system(("cp ~/www/VBS/index.php ~/www/VBS/"+date+"/"+ var+"/Systematics").c_str());
  
  filePath = "../../../";
  if(mad ==1){
    system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Mad/").c_str());  
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Pow.root";
    dataFileName = filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + ".root";
    systFileName =  filePath+"UnfoldFolder"+tightfr+"_Mad/UnfoldData_"+ var + "_"+syst + ".root";
    MCgen = "_Mad";
  }
  else{
     system(("mkdir "+filePath+ "UnfoldFolder"+tightfr+"_Pow/").c_str());     
    matrixFileName = "../" + var + "_test/matrices"+tightfr+"_Mad.root";
    dataFileName = filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + ".root";
    systFileName =  filePath+"UnfoldFolder"+tightfr+"_Pow/UnfoldData_"+ var + "_"+syst + ".root";
    MCgen = "_Pow";
  }
  data = new TFile(dataFileName.c_str());
  unfSyst = new TFile(systFileName.c_str());
  matrix = new TFile(matrixFileName.c_str());
  
  string unfHistoName_4e;
  string recoHistoName_4e;
  string trueHistoName_4e;
  string systHistoName_4e_p;
  string systHistoName_4e_m;
  string histoNameGen_4e;
  string unfHistoName_4m;
  string recoHistoName_4m;
  string trueHistoName_4m;
  string systHistoName_4m_p;
  string systHistoName_4m_m;
  string histoNameGen_4m;
  string unfHistoName_2e2m;
  string recoHistoName_2e2m;
  string trueHistoName_2e2m;
  string systHistoName_2e2m_p;
  string systHistoName_2e2m_m;
  string histoNameGen_2e2m;
  
  unfHistoName_4e = "ZZTo4e_" + var;
  recoHistoName_4e = "ZZTo4e_"+var+"_RECO";
  trueHistoName_4e = "ZZTo4e_"+var+"_GEN";
  systHistoName_4e_p = "ZZTo4e_" + var + "_p";
  systHistoName_4e_m =  "ZZTo4e_" + var + "_m";
  histoNameGen_4e = var + "Gen_qqggJJ_ZZTo4e_st_01"; 
  unfHistoName_4m = "ZZTo4m_" + var;
  recoHistoName_4m = "ZZTo4m_"+var+"_RECO";
  trueHistoName_4m = "ZZTo4m_"+var+"_GEN";
  systHistoName_4m_p = "ZZTo4m_" + var + "_p";
  systHistoName_4m_m =  "ZZTo4m_" + var + "_m";
  histoNameGen_4m = var + "Gen_qqggJJ_ZZTo4e_st_01";
  unfHistoName_2e2m = "ZZTo2e2m_" + var;
  recoHistoName_2e2m = "ZZTo2e2m_"+var+"_RECO";
  trueHistoName_2e2m = "ZZTo2e2m_"+var+"_GEN";
  systHistoName_2e2m_p = "ZZTo2e2m_" + var + "_p";
  systHistoName_2e2m_m =  "ZZTo2e2m_" + var + "_m";
  histoNameGen_2e2m = var + "Gen_qqggJJ_ZZTo2e2m_st_01";

  TH1 * h_true_othMC_4e;
  TH1 *h_true_4e;
  TH1 *h_measured_unf_4e;
  TH1 *h_unfolded_4e;
  TH1 * h_true_othMC_4m;
  TH1 *h_true_4m;
  TH1 *h_measured_unf_4m;
  TH1 *h_unfolded_4m;
  TH1 *h_true_othMC_2e2m;
  TH1 *h_true_2e2m;
  TH1 *h_measured_unf_2e2m;
  TH1 *h_unfolded_2e2m;

  TH1* h_unfolded_4e_p;
  TH1* h_unfolded_4e_m;
  TH1* h_unfolded_4m_p;
  TH1* h_unfolded_4m_m;
  TH1* h_unfolded_2e2m_p;
  TH1* h_unfolded_2e2m_m;

  h_true_othMC_4e = (TH1*) matrix->Get(histoNameGen_4e.c_str());
  h_true_4e = (TH1*) data->Get(trueHistoName_4e.c_str());
  h_measured_unf_4e = (TH1*) data->Get(recoHistoName_4e.c_str());  
  h_unfolded_4e = (TH1*) data->Get(unfHistoName_4e.c_str()); 
  h_true_othMC_4m = (TH1*) matrix->Get(histoNameGen_4m.c_str());
  h_true_4m = (TH1*) data->Get(trueHistoName_4m.c_str());
  h_measured_unf_4m = (TH1*) data->Get(recoHistoName_4m.c_str());  
  h_unfolded_4m = (TH1*) data->Get(unfHistoName_4m.c_str());
  h_true_othMC_2e2m = (TH1*) matrix->Get(histoNameGen_2e2m.c_str());
  h_true_2e2m = (TH1*) data->Get(trueHistoName_2e2m.c_str());
  h_measured_unf_2e2m = (TH1*) data->Get(recoHistoName_2e2m.c_str());  
  h_unfolded_2e2m = (TH1*) data->Get(unfHistoName_2e2m.c_str());

  float dxmax_4e[9];
  float dxmin_4e[9];
  float dxmax_4m[9];
  float dxmin_4m[9];
  float dxmax_2e2m[9];
  float dxmin_2e2m[9];
  TH1* h_max; 
  TH1* h_min; 

  std::cout << "===================================================================================" << std::endl;
  std::cout << "                                 " << syst.c_str() << std::endl;
  std::cout << "===================================================================================" << std::endl;
  if(syst == "MCgen" || syst == "UnfDataOverGenMC"){
    h_unfolded_4e_p = (TH1*) unfSyst->Get(unfHistoName_4e.c_str());
    h_unfolded_4m_p = (TH1*) unfSyst->Get(unfHistoName_4m.c_str());
    h_unfolded_2e2m_p = (TH1*) unfSyst->Get(unfHistoName_2e2m.c_str());
   
    for(int i = 1; i<9; i++){
      dxmax_4e[i-1] = 0;
      dxmin_4e[i-1] = 0;
      dxmax_4m[i-1] = 0;
      dxmin_4m[i-1] = 0; 
      dxmax_2e2m[i-1] = 0;
      dxmin_2e2m[i-1] = 0;
       
      dxmax_4e[i-1] = fabs(h_unfolded_4e_p->GetBinContent(i)- h_unfolded_4e->GetBinContent(i));
      dxmin_4e[i-1] = fabs(h_unfolded_4e_p->GetBinContent(i)- h_unfolded_4e->GetBinContent(i));
      dxmax_4m[i-1] = fabs(h_unfolded_4m_p->GetBinContent(i)- h_unfolded_4m->GetBinContent(i));
      dxmin_4m[i-1] = fabs(h_unfolded_4m_p->GetBinContent(i)- h_unfolded_4m->GetBinContent(i));
      dxmax_2e2m[i-1] = fabs(h_unfolded_2e2m_p->GetBinContent(i)- h_unfolded_2e2m->GetBinContent(i));
      dxmin_2e2m[i-1] = fabs(h_unfolded_2e2m_p->GetBinContent(i)- h_unfolded_2e2m->GetBinContent(i));
      
    }
  }
  else if(syst=="qqgg"||syst=="IrrBkg"||syst=="RedBkg"||syst=="JER"||syst=="JES_ModMat"||syst=="JES_ModData"||syst=="SF"||syst=="SFSq"){
    h_unfolded_4e_p = (TH1*) unfSyst->Get(systHistoName_4e_p.c_str());
    h_unfolded_4m_p = (TH1*) unfSyst->Get(systHistoName_4m_p.c_str());
    h_unfolded_2e2m_p = (TH1*) unfSyst->Get(systHistoName_2e2m_p.c_str());
    h_unfolded_4e_m = (TH1*) unfSyst->Get(systHistoName_4e_m.c_str());
    h_unfolded_4m_m = (TH1*) unfSyst->Get(systHistoName_4m_m.c_str());
    h_unfolded_2e2m_m = (TH1*) unfSyst->Get(systHistoName_2e2m_m.c_str());

    for(int i = 1; i<9; i++){
      dxmax_4e[i-1] = 0;
      dxmin_4e[i-1] = 0;
      dxmax_4m[i-1] = 0;
      dxmin_4m[i-1] = 0; 
      dxmax_2e2m[i-1] = 0;
      dxmin_2e2m[i-1] =0;

      Double_t vec_4e[] ={h_unfolded_4e_p->GetBinContent(i),h_unfolded_4e_m->GetBinContent(i),h_unfolded_4e->GetBinContent(i)};
      Double_t vec_4m[] ={h_unfolded_4m_p->GetBinContent(i),h_unfolded_4m_m->GetBinContent(i),h_unfolded_4m->GetBinContent(i)};
      Double_t vec_2e2m[] ={h_unfolded_2e2m_p->GetBinContent(i),h_unfolded_2e2m_m->GetBinContent(i),h_unfolded_2e2m->GetBinContent(i)};
      float max_4e = 0;
      float min_4e= vec_4e[0];
      float max_4m = 0;
      float min_4m = vec_4m[0]; 
      float max_2e2m= 0;
      float min_2e2m =vec_2e2m[0];
   
      for(int j=0;j<3;j++)
	{
	  if(vec_4e[j]<min_4e) min_4e = vec_4e[j];
	  if(vec_4e[j]>max_4e) max_4e = vec_4e[j];
	  if(vec_4m[j]<min_4m) min_4m = vec_4m[j];
	  if(vec_4m[j]>max_4m) max_4m = vec_4m[j]; 
	  if(vec_2e2m[j]<min_2e2m) min_2e2m = vec_2e2m[j];
	  if(vec_2e2m[j]>max_2e2m) max_2e2m = vec_2e2m[j];
	}
 
      dxmax_4e[i-1] = fabs(max_4e - h_unfolded_4e->GetBinContent(i));
      dxmin_4e[i-1] = fabs(min_4e - h_unfolded_4e->GetBinContent(i));
      dxmax_4m[i-1] = fabs(max_4m - h_unfolded_4m->GetBinContent(i));
      dxmin_4m[i-1] = fabs(min_4m - h_unfolded_4m->GetBinContent(i));
      dxmax_2e2m[i-1] = fabs(max_2e2m - h_unfolded_2e2m->GetBinContent(i));
      dxmin_2e2m[i-1] = fabs(min_2e2m - h_unfolded_2e2m->GetBinContent(i));
    }
  } 
  else std::cout << " Wrong systematic uncertainty!!!!" << std::endl;
 
  h_max = (TH1*)h_unfolded_4e_p->Clone("");
  h_min = (TH1*)h_unfolded_4e_p->Clone(""); 
  h_unfolded = (TH1*)h_unfolded_4e_p->Clone("");
  h_true = (TH1*)h_unfolded_4e_p->Clone(""); 
  h_true_othMC = (TH1*)h_unfolded_4e_p->Clone("");
  float unc_p;
  float unc_m;
  float c_val;
  float c_val_true; 
  float c_val_true_othMC;
  for(int i = 1; i<9; i++){
    c_val = 0;
    c_val_true =0;
    c_val_true_othMC =0;
    unc_p =0;
    unc_m = 0;
    
    c_val =  h_unfolded_4e->GetBinContent(i)+ h_unfolded_4m->GetBinContent(i)+ h_unfolded_2e2m->GetBinContent(i);
    c_val_true =  h_true_4e->GetBinContent(i)+ h_true_4m->GetBinContent(i)+ h_true_2e2m->GetBinContent(i);
    c_val_true_othMC =  h_true_othMC_4e->GetBinContent(i)+ h_true_othMC_4m->GetBinContent(i)+ h_true_othMC_2e2m->GetBinContent(i);
    unc_p = sqrt(dxmax_4e[i-1]*dxmax_4e[i-1]+dxmax_4m[i-1]* dxmax_4m[i-1]+dxmax_2e2m[i-1]*dxmax_2e2m[i-1]);
    unc_m = sqrt(dxmin_4e[i-1]*dxmin_4e[i-1]+dxmin_4m[i-1]* dxmin_4m[i-1]+dxmin_2e2m[i-1]*dxmin_2e2m[i-1]);
    
    h_max->SetBinContent(i,c_val+unc_p); 
    h_min->SetBinContent(i,c_val-unc_m);
    h_unfolded-> SetBinContent(i,c_val);
    h_true->SetBinContent(i,c_val_true);
    h_true_othMC->SetBinContent(i,c_val_true_othMC);

    // cout << "c_val = " << c_val << endl;
    // cout << " c_val+unc_p = " << c_val+ unc_p << endl;
    // cout << " c_val-unc_m = " << c_val - unc_m << endl;
    // cout << "c_true  = " << c_val_true << endl; 
    // cout << "c_true_othMC  = " << c_val_true_othMC << endl;

}
  
  TCanvas *c = new TCanvas ("c","c");
  TLegend *leg = new TLegend(0.55,0.70,0.45,0.90); 
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);

  string XaxisTitle;
  float max_y =0;
  if(var == "Mass"){
    XaxisTitle = "m_{4l} [GeV]";
    YaxisTitle = "Events/GeV";
    // max_y = h_max->GetBinContent(2)+20;//*2+100;
  }
  else if(var == "Jets"){
    XaxisTitle = "N jets (|#eta^{jet}| < 4.7)"; 
    YaxisTitle = "Events";
    //max_y = h_max->GetBinContent(1)+100;//*2+100;
  } 
  else if(var == "Jets_Central"){
    XaxisTitle = "N jets (|#eta^{jet}| < 2.4)";
    YaxisTitle = "Events";
    // max_y = h_max->GetBinContent(1)+100;//*2+100;
  } 
  else if(var == "Mjj"){
    XaxisTitle = "m_{jj} (|#eta^{jet}| < 4.7) [GeV]"; 
    YaxisTitle = "Events/GeV";
    //max_y = h_true->GetBinContent(1)+30;//*2+100;
  } 
 else if(var == "Deta"){
    XaxisTitle = "#Delta#eta_{jj} (|#eta^{jet}| < 4.7)";
    YaxisTitle = "Events/(bin width)";
    //max_y = h_true->GetBinContent(1)+30;//*2+100;
  }   
 else if(var == "Mjj_Central"){
    XaxisTitle = "m_{jj} (|#eta^{jet}| < 2.4) [GeV]"; 
    YaxisTitle = "Events/GeV";
    // max_y = h_true->GetBinContent(1)+30;//*2+100;
  } 
 else if(var == "Deta_Central"){
    XaxisTitle = "#Delta#eta_{jj} (|#eta^{jet}| < 2.4)";
    YaxisTitle = "Events/(bin width)";
    // max_y = h_true->GetBinContent(1)+30;//*2+100;
  }   
 else if(var == "PtJet1"){
    XaxisTitle = "p_{T}^{jet1} [GeV]"; 
    YaxisTitle = "Events/GeV";
    //max_y = h_true->GetBinContent(1)+30;//*2+100;
  } 
 else if(var == "PtJet2"){
    XaxisTitle = "p_{T}^{jet2} [GeV]";
    YaxisTitle = "Events/GeV";
    // max_y = h_true->GetBinContent(1)+20;//*2+100;
  } 
 else if(var == "EtaJet1"){
   XaxisTitle = "|#eta^{jet1}| ";
   YaxisTitle = "Events/(bin width)";
   //max_y = h_true->GetBinContent(1)+60;//*2+100;
  }
 else if(var == "EtaJet2"){
   XaxisTitle = "|#eta^{jet2}| ";
   YaxisTitle = "Events/(bin width)";
   //max_y = h_true->GetBinContent(1)+30;//*2+100;
  }    
  else if(var == "dRZZ"){
   XaxisTitle = "#DeltaR(Z_{1},Z_{2}) "; 
   YaxisTitle = "Events/GeV";
   //max_y = h_true->GetBinContent(1)+30;//*2+100;
  }  

  max_y = h_true->GetMaximum()*1.5; 
  h_max->Scale(1,"width");
  h_min->Scale(1,"width");
  h_unfolded->Scale(1,"width");
  h_true->Scale(1,"width");
  h_true_othMC->Scale(1,"width");
  
 
  int binmax = h_max->GetMaximumBin();
  double binwidth = h_max->GetBinWidth(binmax);  
  
  h_max->Draw("HIST");
  h_max->SetFillColor(2);
  h_max->SetFillStyle(3001);
  h_max->SetLineColor(10);
  h_max->SetTitle("");
  h_max->GetXaxis()->SetRange(0,25);
  h_max->GetXaxis()->SetTitle(XaxisTitle.c_str());  
  h_max->GetYaxis()->SetTitle(YaxisTitle.c_str()); 
  h_max->SetMaximum(max_y/binwidth);
  h_max->SetMinimum(0.);
 
  if(var == "Jets" || var == "Jets_Central"){
   h_max->GetXaxis()->SetBinLabel(1,"0");
   h_max->GetXaxis()->SetBinLabel(2,"1");
   h_max->GetXaxis()->SetBinLabel(3,"2");
   h_max->GetXaxis()->SetBinLabel(4,">2");  
   h_max->GetXaxis()->SetLabelFont(42);
   h_max->GetXaxis()->SetLabelSize(0.055);
  }
  h_max->GetYaxis()->SetTitleOffset(1.6);
  h_max->GetXaxis()->SetTitleOffset(1.3);
  h_min->Draw("HIST SAME");
  h_min->SetFillColor(10);
  h_min->SetLineColor(10);
  
  h_unfolded->Draw("E SAME"); 
  h_unfolded->SetLineColor(1); 
  h_unfolded->SetLineWidth(1);
  h_unfolded->SetMarkerStyle(20);
  h_unfolded->SetMarkerColor(1);
  h_unfolded->SetMarkerSize(1);

  h_true->SetLineColor(8);
  h_true->SetLineStyle(1);  
  h_true->SetLineWidth(1);
  h_true->Draw("HIST SAME");
  h_true_othMC->SetLineColor(9);
  h_true_othMC->SetLineStyle(1);  
  h_true_othMC->SetLineWidth(1);
  h_true_othMC->Draw("HIST SAME");
  h_unfolded->Draw("E SAME"); 
  gPad->RedrawAxis();
  
  lumiTextSize     = 0.7;
  cmsTextSize      = 0.7;
  extraOverCmsTextSize  = 0.80;
  CMS_lumi(c, iPeriod, iPos );
  
  string syst_uncertainty;
  if(syst == "MCgen") syst_uncertainty = "MC generator syst. unc.";
  else if(syst == "UnfDataOverGenMC") syst_uncertainty = "Data/MC syst. unc.";
  else if(syst=="qqgg") syst_uncertainty = "#sigma_{qq}/#sigma_{gg} syst. unc.";
  else if(syst=="IrrBkg") syst_uncertainty = " Irr. Bkg. syst. unc.";
  else if(syst=="RedBkg") syst_uncertainty = " Red. Bkg. syst. unc.";
  else if(syst=="JER") syst_uncertainty = "JER syst. unc.";
  else if(syst=="JES_ModMat") syst_uncertainty = "JES syst. unc. (modified matrix)";
  else if(syst=="JES_ModData") syst_uncertainty = "JES syst. unc. (modified data)";
  else if(syst=="SFSq") syst_uncertainty = "Scale Factor syst. unc.";
    else if(syst=="SF") syst_uncertainty = "Scale Factor syst. unc. (correlated)";
  else std::cout << "wrong systematic uncerainty!!!!" << std::endl;

  leg->AddEntry(h_unfolded,"Unfolded data","lp"); 
  leg->Draw(); 
  if(mad ==1){
    leg->AddEntry(h_true,"MC truth (MadGraph Set)","l"); 
    leg->Draw("SAME"); 
    leg->AddEntry(h_true_othMC,"MC truth (Powheg Set)","l"); 
    leg->Draw("SAME");
  }
  else{ 
    leg->AddEntry(h_true_othMC,"MC truth (MadGraph Set)","l"); 
    leg->Draw("SAME"); 
    leg->AddEntry(h_true,"MC truth (Powheg Set)","l"); 
    leg->Draw("SAME");
  }
  leg->AddEntry(h_max,syst_uncertainty.c_str(),"f"); 
  leg->Draw("SAME");

  string png = " ~/www/VBS/"+date+"/"+ var+"/Systematics/ZZTo4l_"+var+"_"+syst+MCgen+tightfr+".png";
  string pdf =" ~/www/VBS/"+date+"/"+ var+"/Systematics/ZZTo4l_"+var+"_"+syst+MCgen+tightfr+".pdf";
  c->SaveAs(png.c_str());
  c->SaveAs(pdf.c_str());

}

void AllPlots_4l(string var = "Jets", bool mad =1,bool tightregion =0,string date = "test"){
  
  PlotResults_4l(var.c_str(),"MCgen",mad, tightregion, date.c_str());
  PlotResults_4l(var.c_str(),"qqgg",mad, tightregion, date.c_str());
  PlotResults_4l(var.c_str(),"IrrBkg",mad, tightregion, date.c_str());
  PlotResults_4l(var.c_str(),"RedBkg",mad, tightregion, date.c_str());
  //  PlotResults_4l(var.c_str(),"UnfDataOverGenMC",mad,tightregion, date.c_str());
  
  if(var == "Jets"||var == "Mjj"||var == "Deta"||var == "Jets_Central"||var == "Mjj_Central"||var == "Deta_Central"||var == "PtJet1"||var == "PtJet2"||var == "EtaJet1"||var == "EtaJet2"){
    PlotResults_4l(var.c_str(),"JER",mad, tightregion, date.c_str()); 
    PlotResults_4l(var.c_str(),"JES_ModMat",mad, tightregion, date.c_str()); 
    PlotResults_4l(var.c_str(),"JES_ModData",mad, tightregion, date.c_str()); 
  }
  PlotResults_4l(var.c_str(),"SFSq",mad, tightregion, date.c_str()); 
  //  PlotResults_4l(var.c_str(),"SF",mad, tightregion, date.c_str());
}

void PlotAllSyst_4l(string var = "Jets", string date = "test"){
  AllPlots_4l(var.c_str(),1,0, date.c_str());
  AllPlots_4l(var.c_str(),1,1, date.c_str());
  AllPlots_4l(var.c_str(),0,1, date.c_str());
  AllPlots_4l(var.c_str(),0,0, date.c_str());
}

#ifndef __CINT__
int main () { AllYouNeed_Var2(); return 0; }  // Main program when run stand-alone
#endif