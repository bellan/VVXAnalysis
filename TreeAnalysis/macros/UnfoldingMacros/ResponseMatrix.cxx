#if !defined(__CINT__) || defined(__MAKECINT__)
#include "ResponseMatrix.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"
//#include "ZZAnalysis/AnalysisStep/test/Plotter/CMS_lumi.C"
//#include "ZZAnalysis/AnalysisStep/test/Plotter/tdrstyle.C"
#endif

using namespace std;

ClassImp(ResponseMatrix)

ResponseMatrix::ResponseMatrix(): TObject()
{}

ResponseMatrix::ResponseMatrix(bool weight, bool madgraph, bool tightregion): TObject() 
{  

  //  map<string, TH2 > histosGen

  if(tightregion ==1)tightfr = "_fr";
  else tightfr = "";
  
  if(weight==0){
    ggZZTo2e2mu_r  = new TFile("../../results/ZZRecoAnalyzer_SR/ggZZ2e2mu.root");
    ggZZTo4e_r     = new TFile("../../results/ZZRecoAnalyzer_SR/ggZZ4e.root");
    ggZZTo4mu_r    = new TFile("../../results/ZZRecoAnalyzer_SR/ggZZ4mu.root");
    ZZTo4lpow_r    = new TFile("../../results/ZZRecoAnalyzer_SR/ZZTo4l.root");
    ZZTo4lmad_r    = new TFile("../../results/ZZRecoAnalyzer_SR/ZZTo4lamcatnlo.root");

    ZZTo2e2muJJ_r  = new TFile("../../results/ZZRecoAnalyzer_SR/ZZTo2e2muJJ.root");
    ZZTo4eJJ_r     = new TFile("../../results/ZZRecoAnalyzer_SR/ZZTo4eJJ.root");
    ZZTo4muJJ_r    = new TFile("../../results/ZZRecoAnalyzer_SR/ZZTo4muJJ.root"); 
    
    //Truth samples (signal definition distributions) 
    ggZZTo2e2mu_g  = new TFile("../../results/ZZMCAnalyzer_MC/ggZZ2e2mu.root");
    ggZZTo4e_g     = new TFile("../../results/ZZMCAnalyzer_MC/ggZZ4e.root");
    ggZZTo4mu_g    = new TFile("../../results/ZZMCAnalyzer_MC/ggZZ4mu.root");

    ZZTo4lpow_g    = new TFile("../../results/ZZMCAnalyzer_MC/ZZTo4l.root");
    ZZTo4lmad_g    = new TFile("../../results/ZZMCAnalyzer_MC/ZZTo4lamcatnlo.root");

    ZZTo2e2muJJ_g  = new TFile("../../results/ZZMCAnalyzer_MC/ZZTo2e2muJJ.root");
    ZZTo4eJJ_g     = new TFile("../../results/ZZMCAnalyzer_MC/ZZTo4eJJ.root");
    ZZTo4muJJ_g    = new TFile("../../results/ZZMCAnalyzer_MC/ZZTo4muJJ.root"); 

    
    //Files for Powheg and MGatNLO theoretical uncertainties  //To be fix for 13 TeV
    //ZZMCsystNamePow = "../../PowhegSystVar"+tightfr+".root ";
    //ZZMCsystNameMGatNLO = "../../MGatNLOSystVar"+tightfr+".root ";
    //ZZMCsystPow_g = new TFile(ZZMCsystNamePow.c_str());
    //ZZMCsystMGatNLO_g = new TFile(ZZMCsystNameMGatNLO.c_str());
    
    //rescue file
    rescue = new TFile("../../results/ZZRecoAnalyzer_SR/ZZTo4l.root");
    
    //output
    if(madgraph ==1)  {
      fileName    = "matrices" + tightfr+ "_Mad.root";
      fileName_SF = "matrices"+tightfr+ "_SF_Mad.root";
      fileName_JE = "matrices"+tightfr+ "_JESJER_Mad.root";
      mc = "Mad";
    } 
    else {
      fileName    = "matrices" + tightfr+ "_Pow.root"; 
      fileName_SF = "matrices" + tightfr+ "_SF_Pow.root";
      fileName_JE = "matrices" + tightfr+ "_JESJER_Pow.root";
      mc = "Pow";
    }
    W = "";
  }
  
  else{
    
    //Reco samples (response matrices and signal region distributions) 

    ZZTo4lpow_r   = new TFile("../../results/ZZRecoWAnalyzer_SR/ZZTo4l.root");
    ZZTo4lmad_r   = new TFile("../../results/ZZRecoWAnalyzer_SR/ZZTo4lamcatnlo.root");

    ggZZTo2e2mu_r = new TFile("../../results/ZZRecoWAnalyzer_SR/ggZZ2e2mu.root");
    ggZZTo4e_r    = new TFile("../../results/ZZRecoWAnalyzer_SR/ggZZ4e.root");
    ggZZTo4mu_r   = new TFile("../../results/ZZRecoWAnalyzer_SR/ggZZ4mu.root");

    ZZTo2e2muJJ_r = new TFile("../../results/ZZRecoWAnalyzer_SR/ZZTo2e2muJJ.root");
    ZZTo4eJJ_r    = new TFile("../../results/ZZRecoWAnalyzer_SR/ZZTo4eJJ.root");
    ZZTo4muJJ_r   = new TFile("../../results/ZZRecoWAnalyzer_SR/ZZTo4muJJ.root");
    
    //Truth samples (signal definition distributions) 

    ZZTo4lmad_g   = new TFile("../../results/ZZMCWAnalyzer_MC/ZZTo4lamcatnlo.root");
    ZZTo4lpow_g   = new TFile("../../results/ZZMCWAnalyzer_MC/ZZTo4l.root");

    ggZZTo2e2mu_g = new TFile("../../results/ZZMCWAnalyzer_MC/ggZZ2e2mu.root");
    ggZZTo4e_g    = new TFile("../../results/ZZMCWAnalyzer_MC/ggZZ4e.root");
    ggZZTo4mu_g   = new TFile("../../results/ZZMCWAnalyzer_MC/ggZZ4mu.root");

    ZZTo2e2muJJ_g = new TFile("../../results/ZZMCWAnalyzer_MC/ZZTo2e2muJJ.root");
    ZZTo4eJJ_g    = new TFile("../../results/ZZMCWAnalyzer_MC/ZZTo4eJJ.root");
    ZZTo4muJJ_g   = new TFile("../../results/ZZMCWAnalyzer_MC/ZZTo4muJJ.root");

   
    //rescue file
    rescue = new TFile("../../results/ZZRecoWAnalyzer_SR/ZZTo4l.root");
    
    //output 
    if(madgraph ==1)  {
      fileName = "weightedMatrices" + tightfr+ "_Mad.root";
      mc = "Mad";
    } 
    else {
      fileName ="weightedMatrices" + tightfr+ "_Pow.root"; 
      mc = "Pow";
    }
     W = "W_";
  }

  // ggZZTo2e2mu_r->Close();
  // ggZZTo4e_r->Close();
  // ggZZTo4mu_r->Close();  
  // ZZTo4lpow_r->Close();  
  // ZZTo4lmad_r->Close();
  // ZZTo2e2muJJ_r->Close();
  // ZZTo4eJJ_r->Close();
  // ZZTo4muJJ_r->Close();  
  // ggZZTo2e2mu_g->Close();
  // ggZZTo4e_g->Close();
  // ggZZTo4mu_g->Close();
  // ZZTo4lpow_g->Close();
  // ZZTo4lmad_g->Close();
  // ZZTo2e2muJJ_g->Close();
  // ZZTo4eJJ_g->Close();
  // ZZTo4muJJ_g->Close();
}


ResponseMatrix::~ResponseMatrix(){}

//Build the standard response matrix, reco and gen distributions;
//Buind response matrices and distributions varying the gg->ZZ and qq->ZZ cross sections by their uncertainties
void ResponseMatrix::Build(string var, string dataset, string finalstate, int xs_qq, int xs_gg, bool mad)
{ 
  output = new TFile((var+"_test/"+fileName).c_str(), "UPDATE");

  int b = 0;  

  if(var == "Mass") {
    variable = var;
    b=9;
  } 
  else if(var == "dRZZ") {
    variable = var;
    b=7;
  }
  else if(var == "Jets" ||var == "Jets_Central"){
    variable = var + "_JERSmear";
    b=6;
  }
  else if(var == "Mjj" ||var == "Mjj_Central" || var == "Deta" || var == "Deta_Central"){
    variable = var+"_JERSmear";
    b=3;
  }
  else if(var == "PtJet1" ||var == "PtJet2"||var == "EtaJet1" ||var == "EtaJet2" ){
    variable = var + "_JERSmear";
    b=6;
  }
 
  safeMatrixName  = "ResMat_ZZTo2e2m_" + variable +"_"+ W + "01"+ tightfr;
  safeHistoName   = "ZZTo2e2m_" + variable +"_"+ W +"01"; 
  matrixName      = "ResMat_ZZTo" + finalstate + "_" + variable+"_"+ W + dataset+ tightfr;
  histoName_reco  = "ZZTo" + finalstate + "_" + variable+"_"+ W +dataset;
  histoName_gen   =  "ZZTo" + finalstate + "_" + var + "Gen_" + W  +dataset+ tightfr;


  h_Resmat_safe_tmp = (TH2*) rescue->Get(safeMatrixName.c_str()); 
  h_safe_tmp        = (TH1*) rescue->Get(safeHistoName.c_str());  
  h_Resmat_safe     = (TH2*) h_Resmat_safe_tmp->Clone(safeMatrixName.c_str()); 
  h_safe            = (TH1*) h_safe_tmp->Clone(safeHistoName.c_str()); 
 

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
  
  h_Resmat_gg4mu   = (TH2*) ggZZTo4mu_r->Get(matrixName.c_str()); 
  h_Resmat_gg4e    = (TH2*) ggZZTo4e_r->Get(matrixName.c_str()); 
  h_Resmat_gg2e2mu = (TH2*) ggZZTo2e2mu_r->Get(matrixName.c_str()); 

  h_Resmat_4muJJ   = (TH2*) ZZTo4muJJ_r->Get(matrixName.c_str()); 
  h_Resmat_4eJJ    = (TH2*) ZZTo4eJJ_r->Get(matrixName.c_str()); 
  h_Resmat_2e2muJJ = (TH2*) ZZTo2e2muJJ_r->Get(matrixName.c_str()); 

  h_Resmat_4lmad   = (TH2*) ZZTo4lmad_r->Get(matrixName.c_str());
  h_Resmat_4lpow   = (TH2*) ZZTo4lpow_r->Get(matrixName.c_str());

  h_gg4mu          = (TH1*) ggZZTo4mu_r->Get(histoName_reco.c_str());
  h_gg4e           = (TH1*) ggZZTo4e_r->Get(histoName_reco.c_str()); 
  h_gg2e2mu        = (TH1*) ggZZTo2e2mu_r->Get(histoName_reco.c_str()); 

  h_4muJJ          = (TH1*) ZZTo4muJJ_r->Get(histoName_reco.c_str()); 
  h_4eJJ           = (TH1*) ZZTo4eJJ_r->Get(histoName_reco.c_str()); 
  h_2e2muJJ        = (TH1*) ZZTo2e2muJJ_r->Get(histoName_reco.c_str()); 

  h_4lpow          = (TH1*) ZZTo4lpow_r->Get(histoName_reco.c_str()); 
  h_4lmad          = (TH1*) ZZTo4lmad_r->Get(histoName_reco.c_str()); 

  h_gg4mu_gen      = (TH1*) ggZZTo4mu_g->Get(histoName_gen.c_str()); 
  h_gg4e_gen       = (TH1*) ggZZTo4e_g->Get(histoName_gen.c_str()); 
  h_gg2e2mu_gen    = (TH1*) ggZZTo2e2mu_g->Get(histoName_gen.c_str()); 

  h_4muJJ_gen      = (TH1*) ZZTo4muJJ_g->Get(histoName_gen.c_str()); 
  h_4eJJ_gen       = (TH1*) ZZTo4eJJ_g->Get(histoName_gen.c_str()); 
  h_2e2muJJ_gen    = (TH1*) ZZTo2e2muJJ_g->Get(histoName_gen.c_str());

  h_4lpow_gen      = (TH1*) ZZTo4lpow_g->Get(histoName_gen.c_str());
  h_4lmad_gen      = (TH1*) ZZTo4lmad_g->Get(histoName_gen.c_str());


  if(h_Resmat_gg4mu   == NULL)  h_Resmat_gg4mu   = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_gg4e    == NULL)  h_Resmat_gg4e    = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_gg2e2mu == NULL)  h_Resmat_gg2e2mu = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe"); 

  if(h_Resmat_4muJJ   == NULL)  h_Resmat_4muJJ   = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_4eJJ    == NULL)  h_Resmat_4eJJ    = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_2e2muJJ == NULL)  h_Resmat_2e2muJJ = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");

  if(h_gg4mu          == NULL)  h_gg4mu          = (TH1*) h_safe->Clone("h_safe");
  if(h_gg4e           == NULL)  h_gg4e           = (TH1*) h_safe->Clone("h_safe");
  if(h_gg2e2mu        == NULL)  h_gg2e2mu        = (TH1*) h_safe->Clone("h_safe"); 

  if(h_4muJJ          == NULL)  h_4muJJ          = (TH1*) h_safe->Clone("h_safe");
  if(h_4eJJ           == NULL)  h_4eJJ           = (TH1*) h_safe->Clone("h_safe");
  if(h_2e2muJJ        == NULL)  h_2e2muJJ        = (TH1*) h_safe->Clone("h_safe");

  if(h_gg4mu_gen      == NULL)  h_gg4mu_gen      = (TH1*) h_safe->Clone("h_safe");
  if(h_gg4e_gen       == NULL)  h_gg4e_gen       = (TH1*) h_safe->Clone("h_safe");
  if(h_gg2e2mu_gen    == NULL)  h_gg2e2mu_gen    = (TH1*) h_safe->Clone("h_safe");

  if(h_4muJJ_gen      == NULL)  h_4muJJ_gen      = (TH1*) h_safe->Clone("h_safe");
  if(h_4eJJ_gen       == NULL)  h_4eJJ_gen       = (TH1*) h_safe->Clone("h_safe");
  if(h_2e2muJJ_gen    == NULL)  h_2e2muJJ_gen    = (TH1*) h_safe->Clone("h_safe");

  if(h_Resmat_4lpow   == NULL)  h_Resmat_4lpow   = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_4lmad   == NULL)  h_Resmat_4lmad   = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");


  if(h_4lpow          == NULL)  h_4lpow          = (TH2*) h_safe->Clone("h_safe");
  if(h_4lpow_gen      == NULL)  h_4lpow          = (TH2*) h_safe->Clone("h_safe");
  
  //cout << finalstate << " qq4e " <<  h_qq4e->Integral() <<  " qq4m " <<  h_4lpow->Integral() << " qq2e2m " <<  h_qq2e2mu->Integral() << endl;
  //cout << finalstate << " gg4e " <<  h_gg4e->Integral() <<  " gg4m " <<  h_gg4mu->Integral() << " gg2e2m " <<  h_gg2e2mu->Integral() << endl;
  
 TH2 * h_Resmat_gg4mu_cl    = (TH2*) h_Resmat_gg4mu->Clone("h_Resmat_4lpow");
 TH1 * h_gg4mu_cl           = (TH1*) h_gg4mu->Clone("h_4l");
 TH1 * h_gg4mu_gen_cl       = (TH1*) h_gg4mu_gen->Clone("h_4l");
 TH2 * h_Resmat_4muJJ_cl    = (TH2*) h_Resmat_4muJJ->Clone("h_Resmat_4lpow");
 TH1 * h_4muJJ_cl           = (TH1*) h_4muJJ->Clone("h_4l");
 TH1 * h_4muJJ_gen_cl       = (TH1*) h_4muJJ_gen->Clone("h_4l");
 

  h_gg4mu_cl->Add(h_gg4e);
  h_gg4mu_cl->Add(h_gg2e2mu);
  h_gg4mu_gen_cl->Add(h_gg4e_gen);
  h_gg4mu_gen_cl->Add(h_gg2e2mu_gen);
  h_Resmat_gg4mu_cl-> Add(h_Resmat_gg4e, 1);
  h_Resmat_gg4mu_cl-> Add(h_Resmat_gg2e2mu, 1);

  h_4muJJ_cl->Add(h_4eJJ);
  h_4muJJ_cl->Add(h_2e2muJJ);
  h_4muJJ_gen_cl->Add(h_4eJJ_gen);
  h_4muJJ_gen_cl->Add(h_2e2muJJ_gen);
  h_Resmat_4muJJ_cl-> Add(h_Resmat_4eJJ, 1);
  h_Resmat_4muJJ_cl-> Add(h_Resmat_2e2muJJ, 1);
  
  //cout << finalstate << " 4l " <<  h_4lpow_cl->Integral() << " gg " << h_gg4mu_cl->Integral() << " qqJJ " << h_4muJJ_cl->Integral() << endl;
  if(mad == 1){
    h_Resmat_4lTot = (TH2*) h_Resmat_4lmad->Clone("h_Resmat_4l");
    h_4lTot        = (TH1*) h_4lmad->Clone("h_4l");
    h_4lTot_gen    = (TH1*) h_4lmad_gen->Clone("h_4l");
  }
  else{
    h_Resmat_4lTot = (TH2*)h_Resmat_4lpow->Clone("h_Resmat_4l");
    h_4lTot = (TH1*) h_4lpow->Clone("h_4l");
    h_4lTot_gen = (TH1*) h_4lpow_gen->Clone("h_4l");
  }


  h_Resmat_ggTot = (TH2*)h_Resmat_gg4mu_cl->Clone("h_Resmat_gg"); 
  h_Resmat_JJTot = (TH2*)h_Resmat_4muJJ_cl->Clone("h_Resmat_JJ"); 
  
  // cout <<"gg= "<< h_Resmat_gg4mu->Integral(0,1,0,50)<<" " << h_gg4mu->Integral(0,1)<<endl;  
  // cout <<"JJ= " <<h_Resmat_4muJJ->Integral(0,1,0,50)<<" " <<h_4muJJ->Integral(0,1)<<endl;
  // cout <<"4l= "<< h_Resmat_4lTot->Integral(0,1,0,50)<<" " <<h_4l->Integral(0,1)<<endl;
  
  unc_qq = 0.0444; //pdf: 3.4%  scale: 2.85%  //FIXME
  unc_gg = 0.2536; //pdf: 7.10%  scale: 24.35%  
 
  if(xs_qq == 0) { 
    h_Resmat_4lTot->Scale(1);
    h_4lTot->Scale(1);
    h_4lTot_gen->Scale(1); 
  } 
  else if(xs_qq == 1) {
    h_Resmat_4lTot->Scale(1+unc_qq); 
    h_4lTot->Scale(1+unc_qq); 
    h_4lTot_gen->Scale(1+unc_qq);
  }
  else if(xs_qq == -1){
    h_Resmat_4lTot->Scale(1-unc_qq); 
    h_4lTot->Scale(1-unc_qq); 
    h_4lTot_gen->Scale(1-unc_qq);
  }
  else std::cout << "Error: xs_qq must be -1, 0 or 1" << std::endl;

  if(xs_gg == 0) {
    h_Resmat_ggTot->Scale(1);
    h_gg4mu_cl->Scale(1);
    h_gg4mu_gen_cl->Scale(1);
  }
  else if(xs_gg == 1) {
    h_Resmat_ggTot->Scale(1+unc_gg); 
    h_gg4mu_cl->Scale(1+unc_gg); 
    h_gg4mu_gen_cl->Scale(1+unc_gg);
  }
  else if(xs_gg == -1) {
    h_Resmat_ggTot->Scale(1-unc_gg); 
    h_gg4mu_cl->Scale(1-unc_gg); 
    h_gg4mu_gen_cl->Scale(1-unc_gg);
  }
  else std::cout << "Error: xs_gg must be -1, 0 or 1" << std::endl;

  h_Resmat = (TH2*)h_Resmat_4lTot->Clone("h_Resmat_4lTot");
  //To comment only if you do NOT want to use MCFM and Phanton!!
  h_Resmat->Add(h_Resmat_ggTot); 
  h_Resmat->Add(h_Resmat_JJTot);
  
  h_4lTot_c = (TH1*) h_4lTot ->Clone("h_4lTot"); 
    //To comment only if you do NOT want to use MCFM and Phanton!!
  h_4lTot_c->Add(h_gg4mu_cl); 
  h_4lTot_c->Add(h_4muJJ_cl);
 
  //cout << " tot integral " << h_4lTot_c->Integral() << endl;

  h_4lTot_gen_c = (TH1*) h_4lTot_gen ->Clone("h_4lTot_gen");
  //To comment only if you do NOT want to use MCFM and Phanton!!
  h_4lTot_gen_c->Add(h_gg4mu_gen_cl);
  h_4lTot_gen_c->Add(h_4muJJ_gen_cl);
  
  totalint = h_Resmat->Integral();
  //cout << h_Resmat->Integral(0,1,0,50)<<endl;

  h_Resmat_normTot = (TH2*)h_Resmat->Clone("h_Resmat");
  h_Resmat_normTot->Scale(1/totalint);
  
  //std::cout << "total integral " << totalint << std::endl;

  string unc;
  if(xs_qq == 0){
    if(xs_gg == 0) unc = "_st_"; //standard, no variations  
    else if(xs_gg == 1) unc = "_ggp_";
    else if(xs_gg == -1) unc = "_ggm_";
  }
  else if(xs_qq == 1){
    if(xs_gg == 0) unc = "_qqp_"; 
    else if(xs_gg == 1) unc = "_qqp_ggp_";
    else if(xs_gg == -1) unc = "_qqp_ggm_";
  }
  else if(xs_qq == -1){
    if(xs_gg == 0) unc = "_qqm_"; 
    else if(xs_gg == 1) unc = "_qqm_ggp_";
    else if(xs_gg == -1) unc = "_qqm_ggm_";
  }

  string matrixNameFile         = "ResMat_qqggJJ_"+var+"_ZZTo" + finalstate + unc + dataset; 
  string matrixNormTotNameFile  = "ResMat_qqggJJ_"+var+"_normTot_ZZTo" + finalstate + unc + dataset; 
  string histoName_recoFile     = var+"_qqggJJ_ZZTo" + finalstate + unc + dataset; 
  string histoName_genFile      = var+"Gen_qqggJJ_ZZTo" + finalstate + unc + dataset; 
  string histoName_recoFile_err = var+"_statErr_qqggJJ_ZZTo" + finalstate + unc + dataset; 

  h_Resmat->SetTitle(matrixNameFile.c_str());
  h_Resmat_normTot->SetTitle(matrixNormTotNameFile.c_str());
  h_4lTot_c->SetTitle(histoName_recoFile.c_str());
  h_4lTot_gen_c->SetTitle(histoName_genFile.c_str());
  
  output->cd(); 
  h_Resmat->Write(matrixNameFile.c_str());
  h_Resmat_normTot->Write(matrixNormTotNameFile.c_str());
  h_4lTot_c->Write(histoName_recoFile.c_str());
  h_4lTot_gen_c->Write(histoName_genFile.c_str());
  output->Close();
  
}

//Build response matrices and distributions varying lepton scale factors by their uncertainties
void ResponseMatrix::Build_SF(string var, string dataset, string finalstate, string unc, bool mad)
{

  output = new TFile((var+"_test/"+fileName_SF).c_str(), "UPDATE");
   
  int b = 0;  

  if(var == "Jets" ||var == "Jets_Central")b =6;
  else if(var == "PtJet1" ||var =="PtJet2"||var == "EtaJet1" ||var == "EtaJet2") b =6;   
  else if(var == "Mass") b = 9;
  else if(var == "dRZZ") b = 7;  
  else  b=3;
 
  string observable;
  if(var == "Mass" || var == "dRZZ") observable = var;
  else observable = var + "_JERSmear"; 

  safeMatrixName = "ResMat_ZZTo2e2m_" + observable + "_01";
  safeHistoName = "ZZTo2e2m_" + observable + "_01";
  matrixName = "ResMat_ZZTo" + finalstate + "_"+var+"_"+unc+"_"+ dataset + tightfr;
  histoName_reco = "ZZTo" + finalstate + "_"+var+"_"+unc +"_"+ dataset;
  histoName_gen =  "ZZTo" + finalstate + "_"+var+"Gen_"+ dataset + tightfr;

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
 
  float totalint = 0; 
  h_Resmat_gg4mu   = (TH2*) ggZZTo4mu_r->Get(matrixName.c_str()); 
  h_Resmat_gg4e    = (TH2*) ggZZTo4e_r->Get(matrixName.c_str()); 
  h_Resmat_gg2e2mu = (TH2*) ggZZTo2e2mu_r->Get(matrixName.c_str()); 
  h_Resmat_4muJJ   = (TH2*) ZZTo4muJJ_r->Get(matrixName.c_str()); 
  h_Resmat_4eJJ    = (TH2*) ZZTo4eJJ_r->Get(matrixName.c_str()); 
  h_Resmat_2e2muJJ = (TH2*) ZZTo2e2muJJ_r->Get(matrixName.c_str()); 
  h_Resmat_4lmad   = (TH2*) ZZTo4lmad_r->Get(matrixName.c_str());
  h_Resmat_4lpow   = (TH2*) ZZTo4lpow_r->Get(matrixName.c_str()); 

  h_gg4mu          = (TH1*) ggZZTo4mu_r->Get(histoName_reco.c_str()); 
  h_gg4e           = (TH1*) ggZZTo4e_r->Get(histoName_reco.c_str()); 
  h_gg2e2mu        = (TH1*) ggZZTo2e2mu_r->Get(histoName_reco.c_str()); 
  h_4muJJ          = (TH1*) ZZTo4muJJ_r->Get(histoName_reco.c_str()); 
  h_4eJJ           = (TH1*) ZZTo4eJJ_r->Get(histoName_reco.c_str()); 
  h_2e2muJJ        = (TH1*) ZZTo2e2muJJ_r->Get(histoName_reco.c_str()); 
  h_4lmad          = (TH1*) ZZTo4lmad_r->Get(histoName_reco.c_str()); 
  h_4lpow          = (TH1*) ZZTo4lpow_r->Get(histoName_reco.c_str()); 
  
  h_gg4mu_gen      = (TH1*) ggZZTo4mu_g->Get(histoName_gen.c_str()); 
  h_gg4e_gen       = (TH1*) ggZZTo4e_g->Get(histoName_gen.c_str()); 
  h_gg2e2mu_gen    = (TH1*) ggZZTo2e2mu_g->Get(histoName_gen.c_str()); 
  h_4muJJ_gen      = (TH1*) ZZTo4muJJ_g->Get(histoName_gen.c_str()); 
  h_4eJJ_gen       = (TH1*) ZZTo4eJJ_g->Get(histoName_gen.c_str()); 
  h_2e2muJJ_gen    = (TH1*) ZZTo2e2muJJ_g->Get(histoName_gen.c_str()); 
  h_4lmad_gen      = (TH1*) ZZTo4lmad_g->Get(histoName_gen.c_str());
  h_4lpow_gen      = (TH1*) ZZTo4lpow_g->Get(histoName_gen.c_str()); 
  
  if(h_Resmat_gg4mu == NULL)   h_Resmat_gg4mu = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_gg4e == NULL)    h_Resmat_gg4e = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_gg2e2mu == NULL) h_Resmat_gg2e2mu = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe"); 
  if(h_Resmat_4muJJ == NULL)   h_Resmat_4muJJ = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_4eJJ == NULL)    h_Resmat_4eJJ = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_2e2muJJ == NULL) h_Resmat_2e2muJJ = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_gg4mu == NULL)          h_gg4mu = (TH1*) h_safe->Clone("h_safe");
  if(h_gg4e == NULL)           h_gg4e = (TH1*) h_safe->Clone("h_safe");
  if(h_gg2e2mu == NULL)        h_gg2e2mu = (TH1*) h_safe->Clone("h_safe"); 
  if(h_4muJJ == NULL)          h_4muJJ = (TH1*) h_safe->Clone("h_safe");
  if(h_4eJJ == NULL)           h_4eJJ = (TH1*) h_safe->Clone("h_safe");
  if(h_2e2muJJ == NULL)        h_2e2muJJ = (TH1*) h_safe->Clone("h_safe");
  if(h_gg4mu_gen == NULL)      h_gg4mu_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_gg4e_gen == NULL)       h_gg4e_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_gg2e2mu_gen == NULL)    h_gg2e2mu_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_4muJJ_gen == NULL)      h_4muJJ_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_4eJJ_gen == NULL)       h_4eJJ_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_2e2muJJ_gen == NULL)    h_2e2muJJ_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_Resmat_4lpow == NULL)   h_Resmat_4lpow = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_4lpow == NULL)          h_4lpow = (TH2*) h_safe->Clone("h_safe");
  if(h_4lpow_gen == NULL)      h_4lpow_gen = (TH2*) h_safe->Clone("h_safe");
 
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
   

  if(mad ==1){

    h_Resmat_4lTot = (TH2*) h_Resmat_4lmad->Clone("h_Resmat_4lmad");
    h_4lTot        = (TH1*) h_4lmad->Clone("h_4lmad");
    h_4lTot_gen    = (TH1*) h_4lmad_gen->Clone("h_4lmad");
  }
  else{
    h_Resmat_4lTot = (TH2*)h_Resmat_4lpow->Clone("h_Resmat_4lpow");
    h_4lTot        = (TH1*) h_4lpow->Clone("h_4lpow");
    h_4lTot_gen    = (TH1*) h_4lpow_gen->Clone("h_4lpow");
  }

  h_Resmat_ggTot = (TH2*)h_Resmat_gg4mu->Clone("h_Resmat_gg"); 
  h_Resmat_JJTot = (TH2*)h_Resmat_4muJJ->Clone("h_Resmat_JJ"); 
  
  //cout <<"gg= "<< h_Resmat_gg4mu->Integral(0,1,0,50)<<" " << h_gg4mu->Integral(0,1)<<endl;  
  //cout <<"JJ= " <<h_Resmat_4muJJ->Integral(0,1,0,50)<<" " <<h_4muJJ->Integral(0,1)<<endl;
  //cout <<"4l= "<< h_Resmat_4lTot->Integral(0,1,0,50)<<" " <<h_4l->Integral(0,1)<<endl;
  
  h_Resmat = (TH2*)h_Resmat_4lTot->Clone("h_Resmat_4lTot");
  h_Resmat->Add(h_Resmat_ggTot); 
  h_Resmat->Add(h_Resmat_JJTot);
 
  h_4lTot_c = (TH1*) h_4lTot ->Clone("h_4lTot"); 
  h_4lTot_c->Add(h_gg4mu); 
  h_4lTot_c->Add(h_4muJJ);
 
  h_4lTot_gen_c = (TH1*) h_4lTot_gen ->Clone("h_4lTot_gen");
  h_4lTot_gen_c->Add(h_gg4mu_gen);
  h_4lTot_gen_c->Add(h_4muJJ_gen);
 
  totalint = h_Resmat->Integral();
 
  h_Resmat_normTot = (TH2*)h_Resmat->Clone("h_Resmat");
   
  h_Resmat_normTot->Scale(1/totalint);
  
  //std::cout << "total integral " << totalint << std::endl;
  
  string matrixNameFile = "ResMat_qqggJJ_"+var+"_ZZTo" + finalstate + "_" + unc + "_" + dataset; 
  string matrixNormTotNameFile =  "ResMat_qqggJJ_"+var+"_normTot_ZZTo" + finalstate + "_" + unc + "_" + dataset; 
  string histoName_recoFile = var+"_qqggJJ_ZZTo" + finalstate + "_" + unc + "_" + dataset; 
  string histoName_genFile =  var+"Gen_qqggJJ_ZZTo" + finalstate + "_" + unc + "_" + dataset; 
  string histoName_recoFile_err =  var+"_statErr_qqggJJ_ZZTo" + finalstate + "_" + unc + "_" + dataset; 

  h_Resmat->SetTitle(matrixNameFile.c_str());
  h_Resmat_normTot->SetTitle(matrixNormTotNameFile.c_str());
  h_4lTot_c->SetTitle(histoName_recoFile.c_str());
  h_4lTot_gen_c->SetTitle(histoName_genFile.c_str());
 
  TH1 * h_4lTot_err = (TH1*) h_4lTot_c ->Clone("h_4lTot_c");
 
  float err_stat = 0;
  float bin = 0;
  for(int i =1; i <b; i++){
    err_stat = 0;
    bin = 0;
    bin = h_4lTot->GetBinContent(i);
    err_stat = sqrt(bin);
    h_4lTot_err->SetBinError(i,err_stat);
   
  }
  
  output->cd(); 
  h_Resmat->Write(matrixNameFile.c_str());
  h_Resmat_normTot->Write(matrixNormTotNameFile.c_str());
  h_4lTot_c->Write(histoName_recoFile.c_str());
  h_4lTot_gen_c->Write(histoName_genFile.c_str());
  h_4lTot_err->Write(histoName_recoFile_err.c_str());
  output->Close();
  
}

//Build response matrices and distributions varying JER and JES by their uncertainties
void ResponseMatrix::Build_JE(string var, string dataset, string finalstate, string unc, bool mad)
{
  output = new TFile((var+"_test/"+fileName_JE).c_str(), "UPDATE");

  int b = 0;  
  
  if(var == "Jets" ||var == "Jets_Central")b =6;
  else if(var == "PtJet1" ||var =="PtJet2"||var == "EtaJet1" ||var == "EtaJet2") b =6;     
  else  b=3;


  safeMatrixName = "ResMat_ZZTo2e2m_" + var +"_JERSmear_01";
  safeHistoName = "ZZTo2e2m_" + var +"_JERSmear_01"; 
  matrixName = "ResMat_ZZTo" + finalstate + "_"+var+"_"+unc+"Smear_" + dataset + tightfr;
  histoName_reco = "ZZTo" + finalstate + "_"+var+"_"+unc+"Smear_" + dataset;
  histoName_gen =  "ZZTo" + finalstate + "_"+var+"Gen_"+ dataset + tightfr;

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
 
  float totalint = 0; 
  h_Resmat_gg4mu   = (TH2*) ggZZTo4mu_r->Get(matrixName.c_str()); 
  h_Resmat_gg4e    = (TH2*) ggZZTo4e_r->Get(matrixName.c_str()); 
  h_Resmat_gg2e2mu = (TH2*) ggZZTo2e2mu_r->Get(matrixName.c_str()); 
  h_Resmat_4muJJ   = (TH2*) ZZTo4muJJ_r->Get(matrixName.c_str()); 
  h_Resmat_4eJJ    = (TH2*) ZZTo4eJJ_r->Get(matrixName.c_str()); 
  h_Resmat_2e2muJJ = (TH2*) ZZTo2e2muJJ_r->Get(matrixName.c_str()); 
  h_Resmat_4lmad   = (TH2*) ZZTo4lmad_r->Get(matrixName.c_str());
  h_Resmat_4lpow   = (TH2*) ZZTo4lpow_r->Get(matrixName.c_str()); 
  h_gg4mu   = (TH1*) ggZZTo4mu_r->Get(histoName_reco.c_str()); 
  h_gg4e    = (TH1*) ggZZTo4e_r->Get(histoName_reco.c_str()); 
  h_gg2e2mu = (TH1*) ggZZTo2e2mu_r->Get(histoName_reco.c_str()); 
  h_4muJJ   = (TH1*) ZZTo4muJJ_r->Get(histoName_reco.c_str()); 
  h_4eJJ    = (TH1*) ZZTo4eJJ_r->Get(histoName_reco.c_str()); 
  h_2e2muJJ = (TH1*) ZZTo2e2muJJ_r->Get(histoName_reco.c_str()); 
  h_4lmad   = (TH1*) ZZTo4lmad_r->Get(histoName_reco.c_str()); 
  h_4lpow   = (TH1*) ZZTo4lpow_r->Get(histoName_reco.c_str()); 
  h_gg4mu_gen   = (TH1*) ggZZTo4mu_g->Get(histoName_gen.c_str()); 
  h_gg4e_gen    = (TH1*) ggZZTo4e_g->Get(histoName_gen.c_str()); 
  h_gg2e2mu_gen = (TH1*) ggZZTo2e2mu_g->Get(histoName_gen.c_str()); 
  h_4muJJ_gen   = (TH1*) ZZTo4muJJ_g->Get(histoName_gen.c_str()); 
  h_4eJJ_gen    = (TH1*) ZZTo4eJJ_g->Get(histoName_gen.c_str()); 
  h_2e2muJJ_gen = (TH1*) ZZTo2e2muJJ_g->Get(histoName_gen.c_str()); 
  h_4lmad_gen   = (TH1*) ZZTo4lmad_g->Get(histoName_gen.c_str()); 
  h_4lpow_gen   = (TH1*) ZZTo4lpow_g->Get(histoName_gen.c_str()); 
 
  if(h_Resmat_gg4mu == NULL)   h_Resmat_gg4mu = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_gg4e == NULL)    h_Resmat_gg4e = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_gg2e2mu == NULL) h_Resmat_gg2e2mu = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe"); 
  if(h_Resmat_4muJJ == NULL)   h_Resmat_4muJJ = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_4eJJ == NULL)    h_Resmat_4eJJ = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_Resmat_2e2muJJ == NULL) h_Resmat_2e2muJJ = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_gg4mu == NULL)          h_gg4mu = (TH1*) h_safe->Clone("h_safe");
  if(h_gg4e == NULL)           h_gg4e = (TH1*) h_safe->Clone("h_safe");
  if(h_gg2e2mu == NULL)        h_gg2e2mu = (TH1*) h_safe->Clone("h_safe"); 
  if(h_4muJJ == NULL)          h_4muJJ = (TH1*) h_safe->Clone("h_safe");
  if(h_4eJJ == NULL)           h_4eJJ = (TH1*) h_safe->Clone("h_safe");
  if(h_2e2muJJ == NULL)        h_2e2muJJ = (TH1*) h_safe->Clone("h_safe");
  if(h_gg4mu_gen == NULL)      h_gg4mu_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_gg4e_gen == NULL)       h_gg4e_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_gg2e2mu_gen == NULL)    h_gg2e2mu_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_4muJJ_gen == NULL)      h_4muJJ_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_4eJJ_gen == NULL)       h_4eJJ_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_2e2muJJ_gen == NULL)    h_2e2muJJ_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_Resmat_4lpow == NULL)   h_Resmat_4lpow = (TH2*) h_Resmat_safe->Clone("h_Resmat_safe");
  if(h_4lpow == NULL)          h_4lpow = (TH2*) h_safe->Clone("h_safe");
  if(h_4lpow_gen == NULL)      h_4lpow_gen = (TH2*) h_safe->Clone("h_safe");
 
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
   

  if(mad ==1){
    h_Resmat_4lTot = (TH2*)h_Resmat_4lmad->Clone("h_Resmat_4lmad");
    h_4lTot = (TH1*) h_4lmad->Clone("h_4lmad");
    h_4lTot_gen = (TH1*) h_4lmad_gen->Clone("h_4lmad");
  }
  else{
    h_Resmat_4lTot = (TH2*)h_Resmat_4lpow->Clone("h_Resmat_4lpow");
    h_4lTot = (TH1*) h_4lpow->Clone("h_4lpow");
    h_4lTot_gen = (TH1*) h_4lpow_gen->Clone("h_4lpow");
  }

  h_Resmat_ggTot = (TH2*)h_Resmat_gg4mu->Clone("h_Resmat_gg"); 
    h_Resmat_JJTot = (TH2*)h_Resmat_4muJJ->Clone("h_Resmat_JJ"); 
  
  //cout <<"gg= "<< h_Resmat_gg4mu->Integral(0,1,0,50)<<" " << h_gg4mu->Integral(0,1)<<endl;  
  //cout <<"JJ= " <<h_Resmat_4muJJ->Integral(0,1,0,50)<<" " <<h_4muJJ->Integral(0,1)<<endl;
  //cout <<"4l= "<< h_Resmat_4lTot->Integral(0,1,0,50)<<" " <<h_4l->Integral(0,1)<<endl;
  
  h_Resmat = (TH2*)h_Resmat_4lTot->Clone("h_Resmat_4lTot");
  h_Resmat->Add(h_Resmat_ggTot); 
  h_Resmat->Add(h_Resmat_JJTot);
  h_4lTot_c = (TH1*) h_4lTot ->Clone("h_4lTot"); 
  h_4lTot_c->Add(h_gg4mu); 
  h_4lTot_c->Add(h_4muJJ);
 
  h_4lTot_gen_c = (TH1*) h_4lTot_gen ->Clone("h_4lTot_gen");
  h_4lTot_gen_c->Add(h_gg4mu_gen);
  h_4lTot_gen_c->Add(h_4muJJ_gen);
 
  totalint = h_Resmat->Integral();
 
  h_Resmat_normTot = (TH2*)h_Resmat->Clone("h_Resmat");
   
  h_Resmat_normTot->Scale(1/totalint);
  
  //std::cout << "total integral " << totalint << std::endl;
  //  string title = "Response Matrix m_{" + binning + "_" + finalstate + "} " + unc+dataset;   
 
  string matrixNameFile = "ResMat_qqggJJ_"+var+"_ZZTo" + finalstate + "_" + unc + "_" + dataset; 
  string matrixNormTotNameFile =  "ResMat_qqggJJ_"+var+"_normTot_ZZTo" + finalstate + "_" + unc + "_" + dataset; 
  string histoName_recoFile = var+"_qqggJJ_ZZTo" + finalstate + "_" + unc + "_" + dataset; 
  string histoName_genFile =  var+"Gen_qqggJJ_ZZTo" + finalstate + "_" + unc + "_" + dataset; 
  string histoName_recoFile_err =  var+"_statErr_qqggJJ_ZZTo" + finalstate + "_" + unc + "_" + dataset; 

  h_Resmat->SetTitle(matrixNameFile.c_str());
  h_Resmat_normTot->SetTitle(matrixNormTotNameFile.c_str());
  h_4lTot_c->SetTitle(histoName_recoFile.c_str());
  h_4lTot_gen_c->SetTitle(histoName_genFile.c_str());

  TH1 * h_4lTot_err = (TH1*) h_4lTot_c ->Clone("h_4lTot_c");
  // for(int j=1; j<10; j++){
  //   h_4lTot_err->SetBinContent(j,0.); 
  // } 
  float err_stat = 0;
  float bin = 0;
  for(int i =1; i <b; i++){
    err_stat = 0;
    bin = 0;
    bin = h_4lTot->GetBinContent(i);
    err_stat = sqrt(bin);
    h_4lTot_err->SetBinError(i,err_stat);
  }
  output->cd(); 
  h_Resmat->Write(matrixNameFile.c_str());
  h_Resmat_normTot->Write(matrixNormTotNameFile.c_str());
  h_4lTot_c->Write(histoName_recoFile.c_str());
  h_4lTot_gen_c->Write(histoName_genFile.c_str());
  h_4lTot_err->Write(histoName_recoFile_err.c_str());
  output->Close();
}

//Plot distributions
void ResponseMatrix::Plot(string var,string fs, string dataset, string unc, string path)
{
  gROOT->Reset();  
  gROOT->SetStyle("Plain");   
  gStyle->SetOptStat(0);
   
  setTDRStyle(); 
 
  int iPeriod = 4; 
  //int iPos = 11; 
  writeExtraText = true;    
  extraText  = "Simulation";
  extraText2 = "";

  string title;
  string xAxis;
  string yAxis;
  string finalstate;
  double max = 0;  
    
  if(fs == "4m") finalstate = "4#mu";
  else if(fs == "2e2m") finalstate = "2e2#mu";
  else finalstate = fs;
 
  file = new TFile((var+"_test/"+fileName).c_str());
  
  matrixName =  "ResMat_qqggJJ_"+var+"_ZZTo" + fs + "_" + unc+ "_" + dataset;
  TH2D *matrix = (TH2D*) file->Get(matrixName.c_str());
  TCanvas *c = new TCanvas("c","c");
  c->cd();  
  TPad  *pad1 = new TPad("pad1","", 0., 0., 1.0, 1.0);
  pad1->SetLeftMargin(0.15); 
  pad1->SetRightMargin(0.11); 
 
  pad1->Draw();
  pad1->cd();

  if(var =="Mass"){
    xAxis = "reco m_{"+finalstate+"} [GeV]";
    yAxis = "gen  m_{"+finalstate+"} [GeV]";
    max = matrix->GetBinContent(2,2)/2+3;
  } 
  if(var =="dRZZ"){
    xAxis = "reco #DeltaR(Z_1,Z_2)";
    yAxis = "gen  #DeltaR(Z_1,Z_2)";
    max = matrix->GetBinContent(4,4)/2+3;
  }
  else if(var =="Jets"){
    xAxis = "reco N jets (|#eta^{jet}|<4.7)";
    yAxis = "gen N jets (|#eta^{jet}|<4.7)"; 
    max = matrix->GetBinContent(1,1)/2;
  }
  else if(var =="Mjj"){
    xAxis = "reco m_{jj} (|#eta^{jet}|<4.7) [GeV]";
    yAxis = "gen m_{jj} (|#eta^{jet}|<4.7) [GeV]"; 
    max = matrix->GetBinContent(2,2)*1.5;
  }
  else if(var =="Deta"){
    title = "Response Matrix #Delta#eta_{jj} - "+finalstate+"final state (dataset: " + dataset+" unc: "+ unc+")";
    xAxis = "reco #Delta#eta_{jj} (|#eta^{jet}|<4.7)";
    yAxis = "gen #Delta#eta_{jj} (|#eta^{jet}|<4.7)"; 
    max = matrix->GetBinContent(2,2)*1.5;
  }
  else if(var =="Jets_Central"){
    xAxis = "reco N jets (|#eta^{jet}|<2.4)";
    yAxis = "gen N jets (|#eta^{jet}|<2.4)"; 
    max = matrix->GetBinContent(1,1)/3;
  }
  else if(var =="Mjj_Central"){
    title = "Response Matrix m_{jj} with |#eta^{j}| < 2.4 - "+finalstate+"final state (dataset: " + dataset+" unc: "+ unc+")";
    xAxis = "reco m_{jj} (|#eta^{jet}|<2.4) [GeV]";
    yAxis = "gen m_{jj} (|#eta^{jet}|<2.4) [GeV]"; 
    max = matrix->GetBinContent(2,2)*1.5;
  }
  else if(var =="Deta_Central"){
    xAxis = "reco #Delta#eta_{jj} (|#eta^{jet}|<2.4)";
    yAxis = "gen #Delta#eta_{jj} (|#eta^{jet}|<2.4)"; 
    max = matrix->GetBinContent(2,2)*1.5;
  } 
  else if(var =="PtJet1"){
    xAxis = "reco p_{T}^{jet1} [GeV]";
    yAxis = "gen p_{T}^{jet1} [GeV]";
    max = matrix->GetBinContent(1,1);
    matrix->GetXaxis()->SetRangeUser(30,500);
    matrix->GetYaxis()->SetRangeUser(30,500);
  }
  else if(var =="PtJet2"){
    xAxis = "reco p_{T}^{jet2} [GeV]";
    yAxis = "gen p_{T}^{jet2} [GeV]";
    max = matrix->GetBinContent(1,1);
    matrix->GetXaxis()->SetRangeUser(30,500);
    matrix->GetYaxis()->SetRangeUser(30,500);
  }
 else if(var =="EtaJet1"){
   xAxis = "reco |#eta^{jet1}|";
    yAxis = "gen |#eta^{jet1}|";
    max = matrix->GetBinContent(2,2)/2;
 }
 else if(var =="EtaJet2"){
   xAxis = "reco |#eta^{jet2}|";
    yAxis = "gen |#eta^{jet2}|";
    max = matrix->GetBinContent(2,2)/2;
 }
 if(var == "Jets" || var == "Jets_Central"){
    matrix->GetXaxis()->SetBinLabel(1,"0");
    matrix->GetXaxis()->SetBinLabel(2,"1");
    matrix->GetXaxis()->SetBinLabel(3,"2");
    matrix->GetXaxis()->SetBinLabel(4,">2");  
    matrix->GetXaxis()->SetLabelSize(0.05);
 }

 gStyle->SetPaintTextFormat("4.1f");
 // // PrecisionMatrix(matrix);
 // Int_t nbinsx = matrix->GetNbinsX();
 // Int_t nbinsy = matrix->GetNbinsY();
 // for(Int_t bx=1; bx <=nbinsx; bx++){
 //   for(Int_t by=1; by <=nbinsy; by++){
 //     Float_t NewVal = static_cast<int>(matrix->GetBinContent(bx,by) * 100) / 100.0f;
 //     matrix->SetBinContent(bx,by,NewVal);      
 //   }
 // }

 matrix->SetMaximum(max);
 matrix->GetXaxis()->SetTitle(xAxis.c_str());
 matrix->GetYaxis()->SetTitle(yAxis.c_str());
 matrix->SetMarkerColor(kGray+1);
 matrix->SetMarkerSize(1.4);
 matrix->Draw("COLZTEXT");
 matrix->GetXaxis()->SetTitleOffset(1.2);
 matrix->GetYaxis()->SetTitleOffset(1.5);
 
 lumiTextSize     = 0.7;
 cmsTextSize      = 0.7;
  extraOverCmsTextSize  = 0.80;//0.63; 
  CMS_lumi(pad1,iPeriod,0);
 
  string png ="~/www/VBS/"+path+"/"+var+"/"+"ResMat_qqggJJ_"+var+"_ZZTo" + fs + "_" + unc+ "_" + dataset + W + tightfr+ "_"+mc+".png";
  string pdf ="~/www/VBS/"+path+"/"+var+"/"+"ResMat_qqggJJ_"+var+"_ZZTo" + fs + "_" + unc+ "_" + dataset + W + tightfr+ "_"+mc+".pdf";
  
  c->Print(png.c_str());
  c->Print(pdf.c_str());
}

//Build response matrix, reco and gen distributions needed for the theoretical uncertainty on Powheg
void ResponseMatrix::GenMCSystDistributions(string var, string dataset, string finalstate, bool mad)
{
  //MC systematics on Pow
  if(mad == 0)  FolderNameMCSyst = "GenMCUpDownDistributions"+ tightfr+ "_Pow";
  else  FolderNameMCSyst = "GenMCUpDownDistributions"+ tightfr+ "_Mad";
  
  output = new TFile((FolderNameMCSyst+"/MCSystDistributions_"+var+".root").c_str(), "UPDATE");
  int b = 0;  
  
  if(var == "Mass") {
    variable = var;
    b=9;
  } 
  else if(var == "dRZZ") {
    variable = var;
    b=7;
  }
  else if(var == "Jets" ||var == "Jets_Central"){
    variable = var + "_JERSmear";
    b=6;
  }
  else if(var == "Mjj" ||var == "Mjj_Central" || var == "Deta" || var == "Deta_Central"){
    variable = var+"_JERSmear";
    b=3;
  }
  else if(var == "PtJet1" ||var == "PtJet2"||var == "EtaJet1" ||var == "EtaJet2" ){
    variable = var + "_JERSmear";
    b=6;
  }
 
  safeHistoName = "ZZTo2e2m_" + variable +"_"+ W +"01"; 
  histoName_gen =  "ZZTo" + finalstate + "_" + var + "Gen_" + W  +dataset+ tightfr;
  histoMCup = var + "_up_perc";
  histoMCdown = var + "_down_perc";


  h_safe_tmp = (TH1*) rescue->Get(safeHistoName.c_str());  
  h_safe = (TH1*)h_safe_tmp->Clone(safeHistoName.c_str()); 

  for(int k=1; k<b; k++){
    for(int l=1; l<b; l++){
      h_safe->SetBinContent(l,0.);
      h_safe->SetBinError(l,0.);  
    }
  }

  h_gg4mu_gen   = (TH1*) ggZZTo4mu_g->Get(histoName_gen.c_str()); 
  h_gg4e_gen    = (TH1*) ggZZTo4e_g->Get(histoName_gen.c_str()); 
  h_gg2e2mu_gen = (TH1*) ggZZTo2e2mu_g->Get(histoName_gen.c_str()); 
  h_4muJJ_gen   = (TH1*) ZZTo4muJJ_g->Get(histoName_gen.c_str()); 
  h_4eJJ_gen    = (TH1*) ZZTo4eJJ_g->Get(histoName_gen.c_str()); 
  h_2e2muJJ_gen = (TH1*) ZZTo2e2muJJ_g->Get(histoName_gen.c_str()); 
  h_4lmad_gen   = (TH1*) ZZTo4lmad_g->Get(histoName_gen.c_str());
  h_4lpow_gen   = (TH1*) ZZTo4lpow_g->Get(histoName_gen.c_str()); 
  
  if(h_gg4mu_gen   == NULL) h_gg4mu_gen   = (TH1*) h_safe->Clone("h_safe");
  if(h_gg4e_gen    == NULL) h_gg4e_gen    = (TH1*) h_safe->Clone("h_safe");
  if(h_gg2e2mu_gen == NULL) h_gg2e2mu_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_4muJJ_gen   == NULL) h_4muJJ_gen   = (TH1*) h_safe->Clone("h_safe");
  if(h_4eJJ_gen    == NULL) h_4eJJ_gen    = (TH1*) h_safe->Clone("h_safe");
  if(h_2e2muJJ_gen == NULL) h_2e2muJJ_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_4lpow_gen   == NULL) h_4lpow_gen   = (TH2*) h_safe->Clone("h_safe");//
   
  h_4lpow_up_gen   = (TH1*) h_4lpow_gen->Clone("h_4lpow_gen"); //
  h_4lpow_down_gen = (TH1*) h_4lpow_gen->Clone("h_4lpow_gen"); //
    
  if(h_4lpow_up_gen   == NULL) h_4lpow_up_gen   = (TH2*) h_safe->Clone("h_safe");//
  if(h_4lpow_down_gen == NULL) h_4lpow_down_gen = (TH2*) h_safe->Clone("h_safe");//
 
  h_4lpow_up_gen   = (TH1*) ZZMCsystPow_g->Get(histoMCup.c_str());
  h_4lpow_down_gen = (TH1*) ZZMCsystPow_g->Get(histoMCdown.c_str());
  
  int nbins = h_4lpow_gen->GetNbinsX();
  float up_cont = 0;
  float down_cont = 0;
  float def_cont_4lpow = 0;
  
  //reweight powheg distributions from scale variations by the default distribution used in the analysis
  for(int i=1; i <= nbins; i++){
    up_cont = 0;
    down_cont = 0;
    def_cont_4lpow = 0;
    
    up_cont   =    h_4lpowSist_up_gen->GetBinContent(i);  
    down_cont =    h_4lpowSist_down_gen->GetBinContent(i);  

    def_cont_4lpow = h_4lpow_gen->GetBinContent(i);
 
    h_4lpow_up_gen->SetBinContent(i,def_cont_4lpow*up_cont); 
    h_4lpow_down_gen->SetBinContent(i,def_cont_4lpow*down_cont);

    
  }

  h_gg4mu_gen->Add(h_gg4e_gen);
  h_gg4mu_gen->Add(h_gg2e2mu_gen);
  h_4muJJ_gen->Add(h_4eJJ_gen);
  h_4muJJ_gen->Add(h_2e2muJJ_gen);

  
  if(mad ==1){
    h_4lTot_up_gen = (TH1*) h_4lmad_gen->Clone("h_4lmad");
    h_4lTot_down_gen = (TH1*) h_4lmad_gen->Clone("h_4lmad");
  }
  else{
    h_4lTot_up_gen   = (TH1*) h_4lpow_up_gen->Clone("h_4l"); 
    h_4lTot_down_gen = (TH1*) h_4lpow_down_gen->Clone("h_4l");
  }

  h_4lTot_up_gen->Add(h_gg4mu_gen);
  h_4lTot_up_gen->Add(h_4muJJ_gen);
  h_4lTot_down_gen->Add(h_gg4mu_gen);
   h_4lTot_down_gen->Add(h_4muJJ_gen);
 
  string histoName_up =  var+"Gen_qqggJJ_ZZTo" + finalstate +"_up_"+ dataset; 
  string histoName_down =  var+"Gen_qqggJJ_ZZTo" + finalstate  +"_down_"+ dataset; 

  h_4lTot_up_gen->SetTitle(histoName_up.c_str());
  h_4lTot_down_gen->SetTitle(histoName_down.c_str());
  
  output->cd(); 
  h_4lTot_up_gen->Write(histoName_up.c_str());
  h_4lTot_down_gen->Write(histoName_down.c_str());
  output->Close();
  
}

//Build response matrix, reco and gen distributions needed for the theoretical uncertainty on MadGraph5_madatNLO
void ResponseMatrix::GenMGatNLOSystDistributions(string var, string dataset, string finalstate)
{
  
  FolderNameMCSyst = "GenMCUpDownDistributions"+ tightfr+ "_MGatNLO";
  output = new TFile((FolderNameMCSyst+"/MCSystDistributions_"+var+".root").c_str(), "UPDATE");
  
  int b = 0;  
  
  if(var == "Mass") {
    variable = var;
    b=9;
  } 
  else if(var == "dRZZ") {
    variable = var;
    b=7;
  }
  else if(var == "Jets" ||var == "Jets_Central"){
    variable = var + "_JERSmear";
    b=6;
  }
  else if(var == "Mjj" ||var == "Mjj_Central" || var == "Deta" || var == "Deta_Central"){
    variable = var+"_JERSmear";
    b=3;
  }
  else if(var == "PtJet1" ||var == "PtJet2"||var == "EtaJet1" ||var == "EtaJet2" ){
    variable = var + "_JERSmear";
    b=6;
  }
 
  safeHistoName = "ZZTo2e2m_" + variable +"_"+ W +"01"; 
  histoName_gen =  "ZZTo" + finalstate + "_" + var + "Gen_" + W  +dataset+ tightfr;
  histoMCup = var + "_"+ finalstate+"_up";
  histoMCdown = var + "_"+ finalstate+"_down";
  string histoMCcentral = var + "_"+ finalstate+"_default"; 
 
  h_safe_tmp = (TH1*) rescue->Get(safeHistoName.c_str());  
  h_safe = (TH1*)h_safe_tmp->Clone(safeHistoName.c_str()); 

  for(int k=1; k<b; k++){
    for(int l=1; l<b; l++){
      h_safe->SetBinContent(l,0.);
      h_safe->SetBinError(l,0.);  
    }
  }

  h_gg4mu_gen   = (TH1*) ggZZTo4mu_g->Get(histoName_gen.c_str()); 
  h_gg4e_gen    = (TH1*) ggZZTo4e_g->Get(histoName_gen.c_str()); 
  h_gg2e2mu_gen = (TH1*) ggZZTo2e2mu_g->Get(histoName_gen.c_str()); 
  h_4muJJ_gen   = (TH1*) ZZTo4muJJ_g->Get(histoName_gen.c_str()); 
  h_4eJJ_gen    = (TH1*) ZZTo4eJJ_g->Get(histoName_gen.c_str()); 
  h_2e2muJJ_gen = (TH1*) ZZTo2e2muJJ_g->Get(histoName_gen.c_str()); 
  h_4lmad_gen   = (TH1*) ZZTo4lmad_g->Get(histoName_gen.c_str());
  
  //get finalstate distribution for central-up-down
  h_4lpow_gen      = (TH1*)ZZMCsystMGatNLO_g->Get(histoMCcentral.c_str()); 
  h_4lpow_up_gen   = (TH1*)ZZMCsystMGatNLO_g->Get(histoMCup.c_str());
  h_4lpow_down_gen = (TH1*)ZZMCsystMGatNLO_g->Get(histoMCdown.c_str());
    
  if(h_gg4mu_gen      == NULL) h_gg4mu_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_gg4e_gen       == NULL) h_gg4e_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_gg2e2mu_gen    == NULL) h_gg2e2mu_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_4muJJ_gen      == NULL) h_4muJJ_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_4eJJ_gen       == NULL) h_4eJJ_gen = (TH1*) h_safe->Clone("h_safe");
  if(h_2e2muJJ_gen    == NULL) h_2e2muJJ_gen = (TH1*) h_safe->Clone("h_safe");

  if(h_4lpow_gen      == NULL) h_4lpow_gen = (TH2*) h_safe->Clone("h_safe");
  if(h_4lpow_up_gen    == NULL) h_4lpow_up_gen = (TH2*) h_safe->Clone("h_safe");
  if(h_4lpow_down_gen == NULL) h_4lpow_down_gen = (TH2*) h_safe->Clone("h_safe");
  
  h_gg4mu_gen->Add(h_gg4e_gen);
  h_gg4mu_gen->Add(h_gg2e2mu_gen);
  h_4muJJ_gen->Add(h_4eJJ_gen);
  h_4muJJ_gen->Add(h_2e2muJJ_gen);
 
  //scale distributions to the number of events of MGatLO
  TH1 *h_4lpow_cl      = (TH1*) h_4lpow_gen->Clone("h_4lpow_gen");
  TH1 *h_4lpow_up_cl   = (TH1*) h_4lpow_up_gen->Clone("h_4lpow_up_gen");
  TH1 *h_4lpow_down_cl = (TH1*) h_4lpow_down_gen->Clone("h_4lpow_down_gen");
  float norm = h_4lmad_gen->Integral()/h_4lpow_gen->Integral();

  h_4lpow_cl->Scale(norm);
  h_4lpow_up_cl->Scale(norm);
  h_4lpow_down_cl->Scale(norm);

  cout <<"final state: "<< finalstate <<  " MGatNLO = " << h_4lmad_gen->Integral() << " MGatNLO = " << h_4lpow_gen->Integral() << " MGatNLO Up = " << h_4lpow_up_gen->Integral() << " MGatNLO Down = " << h_4lpow_down_gen->Integral() << endl;
  cout << "norm = " << norm << endl;  
    cout <<"final state: "<< finalstate <<  " MGatNLO = " << h_4lpow_cl->Integral() << " MGatNLO Up = " << h_4lpow_up_cl->Integral() << " MGatNLO Down = " << h_4lpow_down_cl->Integral() << endl;
   
    cout << finalstate << " " << h_4lmad_gen->GetBinContent(1) << " "<< h_4lpow_cl->GetBinContent(1) << " " <<h_gg4mu_gen->GetBinContent(1) << " " << /* h_4muJJ_gen->GetBinContent(1) << " " << */h_4lpow_cl->GetBinContent(1) + h_gg4mu_gen->GetBinContent(1) /* +h_4muJJ_gen->GetBinContent(1) */<<endl;
    cout << finalstate << " " << h_4lmad_gen->GetBinContent(2) << " "<< h_4lpow_cl->GetBinContent(2) << " " <<h_gg4mu_gen->GetBinContent(2) << " " << /* h_4muJJ_gen->GetBinContent(2) << " " << */h_4lpow_cl->GetBinContent(2) + h_gg4mu_gen->GetBinContent(2) /*+h_4muJJ_gen->GetBinContent(2) */<<endl;
  h_4lTot_gen = (TH1*) h_4lpow_cl->Clone("h_4lpow_central"); 
  h_4lTot_up_gen = (TH1*) h_4lpow_up_cl->Clone("h_4lpow_up"); 
  h_4lTot_down_gen = (TH1*) h_4lpow_down_cl->Clone("h_4lpow_down");
 
  h_4lTot_gen->Add(h_gg4mu_gen);
  h_4lTot_gen->Add(h_4muJJ_gen);
  h_4lTot_up_gen->Add(h_gg4mu_gen);
  h_4lTot_up_gen->Add(h_4muJJ_gen);
  h_4lTot_down_gen->Add(h_gg4mu_gen);
  h_4lTot_down_gen->Add(h_4muJJ_gen);
  
  cout << finalstate <<  " tot = " <<  h_4lTot_gen->GetBinContent(1) << endl;
  string histoName_up =  var+"Gen_qqggJJ_ZZTo" + finalstate +"_up_"+ dataset; 
  string histoName_down =  var+"Gen_qqggJJ_ZZTo" + finalstate  +"_down_"+ dataset; 
  string histoName_central =  var+"Gen_qqggJJ_ZZTo" + finalstate  +"_central_"+ dataset; 
 
  h_4lTot_gen->SetTitle(histoName_central.c_str());
  h_4lTot_up_gen->SetTitle(histoName_up.c_str());
  h_4lTot_down_gen->SetTitle(histoName_down.c_str());
  
  output->cd(); 
  h_4lTot_gen->Write(histoName_central.c_str());
  h_4lTot_up_gen->Write(histoName_up.c_str());
  h_4lTot_down_gen->Write(histoName_down.c_str());
  output->Close();
  
}
