#if !defined(__CINT__) || defined(__MAKECINT__)
#include "ResponseMatrix.h"
#include "CMS_lumi.C"

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

    gg4l_r       = new TFile("../../results/ZZRecoAnalyzer_SR/gg_4l.root");
    ZZTo4lpow_r  = new TFile("../../results/ZZRecoAnalyzer_SR/ZZTo4l.root");
    ZZTo4lmad_r  = new TFile("../../results/ZZRecoAnalyzer_SR/ZZTo4lamcatnlo.root");
    qq4l2j_r     = new TFile("../../results/ZZRecoAnalyzer_SR/qq_4l2j.root"); 
    
    //Truth samples (signal definition distributions) 
    gg4l_g       = new TFile("../../results/ZZMCAnalyzer_MC/gg_4l.root");
    ZZTo4lpow_g  = new TFile("../../results/ZZMCAnalyzer_MC/ZZTo4l.root");
    ZZTo4lmad_g  = new TFile("../../results/ZZMCAnalyzer_MC/ZZTo4lamcatnlo.root");
    qq4l2j_g     = new TFile("../../results/ZZMCAnalyzer_MC/qq_4l2j.root"); 

    if(madgraph ==1)  {
      fileName    = "matrices" + tightfr+ "_Mad.root";
      fileName_Syst = "matrices"+tightfr+ "_Syst_Mad.root";
      mc = "Mad";
    } 
    else {
      fileName    = "matrices" + tightfr+ "_Pow.root"; 
      fileName_Syst = "matrices" + tightfr+ "_Syst_Pow.root";
      mc = "Pow";
    }
    W = "";
  }
  
  else{

    //Reco samples (response matrices and signal region distributions) 

    ZZTo4lpow_r   = new TFile("../../results/ZZRecoWAnalyzer_SR/ZZTo4l.root");
    ZZTo4lmad_r   = new TFile("../../results/ZZRecoWAnalyzer_SR/ZZTo4lamcatnlo.root");
    gg4l_r        = new TFile("../../results/ZZRecoWAnalyzer_SR/gg_4l.root");
    qq4l2j_r      = new TFile("../../results/ZZRecoWAnalyzer_SR/qq_4l2j.root");
    
    //Truth samples (signal definition distributions) 
    ZZTo4lmad_g   = new TFile("../../results/ZZMCWAnalyzer_MC/ZZTo4lamcatnlo.root");
    ZZTo4lpow_g   = new TFile("../../results/ZZMCWAnalyzer_MC/ZZTo4l.root");
    gg4l_g        = new TFile("../../results/ZZMCWAnalyzer_MC/gg_4l.root");
    qq4l2j_g      = new TFile("../../results/ZZMCWAnalyzer_MC/qq_4l2j.root");

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

  // gg4l_r->Close();  
  // ZZTo4lpow_r->Close();  
  // ZZTo4lmad_r->Close();
  // qq4l2j_r->Close();  
  // gg4l_g->Close();
  // ZZTo4lpow_g->Close();
  // ZZTo4lmad_g->Close();
  // qq4l2j_g->Close();
}



void ResponseMatrix::CloseFiles()
{
  
  gg4l_r->Close();  
  ZZTo4lpow_r->Close();  
  ZZTo4lmad_r->Close();
  qq4l2j_r->Close();  
  gg4l_g->Close();
  ZZTo4lpow_g->Close();
  ZZTo4lmad_g->Close();
  qq4l2j_g->Close();
}


ResponseMatrix::~ResponseMatrix(){}

//Build the standard response matrix, reco and gen distributions;
//Buind response matrices and distributions varying the gg->ZZ and qq->ZZ cross sections by their uncertainties
void ResponseMatrix::Build(string var, string dataset, string finalstate, int xs_qq, int xs_gg, bool mad)
{ 
  output = new TFile((var+"_test/"+fileName).c_str(), "UPDATE");

  variable = var ; 

  matrixName      = "ResMat_ZZTo" + finalstate + "_" + variable+"_"+ W + dataset+ tightfr;
  //if(dataset=="01")  histoName_reco  = "ZZTo" + finalstate + "_" + variable+"_"+ W +dataset+tightfr;
  //else  histoName_reco  = "ZZTo" + finalstate + "_" + variable+"_"+ W +dataset;
  histoName_reco  = "ZZTo" + finalstate + "_" + variable+"_"+ W +dataset;
  histoName_gen   = "ZZTo" + finalstate + "_" + var + "Gen_" + W  +dataset+ tightfr;

  float unc_qq = 0; 
  float unc_gg = 0; 
  float totalint = 0; 
  
  //cout<<" "<<histoName_reco<<endl;

  h_Resmat_gg4l    = (TH2*) gg4l_r->Get(matrixName.c_str()); 
  h_Resmat_qq4l2j  = (TH2*) qq4l2j_r->Get(matrixName.c_str()); 
  h_Resmat_4lmad   = (TH2*) ZZTo4lmad_r->Get(matrixName.c_str());
  h_Resmat_4lpow   = (TH2*) ZZTo4lpow_r->Get(matrixName.c_str());
  h_gg4l           = (TH1*) gg4l_r->Get(histoName_reco.c_str());
  h_qq4l2j         = (TH1*) qq4l2j_r->Get(histoName_reco.c_str()); 
  h_4lpow          = (TH1*) ZZTo4lpow_r->Get(histoName_reco.c_str()); 
  h_4lmad          = (TH1*) ZZTo4lmad_r->Get(histoName_reco.c_str()); 
  h_gg4l_gen       = (TH1*) gg4l_g->Get(histoName_gen.c_str()); 
  h_qq4l2j_gen     = (TH1*) qq4l2j_g->Get(histoName_gen.c_str()); 
  h_4lpow_gen      = (TH1*) ZZTo4lpow_g->Get(histoName_gen.c_str());
  h_4lmad_gen      = (TH1*) ZZTo4lmad_g->Get(histoName_gen.c_str());

  if(h_Resmat_gg4l   == NULL){ cout<<"histogram "<<h_Resmat_gg4l  ->GetName()<<" is null. Abort"<<endl;  abort(); }  
  if(h_Resmat_qq4l2j == NULL){ cout<<"histogram "<<h_Resmat_qq4l2j->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_gg4l          == NULL){ cout<<"histogram "<<h_gg4l         ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_qq4l2j        == NULL){ cout<<"histogram "<<h_qq4l2j       ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_gg4l_gen      == NULL){ cout<<"histogram "<<h_gg4l_gen     ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_qq4l2j_gen    == NULL){ cout<<"histogram "<<h_qq4l2j_gen   ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_Resmat_4lpow  == NULL){ cout<<"histogram "<<h_Resmat_4lpow ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_Resmat_4lmad  == NULL){ cout<<"histogram "<<h_Resmat_4lmad ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_4lpow         == NULL){ cout<<"histogram "<<h_4lpow        ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_4lpow_gen     == NULL){ cout<<"histogram "<<h_4lpow_gen    ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_4lmad         == NULL){ cout<<"histogram "<<h_4lmad        ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_4lmad_gen     == NULL){ cout<<"histogram "<<h_4lmad_gen    ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  
  TH2 * h_Resmat_gg4l_cl    = (TH2*) h_Resmat_gg4l->Clone("h_Resmat_4lpow");
  TH1 * h_gg4l_cl           = (TH1*) h_gg4l->Clone("h_4l");
  TH1 * h_gg4l_gen_cl       = (TH1*) h_gg4l_gen->Clone("h_4l");
  TH2 * h_Resmat_qq4l2j_cl    = (TH2*) h_Resmat_qq4l2j->Clone("h_Resmat_4lpow");
  TH1 * h_qq4l2j_cl           = (TH1*) h_qq4l2j->Clone("h_4l");
  TH1 * h_qq4l2j_gen_cl       = (TH1*) h_qq4l2j_gen->Clone("h_4l");
  
  
  //cout << finalstate << " 4l " <<  h_4lpow_cl->Integral() << " gg " << h_gg4l_cl->Integral() << " qqJJ " << h_qq4l2j_cl->Integral() << endl;
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

  h_Resmat_ggTot = (TH2*)h_Resmat_gg4l_cl->Clone("h_Resmat_gg"); 
  h_Resmat_JJTot = (TH2*)h_Resmat_qq4l2j_cl->Clone("h_Resmat_JJ"); 
  
  // cout <<"gg= "<< h_Resmat_gg4l->Integral(0,1,0,50)<<" " << h_gg4l->Integral(0,1)<<endl;  
  // cout <<"JJ= " <<h_Resmat_qq4l2j->Integral(0,1,0,50)<<" " <<h_qq4l2j->Integral(0,1)<<endl;
  // cout <<"4l= "<< h_Resmat_4lTot->Integral(0,1,0,50)<<" " <<h_4l->Integral(0,1)<<endl;
  
  unc_qq = 0.0285; //pdf: 3.4%  scale: 2.85%  %use only scale because pdf is set somewhere else
  unc_gg = 0.083; //pdf: 3.10%  scale: 8.3%   %use only scale because pdf is set somewhere else
 
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
    h_gg4l_cl->Scale(1);
    h_gg4l_gen_cl->Scale(1);
  }
  else if(xs_gg == 1) {
    h_Resmat_ggTot->Scale(1+unc_gg); 
    h_gg4l_cl->Scale(1+unc_gg); 
    h_gg4l_gen_cl->Scale(1+unc_gg);
  }
  else if(xs_gg == -1) {
    h_Resmat_ggTot->Scale(1-unc_gg); 
    h_gg4l_cl->Scale(1-unc_gg); 
    h_gg4l_gen_cl->Scale(1-unc_gg);
  }
  else std::cout << "Error: xs_gg must be -1, 0 or 1" << std::endl;

  h_Resmat = (TH2*)h_Resmat_4lTot->Clone("h_Resmat_4lTot");
  //To comment only if you do NOT want to use MCFM and Phanton!!
  h_Resmat->Add(h_Resmat_ggTot); 
  h_Resmat->Add(h_Resmat_JJTot);
  
  h_4lTot_c = (TH1*) h_4lTot ->Clone("h_4lTot"); 
    //To comment only if you do NOT want to use MCFM and Phanton!!
  h_4lTot_c->Add(h_gg4l_cl); 
  h_4lTot_c->Add(h_qq4l2j_cl);
 
  //cout << " tot integral " << h_4lTot_c->Integral() << endl;

  h_4lTot_gen_c = (TH1*) h_4lTot_gen ->Clone("h_4lTot_gen");
  //To comment only if you do NOT want to use MCFM and Phanton!!
  h_4lTot_gen_c->Add(h_gg4l_gen_cl);
  h_4lTot_gen_c->Add(h_qq4l2j_gen_cl);
  
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
    if(xs_gg == 0)       unc = "_qqp_"; 
    else if(xs_gg == 1)  unc = "_qqp_ggp_";
    else if(xs_gg == -1) unc = "_qqp_ggm_";
  }
  else if(xs_qq == -1){
    if(xs_gg == 0)       unc = "_qqm_"; 
    else if(xs_gg == 1)  unc = "_qqm_ggp_";
    else if(xs_gg == -1) unc = "_qqm_ggm_";
  }

  string matrixNameFile         = "ResMat_qqggJJ_"+var+"_ZZTo" + finalstate + unc + dataset; 
  string matrixNormTotNameFile  = "ResMat_qqggJJ_"+var+"_normTot_ZZTo" + finalstate + unc + dataset; 
  string histoName_recoFile     = var+"_qqggJJ_ZZTo" + finalstate + unc + dataset; 
  string histoName_genFile      = var+"Gen_qqggJJ_ZZTo" + finalstate + unc + dataset; 
  //  string histoName_recoFile_err = var+"_statErr_qqggJJ_ZZTo" + finalstate + unc + dataset; 


  h_Resmat->SetTitle(matrixNameFile.c_str());
  h_Resmat_normTot->SetTitle(matrixNormTotNameFile.c_str());
  h_4lTot_c->SetTitle(histoName_recoFile.c_str());
  h_4lTot_gen_c->SetTitle(histoName_genFile.c_str());
  output->cd(); 
  h_Resmat->Write(matrixNameFile.c_str(),TObject::kOverwrite);
  h_Resmat_normTot->Write(matrixNormTotNameFile.c_str(),TObject::kOverwrite);
  h_4lTot_c->Write(histoName_recoFile.c_str(),TObject::kOverwrite);
  h_4lTot_gen_c->Write(histoName_genFile.c_str(),TObject::kOverwrite);
  output->Close();

  h_Resmat_gg4l  =NULL;
  h_Resmat_qq4l2j=NULL;
  h_Resmat_4lmad =NULL;
  h_Resmat_4lpow =NULL;
  h_gg4l         =NULL;
  h_qq4l2j       =NULL;
  h_4lpow        =NULL;
  h_4lmad        =NULL;
  h_gg4l_gen     =NULL;
  h_qq4l2j_gen   =NULL;
  h_4lpow_gen    =NULL;
  h_4lmad_gen    =NULL; 

}

void ResponseMatrix::Build_Syst(string var, string dataset, string finalstate, string unc, bool mad)
{

  output = new TFile((var+"_test/"+fileName_Syst).c_str(), "UPDATE");
  string uncGen = "";
  if(unc=="PDFUp") uncGen="_pdfUp";
  else if (unc=="PDFDn") uncGen="_pdfDn";
  else if (unc=="AsUp") uncGen="_asMZUp";
  else if (unc=="AsDn") uncGen="_asMZDn";

  matrixName = "ResMat_ZZTo" + finalstate + "_"+var+"_"+unc+"_"+ dataset + tightfr;
  histoName_reco = "ZZTo" + finalstate + "_"+var+"_"+unc +"_"+ dataset;
  histoName_gen =  "ZZTo" + finalstate + "_"+var+"Gen_"+ dataset + tightfr+ uncGen;


  //  cout<<"unc" <<unc<<" matrixName "<<matrixName<<" histoName_reco "<<histoName_reco<<" histoName_gen "<<histoName_gen<<endl;
  float totalint = 0; 
  h_Resmat_gg4l   = (TH2*) gg4l_r->Get(matrixName.c_str()); 
  h_Resmat_qq4l2j = (TH2*) qq4l2j_r->Get(matrixName.c_str()); 
  h_Resmat_4lmad  = (TH2*) ZZTo4lmad_r->Get(matrixName.c_str());
  h_Resmat_4lpow  = (TH2*) ZZTo4lpow_r->Get(matrixName.c_str()); 
  h_gg4l          = (TH1*) gg4l_r->Get(histoName_reco.c_str()); 
  h_qq4l2j        = (TH1*) qq4l2j_r->Get(histoName_reco.c_str()); 
  h_4lmad         = (TH1*) ZZTo4lmad_r->Get(histoName_reco.c_str()); 
  h_4lpow         = (TH1*) ZZTo4lpow_r->Get(histoName_reco.c_str());   
  h_gg4l_gen      = (TH1*) gg4l_g->Get(histoName_gen.c_str()); 
  h_qq4l2j_gen    = (TH1*) qq4l2j_g->Get(histoName_gen.c_str()); 
  h_4lmad_gen     = (TH1*) ZZTo4lmad_g->Get(histoName_gen.c_str());
  h_4lpow_gen     = (TH1*) ZZTo4lpow_g->Get(histoName_gen.c_str()); 
  
  if(h_Resmat_gg4l   == NULL){ cout<<"histogram "<<h_Resmat_gg4l  ->GetName()<<" is null. Abort"<<endl;  abort(); }  
  if(h_Resmat_qq4l2j == NULL){ cout<<"histogram "<<h_Resmat_qq4l2j->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_gg4l          == NULL){ cout<<"histogram "<<h_gg4l         ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_qq4l2j        == NULL){ cout<<"histogram "<<h_qq4l2j       ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_gg4l_gen      == NULL){ cout<<"histogram "<<h_gg4l_gen     ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_qq4l2j_gen    == NULL){ cout<<"histogram "<<h_qq4l2j_gen   ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_Resmat_4lpow  == NULL){ cout<<"histogram "<<h_Resmat_4lpow ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_Resmat_4lmad  == NULL){ cout<<"histogram "<<h_Resmat_4lmad ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_4lpow         == NULL){ cout<<"histogram "<<h_4lpow        ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_4lpow_gen     == NULL){ cout<<"histogram "<<h_4lpow_gen    ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_4lmad         == NULL){ cout<<"histogram "<<h_4lmad        ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_4lmad_gen     == NULL){ cout<<"histogram "<<h_4lmad_gen    ->GetName()<<" is null. Abort"<<endl;  abort(); } 

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
  h_Resmat_ggTot = (TH2*)h_Resmat_gg4l->Clone("h_Resmat_gg"); 
  h_Resmat_JJTot = (TH2*)h_Resmat_qq4l2j->Clone("h_Resmat_JJ"); 
  
  //cout <<"gg= "<< h_Resmat_gg4l->Integral(0,1,0,50)<<" " << h_gg4l->Integral(0,1)<<endl;  
  //cout <<"JJ= " <<h_Resmat_qq4l2j->Integral(0,1,0,50)<<" " <<h_qq4l2j->Integral(0,1)<<endl;
  //cout <<"4l= "<< h_Resmat_4lTot->Integral(0,1,0,50)<<" " <<h_4l->Integral(0,1)<<endl;
  
  h_Resmat = (TH2*)h_Resmat_4lTot->Clone("h_Resmat_4lTot");
  h_Resmat->Add(h_Resmat_ggTot); 
  h_Resmat->Add(h_Resmat_JJTot);
 
  h_4lTot_c = (TH1*) h_4lTot ->Clone("h_4lTot"); 
  h_4lTot_c->Add(h_gg4l); 
  h_4lTot_c->Add(h_qq4l2j);
 
  h_4lTot_gen_c = (TH1*) h_4lTot_gen ->Clone("h_4lTot_gen");
  h_4lTot_gen_c->Add(h_gg4l_gen);
  h_4lTot_gen_c->Add(h_qq4l2j_gen);
 
  totalint = h_Resmat->Integral();
 
  h_Resmat_normTot = (TH2*)h_Resmat->Clone("h_Resmat");   
  h_Resmat_normTot->Scale(1/totalint);
  
  //std::cout << "total integral " << totalint << std::endl;
  
  string matrixNameFile        = "ResMat_qqggJJ_"+var+"_ZZTo" + finalstate + "_" + unc + "_" + dataset; 
  string matrixNormTotNameFile = "ResMat_qqggJJ_"+var+"_normTot_ZZTo" + finalstate + "_" + unc + "_" + dataset; 
  string histoName_recoFile    = var+"_qqggJJ_ZZTo" + finalstate + "_" + unc + "_" + dataset; 
  string histoName_genFile     = var+"Gen_qqggJJ_ZZTo" + finalstate + "_" + unc + "_" + dataset; 

  h_Resmat->SetTitle(matrixNameFile.c_str());
  h_Resmat_normTot->SetTitle(matrixNormTotNameFile.c_str());
  h_4lTot_c->SetTitle(histoName_recoFile.c_str());
  h_4lTot_gen_c->SetTitle(histoName_genFile.c_str());
 
  output->cd(); 
  h_Resmat->Write(matrixNameFile.c_str(),TObject::kOverwrite);
  h_Resmat_normTot->Write(matrixNormTotNameFile.c_str(),TObject::kOverwrite);
  h_4lTot_c->Write(histoName_recoFile.c_str(),TObject::kOverwrite);
  h_4lTot_gen_c->Write(histoName_genFile.c_str(),TObject::kOverwrite);
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
  else if(var =="dRZZ"){
    xAxis = "reco #DeltaR(Z_1,Z_2)";
    yAxis = "gen  #DeltaR(Z_1,Z_2)";
    max = matrix->GetBinContent(4,4)/2+3;
  }
  else if(var =="PtZZ"){
    xAxis = "reco p_{T}^{4#ell}";
    yAxis = "gen  p_{T}^{4#ell}";
    max = matrix->GetBinContent(4,4)/2+3;
  }
  else if(var =="nJets"){
    xAxis = "reco N jets (|#eta^{jet}|<4.7)";
    yAxis = "gen N jets (|#eta^{jet}|<4.7)"; 
    max = matrix->GetBinContent(1,1)/2;
  }
  else if(var =="nIncJets"){
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
  else if(var =="nJets_Central"){
    xAxis = "reco N jets (|#eta^{jet}|<2.4)";
    yAxis = "gen N jets (|#eta^{jet}|<2.4)"; 
    max = matrix->GetBinContent(1,1)/3;
  }
  else if(var =="nIncJets_Central"){
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
 if(var == "nJets" || var == "nJets_Central"){

    matrix->GetXaxis()->SetBinLabel(1,"0");
    matrix->GetXaxis()->SetBinLabel(2,"1");
    matrix->GetXaxis()->SetBinLabel(3,"2");
    matrix->GetXaxis()->SetBinLabel(4,">2");  
    //    matrix->GetXaxis()->SetBinLabel(5,">3");  
    matrix->GetXaxis()->SetLabelSize(0.05);

    matrix->GetYaxis()->SetBinLabel(1,"0");
    matrix->GetYaxis()->SetBinLabel(2,"1");
    matrix->GetYaxis()->SetBinLabel(3,"2");
    matrix->GetYaxis()->SetBinLabel(4,">2");  
    //    matrix->GetYaxis()->SetBinLabel(5,">3");  
    matrix->GetYaxis()->SetLabelSize(0.05);
 }

 if(var == "nIncJets" || var == "nIncJets_Central"){

    matrix->GetXaxis()->SetBinLabel(1,"#geq0");
    matrix->GetXaxis()->SetBinLabel(2,"#geq1");
    matrix->GetXaxis()->SetBinLabel(3,"#geq2");
    matrix->GetXaxis()->SetBinLabel(4,"#geq3");  
    //    matrix->GetXaxis()->SetBinLabel(5,">3");  
    matrix->GetXaxis()->SetLabelSize(0.05);

    matrix->GetYaxis()->SetBinLabel(1,"#geq0");
    matrix->GetYaxis()->SetBinLabel(2,"#geq1");
    matrix->GetYaxis()->SetBinLabel(3,"#geq2");
    matrix->GetYaxis()->SetBinLabel(4,"#geq3");  
    //    matrix->GetYaxis()->SetBinLabel(5,">3");  
    matrix->GetYaxis()->SetLabelSize(0.05);
 }

 gStyle->SetPaintTextFormat("4.2f");
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
 std::string SavePage = "~/www/PlotsVV/13TeV/";
 string png = SavePage+path+"/"+var+"/"+"ResMat_qqggJJ_"+var+"_ZZTo" + fs + "_" + unc+ "_" + dataset + W + tightfr+ "_"+mc+".png";
 string pdf = SavePage+path+"/"+var+"/"+"ResMat_qqggJJ_"+var+"_ZZTo" + fs + "_" + unc+ "_" + dataset + W + tightfr+ "_"+mc+".pdf";

 c->Print(png.c_str());
 c->Print(pdf.c_str());
 c->Delete();
}

//Build response matrix, reco and gen distributions needed for the theoretical uncertainty on Powheg
void ResponseMatrix::GenMCSystDistributions(string var, string dataset, string finalstate, bool mad)
{
  //MC systematics on Pow
  if(mad == 0)  FolderNameMCSyst = "GenMCUpDownDistributions"+ tightfr+ "_Pow";
  else  FolderNameMCSyst = "GenMCUpDownDistributions"+ tightfr+ "_Mad";
  
  output = new TFile((FolderNameMCSyst+"/MCSystDistributions_"+var+".root").c_str(), "UPDATE");
 
  histoName_gen =  "ZZTo" + finalstate + "_" + var + "Gen_" + W  +dataset+ tightfr;
  histoMCup = var + "_up_perc";
  histoMCdown = var + "_down_perc";



  h_gg4l_gen    = (TH1*) gg4l_g->Get(histoName_gen.c_str()); 
  h_qq4l2j_gen  = (TH1*) qq4l2j_g->Get(histoName_gen.c_str()); 
  h_4lmad_gen   = (TH1*) ZZTo4lmad_g->Get(histoName_gen.c_str());
  h_4lpow_gen   = (TH1*) ZZTo4lpow_g->Get(histoName_gen.c_str()); 

  if(h_gg4l_gen      == NULL){ cout<<"histogram "<<h_gg4l_gen     ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_qq4l2j_gen    == NULL){ cout<<"histogram "<<h_qq4l2j_gen   ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_4lpow_gen     == NULL){ cout<<"histogram "<<h_4lpow        ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_4lmad_gen     == NULL){ cout<<"histogram "<<h_4lmad        ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  

  h_4lpow_up_gen   = (TH1*) h_4lpow_gen->Clone("h_4lpow_gen"); //
  h_4lpow_down_gen = (TH1*) h_4lpow_gen->Clone("h_4lpow_gen"); //

  //check hot
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

  if(mad ==1){
    h_4lTot_up_gen = (TH1*) h_4lmad_gen->Clone("h_4lmad");
    h_4lTot_down_gen = (TH1*) h_4lmad_gen->Clone("h_4lmad");
  }
  else{
    h_4lTot_up_gen   = (TH1*) h_4lpow_up_gen->Clone("h_4l"); 
    h_4lTot_down_gen = (TH1*) h_4lpow_down_gen->Clone("h_4l");
  }

  h_4lTot_up_gen->Add(h_gg4l_gen);
  h_4lTot_up_gen->Add(h_qq4l2j_gen);
  h_4lTot_down_gen->Add(h_gg4l_gen);
   h_4lTot_down_gen->Add(h_qq4l2j_gen);
 
  string histoName_up =  var+"Gen_qqggJJ_ZZTo" + finalstate +"_up_"+ dataset; 
  string histoName_down =  var+"Gen_qqggJJ_ZZTo" + finalstate  +"_down_"+ dataset; 

  h_4lTot_up_gen->SetTitle(histoName_up.c_str());
  h_4lTot_down_gen->SetTitle(histoName_down.c_str());
  
  output->cd(); 
  h_4lTot_up_gen->Write(histoName_up.c_str(),TObject::kOverwrite);
  h_4lTot_down_gen->Write(histoName_down.c_str(),TObject::kOverwrite);
  output->Close();
  
}

//Build response matrix, reco and gen distributions needed for the theoretical uncertainty on MadGraph5_madatNLO
void ResponseMatrix::GenMGatNLOSystDistributions(string var, string dataset, string finalstate)
{
  
  FolderNameMCSyst = "GenMCUpDownDistributions"+ tightfr+ "_MGatNLO";
  output = new TFile((FolderNameMCSyst+"/MCSystDistributions_"+var+".root").c_str(), "UPDATE");

  variable = var; //check hot
 
  histoName_gen =  "ZZTo" + finalstate + "_" + var + "Gen_" + W  +dataset+ tightfr;
  histoMCup = var + "_"+ finalstate+"_up";
  histoMCdown = var + "_"+ finalstate+"_down";
  string histoMCcentral = var + "_"+ finalstate+"_default"; 
 

  h_gg4l_gen   = (TH1*) gg4l_g->Get(histoName_gen.c_str()); 
  h_qq4l2j_gen   = (TH1*) qq4l2j_g->Get(histoName_gen.c_str()); 
  h_4lmad_gen   = (TH1*) ZZTo4lmad_g->Get(histoName_gen.c_str());
  
  //get finalstate distribution for central-up-down
  h_4lpow_gen      = (TH1*)ZZMCsystMGatNLO_g->Get(histoMCcentral.c_str()); 
  h_4lpow_up_gen   = (TH1*)ZZMCsystMGatNLO_g->Get(histoMCup.c_str());
  h_4lpow_down_gen = (TH1*)ZZMCsystMGatNLO_g->Get(histoMCdown.c_str());
    

  if(h_gg4l_gen       == NULL){ cout<<"histogram "<<h_gg4l_gen       ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_qq4l2j_gen     == NULL){ cout<<"histogram "<<h_qq4l2j_gen     ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_4lpow_gen      == NULL){ cout<<"histogram "<<h_4lpow_gen      ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_4lpow_up_gen   == NULL){ cout<<"histogram "<<h_4lpow_up_gen   ->GetName()<<" is null. Abort"<<endl;  abort(); } 
  if(h_4lpow_down_gen == NULL){ cout<<"histogram "<<h_4lpow_down_gen ->GetName()<<" is null. Abort"<<endl;  abort(); } 


 
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
   
    cout << finalstate << " " << h_4lmad_gen->GetBinContent(1) << " "<< h_4lpow_cl->GetBinContent(1) << " " <<h_gg4l_gen->GetBinContent(1) << " " << /* h_qq4l2j_gen->GetBinContent(1) << " " << */h_4lpow_cl->GetBinContent(1) + h_gg4l_gen->GetBinContent(1) /* +h_qq4l2j_gen->GetBinContent(1) */<<endl;
    cout << finalstate << " " << h_4lmad_gen->GetBinContent(2) << " "<< h_4lpow_cl->GetBinContent(2) << " " <<h_gg4l_gen->GetBinContent(2) << " " << /* h_qq4l2j_gen->GetBinContent(2) << " " << */h_4lpow_cl->GetBinContent(2) + h_gg4l_gen->GetBinContent(2) /*+h_qq4l2j_gen->GetBinContent(2) */<<endl;
  h_4lTot_gen = (TH1*) h_4lpow_cl->Clone("h_4lpow_central"); 
  h_4lTot_up_gen = (TH1*) h_4lpow_up_cl->Clone("h_4lpow_up"); 
  h_4lTot_down_gen = (TH1*) h_4lpow_down_cl->Clone("h_4lpow_down");
 
  h_4lTot_gen->Add(h_gg4l_gen);
  h_4lTot_gen->Add(h_qq4l2j_gen);
  h_4lTot_up_gen->Add(h_gg4l_gen);
  h_4lTot_up_gen->Add(h_qq4l2j_gen);
  h_4lTot_down_gen->Add(h_gg4l_gen);
  h_4lTot_down_gen->Add(h_qq4l2j_gen);
  
  cout << finalstate <<  " tot = " <<  h_4lTot_gen->GetBinContent(1) << endl;
  string histoName_up =  var+"Gen_qqggJJ_ZZTo" + finalstate +"_up_"+ dataset; 
  string histoName_down =  var+"Gen_qqggJJ_ZZTo" + finalstate  +"_down_"+ dataset; 
  string histoName_central =  var+"Gen_qqggJJ_ZZTo" + finalstate  +"_central_"+ dataset; 
 
  h_4lTot_gen->SetTitle(histoName_central.c_str());
  h_4lTot_up_gen->SetTitle(histoName_up.c_str());
  h_4lTot_down_gen->SetTitle(histoName_down.c_str());
  
  output->cd(); 
  h_4lTot_gen->Write(histoName_central.c_str(),TObject::kOverwrite);
  h_4lTot_up_gen->Write(histoName_up.c_str(),TObject::kOverwrite);
  h_4lTot_down_gen->Write(histoName_down.c_str(),TObject::kOverwrite);
  output->Close();  
}
