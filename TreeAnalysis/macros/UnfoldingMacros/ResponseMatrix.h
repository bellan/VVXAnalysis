#ifndef RESPONSEMATRIX_h
#define RESPONSEMATRIX_h
/** \class ResponseMatrix
 * 
 *  \author L. Finco - UNITO <linda.finco@cern.ch>
 */
#include <TROOT.h>
#include <TObject.h>
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

class ResponseMatrix : public TObject
{
  
 public:
  
  ResponseMatrix();
  ResponseMatrix(bool weight, bool madgraph, bool tightregion);
  ~ResponseMatrix();
  void Build(string var, string dataset, string finalstate, int xs_qq, int xs_gg, bool mad);
  void Build_SF(string var, string dataset, string finalstate,string unc,bool mad);
  void Build_JE(string var, string dataset, string finalstate,string unc,bool mad);
  void Plot(string var,string fs, string dataset, string unc,string path); 
  void GenMCSystDistributions(string var, string dataset, string finalstate, bool mad);
  void GenMGatNLOSystDistributions(string var, string dataset, string finalstate);
  
  TFile *output; 
  TFile *ggZZTo2e2mu_r;
  TFile *ggZZTo4e_r;
  TFile *ggZZTo4mu_r;
  TFile *ZZJetsTo4l_r;
  TFile *ZZTo2e2muJJ_r;
  TFile *ZZTo4eJJ_r;
  TFile *ZZTo4muJJ_r;
  TFile *ZZTo2e2mu_r;
  TFile *ZZTo4e_r;
  TFile *ZZTo4mu_r;
 
  //Truth samples (signal definition distributions)
  TFile *ggZZTo2e2mu_g;
  TFile *ggZZTo4e_g;
  TFile *ggZZTo4mu_g;
  TFile *ZZJetsTo4l_g;
  TFile *ZZTo2e2muJJ_g;
  TFile *ZZTo4eJJ_g;
  TFile *ZZTo4muJJ_g;
  TFile *ZZTo2e2mu_g;
  TFile *ZZTo4e_g;
  TFile *ZZTo4mu_g;
  TFile *ZZMCsystPow_g;
  TFile *ZZMCsystMGatNLO_g;

  //rescue file
  TFile *rescue;

  TFile *file;
  string fileName;
  string fileName_SF;
  string fileName_JE;
  string FolderNameMCSyst;

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
  TH2 * h_Resmat_qq4mu; 
  TH2 * h_Resmat_qq4e;
  TH2 * h_Resmat_qq2e2mu; 

  TH1 * h_4l;
  TH1 * h_gg4mu;
  TH1 * h_gg4e;
  TH1 * h_gg2e2mu;
  TH1 * h_4muJJ;
  TH1 * h_4eJJ;
  TH1 * h_2e2muJJ;
  TH1 * h_qq4mu;
  TH1 * h_qq4e;
  TH1 * h_qq2e2mu; 
  TH1 * h_4lTot;

  TH1 * h_4l_gen;
  TH1 * h_gg4mu_gen;
  TH1 * h_gg4e_gen;
  TH1 * h_gg2e2mu_gen;
  TH1 * h_4muJJ_gen;
  TH1 * h_4eJJ_gen;
  TH1 * h_2e2muJJ_gen;
  TH1 * h_qq4mu_gen; 
  TH1 * h_qq4e_gen;
  TH1 * h_qq2e2mu_gen;
  TH1 * h_4lTot_gen;
  TH1 * h_4lTot_up_gen; 
  TH1 * h_4lTot_down_gen; 

  TH1 * h_qq4mu_up_gen; 
  TH1 * h_qq4e_up_gen;
  TH1 * h_qq2e2mu_up_gen; 
  TH1 * h_qq4mu_down_gen; 
  TH1 * h_qq4e_down_gen;
  TH1 * h_qq2e2mu_down_gen;

  TH1 * h_qq_down_gen;
  TH1 * h_qq_up_gen;


  TH1 * h_4lTot_c;
  TH1 * h_4lTot_gen_c; 
 
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
  string mc;
  string tightfr; 
  string histoMCup;
  string histoMCdown; 
  string ZZMCsystNamePow;
  string ZZMCsystNameMGatNLO;
 
  ClassDef(ResponseMatrix,1)

};

#endif
