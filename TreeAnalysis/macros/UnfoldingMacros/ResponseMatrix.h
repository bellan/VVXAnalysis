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
#include "tdrstyle.h"
#include "CMS_lumi.h"

class ResponseMatrix : public TObject
{
  
 public:
  
  ResponseMatrix();
  ResponseMatrix(bool weight, bool madgraph, bool tightregion);
  ~ResponseMatrix();
  void Build(string var, string dataset, string finalstate, int xs_qq, int xs_gg, bool mad);
  //  void Build_SF(string var, string dataset, string finalstate,string unc,bool mad);
  // void Build_JE(string var, string dataset, string finalstate,string unc,bool mad);
  void Build_Syst(string var, string dataset, string finalstate,string unc,bool mad);
  void Plot(string var,string fs, string dataset, string unc,string path); 
  void GenMCSystDistributions(string var, string dataset, string finalstate, bool mad);
  void GenMGatNLOSystDistributions(string var, string dataset, string finalstate);
  void CloseFiles();
  
  TFile *output; 
  TFile *gg4l_r;
  TFile *ZZTo4lpow_r;
  TFile *ZZTo4lmad_r;
  TFile *qq4l2j_r;
 
  //Truth samples (signal definition distributions)
  TFile *gg4l_g;
  TFile *ZZTo4lpow_g;
  TFile *ZZTo4lmad_g;
  TFile *qq4l2j_g;
  TFile *ZZMCsystPow_g;
  TFile *ZZMCsystMGatNLO_g;

  TFile *file;
  string fileName;
  string fileName_Syst;
  string FolderNameMCSyst;

  TH2 * h_Resmat;
  TH2 * h_Resmat_4lpow;
  TH2 * h_Resmat_4lmad;

  TH2 * h_Resmat_gg4l;
  TH2 * h_Resmat_qq4l2j;  
  TH2 * h_Resmat_4lTot;
  TH2 * h_Resmat_ggTot;
  TH2 * h_Resmat_JJTot;
    
  TH1 * h_4lpow;
  TH1 * h_4lmad;
  TH1 * h_gg4l;
  TH1 * h_qq4l2j;
  TH1 * h_4lTot;

  TH1 * h_4lpow_gen;
  TH1 * h_4lmad_gen;
  TH1 * h_gg4l_gen;
  TH1 * h_qq4l2j_gen;
   
  TH1 * h_4lTot_gen;
  TH1 * h_4lTot_up_gen; 
  TH1 * h_4lTot_down_gen; 
  TH1 * h_4lpow_up_gen; 
  TH1 * h_4lpow_down_gen; 
  TH1 * h_4lpowSist_up_gen;
  TH1 * h_4lpowSist_down_gen;

  TH1 * h_4lTot_c;
  TH1 * h_4lTot_gen_c;  
  TH2 * h_Resmat_normTot;

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
