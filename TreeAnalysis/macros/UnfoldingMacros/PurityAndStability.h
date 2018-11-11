#ifndef PURITYANDSTABILITY_h
#define PURITYANDSTABILITY_h

#include <TROOT.h>
//#include <TChain.h>
//#include <TFile.h>
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


class PurityAndStability : public TObject
{
  
 public:
  
  PurityAndStability();
  PurityAndStability(bool mad);
  ~PurityAndStability();
  void Build(string var, string finalstate);
  void Plot(string var,  string finalstate,string path); 
  void Plot_PAS(string var,string finalstate,string path); 
  
  TFile *output; 
  TFile *matrixFile;
  TFile *file;
  TFile * madgraph;
  TFile * powheg; 
  
  TH2 * h_Resmat; 
  TH1 * h_purity;
  TH1 * h_stability; 
  TH1 * p_mad;
  TH1 * s_mad;
  TH1 * p_pow;
  TH1 * s_pow;
  string matrixName;
  string histoName;
  string fileName;
  string matrixFileName;
  string mc;
  string fileNameMad;
  string fileNamePow ;
  string pName;
  string sName;

  ClassDef(PurityAndStability,1)

};

#endif
