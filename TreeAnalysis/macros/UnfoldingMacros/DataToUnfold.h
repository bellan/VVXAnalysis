#ifndef DATATOUNFOLD_h
#define DATATOUNFOLD_h

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

class DataToUnfold : public TObject
{
  
 public:
  
  DataToUnfold();
  ~DataToUnfold();
  void Build(string var, string finalstate);
  void Build_JE(string var, string finalstate);
  void Plot(string var,string finalstate,string path); 

  TFile *output; 
  TFile *output_syst; 
  TFile *data;
  TFile *red;
  TFile *ttZ;
  //  TFile *ttWW;
  TFile *WWZ;
  TFile *Irr;
  TFile *file;
  TFile *file_syst;
  //  TFile *rescue;

  TH1 * h_data;
  TH1 * h_data_up; 
  TH1 * h_data_down;
  TH1 * h_totdata;
  TH1 * h_totdata_up;
  TH1 * h_totdata_down;
  TH1 * h_red; 
  TH1 * h_red_up;
  TH1 * h_red_down;
  TH1 * h_ttZ;
  // TH1 * h_ttWW;
  TH1 * h_WWZ;
  TH1 * h_Irr;
  TH1 * h_data_irrp;
  TH1 * h_data_irrm;
  TH1 * h_data_redp;
  TH1 * h_data_redm; 
  //  TH1 * h_safe_tmp;
  //  TH1 * h_safe;
  string variable;
  //string XaxisTitle;
  string fs;  
  string histoName; 
  string histoName_up; 
  string histoName_down; 
  string histoMCName;
  //  string safeHistoName;
  //gstring safeHistoName_mass;
  string dataName;
  string TotdataName;
  string dataName_up;
  string TotdataName_up;
  string dataName_down;
  string TotdataName_down;
  string dataIrrpName;
  string dataIrrmName;
  string dataRedpName;
  string dataRedmName;
  
  double err_red_tot;

  ClassDef(DataToUnfold,1)

};

#endif
