#ifndef H6f_H
#define H6f_H

#include "TH1F.h"
#include "VVXAnalysis/Producers/interface/Boson.h"
#include "TString.h"
#include <vector>

class TFile;


class H6f {

 public:
  H6f(TString name_);
  H6f(TString name_, TFile* file);

  void Fill(Boson *V0, Boson *V1, Boson *V2); 

  void Scale(float w);

  void SetLineColor(Color_t c);

  TString name;

  TH1F* h6fMass;

  TH1F* h2l0Mass;
  TH1F* h2l1Mass;
  TH1F* hjjMass;

  TH1F* h4lMass;

  TH1F* h2l0Pt;
  TH1F* h2l1Pt;
  TH1F* hjjPt;

  TH1F* h4lPt_1;
  TH1F* h4lPt_2;
  TH1F* h4lPt_3;
  TH1F* h4lPt_4;
  TH1F* hjPt_1;
  TH1F* hjPt_2;

};

#endif
