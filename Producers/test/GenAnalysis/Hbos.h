#ifndef Hbos_H
#define Hbos_H

#include "TH1F.h"
#include "TString.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include <vector>

class TFile;

class Hbos {
  
 public:
  Hbos(TString name_);
  Hbos(TString name_, TFile* file);
  
  void FillBos(const phys::Boson<phys::Particle> &Z0, const phys::Boson<phys::Particle> &Z1, const phys::Boson<phys::Particle> &V);

  void Scale(float w);
  
  void SetLineColor(Color_t c);
    
  TString name;
  
  TH1F* hZ0Mass;
  TH1F* hZ1Mass;
  TH1F* hVMass;

  TH1F* hZPt_1;
  TH1F* hZPt_2;
  TH1F* hVPt;
  TH1F* hZZPt;

  TH1F* hZVDR;
  TH1F* hZZDR;
  TH1F* hZZ_VDR;

  TH1F* hZZDeta;
  TH1F* hZZ_VDeta;

};

#endif
