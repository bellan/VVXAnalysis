#ifndef Hjets_H
#define Hjets_H

#include "TH1F.h"
#include "TString.h"
#include <vector>

class TFile;
namespace reco {
  class Candidate;
}


class Hjets {

 public:
  Hjets(TString name_);
  Hjets(TString name_, TFile* file);

  void Filljet(std::vector<const reco::Candidate *> theGenq);
  
  void Scale(float w);
  
  void SetLineColor(Color_t c);

  TString name;
  
  TH1F* hjjDeta;
  TH1F* hjjDtheta;
  TH1F* hjjDphi;
  TH1F* hjjDR;
  
};

#endif
