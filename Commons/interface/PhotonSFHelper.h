#ifndef PHOTONSFHELPER_H
#define PHOTONSFHELPER_H

#include <string>
#include <iostream>
#include <vector>
#include <utility>

#include <cmath>
#include "TString.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH2D.h"

#include <FWCore/ParameterSet/interface/FileInPath.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>


enum SFsyst {central = 0, up = 1, down = 2};

class PhotonSFHelper
{

 public:

  PhotonSFHelper(bool preVFP);
  ~PhotonSFHelper();
  
  float getSF      (int year, /*string phoId,*/ float pt, float eta/*, float SCeta*/) const;
  float getSFError (int year, /*string phoId,*/ float pt, float eta/*, float SCeta*/) const;
   
 private:
   TFile *root_file;
   
   // Photon SF map histograms
   TH2F *h_Pho_2016, *h_Pho_2017, *h_Pho_2018;
      
};

#endif
