#ifndef PHOTONSFHELPER_H
#define PHOTONSFHELPER_H

#include <memory>  // unique_ptr

#include "TH2F.h"

#include <FWCore/ParameterSet/interface/FileInPath.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>


enum SFsyst {central = 0, up = 1, down = 2};

class PhotonSFHelper
{
 public:

  PhotonSFHelper(int year, bool preVFP);
  
  float getSF      (/*string phoId,*/ float pt, float eta/*, float SCeta*/) const;
  float getSFError (/*string phoId,*/ float pt, float eta/*, float SCeta*/) const;
   
 private:
  // Photon SF histogram
  std::unique_ptr<TH2F> h_Pho;
  //std::unique_ptr<TH2F> h_Pho_2016preVFP, h_Pho_2016postVFP, h_Pho_2017, h_Pho_2018;
      
};

#endif
