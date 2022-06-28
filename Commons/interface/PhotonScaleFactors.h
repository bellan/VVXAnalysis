#ifndef VVXAnalysis_Commons_PhotonScaleFactors_H
#define VVXAnalysis_Commons_PhotonScaleFactors_H

/** \class LeptonScaleFactors
 *  This class acts as a wrapper of the efficiencis and scale factors measurements. On top of that, it adds some tools to correct the ZZ event.
 *
 *  $Date: 2022/05/24 13:37:31 $
 *  $Revision: 1.3 $
 *  \author A. Vagnerini - UNITO <antonio.vagnerini@cern.ch>
 */

#include "VVXAnalysis/DataFormats/interface/Photon.h"
//#include "ZZAnalysis/AnalysisStep/interface/LeptonSFHelper.h"
#include "VVXAnalysis/Commons/interface/PhotonSFHelper.h" 

#include <TH2F.h>
#include <TGraphAsymmErrors.h>


class PhotonScaleFactors{
 public:
  PhotonScaleFactors(int year, const std::string& Filename, bool preVFP); //constructor

  //Photon efficiency SF 
  std::pair<double, double> efficiencyScaleFactor(const double& phoPt, const double& phoEta /*, string phoId*/) const{ 
    return std::make_pair(phoSFHelper_.getSF(year_, /*phoId,*/ phoPt, phoEta),phoSFHelper_.getSFError(year_, /*phoId,*/ phoPt, phoEta));
  }

  std::pair<double, double> efficiencyScaleFactor(const phys::Photon& photon) const;

  //Correct weight of event
  double weight(const phys::Photon& photon) const;

  /*
  //Photon fake rate
  std::pair<double,double > fakeRateScaleFactor(const double& phoPt, const double& phoEta, int phoId) const;
  std::pair<double,double > fakeRateScaleFactor(const phys::Photon& photon) const;
  */

 private:

  PhotonSFHelper phoSFHelper_;
  int year_;

  TH2D *f;
  //TH2D *hFRMu_;
  //TH2D *hFREl_;

  std::string year;
};


#endif
