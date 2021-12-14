#ifndef VVXAnalysis_Commons_LeptonScaleFactors_H
#define VVXAnalysis_Commons_LeptonScaleFactors_H

/** \class LeptonScaleFactors
 *  This class acts as a wrapper of the efficiencis and scale factors measurements. On top of that, it adds some tools to correct the ZZ event.
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.3 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "ZZAnalysis/AnalysisStep/interface/LeptonSFHelper.h"

#include <TH2F.h>
#include <TGraphAsymmErrors.h>

class LeptonScaleFactors{
 public:
  LeptonScaleFactors(int year, const std::string& muonFRFilename, const std::string& electronFRFilename, bool preVFP);

  std::pair<double, double> efficiencyScaleFactor(const double& lepPt, const double& lepEta, int lepId, bool isInCracks = false) const{ 
    
    // eta is copied twice because of bad design of LeptonSFHelper class. The eta of the lepton is done beforehand and then it is treated differently in the lepSFHelper for electrons and muons.
    return std::make_pair(lepSFHelper_.getSF(year_, lepId, lepPt, lepEta, lepEta, isInCracks),lepSFHelper_.getSFError(year_, lepId, lepPt, lepEta, lepEta, isInCracks));
  }

  std::pair<double, double> efficiencyScaleFactor(const phys::Lepton& lep) const;

  double weight(const phys::DiBoson<phys::Lepton,phys::Lepton> &ZZ) const;
  double weight(const phys::Boson<phys::Lepton> &Z) const;

  std::pair<double,double > fakeRateScaleFactor(const double& lepPt, const double& lepEta, int lepId) const;
  std::pair<double,double > fakeRateScaleFactor(const phys::Lepton& lep) const;
 private:

  LeptonSFHelper lepSFHelper_;
  int year_;


  TH2D *hFRMu_;
  TH2D *hFREl_;

  std::string year;
};


#endif
