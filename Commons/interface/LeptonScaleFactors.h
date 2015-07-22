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

#include <TH2F.h>

class LeptonScaleFactors{
 public:
  LeptonScaleFactors(const std::string& muonEffFilename, const std::string& electronEFFfilename,
		     const std::string& muonFRFilename, const std::string& electronFRFilename);

  double efficiencyScaleFactor(const double& lepPt, const double& lepEta, int lepId) const;
  double efficiencyScaleFactor(const phys::Lepton& lep) const;

  double efficiencyScaleFactorErr(const double& lepPt, const double& lepEta, int lepId) const; 
  double efficiencyScaleFactorErr(const phys::Lepton& lep) const;

  double weight(const phys::DiBoson<phys::Lepton,phys::Lepton> &ZZ) const;
  double weight(const phys::Boson<phys::Lepton> &Z) const;

  std::pair<double,double> fakeRateScaleFactor(const double& lepPt, const double& lepEta, int lepId) const;
  std::pair<double,double> fakeRateScaleFactor(const phys::Lepton& lep) const;
 private:
  TH2F *hEffMu_;
  TH2F *hEffEl_;

  std::pair<TH1D*,TH1D*> hFRMu_;
  std::pair<TH1D*,TH1D*> hFREl_;
};


#endif
