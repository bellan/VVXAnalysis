#ifndef VVXAnalysis_Commons_LeptonEfficiency_H
#define VVXAnalysis_Commons_LeptonEfficiency_H

/** \class LeptonEfficiency
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

class LeptonEfficiency{
 public:
  LeptonEfficiency();
  double scaleFactor(const double& lepPt, const double& lepEta, int lepId) const;
  double scaleFactor(const phys::Lepton& lep) const;
  double weight(const phys::DiBoson<phys::Lepton,phys::Lepton> &ZZ) const;
  double weight(const phys::Boson<phys::Lepton> &Z) const;

 private:
  TH2F *hMu_;
  TH2F *hEl_;
 
};


#endif
