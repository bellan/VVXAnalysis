#ifndef FrixioneIsolation_h
#define FrixioneIsolation_h

/** \class FrixioneIsolation
*
* Calculates the Frixione isolation for photons; see:
*  - https://journals.aps.org/prd/pdf/10.1103/PhysRevD.105.052003 pag. 12
*  - S. Frixione, Isolated photons in perturbative QCD, Phys. Lett. B 429, 369 (1998).
*
* $Date: 2023/10/16 $
* $Revision: 1.0 $
* \author A. Mecca - UniTo
*/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include <vector>

class FrixioneIsoCalculator {
 public:
  explicit FrixioneIsoCalculator(/* double delta0=0.4,  */double epsilon=1.0/* , int nExponent=1 */)
    /* : delta0_   (delta0) */
    : epsilon_  (epsilon)
    /* , nExponent_(nExponent) */
  {}

  bool isIsolated(const reco::Candidate *photon, const edm::View<reco::Candidate>    &genParticles, double delta0);
  bool isIsolated(const reco::Candidate *photon, double delta0);

  template <class T>
  void cacheVector(const edm::View<T>& genParticles);

 protected:
  double maxEnergyFraction(double delta, double delta0) const; // max (transverse) energy divided by photon pt

 private:
  // Parameters of the equation for E_T^max
  /* double delta0_;  // Isolation cone radius of the photon */
  double epsilon_;
  /* int nExponent_; */

  // Cached vector of pointers (useful in case this is called multiple times per event)
  std::vector<const reco::Candidate*> cachedGenParticles_;
};

#endif
