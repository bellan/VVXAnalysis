#include "VVXAnalysis/Producers/interface/FrixioneIsoCalculator.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <cmath>

double FrixioneIsoCalculator::maxEnergyFraction(double delta, double delta0) const{
  // Case nExponent_ != 1 is not handled
  return epsilon_*(std::cos(1 - delta)/std::cos(1 - delta0));
}

bool FrixioneIsoCalculator::isIsolated(const reco::Candidate* photon, const edm::View<reco::Candidate>& genParticles_cand, double delta0){
  cacheVector(genParticles_cand);
  return isIsolated(photon, delta0);
}

bool FrixioneIsoCalculator::isIsolated(const reco::Candidate* photon, double delta0){
  // Sort gen particles by distance to the photon
  std::sort(cachedGenParticles_.begin(), cachedGenParticles_.end(), [&](const reco::Candidate* part1, const reco::Candidate* part2) {
      return deltaR(*part1, *photon) < deltaR(*part2, *photon);
    });

  // Compute the isolation
  bool excluded_self = false;
  double frixione_sum = 0.;
  double photon_pt = photon->pt();
  for(const reco::Candidate* gp : cachedGenParticles_){
    // Exclude the photon iself from the computation
    if(!excluded_self && gp->p4() == photon->p4()){ // Same particle in different collections
      // Comparing `LorentzVector`s is probably costly, and there is only one genParticle that corresponds
      // to this photon. Therefore I added a guard bool so that the comparison is only done once
      excluded_self = true;
      continue;
    }

    // Compute the isolation
    double dr = deltaR(*gp, *photon);
    if (dr >= delta0)
      break;

    frixione_sum += gp->pt();

    double e_max = photon_pt*maxEnergyFraction(dr, delta0);
    if (frixione_sum > e_max)
      return false;
  }

  return true;
}

template <class T=reco::GenParticle>
void FrixioneIsoCalculator::cacheVector(const edm::View<T>& genParticles_cand){
  cachedGenParticles_.clear();
  cachedGenParticles_.reserve(genParticles_cand.size());
  for(auto& gp : genParticles_cand)
    cachedGenParticles_.emplace_back(static_cast<const reco::Candidate*>(&gp));
}

// Explicit instantiation of template
template void FrixioneIsoCalculator::cacheVector(const edm::View<reco::GenParticle>&);
template void FrixioneIsoCalculator::cacheVector(const edm::View<reco::GenJet     >&);
template void FrixioneIsoCalculator::cacheVector(const edm::View<reco::Candidate  >&);
