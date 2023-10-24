#include "VVXAnalysis/Producers/interface/FrixioneIsoCalculator.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <cmath>
#include <sstream>

double FrixioneIsoCalculator::maxEnergyFraction(double delta, double delta0) const{
  // Case nExponent_ != 1 is not handled
  return epsilon_*(std::cos(1 - delta)/std::cos(1 - delta0));
  // static double cached_part = epsilon_/std::cos(1 - delta0_);
  // return cached_part*std::cos(1 - delta);
}

bool FrixioneIsoCalculator::isIsolated(const reco::GenParticle* photon, const edm::View<reco::GenParticle>& genParticles_cand, double delta0){
  cacheVector(genParticles_cand);
  return isIsolated(photon, delta0);
}

bool FrixioneIsoCalculator::isIsolated(const reco::GenParticle* photon, double delta0){
  // Sort gen particles by distance to the photon
  std::sort(cachedGenParticles_.begin(), cachedGenParticles_.end(), [&](const reco::GenParticle* part1, const reco::GenParticle* part2) {
      return deltaR(*part1, *photon) < deltaR(*part2, *photon);
    });

  // Compute the isolation
  std::stringstream ss;
  ss << "Total particles: " << cachedGenParticles_.size() << " - ph p4: " << photon->p4() << '\n';
  bool excluded_self = false;
  double frixione_sum = 0.;
  double photon_pt = photon->pt();
  for(const reco::GenParticle* gp : cachedGenParticles_){
    ss << "\tgp p4: " << gp->p4();
    // Exclude the photon iself from the computation
    if(!excluded_self && gp->p4() == photon->p4()){ // Same particle in different collections
      ss << " (*)\n";
      // Comparing `LorentzVector`s is probably costly, and there is only one genParticle that corresponds
      // to this photon. Therefore I added a guard bool so that the comparison is only done once
      excluded_self = true;
      continue;
    }

    // Compute the isolation
    double dr = deltaR(*gp, *photon);
    ss << "\t - dr = " << dr;
    if (dr >= delta0){
      ss <<" - breaking out of loop\n";
      break;
    }
    ss << "\t - energy: " << gp->pt() << '\n';

    frixione_sum += gp->pt();

    double e_max = photon_pt*maxEnergyFraction(dr, delta0);
    if (frixione_sum > e_max){
      ss << "FAIL - sum (" << frixione_sum << ") > e_max ("<<e_max<<")\n";
      edm::LogInfo("FrixioneIsoCalculator") << ss.str();
      return false;
    }
  }

  double e_max = photon->pt()*maxEnergyFraction(delta0, delta0);
  ss     << "PASS - sum (" << frixione_sum << ") < e_max ("<<e_max<<")\n";
  edm::LogInfo("FrixioneIsoCalculator") << ss.str();
  return true;
}

void FrixioneIsoCalculator::cacheVector(const edm::View<reco::GenParticle>& genParticles_cand){
  cachedGenParticles_.clear();
  cachedGenParticles_.reserve(genParticles_cand.size());
  for(auto& gp : genParticles_cand)
    cachedGenParticles_.emplace_back(&gp);
}

void FrixioneIsoCalculator::cacheVector(std::vector<const reco::GenParticle*>& genParticles_v){
  cachedGenParticles_.clear();
  cachedGenParticles_.swap(genParticles_v);
}

