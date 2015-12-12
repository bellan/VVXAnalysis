#ifndef VVXAnalysis_Commons_PhysTools_H
#define VVXAnalysis_Commons_PhysTools_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include "VVXAnalysis/DataFormats/interface/GenStatusBit.h"

namespace phys{
  phys::Particle convert(const reco::Candidate &rc);
  phys::Particle convert(const reco::Candidate &rc,std::bitset<15> flags);

}
#endif
