#ifndef VVXAnalysis_Commons_PhysTools_H
#define VVXAnalysis_Commons_PhysTools_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include "VVXAnalysis/DataFormats/interface/GenParticle.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include "VVXAnalysis/DataFormats/interface/GenStatusBit.h"

namespace phys{
  phys::GenParticle convert(const reco::Candidate &rc);
  phys::GenParticle convert(const reco::Candidate &rc,std::bitset<15> flags);

}
#endif
