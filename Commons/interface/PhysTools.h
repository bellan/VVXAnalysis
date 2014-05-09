#ifndef VVXAnalysis_Commons_PhysTools_H
#define VVXAnalysis_Commons_PhysTools_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"

namespace phys{
  phys::Particle convert(const reco::Candidate &rc);
}
#endif
