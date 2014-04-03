#include "DataFormats/Candidate/interface/Candidate.h"

namespace phys{
  phys::Particle convert(const reco::Candidate &rc){
    phys::Particle p(rc.p4(),phys::Particle::computeCharge(rc.pdgId()),rc.pdgId());
    p.setMotherId(rc.mother()->pdgId());
    
    return p;
  }
}
