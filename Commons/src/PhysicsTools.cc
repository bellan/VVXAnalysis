#include "VVXAnalysis/Commons/interface/PhysTools.h"

phys::Particle phys::convert(const reco::Candidate &rc){
  
  phys::Particle p(rc.p4(),phys::Particle::computeCharge(rc.pdgId()),rc.pdgId());
  p.setMotherId(rc.numberOfMothers() > 0 ? rc.mother()->pdgId() : -9999);

  return p;
}
