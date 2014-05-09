#include "VVXAnalysis/Commons/interface/PhysTools.h"

phys::Particle phys::convert(const reco::Candidate &rc){
  
  phys::Particle p(rc.p4(),phys::Particle::computeCharge(rc.pdgId()),rc.pdgId());
  p.setMotherId(rc.mother()->pdgId());
  
  return p;
}
