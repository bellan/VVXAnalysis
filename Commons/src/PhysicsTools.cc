#include "VVXAnalysis/Commons/interface/PhysTools.h"

phys::Particle phys::convert(const reco::Candidate &rc,std::bitset<15> flags){
  
  phys::Particle p(rc.p4(),phys::Particle::computeCharge(rc.pdgId()),rc.pdgId(),flags);

  p.setMotherId(rc.numberOfMothers() > 0 ? rc.mother()->pdgId() : -9999);

  return p;
}

 phys::Particle phys::convert(const reco::Candidate &rc){
  
   std::bitset<15> flags (-99);
   return phys::convert(rc,flags);
  }

