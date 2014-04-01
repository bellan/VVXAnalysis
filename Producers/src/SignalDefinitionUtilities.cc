//-----------FUNCTION: definition of the two ZZ bosons from leptons-------

#include "VVXAnalysis/Producers/interface/SignalDefinitionUtilities.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using std::cout;
using std::endl;

std::pair<phys::Boson<phys::Particle> ,phys::Boson<phys::Particle> > makeZBosonsFromLeptons(const std::vector<const reco::Candidate *>& lm, const std::vector<const reco::Candidate *>& lp, int leptonCode, float mZ){
  std::vector<phys::Particle> plm;
  std::vector<phys::Particle> plp;

  foreach(const reco::Candidate *gp, lm)
    plm.push_back(convert(*gp));
  
  foreach(const reco::Candidate *gp, lp)
    plp.push_back(convert(*gp));
   
  return makeZBosonsFromLeptons(plm, plp, leptonCode, mZ);
}
  
  
std::pair<phys::Boson<phys::Particle> ,phys::Boson<phys::Particle> > makeZBosonsFromLeptons(const std::vector<phys::Particle>& lm, const std::vector<phys::Particle>& lp, int leptonCode, float mZ){
    
    phys::Boson<phys::Particle> Z0;
    phys::Boson<phys::Particle> Z1;
    float minMDiff=99999.;
    if (leptonCode == 4) {
      for (int k=0; k<2; ++k) {
	for (int j=0; j<2; ++j) {
	  float mDiff = fabs((lp[k].p4() + lm[j].p4()).M() - mZ);
	  if ( mDiff < minMDiff ) {
	    minMDiff=mDiff;   
	    
	    Z0.setDaughter(0,lp[k]);
	    Z0.setDaughter(1,lm[j]);

	    Z1.setDaughter(0,lp[(k+1)%2]);
	    Z1.setDaughter(1,lm[(j+1)%2]);
	  }      
	} 	
      }
    }
    else { 
      for (int z=0; z<2; ++z) {
	if ( fabs(lp[z].id()) == fabs(lm[0].id()) ) { 
	  
	  Z0.setDaughter(0,lp[z]);
	  Z0.setDaughter(1,lm[0]);

	  Z1.setDaughter(0,lp[(z+1)%2]);
	  Z1.setDaughter(1,lm[1]);      
	}
      }	
    }

    Z0.setId(23);
    Z1.setId(23);
    
    return std::make_pair(Z0,Z1);
  }



int makeVBosonsFromIds(int j0Id, int j1Id) {    
  
  if ( abs(j0Id) < 6 && abs(j1Id ) < 6) {
    if( (j0Id*j1Id) <0 && (abs(j0Id + j1Id) == 1 || abs(j0Id + j1Id) == 3) ) {
      if( j0Id % 2 == 0 )       return copysign(24,j0Id);  // W
      else if( j1Id % 2 == 0 )  return copysign(24,j1Id);  // W
      else return 0;
    }
    else if( j0Id + j1Id == 0 ) return 23;                 // Z
    else return 0;                             
    
  }
  else return 0;
}

phys::Particle convert(const reco::Candidate &rc){
  
  phys::Particle p(rc.p4(),phys::Particle::computeCharge(rc.pdgId()),rc.pdgId());
  p.setMotherId(rc.mother()->pdgId());
  
  return p;
}
