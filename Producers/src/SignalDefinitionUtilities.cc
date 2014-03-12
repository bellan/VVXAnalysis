//-----------FUNCTION: definition of the two ZZ bosons from leptons-------

#include "VVXAnalysis/Producers/interface/SignalDefinitionUtilities.h"

using std::cout;
using std::endl;

std::pair<phys::Boson<phys::Particle> ,phys::Boson<phys::Particle> > makeZbosonsFromLeptons(const std::vector<const reco::Candidate *>& lm, const std::vector<const reco::Candidate *>& lp, int leptonCode, float mZ){
    
    phys::Boson<phys::Particle> Z0;
    phys::Boson<phys::Particle> Z1;
    float minMDiff=99999.;
    if (leptonCode == 4) {
      for (int k=0; k<2; ++k) {
	for (int j=0; j<2; ++j) {
	  float mDiff = fabs((lp[k]->p4() + lm[j]->p4()).mass() - mZ);
	  if ( mDiff < minMDiff ) {
	    minMDiff=mDiff;   
	    
	    Z0.setDaughter(0,phys::Particle(lp[k]->p4(), phys::Particle::computeCharge(lp[k]->pdgId()), lp[k]->pdgId()));
	    Z0.setDaughter(1,phys::Particle(lm[j]->p4(), phys::Particle::computeCharge(lm[j]->pdgId()), lm[j]->pdgId()));

	    Z1.setDaughter(0,phys::Particle(lp[(k+1)%2]->p4(), phys::Particle::computeCharge(lp[(k+1)%2]->pdgId()), lp[(k+1)%2]->pdgId()));
	    Z1.setDaughter(1,phys::Particle(lm[(j+1)%2]->p4(), phys::Particle::computeCharge(lm[(j+1)%2]->pdgId()), lm[(j+1)%2]->pdgId()));
	  }      
	} 	
      }
    }
    else { 
      for (int z=0; z<2; ++z) {
	if ( fabs(lp[z]->pdgId()) == fabs(lm[0]->pdgId()) ) { 
	  
	  Z0.setDaughter(0,phys::Particle(lp[z]->p4(), phys::Particle::computeCharge(lp[z]->pdgId()), lp[z]->pdgId()));
	  Z0.setDaughter(1,phys::Particle(lm[0]->p4(), phys::Particle::computeCharge(lm[0]->pdgId()), lm[0]->pdgId()));

	  Z1.setDaughter(0,phys::Particle(lp[(z+1)%2]->p4(), phys::Particle::computeCharge(lp[(z+1)%2]->pdgId()), lp[(z+1)%2]->pdgId()));
	  Z1.setDaughter(1,phys::Particle(lm[1]->p4(), phys::Particle::computeCharge(lm[1]->pdgId()), lm[1]->pdgId()));
	}
      }	
    }

    Z0.setId(23);
    Z1.setId(23);
    
    return std::make_pair(Z0,Z1);
  }
