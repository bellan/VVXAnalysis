//-----------FUNCTION: definition of the two ZZ bosons from leptons-------

#include "VVXAnalysis/Producers/interface/SignalDefinitionUtilities.h"

  std::pair<Boson<Particle> ,Boson<Particle> > makeZbosonsFromLeptons(const std::vector<const reco::Candidate *>& lm, const std::vector<const reco::Candidate *>& lp, int leptonCode, float mZ){
    
    Boson<Particle> Z0;
    Boson<Particle> Z1;
    
    float minMDiff=99999.;
    if (leptonCode == 4) {
      
      for (int k=0; k<2; ++k) {
	for (int j=0; j<2; ++j) {
	  float mDiff = fabs((lp[k]->p4() + lm[j]->p4()).mass() - mZ);
	  if ( mDiff < minMDiff ) {
	    minMDiff=mDiff;   
	    
	    Z0.Setdaughter1(lp[k]->p4());            
	    Z0.Setdaughter2(lm[j]->p4());
	    
	    Z1.Setdaughter1(lp[(k+1)%2]->p4());            
	    Z1.Setdaughter2(lm[(j+1)%2]->p4());
	    
	    if ( fabs(lp[0]->pdgId()) == 11 ) {
	      Z0.Setdaughter1Id(11);  //e
	      Z0.Setdaughter2Id(-11); //e
	      Z1.Setdaughter1Id(11);  //e
	      Z1.Setdaughter2Id(-11); //e
	    }
	    if ( fabs(lp[0]->pdgId()) == 13 ) {
	      Z0.Setdaughter1Id(13);  //u   
	      Z0.Setdaughter2Id(-13); //u   
	      Z1.Setdaughter1Id(13);  //u   
	      Z1.Setdaughter2Id(-13); //u   
	    }	
	  }      
	} 	
      }
    }
    else { 
      
      for (int z=0; z<2; ++z) {
	if ( fabs(lp[z]->pdgId()) == fabs(lm[0]->pdgId()) ) { 
	  
	  Z0.Setdaughter1(lp[z]->p4());
	  Z0.Setdaughter2(lm[0]->p4());	  
	  
	  Z1.Setdaughter1(lp[(z+1)%2]->p4());
	  Z1.Setdaughter2(lm[1]->p4());
	  
	  if ( fabs(lm[0]->pdgId()) == 11 ) {
	    Z0.Setdaughter1Id(11);
	    Z0.Setdaughter2Id(-11);
	    Z1.Setdaughter1Id(13); 
	    Z1.Setdaughter2Id(-13);
	  }
	  if ( fabs(lm[0]->pdgId()) == 13 ) {
	    Z0.Setdaughter1Id(13);
	    Z0.Setdaughter2Id(-13);
	    Z1.Setdaughter1Id(11); 
	    Z1.Setdaughter2Id(-11);
	  }
	  
	}
      }	
    }
    
    Z0.setId(23);
    Z1.setId(23);
    
    
    return std::make_pair(Z0,Z1);
    
  }
