//-----------FUNCTION: definition of the two ZZ bosons from leptons-------

#include "VVXAnalysis/Producers/interface/SignalDefinitionUtilities.h"

  std::pair<Boson*,Boson*> makeZbosonsFromLeptons(const std::vector<const reco::Candidate *>& l, const std::vector<const reco::Candidate *>& lm, const std::vector<const reco::Candidate *>& lp, int leptonCode, float mZ){
    
    Boson *Z0 = new Boson();
    Boson *Z1 = new Boson();
    
    float minMDiff=99999.;
    if (leptonCode == 4) {
      
      for (int k=0; k<2; ++k) {
	for (int j=0; j<2; ++j) {
	  float mDiff = fabs((lp[k]->p4() + lm[j]->p4()).mass() - mZ);
	  if ( mDiff < minMDiff ) {
	    minMDiff=mDiff;   
	    
	    Z0->Setdaughter1(lp[k]->p4());            
	    Z0->Setdaughter2(lm[j]->p4());
	    
	    Z1->Setdaughter1(lp[(k+1)%2]->p4());            
	    Z1->Setdaughter2(lm[(j+1)%2]->p4());
	    
	    if ( fabs(l[0]->pdgId()) == 11 ) {
	      Z0->SetdaughtersId(1); //u
	      Z1->SetdaughtersId(1); //u
	    }
	    if ( fabs(l[0]->pdgId()) == 13 ) {
	      Z0->SetdaughtersId(2); //e   
	      Z1->SetdaughtersId(2); //e   
	    }	
	  }      
	} 	
      }
    }
    else { 
      
      for (int z=0; z<2; ++z) {
	if ( fabs(lp[z]->pdgId()) == fabs(lm[0]->pdgId()) ) { 
	  
	  Z0->Setdaughter1(lp[z]->p4());
	  Z0->Setdaughter2(lm[0]->p4());	  
	  
	  Z1->Setdaughter1(lp[(z+1)%2]->p4());
	  Z1->Setdaughter2(lm[1]->p4());
	  
	  if ( fabs(lm[0]->pdgId()) == 11 ) {
	    Z0->SetdaughtersId(1); //u
	    Z1->SetdaughtersId(2); //e
	  }
	  if ( fabs(lm[0]->pdgId()) == 13 ) {
	    Z0->SetdaughtersId(2); //e
	    Z1->SetdaughtersId(1); //u
	  }
	  
	}
      }	
    }
    
    Z0->SetbosonId(23);
    Z1->SetbosonId(23);
    
    
    return std::make_pair(Z0,Z1);
    
  }
