#ifndef ZZWAnalysis_DataFormats_DiBoson_H
#define ZZWAnalysis_DataFormats_DiBoson_H

/** \class DiBoson
 *  
 *  Ids: 4mu = 4*13 = 52; 2e2mu = 2*13 + 2*11 = 48; 4e = 4*11 = 44 
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.3 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include "Particle.h"
#include "Boson.h"

namespace phys {
  
  template <typename P1, typename P2>
    class DiBoson: public Particle {
    
    friend class ::TreePlanter;
  public:
    /// Constructor
    DiBoson(): Particle()
      , triggerWord_(0)
      , regionWord_(0)
      , isBestCand_(false)
      , passFullSel_(false)
      , passTrigger_(false)
      {}
      
    DiBoson(const Boson<P1>& vb1, const Boson<P2>& vb2)
      : Particle(vb1.p4()+vb2.p4(),0,0)
      , daughter0_(vb1)
      , daughter1_(vb2)
      , triggerWord_(0)
      , regionWord_(0)
      , isBestCand_(false)
      , passFullSel_(false)
      , passTrigger_(false)
      {
	for(unsigned int i = 0; i < 2; ++i){
	  id_ += abs(daughter0_.daughter(i).id()) + abs(daughter1_.daughter(i).id()) + daughter0_.daughter(i).id() + daughter1_.daughter(i).id();
	}
	efficiencySF_ = daughter0_.efficiencySF() * daughter1_.efficiencySF();
	
	// compute the fake rate
	// The contribution from non prompt leptons arises from two sources: 3 Passes 1 Fail (3P1F) and 2 Passes 2 Fails (2P2F).
	// Accordingly to e. 16 in AN-2013/108 the number of background events in the search region can be written as two contributions
	// N_bkg_SR = (N_3P1F-N_ZZ_3P1F) * Sum(SFi) + N_2P2F * Sum(SFi*SFj - SFi - SFj)
	// where SFi = fi/(1-fi) is the fake rate scale factor written inside the phys::Lepton data format. 
	// This means that each event, depending on the lepton goodness contribution will get a different weight:
	// 4P0F (SR candidate) = 1
	// 3P1F = SFi
	// 2P2F = SFi*SFj - SFi - SFj
	// The code below shall reproduce this logic.

	// If both bosons are made of good leptons, the the fakeRateSF is 1. If we are in the case 3P1F, then
	// one of the two boson is made of good leptons, i.e., fakeRateSF = 1. Therefore it is safe to take the product of the two VB fake rate, without checking which one is the fake one.
	fakeRateSF_   = daughter0_.fakeRateSF() * daughter1_.fakeRateSF();
	// The case with two leptons failing the selection is a bit more complex, and the fakeRateSF needs to eb overwritten
	if(numberOfGoodGrandDaughters() == 2){  
	  // The formula for the good VB brings a 0, while for the fake one that is not true
	  fakeRateSF_   = daughter0_.daughter(0).fakeRateSF() * daughter0_.daughter(1).fakeRateSF() - daughter0_.daughter(0).fakeRateSF() - daughter0_.daughter(1).fakeRateSF() +
	                  daughter1_.daughter(0).fakeRateSF() * daughter1_.daughter(1).fakeRateSF() - daughter1_.daughter(0).fakeRateSF() - daughter1_.daughter(1).fakeRateSF();
	}
      }
    
    template<typename T1, typename T2>
      DiBoson<T1,T2> clone() const {
      DiBoson<T1,T2> newdiboson(daughter0_.clone<T1>(), daughter1_.clone<T2>());
      newdiboson.setTriggerWord(triggerWord_);
      newdiboson.setRegionWord (regionWord_ );
      newdiboson.setIsBestCand (isBestCand_ );
      newdiboson.setPassFullSel(passFullSel_);
      newdiboson.setPassTrigger(passTrigger_);
      return newdiboson;
    }

    template<typename T>
      Boson<T> daughter(int i) const{
      if(i == 0) return daughter0_;
      else if(i == 1) return daughter1_;
      else { std::cout << "*** DiBoson's daughter not found! ***" << " " << i << std::endl; abort();}
    }

    Boson<P1> first()  const {return daughter0_;}
    Boson<P2> second() const {return daughter1_;}


    // Best candidate in the Control/Search region
    bool isBestCandidate() const {return isBestCand_;}

    // True if pass all requirements on di-boson quantities and on
    // its daughters and grand daughters and...
    bool passFullSelection() const {return passFullSel_;}

    // Type of search/control region
    int region() const {return regionWord_;}

    // Triggers that have been passed
    short trigger() const {return triggerWord_;}
    
    // True if pass the trigger for a given final state
    bool passTrigger() const {return passTrigger_;}

    void setTriggerWord(Short_t tw) {triggerWord_ = tw;}
    void setRegionWord (Int_t   rw) {regionWord_  = rw;} 
    void setIsBestCand (Bool_t  bc) {isBestCand_  = bc;} 
    void setPassFullSel(Bool_t  fs) {passFullSel_ = fs;}
    void setPassTrigger(Bool_t  pt) {passTrigger_ = pt;}

    int numberOfGoodGrandDaughters() const {
      // Put a protection because right now are contemplated only cases where at least one boson is made of good leptons.
      if(daughter0_.numberOfGoodDaughters() < 2 && daughter1_.numberOfGoodDaughters() < 2) abort();
      return daughter0_.numberOfGoodDaughters() + daughter1_.numberOfGoodDaughters();
    }

  private:

    Boson<P1> daughter0_;
    Boson<P2> daughter1_;
    
    Short_t triggerWord_;
    Int_t   regionWord_;
    Bool_t  isBestCand_;
    Bool_t  passFullSel_;
    Bool_t  passTrigger_;

    ClassDef(DiBoson, 1) //
  };
}
#endif
