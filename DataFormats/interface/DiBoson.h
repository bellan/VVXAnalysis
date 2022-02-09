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
      , passFullSel_(false)
      {}
      
    DiBoson(const Boson<P1>& vb1, const Boson<P2>& vb2)
      : Particle(vb1.p4()+vb2.p4(),0,0)
      , daughter0_(vb1)
      , daughter1_(vb2)
      , passFullSel_(false)
      {
	
	for(unsigned int i = 0; i < 2; ++i){
	  id_ += abs(daughter0_.daughter(i).id()) + abs(daughter1_.daughter(i).id()) + daughter0_.daughter(i).id() + daughter1_.daughter(i).id();
	}
	efficiencySF_  = -1;
	efficiencySFUnc_  = -1;
	fakeRateSF_    = -1;
	fakeRateSFUnc_ = -1;
      }
    
    template<typename T1, typename T2>
      DiBoson<T1,T2> clone() const {
      DiBoson<T1,T2> newdiboson(daughter0_.template clone<T1>(), daughter1_.template clone<T2>());
      newdiboson.setPassFullSel(passFullSel_);
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

    Boson<P1> *firstPtr()  {return &daughter0_;}
    Boson<P2> *secondPtr() {return &daughter1_;}


    // True if pass all requirements on di-boson quantities and on
    // its daughters and grand daughters and...
    bool passFullSelection() const {return passFullSel_;}
    void setPassFullSel(Bool_t  fs) {passFullSel_ = fs;}

    int numberOfGoodGrandDaughters() const {
      if(!isValid()) return 0;
      return daughter0_.numberOfGoodDaughters() + daughter1_.numberOfGoodDaughters();
    }

    int numberOfBadGrandDaughters() const {
      if(!isValid()) return 0;
      return daughter0_.numberOfBadDaughters() + daughter1_.numberOfBadDaughters();
    }
    
    double fakeRateSF() const {
      if(id_ == 0) return 1.; // To be checked
      double ifakeRateSF = daughter0_.fakeRateSF() * daughter1_.fakeRateSF();
      return numberOfBadGrandDaughters() == 2 ? -1*ifakeRateSF : ifakeRateSF;
    }

    double efficiencySF() const{return daughter0_.efficiencySF() * daughter1_.efficiencySF();}
    
    double efficiencySFUnc() const {
      return daughter0_.muEffSFUnc()+daughter1_.muEffSFUnc()+daughter0_.eleEffSFUnc()+daughter1_.eleEffSFUnc();
    }
    
    double muEffSFUnc() const {
      return daughter0_.muEffSFUnc()+daughter1_.muEffSFUnc();
    }
    double eleEffSFUnc() const {
      return daughter0_.eleEffSFUnc()+daughter1_.eleEffSFUnc();
    }
    
    
    double fakeRateSFUnc() const{
      return sqrt(pow(daughter0_.fakeRateSF()*daughter1_.fakeRateSFUnc(),2) +
    		  pow(daughter1_.fakeRateSF()*daughter0_.fakeRateSFUnc(),2));
    }
    

  private:
  
    Boson<P1> daughter0_;
    Boson<P2> daughter1_;
    
    Bool_t  passFullSel_;


    ClassDef(DiBoson, 1) //
  };
}
#endif
