#ifndef ZZWAnalysis_DataFormats_Boson_H
#define ZZWAnalysis_DataFormats_Boson_H

/** \class Boson
 *  No description available.
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.3 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include "Particle.h"

namespace phys {

  template <typename P>
  class Boson: public Particle {

    friend class ::TreePlanter;
  public:
    /// Constructor
    Boson(const TLorentzVector& p = TLorentzVector(0.,0.,0.,0.), int pid = 0)
      : Particle(p,0,pid)
      , indexFSR_(-1)
      , hasGoodDaughters_(false){}


    Boson(const P& daughter0, const P& daughter1, int pid = 0)
      : Particle(daughter0.p4()+daughter1.p4(), 0, pid)
      , daughter0_(daughter0)
      , daughter1_(daughter1)
      , indexFSR_(-1)
      , hasGoodDaughters_(false)
      {}
    
      Boson(const Boson<P>& vb)
	: Particle(vb.daughter(0).p4() + vb.daughter(1).p4(), vb.charge(), vb.id())
	, daughter0_(vb.daughter(0))
	, daughter1_(vb.daughter(1))
	, indexFSR_(vb.daughterWithFSR())
	, fsrPhoton_(vb.fsrPhoton())
	, hasGoodDaughters_(vb.hasGoodDaughters()){}

    
    template<typename T>
      Boson<T> clone() const {
      return Boson<T>(daughter0_,daughter1_,id_);
    }

    
    /// Destructor
    virtual ~Boson(){};
    
    void setDaughter(int i, const P& d){
      if(i == 0) daughter0_ = d;
      else if(i == 1) daughter1_ = d;
      else { std::cout << "*** Boson's daughter not found! ***" << " Cannot set it. " << i << std::endl; abort(); }
      p4_ = daughter0_.p4() + daughter1_.p4();
    }

    // Operations
    P daughter(int i) const{
      if(i == 0) return daughter0_;
      else if(i == 1) return daughter1_;
      else { std::cout << "*** Boson's daughter not found! ***" << " " << i << std::endl; abort();}
    }

    void addFSR(int daughter_index, const Particle &photon){
      indexFSR_ = daughter_index;
      fsrPhoton_ = photon;
      p4_ = p4_ + fsrPhoton_.p4();
    }

    

    Particle fsrPhoton() const {return fsrPhoton_;}

    int daughterWithFSR() const {return indexFSR_;}
   
    // the daughters pass the quality criteria
    bool hasGoodDaughters() const {return hasGoodDaughters_;}
    

  protected:
    
  private:
    P daughter0_;
    P daughter1_;

    Int_t indexFSR_;     // daughter with FSR, -1 if no one radiated
    Particle fsrPhoton_;

    Bool_t hasGoodDaughters_;

    ClassDef(Boson, 1) //
  };
}

#endif
