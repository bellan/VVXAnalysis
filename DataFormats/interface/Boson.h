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
#include <bitset>

namespace phys {

  template <typename P>
  class Boson: public Particle {

    friend class ::TreePlanter;
  public:
    /// Constructor
  Boson(const TLorentzVector& pi = TLorentzVector(0.,0.,0.,0.), int pid = 0)
    : Particle(pi,0,pid)
      , indexFSR_(-1)
      , hasGoodDaughters_(false){}

    
  Boson(const P& daughter0, const P& daughter1, int pid = 0)
    : Particle(daughter0.p4()+daughter1.p4(), daughter0.charge()+daughter1.charge(), pid)
      , daughter0_(daughter0)
      , daughter1_(daughter1)
      , indexFSR_(-1)
      , hasGoodDaughters_(false){

      init();

    }
    
  Boson(const Boson<P>& vb)
    : Particle(vb.daughter(0).p4() + vb.daughter(1).p4(), vb.charge(), vb.id())
      , daughter0_(vb.daughter(0))
      , daughter1_(vb.daughter(1))
      , indexFSR_(vb.daughterWithFSR())
      , fsrPhoton0_(vb.fsrPhoton(0))
      , fsrPhoton1_(vb.fsrPhoton(1))
      , hasGoodDaughters_(vb.hasGoodDaughters()){

      init();
    }
    
    
    template<typename T>
      Boson<T> clone() const {
      Boson<T> newboson(daughter0_,daughter1_,id_);
      std::bitset<4> index = std::bitset<4>(indexFSR_);

      if(index.test(0)) newboson.addFSR(0,fsrPhoton0_);
      if(index.test(1)) newboson.addFSR(1,fsrPhoton1_);

      return newboson;
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


    P *daughterPtr(int i) {
      if(i == 0) return &daughter0_;
      else if(i == 1) return &daughter1_;
      else { std::cout << "*** Boson's daughter not found! ***" << " " << i << std::endl; abort();}
    }


    void addFSR(int daughter_index, const Particle &photon){
      std::bitset<4> index = std::bitset<4>(indexFSR_);
      index.set(daughter_index);
      indexFSR_ = index.to_ulong();

      if(daughter_index == 0) fsrPhoton0_ = photon;
      else if(daughter_index == 1) fsrPhoton1_ = photon;
      else { std::cout << "*** (FSR) Boson's daughter not found! ***" << " " << daughter_index << std::endl; abort();}
      p4_ = p4_ + photon.p4();
    }

    

    Particle fsrPhoton(int d) const {
      if(d == 0) return fsrPhoton0_;
      else if(d == 1) return fsrPhoton1_;
      else { std::cout << "*** (FSR) Boson's daughter not found! ***" << " " << d << std::endl; abort();}
    }

    int daughterWithFSR() const {return indexFSR_;}
   
    // the daughters pass the quality criteria
    bool hasGoodDaughters() const {return hasGoodDaughters_;}
    
    // Number of good daughters
    int numberOfGoodDaughters() const {return int(daughter0_.passFullSel()) + int(daughter1_.passFullSel());}

    double fakeRateSF()    const {return daughter0_.fakeRateSF() * daughter1_.fakeRateSF();}
    double fakeRateSFUncHigh() const {return sqrt(pow(daughter0_.fakeRateSF()*daughter1_.fakeRateSFUncHigh(),2) +  
					      pow(daughter1_.fakeRateSF()*daughter0_.fakeRateSFUncHigh(),2));}

    double fakeRateSFUncLow() const {return sqrt(pow(daughter0_.fakeRateSF()*daughter1_.fakeRateSFUncLow(),2) +  
					      pow(daughter1_.fakeRateSF()*daughter0_.fakeRateSFUncLow(),2));}
    
    double efficiencySF() const {return daughter0_.efficiencySF() * daughter1_.efficiencySF();}
    

  protected:
    
  private:
    P daughter0_;
    P daughter1_;

    Int_t indexFSR_;     // daughter with FSR, -1 if no one radiated
    Particle fsrPhoton0_;
    Particle fsrPhoton1_;

    Bool_t hasGoodDaughters_;

    void init(){
      efficiencySF_  = -1;
      fakeRateSF_    = -1;
      fakeRateSFUncHigh_ = -1;
      fakeRateSFUncLow_ = -1;
      
      charge_ = daughter0_.charge() + daughter1_.charge();
      if(indexFSR_ >=0)  { p4_ = p4_ + fsrPhoton0_.p4() + fsrPhoton1_.p4();}
    }

    ClassDef(Boson, 1) //
  };
}

#endif
