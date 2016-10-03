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
      , indexFSR_(0){}

    
  Boson(const P& daughter0, const P& daughter1, int pid = 0)
    : Particle(daughter0.p4()+daughter1.p4(), daughter0.charge()+daughter1.charge(), pid)
      , daughter0_(daughter0)
      , daughter1_(daughter1)
      , indexFSR_(0){

      init();

    }
    
  Boson(const Boson<P>& vb)
    : Particle(vb.daughter(0).p4() + vb.daughter(1).p4(), vb.charge(), vb.id())
      , daughter0_(vb.daughter(0))
      , daughter1_(vb.daughter(1))
      , indexFSR_(vb.daughtersWithFSR())
      , fsrPhoton0_(vb.fsrPhoton(0))
      , fsrPhoton1_(vb.fsrPhoton(1)){

      init();
    }
    
    
    template<typename T>
      Boson<T> clone() const {
      Boson<T> newboson(daughter0_,daughter1_,id_);

      std::bitset<2> index = std::bitset<2>(indexFSR_);
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
      std::bitset<2> index = std::bitset<2>(indexFSR_);

      if(index.test(daughter_index)){
	std::cout<<"This lepton: " << daughter_index << " aleady has a photon associated to it"<<std::endl;
	abort();
      }

      index.set(daughter_index);
      indexFSR_ = index.to_ulong();

      if(daughter_index == 0)      fsrPhoton0_ = photon;
    
      else if(daughter_index == 1) fsrPhoton1_ = photon;
      else { 
	std::cout << "*** (FSR) Boson's daughter not found! ***" << " " << daughter_index << std::endl; 
	abort();}
      
      p4_ = p4_ + photon.p4();
    }

    

    Particle fsrPhoton(int d) const {
      if(d == 0) return fsrPhoton0_;
      else if(d == 1) return fsrPhoton1_;
      else { std::cout << "*** (FSR) Boson's daughter not found! ***" << " " << d << std::endl; abort();}
    }

    int daughtersWithFSR() const {return indexFSR_;}
   
    // Number of good daughters
    int numberOfGoodDaughters() const {return int(daughter0_.passFullSel()) + int(daughter1_.passFullSel());}

    double fakeRateSF()    const {return daughter0_.fakeRateSF() * daughter1_.fakeRateSF();}
    double fakeRateSFUnc() const {return sqrt(pow(daughter0_.fakeRateSF()*daughter1_.fakeRateSFUnc(),2) +
    						  pow(daughter1_.fakeRateSF()*daughter0_.fakeRateSFUnc(),2));}
    
    double efficiencySF() const {return daughter0_.efficiencySF() * daughter1_.efficiencySF();}
    
    double efficiencySFUnc() const {
      double effSF0    = daughter0_.efficiencySF();
      double effSF0Unc = daughter0_.efficiencySFUnc();
      double effSF1    = daughter1_.efficiencySF();
      double effSF1Unc = daughter1_.efficiencySFUnc();      
      
      return sqrt((effSF0Unc*effSF0Unc)/(effSF0*effSF0)+(effSF1Unc*effSF1Unc)/(effSF1*effSF1));
    }
    
  protected:
  private:
    P daughter0_;
    P daughter1_;

    Int_t indexFSR_;     // daughter with FSR, -1 if no one radiated
    Particle fsrPhoton0_;
    Particle fsrPhoton1_;

    void init(){

      efficiencySF_  = -1;
      fakeRateSF_    = -1;
      fakeRateSFUnc_ = -1; 
      
      charge_ = daughter0_.charge() + daughter1_.charge();

      if(indexFSR_ >=0)  {
	p4_ = p4_ + fsrPhoton0_.p4() + fsrPhoton1_.p4();}
      
    }

    ClassDef(Boson, 1) //
  };
}

#endif
