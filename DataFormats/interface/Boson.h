#ifndef ZZWAnalysis_DataFormats_Boson_H
#define ZZWAnalysis_DataFormats_Boson_H

/** \class Boson
 *  No description available.
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.3 $
 *  \author R. Bellan - UCSB <riccardo.bellan@cern.ch>
 */

#include "Particle.h"

namespace phys {

  template <typename P>
  class Boson: public Particle {

  public:
    /// Constructor
    Boson(const TLorentzVector& p = TLorentzVector(0.,0.,0.,0.), int id = 0)
      : Particle(p,0,id){}


    Boson(const P& daughter0, const P& daughter1, int id = 0)
      : Particle(daughter0.p4()+daughter1.p4(), 0, id)
      , daughter0_(daughter0)
      , daughter1_(daughter1)
      {}
        
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

    
  protected:
    
  private:
    P daughter0_;
    P daughter1_;


    ClassDef(Boson, 1) //
  };
}

#endif
