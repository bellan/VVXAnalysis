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
        
    
    /// Destructor
    virtual ~Boson(){};
    
    void setDaughter(int i, const P& d){
      if(i == 0) daughter0_ = d;
      if(i == 1) daughter1_ = d;
      else abort();
    }

    // Operations
    P daughter(int i){
      if(i == 0) return daughter0_;
      if(i == 1) return daughter1_;
      else abort();
    }

    
  protected:
    
  private:
    P daughter0_;
    P daughter1_;


    ClassDef(Boson, 1) //
  };
}

#endif
