#ifndef ZZWAnalysis_DataFormats_Particle_H
#define ZZWAnalysis_DataFormats_Particle_H

/** \class Particle
 *  No description available.
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.3 $
 *  \author R. Bellan - UCSB <riccardo.bellan@cern.ch>
 */

#include <TObject.h>
#include <TLorentzVector.h> 
#include "Math/GenVector/LorentzVector.h"

class TreePlanter;

namespace phys {

  
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

  class Particle: public TObject {
    
    friend class ::TreePlanter;

  public:
    static TLorentzVector convert(const LorentzVector& l) {return TLorentzVector(l.Px(),l.Py(),l.Pz(),l.E());}
    
    /// Constructor
    Particle(const TLorentzVector& p = TLorentzVector(0.,0.,0.,0.), int q =0, int i = 0)
      : p4_(p)
      , charge_(q)
      , id_(i){}

      Particle(const LorentzVector& l, int q =0, int i = 0)
	:p4_(convert(l))
	, charge_(q)
	, id_(i){}

	
    /// Destructor
    virtual ~Particle(){};
    
    // Operations
    TLorentzVector p4() const {return p4_;}
    int id()            const {return id_;}
    int charge()        const {return charge_;}
    double pt()         const {return p4_.Pt();}
    double eta()        const {return p4_.Eta();}
    double phi()        const {return p4_.Phi();}
    
 
   private:
    TLorentzVector p4_;
    Int_t charge_;
    
    
  protected:
    Int_t id_;    


    //Particle genParticle_;

  private:
    ClassDef(Particle, 1) //
  };
}

#endif

