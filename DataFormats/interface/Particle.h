#ifndef VVXAnalysis_DataFormats_Particle_H
#define VVXAnalysis_DataFormats_Particle_H

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
#include <bitset>

#include "GenStatusBit.h"

#include <iostream>
#include <cmath>

class TreePlanter;

namespace phys {
  
  
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
  
  class Particle: public TObject {
    
    friend class ::TreePlanter;
    
  public:
    static TLorentzVector convert(const LorentzVector& l)  {return TLorentzVector(l.Px(),l.Py(),l.Pz(),l.E());}
    static LorentzVector  convert(const TLorentzVector& l) {return LorentzVector(l.Px(),l.Py(),l.Pz(),l.E());}
    
    /// Constructor
    Particle(const TLorentzVector& mom = TLorentzVector(0.,0.,0.,0.), float q = 0, int id = 0, std::bitset<15> flags = 0)
      : p4_(mom)
      , charge_(q)
      , id_(id)
    {
      if (abs(id_)==13) id_= copysign(id, (-1)*q) ;
    }
    
    Particle(const LorentzVector& l, float q = 0, int i = 0, std::bitset<15> flags = 0) : Particle(convert(l), q, i, flags) {}
    
    /// Destructor
    virtual ~Particle(){};
    
    // Operations
    TLorentzVector p4() const {return p4_;}
    int id()            const {return id_;}
    float charge()      const {return charge_;}
    double pt()         const {return p4_.Pt();}
    double eta()        const {return p4_.Eta();}
    double rapidity()   const {return p4_.Rapidity();}
    double phi()        const {return p4_.Phi();}
    double p()          const {return p4_.P();}
    double e()          const {return p4_.E();}
    double mass()       const {return p4_.M();} 
    
    // Method that tries to infer the charge of the particle starting from a pdgId in input
    // to be moved?
    static double computeCharge(int pdgId) {
      double charge = 0;
      if(abs(pdgId) == 1       || abs(pdgId) == 3  || abs(pdgId) == 5) // d, s or b
	charge = -1*copysign(1/3., pdgId);
      else if(abs(pdgId) == 2  || abs(pdgId) == 4  || abs(pdgId) == 6) // u, c or t
	charge = copysign(2/3., pdgId);
      else if(abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15) // e, mu or tau
	charge = -1*copysign(1, pdgId);
      else if(abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16 ||                   // ve, vmu, vtau
	      abs(pdgId) == 21 || abs(pdgId) == 22 || abs(pdgId) == 23 || abs(pdgId) == 25) // gluon, gamma, Z, H
	charge = 0;
      else if(abs(pdgId) == 24) // W
	charge = copysign(1, pdgId);
      
      return charge;
    }
    
    void setId(int pid) {id_ = pid; charge_ = computeCharge(pid);}
    
    void setP4(const TLorentzVector& pi){p4_=pi;}
    
    virtual bool isValid() const {return id_ != 0 && p() > 0;}

  protected:
    TLorentzVector p4_;
    Float_t charge_;
    Int_t id_;    
    
  public:

    virtual bool operator==(const Particle& p1) const{
      return  (*this).id() == p1.id()         &&
      (*this).p4().Px() - p1.p4().Px() < 0.01 &&
      (*this).p4().Py() - p1.p4().Py() < 0.01 &&
      (*this).p4().Pz() - p1.p4().Pz() < 0.01 &&
      (*this).p4().E()  - p1.p4().E()  < 0.01;
    }


    virtual bool operator!=(const Particle& p1) const{
      return !((*this) == p1);
    }

    friend std::ostream&  operator<<(std::ostream& os, const Particle& obj){
      
      os << "ID = " << obj.id() << " p = (" << obj.p4().X() << "," << obj.p4().Py() << "," << obj.p4().Pz() << "," << obj.e() << "), pT = " << obj.pt() << " eta = " << (obj.p() != 0 ? obj.eta() : 0.) << " phi = " << obj.phi() << " mass = " << obj.mass();    
      // write obj to stream
      return os;
    }

  private:
    ClassDef(Particle, 2)
      };
  
}

#endif

