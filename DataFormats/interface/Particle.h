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
    Particle(const TLorentzVector& mom = TLorentzVector(0.,0.,0.,0.), float q =0, int i = 0)
      : p4_(mom)
      , charge_(q)
      , id_(i)
      , motherId_(-99)
      , efficiencySF_(1.)
      , fakeRateSF_(1.)
      , genStatusFlags_(-99)
      {}

      Particle(const LorentzVector& l, float q =0, int i = 0)
	:p4_(convert(l))
	, charge_(q)
	, id_(i)
	, motherId_(-99)
        , efficiencySF_(1.)
        , fakeRateSF_(1.)
        , genStatusFlags_(-99)
      {
	// Correct Id for PF charge change
	if (fabs(id_==13)) id_= fabs(i)*(-1)*q ;
}


      /* Particle(const LorentzVector& l, float q =0, const int &&i = 13) */
      /* 	:p4_(convert(l)) */
      /* 	, charge_(q) */
      /* //, id_(i) */
      /* 	, motherId_(-99) */
      /*   , efficiencySF_(1.) */
      /*   , fakeRateSF_(1.) */
      /*   , genStatusFlags_(-99) */
      /* { */

      /* 	id_= fabs(i)*(-1)*q */
      /* 	  } */

	
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
    
    void setMotherId(int pid) {motherId_ = pid;}

    void setP4(const TLorentzVector& pi){p4_=pi;}
    
    int motherId() const {return motherId_;}

    bool isValid() const {return id_ != 0 && p() > 0;}

    virtual Double_t efficiencySF()  const {return efficiencySF_;}
    virtual Double_t fakeRateSF()    const {return fakeRateSF_;}
    virtual Double_t fakeRateSFUnc() const {return fakeRateSFUnc_;}
    virtual Double_t fakeRateSFVar() const {return fakeRateSFUnc()*fakeRateSFUnc();}

    Bool_t   passFullSel() const {return true;}
 
    // Gen info, in case they are meaningfull
    std::bitset<15> genStatusFlags() const {return genStatusFlags_;}
    void setGenStatusBit(GenStatusBit bit, int val = 1) {genStatusFlags_.set(bit, val);}
    

  protected:
    TLorentzVector p4_;
    Float_t charge_;
    Int_t id_;    
    Int_t motherId_;
    Double_t efficiencySF_;
    Double_t fakeRateSF_;
    Double_t fakeRateSFUnc_;
    std::bitset<15> genStatusFlags_;

  public:
    friend std::ostream&  operator<<(std::ostream& os, const Particle& obj){
      
      os << "ID = " << obj.id() << " p = (" << obj.p4().X() << "," << obj.p4().Py() << "," << obj.p4().Pz() << "," << obj.e() << "), pT = " << obj.pt() << " eta = " << (obj.p() != 0 ? obj.eta() : 0.) << " phi = " << obj.phi() << " mass = " << obj.mass();    
      // write obj to stream
      return os;
    }
    
  private:
    ClassDef(Particle, 1) //     
  };

}

#endif

