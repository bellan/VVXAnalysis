#ifndef VVXAnalysis_DataFormats_Proton_H
#define VVXAnalysis_DataFormats_Proton_H

/** \class Proton
 *  No description available.
 *
 *  $Date: 2022/03/21 13:37:31 $
 *  $Revision: ??? $
 *  \author G. Marozzo - UNITO <giovanni.marozzo@edu.unito.it>
 */

#include <TObject.h>
#include <TLorentzVector.h> 

#include <iostream>
#include <cmath>

class TreePlanter;

namespace phys {
  
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
  
  class Proton: public TObject {
    
    friend class ::TreePlanter;
    
  public:
 
    /// Constructor

    Proton(double xi = 0., double x0 = 0., double y0 = 0., double thetax0 = 0., double thetay0 = 0., double Ebeam=6500)
      : xi_(xi)
      , x0_(x0)
      , y0_(y0)
      , thetax0_(thetax0)
      , thetay0_(thetay0)
      , Ebeam_(Ebeam)
      , E_((1-xi)*Ebeam)
      , efficiencySF_(1.)
      , efficiencySFUnc_(0.)
      , fakeRateSF_(1.){
    }
    
    
    /// Destructor
    virtual ~Proton(){};
    
    // Operations
    Double_t xi() const {return xi_;}
    Double_t x0() const {return x0_;}
    Double_t y0() const {return y0_;}
    Double_t thetax0() const {return thetax0_;}
    Double_t thetay0() const {return thetay0_;}
    Double_t Ebeam() const {return Ebeam_;}
    
    TLorentzVector p4() const {return TLorentzVector(E_*thetax0_,E_*thetay0_,E_*sqrt(1-thetax0_*thetax0_-thetay0_*thetay0_),E_);}
    Double_t eta()        const {return this->p4().Eta();}
    Double_t rapidity()   const {return this->p4().Rapidity();}
    Double_t phi()        const {return this->p4().Phi();}
    Double_t p()          const {return this->p4().P();}
    Double_t e()          const {return this->p4().E();}
    
    
    
    bool isValid() const {if(xi_>0.04) return 1;
                          else return 0;};

    void setEfficenySFUnc(float effSfUnc ) {efficiencySFUnc_ = effSfUnc;}
    
    Double_t efficiencySF()  const {return efficiencySF_;}
    Double_t efficiencySFUnc()  const {return efficiencySFUnc_;}
    Double_t fakeRateSF()    const {return fakeRateSF_;}
    Double_t fakeRateSFUnc() const {return fakeRateSFUnc_;} 
    Double_t fakeRateSFVar() const {return fakeRateSFUnc()*fakeRateSFUnc();}
    Bool_t   passFullSel() const {return true;}
    

  protected:
 
    //kinematic variables after scattering
    Double_t xi_;
    Double_t E_;
    Double_t x0_;
    Double_t y0_;
    Double_t thetax0_;
    Double_t thetay0_;
    Double_t Ebeam_;
    
    Double_t efficiencySF_;
    Double_t efficiencySFUnc_;
    Double_t fakeRateSF_;
    Double_t fakeRateSFUnc_;
    
  public:

    //to be changed
    bool operator==(const Proton& p1) const{
      return  (*this).p4().Px() - p1.p4().Px() < 0.01 &&
      (*this).p4().Py() - p1.p4().Py() < 0.01 &&
      (*this).p4().Pz() - p1.p4().Pz() < 0.01 &&
      (*this).p4().E()  - p1.p4().E()  < 0.01;
    }


    bool operator!=(Proton p1) const{
      return !((*this) == p1);
    }


    
  private:
    ClassDef(Proton, 1) //     
      };
  
}

#endif
