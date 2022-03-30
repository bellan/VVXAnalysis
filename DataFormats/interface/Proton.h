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

    Proton(bool multiRP=0, double xi=0., double vx=0., double vy=0., double thetaX=0., double thetaY=0., double time=0., double Ebeam=6500)
      : xi_(xi)
      , vx_(vx)
      , vy_(vy)
      , thetaX_(thetaX)
      , thetaY_(thetaY)
      , time_(time)
      , E_((1-xi)*Ebeam)
      , Ebeam_(Ebeam)
      , efficiencySF_(1.)
      , efficiencySFUnc_(0.)
      , ismultiRP_(multiRP){
    }
    
    
    /// Destructor
    virtual ~Proton(){}
    
    // Operations
    Double_t xi() const {return xi_;}
    Double_t vx() const {return vx_;}
    Double_t vy() const {return vy_;}
    Double_t thetaX() const {return thetaX_;}
    Double_t thetaY() const {return thetaY_;}
    Double_t time() const {return time_;}
    
    Double_t Ebeam() const {return Ebeam_;}
    
    Double_t xiError() const {return xiError_;}
    Double_t vxError() const {return vxError_;}
    Double_t vyError() const {return vyError_;}
    Double_t thetaXError() const {return thetaXError_;}
    Double_t thetaYError() const {return thetaYError_;}
    Double_t timeError() const {return timeError_;}
    
    TLorentzVector p4()   const {return TLorentzVector(E_*thetaX_,E_*thetaY_,E_*sqrt(1-thetaX_*thetaX_-thetaY_*thetaY_),E_);}
    Double_t eta()        const {return this->p4().Eta();}
    Double_t rapidity()   const {return this->p4().Rapidity();}
    Double_t phi()        const {return this->p4().Phi();}
    Double_t p()          const {return this->p4().P();}
    Double_t e()          const {return this->p4().E();}
    
    Bool_t ismultiRP()    const {return ismultiRP_;}
    
    bool appropriatexi() const {if(xi_>0.04) return 1;
                                else return 0;}

    void setEfficenySFUnc(float effSfUnc ) {efficiencySFUnc_ = effSfUnc;}
    
    Double_t efficiencySF()  const {return efficiencySF_;}
    Double_t efficiencySFUnc()  const {return efficiencySFUnc_;}
    Bool_t   passFullSel() const {return true;}
    

  protected:
 
    //kinematic variables after scattering
    Double_t xi_;
    Double_t vx_;
    Double_t vy_;
    Double_t thetaX_;
    Double_t thetaY_;
    Double_t time_;

    Double_t xiError_;
    Double_t vxError_;
    Double_t vyError_;
    Double_t thetaXError_;
    Double_t thetaYError_;
    Double_t timeError_;
    
    Double_t E_;
    Double_t Ebeam_;
    
    Double_t efficiencySF_;
    Double_t efficiencySFUnc_;
    
    Bool_t valid_fit_;
    
    Bool_t ismultiRP_;
 
    
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
