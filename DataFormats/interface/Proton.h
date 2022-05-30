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

#include "Particle.h"

#include <iostream>
#include <cmath>

class TreePlanter;

namespace phys {
  
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
  
  class Proton: public TObject {
    
    friend class ::TreePlanter;
    
  public:
 
    /// Constructor

    Proton(bool multiRP=0, double xi=0., double vx=0., double vy=0., double thetaX=0., double thetaY=0., double Ebeam=6500)
      : xi_(xi)
      , vx_(vx)
      , vy_(vy)
      , thetaX_(thetaX)
      , thetaY_(thetaY)
      , Ebeam_(Ebeam)
      , ismultiRP_(multiRP){
    }
    
    Proton(Particle genProton)
      : xi_(1-genProton.e()/6500){
    }
    
    Proton(double xi, bool LHCSector)
      : xi_(xi)
      , LHCSector_(LHCSector){
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
    Double_t E()     const {return 6500*(1-xi_);}
    
    Bool_t LHCSector() const {return LHCSector_;}
    
    Double_t xiError() const {return xiError_;}
    Double_t vxError() const {return vxError_;}
    Double_t vyError() const {return vyError_;}
    Double_t thetaXError() const {return thetaXError_;}
    Double_t thetaYError() const {return thetaYError_;}
    Double_t timeError() const {return timeError_;}
    
    TLorentzVector p4()   const {
    double pz;
    if(!LHCSector_) pz=this->E();
    else pz=(-1)*this->E();
    return TLorentzVector(pz*TMath::Sin(thetaX_),pz*TMath::Sin(thetaY_),pz,this->E());}
    
    Double_t pt()         const {return this->p4().Pt();}
    Double_t eta()        const {return this->p4().Eta();}
    Double_t rapidity()   const {return this->p4().Rapidity();}
    Double_t phi()        const {return this->p4().Phi();}
    Double_t p()          const {return this->p4().P();}
    
    Bool_t ismultiRP()    const {return ismultiRP_;}
    Bool_t valid_fit()    const {return valid_fit_;}
    
    bool appropriatexi() const {if(xi_>0.05) return 1;
                                else return 0;}
                                

  protected:
 
    //kinematic variables after scattering
    Double_t xi_;
    Double_t vx_;
    Double_t vy_;
    Double_t thetaX_;
    Double_t thetaY_;
    Double_t time_;
    
    Bool_t LHCSector_;       //  0 -> sector 45, 1 -> sector 56

    Double_t xiError_;
    Double_t vxError_;
    Double_t vyError_;
    Double_t thetaXError_;
    Double_t thetaYError_;
    Double_t timeError_;
   
    Double_t Ebeam_;
    
    Bool_t valid_fit_;
    
    Bool_t ismultiRP_;


    
  private:
    ClassDef(Proton, 4) //     
      };
  
}

#endif
