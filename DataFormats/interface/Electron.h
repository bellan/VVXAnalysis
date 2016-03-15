#ifndef VVXAnalysis_DataFormats_Electron_H
#define VVXAnalysis_DataFormats_Electron_H

/** \class Electron
 *  No description available.
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.3 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include "Lepton.h"

namespace phys {

  class Electron: public Lepton {

    friend class ::TreePlanter;
  public:

    /// Constructor
    Electron(const TLorentzVector& p = TLorentzVector(0.,0.,0.,0.), float q =0, int pid = 0)
      : Lepton(p,q,pid)
      , energy_(-9999.)
      , phiWidth_(-9999.)
      , etaWidth_(-9999.)
      , BDT_(-9999.)
      , isBDT_(false)
      , missingHit_(-1)
      , nCrystals_(-1)
      {}
      
      Electron(const Lepton& lep)
	: Lepton(lep)
	, energy_(-9999.)
	, phiWidth_(-9999.)
	, etaWidth_(-9999.)
	, BDT_(-9999.)
	, isBDT_(false)
	, missingHit_(-1)
	, nCrystals_(-1)
	{
	  if(abs(id()) != 11) Particle::id_ = 111;
	}
    


    /// Destructor
    virtual ~Electron(){};
    
    // Operations
    Double_t energy()     const {return energy_;}	 
    Double_t phiWidth()   const {return phiWidth_;}	 
    Double_t etaWidth()   const {return etaWidth_;}	 
    Double_t BDT()        const {return BDT_;}  	 
    Bool_t   isBDT()      const {return isBDT_;}	 
    Int_t    missingHit() const {return missingHit_;}
    Int_t    nCrystals()  const {return nCrystals_;} 
    
    Bool_t passBDT() const {

      bool lowPt   = pt() > 7 && pt() < 10;
      bool highPt  = pt() >= 10;
      if(!lowPt && !highPt) return false;

      bool lowEta  = fabs(eta()) < 0.8;
      bool midEta  = fabs(eta()) >= 0.8 && fabs(eta()) < 1.479;
      bool highEta = fabs(eta()) >= 1.479 && fabs(eta()) < 2.5; 
      if(!lowEta && !midEta && ! highEta) return false;

      if(pfCombRelIso() > 0.4 || missingHit() > 1 || sip() >=  4) return false;

      if(lowPt){
	if(lowEta)  return BDT() > 0.47;
	if(midEta)  return BDT() > 0.004;
	if(highEta) return BDT() > 0.295;
      }
      if(highPt){
	if(lowEta)  return BDT() > -0.34;
	if(midEta)  return BDT() > -0.65;
	if(highEta) return BDT() >  0.60;
      }
      return false;
    }


  protected:
    
  private:
    Double_t energy_;
    Double_t phiWidth_;
    Double_t etaWidth_;
    Double_t BDT_;  
    Bool_t isBDT_;
    Int_t  missingHit_; 
    Int_t  nCrystals_;
    
    ClassDef(Electron, 1) //
  };
}

#endif

