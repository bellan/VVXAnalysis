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

  public:
    
    /// Constructor
    Electron(const TLorentzVector& p = TLorentzVector(0.,0.,0.,0.), int q =0, int id = 0)
      : Lepton(p,q,id)
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
	  Particle::id_ = 11;
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
    

  protected:
    
  public:
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

