#ifndef VVXAnalysis_DataFormats_Lepton_H
#define VVXAnalysis_DataFormats_Lepton_H

/** \class Lepton
 *  No description available.
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.5 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include "Particle.h"
#include "Jet.h"

class TreePlanter;

namespace phys {

  class Lepton: public Particle {

    friend class ::TreePlanter;

  public:
    
    /// Constructor
    Lepton(const TLorentzVector& p = TLorentzVector(0.,0.,0.,0.), float q =0, int pid = 0)
      : Particle(p,q, pid)
      , dxy_(-9999.)               
      , dz_(-9999.)                
      , sip_(-9999.)
      , combRelIso_(-9999.)
      , pfChargedHadIso_(-9999.)
      , pfNeutralHadIso_(-9999.)
      , pfPhotonIso_(-9999.)
      , pfCombRelIso_(-9999.)
      , pfCombRelIsoFSRCorr_(-9999.)
      , rho_(-9999.) 
      , isPF_(false)
      , matchHLT_(false)
      , isGood_(false)
      //, nearestjet(TLorentzVector(0.,0.,0.,0.), 0)
      {}
    
    /// Destructor
    virtual ~Lepton(){};
    
    // Operations
    Double_t dxy()                 const {return dxy_;}            
    Double_t dz()                  const {return dz_;}             
    Double_t sip()                 const {return sip_;}	      
    Double_t combRelIso()          const {return combRelIso_;}     
    Double_t pfChargedHadIso()     const {return pfChargedHadIso_;}
    Double_t pfNeutralHadIso()     const {return pfNeutralHadIso_;}
    Double_t pfPhotonIso()         const {return pfPhotonIso_;}    
    Double_t pfCombRelIso()        const {return pfCombRelIso_;}   
    Double_t pfCombRelIsoFSRCorr() const {return pfCombRelIsoFSRCorr_;}   
    Double_t rho()                 const {return rho_;}             
    Bool_t   isPF()                const {return isPF_;}
    Bool_t   matchHLT()            const {return matchHLT_;}
    Bool_t   isGood()              const {return isGood_;}

  protected:
    

  private:
    Double_t dxy_;               
    Double_t dz_;                
    Double_t sip_;
    Double_t combRelIso_;
    Double_t pfChargedHadIso_;
    Double_t pfNeutralHadIso_;
    Double_t pfPhotonIso_;
    Double_t pfCombRelIso_;
    Double_t pfCombRelIsoFSRCorr_;
    Double_t rho_;

    Bool_t isPF_;
    Bool_t matchHLT_;
    Bool_t isGood_;

    //Jet nearestjet;    

    ClassDef(Lepton, 1) //
  };
}

#endif

