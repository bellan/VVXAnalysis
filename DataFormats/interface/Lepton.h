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
    Lepton(const TLorentzVector& pin = TLorentzVector(0.,0.,0.,0.), float q =0, int pid = 0)
      : Particle(pin,q, pid)
      , dxy_(-9999.)               
      , dz_(-9999.)                
      , sip_(-9999.)
       , pfCombRelIso_(-9999.)
      , pfCombRelIsoFSRCorr_(-9999.)
       , scEta_(-9999.)
      , matchHLT_(false)
      , isGood_(false)
      , isInCracks_(false)

      //, nearestjet(TLorentzVector(0.,0.,0.,0.), 0)
      {}
    
    /// Destructor
    virtual ~Lepton(){};
    
    // Operations
    Double_t dxy()                 const {return dxy_;}            
    Double_t dz()                  const {return dz_;}             
    Double_t sip()                 const {return sip_;}	      
    Double_t pfCombRelIso()        const {return pfCombRelIso_;}   
    Double_t pfCombRelIsoFSRCorr() const {return pfCombRelIsoFSRCorr_;}   
    Double_t scEta()               const {return scEta_;}             
    Bool_t   matchHLT()            const {return matchHLT_;}
    Bool_t   isGood()              const {return isGood_;}

    Bool_t  isInCracks()           const {return isInCracks_;}

    Bool_t   passFullSelNoFSRCorr()const {return isGood_ && pfCombRelIso_ < (abs(id_) == 13 ? 0.35 : 0.35);}
    //Bool_t   passFullSel()         const {return isGood_ && pfCombRelIsoFSRCorr_ < (abs(id_) == 13 ? 0.35 : 0.35);}
    Bool_t   passFullSel()         const {return isGood_;} // In ZZAnalysis:Run2Legacy the iso is included in the ID

    // The fake rate is set to a value different from 1 even for true leptons.
    void setFakeRateSF(const std::pair<double,double> & sf) {
      fakeRateSF_    =  sf.first; 
      fakeRateSFUnc_ =  sf.second;
    }
    
    Double_t fakeRateSF()        const {return passFullSel() ? 1. : fakeRateSF_;}
    Double_t fakeRateSFUnc()     const {return passFullSel() ? 1. : fakeRateSFUnc_;}

  protected:
    

  private:
    Double_t dxy_;               
    Double_t dz_;                
    Double_t sip_;
    Double_t pfCombRelIso_;
    Double_t pfCombRelIsoFSRCorr_;
    Double_t scEta_;

    Bool_t matchHLT_;
    Bool_t isGood_;

    Bool_t isInCracks_;

    ClassDef(Lepton, 1) //
  };
}

#endif

