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
    enum class IdWp  : UInt_t { Veto=0     , Loose, Medium, Tight };  // Working points for cut-based ID
    enum class IsoWp : UInt_t { VeryLoose=0, Loose, Medium, Tight, VeryTight, VeryVeryTight };  // Isolations
    
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
      , sigmaIetaIeta_(-1.)
      , HoverE_(-1.)
      , EoverP_(-1.)
      , missingHits_(-1)
      , passConversionVeto_(false)
      , mvaValue_(-1.)
      , pogID_(0)
      , isoPF_(0)

      //, nearestjet(TLorentzVector(0.,0.,0.,0.), 0)
      {}
    
    /// Destructor
    virtual ~Lepton();
    
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
    
    Bool_t passPogID()             const {return PogID_;}
    
    Double_t sigmaIetaIeta()       const {return sigmaIetaIeta_;}
    Double_t HoverE()              const {return HoverE_;}
    Double_t EoverP()              const {return EoverP_;}
    Int_t missingHits()            const {return missingHits_;}
    Bool_t passConversionVeto()    const {return passConversionVeto_;}
    
    Double_t mvaValue()            const {return mvaValue_;}
    Bool_t passPogID(IdWp  wp)     const {return pogID_.test(static_cast<UInt_t>(wp));}
    Bool_t passIsoPF(IsoWp wp)     const {return isoPF_.test(static_cast<UInt_t>(wp));}
    std::bitset<4> pogID()         const {return pogID_;}  // For test purposes only; analyzers should use the above two functions
    std::bitset<6> isoPF()         const {return isoPF_;}

    // The fake rate is set to a value different from 1 even for true leptons.
    void setFakeRateSF(const std::pair<double,double> & sf) {
      fakeRateSF_    =  sf.first; 
      fakeRateSFUnc_ =  sf.second;
    }
    
    Double_t fakeRateSF()        const {return passFullSel() ? 1. : fakeRateSF_;}
    Double_t fakeRateSFUnc()     const {return passFullSel() ? 0. : fakeRateSFUnc_;}
    
    Double_t Em1_Pm1()           const {return (1/EoverP_-1)/p4_.P();}

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
    
    Bool_t PogID_;
    
    // Electron-specific variables
    float sigmaIetaIeta_;
    float HoverE_;
    float EoverP_;
    Int_t missingHits_;
    Bool_t passConversionVeto_;

    float mvaValue_;
    std::bitset<4> pogID_;
    std::bitset<6> isoPF_; // Only for muons

    ClassDef(Lepton, 2) //
  };
}

#endif

