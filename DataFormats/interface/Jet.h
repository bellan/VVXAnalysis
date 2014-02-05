#ifndef ZZWAnalysis_DataFormats_Jet_H
#define ZZWAnalysis_DataFormats_Jet_H

/** \class Jet
 *  No description available.
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UCSB <riccardo.bellan@cern.ch>
 */

#include "Particle.h"

class TreePlanter;

namespace phys {

  class Jet: public Particle {
    // Friends
    friend class ::TreePlanter;
    
  public:
    
    /// Constructor
    Jet(const TLorentzVector& p = TLorentzVector(0.,0.,0.,0.), int q =0, int id = 0)
      : Particle(p,q, id)
      , nConstituents_(-1)
      , nCharged_(-1)
      , nNeutral_(-1)
      , neutralHadronEnergyFraction_(-9999.)
      , chargedHadronEnergyFraction_(-9999.)
      , chargedEmEnergyFraction_(-9999.)    
      , neutralEmEnergyFraction_(-9999.)    
      , muonEnergyFraction_(-9999.)
      , csvtagger_(-2)
      , girth_(-9999.)
      , girth_charged_(-9999.)
      , ptd_(-9999.)
      , rms_(-9999.)
      , beta_(-9999.)
      , jetArea_(-9999.)
      , secvtxMass_(-9999.)
      , Lxy_(-9999.)
      , LxyErr_(-9999.)
      , rawFactor_(-9999.)
      , uncOnFourVectorScale_(-9999.)
      , puMVAFull_(-9999.)
      , puMVASimple_(-9999.)
      , puCutBased_(-9999.)
      , pass_puMVAFull_loose_(false)
      , pass_pUMVAFull_medium_(false)
      , pass_pUMVAFull_tight_(false) 
      , pass_puMVASimple_loose_(false) 
      , pass_puMVASimple_medium_(false)
      , pass_puMVASimple_tight_(false) 
      , pass_puCutBased_loose_(false) 
      , pass_puCutBased_medium_(false)
      , pass_puCutBased_tight_(false) 
      , mcPartonFlavour_(-1)
    {}           
    
    /// Destructor
    virtual ~Jet(){};
    
    // Operations

    Int_t nConstituents() const {return nConstituents_;} 
    Int_t nCharged()      const {return nCharged_;}	  
    Int_t nNeutral()      const {return nNeutral_;}      

    Double_t neutralHadronEnergyFraction() const {return neutralHadronEnergyFraction_;}    
    Double_t chargedHadronEnergyFraction() const {return chargedHadronEnergyFraction_;} 
    Double_t chargedEmEnergyFraction()     const {return chargedEmEnergyFraction_;}     
    Double_t neutralEmEnergyFraction()     const {return neutralEmEnergyFraction_;}     
    Double_t muonEnergyFraction()          const {return muonEnergyFraction_;}          

    // B-tagging info
    Double_t csvtagger()     const {return csvtagger_;}         
    
    // Quark-Gluon discrimination variables
    Double_t girth()         const {return girth_;}
    Double_t girth_charged() const {return girth_charged_;}
    // sum pt^2 / (sum pt)^2
    Double_t ptd()           const {return ptd_;}
    // jet width
    Double_t rms()           const {return rms_;}
    // fraction of charged hadron sum pt from charged hadrons from the primary vertex
    Double_t beta()          const {return beta_;}
    // return the jet area 
    Double_t jetArea()       const {return jetArea_;}

    // return secondary vertex b-tagging information
    Double_t secvtxMass()    const {return secvtxMass_;}
    Double_t Lxy()           const {return Lxy_;}	 
    Double_t LxyErr()        const {return LxyErr_;}    
      
    // return a correction factor that can be applied to the jet energy or pT to bring it back to the uncorrected value
    Double_t rawFactor()     const {return rawFactor_;}
    
    // Uncertainty on four vector energy scale
    Double_t uncOnFourVectorScale() const {return  uncOnFourVectorScale_;}

    // PU ID:
    Double_t puMVAFull()     const {return puMVAFull_;}	      
    Double_t puMVASimple()   const {return puMVASimple_;}	      
    Double_t puCutBased()    const {return puCutBased_;}	      
					 				      
    Bool_t   pass_puMVAFull_loose()    const {return pass_puMVAFull_loose_;}   
    Bool_t   pass_pUMVAFull_medium()   const {return pass_pUMVAFull_medium_;}  
    Bool_t   pass_pUMVAFull_tight()    const {return pass_pUMVAFull_tight_;}   
					 				      
    Bool_t   pass_puMVASimple_loose()  const {return pass_puMVASimple_loose_;} 
    Bool_t   pass_puMVASimple_medium() const {return pass_puMVASimple_medium_;}
    Bool_t   pass_puMVASimple_tight()  const {return pass_puMVASimple_tight_;} 
					 	    		      
    Bool_t   pass_puCutBased_loose()   const {return pass_puCutBased_loose_;}  
    Bool_t   pass_puCutBased_medium()  const {return pass_puCutBased_medium_;} 
    Bool_t   pass_puCutBased_tight()   const {return pass_puCutBased_tight_;}  

    // return the matched MC parton flavour
    Int_t mcPartonFlavour() const {return mcPartonFlavour_;}

    
  protected:
    
  private:

    Int_t nConstituents_;
    Int_t nCharged_;
    Int_t nNeutral_;

    Double_t neutralHadronEnergyFraction_;    
    Double_t chargedHadronEnergyFraction_;
    Double_t chargedEmEnergyFraction_;    
    Double_t neutralEmEnergyFraction_;    
    Double_t muonEnergyFraction_;

    // B-tagging info
    Double_t csvtagger_;         
    
    // Quark-Gluon discrimination variables
    Double_t girth_;
    Double_t girth_charged_;
    // sum pt^2 / (sum pt)^2
    Double_t ptd_;
    // jet width
    Double_t rms_;
    // fraction of charged hadron sum pt from charged hadrons from the primary vertex
    Double_t beta_;
    // return the jet area 
    Double_t jetArea_;

    // return secondary vertex b-tagging information
    Double_t secvtxMass_;
    Double_t Lxy_;
    Double_t LxyErr_;
      
    // return a correction factor that can be applied to the jet energy or pT to bring it back to the uncorrected value
    Double_t rawFactor_;
    
    // Uncertainty on four vector energy scale
    Double_t uncOnFourVectorScale_;

    // PU ID:
    Double_t puMVAFull_;
    Double_t puMVASimple_;
    Double_t puCutBased_;

    Bool_t   pass_puMVAFull_loose_;
    Bool_t   pass_pUMVAFull_medium_;
    Bool_t   pass_pUMVAFull_tight_; 

    Bool_t   pass_puMVASimple_loose_; 
    Bool_t   pass_puMVASimple_medium_;
    Bool_t   pass_puMVASimple_tight_; 

    Bool_t   pass_puCutBased_loose_; 
    Bool_t   pass_puCutBased_medium_;
    Bool_t   pass_puCutBased_tight_; 

    // return the matched MC parton flavour
    Int_t mcPartonFlavour_;
    
    ClassDef(Jet, 1) //
  };
}

#endif

