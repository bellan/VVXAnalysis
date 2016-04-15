#ifndef VVXAnalysis_DataFormats_Jet_H
#define VVXAnalysis_DataFormats_Jet_H

/** \class Jet
 *  No description available.
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include "Particle.h"

class TreePlanter;

namespace phys {

  class Jet: public Particle {
    // Friends
    friend class ::TreePlanter;
    
  public:
    
    enum JERVariations{central,up,down};

    /// Constructor
    Jet(const TLorentzVector& p = TLorentzVector(0.,0.,0.,0.), float q =0, int pid = 0)
      : Particle(p, q, pid)
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
      , jetArea_(-9999.)
      , secvtxMass_(-9999.)
      , qgLikelihood_ (-99.)
      , rawFactor_(-9999.)
      , jecUnc_(-9999.)
      , mcPartonFlavour_(-1)
      , sigma_MC_pt_(-9999.)
      , sigma_MC_phi_(-9999.)
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
    // return the jet area 
    Double_t jetArea()       const {return jetArea_;}

    // return secondary vertex b-tagging information
    Double_t secvtxMass()    const {return secvtxMass_;}
      
    // return a correction factor that can be applied to the jet energy or pT to bring it back to the uncorrected value
    Double_t rawFactor()     const {return rawFactor_;}
    
    // Uncertainty on four vector energy scale
    Double_t jecUncertainty() const {return  jecUnc_;}

    // JER
    Double_t sigma_MC_pt()  const {return sigma_MC_pt_;}
    Double_t sigma_MC_phi()  const {return sigma_MC_phi_;}

    // Values from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    Double_t jer_c(JERVariations jervar)     const {
      switch(jervar){
      case(central): return jer_c_;
      case(up):      return jer_cup_;
      case(down):    return jer_cdown_;
      default:       return -9999.;
      }
    }

    Double_t jer_width(JERVariations jervar) const {return sqrt(pow(jer_c(jervar),2)-1)*sigma_MC_pt();}
    

    // return the matched MC parton flavour
    Int_t mcPartonFlavour() const {return mcPartonFlavour_;}

    bool    passLooseJetID() const {
      return 
	( fabs(eta()) <= 3 &&
	 (neutralHadronEnergyFraction_   < 0.99                        &&
	  neutralEmEnergyFraction_       < 0.99                        &&
	  nConstituents_                 > 1                           &&
	  (chargedHadronEnergyFraction_  > 0    || fabs(eta()) > 2.4)  &&
	  (nCharged_                     > 0    || fabs(eta()) > 2.4)  &&
	  (chargedEmEnergyFraction_      < 0.99 || fabs(eta()) > 2.4)))
	||
	(fabs(eta()) > 3 && neutralEmEnergyFraction_ < 0.90 && nNeutral_ > 10);
    }

    Double_t qgLikelihood() const {return qgLikelihood_;}


    bool  passPUID() const {return true;}

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
    // return the jet area 
    Double_t jetArea_;

    // return secondary vertex b-tagging information
    Double_t secvtxMass_;

    Double_t qgLikelihood_;
      
    // return a correction factor that can be applied to the jet energy or pT to bring it back to the uncorrected value
    Double_t rawFactor_;
    
    // Uncertainty on four vector energy scale
    Double_t jecUnc_;

    // return the matched MC parton flavour
    Int_t mcPartonFlavour_;
    
    // Jet MC resolution, for JER determination. 
    Double_t sigma_MC_pt_;
    Double_t sigma_MC_phi_;
    Double_t jer_c_;
    Double_t jer_cup_;
    Double_t jer_cdown_;

    ClassDef(Jet, 1) //
  };
}

#endif

