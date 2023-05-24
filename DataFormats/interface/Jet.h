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
    
    //    enum JERVariations{central,up,down}; /DEL
    
    struct JetScores {
      Double_t TvsQCD = -2;
      Double_t WvsQCD = -2;
      Double_t ZvsQCD = -2;
      Double_t ZbbvsQCD = -2;
      Double_t HbbvsQCD = -2;
      Double_t H4qvsQCD = -2;
    };
    
    struct DeepFlavourScores {
      Double_t probb   = -2.;
      Double_t probc   = -2.;
      Double_t probg   = -2.;
      Double_t problepb= -2.;
      Double_t probbb  = -2.;
      Double_t probuds = -2.;
    };
    
    /// Constructor
    Jet(const TLorentzVector& p = TLorentzVector(0.,0.,0.,0.), float q =0, int pid = 0)
      : Particle(p, q, pid)
      , chargedMultiplicity_(-1)
      , neutralMultiplicity_(-1)
      , neutralHadronEnergyFraction_(-9999.)
      , chargedHadronEnergyFraction_(-9999.)
      , chargedEmEnergyFraction_(-9999.)
      , neutralEmEnergyFraction_(-9999.)
      , muonEnergyFraction_(-9999.)
      , electronEnergyFraction_(-9999.)
      , photonEnergyFraction_(-9999.)
      , csvtagger_(-2)
      , deepAK8_()
      , deepAK8_MD_()
      , particleNet_()
      //, particleNet_MD_()
      , girth_(-9999.)
      , girth_charged_(-9999.)
      , ptd_(-9999.)
      , jetArea_(-9999.)
      , secvtxMass_(-9999.)
      , qgLikelihood_ (-99.)
      , rawFactor_(-9999.)
      , jecUnc_(-9999.)
      , mcPartonFlavour_(-1)
      , passLooseId_(false)
      , fullPuId_(-1)
      , tau1_(-999)
      , tau2_(-999)
      , tau3_(-999)
      , corrPrunedMass_(-999)
      , prunedMass_(-999)
      , softDropMass_(-999)
      , puppiTau1_(-999)
      , puppiTau2_(-999)
      , puppiTau3_(-999)
      , puppiMass_(-999)

    {}           
    
    /// Destructor
    virtual ~Jet(){};
    
    // Operations

    Int_t chargedMultiplicity() const { return chargedMultiplicity_; }
    Int_t neutralMultiplicity() const { return neutralMultiplicity_; }

    float neutralHadronEnergyFraction() const { return neutralHadronEnergyFraction_; }
    float chargedHadronEnergyFraction() const { return chargedHadronEnergyFraction_; }
    float chargedEmEnergyFraction()     const { return chargedEmEnergyFraction_; }
    float neutralEmEnergyFraction()     const { return neutralEmEnergyFraction_; }
    float muonEnergyFraction()          const { return muonEnergyFraction_; }
    float electronEnergyFraction()      const { return electronEnergyFraction_; }
    float photonEnergyFraction()        const { return photonEnergyFraction_; }

    // B-tagging info
    Double_t csvtagger()     const {return csvtagger_;}         
    
    // DeepAK8 and ParticleNet score getters
    const JetScores& deepAK8()        const { return deepAK8_;     }
    const JetScores& deepAK8_MD()     const { return deepAK8_MD_;  }
    const JetScores& particleNet()    const { return particleNet_; }
    // const JetScores& particleNet_MD() const { return particleNet_MD_; }
    const DeepFlavourScores& deepFlavour() const { return deepFlavour_; }

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

    // Pile-up full-id
    Bool_t fullPuId(int level = 1) const {
      if(level == 1)       return (fullPuId_ & (1 << 1));
      else if(level == 0)  return (fullPuId_ & (1 << 2));
      else if(level == 2)  return (fullPuId_ & (1 << 0));
      else { 
	std::cout<<"Error in Pu-id request. the level value "<<level<<" is uncorrect. Choose between 0 for loose, 1 for medium or 2 fot tight"<<std::endl;  
	abort();
      } 
    } 
    
    // JER: Values from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    Double_t ptJerUp() const {return  pt_jerup_;}
    Double_t ptJerDn() const {return  pt_jerdn_;}
    Double_t ptNoJer() const {return  pt_nojer_;}

    // return the matched MC parton flavour
    Int_t mcPartonFlavour() const {return mcPartonFlavour_;}

    bool passLooseJetID() const {
      return passLooseId_;
    }

    Double_t qgLikelihood() const {return qgLikelihood_;}

    bool  passPUID() const {return true;}

    // AK8 methods
    double tau1()           const {return tau1_;}
    double tau2()           const {return tau2_;}
    double tau3()           const {return tau3_;}
    double corrPrunedMass() const {return corrPrunedMass_;}
    double prunedMass()     const {return prunedMass_;}
    double softDropMass()   const {return softDropMass_;}
    double puppiTau1()      const {return puppiTau1_;}
    double puppiTau2()      const {return puppiTau2_;}
    double puppiTau3()      const {return puppiTau3_;}
    double puppiMass()      const {return puppiMass_;}
    
    inline double chosenAlgoMass() const {return softDropMass_;} //puppiMass_ + 11.85;}

  protected:
    
  private:
    Int_t chargedMultiplicity_;
    Int_t neutralMultiplicity_;
    Int_t pileUpId_;
    Int_t partonId_;
    Int_t partonFlavour_;
    Int_t hadronFlavour_;

    float neutralHadronEnergyFraction_;
    float chargedHadronEnergyFraction_;
    float chargedEmEnergyFraction_;
    float neutralEmEnergyFraction_;
    float muonEnergyFraction_;
    float electronEnergyFraction_;
    float photonEnergyFraction_;

    // B-tagging info
    Double_t csvtagger_;
    
    // DeepAK8 and ParticleNet scores
    JetScores deepAK8_, deepAK8_MD_, particleNet_; //, particleNet_MD_;
    
    // Deep Flavour
    DeepFlavourScores deepFlavour_;
    
    // Quark-Gluon discrimination variables
    Double_t girth_;
    Double_t girth_charged_;
    // sum pt^2 / (sum pt)^2
    Double_t ptd_;
    // jet width (area) 
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

    // loose ID
    bool passLooseId_;

    // full Pile up ID
    Int_t fullPuId_;
    
    // Jet MC resolution, for JER determination. 
    Double_t pt_nojer_;
    Double_t pt_jerup_;
    Double_t pt_jerdn_;


    // AK8
    Double_t tau1_;
    Double_t tau2_;
    Double_t tau3_;
    Double_t corrPrunedMass_;
    Double_t prunedMass_;
    Double_t softDropMass_;
    Double_t puppiTau1_;
    Double_t puppiTau2_;
    Double_t puppiTau3_;
    Double_t puppiMass_;


    ClassDef(Jet, 2)
  };
}

#endif

