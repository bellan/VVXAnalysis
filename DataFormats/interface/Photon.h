#ifndef VVXAnalysis_DataFormats_Photon_H
#define VVXAnalysis_DataFormats_Photon_H

/** \class Photon
 *  No description available.
 *
 *  $Date: 2021/09/14 12:34:31 $
 *  $Revision: 1.5 $
 *  \author A. Mecca - UNITO <alberto.mecca@cern.ch>
 */

#include "RecoParticle.h"

class TreePlanter;

namespace phys {

  class Photon: public RecoParticle {

    friend class ::TreePlanter;

  public:
    enum IdWp {None, VeryLoose, Loose, Medium, Tight};  // Working points for ID
    
    enum class IDcut : UInt_t { // https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2#Applying_Individual_Cuts_of_a_Se
      MinPt  = 0x1,
      SCEta  = 0x2,
      HoverE = 0x4,
      sieie  = 0x8,
      chIso  = 0x10,
      neIso  = 0x20,
      phIso  = 0x40
    };

    enum class MVAwp {
      wp80,
      wp90
    };

    static constexpr UInt_t ID_BITMASK =
      static_cast<UInt_t>(IDcut::HoverE) |
      static_cast<UInt_t>(IDcut::sieie ) |
      static_cast<UInt_t>(IDcut::chIso ) |
      static_cast<UInt_t>(IDcut::neIso ) |
      static_cast<UInt_t>(IDcut::phIso ) ;

    static constexpr float TRANSITION_BARREL_ENDCAP = 1.479;
    
    // Constructor
    Photon(const TLorentzVector& pin = TLorentzVector(0.,0.,0.,0.), float q = 0, int id = 22):
        RecoParticle(pin, q, id)
	  {};
    
    virtual ~Photon(){};
    
    // Getters
    bool cutBasedIDTight()  const { return cutBasedIDTight_ ; }
    bool cutBasedIDMedium() const { return cutBasedIDMedium_; }
    bool cutBasedIDLoose()  const { return cutBasedIDLoose_ ; }
    
    float MVAvalue() const { return MVAvalue_; }
    
    float chargedIsolation()      const { return chargedIsolation_      ; } 
    float neutralHadronIsolation()const { return neutralHadronIsolation_; }
    float photonIsolation()       const { return photonIsolation_       ; }
    
    float sigmaIetaIeta() const { return sigmaIetaIeta_; }
    float HoverE()        const { return HoverE_       ; }
    
    bool passElectronVeto() const { return passElectronVeto_; }
    bool hasPixelSeed()     const { return hasPixelSeed_    ; }
    
    float seedEnergy()      const { return seedEnergy_; }
    float scEta() const { return scEta_; }
    
    float puppiChargedHadronIso() const { return puppiChargedHadronIso_; }
    float puppiNeutralHadronIso() const { return puppiNeutralHadronIso_; }
    float puppiPhotonIso()        const { return puppiPhotonIso_       ; }

    // Legacy methods for energy scale/smear systematics
    Float_t energyScaleUp()         const { return scale_total_up_*e(); };
    Float_t energyScaleDown()       const { return scale_total_dn_*e(); };
    Float_t energySigmaUp()         const { return sigma_total_up_*e(); };
    Float_t energySigmaDown()       const { return sigma_total_dn_*e(); };

  private:
    
    // ---- photon ID's ----
    bool cutBasedIDTight_;
    bool cutBasedIDMedium_;
    bool cutBasedIDLoose_;
    unsigned int cutIDbitsLoose_ ;  // An integer in which each bit represents a cut;
    unsigned int cutIDbitsMedium_;  // actually only the last 7 are used (see the enum IDcut).
    unsigned int cutIDbitsTight_ ;  // Let's trust ROOT's compression to save on the disk usage.
    
    float MVAvalue_;  // PhotonMVAEstimatorRunIIFall17v2Values
    
    // ---- Isolations ----
    float chargedIsolation_;
    float neutralHadronIsolation_;
    float photonIsolation_;
    
    // ---- Variables used in the ID ----
    float sigmaIetaIeta_;
    float HoverE_;
    
    // ---- conversion veto ----
    bool passElectronVeto_;
    bool hasPixelSeed_;
    
    // ---- input variables for regression energy corrections ----
    float seedEnergy_;
    float scEta_;  // superCluster Eta
    
    // PUPPI isolations
    float puppiChargedHadronIso_;
    float puppiNeutralHadronIso_;
    float puppiPhotonIso_;
    
    // Energy corrections
    /* float energyScaleUp_; // energy with the ecal energy scale shifted 1 sigma up (adding gain/stat/syst in quadrature) */
    /* float energyScaleDown_; // energy with the ecal energy scale shifted 1 sigma down (adding gain/stat/syst in quadrature)  */
    /* float energySigmaUp_; // energy with the ecal energy smearing value shifted 1 sigma up */
    /* float energySigmaDown_; // energy with the ecal energy smearing value shifted 1 sigma down */
    
  public:
    // Utils
    bool cutBasedID(IdWp wp) const{
      switch(wp){
      case VeryLoose:
	return (cutIDbitsLoose_ & static_cast<UInt_t>(IDcut::HoverE))
	  &&   (cutIDbitsLoose_ & static_cast<UInt_t>(IDcut::neIso ))
	  &&   (cutIDbitsLoose_ & static_cast<UInt_t>(IDcut::phIso ));
      case Loose:
	return cutBasedIDLoose() ; break;
      case Medium:
	return cutBasedIDMedium(); break;
      case Tight:
	return cutBasedIDTight() ; break;
      default:
	return true;
      }
    }

    UInt_t nCutsPass(IdWp wp) const{
      return std::bitset<32>(ID_BITMASK & getCutIdBits(wp)).count();
    }

    bool cutBasedID(IdWp wp, IDcut cut) const{
	return getCutIdBits(wp) & static_cast<UInt_t>(cut);  // bitwise AND
    }

    UInt_t getCutIdBits(IdWp wp) const{
      switch(wp){
      case VeryLoose:
	return cutIDbitsLoose_ | (static_cast<UInt_t>(IDcut::chIso) | static_cast<UInt_t>(IDcut::sieie)) ; break;
      case Loose:
	return cutIDbitsLoose_ ; break;
      case Medium:
	return cutIDbitsMedium_; break;
      case Tight:
	return cutIDbitsLoose_ ; break;
      default:
	return ID_BITMASK; // equivalent to return true
      }
    }

    bool passMVA(MVAwp wp) const{
      // See https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_X/RecoEgamma/PhotonIdentification/python/Identification/mvaPhotonID_Fall17_94X_V2_cff.py
      switch(wp){
      case MVAwp::wp90:
	return MVAvalue_ > (isBarrel() ? -0.02 : -0.26);
      case MVAwp::wp80:
	return MVAvalue_ > (isBarrel() ? +0.42 : +0.14);
      default:
	std::cerr << "WARN: unknown MVAwp in Photon::passMVA()\n";
	return true;
      }
    }

    inline bool isBarrel() const { return fabs(eta()) < TRANSITION_BARREL_ENDCAP; }
    
    bool passHoverE(IdWp wp) const{
      if(isBarrel())
	switch(wp){
	case IdWp::Tight:  return HoverE() < 0.02148;
	case IdWp::Medium: return HoverE() < 0.02197;
	case IdWp::Loose:  return HoverE() < 0.04596;
	default: return true;
	}
      else
	switch(wp){
	case IdWp::Tight:  return HoverE() < 0.0321;
	case IdWp::Medium: return HoverE() < 0.0326;
	case IdWp::Loose:  return HoverE() < 0.0590;
	default: return true;
	}
    }

    bool passSigmaiEtaiEta(IdWp wp) const{
      if(isBarrel())
	switch(wp){
	case IdWp::Tight:  return sigmaIetaIeta() < 0.00996;
	case IdWp::Medium: return sigmaIetaIeta() < 0.01015;
	case IdWp::Loose:  return sigmaIetaIeta() < 0.0106 ;
	default: return true;
	}
      else
	switch(wp){
	case IdWp::Tight:  return sigmaIetaIeta() < 0.0271;
	case IdWp::Medium: return sigmaIetaIeta() < 0.0272;
	case IdWp::Loose:  return sigmaIetaIeta() < 0.0272;
	default: return true;
	}
    }

    bool passChargedIsolation(IdWp wp) const{
      if(isBarrel())
	switch(wp){
	case IdWp::Tight:  return chargedIsolation() < 0.65 ;
	case IdWp::Medium: return chargedIsolation() < 1.141;
	case IdWp::Loose:  return chargedIsolation() < 1.694;
	default: return true;
	}
      else
	switch(wp){
	case IdWp::Tight:  return chargedIsolation() < 0.517;
	case IdWp::Medium: return chargedIsolation() < 1.051;
	case IdWp::Loose:  return chargedIsolation() < 2.089;
	default: return true;
	}
    }

    bool passNeutralIsolation(IdWp wp) const{
      if(isBarrel())
	switch(wp){
	case IdWp::Tight:  return neutralHadronIsolation() < 0.317  + 0.01512 *pt() + 2.259e-05 *pt()*pt();
	case IdWp::Medium: return neutralHadronIsolation() < 1.189  + 0.01512 *pt() + 2.259e-05 *pt()*pt();
	case IdWp::Loose:  return neutralHadronIsolation() < 24.032 + 0.01512 *pt() + 2.259e-05 *pt()*pt();
	default: return true;
	}
      else
	switch(wp){
	case IdWp::Tight:  return neutralHadronIsolation() < 2.716  + 0.0117 *pt() + 2.3e-05 *pt()*pt();
	case IdWp::Medium: return neutralHadronIsolation() < 2.718  + 0.0117 *pt() + 2.3e-05 *pt()*pt();
	case IdWp::Loose:  return neutralHadronIsolation() < 19.722 + 0.0117 *pt() + 2.3e-05 *pt()*pt();
	default: return true;
	}
    }
    
    bool passPhotonIsolation(IdWp wp) const{
      if(isBarrel())
	switch(wp){
	case IdWp::Tight:  return photonIsolation() < 2.044 + 0.004017 *pt();
	case IdWp::Medium: return photonIsolation() < 2.08  + 0.004017 *pt();
	case IdWp::Loose:  return photonIsolation() < 2.876 + 0.004017 *pt();
	default: return true;
	}
      else
	switch(wp){
	case IdWp::Tight:  return photonIsolation() < 3.032 + 0.0037 *pt();
	case IdWp::Medium: return photonIsolation() < 3.867 + 0.0037 *pt();
	case IdWp::Loose:  return photonIsolation() < 4.162 + 0.0037 *pt();
	default: return true;
	}
    }

    ClassDef(Photon, 2)
  };
}


#endif

