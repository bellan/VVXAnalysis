#ifndef VVXAnalysis_DataFormats_Photon_H
#define VVXAnalysis_DataFormats_Photon_H

/** \class Photon
 *  No description available.
 *
 *  $Date: 2021/09/14 12:34:31 $
 *  $Revision: 1.5 $
 *  \author A. Mecca - UNITO <alberto.mecca@cern.ch>
 */

#include "Particle.h"

class TreePlanter;

namespace phys {

  class Photon: public Particle {

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
    Photon(const TLorentzVector& pin = TLorentzVector(0.,0.,0.,0.), float q =0, int id = 0):
        Particle(pin, q, id)
	  {};
    
    virtual ~Photon();
    
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
    /* float eMax() const { return eMax_; } */
    /* float e2nd() const { return e2nd_; } */
    /* float e3x3() const { return e3x3_; } */
    /* float eTop() const { return eTop_; } */
    /* float eBottom() const { return eBottom_; } */
    /* float eLeft() const { return eLeft_; } */
    /* float eRight() const { return eRight_; } */
    /* float see() const { return see_; } */
    /* float spp() const { return spp_; } */
    /* float sep() const { return sep_; } */
    /* float maxDR() const { return maxDR_; } */
    /* float maxDRDPhi() const { return maxDRDPhi_; } */
    /* float maxDRDEta() const { return maxDRDEta_; } */
    /* float maxDRRawEnergy() const { return maxDRRawEnergy_; } */
    /* float subClusRawE1() const { return subClusRawE1_; } */
    /* float subClusRawE2() const { return subClusRawE2_; } */
    /* float subClusRawE3() const { return subClusRawE3_; } */
    /* float subClusDPhi1() const { return subClusDPhi1_; } */
    /* float subClusDPhi2() const { return subClusDPhi2_; } */
    /* float subClusDPhi3() const { return subClusDPhi3_; } */
    /* float subClusDEta1() const { return subClusDEta1_; } */
    /* float subClusDEta2() const { return subClusDEta2_; } */
    /* float subClusDEta3() const { return subClusDEta3_; } */
    /* float cryEta() const { return cryEta_; } */
    /* float cryPhi() const { return cryPhi_; } */
    /* float iEta()   const { return iEta_  ; } */
    /* float iPhi()   const { return iPhi_  ; } */
    float scEta() const { return scEta_; }
    
    float puppiChargedHadronIso() const { return puppiChargedHadronIso_; }
    float puppiNeutralHadronIso() const { return puppiNeutralHadronIso_; }
    float puppiPhotonIso()        const { return puppiPhotonIso_       ; }
    
    float ecalEnergyPreCorr()     const { return ecalEnergyPreCorr_; };
    float ecalEnergyErrPreCorr()  const { return ecalEnergyErrPreCorr_; };
    float ecalEnergyPostCorr()    const { return ecalEnergyPostCorr_; };
    float ecalEnergyErrPostCorr() const { return ecalEnergyErrPostCorr_; };
    /* float ecalTrkEnergyPreCorr() const { return ecalTrkEnergyPreCorr_; }; */
    /* float ecalTrkEnergyErrPreCorr() const { return ecalTrkEnergyErrPreCorr_; }; */
    /* float ecalTrkEnergyPostCorr() const { return ecalTrkEnergyPostCorr_; }; */
    /* float ecalTrkEnergyErrPostCorr() const { return ecalTrkEnergyErrPostCorr_; }; */
    float energyScaleValue()      const { return energyScaleValue_   ; };
    float energySigmaValue()      const { return energySigmaValue_   ; };
    float energySmearNrSigma()    const { return energySmearNrSigma_ ; };
    float energyScaleUp()         const { return energyScaleUp_      ; };
    float energyScaleDown()       const { return energyScaleDown_    ; };
    float energyScaleStatUp()     const { return energyScaleStatUp_  ; };
    float energyScaleStatDown()   const { return energyScaleStatDown_; };
    float energyScaleSystUp()     const { return energyScaleSystUp_  ; };
    float energyScaleSystDown()   const { return energyScaleSystDown_; };
    float energyScaleGainUp()     const { return energyScaleGainUp_  ; };
    float energyScaleGainDown()   const { return energyScaleGainDown_; };
    float energyScaleEtUp()       const { return energyScaleEtUp_    ; };
    float energyScaleEtDown()     const { return energyScaleEtDown_  ; };
    float energySigmaUp()         const { return energySigmaUp_      ; };
    float energySigmaDown()       const { return energySigmaDown_    ; };
    float energySigmaPhiUp()      const { return energySigmaPhiUp_   ; };
    float energySigmaPhiDown()    const { return energySigmaPhiDown_ ; };
    float energySigmaRhoUp()      const { return energySigmaRhoUp_   ; };
    float energySigmaRhoDown()    const { return energySigmaRhoDown_ ; };

    
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
    /* float eMax_; */
    /* float e2nd_; */
    /* float e3x3_; */
    /* float eTop_; */
    /* float eBottom_; */
    /* float eLeft_; */
    /* float eRight_; */
    
    /* float see_; */
    /* float spp_; */
    /* float sep_; */
    
    /* float maxDR_; */
    /* float maxDRDPhi_; */
    /* float maxDRDEta_; */
    /* float maxDRRawEnergy_; */
    
    /* float subClusRawE1_; */
    /* float subClusRawE2_; */
    /* float subClusRawE3_; */
    
    /* float subClusDPhi1_; */
    /* float subClusDPhi2_; */
    /* float subClusDPhi3_; */
    
    /* float subClusDEta1_; */
    /* float subClusDEta2_; */
    /* float subClusDEta3_; */
    
    /* float cryEta_; */
    /* float cryPhi_; */
    /* float iEta_; */
    /* float iPhi_; */
    float scEta_;  // superCluster Eta
    
    // PUPPI isolations
    float puppiChargedHadronIso_;
    float puppiNeutralHadronIso_;
    float puppiPhotonIso_;
    
    // Energy corrections
    float ecalEnergyPreCorr_; // ecalEnergy before scale & smearing corrections
    float ecalEnergyErrPreCorr_; // resolution estimate on the ecalEnergy before scale & smearing corrections
    float ecalEnergyPostCorr_; // ecalEnergy of electron after scale & smearing corrections
    float ecalEnergyErrPostCorr_; // resolution estimate on the ecalEnergy after scale & smearing corrections
    /* float ecalTrkEnergyPreCorr_; // ECAL-Trk combined electron energy before scale & smearing corrections */
    /* float ecalTrkEnergyErrPreCorr_; // resolution estimate of the ECAL-Trk combined electron energy before scale & smearing corrections */
    /* float ecalTrkEnergyPostCorr_; // ECAL-Trk combined electron energy after scale & smearing corrections */
    /* float ecalTrkEnergyErrPostCorr_; // resolution estimate of the ECAL-Trk combined electron energy after scale & smearing corrections */
    float energyScaleValue_; // value of the scale correction, MC ignores this value and takes 1
    float energySigmaValue_; // value of the resolution correction
    float energySmearNrSigma_; // a Gaussian random number to smear by (deterministic based on supercluster), data ignores this value and takes 0
    float energyScaleUp_; // energy with the ecal energy scale shifted 1 sigma up (adding gain/stat/syst in quadrature)
    float energyScaleDown_; // energy with the ecal energy scale shifted 1 sigma down (adding gain/stat/syst in quadrature) 
    float energyScaleStatUp_; // energy with the ecal energy scale shifted 1 sigma(stat) up
    float energyScaleStatDown_; // energy with the ecal energy scale shifted 1 sigma(stat) down
    float energyScaleSystUp_; // energy with the ecal energy scale shifted 1 sigma(syst) up
    float energyScaleSystDown_; // energy with the ecal energy scale shifted 1 sigma(syst) down
    float energyScaleGainUp_; // energy with the ecal energy scale shifted 1 sigma(gain) up
    float energyScaleGainDown_; // energy with the ecal energy scale shifted 1 sigma(gain) down
    float energyScaleEtUp_; // 2016 legacy only: adhoc error derived from closure vs Et
    float energyScaleEtDown_; // 2016 legacy only: adhoc error dervied from closure vs Et
    float energySigmaUp_; // energy with the ecal energy smearing value shifted 1 sigma up
    float energySigmaDown_; // energy with the ecal energy smearing value shifted 1 sigma down
    float energySigmaPhiUp_; // energy with the ecal energy smearing value shifted 1 sigma(phi) up
    float energySigmaPhiDown_; // energy with the ecal energy smearing value shifted 1 sigma(phi) down
    float energySigmaRhoUp_; // energy with the ecal energy smearing value shifted 1 sigma(rho) up
    float energySigmaRhoDown_; // energy with the ecal energy smearing value shifted 1 sigma(rho) down
    
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

    ClassDef(Photon, 1)
  };
}


#endif

