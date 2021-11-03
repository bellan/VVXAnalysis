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
    typedef std::pair<std::string, Bool_t> IdPair;
    
    // Constructor
    Photon(const TLorentzVector& pin = TLorentzVector(0.,0.,0.,0.), float q =0, int id = 0):
        Particle(pin, q, id)
	  {};
    
    virtual ~Photon(){};
    
    // Getters
    bool cutBasedIDTight()  const { return cutBasedIDTight_ ; }
    bool cutBasedIDMedium() const { return cutBasedIDMedium_; }
    bool cutBasedIDLoose()  const { return cutBasedIDLoose_ ; }
    
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
    
    ClassDef(Photon, 1)
  };
}
#endif

