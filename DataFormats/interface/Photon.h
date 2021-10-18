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
    bool cutBasedIDTight()  { return cutBasedIDTight_ ; }
    bool cutBasedIDMedium() { return cutBasedIDMedium_; }
    bool cutBasedIDLoose()  { return cutBasedIDLoose_ ; }
    
    float chargedIsolation()      { return chargedIsolation_      ; } 
    float neutralHadronIsolation(){ return neutralHadronIsolation_; }
    float photonIsolation()       { return photonIsolation_       ; }
    
    bool passElectronVeto() { return passElectronVeto_; }
    bool hasPixelSeed()     { return hasPixelSeed_    ; }
    
    float seedEnergy() { return seedEnergy_; }
    /* float eMax() { return eMax_; } */
    /* float e2nd() { return e2nd_; } */
    /* float e3x3() { return e3x3_; } */
    /* float eTop() { return eTop_; } */
    /* float eBottom() { return eBottom_; } */
    /* float eLeft() { return eLeft_; } */
    /* float eRight() { return eRight_; } */
    /* float see() { return see_; } */
    /* float spp() { return spp_; } */
    /* float sep() { return sep_; } */
    /* float maxDR() { return maxDR_; } */
    /* float maxDRDPhi() { return maxDRDPhi_; } */
    /* float maxDRDEta() { return maxDRDEta_; } */
    /* float maxDRRawEnergy() { return maxDRRawEnergy_; } */
    /* float subClusRawE1() { return subClusRawE1_; } */
    /* float subClusRawE2() { return subClusRawE2_; } */
    /* float subClusRawE3() { return subClusRawE3_; } */
    /* float subClusDPhi1() { return subClusDPhi1_; } */
    /* float subClusDPhi2() { return subClusDPhi2_; } */
    /* float subClusDPhi3() { return subClusDPhi3_; } */
    /* float subClusDEta1() { return subClusDEta1_; } */
    /* float subClusDEta2() { return subClusDEta2_; } */
    /* float subClusDEta3() { return subClusDEta3_; } */
    float cryEta() { return cryEta_; }
    float cryPhi() { return cryPhi_; }
    float iEta() { return iEta_; }
    float iPhi() { return iPhi_; }
    float puppiChargedHadronIso() { return puppiChargedHadronIso_; }
    float puppiNeutralHadronIso() { return puppiNeutralHadronIso_; }
    float puppiPhotonIso() { return puppiPhotonIso_; }
    
    float ecalEnergyPreCorr() { return ecalEnergyPreCorr_; };
    float ecalEnergyErrPreCorr() { return ecalEnergyErrPreCorr_; };
    float ecalEnergyPostCorr() { return ecalEnergyPostCorr_; };
    float ecalEnergyErrPostCorr() { return ecalEnergyErrPostCorr_; };
    /* float ecalTrkEnergyPreCorr() { return ecalTrkEnergyPreCorr_; }; */
    /* float ecalTrkEnergyErrPreCorr() { return ecalTrkEnergyErrPreCorr_; }; */
    /* float ecalTrkEnergyPostCorr() { return ecalTrkEnergyPostCorr_; }; */
    /* float ecalTrkEnergyErrPostCorr() { return ecalTrkEnergyErrPostCorr_; }; */
    float energyScaleValue() { return energyScaleValue_; };
    float energySigmaValue() { return energySigmaValue_; };
    float energySmearNrSigma() { return energySmearNrSigma_; };
    float energyScaleUp() { return energyScaleUp_; };
    float energyScaleDown() { return energyScaleDown_; };
    float energyScaleStatUp() { return energyScaleStatUp_; };
    float energyScaleStatDown() { return energyScaleStatDown_; };
    float energyScaleSystUp() { return energyScaleSystUp_; };
    float energyScaleSystDown() { return energyScaleSystDown_; };
    float energyScaleGainUp() { return energyScaleGainUp_; };
    float energyScaleGainDown() { return energyScaleGainDown_; };
    float energyScaleEtUp() { return energyScaleEtUp_; };
    float energyScaleEtDown() { return energyScaleEtDown_; };
    float energySigmaUp() { return energySigmaUp_; };
    float energySigmaDown() { return energySigmaDown_; };
    float energySigmaPhiUp() { return energySigmaPhiUp_; };
    float energySigmaPhiDown() { return energySigmaPhiDown_; };
    float energySigmaRhoUp() { return energySigmaRhoUp_; };
    float energySigmaRhoDown() { return energySigmaRhoDown_; };

    
  private:
    
    // ---- photon ID's ----
    bool cutBasedIDTight_;
    bool cutBasedIDMedium_;
    bool cutBasedIDLoose_;
    
    // ---- Isolations ----
    float chargedIsolation_;
    float neutralHadronIsolation_;
    float photonIsolation_;
    
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
    
    float cryEta_;
    float cryPhi_;
    float iEta_;
    float iPhi_;
    
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

