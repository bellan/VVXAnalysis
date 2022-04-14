#ifndef VVGammaAnalyzer_h
#define VVGammaAnalyzer_h

/** \class VVGammaAnalyzer
 *  ZZG --> 4l G,   WZG --> 3l nu G,   VZG --> 2l 2j G
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author A. Mecca - UNITO <alberto.mecca@cern.ch>
 */

#include <set>
#include <map>

#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"

class VVGammaAnalyzer: public EventAnalyzer, RegistrableAnalysis<VVGammaAnalyzer>{

public:

  //, const std::string& filename, const double& lumi = 1., const double& externalXSection = -1., bool doBasicPlots = false

 VVGammaAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<VVGammaAnalyzer>(*this)), 
		   configuration){
    //theHistograms.profile(genCategory);
    // Memory allocation
    leptons_      = new std::vector<const phys::Lepton*>;
    genQuarks_    = new std::vector<phys::Particle>;
    genChLeptons_ = new std::vector<phys::Particle>;
    genNeutrinos_ = new std::vector<phys::Particle>;
    genPhotons_   = new std::vector<phys::Particle>;
    
    genZlepCandidates_ = new std::vector<phys::Boson<phys::Particle>>;
    genWlepCandidates_ = new std::vector<phys::Boson<phys::Particle>>;
    genZhadCandidates_ = new std::vector<phys::Boson<phys::Particle>>;
    genWhadCandidates_ = new std::vector<phys::Boson<phys::Particle>>;
    kinPhotons_  = new std::vector<phys::Photon>;
    goodPhotons_ = new std::vector<phys::Photon>;
  }

  virtual ~VVGammaAnalyzer(){
    delete leptons_;
    delete genQuarks_;
    delete genChLeptons_;
    delete genNeutrinos_;
    delete genPhotons_;
    
    delete genZlepCandidates_;
    delete genWlepCandidates_;
    delete genZhadCandidates_;
    delete genWhadCandidates_;
    
    delete kinPhotons_;
    delete goodPhotons_;
    
    //delete photonSFhist;
  }
	
	virtual void begin();
  	
	virtual Int_t cut();

	virtual void analyze();

	virtual void end(TFile &);
	
	virtual void finish();

 private:
 	
	std::vector<const phys::Lepton*>* leptons_;
 	// Photons with pt > 20 (already in ntuples), at least loose (cut-based) ID
	std::vector<phys::Photon>* kinPhotons_;
 	std::vector<phys::Photon>* goodPhotons_;
 	
 	// Vectors of gen particles
 	std::vector<phys::Particle>* genQuarks_;
	std::vector<phys::Particle>* genChLeptons_;
	std::vector<phys::Particle>* genNeutrinos_;
	std::vector<phys::Particle>* genPhotons_;
	// Vectors of gen Bosons
	std::vector<phys::Boson<phys::Particle>>* genZlepCandidates_;
	std::vector<phys::Boson<phys::Particle>>* genWlepCandidates_;
	std::vector<phys::Boson<phys::Particle>>* genZhadCandidates_;
	std::vector<phys::Boson<phys::Particle>>* genWhadCandidates_;
	// Gen objects
	phys::DiBoson<phys::Particle, phys::Particle> genZZ_;
	phys::DiBoson<phys::Particle, phys::Particle> genWZ_;
	
	//TH2F* photonSFhist = nullptr;
 	
 	// Objects reconstruction for each event
 	void genEventSetup();
 	
 	// Basic histograms
 	void genEventHistos();
	void baseHistos_cut();
	//void baseHistos_analyze();
	void PKU_comparison();
 	
 	// Sub analyses
	void LeptonFakeRate();
	void PhotonFakeRate();
 	void effPhotons(); // uses goodPhotons_ and genPhotons_
 	
 	
 	void endNameHistos();
 	
 	// Utilities
	void initCherryPick();
	bool cherrypickEvt() const;
	std::map<phys::RegionTypes,
	  std::map<unsigned long,        // run
	    std::map<unsigned long,      // lumi block
	      std::set<unsigned long>>>> // event
	        cherryEvents;
 	/*
 	double getSF(const phys::Photon& ph){ 
 		if(photonSFhist)
 			return photonSFhist->GetBinContent(photonSFhist->FindBin(ph.eta(), ph.pt()));
 		else return 1.;
	}*/
 	
 	template<class T, class V>
	bool haveCommonDaughter(const phys::Boson<T>& a, const phys::Boson<V>& b, const float tol=0.001){
		return (
			(a.daughter(0).p4() - b.daughter(0).p4()).P() < tol ||
		  (a.daughter(0).p4() - b.daughter(1).p4()).P() < tol ||
		  (a.daughter(1).p4() - b.daughter(0).p4()).P() < tol ||
		  (a.daughter(1).p4() - b.daughter(1).p4()).P() < tol   );
	}
	
	template<class T1, class T2, class V1, class V2>
	bool haveCommonDaughter(const phys::DiBoson<T1, T2>& a, const phys::DiBoson<V1, V2>& b, const float tol=0.001){
		return (
			haveCommonDaughter(a.daughter(0), b.daughter(0), tol) ||
		  haveCommonDaughter(a.daughter(0), b.daughter(1), tol) ||
		  haveCommonDaughter(a.daughter(1), b.daughter(0), tol) ||
		  haveCommonDaughter(a.daughter(1), b.daughter(1), tol)   );
	}
 	
 	
  friend class Selector<VVGammaAnalyzer>;
  template< class PAR >
    bool ZBosonDefinition(phys::Boson<PAR> cand) const{
    bool checkCharge = cand.daughter(0).charge() + cand.daughter(1).charge() == 0;
    return checkCharge && fabs(cand.p4().M() - phys::ZMASS) < 30;
  }


  template< class PAR >
    bool WBosonDefinition(phys::Boson<PAR> cand) {
		  bool gooddaughters = false;
		  if(fabs(cand.daughter(0).eta()) < 2.5 && cand.daughter(0).pt() > 30 &&
		     cand.daughter(0).passPUID() && cand.daughter(0).passLooseJetID() &&
		     fabs(cand.daughter(1).eta()) < 2.5 && cand.daughter(1).pt() > 30 &&
		     cand.daughter(1).passPUID() && cand.daughter(1).passLooseJetID())
		    gooddaughters = true;

		  if(fabs(cand.p4().M() - phys::WMASS) < 150 && gooddaughters)
		    return true;
		  return false;

		}
  
    bool GenWBosonDefinition(phys::Boson<phys::Particle> cand) {
		  bool gooddaughters = false;
		  if(fabs(cand.daughter(0).eta()) < 2.5 && cand.daughter(0).pt() > 30 &&
		     fabs(cand.daughter(1).eta()) < 2.5 && cand.daughter(1).pt() > 30)
		    gooddaughters = true;

		  if(fabs(cand.p4().M() - phys::WMASS) < 150 && gooddaughters)
		    return true;
		  return false;
  	}
		
		static const std::vector<double> pt_bins;
		static const std::vector<double> pt_bins_LFR;
		static const std::vector<double> eta_bins;
		static const std::vector<double> aeta_bins;
  
  protected:
		unsigned long evtN_;  // Used in the print in cut()
		std::map<phys::RegionTypes, unsigned long> evtNInReg_, analyzedNInReg_;  // Used to count processed events.
		clock_t startTime_;  // Used to calculate elapsed time
		std::map<phys::RegionTypes, float> evtWInReg_, analyzedWInReg_;  // Weighted events passing cut and total
};
#endif

