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
    leptons_      = new std::vector<phys::Lepton>;
    genQuarks_    = new std::vector<phys::Particle>;
    genChLeptons_ = new std::vector<phys::Particle>;
    genNeutrinos_ = new std::vector<phys::Particle>;
    genPhotons_   = new std::vector<phys::Particle>;
    
    genZlepCandidates_ = new std::vector<phys::Boson<phys::Particle>>;
    genWlepCandidates_ = new std::vector<phys::Boson<phys::Particle>>;
    genZhadCandidates_ = new std::vector<phys::Boson<phys::Particle>>;
    genWhadCandidates_ = new std::vector<phys::Boson<phys::Particle>>;
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
  }
	
	virtual void begin();
  	
	void initEvent();
	virtual Int_t cut();

	virtual void analyze();

	virtual void end(TFile &);
	
	virtual void finish();


 private:
	
	std::vector<phys::Lepton>* leptons_;
	
	// Systematics: photons {EScale, ESigma} x {Up, Down} + {central}
	std::map<const char*, std::unique_ptr<std::vector<phys::Photon>>> kinPhotons_;    // Only kinematic selection
	std::map<const char*, std::unique_ptr<std::vector<phys::Photon>>> loosePhotons_;  // Loose ID: currently 3/5 cuts of ID
	std::map<const char*, std::unique_ptr<std::vector<phys::Photon>>> goodPhotons_;   // Tight ID: currently Loose WP of POG cut-based ID
 	
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
	// V --> j (j)
	phys::Boson<phys::Jet> candVTojj_, fakeVTojj_;
	phys::Jet              candVToJ_ , fakeVToJ_ ;
	
	std::unique_ptr<TH2F> hPhotonFR_;
	std::unique_ptr<TH2F> hPhotonEff_;
        std::string channelReco_;
 	
 	// Objects reconstruction for each event
        void makeChannelReco();  // sets channelReco_
	void genEventSetup();
	/* const phys::Jet* candAK8(const std::vector<phys::Jet>*) const; */
	void hadronicObjectsReconstruction();
	static std::vector<phys::Boson<phys::Jet>> candidatesVTojj(const std::vector<phys::Jet>&);
 	
 	// Basic histograms
 	void genEventHistos();
	void baseHistos_cut();
	void photonHistos();
	void jetHistos();
	//void baseHistos_analyze();
	void PKU_comparison();
 	
 	// Sub analyses
	void plotsVVGstatus(const char* name, const char* title, const TLorentzVector& p4_VV, const char* mType="mass");
	void leptonFakeRate();
	void photonFakeRate();
 	void photonEfficiency(const std::vector<phys::Photon>&, const char*);
 	void photonIsolation(const std::vector<phys::Photon>&, const char*);
	void systematicsStudy();
	void SYSplots(const char* syst, const double& weight, const phys::Photon*);
        void debug3Lregion();
	void photonGenStudy();

 	void endNameHistos();
 	
 	// Utilities
	template <class T, class UnaryPredicate>
	  static std::vector<phys::Boson<T>> makeBosons(const std::vector<T>&, UnaryPredicate);
	void initCherryPick();
	bool cherrypickEvt() const;
	std::map<phys::RegionTypes,
	  std::map<unsigned long,        // run
	    std::map<unsigned long,      // lumi block
	      std::set<unsigned long>>>> // event
	        cherryEvents;
 	
        double getPhotonFR   (const phys::Photon& ph) const;
        double getPhotonFRUnc(const phys::Photon& ph) const;

	double getPhotonEff   (const phys::Photon& ph) const;
	double getPhotonEffUnc(const phys::Photon& ph) const;
	
	static bool is4Lregion(const phys::RegionTypes reg){
	  return reg == phys::SR4P || reg == phys::CR3P1F || reg == phys::CR2P2F ||
	    reg == phys::SR4P_1L || reg == phys::SR4P_1P || reg == phys::CR4P_1F;
	}
	static bool is3Lregion(const phys::RegionTypes reg){
	  return reg == phys::SR3P || reg == phys::CR110 || reg == phys::CR101 || reg == phys::CR011 ||
	    reg == phys::CR100 || reg == phys::CR010 || reg == phys::CR001 || reg == phys::CR000 ||
	    reg == phys::SR3P_1L || reg == phys::SR3P_1P || reg == phys::CR3P_1F;
	}
	static bool is2Lregion(const phys::RegionTypes reg){
	  return reg == phys::SR2P || reg == phys::SR2P_1L || reg == phys::SR2P_1P || reg == phys::CR2P_1F;
	}
	
	static bool passVeryLoose(const phys::Photon& ph);
	
	static char phABCD(const phys::Photon&, const phys::Photon::IDwp);
	static char phABCD_study(const phys::Photon&, const double& barrel_thr, const double& endcap_thr);
 	
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
		static const std::vector<double> ph_pt_bins;
		static const std::vector<double> ph_aeta_bins;
		static const std::vector<double> mVV_bins;
		static const std::vector<double> mVVG_bins;
  
  protected:
		unsigned long evtN_;  // Used in the print in cut()
		std::map<phys::RegionTypes, unsigned long> evtNInReg_, analyzedNInReg_;  // Used to count processed events.
		clock_t startTime_;  // Used to calculate elapsed time
		std::map<phys::RegionTypes, float> evtWInReg_, analyzedWInReg_;  // Weighted events passing cut and total
};
#endif

