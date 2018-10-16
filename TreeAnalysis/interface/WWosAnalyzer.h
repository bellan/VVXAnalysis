#ifndef WWosAnalyzer_h
#define WWosAnalyzer_h

/** \Class WWosAnalyzer
 *  Concrete class for WW with opposite sign analysis
 *
 *  $Date: 2018/08/04 13:37:31 $
 *  $Revision: 0.3 $
 *  \author A. Mecca alberto.mecca@edu.unito.it
 */

#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"

// TTJets_SMALL      9055 (40 Mb sample) 
// WWEW            499500
// WWQCD           499600
// WWEWQCD         467264
// TTTo2L2Nu     18539890
// DYJetsToLL_M50 3014382

//##### New samples #####
// WWEW            499500 	~1'15"
// WWQCD           499600 	~1'22"
// WWEWQCD         467264		~1'17"
// WZ              163611 	~0'31"
// TTJets          166862 	~0'34"
// TTTo2L2Nu      2321521 	~7'44"
// DYJetsToLL_M50   22199 	~0'04"

//#define DO_GEN_PARTICLES_ANALYSIS
#ifdef DO_GEN_PARTICLES_ANALYSIS
	//#define DO_STATISTICS_ON_PARTICLES
	//#define DO_STATISTICS_ON_EVENTS
	//#define DO_EFFICIENCY_ANALYSIS
#endif

#define LEPTON_CUT	//Cuts events in which there are not exactly 2 leptons
//#define Jet_CUT			//Cuts events in which there are not less than 2 jets

class WWosAnalyzer: public EventAnalyzer, RegistrableAnalysis<WWosAnalyzer>{
	public:
		WWosAnalyzer(const AnalysisConfiguration& configuration) 
				: EventAnalyzer(*(new Selector<WWosAnalyzer>(*this)), configuration){
    	//theHistograms.profile(genCategory);
 	 	}

		virtual ~WWosAnalyzer(){}
  	
		virtual void begin();
  	
		virtual Int_t cut();

		virtual void analyze();
  
		virtual void end(TFile &);

	private:
		friend class Selector<WWosAnalyzer>;
		template< class PAR >
		bool ZBosonDefinition(phys::Boson<PAR> cand) const{  //candidate
			bool checkCharge = cand.daughter(0).charge() + cand.daughter(1).charge() == 0;
			bool checkMass = fabs(cand.p4().M() - phys::ZMASS) < 30;
			return checkCharge && checkMass;
		}

		template< class PAR >	//Apparently, PAR must inherit from Jet
		bool WBosonDefinition(phys::Boson<PAR> cand ) {	//candidate
			bool gooddaughters = false;
			if(fabs(cand.daughter(0).eta()) < 2.5 && cand.daughter(0).pt() > 30 &&
					cand.daughter(0).passPUID() && cand.daughter(0).passLooseJetID() &&
					fabs(cand.daughter(1).eta()) < 2.5 && cand.daughter(1).pt() > 30 &&
					cand.daughter(1).passPUID() && cand.daughter(1).passLooseJetID()
				) //passLooseJetID and passPUID are inherited from VVXAnalysis/DataFormats/interface/Jet.h
				gooddaughters = true;

			if(fabs(cand.p4().M() - phys::WMASS) < 150 && gooddaughters)
				return true;
			return false;
		}
		
		bool leptonCut();
		bool jetCsvtagCut(const float &max = 1.){ return jetCsvtagCut(max,max); };
		bool jetCsvtagCut(const float &leadMax, const float &trailMax);
		bool newCuts();
				
		void leptonPlots(const phys::Particle* lead, const phys::Particle* trail, const std::string& type = std::string(""), bool useWeight = true);
		void leptonCutAnalysis(const phys::Particle* lead, const phys::Particle* trail, const std::string& type = std::string(""), bool useWeight = true);
		void lepton2DGraph(const phys::Particle* lead, const phys::Particle* trail, const std::string& type = std::string(""), bool useWeight = true);
		void nestedPtCutHistogram(const phys::Particle* lead, const phys::Particle* trail, const std::string& type, float ptCut, float weight);
		
		void jetPlots(const std::string& type = "");
		void miscPlots(const phys::Particle* lead, const phys::Particle* trail, const std::string& type = std::string(""), bool useWeight = true);
		
		void genParticlesAnalysis();	//All the work realate to efficiency/resolution analysis
		void genParticlesCategorization();	//Divides genParticle by the id
		
		void endGenParticleAnalysis(); //stuff from end();
		void nameCutGraph();				//Sets the name of the bins in the graph "Cuts"
		
		//Function for efficiency analysis
		template <class T, class P, typename C>
		void analyzeEfficiency(std::vector<T>* gen, std::vector<P>* rec, std::string name, C& counter, Float_t maxDeltaR = 0.2);
		
		template <class P, class T>
		bool checkMatch(const /*phys::Particle&*/P& reconstructed, const /*phys::Particle&*/ T& generated,  const float& tolerance);	//Checks deltaR
		
		template <class P>
		phys::Particle* findMatchingParticle(const phys::Particle& rec, std::vector<P>* candidates);
		
		#ifdef DO_EFFICIENCY_ANALYSIS
		template <class P, class T>
		void resolutionAnalysis(const T& rec, const P& gen, std::string name);
		void fitResolutionPt(std::string name);
		void fitResolutionE(std::string name);
		
		void normalizeHistograms(std::string name);
		void normalizeEta(std::string name);
		void normalizePhi(std::string name);
		void normalizePt(std::string name);
		TGraphAsymmErrors* myTGraphAsymmErrors(TH1* num, TH1* denom);
		#endif
		
		void fillBasicPlots();
		void fillParticlePlots(const std::string &, const phys::Particle & );
		
		void doSomeFits();
		
		#ifdef DO_STATISTICS_ON_PARTICLES
		void initStatistics();
		void tempStatisticParticles(const phys::Particle&);
		#endif
		#ifdef DO_STATISTICS_ON_EVENTS
		void tempStatisticEvents();
		#endif
		
		phys::Lepton leadLepton;	//electron or muon with the largest pt among the leptons in this event
		phys::Lepton trailLepton;	//        ...               second largest pt          ...
		enum WWosEventTypes {ee, em, mm};
		WWosEventTypes thisEventType;		//{"ee", "em", "mm"}
		
		//Generated particles (prompt)
		std::vector<phys::Particle>* genElectrons;		//genElectrons in this event
		std::vector<phys::Particle>* genMuons;				//genMuons in this event
		std::vector<phys::Particle>* genLeptons;			//every genLepton in this event
		std::vector<phys::Particle>* genNeutrinos;		
		
		std::vector<phys::Particle>* genCleanedJets;	//every prompt genJet in this event	
		std::vector<phys::Particle>* fakeJets; 				//isolated particles that appear in genJets
		
		#ifdef DO_STATISTICS_ON_PARTICLES
		//Various counters
		int particleCounter;
		int eCounter;
		int mCounter;
	
		int promptCounter;
		int peCounter;
		int pmCounter;
		#endif
		
		int lostEvents;
		
		int eeEvents; //electronEvents;
		int mmEvents; //muonEvents;
		int emEvents; //mixedEvents;
		int multiSignEvents;
		int passingSelection;
		int passingCut;
		
		#ifdef DO_GEN_PARTICLES_ANALYSIS
		long matchedElectrons;
		long matchedMuons;
		int matchedJets;
		int wrongChargeE;
		int wrongChargeM;
		long totalElectrons;
		long totalMuons;
		
		long totjets;
		long totgenJets;
		long totgenCleanedJets;
		long totgenParticles;
		#endif
		
		clock_t startTime;
		unsigned long evtN;
};

void getFitInfo(TF1*);
//void confidence_binomial_clopper_pearson(int n, int k, double &xlow, double &xhigh, double level=.68540158589942957);
#endif


/*  LocalWords:  WWosAnalyzer
 */
