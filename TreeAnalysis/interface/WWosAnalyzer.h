#ifndef WWosAnalyzer_h
#define WWosAnalyzer_h

/** \Class WWosAnalyzer
 *  Concrete class for WW with opposite sign analysis
 *
 *  $Date: 2018/08/04 13:37:31 $
 *  $Revision: 0.1 $
 *  \author A. Mecca alberto.mecca@edu.unito.it
 */

#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"

//#define TTJets_SMALL
#define WWEW
//#define WWQCD
//#define WWEWQCD

#ifdef TTJets_SMALL	//40 Mb sample
	#define NUMBER_OF_EVENTS 9055
#elif defined WWEW
	#define NUMBER_OF_EVENTS 499500
#elif defined WWQCD
	#define NUMBER_OF_EVENTS 499600
#elif defined WWEWQCD
	#define NUMBER_OF_EVENTS 467264
#elif defined TTTo2L2Nu
	#define NUMBER_OF_EVENTS 18539890
#else
	#include <climits>
	#define NUMBER_OF_EVENTS ULONG_MAX
#endif

//#define DO_STATISTICS_ON_PARTICLES
//#define DO_STATISTICS_ON_EVENTS

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

		template< class PAR >
		bool WBosonDefinition(phys::Boson<PAR> cand ) {	//candidate

			bool gooddaughters = false;
    	
			if(fabs(cand.daughter(0).eta()) < 2.5 && cand.daughter(0).pt() > 30 &&
					cand.daughter(0).passPUID() && cand.daughter(0).passLooseJetID() &&
					fabs(cand.daughter(1).eta()) < 2.5 && cand.daughter(1).pt() > 30 &&
					cand.daughter(1).passPUID() && cand.daughter(1).passLooseJetID()
				)
				gooddaughters = true;

			if(fabs(cand.p4().M() - phys::WMASS) < 150 && gooddaughters)
				return true;
			return false;
		}
		
		//Function for efficiency analysis
		template <class T, class P, typename C>
		void analyzeEfficiency(std::vector<T>* gen, std::vector<P>* rec, std::string name, C& counter);
		
		template <class P, class T>
		bool checkMatch(const /*phys::Particle&*/P& reconstructed, const /*phys::Particle&*/ T& generated,  const float& tolerance);	//Checks deltaR
		phys::Particle* findMatchingParticle(const phys::Particle& rec, std::vector<phys::Lepton>* candidates);
		
		void normalizeHistograms(std::string name);
		
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
		
		//Generated particles (prompt)
		std::vector<phys::Particle>* genElectrons;	//genElectrons in this event
		std::vector<phys::Particle>* genMuons;	//muons in this event
		std::vector<phys::Particle>* genLeptons;	//every lepton in this event
		
		#ifdef DO_STATISTICS_ON_PARTICLES
		//Various counters
		int particleCounter;
		int eCounter;
		int mCounter;
	
		int promptCounter;
		int peCounter;
		int pmCounter;
		#endif
		
		int electronEvents;
		int muonEvents;
		int passingSelection;
		
		long matchedElectrons;
		long matchedMuons;
		int wrongChargeE;
		int wrongChargeM;
		long totalElectrons;
		long totalMuons;
		
		clock_t startTime;
};

void getFitInfo(TF1*);
#endif


/*  LocalWords:  WWosAnalyzer
 */
