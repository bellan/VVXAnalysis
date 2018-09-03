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
  
		virtual void end();

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
		
		bool electronMatch(const phys::Particle &, const phys::Particle &);
		void findElectronMatch(const phys::Particle &, std::vector<phys::Particle> *);
		
		void inTheLastEvent();		//I had problems with the execution of end()
		
		void fillBasicPlots();
		void fillParticlePlots(const std::string &, const phys::Particle & );
		
		void doSomeFits();
		void getFitInfo(TF1*);
		
		void initStatistics();
		void tempStatisticParticles(const phys::Particle&);
		void tempStatisticEvents();
		
		//Various counters
		int counter;
		int eCounter;
		int mCounter;
	
		int promptCounter;
		int peCounter;
		int pmCounter;
		
		int electronEvents;
		int muonEvents;
		int passingSelection;
		
		long long matchedElectrons;
		long long totalElectrons;
};
#endif


/*  LocalWords:  WWosAnalyzer
 */
