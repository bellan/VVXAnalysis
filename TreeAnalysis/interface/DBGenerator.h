#ifndef DBGenerator_h
#define DBGenerator_h

/** \Class DBGenerator
 *  Used to generate a database (.csv for now) for scikit Machine Learning, 
 *	therefore it deviates somewhat from the intended purpose of an EventAnalyzer.
 *
 *  $Date: 2020/01/23 17:00:00 $
 *  $Revision: 1.0 $
 *  \author A. Mecca alberto.mecca@edu.unito.it
 */

#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VZZAnalyzer.h"

#include <iostream>
#include <fstream>

//class DBGenerator: public EventAnalyzer, RegistrableAnalysis<DBGenerator>{

class DBGenerator: public VZZAnalyzer, RegistrableAnalysis<DBGenerator>{
	public:
		//enum VCandType {None, W, Z}; //VCandidate is closest to Wmass or Zmass?
		const std::vector<std::string> signalNames_ = {"WZZ", "ZZZ"};
		
		DBGenerator(const AnalysisConfiguration& configuration) 
				//: EventAnalyzer(*(new Selector<DBGenerator>(*this)), configuration){
				: VZZAnalyzer(configuration) {
    	//theHistograms.profile(genCategory);
 	 	}

		virtual ~DBGenerator(){}
  	
		virtual void begin();
  	
		virtual Int_t cut();

		virtual void analyze();
  
		virtual void end(TFile &);
		
		void printZeroes(size_t nzeros);
		void printVars(size_t n, ...);
		
		// ----- ----- Event-specific variables calculation ----- ----- 
		//void fillGenHadVBs(); // see VZZAnalyzer.h
		
		template <class J = phys::Jet>
		//phys::Boson<J>* findBestVFromPair(const std::vector<J>*, VCandType& candType);
		
		int isSignal();
	private:
		//std::vector<phys::Boson<phys::Particle>>* genHadVBs_ = nullptr;  // genVBParticles with hadronic daugthers
		
		static const char SEP_CHAR = ',';	//separatory char used in the .csv
		clock_t startTime; //Used to calculate elapsed time
		unsigned long evtN; //Used to count processed events
		bool isSigFile_ = false;  // Set during begin()
		
		std::ofstream outputFile; //a .csv file the data in the tree is written to
		
		friend class Selector<DBGenerator>;
		
		
		//Miscellanous stuff required by EventAnalyzer
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
};

#endif
