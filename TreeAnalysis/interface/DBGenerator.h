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

#include <iostream>
#include <fstream>

class DBGenerator: public EventAnalyzer, RegistrableAnalysis<DBGenerator>{
	public:
		enum VCandType {None, W, Z}; //VCandidate is closest to Wmass or Zmass?
		
		DBGenerator(const AnalysisConfiguration& configuration) 
				: EventAnalyzer(*(new Selector<DBGenerator>(*this)), configuration){
    	//theHistograms.profile(genCategory);
 	 	}

		virtual ~DBGenerator(){}
  	
		virtual void begin();
  	
		virtual Int_t cut();

		virtual void analyze();
  
		virtual void end(TFile &);
		
		void printZeros(size_t nzeros);
		void printVars(int n, ...);
		bool findBestVCandidate(const std::vector<phys::Jet>*, phys::Boson<phys::Jet>& VCandidate, VCandType& candType);
	private:
		static const char SEP_CHAR = ',';	//separatory char used in the .csv
		clock_t startTime; //Used to calculate elapsed time
		unsigned long evtN; //Used to count processed events
		
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
