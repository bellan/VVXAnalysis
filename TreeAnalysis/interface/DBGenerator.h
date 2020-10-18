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
		
		// ----- ----- Write to databases ----- -----
		void mainEvtRec(int sigRecType, const phys::Particle* recV);
		
		void writeTagger(std::ofstream& outFile4, std::ofstream& outFile8);  // We want to choose the best AK8 (pair of AK4), using info from the original quarks
		void writeInfoAK4(const phys::Boson<phys::Jet>* bestAK4, std::ofstream& outFile);
		void writeInfoAK8(const phys::Jet* bestAK8, std::ofstream& outFile);
		
		// ----- ----- Helpers ----- -----
		static void printZeroes(std::ofstream& outFile, size_t nzeros);
		static void printVars(std::ofstream& outFile, size_t n, ...);
		static void printVars(std::ofstream& outFile, size_t n, const double* buf);
		static void printVars(std::ofstream& outFile, const std::vector<double>& vect);
		
		// This will be replaced by VZZAnalyzer::getAK4Features() once a new predictor is trained
		//std::vector<double> getAK4FeaturesTEMP(const phys::Boson<phys::Jet>& jj);
		
	private:
		static const char SEP_CHAR = ',';	//separatory char used in the .csv
		//clock_t startTime_; //Used to calculate elapsed time
		//unsigned long evtN_; //Used to count processed events
		bool isSigFile_ = false;  // Set during begin()
		
		std::ofstream outputFile_; //a .csv file the data in the tree is written to
		std::ofstream outputFileAK4_; //a .csv file; this one is for the Jet Tagger
		std::ofstream outputFileAK8_; //a .csv file; this one is for the Jet Tagger
		
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
