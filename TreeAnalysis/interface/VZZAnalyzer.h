#ifndef VZZAnalyzer_h
#define VZZAnalyzer_h

/** \Class VZZAnalyzer
 *  Looks for ZZZ or WZZ, with ZZ -> 4l and W,Z -> JJ.
 *
 *  $Date: 2020/01/27 12:00:00 $
 *  $Revision: 1.0 $
 *  \author A. Mecca alberto.mecca@edu.unito.it
 */

#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"


#include <iostream>
#include <fstream>

class VZZAnalyzer: public EventAnalyzer, RegistrableAnalysis<VZZAnalyzer>{
	public:
		enum VCandType {None, W, Z}; //VCandidate is closest to Wmass or Zmass?
		VZZAnalyzer(const AnalysisConfiguration& configuration) 
				: EventAnalyzer(*(new Selector<VZZAnalyzer>(*this)), configuration){
    	//theHistograms.profile(genCategory);
    	candType = VCandType::None;
 	 	}

		virtual ~VZZAnalyzer(){}
  	
		virtual void begin();
  	
		virtual Int_t cut();

		virtual void analyze();
  
		virtual void end(TFile &);
		
		void simpleGraphs();
		void jetRecoGraphs();
		bool findBestVCandidate(); //true: candidate fits the (W/Z)BosonDefinition
		
	private:
		//Counters, etc.
		clock_t startTime_; //Used to calculate elapsed time
		unsigned long evtN_; //Used to count processed events
		phys::Boson<phys::Jet> VCandidate;
		VCandType candType;
		
		friend class Selector<VZZAnalyzer>;
		
		
		//Miscellanous stuff required by EventAnalyzer
		template< class PAR >
		bool ZBosonDefinition(phys::Boson<PAR> cand) const{  //candidate
			bool checkMass = fabs(cand.p4().M() - phys::ZMASS) < 50;
			return checkMass;
		}

		template< class PAR >	//Apparently, PAR must inherit from Jet
		bool WBosonDefinition(phys::Boson<PAR> cand ) const{	//candidate
			bool checkMass = fabs(cand.p4().M() - phys::WMASS) < 50;
			return checkMass;
		}
};

#endif
