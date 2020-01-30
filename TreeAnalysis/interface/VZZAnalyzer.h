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
		VZZAnalyzer(const AnalysisConfiguration& configuration) 
				: EventAnalyzer(*(new Selector<VZZAnalyzer>(*this)), configuration){
    	//theHistograms.profile(genCategory);
 	 	}

		virtual ~VZZAnalyzer(){}
  	
		virtual void begin();
  	
		virtual Int_t cut();

		virtual void analyze();
  
		virtual void end(TFile &);
		
	private:
		clock_t startTime; //Used to calculate elapsed time
		unsigned long evtN; //Used to count processed events
		
		friend class Selector<VZZAnalyzer>;
		
		
		
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
