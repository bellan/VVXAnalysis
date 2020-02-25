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
    	//candType = VCandType::None;
 	 	}

		virtual ~VZZAnalyzer(){}
  	
		virtual void begin();
  	
		virtual Int_t cut();

		virtual void analyze();
  
		virtual void end(TFile &);
		
		void simpleGraphs();
		void jetRecoGraphs();
		template<class P = phys::Jet>
		bool findBestVFromSing(/*const*/ std::vector<P>*, const P*& VCandidate, VCandType& candType);
		//Passing pointer by reference so that it gets updated
		template <class J = phys::Jet>
		bool findBestVFromPair(const std::vector<J>*, phys::Boson<J>*& VCandidate, VCandType& candType);
		//Searches among the Jets in the vector (which is likely either "jets" or "jetsAK8") and finds the pair candidate with mass closest to W or Z (modifying "candType"), and stores it in "VCandidate". true: candidate fits the (W/Z) BosonDefinition
		template <class P = phys::Particle>
		bool findBestVPoint(std::vector<const P*>* js, const P*& thisCandidate, VCandType& thisCandType); //Uses a vector<P*> instead of a vector<P>
		
		template <class P, class R = phys::Boson<phys::Particle>> // P = Jet or Particle
		P*              closestSing(std::vector<P>* cands, const R& reference);
		template <class P, class R = phys::Boson<phys::Particle>> // P = Jet or Particle
		phys::Boson<P>* closestPair(std::vector<P>* cands, const R& reference);
		
		void bestCandidateAnalisys();
		void ptCutMVA(); //select the best cut in pt for the jetsAK8 --> corrPrunedMass
		void jetAnalisys(); //how often ak8 are better than ak4? Does it depend on pt? |p|? Mass?
		
		
	private:
		//Counters, etc.
		clock_t startTime_; //Used to calculate elapsed time
		unsigned long evtN_; //Used to count processed events
		unsigned int singWFromJets_, pairWFromJets_, singWFromJetsAK8_, pairWFromJetsAK8_;
		unsigned int singZFromJets_, pairZFromJets_, singZFromJetsAK8_, pairZFromJetsAK8_;
		unsigned int recVBtot_; //Evts where the reconstructed VB is acceptable (VBosonDefinition)
		unsigned int goodRec_, withGenVB_; //Evts wit a recVB near a genVB / Evts with a genVB
		//unsigned int finalFromJetsSing_, finalFromJetsPair_, finalFromJetsAK8sing_, finalFromJetsAK8Pair_;
		//phys::Boson<phys::Jet> VCandidate_;
		//VCandType candType_;
		
		friend class Selector<VZZAnalyzer>;
		
		
		//Signal definition
		template <class PAR>
		bool ZBosonDefinition(phys::Boson<PAR>& cand) const{  //candidate
			bool checkMass = fabs(cand.p4().M() - phys::ZMASS) < 40; //temp
			return checkMass;
		}
		template <class PAR>	 //Apparently, PAR must inherit from Jet
		bool WBosonDefinition(phys::Boson<PAR>& cand) const{	//candidate
			bool checkMass = fabs(cand.p4().M() - phys::WMASS) < 40;
			return checkMass;
		}
		
		template <class PAR>  //usually PAR is Jet
		bool ZBosonDefinition(PAR& cand) const{
			bool checkMass = fabs(cand.p4().M() - phys::ZMASS) < 40;
			return checkMass;
		}
		template <class PAR>	
		bool WBosonDefinition(PAR& cand) const{	
			bool checkMass = fabs(cand.p4().M() - phys::WMASS) < 40;
			return checkMass;
		}
		/*
		bool isGoodW(const phys::Boson<phys::Particle>& genV) const{
			bool goodMass = fabs(genV.p4().M() - phys::WMASS) < 40;
			bool goodId = abs(genV.id()) == 24;
			bool goodDaugthers = (genV.daughter(0)).pt()>30 &&  (genV.daughter(1)).pt()>30;
			return (goodId && goodMass && goodDaugthers);
		}
		bool isGoodZ(const phys::Boson<phys::Particle>& genV){
			bool goodMass = fabs(genV.p4().M() - phys::ZMASS) < 40;
			bool goodId = genV.id() == 23;
			bool goodDaugthers = (genV.daughter(0)).pt()>30 &&  (genV.daughter(1)).pt()>30;
			return (goodId && goodMass && goodDaugthers);
		}*/
};

#endif
