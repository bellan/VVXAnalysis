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
		
		// ----- ----- Helper functions ----- ----- 
		template<class P = phys::Jet>
		bool findBestVFromSing(/*const*/ std::vector<P>*, const P*& VCandidate, VCandType& candType);
		//Passing pointer by reference so that it gets updated
		template <class J = phys::Jet>
		bool findBestVFromPair(const std::vector<J>*, phys::Boson<J>*& VCandidate, VCandType& candType);
		//Searches among the Jets in the vector (which is likely either "jets" or "jetsAK8") and finds the pair candidate with mass closest to W or Z (modifying "candType"), and stores it in "VCandidate". true: candidate fits the (W/Z) BosonDefinition
		template <class P = phys::Particle>
		bool findBestVPoint(std::vector<const P*>* js, const P*& thisCandidate, VCandType& thisCandType); //Uses a vector<P*> instead of a vector<P>
		
		template <class P, class R = phys::Boson<phys::Particle>> // P = Jet or Particle
		P*              closestSing(std::vector<P>* cands, const R& reference); //max dR = 0.4
		template <class P, class R = phys::Boson<phys::Particle>> // P = Jet or Particle
		phys::Boson<P>* closestPair(std::vector<P>* cands, const R& reference); //max dR = 0.4
		template <class P, class R = phys::DiBoson<phys::Lepton, phys::Lepton>>
		P* furthestSing(std::vector<P>* cands, const R& reference, const float& minDR = 2., const std::pair<float,float>& massLimits = std::make_pair(1.,1000.));
		template <class P, class R = phys::DiBoson<phys::Lepton, phys::Lepton>>
		phys::Boson<P>* furthestPair(std::vector<P>* cands, const R& reference, const float& minDR = 2., const std::pair<float,float>& massLimits = std::make_pair(1.,1000.));
		
		template <class P = phys::Particle>
		inline double getRefinedMass(const P& p) const{ return p.mass(); } 
		// If I used "const P&" it would have precedence on the template specialization when the object is a pointer to non-const, since "const P&" matches const P* & (const reference to non-const pointer to non-const)
		//template <class P = phys::Particle*> 
		//inline double getRefinedMass(P* p) const{ return p->mass();}
		
		//I give up, c++ templates are too complicated when mixed with cv-qualifiers. It's easier to write the function overloads
		inline double getRefinedMass(const phys::Particle* p) const{ return p->mass(); }
		template <class P = phys::Particle>
		inline double getRefinedMass(const phys::Boson<P>* p) const{ return p->mass(); }
		
		//Function overload of the template for Jets 
		inline double getRefinedMass(const phys::Jet& j) const{ return j.corrPrunedMass(); }
		inline double getRefinedMass(const phys::Jet* j) const{ return j->corrPrunedMass();}
		
		// ----- ----- Event-specific varibles calculation ----- ----- 
		void fillGenHadVBs(); //old: Fills the vector only if it is empty
		void fillRecHadVBs();  //old: Fills the vector only if it is empty
		void calcS();
		
		// ----- ----- Large sub-analisys ----- ----- 
		void simpleGraphs();
		void jetRecoGraphs();
		
		void bestCandidateAnalysis();
		void endBestCandAnalysis(TFile& fout); //Divide histograms to obtain efficiency
		
		void ptCutMVA(); //select the best cut in pt for the jetsAK8 --> corrPrunedMass
		
		void closestJetAnalisys(); //how often ak8 are better than ak4? Does it depend on pt? |p|? Mass?
		void endClosestJetAn(); //Histogram normalization
		
		void furthestJetMVA(); //Looking for a variable that lets us choose whether to use an AK4 or an AK8
		void minPtJetMVA(); //Still looking. Maybe minimizing total pt can help us choose AK4/8
		
		void bestZMassJetMVA();
		void specialPeakAnalisys(const phys::Particle& theGenAK8); //Is there a pair of AK4 that reconstructs theese events with an AK8 but low-pt ZZ?
		
		void endResolutionAnalisys(TFile& fout); //Calculates AK4-AK8 resolution per bin of ZZ pt
		
	private:
		std::vector<phys::Boson<phys::Particle>>* genHadVBs_ = nullptr; //genVBParticles with hadronic daugthers
		std::vector<phys::Boson<phys::Jet>>* AK4pairs_ = nullptr; //pairs of all the reconstructed AK4 jets
		
		float sAK4g_; //invariant mass of ZZ and all the AK4s gen
		float sAK8g_; //                                 AK8s
		float sAK4r_; //                                 AK4s rec
		float sAK8r_; //                                 AK8s
		// We will use #hat{s} for proper energy of ZZ+best AK8/pair pf AK4
		
		// ----- ----- Counters, ecc. ----- ----- 
		clock_t startTime_; //Used to calculate elapsed time
		unsigned long evtN_, analyzedN_; //Used to count processed events
		unsigned int singWFromJets_, pairWFromJets_, singWFromJetsAK8_, pairWFromJetsAK8_;
		unsigned int singZFromJets_, pairZFromJets_, singZFromJetsAK8_, pairZFromJetsAK8_;
		unsigned int recVBtot_; //Evts where the reconstructed VB is acceptable (VBosonDefinition)
		unsigned int goodRec_, withGenVB_; //Evts wit a recVB near a genVB / Evts with a genVB
		unsigned int win4_,win8_;//How often an AK4/AK8 reconstructs better an hadronic decaying VB
		
		friend class Selector<VZZAnalyzer>;
		
		
		// ----- ----- Signal definition ----- ----- 
		int isSignal() const;  // 0-->no  1-->AK4  2-->AK8
		
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
			bool checkMass = fabs(getRefinedMass(cand) - phys::ZMASS) < 40;
			return checkMass;
		}
		template <class PAR>	
		bool WBosonDefinition(PAR& cand) const{	
			bool checkMass = fabs(getRefinedMass(cand) - phys::WMASS) < 40;
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
		inline bool tauCut(const phys::Jet& jet8, const double& thr = 0.6){
			return (jet8.tau2()/jet8.tau1() < thr);
		}
		inline bool tauCut(const phys::Jet* jet8, const double& thr = 0.6){
			return tauCut(*jet8, thr);
		}
		inline bool tauCut(std::vector<phys::Jet>::iterator jet8, const double& thr = 0.6){
			return tauCut(*jet8, thr);
		}
};

#endif
