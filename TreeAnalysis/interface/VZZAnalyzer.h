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
#include "VVXAnalysis/Commons/interface/Comparators.h"

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <iostream>
#include <fstream>
#include <algorithm>  // std::min


class VZZAnalyzer: public EventAnalyzer, RegistrableAnalysis<VZZAnalyzer>{
	public:
		enum VCandType {None, W, Z}; //VCandidate is closest to Wmass or Zmass?
		const char* massAlgsNames_[6] = {"mass", "secvtxMass", "corrPrunedMass", "prunedMass", "softDropMass", "puppiMass"};
		
		VZZAnalyzer(const AnalysisConfiguration& configuration) 
				: EventAnalyzer(*(new Selector<VZZAnalyzer>(*this)), configuration){
    	//theHistograms.profile(genCategory);
    	//std::cout<<"Creating VZZAnalyzer"<<std::endl;
 	 	}

		virtual ~VZZAnalyzer(){
			//std::cout<<"Destroying VZZAnalyzer"<<std::endl;
			if(genHadVBs_) delete genHadVBs_; //Deallocates memory
			if(AK4pairs_)  delete AK4pairs_;
			if(AllGenVBjj_)delete AllGenVBjj_;
			if(AK4GenRec_) delete AK4GenRec_;
			if(genZZ_)     delete genZZ_;
			//if(qq_)        delete qq_;
			if(sigVB_)     delete sigVB_;
			
			Py_XDECREF(AK4_classifier_);  // Free memory in event of crash (end() is not called)
			Py_XDECREF(AK8_classifier_);
			Py_XDECREF(helper_module_);  // XDECREF: checks if the reference count is >=0
			//Py_FinalizeEx();  // Close Python interpreter
		}
  	
		virtual void begin();
  	
		virtual Int_t cut();

		virtual void analyze();
  
		virtual void end(TFile &);
		
		// ----- ----- Helper functions ----- ----- 
		template<class P = phys::Jet>
		const P* findBestVFromSing(std::vector<P>*, VCandType&);
		//The candidate is NOT copied, and the pointer returned points to the original
		template<class P = phys::Jet> const P* findBestVFromSing(std::vector<P>* v){
			VCandType temp = VCandType::None;
			return findBestVFromSing(v, temp);
		}
		template <class J = phys::Jet>
		phys::Boson<J>* findBestVFromPair(const std::vector<J>*, VCandType&);
		//Searches among the Jets in the vector and finds the pair candidate with mass closest to W or Z (modifying "candType"). Returns the candidate only if it fits the (W/Z)BosonDefinition
		template <class J = phys::Jet> phys::Boson<J>* findBestVFromPair(const std::vector<J>* v){
			VCandType temp = VCandType::None;
			return findBestVFromPair(v, temp);
		}
		template <class P = phys::Particle>
		const P* findBestVPoint(std::vector<const P*>& js, VCandType& thisCandType); //Uses a vector<P*> instead of a vector<P>
		
		//Implemented in VZZAnalyzer_impl.cc
		template <class P, class R = phys::Boson<phys::Particle>> // P = Jet or Particle
		const P*        closestSing(std::vector<P>* cands, const R& ref, size_t& k); //max dR=0.4
		template <class P, class R = phys::Boson<phys::Particle>> // P = Jet or Particle
		const P*        closestSing(std::vector<P>* cands, const R& ref){
			size_t k = 0;
			return(closestSing(cands, ref, k));
		};
		template <class P, class R = phys::Boson<phys::Particle>> // P = Jet or Particle
		phys::Boson<P>* closestPair(std::vector<P>* cands, const R& ref); //max dR = 0.4
		template <class P, class R = phys::DiBoson<phys::Lepton, phys::Lepton>>
		P* furthestSing(std::vector<P>* cands, const R& reference, const float& minDR = 2., const std::pair<float,float>& massLimits = std::make_pair(1.,1000.));
		template <class P, class R = phys::DiBoson<phys::Lepton, phys::Lepton>>
		phys::Boson<P>* furthestPair(std::vector<P>* cands, const R& reference, const float& minDR = 2., const std::pair<float,float>& massLimits = std::make_pair(1.,1000.));
		
		template <class P = phys::Particle>
		static inline double getRefinedMass(const P& p){ return p.mass(); } 
		// If I used "const P&" it would have precedence on the template specialization when the object is a pointer to non-const, since "const P&" matches const P* & (const reference to non-const pointer to non-const)
		//template <class P = phys::Particle*> 
		//inline double getRefinedMass(P* p) const{ return p->mass();}
		
		//I give up, c++ templates are too complicated when mixed with cv-qualifiers. It's easier to write the function overloads
		static inline double getRefinedMass(const phys::Particle* p){ return p->mass(); }
		template <class P = phys::Particle>
		static inline double getRefinedMass(const phys::Boson<P>* p){ return p->mass(); }
		
		//Function overload of the template for Jets 
		static inline double getRefinedMass(const phys::Jet& j) { return j.chosenAlgoMass(); }
		static inline double getRefinedMass(const phys::Jet* j) { return j->chosenAlgoMass();}
		
		static inline double minDM(const double& mass, const double& r1 = phys::ZMASS, const double& r2 = phys::WMASS) { return std::min( fabs(mass-r1), fabs(mass-r2) ); }
		
		// ----- ----- Predictions from scikit classifiers ----- ----- 
		double getPyPrediction(const std::vector<double>&, const PyObject*) const; // uses module_helper
		static std::vector<double> getAK4features(const phys::Boson<phys::Jet>&);
		static std::vector<double> getAK8features(const phys::Jet&);
		
		// ----- ----- Event-specific variables calculation ----- ----- 
		void fillGenHadVBs(); //old: Fills the vector only if it is empty
		void fillRecHadVBs(); //old: Fills the vector only if it is empty
		void calcS();
		void fillAK4GenRec(bool doGraphs = true); // Fills AK4GenRec_
		void fillGenVBtoAK4();
		void makeGenZZ();
		
		// ----- ----- Large sub-analisys ----- ----- 
		void baseHistos();  // called in cut(), so runs every event
		void simpleGraphs();
		void jetRecoGraphs();
		
		void bestCandidateAnalysis();
		void endBestCandAnalysis(TFile& f_out); //Divide histograms to obtain efficiency
		
		void ptCutMVA(); //select the best cut in pt for the jetsAK8 --> corrPrunedMass
		
		void closestJetAnalisys(); //how often ak8 are better than ak4? Does it depend on pt? |p|? Mass?
		void endClosestJetAn(); //Histogram normalization
		
		void furthestJetMVA(); //Looking for a variable that lets us choose whether to use an AK4 or an AK8
		void minPtJetMVA(); //Still looking. Maybe minimizing total pt can help us choose AK4/8
		
		void bestZMassJetMVA();
		void specialPeakAnalisys(const phys::Particle& theGenAK8); //Is there a pair of AK4 that reconstructs theese events with an AK8 but low-pt ZZ?
		void endResolutionAnalisys(TFile& f_out); //Calculates AK4-AK8 resolution per bin of ZZ pt
		
		void resolutionZmass();  // same as bestZMassJetMVA() but without dR(jet, ZZ) > 2.
		void endResolutionZmass(TFile& f_out);
		//Permanently moved into an external macro
		
		void AK8nearGenHadVB();  // Is there an AK8 near genHadVBs_->front() ?
		
		void makeQQ(bool doGraphs = false); //phys::Boson<phys::Particle>* genQuarksID(); // Returns a boson created wit the 2 quarks only if the final state is 4l 2q
		void genQuarksAnalisys(/*const phys::Boson<phys::Particle>* qq_*/);  // gen quarks from VB decay
		//void endGenQuarksAnalisys(TFile& f_out);  //better implementation in macro/Efficiency
		
		//Analisys on the resolution/efficiency of reconstruction of AK4 and AK8
		void reconstructionAK();  //runs reconstructionAK4() and 8() which calculate efficiency
		std::pair<const phys::Boson<phys::Particle>*, phys::Boson<phys::Jet>*> reconstructionAK4(); // returns nullptr(s) if there was no reconstructed/generated candidate
		std::pair<const phys::Particle*, phys::Jet*> reconstructionAK8(/*bool doGraphs=false*/);
		void AKrace(std::pair<const phys::Boson<phys::Particle>*, phys::Boson<phys::Jet>*>, std::pair<const phys::Particle*, phys::Jet*>);
		void AK8recover();  // Is it possible to recover some signal by using AK8?
		void endReconstructionAK();
		
		void AK8MassAlgorithms();
		
		void genTauAnalisys();
		
		void genSignalGraphs();  // Reads sigType_
		void recSignalGraphs();
		void recSignalAnalysis();
		void endSignalEff(TFile&);
		
		void genHadVBCategorization();  // Two genHadVB (W2, Z3) can be formed from the same jets. How often does that happen?
		void endGenHadVbCateg();
		
		void endNameCuts();  // Gives names to the xAxis labels of "Cuts"
		
		void endVHad_vs_ZZpt(TFile&);  // Stack plot of "Vjj_ZZpt" and "VJ_ZZpt"
		
	private:
		float sAK4g_; //invariant mass of ZZ and all the AK4s gen
		float sAK8g_; //                                 AK8s
		float sAK4r_; //                                 AK4s rec
		float sAK8r_; //                                 AK8s
		// We will use #hat{s} for proper energy of ZZ+best AK8/pair pf AK4
		
		// ----- ----- Counters, ecc. ----- ----- 
		unsigned int singWFromJets_, pairWFromJets_, singWFromJetsAK8_, pairWFromJetsAK8_;
		unsigned int singZFromJets_, pairZFromJets_, singZFromJetsAK8_, pairZFromJetsAK8_;
		unsigned int recVBtot_; //Evts where the reconstructed VB is acceptable (VBosonDefinition)
		unsigned int goodRec_, withGenVB_; //Evts wit a recVB near a genVB / Evts with a genVB
		unsigned int win4_,win8_;//How often an AK4/AK8 reconstructs better an hadronic decaying VB
		
		unsigned int Nhad_genVB_;  // Generated VB with hadronic decay
		unsigned int Ncms_recVB_;  // Reconstructed by the detector
		unsigned int Nall_recVB_60_120_;  // Reconstructed by the algorithm (60<mrec<120)
		unsigned int Ngood_recVB_60_120_; // Rec. by the algorithm that match the rec. from the detector 
		//in SignalDefinitions.cc for generated VB (bestJetPairW, and bestJetPairZ): (|m-mV|<30) 
		unsigned int Nall_recVB_m20_;  // Reconstructed by the algorithm (|m-mV|<20) 
		unsigned int Ngood_recVB_m20_;
		
		unsigned int N_gen8_, Ncms_rec8_;  // Implicit cut 60<m_gen<120
		unsigned int Ncms_rec8_pTau35_, Nall_rec8_pTau35_; //Reconstructed with a puppiTau21 < 0.35
		unsigned int Ngood_rec8_pTau35_;  //Reconstructed near a recAK8 that is close to a genAK8
		
		unsigned int NnoAK4_, N8genVB_, N8recover_;  // See AK8recover()
		
		friend class Selector<VZZAnalyzer>;
		
	protected:
		unsigned long evtN_, analyzedN_; //Used to count processed events.
		clock_t startTime_; //Used to calculate elapsed time
		float analyzedW_;  // Weighted events passing cut
		
		PyObject* helper_module_;  // Python mini-module to unpickle AK4_classifier_ and obtain predictions from it
		PyObject* AK4_classifier_;  // A scikit-learn classifier, trained on pairs of AK4 to distinguish from W/Z induced jets and background QCD processes. Must implement the method "predict_proba(<2D matrix>)"
		PyObject* AK8_classifier_;
		
		//phys::DiBoson<phys::Particle, phys::Particle>* ZZ_gen = nullptr;
		std::vector<phys::Boson<phys::Particle>>* genHadVBs_ = nullptr;  // genVBParticles with hadronic daugthers
		std::vector<phys::Boson<phys::Particle>>* AllGenVBjj_ = nullptr;  // All the unique pairs of genAK4 with |m - M_{Z,W}| < 30.  genVB is limited to 2 (Z3 and W2 in SignalDefinitions)
		std::vector<phys::Boson<phys::Jet>>* AK4pairs_ = nullptr;  // All the pairs of all the reconstructed AK4 jets
		std::vector<std::pair<phys::Boson<phys::Particle>, phys::Boson<phys::Jet>>>* AK4GenRec_ = nullptr;  // Pairs of gen VB-->jj succesfully reconstructed by the detector
		
		phys::DiBoson<phys::Particle, phys::Particle>* genZZ_; //Always != nullptr from begin() to end().
		phys::Boson<phys::Particle> qq_;
		
		// ----- ----- Signal definition ----- ----- 
		int sigType_;            // 0-->no  |     1-->AK4     |  2-->AK8
		phys::Particle* sigVB_;  // nullptr | Boson<Particle> | Particle
		int isSignal(bool doGraphs = false);  
		
		template <class PAR>
		bool ZBosonDefinition(phys::Boson<PAR>& cand) const{  //candidate
			//bool checkMass = fabs(cand.p4().M() - phys::ZMASS) < 40; //temp
			bool checkMass = (60. < cand.p4().M() && cand.p4().M() < 120.);
			return checkMass;
		}
		template <class PAR>
		bool WBosonDefinition(phys::Boson<PAR>& cand) const{	//candidate
			//bool checkMass = fabs(cand.p4().M() - phys::WMASS) < 40;
			bool checkMass = (60. < cand.p4().M() && cand.p4().M() < 120.);
			return checkMass;
		}
		
		template <class PAR>  //usually PAR is Jet
		bool ZBosonDefinition(PAR& cand) const{
			//bool checkMass = fabs(getRefinedMass(cand) - phys::ZMASS) < 40;
			float mass = getRefinedMass(cand);
			bool checkMass = (60. < mass && mass < 120.);
			return checkMass;
		}
		template <class PAR>	
		bool WBosonDefinition(PAR& cand) const{	
			//bool checkMass = fabs(getRefinedMass(cand) - phys::WMASS) < 40;
			float mass = getRefinedMass(cand);
			bool checkMass = (60. < mass && mass < 120.);
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


/*template <class T>
class SortableMap{
	public:
		SortableMap() {};
		SortableMap(const std::vector<T>& vect){
			for(size_t i = 0; i < vect.size(); ++i)
				map.push_back(make_pair(i, vect.at(i)));
		}
		
		void sort_by_val(bool comp(const T&, const T&)){
			auto lambda = [comp](const std::pair<int, T>& a, const std::pair<int, T>& b) { return comp(a.second, b.second); };
			sort(map.begin(), map.end(), lambda);
		}
		
		size_t find_by_key(int k) const{
			for(size_t i = 0; i < map.size(); ++i)
				if(map.at(i).first == k)
					return i;
		}
		
		template <class P = phys::Particle>
		size_t find_closest(const P& reference) const{
			//closestSing
			if(map.size() < 1)  return -1;
			if(map.size() == 1) return 0;
			
			int index = -1;
			float minDR = 0.4; //starting value = the threshold we use
			for(int i = 0; i < map.size(); ++i){
			float DR = physmath::deltaR( map.at(i).second, reference );
				if(DR < minDR){
					minDR = DR;
					index = i;
				}
			}
			if(physmath::deltaR(reference, map.at(index).second) < 0.4)
				return index;
			else return nullptr;
		}
		
		std::vector<std::pair<int, T>> map;
};*/


#endif
