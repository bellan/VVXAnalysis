#ifndef VVXAnalysis_Commons_Comparators_H
#define VVXAnalysis_Commons_Comparators_H

#include <algorithm>  // std::min()
#include "VVXAnalysis/Commons/interface/Utilities.h"

namespace phys{
	
	// ----- Energy
	struct EComparator{
		template<typename PAR>
		bool operator()(const PAR& a, const PAR& b) const{
			return a.e() > b.e();
		}
	};
  
  	
  	// ----- Pt
	struct PtComparator{
		template<typename LEP>
		bool operator()( const LEP & a, const LEP & b) const{ 
			return a.pt() > b.pt(); 
		}
		
		template<typename DiBOS1, typename DiBOS2>
		bool operator()(const DiBOS1 & a , const DiBOS2 & b) const{
			return a.first().daughter(0).pt() > b.first().daughter(0).pt();
		}
		
		template<typename DiBOS1, typename DiBOS2>
    		bool operator()(const DiBOS1 * a , const DiBOS2 * b) const{
			return a->first().daughter(0).pt() > b->first().daughter(0).pt();
		}
	};
  
	struct ScalarSumPtComparator{
		template<typename PAR>
		bool operator()(const PAR & a, const PAR & b) const{ 
			double sumPta = a.daughter(0).pt() + a.daughter(1).pt();
			double sumPtb = b.daughter(0).pt() + b.daughter(1).pt();
			return sumPta > sumPtb; 
		}
	};
  
	struct WJetPtComparator{
		template<typename BOS>
		bool operator()( const BOS & w1, const BOS & w2) const{
			double w1pt1 = w1.daughter(0).pt();
			double w1pt2 = w1.daughter(1).pt();
			if (w1pt2>w1pt1) std::swap(w1pt1,w1pt2);
			double w2pt1 = w2.daughter(0).pt();
			double w2pt2 = w2.daughter(1).pt();
			if (w2pt2>w2pt1) std::swap(w2pt1,w2pt2);
      
			if (w1pt2==w2pt2) {
				return w1pt1>w2pt1;
			} else return (w1pt2>w2pt2);
		}
	};
	
	struct PtTotRefComparator{ //Compares the pt of the sum of a particle's momentum and the one of a reference particle. Used to find the particle that minimizes total pt in an event
  	template <class P>
  	PtTotRefComparator(const P& refParticle) : refp4_(refParticle.p4()) {}
  	PtTotRefComparator(const TLorentzVector p4) : refp4_(p4) {}
  	
  	template <class P1, class P2>
  	bool operator()(const P1& p1, const P2& p2){
  		return (p1.p4() + refp4_).Pt() < (p2.p4() + refp4_).Pt();
  	}
  	
  	template <class P1, class P2>
  	bool operator()(const P1* p1, const P2* p2){
  		return (p1->p4() + refp4_).Pt() < (p1->p4() + refp4_).Pt();
  	}
  	
  	TLorentzVector refp4_;
  };
  
	
	// ----- Mass
	struct MassComparator{
		MassComparator(const double& ref): ref_(ref){}
		template<typename PAR>
		bool operator()(const PAR & a , const PAR & b) const{ 
			return fabs(a.p4().M()-ref_) < fabs(b.p4().M()-ref_); 
		}
		
		template<typename PAR>
		bool operator()(const PAR * a , const PAR * b) const{ 
			return fabs(a->p4().M()-ref_) < fabs(b->p4().M()-ref_); 
		}
    
    
		MassComparator(int id, const double& ref): id_(id), ref_(ref){}
    
		typedef std::tuple<uint,uint,int,double> QuarkPairFeatures;
		typedef std::vector<QuarkPairFeatures>   QuarkPairsFeatures;
        
		bool operator()(const QuarkPairFeatures & a , const QuarkPairFeatures & b) const{ 
			if       (abs(std::get<2>(a)) == id_ && abs(std::get<2>(b)) != id_) 
				return true;
			else if  (abs(std::get<2>(a)) != id_ && abs(std::get<2>(b)) == id_) 
				return false;
			else if  (abs(std::get<2>(a)) == id_ && abs(std::get<2>(b)) == id_)
				return fabs(std::get<3>(a) - ref_) < fabs(std::get<3>(b) - ref_); 
			else    
				return true;
		}
		
		int id_;
		double ref_;
  	};
  
  struct Mass2Comparator{ //Used to order Particles by the closest to any of 2 different mass values (e.g. the closest to ZMASS or WMASS)
		double ref1_, ref2_;
		Mass2Comparator(const double& ref1, const double& ref2): ref1_(ref1), ref2_(ref2){}
		
		template<typename PAR>
		bool operator()(const PAR & a, const PAR & b) const{ 
			double diffA = std::min(fabs(a.mass() - ref1_), fabs(a.mass() - ref2_));
			double diffB = std::min(fabs(b.mass() - ref1_), fabs(b.mass() - ref2_));
			return diffA < diffB;
		}
		
		bool operator()(const phys::Jet & a, const phys::Jet & b) const{ 
			double diffA = std::min(fabs(a.chosenAlgoMass() - ref1_), fabs(a.chosenAlgoMass() - ref2_));
			double diffB = std::min(fabs(b.chosenAlgoMass() - ref1_), fabs(b.chosenAlgoMass() - ref2_));
			return diffA < diffB;
		}
		
		template<typename PAR>
		bool operator()(PAR* a, PAR* b) const{ 
			return operator()(*a, *b);
		}
		
		bool operator()(const std::tuple<size_t,size_t,double> &a, const std::tuple<size_t,size_t,double>& b) const{
			double diffA = std::min(fabs(std::get<2>(a) - ref1_), fabs(std::get<2>(a) - ref2_));
			double diffB = std::min(fabs(std::get<2>(b) - ref1_), fabs(std::get<2>(b) - ref2_));
			return diffA < diffB;
		}
	};
	
	struct pairMassComparator{
		pairMassComparator(const int element, const double& ref): element_(element), ref_(ref){}
  		template<typename PAR, typename BOS>
    		bool operator()(const std::pair<BOS, PAR> & a , const std::pair<BOS, PAR> & b) const{
    			if(element_ == 0){
      				return fabs(a.first.p4().M() - ref_) < fabs(b.first.p4().M() - ref_); 
    			}
    			else{
      				return fabs(a.second.p4().M() - ref_) < fabs(b.second.p4().M() - ref_); 
    			}
  		}
  
  		template<typename PAR, typename BOS>
    		bool operator()(const std::pair<BOS, PAR> * a , const std::pair<BOS, PAR> * b) const{ 
    			if(element_ == 0){
      				return fabs(a->first.p4().M() - ref_) < fabs(b->first.p4().M() - ref_); 
    			}
    			else{
      				return fabs(a->second.p4().M() - ref_) < fabs(b->second.p4().M() - ref_); 
    			}
  		}
  
  		template<typename DiBOS>
    		bool operator()(const DiBOS & a , const DiBOS & b) const{ 
    			if(element_ == 0){
      				return fabs(a.first().p4().M() - ref_) < fabs(b.first().p4().M() - ref_); 
    			}
    			else{
      				return fabs(a.second().p4().M() - ref_) < fabs(b.second().p4().M() - ref_); 
    			}
  		}
  
  		template<typename DiBOS>
    			bool operator()(const DiBOS * a , const DiBOS * b) const{ 
    			if(element_ == 0){
      				return fabs(a->first().p4().M() - ref_) < fabs(b->first().p4().M() - ref_); 
    			}
    			else{
      				return fabs(a->second().p4().M() - ref_) < fabs(b->second().p4().M() - ref_); 
    			}
  		}
  
  		int element_;
  		double ref_;
	};
	
	
	// ----- Transverse mass
	struct mTComparator{
		mTComparator(const double& ref): ref_(ref){}
		template<typename BOS>
		bool operator()(const BOS & a , const BOS & b) const{ 
    			return fabs(physmath::mT(a.daughter(0), a.daughter(1))-ref_) < fabs(physmath::mT(b.daughter(0), b.daughter(1))-ref_); 
  		}
  
  		template<typename PAR>
  		bool operator()(const std::pair<PAR, PAR> & a , const std::pair<PAR, PAR> & b) const{ 
    			return fabs(physmath::mT(a.first, a.second)-ref_) < fabs(physmath::mT(b.first, b.second)-ref_); 
  		}
  
		double ref_;
	};
	
	struct pairTrmassComparator{
		pairTrmassComparator(const int element, const double& ref): element_(element), ref_(ref){}
  		template<typename PAR, typename BOS>
    		bool operator()(const std::pair<BOS, PAR> & a , const std::pair<BOS, PAR> & b) const{
    			if(element_ == 0){
      				return fabs(a.first.p4().Mt() - ref_) < fabs(b.first.p4().Mt() - ref_); 
    			}
    			else{
      				return fabs(a.second.p4().Mt() - ref_) < fabs(b.second.p4().Mt() - ref_); 
    			}
  		}
  
  		template<typename PAR, typename BOS>
    		bool operator()(const std::pair<BOS, PAR> * a , const std::pair<BOS, PAR> * b) const{ 
    			if(element_ == 0){
      				return fabs(a->first.p4().Mt() - ref_) < fabs(b->first.p4().Mt() - ref_); 
    			}
    			else{
      				return fabs(a->second.p4().Mt() - ref_) < fabs(b->second.p4().Mt() - ref_); 
    			}
  		}
  
  		template<typename DiBOS>
    		bool operator()(const DiBOS & a , const DiBOS & b) const{ 
    			if(element_ == 0){
      				return fabs(a.first().p4().Mt() - ref_) < fabs(b.first().p4().Mt() - ref_); 
    			}
    			else{
      				return fabs(a.second().p4().Mt() - ref_) < fabs(b.second().p4().Mt() - ref_); 
    			}
  		}
  
  		template<typename DiBOS>
    		bool operator()(const DiBOS * a , const DiBOS * b) const{ 
    			if(element_ == 0){
      				return fabs(a->first().p4().Mt() - ref_) < fabs(b->first().p4().Mt() - ref_); 
    			}
    			else{
      				return fabs(a->second().p4().Mt() - ref_) < fabs(b->second().p4().Mt() - ref_); 
    			}
  		}
  
  		int element_;
  		double ref_;
	};
  
	
	// ----- Delta R
	struct VdeltaRComparator{
		template<typename PAIR>
		bool operator()(const PAIR & a, const PAIR & b) const{
			return physmath::deltaR(a.first, a.second) < physmath::deltaR(b.first, b.second);
		}
		
		template<typename P>
		bool operator()(const Boson<P> & a, const Boson<P> & b) const{
			return physmath::deltaR(a.daughter(0), a.daughter(1)) < physmath::deltaR(b.daughter(0), b.daughter(1));
		}
	};
	
	struct deltaRComparator{		
		template<typename VV>
		bool operator()(const VV & a, const VV & b) const{
    			return physmath::deltaR(a.first.daughter(0), a.first.daughter(1)) < physmath::deltaR(b.first.daughter(0), b.first.daughter(1));
  		}
	};
	
	struct DeltaRComparator{ //Used to sort particles in order of deltaR from a reference momentum
  		template <class P>
  		DeltaRComparator(const P& refParticle) : refp4_(refParticle.p4()) {}
  		DeltaRComparator(const TLorentzVector p4) : refp4_(p4) {}
  	
  		template <class P1, class P2>
  		bool operator()(const P1& p1, const P2& p2){
  			return physmath::deltaR(p1.p4(), refp4_) < physmath::deltaR(p2.p4(), refp4_);
  		}
  		
  		template <class P1, class P2>
  		bool operator()(const P1* p1, const P2* p2){
  			return physmath::deltaR(p1->p4(), refp4_) < physmath::deltaR(p2->p4(), refp4_);
  		}
  	
  		TLorentzVector refp4_;
	};
	
	struct PairComparator{
		template <class T, class U>
		bool operator()(const std::pair<T, U>& a, const std::pair<T, U>& b){
			return a.second < b.second;
		}
	};
  
}
#endif
