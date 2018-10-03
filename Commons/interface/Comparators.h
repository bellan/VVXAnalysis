#ifndef VVXAnalysis_Commons_Comparators_H
#define VVXAnalysis_Commons_Comparators_H

#include "VVXAnalysis/Commons/interface/Utilities.h"

namespace phys{
	
	struct EComparator{
		template<typename PAR>
		bool operator()(const PAR& a, const PAR& b) const{
			return a.e() > b.e();
		}
	};
  
	struct PtComparator{
		template<typename LEP>
		bool operator()( const LEP & a, const LEP & b) const{ 
			return a.pt() > b.pt(); 
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
			if       (abs(std::get<2>(a)) == id_ && abs(std::get<2>(b)) != id_) return true;
			else if  (abs(std::get<2>(a)) != id_ && abs(std::get<2>(b)) == id_) return false;
			else if  (abs(std::get<2>(a)) == id_ && abs(std::get<2>(b)) == id_)
				return fabs(std::get<3>(a) - ref_) < fabs(std::get<3>(b) - ref_); 
			else    
				return true;
		}
		
		int id_;
		double ref_;
  };
  
	struct VdeltaRComparator{
		template<typename PAIR>
		bool operator()(const PAIR & a,
                    const PAIR & b) const{
		return physmath::deltaR(a.first, a.second) < physmath::deltaR(b.first, b.second);
		}
	};
  
}
#endif
