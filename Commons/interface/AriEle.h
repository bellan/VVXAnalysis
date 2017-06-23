#ifndef VVXAnalysis_Commons_AriEle_h
#define VVXAnalysis_Commons_AriEle_h

#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

using namespace phys;
using namespace std;

typedef std::pair<phys::Boson<phys::Particle>, phys::Particle> Zltype;
typedef DiBoson<Particle, Particle> ZZtype;
typedef Boson<Particle> Ztype;
typedef pair<Particle, Particle> pairParticle;

template<typename T> double mT(const T& p1, const T& p2){
  return sqrt( 2*p1.pt()*p2.pt()*(1-TMath::Cos(physmath::deltaPhi(p1.phi(), p2.phi()))) );
}
template<typename T> double mT(const T& p1, const T& p2, const T& p3){
  return sqrt(mT(p1, p2)*mT(p1, p2) + mT(p1, p3)*mT(p1, p3) + mT(p2, p3)*mT(p2, p3));
}


struct mTComparator{
mTComparator(const double& ref): ref_(ref){}
  template<typename BOS>
  bool operator()(const BOS & a , 
		  const BOS & b) const{ 
    return fabs(mT(a.daughter(0), a.daughter(1))-ref_) < fabs(mT(b.daughter(0), b.daughter(1))-ref_); 
  }
  
  template<typename PAR>
  bool operator()(const std::pair<PAR, PAR> & a , 
		  const std::pair<PAR, PAR> & b) const{ 
    return fabs(mT(a.first, a.second)-ref_) < fabs(mT(b.first, b.second)-ref_); 
  }
  /*    template<typename PAR>
	bool operator()(const PAR * a , 
	const PAR * b) const{ 
	return fabs(a->p4().M()-ref_) < fabs(b->p4().M()-ref_); 
	}*/
  double ref_;
};
  
struct massComparator{
massComparator(bool element, const double& ref): element_(element), ref_(ref){}
  template<typename T1, typename T2>
    bool operator()(const std::pair<T1, T2> & a, const std::pair<T1, T2> & b) const{    
    if(element_ == 0) return fabs(a.first.p4().M()-ref_) < fabs(b.first.p4.M()-ref_);
    else if (element_ == 1) return fabs(a.second.p4().M()-ref_) < fabs(b.second.p4.M()-ref_);    
  }

  bool element_;
  double ref_;
};
  
template<typename T> bool isTheSame(const T& p1, const T& p2){
  return !(p1.pt() != p2.pt());
  //return abs(p1.pt()- p2.pt()) < 0.1;
}


#endif
