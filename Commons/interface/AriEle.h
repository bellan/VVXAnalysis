#ifndef VVXAnalysis_Commons_AriEle_h
#define VVXAnalysis_Commons_AriEle_h

#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

using namespace phys;
using namespace std;

typedef std::pair<phys::Boson<phys::Particle>, phys::Particle> Zltype;
typedef DiBoson<Particle, Particle> ZZtype;
typedef Boson<Particle> Ztype;

template<typename T> double mT(const T& p1, const T& p2){
  return sqrt( 2*p1.pt()*p2.pt()*(1-TMath::Cos(physmath::deltaPhi(p1.phi(), p2.phi()))) );
}

struct mTComparator{
  mTComparator(const double& ref): ref_(ref){}
    template<typename PAR>
    bool operator()(const PAR & a , 
		    const PAR & b) const{ 
      return fabs(mT(a.daughter(0), a.daughter(1))-ref_) < fabs(mT(b.daughter(0), b.daughter(1))-ref_); 
    }
  /*    template<typename PAR>
	bool operator()(const PAR * a , 
	const PAR * b) const{ 
	return fabs(a->p4().M()-ref_) < fabs(b->p4().M()-ref_); 
	}*/
   double ref_;
};

template<typename T> bool isTheSame(const T& p1, const T& p2){
  return !(p1.pt() != p2.pt());
  //return abs(p1.pt()- p2.pt()) < 0.1;
}


#endif
