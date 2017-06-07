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


#endif
