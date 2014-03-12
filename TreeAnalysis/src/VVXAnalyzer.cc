#include "VVXAnalysis/TreeAnalysis/interface/VVXAnalyzer.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using std::cout;
using std::endl;


using namespace phys;

void VVXAnalyzer::analyze(){

  std::vector<const Particle* > Z;
  foreach(const Boson<Lepton>& z, *Zmm)
    Z.push_back(&z);

  foreach(const Boson<Electron>& z, *Zee)
    Z.push_back(&z);

  std::stable_sort(Z.begin(),Z.end(),MassComparator(ZMASS));

}
