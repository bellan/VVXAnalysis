#include "VVXAnalysis/TreeAnalysis/interface/WlllnuAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp>
using namespace boost::assign;

using std::cout;
using std::endl;


using namespace phys;

Int_t WlllnuAnalyzer::cut() {
  return 1;
}


void WlllnuAnalyzer::analyze(){
  
  
  cout << "------------------------------------------------------------------"<<endl;
  cout << "Run: " << run << " event: " << event << endl;
  
  foreach(const phys::Particle &gen, *genParticles){
    if(abs(gen.id()) != 11 && abs(gen.id()) != 13 || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) continue;
       cout << "id: " << gen.id() << " pt: " << gen.pt() << endl;
       theHistograms.fill("ptAllGenParticle","pt ", 100, 0, 100, gen.pt());
  }
  
}
  
