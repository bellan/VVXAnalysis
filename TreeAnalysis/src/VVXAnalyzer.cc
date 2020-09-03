#include "VVXAnalysis/TreeAnalysis/interface/VVXAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;


using namespace phys;

Int_t VVXAnalyzer::cut() {
  
  return 1;
}



void VVXAnalyzer::analyze(){

  cout << "------------------------------------------------------------------"<<endl;
  cout << "Run: " << run << " event: " << event << endl;

  cout << "# muons: " << muons->size() <<endl;

  theHistograms.fill("nmuons", "Number of muons per event", 6,0,6, muons->size());


  foreach(const phys::Lepton muon, *muons)
    cout << "pt muon: " << muon.pt() << endl;


  
  
  
}
