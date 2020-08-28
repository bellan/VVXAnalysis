#include "VVXAnalysis/TreeAnalysis/interface/VZZaQGCAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;


using namespace phys;

Int_t VZZaQGCAnalyzer::cut() {
  
  return 1;
}



void VZZaQGCAnalyzer::analyze(){

  //cout << "------------------------------------------------------------------"<<endl;
  //cout << "Run: " << run << " event: " << event << endl;
  
  
}
