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
  
 
}
  
