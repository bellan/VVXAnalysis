#include "VVXAnalysis/TreeAnalysis/interface/VVXAnalyzer.h"

#include "VVXAnalysis/TreeAnalysis/interface/Colours.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace std;

void VVXAnalyzer::analyze(){
  Colour cc("\033[0;41m");
  cout<<cc(3)<<45<<cc("TROTA")<<" " <<cc(9.94)<<endl;
  
}
