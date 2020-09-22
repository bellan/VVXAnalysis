#include "VVXAnalysis/TreeAnalysis/interface/VVXAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Colours.h"

#include "VVXAnalysis/Commons/interface/GenVBHelper.h"


#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;
using namespace colour;

using namespace phys;


Int_t VVXAnalyzer::cut() {
  
  return 1;
}

void VVXAnalyzer::analyze(){

  theHistograms.fill("nZtoChLep"    , "Number of Z->ll per event" , 7,0,7, genVBHelper_.ZtoChLep().size());
  theHistograms.fill("nZtoNeutrinos", "Number of Z->nn per event" , 7,0,7, genVBHelper_.ZtoNeutrinos().size());
  theHistograms.fill("nWtoLep"      , "Number of W->lnu per event", 7,0,7, genVBHelper_.WtoLep().size());
  theHistograms.fill("nZtoQ"        , "Number of Z->qq per event" , 7,0,7, genVBHelper_.ZtoQ().size());
  theHistograms.fill("nWtoQ"        , "Number of W->qq' per event", 7,0,7, genVBHelper_.WtoQ().size());

  int nVBs = genVBHelper_.ZtoChLep().size() + genVBHelper_.ZtoNeutrinos().size() + genVBHelper_.WtoLep().size() + genVBHelper_.ZtoQ().size() + genVBHelper_.WtoQ().size();
  theHistograms.fill("nVBs", "Number of VB per event", 7,0,7, nVBs);
}















