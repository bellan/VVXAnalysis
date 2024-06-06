#include "VVXAnalysis/TreeAnalysis/interface/WlllnuAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include "VVXAnalysis/DataFormats/interface/TypeDefs.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"

#include <TCanvas.h>
#include <TH1F.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp>
using namespace boost::assign;
using namespace colour;

using std::cout;
using std::endl;


using namespace phys;
using namespace physmath;

void WlllnuAnalyzer::begin(){
}

Int_t WlllnuAnalyzer::cut() {
  return 1;
}

void WlllnuAnalyzer::analyze(){

  cout << "------------------------------------------------------------------"<<endl;
  cout << "Run: " << run << " event: " << event << endl;
  
  foreach(const phys::Particle &gen, *genParticles)
    theHistograms->fill("gen_ID"      , "gen ID"  , 20 , 0, 20, gen.id(), theWeight);
  
}


void WlllnuAnalyzer::end(TFile &){
}
