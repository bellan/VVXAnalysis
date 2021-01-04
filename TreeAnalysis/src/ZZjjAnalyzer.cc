#include "VVXAnalysis/TreeAnalysis/interface/ZZjjAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Colours.h"

#include "VVXAnalysis/Commons/interface/GenVBHelper.h"


#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using namespace colour;
using namespace physmath;
using namespace phys;
using namespace std;


void ZZjjAnalyzer::GenAnalysis(DiBosonParticle &ZZ, Particle &Jet0, Particle &Jet1){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~ Begin of gen Analysis ~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Cuts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~ Histograms ~~~~~~~~~~~~~~~~~~~~~~~~~
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~ End of gen Analysis ~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}

void ZZjjAnalyzer::RecoAnalysis(DiBosonLepton &ZZ, Particle &Jet0, Particle &Jet1){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~ Begin of reco Analysis ~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Cuts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~ Histograms ~~~~~~~~~~~~~~~~~~~~~~~~~

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~ End of reco Analysis ~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}


void ZZjjAnalyzer::GenRecoAnalysis(const DiBosonParticle genZZ, const Particle genJet0, const Particle genJet1, const DiBosonLepton recoZZ, const Particle recoJet0, const Particle recoJet1){
  //eventGenReco++;
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~ Begin of Reco vs Gen ~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~ Histograms ~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~ End of Reco vs Gen ~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}


void ZZjjAnalyzer::begin(){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~ Begin function ~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  cout << "\n--------------------------------------------------------------------------" << endl;
  cout << "\n Starting ZZAnalysis on sample in " << fileName << endl;

  // event counters
  //eventSampe = 0;
  //eventGen = 0;
  //eventReco = 0;
  //eventGenReco = 0;
  
  // free counters
  
  
  // time begins
  begintime = ((float)clock())/CLOCKS_PER_SEC;
}


void ZZjjAnalyzer::analyze(){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~ Main analysis ~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  //eventSample++;
  
  //gen variables	
  DiBosonParticle genZZ;
  Particle genJet0;
  Particle genJet1;
  
  //reco variables
  DiBosonLepton recoZZ;
  Particle recoJet0;
  Particle recoJet1;
  
  //Gen analysis
  ZZjjAnalyzer::GenAnalysis(genZZ, genJet0, genJet1);
  
  //Reco analysis
  //if(genZZ.pt() != 0){
    ZZjjAnalyzer::RecoAnalysis(recoZZ, recoJet0, recoJet1);
  //}

  ///*
  if(genZZ.pt() != 0. && recoZZ.pt() != 0.){
    //Reco vs Gen analysis
    ZZjjAnalyzer::GenRecoAnalysis(genZZ, genJet0, genJet1, recoZZ, recoJet0, recoJet1);
  }
  //*/
}


Int_t ZZjjAnalyzer::cut(){  
  return 1;
}


void ZZjjAnalyzer::end(TFile &){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~ End banner ~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout << "\n--------------------------------------------------------------------------" << endl;
  
  //cout << "\nEvents of the sample analyzed:                       " << setw(9) << eventSample << endl;
  //cout << "Gen events analyzed:                                 " << setw(9) << eventGen << endl;
  //cout << "Reco events analyzed:                                " << setw(9) << eventReco << endl;
  //cout << "Gen&Reco events analyzed:                            " << setw(9) << eventGenReco << endl;

    
  // execution time
  endtime = ((float)clock())/CLOCKS_PER_SEC;
  helper_->printTime(begintime, endtime);
  cout << "\n--------------------------------------------------------------------------" << endl;
}
