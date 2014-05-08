#include "VVXAnalysis/TreeAnalysis/interface/VVXAnalyzer.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using std::cout;
using std::endl;


using namespace phys;

Int_t VVXAnalyzer::cut() {
  
  int nCentralJets30 = 0;
  foreach(const Jet &jet, *jets){
    if(fabs(jet.eta()) < 2.5 && jet.pt() > 30) ++nCentralJets30;
  }
  
  if(Zmm->size()+Zee->size() < 2) return -1;

  theHistograms.fill("nCentralJets30", "nCentralJets30",  10, 0, 10, nCentralJets30, theWeight); 
  if(nCentralJets30 >= 2){
    theHistograms.fill("deltaEtaJJ", "#Delta #eta(j,j)",  50, 0, 8, fabs(jets->at(0).eta() - jets->at(1).eta()), theWeight); 
    theHistograms.fill("mJJ", "m_{jj}",  100, 0, 1000, (jets->at(0).p4() + jets->at(1).p4()).M(), theWeight); 
  }
  if(nCentralJets30 >= 1) return 1;
  else return -1;

}

void VVXAnalyzer::analyze(){

   foreach(const Boson<Lepton>& z, *Zmm)
     theHistograms.fill("ZCandMass"    , "ZCandMass", 200, 0, 500, z.p4().M(), theWeight);
   foreach(const Boson<Electron>& z, *Zee)
     theHistograms.fill("ZCandMass"    , "ZCandMass", 200, 0, 500, z.p4().M(), theWeight);

   
   
}



