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
  
  foreach(const Jet& jet1, *jets)
    foreach(const Jet& jet2, *jets){
    double mjj = (jet1.p4() + jet2.p4()).M();
    if (mjj > 150) return 1;
  }
  int countMu20 = 0;
  int countEle20 = 0;
  foreach(const Lepton& mu, *muons)    countMu20  += mu.pt() > 20 ? 1: 0;
  foreach(const Lepton& e, *electrons) countEle20 += e.pt() > 20 ? 1: 0; 
  if(countEle20 + countMu20 >= 2) return 1;
  
  return -1;
}



void VVXAnalyzer::analyze(){

  cout << "------------------------------------------------------------------"<<endl;
  cout << "Run: " << run << " event: " << event << endl;
  
  cout << "Gen Particles"<<endl;
  std::vector<phys::Particle> leptons;
  
  foreach(const phys::Particle &gen, *genParticles){
    bool isPrompt = gen.genStatusFlags().test(phys::GenStatusBit::isPrompt) && gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess);
    if(!isPrompt) cout << "Non prompt - ";
    // cout << "id: " << gen.id() << " pt: " << gen.pt() << " eta: " << gen.eta() << " phi: " << gen.phi() << endl;   
    if((abs(gen.id()) >= 11 && abs(gen.id()) <= 16 || gen.id() == 22) && gen.pt() > 20)
      leptons.push_back(gen);
    //theHistogram.Fill("GenJetsEta", "GenJetsEta", 100, -5, 5, gen.eta(), 1);
  }
  
  cout << "Gen Jets"<<endl;
  foreach(const phys::Particle &gen, *genJets){
    cout << "id: " << gen.id() << " pt: " << gen.pt() << " eta: " << gen.eta() << " phi: " << gen.phi() << endl;   
    theHistograms.fill("GenJetsEta", "GenJetsEta", 100, -5, 5, gen.eta(), theWeight);
  }

  //cout << "Gen Jets"<<endl;

  foreach(const phys::Particle &gen, *genJets){
    bool matched = false;
    foreach(const phys::Particle &lep, leptons)
      if(physmath::deltaR(lep,gen) < 0.4) matched = true;
    
    if(!matched && gen.pt() > 30)
      theHistograms.fill("GenCleanedJetsEta", "GenCleanedJetsEta", 100, -5, 5, gen.eta(), theWeight);
  }

  theHistograms.fill("njets", "njets", 10, 0, 10, jets->size(), theWeight);
  foreach(const phys::Jet &jet, *jets){
    cout << "pt: " << jet.pt() << " eta: " << jet.eta() << " phi: " << jet.phi() << endl;   
    theHistograms.fill("jetsEta", "jetsEta", 100, -5, 5, jet.eta(), theWeight);
    theHistograms.fill("jetsPt", "jetsPt", 100, 0, 500, jet.pt(), theWeight);
  }

  //cout << "Gen Jets"<<endl;

  
  
}
