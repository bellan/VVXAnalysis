#include "VVXAnalysis/TreeAnalysis/interface/WZZAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;


using namespace phys;

Int_t WZZAnalyzer::cut() {
  
  return 1;
}



void WZZAnalyzer::analyze(){

  cout << "-----------------------------------------------------------------"<<endl;
  cout << "Run: " << run << " event: " << event << endl;
  
  /*
   *
   *------------INFORMATIONS ABOUT MUONS------------------
   *
*/  
  foreach(const phys::Particle& gen, *genParticles){
    if(abs(gen.id()) == 13 )//for muon: mc pdg id = 13 by notation
      theHistograms.fill("genMuonPt","pt of gen muons",20,0,200,gen.pt(),theWeight);
      theHistograms.fill("genMuonEta","eta of gen muons",30,0,2.5,fabs(gen.eta()),theWeight);
      theHistograms.fill("genMuonPhi","phi of gen muons",30,-3.2,3.2,gen.phi(),theWeight);
      theHistograms.fill("genMuonRapidity","rapidity of gen muons",30,0,2.5,fabs(gen.rapidity()),theWeight);

      
  }
   
  foreach(const phys::Lepton& mu, *muons){
    theHistograms.fill("muonsPt","pt of muons",20,0,200,mu.pt());
    theHistograms.fill("muonsPtEta","eta vs pt of muons",60,0,200,30,0,2.5, mu.pt(),fabs(mu.eta()));
    theHistograms.fill("muonsPtRapidity","rapidity vs pt of muons",60,0,200,30,0,2.5, mu.pt(),fabs(mu.rapidity()));
    
  }


   
  /*
   *
   *------------INFORMATIONS ABOUT ELECTRONS------------------
   *  
*/  

   foreach(const phys::Particle& gen, *genParticles){
    if(abs(gen.id()) == 11 )//for electron: mc pdg id = 11 by notation
      theHistograms.fill("genElectronPt","pt of gen electrons",20,0,200,gen.pt(),theWeight);
      theHistograms.fill("genElectronEta","eta of gen electrons",30,0,2.5,fabs(gen.eta()),theWeight);
      theHistograms.fill("genElectronPhi","phi of gen electrons",30,-3.2,3.2,gen.phi(),theWeight);
      theHistograms.fill("genElectronRapidity","rapidity of gen electrons",30,0,2.5,fabs(gen.rapidity()),theWeight);

  }
   
   foreach(const phys::Lepton& e, *electrons){
    theHistograms.fill("electronsPt","pt of electrons",20,0,200,e.pt());
    theHistograms.fill("electronsPtEta","eta vs pt of electrons",60,0,200,30,0,2.5, e.pt(),fabs(e.eta()));
    theHistograms.fill("electronsPtRapidity","rapidity vs pt of electrons",60,0,200,30,0,2.5, e.pt(),fabs(e.rapidity()));

  }
/*
   *
   *------------INFORMATIONS ABOUT JETS------------------
   *  
*/  


   /* foreach(const phys::Particle& gen, *pgenJets){
    if(abs(gen.id()) == 13 )//for muon: mc pdg id = 13 by notation
      theHistograms.fill("genMuonPt","pt of gen muons",20,0,200,gen.pt(),theWeight);
      theHistograms.fill("genMuonEta","eta of gen muons",30,0,2.5,fabs(gen.eta()),theWeight);

  }
   *//*
    
  foreach(const phys::Jet& mu, *jets){
    theHistograms.fill("muonsPt","pt of muons",20,0,200,mu.pt());
    theHistograms.fill("muonsPtEta","eta vs pt of muons",60,0,200,30,0,2.5, mu.pt(),fabs(mu.eta()));
  }

*/



}
