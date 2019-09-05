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
   *-----------------GenParticles------------------
   *  
*/  

  
  int ngenElectrons=0;
  int ngenMuons=0;
  
  foreach(const phys::Particle& gen, *genParticles){
    if(abs(gen.id()) == 13 ){//for muons: mc pdg id = 13 by notation
      theHistograms.fill("genMuonPt","pt of gen muons",20,0,200,gen.pt(),theWeight);
      theHistograms.fill("genMuonEta","eta of gen muons",30,0,2.5,fabs(gen.eta()),theWeight);
      theHistograms.fill("genMuonPhi","phi of gen muons",30,-3.2,3.2,gen.phi(),theWeight);
      theHistograms.fill("genMuonRapidity","rapidity of gen muons",30,0,2.5,fabs(gen.rapidity()),theWeight);
      ngenMuons++;
    }


    
    if(abs(gen.id()) == 11 ){//for electrons: mc pdg id = 11 by notation
      theHistograms.fill("genElectronPt","pt of gen electrons",20,0,200,gen.pt(),theWeight);
      theHistograms.fill("genElectronEta","eta of gen electrons",30,0,2.5,fabs(gen.eta()),theWeight);
      theHistograms.fill("genElectronPhi","phi of gen electrons",30,-3.2,3.2,gen.phi(),theWeight);
      theHistograms.fill("genElectronRapidity","rapidity of gen electrons",30,0,2.5,fabs(gen.rapidity()),theWeight);
      ngenElectrons++;
    }
   
  }
  theHistograms.fill("ngenElectrons","Number of genElectrons",  10, 0, 10 , ngenElectrons, theWeight);
  theHistograms.fill("ngenMuons","Number of genMuons",  10, 0, 10 , ngenMuons, theWeight);  






  /*
   *
   *---------------MUONS------------------
   *
*/

  
  theHistograms.fill("nmuons","Number of muons",  10, 0, 10 , muons->size(), theWeight);

  
  foreach(const phys::Lepton& mu, *muons){
    theHistograms.fill("muonsPt","pt of muons",20,0,200,mu.pt());
    theHistograms.fill("muonsPtEta","eta vs pt of muons",60,0,200,30,0,2.5, mu.pt(),fabs(mu.eta()));
    theHistograms.fill("muonsPtRapidity","rapidity vs pt of muons",60,0,200,30,0,2.5, mu.pt(),fabs(mu.rapidity()));
    theHistograms.fill("muonsPhiEta","eta vs phi of muons",60,-3.2,3.2,30,0,2.5, mu.phi(),fabs(mu.eta()));
    theHistograms.fill("muonsPhiPt","pt vs phi of muons",60,-3.2,3.2,60,0,200, mu.phi(),mu.pt());

    
  }

  /*
  loop(const phys::Particle& gen, *genParticles){
    if(abs(gen.id()) == 13)

  } 

  */
  /*
   *
   *-----------------ELECTRONS------------------
   *  
*/  

  
  theHistograms.fill("nelectrons","Number of electrons",  10, 0, 10 , electrons->size(), theWeight);

  
   foreach(const phys::Lepton& e, *electrons){
    theHistograms.fill("electronsPt","pt of electrons",20,0,200,e.pt());
    theHistograms.fill("electronsPtEta","eta vs pt of electrons",60,0,200,30,0,2.5, e.pt(),fabs(e.eta()));
    theHistograms.fill("electronsPtRapidity","rapidity vs pt of electrons",60,0,200,30,0,2.5, e.pt(),fabs(e.rapidity()));
    theHistograms.fill("electronsPhiEta","eta vs phi of electrons",60,-3.2,3.2,30,0,2.5, e.phi(),fabs(e.eta()));
    theHistograms.fill("electronsPhiPt","pt vs phi of electrons",60,-3.2,3.2,60,0,200, e.phi(),e.pt());

    
  }

   theHistograms.fill("nelectronsmuons","Number of electrons and muons",  10, 0, 10 , 10, 0 , 10 , electrons->size(), muons->size(), theWeight);

   
   /*
   *
   *-----------------JETS------------------
   *  
*/  


  foreach(const phys::Particle& gen, *genJets){
      theHistograms.fill("genJetPt","pt of gen jets",20,0,200,gen.pt(),theWeight);
      theHistograms.fill("genJetEta","eta of gen jets",30,0,2.5,fabs(gen.eta()),theWeight);
      theHistograms.fill("genJetPhi","phi of gen jets",30,-3.2,3.2,gen.phi(),theWeight);
      theHistograms.fill("genJetRapidity","rapidity of gen jets",30,0,2.5,fabs(gen.rapidity()),theWeight);

      
  }

  theHistograms.fill("njets","Number of jets",  10, 0, 10 , jets->size(), theWeight);

  
  foreach(const phys::Jet& j, *jets){
    theHistograms.fill("jetsPt","pt of jets",20,30,200,j.pt());
    theHistograms.fill("jetsEta","eta of jets",30,0,4.7,fabs(j.eta()));
    theHistograms.fill("jetsPtEta","eta vs pt of jets",20,30,200,30,0,4.7, j.pt(),fabs(j.eta()));
    theHistograms.fill("jetsRapidity","rapidity of jets",30,0,5,fabs(j.rapidity()));
    theHistograms.fill("jetsPtRapidity","rapidity vs pt of jets",20,30,200,30,0,5, j.pt(),fabs(j.rapidity()));
    theHistograms.fill("jetsPhi","phi of jets",60,-3.2,3.2,j.phi());
    theHistograms.fill("jetsPtPhi","phi vs pt of jets",20,30,200,60,-3.2,3.2, j.pt(),j.phi());
    theHistograms.fill("jetsEtaPhi","phi vs eta of jets",30,0,4.7,60,-3.2,3.2,fabs(j.eta()),j.phi());




  }


  /*
   *
   *-------------------JETS(AK8)------------------
   *  
*/  

    foreach(const phys::Particle& gen, *genJetsAK8){
      theHistograms.fill("genJetAK8Pt","pt of gen jetsAK8",20,0,200,gen.pt(),theWeight);
      theHistograms.fill("genJetAK8Eta","eta of gen jetsAK8",30,0,2.5,fabs(gen.eta()),theWeight);
      theHistograms.fill("genJetAK8Phi","phi of gen jetsAK8",30,-3.2,3.2,gen.phi(),theWeight);
      theHistograms.fill("genJetAK8Rapidity","rapidity of gen jetsAK8",30,0,2.5,fabs(gen.rapidity()),theWeight);

      
  }


 foreach(const phys::Jet& j, *jetsAK8){
    theHistograms.fill("jetsAK8Pt","pt of jetsAK8",20,30,200,j.pt());
    theHistograms.fill("jetsAK8Eta","eta of jetsAK8",30,0,4.7,fabs(j.eta()));
    theHistograms.fill("jetsAK8PtEta","eta vs pt of jetsAK8",20,30,200,30,0,4.7, j.pt(),fabs(j.eta()));
    theHistograms.fill("jetsAK8Rapidity","rapidity of jetsAK8",30,0,5,fabs(j.rapidity()));
    theHistograms.fill("jetsAK8PtRapidity","rapidity vs pt of jetsAK8",20,30,200,30,0,5, j.pt(),fabs(j.rapidity()));
    theHistograms.fill("jetsAK8Phi","phi of jetsAK8",60,-3.2,3.2,j.phi());
    theHistograms.fill("jetsAK8PtPhi","phi vs pt of jetsAK8",20,30,200,60,-3.2,3.2, j.pt(),j.phi());
    theHistograms.fill("jetsAK8EtaPhi","phi vs eta of jetsAK8",30,0,4.7,60,-3.2,3.2,fabs(j.eta()),j.phi());

 }

}
