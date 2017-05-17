#include "VVXAnalysis/TreeAnalysis/interface/WlllnuAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp>
using namespace boost::assign;

using std::cout;
using std::endl;


using namespace phys;

Int_t WlllnuAnalyzer::cut() {
  return 1;
}


void WlllnuAnalyzer::analyze(){


  cout << "------------------------------------------------------------------"<<endl;
  cout << "Run: " << run << " event: " << event << endl;
  
  //Leptons
  
  std::vector<phys::Particle>  leptons;
  //  phys::Lepton<phys::Particle> lep; //why doesn't work? error: ‘phys::Lepton’ is not a template
  
  foreach(const phys::Particle &gen, *genParticles){
    if( (abs(gen.id()) != 11 && abs(gen.id()) != 13) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) continue;
    cout << "id: " << gen.id() << " pt: " << gen.pt() << " mass: " << gen.mass() << endl;
    theHistograms.fill("ptAllGenParticle","pt ", 100, 0, 100, gen.pt());
    theHistograms.fill("etaAllGenParticle","eta ", 100, 0, 100, gen.eta());
    theHistograms.fill("massAllGenParticle","mass ", 100, 0, 100, gen.mass());
    leptons.push_back(gen);
  }
  //cout << "lepton 0: mass" << leptons.at(0).mass() << endl;

  
  //----------------------------------------------------------------//  
  
  //Bosons
    
  std::vector<Boson<phys::Particle> > Zcandidates;
  
  phys::Boson<phys::Particle> z1;
  phys::Boson<phys::Particle> z2;
  double zMass = 91.1876; //Gev
  
  foreach(const phys::Boson<phys::Particle> &gen, *genVBParticles){
    if(gen.id() == 23 && (abs(gen.daughter(0).id()) == 11 || abs(gen.daughter(0).id()) == 13)){
      if(z1.id() != 23) z1 = gen;
      else if ( abs(z1.mass() - zMass) < abs(gen.mass() - zMass)) z2 = gen;
      else {
	z2 = z1;
	z1 = gen;
      }
      theHistograms.fill("ptAllZ","pt ", 100, 0, 100, gen.pt());  
    }
  }
  cout << "z1: " << z1.id() << " pt: " << z1.pt() << " mass: " << z1.mass() << endl;
  cout << "z2: " << z2.id() << " pt: " << z2.pt() << " mass: " << z2.mass() << endl;
  
  
  
}
