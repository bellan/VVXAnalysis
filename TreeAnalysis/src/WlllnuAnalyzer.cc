#include "VVXAnalysis/TreeAnalysis/interface/WlllnuAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"

#include <TCanvas.h>
#include <TH1F.h>

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


template<typename T> double mT(const T& p1, const T& p2){

  return sqrt( 2*p1.pt()*p2.pt()*(1-TMath::Cos(physmath::deltaPhi(p1.phi(), p2.phi()))) );

}


void WlllnuAnalyzer::analyze(){
  
  cout << "------------------------------------------------------------------"<<endl;
  cout << "Run: " << run << " event: " << event << endl;
  
  //leptons
  
  std::vector<phys::Particle>  leptons;
  std::vector<Particle> electrons; // Z daughters
  std::vector<Particle> muons; // Z daughters  
  int finalid = 0;

   //find leptons
  foreach(const phys::Particle &gen, *genParticles){    
    if( (abs(gen.id()) != 11 && abs(gen.id()) != 13) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) continue;
    finalid += abs(gen.id());
    cout << " genLepton: " << gen << endl;
    //cout << "id: " << gen.id() << " pt: " << gen.pt() << " mass: " << gen.mass() << endl;
    theHistograms.fill("ptAllGenParticle","pt ", 100, 0, 200, gen.pt());
    theHistograms.fill("etaAllGenParticle","eta ", 100, -10, 10, gen.eta());
    theHistograms.fill("massAllGenParticle","mass ", 100, 0, 0.2, gen.mass());
    theHistograms.fill("YAllGenParticle","Y ", 100, 0, 100, gen.rapidity());
    leptons.push_back(gen);
    if (gen.id() == 11) electrons.insert(electrons.begin(),gen); //fist all e-, then all e+
    else if (gen.id() == -11) electrons.push_back(gen);
    else if (gen.id() == 13) muons.insert(muons.begin(),gen); //first all mu-, then all mu+
    else if (gen.id() == -13) muons.push_back(gen);
  }
  
  //check 4l in final state
  if( (electrons.size() + muons.size() != 4) || (finalid != 44 && finalid != 48 && finalid != 52) ){
    cout << "\n\tEvent without 4l" << endl;
    return;    
  }
  
  //build z0, z1
  std::vector<Boson<Particle> > Zcandidates;
  phys::Boson<phys::Particle> z0;
  phys::Boson<phys::Particle> z1;
    
  if (finalid == 44){ //4e
    Zcandidates.push_back(phys::Boson<Particle>(electrons[0],electrons[2], 23));
    Zcandidates.push_back(phys::Boson<Particle>(electrons[1],electrons[3], 23));
    Zcandidates.push_back(phys::Boson<Particle>(electrons[0],electrons[3], 23));
    Zcandidates.push_back(phys::Boson<Particle>(electrons[1],electrons[2], 23));
  }
  if (finalid == 52){ //4mu
    Zcandidates.push_back(phys::Boson<Particle>(muons[0],muons[2], 23));
    Zcandidates.push_back(phys::Boson<Particle>(muons[1],muons[3], 23));
    Zcandidates.push_back(phys::Boson<Particle>(muons[0],muons[3], 23));
    Zcandidates.push_back(phys::Boson<Particle>(muons[1],muons[2], 23));
  }
  if (finalid == 48){ //2e2mu
    Zcandidates.push_back(phys::Boson<Particle>(electrons[0],electrons[1], 23));
    Zcandidates.push_back(phys::Boson<Particle>(muons[0],muons[1], 23));
  }
  z0 = Zcandidates[0];
  foreach(const phys::Boson<phys::Particle> cand, Zcandidates)
    if ( abs(cand.mass() - ZMASS) < abs(z0.mass() - ZMASS)) z0 = cand;
  foreach(const phys::Boson<phys::Particle> cand, Zcandidates)
    if ( cand.daughter(0).pt() != z0.daughter(0).pt() &&  cand.daughter(1).pt() != z0.daughter(1).pt() ) z1 = cand;//camparing pt is not enough, need function to compare particles: maybe isAlmostEqual ??
  //!!!need to find a way to check if 2 particles are the same!!!!

  theHistograms.fill("massZParticles","Z mass ", 1000, 50, 150, z0.mass());
  theHistograms.fill("massZparticles","Z mass ", 1000, 50, 150, z1.mass());
  
  // cout << "\n z0: " << z0 << endl;
  cout << "\n z0.daughter(0).pt(): " << z0.daughter(0).pt() << " z0.daughter(1).pt(): " << z0.daughter(1).pt() << " z0 mass: " << z0.mass() << endl;
  //cout << " z1: " << z1 << endl;
  cout << " z1.daughter(0).pt(): " << z1.daughter(0).pt() << " z1.daughter(1).pt(): " << z1.daughter(1).pt() << " z1 mass: " << z1.mass() << endl;
    
  /*
  //comparing leptons with z0, z1: change it!!

  cout << "\n\t Comparing mT leptons and mass z0/z1" << endl;
  for(unsigned int i=0; i<leptons.size(); i++){
    for(unsigned int j=i++; j<leptons.size(); j++){
      if(leptons[i].id() == -leptons[j].id()){
	if(abs(leptons[i].id()) == abs(z0.daughter(0).id()) )
	  cout << i << " " << j << " z0: " << mT(leptons[i], leptons[j]) - z0.mass() << endl;
	else if(abs(leptons[i].id()) == abs(z1.daughter(0).id()) )
	  cout << i << " " << j << " z1: " << mT(leptons[i], leptons[j]) - z1.mass() << endl;
      }
      else cout << i << " " << j << " no Z daughter candididates" << endl;
    }  
  }
  */
  //deltaR, deltaEta, deltaPhi
  
  double deltaEta0 = z0.daughter(0).eta() - z1.eta(); //between z1 and z0's first daughter
  double deltaEta1 = z0.daughter(1).eta() - z1.eta(); //between z1 and z0's first daughter
  /* 
  cout << "\n\ndeltaR daughter 0: " << deltaR(z1.p4(), z0.daughter(0).p4()) << endl;
  cout << "deltaR daughter 1: " << deltaR(z1.p4(), z0.daughter(1).p4()) << endl; 
  cout << "deltaPhi daughter 0: " << deltaPhi(z0.daughter(0).p4(), z1.p4()) << endl;
  cout << "deltaPhi daughter 1: " << deltaPhi(z0.daughter(1).p4(), z1.p4()) << endl;
  cout << "deltaEta daughter 0: " << deltaEta0 << endl;
  cout << "deltaEta daughter 1: " << deltaEta1 << endl;
  */
  // ZZ
  
  DiBoson<phys::Particle,phys::Particle>  ZZ(z0,z1);
  cout << "\n\n ZZ: " << ZZ.id() << " pt: " << ZZ.pt() << " mass: " << ZZ.mass() << "\n daughters: " << ZZ.first().id() << "\t" << ZZ.second().id() << " Y: " << ZZ.rapidity()
       << endl; //why daughter(0) instead of first() is not working?? daughter<Particle>(0)
  
  theHistograms.fill("ptZZ","pt ", 100, 0, 100, ZZ.pt());
  theHistograms.fill("etaZZ","eta ", 100, 0, 100, ZZ.eta());
  theHistograms.fill("massZZ","mass ", 1000, 0, 400, ZZ.pt());
  theHistograms.fill("YZZ","Y ", 100, 0, 100, ZZ.rapidity());
  //<phys::Boson<phys::Particle>,phys::Boson<phys::Particle> >

  //----------------------------------------------------------------//
  
}
