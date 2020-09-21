#include "VVXAnalysis/TreeAnalysis/interface/VVXAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Colours.h"


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

template<class T>
std::vector<phys::Boson<phys::Particle> > VVXAnalyzer::getZtoX(const std::vector<phys::Particle> & collectionX, T condition){

  std::vector<phys::Boson<phys::Particle> > ZtoX;
  std::vector<phys::Boson<phys::Particle> > ZtoXcand;
   
  foreach(const phys::Particle &p1, collectionX) foreach(const phys::Particle &p2, collectionX)
    if(condition(p1,p2)) ZtoXcand.push_back(phys::Boson<phys::Particle>(p1,p2,23));
  std::stable_sort(ZtoXcand.begin(), ZtoXcand.end(), phys::MassComparator(phys::ZMASS));
  if(!ZtoXcand.empty()){

    phys::Boson<phys::Particle> Z0 = ZtoXcand.at(0);
    phys::Boson<phys::Particle> Z1;

    std::stable_sort(ZtoXcand.begin(), ZtoXcand.end(), phys::ScalarSumPtComparator());
    
    foreach(const phys::Boson<phys::Particle> &z, ZtoXcand){
      if(!Z0.overlapWithDaughters(z.daughter(0)) && !Z0.overlapWithDaughters(z.daughter(1))){
	Z1 = z;
	break;
      }
    }  
    if(Z0.mass() >= 60 && Z0.mass() <= 120) ZtoX.push_back(Z0);
    if(Z1.mass() >= 60 && Z1.mass() <= 120) ZtoX.push_back(Z1);
  }

  foreach(const phys::Boson<phys::Particle> &zcan, ZtoXcand){
    bool match = false;
    foreach(const phys::Boson<phys::Particle> &zxx, ZtoX)
      if(zcan.overlapWithDaughters(zxx.daughter(0)) || zcan.overlapWithDaughters(zxx.daughter(1))) match = true;
    if(!match && zcan.mass() >= 60 && zcan.mass() <= 120) ZtoX.push_back(zcan);
  }
  return ZtoX;
}



void VVXAnalyzer::analyze(){

  // Warning. Inside genVBParticles there are ONLY bosons that satisfy the following criterions 
  // 1- build a ZZ and pass certain mass cuts (fully leptonic decay only!). 
  // 2- build a WZ and pass certain mass cuts (fully leptonic decay only!).
  // 3- an additional hadronically decaying VB is reconstructed using genJets
  // Taus are not stable particles therefore are inside a separate collection!
  
  // Prepare the containers
  std::vector<phys::Particle> quarks;
  std::vector<phys::Particle> chleptons, posleptons, negleptons; // chleptons is an UNCLEANED collection
  std::vector<phys::Particle> neutrinos, photons;
  std::vector<phys::Boson<phys::Particle> > ZtoChLep;
  std::vector<phys::Boson<phys::Particle> > ZtoNeutrinos;
  std::vector<phys::Boson<phys::Particle> > WtoLep;
  std::vector<phys::Boson<phys::Particle> > ZtoQ;
  std::vector<phys::Boson<phys::Particle> > WtoQ;
  

  // Divide Z from W that decay leptonically
  foreach(const phys::Boson<phys::Particle> &vb, *genVBParticles)
    if(vb.decayType() == 1)
      if     (    vb.id()  == 23) ZtoChLep.push_back(vb);
      else if(abs(vb.id()) == 24) WtoLep.push_back(vb);     // it should be empty here
  
  // categorize basic particles
  foreach(const phys::Particle gp, *genParticles){   
    if(abs(gp.id()) == 11 || abs(gp.id()) == 13)                            chleptons.push_back(gp);
    else if(abs(gp.id()) <= 6)                                              quarks.push_back(gp);
    else if(abs(gp.id()) == 12 || abs(gp.id()) == 14 || abs(gp.id()) == 16) neutrinos.push_back(gp);
    else if(gp.id() == 22)                                                  photons.push_back(gp);
  }

  // Clean the chlepton collections from leptons already used to build up the Z and W that decays leptonically
  foreach(phys::Particle &lep, chleptons){
    // Add FSR photons to the lepton 4-momentum
    foreach(const phys::Particle &pho, photons) if(physmath::deltaR(lep,pho) < 0.1) lep.setP4(lep.p4()+pho.p4());

    // Check if the lepton was already used to build up the VB in the genVB collection
    bool matching=false;
    foreach(const phys::Boson<phys::Particle> &vb, ZtoChLep) if(vb.overlapWithDaughters(lep)) matching=true;
    foreach(const phys::Boson<phys::Particle> &vb, WtoLep  ) if(vb.overlapWithDaughters(lep)) matching=true;
    if(matching) continue;

    if(lep.id() == -11 || lep.id() == -13) posleptons.push_back(lep);
    if(lep.id() ==  11 || lep.id() ==  13) negleptons.push_back(lep);
  }
  

  
  // -------------------- Search for additional Z->l+l- candidates ----------------------------
  std::vector<phys::Boson<phys::Particle> > Zllcand;
  foreach(const phys::Particle &p, posleptons) foreach(const phys::Particle &m, negleptons)
    if(abs(p.id()) == abs(m.id())) Zllcand.push_back(phys::Boson<phys::Particle>(m,p,23));

  // This should be the best choice, because either there are 2 Z bosons already in the events, then at most there will be an additional
  // Z, so any ordering fits, or there are 0 Z boson, therefore the best additional boson should come with the mass sorting.
  std::stable_sort(Zllcand.begin(), Zllcand.end(), phys::MassComparator(phys::ZMASS));
  
  if(!Zllcand.empty() && Zllcand.at(0).mass() >= 60 && Zllcand.at(0).mass() <= 120) ZtoChLep.push_back(Zllcand.at(0));
  // ------------------------------------------------------------------------------------------




  
  // ---------------------- Search for Z->nunu in the event -----------------------------------
  ZtoNeutrinos = getZtoX(neutrinos,ZnnCondition());

  
  // std::vector<phys::Boson<phys::Particle> > Znncand;
   
  // foreach(const phys::Particle &n1, neutrinos) foreach(const phys::Particle &n2, neutrinos)
  //   if(n1.id() + n2.id() == 0 && n1.id() > 0) Znncand.push_back(phys::Boson<phys::Particle>(n1,n2,23));
  // std::stable_sort(Znncand.begin(), Znncand.end(), phys::MassComparator(phys::ZMASS));
  // if(!Znncand.empty()){

  //   phys::Boson<phys::Particle> Z0 = Znncand.at(0);
  //   phys::Boson<phys::Particle> Z1;

  //   std::stable_sort(Znncand.begin(), Znncand.end(), phys::ScalarSumPtComparator());
    
  //   foreach(const phys::Boson<phys::Particle> &z, Znncand){
  //     if(!Z0.overlapWithDaughters(z.daughter(0)) && !Z0.overlapWithDaughters(z.daughter(1))){
  // 	Z1 = z;
  // 	break;
  //     }
  //   }  
  //   if(Z0.mass() >= 60 && Z0.mass() <= 120) ZtoNeutrinos.push_back(Z0);
  //   if(Z1.mass() >= 60 && Z1.mass() <= 120) ZtoNeutrinos.push_back(Z1);
  // }

  // foreach(const phys::Boson<phys::Particle> &zcan, Znncand){
  //   bool match = false;
  //   foreach(const phys::Boson<phys::Particle> &znn, ZtoNeutrinos)
  //     if(zcan.overlapWithDaughters(znn.daughter(0)) || zcan.overlapWithDaughters(znn.daughter(1))) match = true;
  //   if(!match && zcan.mass() >= 60 && zcan.mass() <= 120) ZtoNeutrinos.push_back(zcan);
  // }
  // ---------------------------------------------------------------------------------------------


  
  // ----------------------------- Search for W->ln in the event ---------------------------------------
  // at most 1 is expected...
  std::vector<phys::Boson<phys::Particle> > Wlncand;

  // Clean the neutrinos collections from neutrinos already used to build up the Z
  foreach(phys::Particle &n, neutrinos){
    bool matching=false;
    foreach(const phys::Boson<phys::Particle> &vb, ZtoNeutrinos)
      if(vb.overlapWithDaughters(n)) matching=true;
    if(matching) continue;

    foreach(const phys::Particle &l, posleptons)
      if(abs(l.id() + n.id()) == 1 && abs(n.id()) > abs(l.id())) Wlncand.push_back(phys::Boson<phys::Particle>(l,n,copysign(24,n.id())));
    
    foreach(const phys::Particle &l, negleptons)
      if(abs(l.id() + n.id()) == 1 && abs(n.id()) > abs(l.id())) Wlncand.push_back(phys::Boson<phys::Particle>(l,n,copysign(24,n.id())));
  }

  std::stable_sort(Wlncand.begin(), Wlncand.end(), phys::MassComparator(phys::WMASS));
  if(!Wlncand.empty()){
    phys::Boson<phys::Particle> W = Wlncand.at(0);

    if(W.mass() >= 50 && W.mass() <= 110) WtoLep.push_back(W);
  }

  foreach(const phys::Boson<phys::Particle> &wcan, Wlncand){
    bool match = false;
    foreach(const phys::Boson<phys::Particle> &wln, WtoLep)
      if(wcan.overlapWithDaughters(wln.daughter(0)) || wcan.overlapWithDaughters(wln.daughter(1))) match = true;
    if(!match && wcan.mass() >= 50 && wcan.mass() <= 110) WtoLep.push_back(wcan);
  }
  // ---------------------------------------------------------------------------------------------





  // ---------------------- Search for Z->qq in the event -----------------------------------
  std::vector<phys::Boson<phys::Particle> > Zqqcand;
  foreach(const phys::Particle &q1, quarks) foreach(const phys::Particle &q2, quarks)
    if(q1.id() + q2.id() == 0 && q1.id() > 0) Zqqcand.push_back(phys::Boson<phys::Particle>(q1,q2,23));

  std::stable_sort(Zqqcand.begin(), Zqqcand.end(), phys::MassComparator(phys::ZMASS));
  
  if(!Zqqcand.empty()){

    phys::Boson<phys::Particle> Z0 = Zqqcand.at(0);
    phys::Boson<phys::Particle> Z1;

    // attenzione! da controllare. Se il falvour è diverso, usare MassComparator?
    std::stable_sort(Zqqcand.begin(), Zqqcand.end(), phys::ScalarSumPtComparator());
    
    foreach(const phys::Boson<phys::Particle> &z, Zqqcand){
      if(!Z0.overlapWithDaughters(z.daughter(0)) && !Z0.overlapWithDaughters(z.daughter(1))){
	Z1 = z;
	break;
      }
    }  
    if(Z0.mass() >= 60 && Z0.mass() <= 120) ZtoQ.push_back(Z0);
    if(Z1.mass() >= 60 && Z1.mass() <= 120) ZtoQ.push_back(Z1);
  }

  foreach(const phys::Boson<phys::Particle> &zcan, Zqqcand){
    bool match = false;
    foreach(const phys::Boson<phys::Particle> &zqq, ZtoQ)
      if(zcan.overlapWithDaughters(zqq.daughter(0)) || zcan.overlapWithDaughters(zqq.daughter(1))) match = true;
    if(!match && zcan.mass() >= 60 && zcan.mass() <= 120) ZtoQ.push_back(zcan);
  }
  // ---------------------------------------------------------------------------------------------






  // ---------------------- Search for W->qq' in the event -----------------------------------

  // Clean the  quark collection from quarks already used to build up the Z
  std::vector<phys::Particle> nonoverappingq;    
  foreach(const phys::Particle &q, quarks){
    bool matching=false;
    foreach(const phys::Boson<phys::Particle> &vb, ZtoQ)
      if(vb.overlapWithDaughters(q)) matching=true;
    if(matching) continue;
    nonoverappingq.push_back(q);
  }
  std::vector<phys::Boson<phys::Particle> > Wqqcand;    
  foreach(const phys::Particle &q1, nonoverappingq) foreach(const phys::Particle &q2, nonoverappingq)
    if( abs((q1.id() + q2.id()))%2 == 1 && q1.id() > 0) Wqqcand.push_back(phys::Boson<phys::Particle>(q1,q2, copysign(24,-1))); // FIXME
    
  
  std::stable_sort(Wqqcand.begin(), Wqqcand.end(), phys::MassComparator(phys::WMASS));
  
  if(!Wqqcand.empty()){

    phys::Boson<phys::Particle> W = Wqqcand.at(0);
    if(W.mass() >= 50 && W.mass() <= 110) WtoQ.push_back(W);
  }

  foreach(const phys::Boson<phys::Particle> &wcan, Wqqcand){
    bool match = false;
    foreach(const phys::Boson<phys::Particle> &wqq, WtoQ)
      if(wcan.overlapWithDaughters(wqq.daughter(0)) || wcan.overlapWithDaughters(wqq.daughter(1))) match = true;
    if(!match && wcan.mass() >= 50 && wcan.mass() <= 110) WtoQ.push_back(wcan);
  }
  // ---------------------------------------------------------------------------------------------

  
  theHistograms.fill("nZtoChLep", "Number of Z->ll per event", 7,0,7, ZtoChLep.size());
  theHistograms.fill("nZtoNeutrinos", "Number of Z->nn per event", 7,0,7, ZtoNeutrinos.size());
  theHistograms.fill("nWtoLep", "Number of W->lnu per event", 7,0,7, WtoLep.size());
  theHistograms.fill("nZtoQ", "Number of Z->qq per event", 7,0,7, ZtoQ.size());
  theHistograms.fill("nWtoQ", "Number of W->qq' per event", 7,0,7, WtoQ.size());

  int nVBs = ZtoChLep.size() + ZtoNeutrinos.size() + WtoLep.size() + ZtoQ.size() + WtoQ.size();
  theHistograms.fill("nVBs", "Number of VB per event", 7,0,7, nVBs);
  
  //theHistograms.fill("nZnncand", "Number of Z->nn per event", 7,0,7, Znncand.size());
  theHistograms.fill("nneutrinos", "Number of neutrinos per event", 7,0,7, neutrinos.size());
  // // theHistograms.fill("nlepVB", "Number of lepVB per event", 7,0,7, lepVBs.size());
  theHistograms.fill("nquarks", "Number of quarks per event", 7,0,7, quarks.size());
  //   //  
 

  
}















