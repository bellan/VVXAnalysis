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

template<class T1, class T2>
std::vector<phys::Boson<phys::Particle> > VVXAnalyzer::getVtoX(const std::vector<phys::Particle> & collectionX1,
							       const std::vector<phys::Particle> & collectionX2,
							       T1 idcondition, T2 masswindow, const double& referenceMass){

  std::vector<phys::Boson<phys::Particle> > VtoX;
  std::vector<phys::Boson<phys::Particle> > VtoXcand;
   
  foreach(const phys::Particle &p1, collectionX1) foreach(const phys::Particle &p2, collectionX2)
    if(idcondition(p1,p2)) VtoXcand.push_back(phys::Boson<phys::Particle>(p1,p2,
									  p1.id()+p2.id()==0 ? 23 : copysign(24,p1.charge()+p2.charge()) ));
  std::stable_sort(VtoXcand.begin(), VtoXcand.end(), phys::MassComparator(referenceMass));
  if(!VtoXcand.empty()){

    phys::Boson<phys::Particle> V0 = VtoXcand.at(0);
    phys::Boson<phys::Particle> V1;

    std::stable_sort(VtoXcand.begin(), VtoXcand.end(), phys::ScalarSumPtComparator());
    
    foreach(const phys::Boson<phys::Particle> &v, VtoXcand){
      if(!V0.overlapWithDaughters(v.daughter(0)) && !V0.overlapWithDaughters(v.daughter(1))){
	V1 = v;
	break;
      }
    }  
    if(masswindow(V0)) VtoX.push_back(V0);
    if(masswindow(V1)) VtoX.push_back(V1);
  }

  foreach(const phys::Boson<phys::Particle> &vcan, VtoXcand){
    bool match = false;
    foreach(const phys::Boson<phys::Particle> &vxx, VtoX)
      if(vcan.overlapWithDaughters(vxx.daughter(0)) || vcan.overlapWithDaughters(vxx.daughter(1))) match = true;
    if(!match && masswindow(vcan)) VtoX.push_back(vcan);
  }
  return VtoX;
}

std::vector<phys::Particle> VVXAnalyzer::removeOverlaps(const std::vector<phys::Particle> &collectionX, const std::vector<phys::Boson<phys::Particle> >& collectionVB){
  std::vector<phys::Particle> nonoverapping;    
  foreach(const phys::Particle &p, collectionX){
    bool matching=false;
    foreach(const phys::Boson<phys::Particle> &vb, collectionVB)
      if(vb.overlapWithDaughters(p)) matching=true;
    if(matching) continue;
    nonoverapping.push_back(p);
  }
  return nonoverapping;
}




void VVXAnalyzer::analyze(){

  // Warning. Inside genVBParticles there are ONLY bosons that satisfy the following criterions 
  // 1- build a ZZ and pass certain mass cuts (fully leptonic decay only!). 
  // 2- build a WZ and pass certain mass cuts (fully leptonic decay only!).
  // 3- an additional hadronically decaying VB is reconstructed using genJets
  // Taus are not stable particles therefore are inside a separate collection!
  
  // Prepare the containers
  std::vector<phys::Particle> quarks, antiquarks;
  std::vector<phys::Particle> chleptons, posleptons, negleptons, photons; // chleptons is an UNCLEANED collection
  std::vector<phys::Particle> neutrinos, antineutrinos;
  
  std::vector<phys::Boson<phys::Particle> > ZtoChLep;
    

  // Take Z from the event topology (there should not be W bosons...)
  foreach(const phys::Boson<phys::Particle> &vb, *genVBParticles)
    if(vb.decayType() == 1 && vb.id()  == 23) ZtoChLep.push_back(vb);
  
  // categorize basic particles
  foreach(const phys::Particle gp, *genParticles){   
    if(abs(gp.id()) == 11 || abs(gp.id()) == 13)                            chleptons.push_back(gp);
    else if(abs(gp.id()) <= 6)
      if(gp.id() > 0)                                                       quarks.push_back(gp);
      else                                                                  antiquarks.push_back(gp);
    else if(abs(gp.id()) == 12 || abs(gp.id()) == 14 || abs(gp.id()) == 16)
      if(gp.id() > 0)                                                       neutrinos.push_back(gp);
      else                                                                  antineutrinos.push_back(gp);
    else if(gp.id() == 22)                                                  photons.push_back(gp);
  }

  // Clean the chlepton collections from leptons already used to build up the Z and W that decays leptonically
  foreach(phys::Particle &lep, chleptons){
    // Add FSR photons to the lepton 4-momentum
    foreach(const phys::Particle &pho, photons) if(physmath::deltaR(lep,pho) < 0.1) lep.setP4(lep.p4()+pho.p4());

    // Check if the lepton was already used to build up the VB in the genVB collection
    bool matching=false;
    foreach(const phys::Boson<phys::Particle> &vb, ZtoChLep) if(vb.overlapWithDaughters(lep)) matching=true;
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
  
  if(!Zllcand.empty() && ZMassWindow()(Zllcand.at(0))) ZtoChLep.push_back(Zllcand.at(0));
  // ------------------------------------------------------------------------------------------




  
  // ---------------------- Search for Z->nunu in the event -----------------------------------------
  std::vector<phys::Boson<phys::Particle> > ZtoNeutrinos = getVtoX(neutrinos, antineutrinos, 
								   ZDaughtersIdCondition(), ZMassWindow(), phys::ZMASS);
  // ------------------------------------------------------------------------------------------------
  
  // ---------------------- Search for Z->qq in the event -------------------------------------------
  std::vector<phys::Boson<phys::Particle> > ZtoQ          = getVtoX(quarks, antiquarks,
								    ZDaughtersIdCondition(), ZMassWindow(), phys::ZMASS);
  // ------------------------------------------------------------------------------------------------

  // ---------------------- Search for W->qq' in the event -----------------------------------
  std::vector<phys::Boson<phys::Particle> > WtoQ          = getVtoX(removeOverlaps(quarks,ZtoQ), removeOverlaps(antiquarks,ZtoQ),
								    WqqDaughtersIdCondition(), WMassWindow(), phys::WMASS);
  // ------------------------------------------------------------------------------------------------

  // ----------------------------- Search for W->ln in the event ---------------------------------------
  // at most 1 is expected...
  std::vector<phys::Boson<phys::Particle> > WtoLep        = getVtoX(posleptons, removeOverlaps(neutrinos,ZtoNeutrinos),
								    WlnDaughtersIdCondition(), WMassWindow(), phys::WMASS);
  std::vector<phys::Boson<phys::Particle> > WmtoLep       = getVtoX(negleptons, removeOverlaps(antineutrinos,ZtoNeutrinos),
								    WlnDaughtersIdCondition(), WMassWindow(), phys::WMASS);
  WtoLep.insert(WtoLep.end(), WmtoLep.begin(), WmtoLep.end());



  
  theHistograms.fill("nZtoChLep", "Number of Z->ll per event", 7,0,7, ZtoChLep.size());
  theHistograms.fill("nZtoNeutrinos", "Number of Z->nn per event", 7,0,7, ZtoNeutrinos.size());
  theHistograms.fill("nWtoLep", "Number of W->lnu per event", 7,0,7, WtoLep.size());
  theHistograms.fill("nZtoQ", "Number of Z->qq per event", 7,0,7, ZtoQ.size());
  theHistograms.fill("nWtoQ", "Number of W->qq' per event", 7,0,7, WtoQ.size());

  int nVBs = ZtoChLep.size() + ZtoNeutrinos.size() + WtoLep.size() + ZtoQ.size() + WtoQ.size();
  theHistograms.fill("nVBs", "Number of VB per event", 7,0,7, nVBs);
  
  theHistograms.fill("nneutrinos", "Number of neutrinos per event", 7,0,7, neutrinos.size());
  theHistograms.fill("nquarks", "Number of quarks per event", 7,0,7, quarks.size());
  //   //  
 

  
}















