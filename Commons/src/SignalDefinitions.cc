//-----------FUNCTION: definition of the two ZZ bosons from leptons-------

#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
//#include "VVXAnalysis/Commons/interface/PhysTools.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/DataFormats/interface/GenStatusBit.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include <bitset>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/assign/std/vector.hpp>
using namespace boost::assign;

using namespace std;

typedef std::tuple<uint,uint,int,double> QuarkPairFeatures; //0,1-->indexes of jets;  2-->Id of boson (if known);  3--> mass of boson
typedef std::vector<QuarkPairFeatures> QuarkPairsFeatures;


// Gen Category:

// 0) ZZ->4l
// 1) WZ->3lnu
// 2) ssWW->2l2nu (not set)
// 3) osWW->2l2nu (not set)
// 4)  



//==============================================================//
//                         ZZ analysis                          //
//==============================================================//


///////-------- Getzz: given a vector of Z Boson<Particle> returns a tuple of the 2 Z chosen-------------------------


std::tuple<bool, phys::Boson<phys::Particle>, phys::Boson<phys::Particle> > zz::getZZ(const std::vector<phys::Boson<phys::Particle> >  &ZLL){
  
  
  if(ZLL.size() < 2) return std::make_tuple(false, phys::Boson<phys::Particle>(),phys::Boson<phys::Particle>());
  
  std::vector<phys::Boson<phys::Particle> > Zll;  
  
  foreach(const phys::Boson<phys::Particle>& z,   ZLL) Zll.push_back(z.clone<phys::Particle>());

  // Search for the first Z that must have the mass closest to the nominal Z mass
 
  std::stable_sort(Zll.begin(), Zll.end(), phys::MassComparator(phys::ZMASS));
  phys::Boson<phys::Particle> Z0 = Zll.at(0);
  
  // Search for the second Z that must not be composed by the same daughters as the first Z and is chosen by the largest sum of Z daughters'  as in H->ZZ and ZZ analysies. 
  phys::Boson<phys::Particle> Z1;
  
  std::stable_sort(Zll.begin(), Zll.end(), phys::ScalarSumPtComparator());
  
  foreach(const phys::Boson<phys::Particle> &z, Zll){
    
    
    double dP00 =(Z0.daughter(0).p4() - z.daughter(0).p4()).P();
    double dP01 =(Z0.daughter(0).p4() - z.daughter(1).p4()).P();
    double dP10 =(Z0.daughter(1).p4() - z.daughter(0).p4()).P();
    double dP11 =(Z0.daughter(1).p4() - z.daughter(1).p4()).P();   

    
    if (dP00 > 1e-5 && dP01 > 1e-5 && dP10 > 1e-5 && dP11 > 1e-5){
      Z1 = z;
      break;
    }
  }  
  
  if(Z1.id() == 0) return std::make_tuple(false, Z0, phys::Boson<phys::Particle>());
  
  // Now check that the 4 leptons are not mismatched due to the presence of low mass resonances 
  
  bool passllLowMass = true;


  for(int i = 0; i <=1; ++i) {
  //not very efficent or elegant . To be fixed  

    if(  (( (Z0.daughter(i).p4()+ Z0.daughter((i+1)%2).p4()).M() < 4 ) || ( (Z1.daughter(i).p4()+ Z1.daughter((i+1)%2).p4()).M() < 4)) || 
	 (  (abs(Z0.daughter(0).id())==abs(Z1.daughter(0).id()))  && 
	    ( ( (Z0.daughter(0).id() ==  -1*(Z1.daughter(i).id())  ) && ( (Z0.daughter(0).p4()+ Z1.daughter(i).p4()).M()<4. ) ) ||
	      ( (Z0.daughter(1).id() ==  -1*(Z1.daughter(i).id())   ) && ( (Z0.daughter(1).p4()+ Z1.daughter(i).p4()).M()<4. ) ) ) )
	 ) passllLowMass = false;
            
    }
    

  bool inZMassWindow = true;

  if(Z0.mass() > 120. || Z0.mass() < 60. || Z1.mass() > 120. || Z1.mass() < 60.)
    inZMassWindow = false;


  if(!passllLowMass || !inZMassWindow ) return std::make_tuple(false, Z0, Z1);
  else return std::make_tuple(true, Z0, Z1);
}



zz::SignalTopology zz::getSignalTopology(const std::vector<phys::Particle> &theGenl, std::vector<phys::Particle> &theGenj, std::vector<phys::Particle> &theGenjAK8){
  
  bitset<16>  topology(0);   
  
  std::vector<phys::Particle> theGenlm, theGenlp, genNu;

  foreach(const phys::Particle &p, theGenl){    


    if(abs(p.id()) != 11 && abs(p.id()) != 13 && abs(p.id()) != 12 && abs(p.id()) != 14) continue;

    if(abs(p.id()) == 12 || abs(p.id()) == 14) {genNu.push_back(p); continue;}

    if (p.id() > 0) theGenlm.push_back(p); // negative leptons                                          
    else            theGenlp.push_back(p); // positive leptons 
  }
  
  // Creating and filling the vector of Z boson candidates
  std::vector<phys::Boson<phys::Particle> > Z;
  
  foreach(const phys::Particle &p, theGenlp) foreach(const phys::Particle &m, theGenlm)
    if(abs(p.id()) == abs(m.id())) Z.push_back(phys::Boson<phys::Particle>(m,p,23));
  
  // ---------------------------- Creation of the bosons ----------------------------

  phys::Boson<phys::Particle> Z0, Z1, Z2, Z3, W0, W1, W2;
  

  int numberOfChargedLeptons = theGenlp.size() + theGenlm.size();

  // Not enough lepton to form a WZ or a ZZ boson pair. FIXME! ssWW not contemplated yet!
  if(theGenlp.size() == 0 || theGenlm.size() == 0) return std::make_tuple(topology.to_ulong(), Z0, Z1, Z2, Z3, W0, W1, W2);
    

  // ZZ --> 4l FS
  if(Z.size() >= 2 && numberOfChargedLeptons >=4){

    std::tuple<bool, phys::Boson<phys::Particle>,phys::Boson<phys::Particle> > Zpair = zz::getZZ(Z);
    // Not enough Zs of good quality 
    if(!std::get<0>(Zpair)) return std::make_tuple(topology.to_ulong(), Z0, Z1, Z2, Z3, W0, W1, W2); 
    
    Z0 = std::get<1>(Zpair); 
    Z1 = std::get<2>(Zpair); 
    topology.set(0);
  }

  // WZ --> 3l2nu FS
  else if(Z.size() >= 1 && numberOfChargedLeptons >= 3 && genNu.size() > 0){
    std::tuple<bool, phys::Boson<phys::Particle>,phys::Boson<phys::Particle> > WZpair = wz::getWZ(theGenlm, theGenlp, genNu);
    // Not enough bosons of good quality 
    if(!std::get<0>(WZpair)) return std::make_tuple(topology.to_ulong(), Z0, Z1, Z2, Z3, W0, W1, W2); 

    W0 = std::get<1>(WZpair); 
    Z0 = std::get<2>(WZpair); 
    topology.set(1);
  }

  else return std::make_tuple(topology.to_ulong(), Z0, Z1, Z2, Z3, W0, W1, W2);
  
  // ssWW topology.set(2);
  // osWW topology.set(3);

   

  phys::DiBoson<phys::Particle,phys::Particle> ZZ(Z0,Z1);

  if(vv::inLeptonAcceptance(concatenate(theGenlp,theGenlm),5,2.5,5,2.5))  topology.set(4);   // Fiducial acceptance
  if(vv::inLeptonAcceptance(concatenate(theGenlp,theGenlm),7,2.5,5,2.4))  topology.set(5);   // Detector acceptance
  if(inTriggerPlateau(concatenate(theGenlp,theGenlm)))                    topology.set(6);   // is on trigger plateau    
  if( ZZ.mass() > 100)                                                    topology.set(7);   // Remove additional VB decaying hadronically


  //int Z0DaugID = Z0.daughter(0).id();  
  //int Z1DaugID = Z1.daughter(1).id();
  //if(abs(Z0DaugID) == 13 || abs(Z1DaugID) == 13) topology.set(8);
  //if(abs(Z0DaugID) == 11 || abs(Z1DaugID) == 11) topology.set(9);


  std::vector<phys::Boson<phys::Particle> > lepBosons;
  lepBosons += Z0,Z1,W0,W1;
  std::vector<phys::Boson<phys::Particle> > hadBosons = vv::categorizeHadronicPartOftheEvent(theGenj, theGenjAK8, lepBosons, topology);
  Z3 = hadBosons.at(0);
  W2 = hadBosons.at(1);

  return std::make_tuple(topology.to_ulong(), Z0, Z1, Z2, Z3, W0, W1, W2);
}


// This function cleans the jets from leptons!
std::vector<phys::Boson<phys::Particle> > vv::categorizeHadronicPartOftheEvent(std::vector<phys::Particle> &theGenj,
									       std::vector<phys::Particle> &theGenjAK8,
									       const std::vector<phys::Boson<phys::Particle> >& bosonsToLeptons, 
									       bitset<16>& topology){
  
  phys::Boson<phys::Particle> Z3, W2;

  // Clean the gen jet collection properly
  std::vector<phys::Particle> tmp;
  foreach(const phys::Particle& jet, theGenj){
    bool match = false;
    foreach(const phys::Boson<phys::Particle>& boson, bosonsToLeptons)
      if(boson.daughter(0).id()*boson.daughter(1).id() != 0 && 
	 (physmath::deltaR(boson.daughter(0),jet) < 0.4 ||
	  physmath::deltaR(boson.daughter(1),jet) < 0.4))
	match = true;
    if(!match) tmp.push_back(jet);
  }
  theGenj = tmp;


  // Clean the gen jet collection properly, 0.8 to be checked!
  tmp = std::vector<phys::Particle>();
  foreach(const phys::Particle& jet, theGenjAK8){
    bool match = false;
    foreach(const phys::Boson<phys::Particle>& boson, bosonsToLeptons)
      if(boson.daughter(0).id()*boson.daughter(1).id() != 0 && 
	 (physmath::deltaR(boson.daughter(0),jet) < 0.8 ||
	  physmath::deltaR(boson.daughter(1),jet) < 0.8))
	match = true;
    if(!match) tmp.push_back(jet);
  }
  theGenjAK8 = tmp;

  // FIXME: AK8 numerology still to be added

  
  bool foundWjj(false), foundZjj(false);
 

  // Now check if there are  additional W/Z boson candidates which decays into hadrons
  if(theGenj.size() >= 2) {
    
    // ---- Search for W/Z  ----- 
    QuarkPairsFeatures jetPairsFeatures;

    // make all jet-jet pairs
    for(uint i = 0;  i < theGenj.size()-1; ++i) for(uint j = i+1;  j < theGenj.size(); ++j)
						  jetPairsFeatures.push_back(std::make_tuple(i, j, 0, (theGenj[i].p4() + theGenj[j].p4()).M()));

    // Search for W-->jj
    QuarkPairFeatures bestJetPairW;
    std::stable_sort(jetPairsFeatures.begin(), jetPairsFeatures.end(), phys::MassComparator(0, phys::WMASS));
    bestJetPairW = jetPairsFeatures.front();

    if ( fabs(std::get<3>(bestJetPairW) - phys::WMASS) < 30. ){
      phys::Particle q0 = theGenj[std::get<0>(bestJetPairW)];
      phys::Particle q1 = theGenj[std::get<1>(bestJetPairW)];
      if (q0.pt() < q1.pt()) std::swap(q0,q1); 
      W2 = phys::Boson<phys::Particle>(q0, q1, 24);
      foundWjj = true;
    }
    
    // Search for Z-->jj
    QuarkPairFeatures bestJetPairZ;
    std::stable_sort(jetPairsFeatures.begin(), jetPairsFeatures.end(), phys::MassComparator(0, phys::ZMASS));
    bestJetPairZ = jetPairsFeatures.front();
    
    //
    if ( fabs(std::get<3>(bestJetPairZ) - phys::ZMASS) < 30. ){
      phys::Particle q0 = theGenj[std::get<0>(bestJetPairZ)];
      phys::Particle q1 = theGenj[std::get<1>(bestJetPairZ)];
      if (q0.pt() < q1.pt()) std::swap(q0,q1); 
      Z3 = phys::Boson<phys::Particle>(q0, q1, 23);
      foundZjj = true;
    }
  } 
  
  
  bool foundWj = false;
  bool foundZj = false;
  if(theGenjAK8.size() > 0){
  	// Search for W-->j
		std::stable_sort(theGenjAK8.begin(), theGenjAK8.end(), phys::MassComparator(phys::WMASS));
		foundWj = fabs(theGenjAK8.front().mass() - phys::WMASS) < 30.;
		
		// Search for Z-->j
		std::stable_sort(theGenjAK8.begin(), theGenjAK8.end(), phys::MassComparator(phys::ZMASS));
		foundZj = fabs(theGenjAK8.front().mass() - phys::ZMASS) < 30.;
	}
  
  //------------------------- end searching for the bosons ----------------------------------------

  
  // // Some usefull counters for topology characterization
  // int countJets = 0;
  // int countCentralJets = 0;
  // foreach(const phys::Particle& jet, theGenj)
  //   if(jet.pt() > 30.){
  //     if(abs(jet.eta()) < 4.7) ++countJets;
  //     if(abs(jet.eta()) < 2.4) ++countCentralJets;
  //   }


  //bool hasJets = countJets > 0;
  
  //bool hasAtLeast2jets  =  countJets > 1;
  
  //bool hasCentralJets   = countCentralJets > 0;

  //if(hasJets)          topology.set(4);       //jets (pT>30 GeV and |eta| < 4.7)
  
  //if(hasAtLeast2jets)  topology.set(5);       //2jets
  
  //if(hasCentralJets)   topology.set(6);       //jets (pT>30 GeV and |eta| < 2.4)
  
  if(foundWjj || foundZjj || foundWj || foundZj) topology.set(8); //Hadronic signal
  
  if(foundWjj)         topology.set(9);       //hadronic W
  
  if(foundZjj)         topology.set(10);      //hadronic Z
  
  if(foundWj)          topology.set(11);      //W->j (AK8)
  
  if(foundZj)          topology.set(12);      //Z->j (AK8)
  


  std::vector<phys::Boson<phys::Particle> > hadBosons;
  hadBosons += Z3, W2;
  return hadBosons;

}



bool vv::checkLeptonAcceptance(const phys::Particle &lepton, 
			       const double& pt_e, const double&eta_e, 
			       const double& pt_mu, const double& eta_mu){
  
  
  if     (abs(lepton.id()) == 11 && lepton.pt() > pt_e  && fabs(lepton.eta()) < eta_e ) return true;
  else if(abs(lepton.id()) == 13 && lepton.pt() > pt_mu && fabs(lepton.eta()) < eta_mu) return true;
  else return false;

  //  if(lepton.pt() > 5. && fabs(lepton.eta()) < 2.5) return true;
}

 

bool vv::inLeptonAcceptance(const std::vector<phys::Particle> &leptons,
			    const double& pt_e , const double&eta_e   , 
			    const double& pt_mu, const double& eta_mu ){
  
  bool inAcceptance = true;

  foreach(const phys::Particle& lep, leptons) inAcceptance = inAcceptance && checkLeptonAcceptance(lep, pt_e, eta_e, pt_mu, eta_mu);

  return inAcceptance;
}



bool zz::inTriggerPlateau(const std::vector<phys::Particle> &leptons){
  
  int pt10 = 0; 
  int pt20 = 0;
  
  foreach(const phys::Particle& lep, leptons){
    if(lep.pt() > 10) ++pt10;
    if(lep.pt() > 20) ++pt20;
  }
  if(pt10 > 1 && pt20 > 0) return true;   
  else return false;
}



bool zz::inHiggsFiducialRegion(const zz::SignalTopology &topology){

  if(!std::get<1>(topology).isValid() || !std::get<2>(topology).isValid()) return false;

  phys::Boson<phys::Particle> Z0 = std::get<1>(topology);
  phys::Boson<phys::Particle> Z1 = std::get<2>(topology);
  phys::DiBoson<phys::Particle,phys::Particle> ZZ(Z0,Z1);

  bool massrequirements =  Z0.mass() > 40 && Z0.mass() < 120 && Z1.mass() > 12 && Z1.mass() < 120 &&  ZZ.mass() > 100;
  
  int pt10 = 0; 
  int pt20 = 0;

  bool acceptance = false;

  // Ask for leptons within eta and pt acceptance
  if(std::get<1>(topology).isValid() && std::get<2>(topology).isValid()) 
    acceptance = 
      vv::checkLeptonAcceptance(std::get<1>(topology).daughter(0),7,2.5,5,2.4) && vv::checkLeptonAcceptance(std::get<1>(topology).daughter(1),7,2.5,5,2.4) &&
      vv::checkLeptonAcceptance(std::get<2>(topology).daughter(0),7,2.5,5,2.4) && vv::checkLeptonAcceptance(std::get<2>(topology).daughter(1),7,2.5,5,2.4);
  
  for(int i = 0; i < 2; ++i){
      if(std::get<1>(topology).daughter(i).pt() > 10) ++pt10; 
      if(std::get<1>(topology).daughter(i).pt() > 20) ++pt20;
      if(std::get<2>(topology).daughter(i).pt() > 10) ++pt10; 
      if(std::get<2>(topology).daughter(i).pt() > 20) ++pt20;
     }
  bool trigger = false;
  if(pt10 > 1 && pt20 > 0) trigger = true;
  return massrequirements && acceptance && trigger;
}


bool zz::inTightFiducialRegion(const zz::SignalTopology &topology){

  if(std::get<0>(topology) <= 0) return false;

  int pt10 = 0; 
  int pt20 = 0;

  bool acceptance = false;

  // Ask for leptons within eta and pt acceptance
  if(std::get<1>(topology).isValid() && std::get<2>(topology).isValid()) 
    acceptance = 
      vv::checkLeptonAcceptance(std::get<1>(topology).daughter(0),7,2.5,5,2.4) && vv::checkLeptonAcceptance(std::get<1>(topology).daughter(1),7,2.5,5,2.4) &&
      vv::checkLeptonAcceptance(std::get<2>(topology).daughter(0),7,2.5,5,2.4) && vv::checkLeptonAcceptance(std::get<2>(topology).daughter(1),7,2.5,5,2.4);
  
  for(int i = 0; i < 2; ++i){
      if(std::get<1>(topology).daughter(i).pt() > 10) ++pt10; 
      if(std::get<1>(topology).daughter(i).pt() > 20) ++pt20;
      if(std::get<2>(topology).daughter(i).pt() > 10) ++pt10; 
      if(std::get<2>(topology).daughter(i).pt() > 20) ++pt20;
     }
  bool trigger = false;
  if(pt10 > 1 && pt20 > 0) trigger = true;
  return acceptance && trigger;
}




std::tuple<bool, phys::Boson<phys::Particle>, phys::Boson<phys::Particle> > wz::getWZ(const std::vector<phys::Particle>& lepMinus,
										      const std::vector<phys::Particle>& lepPlus,
										      const std::vector<phys::Particle>& neutrinos){
  // Cases of WZ:
  // H) --> e+  nu_e     , e+e-     -11  12 -11 11 =  1
  // H) --> e-  antinu_e , e+e-      11 -12 -11 11 = -1
  // M) --> mu+ nu_mu    , e+e-     -13  14 -11 11 =  1
  // M) --> mu- antinu_mu, e+e-      13 -14 -11 11 = -1
  // M) --> e+  nu_e     , mu+mu-   -11  12 -13 13 =  1
  // M) --> e-  antie    , mu+mu-    11 -12 -13 13 = -1
  // H) --> mu+ nu_mu    , mu+mu-   -13  14 -13 13 =  1
  // H) --> mu- antinu_mu, mu+mu-    13 -14 -13 13 = -1


  int sumId = 0;
  foreach(const phys::Particle& l, lepMinus)  sumId += l.id();
  foreach(const phys::Particle& l, lepPlus)   sumId += l.id();
  foreach(const phys::Particle& l, neutrinos) sumId += l.id();

  if(abs(sumId) != 1) return std::make_tuple(false, phys::Boson<phys::Particle>(),phys::Boson<phys::Particle>());

  phys::Boson<phys::Particle> W,Z;

  
  // Creation and filling of the vector of Z boson candidates
  std::vector<phys::Boson<phys::Particle> > Zs;
  
  foreach(const phys::Particle &p, lepPlus) foreach(const phys::Particle &m, lepMinus)
    if(abs(p.id()) == abs(m.id())) Zs.push_back(phys::Boson<phys::Particle>(m,p,23));

  // Mixed (M) cases presents only 1 Z candidate
  if(Zs.size() == 1){
    Z = Zs.at(0);
    phys::Particle lepFromW, neuFromW;
    neuFromW = neutrinos.at(0);
    if(neuFromW.id() > 0)
      foreach(const phys::Particle &p, lepPlus){
	if(abs(p.id()) == neuFromW.id()){
	  lepFromW = p;
	  break;
	}
      }
    else
      foreach(const phys::Particle &p, lepMinus){
	if(p.id() == abs(neuFromW.id())){
	  lepFromW = p;
	  break;
	}
      }
    W = phys::Boson<phys::Particle>(lepFromW,neuFromW, copysign(24,neuFromW.id()));
  }
  // Homogeneous (H) cases presents 2 Z candidates
  if(Zs.size() == 2){
    Z = fabs((Zs.at(0).mass() - phys::ZMASS)) < fabs((Zs.at(1).mass() - phys::ZMASS)) ? Zs.at(0) : Zs.at(1);
    phys::Particle lepFromW, neuFromW;
    neuFromW = neutrinos.at(0);
    if(neuFromW.id() > 0)
      foreach(const phys::Particle &p, lepPlus){
	if(abs(p.id()) == neuFromW.id() && 
	   p != Z.daughter(0) && p != Z.daughter(1)){
	  lepFromW = p;
	  break;
	}
      }
    else
      foreach(const phys::Particle &p, lepMinus){
	if(p.id() == abs(neuFromW.id()) &&
	   p != Z.daughter(0) && p != Z.daughter(1)){
	  lepFromW = p;
	  break;
	}
      }
    W = phys::Boson<phys::Particle>(lepFromW,neuFromW, copysign(24,neuFromW.id()));    
  }
  else{return std::make_tuple(false, phys::Boson<phys::Particle>(),phys::Boson<phys::Particle>());}

  // FIXME: Need to add mass requirements here

  return std::make_tuple(true, W, Z);

}










////////////////////////////////////////////////////////////////////////////////
// OLD stuff //
////////////////////////////////////////////////////////////////////////////////


std::pair<phys::Boson<phys::Particle> ,phys::Boson<phys::Particle> > zzw::makeZBosonsFromLeptons(const std::vector<phys::Particle>& lm, const std::vector<phys::Particle>& lp, int leptonCode, float mZ){
    
    phys::Boson<phys::Particle> Z0;
    phys::Boson<phys::Particle> Z1;
    float minMDiff=99999.;
    if (leptonCode == 4) {
      for (int k=0; k<2; ++k) {
	for (int j=0; j<2; ++j) {
	  float mDiff = fabs((lp[k].p4() + lm[j].p4()).M() - mZ);
	  if ( mDiff < minMDiff ) {
	    minMDiff=mDiff;   
	    
	    Z0.setDaughter(0,lp[k]);
	    Z0.setDaughter(1,lm[j]);

	    Z1.setDaughter(0,lp[(k+1)%2]);
	    Z1.setDaughter(1,lm[(j+1)%2]);
	  }      
	} 	
      }
    }
    else { 
      for (int z=0; z<2; ++z) {
	if ( fabs(lp[z].id()) == fabs(lm[0].id()) ) { 
	  
	  Z0.setDaughter(0,lp[z]);
	  Z0.setDaughter(1,lm[0]);

	  Z1.setDaughter(0,lp[(z+1)%2]);
	  Z1.setDaughter(1,lm[1]);      
	}
      }	
    }

    Z0.setId(23);
    Z1.setId(23);
    
    return std::make_pair(Z0,Z1);
  }





zzw::GenTopology zzw::getGenTopology(int signalDefinition, 
				     const std::vector<phys::Particle> &theGenl, const std::vector<phys::Particle> &theGenj, 
				     const std::vector<phys::Particle> &theGenZ, const std::vector<phys::Particle> &theGenW){
  int categoryNum = 999; 
  
  int numMu       = 0, numE     = 0;
  std::vector<phys::Particle> theGenlm, theGenlp, theGenq;

  foreach(const phys::Particle &p, theGenj) if (abs(p.id()) < 7) theGenq.push_back(p);  // quarks
  
  foreach(const phys::Particle &p, theGenl){
    numE  = abs(p.id()) == 11 ? numE+1  : numE; 
    numMu = abs(p.id()) == 13 ? numMu+1 : numMu;
   
    if (p.id() > 0)                   theGenlm.push_back(p); // negative leptons                                          
    else                              theGenlp.push_back(p); // positive leptons 
  }

  int leptonCode = 0;
  if(  numMu == 2 && numE == 2)                                leptonCode = 2;
  if( (numMu == 4 && numE == 0) || (numMu == 0 && numE == 4) ) leptonCode = 4;

  phys::Boson<phys::Particle> Z0, Z1, Z2, W;

  if ( leptonCode == 2 || leptonCode == 4) {

    bool isWloose  = false, isZloose  = false, isWtight  = false, isZtight  = false;

    phys::Particle q0, q1;
    int bosonId = -99;

    if(theGenq.size() >= 2) {
      QuarkPairsFeatures quarkPairsFeatures;
    
      for(uint i = 0;  i < theGenq.size()-1; ++i) for(uint j = i+1;  j < theGenq.size(); ++j)
	quarkPairsFeatures.push_back(std::make_tuple(i, j, vvx::makeVBosonsFromIds(theGenq[i].id(), theGenq[j].id()), (theGenq[i].p4() + theGenq[j].p4()).M()));
      
      // ----- Search for a true W in the event -----
      std::stable_sort(quarkPairsFeatures.begin(), quarkPairsFeatures.end(), phys::MassComparator(24, phys::WMASS));
      QuarkPairFeatures bestQuarkPair = quarkPairsFeatures.front();
      if(abs(std::get<2>(bestQuarkPair)) == 24 and fabs(std::get<3>(bestQuarkPair) - phys::WMASS) < 10){
	q0 = theGenq[std::get<0>(bestQuarkPair)];
	q1 = theGenq[std::get<1>(bestQuarkPair)];
	if ( q0.pt() < q1.pt() ) { q0 = theGenq[std::get<1>(bestQuarkPair)];  q1 = theGenq[std::get<0>(bestQuarkPair)]; }
	bosonId =  std::get<2>(bestQuarkPair);
	isWtight = true;
      }
      // --------------------------------------------
      
      // ----- Search for a true hadronic Z in the event -----
      if(!isWtight){
	std::stable_sort(quarkPairsFeatures.begin(), quarkPairsFeatures.end(), phys::MassComparator(23, phys::ZMASS));
	QuarkPairFeatures bestQuarkPair = quarkPairsFeatures.front();
	if(abs(std::get<2>(bestQuarkPair)) == 23 and fabs(std::get<3>(bestQuarkPair) - phys::ZMASS) < 10){
	  q0 = theGenq[std::get<0>(bestQuarkPair)];
	  q1 = theGenq[std::get<1>(bestQuarkPair)];
	  if ( q0.pt() < q1.pt() ) { q0 = theGenq[std::get<1>(bestQuarkPair)];  q1 = theGenq[std::get<0>(bestQuarkPair)]; }
	  bosonId =  std::get<2>(bestQuarkPair);
	  isZtight = true;
	}
      }
    // --------------------------------------------
    }


    if(theGenj.size() >= 2) {

      // ---- Search for loose W/Z (i.e., not true boson, but rather combinations of partons that resemble a boson ----- 
      if(!isWtight && !isZtight){
	QuarkPairsFeatures jetPairsFeatures;
	for(uint i = 0;  i < theGenj.size()-1; ++i) for(uint j = i+1;  j < theGenj.size(); ++j)
	  jetPairsFeatures.push_back(std::make_tuple(i, j, 0, (theGenj[i].p4() + theGenj[j].p4()).M()));

	QuarkPairFeatures bestJetPairW;
	std::stable_sort(jetPairsFeatures.begin(), jetPairsFeatures.end(), phys::MassComparator(0, phys::WMASS));
	bestJetPairW = jetPairsFeatures.front();
      
	QuarkPairFeatures bestJetPairZ;
	std::stable_sort(jetPairsFeatures.begin(), jetPairsFeatures.end(), phys::MassComparator(0, phys::ZMASS));
	bestJetPairZ = jetPairsFeatures.front();

	if ( fabs(std::get<3>(bestJetPairZ) - phys::ZMASS) < 10. ){
	  isZloose = true; 
	  bosonId = 123;
	  q0 = theGenj[std::get<0>(bestJetPairZ)];
	  q1 = theGenj[std::get<1>(bestJetPairZ)];
	  if ( q0.pt() < q1.pt() ) { q0 = theGenj[std::get<1>(bestJetPairZ)];  q1 = theGenj[std::get<0>(bestJetPairZ)]; }
	}

	// Give priority to W loose, accordingly with background categorization
	if ( fabs(std::get<3>(bestJetPairW) - phys::WMASS) < 10. ){
	  isWloose = true;
	  bosonId = 124;  
	  q0 = theGenj[std::get<0>(bestJetPairW)];
	  q1 = theGenj[std::get<1>(bestJetPairW)];
	  if ( q0.pt() < q1.pt() ) { q0 = theGenj[std::get<1>(bestJetPairW)];  q1 = theGenj[std::get<0>(bestJetPairW)]; }
	}
      }
      // --------------------------------------------
    }
    
    //--------------------1: MC history------------------------------------
    if ( signalDefinition==1 ) {              
      
      bool LeptonsMotherSelec = true;   
      for(int t=0; t<4; ++t) LeptonsMotherSelec = LeptonsMotherSelec && theGenl[t].motherId() == 23;
      
      
      if ( theGenW.size() == 1) isWtight = true;      //definition of tight W (mass + cat)
      if ( theGenZ.size() == 3) isZtight = true;      //definition of tight Z (mass + cat)
      
      if ( theGenZ.size() >= 2 && LeptonsMotherSelec ) {

	Z0.setDaughter(0, theGenl[0]);
	Z0.setDaughter(1, theGenl[1]);
	Z0.setId(theGenZ[0].id());

 	Z1.setDaughter(0, theGenl[2]);
	Z1.setDaughter(1, theGenl[3]);
	Z1.setId(theGenZ[1].id());


	if ( isWtight ) {       
	  // FIXME... it should be done using parentage
	  W.setDaughter(0, q0);
	  W.setDaughter(1, q1);
	  W.setId(theGenW[0].id());  
	}
	 
	else if ( isZtight ) {
	  // FIXME... it should be done using parentage
	  Z2.setDaughter(0, q0);
	  Z2.setDaughter(1, q1);
	  Z2.setId(theGenZ[2].id());
	}
      }
    }
    

    // -----------------2: Real signal, MadGraph pairing------------------
    else if ( signalDefinition==2 ) {         
      
      Z0.setDaughter(0, theGenl[0]);
      Z0.setDaughter(1, theGenl[1]);
      Z0.setId(theGenZ[0].id());
      
      Z1.setDaughter(0, theGenl[2]);
      Z1.setDaughter(1, theGenl[3]);
      Z1.setId(theGenZ[1].id());
      
      if (isWloose || isWtight) {      //definition of tight W (mass + cat)
	
    	W.setDaughter(0, q0);
	W.setDaughter(1, q1);
	W.setId(bosonId);
	
      } else if (isZloose || isZtight) {     //definition of tight Z (mass + cat)
	
	Z2.setDaughter(0, q0);
	Z2.setDaughter(1, q1);
	Z2.setId(bosonId);
      } 	
    }
 
   
    //-----------------3: Real signal, real pairing-----------------------
    else if ( signalDefinition==3 ) {         
      std::pair<phys::Boson<phys::Particle>, phys::Boson<phys::Particle> > ZZ = makeZBosonsFromLeptons(theGenlm, theGenlp, leptonCode, phys::ZMASS);

      Z0 = ZZ.first;
      Z1 = ZZ.second;
      
      if (isWloose || isWtight) {    //definition of tight W (mass + cat)
	
	W.setDaughter(0, q0);
	W.setDaughter(1, q1);
	W.setId(bosonId);
	
      } else if (isZloose || isZtight) {   //definition of tight Z (mass + cat)
	
	Z2.setDaughter(0, q0);
	Z2.setDaughter(1, q1);
	Z2.setId(bosonId);
      }     
    } 
    
    else { cout << "*** Signal definition not found! ***" << endl; abort(); }

    //=====================================================================================

    bool hasZZ4l    = fabs(Z0.p4().M()-phys::ZMASS) < 10. && fabs(Z1.p4().M()-phys::ZMASS) < 10.;    
    bool isMySignal = hasZZ4l && isWtight;
    bool has3Z      = hasZZ4l && isZtight;
      
    bool passEtaAccLep = true;
      
    //for(int i=0; i<4; ++i) passEtaAccLep = passEtaAccLep && fabs(theGenl[i].eta()) < 2.5;
     
    // eta cut for all leptons
    
    if ( passEtaAccLep) {

      // ========== Signal: ZZW ==========

      if ( isMySignal ){
	categoryNum = 0;
	//bool passZZacc = fabs(Z0.daughter(0).eta()) < 2.5 && fabs(Z0.daughter(1).eta()) < 2.5 &&
	//  fabs(Z1.daughter(0).eta()) < 2.5 && fabs(Z1.daughter(1).eta()) < 2.5;
	
	//	bool passWacc = fabs(W.daughter(0).eta()) < 2.5 && fabs(W.daughter(1).eta()) < 2.5 && 
	//  (fabs(W.daughter(0).pt()) > 20 || fabs(W.daughter(1).pt()) > 20); 
	
	//if(passZZacc && passWacc) categoryNum = 0;
	//if(passZZacc) categoryNum = 101;
	//if(passWacc) categoryNum = 102;
      }

      // ========== Background ==========

      else {     
	
	if( hasZZ4l ){

	  // ---------- ZZZ ---------- 

	  if ( has3Z )                                   categoryNum = 1;
	  
	  // ---------- ZZWloose ----------

	  else if ( !has3Z && isWloose )                 categoryNum = 2;
	  
	  // ---------- ZZZloose ----------

	  else if ( !has3Z && !isWloose && isZloose )    categoryNum = 3;
	  
	  // ---------- ZZ+X ----------

	  else                                           categoryNum = 4;
	}

	else {
	  
	  // ---------- WZ+X ----------

	  if ( isWtight )                                categoryNum = 5;
	  	  
	  // ---------- ZZjj+X ----------
	  
	  else if ( isZtight )                           categoryNum = 6;
	  
	  // ---------- ZWloose+X ----------
	  
	  else if ( !isWtight && !isZtight && isWloose ) categoryNum = 7;
	  	  
	  // ---------- ZZjj+X ----------

	  else if ( !isZtight && !isWloose && isZloose ) categoryNum = 8;
	  
	  // ---------- Z+X+Y----------

	  else              	                         categoryNum = 9;
  	}		
      }      
    }
  }
 

 //  if(categoryNum != 0){
//     cout << "----------------------------------------------------------------------------------" << endl;
//     cout << "Category: " << categoryNum << endl;
//     // cout<<"Run: " << event.run() << " event: " << event.id().event() <<endl;
//     cout << "----------------------------------------------------------------------------------" << endl;

//     foreach(const phys::Particle& p, theGenl) cout<< p.id() << " " << p.eta() << endl;
//     foreach(const phys::Particle& p, theGenj) cout<< p.id() << endl;

//     cout<<"Z0: " << Z0.id() << " " << Z0.p4().M() << " daughters id: "  << Z0.daughter(0).id() << " " << Z0.daughter(1).id() << endl;
//     cout<<"Z1: " << Z1.id() << " " << Z1.p4().M() << " daughters id: "  << Z1.daughter(0).id() << " " << Z1.daughter(1).id() << endl;
//     cout<<"Z2: " << Z2.id() << " " << Z2.p4().M() << " daughters id: "  << Z2.daughter(0).id() << " " << Z2.daughter(1).id() << endl;
//     cout<<"W: "  << W.id()  << " " << W.p4().M()  << " daughters id: "  << W.daughter(0).id()  << " " << W.daughter(1).id()  << endl;
    
//     foreach(const phys::Particle& p, theGenZ) cout<<"Z true: " << p.id() << " " << p.p4().M() << endl;
//     foreach(const phys::Particle& p, theGenW) cout<<"W true: " << p.id() << " " << p.p4().M() << endl;

//   }

  return std::make_tuple(categoryNum, Z0, Z1, Z2, W);
}





int vvx::makeVBosonsFromIds(int j0Id, int j1Id) {    
  
  if ( abs(j0Id) < 6 && abs(j1Id ) < 6) {
    if( (j0Id*j1Id) <0 && (abs(j0Id + j1Id) == 1 || abs(j0Id + j1Id) == 3) ) {
      if( j0Id % 2 == 0 )       return copysign(24,j0Id);  // W
      else if( j1Id % 2 == 0 )  return copysign(24,j1Id);  // W
      else return 0;
    }
    else if( j0Id + j1Id == 0 ) return 23;                 // Z
    else return 0;                             
    
  }
  else return 0;
}
