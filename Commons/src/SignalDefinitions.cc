//-----------FUNCTION: definition of the two ZZ bosons from leptons-------

#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/PhysTools.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using std::cout;
using std::endl;

const float MZ = 91.19;
const float MW = 80.39;

typedef std::tuple<uint,uint,int,double> QuarkPairFeatures;
typedef std::vector<QuarkPairFeatures>   QuarkPairsFeatures;

struct MassComparator{
  MassComparator(int id, const double& ref): id_(id), ref_(ref){}
  
  bool operator()(const QuarkPairFeatures & a , 
		  const QuarkPairFeatures & b) const{ 
    if       (abs(std::get<2>(a)) == id_ && abs(std::get<2>(b)) != id_) return true;
    else if  (abs(std::get<2>(a)) != id_ && abs(std::get<2>(b)) == id_) return false;
    else if  (abs(std::get<2>(a)) == id_ && abs(std::get<2>(b)) == id_)
      return fabs(std::get<3>(a) - ref_) < fabs(std::get<3>(b) - ref_); 
    else    
      return true;
  }
  
  int id_;
  double ref_;
};


std::pair<phys::Boson<phys::Particle> ,phys::Boson<phys::Particle> > zzw::makeZBosonsFromLeptons(const std::vector<const reco::Candidate *>& lm, const std::vector<const reco::Candidate *>& lp, int leptonCode, float mZ){
  std::vector<phys::Particle> plm;
  std::vector<phys::Particle> plp;

  foreach(const reco::Candidate *gp, lm)
    plm.push_back(phys::convert(*gp));
  
  foreach(const reco::Candidate *gp, lp)
    plp.push_back(phys::convert(*gp));
   
  return zzw::makeZBosonsFromLeptons(plm, plp, leptonCode, mZ);
}
  
  
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



int zzw::makeVBosonsFromIds(int j0Id, int j1Id) {    
  
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
	quarkPairsFeatures.push_back(std::make_tuple(i, j, makeVBosonsFromIds(theGenq[i].id(), theGenq[j].id()), (theGenq[i].p4() + theGenq[j].p4()).M()));
      
      // ----- Search for a true W in the event -----
      std::stable_sort(quarkPairsFeatures.begin(), quarkPairsFeatures.end(), MassComparator(24, MW));
      QuarkPairFeatures bestQuarkPair = quarkPairsFeatures.front();
      if(abs(std::get<2>(bestQuarkPair)) == 24 and fabs(std::get<3>(bestQuarkPair) - MW) < 10){
	q0 = theGenq[std::get<0>(bestQuarkPair)];
	q1 = theGenq[std::get<1>(bestQuarkPair)];
	if ( q0.pt() < q1.pt() ) { q0 = theGenq[std::get<1>(bestQuarkPair)];  q1 = theGenq[std::get<0>(bestQuarkPair)]; }
	bosonId =  std::get<2>(bestQuarkPair);
	isWtight = true;
      }
      // --------------------------------------------
      
      // ----- Search for a true Z in the event -----
      if(!isWtight){
	std::stable_sort(quarkPairsFeatures.begin(), quarkPairsFeatures.end(), MassComparator(23, MZ));
	QuarkPairFeatures bestQuarkPair = quarkPairsFeatures.front();
	if(abs(std::get<2>(bestQuarkPair)) == 23 and fabs(std::get<3>(bestQuarkPair) - MZ) < 10){
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
	std::stable_sort(jetPairsFeatures.begin(), jetPairsFeatures.end(), MassComparator(0, MW));
	bestJetPairW = jetPairsFeatures.front();
      
	QuarkPairFeatures bestJetPairZ;
	std::stable_sort(jetPairsFeatures.begin(), jetPairsFeatures.end(), MassComparator(0, MZ));
	bestJetPairZ = jetPairsFeatures.front();

	if ( fabs(std::get<3>(bestJetPairZ) - MZ) < 10. ){
	  isZloose = true; 
	  bosonId = 123;
	  q0 = theGenj[std::get<0>(bestJetPairZ)];
	  q1 = theGenj[std::get<1>(bestJetPairZ)];
	  if ( q0.pt() < q1.pt() ) { q0 = theGenj[std::get<1>(bestJetPairZ)];  q1 = theGenj[std::get<0>(bestJetPairZ)]; }
	}

	// Give priority to W loose, accordingly with background categorization
	if ( fabs(std::get<3>(bestJetPairW) - MW) < 10. ){
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
      std::pair<phys::Boson<phys::Particle>, phys::Boson<phys::Particle> > ZZ = makeZBosonsFromLeptons(theGenlm, theGenlp, leptonCode, MZ);

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

    bool hasZZ4l    = fabs(Z0.p4().M()-MZ) < 10. && fabs(Z1.p4().M()-MZ) < 10.;    
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
