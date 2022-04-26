#include "VVXAnalysis/TreeAnalysis/interface/PPZZAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/assign/std/vector.hpp>
using namespace boost::assign;

using std::cout;
using std::endl;


using namespace phys;

Int_t PPZZAnalyzer::cut() {

  return 1;
}



void PPZZAnalyzer::analyze(){

  //cout << "------------------------------------------------------------------"<<endl;
  //cout << "Run: " << run << " event: " << event << endl;

  std::vector<phys::Proton> eventprotons;
  std::vector<phys::Boson<phys::Lepton>> eventZ;

  foreach(const phys::Proton sproton, *singleRPprotons){
    if(sproton.xi()>0.04){
      eventprotons.push_back(sproton);}}
  foreach(const phys::Proton mproton, *multiRPprotons){
    if(mproton.xi()>0.04){
      eventprotons.push_back(mproton);}}

  phys::Proton p1, p2;
  phys::ProtonPair pp;
  phys::Lepton muon1, muon2;
  phys::Electron electron1, electron2;


  //constructing Z bosons
  
  if(muons->size()>=2){
    foreach(const phys::Lepton muon, *muons){
      if(muon.id()==13){
	muon1=muon;}
      foreach(const phys::Lepton muonp, *muons){
	if(muon.id()==-13){
	  muon2=muonp;}
	phys::Boson<phys::Lepton> Z1 = phys::Boson<phys::Lepton>(muon1,muon2);
	if(Z1.mass()>70 && Z1.mass()<110){
	  theHistograms->fill("mZmu","mass of mumu pair",25,50,150,Z1.mass(),theWeight);
	  eventZ.push_back(Z1);
	}
      }
    }
  }

  if(electrons->size()>=2){
    foreach(const phys::Lepton electron, *electrons){
      if(electron.id()==11){
	electron1=electron;}
      foreach(const phys::Lepton electronp, *electrons){
	if(electron.id()==-11){
	  electron2=electronp;}
	phys::Boson<phys::Lepton> Z2 = phys::Boson<phys::Lepton>(electron1,electron2);
	if(Z2.mass()>70 && Z2.mass()<110){
	  theHistograms->fill("mZee","mass of ee pair",25,50,150,Z2.mass(),theWeight);
	  eventZ.push_back(Z2);
	}
      }
    }
  }

  double bestmass=9999.;
  double maxpt=0;
  double secondmaxpt=0;
  int bestmassindex=0;
  int maxptindex=0;
  int secondmaxptindex=0;
  
  for(int i=0; i<eventZ.size();i++){
    if(eventZ[i].mass()-phys::ZMASS<bestmass){
      bestmassindex=i;
      bestmass=eventZ[i].mass()-phys::ZMASS;}
    if((eventZ[i].daughter(0).pt()+eventZ[i].daughter(1).pt())>maxpt){
      secondmaxpt=maxpt;
      secondmaxptindex=maxptindex;
      maxpt=eventZ[i].daughter(0).pt()+eventZ[i].daughter(1).pt();
      maxptindex=i;
    }   
  }

  TLorentzVector p4;
  if(eventZ.size()>=2){
    phys::Boson<phys::Lepton> Z1 = eventZ[bestmassindex];
    phys::Boson<phys::Lepton> Z2;
  
    if(bestmassindex!=maxptindex){
      Z2=eventZ[maxptindex];}
    else {Z2=eventZ[secondmaxptindex];}
  
    p4=Z1.p4()+Z2.p4();
    theHistograms->fill("mZZ","mass of ZZ pair",50,0,2000,p4.M(),theWeight);
    theHistograms->fill("yZZ","Expected rapidity of ZZ pair",25,0,1,p4.Rapidity(),theWeight);
  }
  

  TH2F *th2 = new TH2F("th2","2D matching map",100,-5,5,100,-1,1);

  for(int i=0; i<eventprotons.size();i++){
    p1=eventprotons[i];
    for(int j=i+1; j<eventprotons.size();j++){
      p2=eventprotons[j];
      pp=phys::ProtonPair(p1,p2);
      theHistograms->fill("mpp","Expected mass of ZZ pair",25,0,650,pp.mpp(),theWeight);
      theHistograms->fill("ypp","Expected rapidity of ZZ pair",25,0,1,pp.ypp(),theWeight);

      if(eventZ.size()>=2){
	theHistograms->fill("massdiff","Difference between pp and ZZ mass",50,-5,5,1-p4.M()/pp.mpp(),theWeight);
	//theHistograms->fill("massdiff2","Difference between pp and ZZ mass",25,-50,50,p4.M()-pp.mpp(),theWeight);
	theHistograms->fill("ydiff","Difference between pp and ZZ rapidity",25,-1,1,pp.ypp()-p4.Rapidity(),theWeight);

	th2->Fill(1-p4.M()/pp.mpp(),pp.ypp()-p4.Rapidity());
      }
    }
  }

  //theHistograms->fill("mZZ2","mass of ZZ pair",50,0,2000,ZZ->mass(),theWeight);
  //TFile *file1= new TFile("th2.root","UPDATE");
  //th2->Write();
  //file1->Write("th2");
  //file1->Close();

  eventprotons.clear();
  eventZ.clear();
}
