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

  foreach(const phys::Proton sproton, *singleRPprotons){
    if(sproton.xi()>0.04){
      eventprotons.push_back(sproton);}}
  foreach(const phys::Proton mproton, *multiRPprotons){
    if(mproton.xi()>0.04){
      eventprotons.push_back(mproton);}}

  phys::Proton p1, p2;
  phys::ProtonPair pp;
  
  if(ZZ->passFullSelection()){
    theHistograms->fill("mZZ","mass of ZZ pair",50,300,10000,ZZ->mass(),theWeight);
    theHistograms->fill("yZZ","rapidity of ZZ pair",25,0,1,ZZ->rapidity(),theWeight);
  }
  

  for(unsigned int i=0; i<eventprotons.size();i++){
    p1=eventprotons[i];
    for(unsigned int j=i+1; j<eventprotons.size();j++){
      p2=eventprotons[j];
      pp=phys::ProtonPair(p1,p2);
      theHistograms->fill("mpp","Expected mass of ZZ pair",25,300,2800,pp.mpp(),theWeight);
      theHistograms->fill("ypp","Expected rapidity of ZZ pair",25,0,1,pp.ypp(),theWeight);

      if(ZZ->passFullSelection()){
	double massdif = 1-ZZ->mass()/pp.mpp();
	double ydif = pp.ypp()-ZZ->rapidity();
	theHistograms->fill("massdiff","Difference between pp and ZZ mass",50,-5,5,massdif,theWeight);
	theHistograms->fill("ydiff","Difference between pp and ZZ rapidity",25,-1,1,ydif,theWeight);
	theHistograms->fill("th2","2D matching distribution",50,-1.5,0.5,50,-1,1,massdif,ydif,theWeight);
	theHistograms->fill("th2xi","2D xi distribution",50,0,0.25,50,0,0.25,p1.xi(),p2.xi(),theWeight);
	
	if(abs(ydif)<0.2 && abs(massdif)<0.5){
	  theHistograms->fill("goodmZZ","mass of matched ZZ pair",10,300,2800,ZZ->mass(),theWeight);
	}
      }
    }
  }

  eventprotons.clear();
}
/*
void PPZZAnalyzer::finish(){
  TFile *file1= new TFile("th2.root","UPDATE");
  file1->Write("th2");
  file1->Close();
}
*/
