#include "VVXAnalysis/TreeAnalysis/interface/ZZMCAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"

#include <boost/assign/std/vector.hpp> 
#include <boost/assert.hpp> 
using namespace std;
using namespace boost::assign; // bring 'operator+=()' into scope

#include <boost/foreach.hpp>
#include <sstream> 
#include <string> 
#define foreach BOOST_FOREACH

using std::cout;
using std::endl;


using namespace phys;

void ZZMCAnalyzer::ZZplots(int id){

  if(ZZ->id() != id && id != -1) return; // -1 here means generic 4l final state

  string decay  = "4l";
  string sample = "01";

  if      (id == 52) {decay = "4m";}
  else if (id == 48) {decay = "2e2m";}
  else if (id == 44) {decay = "4e";}

  // if(e < nentries/2) {sample = "0";} 
  // else {sample = "1";}

  m4L_gen = sqrt((genVBParticles->at(0).p4()+genVBParticles->at(1).p4())*(genVBParticles->at(0).p4()+genVBParticles->at(1).p4()));

 theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theWeight);
 theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_"+sample, std::string("Response matrix invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, Xbins, ZZ->mass() ,m4L_gen , theWeight);
 
}


void ZZMCAnalyzer::analyze(){
 
   e++; 
  //cout << "n event: " << e << endl;
  nentries =  tree()->GetEntries();

  if(topology.test(0)){
    ZZplots();   // ZZ --> 4l
    ZZplots(52); // ZZ --> 4m
    ZZplots(48); // ZZ --> 2e2m
    ZZplots(44); // ZZ --> 4e
  }
}


void ZZMCAnalyzer::begin() {

  Xbins += 100,200,250,300,350,400,500,600,800;
  m4L_gen =0;
 
}

void ZZMCAnalyzer::end( TFile &) {
  cout <<nentries<< endl;
  cout << e << endl;
}  
