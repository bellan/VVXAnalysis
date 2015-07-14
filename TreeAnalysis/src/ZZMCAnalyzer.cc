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

void ZZMCAnalyzer::ZZplots(string decay){

 if(decay == "None"){
   std::cout<<"Check decay channel"<<std::endl;
   return;
 }
 
 string sample = "01";
 if(PreCounter < nentries/2) {sample = "0";} 
 else {sample = "1";}
 
 m4L_gen = sqrt((genVBParticles->at(0).p4()+genVBParticles->at(1).p4())*(genVBParticles->at(0).p4()+genVBParticles->at(1).p4()));
 
 Int_t njets = genJets->size();
 if (njets>3) njets=3;
 

 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight());  
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight());

 theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight());
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight());
 

 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenPU_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.weight());  
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGenPU_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.weight());
  
  if(regionWord.test(3)) {
    theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenReco_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight);  
    theHistograms.fill(std::string("ZZTo")+decay+"_MassGenReco_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen,theWeight);      
}

}
void ZZMCAnalyzer::analyze(){
 
  PreCounter+=1;


if (topology.test(0)){

    bool Ele  = 0;
    bool Muon = 0;


    foreach(const phys::Particle &gen, *genParticles)
    
      if(abs(gen.id())==13) Muon = 1;
      else if(abs(gen.id())==11) Ele = 1; 

 
 
    std::string decay="None";
    
    if(Ele&Muon)       {decay = "2e2m";} 
    else if(!Ele&Muon) {decay = "4m";}    
    else if(Ele&!Muon) {decay = "4e";}   
    else std::cout<<"NO FINALSTATE"<<std::endl;
   
	ZZplots("4l");
	ZZplots(decay);
 
 }  
}

void ZZMCAnalyzer::begin() {

  nentries =  tree()->GetEntries();
  PreCounter = 0;
  Xbins += 100,200,250,300,350,400,500,600,800;
  m4L_gen = 0;

}

void ZZMCAnalyzer::end( TFile &) {
  cout <<"Tree Entries"<<nentries<< endl;
 
}  
