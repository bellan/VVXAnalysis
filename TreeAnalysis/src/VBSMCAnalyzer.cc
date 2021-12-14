#include "VVXAnalysis/TreeAnalysis/interface/VBSMCAnalyzer.h"
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

void VBSMCAnalyzer::ZZplots(string decay){
 
  
  if(decay == "None"){
    std::cout<<"Check decay channel"<<std::endl;
    return;
  }
  m4L_gen  = sqrt((genVBParticles->at(0).p4()+genVBParticles->at(1).p4())*(genVBParticles->at(0).p4()+genVBParticles->at(1).p4()));
  
  if (m4L_gen>=800) m4L_gen = 799;
   
  njets = genJets->size(); 

  Float_t w_kf = 1.;
  
  if((theMCInfo.fileName()=="ggZZ2e2mu") || (theMCInfo.fileName()=="ggZZ4e") || (theMCInfo.fileName()=="ggZZ4mu") || (theMCInfo.fileName()=="ggTo2e2mu_Contin_MCFM701") || (theMCInfo.fileName()=="ggTo4e_Contin_MCFM701") || (theMCInfo.fileName()=="ggTo4mu_Contin_MCFM701"))  w_kf = 1.7 ; 
  else if((theMCInfo.fileName()=="ZZTo4l")) w_kf = 1.1; 
  
  
  if(njets>=2){  
    
    deta_gen = fabs(genJets->at(0).eta() - genJets->at(1).eta());
    mjj_gen =  (genJets->at(0).p4() + genJets->at(1).p4()).M();
    
    if(mjj_gen>400 && deta_gen>2.4){   
      
      theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_01_fr", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight()*w_kf);
      
      
      if((region_ == phys::MC && regionWord.test(phys::SR4P)) || ((region_ == phys::MC_HZZ) && regionWord.test(phys::SR_HZZ))){
	
	
	theHistograms.fill(std::string("ZZTo")+decay+"_MassGenReco_01_fr", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen,theWeight*w_kf); 
	
	//SCALE FACTOR HISTOGRAMS
	
	Float_t scaleFacErrSq = ZZ->efficiencySFUnc();
	
	theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqMinus_01_fr", "", Xbins , m4L_gen, theWeight*w_kf*(1-scaleFacErrSq));   
	theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqPlus_01_fr","", Xbins , m4L_gen, theWeight*w_kf*(1+scaleFacErrSq));
	
      }
    } 
  }
}

void VBSMCAnalyzer::analyze(){
  

  if( ((region_ == phys::MC && topology.test(2)) && (region_ == phys::MC && topology.test(3)) ) || ( (region_ == phys::MC_HZZ && topology.test(0) ) &&  (region_ == phys::MC_HZZ && topology.test(1) ) ) ){       

  int z1 = abs(genVBParticles->at(0).daughter(0).id()); 
  int z2 = abs(genVBParticles->at(1).daughter(0).id()); 

  std::string decay="None";                                                                                                                                                                                

  if((z1==11 && z2==13) || (z1==13 && z2==11)) {decay = "2e2m";} 
  else if(z1==13 && z2==13) {decay = "4m";} 
  else if(z1==11 && z2==11) {decay = "4e";}
  else {cout<<"Wrong decay, check z doughters: Z0 l0"<<genVBParticles->at(0).daughter(0).id()<<" Z0 l1 "<<genVBParticles->at(0).daughter(1).id()<<" Z1 l0 "<<genVBParticles->at(1).daughter(0).id()<<" Z1 l1 "<<genVBParticles->at(1).daughter(1).id()<<endl; abort();} 
 
    ZZplots(decay);   
  }
 }  


void VBSMCAnalyzer::begin() {

  Xbins += 100,200,250,300,350,400,500,600,800; 
  m4L_gen = 0;
  mjj_gen = 0;
  deta_gen =0;
}

void VBSMCAnalyzer::end( TFile &) {

}  
