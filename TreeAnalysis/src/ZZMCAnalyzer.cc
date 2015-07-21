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

 njets = genJets->size();
 if (njets>3) njets=3;
 
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight());  
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight());
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight());
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight());
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenPU_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.weight());  
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGenPU_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.weight());


 //To calculate distributions weighted for the ratio between the unfolded data and the generator MC (an early unfolding is required)
 string UnfOverMC_Jets = "ZZTo"+decay+"_Jets_Ratio"; 
 h_UnfOverMC_Jets = (TH1*) UnfOverMC->Get(UnfOverMC_Jets.c_str()); 
 int bin_Jets = h_UnfOverMC_Jets->FindBin(njets); 
 float w_Jets = h_UnfOverMC_Jets->GetBinContent(bin_Jets);
  theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_W_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight()*w_Jets);  
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_W_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight()*w_Jets);
 
 string UnfOverMC_Mass = "ZZTo"+decay+"_Mass_Ratio"; 
 h_UnfOverMC_Mass = (TH1*) UnfOverMC->Get(UnfOverMC_Mass.c_str()); 
 int bin_Mass = h_UnfOverMC_Mass->FindBin(m4L_gen); 
 float w_Mass = h_UnfOverMC_Mass->GetBinContent(bin_Mass); 
 cout  << m4L_gen << " " << bin_Mass << " " << w_Mass << " " <<  theMCInfo.sampleWeight()<< " " << theMCInfo.sampleWeight()*w_Mass << endl;
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_W_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight()*w_Mass);
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_W_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight()*w_Mass);

  
 if(regionWord.test(3)) {

   Float_t errSFLep1 =0;  Float_t errSFLep2 = 0;  Float_t errSFLep3 =0 ;  Float_t errSFLep4 = 0;   

   Float_t Lep1ID = ZZ->first().daughter(0).id();
   Float_t Lep1Pt = ZZ->first().daughter(0).pt();
   Float_t Lep1Eta = ZZ->first().daughter(0).eta();
   
   Float_t Lep2ID = ZZ->first().daughter(1).id();
   Float_t Lep2Pt = ZZ->first().daughter(1).pt();
   Float_t Lep2Eta = ZZ->first().daughter(1).eta();
   
   Float_t Lep3ID = ZZ->second().daughter(0).id();
   Float_t Lep3Pt = ZZ->second().daughter(0).pt();
   Float_t Lep3Eta = ZZ->second().daughter(0).eta();
   
   Float_t Lep4ID = ZZ->second().daughter(1).id();
   Float_t Lep4Pt = ZZ->second().daughter(1).pt();
   Float_t Lep4Eta = ZZ->second().daughter(1).eta();
   
   Float_t  SFLep1 =  lepSF.efficiencyScaleFactorErr(Lep1Pt,Lep1Eta,Lep1ID,errSFLep1);
   Float_t  SFLep2 =  lepSF.efficiencyScaleFactorErr(Lep2Pt,Lep2Eta,Lep2ID,errSFLep2);
   Float_t  SFLep3 =  lepSF.efficiencyScaleFactorErr(Lep3Pt,Lep3Eta,Lep3ID,errSFLep3);
   Float_t  SFLep4 =  lepSF.efficiencyScaleFactorErr(Lep4Pt,Lep4Eta,Lep4ID,errSFLep4);

    Float_t scaleFacErrSq = 0;

 Float_t errCorrSyst = 0;
      if(Lep1ID == 13){
	  if(Lep1Pt >= 15.) errCorrSyst = 0.000025;
	  else errCorrSyst = 0.000225;
      	  errSFLep1 = sqrt(errCorrSyst+errSFLep1*errSFLep1);
	 }

      if(Lep2ID == 13){
	if(Lep2Pt >= 15.) errCorrSyst = 0.000025;
	else errCorrSyst = 0.000225;
	errSFLep2 = sqrt(errCorrSyst+errSFLep2*errSFLep2);
	  }
      
      if(Lep3ID == 13){
	if(Lep3Pt >= 15.) errCorrSyst = 0.000025;
	else errCorrSyst = 0.000225;
	errSFLep3 = sqrt(errCorrSyst+errSFLep3*errSFLep3);
	  }
      
      if(Lep4ID == 13){
	if(Lep4Pt >= 15.) errCorrSyst = 0.000025;
	else errCorrSyst = 0.000225;
	errSFLep4 = sqrt(errCorrSyst+errSFLep4*errSFLep4);
	  }
      
      scaleFacErrSq = sqrt((errSFLep1*errSFLep1)/(SFLep1*SFLep1)+(errSFLep2*errSFLep2)/(SFLep2*SFLep2)+(errSFLep3*errSFLep3)/(SFLep3*SFLep3)+(errSFLep4*errSFLep4)/(SFLep4*SFLep4));


  theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenReco_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight);  
    theHistograms.fill(std::string("ZZTo")+decay+"_MassGenReco_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen,theWeight);      

    theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqMinus_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqMinus_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen, theWeight*(1-scaleFacErrSq));   
   

   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqPlus_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqPlus_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen, theWeight*(1+scaleFacErrSq));   


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
  //UnfOverMC = new TFile("/macros/UnfoldingMacros/UnfoldFolder/Ratio_UnfoldedDataOverGenMC.root");
  nentries =  tree()->GetEntries();
  PreCounter = 0;
  Xbins += 100,200,250,300,350,400,500,600,800;
  m4L_gen = 0;
  njets = 0;
}

void ZZMCAnalyzer::end( TFile &) {
  cout <<"Tree Entries"<<nentries<< endl;
 
}  
