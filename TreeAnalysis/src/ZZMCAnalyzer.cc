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
 // if(PreCounter < nentries/2) {sample = "0";} 
 //else {sample = "1";}

 m4L_gen = sqrt((genVBParticles->at(0).p4()+genVBParticles->at(1).p4())*(genVBParticles->at(0).p4()+genVBParticles->at(1).p4()));
 if (m4L_gen>=800) m4L_gen = 799;
 

 // //Unfolding and jets part commented for now

 // njets = genJets->size();
 //if (njets>3) njets=3;
 // ncentraljets = centralGenJets->size();
 //if (ncentraljets>3) ncentraljets=3; 
//  theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight());  
//  theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight());
// theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGen_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theMCInfo.sampleWeight());  
//  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGen_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theMCInfo.sampleWeight());
//  theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight());

 theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight());

 

 // //Unfolding and jets part commented for now

 //To calculate distributions weighted for the ratio between the unfolded data and the generator MC (an early unfolding is required)
 // string UnfOverMC_Jets = "ZZTo"+decay+"_Jets_Ratio"; 
 // h_UnfOverMC_Jets = (TH1*) UnfOverMC->Get(UnfOverMC_Jets.c_str());  
 // h_UnfOverMC_Jets->Draw();
 // int bin_Jets = h_UnfOverMC_Jets->FindBin(njets);
 //  float w_Jets = h_UnfOverMC_Jets->GetBinContent(bin_Jets);
 //  theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_W_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight()*w_Jets);  
 // theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_W_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight()*w_Jets);
 
 // string UnfOverMC_Mass = "ZZTo"+decay+"_Mass_Ratio"; 
 // h_UnfOverMC_Mass = (TH1*) UnfOverMC_Pow->Get(UnfOverMC_Mass.c_str()); 
 // int bin_Mass = h_UnfOverMC_Mass->FindBin(m4L_gen); 
 // float w_Mass = h_UnfOverMC_Mass->GetBinContent(bin_Mass); 
 // cout  << m4L_gen << " " << bin_Mass << " " << w_Mass << " " <<  theMCInfo.sampleWeight()<< " " << theMCInfo.sampleWeight()*w_Mass << endl;
 // theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_W_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight()*w_Mass);
 // theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_W_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight()*w_Mass);

 // if(njets>=2){  
   
 //   deta_gen = fabs(genJets->at(0).eta() - genJets->at(1).eta());
 //   mjj_gen =  (genJets->at(0).p4() + genJets->at(1).p4()).M();
 //   if (deta_gen>=4.7) deta_gen = 4.6;
 //   if (mjj_gen>=800) mjj_gen = 799;
   
 //   theHistograms.fill(std::string("ZZTo")+decay+"_MjjGen_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,theMCInfo.sampleWeight());  
 //   theHistograms.fill(std::string("ZZTo")+decay+"_MjjGen_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,theMCInfo.sampleWeight());  
   
 //   theHistograms.fill(std::string("ZZTo")+decay+"_DetaGen_"+sample, std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,theMCInfo.sampleWeight());  
 //   theHistograms.fill(std::string("ZZTo")+decay+"_DetaGen_01", std::string("#Delta#eta__{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,theMCInfo.sampleWeight());
 
   // string UnfOverMC_Mjj = "ZZTo"+decay+"_Mjj_Ratio"; 
   // h_UnfOverMC_Mjj = (TH1*) UnfOverMC->Get(UnfOverMC_Mjj.c_str());  
   // h_UnfOverMC_Mjj->Draw();
   // int bin_Mjj = h_UnfOverMC_Mjj->FindBin(mjj_gen);
   // float w_Mjj = h_UnfOverMC_Mjj->GetBinContent(bin_Mjj);
   // theHistograms.fill(std::string("ZZTo")+decay+"_MjjGen_W_"+sample, std::string("m_{jj}of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,theMCInfo.sampleWeight()*w_Mjj);  
   // theHistograms.fill(std::string("ZZTo")+decay+"_MjjGen_W_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,theMCInfo.sampleWeight()*w_Mjj);  
   
   // string UnfOverMC_Deta = "ZZTo"+decay+"_Deta_Ratio"; 
   // h_UnfOverMC_Deta = (TH1*) UnfOverMC->Get(UnfOverMC_Deta.c_str());  
   // h_UnfOverMC_Deta->Draw();
   // int bin_Deta = h_UnfOverMC_Deta->FindBin(deta_gen);
   // float w_Deta = h_UnfOverMC_Deta->GetBinContent(bin_Deta);
   // theHistograms.fill(std::string("ZZTo")+decay+"_DetaGen_W_"+sample, std::string("#Delta#eta_{jj}of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,theMCInfo.sampleWeight()*w_Deta);  
   // theHistograms.fill(std::string("ZZTo")+decay+"_DetaGen_W_01", std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,theMCInfo.sampleWeight()*w_Deta);  
 
 //}

 // if(ncentraljets>=2){  
   
 //   deta_gen_cj = fabs(centralGenJets->at(0).eta() - centralGenJets->at(1).eta());
 //   mjj_gen_cj =  (centralGenJets->at(0).p4() + centralGenJets->at(1).p4()).M();
 //   if (deta_gen_cj>=4.7) deta_gen_cj = 4.6;
 //   if (mjj_gen_cj>=800) mjj_gen_cj = 799;
   
 //   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGen_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theMCInfo.sampleWeight());  
 //   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGen_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theMCInfo.sampleWeight());  
   
 //   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGen_"+sample, std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theMCInfo.sampleWeight());  
 //   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGen_01", std::string("#Delta#eta__{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theMCInfo.sampleWeight());
 
 // }

 // if(regionWord.test(3)) { //HZZ
  if(regionWord.test(26)) { //ZZ


    Float_t  scaleFacErrSq = ZZ->efficiencySFUnc() ;
      
      // theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenReco_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight);  
      // theHistograms.fill(std::string("ZZTo")+decay+"_MassGenReco_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen,theWeight);      
      // theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqMinus_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*(1-scaleFacErrSq));  
      // theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqMinus_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen, theWeight*(1-scaleFacErrSq));   
      
      // theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqPlus_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*(1+scaleFacErrSq));  
      // theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqPlus_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen, theWeight*(1+scaleFacErrSq));   


      theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenReco_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight);  
      theHistograms.fill(std::string("ZZTo")+decay+"_MassGenReco_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen,theWeight);      

      theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqMinus_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*(1-scaleFacErrSq));  
      theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqMinus_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen, theWeight*(1-scaleFacErrSq));   
      
      theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqPlus_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*(1+scaleFacErrSq));  
      theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqPlus_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen, theWeight*(1+scaleFacErrSq));   

      // Histograms for scale factor correlated errors
			     
      // theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrMinus_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*(1-scaleFacErrSq));  
      // theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrMinus_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen, theWeight*(1-scaleFacErr));   
      
      // theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrPlus_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*(1+scaleFacErrSq));  
      // theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrPlus_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen, theWeight*(1+scaleFacErr));   
            
      // theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrMinus_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*(1-scaleFacErrSq));  
      // theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrMinus_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen, theWeight*(1-scaleFacErr));   
      
      // theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrPlus_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*(1+scaleFacErrSq));  
      // theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrPlus_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen, theWeight*(1+scaleFacErr));  
      
  }

}
void ZZMCAnalyzer::analyze(){
 
  PreCounter+=1;

  Float_t Lep1Pt = ZZ->first().daughter(0).pt();
  Float_t Lep2Pt = ZZ->first().daughter(1).pt();
  Float_t Lep3Pt = ZZ->second().daughter(0).pt();
  Float_t Lep4Pt = ZZ->second().daughter(1).pt();

  Float_t Lep1Eta = ZZ->first().daughter(0).eta();
  Float_t Lep2Eta = ZZ->first().daughter(1).eta();
  Float_t Lep3Eta = ZZ->second().daughter(0).eta();
  Float_t Lep4Eta = ZZ->second().daughter(1).eta();
  
  bool isFR;

  if( Lep1Pt>10. && Lep2Pt>10. && Lep3Pt>10. && Lep4Pt>10. && Lep1Eta<2.5 && Lep2Eta<2.5 && Lep3Eta < 2.5 && Lep4Eta < 2.5 )  isFR=1;

  if (topology.test(0) && isFR){

    int Ele  = 0;
    int Muon = 0; 
    int lep = 0;

    //bool bad1 = 0;
    //bool bad2 = 0;

    std::string eventstr=std::to_string(run)+":"+std::to_string(lumiBlock)+":"+std::to_string(event);
    //    std::cout<<"\nnew ev "<<std::endl;
    //    std::cout<<"\n"<<eventstr<<std::endl;

    foreach(const phys::Particle &gen, *genParticles){


      if((!gen.genStatusFlags().test(GenStatusBit::isPrompt)) || (!gen.genStatusFlags().test(GenStatusBit::isHardProcess)))  continue;
   
      //if(!gen.genStatusFlags().test(GenStatusBit::isPrompt))  continue;
 
 
      if(abs(gen.id())==13){
	Muon += 1; 
	lep +=1;
	//if(!gen.genStatusFlags().test(GenStatusBit::isHardProcess)) bad1=1;//std::cout<<"muon hard\n";  //std::cout<<"gen id "<< gen.id()<<" status "<<gen.genStatusFlags()<<" mom "<<gen.motherId()<<std::endl;   
      }
      else if(abs(gen.id())==11){
	Ele += 1;
	lep+=1; 
	//if(gen.genStatusFlags().test(GenStatusBit::isHardProcess)) bad2=1; //std::cout<<"ele hard\n";//bad1=1;// std::cout<<" is not hard "<<std::endl;  //std::cout<<"gen id "<< gen.id()<<" status "<<gen.genStatusFlags()<<" mom "<<gen.motherId()<<std::endl;   
      }
      // if(lep>4) std::cout<<"more leptons\n";
    } 
    std::string decay="None";
    
    if(Ele==Muon)       {decay = "2e2m";} 
    else if(Ele<Muon) {decay = "4m";}    
    else if(Ele>Muon) {decay = "4e";}   
    else std::cout<<"NO FINALSTATE"<<std::endl;
    //if(lep>4) std::cout<<"more leptons\n "<<eventstr<<" decay "<<decay<<std::endl;

    //ZZplots("4l");
    ZZplots(decay);    
 }  
}

void ZZMCAnalyzer::begin() {
  //  UnfOverMC = new TFile("macros/UnfoldingMacros/UnfoldFolder/Ratio_UnfoldedDataOverGenMC.root");
  // UnfOverMC_Pow = new TFile("macros/UnfoldingMacros/UnfoldFolder_Pow/Ratio_UnfoldedDataOverGenMC.root");
  nentries =  tree()->GetEntries();
  PreCounter = 0;
  Xbins += 100,200,250,300,350,400,500,600,800; 
 
  // Xbins_deta += 0,2.4,4.7;
  // Xbins_mjj += 0.,200,800;
  // m4L_gen = 0;
  // njets = 0;
  // mjj_gen = 0;
  // deta_gen = 0;
  // ncentraljets = 0;
  // mjj_gen_cj = 0;
  // deta_gen_cj = 0;
}

void ZZMCAnalyzer::end( TFile &) {
  cout <<"Tree Entries"<<nentries<< endl;
 
}  
