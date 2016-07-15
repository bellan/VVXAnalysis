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


  bool isTightFr =kFALSE;

  nEvent ++;

 if(decay == "None"){
   std::cout<<"Check decay channel"<<std::endl;
   return;
 }
 string region;
 string sample = "01";
 if(PreCounter < nentries/2) {sample = "0";} 
 else {sample = "1";}
 
 m4L_gen = sqrt((genVBParticles->at(0).p4()+genVBParticles->at(1).p4())*(genVBParticles->at(0).p4()+genVBParticles->at(1).p4()));


 if (m4L_gen>=800) m4L_gen = 799;
 
 drzz_gen =physmath::deltaR(genVBParticles->at(0),genVBParticles->at(1));
 //if(drzz_gen>6) drzz_gen = 5.9; //overflow bin

 njets = genJets->size();
 if (njets>3) njets=3;

 ncentraljets = centralGenJets->size();
 if (ncentraljets>3) ncentraljets=3;
 
 ptzz_gen =  (genVBParticles->at(0).p4()+genVBParticles->at(1).p4()).Pt();
 if (ptzz_gen>300) ptzz_gen=299;

 dphizz_gen = 0;
 dphizz_gen = fabs(physmath::deltaPhi(genVBParticles->at(0).phi(),genVBParticles->at(1).phi())); 

 w_kf = 1; 
 if((theMCInfo.fileName()=="ggZZ2e2mu") || (theMCInfo.fileName()=="ggZZ4e") || (theMCInfo.fileName()=="ggZZ4mu"))   w_kf = kFactor_ggZZ ; 
 else if((theMCInfo.fileName()=="ZZTo4l") || (theMCInfo.fileName()=="ZZTo4lamcatnlo")) w_kf = kFactor_qqZZM * kFactor_EWKqqZZ ; 
 
 
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight()*w_kf);  
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight()*w_kf);
 theHistograms.fill(std::string("ZZTo")+decay+"_Jets_CentralGen_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theMCInfo.sampleWeight()*w_kf);  
 theHistograms.fill(std::string("ZZTo")+decay+"_Jets_CentralGen_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theMCInfo.sampleWeight()*w_kf);
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight()*w_kf);
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight()*w_kf);
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGen_"+sample,"", Xbins_drzz, drzz_gen,theMCInfo.sampleWeight()*w_kf);
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGen_01","", Xbins_drzz, drzz_gen,theMCInfo.sampleWeight()*w_kf);
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGen_"+sample,"", Xbins_ptzz, ptzz_gen,theMCInfo.sampleWeight()*w_kf);
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGen_01","", Xbins_ptzz, ptzz_gen,theMCInfo.sampleWeight()*w_kf);
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGen_"+sample,"", Xbins_dphizz, dphizz_gen,theMCInfo.sampleWeight()*w_kf);
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGen_01","", Xbins_dphizz, dphizz_gen,theMCInfo.sampleWeight()*w_kf);

 if(njets >=1){
   ptjet1_gen = genJets->at(0).pt();
   if(ptjet1_gen>=500) ptjet1_gen=499;
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1Gen_"+sample,"",Xbins_ptjet1,ptjet1_gen,theMCInfo.sampleWeight()*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1Gen_01","",Xbins_ptjet1,ptjet1_gen,theMCInfo.sampleWeight()*w_kf);
   etajet1_gen = fabs(genJets->at(0).eta());
   if(etajet1_gen>=4.7) etajet1_gen=4.6;
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1Gen_"+sample,"",Xbins_etajet1,etajet1_gen,theMCInfo.sampleWeight()*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1Gen_01","",Xbins_etajet1,etajet1_gen,theMCInfo.sampleWeight()*w_kf);
 } 

 if(njets>=2){  
  
   deta_gen = fabs(genJets->at(0).eta() - genJets->at(1).eta());
   mjj_gen =  (genJets->at(0).p4() + genJets->at(1).p4()).M();
   ptjet2_gen = genJets->at(1).pt();
   etajet2_gen = fabs(genJets->at(1).eta());
   
   if (deta_gen>=4.7) deta_gen = 4.6;
   if (mjj_gen>=800) mjj_gen = 799;
   if(ptjet2_gen>=500) ptjet2_gen=499; 
   if(etajet2_gen>=4.7) etajet2_gen=4.6;

   theHistograms.fill(std::string("ZZTo")+decay+"_MjjGen_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,theMCInfo.sampleWeight()*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjGen_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,theMCInfo.sampleWeight()*w_kf);  
   
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaGen_"+sample, std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,theMCInfo.sampleWeight()*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaGen_01", std::string("#Delta#eta__{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,theMCInfo.sampleWeight()*w_kf);

   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2Gen_"+sample,"",Xbins_ptjet2,ptjet2_gen,theMCInfo.sampleWeight()*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2Gen_01","",Xbins_ptjet2,ptjet2_gen,theMCInfo.sampleWeight()*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2Gen_"+sample,"",Xbins_etajet2,etajet2_gen,theMCInfo.sampleWeight()*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2Gen_01","",Xbins_etajet2,etajet2_gen,theMCInfo.sampleWeight()*w_kf);
}

 if(ncentraljets>=2){  
   
   deta_gen_cj = fabs(centralGenJets->at(0).eta() - centralGenJets->at(1).eta());
   mjj_gen_cj =  (centralGenJets->at(0).p4() + centralGenJets->at(1).p4()).M();
   if (deta_gen_cj>=4.7) deta_gen_cj = 4.6;
   if (mjj_gen_cj>=800) mjj_gen_cj = 799;
   
   theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_CentralGen_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theMCInfo.sampleWeight()*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_CentralGen_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theMCInfo.sampleWeight()*w_kf);  
   
   theHistograms.fill(std::string("ZZTo")+decay+"_Deta_CentralGen_"+sample, std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theMCInfo.sampleWeight()*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_Deta_CentralGen_01", std::string("#Delta#eta__{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theMCInfo.sampleWeight()*w_kf);
 
 }
 ///////////////////////////////////////////////////////////////////////IN FIDUCIAL REGION/////////////////////////////////////////////////////////

 
 //To run in the fiducial region (otherwise to be commented): 


 if(region_ == phys::MC && topology.test(3))              isTightFr = kTRUE; 
 else if (region_ == phys::MC_HZZ && topology.test(2))    isTightFr = kTRUE;        
 
 
 if(isTightFr){
  region = "_fr";
  inFiducialRegion ++;
  
  theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_"+sample+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight()*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_01"+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight()*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_Jets_CentralGen_"+sample+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theMCInfo.sampleWeight()*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_Jets_CentralGen_01"+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theMCInfo.sampleWeight()*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_"+sample+region, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight()*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_01"+region, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight()*w_kf);
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGen_"+sample+region,"", Xbins_drzz, drzz_gen,theMCInfo.sampleWeight()*w_kf);
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGen_01"+region,"", Xbins_drzz, drzz_gen,theMCInfo.sampleWeight()*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGen_"+sample+region,"", Xbins_ptzz, ptzz_gen,theMCInfo.sampleWeight()*w_kf);
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGen_01"+region,"", Xbins_ptzz, ptzz_gen,theMCInfo.sampleWeight()*w_kf);
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGen_"+sample+region,"", Xbins_dphizz, dphizz_gen,theMCInfo.sampleWeight()*w_kf);
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGen_01"+region,"", Xbins_dphizz, dphizz_gen,theMCInfo.sampleWeight()*w_kf);

   if(njets >=1){
     ptjet1_gen = genJets->at(0).pt();
     if(ptjet1_gen>=500) ptjet1_gen=499;
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1Gen_"+sample+region,"",Xbins_ptjet1,ptjet1_gen,theMCInfo.sampleWeight()*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1Gen_01"+region,"",Xbins_ptjet1,ptjet1_gen,theMCInfo.sampleWeight()*w_kf);
     etajet1_gen = fabs(genJets->at(0).eta()); 
     
     if(etajet1_gen>=4.7) etajet1_gen=4.6;
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1Gen_"+sample+region,"",Xbins_etajet1,etajet1_gen,theMCInfo.sampleWeight()*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1Gen_01"+region,"",Xbins_etajet1,etajet1_gen,theMCInfo.sampleWeight()*w_kf);
   } 
   
   if(njets>=2){  
     
     deta_gen = fabs(genJets->at(0).eta() - genJets->at(1).eta());
     mjj_gen =  (genJets->at(0).p4() + genJets->at(1).p4()).M();
     ptjet2_gen = genJets->at(1).pt();
     etajet2_gen = fabs(genJets->at(1).eta());
     
     if (deta_gen>=4.7) deta_gen = 4.6;
     if (mjj_gen>=800) mjj_gen = 799;
     if(ptjet2_gen>=500) ptjet2_gen=499; 
     if(etajet2_gen>=4.7) etajet2_gen=4.6;
 
     theHistograms.fill(std::string("ZZTo")+decay+"_MjjGen_"+sample+region, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,theMCInfo.sampleWeight()*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_MjjGen_01"+region, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,theMCInfo.sampleWeight()*w_kf);  
     
     theHistograms.fill(std::string("ZZTo")+decay+"_DetaGen_"+sample+region, std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,theMCInfo.sampleWeight()*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_DetaGen_01"+region, std::string("#Delta#eta__{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,theMCInfo.sampleWeight()*w_kf);
     
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2Gen_"+sample+region,"",Xbins_ptjet2,ptjet2_gen,theMCInfo.sampleWeight()*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2Gen_01"+region,"",Xbins_ptjet2,ptjet2_gen,theMCInfo.sampleWeight()*w_kf);
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2Gen_"+sample+region,"",Xbins_etajet2,etajet2_gen,theMCInfo.sampleWeight()*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2Gen_01"+region,"",Xbins_etajet2,etajet2_gen,theMCInfo.sampleWeight()*w_kf);
   }
   
   if(ncentraljets>=2){  
     
     deta_gen_cj = fabs(centralGenJets->at(0).eta() - centralGenJets->at(1).eta());
     mjj_gen_cj =  (centralGenJets->at(0).p4() + centralGenJets->at(1).p4()).M();
     if (deta_gen_cj>=4.7) deta_gen_cj = 4.6;
     if (mjj_gen_cj>=800) mjj_gen_cj = 799;
     
     theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_CentralGen_"+sample+region, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theMCInfo.sampleWeight()*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_CentralGen_01"+region, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theMCInfo.sampleWeight()*w_kf);  
     
     theHistograms.fill(std::string("ZZTo")+decay+"_Deta_CentralGen_"+sample+region, std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theMCInfo.sampleWeight()*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_Deta_CentralGen_01"+region, std::string("#Delta#eta__{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theMCInfo.sampleWeight()*w_kf);
     
   }
 }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


if((region_ == phys::MC && regionWord.test(26)) || ((region_ == phys::MC_HZZ) && regionWord.test(3))){


   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenReco_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*w_kf); 
   theHistograms.fill(std::string("ZZTo")+decay+"_Jets_CentralGenReco_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theWeight*w_kf);   
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGenReco_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen,theWeight*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenReco_"+sample,"", Xbins_drzz , drzz_gen,theWeight*w_kf);      
   theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenReco_"+sample,"", Xbins_ptzz , ptzz_gen,theWeight*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenReco_"+sample,"", Xbins_dphizz , dphizz_gen,theWeight*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenReco_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGenReco_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen,theWeight*w_kf); 
   theHistograms.fill(std::string("ZZTo")+decay+"_Jets_CentralGenReco_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theWeight*w_kf);      
   theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenReco_01","", Xbins_drzz , drzz_gen,theWeight*w_kf);        
   theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenReco_01","", Xbins_ptzz , ptzz_gen,theWeight*w_kf);       
   theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenReco_01","", Xbins_dphizz , dphizz_gen,theWeight*w_kf);        

   if(njets >=1){
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenReco_"+sample,"",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenReco_01","",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf);
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenReco_"+sample,"",Xbins_etajet1,etajet1_gen,theWeight*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenReco_01","",Xbins_etajet1,etajet1_gen,theWeight*w_kf);
   } 

 if(njets>=2){  
  
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenReco_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,theWeight*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenReco_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,theWeight*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenReco_"+sample, std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,theWeight*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenReco_01", std::string("#Delta#eta__{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,theWeight*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenReco_"+sample,"",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenReco_01","",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenReco_"+sample,"",Xbins_etajet2,etajet2_gen,theWeight*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenReco_01","",Xbins_etajet2,etajet2_gen,theWeight*w_kf);
}

 if(ncentraljets>=2){  
   
   theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_CentralGenReco_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theWeight*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_CentralGenReco_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theWeight*w_kf);  
   
   theHistograms.fill(std::string("ZZTo")+decay+"_Deta_CentralGenReco_"+sample, std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theWeight*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_Deta_CentralGenReco_01", std::string("#Delta#eta__{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theWeight*w_kf);
 
 }

 
 //SCALE FACTOR HISTOGRAMS

 Float_t scaleFacErrSq = ZZ->efficiencySFUnc();

 theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqMinus_"+sample, "", Xbins , m4L_gen, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqMinus_01", "", Xbins , m4L_gen, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqPlus_"+sample, "", Xbins , m4L_gen, theWeight*w_kf*(1+scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqPlus_01","", Xbins , m4L_gen, theWeight*w_kf*(1+scaleFacErrSq));
 
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqMinus_"+sample,"",4,0,4,njets,theWeight*w_kf*(1-scaleFacErrSq)); 
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqMinus_01","",4,0,4,njets,theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqPlus_"+sample,"",4,0,4,njets,theWeight*w_kf*(1+scaleFacErrSq));  
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqPlus_01","",4,0,4,njets,theWeight*w_kf*(1+scaleFacErrSq)); 
 
 theHistograms.fill(std::string("ZZTo")+decay+"_Jets_CentralGenRecoSFErrSqMinus_"+sample,"",4,0,4,ncentraljets,theWeight*w_kf*(1-scaleFacErrSq)); 
 theHistograms.fill(std::string("ZZTo")+decay+"_Jets_CentralGenRecoSFErrSqMinus_01","",4,0,4,ncentraljets,theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_Jets_CentralGenRecoSFErrSqPlus_"+sample,"",4,0,4,ncentraljets,theWeight*w_kf*(1+scaleFacErrSq));  
 theHistograms.fill(std::string("ZZTo")+decay+"_Jets_CentralGenRecoSFErrSqPlus_01", "",4,0,4,ncentraljets,theWeight*w_kf*(1+scaleFacErrSq));  
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenRecoSFErrSqMinus_"+sample, "", Xbins_drzz , drzz_gen, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenRecoSFErrSqMinus_01", "", Xbins_drzz , drzz_gen, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenRecoSFErrSqPlus_"+sample, "", Xbins_drzz , drzz_gen, theWeight*w_kf*(1+scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenRecoSFErrSqPlus_01","", Xbins_drzz , drzz_gen, theWeight*w_kf*(1+scaleFacErrSq));

 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenRecoSFErrSqMinus_"+sample, "", Xbins_ptzz , ptzz_gen, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenRecoSFErrSqMinus_01", "", Xbins_ptzz , ptzz_gen, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenRecoSFErrSqPlus_"+sample, "", Xbins_ptzz , ptzz_gen, theWeight*w_kf*(1+scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenRecoSFErrSqPlus_01","", Xbins_ptzz , ptzz_gen, theWeight*w_kf*(1+scaleFacErrSq));
   
  theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenRecoSFErrSqMinus_"+sample, "", Xbins_dphizz , dphizz_gen, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenRecoSFErrSqMinus_01", "", Xbins_dphizz , dphizz_gen, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenRecoSFErrSqPlus_"+sample, "", Xbins_dphizz , dphizz_gen, theWeight*w_kf*(1+scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenRecoSFErrSqPlus_01","", Xbins_dphizz , dphizz_gen, theWeight*w_kf*(1+scaleFacErrSq)); 



 if(njets >=1){
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrSqMinus_"+sample,"",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrSqMinus_01","",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErrSq));
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrSqPlus_"+sample,"",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrSqPlus_01","",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErrSq));
   
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrSqMinus_"+sample,"",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrSqMinus_01","",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErrSq));
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrSqPlus_"+sample,"",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrSqPlus_01","",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErrSq));
   
 } 
 
 if(njets>=2){  
   // cout <<  genJets->at(0).motherId() << " " << genJets->at(1).motherId() << endl;
   
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenRecoSFErrSqMinus_"+sample,"",Xbins_mjj,mjj_gen,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenRecoSFErrSqMinus_01","",Xbins_mjj,mjj_gen,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenRecoSFErrSqPlus_"+sample,"",Xbins_mjj,mjj_gen,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenRecoSFErrSqPlus_01","",Xbins_mjj,mjj_gen,theWeight*w_kf*(1+scaleFacErrSq));  
   
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenRecoSFErrSqMinus_"+sample,"",Xbins_deta,deta_gen,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenRecoSFErrSqMinus_01","",Xbins_deta,deta_gen,theWeight*w_kf*(1-scaleFacErrSq));
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenRecoSFErrSqPlus_"+sample,"",Xbins_deta,deta_gen,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenRecoSFErrSqPlus_01","",Xbins_deta,deta_gen,theWeight*w_kf*(1+scaleFacErrSq));
   
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenRecoSFErrSqMinus_"+sample,"",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenRecoSFErrSqMinus_01","",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf*(1-scaleFacErrSq));
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenRecoSFErrSqPlus_"+sample,"",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenRecoSFErrSqPlus_01","",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf*(1+scaleFacErrSq));
   
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenRecoSFErrSqMinus_"+sample,"",Xbins_etajet2,etajet2_gen,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenRecoSFErrSqMinus_01","",Xbins_etajet2,etajet2_gen,theWeight*w_kf*(1-scaleFacErrSq));
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenRecoSFErrSqPlus_"+sample,"",Xbins_etajet2,etajet2_gen,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenRecoSFErrSqPlus_01","",Xbins_etajet2,etajet2_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
   
    
 }
 
 if(ncentraljets>=2){  
   
   theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_CentralGenRecoSFErrSqMinus_"+sample,"",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_CentralGenRecoSFErrSqMinus_01","",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_CentralGenRecoSFErrSqPlus_"+sample,"",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_CentralGenRecoSFErrSqPlus_01","",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
   
   theHistograms.fill(std::string("ZZTo")+decay+"_Deta_CentralGenRecoSFErrSqMinus_"+sample,"",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_Deta_CentralGenRecoSFErrSqMinus_01","",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1-scaleFacErrSq));
   theHistograms.fill(std::string("ZZTo")+decay+"_Deta_CentralGenRecoSFErrSqPlus_"+sample,"",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_Deta_CentralGenRecoSFErrSqPlus_01","",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
   
 }
 
 ///////////////////////////////////////////////////////////////////////IN FIDUCIAL REGION/////////////////////////////////////////////////////////
 
 //To run in the fiducial region (otherwise to be commented): 


  if(isTightFr){
   
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenReco_"+sample+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenReco_01"+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_Jets_CentralGenReco_"+sample+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theWeight*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_Jets_CentralGenReco_01"+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theWeight*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGenReco_"+sample+region, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theWeight*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGenReco_01"+region, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theWeight*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenReco_"+sample+region, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins_drzz, drzz_gen,theWeight*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenReco_01"+region, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins_drzz, drzz_gen,theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenReco_"+sample+region, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins_ptzz, ptzz_gen,theWeight*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenReco_01"+region, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins_ptzz, ptzz_gen,theWeight*w_kf); 
   theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenReco_"+sample+region, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins_dphizz, dphizz_gen,theWeight*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenReco_01"+region, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins_dphizz, dphizz_gen,theWeight*w_kf);
   
   if(njets >=1){
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenReco_"+sample+region,"",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenReco_01"+region,"",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf);
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenReco_"+sample+region,"",Xbins_etajet1,etajet1_gen,theWeight*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenReco_01"+region,"",Xbins_etajet1,etajet1_gen,theWeight*w_kf);
   } 
   
   if(njets>=2){  
     
     theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenReco_"+sample+region, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,theWeight*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenReco_01"+region, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,theWeight*w_kf);  
     
     theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenReco_"+sample+region, std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,theWeight*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenReco_01"+region, std::string("#Delta#eta__{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,theWeight*w_kf);
     
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenReco_"+sample+region,"",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenReco_01"+region,"",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf);
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenReco_"+sample+region,"",Xbins_etajet2,etajet2_gen,theWeight*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenReco_01"+region,"",Xbins_etajet2,etajet2_gen,theWeight*w_kf);
   }
   
   if(ncentraljets>=2){  
     
     theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_CentralGenReco_"+sample+region, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theWeight*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_CentralGenReco_01"+region, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theWeight*w_kf);  
     
     theHistograms.fill(std::string("ZZTo")+decay+"_Deta_CentralGenReco_"+sample+region, std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theWeight*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_Deta_CentralGenReco_01"+region, std::string("#Delta#eta__{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theWeight*w_kf);
     
   }
   
   //SCALE FACTOR HISTOGRAMS
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqMinus_"+sample+region, "", Xbins , m4L_gen, theWeight*w_kf*(1-scaleFacErrSq));   
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqMinus_01_fr", "", Xbins , m4L_gen, theWeight*w_kf*(1-scaleFacErrSq));   
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqPlus_"+sample+region, "", Xbins , m4L_gen, theWeight*w_kf*(1+scaleFacErrSq));   
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqPlus_01_fr","", Xbins , m4L_gen, theWeight*w_kf*(1+scaleFacErrSq));
   
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqMinus_"+sample+region,"",4,0,4,njets,theWeight*w_kf*(1-scaleFacErrSq)); 
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqMinus_01_fr","",4,0,4,njets,theWeight*w_kf*(1-scaleFacErrSq));   
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqPlus_"+sample+region,"",4,0,4,njets,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqPlus_01_fr","",4,0,4,njets,theWeight*w_kf*(1+scaleFacErrSq)); 
   
   theHistograms.fill(std::string("ZZTo")+decay+"_Jets_CentralGenRecoSFErrSqMinus_"+sample+region,"",4,0,4,ncentraljets,theWeight*w_kf*(1-scaleFacErrSq)); 
   theHistograms.fill(std::string("ZZTo")+decay+"_Jets_CentralGenRecoSFErrSqMinus_01_fr","",4,0,4,ncentraljets,theWeight*w_kf*(1-scaleFacErrSq));   
   theHistograms.fill(std::string("ZZTo")+decay+"_Jets_CentralGenRecoSFErrSqPlus_"+sample+region,"",4,0,4,ncentraljets,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_Jets_CentralGenRecoSFErrSqPlus_01_fr", "",4,0,4,ncentraljets,theWeight*w_kf*(1+scaleFacErrSq));  

   theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenRecoSFErrSqMinus_"+sample+region, "", Xbins_drzz , drzz_gen, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenRecoSFErrSqMinus_01_fr", "", Xbins_drzz , drzz_gen, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenRecoSFErrSqPlus_"+sample+region, "", Xbins_drzz , drzz_gen, theWeight*w_kf*(1+scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenRecoSFErrSqPlus_01_fr","", Xbins_drzz , drzz_gen, theWeight*w_kf*(1+scaleFacErrSq));   

 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenRecoSFErrSqMinus_"+sample+region, "", Xbins_ptzz , ptzz_gen, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenRecoSFErrSqMinus_01_fr", "", Xbins_ptzz , ptzz_gen, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenRecoSFErrSqPlus_"+sample+region, "", Xbins_ptzz , ptzz_gen, theWeight*w_kf*(1+scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenRecoSFErrSqPlus_01_fr","", Xbins_ptzz , ptzz_gen, theWeight*w_kf*(1+scaleFacErrSq)); 

 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenRecoSFErrSqMinus_"+sample+region, "", Xbins_dphizz , dphizz_gen, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenRecoSFErrSqMinus_01_fr", "", Xbins_dphizz , dphizz_gen, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenRecoSFErrSqPlus_"+sample+region, "", Xbins_dphizz , dphizz_gen, theWeight*w_kf*(1+scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenRecoSFErrSqPlus_01_fr","", Xbins_dphizz , dphizz_gen, theWeight*w_kf*(1+scaleFacErrSq)); 


   if(njets >=1){
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrSqMinus_"+sample+region,"",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrSqMinus_01_fr","",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErrSq));
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrSqPlus_"+sample+region,"",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrSqPlus_01_fr","",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErrSq));
     
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrSqMinus_"+sample+region,"",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrSqMinus_01_fr","",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErrSq));
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrSqPlus_"+sample+region,"",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrSqPlus_01_fr","",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErrSq));
     

   } 
   
   if(njets>=2){  
     // cout <<  genJets->at(0).motherId() << " " << genJets->at(1).motherId() << endl;
     
     theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenRecoSFErrSqMinus_"+sample+region,"",Xbins_mjj,mjj_gen,theWeight*w_kf*(1-scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenRecoSFErrSqMinus_01_fr","",Xbins_mjj,mjj_gen,theWeight*w_kf*(1-scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenRecoSFErrSqPlus_"+sample+region,"",Xbins_mjj,mjj_gen,theWeight*w_kf*(1+scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenRecoSFErrSqPlus_01_fr","",Xbins_mjj,mjj_gen,theWeight*w_kf*(1+scaleFacErrSq));  
     
     theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenRecoSFErrSqMinus_"+sample+region,"",Xbins_deta,deta_gen,theWeight*w_kf*(1-scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenRecoSFErrSqMinus_01_fr","",Xbins_deta,deta_gen,theWeight*w_kf*(1-scaleFacErrSq));
     theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenRecoSFErrSqPlus_"+sample+region,"",Xbins_deta,deta_gen,theWeight*w_kf*(1+scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenRecoSFErrSqPlus_01_fr","",Xbins_deta,deta_gen,theWeight*w_kf*(1+scaleFacErrSq));
     
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenRecoSFErrSqMinus_"+sample+region,"",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf*(1-scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenRecoSFErrSqMinus_01_fr","",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf*(1-scaleFacErrSq));
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenRecoSFErrSqPlus_"+sample+region,"",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf*(1+scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenRecoSFErrSqPlus_01_fr","",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf*(1+scaleFacErrSq));
     
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenRecoSFErrSqMinus_"+sample+region,"",Xbins_etajet2,etajet2_gen,theWeight*w_kf*(1-scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenRecoSFErrSqMinus_01_fr","",Xbins_etajet2,etajet2_gen,theWeight*w_kf*(1-scaleFacErrSq));
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenRecoSFErrSqPlus_"+sample+region,"",Xbins_etajet2,etajet2_gen,theWeight*w_kf*(1+scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenRecoSFErrSqPlus_01_fr","",Xbins_etajet2,etajet2_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
     
 
   }
   
   if(ncentraljets>=2){  
     
     theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_CentralGenRecoSFErrSqMinus_"+sample+region,"",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_CentralGenRecoSFErrSqMinus_01_fr","",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_CentralGenRecoSFErrSqPlus_"+sample+region,"",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_CentralGenRecoSFErrSqPlus_01_fr","",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
     
     theHistograms.fill(std::string("ZZTo")+decay+"_Deta_CentralGenRecoSFErrSqMinus_"+sample+region,"",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1-scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_Deta_CentralGenRecoSFErrSqMinus_01_fr","",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1-scaleFacErrSq));
     theHistograms.fill(std::string("ZZTo")+decay+"_Deta_CentralGenRecoSFErrSqPlus_"+sample+region,"",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1+scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_Deta_CentralGenRecoSFErrSqPlus_01_fr","",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
     
   }
  }
 }
}

void ZZMCAnalyzer::analyze(){
  
  PreCounter+=1;
   
  if((region_ == phys::MC && topology.test(2)) || (region_ == phys::MC_HZZ && topology.test(0) ) ){       

  int Ele  = 0;
  int Muon = 0; 
  int lep = 0;
  
  foreach(const phys::Particle &gen, *genParticles){
    if((!gen.genStatusFlags().test(GenStatusBit::isPrompt)) || (!gen.genStatusFlags().test(GenStatusBit::fromHardProcess)))  continue;
    if(abs(gen.id())==13){	
      Muon += 1; 
      lep +=1;
      }
      else if(abs(gen.id())==11){
	Ele += 1;
	lep+=1; 
      }
     }

    std::string decay="None";
    if(Ele==Muon && Ele!=0)  {decay = "2e2m";} 
    else if(Ele<Muon) {decay = "4m";}    
    else if(Ele>Muon) {decay = "4e";}   
    else if(Ele==Muon && Ele == 0) std::cout<<"NO FINALSTATE"<<std::endl;
    else std::cout<<"NO FINALSTATE"<<std::endl;   

    if(lep>4) std::cout<<"To many leptons: "<<lep<<std::endl; 
    else if(lep<4)  std::cout<<"To few leptons: "<<lep<<std::endl; 

    //ZZplots("4l");
    ZZplots(decay);   
 }  
}

void ZZMCAnalyzer::begin() {
  nentries =  tree()->GetEntries();
  PreCounter = 0;
  nEvent = 0;
  inFiducialRegion=0;
  Xbins += 100,200,250,300,350,400,500,600,800; 
  Xbins_deta += 0,2.4,4.7;
  Xbins_mjj += 0.,200,800;
  Xbins_ptjet1 += 30,50,100,200,300,500;
  Xbins_ptjet2 += 30,100,200,500;
  Xbins_etajet1 += 0,1.5,3,4.7;
  Xbins_etajet2 += 0,1.5,3,4.7; 
  Xbins_drzz += 0,1,2,3,4,5,6;
  Xbins_ptzz += 0,25,50,75,100,150,200,300;
  Xbins_dphizz += 0,1.5,2.,2.25,2.5,2.75,3,3.25;
  
  m4L_gen = 0;
  njets = 0;
  mjj_gen = 0;
  deta_gen = 0;
  ncentraljets = 0;
  mjj_gen_cj = 0;
  deta_gen_cj = 0;
  ptjet1_gen = 0;
  ptjet2_gen = 0; 
  etajet1_gen = 0;
  etajet2_gen = 0; 
  drzz_gen =0;
  ptzz_gen =0;
  dphizz_gen =0;
}

void ZZMCAnalyzer::end( TFile &) {
  cout <<"Tree Entries"<<nentries<< endl;
  cout <<"events in the fiducial region"<<inFiducialRegion<< endl;
}  
