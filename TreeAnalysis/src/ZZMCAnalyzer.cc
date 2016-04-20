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
 
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight()*w_kf);  
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight()*w_kf);
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGen_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theMCInfo.sampleWeight()*w_kf);  
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGen_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theMCInfo.sampleWeight()*w_kf);
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
   
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGen_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theMCInfo.sampleWeight()*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGen_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theMCInfo.sampleWeight()*w_kf);  
   
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGen_"+sample, std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theMCInfo.sampleWeight()*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGen_01", std::string("#Delta#eta__{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theMCInfo.sampleWeight()*w_kf);
 
 }
 ///////////////////////////////////////////////////////////////////////IN FIDUCIAL REGION/////////////////////////////////////////////////////////
 zz::SignalTopology zzSignalTopology = zz::getSignalTopology(*genParticles, *genJets);
 
 //To run in the fiducial region (otherwise to be commented): 
 if(zz::inTightFiducialRegion(zzSignalTopology)){
   region = "_fr";
   inFiducialRegion ++;
   
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_"+sample+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight()*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_01"+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight()*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGen_"+sample+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theMCInfo.sampleWeight()*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGen_01"+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theMCInfo.sampleWeight()*w_kf);
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
     
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGen_"+sample+region, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theMCInfo.sampleWeight()*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGen_01"+region, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theMCInfo.sampleWeight()*w_kf);  
     
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGen_"+sample+region, std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theMCInfo.sampleWeight()*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGen_01"+region, std::string("#Delta#eta__{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theMCInfo.sampleWeight()*w_kf);
     
   }
 }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


 if(regionWord.test(3)) {
   
   //Float_t errSFLep1 =0;  Float_t errSFLep2 = 0;  Float_t errSFLep3 =0 ;  Float_t errSFLep4 = 0;   
   
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
   
   // Float_t  SFLep1 =  lepSF.efficiencyScaleFactorErr(Lep1Pt,Lep1Eta,Lep1ID,errSFLep1);
   // Float_t  SFLep2 =  lepSF.efficiencyScaleFactorErr(Lep2Pt,Lep2Eta,Lep2ID,errSFLep2);
   // Float_t  SFLep3 =  lepSF.efficiencyScaleFactorErr(Lep3Pt,Lep3Eta,Lep3ID,errSFLep3);
   // Float_t  SFLep4 =  lepSF.efficiencyScaleFactorErr(Lep4Pt,Lep4Eta,Lep4ID,errSFLep4);
   
   Float_t  SFLep1 =  lepSF.efficiencyScaleFactor(Lep1Pt,Lep1Eta,Lep1ID);
   Float_t  SFLep2 =  lepSF.efficiencyScaleFactor(Lep2Pt,Lep2Eta,Lep2ID);
   Float_t  SFLep3 =  lepSF.efficiencyScaleFactor(Lep3Pt,Lep3Eta,Lep3ID);
   Float_t  SFLep4 =  lepSF.efficiencyScaleFactor(Lep4Pt,Lep4Eta,Lep4ID);
   
   Float_t errSFLep1 =  lepSF.efficiencyScaleFactorErr(Lep1Pt,Lep1Eta,Lep1ID);
   Float_t errSFLep2 =  lepSF.efficiencyScaleFactorErr(Lep2Pt,Lep2Eta,Lep2ID);
   Float_t errSFLep3 =  lepSF.efficiencyScaleFactorErr(Lep3Pt,Lep3Eta,Lep3ID);
   Float_t errSFLep4 =  lepSF.efficiencyScaleFactorErr(Lep4Pt,Lep4Eta,Lep4ID);
   
   Float_t scaleFacErr = 0;
   
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
   scaleFacErr = errSFLep1/SFLep1+errSFLep2/SFLep2+errSFLep3/SFLep3+errSFLep4/(SFLep4);
   
   
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenReco_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*w_kf); 
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenReco_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theWeight*w_kf);   
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGenReco_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen,theWeight*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenReco_"+sample,"", Xbins_drzz , drzz_gen,theWeight*w_kf);      
   theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenReco_"+sample,"", Xbins_ptzz , ptzz_gen,theWeight*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenReco_"+sample,"", Xbins_dphizz , dphizz_gen,theWeight*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenReco_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGenReco_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay+"of reco events" , Xbins , m4L_gen,theWeight*w_kf); 
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenReco_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theWeight*w_kf);      
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
   
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenReco_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theWeight*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenReco_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theWeight*w_kf);  
   
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenReco_"+sample, std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theWeight*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenReco_01", std::string("#Delta#eta__{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theWeight*w_kf);
 
 }
 
 //SCALE FACTOR HISTOGRAMS
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqMinus_"+sample, "", Xbins , m4L_gen, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqMinus_01", "", Xbins , m4L_gen, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqPlus_"+sample, "", Xbins , m4L_gen, theWeight*w_kf*(1+scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrSqPlus_01","", Xbins , m4L_gen, theWeight*w_kf*(1+scaleFacErrSq));
 
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqMinus_"+sample,"",4,0,4,njets,theWeight*w_kf*(1-scaleFacErrSq)); 
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqMinus_01","",4,0,4,njets,theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqPlus_"+sample,"",4,0,4,njets,theWeight*w_kf*(1+scaleFacErrSq));  
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqPlus_01","",4,0,4,njets,theWeight*w_kf*(1+scaleFacErrSq)); 
 
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenRecoSFErrSqMinus_"+sample,"",4,0,4,ncentraljets,theWeight*w_kf*(1-scaleFacErrSq)); 
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenRecoSFErrSqMinus_01","",4,0,4,ncentraljets,theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenRecoSFErrSqPlus_"+sample,"",4,0,4,ncentraljets,theWeight*w_kf*(1+scaleFacErrSq));  
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenRecoSFErrSqPlus_01", "",4,0,4,ncentraljets,theWeight*w_kf*(1+scaleFacErrSq));  
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

 // Histograms for scale factor correlated errors
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrMinus_"+sample,"", Xbins , m4L_gen, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrMinus_01","", Xbins , m4L_gen, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrPlus_"+sample,"", Xbins , m4L_gen, theWeight*w_kf*(1+scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrPlus_01","", Xbins , m4L_gen, theWeight*w_kf*(1+scaleFacErr)); 
 
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrMinus_"+sample,"",4,0,4,njets,theWeight*w_kf*(1-scaleFacErr));  
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrMinus_01","",4,0,4,njets,theWeight*w_kf*(1-scaleFacErr));  
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrPlus_"+sample,"",4,0,4,njets,theWeight*w_kf*(1+scaleFacErr));  
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrPlus_01","",4,0,4,njets,theWeight*w_kf*(1+scaleFacErr));   
 
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenRecoSFErrMinus_"+sample,"",4,0,4,ncentraljets,theWeight*w_kf*(1-scaleFacErr));  
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenRecoSFErrMinus_01","",4,0,4,ncentraljets,theWeight*w_kf*(1-scaleFacErr));  
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenRecoSFErrPlus_"+sample,"",4,0,4,ncentraljets,theWeight*w_kf*(1+scaleFacErr));  
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenRecoSFErrPlus_01","",4,0,4,ncentraljets,theWeight*w_kf*(1+scaleFacErr));  
 
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenRecoSFErrMinus_"+sample,"", Xbins_drzz , drzz_gen, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenRecoSFErrMinus_01","", Xbins_drzz , drzz_gen, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenRecoSFErrPlus_"+sample,"", Xbins_drzz , drzz_gen, theWeight*w_kf*(1+scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenRecoSFErrPlus_01","", Xbins_drzz , drzz_gen, theWeight*w_kf*(1+scaleFacErr)); 

  theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenRecoSFErrMinus_"+sample,"", Xbins_ptzz , ptzz_gen, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenRecoSFErrMinus_01","", Xbins_ptzz , ptzz_gen, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenRecoSFErrPlus_"+sample,"", Xbins_ptzz , ptzz_gen, theWeight*w_kf*(1+scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenRecoSFErrPlus_01","", Xbins_ptzz , ptzz_gen, theWeight*w_kf*(1+scaleFacErr));

 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenRecoSFErrMinus_"+sample,"", Xbins_dphizz , dphizz_gen, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenRecoSFErrMinus_01","", Xbins_dphizz , dphizz_gen, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenRecoSFErrPlus_"+sample,"", Xbins_dphizz , dphizz_gen, theWeight*w_kf*(1+scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenRecoSFErrPlus_01","", Xbins_dphizz , dphizz_gen, theWeight*w_kf*(1+scaleFacErr));

 if(njets >=1){
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrSqMinus_"+sample,"",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrSqMinus_01","",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErrSq));
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrSqPlus_"+sample,"",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrSqPlus_01","",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErrSq));
   
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrSqMinus_"+sample,"",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrSqMinus_01","",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErrSq));
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrSqPlus_"+sample,"",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrSqPlus_01","",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErrSq));
   
   // Histograms for scale factor correlated errors
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrMinus_"+sample,"",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrMinus_01","",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErr));
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrPlus_"+sample,"",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrPlus_01","",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErr));
   
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrMinus_"+sample,"",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrMinus_01","",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErr));
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrPlus_"+sample,"",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrPlus_01","",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErr));
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
   
   // Histograms for scale factor correlated errors
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenRecoSFErrMinus_"+sample,"",Xbins_mjj,mjj_gen,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenRecoSFErrMinus_01","",Xbins_mjj,mjj_gen,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenRecoSFErrPlus_"+sample,"",Xbins_mjj,mjj_gen,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenRecoSFErrPlus_01","",Xbins_mjj,mjj_gen,theWeight*w_kf*(1+scaleFacErr));  
   
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenRecoSFErrMinus_"+sample,"",Xbins_deta,deta_gen,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenRecoSFErrMinus_01","",Xbins_deta,deta_gen,theWeight*w_kf*(1-scaleFacErr));
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenRecoSFErrPlus_"+sample,"",Xbins_deta,deta_gen,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenRecoSFErrPlus_01","",Xbins_deta,deta_gen,theWeight*w_kf*(1+scaleFacErr));
   
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenRecoSFErrMinus_"+sample,"",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenRecoSFErrMinus_01","",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf*(1-scaleFacErr));
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenRecoSFErrPlus_"+sample,"",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenRecoSFErrPlus_01","",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf*(1+scaleFacErr));
   
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenRecoSFErrMinus_"+sample,"",Xbins_etajet2,etajet2_gen,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenRecoSFErrMinus_01","",Xbins_etajet2,etajet2_gen,theWeight*w_kf*(1-scaleFacErr));
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenRecoSFErrPlus_"+sample,"",Xbins_etajet2,etajet2_gen,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenRecoSFErrPlus_01","",Xbins_etajet2,etajet2_gen,theWeight*w_kf*(1+scaleFacErr)); 
 }
 
 if(ncentraljets>=2){  
   
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenRecoSFErrSqMinus_"+sample,"",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenRecoSFErrSqMinus_01","",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenRecoSFErrSqPlus_"+sample,"",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenRecoSFErrSqPlus_01","",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
   
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenRecoSFErrSqMinus_"+sample,"",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenRecoSFErrSqMinus_01","",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1-scaleFacErrSq));
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenRecoSFErrSqPlus_"+sample,"",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenRecoSFErrSqPlus_01","",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
   
   // Histograms for scale factor correlated errors
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenRecoSFErrMinus_"+sample,"",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenRecoSFErrMinus_01","",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenRecoSFErrPlus_"+sample,"",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenRecoSFErrPlus_01","",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErr)); 
   
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenRecoSFErrMinus_"+sample,"",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenRecoSFErrMinus_01","",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1-scaleFacErr));
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenRecoSFErrPlus_"+sample,"",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenRecoSFErrPlus_01","",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1+scaleFacErr));
 }
 
 ///////////////////////////////////////////////////////////////////////IN FIDUCIAL REGION/////////////////////////////////////////////////////////
 
 //To run in the fiducial region (otherwise to be commented): 
 if(zz::inTightFiducialRegion(zzSignalTopology)){
   
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenReco_"+sample+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenReco_01"+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*w_kf);
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenReco_"+sample+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theWeight*w_kf);  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenReco_01"+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theWeight*w_kf);
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
     
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenReco_"+sample+region, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theWeight*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenReco_01"+region, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,theWeight*w_kf);  
     
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenReco_"+sample+region, std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theWeight*w_kf);  
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenReco_01"+region, std::string("#Delta#eta__{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,theWeight*w_kf);
     
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
   
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenRecoSFErrSqMinus_"+sample+region,"",4,0,4,ncentraljets,theWeight*w_kf*(1-scaleFacErrSq)); 
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenRecoSFErrSqMinus_01_fr","",4,0,4,ncentraljets,theWeight*w_kf*(1-scaleFacErrSq));   
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenRecoSFErrSqPlus_"+sample+region,"",4,0,4,ncentraljets,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenRecoSFErrSqPlus_01_fr", "",4,0,4,ncentraljets,theWeight*w_kf*(1+scaleFacErrSq));  

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

   // Histograms for scale factor correlated errors
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrMinus_"+sample+region,"", Xbins , m4L_gen, theWeight*w_kf*(1-scaleFacErr));   
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrMinus_01_fr","", Xbins , m4L_gen, theWeight*w_kf*(1-scaleFacErr));   
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrPlus_"+sample+region,"", Xbins , m4L_gen, theWeight*w_kf*(1+scaleFacErr));   
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGenRecoSFErrPlus_01_fr","", Xbins , m4L_gen, theWeight*w_kf*(1+scaleFacErr)); 
   
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrMinus_"+sample+region,"",4,0,4,njets,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrMinus_01_fr","",4,0,4,njets,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrPlus_"+sample+region,"",4,0,4,njets,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrPlus_01_fr","",4,0,4,njets,theWeight*w_kf*(1+scaleFacErr));   
   
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenRecoSFErrMinus_"+sample+region,"",4,0,4,ncentraljets,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenRecoSFErrMinus_01_fr","",4,0,4,ncentraljets,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenRecoSFErrPlus_"+sample+region,"",4,0,4,ncentraljets,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGenRecoSFErrPlus_01_fr","",4,0,4,ncentraljets,theWeight*w_kf*(1+scaleFacErr));  
   
    theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenRecoSFErrMinus_"+sample+region,"", Xbins_drzz , drzz_gen, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenRecoSFErrMinus_01_fr","", Xbins_drzz , drzz_gen, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenRecoSFErrPlus_"+sample+region,"", Xbins_drzz , drzz_gen, theWeight*w_kf*(1+scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZGenRecoSFErrPlus_01_fr","", Xbins_drzz , drzz_gen, theWeight*w_kf*(1+scaleFacErr)); 

    theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenRecoSFErrMinus_"+sample+region,"", Xbins_ptzz , ptzz_gen, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenRecoSFErrMinus_01_fr","", Xbins_ptzz , ptzz_gen, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenRecoSFErrPlus_"+sample+region,"", Xbins_ptzz , ptzz_gen, theWeight*w_kf*(1+scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZGenRecoSFErrPlus_01_fr","", Xbins_ptzz , ptzz_gen, theWeight*w_kf*(1+scaleFacErr)); 

  theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenRecoSFErrMinus_"+sample+region,"", Xbins_dphizz , dphizz_gen, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenRecoSFErrMinus_01_fr","", Xbins_dphizz , dphizz_gen, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenRecoSFErrPlus_"+sample+region,"", Xbins_dphizz , dphizz_gen, theWeight*w_kf*(1+scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZGenRecoSFErrPlus_01_fr","", Xbins_dphizz , dphizz_gen, theWeight*w_kf*(1+scaleFacErr)); 

   if(njets >=1){
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrSqMinus_"+sample+region,"",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrSqMinus_01_fr","",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErrSq));
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrSqPlus_"+sample+region,"",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrSqPlus_01_fr","",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErrSq));
     
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrSqMinus_"+sample+region,"",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrSqMinus_01_fr","",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErrSq));
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrSqPlus_"+sample+region,"",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrSqPlus_01_fr","",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErrSq));
     
     // Histograms for scale factor correlated errors
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrMinus_"+sample+region,"",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrMinus_01_fr","",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErr));
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrPlus_"+sample+region,"",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1GenRecoSFErrPlus_01_fr","",Xbins_ptjet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErr));
     
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrMinus_"+sample+region,"",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrMinus_01_fr","",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErr));
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrPlus_"+sample+region,"",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1GenRecoSFErrPlus_01_fr","",Xbins_etajet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErr));
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
     
     // Histograms for scale factor correlated errors
     theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenRecoSFErrMinus_"+sample+region,"",Xbins_mjj,mjj_gen,theWeight*w_kf*(1-scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenRecoSFErrMinus_01_fr","",Xbins_mjj,mjj_gen,theWeight*w_kf*(1-scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenRecoSFErrPlus_"+sample+region,"",Xbins_mjj,mjj_gen,theWeight*w_kf*(1+scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_MjjGenRecoSFErrPlus_01_fr","",Xbins_mjj,mjj_gen,theWeight*w_kf*(1+scaleFacErr));  
     
     theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenRecoSFErrMinus_"+sample+region,"",Xbins_deta,deta_gen,theWeight*w_kf*(1-scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenRecoSFErrMinus_01_fr","",Xbins_deta,deta_gen,theWeight*w_kf*(1-scaleFacErr));
     theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenRecoSFErrPlus_"+sample+region,"",Xbins_deta,deta_gen,theWeight*w_kf*(1+scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_DetaGenRecoSFErrPlus_01_fr","",Xbins_deta,deta_gen,theWeight*w_kf*(1+scaleFacErr));
     
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenRecoSFErrMinus_"+sample+region,"",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf*(1-scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenRecoSFErrMinus_01_fr","",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf*(1-scaleFacErr));
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenRecoSFErrPlus_"+sample+region,"",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf*(1+scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2GenRecoSFErrPlus_01_fr","",Xbins_ptjet2,ptjet2_gen,theWeight*w_kf*(1+scaleFacErr));
     
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenRecoSFErrMinus_"+sample+region,"",Xbins_etajet2,etajet2_gen,theWeight*w_kf*(1-scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenRecoSFErrMinus_01_fr","",Xbins_etajet2,etajet2_gen,theWeight*w_kf*(1-scaleFacErr));
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenRecoSFErrPlus_"+sample+region,"",Xbins_etajet2,etajet2_gen,theWeight*w_kf*(1+scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2GenRecoSFErrPlus_01_fr","",Xbins_etajet2,etajet2_gen,theWeight*w_kf*(1+scaleFacErr)); 
   }
   
   if(ncentraljets>=2){  
     
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenRecoSFErrSqMinus_"+sample+region,"",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenRecoSFErrSqMinus_01_fr","",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenRecoSFErrSqPlus_"+sample+region,"",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenRecoSFErrSqPlus_01_fr","",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
     
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenRecoSFErrSqMinus_"+sample+region,"",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1-scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenRecoSFErrSqMinus_01_fr","",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1-scaleFacErrSq));
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenRecoSFErrSqPlus_"+sample+region,"",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1+scaleFacErrSq));  
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenRecoSFErrSqPlus_01_fr","",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
     
     // Histograms for scale factor correlated errors
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenRecoSFErrMinus_"+sample+region,"",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenRecoSFErrMinus_01_fr","",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenRecoSFErrPlus_"+sample+region,"",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGenRecoSFErrPlus_01_fr","",Xbins_mjj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErr)); 
     
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenRecoSFErrMinus_"+sample+region,"",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1-scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenRecoSFErrMinus_01_fr","",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1-scaleFacErr));
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenRecoSFErrPlus_"+sample+region,"",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1+scaleFacErr));  
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGenRecoSFErrPlus_01_fr","",Xbins_deta,deta_gen_cj,theWeight*w_kf*(1+scaleFacErr));
   }
 }
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
