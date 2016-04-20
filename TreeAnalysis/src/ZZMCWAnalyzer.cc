#include "VVXAnalysis/TreeAnalysis/interface/ZZMCWAnalyzer.h"
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

void ZZMCWAnalyzer::ZZplots(string decay){

 if(decay == "None"){
   std::cout<<"Check decay channel"<<std::endl;
   return;
 }
 
 string sample = "01";
 if(PreCounter < nentries/2) {sample = "0";} 
 else {sample = "1";}
 
 m4L_gen = sqrt((genVBParticles->at(0).p4()+genVBParticles->at(1).p4())*(genVBParticles->at(0).p4()+genVBParticles->at(1).p4()));
 if (m4L_gen>=800) m4L_gen = 799;
 
 njets = genJets->size();
 if (njets>3) njets=3;

 ncentraljets = centralGenJets->size();
 if (ncentraljets>3) ncentraljets=3;
 
 if(decay == "4m"){
   bin_Jets = h_UnfOverMC_Jets_4m->FindBin(njets);  
   w_Jets =  h_UnfOverMC_Jets_4m->GetBinContent(bin_Jets);
   bin_Mass = h_UnfOverMC_Mass_4m->FindBin(m4L_gen);  
   w_Mass =  h_UnfOverMC_Mass_4m->GetBinContent(bin_Mass);
   bin_CentralJets = h_UnfOverMC_CentralJets_4m->FindBin(ncentraljets);  
   w_CentralJets =  h_UnfOverMC_CentralJets_4m->GetBinContent(bin_CentralJets);
 }  
 if(decay == "4e"){
   bin_Jets = h_UnfOverMC_Jets_4e->FindBin(njets);  
   w_Jets =  h_UnfOverMC_Jets_4e->GetBinContent(bin_Jets);
   bin_Mass = h_UnfOverMC_Mass_4e->FindBin(m4L_gen);  
   w_Mass =  h_UnfOverMC_Mass_4e->GetBinContent(bin_Mass);
   bin_CentralJets = h_UnfOverMC_CentralJets_4e->FindBin(ncentraljets);  
   w_CentralJets =  h_UnfOverMC_CentralJets_4e->GetBinContent(bin_CentralJets);
 }  
 if(decay == "2e2m"){
   bin_Jets = h_UnfOverMC_Jets_2e2m->FindBin(njets);  
   w_Jets =  h_UnfOverMC_Jets_2e2m->GetBinContent(bin_Jets);
   bin_Mass = h_UnfOverMC_Mass_2e2m->FindBin(m4L_gen);  
   w_Mass =  h_UnfOverMC_Mass_2e2m->GetBinContent(bin_Mass);
   bin_CentralJets = h_UnfOverMC_CentralJets_2e2m->FindBin(ncentraljets);  
   w_CentralJets =  h_UnfOverMC_CentralJets_2e2m->GetBinContent(bin_CentralJets);
 }  
 
 //To calculate distributions weighted for the ratio between the unfolded data and the generator MC (an earlier unfolding is required)
 //JETS
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_W_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight()*w_Jets);  
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_W_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight()*w_Jets);
 
 //MASS
 // cout  << m4L_gen << " " << bin_Mass << " " << w_Mass << " " <<  theMCInfo.sampleWeight()<< " " << theMCInfo.sampleWeight()*w_Mass << endl;
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_W_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight()*w_Mass);
 theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_W_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight()*w_Mass);

 //CENTRALJETS
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGen_W_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theMCInfo.sampleWeight()*w_CentralJets);  
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGen_W_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theMCInfo.sampleWeight()*w_CentralJets);
 
 //PTJET1 and ETAJET1
 if(njets>=1){ 
   ptjet1_gen = genJets->at(0).pt();
   if(ptjet1_gen >=500) ptjet1_gen =499;
   etajet1_gen = fabs(genJets->at(0).eta());
   if(etajet1_gen >=6) etajet1_gen =5;
   
   if(decay == "4m"){
     bin_PtJet1 = h_UnfOverMC_PtJet1_4m->FindBin(ptjet1_gen);  
     w_PtJet1 =  h_UnfOverMC_PtJet1_4m->GetBinContent(bin_PtJet1);
     bin_EtaJet1 = h_UnfOverMC_EtaJet1_4m->FindBin(etajet1_gen);  
     w_EtaJet1 =  h_UnfOverMC_EtaJet1_4m->GetBinContent(bin_EtaJet1);
   } 
   if(decay == "4e"){
     bin_PtJet1 = h_UnfOverMC_PtJet1_4e->FindBin(ptjet1_gen);  
     w_PtJet1 =  h_UnfOverMC_PtJet1_4e->GetBinContent(bin_PtJet1);
     bin_EtaJet1 = h_UnfOverMC_EtaJet1_4e->FindBin(etajet1_gen);  
     w_EtaJet1 =  h_UnfOverMC_EtaJet1_4e->GetBinContent(bin_EtaJet1);
   } 
   if(decay == "2e2m"){
     bin_PtJet1 = h_UnfOverMC_PtJet1_2e2m->FindBin(ptjet1_gen);  
     w_PtJet1 =  h_UnfOverMC_PtJet1_2e2m->GetBinContent(bin_PtJet1);
     bin_EtaJet1 = h_UnfOverMC_EtaJet1_2e2m->FindBin(etajet1_gen);  
     w_EtaJet1 =  h_UnfOverMC_EtaJet1_2e2m->GetBinContent(bin_EtaJet1);
   } 
   
   //PTJET1
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1Gen_W_"+sample, "p_T^{jet1}", Xbins_ptjet1,ptjet1_gen,theWeight*w_PtJet1);
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1Gen_W_01", "p_T^{jet1}", Xbins_ptjet1,ptjet1_gen,theWeight*w_PtJet1);
   
   //ETAJET1
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1Gen_W_"+sample, "#eta^{jet1}", Xbins_etajet1,etajet1_gen,theWeight*w_EtaJet1);
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1Gen_W_01", "#eta^{jet1}", Xbins_etajet1,etajet1_gen,theWeight*w_EtaJet1);
 }


 //MJJ AND DETA
 if(njets>=2){  
   
   deta_gen = fabs(genJets->at(0).eta() - genJets->at(1).eta());
   mjj_gen =  (genJets->at(0).p4() + genJets->at(1).p4()).M();
   ptjet2_gen = genJets->at(1).pt(); 
   etajet2_gen = fabs(genJets->at(1).eta());
   if (deta_gen>=4.7) deta_gen = 4.6;
   if (mjj_gen>=800) mjj_gen = 799;
   if(ptjet2_gen >=500) ptjet2_gen =499;
   if(etajet2_gen >=6) etajet2_gen =5;

   if(decay == "4m"){
     bin_Mjj = h_UnfOverMC_Mjj_4m->FindBin(mjj_gen);  
     w_Mjj =  h_UnfOverMC_Mjj_4m->GetBinContent(bin_Mjj);
     bin_Deta = h_UnfOverMC_Deta_4m->FindBin(deta_gen);  
     w_Deta =  h_UnfOverMC_Deta_4m->GetBinContent(bin_Deta);
     bin_PtJet2 = h_UnfOverMC_PtJet2_4m->FindBin(ptjet2_gen);  
     w_PtJet2 =  h_UnfOverMC_PtJet2_4m->GetBinContent(bin_PtJet2);
     bin_EtaJet2 = h_UnfOverMC_EtaJet2_4m->FindBin(etajet2_gen);  
     w_EtaJet2 =  h_UnfOverMC_EtaJet2_4m->GetBinContent(bin_EtaJet2);
   } 
   if(decay == "4e"){
     bin_Mjj = h_UnfOverMC_Mjj_4e->FindBin(mjj_gen);  
     w_Mjj =  h_UnfOverMC_Mjj_4e->GetBinContent(bin_Mjj);
     bin_Deta = h_UnfOverMC_Deta_4e->FindBin(deta_gen);  
     w_Deta =  h_UnfOverMC_Deta_4e->GetBinContent(bin_Deta);
     bin_PtJet2 = h_UnfOverMC_PtJet2_4e->FindBin(ptjet2_gen);  
     w_PtJet2 =  h_UnfOverMC_PtJet2_4e->GetBinContent(bin_PtJet2);
     bin_EtaJet2 = h_UnfOverMC_EtaJet2_4e->FindBin(etajet2_gen);  
     w_EtaJet2 =  h_UnfOverMC_EtaJet2_4e->GetBinContent(bin_EtaJet2);
   } 
   if(decay == "2e2m"){
     bin_Mjj = h_UnfOverMC_Mjj_2e2m->FindBin(mjj_gen);  
     w_Mjj =  h_UnfOverMC_Mjj_2e2m->GetBinContent(bin_Mjj);
     bin_Deta = h_UnfOverMC_Deta_2e2m->FindBin(deta_gen);  
     w_Deta =  h_UnfOverMC_Deta_2e2m->GetBinContent(bin_Deta);
     bin_PtJet2 = h_UnfOverMC_PtJet2_2e2m->FindBin(ptjet2_gen);  
     w_PtJet2 =  h_UnfOverMC_PtJet2_2e2m->GetBinContent(bin_PtJet2);
     bin_EtaJet2 = h_UnfOverMC_EtaJet2_2e2m->FindBin(etajet2_gen);  
     w_EtaJet2 =  h_UnfOverMC_EtaJet2_2e2m->GetBinContent(bin_EtaJet2);
   } 

   //MJJ
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjGen_W_"+sample, std::string("m_{jj}of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,theMCInfo.sampleWeight()*w_Mjj);  
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjGen_W_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,theMCInfo.sampleWeight()*w_Mjj);  
   
   //DETA
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaGen_W_"+sample, std::string("#Delta#eta_{jj}of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,theMCInfo.sampleWeight()*w_Deta);  
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaGen_W_01", std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,theMCInfo.sampleWeight()*w_Deta);  

   //PTJET2
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2Gen_W_"+sample, "p_T^{jet2}", Xbins_ptjet2,ptjet2_gen,theWeight*w_PtJet2);
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2Gen_W_01", "p_T^{jet2}", Xbins_ptjet2,ptjet2_gen,theWeight*w_PtJet2);
   
   //ETAJET2
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2Gen_W_"+sample, "#eta^{jet2}", Xbins_etajet2,etajet2_gen,theWeight*w_EtaJet2);
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2Gen_W_01", "#eta^{jet2}", Xbins_etajet2,etajet2_gen,theWeight*w_EtaJet2);
 
 }

 //CENTRALMJJ AND DETA
 if(ncentraljets>=2){  
   
   deta_gen_cj = fabs(centralGenJets->at(0).eta() - centralGenJets->at(1).eta());
   mjj_gen_cj =  (centralGenJets->at(0).p4() + centralGenJets->at(1).p4()).M();
   if (deta_gen_cj>=4.7) deta_gen_cj = 4.6;
   if (mjj_gen_cj>=800) mjj_gen_cj = 799;
   
   if(decay == "4m"){
     bin_CentralMjj = h_UnfOverMC_CentralMjj_4m->FindBin(mjj_gen_cj);  
     w_CentralMjj =  h_UnfOverMC_CentralMjj_4m->GetBinContent(bin_CentralMjj);
     bin_CentralDeta = h_UnfOverMC_CentralDeta_4m->FindBin(deta_gen_cj);  
     w_CentralDeta =  h_UnfOverMC_CentralDeta_4m->GetBinContent(bin_CentralDeta);
   } 
   if(decay == "4e"){
     bin_CentralMjj = h_UnfOverMC_CentralMjj_4e->FindBin(mjj_gen_cj);  
     w_CentralMjj =  h_UnfOverMC_CentralMjj_4e->GetBinContent(bin_CentralMjj);
     bin_CentralDeta = h_UnfOverMC_CentralDeta_4e->FindBin(deta_gen_cj);  
     w_CentralDeta =  h_UnfOverMC_CentralDeta_4e->GetBinContent(bin_CentralDeta);
   } 
   if(decay == "2e2m"){
     bin_CentralMjj = h_UnfOverMC_CentralMjj_2e2m->FindBin(mjj_gen_cj);  
     w_CentralMjj =  h_UnfOverMC_CentralMjj_2e2m->GetBinContent(bin_CentralMjj);
     bin_CentralDeta = h_UnfOverMC_CentralDeta_2e2m->FindBin(deta_gen_cj);  
     w_CentralDeta =  h_UnfOverMC_CentralDeta_2e2m->GetBinContent(bin_CentralDeta);
   }  
   
   //CENTRALMJJ
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGen_W_"+sample, "m_{jj}", Xbins_mjj,mjj_gen_cj,theMCInfo.sampleWeight()*w_CentralMjj);
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGen_W_01", "m_{jj}", Xbins_mjj, mjj_gen_cj,theMCInfo.sampleWeight()*w_CentralMjj);
   
   //CENTRALDETA
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGen_W_"+sample, "#Delta#eta_{jj}", Xbins_deta,deta_gen_cj,theMCInfo.sampleWeight()*w_CentralDeta);
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGen_W_01", "#Delta#eta_{jj}", Xbins_deta,deta_gen_cj,theMCInfo.sampleWeight()*w_CentralDeta);
 }

 //////////////////////////////////////////////IN TIGHT FIDUCIAL REGION///////////////////////////////////////////
  zz::SignalTopology zzSignalTopology = zz::getSignalTopology(*genParticles, *genJets);
  string region;
 if(zz::inTightFiducialRegion(zzSignalTopology)){
   region = "_fr";
   
   if(decay == "4m"){
     bin_Jets_fr = h_UnfOverMC_fr_Jets_4m->FindBin(njets);  
     w_Jets_fr =  h_UnfOverMC_fr_Jets_4m->GetBinContent(bin_Jets_fr);
     bin_Mass_fr = h_UnfOverMC_fr_Mass_4m->FindBin(m4L_gen);  
     w_Mass_fr =  h_UnfOverMC_fr_Mass_4m->GetBinContent(bin_Mass_fr);
     bin_CentralJets_fr = h_UnfOverMC_fr_CentralJets_4m->FindBin(ncentraljets);  
     w_CentralJets_fr =  h_UnfOverMC_fr_CentralJets_4m->GetBinContent(bin_CentralJets_fr);
   }  
   if(decay == "4e"){
     bin_Jets_fr = h_UnfOverMC_fr_Jets_4e->FindBin(njets);  
     w_Jets_fr =  h_UnfOverMC_fr_Jets_4e->GetBinContent(bin_Jets_fr);
     bin_Mass_fr = h_UnfOverMC_fr_Mass_4e->FindBin(m4L_gen);  
     w_Mass_fr =  h_UnfOverMC_fr_Mass_4e->GetBinContent(bin_Mass_fr);
     bin_CentralJets_fr = h_UnfOverMC_fr_CentralJets_4e->FindBin(ncentraljets);  
     w_CentralJets_fr =  h_UnfOverMC_fr_CentralJets_4e->GetBinContent(bin_CentralJets_fr);
   }  
   if(decay == "2e2m"){
     bin_Jets_fr = h_UnfOverMC_fr_Jets_2e2m->FindBin(njets);  
     w_Jets_fr =  h_UnfOverMC_fr_Jets_2e2m->GetBinContent(bin_Jets_fr);
     bin_Mass_fr = h_UnfOverMC_fr_Mass_2e2m->FindBin(m4L_gen);  
     w_Mass_fr =  h_UnfOverMC_fr_Mass_2e2m->GetBinContent(bin_Mass_fr);
     bin_CentralJets_fr = h_UnfOverMC_fr_CentralJets_2e2m->FindBin(ncentraljets);  
     w_CentralJets_fr =  h_UnfOverMC_fr_CentralJets_2e2m->GetBinContent(bin_CentralJets_fr);
   }  
   
   //To calculate distributions weighted for the ratio between the unfolded data and the generator MC (an early unfolding is required)
   //JETS
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_W_"+sample+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight()*w_Jets_fr);  
   theHistograms.fill(std::string("ZZTo")+decay+"_JetsGen_W_01"+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theMCInfo.sampleWeight()*w_Jets_fr);
   
   //MASS
   // cout  << m4L_gen << " " << bin_Mass << " " << w_Mass << " " <<  theMCInfo.sampleWeight()<< " " << theMCInfo.sampleWeight()*w_Mass << endl;
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_W_"+sample+region, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight()*w_Mass_fr);
   theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_W_01"+region, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight()*w_Mass_fr);
   
   //CENTRALJETS
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGen_W_"+sample+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theMCInfo.sampleWeight()*w_CentralJets_fr);  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsGen_W_01"+region, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theMCInfo.sampleWeight()*w_CentralJets_fr);
 
   //PTJET1 and ETAJET1
   if(njets>=1){ 
       
     if(decay == "4m"){
       bin_PtJet1_fr = h_UnfOverMC_fr_PtJet1_4m->FindBin(ptjet1_gen);  
       w_PtJet1_fr =  h_UnfOverMC_fr_PtJet1_4m->GetBinContent(bin_PtJet1_fr);
       bin_EtaJet1_fr = h_UnfOverMC_fr_EtaJet1_4m->FindBin(etajet1_gen);  
       w_EtaJet1_fr =  h_UnfOverMC_fr_EtaJet1_4m->GetBinContent(bin_EtaJet1_fr);
     } 
     if(decay == "4e"){
       bin_PtJet1_fr = h_UnfOverMC_fr_PtJet1_4e->FindBin(ptjet1_gen);  
       w_PtJet1_fr =  h_UnfOverMC_fr_PtJet1_4e->GetBinContent(bin_PtJet1_fr);
       bin_EtaJet1_fr = h_UnfOverMC_fr_EtaJet1_4e->FindBin(etajet1_gen);  
       w_EtaJet1_fr =  h_UnfOverMC_fr_EtaJet1_4e->GetBinContent(bin_EtaJet1_fr);
     } 
     if(decay == "2e2m"){
       bin_PtJet1_fr = h_UnfOverMC_fr_PtJet1_2e2m->FindBin(ptjet1_gen);  
       w_PtJet1_fr =  h_UnfOverMC_fr_PtJet1_2e2m->GetBinContent(bin_PtJet1_fr);
       bin_EtaJet1_fr = h_UnfOverMC_fr_EtaJet1_2e2m->FindBin(etajet1_gen);  
       w_EtaJet1_fr =  h_UnfOverMC_fr_EtaJet1_2e2m->GetBinContent(bin_EtaJet1_fr);
     } 
     
     //PTJET1
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1Gen_W_"+sample+region, "p_T^{jet1}", Xbins_ptjet1,ptjet1_gen,theWeight*w_PtJet1_fr);
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1Gen_W_01"+region, "p_T^{jet1}", Xbins_ptjet1,ptjet1_gen,theWeight*w_PtJet1_fr);
     
   //ETAJET1
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1Gen_W_"+sample+region, "#eta^{jet1}", Xbins_etajet1,etajet1_gen,theWeight*w_EtaJet1_fr);
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1Gen_W_01"+region, "#eta^{jet1}", Xbins_etajet1,etajet1_gen,theWeight*w_EtaJet1_fr);
   }
   
   
   //MJJ AND DETA
   if(njets>=2){  
     
     if(decay == "4m"){
       bin_Mjj_fr = h_UnfOverMC_fr_Mjj_4m->FindBin(mjj_gen);  
       w_Mjj_fr =  h_UnfOverMC_fr_Mjj_4m->GetBinContent(bin_Mjj_fr);
       bin_Deta_fr = h_UnfOverMC_fr_Deta_4m->FindBin(deta_gen);  
       w_Deta_fr =  h_UnfOverMC_fr_Deta_4m->GetBinContent(bin_Deta_fr);
       bin_PtJet2_fr = h_UnfOverMC_fr_PtJet2_4m->FindBin(ptjet2_gen);  
       w_PtJet2_fr =  h_UnfOverMC_fr_PtJet2_4m->GetBinContent(bin_PtJet2_fr);
       bin_EtaJet2_fr = h_UnfOverMC_fr_EtaJet2_4m->FindBin(etajet2_gen);  
       w_EtaJet2_fr =  h_UnfOverMC_fr_EtaJet2_4m->GetBinContent(bin_EtaJet2_fr);
   } 
     if(decay == "4e"){
       bin_Mjj_fr = h_UnfOverMC_fr_Mjj_4e->FindBin(mjj_gen);  
       w_Mjj_fr =  h_UnfOverMC_fr_Mjj_4e->GetBinContent(bin_Mjj_fr);
       bin_Deta_fr = h_UnfOverMC_fr_Deta_4e->FindBin(deta_gen);  
       w_Deta_fr =  h_UnfOverMC_fr_Deta_4e->GetBinContent(bin_Deta_fr);
       bin_PtJet2_fr = h_UnfOverMC_fr_PtJet2_4e->FindBin(ptjet2_gen);  
       w_PtJet2_fr =  h_UnfOverMC_fr_PtJet2_4e->GetBinContent(bin_PtJet2_fr);
       bin_EtaJet2_fr = h_UnfOverMC_fr_EtaJet2_4e->FindBin(etajet2_gen);  
       w_EtaJet2_fr =  h_UnfOverMC_fr_EtaJet2_4e->GetBinContent(bin_EtaJet2_fr);
     } 
     if(decay == "2e2m"){
       bin_Mjj_fr = h_UnfOverMC_fr_Mjj_2e2m->FindBin(mjj_gen);  
       w_Mjj_fr =  h_UnfOverMC_fr_Mjj_2e2m->GetBinContent(bin_Mjj_fr);
       bin_Deta_fr = h_UnfOverMC_fr_Deta_2e2m->FindBin(deta_gen);  
       w_Deta_fr =  h_UnfOverMC_fr_Deta_2e2m->GetBinContent(bin_Deta_fr);
       bin_PtJet2_fr = h_UnfOverMC_fr_PtJet2_2e2m->FindBin(ptjet2_gen);  
       w_PtJet2_fr =  h_UnfOverMC_fr_PtJet2_2e2m->GetBinContent(bin_PtJet2_fr);
       bin_EtaJet2_fr = h_UnfOverMC_fr_EtaJet2_2e2m->FindBin(etajet2_gen);  
       w_EtaJet2_fr =  h_UnfOverMC_fr_EtaJet2_2e2m->GetBinContent(bin_EtaJet2_fr);
     } 
     
     //MJJ
     theHistograms.fill(std::string("ZZTo")+decay+"_MjjGen_W_"+sample+region, std::string("m_{jj}of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,theMCInfo.sampleWeight()*w_Mjj_fr);  
     theHistograms.fill(std::string("ZZTo")+decay+"_MjjGen_W_01"+region, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,theMCInfo.sampleWeight()*w_Mjj_fr);  
     
     //DETA
     theHistograms.fill(std::string("ZZTo")+decay+"_DetaGen_W_"+sample+region, std::string("#Delta#eta_{jj}of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,theMCInfo.sampleWeight()*w_Deta_fr);  
     theHistograms.fill(std::string("ZZTo")+decay+"_DetaGen_W_01"+region, std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,theMCInfo.sampleWeight()*w_Deta_fr);  
     
     //PTJET2
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2Gen_W_"+sample+region, "p_T^{jet2}", Xbins_ptjet2,ptjet2_gen,theWeight*w_PtJet2_fr);
     theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2Gen_W_01"+region, "p_T^{jet2}", Xbins_ptjet2,ptjet2_gen,theWeight*w_PtJet2_fr);
     
     //ETAJET2
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2Gen_W_"+sample+region, "#eta^{jet2}", Xbins_etajet2,etajet2_gen,theWeight*w_EtaJet2_fr);
     theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2Gen_W_01"+region, "#eta^{jet2}", Xbins_etajet2,etajet2_gen,theWeight*w_EtaJet2_fr);
     
   }
   
   //CENTRALMJJ AND DETA
   if(ncentraljets>=2){  
     
     if(decay == "4m"){
       bin_CentralMjj_fr = h_UnfOverMC_fr_CentralMjj_4m->FindBin(mjj_gen_cj);  
       w_CentralMjj_fr =  h_UnfOverMC_fr_CentralMjj_4m->GetBinContent(bin_CentralMjj_fr);
       bin_CentralDeta_fr = h_UnfOverMC_fr_CentralDeta_4m->FindBin(deta_gen_cj);  
       w_CentralDeta_fr =  h_UnfOverMC_fr_CentralDeta_4m->GetBinContent(bin_CentralDeta_fr);
     } 
     if(decay == "4e"){
       bin_CentralMjj_fr = h_UnfOverMC_fr_CentralMjj_4e->FindBin(mjj_gen_cj);  
       w_CentralMjj_fr =  h_UnfOverMC_fr_CentralMjj_4e->GetBinContent(bin_CentralMjj_fr);
       bin_CentralDeta_fr = h_UnfOverMC_fr_CentralDeta_4e->FindBin(deta_gen_cj);  
       w_CentralDeta_fr =  h_UnfOverMC_fr_CentralDeta_4e->GetBinContent(bin_CentralDeta_fr);
     } 
     if(decay == "2e2m"){
       bin_CentralMjj_fr = h_UnfOverMC_fr_CentralMjj_2e2m->FindBin(mjj_gen_cj);  
       w_CentralMjj_fr =  h_UnfOverMC_fr_CentralMjj_2e2m->GetBinContent(bin_CentralMjj_fr);
       bin_CentralDeta_fr = h_UnfOverMC_fr_CentralDeta_2e2m->FindBin(deta_gen_cj);  
       w_CentralDeta_fr =  h_UnfOverMC_fr_CentralDeta_2e2m->GetBinContent(bin_CentralDeta_fr);
     }  
     
     //CENTRALMJJ
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGen_W_"+sample+region, "m_{jj}", Xbins_mjj,mjj_gen_cj,theMCInfo.sampleWeight()*w_CentralMjj_fr);
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjGen_W_01"+region, "m_{jj}", Xbins_mjj, mjj_gen_cj,theMCInfo.sampleWeight()*w_CentralMjj_fr);
     
     //CENTRALDETA
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGen_W_"+sample+region, "#Delta#eta_{jj}", Xbins_deta,deta_gen_cj,theMCInfo.sampleWeight()*w_CentralDeta_fr);
     theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaGen_W_01"+region, "#Delta#eta_{jj}", Xbins_deta,deta_gen_cj,theMCInfo.sampleWeight()*w_CentralDeta_fr);
     
   }
 }
}
void ZZMCWAnalyzer::analyze(){
 
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

void ZZMCWAnalyzer::begin() {
  UnfOverMC = new TFile("macros/test/UnfoldFolder_Mad/Ratio_UnfoldedDataOverGenMC.root");
  UnfOverMC_Pow = new TFile("macros/test/UnfoldFolder_Pow/Ratio_UnfoldedDataOverGenMC.root");
  UnfOverMC_fr = new TFile("macros/test/UnfoldFolder_fr_Mad/Ratio_UnfoldedDataOverGenMC.root");
  UnfOverMC_fr_Pow = new TFile("macros/test/UnfoldFolder_fr_Pow/Ratio_UnfoldedDataOverGenMC.root");

  h_UnfOverMC_Mass_4e = (TH1*) UnfOverMC_Pow->Get("ZZTo4e_Mass_Ratio"); 
  h_UnfOverMC_Mass_4m = (TH1*) UnfOverMC_Pow->Get("ZZTo4m_Mass_Ratio"); 
  h_UnfOverMC_Mass_2e2m = (TH1*) UnfOverMC_Pow->Get("ZZTo2e2m_Mass_Ratio");
  h_UnfOverMC_Jets_4e = (TH1*) UnfOverMC->Get("ZZTo4e_Jets_Ratio"); 
  h_UnfOverMC_Jets_4m = (TH1*) UnfOverMC->Get("ZZTo4m_Jets_Ratio"); 
  h_UnfOverMC_Jets_2e2m = (TH1*) UnfOverMC->Get("ZZTo2e2m_Jets_Ratio");
  h_UnfOverMC_Mjj_4e = (TH1*) UnfOverMC->Get("ZZTo4e_Mjj_Ratio"); 
  h_UnfOverMC_Mjj_4m = (TH1*) UnfOverMC->Get("ZZTo4m_Mjj_Ratio"); 
  h_UnfOverMC_Mjj_2e2m = (TH1*) UnfOverMC->Get("ZZTo2e2m_Mjj_Ratio");
  h_UnfOverMC_Deta_4e = (TH1*) UnfOverMC->Get("ZZTo4e_Deta_Ratio"); 
  h_UnfOverMC_Deta_4m = (TH1*) UnfOverMC->Get("ZZTo4m_Deta_Ratio"); 
  h_UnfOverMC_Deta_2e2m = (TH1*) UnfOverMC->Get("ZZTo2e2m_Deta_Ratio");
  h_UnfOverMC_CentralJets_4e = (TH1*) UnfOverMC->Get("ZZTo4e_CentralJets_Ratio"); 
  h_UnfOverMC_CentralJets_4m = (TH1*) UnfOverMC->Get("ZZTo4m_CentralJets_Ratio"); 
  h_UnfOverMC_CentralJets_2e2m = (TH1*) UnfOverMC->Get("ZZTo2e2m_CentralJets_Ratio");
  h_UnfOverMC_CentralMjj_4e = (TH1*) UnfOverMC->Get("ZZTo4e_CentralMjj_Ratio"); 
  h_UnfOverMC_CentralMjj_4m = (TH1*) UnfOverMC->Get("ZZTo4m_CentralMjj_Ratio"); 
  h_UnfOverMC_CentralMjj_2e2m = (TH1*) UnfOverMC->Get("ZZTo2e2m_CentralMjj_Ratio");
  h_UnfOverMC_CentralDeta_4e = (TH1*) UnfOverMC->Get("ZZTo4e_CentralDeta_Ratio"); 
  h_UnfOverMC_CentralDeta_4m = (TH1*) UnfOverMC->Get("ZZTo4m_CentralDeta_Ratio"); 
  h_UnfOverMC_CentralDeta_2e2m = (TH1*) UnfOverMC->Get("ZZTo2e2m_CentralDeta_Ratio"); 
  h_UnfOverMC_PtJet1_4e = (TH1*) UnfOverMC->Get("ZZTo4e_PtJet1_Ratio"); 
  h_UnfOverMC_PtJet1_4m = (TH1*) UnfOverMC->Get("ZZTo4m_PtJet1_Ratio"); 
  h_UnfOverMC_PtJet1_2e2m = (TH1*) UnfOverMC->Get("ZZTo2e2m_PtJet1_Ratio"); 
  h_UnfOverMC_PtJet2_4e = (TH1*) UnfOverMC->Get("ZZTo4e_PtJet2_Ratio"); 
  h_UnfOverMC_PtJet2_4m = (TH1*) UnfOverMC->Get("ZZTo4m_PtJet2_Ratio"); 
  h_UnfOverMC_PtJet2_2e2m = (TH1*) UnfOverMC->Get("ZZTo2e2m_PtJet2_Ratio");
  h_UnfOverMC_EtaJet1_4e = (TH1*) UnfOverMC->Get("ZZTo4e_EtaJet1_Ratio"); 
  h_UnfOverMC_EtaJet1_4m = (TH1*) UnfOverMC->Get("ZZTo4m_EtaJet1_Ratio"); 
  h_UnfOverMC_EtaJet1_2e2m = (TH1*) UnfOverMC->Get("ZZTo2e2m_EtaJet1_Ratio"); 
  h_UnfOverMC_EtaJet2_4e = (TH1*) UnfOverMC->Get("ZZTo4e_EtaJet2_Ratio"); 
  h_UnfOverMC_EtaJet2_4m = (TH1*) UnfOverMC->Get("ZZTo4m_EtaJet2_Ratio"); 
  h_UnfOverMC_EtaJet2_2e2m = (TH1*) UnfOverMC->Get("ZZTo2e2m_EtaJet2_Ratio");  

  h_UnfOverMC_fr_Mass_4e = (TH1*) UnfOverMC_fr_Pow->Get("ZZTo4e_Mass_Ratio"); 
  h_UnfOverMC_fr_Mass_4m = (TH1*) UnfOverMC_fr_Pow->Get("ZZTo4m_Mass_Ratio"); 
  h_UnfOverMC_fr_Mass_2e2m = (TH1*) UnfOverMC_fr_Pow->Get("ZZTo2e2m_Mass_Ratio");
  h_UnfOverMC_fr_Jets_4e = (TH1*) UnfOverMC_fr->Get("ZZTo4e_Jets_Ratio"); 
  h_UnfOverMC_fr_Jets_4m = (TH1*) UnfOverMC_fr->Get("ZZTo4m_Jets_Ratio"); 
  h_UnfOverMC_fr_Jets_2e2m = (TH1*) UnfOverMC_fr->Get("ZZTo2e2m_Jets_Ratio");
  h_UnfOverMC_fr_Mjj_4e = (TH1*) UnfOverMC_fr->Get("ZZTo4e_Mjj_Ratio"); 
  h_UnfOverMC_fr_Mjj_4m = (TH1*) UnfOverMC_fr->Get("ZZTo4m_Mjj_Ratio"); 
  h_UnfOverMC_fr_Mjj_2e2m = (TH1*) UnfOverMC_fr->Get("ZZTo2e2m_Mjj_Ratio");
  h_UnfOverMC_fr_Deta_4e = (TH1*) UnfOverMC_fr->Get("ZZTo4e_Deta_Ratio"); 
  h_UnfOverMC_fr_Deta_4m = (TH1*) UnfOverMC_fr->Get("ZZTo4m_Deta_Ratio"); 
  h_UnfOverMC_fr_Deta_2e2m = (TH1*) UnfOverMC_fr->Get("ZZTo2e2m_Deta_Ratio");
  h_UnfOverMC_fr_CentralJets_4e = (TH1*) UnfOverMC_fr->Get("ZZTo4e_CentralJets_Ratio"); 
  h_UnfOverMC_fr_CentralJets_4m = (TH1*) UnfOverMC_fr->Get("ZZTo4m_CentralJets_Ratio"); 
  h_UnfOverMC_fr_CentralJets_2e2m = (TH1*) UnfOverMC_fr->Get("ZZTo2e2m_CentralJets_Ratio");
  h_UnfOverMC_fr_CentralMjj_4e = (TH1*) UnfOverMC_fr->Get("ZZTo4e_CentralMjj_Ratio"); 
  h_UnfOverMC_fr_CentralMjj_4m = (TH1*) UnfOverMC_fr->Get("ZZTo4m_CentralMjj_Ratio"); 
  h_UnfOverMC_fr_CentralMjj_2e2m = (TH1*) UnfOverMC_fr->Get("ZZTo2e2m_CentralMjj_Ratio");
  h_UnfOverMC_fr_CentralDeta_4e = (TH1*) UnfOverMC_fr->Get("ZZTo4e_CentralDeta_Ratio"); 
  h_UnfOverMC_fr_CentralDeta_4m = (TH1*) UnfOverMC_fr->Get("ZZTo4m_CentralDeta_Ratio"); 
  h_UnfOverMC_fr_CentralDeta_2e2m = (TH1*) UnfOverMC_fr->Get("ZZTo2e2m_CentralDeta_Ratio"); 
  h_UnfOverMC_fr_PtJet1_4e = (TH1*) UnfOverMC_fr->Get("ZZTo4e_PtJet1_Ratio"); 
  h_UnfOverMC_fr_PtJet1_4m = (TH1*) UnfOverMC_fr->Get("ZZTo4m_PtJet1_Ratio"); 
  h_UnfOverMC_fr_PtJet1_2e2m = (TH1*) UnfOverMC_fr->Get("ZZTo2e2m_PtJet1_Ratio"); 
  h_UnfOverMC_fr_PtJet2_4e = (TH1*) UnfOverMC_fr->Get("ZZTo4e_PtJet2_Ratio"); 
  h_UnfOverMC_fr_PtJet2_4m = (TH1*) UnfOverMC_fr->Get("ZZTo4m_PtJet2_Ratio"); 
  h_UnfOverMC_fr_PtJet2_2e2m = (TH1*) UnfOverMC_fr->Get("ZZTo2e2m_PtJet2_Ratio");
  h_UnfOverMC_fr_EtaJet1_4e = (TH1*) UnfOverMC_fr->Get("ZZTo4e_EtaJet1_Ratio"); 
  h_UnfOverMC_fr_EtaJet1_4m = (TH1*) UnfOverMC_fr->Get("ZZTo4m_EtaJet1_Ratio"); 
  h_UnfOverMC_fr_EtaJet1_2e2m = (TH1*) UnfOverMC_fr->Get("ZZTo2e2m_EtaJet1_Ratio"); 
  h_UnfOverMC_fr_EtaJet2_4e = (TH1*) UnfOverMC_fr->Get("ZZTo4e_EtaJet2_Ratio"); 
  h_UnfOverMC_fr_EtaJet2_4m = (TH1*) UnfOverMC_fr->Get("ZZTo4m_EtaJet2_Ratio"); 
  h_UnfOverMC_fr_EtaJet2_2e2m = (TH1*) UnfOverMC_fr->Get("ZZTo2e2m_EtaJet2_Ratio");  

  nentries =  tree()->GetEntries();
  PreCounter = 0;
  Xbins += 100,200,250,300,350,400,500,600,800; 
  Xbins_deta += 0,2.4,4.7;
  Xbins_mjj += 0.,200,800; 
  Xbins_ptjet1 += 30,50,100,200,300,500;
  Xbins_ptjet2 += 30,100,200,500;
  Xbins_etajet1 += 0,1.5,3,4.7;
  Xbins_etajet2 +=0,1.5,3,4.7;
  
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

  bin_Jets=0;
  bin_Mass=0;
  bin_Mjj=0;
  bin_Deta=0;
  bin_CentralJets=0;
  bin_CentralMjj=0;
  bin_CentralDeta=0; 
  bin_PtJet1 = 0;
  bin_PtJet2 = 0;
  bin_EtaJet1 = 0;
  bin_EtaJet2 = 0;
  w_Jets=0;
  w_Mass=0;
  w_Mjj=0;
  w_Deta=0;
  w_CentralJets=0;
  w_CentralMjj=0;
  w_CentralDeta=0;
  w_PtJet1 = 0;
  w_PtJet2 = 0;
  w_EtaJet1 = 0;
  w_EtaJet2 = 0;
  
  bin_Mass_fr = 0;
  bin_Mjj_fr = 0;
  bin_Deta_fr = 0;
  bin_CentralJets_fr = 0;
  bin_CentralMjj_fr = 0;
  bin_CentralDeta_fr = 0;
  bin_PtJet1_fr = 0;
  bin_PtJet2_fr = 0;
  bin_EtaJet1_fr = 0;
  bin_EtaJet2_fr = 0;
  w_Jets_fr = 0;
  w_Mass_fr = 0;
  w_Mjj_fr = 0;
  w_Deta_fr = 0;
  w_CentralJets_fr = 0;
  w_CentralMjj_fr = 0;
  w_CentralDeta_fr = 0;
  w_PtJet1_fr = 0;
  w_PtJet2_fr = 0;
  w_EtaJet1_fr = 0;
  w_EtaJet2_fr = 0;
}

void ZZMCWAnalyzer::end( TFile &) {
  cout <<"Tree Entries"<<nentries<< endl;
 
}  
