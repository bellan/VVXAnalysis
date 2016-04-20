#include "VVXAnalysis/TreeAnalysis/interface/ZZRecoWAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"

#include "TRandom.h"
#include "TTree.h"
#include <boost/assign/std/vector.hpp> 
#include <boost/assert.hpp> 
using namespace std;
using namespace boost::assign; // bring 'operator+=()' into scope

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using std::cout;
using std::endl;

using namespace phys;

void ZZRecoWAnalyzer::ZZplots(int id, int e){

  if(ZZ->id() != id && id != -1) return; // -1 here means generic 4l final state
  
  std::string decay  = "4l";
  std::string sample;
 
  if(id == 52) decay = "4m";
  else if (id == 48) decay = "2e2m";
  else if (id == 44) decay = "4e";
    
  if(e < nentries/2){sample = "0";}
  else {sample = "1";}
  
  Int_t njets = jets->size();
  if (njets>3) njets=3;
  Int_t ncentraljets = centralJets->size();
  if (ncentraljets>3) ncentraljets=3;
  Int_t npjets = pjets->size();
  if (npjets>3) npjets=3;
 
  //1D Reco Mass Distributions
  m4L = ZZ->mass();  
  if(m4L > 800) m4L = 799;

  //******************************************************************RECO JETS****************************************************************************************
  float centralDeta = 0;
  float centralMjj = 0;
  float centralPtJet1 = 0;
  float centralPtJet2 = 0;
  float centralEtaJet1 = 0;
  float centralEtaJet2 = 0;
  // cout << "1" << endl; 
  if(nCentralJERjets>=1){
    centralPtJet1 = CentralJER_jetPt->at(0);
    if (centralPtJet1>=500) centralPtJet1 = 499; 
    centralEtaJet1 = fabs(CentralJER_jets->at(0).eta());
    if (centralEtaJet1>=6) centralEtaJet1 = 5;
  }
  
 
  //1D Reco DeltaEta and mJJ Distributions - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)
   if(nCentralJERjets>=2){
     
     centralDeta = fabs(CentralJER_jets->at(0).eta() - CentralJER_jets->at(1).eta());
     centralMjj =  (CentralJER_jets->at(0).p4() + CentralJER_jets->at(1).p4()).M();
     centralPtJet2 = CentralJER_jetPt->at(1);
     centralEtaJet2 = fabs(CentralJER_jets->at(1).eta());
     
     if (centralDeta>=4.7) centralDeta = 4.6;
     if (centralMjj>=800) centralMjj = 799;
     if (centralPtJet2>=500) centralPtJet2 = 499; 
     if (centralEtaJet2>=6) centralEtaJet2 = 5;
   }
  
   //******************************************************************RECO CENTRAL JETS****************************************************************************************

  float centralDeta_cj = 0;
  float centralMjj_cj = 0;
 
  //1D Reco DeltaEta and mJJ Distributions - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)
  if(nCentralJERcentraljets>=2){
  
    centralDeta_cj = fabs(CentralJER_centraljets->at(0).eta() - CentralJER_centraljets->at(1).eta());
    centralMjj_cj =  (CentralJER_centraljets->at(0).p4() + CentralJER_centraljets->at(1).p4()).M();
    
    if (centralDeta_cj>=4.7) centralDeta_cj = 4.6;
    if (centralMjj_cj>=800) centralMjj_cj = 799;
    
  }
 
  //if MC gen (for response matrices only)
  if (genCategory !=-1){
    if(topology.test(0)){
      
      ngenjets =  genJets->size(); 
      if (ngenjets>3) ngenjets=3;
      
      ngencentraljets =  centralGenJets->size(); 
      if (ngencentraljets>3) ngencentraljets=3;

      m4L_gen = sqrt((genVBParticles->at(0).p4()+genVBParticles->at(1).p4())*(genVBParticles->at(0).p4()+genVBParticles->at(1).p4()));
      if (m4L_gen>=800) m4L_gen = 799;
      
      //****************************************************************************RECO&GEN****************************************************************************************
      //Matrix and reco distribution weighted for the ratio between the unfolded data and the generator level information in order to 
      //compute the relative systematic uncertainty
      //distributions to evaluate data/MC systematic uncertainty
      
      if(decay == "4m"){
 	bin_Jets = h_UnfOverMC_Jets_4m->FindBin(ngenjets);  
 	w_Jets =  h_UnfOverMC_Jets_4m->GetBinContent(bin_Jets);
	bin_Mass = h_UnfOverMC_Mass_4m->FindBin(m4L_gen);  
 	w_Mass =  h_UnfOverMC_Mass_4m->GetBinContent(bin_Mass);
	bin_CentralJets = h_UnfOverMC_CentralJets_4m->FindBin(ngencentraljets);  
 	w_CentralJets =  h_UnfOverMC_CentralJets_4m->GetBinContent(bin_CentralJets);
      }  
      if(decay == "4e"){
 	bin_Jets = h_UnfOverMC_Jets_4e->FindBin(ngenjets);  
 	w_Jets =  h_UnfOverMC_Jets_4e->GetBinContent(bin_Jets);
	bin_Mass = h_UnfOverMC_Mass_4e->FindBin(m4L_gen);  
 	w_Mass =  h_UnfOverMC_Mass_4e->GetBinContent(bin_Mass);
	bin_CentralJets = h_UnfOverMC_CentralJets_4e->FindBin(ngencentraljets);  
 	w_CentralJets =  h_UnfOverMC_CentralJets_4e->GetBinContent(bin_CentralJets);
      }  
      if(decay == "2e2m"){
 	bin_Jets = h_UnfOverMC_Jets_2e2m->FindBin(ngenjets);  
 	w_Jets =  h_UnfOverMC_Jets_2e2m->GetBinContent(bin_Jets);
	bin_Mass = h_UnfOverMC_Mass_2e2m->FindBin(m4L_gen);  
 	w_Mass =  h_UnfOverMC_Mass_2e2m->GetBinContent(bin_Mass);
	bin_CentralJets = h_UnfOverMC_CentralJets_2e2m->FindBin(ngencentraljets);  
 	w_CentralJets =  h_UnfOverMC_CentralJets_2e2m->GetBinContent(bin_CentralJets);
      }  
      
      //JETS
      theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_W_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERjets , theWeight*w_Jets);
      theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_W_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4, nCentralJERjets, theWeight*w_Jets);
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERCentralSmear_W_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_Jets); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERCentralSmear_W_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_Jets);
      
      //MASS
      theHistograms.fill(std::string("ZZTo")+decay+"_Mass_W_01", std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L,theWeight*w_Mass);
      theHistograms.fill(std::string("ZZTo")+decay+"_Mass_W_"+sample, std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L,theWeight*w_Mass);
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_W_"+sample, std::string("Response matrix invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, Xbins, m4L ,m4L_gen , theWeight*w_Mass); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_W_01", std::string("Response matrix invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, Xbins, m4L ,m4L_gen , theWeight*w_Mass);
      
      //CENTRALJETS
       theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERCentralSmear_W_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERcentraljets , theWeight*w_CentralJets);
      theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERCentralSmear_W_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4, nCentralJERcentraljets, theWeight*w_CentralJets);
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERCentralSmear_W_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_CentralJets); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERCentralSmear_W_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_CentralJets);     

      //      cout << ngencentraljets <<" "<< nCentralJERcentraljets <<" " <<w_CentralJets <<endl;

      //PTJET1 and ETAJET1
      if(ngenjets>=1){ 

	ptjet1_gen = genJets->at(0).pt();
	if(ptjet1_gen >=500) ptjet1_gen =499;
	etajet1_gen = fabs(genJets->at(0).eta());
	if(etajet1_gen >=6) etajet1_gen =5;

	if(nCentralJERjets>=1){
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
	  theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JERCentralSmear_W_"+sample, "p_T^{jet1}", Xbins_ptjet1,centralPtJet1,theWeight*w_PtJet1);
	  theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JERCentralSmear_W_01", "p_T^{jet1}", Xbins_ptjet1, centralPtJet1,theWeight*w_PtJet1);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JERCentralSmear_W_"+sample, std::string("Response matrix p_T^{jet1} of ZZ_{1}#rightarrow ")+decay, Xbins_ptjet1,Xbins_ptjet1, centralPtJet1,ptjet1_gen,theWeight*w_PtJet1);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JERCentralSmear_W_01", std::string("Response matrix  p_T^{jet1} of ZZ_{1}#rightarrow ")+decay,Xbins_ptjet1,Xbins_ptjet1, centralPtJet1,ptjet1_gen,theWeight*w_PtJet1);
	 
	  //ETAJET1
	  theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JERCentralSmear_W_"+sample, "#eta^{jet1}", Xbins_etajet1,centralEtaJet1,theWeight*w_EtaJet1);
	  theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JERCentralSmear_W_01", "#eta^{jet1}", Xbins_etajet1, centralEtaJet1,theWeight*w_EtaJet1);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JERCentralSmear_W_"+sample, std::string("Response matrix #eta^{jet1} of ZZ_{1}#rightarrow ")+decay, Xbins_etajet1,Xbins_etajet1, centralEtaJet1,etajet1_gen,theWeight*w_EtaJet1);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JERCentralSmear_W_01", std::string("Response matrix  #eta^{jet1} of ZZ_{1}#rightarrow ")+decay,Xbins_etajet1,Xbins_etajet1, centralEtaJet1,etajet1_gen,theWeight*w_EtaJet1);
	}
      }

      //DETA AND MJJ
      if(ngenjets>=2){  
 	deta_gen = fabs(genJets->at(0).eta() - genJets->at(1).eta());
 	mjj_gen =  (genJets->at(0).p4() + genJets->at(1).p4()).M();
	ptjet2_gen = genJets->at(1).pt(); 
	etajet2_gen = fabs(genJets->at(1).eta());

 	if (deta_gen>=4.7) deta_gen = 4.6;
 	if (mjj_gen>=800) mjj_gen = 799;
	if(ptjet2_gen >=500) ptjet2_gen =499;
	if(etajet2_gen >=6) etajet2_gen =5;
	
	if(nCentralJERjets>=2){
	  
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
	  theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERCentralSmear_W_"+sample, "m_{jj}", Xbins_mjj,centralMjj,theWeight*w_Mjj);
	  theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERCentralSmear_W_01", "m_{jj}", Xbins_mjj, centralMjj,theWeight*w_Mjj);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERCentralSmear_W_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_mjj,Xbins_mjj, centralMjj,mjj_gen,theWeight*w_Mjj);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERCentralSmear_W_01", std::string("Response matrix  m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,Xbins_mjj, centralMjj,mjj_gen,theWeight*w_Mjj);
	 
	  //DETA
	  theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERCentralSmear_W_"+sample, "#Delta#eta_{jj}", Xbins_deta,centralDeta,theWeight*w_Deta);
	  theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERCentralSmear_W_01", "#Delta#eta_{jj}", Xbins_deta, centralDeta,theWeight*w_Deta);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERCentralSmear_W_"+sample, std::string("Response matrix #Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_deta,Xbins_deta, centralDeta,deta_gen,theWeight*w_Deta);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERCentralSmear_W_01", std::string("Response matrix  #Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,Xbins_deta, centralDeta,deta_gen,theWeight*w_Deta);

	  //PTJET2
	  theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JERCentralSmear_W_"+sample, "p_{T}^{jet2}", Xbins_ptjet2,centralPtJet2,theWeight*w_PtJet2);
	  theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JERCentralSmear_W_01", "p_{T}^{jet2}", Xbins_ptjet2, centralPtJet2,theWeight*w_PtJet2);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JERCentralSmear_W_"+sample, std::string("Response matrix p_{T}^{jet2} of ZZ_{1}#rightarrow ")+decay, Xbins_ptjet2,Xbins_ptjet2, centralPtJet2,ptjet2_gen,theWeight*w_PtJet2);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JERCentralSmear_W_01", std::string("Response matrix p_{T}^{jet2} of ZZ_{1}#rightarrow ")+decay,Xbins_ptjet2,Xbins_ptjet2, centralPtJet2,ptjet2_gen,theWeight*w_PtJet2);
	 
	  //ETAJET2
	  theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JERCentralSmear_W_"+sample, "#eta^{jet2}", Xbins_etajet2,centralEtaJet2,theWeight*w_EtaJet2);
	  theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JERCentralSmear_W_01", "#eta^{jet2}", Xbins_etajet2, centralEtaJet2,theWeight*w_EtaJet2);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JERCentralSmear_W_"+sample, std::string("Response matrix #eta^{jet2} of ZZ_{1}#rightarrow ")+decay, Xbins_etajet2,Xbins_etajet2, centralEtaJet2,etajet2_gen,theWeight*w_EtaJet2);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JERCentralSmear_W_01", std::string("Response matrix  #eta^{jet2} of ZZ_{1}#rightarrow ")+decay,Xbins_etajet2,Xbins_etajet2, centralEtaJet2,etajet2_gen,theWeight*w_EtaJet2);
	} 
      }
                
      //CENTRALMJJ AND CENTRALDETA
      if(ngencentraljets>=2){  

 	deta_gen_cj = fabs(centralGenJets->at(0).eta() - centralGenJets->at(1).eta());
 	mjj_gen_cj =  (centralGenJets->at(0).p4() + centralGenJets->at(1).p4()).M();
	
 	if (deta_gen_cj>=4.7) deta_gen_cj = 4.6;
 	if (mjj_gen_cj>=800) mjj_gen_cj = 799;
	
	if(nCentralJERcentraljets>=2){

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
	  theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERCentralSmear_W_"+sample, "m_{jj}", Xbins_mjj,centralMjj_cj,theWeight*w_CentralMjj);
	  theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERCentralSmear_W_01", "m_{jj}", Xbins_mjj, centralMjj_cj,theWeight*w_CentralMjj);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERCentralSmear_W_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_mjj,Xbins_mjj, centralMjj_cj,mjj_gen_cj,theWeight*w_CentralMjj);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERCentralSmear_W_01", std::string("Response matrix  m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,Xbins_mjj, centralMjj_cj,mjj_gen_cj,theWeight*w_CentralMjj);
	  
	  //CENTRALDETA
	  theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERCentralSmear_W_"+sample, "#Delta#eta_{jj}", Xbins_deta,centralDeta_cj,theWeight*w_CentralDeta);
	  theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERCentralSmear_W_01", "#Delta#eta_{jj}", Xbins_deta, centralDeta_cj,theWeight*w_CentralDeta);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERCentralSmear_W_"+sample, std::string("Response matrix #Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_deta,Xbins_deta, centralDeta_cj,deta_gen_cj,theWeight*w_CentralDeta);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERCentralSmear_W_01", std::string("Response matrix  #Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,Xbins_deta, centralDeta_cj,deta_gen_cj,theWeight*w_CentralDeta);
	}
      }
    } //end if(topology.test(0))
    
      // if the event is reconstructed but not generated as signal, put the w_Mass=1
    else{
      theHistograms.fill(std::string("ZZTo")+decay+"_Mass_W_01", std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L,theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_Mass_W_"+sample, std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L,theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_W_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERjets , theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_W_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4, nCentralJERjets, theWeight); 
      theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERCentralSmear_W_"+sample, "m_{jj}", Xbins_mjj,centralMjj,theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERCentralSmear_W_01", "m_{jj}", Xbins_mjj, centralMjj,theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERCentralSmear_W_"+sample, "#Delta#eta_{jj}", Xbins_deta,centralDeta,theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERCentralSmear_W_01", "#Delta#eta_{jj}", Xbins_deta, centralDeta,theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERCentralSmear_W_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERcentraljets , theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERCentralSmear_W_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4, nCentralJERcentraljets, theWeight); 
      theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERCentralSmear_W_"+sample, "m_{jj}", Xbins_mjj,centralMjj_cj,theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERCentralSmear_W_01", "m_{jj}", Xbins_mjj, centralMjj_cj,theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERCentralSmear_W_"+sample, "#Delta#eta_{jj}", Xbins_deta,centralDeta_cj,theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERCentralSmear_W_01", "#Delta#eta_{jj}", Xbins_deta, centralDeta_cj,theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JERCentralSmear_W_"+sample, "p_{T}^{jet1}", Xbins_ptjet1,centralPtJet1,theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JERCentralSmear_W_01", "p_{T}^{jet1}", Xbins_ptjet1, centralPtJet1,theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JERCentralSmear_W_"+sample, "p_{T}^{jet2}", Xbins_ptjet2,centralPtJet2,theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JERCentralSmear_W_01", "p_{T}^{jet2}", Xbins_ptjet2, centralPtJet2,theWeight); 
      theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JERCentralSmear_W_"+sample, "#eta^{jet1}", Xbins_etajet1,centralEtaJet1,theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JERCentralSmear_W_01", "#eta^{jet1}", Xbins_etajet1, centralEtaJet1,theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JERCentralSmear_W_"+sample, "#eta^{jet2}", Xbins_etajet2,centralEtaJet2,theWeight);
      theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JERCentralSmear_W_01", "#eta^{jet2}", Xbins_etajet2, centralEtaJet2,theWeight);
    }

    /////////////////////////////////////////////////IN TIGHTER FIDUCIAL REGION////////////////////////////////////////////////
  
    zz::SignalTopology zzSignalTopology = zz::getSignalTopology(*genParticles, *genJets);
    string region; 
    if(zz::inTightFiducialRegion(zzSignalTopology)){ 
      region = "_fr";
      if(topology.test(0)){
	if(decay == "4m"){
	  bin_Jets_fr = h_UnfOverMC_fr_Jets_4m->FindBin(ngenjets);
	  w_Jets_fr =  h_UnfOverMC_fr_Jets_4m->GetBinContent(bin_Jets_fr);
	  bin_Mass_fr = h_UnfOverMC_fr_Mass_4m->FindBin(m4L_gen);  
	  w_Mass_fr =  h_UnfOverMC_fr_Mass_4m->GetBinContent(bin_Mass_fr); 
	  bin_CentralJets_fr = h_UnfOverMC_fr_CentralJets_4m->FindBin(ngencentraljets); 
	  w_CentralJets_fr =  h_UnfOverMC_fr_CentralJets_4m->GetBinContent(bin_CentralJets_fr);
	}  
	if(decay == "4e"){
	  bin_Jets_fr = h_UnfOverMC_fr_Jets_4e->FindBin(ngenjets); 
	  w_Jets_fr =  h_UnfOverMC_fr_Jets_4e->GetBinContent(bin_Jets_fr);
	  bin_Mass_fr = h_UnfOverMC_fr_Mass_4e->FindBin(m4L_gen); 
	  w_Mass_fr =  h_UnfOverMC_fr_Mass_4e->GetBinContent(bin_Mass_fr);
	  bin_CentralJets_fr = h_UnfOverMC_fr_CentralJets_4e->FindBin(ngencentraljets);
	  w_CentralJets_fr =  h_UnfOverMC_fr_CentralJets_4e->GetBinContent(bin_CentralJets_fr);
	}  
	if(decay == "2e2m"){
	  bin_Jets_fr = h_UnfOverMC_fr_Jets_2e2m->FindBin(ngenjets);
	  w_Jets_fr =  h_UnfOverMC_fr_Jets_2e2m->GetBinContent(bin_Jets_fr);
	  bin_Mass_fr = h_UnfOverMC_fr_Mass_2e2m->FindBin(m4L_gen);  
	  w_Mass_fr =  h_UnfOverMC_fr_Mass_2e2m->GetBinContent(bin_Mass_fr);
	  bin_CentralJets_fr = h_UnfOverMC_fr_CentralJets_2e2m->FindBin(ngencentraljets);
	  w_CentralJets_fr =  h_UnfOverMC_fr_CentralJets_2e2m->GetBinContent(bin_CentralJets_fr);
	}  

	//JETS
	theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_W_"+sample+region, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERjets , theWeight*w_Jets_fr);
	theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_W_01"+region, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4, nCentralJERjets, theWeight*w_Jets_fr);
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERCentralSmear_W_"+sample+region, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_Jets_fr); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERCentralSmear_W_01"+region, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_Jets_fr);

	//MASS
	theHistograms.fill(std::string("ZZTo")+decay+"_Mass_W_01"+region, std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L,theWeight*w_Mass_fr);
	theHistograms.fill(std::string("ZZTo")+decay+"_Mass_W_"+sample+region, std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L,theWeight*w_Mass_fr);
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_W_"+sample+region, std::string("Response matrix invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, Xbins, m4L ,m4L_gen , theWeight*w_Mass_fr); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_W_01"+region, std::string("Response matrix invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, Xbins, m4L ,m4L_gen , theWeight*w_Mass_fr);
	
	//CENTRALJETS
	theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERCentralSmear_W_"+sample+region, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERcentraljets , theWeight*w_CentralJets_fr);
	theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERCentralSmear_W_01"+region, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4, nCentralJERcentraljets, theWeight*w_CentralJets_fr);
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERCentralSmear_W_"+sample+region, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_CentralJets_fr); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERCentralSmear_W_01"+region, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_CentralJets_fr);     
	
	//      cout << ngencentraljets <<" "<< nCentralJERcentraljets <<" " <<w_CentralJets <<endl;

	//PTJET1 and ETAJET1
	if(ngenjets>=1){
	  if(nCentralJERjets>=1){ 
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
	    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JERCentralSmear_W_"+sample+region, "p_T^{jet1}", Xbins_ptjet1,centralPtJet1,theWeight*w_PtJet1_fr);
	    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JERCentralSmear_W_01"+region, "p_T^{jet1}", Xbins_ptjet1, centralPtJet1,theWeight*w_PtJet1_fr);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JERCentralSmear_W_"+sample+region, std::string("Response matrix p_T^{jet1} of ZZ_{1}#rightarrow ")+decay, Xbins_ptjet1,Xbins_ptjet1, centralPtJet1,ptjet1_gen,theWeight*w_PtJet1_fr);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JERCentralSmear_W_01"+region, std::string("Response matrix  p_T^{jet1} of ZZ_{1}#rightarrow ")+decay,Xbins_ptjet1,Xbins_ptjet1, centralPtJet1,ptjet1_gen,theWeight*w_PtJet1_fr);
	   
	    //ETAJET1
	    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JERCentralSmear_W_"+sample+region, "#eta^{jet1}", Xbins_etajet1,centralEtaJet1,theWeight*w_EtaJet1_fr);
	    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JERCentralSmear_W_01"+region, "#eta^{jet1}", Xbins_etajet1, centralEtaJet1,theWeight*w_EtaJet1_fr);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JERCentralSmear_W_"+sample+region, std::string("Response matrix #eta^{jet1} of ZZ_{1}#rightarrow ")+decay, Xbins_etajet1,Xbins_etajet1, centralEtaJet1,etajet1_gen,theWeight*w_EtaJet1_fr);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JERCentralSmear_W_01"+region, std::string("Response matrix  #eta^{jet1} of ZZ_{1}#rightarrow ")+decay,Xbins_etajet1,Xbins_etajet1, centralEtaJet1,etajet1_gen,theWeight*w_EtaJet1_fr);
	  }
	}
	
	//DETA AND MJJ
	if(ngenjets>=2){  
	  if(nCentralJERjets>=2){
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
	    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERCentralSmear_W_"+sample+region, "m_{jj}", Xbins_mjj,centralMjj,theWeight*w_Mjj_fr);
	    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERCentralSmear_W_01"+region, "m_{jj}", Xbins_mjj, centralMjj,theWeight*w_Mjj_fr);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERCentralSmear_W_"+sample+region, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_mjj,Xbins_mjj, centralMjj,mjj_gen,theWeight*w_Mjj_fr);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERCentralSmear_W_01"+region, std::string("Response matrix  m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,Xbins_mjj, centralMjj,mjj_gen,theWeight*w_Mjj_fr);
	    
	    //DETA
	    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERCentralSmear_W_"+sample+region, "#Delta#eta_{jj}", Xbins_deta,centralDeta,theWeight*w_Deta_fr);
	    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERCentralSmear_W_01"+region, "#Delta#eta_{jj}", Xbins_deta, centralDeta,theWeight*w_Deta_fr);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERCentralSmear_W_"+sample+region, std::string("Response matrix #Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_deta,Xbins_deta, centralDeta,deta_gen,theWeight*w_Deta_fr);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERCentralSmear_W_01"+region, std::string("Response matrix  #Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,Xbins_deta, centralDeta,deta_gen,theWeight*w_Deta_fr);
	    
	    //PTJET2
	    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JERCentralSmear_W_"+sample+region, "p_{T}^{jet2}", Xbins_ptjet2,centralPtJet2,theWeight*w_PtJet2_fr);
	    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JERCentralSmear_W_01"+region, "p_{T}^{jet2}", Xbins_ptjet2, centralPtJet2,theWeight*w_PtJet2_fr);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JERCentralSmear_W_"+sample+region, std::string("Response matrix p_{T}^{jet2} of ZZ_{1}#rightarrow ")+decay, Xbins_ptjet2,Xbins_ptjet2, centralPtJet2,ptjet2_gen,theWeight*w_PtJet2_fr);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JERCentralSmear_W_01"+region, std::string("Response matrix p_{T}^{jet2} of ZZ_{1}#rightarrow ")+decay,Xbins_ptjet2,Xbins_ptjet2, centralPtJet2,ptjet2_gen,theWeight*w_PtJet2_fr);
	    
	    //ETAJET2
	    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JERCentralSmear_W_"+sample+region, "#eta^{jet2}", Xbins_etajet2,centralEtaJet2,theWeight*w_EtaJet2_fr);
	    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JERCentralSmear_W_01"+region, "#eta^{jet2}", Xbins_etajet2, centralEtaJet2,theWeight*w_EtaJet2_fr);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JERCentralSmear_W_"+sample+region, std::string("Response matrix #eta^{jet2} of ZZ_{1}#rightarrow ")+decay, Xbins_etajet2,Xbins_etajet2, centralEtaJet2,etajet2_gen,theWeight*w_EtaJet2_fr);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JERCentralSmear_W_01"+region, std::string("Response matrix  #eta^{jet2} of ZZ_{1}#rightarrow ")+decay,Xbins_etajet2,Xbins_etajet2, centralEtaJet2,etajet2_gen,theWeight*w_EtaJet2_fr);
	  } 
	}
	
	//CENTRALMJJ AND CENTRALDETA
	if(ngencentraljets>=2){  
	  if(nCentralJERcentraljets>=2){
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
	    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERCentralSmear_W_"+sample+region, "m_{jj}", Xbins_mjj,centralMjj_cj,theWeight*w_CentralMjj_fr);
	    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERCentralSmear_W_01"+region, "m_{jj}", Xbins_mjj, centralMjj_cj,theWeight*w_CentralMjj_fr);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERCentralSmear_W_"+sample+region, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_mjj,Xbins_mjj, centralMjj_cj,mjj_gen_cj,theWeight*w_CentralMjj_fr);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERCentralSmear_W_01"+region, std::string("Response matrix  m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,Xbins_mjj, centralMjj_cj,mjj_gen_cj,theWeight*w_CentralMjj_fr);
	    
	    //CENTRALDETA
	    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERCentralSmear_W_"+sample+region, "#Delta#eta_{jj}", Xbins_deta,centralDeta_cj,theWeight*w_CentralDeta_fr);
	    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERCentralSmear_W_01"+region, "#Delta#eta_{jj}", Xbins_deta, centralDeta_cj,theWeight*w_CentralDeta_fr);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERCentralSmear_W_"+sample+region, std::string("Response matrix #Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_deta,Xbins_deta, centralDeta_cj,deta_gen_cj,theWeight*w_CentralDeta_fr);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERCentralSmear_W_01"+region, std::string("Response matrix  #Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,Xbins_deta, centralDeta_cj,deta_gen_cj,theWeight*w_CentralDeta_fr);
	  }
	}
      } //end if(topology.test(0))
   
      // if the event is reconstructed but not generated as signal, put the w_Mass=1
      else{
	theHistograms.fill(std::string("ZZTo")+decay+"_Mass_W_01"+region, std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L,theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_Mass_W_"+sample+region, std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L,theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_W_"+sample+region, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERjets , theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_W_01"+region, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4, nCentralJERjets, theWeight); 
	theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERCentralSmear_W_"+sample+region, "m_{jj}", Xbins_mjj,centralMjj,theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERCentralSmear_W_01"+region, "m_{jj}", Xbins_mjj, centralMjj,theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERCentralSmear_W_"+sample+region, "#Delta#eta_{jj}", Xbins_deta,centralDeta,theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERCentralSmear_W_01"+region, "#Delta#eta_{jj}", Xbins_deta, centralDeta,theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERCentralSmear_W_"+sample+region, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERcentraljets , theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERCentralSmear_W_01"+region, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4, nCentralJERcentraljets, theWeight); 
	theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERCentralSmear_W_"+sample+region, "m_{jj}", Xbins_mjj,centralMjj_cj,theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERCentralSmear_W_01"+region, "m_{jj}", Xbins_mjj, centralMjj_cj,theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERCentralSmear_W_"+sample+region, "#Delta#eta_{jj}", Xbins_deta,centralDeta_cj,theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERCentralSmear_W_01"+region, "#Delta#eta_{jj}", Xbins_deta, centralDeta_cj,theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JERCentralSmear_W_"+sample+region, "p_{T}^{jet1}", Xbins_ptjet1,centralPtJet1,theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JERCentralSmear_W_01"+region, "p_{T}^{jet1}", Xbins_ptjet1, centralPtJet1,theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JERCentralSmear_W_"+sample+region, "p_{T}^{jet2}", Xbins_ptjet2,centralPtJet2,theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JERCentralSmear_W_01"+region, "p_{T}^{jet2}", Xbins_ptjet2, centralPtJet2,theWeight); 
	theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JERCentralSmear_W_"+sample+region, "#eta^{jet1}", Xbins_etajet1,centralEtaJet1,theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JERCentralSmear_W_01"+region, "#eta^{jet1}", Xbins_etajet1, centralEtaJet1,theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JERCentralSmear_W_"+sample+region, "#eta^{jet2}", Xbins_etajet2,centralEtaJet2,theWeight);
	theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JERCentralSmear_W_01"+region, "#eta^{jet2}", Xbins_etajet2, centralEtaJet2,theWeight);
      }
      
    }//end if(zz::inTightFiducialRegion(zzSignalTopology)) 
    
  }//end if(genCategory != -1)    
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if(region_ == phys::CR3P1F || region_ == phys::CR2P2F) {
    theHistograms.fill(std::string("ZZTo")+decay+"_Mass"+"_FRVar", std::string("Var From FR Invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, m4L,ZZ->fakeRateSFVar());
    theHistograms.fill(std::string("ZZTo")+decay+"_Jets"+"_FRVar", std::string("Var From FR Invariant mass of ZZ_{1}#rightarrow ")+decay ,4,0,4,njets,ZZ->fakeRateSFVar());
  }
}// end ZZplots()

//Smearing function
double ZZRecoWAnalyzer::JER_PtSmear(double pt, double width)
{
  double ptsmear= gRandom->Gaus(pt,width);    
  return ptsmear;
}

void ZZRecoWAnalyzer::analyze(){
  
  e++;
  
//Smearing jet pt (Central for the standard analysis) without requiring it is signal (1D distributions made of ALL reco events)
    CentralJER_jets->clear();
    CentralJER_centraljets->clear();
    CentralJER_jetPt->clear();
    
    //Loop on all reco jets
    foreach(const phys::Jet &jet, *pjets){
      
      double jetPt = 0;
      jetPt = jet.pt();
      
      //JER correction (applied only on MC reco). Up and down velues used to assess systematic uncertainty on jet energy resolution
      double newJetPtJER =0; 
      double width = 0;
           
      width = jet.jer_width(phys::Jet::central);
           
      newJetPtJER = JER_PtSmear(jetPt, width);
            
      if(newJetPtJER > 30) {
	CentralJER_jets->push_back(jet);
	CentralJER_jetPt->push_back(newJetPtJER);
      }
      if(newJetPtJER > 30 && fabs(jet.eta()<2.4)) CentralJER_centraljets->push_back(jet);
       
    }
    
    nCentralJERjets = CentralJER_jets->size();
    nCentralJERcentraljets = CentralJER_centraljets->size();
        
    if (nCentralJERjets>3) nCentralJERjets=3;
    if (nCentralJERcentraljets>3) nCentralJERcentraljets=3;
    
    stable_sort(CentralJER_jets->begin(), CentralJER_jets->end(), PtComparator());
    stable_sort(CentralJER_centraljets->begin(), CentralJER_centraljets->end(), PtComparator());
    stable_sort(CentralJER_jetPt->begin(), CentralJER_jetPt->end(),std::greater<float>());
    
    // ZZplots();   // ZZ --> 4l 
    ZZplots(52,e); // ZZ --> 4m
    ZZplots(48,e); // ZZ --> 2e2m
    ZZplots(44,e); // ZZ --> 4e
    
  }
  
void ZZRecoWAnalyzer::begin() {
  
  CentralJER_jets  = new std::vector<phys::Jet>();
  CentralJER_centraljets  = new std::vector<phys::Jet>();
  CentralJER_jetPt  = new std::vector<double>();
  
  nCentralJERjets=0;
  nCentralJERcentraljets=0;
  
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

  nentries = tree()->GetEntries("ZZCand.passFullSel_"); 
 
  Xbins += 100,200,250,300,350,400,500,600,800;
  Xbins_deta += 0,2.4,4.7;
  Xbins_mjj += 0.,200,800;
  Xbins_ptjet1 += 30,50,100,200,300,500;
  Xbins_ptjet2 += 30,100,200,500;  
  Xbins_etajet1 += 0,1.5,3,4.7;
  Xbins_etajet2 +=0,1.5,3,4.7;
 
  m4L = 0;
  m4L_gen = 0; 
  ngenjets = 0;
  ngencentraljets =0; 
  mjj_gen = 0;
  deta_gen = 0; 
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
  w_EtaJet2 = 0; bin_Jets=0;
 
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

void ZZRecoWAnalyzer::end( TFile &) {
   cout<<theMCInfo.analyzedEvents()<< " " << e <<endl;
  
  if(region_ == phys::CR3P1F || region_ == phys::CR2P2F) {
    vector<std::string>  FinalState = {"4m","4e","2e2m"};
    
    for (std::vector<std::string>::iterator it = FinalState.begin() ; it != FinalState.end(); ++it){
      
      TH1 *hvar =  theHistograms.get(("ZZTo"+*it+"_Mass_FRVar").c_str());
      
      TH1 *h = theHistograms.get(("ZZTo"+*it+"_Mass_01").c_str());
      
      if(!h) continue;
      for(int i = 1; i<=h->GetNbinsX();i++){
	
	Float_t Err = h->GetBinError(i);
	h->SetBinError(i,sqrt(Err*Err+hvar->GetBinContent(i)));
      }
    }
  }  
}
