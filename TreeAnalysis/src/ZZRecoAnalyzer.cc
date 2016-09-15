#include "VVXAnalysis/TreeAnalysis/interface/ZZRecoAnalyzer.h"
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

void ZZRecoAnalyzer::ZZplots(int id, int e){
  nEvent++;

  if(ZZ->id() != id && id != -1) return; // -1 here means generic 4l final state
  
  std::string decay  = "4l";
  std::string sample;
  std::string region = "";
  if      (id == 52) {decay = "4m";}
  else if (id == 48) {decay = "2e2m";}
  else if (id == 44) {decay = "4e";}
  
  if(e < nentries/2){sample = "0";}
  else {sample = "1";}
  
  Int_t njets = jets->size();
  if (njets>3) njets=3;
  Int_t ncentraljets = centralJets->size();
  if (ncentraljets>3) ncentraljets=3;
  Int_t npjets = pjets->size();
  if (npjets>3) npjets=3;
 


 // cout << npjets << endl;
  
  //////////////////////////////////////////////DATA ONLY///////////////////////////////////////////////////////
  
  //**********************************************JETS*********************************************************
  
  // //1D Reco nJet Distributions (no JER smearing, for data)
  // theHistograms.fill(std::string("ZZTo")+decay+"_Jets_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight); 
  // theHistograms.fill(std::string("ZZTo")+decay+"_Jets_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight); 
 
  float ptJet1 =0;
  float upPtJet1JESData=0;
  float downPtJet1JESData=0;
  float ptJet2 =0;
  float upPtJet2JESData=0;
  float downPtJet2JESData=0;
  float etaJet1 =0;
  float upEtaJet1JESData=0;
  float downEtaJet1JESData=0;
  float etaJet2 =0;
  float upEtaJet2JESData=0;
  float downEtaJet2JESData=0;
  float deta=0;
  float mjj=0;
  float centraldeta=0;
  float centralmjj=0;

  //1D Reco PtJet1 distributions  - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)
  if(njets>=1){
    ptJet1 = jets->at(0).pt();
    if (ptJet1>=300) ptJet1 = 299; //if (ptJet1>=500) ptJet1 = 499;
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_"+sample, "", Xbins_ptjet1, ptJet1, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_01", "", Xbins_ptjet1, ptJet1, theWeight);
    // theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_"+sample, "", 100,0,500, ptJet1, theWeight);
    // theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_01", "", 100,0,500, ptJet1, theWeight); 
    etaJet1 = fabs(jets->at(0).eta());
    if (etaJet1>=4.7) etaJet1 = 4.6;
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_"+sample, "", Xbins_etajet1, etaJet1, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_01", "", Xbins_etajet1, etaJet1, theWeight);
  }

  //1D Reco PtJet1 Distributions - JER/JES smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
  if(nUpJESDatajets>=1){
    upPtJet1JESData = UpJESData_jetPt->at(0);
    if (upPtJet1JESData>=300) upPtJet1JESData = 299;  //if (upPtJet1JESData>=500) upPtJet1JESData = 499;
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JESDataUpSmear_"+sample, "", Xbins_ptjet1, upPtJet1JESData, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JESDataUpSmear_01", "", Xbins_ptjet1, upPtJet1JESData, theWeight); 
    upEtaJet1JESData = fabs(UpJESData_jets->at(0).eta());
    if (upEtaJet1JESData>=4.7) upEtaJet1JESData = 4.6;
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JESDataUpSmear_"+sample, "", Xbins_etajet1, upEtaJet1JESData, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JESDataUpSmear_01", "", Xbins_etajet1, upEtaJet1JESData, theWeight);
  }
  if(nDownJESDatajets>=1){
    downPtJet1JESData = DownJESData_jetPt->at(0);
    if (downPtJet1JESData>=300) downPtJet1JESData = 299; //if (downPtJet1JESData>=500) downPtJet1JESData = 499;
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JESDataDownSmear_"+sample, "", Xbins_ptjet1, downPtJet1JESData, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JESDataDownSmear_01", "", Xbins_ptjet1, downPtJet1JESData, theWeight); 
    downEtaJet1JESData = fabs(DownJESData_jets->at(0).eta());
    if (downEtaJet1JESData>=4.7) downEtaJet1JESData = 4.6;
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JESDataDownSmear_"+sample, "", Xbins_etajet1, downEtaJet1JESData, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JESDataDownSmear_01", "", Xbins_etajet1, downEtaJet1JESData, theWeight);
  }
  
  //1D Reco Deta and Mjj Distributions (no JER smearing, for data)
  if(njets >=2){
    deta = fabs(jets->at(0).eta() - jets->at(1).eta());
    if (deta>=4.7) deta = 4.6;
    mjj =  (jets->at(0).p4() + jets->at(1).p4()).M();
    if (mjj>=800) mjj = 799;
    ptJet2 = jets->at(1).pt();
    if (ptJet2>=200) ptJet2 = 199; //if (ptJet2>=500) ptJet2 = 499;
    float dphi = deltaPhi(jets->at(0).phi(),jets->at(1).phi());
    if (dphi>=6) dphi = 5;
    
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_"+sample, "", Xbins_ptjet2, ptJet2, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_01", "", Xbins_ptjet2, ptJet2, theWeight);
    // theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_"+sample, "", 100,0,500, ptJet2, theWeight);
    // theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_01", "",100,0,500 , ptJet2, theWeight);
    etaJet2 = fabs(jets->at(1).eta());
    if (etaJet2>=4.7) etaJet2 = 4.6;
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_"+sample, "", Xbins_etajet2, etaJet2, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_01", "", Xbins_etajet2, etaJet2, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_"+sample,"",Xbins_mjj,mjj,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_01","",Xbins_mjj,mjj,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_"+sample,"",Xbins_deta,deta,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_01","",Xbins_deta,deta,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_Phi_01","",Xbins_dphi,dphi,theWeight); 
  }
  
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESDataUpSmear_"+sample, "", 4,0,4,nUpJESDatajets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESDataDownSmear_"+sample, "", 4,0,4,nDownJESDatajets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESDataUpSmear_01", "", 4,0,4,nUpJESDatajets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESDataDownSmear_01", "", 4,0,4,nDownJESDatajets , theWeight);
  
  if(nUpJESDatajets >=2){
    float upDetaJESData = fabs(UpJESData_jets->at(0).eta() - UpJESData_jets->at(1).eta());
    float upMjjJESData =  (UpJESData_jets->at(0).p4() + UpJESData_jets->at(1).p4()).M();
    upPtJet2JESData = UpJESData_jetPt->at(1); 
    upEtaJet2JESData = fabs(UpJESData_jets->at(1).eta());
    if(upDetaJESData>=4.7)upDetaJESData= 4.6;
    if(upMjjJESData>=800) upMjjJESData = 799;
    if (upPtJet2JESData>=200) upPtJet2JESData = 199;// if (upPtJet2JESData>=500) upPtJet2JESData = 499;
    if (upEtaJet2JESData>=4.7) upEtaJet2JESData = 4.6;
    
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JESDataUpSmear_"+sample,"",Xbins_mjj,upMjjJESData,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JESDataUpSmear_01","",Xbins_mjj,upMjjJESData,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JESDataUpSmear_"+sample,"",Xbins_deta,upDetaJESData,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JESDataUpSmear_01","",Xbins_deta,upDetaJESData,theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JESDataUpSmear_"+sample, "", Xbins_ptjet2, upPtJet2JESData, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JESDataUpSmear_01", "", Xbins_ptjet2, upPtJet2JESData, theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JESDataUpSmear_"+sample, "", Xbins_etajet2, upEtaJet2JESData, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JESDataUpSmear_01", "", Xbins_etajet2, upEtaJet2JESData, theWeight);
  }
  if(nDownJESDatajets >=2){
    float downDetaJESData = fabs(DownJESData_jets->at(0).eta() - DownJESData_jets->at(1).eta());
    float downMjjJESData =  (DownJESData_jets->at(0).p4() + DownJESData_jets->at(1).p4()).M();
    downPtJet2JESData = DownJESData_jetPt->at(1);
    downEtaJet2JESData = fabs(DownJESData_jets->at(1).eta());
    if(downDetaJESData>=4.7)downDetaJESData= 4.6;
    if(downMjjJESData>=800) downMjjJESData = 799;
    if (downPtJet2JESData>=200) downPtJet2JESData = 199; // if (downPtJet2JESData>=500) downPtJet2JESData = 499; 
    if (downEtaJet2JESData>=4.7) downEtaJet2JESData = 4.6;
    
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JESDataDownSmear_"+sample,"",Xbins_mjj,downMjjJESData,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JESDataDownSmear_01","",Xbins_mjj,downMjjJESData,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JESDataDownSmear_"+sample,"",Xbins_deta,downDetaJESData,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JESDataDownSmear_01","",Xbins_deta,downDetaJESData,theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JESDataDownSmear_"+sample, "", Xbins_ptjet2, downPtJet2JESData, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JESDataDownSmear_01", "", Xbins_ptjet2, downPtJet2JESData, theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JESDataDownSmear_"+sample, "", Xbins_etajet2, downEtaJet2JESData, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JESDataDownSmear_01", "", Xbins_etajet2, downEtaJet2JESData, theWeight); 
  }
  
  
  
  //**********************************************CENTRAL JETS*********************************************************
  
  // //1D Reco nJet Distributions (no JER smearing, for data)
  // theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_"+sample, std::string("Number of centraljets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theWeight); 
  // theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_01", std::string("Number of centraljets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theWeight); 
  
  //1D Reco Deta and Mjj Distributions (no JER smearing, for data)
  if(ncentraljets >=2){
    centraldeta = fabs(centralJets->at(0).eta() - centralJets->at(1).eta());
    if (centraldeta >=4.7) centraldeta = 4.6;
    centralmjj =  (centralJets->at(0).p4() + centralJets->at(1).p4()).M();
    if (centralmjj>=800) centralmjj = 799;
    
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_"+sample,"",Xbins_mjj,centralmjj,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_01","",Xbins_mjj,centralmjj,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_"+sample,"",Xbins_deta,centraldeta,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_01","",Xbins_deta,centraldeta,theWeight); 
    
  }
  
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JESDataUpSmear_"+sample, "", 4,0,4,nUpJESDatacentraljets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JESDataDownSmear_"+sample, "", 4,0,4,nDownJESDatacentraljets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JESDataUpSmear_01", "", 4,0,4,nUpJESDatacentraljets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JESDataDownSmear_01", "", 4,0,4,nDownJESDatacentraljets , theWeight);
  
  if(nUpJESDatacentraljets >=2){
    float upCentraldetaJESData = fabs(UpJESData_centraljets->at(0).eta() - UpJESData_centraljets->at(1).eta());
    float upCentralmjjJESData =  (UpJESData_centraljets->at(0).p4() + UpJESData_centraljets->at(1).p4()).M();
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JESDataUpSmear_"+sample,"",Xbins_mjj,upCentralmjjJESData,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JESDataUpSmear_01","",Xbins_mjj,upCentralmjjJESData,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JESDataUpSmear_"+sample,"",Xbins_deta,upCentraldetaJESData,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JESDataUpSmear_01","",Xbins_deta,upCentraldetaJESData,theWeight); 
  }
  if(nDownJESDatacentraljets >=2){
    float downCentraldetaJESData = fabs(DownJESData_centraljets->at(0).eta() - DownJESData_centraljets->at(1).eta());
    float downCentralmjjJESData =  (DownJESData_centraljets->at(0).p4() + DownJESData_centraljets->at(1).p4()).M();
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JESDataDownSmear_"+sample,"",Xbins_mjj,downCentralmjjJESData,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JESDataDownSmear_01","",Xbins_mjj,downCentralmjjJESData,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JESDataDownSmear_"+sample,"",Xbins_deta,downCentraldetaJESData,theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JESDataDownSmear_01","",Xbins_deta,downCentraldetaJESData,theWeight); 
  }
  
  // ***********************************************************************************************************************************************************
  
  ////////////////////////////////////////////////////DATA AND MONTECARLO/////////////////////////////////////////////////
  
  //NNLO/NLO k_factor w_kf = 0;
 
  string sampleName;
  //if powheg samples
  sampleName = "ZZTo4mu";
  //else FIXME
  //sampleName = "oth";
  
  w_kf = 0;
  dphizz_gen = 0;
  dphizz = 0;
  dphizz = fabs(physmath::deltaPhi(ZZ->first().phi(),ZZ->second().phi()));
  w_kf = 1;

  // if (genCategory !=-1){
  //   if(topology.test(0)){  
  //     dphizz_gen = fabs(physmath::deltaPhi(genVBParticles->at(0).phi(),genVBParticles->at(1).phi())); 
  //     //w_kf = kfNNLO.weight(sampleName.c_str(), decay.c_str(), dphizz_gen,nEvent); 
  //     w_kf = kfNNLO.weight(sampleName.c_str(), dphizz_gen, nEvent); 

  //     //if(dphizz_gen >=3.1 && dphizz_gen <= TMath::Pi()) cout << dphizz_gen << endl;
  //     // cout << "GEN OK: "<< ZZ->first().phi()<< " " << ZZ->second().phi() << " " <<dphizz << " " << genVBParticles->at(0).phi() << " " << genVBParticles->at(1).phi() << " " << dphizz_gen << " " << w_kf << " " << endl;
  //   }
  //   else {
  //     w_kf =1;
  //     // cout << "NO GEN: " <<ZZ->first().phi()<< " " << ZZ->second().phi() << " " <<dphizz << " " << w_kf << endl; 
  //   }
  // }
  // else  {
  //   w_kf =1;
  //   //cout << w_kf << endl;
  // }

  // //if(w_kf==0 || w_kf < 0) cout << "******************** " << w_kf << endl;

  // w_kf = 1.15;



  //1D Reco Mass Distributions
  m4L = ZZ->mass();  
  ptzz = ZZ->pt();
  drzz =physmath::deltaR(ZZ->first(),ZZ->second());
  //if(drzz > 6) drzz = 5.9; //overflow bin
  if(m4L > 800) m4L = 799;
  theHistograms.fill(std::string("ZZTo")+decay+"_Mass_01","", Xbins, m4L,theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_Mass_"+sample,"", Xbins, m4L,theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_dRZZ_01","", Xbins_drzz, drzz,theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_dRZZ_"+sample,"", Xbins_drzz, drzz,theWeight*w_kf);
  if(njets >=2){
    theHistograms.fill(std::string("ZZTo")+decay+"_ZZMass2j_01","", Xbins, m4L,theWeight*w_kf);
  }

  if(ptzz>=300) ptzz = 299;
  theHistograms.fill(std::string("ZZTo")+decay+"_PtZZ_01", Xbins_ptzz, ptzz,theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_PtZZ_"+sample, Xbins_ptzz, ptzz,theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZ_01", Xbins_dphizz, dphizz,theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZ_"+sample, Xbins_dphizz, dphizz,theWeight*w_kf);

  // //Pt and eta distribution for the leading, sub-leading and sub-sub-leading jets (if they exist)
 //  if(npjets>0){  
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Pt0_01", std::string("p_{t}^{jet0} of ZZ_{1}#rightarrow ")+decay,50,0,350,pjets->at(0).pt(), theWeight*w_kf);  
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Eta0_01", std::string("eta^{jet0} of ZZ_{1}#rightarrow ")+decay,15,-6,6,pjets->at(0).eta(), theWeight*w_kf);  
 //  }
 //  if(npjets>1){  
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Pt1_01", std::string("p_{t}^{jet1} of ZZ_{1}#rightarrow ")+decay,50,0,350,pjets->at(1).pt(), theWeight*w_kf);  
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Eta1_01", std::string("eta^{jet1} of ZZ_{1}#rightarrow ")+decay,15,-6,6,pjets->at(1).eta(), theWeight*w_kf);  
 //  }
 //  if(npjets>2){  
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Pt2_01", std::string("p_{t}^{jet2} of ZZ_{1}#rightarrow ")+decay,50,0,350,pjets->at(2).pt(), theWeight*w_kf);  
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Eta2_01", std::string("eta^{jet2} of ZZ_{1}#rightarrow ")+decay,15,-6,6,pjets->at(2).eta(), theWeight*w_kf);  
 //  }
  

  //1D Reco nJet Distributions (no JER smearing, for data)
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_"+sample,"",4,0,4,njets,theWeight*w_kf); 
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_01","",4,0,4,njets,theWeight*w_kf); 

 //1D Reco ncentralJet Distributions (no JER smearing, for data)
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_"+sample,"",4,0,4,ncentraljets,theWeight*w_kf); 
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_01","",4,0,4,ncentraljets,theWeight*w_kf); 
  
  //1D Reco nPJet Distributions (no JER smearing, for data)
  theHistograms.fill(std::string("ZZTo")+decay+"_PJets_"+sample,"",4,0,4,npjets,theWeight*w_kf); 
  theHistograms.fill(std::string("ZZTo")+decay+"_PJets_01","",4,0,4,npjets,theWeight*w_kf); 
  
 
  /////////////////////////////////////////////////////MONTECARLO ONLY//////////////////////////////////////////////////////


  //******************************************************************RECO JETS****************************************************************************************

  //1D Reco nJet Distributions - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_"+sample, "", 4,0,4,nCentralJERjets , theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_01", "", 4,0,4,nCentralJERjets , theWeight*w_kf);
  
  //1D Reco nJet Distributions - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERUpSmear_"+sample, "", 4,0,4,nUpJERjets , theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERDownSmear_"+sample, "", 4,0,4,nDownJERjets , theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERUpSmear_01", "", 4,0,4,nUpJERjets , theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERDownSmear_01", "", 4,0,4,nDownJERjets , theWeight*w_kf);

 //1D Reco nJet Distributions - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESUpSmear_"+sample, "", 4,0,4,nUpJESjets , theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESDownSmear_"+sample, "", 4,0,4,nDownJESjets , theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESUpSmear_01", "", 4,0,4,nUpJESjets , theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESDownSmear_01", "", 4,0,4,nDownJESjets , theWeight*w_kf);
  
  // //Pt and eta distribution for the leading, sub-leading and sub-sub-leading jets (if they exist)
  // if(nCentralJERjets>0){  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Pt0_JERCentralSmear_01","",50,0,350,pjets->at(0).pt(), theWeight*w_kf);  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Eta0_JERCentralSmear_01","",30,0,6,pjets->at(0).eta(), theWeight*w_kf);  
  // }
  // if(nCentralJERjets>1){  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Pt1_JERCentralSmear_01","",50,0,350,pjets->at(1).pt(), theWeight*w_kf);  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Eta1_JERCentralSmear_01","",30,0,6,pjets->at(1).eta(), theWeight*w_kf);  
  // }
  // if(nCentralJERjets>2){  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Pt2_JERCentralSmear_01","",50,0,350,pjets->at(2).pt(), theWeight*w_kf);  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Eta2_JERCentralSmear_01","",30,0,6,pjets->at(2).eta(), theWeight*w_kf);  
  // }

  float centralDeta = 0;
  float centralMjj = 0;
  float upDetaJER  = 0;
  float upMjjJER = 0;  
  float downDetaJER  = 0;
  float downMjjJER  = 0;
  float upDetaJES  = 0;
  float upMjjJES = 0;
  float downDetaJES = 0; 
  float downMjjJES = 0;

  float centralPtJet1 = 0;
  float centralPtJet2 = 0;
  float upPtJet1JER  = 0;
  float upPtJet2JER = 0;  
  float downPtJet1JER  = 0;
  float downPtJet2JER  = 0;
  float upPtJet1JES  = 0;
  float upPtJet2JES = 0;
  float downPtJet1JES = 0; 
  float downPtJet2JES = 0;

  float centralEtaJet1 = 0;
  float centralEtaJet2 = 0;
  float upEtaJet1JER  = 0;
  float upEtaJet2JER = 0;  
  float downEtaJet1JER  = 0;
  float downEtaJet2JER  = 0;
  float upEtaJet1JES  = 0;
  float upEtaJet2JES = 0;
  float downEtaJet1JES = 0; 
  float downEtaJet2JES = 0;

  float centralDphi =0;

  //1D Reco PtJet1 distributions  - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)
  if(nCentralJERjets>=1){
    centralPtJet1 = CentralJER_jetPt->at(0);
    if (centralPtJet1>=300) centralPtJet1 = 299; //if (centralPtJet1>=500) centralPtJet1 = 499;
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JERCentralSmear_"+sample, "", Xbins_ptjet1, centralPtJet1, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JERCentralSmear_01", "", Xbins_ptjet1, centralPtJet1, theWeight*w_kf); 
    centralEtaJet1 = fabs(CentralJER_jets->at(0).eta());  
   
    if (centralEtaJet1>=4.7) centralEtaJet1 = 4.6;
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JERCentralSmear_"+sample, "", Xbins_etajet1, centralEtaJet1, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JERCentralSmear_01", "", Xbins_etajet1, centralEtaJet1, theWeight*w_kf);
  }

  //1D Reco DeltaEta and mJJ Distributions - JER/JES smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
  if(nUpJERjets>=1){
    upPtJet1JER = UpJER_jetPt->at(0);
    if (upPtJet1JER>=300) upPtJet1JER = 299;//if (upPtJet1JER>=500) upPtJet1JER = 499;
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JERUpSmear_"+sample, "", Xbins_ptjet1, upPtJet1JER, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JERUpSmear_01", "", Xbins_ptjet1, upPtJet1JER, theWeight*w_kf); 
    upEtaJet1JER = fabs(UpJER_jets->at(0).eta());
    if (upEtaJet1JER>=4.7) upEtaJet1JER = 4.6;
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JERUpSmear_"+sample, "", Xbins_etajet1, upEtaJet1JER, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JERUpSmear_01", "", Xbins_etajet1, upEtaJet1JER, theWeight*w_kf);
  }
  if(nDownJERjets>=1){
    downPtJet1JER = DownJER_jetPt->at(0);
    if (downPtJet1JER>=300) downPtJet1JER = 299; //if (downPtJet1JER>=500) downPtJet1JER = 499;
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JERDownSmear_"+sample, "", Xbins_ptjet1, downPtJet1JER, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JERDownSmear_01", "", Xbins_ptjet1, downPtJet1JER, theWeight*w_kf); 
    downEtaJet1JER = fabs( DownJER_jets->at(0).eta());
    if (downEtaJet1JER>=4.7) downEtaJet1JER = 4.6;
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JERDownSmear_"+sample, "", Xbins_etajet1, downEtaJet1JER, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JERDownSmear_01", "", Xbins_etajet1, downEtaJet1JER, theWeight*w_kf);
  }
  if(nUpJESjets>=1){
    upPtJet1JES = UpJES_jetPt->at(0);
    if (upPtJet1JES>=300) upPtJet1JES= 299;
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JESUpSmear_"+sample, "", Xbins_ptjet1, upPtJet1JES, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JESUpSmear_01", "", Xbins_ptjet1, upPtJet1JES, theWeight*w_kf); 
    upEtaJet1JES =  fabs(UpJES_jets->at(0).eta());
    if (upEtaJet1JES>=4.7) upEtaJet1JES= 4.6;
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JESUpSmear_"+sample, "", Xbins_etajet1, upEtaJet1JES, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JESUpSmear_01", "", Xbins_etajet1, upEtaJet1JES, theWeight*w_kf);
  }
  if(nDownJESjets>=1){
    downPtJet1JES = DownJES_jetPt->at(0);
    if (downPtJet1JES>=300) downPtJet1JES = 299;
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JESDownSmear_"+sample, "", Xbins_ptjet1, downPtJet1JES, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1_JESDownSmear_01", "", Xbins_ptjet1, downPtJet1JES, theWeight*w_kf);
    downEtaJet1JES = fabs( DownJES_jets->at(0).eta());
    if (downEtaJet1JES>=4.7) downEtaJet1JES = 4.6;
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JESDownSmear_"+sample, "", Xbins_etajet1, downEtaJet1JES, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1_JESDownSmear_01", "", Xbins_etajet1, downEtaJet1JES, theWeight*w_kf);
  }

  //1D Reco DeltaEta and mJJ Distributions - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)
   if(nCentralJERjets>=2){
  
    centralDeta = fabs(CentralJER_jets->at(0).eta() - CentralJER_jets->at(1).eta());
    centralMjj =  (CentralJER_jets->at(0).p4() + CentralJER_jets->at(1).p4()).M();
    centralPtJet2 = CentralJER_jetPt->at(1);
    centralEtaJet2 = fabs( CentralJER_jets->at(1).eta()); 
    centralDphi = deltaPhi(CentralJER_jets->at(0).phi(),CentralJER_jets->at(1).phi());
   
    if (centralDeta>=4.7) centralDeta = 4.6;
    if (centralMjj>=800) centralMjj = 799;
    if (centralPtJet2>=200) centralPtJet2 = 199; 
    if (centralEtaJet2>=4.7) centralEtaJet2 = 4.6;
    if (centralDphi>=6)  centralDphi= 5;

    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERCentralSmear_"+sample, "reco m_{jj}", Xbins_mjj, centralMjj, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERCentralSmear_01", "reco m_{jj}", Xbins_mjj, centralMjj, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERCentralSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, centralDeta, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERCentralSmear_01", "reco  #Delta#eta_{jj}", Xbins_deta, centralDeta, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JERCentralSmear_"+sample, "", Xbins_ptjet2, centralPtJet2, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JERCentralSmear_01", "", Xbins_ptjet2, centralPtJet2, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JERCentralSmear_"+sample, "", Xbins_etajet2, centralEtaJet2, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JERCentralSmear_01", "", Xbins_etajet2, centralEtaJet2, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_Phi_JERCentralSmear_01", "", Xbins_dphi, centralDphi, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_ZZMass2j_JERCentralSmear_01", std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L,theWeight*w_kf);
  }
  
  //1D Reco DeltaEta and mJJ Distributions - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
  if(nUpJERjets>=2){

    upDetaJER = fabs(UpJER_jets->at(0).eta() - UpJER_jets->at(1).eta());
    upMjjJER =  (UpJER_jets->at(0).p4() + UpJER_jets->at(1).p4()).M();
    upPtJet2JER = UpJER_jetPt->at(1); 
    upEtaJet2JER =  fabs(UpJER_jets->at(1).eta());

    if (upDetaJER>=4.7) upDetaJER = 4.6;
    if (upMjjJER>=800) upMjjJER = 799;
    if (upPtJet2JER>=200) upPtJet2JER = 199; 
    if (upEtaJet2JER>=4.7) upEtaJet2JER = 4.6;

    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERUpSmear_"+sample, "reco m_{jj}", Xbins_mjj, upMjjJER, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERUpSmear_01", "reco m_{jj}", Xbins_mjj, upMjjJER, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERUpSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, upDetaJER, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERUpSmear_01", "reco #Delta#eta_{jj}", Xbins_deta, upDetaJER, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JERUpSmear_"+sample, "", Xbins_ptjet2, upPtJet2JER, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JERUpSmear_01", "", Xbins_ptjet2, upPtJet2JER, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JERUpSmear_"+sample, "", Xbins_etajet2, upEtaJet2JER, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JERUpSmear_01", "", Xbins_etajet2, upEtaJet2JER, theWeight*w_kf);
  }
 
  if(nDownJERjets>=2){

    downDetaJER = fabs(DownJER_jets->at(0).eta() - DownJER_jets->at(1).eta());
    downMjjJER =  (DownJER_jets->at(0).p4() + DownJER_jets->at(1).p4()).M();
    downPtJet2JER = DownJER_jetPt->at(1);
    downEtaJet2JER =  fabs(DownJER_jets->at(1).eta()); 

    if (downDetaJER>=4.7) downDetaJER = 4.6;
    if (downMjjJER>=800) downMjjJER = 799;
    if (downPtJet2JER>=200) downPtJet2JER = 199; 
    if (downEtaJet2JER>=4.7) downEtaJet2JER = 4.6;

    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERDownSmear_"+sample, "reco m_{jj}", Xbins_mjj, downMjjJER, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERDownSmear_01", "reco m_{jj}", Xbins_mjj, downMjjJER, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERDownSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, downDetaJER, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERDownSmear_01", "reco #Delta#eta_{jj}", Xbins_deta, downDetaJER, theWeight*w_kf);  
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JERDownSmear_"+sample, "", Xbins_ptjet2, downPtJet2JER, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JERDownSmear_01", "", Xbins_ptjet2, downPtJet2JER, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JERDownSmear_"+sample, "", Xbins_etajet2, downEtaJet2JER, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JERDownSmear_01", "", Xbins_etajet2, downEtaJet2JER, theWeight*w_kf);

  }
 
  //1D Reco DeltaEta and mJJ Distributions - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
  if(nUpJESjets>=2){

    upDetaJES = fabs(UpJES_jets->at(0).eta() - UpJES_jets->at(1).eta());
    upMjjJES =  (UpJES_jets->at(0).p4() + UpJES_jets->at(1).p4()).M();
    upPtJet2JES = UpJES_jetPt->at(1); 
    upEtaJet2JES =  fabs(UpJES_jets->at(1).eta());

    if (upDetaJES>=4.7) upDetaJES = 4.6;
    if (upMjjJES>=800) upMjjJES = 799;
    if (upPtJet2JES>=200) upPtJet2JES = 199; 
    if (upEtaJet2JES>=4.7) upEtaJet2JES = 4.6;

    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JESUpSmear_"+sample, "reco m_{jj}", Xbins_mjj, upMjjJES, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JESUpSmear_01", "reco m_{jj}", Xbins_mjj, upMjjJES, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JESUpSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, upDetaJES, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JESUpSmear_01", "reco #Delta#eta_{jj}", Xbins_deta, upDetaJES, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JESUpSmear_"+sample, "", Xbins_ptjet2, upPtJet2JES, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JESUpSmear_01", "", Xbins_ptjet2, upPtJet2JES, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JESUpSmear_"+sample, "", Xbins_etajet2, upEtaJet2JES, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JESUpSmear_01", "", Xbins_etajet2, upEtaJet2JES, theWeight*w_kf);
  }
 
  if(nDownJESjets>=2){

    downDetaJES = fabs(DownJES_jets->at(0).eta() - DownJES_jets->at(1).eta());
    downMjjJES =  (DownJES_jets->at(0).p4() + DownJES_jets->at(1).p4()).M();
    downPtJet2JES = DownJES_jetPt->at(1); 
    downEtaJet2JES =  fabs(DownJES_jets->at(1).eta());

    if (downDetaJES>=4.7) downDetaJES = 4.6;
    if (downMjjJES>=800) downMjjJES = 799;
    if (downPtJet2JES>=200) downPtJet2JES = 199; 
    if (downEtaJet2JES>=4.7) downEtaJet2JES = 4.6;

    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JESDownSmear_"+sample, "reco m_{jj}", Xbins_mjj, downMjjJES, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JESDownSmear_01", "reco m_{jj}", Xbins_mjj, downMjjJES, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JESDownSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, downDetaJES, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JESDownSmear_01", "reco #Delta#eta_{jj}", Xbins_deta, downDetaJES, theWeight*w_kf);  
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JESDownSmear_"+sample, "", Xbins_ptjet2, downPtJet2JES, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2_JESDownSmear_01", "", Xbins_ptjet2, downPtJet2JES, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JESDownSmear_"+sample, "", Xbins_etajet2, downEtaJet2JES, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2_JESDownSmear_01", "", Xbins_etajet2, downEtaJet2JES, theWeight*w_kf);
  } 


 //******************************************************************RECO CENTRAL JETS****************************************************************************************

  //1D Reco nJet Distributions - JER smearing (CentralJets_JERCentralSmear to be used in the standard analysis)
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERCentralSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERcentraljets , theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERCentralSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERcentraljets , theWeight*w_kf);
  
  //1D Reco nJet Distributions - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERUpSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJERcentraljets , theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERDownSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJERcentraljets , theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERUpSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJERcentraljets , theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERDownSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJERcentraljets , theWeight*w_kf);

 //1D Reco nJet Distributions - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JESUpSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJEScentraljets , theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JESDownSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJEScentraljets , theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JESUpSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJEScentraljets , theWeight*w_kf);
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JESDownSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJEScentraljets , theWeight*w_kf);
  
  // //Pt and eta distribution for the leading, sub-leading and sub-sub-leading centraljets (if they exist)
  // if(nCentralJERcentraljets>0){  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Pt0_JERCentralSmear_01", std::string("p_{t}^{jet0} of ZZ_{1}#rightarrow ")+decay,50,0,350,FIXMEpjets->at(0).pt(), theWeight*w_kf);  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Eta0_JERCentralSmear_01", std::string("eta^{jet0} of ZZ_{1}#rightarrow ")+decay,30,0,4.7,pjets->at(0).eta(), theWeight*w_kf);  
  // }
  // if(nCentralJERcentraljets>1){  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Pt1_JERCentralSmear_01", std::string("p_{t}^{jet1} of ZZ_{1}#rightarrow ")+decay,50,0,350,pjets->at(1).pt(), theWeight*w_kf);  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Eta1_JERCentralSmear_01", std::string("eta^{jet1} of ZZ_{1}#rightarrow ")+decay,30,0,4.7,pjets->at(1).eta(), theWeight*w_kf);  
  // }
  // if(nCentralJERcentraljets>2){  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Pt2_JERCentralSmear_01", std::string("p_{t}^{jet2} of ZZ_{1}#rightarrow ")+decay,50,0,350,pjets->at(2).pt(), theWeight*w_kf);  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Eta2_JERCentralSmear_01", std::string("eta^{jet2} of ZZ_{1}#rightarrow ")+decay,30,0,4.7,pjets->at(2).eta(), theWeight*w_kf);  
  // }

  float centralDeta_cj = 0;
  float centralMjj_cj = 0;
  float upDetaJER_cj  = 0;
  float upMjjJER_cj = 0;  
  float downDetaJER_cj  = 0;
  float downMjjJER_cj  = 0;
  float upDetaJES_cj  = 0;
  float upMjjJES_cj = 0;
  float downDetaJES_cj = 0; 
  float downMjjJES_cj = 0;

  //1D Reco DeltaEta and mJJ Distributions - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)

  if(nCentralJERcentraljets>=2){
  
    centralDeta_cj = fabs(CentralJER_centraljets->at(0).eta() - CentralJER_centraljets->at(1).eta());
    centralMjj_cj =  (CentralJER_centraljets->at(0).p4() + CentralJER_centraljets->at(1).p4()).M();
    
    if (centralDeta_cj>=4.7) centralDeta_cj = 4.6;
    if (centralMjj_cj>=800) centralMjj_cj = 799;

    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERCentralSmear_"+sample, "reco m_{jj}", Xbins_mjj, centralMjj_cj, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERCentralSmear_01", "reco m_{jj}", Xbins_mjj, centralMjj_cj, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERCentralSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, centralDeta_cj, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERCentralSmear_01", "reco  #Delta#eta_{jj}", Xbins_deta, centralDeta_cj, theWeight*w_kf);
  }
  
  //1D Reco DeltaEta and mJJ Distributions - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
  if(nUpJERcentraljets>=2){
    
    upDetaJER_cj = fabs(UpJER_centraljets->at(0).eta() - UpJER_centraljets->at(1).eta());
    upMjjJER_cj =  (UpJER_centraljets->at(0).p4() + UpJER_centraljets->at(1).p4()).M();
    
    if (upDetaJER_cj>=4.7) upDetaJER_cj = 4.6;
    if (upMjjJER_cj>=800) upMjjJER_cj = 799;
    
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERUpSmear_"+sample, "reco m_{jj}", Xbins_mjj, upMjjJER_cj, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERUpSmear_01", "reco m_{jj}", Xbins_mjj, upMjjJER_cj, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERUpSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, upDetaJER_cj, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERUpSmear_01", "reco #Delta#eta_{jj}", Xbins_deta, upDetaJER_cj, theWeight*w_kf);
  }

  if(nDownJERcentraljets>=2){

    downDetaJER_cj = fabs(DownJER_centraljets->at(0).eta() - DownJER_centraljets->at(1).eta());
    downMjjJER_cj =  (DownJER_centraljets->at(0).p4() + DownJER_centraljets->at(1).p4()).M();
 
    if (downDetaJER_cj>=4.7) downDetaJER_cj = 4.6;
    if (downMjjJER_cj>=800) downMjjJER_cj = 799;
    
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERDownSmear_"+sample, "reco m_{jj}", Xbins_mjj, downMjjJER_cj, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERDownSmear_01", "reco m_{jj}", Xbins_mjj, downMjjJER_cj, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERDownSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, downDetaJER_cj, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERDownSmear_01", "reco #Delta#eta_{jj}", Xbins_deta, downDetaJER_cj, theWeight*w_kf);
  }
 
  //1D Reco DeltaEta and mJJ Distributions - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
  if(nUpJEScentraljets>=2){

    upDetaJES_cj = fabs(UpJES_centraljets->at(0).eta() - UpJES_centraljets->at(1).eta());
    upMjjJES_cj =  (UpJES_centraljets->at(0).p4() + UpJES_centraljets->at(1).p4()).M();

    if (upDetaJES_cj>=4.7) upDetaJES_cj = 4.6;
    if (upMjjJES_cj>=800) upMjjJES_cj = 799;

    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JESUpSmear_"+sample, "reco m_{jj}", Xbins_mjj, upMjjJES_cj, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JESUpSmear_01", "reco m_{jj}", Xbins_mjj, upMjjJES_cj, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JESUpSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, upDetaJES_cj, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JESUpSmear_01", "reco #Delta#eta_{jj}", Xbins_deta, upDetaJES_cj, theWeight*w_kf);
  }
  
  if(nDownJEScentraljets>=2){

    downDetaJES_cj = fabs(DownJES_centraljets->at(0).eta() - DownJES_centraljets->at(1).eta());
    downMjjJES_cj =  (DownJES_centraljets->at(0).p4() + DownJES_centraljets->at(1).p4()).M();

    if (downDetaJES_cj>=4.7) downDetaJES_cj = 4.6;
    if (downMjjJES_cj>=800) downMjjJES_cj = 799;

    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JESDownSmear_"+sample, "reco m_{jj}", Xbins_mjj, downMjjJES_cj, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JESDownSmear_01", "reco m_{jj}", Xbins_mjj, downMjjJES_cj, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JESDownSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, downDetaJES_cj, theWeight*w_kf);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JESDownSmear_01", "reco #Delta#eta_{jj}", Xbins_deta, downDetaJES_cj, theWeight*w_kf);
  }
    
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
 
  /////////////////////////////////////////////////////SCALE FACTOR HISTOGRAMS/////////////////////////////////////////////////////////////////////////////////////
 theHistograms.fill(std::string("ZZTo")+decay+"_MassSFErrSqMinus_"+sample, "", Xbins , m4L, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_MassSFErrSqMinus_01", "", Xbins , m4L, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_MassSFErrSqPlus_"+sample, "", Xbins , m4L, theWeight*w_kf*(1+scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_MassSFErrSqPlus_01","", Xbins , m4L, theWeight*w_kf*(1+scaleFacErrSq));

 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZSFErrSqMinus_"+sample, "", Xbins_drzz , drzz, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZSFErrSqMinus_01", "", Xbins_drzz , drzz, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZSFErrSqPlus_"+sample, "", Xbins_drzz , drzz, theWeight*w_kf*(1+scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZSFErrSqPlus_01","", Xbins_drzz , drzz, theWeight*w_kf*(1+scaleFacErrSq));

 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZSFErrSqMinus_"+sample, "", Xbins_ptzz , ptzz, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZSFErrSqMinus_01", "", Xbins_ptzz , ptzz, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZSFErrSqPlus_"+sample, "", Xbins_ptzz , ptzz, theWeight*w_kf*(1+scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZSFErrSqPlus_01","", Xbins_ptzz , ptzz, theWeight*w_kf*(1+scaleFacErrSq)); 

 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZSFErrSqMinus_"+sample, "", Xbins_dphizz , dphizz, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZSFErrSqMinus_01", "", Xbins_dphizz , dphizz, theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZSFErrSqPlus_"+sample, "", Xbins_dphizz , dphizz, theWeight*w_kf*(1+scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZSFErrSqPlus_01","", Xbins_dphizz , dphizz, theWeight*w_kf*(1+scaleFacErrSq));

 theHistograms.fill(std::string("ZZTo")+decay+"_JetsSFErrSqMinus_"+sample,"",4,0,4,nCentralJERjets,theWeight*w_kf*(1-scaleFacErrSq)); 
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsSFErrSqMinus_01","",4,0,4,nCentralJERjets,theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsSFErrSqPlus_"+sample,"",4,0,4,nCentralJERjets,theWeight*w_kf*(1+scaleFacErrSq));  
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsSFErrSqPlus_01","",4,0,4,nCentralJERjets,theWeight*w_kf*(1+scaleFacErrSq)); 
 
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsSFErrSqMinus_"+sample,"",4,0,4,nCentralJERcentraljets,theWeight*w_kf*(1-scaleFacErrSq)); 
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsSFErrSqMinus_01","",4,0,4,nCentralJERcentraljets,theWeight*w_kf*(1-scaleFacErrSq));   
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsSFErrSqPlus_"+sample,"",4,0,4,nCentralJERcentraljets,theWeight*w_kf*(1+scaleFacErrSq));  
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsSFErrSqPlus_01", "",4,0,4,nCentralJERcentraljets,theWeight*w_kf*(1+scaleFacErrSq));     
 
 // Histograms for scale factor correlated errors
 theHistograms.fill(std::string("ZZTo")+decay+"_MassSFErrMinus_"+sample,"", Xbins , m4L, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_MassSFErrMinus_01","", Xbins , m4L, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_MassSFErrPlus_"+sample,"", Xbins , m4L, theWeight*w_kf*(1+scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_MassSFErrPlus_01","", Xbins , m4L, theWeight*w_kf*(1+scaleFacErr)); 

 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZSFErrMinus_"+sample,"", Xbins_drzz , drzz, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZSFErrMinus_01","", Xbins_drzz , drzz, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZSFErrPlus_"+sample,"", Xbins_drzz , drzz, theWeight*w_kf*(1+scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_dRZZSFErrPlus_01","", Xbins_drzz , drzz, theWeight*w_kf*(1+scaleFacErr)); 

 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZSFErrMinus_"+sample,"", Xbins_ptzz , ptzz, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZSFErrMinus_01","", Xbins_ptzz , ptzz, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZSFErrPlus_"+sample,"", Xbins_ptzz , ptzz, theWeight*w_kf*(1+scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_PtZZSFErrPlus_01","", Xbins_ptzz , ptzz, theWeight*w_kf*(1+scaleFacErr)); 

 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZSFErrMinus_"+sample,"", Xbins_dphizz , dphizz, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZSFErrMinus_01","", Xbins_dphizz , dphizz, theWeight*w_kf*(1-scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZSFErrPlus_"+sample,"", Xbins_dphizz , dphizz, theWeight*w_kf*(1+scaleFacErr));   
 theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZSFErrPlus_01","", Xbins_dphizz , dphizz, theWeight*w_kf*(1+scaleFacErr));

 theHistograms.fill(std::string("ZZTo")+decay+"_JetsSFErrMinus_"+sample,"",4,0,4,nCentralJERjets,theWeight*w_kf*(1-scaleFacErr));  
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsSFErrMinus_01","",4,0,4,nCentralJERjets,theWeight*w_kf*(1-scaleFacErr));  
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsSFErrPlus_"+sample,"",4,0,4,nCentralJERjets,theWeight*w_kf*(1+scaleFacErr));  
 theHistograms.fill(std::string("ZZTo")+decay+"_JetsSFErrPlus_01","",4,0,4,nCentralJERjets,theWeight*w_kf*(1+scaleFacErr));   
 
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsSFErrMinus_"+sample,"",4,0,4,nCentralJERcentraljets,theWeight*w_kf*(1-scaleFacErr));  
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsSFErrMinus_01","",4,0,4,nCentralJERcentraljets,theWeight*w_kf*(1-scaleFacErr));  
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsSFErrPlus_"+sample,"",4,0,4,nCentralJERcentraljets,theWeight*w_kf*(1+scaleFacErr));  
 theHistograms.fill(std::string("ZZTo")+decay+"_CentralJetsSFErrPlus_01","",4,0,4,nCentralJERcentraljets,theWeight*w_kf*(1+scaleFacErr));  
 
 
 
 if(nCentralJERjets >=1){
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1SFErrSqMinus_"+sample,"",Xbins_ptjet1,centralPtJet1,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1SFErrSqMinus_01","",Xbins_ptjet1,centralPtJet1,theWeight*w_kf*(1-scaleFacErrSq));
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1SFErrSqPlus_"+sample,"",Xbins_ptjet1,centralPtJet1,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1SFErrSqPlus_01","",Xbins_ptjet1,centralPtJet1,theWeight*w_kf*(1+scaleFacErrSq));
   
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1SFErrSqMinus_"+sample,"",Xbins_etajet1,centralEtaJet1,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1SFErrSqMinus_01","",Xbins_etajet1,centralEtaJet1,theWeight*w_kf*(1-scaleFacErrSq));
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1SFErrSqPlus_"+sample,"",Xbins_etajet1,centralEtaJet1,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1SFErrSqPlus_01","",Xbins_etajet1,centralEtaJet1,theWeight*w_kf*(1+scaleFacErrSq));
   
   // Histograms for scale factor correlated errors
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1SFErrMinus_"+sample,"",Xbins_ptjet1,centralPtJet1,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1SFErrMinus_01","",Xbins_ptjet1,centralPtJet1,theWeight*w_kf*(1-scaleFacErr));
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1SFErrPlus_"+sample,"",Xbins_ptjet1,centralPtJet1,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1SFErrPlus_01","",Xbins_ptjet1,centralPtJet1,theWeight*w_kf*(1+scaleFacErr));
   
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1SFErrMinus_"+sample,"",Xbins_etajet1,centralEtaJet1,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1SFErrMinus_01","",Xbins_etajet1,centralEtaJet1,theWeight*w_kf*(1-scaleFacErr));
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1SFErrPlus_"+sample,"",Xbins_etajet1,centralEtaJet1,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1SFErrPlus_01","",Xbins_etajet1,centralEtaJet1,theWeight*w_kf*(1+scaleFacErr));
 } 
 
 if(nCentralJERjets>=2){  
   // cout <<  genCentralJERjets->at(0).motherId() << " " << genCentralJERjets->at(1).motherId() << endl;
   
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjSFErrSqMinus_"+sample,"",Xbins_mjj,centralMjj,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjSFErrSqMinus_01","",Xbins_mjj,centralMjj,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjSFErrSqPlus_"+sample,"",Xbins_mjj,centralMjj,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjSFErrSqPlus_01","",Xbins_mjj,centralMjj,theWeight*w_kf*(1+scaleFacErrSq));  
   
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaSFErrSqMinus_"+sample,"",Xbins_deta,centralDeta,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaSFErrSqMinus_01","",Xbins_deta,centralDeta,theWeight*w_kf*(1-scaleFacErrSq));
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaSFErrSqPlus_"+sample,"",Xbins_deta,centralDeta,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaSFErrSqPlus_01","",Xbins_deta,centralDeta,theWeight*w_kf*(1+scaleFacErrSq));
   
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2SFErrSqMinus_"+sample,"",Xbins_ptjet2,centralPtJet2,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2SFErrSqMinus_01","",Xbins_ptjet2,centralPtJet2,theWeight*w_kf*(1-scaleFacErrSq));
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2SFErrSqPlus_"+sample,"",Xbins_ptjet2,centralPtJet2,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2SFErrSqPlus_01","",Xbins_ptjet2,centralPtJet2,theWeight*w_kf*(1+scaleFacErrSq));
   
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2SFErrSqMinus_"+sample,"",Xbins_etajet2,centralEtaJet2,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2SFErrSqMinus_01","",Xbins_etajet2,centralEtaJet2,theWeight*w_kf*(1-scaleFacErrSq));
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2SFErrSqPlus_"+sample,"",Xbins_etajet2,centralEtaJet2,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2SFErrSqPlus_01","",Xbins_etajet2,centralEtaJet2,theWeight*w_kf*(1+scaleFacErrSq)); 
   
   // Histograms for scale factor correlated errors
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjSFErrMinus_"+sample,"",Xbins_mjj,centralMjj,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjSFErrMinus_01","",Xbins_mjj,centralMjj,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjSFErrPlus_"+sample,"",Xbins_mjj,centralMjj,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_MjjSFErrPlus_01","",Xbins_mjj,centralMjj,theWeight*w_kf*(1+scaleFacErr));  
   
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaSFErrMinus_"+sample,"",Xbins_deta,centralDeta,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaSFErrMinus_01","",Xbins_deta,centralDeta,theWeight*w_kf*(1-scaleFacErr));
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaSFErrPlus_"+sample,"",Xbins_deta,centralDeta,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_DetaSFErrPlus_01","",Xbins_deta,centralDeta,theWeight*w_kf*(1+scaleFacErr));
   
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2SFErrMinus_"+sample,"",Xbins_ptjet2,centralPtJet2,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2SFErrMinus_01","",Xbins_ptjet2,centralPtJet2,theWeight*w_kf*(1-scaleFacErr));
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2SFErrPlus_"+sample,"",Xbins_ptjet2,centralPtJet2,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2SFErrPlus_01","",Xbins_ptjet2,centralPtJet2,theWeight*w_kf*(1+scaleFacErr));
   
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2SFErrMinus_"+sample,"",Xbins_etajet2,centralEtaJet2,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2SFErrMinus_01","",Xbins_etajet2,centralEtaJet2,theWeight*w_kf*(1-scaleFacErr));
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2SFErrPlus_"+sample,"",Xbins_etajet2,centralEtaJet2,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2SFErrPlus_01","",Xbins_etajet2,centralEtaJet2,theWeight*w_kf*(1+scaleFacErr)); 
 }
 
 if(nCentralJERcentraljets>=2){  
   
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjSFErrSqMinus_"+sample,"",Xbins_mjj,centralMjj_cj,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjSFErrSqMinus_01","",Xbins_mjj,centralMjj_cj,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjSFErrSqPlus_"+sample,"",Xbins_mjj,centralMjj_cj,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjSFErrSqPlus_01","",Xbins_mjj,centralMjj_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
   
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaSFErrSqMinus_"+sample,"",Xbins_deta,centralDeta_cj,theWeight*w_kf*(1-scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaSFErrSqMinus_01","",Xbins_deta,centralDeta_cj,theWeight*w_kf*(1-scaleFacErrSq));
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaSFErrSqPlus_"+sample,"",Xbins_deta,centralDeta_cj,theWeight*w_kf*(1+scaleFacErrSq));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaSFErrSqPlus_01","",Xbins_deta,centralDeta_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
   
   // Histograms for scale factor correlated errors
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjSFErrMinus_"+sample,"",Xbins_mjj,centralMjj_cj,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjSFErrMinus_01","",Xbins_mjj,centralMjj_cj,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjSFErrPlus_"+sample,"",Xbins_mjj,centralMjj_cj,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjjSFErrPlus_01","",Xbins_mjj,centralMjj_cj,theWeight*w_kf*(1+scaleFacErr)); 
   
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaSFErrMinus_"+sample,"",Xbins_deta,centralDeta_cj,theWeight*w_kf*(1-scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaSFErrMinus_01","",Xbins_deta,centralDeta_cj,theWeight*w_kf*(1-scaleFacErr));
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaSFErrPlus_"+sample,"",Xbins_deta,centralDeta_cj,theWeight*w_kf*(1+scaleFacErr));  
   theHistograms.fill(std::string("ZZTo")+decay+"_CentralDetaSFErrPlus_01","",Xbins_deta,centralDeta_cj,theWeight*w_kf*(1+scaleFacErr));
 }
 
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // stable_sort(genJets->begin(), genJets->end(), PtComparator());
  //if MC gen (for response matrices only)    
  if (genCategory !=-1){
    if(topology.test(0)){  
      zz::SignalTopology zzSignalTopology = zz::getSignalTopology(*genParticles, *genJets);
      
      ngenjets =  genJets->size(); 
      if (ngenjets>3) ngenjets=3;
      
      ngencentraljets =  centralGenJets->size(); 
      if (ngencentraljets>3) ngencentraljets=3;
      
      m4L_gen = sqrt((genVBParticles->at(0).p4()+genVBParticles->at(1).p4())*(genVBParticles->at(0).p4()+genVBParticles->at(1).p4()));
      if (m4L_gen>=800) m4L_gen = 799;
      
      drzz_gen =physmath::deltaR(genVBParticles->at(0),genVBParticles->at(1));
      ptzz_gen =  (genVBParticles->at(0).p4()+genVBParticles->at(1).p4()).Pt();
      if (ptzz_gen>300) ptzz_gen=299;
      // dphizz_gen = fabs(genVBParticles->at(0).phi()-genVBParticles->at(1).phi());
      //if(dphizz_gen > TMath::Pi()) dphizz_gen = (TMath::TwoPi()) - dphizz_gen;
      dphizz_gen = fabs(physmath::deltaPhi(genVBParticles->at(0).phi(),genVBParticles->at(1).phi())); 
      //Response Matrix Mass (Reco&Gen) 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_"+sample,"", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_01","", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZ_"+sample,"", Xbins_drzz, Xbins_drzz,drzz ,drzz_gen , theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZ_01","", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf);
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZ_"+sample,"", Xbins_ptzz, Xbins_ptzz,ptzz ,ptzz_gen , theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZ_01","", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZ_"+sample,"", Xbins_dphizz, Xbins_dphizz,dphizz ,dphizz_gen , theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZ_01","", Xbins_dphizz, Xbins_dphizz, dphizz ,dphizz_gen , theWeight*w_kf);
      //********************************************************************RECO&GEN JETS*********************************************************************************
      
      //Response Matrix nJets (Reco&Gen) - No JER smearing
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_"+sample,"", 4,0,4,4,0,4, njets,ngenjets, theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_01","", 4,0,4,4,0,4, njets,ngenjets, theWeight*w_kf); 
      
      //Response Matrix nJets (Reco&Gen) - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERCentralSmear_"+sample,"", 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERCentralSmear_01","", 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_kf);
      
      // Response Matrix nJets (Reco&Gen) - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERUpSmear_"+sample,"", 4,0,4,4,0,4,nUpJERjets,ngenjets, theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERUpSmear_01","", 4,0,4,4,0,4,nUpJERjets,ngenjets, theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERDownSmear_"+sample,"", 4,0,4,4,0,4,nDownJERjets,ngenjets, theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERDownSmear_01","", 4,0,4,4,0,4,nDownJERjets,ngenjets, theWeight*w_kf);
      
      // Response Matrix nJets (Reco&Gen) - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JESUpSmear_"+sample,"", 4,0,4,4,0,4,nUpJESjets,ngenjets, theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JESUpSmear_01","", 4,0,4,4,0,4,nUpJESjets,ngenjets, theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JESDownSmear_"+sample,"", 4,0,4,4,0,4,nDownJESjets,ngenjets, theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JESDownSmear_01","",4,0,4,4,0,4,nDownJESjets,ngenjets, theWeight*w_kf);
      
      if(ngenjets >=1){
	ptjet1_gen = genJets->at(0).pt();
	if(ptjet1_gen >=300) ptjet1_gen =299;
	etajet1_gen =  fabs(genJets->at(0).eta());
	if(etajet1_gen >=4.7) etajet1_gen =4.6;
	
	if(nCentralJERjets>=1){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JERCentralSmear_"+sample,"",Xbins_ptjet1,Xbins_ptjet1, centralPtJet1,ptjet1_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JERCentralSmear_01","", Xbins_ptjet1,Xbins_ptjet1,centralPtJet1,ptjet1_gen,theWeight*w_kf);  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JERCentralSmear_"+sample,"",Xbins_etajet1,Xbins_etajet1, centralEtaJet1,etajet1_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JERCentralSmear_01","", Xbins_etajet1,Xbins_etajet1,centralEtaJet1,etajet1_gen,theWeight*w_kf);  
	}
	if(nUpJERjets>=1){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JERUpSmear_"+sample,"",Xbins_ptjet1,Xbins_ptjet1, upPtJet1JER,ptjet1_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JERUpSmear_01","", Xbins_ptjet1,Xbins_ptjet1,upPtJet1JER,ptjet1_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JERUpSmear_"+sample,"",Xbins_etajet1,Xbins_etajet1, upEtaJet1JER,etajet1_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JERUpSmear_01","", Xbins_etajet1,Xbins_etajet1,upEtaJet1JER,etajet1_gen,theWeight*w_kf);  
	}
	if(nDownJERjets>=1){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JERDownSmear_"+sample,"",Xbins_ptjet1,Xbins_ptjet1, downPtJet1JER,ptjet1_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JERDownSmear_01","", Xbins_ptjet1,Xbins_ptjet1,downPtJet1JER,ptjet1_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JERDownSmear_"+sample,"",Xbins_etajet1,Xbins_etajet1, downEtaJet1JER,etajet1_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JERDownSmear_01","", Xbins_etajet1,Xbins_etajet1,downEtaJet1JER,etajet1_gen,theWeight*w_kf);   
	}
	if(nUpJESjets>=1){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JESUpSmear_"+sample,"",Xbins_ptjet1,Xbins_ptjet1, upPtJet1JES,ptjet1_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JESUpSmear_01","", Xbins_ptjet1,Xbins_ptjet1,upPtJet1JES,ptjet1_gen,theWeight*w_kf);   
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JESUpSmear_"+sample,"",Xbins_etajet1,Xbins_etajet1, upEtaJet1JES,etajet1_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JESUpSmear_01","", Xbins_etajet1,Xbins_etajet1,upEtaJet1JES,etajet1_gen,theWeight*w_kf);  
	}
	if(nDownJESjets>=1){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JESDownSmear_"+sample,"",Xbins_ptjet1,Xbins_ptjet1, downPtJet1JES,ptjet1_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JESDownSmear_01","", Xbins_ptjet1,Xbins_ptjet1,downPtJet1JES,ptjet1_gen,theWeight*w_kf);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JESDownSmear_"+sample,"",Xbins_etajet1,Xbins_etajet1, downEtaJet1JES,etajet1_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JESDownSmear_01","", Xbins_etajet1,Xbins_etajet1,downEtaJet1JES,etajet1_gen,theWeight*w_kf);   
	}
      }//ngenjets>=1
      
      if(ngenjets>=2){  
	deta_gen = fabs(genJets->at(0).eta() - genJets->at(1).eta());
	mjj_gen =  (genJets->at(0).p4() + genJets->at(1).p4()).M();
	ptjet2_gen = genJets->at(1).pt(); 
	etajet2_gen =  fabs(genJets->at(1).eta());
	
	if(ptjet2_gen >=200) ptjet2_gen =199;
	if(etajet2_gen >=4.7) etajet2_gen =4.6;
	if (deta_gen>=4.7) deta_gen = 4.6;
	if (mjj_gen>=800) mjj_gen = 799;
	
	//Response Matrix  DeltaEta and mJJ (Reco&Gen) - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)
	if(nCentralJERjets>=2){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERCentralSmear_"+sample,"",Xbins_mjj,Xbins_mjj, centralMjj,mjj_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERCentralSmear_01","", Xbins_mjj,Xbins_mjj,centralMjj,mjj_gen,theWeight*w_kf);  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERCentralSmear_"+sample,"",Xbins_deta,Xbins_deta, centralDeta,deta_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERCentralSmear_01","", Xbins_deta,Xbins_deta,centralDeta,deta_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JERCentralSmear_"+sample,"",Xbins_ptjet2,Xbins_ptjet2, centralPtJet2,ptjet2_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JERCentralSmear_01","", Xbins_ptjet2,Xbins_ptjet2,centralPtJet2,ptjet2_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JERCentralSmear_"+sample,"",Xbins_etajet2,Xbins_etajet2, centralEtaJet2,etajet2_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JERCentralSmear_01","", Xbins_etajet2,Xbins_etajet2,centralEtaJet2,etajet2_gen,theWeight*w_kf);   
	}
	
	// Response Matrix  DeltaEta and mJJ (Reco&Gen) - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
	if(nUpJERjets>=2){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERUpSmear_"+sample,"",Xbins_mjj,Xbins_mjj, upMjjJER,mjj_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERUpSmear_01","", Xbins_mjj,Xbins_mjj,upMjjJER,mjj_gen,theWeight*w_kf);  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERUpSmear_"+sample,"",Xbins_deta,Xbins_deta, upDetaJER,deta_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERUpSmear_01","", Xbins_deta,Xbins_deta,upDetaJER,deta_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JERUpSmear_"+sample,"",Xbins_ptjet2,Xbins_ptjet2, upPtJet2JER,ptjet2_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JERUpSmear_01","", Xbins_ptjet2,Xbins_ptjet2,upPtJet2JER,ptjet2_gen,theWeight*w_kf);  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JERUpSmear_"+sample,"",Xbins_etajet2,Xbins_etajet2, upEtaJet2JER,etajet2_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JERUpSmear_01","", Xbins_etajet2,Xbins_etajet2,upEtaJet2JER,etajet2_gen,theWeight*w_kf); 
	}
	
	if(nDownJERjets>=2){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERDownSmear_"+sample,"",Xbins_mjj,Xbins_mjj, downMjjJER,mjj_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERDownSmear_01","", Xbins_mjj,Xbins_mjj,downMjjJER,mjj_gen,theWeight*w_kf);  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERDownSmear_"+sample,"",Xbins_deta,Xbins_deta, downDetaJER,deta_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERDownSmear_01","", Xbins_deta,Xbins_deta,downDetaJER,deta_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JERDownSmear_"+sample,"",Xbins_ptjet2,Xbins_ptjet2, downPtJet2JER,ptjet2_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JERDownSmear_01","", Xbins_ptjet2,Xbins_ptjet2,downPtJet2JER,ptjet2_gen,theWeight*w_kf);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JERDownSmear_"+sample,"",Xbins_etajet2,Xbins_etajet2, downEtaJet2JER,etajet2_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JERDownSmear_01","", Xbins_etajet2,Xbins_etajet2,downEtaJet2JER,etajet2_gen,theWeight*w_kf);  
	}
	
	// Response Matrix  DeltaEta and mJJ (Reco&Gen) - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
	if(nUpJESjets>=2){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JESUpSmear_"+sample,"",Xbins_mjj,Xbins_mjj, upMjjJES,mjj_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JESUpSmear_01","", Xbins_mjj,Xbins_mjj,upMjjJES,mjj_gen,theWeight*w_kf);  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JESUpSmear_"+sample,"",Xbins_deta,Xbins_deta, upDetaJES,deta_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JESUpSmear_01","", Xbins_deta,Xbins_deta,upDetaJES,deta_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JESUpSmear_"+sample,"",Xbins_ptjet2,Xbins_ptjet2, upPtJet2JES,ptjet2_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JESUpSmear_01","", Xbins_ptjet2,Xbins_ptjet2,upPtJet2JES,ptjet2_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JESUpSmear_"+sample,"",Xbins_etajet2,Xbins_etajet2, upEtaJet2JES,etajet2_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JESUpSmear_01","", Xbins_etajet2,Xbins_etajet2,upEtaJet2JES,etajet2_gen,theWeight*w_kf); 
	}
	
	if(nDownJESjets>=2){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JESDownSmear_"+sample,"",Xbins_mjj,Xbins_mjj, downMjjJES,mjj_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JESDownSmear_01","", Xbins_mjj,Xbins_mjj,downMjjJES,mjj_gen,theWeight*w_kf);  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JESDownSmear_"+sample,"",Xbins_deta,Xbins_deta, downDetaJES,deta_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JESDownSmear_01","", Xbins_deta,Xbins_deta,downDetaJES,deta_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JESDownSmear_"+sample,"",Xbins_ptjet2,Xbins_ptjet2, downPtJet2JES,ptjet2_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JESDownSmear_01","", Xbins_ptjet2,Xbins_ptjet2,downPtJet2JES,ptjet2_gen,theWeight*w_kf);
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JESDownSmear_"+sample,"",Xbins_etajet2,Xbins_etajet2, downEtaJet2JES,etajet2_gen,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JESDownSmear_01","", Xbins_etajet2,Xbins_etajet2,downEtaJet2JES,etajet2_gen,theWeight*w_kf); 
	}
      }//ngenjets>=2
      
      
      //***************************************************RECO&GEN CENTRAL JETS************************************************************************
      
      //Response Matrix nCentralJets (Reco&Gen) - No JER smearing
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_"+sample,"", 4,0,4,4,0,4, ncentraljets,ngencentraljets, theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_01","", 4,0,4,4,0,4, ncentraljets,ngencentraljets, theWeight*w_kf); 
      
      //Response Matrix nCentralJets (Reco&Gen) - JER smearing (CentralJets_JERCentralSmear to be used in the standard analysis)
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERCentralSmear_"+sample,"", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERCentralSmear_01","", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf);
      
      // Response Matrix nCentralJets (Reco&Gen) - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERUpSmear_"+sample,"", 4,0,4,4,0,4,nUpJERcentraljets,ngencentraljets, theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERUpSmear_01","", 4,0,4,4,0,4,nUpJERcentraljets,ngencentraljets, theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERDownSmear_"+sample,"", 4,0,4,4,0,4,nDownJERcentraljets,ngencentraljets, theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERDownSmear_01","", 4,0,4,4,0,4,nDownJERcentraljets,ngencentraljets, theWeight*w_kf);
      
      // Response Matrix nCentralJets (Reco&Gen) - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JESUpSmear_"+sample,"", 4,0,4,4,0,4,nUpJEScentraljets,ngencentraljets, theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JESUpSmear_01","", 4,0,4,4,0,4,nUpJEScentraljets,ngencentraljets, theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JESDownSmear_"+sample,"", 4,0,4,4,0,4,nDownJEScentraljets,ngencentraljets, theWeight*w_kf); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JESDownSmear_01","", 4,0,4,4,0,4,nDownJEScentraljets,ngencentraljets, theWeight*w_kf);
      
      if(ngencentraljets>=2){  
	deta_gen_cj = fabs(centralGenJets->at(0).eta() - centralGenJets->at(1).eta());
	mjj_gen_cj =  (centralGenJets->at(0).p4() + centralGenJets->at(1).p4()).M();
	
	if (deta_gen_cj>=4.7) deta_gen_cj = 4.6;
	if (mjj_gen_cj>=800) mjj_gen_cj = 799;
	
	//Response Matrix  DeltaEta and mJJ (Reco&Gen) - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)
	if(nCentralJERcentraljets>=2){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERCentralSmear_"+sample,"",Xbins_mjj,Xbins_mjj, centralMjj_cj,mjj_gen_cj,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERCentralSmear_01","", Xbins_mjj,Xbins_mjj,centralMjj_cj,mjj_gen_cj,theWeight*w_kf);  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERCentralSmear_"+sample,"",Xbins_deta,Xbins_deta, centralDeta_cj,deta_gen_cj,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERCentralSmear_01","", Xbins_deta,Xbins_deta,centralDeta_cj,deta_gen_cj,theWeight*w_kf); 
	  
	}
	
	// Response Matrix  DeltaEta and mJJ (Reco&Gen) - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
	if(nUpJERcentraljets>=2){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERUpSmear_"+sample,"",Xbins_mjj,Xbins_mjj, upMjjJER_cj,mjj_gen_cj,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERUpSmear_01","", Xbins_mjj,Xbins_mjj,upMjjJER_cj,mjj_gen_cj,theWeight*w_kf);  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERUpSmear_"+sample,"",Xbins_deta,Xbins_deta, upDetaJER_cj,deta_gen_cj,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERUpSmear_01","", Xbins_deta,Xbins_deta,upDetaJER_cj,deta_gen_cj,theWeight*w_kf); 
	}
	
	if(nDownJERcentraljets>=2){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERDownSmear_"+sample,"",Xbins_mjj,Xbins_mjj, downMjjJER_cj,mjj_gen_cj,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERDownSmear_01","", Xbins_mjj,Xbins_mjj,downMjjJER_cj,mjj_gen_cj,theWeight*w_kf);  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERDownSmear_"+sample,"",Xbins_deta,Xbins_deta, downDetaJER_cj,deta_gen_cj,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERDownSmear_01","", Xbins_deta,Xbins_deta,downDetaJER_cj,deta_gen_cj,theWeight*w_kf); 
	}
	
	// Response Matrix  DeltaEta and mJJ (Reco&Gen) - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
	if(nUpJEScentraljets>=2){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JESUpSmear_"+sample,"",Xbins_mjj,Xbins_mjj, upMjjJES_cj,mjj_gen_cj,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JESUpSmear_01","", Xbins_mjj,Xbins_mjj,upMjjJES_cj,mjj_gen_cj,theWeight*w_kf);  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JESUpSmear_"+sample,"",Xbins_deta,Xbins_deta, upDetaJES_cj,deta_gen_cj,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JESUpSmear_01","", Xbins_deta,Xbins_deta,upDetaJES_cj,deta_gen_cj,theWeight*w_kf); 
	}
	
	if(nDownJEScentraljets>=2){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JESDownSmear_"+sample,"",Xbins_mjj,Xbins_mjj, downMjjJES_cj,mjj_gen_cj,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JESDownSmear_01","", Xbins_mjj,Xbins_mjj,downMjjJES_cj,mjj_gen_cj,theWeight*w_kf);  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JESDownSmear_"+sample,"",Xbins_deta,Xbins_deta, downDetaJES_cj,deta_gen_cj,theWeight*w_kf); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JESDownSmear_01","", Xbins_deta,Xbins_deta,downDetaJES_cj,deta_gen_cj,theWeight*w_kf); 
	}
      }
      
      ////////////////////////////////////////////////SCALE FACTOR MATRICES///////////////////////////////////////////////////////
     
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MassSFErrSqMinus_"+sample,"", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf*(1-scaleFacErrSq)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MassSFErrSqMinus_01","", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf*(1-scaleFacErrSq));
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MassSFErrSqPlus_"+sample,"", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf*(1+scaleFacErrSq)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MassSFErrSqPlus_01","", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf*(1+scaleFacErrSq));
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZSFErrSqMinus_"+sample,"", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf*(1-scaleFacErrSq)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZSFErrSqMinus_01","", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf*(1-scaleFacErrSq));
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZSFErrSqPlus_"+sample,"", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf*(1+scaleFacErrSq)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZSFErrSqPlus_01","", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf*(1+scaleFacErrSq)); 

theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZSFErrSqMinus_"+sample,"", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf*(1-scaleFacErrSq)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZSFErrSqMinus_01","", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf*(1-scaleFacErrSq));
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZSFErrSqPlus_"+sample,"", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf*(1+scaleFacErrSq)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZSFErrSqPlus_01","", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf*(1+scaleFacErrSq));

      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZSFErrSqMinus_"+sample,"", Xbins_dphizz, Xbins_dphizz, dphizz ,dphizz_gen , theWeight*w_kf*(1-scaleFacErrSq)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZSFErrSqMinus_01","", Xbins_dphizz, Xbins_dphizz, dphizz ,dphizz_gen , theWeight*w_kf*(1-scaleFacErrSq));
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZSFErrSqPlus_"+sample,"", Xbins_dphizz, Xbins_dphizz, dphizz ,dphizz_gen , theWeight*w_kf*(1+scaleFacErrSq)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZSFErrSqPlus_01","", Xbins_dphizz, Xbins_dphizz, dphizz ,dphizz_gen , theWeight*w_kf*(1+scaleFacErrSq));

      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_JetsSFErrSqMinus_"+sample,"", 4,0,4,4,0,4,nCentralJERjets,ngenjets,theWeight*w_kf*(1-scaleFacErrSq)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_JetsSFErrSqMinus_01","", 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_kf*(1-scaleFacErrSq));
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_JetsSFErrSqPlus_"+sample,"", 4,0,4,4,0,4,nCentralJERjets,ngenjets,theWeight*w_kf*(1+scaleFacErrSq)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_JetsSFErrSqPlus_01","", 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_kf*(1+scaleFacErrSq));
 
      // Histograms for scale factor correlated errors
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MassSFErrMinus_"+sample,"", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf*(1-scaleFacErr)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MassSFErrMinus_01","", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf*(1-scaleFacErr));
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MassSFErrPlus_"+sample,"", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf*(1+scaleFacErr)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MassSFErrPlus_01","", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf*(1+scaleFacErr));
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZSFErrMinus_"+sample,"", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf*(1-scaleFacErr)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZSFErrMinus_01","", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf*(1-scaleFacErr));
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZSFErrPlus_"+sample,"", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf*(1+scaleFacErr)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZSFErrPlus_01","", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf*(1+scaleFacErr));
      
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZSFErrMinus_"+sample,"", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf*(1-scaleFacErr)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZSFErrMinus_01","", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf*(1-scaleFacErr));
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZSFErrPlus_"+sample,"", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf*(1+scaleFacErr)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZSFErrPlus_01","", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf*(1+scaleFacErr)); 

      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZSFErrMinus_"+sample,"", Xbins_dphizz, Xbins_dphizz, dphizz ,dphizz_gen , theWeight*w_kf*(1-scaleFacErr)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZSFErrMinus_01","", Xbins_dphizz, Xbins_dphizz, dphizz ,dphizz_gen , theWeight*w_kf*(1-scaleFacErr));
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZSFErrPlus_"+sample,"", Xbins_dphizz, Xbins_dphizz, dphizz ,dphizz_gen , theWeight*w_kf*(1+scaleFacErr)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZSFErrPlus_01","", Xbins_dphizz, Xbins_dphizz, dphizz ,drzz_gen , theWeight*w_kf*(1+scaleFacErr));
      
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_JetsSFErrMinus_"+sample,"", 4,0,4,4,0,4,nCentralJERjets,ngenjets,theWeight*w_kf*(1-scaleFacErr)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_JetsSFErrMinus_01","", 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_kf*(1-scaleFacErr));
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_JetsSFErrPlus_"+sample,"", 4,0,4,4,0,4,nCentralJERjets,ngenjets,theWeight*w_kf*(1+scaleFacErr)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_JetsSFErrPlus_01","", 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_kf*(1+scaleFacErr));

      if(ngenjets >=1){
	if(nCentralJERjets>=1){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1SFErrSqMinus_"+sample,"",Xbins_ptjet1,Xbins_ptjet1,centralPtJet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1SFErrSqMinus_01","", Xbins_ptjet1,Xbins_ptjet1,centralPtJet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErrSq));  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1SFErrSqPlus_"+sample,"",Xbins_ptjet1,Xbins_ptjet1, centralPtJet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1SFErrSqPlus_01","", Xbins_ptjet1,Xbins_ptjet1,centralPtJet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
	  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1SFErrSqMinus_"+sample,"",Xbins_etajet1,Xbins_etajet1, centralEtaJet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1SFErrSqMinus_01","", Xbins_etajet1,Xbins_etajet1,centralEtaJet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErrSq));   
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1SFErrSqPlus_"+sample,"",Xbins_etajet1,Xbins_etajet1, centralEtaJet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1SFErrSqPlus_01","", Xbins_etajet1,Xbins_etajet1,centralEtaJet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErrSq));   
	  // Histograms for scale factor correlated errors
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1SFErrMinus_"+sample,"",Xbins_ptjet1,Xbins_ptjet1,centralPtJet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1SFErrMinus_01","", Xbins_ptjet1,Xbins_ptjet1,centralPtJet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErr));  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1SFErrPlus_"+sample,"",Xbins_ptjet1,Xbins_ptjet1, centralPtJet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1SFErrPlus_01","", Xbins_ptjet1,Xbins_ptjet1,centralPtJet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErr)); 
	  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1SFErrMinus_"+sample,"",Xbins_etajet1,Xbins_etajet1, centralEtaJet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1SFErrMinus_01","", Xbins_etajet1,Xbins_etajet1,centralEtaJet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErr));   
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1SFErrPlus_"+sample,"",Xbins_etajet1,Xbins_etajet1, centralEtaJet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1SFErrPlus_01","", Xbins_etajet1,Xbins_etajet1,centralEtaJet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErr)); 
	}

      }//ngenjets>=1
      
      if(ngenjets>=2){ 
	
	//Response Matrix  DeltaEta and mJJ (Reco&Gen) - JER smearing (JetsSFErrSqMinus to be used in the standard analysis)
	if(nCentralJERjets>=2){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MjjSFErrSqMinus_"+sample,"",Xbins_mjj,Xbins_mjj, centralMjj,mjj_gen,theWeight*w_kf*(1-scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MjjSFErrSqMinus_01","", Xbins_mjj,Xbins_mjj,centralMjj,mjj_gen,theWeight*w_kf*(1-scaleFacErrSq));  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MjjSFErrSqPlus_"+sample,"",Xbins_mjj,Xbins_mjj, centralMjj,mjj_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MjjSFErrSqPlus_01","", Xbins_mjj,Xbins_mjj,centralMjj,mjj_gen,theWeight*w_kf*(1+scaleFacErrSq));   

	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DetaSFErrSqMinus_"+sample,"",Xbins_deta,Xbins_deta, centralDeta,deta_gen,theWeight*w_kf*(1-scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DetaSFErrSqMinus_01","", Xbins_deta,Xbins_deta,centralDeta,deta_gen,theWeight*w_kf*(1-scaleFacErrSq));  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DetaSFErrSqPlus_"+sample,"",Xbins_deta,Xbins_deta, centralDeta,deta_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DetaSFErrSqPlus_01","", Xbins_deta,Xbins_deta,centralDeta,deta_gen,theWeight*w_kf*(1+scaleFacErrSq)); 

	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2SFErrSqMinus_"+sample,"",Xbins_ptjet2,Xbins_ptjet2, centralPtJet2,ptjet2_gen,theWeight*w_kf*(1-scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2SFErrSqMinus_01","", Xbins_ptjet2,Xbins_ptjet2,centralPtJet2,ptjet2_gen,theWeight*w_kf*(1-scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2SFErrSqPlus_"+sample,"",Xbins_ptjet2,Xbins_ptjet2, centralPtJet2,ptjet2_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2SFErrSqPlus_01","", Xbins_ptjet2,Xbins_ptjet2,centralPtJet2,ptjet2_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2SFErrSqMinus_"+sample,"",Xbins_etajet2,Xbins_etajet2, centralEtaJet2,etajet2_gen,theWeight*w_kf*(1-scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2SFErrSqMinus_01","", Xbins_etajet2,Xbins_etajet2,centralEtaJet2,etajet2_gen,theWeight*w_kf*(1-scaleFacErrSq));  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2SFErrSqPlus_"+sample,"",Xbins_etajet2,Xbins_etajet2, centralEtaJet2,etajet2_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2SFErrSqPlus_01","", Xbins_etajet2,Xbins_etajet2,centralEtaJet2,etajet2_gen,theWeight*w_kf*(1+scaleFacErrSq)) ;
	  // Histograms for scale factor correlated errors 	  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MjjSFErrMinus_"+sample,"",Xbins_mjj,Xbins_mjj, centralMjj,mjj_gen,theWeight*w_kf*(1-scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MjjSFErrMinus_01","", Xbins_mjj,Xbins_mjj,centralMjj,mjj_gen,theWeight*w_kf*(1-scaleFacErr));  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MjjSFErrPlus_"+sample,"",Xbins_mjj,Xbins_mjj, centralMjj,mjj_gen,theWeight*w_kf*(1+scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MjjSFErrPlus_01","", Xbins_mjj,Xbins_mjj,centralMjj,mjj_gen,theWeight*w_kf*(1+scaleFacErr));   

	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DetaSFErrMinus_"+sample,"",Xbins_deta,Xbins_deta, centralDeta,deta_gen,theWeight*w_kf*(1-scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DetaSFErrMinus_01","", Xbins_deta,Xbins_deta,centralDeta,deta_gen,theWeight*w_kf*(1-scaleFacErr));  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DetaSFErrPlus_"+sample,"",Xbins_deta,Xbins_deta, centralDeta,deta_gen,theWeight*w_kf*(1+scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DetaSFErrPlus_01","", Xbins_deta,Xbins_deta,centralDeta,deta_gen,theWeight*w_kf*(1+scaleFacErr)); 

	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2SFErrMinus_"+sample,"",Xbins_ptjet2,Xbins_ptjet2, centralPtJet2,ptjet2_gen,theWeight*w_kf*(1-scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2SFErrMinus_01","", Xbins_ptjet2,Xbins_ptjet2,centralPtJet2,ptjet2_gen,theWeight*w_kf*(1-scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2SFErrPlus_"+sample,"",Xbins_ptjet2,Xbins_ptjet2, centralPtJet2,ptjet2_gen,theWeight*w_kf*(1+scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2SFErrPlus_01","", Xbins_ptjet2,Xbins_ptjet2,centralPtJet2,ptjet2_gen,theWeight*w_kf*(1+scaleFacErr)); 
 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2SFErrMinus_"+sample,"",Xbins_etajet2,Xbins_etajet2, centralEtaJet2,etajet2_gen,theWeight*w_kf*(1-scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2SFErrMinus_01","", Xbins_etajet2,Xbins_etajet2,centralEtaJet2,etajet2_gen,theWeight*w_kf*(1-scaleFacErr));  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2SFErrPlus_"+sample,"",Xbins_etajet2,Xbins_etajet2, centralEtaJet2,etajet2_gen,theWeight*w_kf*(1+scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2SFErrPlus_01","", Xbins_etajet2,Xbins_etajet2,centralEtaJet2,etajet2_gen,theWeight*w_kf*(1+scaleFacErr)) ;
	}
	
      }//ngenjets>=2
      
      
      //***************************************************RECO&GEN CENTRAL JETS************************************************************************
      
      //Response Matrix nCentralJets (Reco&Gen) - JER smearing (CentralJetsSFErrMinus to be used in the standard analysis)
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJetsSFErrSqMinus_"+sample,"", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf*(1-scaleFacErrSq)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJetsSFErrSqMinus_01","", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf*(1-scaleFacErrSq));
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJetsSFErrSqPlus_"+sample,"", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf*(1+scaleFacErrSq)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJetsSFErrSqPlus_01","", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf*(1+scaleFacErrSq));
     
      // Histograms for scale factor correlated errors
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJetsSFErrMinus_"+sample,"", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf*(1-scaleFacErr)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJetsSFErrMinus_01","", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf*(1-scaleFacErr));
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJetsSFErrPlus_"+sample,"", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf*(1+scaleFacErr)); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJetsSFErrPlus_01","", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf*(1+scaleFacErr));

      if(ngencentraljets>=2){  
	
	//Response Matrix  DeltaEta and mJJ (Reco&Gen) - JER smearing (JetsSFErrSqMinus to be used in the standard analysis)
	if(nCentralJERcentraljets>=2){
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjjSFErrSqMinus_"+sample,"",Xbins_mjj,Xbins_mjj, centralMjj_cj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjjSFErrSqMinus_01","", Xbins_mjj,Xbins_mjj,centralMjj_cj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjjSFErrSqPlus_"+sample,"",Xbins_mjj,Xbins_mjj, centralMjj_cj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjjSFErrSqPlus_01","", Xbins_mjj,Xbins_mjj,centralMjj_cj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
	  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDetaSFErrSqMinus_"+sample,"",Xbins_deta,Xbins_deta, centralDeta_cj,deta_gen_cj,theWeight*w_kf*(1-scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDetaSFErrSqMinus_01","", Xbins_deta,Xbins_deta,centralDeta_cj,deta_gen_cj,theWeight*w_kf*(1-scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDetaSFErrSqPlus_"+sample,"",Xbins_deta,Xbins_deta, centralDeta_cj,deta_gen_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDetaSFErrSqPlus_01","", Xbins_deta,Xbins_deta,centralDeta_cj,deta_gen_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
 
	  // Histograms for scale factor correlated errors 	  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjjSFErrMinus_"+sample,"",Xbins_mjj,Xbins_mjj, centralMjj_cj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjjSFErrMinus_01","", Xbins_mjj,Xbins_mjj,centralMjj_cj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjjSFErrPlus_"+sample,"",Xbins_mjj,Xbins_mjj, centralMjj_cj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjjSFErrPlus_01","", Xbins_mjj,Xbins_mjj,centralMjj_cj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErr)); 
	  
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDetaSFErrMinus_"+sample,"",Xbins_deta,Xbins_deta, centralDeta_cj,deta_gen_cj,theWeight*w_kf*(1-scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDetaSFErrMinus_01","", Xbins_deta,Xbins_deta,centralDeta_cj,deta_gen_cj,theWeight*w_kf*(1-scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDetaSFErrPlus_"+sample,"",Xbins_deta,Xbins_deta, centralDeta_cj,deta_gen_cj,theWeight*w_kf*(1+scaleFacErr)); 
	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDetaSFErrPlus_01","", Xbins_deta,Xbins_deta,centralDeta_cj,deta_gen_cj,theWeight*w_kf*(1+scaleFacErr));
	}
      }

      
      ///////////////////////////////////////////////////MATRICES IN FIDUCIAL REGION///////////////////////////////////////////////////////////////	
      // To run in the fiducial region (otherwise to be commented): 
      
      if(zz::inTightFiducialRegion(zzSignalTopology)){
	region = "_fr";
	inFiducialRegion ++;
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_"+sample+region,"", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_01"+region,"", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf);
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZ_"+sample+region,"", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZ_01"+region,"", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf);
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZ_"+sample+region,"", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZ_01"+region,"", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf);	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZ_"+sample+region,"", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf); 

	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZ_"+sample+region,"", Xbins_dphizz, Xbins_dphizz, dphizz ,dphizz_gen , theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZ_01"+region,"", Xbins_dphizz, Xbins_dphizz, dphizz ,dphizz_gen , theWeight*w_kf);
	//****************************************************************************RECO&GEN JETS****************************************************************************************
	
	//Response Matrix nJets (Reco&Gen) - No JER smearing
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_"+sample+region,"", 4,0,4,4,0,4, njets,ngenjets, theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_01"+region,"", 4,0,4,4,0,4, njets,ngenjets, theWeight*w_kf); 
	
	//Response Matrix nJets (Reco&Gen) - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERCentralSmear_"+sample+region,"", 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERCentralSmear_01"+region,"", 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_kf);
	
	// Response Matrix nJets (Reco&Gen) - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERUpSmear_"+sample+region,"", 4,0,4,4,0,4,nUpJERjets,ngenjets, theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERUpSmear_01"+region,"", 4,0,4,4,0,4,nUpJERjets,ngenjets, theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERDownSmear_"+sample+region,"", 4,0,4,4,0,4,nDownJERjets,ngenjets, theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERDownSmear_01"+region,"", 4,0,4,4,0,4,nDownJERjets,ngenjets, theWeight*w_kf);
	
	// Response Matrix nJets (Reco&Gen) - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JESUpSmear_"+sample+region,"", 4,0,4,4,0,4,nUpJESjets,ngenjets, theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JESUpSmear_01"+region,"", 4,0,4,4,0,4,nUpJESjets,ngenjets, theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JESDownSmear_"+sample+region,"", 4,0,4,4,0,4,nDownJESjets,ngenjets, theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JESDownSmear_01"+region,"", 4,0,4,4,0,4,nDownJESjets,ngenjets, theWeight*w_kf);
	
	if(ngenjets >=1){
	  ptjet1_gen = genJets->at(0).pt();
	  if(ptjet1_gen >=300) ptjet1_gen =299;
	  etajet1_gen = fabs(genJets->at(0).eta());
	  if(etajet1_gen >=4.7) etajet1_gen =4.6;
	  
	  if(nCentralJERjets>=1){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JERCentralSmear_"+sample+region,"",Xbins_ptjet1,Xbins_ptjet1, centralPtJet1,ptjet1_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JERCentralSmear_01"+region,"", Xbins_ptjet1,Xbins_ptjet1,centralPtJet1,ptjet1_gen,theWeight*w_kf);  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JERCentralSmear_"+sample+region,"",Xbins_etajet1,Xbins_etajet1, centralEtaJet1,etajet1_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JERCentralSmear_01"+region,"", Xbins_etajet1,Xbins_etajet1,centralEtaJet1,etajet1_gen,theWeight*w_kf);  
	  }
	  if(nUpJERjets>=1){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JERUpSmear_"+sample+region,"",Xbins_ptjet1,Xbins_ptjet1, upPtJet1JER,ptjet1_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JERUpSmear_01"+region,"", Xbins_ptjet1,Xbins_ptjet1,upPtJet1JER,ptjet1_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JERUpSmear_"+sample+region,"",Xbins_etajet1,Xbins_etajet1, upEtaJet1JER,etajet1_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JERUpSmear_01"+region,"", Xbins_etajet1,Xbins_etajet1,upEtaJet1JER,etajet1_gen,theWeight*w_kf);  
	  }
	  if(nDownJERjets>=1){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JERDownSmear_"+sample+region,"",Xbins_ptjet1,Xbins_ptjet1, downPtJet1JER,ptjet1_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JERDownSmear_01"+region,"", Xbins_ptjet1,Xbins_ptjet1,downPtJet1JER,ptjet1_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JERDownSmear_"+sample+region,"",Xbins_etajet1,Xbins_etajet1, downEtaJet1JER,etajet1_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JERDownSmear_01"+region,"", Xbins_etajet1,Xbins_etajet1,downEtaJet1JER,etajet1_gen,theWeight*w_kf);   
	  }
	  if(nUpJESjets>=1){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JESUpSmear_"+sample+region,"",Xbins_ptjet1,Xbins_ptjet1, upPtJet1JES,ptjet1_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JESUpSmear_01"+region,"", Xbins_ptjet1,Xbins_ptjet1,upPtJet1JES,ptjet1_gen,theWeight*w_kf);   
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JESUpSmear_"+sample+region,"",Xbins_etajet1,Xbins_etajet1, upEtaJet1JES,etajet1_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JESUpSmear_01"+region,"", Xbins_etajet1,Xbins_etajet1,upEtaJet1JES,etajet1_gen,theWeight*w_kf);  
	  }
	  if(nDownJESjets>=1){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JESDownSmear_"+sample+region,"",Xbins_ptjet1,Xbins_ptjet1, downPtJet1JES,ptjet1_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1_JESDownSmear_01"+region,"", Xbins_ptjet1,Xbins_ptjet1,downPtJet1JES,ptjet1_gen,theWeight*w_kf);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JESDownSmear_"+sample+region,"",Xbins_etajet1,Xbins_etajet1, downEtaJet1JES,etajet1_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1_JESDownSmear_01"+region,"", Xbins_etajet1,Xbins_etajet1,downEtaJet1JES,etajet1_gen,theWeight*w_kf);   
	  }
	}//ngenjets>=1
	
	if(ngenjets>=2){  
	  deta_gen = fabs(genJets->at(0).eta() - genJets->at(1).eta());
	  mjj_gen =  (genJets->at(0).p4() + genJets->at(1).p4()).M();
	  ptjet2_gen = genJets->at(1).pt(); 
	  etajet2_gen = fabs(genJets->at(1).eta());
	  
	  if(ptjet2_gen >=200) ptjet2_gen =199;
	  if(etajet2_gen >=4.7) etajet2_gen =4.6;
	  if (deta_gen>=4.7) deta_gen = 4.6;
	  if (mjj_gen>=800) mjj_gen = 799;
	  
	  //Response Matrix  DeltaEta and mJJ (Reco&Gen) - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)
	  if(nCentralJERjets>=2){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERCentralSmear_"+sample+region,"",Xbins_mjj,Xbins_mjj, centralMjj,mjj_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERCentralSmear_01"+region,"", Xbins_mjj,Xbins_mjj,centralMjj,mjj_gen,theWeight*w_kf);  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERCentralSmear_"+sample+region,"",Xbins_deta,Xbins_deta, centralDeta,deta_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERCentralSmear_01"+region,"", Xbins_deta,Xbins_deta,centralDeta,deta_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JERCentralSmear_"+sample+region,"",Xbins_ptjet2,Xbins_ptjet2, centralPtJet2,ptjet2_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JERCentralSmear_01"+region,"", Xbins_ptjet2,Xbins_ptjet2,centralPtJet2,ptjet2_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JERCentralSmear_"+sample+region,"",Xbins_etajet2,Xbins_etajet2, centralEtaJet2,etajet2_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JERCentralSmear_01"+region,"", Xbins_etajet2,Xbins_etajet2,centralEtaJet2,etajet2_gen,theWeight*w_kf);   
	  }
	  
	  // Response Matrix  DeltaEta and mJJ (Reco&Gen) - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
	  if(nUpJERjets>=2){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERUpSmear_"+sample+region,"",Xbins_mjj,Xbins_mjj, upMjjJER,mjj_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERUpSmear_01"+region,"", Xbins_mjj,Xbins_mjj,upMjjJER,mjj_gen,theWeight*w_kf);  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERUpSmear_"+sample+region,"",Xbins_deta,Xbins_deta, upDetaJER,deta_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERUpSmear_01"+region,"", Xbins_deta,Xbins_deta,upDetaJER,deta_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JERUpSmear_"+sample+region,"",Xbins_ptjet2,Xbins_ptjet2, upPtJet2JER,ptjet2_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JERUpSmear_01"+region,"", Xbins_ptjet2,Xbins_ptjet2,upPtJet2JER,ptjet2_gen,theWeight*w_kf);  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JERUpSmear_"+sample+region,"",Xbins_etajet2,Xbins_etajet2, upEtaJet2JER,etajet2_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JERUpSmear_01"+region,"", Xbins_etajet2,Xbins_etajet2,upEtaJet2JER,etajet2_gen,theWeight*w_kf); 
	  }
	  
	  if(nDownJERjets>=2){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERDownSmear_"+sample+region,"",Xbins_mjj,Xbins_mjj, downMjjJER,mjj_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERDownSmear_01"+region,"", Xbins_mjj,Xbins_mjj,downMjjJER,mjj_gen,theWeight*w_kf);  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERDownSmear_"+sample+region,"",Xbins_deta,Xbins_deta, downDetaJER,deta_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERDownSmear_01"+region,"", Xbins_deta,Xbins_deta,downDetaJER,deta_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JERDownSmear_"+sample+region,"",Xbins_ptjet2,Xbins_ptjet2, downPtJet2JER,ptjet2_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JERDownSmear_01"+region,"", Xbins_ptjet2,Xbins_ptjet2,downPtJet2JER,ptjet2_gen,theWeight*w_kf);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JERDownSmear_"+sample+region,"",Xbins_etajet2,Xbins_etajet2, downEtaJet2JER,etajet2_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JERDownSmear_01"+region,"", Xbins_etajet2,Xbins_etajet2,downEtaJet2JER,etajet2_gen,theWeight*w_kf);  
	  }
	  
	  // Response Matrix  DeltaEta and mJJ (Reco&Gen) - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
	  if(nUpJESjets>=2){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JESUpSmear_"+sample+region,"",Xbins_mjj,Xbins_mjj, upMjjJES,mjj_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JESUpSmear_01"+region,"", Xbins_mjj,Xbins_mjj,upMjjJES,mjj_gen,theWeight*w_kf);  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JESUpSmear_"+sample+region,"",Xbins_deta,Xbins_deta, upDetaJES,deta_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JESUpSmear_01"+region,"", Xbins_deta,Xbins_deta,upDetaJES,deta_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JESUpSmear_"+sample+region,"",Xbins_ptjet2,Xbins_ptjet2, upPtJet2JES,ptjet2_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JESUpSmear_01"+region,"", Xbins_ptjet2,Xbins_ptjet2,upPtJet2JES,ptjet2_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JESUpSmear_"+sample+region,"",Xbins_etajet2,Xbins_etajet2, upEtaJet2JES,etajet2_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JESUpSmear_01"+region,"", Xbins_etajet2,Xbins_etajet2,upEtaJet2JES,etajet2_gen,theWeight*w_kf); 
	  }
	  
	  if(nDownJESjets>=2){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JESDownSmear_"+sample+region,"",Xbins_mjj,Xbins_mjj, downMjjJES,mjj_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JESDownSmear_01"+region,"", Xbins_mjj,Xbins_mjj,downMjjJES,mjj_gen,theWeight*w_kf);  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JESDownSmear_"+sample+region,"",Xbins_deta,Xbins_deta, downDetaJES,deta_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JESDownSmear_01"+region,"", Xbins_deta,Xbins_deta,downDetaJES,deta_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JESDownSmear_"+sample+region,"",Xbins_ptjet2,Xbins_ptjet2, downPtJet2JES,ptjet2_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2_JESDownSmear_01"+region,"", Xbins_ptjet2,Xbins_ptjet2,downPtJet2JES,ptjet2_gen,theWeight*w_kf);
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JESDownSmear_"+sample+region,"",Xbins_etajet2,Xbins_etajet2, downEtaJet2JES,etajet2_gen,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2_JESDownSmear_01"+region,"", Xbins_etajet2,Xbins_etajet2,downEtaJet2JES,etajet2_gen,theWeight*w_kf); 
	  }
	}//ngenjets>=2
	
	
	//**************************************************RECO&GEN CENTRAL JETS******************************************************************************
	
	//Response Matrix nCentralJets (Reco&Gen) - No JER smearing
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_"+sample+region,"", 4,0,4,4,0,4, ncentraljets,ngencentraljets, theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_01"+region,"", 4,0,4,4,0,4, ncentraljets,ngencentraljets, theWeight*w_kf); 
	
	//Response Matrix nCentralJets (Reco&Gen) - JER smearing (CentralJets_JERCentralSmear to be used in the standard analysis)
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERCentralSmear_"+sample+region,"", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERCentralSmear_01"+region,"", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf);
	
	// Response Matrix nCentralJets (Reco&Gen) - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERUpSmear_"+sample+region,"", 4,0,4,4,0,4,nUpJERcentraljets,ngencentraljets, theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERUpSmear_01"+region,"", 4,0,4,4,0,4,nUpJERcentraljets,ngencentraljets, theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERDownSmear_"+sample+region,"", 4,0,4,4,0,4,nDownJERcentraljets,ngencentraljets, theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERDownSmear_01"+region,"", 4,0,4,4,0,4,nDownJERcentraljets,ngencentraljets, theWeight*w_kf);
	
	// Response Matrix nCentralJets (Reco&Gen) - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JESUpSmear_"+sample+region,"", 4,0,4,4,0,4,nUpJEScentraljets,ngencentraljets, theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JESUpSmear_01"+region,"", 4,0,4,4,0,4,nUpJEScentraljets,ngencentraljets, theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JESDownSmear_"+sample+region,"", 4,0,4,4,0,4,nDownJEScentraljets,ngencentraljets, theWeight*w_kf); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JESDownSmear_01"+region,"", 4,0,4,4,0,4,nDownJEScentraljets,ngencentraljets, theWeight*w_kf);
	
	if(ngencentraljets>=2){  
	  deta_gen_cj = fabs(centralGenJets->at(0).eta() - centralGenJets->at(1).eta());
	  mjj_gen_cj =  (centralGenJets->at(0).p4() + centralGenJets->at(1).p4()).M();
	  
	  if (deta_gen_cj>=4.7) deta_gen_cj = 4.6;
	  if (mjj_gen_cj>=800) mjj_gen_cj = 799;
	  
	  //Response Matrix  DeltaEta and mJJ (Reco&Gen) - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)
	  if(nCentralJERcentraljets>=2){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERCentralSmear_"+sample+region,"",Xbins_mjj,Xbins_mjj, centralMjj_cj,mjj_gen_cj,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERCentralSmear_01"+region,"", Xbins_mjj,Xbins_mjj,centralMjj_cj,mjj_gen_cj,theWeight*w_kf);  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERCentralSmear_"+sample+region,"",Xbins_deta,Xbins_deta, centralDeta_cj,deta_gen_cj,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERCentralSmear_01"+region,"", Xbins_deta,Xbins_deta,centralDeta_cj,deta_gen_cj,theWeight*w_kf); 
	    
	  }
	  
	  // Response Matrix  DeltaEta and mJJ (Reco&Gen) - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
	  if(nUpJERcentraljets>=2){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERUpSmear_"+sample+region,"",Xbins_mjj,Xbins_mjj, upMjjJER_cj,mjj_gen_cj,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERUpSmear_01"+region,"", Xbins_mjj,Xbins_mjj,upMjjJER_cj,mjj_gen_cj,theWeight*w_kf);  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERUpSmear_"+sample+region,"",Xbins_deta,Xbins_deta, upDetaJER_cj,deta_gen_cj,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERUpSmear_01"+region,"", Xbins_deta,Xbins_deta,upDetaJER_cj,deta_gen_cj,theWeight*w_kf); 
	  }
	  
	  if(nDownJERcentraljets>=2){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERDownSmear_"+sample+region,"",Xbins_mjj,Xbins_mjj, downMjjJER_cj,mjj_gen_cj,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERDownSmear_01"+region,"", Xbins_mjj,Xbins_mjj,downMjjJER_cj,mjj_gen_cj,theWeight*w_kf);  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERDownSmear_"+sample+region,"",Xbins_deta,Xbins_deta, downDetaJER_cj,deta_gen_cj,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERDownSmear_01"+region,"", Xbins_deta,Xbins_deta,downDetaJER_cj,deta_gen_cj,theWeight*w_kf); 
	  }
	  
	  // Response Matrix  DeltaEta and mJJ (Reco&Gen) - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
	  if(nUpJEScentraljets>=2){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JESUpSmear_"+sample+region,"",Xbins_mjj,Xbins_mjj, upMjjJES_cj,mjj_gen_cj,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JESUpSmear_01"+region,"", Xbins_mjj,Xbins_mjj,upMjjJES_cj,mjj_gen_cj,theWeight*w_kf);  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JESUpSmear_"+sample+region,"",Xbins_deta,Xbins_deta, upDetaJES_cj,deta_gen_cj,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JESUpSmear_01"+region,"", Xbins_deta,Xbins_deta,upDetaJES_cj,deta_gen_cj,theWeight*w_kf); 
	  }
	  
	  if(nDownJEScentraljets>=2){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JESDownSmear_"+sample+region,"",Xbins_mjj,Xbins_mjj, downMjjJES_cj,mjj_gen_cj,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JESDownSmear_01"+region,"", Xbins_mjj,Xbins_mjj,downMjjJES_cj,mjj_gen_cj,theWeight*w_kf);  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JESDownSmear_"+sample+region,"",Xbins_deta,Xbins_deta, downDetaJES_cj,deta_gen_cj,theWeight*w_kf); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JESDownSmear_01"+region,"", Xbins_deta,Xbins_deta,downDetaJES_cj,deta_gen_cj,theWeight*w_kf); 
	  }
	}
	
	////////////////////////////////////////////////SCALE FACTOR MATRICES///////////////////////////////////////////////////////
	
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MassSFErrSqMinus_"+sample+region,"", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf*(1-scaleFacErrSq)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MassSFErrSqMinus_01_fr","", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf*(1-scaleFacErrSq));
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MassSFErrSqPlus_"+sample+region,"", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf*(1+scaleFacErrSq)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MassSFErrSqPlus_01_fr","", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf*(1+scaleFacErrSq));

	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZSFErrSqMinus_"+sample+region,"", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf*(1-scaleFacErrSq)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZSFErrSqMinus_01_fr","", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf*(1-scaleFacErrSq));
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZSFErrSqPlus_"+sample+region,"", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf*(1+scaleFacErrSq)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZSFErrSqPlus_01_fr","", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf*(1+scaleFacErrSq));

	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZSFErrSqMinus_"+sample+region,"", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf*(1-scaleFacErrSq)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZSFErrSqMinus_01_fr","", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf*(1-scaleFacErrSq));
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZSFErrSqPlus_"+sample+region,"", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf*(1+scaleFacErrSq)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZSFErrSqPlus_01_fr","", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf*(1+scaleFacErrSq));

	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZSFErrSqMinus_"+sample+region,"", Xbins_ptzz, Xbins_ptzz, dphizz ,dphizz_gen , theWeight*w_kf*(1-scaleFacErrSq)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZSFErrSqMinus_01_fr","", Xbins_dphizz, Xbins_dphizz, dphizz ,dphizz_gen , theWeight*w_kf*(1-scaleFacErrSq));
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZSFErrSqPlus_"+sample+region,"", Xbins_dphizz, Xbins_dphizz, dphizz ,dphizz_gen , theWeight*w_kf*(1+scaleFacErrSq)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZSFErrSqPlus_01_fr","", Xbins_dphizz, Xbins_dphizz, dphizz ,dphizz_gen , theWeight*w_kf*(1+scaleFacErrSq));
	
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_JetsSFErrSqMinus_"+sample+region,"", 4,0,4,4,0,4,nCentralJERjets,ngenjets,theWeight*w_kf*(1-scaleFacErrSq)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_JetsSFErrSqMinus_01_fr","", 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_kf*(1-scaleFacErrSq));
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_JetsSFErrSqPlus_"+sample+region,"", 4,0,4,4,0,4,nCentralJERjets,ngenjets,theWeight*w_kf*(1+scaleFacErrSq)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_JetsSFErrSqPlus_01_fr","", 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_kf*(1+scaleFacErrSq));
	
	// Histograms for scale factor correlated errors
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MassSFErrMinus_"+sample+region,"", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf*(1-scaleFacErr)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MassSFErrMinus_01_fr","", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf*(1-scaleFacErr));
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MassSFErrPlus_"+sample+region,"", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf*(1+scaleFacErr)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MassSFErrPlus_01_fr","", Xbins, Xbins, m4L ,m4L_gen , theWeight*w_kf*(1+scaleFacErr));
	
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZSFErrMinus_"+sample+region,"", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf*(1-scaleFacErr)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZSFErrMinus_01_fr","", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf*(1-scaleFacErr));
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZSFErrPlus_"+sample+region,"", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf*(1+scaleFacErr)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_dRZZSFErrPlus_01_fr","", Xbins_drzz, Xbins_drzz, drzz ,drzz_gen , theWeight*w_kf*(1+scaleFacErr));
	
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZSFErrMinus_"+sample+region,"", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf*(1-scaleFacErr)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZSFErrMinus_01_fr","", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf*(1-scaleFacErr));
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZSFErrPlus_"+sample+region,"", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf*(1+scaleFacErr)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtZZSFErrPlus_01_fr","", Xbins_ptzz, Xbins_ptzz, ptzz ,ptzz_gen , theWeight*w_kf*(1+scaleFacErr));
	
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZSFErrMinus_"+sample+region,"", Xbins_ptzz, Xbins_dphizz, dphizz ,dphizz_gen , theWeight*w_kf*(1-scaleFacErr)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZSFErrMinus_01_fr","", Xbins_dphizz, Xbins_dphizz, dphizz ,dphizz_gen , theWeight*w_kf*(1-scaleFacErr));
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZSFErrPlus_"+sample+region,"", Xbins_dphizz, Xbins_dphizz, dphizz ,dphizz_gen , theWeight*w_kf*(1+scaleFacErr)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DphiZZSFErrPlus_01_fr","", Xbins_dphizz, Xbins_dphizz, dphizz ,dphizz_gen , theWeight*w_kf*(1+scaleFacErr));

	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_JetsSFErrMinus_"+sample+region,"", 4,0,4,4,0,4,nCentralJERjets,ngenjets,theWeight*w_kf*(1-scaleFacErr)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_JetsSFErrMinus_01_fr","", 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_kf*(1-scaleFacErr));
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_JetsSFErrPlus_"+sample+region,"", 4,0,4,4,0,4,nCentralJERjets,ngenjets,theWeight*w_kf*(1+scaleFacErr)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_JetsSFErrPlus_01_fr","", 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_kf*(1+scaleFacErr));
	
	if(ngenjets >=1){
	  if(nCentralJERjets>=1){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1SFErrSqMinus_"+sample+region,"",Xbins_ptjet1,Xbins_ptjet1,centralPtJet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1SFErrSqMinus_01_fr","", Xbins_ptjet1,Xbins_ptjet1,centralPtJet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErrSq));  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1SFErrSqPlus_"+sample+region,"",Xbins_ptjet1,Xbins_ptjet1, centralPtJet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1SFErrSqPlus_01_fr","", Xbins_ptjet1,Xbins_ptjet1,centralPtJet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
	    
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1SFErrSqMinus_"+sample+region,"",Xbins_etajet1,Xbins_etajet1, centralEtaJet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1SFErrSqMinus_01_fr","", Xbins_etajet1,Xbins_etajet1,centralEtaJet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErrSq));   
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1SFErrSqPlus_"+sample+region,"",Xbins_etajet1,Xbins_etajet1, centralEtaJet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1SFErrSqPlus_01_fr","", Xbins_etajet1,Xbins_etajet1,centralEtaJet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErrSq));   
	    
	    // Histograms for scale factor correlated errors
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1SFErrMinus_"+sample+region,"",Xbins_ptjet1,Xbins_ptjet1,centralPtJet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1SFErrMinus_01_fr","", Xbins_ptjet1,Xbins_ptjet1,centralPtJet1,ptjet1_gen,theWeight*w_kf*(1-scaleFacErr));  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1SFErrPlus_"+sample+region,"",Xbins_ptjet1,Xbins_ptjet1, centralPtJet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet1SFErrPlus_01_fr","", Xbins_ptjet1,Xbins_ptjet1,centralPtJet1,ptjet1_gen,theWeight*w_kf*(1+scaleFacErr)); 
	    
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1SFErrMinus_"+sample+region,"",Xbins_etajet1,Xbins_etajet1, centralEtaJet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1SFErrMinus_01_fr","", Xbins_etajet1,Xbins_etajet1,centralEtaJet1,etajet1_gen,theWeight*w_kf*(1-scaleFacErr));   
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1SFErrPlus_"+sample+region,"",Xbins_etajet1,Xbins_etajet1, centralEtaJet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet1SFErrPlus_01_fr","", Xbins_etajet1,Xbins_etajet1,centralEtaJet1,etajet1_gen,theWeight*w_kf*(1+scaleFacErr)); 
	  }
	  
	}//ngenjets>=1
	
	if(ngenjets>=2){ 
	  
	  //Response Matrix  DeltaEta and mJJ (Reco&Gen) - JER smearing (JetsSFErrSqMinus to be used in the standard analysis)
	  if(nCentralJERjets>=2){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MjjSFErrSqMinus_"+sample+region,"",Xbins_mjj,Xbins_mjj, centralMjj,mjj_gen,theWeight*w_kf*(1-scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MjjSFErrSqMinus_01_fr","", Xbins_mjj,Xbins_mjj,centralMjj,mjj_gen,theWeight*w_kf*(1-scaleFacErrSq));  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MjjSFErrSqPlus_"+sample+region,"",Xbins_mjj,Xbins_mjj, centralMjj,mjj_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MjjSFErrSqPlus_01_fr","", Xbins_mjj,Xbins_mjj,centralMjj,mjj_gen,theWeight*w_kf*(1+scaleFacErrSq));   
	    
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DetaSFErrSqMinus_"+sample+region,"",Xbins_deta,Xbins_deta, centralDeta,deta_gen,theWeight*w_kf*(1-scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DetaSFErrSqMinus_01_fr","", Xbins_deta,Xbins_deta,centralDeta,deta_gen,theWeight*w_kf*(1-scaleFacErrSq));  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DetaSFErrSqPlus_"+sample+region,"",Xbins_deta,Xbins_deta, centralDeta,deta_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DetaSFErrSqPlus_01_fr","", Xbins_deta,Xbins_deta,centralDeta,deta_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
	    
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2SFErrSqMinus_"+sample+region,"",Xbins_ptjet2,Xbins_ptjet2, centralPtJet2,ptjet2_gen,theWeight*w_kf*(1-scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2SFErrSqMinus_01_fr","", Xbins_ptjet2,Xbins_ptjet2,centralPtJet2,ptjet2_gen,theWeight*w_kf*(1-scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2SFErrSqPlus_"+sample+region,"",Xbins_ptjet2,Xbins_ptjet2, centralPtJet2,ptjet2_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2SFErrSqPlus_01_fr","", Xbins_ptjet2,Xbins_ptjet2,centralPtJet2,ptjet2_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
	    
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2SFErrSqMinus_"+sample+region,"",Xbins_etajet2,Xbins_etajet2, centralEtaJet2,etajet2_gen,theWeight*w_kf*(1-scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2SFErrSqMinus_01_fr","", Xbins_etajet2,Xbins_etajet2,centralEtaJet2,etajet2_gen,theWeight*w_kf*(1-scaleFacErrSq));  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2SFErrSqPlus_"+sample+region,"",Xbins_etajet2,Xbins_etajet2, centralEtaJet2,etajet2_gen,theWeight*w_kf*(1+scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2SFErrSqPlus_01_fr","", Xbins_etajet2,Xbins_etajet2,centralEtaJet2,etajet2_gen,theWeight*w_kf*(1+scaleFacErrSq)) ;
	    
	    // Histograms for scale factor correlated errors 	  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MjjSFErrMinus_"+sample+region,"",Xbins_mjj,Xbins_mjj, centralMjj,mjj_gen,theWeight*w_kf*(1-scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MjjSFErrMinus_01_fr","", Xbins_mjj,Xbins_mjj,centralMjj,mjj_gen,theWeight*w_kf*(1-scaleFacErr));  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MjjSFErrPlus_"+sample+region,"",Xbins_mjj,Xbins_mjj, centralMjj,mjj_gen,theWeight*w_kf*(1+scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_MjjSFErrPlus_01_fr","", Xbins_mjj,Xbins_mjj,centralMjj,mjj_gen,theWeight*w_kf*(1+scaleFacErr));   
	    
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DetaSFErrMinus_"+sample+region,"",Xbins_deta,Xbins_deta, centralDeta,deta_gen,theWeight*w_kf*(1-scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DetaSFErrMinus_01_fr","", Xbins_deta,Xbins_deta,centralDeta,deta_gen,theWeight*w_kf*(1-scaleFacErr));  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DetaSFErrPlus_"+sample+region,"",Xbins_deta,Xbins_deta, centralDeta,deta_gen,theWeight*w_kf*(1+scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_DetaSFErrPlus_01_fr","", Xbins_deta,Xbins_deta,centralDeta,deta_gen,theWeight*w_kf*(1+scaleFacErr)); 
	    
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2SFErrMinus_"+sample+region,"",Xbins_ptjet2,Xbins_ptjet2, centralPtJet2,ptjet2_gen,theWeight*w_kf*(1-scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2SFErrMinus_01_fr","", Xbins_ptjet2,Xbins_ptjet2,centralPtJet2,ptjet2_gen,theWeight*w_kf*(1-scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2SFErrPlus_"+sample+region,"",Xbins_ptjet2,Xbins_ptjet2, centralPtJet2,ptjet2_gen,theWeight*w_kf*(1+scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_PtJet2SFErrPlus_01_fr","", Xbins_ptjet2,Xbins_ptjet2,centralPtJet2,ptjet2_gen,theWeight*w_kf*(1+scaleFacErr)); 
	    
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2SFErrMinus_"+sample+region,"",Xbins_etajet2,Xbins_etajet2, centralEtaJet2,etajet2_gen,theWeight*w_kf*(1-scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2SFErrMinus_01_fr","", Xbins_etajet2,Xbins_etajet2,centralEtaJet2,etajet2_gen,theWeight*w_kf*(1-scaleFacErr));  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2SFErrPlus_"+sample+region,"",Xbins_etajet2,Xbins_etajet2, centralEtaJet2,etajet2_gen,theWeight*w_kf*(1+scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_EtaJet2SFErrPlus_01_fr","", Xbins_etajet2,Xbins_etajet2,centralEtaJet2,etajet2_gen,theWeight*w_kf*(1+scaleFacErr)) ;
	  }
	  
	}//ngenjets>=2
	
	
	//***************************************************RECO&GEN CENTRAL JETS************************************************************************
	
	//Response Matrix nCentralJets (Reco&Gen) - JER smearing (CentralJetsSFErrMinus to be used in the standard analysis)
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJetsSFErrSqMinus_"+sample+region,"", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf*(1-scaleFacErrSq)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJetsSFErrSqMinus_01_fr","", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf*(1-scaleFacErrSq));
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJetsSFErrSqPlus_"+sample+region,"", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf*(1+scaleFacErrSq)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJetsSFErrSqPlus_01_fr","", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf*(1+scaleFacErrSq));
	
	// Histograms for scale factor correlated errors
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJetsSFErrMinus_"+sample+region,"", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf*(1-scaleFacErr)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJetsSFErrMinus_01_fr","", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf*(1-scaleFacErr));
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJetsSFErrPlus_"+sample+region,"", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf*(1+scaleFacErr)); 
	theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJetsSFErrPlus_01_fr","", 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight*w_kf*(1+scaleFacErr));
	
	if(ngencentraljets>=2){  
	  
	  //Response Matrix  DeltaEta and mJJ (Reco&Gen) - JER smearing (JetsSFErrSqMinus to be used in the standard analysis)
	  if(nCentralJERcentraljets>=2){
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjjSFErrSqMinus_"+sample+region,"",Xbins_mjj,Xbins_mjj, centralMjj_cj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjjSFErrSqMinus_01_fr","", Xbins_mjj,Xbins_mjj,centralMjj_cj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjjSFErrSqPlus_"+sample+region,"",Xbins_mjj,Xbins_mjj, centralMjj_cj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjjSFErrSqPlus_01_fr","", Xbins_mjj,Xbins_mjj,centralMjj_cj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
	    
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDetaSFErrSqMinus_"+sample+region,"",Xbins_deta,Xbins_deta, centralDeta_cj,deta_gen_cj,theWeight*w_kf*(1-scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDetaSFErrSqMinus_01_fr","", Xbins_deta,Xbins_deta,centralDeta_cj,deta_gen_cj,theWeight*w_kf*(1-scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDetaSFErrSqPlus_"+sample+region,"",Xbins_deta,Xbins_deta, centralDeta_cj,deta_gen_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDetaSFErrSqPlus_01_fr","", Xbins_deta,Xbins_deta,centralDeta_cj,deta_gen_cj,theWeight*w_kf*(1+scaleFacErrSq)); 
	    
	    // Histograms for scale factor correlated errors 	  
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjjSFErrMinus_"+sample+region,"",Xbins_mjj,Xbins_mjj, centralMjj_cj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjjSFErrMinus_01_fr","", Xbins_mjj,Xbins_mjj,centralMjj_cj,mjj_gen_cj,theWeight*w_kf*(1-scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjjSFErrPlus_"+sample+region,"",Xbins_mjj,Xbins_mjj, centralMjj_cj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjjSFErrPlus_01_fr","", Xbins_mjj,Xbins_mjj,centralMjj_cj,mjj_gen_cj,theWeight*w_kf*(1+scaleFacErr)); 
	    
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDetaSFErrMinus_"+sample+region,"",Xbins_deta,Xbins_deta, centralDeta_cj,deta_gen_cj,theWeight*w_kf*(1-scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDetaSFErrMinus_01_fr","", Xbins_deta,Xbins_deta,centralDeta_cj,deta_gen_cj,theWeight*w_kf*(1-scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDetaSFErrPlus_"+sample+region,"",Xbins_deta,Xbins_deta, centralDeta_cj,deta_gen_cj,theWeight*w_kf*(1+scaleFacErr)); 
	    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDetaSFErrPlus_01_fr","", Xbins_deta,Xbins_deta,centralDeta_cj,deta_gen_cj,theWeight*w_kf*(1+scaleFacErr));
	  }
	}
	
      }
      
      //*************************************************************************************************************************************************************************************************
      
    } //end if(topology.test(0))
    
  }//end if(genCategory != -1)  
  
  theHistograms.fill(std::string("ZZTo")+decay+"_MinDeltaR_01", "",100,0,10, min_dR, theWeight*w_kf);
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if(region_ == phys::CR3P1F || region_ == phys::CR2P2F) {
    theHistograms.fill(std::string("ZZTo")+decay+"_Mass"+"_FRVar","", Xbins, m4L,ZZ->fakeRateSFVar()); 
    theHistograms.fill(std::string("ZZTo")+decay+"_dRZZ"+"_FRVar","", Xbins_drzz, drzz,ZZ->fakeRateSFVar());
    theHistograms.fill(std::string("ZZTo")+decay+"_PtZZ"+"_FRVar","", Xbins_ptzz, ptzz,ZZ->fakeRateSFVar());
    theHistograms.fill(std::string("ZZTo")+decay+"_DphiZZ"+"_FRVar","", Xbins_dphizz, dphizz,ZZ->fakeRateSFVar());
    theHistograms.fill(std::string("ZZTo")+decay+"_Jets"+"_FRVar", "",4,0,4,njets,ZZ->fakeRateSFVar());
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets"+"_FRVar","",4,0,4,ncentraljets,ZZ->fakeRateSFVar());
    if(njets>=1){ 
      theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1"+"_FRVar", "",Xbins_ptjet1,ptJet1,ZZ->fakeRateSFVar());
      theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1"+"_FRVar", "",Xbins_etajet1,etaJet1,ZZ->fakeRateSFVar());
    } 
    if(njets>=2){ 
      theHistograms.fill(std::string("ZZTo")+decay+"_Deta"+"_FRVar", "",Xbins_deta,deta,ZZ->fakeRateSFVar());
      theHistograms.fill(std::string("ZZTo")+decay+"_Mjj"+"_FRVar", "",Xbins_mjj,mjj,ZZ->fakeRateSFVar());
      theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2"+"_FRVar", "",Xbins_ptjet2,ptJet2,ZZ->fakeRateSFVar());
      theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2"+"_FRVar", "",Xbins_etajet2,etaJet2,ZZ->fakeRateSFVar());
    }  
    if(njets>=2){ 
      theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta"+"_FRVar", "",Xbins_deta,centraldeta,ZZ->fakeRateSFVar());
      theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj"+"_FRVar", "",Xbins_mjj,centralmjj,ZZ->fakeRateSFVar());
    }
  }
}// end ZZplots()

//Smearing function
double ZZRecoAnalyzer::JER_PtSmear(double pt, double width)
{
  double ptsmear= gRandom->Gaus(pt,width);    
  return ptsmear;
}

void ZZRecoAnalyzer::analyze(){
  
  e++;
  
  // ///////////////////////////////////////////////DATA ONLY/////////////////////////////////////////////////
  // //for data only (JES uncertainty) 
  
  UpJESData_jets->clear();
  DownJESData_jets->clear();
  UpJESData_centraljets->clear();
  DownJESData_centraljets->clear();
  
  UpJESData_jetPt->clear();
  DownJESData_jetPt->clear();
  
  foreach(const phys::Jet &dataJet, *pjets){
    double dataJetPt = 0;
    dataJetPt = dataJet.pt();
    
    //JES correction: Up and down velues used to assess systematic uncertainty on jet energy resolution
    double newJetPtJESData_up =0;  
    double newJetPtJESData_down =0;
    
    newJetPtJESData_up = dataJetPt*(1+dataJet.uncOnFourVectorScale());
    newJetPtJESData_down = dataJetPt*(1-dataJet.uncOnFourVectorScale());
    
    if(newJetPtJESData_up > 30) {
      UpJESData_jets->push_back(dataJet); 
      UpJESData_jetPt->push_back(newJetPtJESData_up);
    }
    if(newJetPtJESData_down > 30) {
      DownJESData_jets->push_back(dataJet);
      DownJESData_jetPt->push_back(newJetPtJESData_down);
    }    
    if(newJetPtJESData_up > 30 && fabs(dataJet.eta())<2.4) UpJESData_centraljets->push_back(dataJet); 
    if(newJetPtJESData_down > 30 && fabs(dataJet.eta())<2.4) DownJESData_centraljets->push_back(dataJet);

    //cout << decay << " " << dataJetPt << " " <<   newJetPtJESData_up << " " << newJetPtJESData_down << " " << njets << endl;
    
  }
  
  nUpJESDatajets = UpJESData_jets->size();
  nDownJESDatajets = DownJESData_jets->size();
  nUpJESDatacentraljets = UpJESData_centraljets->size();
  nDownJESDatacentraljets = DownJESData_centraljets->size();
  if (nUpJESDatajets>3) nUpJESDatajets=3;
  if (nDownJESDatajets>3) nDownJESDatajets=3;
  if (nUpJESDatacentraljets>3) nUpJESDatacentraljets=3;
  if (nDownJESDatacentraljets>3) nDownJESDatacentraljets=3;
  
  stable_sort(UpJESData_jets->begin(), UpJESData_jets->end(), PtComparator());
  stable_sort(DownJESData_jets->begin(), DownJESData_jets->end(), PtComparator());
  stable_sort(UpJESData_centraljets->begin(), UpJESData_centraljets->end(), PtComparator());
  stable_sort(DownJESData_centraljets->begin(), DownJESData_centraljets->end(), PtComparator());
  
  stable_sort(UpJESData_jetPt->begin(), UpJESData_jetPt->end(),std::greater<float>());
  stable_sort(DownJESData_jetPt->begin(), DownJESData_jetPt->end(),std::greater<float>());
  
  ///////////////////////////////////////////////// MC ONLY/////////////////////////////////////////////////
//Smearing jet pt (Up and Down distributions just for systematic uncertainty estimate, Central for the standard analysis) without requiring it is signal (1D distributions made of ALL reco events)
    CentralJER_jets->clear();
    UpJER_jets->clear();
    DownJER_jets->clear();
    CentralJER_centraljets->clear();
    UpJER_centraljets->clear();
    DownJER_centraljets->clear();

    UpJES_jets->clear();
    DownJES_jets->clear();
    UpJES_centraljets->clear();
    DownJES_centraljets->clear();

    CentralJER_jetPt->clear();
    UpJER_jetPt->clear();
    DownJER_jetPt->clear();
    UpJES_jetPt->clear();
    DownJES_jetPt->clear();
       
    //Loop on all reco jets
    foreach(const phys::Jet &jet, *pjets){
      //jet_count ++;      

      double jetPt = 0;
      jetPt = jet.pt();
      
      //JER correction (applied only on MC reco). Up and down velues used to assess systematic uncertainty on jet energy resolution
      double newJetPtJER =0; 
      double newJetPtJER_up =0;  
      double newJetPtJER_down =0;
      double width = 0;
      double width_up = 0; 
      double width_down = 0;
      
      width = jet.jer_width(phys::Jet::central);
      width_up = jet.jer_width(phys::Jet::up); 
      width_down = jet.jer_width(phys::Jet::down);
      
      newJetPtJER = JER_PtSmear(jetPt, width);
      newJetPtJER_up = JER_PtSmear(jetPt, width_up);  
      newJetPtJER_down = JER_PtSmear(jetPt, width_down);
      
      if(newJetPtJER > 30){
	CentralJER_jets->push_back(jet);
	CentralJER_jetPt->push_back(newJetPtJER);
      }
      if(newJetPtJER_up > 30) {
	UpJER_jets->push_back(jet); 
      	UpJER_jetPt->push_back(newJetPtJER_up);
      }
      if(newJetPtJER_down > 30){
	DownJER_jets->push_back(jet);
	DownJER_jetPt->push_back(newJetPtJER_down);
      }
      if(newJetPtJER > 30 && fabs(jet.eta())<2.4) CentralJER_centraljets->push_back(jet);
      if(newJetPtJER_up > 30 && fabs(jet.eta())<2.4) UpJER_centraljets->push_back(jet); 
      if(newJetPtJER_down > 30 && fabs(jet.eta())<2.4) DownJER_centraljets->push_back(jet);
      
      //JES correction: Up and down velues used to assess systematic uncertainty on jet energy resolution
      double newJetPtJES_up =0;  
      double newJetPtJES_down =0;
      
      //cout << jet.uncOnFourVectorScale() << endl;
      
      newJetPtJES_up = newJetPtJER*(1+jet.uncOnFourVectorScale());
      newJetPtJES_down = newJetPtJER*(1-jet.uncOnFourVectorScale());
      
      if(newJetPtJES_up > 30) {
	UpJES_jets->push_back(jet);
	UpJES_jetPt->push_back(newJetPtJES_up);
      } 
      if(newJetPtJES_down > 30) {
	DownJES_jets->push_back(jet);
	DownJES_jetPt->push_back(newJetPtJES_down);
      }
      if(newJetPtJES_up > 30 && fabs(jet.eta())<2.4 ) UpJES_centraljets->push_back(jet); 
      if(newJetPtJES_down > 30 &&  fabs(jet.eta())<2.4) DownJES_centraljets->push_back(jet);
      //cout << newJetPtJER << " " <<newJetPtJES_up << " " << newJetPtJES_down << endl;
    }
    
    nCentralJERjets = CentralJER_jets->size();
    nUpJERjets = UpJER_jets->size();
    nDownJERjets = DownJER_jets->size();
    nCentralJERcentraljets = CentralJER_centraljets->size();
    nUpJERcentraljets = UpJER_centraljets->size();
    nDownJERcentraljets = DownJER_centraljets->size();    
    
    if (nCentralJERjets>3) nCentralJERjets=3;
    if (nUpJERjets>3) nUpJERjets=3;
    if (nDownJERjets>3) nDownJERjets=3;
    if (nCentralJERcentraljets>3) nCentralJERcentraljets=3;
    if (nUpJERcentraljets>3) nUpJERcentraljets=3;
    if (nDownJERcentraljets>3) nDownJERcentraljets=3;    


    nUpJESjets = UpJES_jets->size();
    nDownJESjets = DownJES_jets->size();
    nUpJEScentraljets = UpJES_centraljets->size();
    nDownJEScentraljets = DownJES_centraljets->size();
    
    if (nUpJESjets>3) nUpJESjets=3;
    if (nDownJESjets>3) nDownJESjets=3;
    if (nUpJEScentraljets>3) nUpJEScentraljets=3;
    if (nDownJEScentraljets>3) nDownJEScentraljets=3;
    
    stable_sort(CentralJER_jets->begin(), CentralJER_jets->end(), PtComparator());
    stable_sort(UpJER_jets->begin(), UpJER_jets->end(), PtComparator());
    stable_sort(DownJER_jets->begin(), DownJER_jets->end(), PtComparator());
    stable_sort(UpJES_jets->begin(), UpJES_jets->end(), PtComparator());
    stable_sort(DownJES_jets->begin(), DownJES_jets->end(), PtComparator());
    stable_sort(CentralJER_centraljets->begin(), CentralJER_centraljets->end(), PtComparator());
    stable_sort(UpJER_centraljets->begin(), UpJER_centraljets->end(), PtComparator());
    stable_sort(DownJER_centraljets->begin(), DownJER_centraljets->end(), PtComparator());
    stable_sort(UpJES_centraljets->begin(), UpJES_centraljets->end(), PtComparator());
    stable_sort(DownJES_centraljets->begin(), DownJES_centraljets->end(), PtComparator());
  
    stable_sort(CentralJER_jetPt->begin(), CentralJER_jetPt->end(),std::greater<float>());
    stable_sort(UpJER_jetPt->begin(), UpJER_jetPt->end(),std::greater<float>());
    stable_sort(DownJER_jetPt->begin(), DownJER_jetPt->end(),std::greater<float>());
    stable_sort(UpJES_jetPt->begin(), UpJES_jetPt->end(),std::greater<float>());
    stable_sort(DownJES_jetPt->begin(), DownJES_jetPt->end(),std::greater<float>());


    //to check if jets associated with leptons are correctly deleted
    float dRl10 =0;   
    float dRl20 =0;   
    float dRl11 =0;   
    float dRl21 =0;
    
    
    foreach(const phys::Jet &jet, *jets){
      theHistograms.fill("uncOnFourVectorScale", 50, 0, 0.1, jet.uncOnFourVectorScale(), theWeight*w_kf);
      if(jets->size()>=2){
	dRl10 =  deltaR(jet.eta(),jet.phi(),ZZ->first().daughter(0).eta(),ZZ->first().daughter(0).phi());
	dRl11 =  deltaR(jet.eta(),jet.phi(),ZZ->first().daughter(1).eta(),ZZ->first().daughter(1).phi());
	dRl20 =  deltaR(jet.eta(),jet.phi(),ZZ->second().daughter(0).eta(),ZZ->second().daughter(0).phi());
	dRl21 =  deltaR(jet.eta(),jet.phi(),ZZ->second().daughter(1).eta(),ZZ->second().daughter(1).phi());
	float dR[] ={dRl10,dRl11,dRl20,dRl21};
	min_dR = *std::min_element(dR,dR+4);
      }
    
    }
    //     // ZZplots();   // ZZ --> 4l 
    ZZplots(52,e); // ZZ --> 4m
    ZZplots(48,e); // ZZ --> 2e2m
    ZZplots(44,e); // ZZ --> 4e
    
}
  
void ZZRecoAnalyzer::begin() {
  
  UpJESData_jets  = new std::vector<phys::Jet>();
  DownJESData_jets  = new std::vector<phys::Jet>();
  UpJESData_centraljets  = new std::vector<phys::Jet>();
  DownJESData_centraljets  = new std::vector<phys::Jet>();

  UpJESData_jetPt  = new std::vector<double>();
  DownJESData_jetPt  = new std::vector<double>();
  
  CentralJER_jets  = new std::vector<phys::Jet>();
  UpJER_jets  = new std::vector<phys::Jet>();
  DownJER_jets  = new std::vector<phys::Jet>();
  UpJES_jets  = new std::vector<phys::Jet>();
  DownJES_jets  = new std::vector<phys::Jet>();
  CentralJER_centraljets  = new std::vector<phys::Jet>();
  UpJER_centraljets  = new std::vector<phys::Jet>();
  DownJER_centraljets  = new std::vector<phys::Jet>();
  UpJES_centraljets  = new std::vector<phys::Jet>();
  DownJES_centraljets  = new std::vector<phys::Jet>();

  CentralJER_jetPt  = new std::vector<double>();
  UpJER_jetPt  = new std::vector<double>();
  DownJER_jetPt  = new std::vector<double>();
  UpJES_jetPt  = new std::vector<double>();
  DownJES_jetPt  = new std::vector<double>();

  jet_count = 0;
  inFiducialRegion = 0;
  inGenAndReco =0;
  nEvent = 0;

  nUpJESDatajets = 0;
  nDownJESDatajets = 0;
  nUpJESDatacentraljets = 0;
  nDownJESDatacentraljets = 0;

  nCentralJERjets=0;
  nUpJERjets=0;
  nDownJERjets=0; 
  nUpJESjets=0; 
  nDownJESjets=0;
  nCentralJERcentraljets=0;
  nUpJERcentraljets=0;
  nDownJERcentraljets=0; 
  nUpJEScentraljets=0; 
  nDownJEScentraljets=0;
 
  nentries = tree()->GetEntries("ZZCand.passFullSel_"); 
  Xbins += 100,200,250,300,350,400,500,600,800;
  Xbins_deta += 0,2.4,4.7;
  Xbins_mjj += 0.,200,800;
  Xbins_ptjet1 += 30,50,100,200,300; //30,50,100,200,300,500;
  Xbins_ptjet2 += 30,100,200; //30,100,200,500;  
  Xbins_etajet1 += 0,1.5,3,4.7;
  Xbins_etajet2 +=0,1.5,3,4.7;
  Xbins_dphi +=  0,2,3,4; 
  Xbins_drzz += 0,1,2,3,4,5,6;
  Xbins_ptzz += 0,25,50,75,100,150,200,300;
  Xbins_dphizz += 0,1.5,2.,2.25,2.5,2.75,3,3.25;

  m4L = 0; 
  drzz= 0;
  m4L_gen = 0; 
  drzz_gen= 0;
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
  ptzz = 0;
  ptzz_gen = 0;
  dphizz = 0; 
  dphizz_gen = 0;
  min_dR=0;   
}

void ZZRecoAnalyzer::end( TFile &) {
  cout<<theMCInfo.analyzedEvents()<< " " << e <<endl;
  cout << "Events Reconstructed in the fiducial region  and generated with the signal definition " << inGenAndReco << endl;   
  cout << "Events in the Fiducial Region " << inFiducialRegion << endl;
  
  if(region_ == phys::CR3P1F || region_ == phys::CR2P2F) {
    vector<std::string>  FinalState = {"4m","4e","2e2m"};
    vector<std::string>  Variable = {"Mass","Jets","CentralJets","PtJet1","EtaJet1","PtJet2","EtaJet2","Deta","Mjj","CentralDeta","CentralMjj","dRZZ", "PtZZ","DphiZZ"};
   
    for (std::vector<std::string>::iterator var = Variable.begin() ; var != Variable.end(); ++var){
      for (std::vector<std::string>::iterator it = FinalState.begin() ; it != FinalState.end(); ++it){
	
  	TH1 *hvar =  theHistograms.get(("ZZTo"+*it+"_"+*var+"_FRVar").c_str());
	
  	TH1 *h = theHistograms.get(("ZZTo"+*it+"_"+*var+"_01").c_str());
	
  	if(!h) continue;
  	for(int i = 1; i<=h->GetNbinsX();i++){
	  
  	  Float_t Err = h->GetBinError(i);
  	  h->SetBinError(i,sqrt(Err*Err+hvar->GetBinContent(i)));
  	}
      }
    }  
  }
}
