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

  if(ZZ->id() != id && id != -1) return; // -1 here means generic 4l final state
  
  std::string decay  = "4l";
  std::string sample;
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
  
  // //////////////////////////////////////////////DATA ONLY///////////////////////////////////////////////////////
  
 //  //**********************************************JETS*********************************************************
  
 //  // //1D Reco nJet Distributions (no JER smearing, for data)
 //  // theHistograms.fill(std::string("ZZTo")+decay+"_Jets_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight); 
 //  // theHistograms.fill(std::string("ZZTo")+decay+"_Jets_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight); 
  
 //  //1D Reco Deta and Mjj Distributions (no JER smearing, for data)
 //  if(njets >=2){
 //    float deta = fabs(jets->at(0).eta() - jets->at(1).eta());
 //    if (deta>=4.7) deta = 4.6;
 //    float  mjj =  (jets->at(0).p4() + jets->at(1).p4()).M();
 //    if (mjj>=800) mjj = 799;
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta,theWeight); 
 //  }

 //  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESDataUpSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJESDatajets , theWeight);
 //  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESDataDownSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJESDatajets , theWeight);
 //  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESDataUpSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJESDatajets , theWeight);
 //  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESDataDownSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJESDatajets , theWeight);
  
 //  if(nUpJESDatajets >=2){
 //    float upDetaJESData = fabs(UpJESData_jets->at(0).eta() - UpJESData_jets->at(1).eta());
 //    float upMjjJESData =  (UpJESData_jets->at(0).p4() + UpJESData_jets->at(1).p4()).M();
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JESDataUpSmear_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,upMjjJESData,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JESDataUpSmear_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,upMjjJESData,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JESDataUpSmear_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,upDetaJESData,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JESDataUpSmear_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,upDetaJESData,theWeight); 
 //  }
 //  if(nDownJESDatajets >=2){
 //    float downDetaJESData = fabs(DownJESData_jets->at(0).eta() - DownJESData_jets->at(1).eta());
 //    float downMjjJESData =  (DownJESData_jets->at(0).p4() + DownJESData_jets->at(1).p4()).M();
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JESDataDownSmear_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,downMjjJESData,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JESDataDownSmear_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,downMjjJESData,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JESDataDownSmear_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,downDetaJESData,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JESDataDownSmear_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,downDetaJESData,theWeight); 
 // }
  
  

 //  //**********************************************CENTRAL JETS*********************************************************
  
 //  // //1D Reco nJet Distributions (no JER smearing, for data)
 //  // theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_"+sample, std::string("Number of centraljets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theWeight); 
 //  // theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_01", std::string("Number of centraljets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theWeight); 
  
 //  //1D Reco Deta and Mjj Distributions (no JER smearing, for data)
 //  if(ncentraljets >=2){
 //    float centraldeta = fabs(centralJets->at(0).eta() - centralJets->at(1).eta());
 //    if (centraldeta >=4.7) centraldeta = 4.6;
 //    float  centralmjj =  (centralJets->at(0).p4() + centralJets->at(1).p4()).M();
 //    if (centralmjj>=800) centralmjj = 799;
 //    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,centralmjj,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,centralmjj,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,centraldeta,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,centraldeta,theWeight); 
 //  }
  
 //  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JESDataUpSmear_"+sample, "Number of reco centralJets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJESDatacentraljets , theWeight);
 //  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JESDataDownSmear_"+sample, "Number of reco centralJets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJESDatacentraljets , theWeight);
 //  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JESDataUpSmear_01", "Number of reco centralJets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJESDatacentraljets , theWeight);
 //  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JESDataDownSmear_01", "Number of reco centralJets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJESDatacentraljets , theWeight);
  
 //  if(nUpJESDatacentraljets >=2){
 //    float upCentraldetaJESData = fabs(UpJESData_centraljets->at(0).eta() - UpJESData_centraljets->at(1).eta());
 //    float upCentralmjjJESData =  (UpJESData_centraljets->at(0).p4() + UpJESData_centraljets->at(1).p4()).M();
 //    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JESDataUpSmear_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,upCentralmjjJESData,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JESDataUpSmear_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,upCentralmjjJESData,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JESDataUpSmear_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,upCentraldetaJESData,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JESDataUpSmear_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,upCentraldetaJESData,theWeight); 
 //  }
 // if(nDownJESDatacentraljets >=2){
 //    float downCentraldetaJESData = fabs(DownJESData_centraljets->at(0).eta() - DownJESData_centraljets->at(1).eta());
 //    float downCentralmjjJESData =  (DownJESData_centraljets->at(0).p4() + DownJESData_centraljets->at(1).p4()).M();
 //    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JESDataDownSmear_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,downCentralmjjJESData,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JESDataDownSmear_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,downCentralmjjJESData,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JESDataDownSmear_"+sample, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,downCentraldetaJESData,theWeight); 
 //    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JESDataDownSmear_01", std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,downCentraldetaJESData,theWeight); 
 //  }
 
 //***********************************************************************************************************************************************************

  ////////////////////////////////////////////////////DATA AND MONTECARLO/////////////////////////////////////////////////

  //1D Reco Mass Distributions
  m4L = ZZ->mass();  
  if(m4L > 800) m4L = 799;
  theHistograms.fill(std::string("ZZTo")+decay+"_Mass_01", std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L,theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_Mass_"+sample, std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L,theWeight);
 
  
 // //Pt and eta distribution for the leading, sub-leading and sub-sub-leading jets (if they exist)
 //  if(npjets>0){  
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Pt0_01", std::string("p_{t}^{jet0} of ZZ_{1}#rightarrow ")+decay,50,0,350,pjets->at(0).pt(), theWeight);  
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Eta0_01", std::string("eta^{jet0} of ZZ_{1}#rightarrow ")+decay,60,-6,6,pjets->at(0).eta(), theWeight);  
 //  }
 //  if(npjets>1){  
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Pt1_01", std::string("p_{t}^{jet1} of ZZ_{1}#rightarrow ")+decay,50,0,350,pjets->at(1).pt(), theWeight);  
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Eta1_01", std::string("eta^{jet1} of ZZ_{1}#rightarrow ")+decay,60,-6,6,pjets->at(1).eta(), theWeight);  
 //  }
 //  if(npjets>2){  
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Pt2_01", std::string("p_{t}^{jet2} of ZZ_{1}#rightarrow ")+decay,50,0,350,pjets->at(2).pt(), theWeight);  
 //    theHistograms.fill(std::string("ZZTo")+decay+"_Eta2_01", std::string("eta^{jet2} of ZZ_{1}#rightarrow ")+decay,60,-6,6,pjets->at(2).eta(), theWeight);  
 //  }
  

  //1D Reco nJet Distributions (no JER smearing, for data)
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight); 

 //1D Reco ncentralJet Distributions (no JER smearing, for data)
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,ncentraljets,theWeight); 
  
  //1D Reco nPJet Distributions (no JER smearing, for data)
  theHistograms.fill(std::string("ZZTo")+decay+"_PJets_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,npjets,theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_PJets_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,npjets,theWeight); 

  /////////////////////////////////////////////////////MONTECARLO ONLY//////////////////////////////////////////////////////


  //******************************************************************RECO JETS****************************************************************************************

  //1D Reco nJet Distributions - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERjets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERjets , theWeight);
  
  //1D Reco nJet Distributions - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERUpSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJERjets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERDownSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJERjets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERUpSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJERjets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERDownSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJERjets , theWeight);

 //1D Reco nJet Distributions - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESUpSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJESjets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESDownSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJESjets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESUpSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJESjets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESDownSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJESjets , theWeight);
  
  // //Pt and eta distribution for the leading, sub-leading and sub-sub-leading jets (if they exist)
  // if(nCentralJERjets>0){  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Pt0_JERCentralSmear_01", std::string("p_{t}^{jet0} of ZZ_{1}#rightarrow ")+decay,50,0,350,pjets->at(0).pt(), theWeight);  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Eta0_JERCentralSmear_01", std::string("eta^{jet0} of ZZ_{1}#rightarrow ")+decay,30,0,6,pjets->at(0).eta(), theWeight);  
  // }
  // if(nCentralJERjets>1){  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Pt1_JERCentralSmear_01", std::string("p_{t}^{jet1} of ZZ_{1}#rightarrow ")+decay,50,0,350,pjets->at(1).pt(), theWeight);  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Eta1_JERCentralSmear_01", std::string("eta^{jet1} of ZZ_{1}#rightarrow ")+decay,30,0,6,pjets->at(1).eta(), theWeight);  
  // }
  // if(nCentralJERjets>2){  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Pt2_JERCentralSmear_01", std::string("p_{t}^{jet2} of ZZ_{1}#rightarrow ")+decay,50,0,350,pjets->at(2).pt(), theWeight);  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Eta2_JERCentralSmear_01", std::string("eta^{jet2} of ZZ_{1}#rightarrow ")+decay,30,0,6,pjets->at(2).eta(), theWeight);  
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
  
  //1D Reco DeltaEta and mJJ Distributions - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)
   if(nCentralJERjets>=2){
  
    centralDeta = fabs(CentralJER_jets->at(0).eta() - CentralJER_jets->at(1).eta());
    centralMjj =  (CentralJER_jets->at(0).p4() + CentralJER_jets->at(1).p4()).M();
    
    if (centralDeta>=4.7) centralDeta = 4.6;
    if (centralMjj>=800) centralMjj = 799;

    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERCentralSmear_"+sample, "reco m_{jj}", Xbins_mjj, centralMjj, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERCentralSmear_01", "reco m_{jj}", Xbins_mjj, centralMjj, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERCentralSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, centralDeta, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERCentralSmear_01", "reco  #Delta#eta_{jj}", Xbins_deta, centralDeta, theWeight);
  }
  
  //1D Reco DeltaEta and mJJ Distributions - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
  if(nUpJERjets>=2){

    upDetaJER = fabs(UpJER_jets->at(0).eta() - UpJER_jets->at(1).eta());
    upMjjJER =  (UpJER_jets->at(0).p4() + UpJER_jets->at(1).p4()).M();

    if (upDetaJER>=4.7) upDetaJER = 4.6;
    if (upMjjJER>=800) upMjjJER = 799;

    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERUpSmear_"+sample, "reco m_{jj}", Xbins_mjj, upMjjJER, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERUpSmear_01", "reco m_{jj}", Xbins_mjj, upMjjJER, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERUpSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, upDetaJER, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERUpSmear_01", "reco #Delta#eta_{jj}", Xbins_deta, upDetaJER, theWeight);
  }
 
  if(nDownJERjets>=2){

    downDetaJER = fabs(DownJER_jets->at(0).eta() - DownJER_jets->at(1).eta());
    downMjjJER =  (DownJER_jets->at(0).p4() + DownJER_jets->at(1).p4()).M();
 
    if (downDetaJER>=4.7) downDetaJER = 4.6;
    if (downMjjJER>=800) downMjjJER = 799;
    
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERDownSmear_"+sample, "reco m_{jj}", Xbins_mjj, downMjjJER, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERDownSmear_01", "reco m_{jj}", Xbins_mjj, downMjjJER, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERDownSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, downDetaJER, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERDownSmear_01", "reco #Delta#eta_{jj}", Xbins_deta, downDetaJER, theWeight);
  }
 
  //1D Reco DeltaEta and mJJ Distributions - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
  if(nUpJESjets>=2){

    upDetaJES = fabs(UpJES_jets->at(0).eta() - UpJES_jets->at(1).eta());
    upMjjJES =  (UpJES_jets->at(0).p4() + UpJES_jets->at(1).p4()).M();

    if (upDetaJES>=4.7) upDetaJES = 4.6;
    if (upMjjJES>=800) upMjjJES = 799;

    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JESUpSmear_"+sample, "reco m_{jj}", Xbins_mjj, upMjjJES, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JESUpSmear_01", "reco m_{jj}", Xbins_mjj, upMjjJES, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JESUpSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, upDetaJES, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JESUpSmear_01", "reco #Delta#eta_{jj}", Xbins_deta, upDetaJES, theWeight);
  }
  
  if(nDownJESjets>=2){

    downDetaJES = fabs(DownJES_jets->at(0).eta() - DownJES_jets->at(1).eta());
    downMjjJES =  (DownJES_jets->at(0).p4() + DownJES_jets->at(1).p4()).M();

    if (downDetaJES>=4.7) downDetaJES = 4.6;
    if (downMjjJES>=800) downMjjJES = 799;

    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JESDownSmear_"+sample, "reco m_{jj}", Xbins_mjj, downMjjJES, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JESDownSmear_01", "reco m_{jj}", Xbins_mjj, downMjjJES, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JESDownSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, downDetaJES, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JESDownSmear_01", "reco #Delta#eta_{jj}", Xbins_deta, downDetaJES, theWeight);
  }


 //******************************************************************RECO CENTRAL JETS****************************************************************************************

  //1D Reco nJet Distributions - JER smearing (CentralJets_JERCentralSmear to be used in the standard analysis)
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERCentralSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERcentraljets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERCentralSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERcentraljets , theWeight);
  
  //1D Reco nJet Distributions - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERUpSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJERcentraljets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERDownSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJERcentraljets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERUpSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJERcentraljets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JERDownSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJERcentraljets , theWeight);

 //1D Reco nJet Distributions - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JESUpSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJEScentraljets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JESDownSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJEScentraljets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JESUpSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJEScentraljets , theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_CentralJets_JESDownSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJEScentraljets , theWeight);
  
  // //Pt and eta distribution for the leading, sub-leading and sub-sub-leading centraljets (if they exist)
  // if(nCentralJERcentraljets>0){  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Pt0_JERCentralSmear_01", std::string("p_{t}^{jet0} of ZZ_{1}#rightarrow ")+decay,50,0,350,FIXMEpjets->at(0).pt(), theWeight);  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Eta0_JERCentralSmear_01", std::string("eta^{jet0} of ZZ_{1}#rightarrow ")+decay,30,0,6,pjets->at(0).eta(), theWeight);  
  // }
  // if(nCentralJERcentraljets>1){  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Pt1_JERCentralSmear_01", std::string("p_{t}^{jet1} of ZZ_{1}#rightarrow ")+decay,50,0,350,pjets->at(1).pt(), theWeight);  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Eta1_JERCentralSmear_01", std::string("eta^{jet1} of ZZ_{1}#rightarrow ")+decay,30,0,6,pjets->at(1).eta(), theWeight);  
  // }
  // if(nCentralJERcentraljets>2){  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Pt2_JERCentralSmear_01", std::string("p_{t}^{jet2} of ZZ_{1}#rightarrow ")+decay,50,0,350,pjets->at(2).pt(), theWeight);  
  //   theHistograms.fill(std::string("ZZTo")+decay+"_Eta2_JERCentralSmear_01", std::string("eta^{jet2} of ZZ_{1}#rightarrow ")+decay,30,0,6,pjets->at(2).eta(), theWeight);  
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

    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERCentralSmear_"+sample, "reco m_{jj}", Xbins_mjj, centralMjj_cj, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERCentralSmear_01", "reco m_{jj}", Xbins_mjj, centralMjj_cj, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERCentralSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, centralDeta_cj, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERCentralSmear_01", "reco  #Delta#eta_{jj}", Xbins_deta, centralDeta_cj, theWeight);
  }
  
  //1D Reco DeltaEta and mJJ Distributions - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
  if(nUpJERcentraljets>=2){
    
    upDetaJER_cj = fabs(UpJER_centraljets->at(0).eta() - UpJER_centraljets->at(1).eta());
    upMjjJER_cj =  (UpJER_centraljets->at(0).p4() + UpJER_centraljets->at(1).p4()).M();
    
    if (upDetaJER_cj>=4.7) upDetaJER_cj = 4.6;
    if (upMjjJER_cj>=800) upMjjJER_cj = 799;
    
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERUpSmear_"+sample, "reco m_{jj}", Xbins_mjj, upMjjJER_cj, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERUpSmear_01", "reco m_{jj}", Xbins_mjj, upMjjJER_cj, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERUpSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, upDetaJER_cj, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERUpSmear_01", "reco #Delta#eta_{jj}", Xbins_deta, upDetaJER_cj, theWeight);
  }

  if(nDownJERcentraljets>=2){

    downDetaJER_cj = fabs(DownJER_centraljets->at(0).eta() - DownJER_centraljets->at(1).eta());
    downMjjJER_cj =  (DownJER_centraljets->at(0).p4() + DownJER_centraljets->at(1).p4()).M();
 
    if (downDetaJER_cj>=4.7) downDetaJER_cj = 4.6;
    if (downMjjJER_cj>=800) downMjjJER_cj = 799;
    
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERDownSmear_"+sample, "reco m_{jj}", Xbins_mjj, downMjjJER_cj, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JERDownSmear_01", "reco m_{jj}", Xbins_mjj, downMjjJER_cj, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERDownSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, downDetaJER_cj, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JERDownSmear_01", "reco #Delta#eta_{jj}", Xbins_deta, downDetaJER_cj, theWeight);
  }
 
  //1D Reco DeltaEta and mJJ Distributions - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
  if(nUpJEScentraljets>=2){

    upDetaJES_cj = fabs(UpJES_centraljets->at(0).eta() - UpJES_centraljets->at(1).eta());
    upMjjJES_cj =  (UpJES_centraljets->at(0).p4() + UpJES_centraljets->at(1).p4()).M();

    if (upDetaJES_cj>=4.7) upDetaJES_cj = 4.6;
    if (upMjjJES_cj>=800) upMjjJES_cj = 799;

    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JESUpSmear_"+sample, "reco m_{jj}", Xbins_mjj, upMjjJES_cj, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JESUpSmear_01", "reco m_{jj}", Xbins_mjj, upMjjJES_cj, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JESUpSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, upDetaJES_cj, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JESUpSmear_01", "reco #Delta#eta_{jj}", Xbins_deta, upDetaJES_cj, theWeight);
  }
  
  if(nDownJEScentraljets>=2){

    downDetaJES_cj = fabs(DownJES_centraljets->at(0).eta() - DownJES_centraljets->at(1).eta());
    downMjjJES_cj =  (DownJES_centraljets->at(0).p4() + DownJES_centraljets->at(1).p4()).M();

    if (downDetaJES_cj>=4.7) downDetaJES_cj = 4.6;
    if (downMjjJES_cj>=800) downMjjJES_cj = 799;

    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JESDownSmear_"+sample, "reco m_{jj}", Xbins_mjj, downMjjJES_cj, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralMjj_JESDownSmear_01", "reco m_{jj}", Xbins_mjj, downMjjJES_cj, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JESDownSmear_"+sample, "reco #Delta#eta_{jj}", Xbins_deta, downDetaJES_cj, theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_CentralDeta_JESDownSmear_01", "reco #Delta#eta_{jj}", Xbins_deta, downDetaJES_cj, theWeight);
  }

  
  //if MC gen (for response matrices only)
  if (genCategory !=-1){
    if(topology.test(0)){
      

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
   
//   Float_t scaleFacErr = 0;

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
 
      theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqPlus_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*(1+scaleFacErrSq));  

      theHistograms.fill(std::string("ZZTo")+decay+"_JetsGenRecoSFErrSqMinus_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight*(1-scaleFacErrSq));   


     // stable_sort(genJets->begin(), genJets->end(), PtComparator());
      
      ngenjets =  genJets->size(); 
      if (ngenjets>3) ngenjets=3;
      
      ngencentraljets =  centralGenJets->size(); 
      if (ngencentraljets>3) ngencentraljets=3;
      
      m4L_gen = sqrt((genVBParticles->at(0).p4()+genVBParticles->at(1).p4())*(genVBParticles->at(0).p4()+genVBParticles->at(1).p4()));
      if (m4L_gen>=800) m4L_gen = 799;
      
      //Response Matrix Mass (Reco&Gen) 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_"+sample, std::string("Response matrix invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, Xbins, m4L ,m4L_gen , theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_01", std::string("Response matrix invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, Xbins, m4L ,m4L_gen , theWeight);
      
      //****************************************************************************RECO&GEN JETS****************************************************************************************

      //Response Matrix nJets (Reco&Gen) - No JER smearing
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4, njets,ngenjets, theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4, njets,ngenjets, theWeight); 
      
      //Response Matrix nJets (Reco&Gen) - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERCentralSmear_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERCentralSmear_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight);
      
      // Response Matrix nJets (Reco&Gen) - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERUpSmear_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nUpJERjets,ngenjets, theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERUpSmear_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nUpJERjets,ngenjets, theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERDownSmear_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nDownJERjets,ngenjets, theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERDownSmear_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nDownJERjets,ngenjets, theWeight);
      
      // Response Matrix nJets (Reco&Gen) - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JESUpSmear_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nUpJESjets,ngenjets, theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JESUpSmear_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nUpJESjets,ngenjets, theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JESDownSmear_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nDownJESjets,ngenjets, theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JESDownSmear_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nDownJESjets,ngenjets, theWeight);
    
      if(ngenjets>=2){  
 	deta_gen = fabs(genJets->at(0).eta() - genJets->at(1).eta());
 	mjj_gen =  (genJets->at(0).p4() + genJets->at(1).p4()).M();
	 
 	if (deta_gen>=4.7) deta_gen = 4.6;
 	if (mjj_gen>=800) mjj_gen = 799;
 	
 	//Response Matrix  DeltaEta and mJJ (Reco&Gen) - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)
      	if(nCentralJERjets>=2){
      	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERCentralSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,Xbins_mjj, centralMjj,mjj_gen,theWeight); 
      	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERCentralSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_mjj,Xbins_mjj,centralMjj,mjj_gen,theWeight);  
      	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERCentralSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,Xbins_deta, centralDeta,deta_gen,theWeight); 
      	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERCentralSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_deta,Xbins_deta,centralDeta,deta_gen,theWeight); 
	  
 	// //Matrix and reco distribution weighted for the ratio between the unfolded data and the generator level information in order to 
 	// //compute the relative systematic uncertainty
 	//   string UnfOverMC_Mjj = "ZZTo"+decay+"_Mjj_Ratio";  
 	//   h_UnfOverMC_Mjj = (TH1*) UnfOverMC->Get(UnfOverMC_Mjj.c_str()); 
 	//   int bin_Mjj = h_UnfOverMC_Mjj->FindBin(mjj_gen);  
 	//   float w_Mjj =  h_UnfOverMC_Mjj->GetBinContent(bin_Mjj);  
	  
 	//   //1D reco distribution built not for all reco events, but only events gen&reco
 	//   theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERCentralSmear_W_"+sample, "m_{jj}", Xbins_mjj,centralMjj,theWeight*w_Mjj);
 	//   theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERCentralSmear_W_01", "m_{jj}", Xbins_mjj, centralMjj,theWeight*w_Mjj);
 	//   theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERCentralSmear_W_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_mjj,Xbins_mjj, centralMjj,mjj_gen,theWeight*w_Mjj);
 	//   theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERCentralSmear_W_01", std::string("Response matrix  m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,Xbins_mjj, centralMjj,mjj_gen,theWeight*w_Mjj);
	  
 	//   string UnfOverMC_Deta = "ZZTo"+decay+"_Deta_Ratio";  
 	//   h_UnfOverMC_Deta = (TH1*) UnfOverMC->Get(UnfOverMC_Deta.c_str()); 
 	//   int bin_Deta = h_UnfOverMC_Deta->FindBin(deta_gen);  
 	//   float w_Deta =  h_UnfOverMC_Deta->GetBinContent(bin_Deta);  
	  
 	//   //1D reco distribution built not for all reco events, but only events gen&reco
 	//   theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERCentralSmear_W_"+sample, "#Delta#eta_{jj}", Xbins_deta,centralDeta,theWeight*w_Deta);
 	//   theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERCentralSmear_W_01", "#Delta#eta_{jj}", Xbins_deta, centralDeta,theWeight*w_Deta);
	  
 	//   //Matrix and reco distribution weighted for the ratio between the unfolded data and the generator level information in order to 
 	//   //compute the relative systematic uncertainty
 	//   theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERCentralSmear_W_"+sample, std::string("Response matrix #Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_deta,Xbins_deta, centralDeta,deta_gen,theWeight*w_Deta);
 	//   theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERCentralSmear_W_01", std::string("Response matrix  #Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,Xbins_deta, centralDeta,deta_gen,theWeight*w_Deta);
 	}
	
 	// Response Matrix  DeltaEta and mJJ (Reco&Gen) - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
 	if(nUpJERjets>=2){
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERUpSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,Xbins_mjj, upMjjJER,mjj_gen,theWeight); 
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERUpSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_mjj,Xbins_mjj,upMjjJER,mjj_gen,theWeight);  
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERUpSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,Xbins_deta, upDetaJER,deta_gen,theWeight); 
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERUpSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_deta,Xbins_deta,upDetaJER,deta_gen,theWeight); 
 	}
	
 	if(nDownJERjets>=2){
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERDownSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,Xbins_mjj, downMjjJER,mjj_gen,theWeight); 
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JERDownSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_mjj,Xbins_mjj,downMjjJER,mjj_gen,theWeight);  
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERDownSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,Xbins_deta, downDetaJER,deta_gen,theWeight); 
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JERDownSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_deta,Xbins_deta,downDetaJER,deta_gen,theWeight); 
 	}
	
 	// Response Matrix  DeltaEta and mJJ (Reco&Gen) - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
 	if(nUpJESjets>=2){
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JESUpSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,Xbins_mjj, upMjjJES,mjj_gen,theWeight); 
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JESUpSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_mjj,Xbins_mjj,upMjjJES,mjj_gen,theWeight);  
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JESUpSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,Xbins_deta, upDetaJES,deta_gen,theWeight); 
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JESUpSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_deta,Xbins_deta,upDetaJES,deta_gen,theWeight); 
 	}
	
 	if(nDownJESjets>=2){
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JESDownSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,Xbins_mjj, downMjjJES,mjj_gen,theWeight); 
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mjj_JESDownSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_mjj,Xbins_mjj,downMjjJES,mjj_gen,theWeight);  
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JESDownSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,Xbins_deta, downDetaJES,deta_gen,theWeight); 
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Deta_JESDownSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_deta,Xbins_deta,downDetaJES,deta_gen,theWeight); 
 	}
      }
      

      //****************************************************************************RECO&GEN CENTRAL JETS****************************************************************************************

      //Response Matrix nCentralJets (Reco&Gen) - No JER smearing
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4, ncentraljets,ngencentraljets, theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4, ncentraljets,ngencentraljets, theWeight); 
      
      //Response Matrix nCentralJets (Reco&Gen) - JER smearing (CentralJets_JERCentralSmear to be used in the standard analysis)
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERCentralSmear_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERCentralSmear_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nCentralJERcentraljets,ngencentraljets, theWeight);
      
      // Response Matrix nCentralJets (Reco&Gen) - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERUpSmear_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nUpJERcentraljets,ngencentraljets, theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERUpSmear_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nUpJERcentraljets,ngencentraljets, theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERDownSmear_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nDownJERcentraljets,ngencentraljets, theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JERDownSmear_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nDownJERcentraljets,ngencentraljets, theWeight);
      
      // Response Matrix nCentralJets (Reco&Gen) - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JESUpSmear_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nUpJEScentraljets,ngencentraljets, theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JESUpSmear_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nUpJEScentraljets,ngencentraljets, theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JESDownSmear_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nDownJEScentraljets,ngencentraljets, theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralJets_JESDownSmear_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nDownJEScentraljets,ngencentraljets, theWeight);
    
      if(ngencentraljets>=2){  
 	deta_gen_cj = fabs(centralGenJets->at(0).eta() - centralGenJets->at(1).eta());
 	mjj_gen_cj =  (centralGenJets->at(0).p4() + centralGenJets->at(1).p4()).M();
	 
 	if (deta_gen_cj>=4.7) deta_gen_cj = 4.6;
 	if (mjj_gen_cj>=800) mjj_gen_cj = 799;
 	
 	//Response Matrix  DeltaEta and mJJ (Reco&Gen) - JER smearing (Jets_JERCentralSmear to be used in the standard analysis)
      	if(nCentralJERcentraljets>=2){
      	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERCentralSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,Xbins_mjj, centralMjj_cj,mjj_gen_cj,theWeight); 
      	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERCentralSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_mjj,Xbins_mjj,centralMjj_cj,mjj_gen_cj,theWeight);  
      	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERCentralSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,Xbins_deta, centralDeta_cj,deta_gen_cj,theWeight); 
      	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERCentralSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_deta,Xbins_deta,centralDeta_cj,deta_gen_cj,theWeight); 
	  
 	}
	
 	// Response Matrix  DeltaEta and mJJ (Reco&Gen) - JER smearing Up and Down (To be used in the evaluation of JER systematics uncertainties)
 	if(nUpJERcentraljets>=2){
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERUpSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,Xbins_mjj, upMjjJER_cj,mjj_gen_cj,theWeight); 
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERUpSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_mjj,Xbins_mjj,upMjjJER_cj,mjj_gen_cj,theWeight);  
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERUpSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,Xbins_deta, upDetaJER_cj,deta_gen_cj,theWeight); 
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERUpSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_deta,Xbins_deta,upDetaJER_cj,deta_gen_cj,theWeight); 
 	}
	
 	if(nDownJERcentraljets>=2){
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERDownSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,Xbins_mjj, downMjjJER_cj,mjj_gen_cj,theWeight); 
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JERDownSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_mjj,Xbins_mjj,downMjjJER_cj,mjj_gen_cj,theWeight);  
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERDownSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,Xbins_deta, downDetaJER_cj,deta_gen_cj,theWeight); 
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JERDownSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_deta,Xbins_deta,downDetaJER_cj,deta_gen_cj,theWeight); 
 	}
	
 	// Response Matrix  DeltaEta and mJJ (Reco&Gen) - JES smearing Up and Down (To be used in the evaluation of JES systematics uncertainties)
 	if(nUpJEScentraljets>=2){
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JESUpSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,Xbins_mjj, upMjjJES_cj,mjj_gen_cj,theWeight); 
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JESUpSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_mjj,Xbins_mjj,upMjjJES_cj,mjj_gen_cj,theWeight);  
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JESUpSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,Xbins_deta, upDetaJES_cj,deta_gen_cj,theWeight); 
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JESUpSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_deta,Xbins_deta,upDetaJES_cj,deta_gen_cj,theWeight); 
 	}
	
 	if(nDownJEScentraljets>=2){
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JESDownSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,Xbins_mjj, downMjjJES_cj,mjj_gen_cj,theWeight); 
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralMjj_JESDownSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_mjj,Xbins_mjj,downMjjJES_cj,mjj_gen_cj,theWeight);  
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JESDownSmear_"+sample, std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,Xbins_deta, downDetaJES_cj,deta_gen_cj,theWeight); 
 	  theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_CentralDeta_JESDownSmear_01", std::string("Response matrix m_{jj} of ZZ_{1}#rightarrow ")+decay, Xbins_deta,Xbins_deta,downDetaJES_cj,deta_gen_cj,theWeight); 
 	}
      }
      
      //*************************************************************************************************************************************************************************************************

	
 	// //Matrix and reco distribution weighted for the ratio between the unfolded data and the generator level information in order to 
 	// //compute the relative systematic uncertainty
 	// //distributions to evaluate data/MC systematic uncertainty
 	// string UnfOverMC_Jets = "ZZTo"+decay+"_Jets_Ratio";  
 	// h_UnfOverMC_Jets = (TH1*) UnfOverMC->Get(UnfOverMC_Jets.c_str()); 
 	// int bin_Jets = h_UnfOverMC_Jets->FindBin(ngencentraljets);  
 	// float w_Jets =  h_UnfOverMC_Jets->GetBinContent(bin_Jets);  
	
 	// //1D reco distribution built not for all reco events, but only events gen&reco
 	// theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_W_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERjets , theWeight*w_Jets);
 	// theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_W_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4, nCentralJERjets, theWeight*w_Jets);
 	// theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERCentralSmear_W_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nCentralJERjets,ngencentraljets, theWeight*w_Jets); 
 	// theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERCentralSmear_W_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nCentralJERjets,ngencentraljets, theWeight*w_Jets);
	
 	// //Matrix and reco distribution weighted for the ratio between the unfolded data and the generator level information in order to 
 	// //compute the relative systematic uncertainty (an early unfolding is required)
 	// string UnfOverMC_Mass = "ZZTo"+decay+"_Mass_Ratio";
 	// h_UnfOverMC_Mass = (TH1*) UnfOverMC_Pow->Get(UnfOverMC_Mass.c_str()); 
 	// int bin_Mass = h_UnfOverMC_Mass->FindBin(m4L_gen);
 	// float w_Mass =  h_UnfOverMC_Mass->GetBinContent(bin_Mass);
	
 	// //1D reco distribution built not for all reco events, but only events gen&reco.
 	// theHistograms.fill(std::string("ZZTo")+decay+"_Mass_W_01", std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L,theWeight*w_Mass);
 	// theHistograms.fill(std::string("ZZTo")+decay+"_Mass_W_"+sample, std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L,theWeight*w_Mass);
 	// theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_W_"+sample, std::string("Response matrix invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, Xbins, m4L ,m4L_gen , theWeight*w_Mass); 
 	// theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_W_01", std::string("Response matrix invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, Xbins, m4L ,m4L_gen , theWeight*w_Mass);
	
 	// // cout << decay.c_str() << " " << m4L << " " << m4L_gen << " " << w_Mass << endl;
	
    } //end if(topology.test(0))
  
      // // if the event is reconstructed but not generated as signal, put the w_Mass=1
      // else{
      // 	theHistograms.fill(std::string("ZZTo")+decay+"_Mass_W_01", std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L,theWeight);
      // 	theHistograms.fill(std::string("ZZTo")+decay+"_Mass_W_"+sample, std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L,theWeight);
      // 	theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_W_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERjets , theWeight);
      // 	theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_W_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4, nCentralJERjets, theWeight); 
      // 	theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERCentralSmear_W_"+sample, "m_{jj}", Xbins_mjj,centralMjj,theWeight);
      // 	theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_JERCentralSmear_W_01", "m_{jj}", Xbins_mjj, centralMjj,theWeight);
      // 	theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERCentralSmear_W_"+sample, "#Delta#eta_{jj}", Xbins_deta,centralDeta,theWeight);
      // 	theHistograms.fill(std::string("ZZTo")+decay+"_Deta_JERCentralSmear_W_01", "#Delta#eta_{jj}", Xbins_deta, centralDeta,theWeight);
      // }
  
  
       
  // //**************************************************For VBS Analysis*******************************************************************************************

  // if(nCentralJERjets>=2){
    
  //   float deta = fabs(CentralJER_jets->at(0).eta() - CentralJER_jets->at(1).eta());
  //   float mjj =  (CentralJER_jets->at(0).p4() + CentralJER_jets->at(1).p4()).M();
  // //theHistograms.fill(std::string("DeltaEtaJJ_ZZTo") +decay+"_01", std::string("#Delta #eta(j,j) between the two most energetic jets ZZ#rightarrow ") +decay,  50, 0, 8, fabs(CentralJER_jets->at(0).eta() - CentralJER_jets->at(1).eta()), theWeight);
  // //  theHistograms.fill(std::string("mJJ")+decay+"_01", "m_{jj}", 100,0,3000,mjj, theWeight);  
    
  //   theHistograms.fill(std::string("DeltaEtaJJ_ZZTo") +decay+"_01", std::string("#Delta #eta(j,j) between the two most energetic jets ZZ#rightarrow ") +decay,  Xbins_deta, deta, theWeight); 
  //   theHistograms.fill(std::string("mJJ_ZZTo")+decay+"_01", "m_{jj}",  Xbins_mjj,mjj, theWeight); 
    
    
  //   if(CentralJER_jets->at(0).pt()>100){
  //     theHistograms.fill(std::string("DeltaEtaJJ_pt100_ZZTo") +decay+"_01", std::string("#Delta #eta(j,j) between the two most energetic jets ZZ#rightarrow ") +decay,  Xbins_deta, deta, theWeight); 
  //     theHistograms.fill(std::string("mJJ_pt100_ZZTo")+decay+"_01", "m_{jj}",  Xbins_mjj,mjj, theWeight); 
  //   }
    
  //   if(met->pt()<60){
  //     theHistograms.fill(std::string("DeltaEtaJJ_met60_ZZTo") +decay+"_01", std::string("#Delta #eta(j,j) between the two most energetic jets ZZ#rightarrow ") +decay,  Xbins_deta, deta, theWeight); 
  //     theHistograms.fill(std::string("mJJ_met60_ZZTo")+decay+"_01", "m_{jj}",  Xbins_mjj,mjj, theWeight); 
    
  //     if(CentralJER_jets->at(0).pt()>100){
  // 	theHistograms.fill(std::string("DeltaEtaJJ_pt100_met60_ZZTo") +decay+"_01", std::string("#Delta #eta(j,j) between the two most energetic jets ZZ#rightarrow ") +decay,  Xbins_deta, deta, theWeight); 
  // 	theHistograms.fill(std::string("mJJ_pt100_met60_ZZTo")+decay+"_01", "m_{jj}",  Xbins_mjj,mjj, theWeight); 
  //     }
  //   }
    
  //   if(CentralJER_jets->at(1).pt()>100 && CentralJER_jets->at(1).pt()>70 && met->pt()<60 ){
    
  //     theHistograms.fill(std::string("DeltaEtaJJ_pt100_met60_pt70_ZZTo") +decay+"_01", std::string("#Delta #eta(j,j) between the two most energetic jets ZZ#rightarrow ") +decay,  Xbins_deta, deta, theWeight); 
  //     theHistograms.fill(std::string("mJJ_pt100_met60_pt70_ZZTo")+decay+"_01", "m_{jj}",  Xbins_mjj,mjj, theWeight);
    
  //     if(mjj > 500)   theHistograms.fill(std::string("DeltaEtaJJ_pt100_met60_pt70_mjj500_ZZTo") +decay+"_01", std::string("#Delta #eta(j,j) between the two most energetic jets ZZ#rightarrow ") +decay,  Xbins_deta, deta, theWeight); 
  //     if(mjj > 200)   theHistograms.fill(std::string("DeltaEtaJJ_pt100_met60_pt70_mjj200_ZZTo") +decay+"_01", std::string("#Delta #eta(j,j) between the two most energetic jets ZZ#rightarrow ") +decay,  Xbins_deta, deta, theWeight); 
  //     if(deta > 2) theHistograms.fill(std::string("mJJ_pt100_met60_pt70_deta2_ZZTo")+decay+"_01", "m_{jj}",  Xbins_mjj,mjj, theWeight);
    
  //   }
    
  // }

    //*************************************************************************************************************************************************************
  }//end if(genCategory != -1)    
 
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  if(region_ == phys::CR3P1F || region_ == phys::CR2P2F) {
    theHistograms.fill(std::string("ZZTo")+decay+"_Mass"+"_FRVar", std::string("Var From FR Invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, m4L,ZZ->fakeRateSFVar());
    theHistograms.fill(std::string("ZZTo")+decay+"_Jets"+"_FRVar", std::string("Var From FR Invariant mass of ZZ_{1}#rightarrow ")+decay ,4,0,4,njets,ZZ->fakeRateSFVar());
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
 
  // UpJESData_jets->clear();
  // DownJESData_jets->clear();
  // UpJESData_centraljets->clear();
  // DownJESData_centraljets->clear();
  
  // foreach(const phys::Jet &jet, *pjets){
  //   double jetPt = 0;
  //   jetPt = jet.pt();
    
  //   //JES correction: Up and down velues used to assess systematic uncertainty on jet energy resolution
  //   double newJetPtJESData_up =0;  
  //   double newJetPtJESData_down =0;
    
  //   newJetPtJESData_up = jetPt*(1+jet.uncOnFourVectorScale());
  //   newJetPtJESData_down = jetPt*(1-jet.uncOnFourVectorScale());
    
  //   if(newJetPtJESData_up > 30) UpJESData_jets->push_back(jet); 
  //   if(newJetPtJESData_down > 30) DownJESData_jets->push_back(jet);
  //   if(newJetPtJESData_up > 30 && jet.eta()<2.4) UpJESData_centraljets->push_back(jet); 
  //   if(newJetPtJESData_down > 30 && jet.eta()<2.4) DownJESData_centraljets->push_back(jet);

  //   //cout << decay << " " << jetPt << " " <<   newJetPtJESData_up << " " << newJetPtJESData_down << " " << njets << endl;
    
  // }
  
  // nUpJESDatajets = UpJESData_jets->size();
  // nDownJESDatajets = DownJESData_jets->size();
  // nUpJESDatacentraljets = UpJESData_centraljets->size();
  // nDownJESDatacentraljets = DownJESData_centraljets->size();
  // if (nUpJESDatajets>3) nUpJESDatajets=3;
  // if (nDownJESDatajets>3) nDownJESDatajets=3;
  // if (nUpJESDatacentraljets>3) nUpJESDatacentraljets=3;
  // if (nDownJESDatacentraljets>3) nDownJESDatacentraljets=3;
  
  // stable_sort(UpJESData_jets->begin(), UpJESData_jets->end(), PtComparator());
  // stable_sort(DownJESData_jets->begin(), DownJESData_jets->end(), PtComparator());
  // stable_sort(UpJESData_centraljets->begin(), UpJESData_centraljets->end(), PtComparator());
  // stable_sort(DownJESData_centraljets->begin(), DownJESData_centraljets->end(), PtComparator());
  

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
    
    
    //Loop on all reco jets
    foreach(const phys::Jet &jet, *pjets){
      
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
      
      if(newJetPtJER > 30) CentralJER_jets->push_back(jet);
      if(newJetPtJER_up > 30) UpJER_jets->push_back(jet); 
      if(newJetPtJER_down > 30) DownJER_jets->push_back(jet);
      if(newJetPtJER > 30 && jet.eta()<2.4) CentralJER_centraljets->push_back(jet);
      if(newJetPtJER_up > 30 && jet.eta()<2.4) UpJER_centraljets->push_back(jet); 
      if(newJetPtJER_down > 30 && jet.eta()<2.4) DownJER_centraljets->push_back(jet);
      
      //JES correction: Up and down velues used to assess systematic uncertainty on jet energy resolution
      double newJetPtJES_up =0;  
      double newJetPtJES_down =0;
      
      //cout << jet.uncOnFourVectorScale() << endl;
      
      newJetPtJES_up = newJetPtJER*(1+jet.uncOnFourVectorScale());
      newJetPtJES_down = newJetPtJER*(1-jet.uncOnFourVectorScale());
      
      if(newJetPtJES_up > 30) UpJES_jets->push_back(jet); 
      if(newJetPtJES_down > 30) DownJES_jets->push_back(jet);
      if(newJetPtJES_up > 30 && jet.eta()<2.4 ) UpJES_centraljets->push_back(jet); 
      if(newJetPtJES_down > 30 &&  jet.eta()<2.4) DownJES_centraljets->push_back(jet);
      
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
    
    // ZZplots();   // ZZ --> 4l 
    ZZplots(52,e); // ZZ --> 4m
    ZZplots(48,e); // ZZ --> 2e2m
    ZZplots(44,e); // ZZ --> 4e
    
  }
  
void ZZRecoAnalyzer::begin() {
  
  UpJESData_jets  = new std::vector<phys::Jet>();
  DownJESData_jets  = new std::vector<phys::Jet>();
  UpJESData_centraljets  = new std::vector<phys::Jet>();
  DownJESData_centraljets  = new std::vector<phys::Jet>();
  
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
 
  UnfOverMC = new TFile("macros/UnfoldingMacros/UnfoldFolder/Ratio_UnfoldedDataOverGenMC.root"); 
  UnfOverMC_Pow = new TFile("macros/UnfoldingMacros/UnfoldFolder_Pow/Ratio_UnfoldedDataOverGenMC.root"); 
  nentries = tree()->GetEntries("ZZCand.passFullSel_"); 
  Xbins += 100,200,250,300,350,400,500,600,800;
  Xbins_deta += 0,2.4,4.7;
  Xbins_mjj += 0.,200,800;
  m4L = 0;
  m4L_gen = 0; 
  ngenjets = 0;
  ngencentraljets =0; 
  mjj_gen = 0;
  deta_gen = 0; 
  mjj_gen_cj = 0;
  deta_gen_cj = 0; 
}

void ZZRecoAnalyzer::end( TFile &) {
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
