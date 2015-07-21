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
  
  
  //////////////////////////////////////////////DATA ONLY///////////////////////////////////////////////////////


  // //1D Reco nJet Distributions (no JER smearing, for data)
  // theHistograms.fill(std::string("ZZTo")+decay+"_Jets_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight); 
  // theHistograms.fill(std::string("ZZTo")+decay+"_Jets_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight); 
  
  // //for data only (JES uncertainty) 
  // UpJESData_jets  = new std::vector<phys::Jet>();
  // DownJESData_jets  = new std::vector<phys::Jet>();
  // UpJESData_jets->clear();
  // DownJESData_jets->clear();
  
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

  //   //cout << decay << " " << jetPt << " " <<   newJetPtJESData_up << " " << newJetPtJESData_down << " " << njets << endl;
    
  // }
  
  // Int_t nUpJESDatajets = UpJESData_jets->size();
  // Int_t nDownJESDatajets = DownJESData_jets->size();
  // if (nUpJESDatajets>3) nUpJESDatajets=3;
  // if (nDownJESDatajets>3) nDownJESDatajets=3;
  
  // theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESDataUpSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJESDatajets , theWeight);
  // theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESDataDownSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJESDatajets , theWeight);
  // theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESDataUpSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nUpJESDatajets , theWeight);
  // theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JESDataDownSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nDownJESDatajets , theWeight);
  
  ////////////////////////////////////////////////////DATA AND MONTECARLO/////////////////////////////////////////////////

  //1D Reco Mass Distributions
  theHistograms.fill(std::string("ZZTo")+decay+"_Mass_01", std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, ZZ->mass(),theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_Mass_"+sample, std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, ZZ->mass(),theWeight);
  

  /////////////////////////////////////////////////////MONTECARLO ONLY//////////////////////////////////////////////////////

  //Smearing jet pt (Up and Down distributions just for systematic uncertainty estimate, Central for the standard analysis) without requiring it is signal (1D distributions made of ALL reco events)
  CentralJER_jets  = new std::vector<phys::Jet>();
  UpJER_jets  = new std::vector<phys::Jet>();
  DownJER_jets  = new std::vector<phys::Jet>();
  
  UpJES_jets  = new std::vector<phys::Jet>();
  DownJES_jets  = new std::vector<phys::Jet>();
  
  CentralJER_jets->clear();
  UpJER_jets->clear();
  DownJER_jets->clear();
  
  UpJES_jets->clear();
  DownJES_jets->clear();
  
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

    //JES correction: Up and down velues used to assess systematic uncertainty on jet energy resolution
    double newJetPtJES_up =0;  
    double newJetPtJES_down =0;

    //cout << jet.uncOnFourVectorScale() << endl;

    newJetPtJES_up = newJetPtJER*(1+jet.uncOnFourVectorScale());
    newJetPtJES_down = newJetPtJER*(1-jet.uncOnFourVectorScale());
  
    if(newJetPtJES_up > 30) UpJES_jets->push_back(jet); 
    if(newJetPtJES_down > 30) DownJES_jets->push_back(jet);
    
  }
  
  Int_t nCentralJERjets = CentralJER_jets->size();
  Int_t nUpJERjets = UpJER_jets->size();
  Int_t nDownJERjets = DownJER_jets->size();
  
  if (nCentralJERjets>3) nCentralJERjets=3;
  if (nUpJERjets>3) nUpJERjets=3;
  if (nDownJERjets>3) nDownJERjets=3;

  Int_t nUpJESjets = UpJES_jets->size();
  Int_t nDownJESjets = DownJES_jets->size();
  if (nUpJESjets>3) nUpJESjets=3;
  if (nDownJESjets>3) nDownJESjets=3;

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

  //if MC gen (for response matrices only)
  if (genCategory !=-1){
    if(topology.test(0)){
      
      ngenjets =  genJets->size(); 
      if (ngenjets>3) ngenjets=3;
      m4L_gen = sqrt((genVBParticles->at(0).p4()+genVBParticles->at(1).p4())*(genVBParticles->at(0).p4()+genVBParticles->at(1).p4()));

      //Response Matrix Mass (Reco&Gen) 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_"+sample, std::string("Response matrix invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, Xbins, ZZ->mass() ,m4L_gen , theWeight); 
      theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_01", std::string("Response matrix invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, Xbins, ZZ->mass() ,m4L_gen , theWeight);
      
      // //Response Matrix nJets (Reco&Gen) - No JER smearing
      // theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4, njets,ngenjets, theWeight); 
      // theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4, njets,ngenjets, theWeight); 
      
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
      
      // //Matrix and reco distribution weighted for the ratio between the unfolded data and the generator level information in order to compute the relative systematic uncertainty
      // //distributions to evaluate data/MC systematic uncertainty
      // string UnfOverMC_Jets = "ZZTo"+decay+"_Jets_Ratio";
      // h_UnfOverMC_Jets = (TH1*) UnfOverMC->Get(UnfOverMC_Jets.c_str()); 
      // int bin_Jets = h_UnfOverMC_Jets->FindBin(nCentralJERjets);
      // float w_Jets =  h_UnfOverMC_Jets->GetBinContent(bin_Jets);
       
      // //1D reco distribution built not for all reco events, but only events gen&reco
      // theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_W_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERjets , theWeight*w_Jets);
      // theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_W_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4, nCentralJERjets, theWeight*w_Jets);
      // theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERCentralSmear_W_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_Jets); 
      // theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_JERCentralSmear_W_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4,nCentralJERjets,ngenjets, theWeight*w_Jets);
      
      // //Matrix and reco distribution weighted for the ratio between the unfolded data and the generator level information in order to compute the relative systematic uncertainty (an early unfolding is required)
      // string UnfOverMC_Mass = "ZZTo"+decay+"_Mass_Ratio";
      // h_UnfOverMC_Mass = (TH1*) UnfOverMC->Get(UnfOverMC_Mass.c_str()); 
      // int bin_Mass = h_UnfOverMC_Mass->FindBin(m4L_gen);
      // float w_Mass =  h_UnfOverMC_Mass->GetBinContent(bin_Mass);
      
      // //1D reco distribution built not for all reco events, but only events gen&reco.
      // theHistograms.fill(std::string("ZZTo")+decay+"_Mass_W_01", std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, ZZ->mass(),theWeight*w_Mass);
      // theHistograms.fill(std::string("ZZTo")+decay+"_Mass_W_"+sample, std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, ZZ->mass(),theWeight*w_Mass);
      // theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_W_"+sample, std::string("Response matrix invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, Xbins, ZZ->mass() ,m4L_gen , theWeight*w_Mass); 
      // theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_W_01", std::string("Response matrix invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, Xbins, ZZ->mass() ,m4L_gen , theWeight*w_Mass);
      
      // cout << decay.c_str() << " " << ZZ->mass() << " " << m4L_gen << " " << w_Mass << endl;
      
    }

    // if the event is reconstructed but not generated as signal, put the w_Mass=1
    // else{
    //  theHistograms.fill(std::string("ZZTo")+decay+"_Mass_W_01", std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, ZZ->mass(),theWeight);
    //  theHistograms.fill(std::string("ZZTo")+decay+"_Mass_W_"+sample, std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, ZZ->mass(),theWeight);
    //   theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_W_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,nCentralJERjets , theWeight);
    //   theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_W_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4, nCentralJERjets, theWeight);
    // }
  }
  
  // stable_sort(CentralJER_jets->begin(), CentralJER_jets->end(), PtComparator());
  // if(nCentralJERjets>=2){
  //   theHistograms.fill(std::string("DeltaEtaJJ_ZZTo") +decay+"_01", std::string("#Delta #eta(j,j) between the two most energetic jets ZZ#rightarrow ") +decay,  50, 0, 8, fabs(CentralJER_jets->at(0).eta() - CentralJER_jets->at(1).eta()), theWeight); 
  //   theHistograms.fill(std::string("mJJ")+decay+"_01", "m_{jj}",  100, 0, 3000, (CentralJER_jets->at(0).p4() + CentralJER_jets->at(1).p4()).M(), theWeight); 
  // }
  
  
  if(region_ == phys::CR3P1F || region_ == phys::CR2P2F) {
    theHistograms.fill(std::string("ZZTo")+decay+"_Mass"+"_FRVar", std::string("Var From FR Invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, ZZ->mass(),ZZ->fakeRateSFVar());
    theHistograms.fill(std::string("ZZTo")+decay+"_Jets"+"_FRVar", std::string("Var From FR Invariant mass of ZZ_{1}#rightarrow ")+decay ,4,0,4,njets,ZZ->fakeRateSFVar());
  }
}

//Smearing function
double ZZRecoAnalyzer::JER_PtSmear(double pt, double width)
{
  double ptsmear= gRandom->Gaus(pt,width);    
  return ptsmear;
}

void ZZRecoAnalyzer::analyze(){
  
  e++;
  
  // ZZplots();   // ZZ --> 4l
  ZZplots(52,e); // ZZ --> 4m
  ZZplots(48,e); // ZZ --> 2e2m
  ZZplots(44,e); // ZZ --> 4e
  
 }

void ZZRecoAnalyzer::begin() {
  //UnfOverMC = new TFile("macros/UnfoldingMacros/UnfoldFolderRatio_UnfoldedDataOverGenMC.root");
  nentries = tree()->GetEntries("ZZCand.passFullSel_");
  Xbins += 100,200,250,300,350,400,500,600,800;
  m4L_gen = 0;
  ngenjets = 0;
}

void ZZRecoAnalyzer::end( TFile &) {
   cout<<theMCInfo.analyzedEvents()<< " " << e <<endl;
  
  if(region_ == phys::CR3P1F || region_ == phys::CR2P2F) {
    vector<std::string>  FinalState = {"4m","4e","2e2m"};
    
    for (std::vector<std::string>::iterator it = FinalState.begin() ; it != FinalState.end(); ++it){
      
      TH1 *hvar =  new TH1F();
      hvar =  theHistograms.get(("ZZTo"+*it+"_Mass_FRVar").c_str());
      
      TH1 *h =  new TH1F();
      h =  theHistograms.get(("ZZTo"+*it+"_Mass_01").c_str());
      
      if(!h) continue;
      for(int i = 1; i<=h->GetNbinsX();i++){
	
	Float_t Err = h->GetBinError(i);
	h->SetBinError(i,sqrt(Err*Err+hvar->GetBinContent(i)));
      }
    }
  }  
}
