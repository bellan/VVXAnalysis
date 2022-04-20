#include "VVXAnalysis/TreeAnalysis/interface/VVGammaAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include "VVXAnalysis/Commons/interface/GenVBHelper.h"

#include "TTree.h"

#include <algorithm>  // std::move
#include <fstream>

#include <boost/foreach.hpp>
// #include <boost/range/join.hpp>  // boost::join
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;
using std::vector;
using namespace colour;

using namespace phys;
using namespace physmath;

#define CUT_LAYOUT 6,0.5,6.5
#define DEBUG

void VVGammaAnalyzer::begin(){
  cout<<'\n';
  for(char i=0; i<25; ++i) cout<<'-';
  cout<<" Start of VVGammaAnalyzer ";
  for(char i=0; i<25; ++i) cout<<'-';
  cout<<'\n';
  /*
    TFile* photonSFfile = TFile::Open("../data/2018_PhotonsMedium.root", "r");
    if(! (photonSFfile && photonSFfile->IsOpen()) ){
    cout<<"Warning: photon SF file not found\n";
    } else{
    photonSFhist = (TH2F*) photonSFfile->Get("EGamma_SF2D")->Clone();
    photonSFfile->Close();
    }
    delete photonSFfile;*/

  initCherryPick();

  size_t digits = std::to_string(tree()->GetEntries()).length();
  std::string spaces( digits, ' ' );
  // cout<<"Analyzed:\t"<<spaces<<'/'<<tree()->GetEntries()<<std::flush;
  cout<<'\n'; //TEMP

  return;
}


Int_t VVGammaAnalyzer::cut() {
  evtN_++; evtNInReg_[region_]++; evtWInReg_[region_] += theWeight;
  // cout<<"\r\t\t"<<evtN_;  //TEMP
  cout << regionType(region_) << ':' << run << ':' << lumiBlock << ':' << event << '\n';
  if(cherrypickEvt()){
    theHistograms->fill("cherry_ZZ_mass"    , "Events in {UL#backslashLegacy}: m_{ZZ};m_{ZZ} [GeV/c^{2}]", 25,0.,500., ZZ->mass()                      , theWeight);
    theHistograms->fill("cherry_ZZ_goodLept", "Events in {UL#backslashLegacy}: # good leptons"           , 5,-0.5,4.5, ZZ->numberOfGoodGrandDaughters(), theWeight);
    theHistograms->fill("cherry_ZZ_badLept" , "Events in {UL#backslashLegacy}: # bad leptons"            , 5,-0.5,4.5, ZZ->numberOfBadGrandDaughters() , theWeight);
  }

  theHistograms->fill("AAA_cuts"  , "Cuts weighted"  , CUT_LAYOUT, 1, theWeight);
  theHistograms->fill("AAA_cuts_u", "Cuts unweighted", CUT_LAYOUT, 1);
	
  // Cleanup
  leptons_    ->clear();
  kinPhotons_ ->clear();
  goodPhotons_->clear();
  
  bool haveZVlep = false;
  bool have2l2j = false;
  bool haveGoodPhoton = false;
	
  theHistograms->fill("POG_leptons", "n leptons", 7,-0.5,6.5, electrons->size()+muons->size(), theWeight);
  
  
  // Contruct vector with leptons from dibosons
  if(ZZ && ZZ->pt() > 1.)
    leptons_->insert(leptons_->end(), {
	  ZZ->first().daughterPtr(0), 
	  ZZ->first().daughterPtr(1), 
	  ZZ->second().daughterPtr(0), 
	  ZZ->second().daughterPtr(1)
	  });
  else if(ZW && ZW->pt() > 1.)
    leptons_->insert(leptons_->end(), {
	  ZW->first().daughterPtr(0), 
	  ZW->first().daughterPtr(1), 
	  ZW->second().daughterPtr(0)
	  });
  else if(ZL && ZL->first.pt() > 1.){
    leptons_->insert(leptons_->end(), {
	  ZL->first.daughterPtr(0),
	  ZL->first.daughterPtr(1),
	  &(ZL->second)
      });
  }
  
  unsigned int nGoodLeptons = 0;
  for(const Lepton* l : *leptons_)
    if(l->passFullSel()) nGoodLeptons++;
  theHistograms->fill("goodLeptons_N", "Number of leptons passing full selection", 5, -0.5, 4.5, nGoodLeptons, theWeight);
	
  vector<Photon> kinPhotons_EScale_Up;
  vector<Photon> kinPhotons_EScale_Down;
  vector<Photon> goodPhotons_EScale_Up;
  vector<Photon> goodPhotons_EScale_Down;
  for(auto ph : *photons){
    //Pixel seed and electron veto
    if(ph.hasPixelSeed() || !ph.passElectronVeto()) continue;
		
    //Kinematic selection
    // if(ph.pt() < 20) continue;
    float ph_aeta = fabs(ph.eta());
    if(ph_aeta > 2.4) continue;
    if(ph_aeta > 1.4442 && ph_aeta < 1.566) continue;
		
    //Electrons and muons matching
    bool match = false;
    for(const Lepton* lep : *leptons_){
      if(deltaR(ph,*lep) < 0.3){
	match = true;
	break;
      }
    }
    if(match) continue;
    
    if(ph.pt() * ph.energyScaleUp()/ph.e() > 20){
      kinPhotons_EScale_Up.push_back(ph);
      if(ph.cutBasedIDLoose())
	goodPhotons_EScale_Up.push_back(ph);
    }
    if(ph.pt() * ph.energyScaleDown()/ph.e() > 20){
      kinPhotons_EScale_Down.push_back(ph);
      if(ph.cutBasedIDLoose())
	goodPhotons_EScale_Down.push_back(ph);
    }

    if(ph.pt() < 20) continue;
    
    kinPhotons_->push_back(ph);
    if(ph.cutBasedIDLoose())
      goodPhotons_->push_back(ph);
    
  }
  
  if(kinPhotons_->size() > 1)            theHistograms->fill("ph_eScale_count", "Num photons pass selection", 6,-0.5,5.5, 0, theWeight);
  if(kinPhotons_EScale_Down.size() > 1)  theHistograms->fill("ph_eScale_count", "Num photons pass selection", 6,-0.5,5.5, 1, theWeight);
  if(kinPhotons_EScale_Up.size() > 1)    theHistograms->fill("ph_eScale_count", "Num photons pass selection", 6,-0.5,5.5, 2, theWeight);
  if(goodPhotons_->size() > 1)           theHistograms->fill("ph_eScale_count", "Num photons pass selection", 6,-0.5,5.5, 3, theWeight);
  if(goodPhotons_EScale_Down.size() > 1) theHistograms->fill("ph_eScale_count", "Num photons pass selection", 6,-0.5,5.5, 4, theWeight);
  if(goodPhotons_EScale_Up.size() > 1)   theHistograms->fill("ph_eScale_count", "Num photons pass selection", 6,-0.5,5.5, 5, theWeight);
  
  // Basic histos
  
  for(size_t i = 0; i < photons->size() ; i++){
    const Photon& ph = photons->at(i);
    theHistograms->fill("ph_E"                , "E;[GeV]"                             , 50,0.,500., ph.e()                           , theWeight);
    theHistograms->fill("ph_E_m_ScaleUp"      , "energyScaleUp - E;#DeltaE[GeV]"      , 60,-1.,5.,  ph.energyScaleUp() - ph.e()      , theWeight);
    theHistograms->fill("ph_E_m_ScaleDown"    , "E - energyScaleDown;#DeltaE[GeV]"    , 60,-1.,5.,  ph.e() - ph.energyScaleDown()    , theWeight);
    theHistograms->fill("ph_E_m_ScaleStatUp"  , "energyScaleStatUp - E;#DeltaE[GeV]"  , 60,-1.,5.,  ph.energyScaleStatUp() - ph.e()  , theWeight);
    theHistograms->fill("ph_E_m_ScaleStatDown", "E - energyScaleStatDown;#DeltaE[GeV]", 60,-1.,5.,  ph.e() - ph.energyScaleStatDown(), theWeight);
    theHistograms->fill("ph_E_m_ScaleSystUp"  , "energyScaleSystUp - E;#DeltaE[GeV]"  , 60,-1.,5.,  ph.energyScaleSystUp() - ph.e()  , theWeight);
    theHistograms->fill("ph_E_m_ScaleSystDown", "E - energyScaleSystDown;#DeltaE[GeV]", 60,-1.,5.,  ph.e() - ph.energyScaleSystDown(), theWeight);
    theHistograms->fill("ph_E_m_ScaleGainUp"  , "energyScaleGainUp - E#;DeltaE[GeV]"  , 60,-1.,5.,  ph.energyScaleGainUp() - ph.e()  , theWeight);
    theHistograms->fill("ph_E_m_ScaleGainDown", "E - energyScaleGainDown;#DeltaE[GeV]", 60,-1.,5.,  ph.e() - ph.energyScaleGainDown(), theWeight);
    theHistograms->fill("ph_E_m_ScaleEtUp"    , "energyScaleEtUp - E;#DeltaE[GeV]"    , 60,-1.,5.,  ph.energyScaleEtUp() - ph.e()    , theWeight);
    theHistograms->fill("ph_E_m_ScaleEtDown"  , "E - energyScaleEtDown;#DeltaE[GeV]"  , 60,-1.,5.,  ph.e() - ph.energyScaleEtDown()  , theWeight);
  }
  
  baseHistos_cut();
  // PKU_comparison();
  theHistograms->fill("goodPhoton_N", "n good #gamma | ZZ/WZ exists", 7,-0.5,6.5, goodPhotons_->size(), theWeight);
  
  
  // ----- BASELINE SELECTION -----
  // ----- Cut1: require at least a ZZ or a WZ candidate
  haveZVlep = (ZZ && ZZ->pt() > 1.) || (ZW && ZW->pt() > 1.);	
  if(haveZVlep){
    theHistograms->fill("AAA_cuts"  , "Cuts weighted"  , CUT_LAYOUT, 2, theWeight);
    theHistograms->fill("AAA_cuts_u", "Cuts unweighted", CUT_LAYOUT, 2);
  }
	
  have2l2j = (muons->size()+electrons->size()==2) && (jets->size()==2 || jetsAK8->size()>=1);
  if(have2l2j){
    theHistograms->fill("AAA_cuts"  , "Cuts weighted"  , CUT_LAYOUT, 3, theWeight);
    theHistograms->fill("AAA_cuts_u", "Cuts unweighted", CUT_LAYOUT, 3);
  }
  
  // ----- Cut2: Require at least 1 loose photon with pt > 20 GeV
  haveGoodPhoton = goodPhotons_->size() >= 1;
  if(haveGoodPhoton){
    theHistograms->fill("AAA_cuts"  , "Cuts weighted"  , CUT_LAYOUT, 4, theWeight);
    theHistograms->fill("AAA_cuts_u", "Cuts unweighted", CUT_LAYOUT, 4);
  }
  
  //return (haveGoodPhoton ? 1 : -1);
  return (haveZVlep ? 1 : -1); //TEMP
  
  if(!haveZVlep) 
    return -1;
  else if(!haveGoodPhoton)
    return -1;
  else
    return 1;
}

void VVGammaAnalyzer::analyze(){
  analyzedNInReg_[region_]++; analyzedWInReg_[region_] += theWeight;
  theHistograms->fill("AAA_cuts"  , "Cuts weighted"  , CUT_LAYOUT, 5, theWeight);
  theHistograms->fill("AAA_cuts_u", "Cuts unweighted", CUT_LAYOUT, 5);
  
  if(goodPhotons_->size() >= 1 && goodPhotons_->front().cutBasedIDLoose()){
    theHistograms->fill("AAA_cuts"  , "Cuts weighted"  , CUT_LAYOUT, 6, theWeight);
    theHistograms->fill("AAA_cuts_u", "Cuts unweighted", CUT_LAYOUT, 6);
  }
  
  unsigned int nEl(0), nMu(0);
  for(const Lepton* l : *leptons_){
    switch( abs(l->id()) ){
    case 11:
      nEl++;
      break;
    case 13:
      nMu++;
      break;
    default:
      cout<<"Error: found lept from ZZ/ZW/ZL with ID: "<<l->id()<<'\n';
    }
  }
  std::string channel_reco = Form("%ue%um", nEl, nMu);
  
  LeptonFakeRate();
  PhotonFakeRate();
  
  Photon* thePh = nullptr;
  if(goodPhotons_->size() >= 1) thePh = &goodPhotons_->at(0);
  
  // SR: loose photon ID - CR: kin selection && !loose photon ID
  // Photon& thePhoton = goodPhotons_->front();
  // theHistograms->fill("photon_pt", "p_{t}^{#gamma};GeV/c", 50,0.,200., thePhoton.pt(), theWeight);
  bool four_lep  = (region_ >= SR4P && region_ <= CR4P_1F) && ZZ && ZZ->pt() > 1.;
  bool three_lep = ((region_>=CR110 && region_<= CR000) || region_ == SR3P ) && ZW && ZW->pt() > 1.;
  bool LFR_lep   = region_ == CRLFR && ZL && ZL->first.pt() > 1.;
  
  if(four_lep){
    if(thePh)
      theHistograms->fill("ZZG_mass_"+channel_reco,"m_{4l#gamma};[GeV/c^{2}]",25,0.,500,(ZZ->p4()+thePh->p4()).M(),theWeight);
    
    theHistograms->fill("ZZ_mass"               , "m_{4l};GeV/c^{2}", 25,0.,500. , ZZ->mass()                   , theWeight);
    theHistograms->fill("Z0_mass"               , "m_{Z0};GeV/c^{2}", 35,55.,125., ZZ->first().mass()           , theWeight);
    theHistograms->fill("Z1_mass"               , "m_{Z1};GeV/c^{2}", 35,55.,125., ZZ->second().mass()          , theWeight);
    theHistograms->fill("ZZ_pt"                 , "p_{t,ZZ};GeV/c"  , 20,0.,400. , ZZ->pt()                     , theWeight);
    theHistograms->fill("Z0_l0_pt"              , "p_{t,l00};GeV/c" , 20,0.,400. , ZZ->first().daughter(0).pt() , theWeight);
    theHistograms->fill("Z0_l1_pt"              , "p_{t,l01};GeV/c" , 20,0.,400. , ZZ->first().daughter(1).pt() , theWeight);
    theHistograms->fill("Z1_l0_pt"              , "p_{t,l10};GeV/c" , 20,0.,400. , ZZ->second().daughter(0).pt(), theWeight);
    theHistograms->fill("Z1_l1_pt"              , "p_{t,l11};GeV/c" , 20,0.,400. , ZZ->second().daughter(1).pt(), theWeight);

    theHistograms->fill("ZZ_mass_" +channel_reco, "m_{4l};GeV/c^{2}", 25,0.,500. , ZZ->mass()                   , theWeight);
    theHistograms->fill("Z0_mass"  +channel_reco, "m_{Z0};GeV/c^{2}", 35,55.,125., ZZ->first().mass()           , theWeight);
    theHistograms->fill("Z1_mass_" +channel_reco, "m_{Z1};GeV/c^{2}", 35,55.,125., ZZ->second().mass()          , theWeight);
    theHistograms->fill("ZZ_pt_"   +channel_reco, "p_{t,ZZ};GeV/c"  , 20,0.,400. , ZZ->pt()                     , theWeight);
    theHistograms->fill("Z0_l0_pt_"+channel_reco, "p_{t,l00};GeV/c" , 20,0.,400. , ZZ->first().daughter(0).pt() , theWeight);
    theHistograms->fill("Z0_l1_pt_"+channel_reco, "p_{t,l01};GeV/c" , 20,0.,400. , ZZ->first().daughter(1).pt() , theWeight);    
    theHistograms->fill("Z1_l0_pt_"+channel_reco, "p_{t,l10};GeV/c" , 20,0.,400. , ZZ->second().daughter(0).pt(), theWeight);    
    theHistograms->fill("Z1_l1_pt_"+channel_reco, "p_{t,l11};GeV/c" , 20,0.,400. , ZZ->second().daughter(1).pt(), theWeight);
  }
  else if(three_lep){
    if(thePh)
      theHistograms->fill("ZWG_massT_"+channel_reco,"m_{T,3l#gamma};[GeV/c^{2}]",25,0.,500,(ZZ->p4()+thePh->p4()).Mt(),theWeight);
    
    theHistograms->fill("ZW_massT"              , "m_{T,3l};GeV/c^{2}", 25,0.,500. , ZW->p4().Mt()                , theWeight);
    theHistograms->fill("Z_mass"                , "m_{Z};GeV/c^{2}"   , 35,55.,125., ZW->first().mass()           , theWeight);
    theHistograms->fill("W_massT"               , "m_{T,W};GeV/c^{2}" , 35,55.,125., ZW->second().p4().Mt()       , theWeight);
    theHistograms->fill("ZW_pt"                 , "p_{t,ZW};GeV/c"    , 20,0.,400. , ZW->pt()                     , theWeight);
    theHistograms->fill("Z_l0_pt"               , "p_{t,l00};GeV/c"   , 20,0.,400. , ZW->first().daughter(0).pt() , theWeight);
    theHistograms->fill("Z_l1_pt"               , "p_{t,l01};GeV/c"   , 20,0.,400. , ZW->first().daughter(1).pt() , theWeight);
    theHistograms->fill("W_l_pt"                , "p_{t,l10};GeV/c"   , 20,0.,400. , ZW->second().daughter(0).pt(), theWeight);
    theHistograms->fill("W_MET_pt"              , "p_{t,MET};GeV/c"   , 20,0.,400. , ZW->second().daughter(1).pt(), theWeight);

    theHistograms->fill("ZW_massT_"+channel_reco, "m_{T,3l};GeV/c^{2}", 25,0.,500. , ZW->p4().Mt()                , theWeight);
    theHistograms->fill("Z_mass_"  +channel_reco, "m_{Z};GeV/c^{2}"   , 35,55.,125., ZW->first().mass()           , theWeight);
    theHistograms->fill("W_massT_" +channel_reco, "m_{T,W};GeV/c^{2}" , 35,55.,125., ZW->second().p4().Mt()       , theWeight);
    theHistograms->fill("ZW_pt_"   +channel_reco, "p_{t,ZW};GeV/c"    , 20,0.,400. , ZW->pt()                     , theWeight);
    theHistograms->fill("Z_l0_pt_" +channel_reco, "p_{t,l00};GeV/c"   , 20,0.,400. , ZW->first().daughter(0).pt() , theWeight);
    theHistograms->fill("Z_l1_pt_" +channel_reco, "p_{t,l01};GeV/c"   , 20,0.,400. , ZW->first().daughter(1).pt() , theWeight);
    theHistograms->fill("W_l_pt_"  +channel_reco, "p_{t,l10};GeV/c"   , 20,0.,400. , ZW->second().daughter(0).pt(), theWeight);
    theHistograms->fill("W_MET_pt_"+channel_reco, "p_{t,MET};GeV/c"   , 20,0.,400. , ZW->second().daughter(1).pt(), theWeight);
  }
  else if(LFR_lep){
    theHistograms->fill("ZL_mass_"              , "m_{3l};GeV/c^{2}", 25,0.,500. , (ZL->first.p4()+ZL->second.p4()).M(), theWeight);
    theHistograms->fill("Z_mass_"               , "m_{Z};GeV/c^{2}" , 35,55.,125., ZL->first.mass()                    , theWeight);
    theHistograms->fill("Z_l0_pt_"              , "p_{t,l00};GeV/c" , 20,0.,400. , ZL->first.daughter(0).pt()          , theWeight);
    theHistograms->fill("Z_l1_pt_"              , "p_{t,l00};GeV/c" , 20,0.,400. , ZL->first.daughter(1).pt()          , theWeight);
    theHistograms->fill("L_pt_"                 , "p_{t,l3};GeV/c"  , 20,0.,400. , ZL->second.pt()                     , theWeight);

    theHistograms->fill("ZL_mass_" +channel_reco, "m_{3l};GeV/c^{2}", 25,0.,500. , (ZL->first.p4()+ZL->second.p4()).M(), theWeight);
    theHistograms->fill("Z_mass_"  +channel_reco, "m_{Z};GeV/c^{2}" , 35,55.,125., ZL->first.mass()                    , theWeight);
    theHistograms->fill("Z_l0_pt_" +channel_reco, "p_{t,l00};GeV/c" , 20,0.,400. , ZL->first.daughter(0).pt()          , theWeight);
    theHistograms->fill("Z_l1_pt_" +channel_reco, "p_{t,l00};GeV/c" , 20,0.,400. , ZL->first.daughter(1).pt()          , theWeight);
    theHistograms->fill("L_pt_"    +channel_reco, "p_{t,l3};GeV/c"  , 20,0.,400. , ZL->second.pt()                     , theWeight);
  }
  	
  /*
    theHistograms->fill("nZtoChLep"    , "Number of Z->ll per event" , 7,0,7, genVBHelper_.ZtoChLep().size());
    theHistograms->fill("nZtoNeutrinos", "Number of Z->nn per event" , 7,0,7, genVBHelper_.ZtoNeutrinos().size());
    theHistograms->fill("nWtoLep"      , "Number of W->lnu per event", 7,0,7, genVBHelper_.WtoLep().size());
    theHistograms->fill("nZtoQ"        , "Number of Z->qq per event" , 7,0,7, genVBHelper_.ZtoQ().size());
    theHistograms->fill("nWtoQ"        , "Number of W->qq' per event", 7,0,7, genVBHelper_.WtoQ().size());
  */

  //int nVBs = genVBHelper_.ZtoChLep().size() + genVBHelper_.ZtoNeutrinos().size() + genVBHelper_.WtoLep().size() + genVBHelper_.ZtoQ().size() + genVBHelper_.WtoQ().size();
  //theHistograms->fill("nVBs", "Number of VB per event", 7,0,7, nVBs);
}


void VVGammaAnalyzer::end(TFile& fout){
  unsigned long evtN(0), analyzedN(0);
  float evtW(0.), analyzedW(0.);
  if(evtNInReg_.find(region_) != evtNInReg_.end()){
    evtN      = evtNInReg_.at(region_);
    evtW      = evtWInReg_.at(region_);
  }
  if(analyzedNInReg_.find(region_) != analyzedNInReg_.end()){
    analyzedN = analyzedNInReg_.at(region_);
    analyzedW = analyzedWInReg_.at(region_);
  }  // Otherwise there's not a single event in this region
  
  cout<<"\n\t----- "<<phys::regionType(region_)<<"-----\n";
  cout<<Form("Total events: %lu (weighted: %.3g)\n", evtN, evtW);
  cout<<Form("Passing cut:  %lu (weighted: %.3g)\n", analyzedN, analyzedW);
  cout<<Form("Fraction:     %.1f %% (weighted: %.1f %%)\n", 100.*analyzedN/evtN, 100.*analyzedW/evtW);
  
  // Label names
  endNameHistos();
}


void VVGammaAnalyzer::endNameHistos(){
  TH1* cuts = theHistograms->get("AAA_cuts");
  TH1* cuts_u = theHistograms->get("AAA_cuts_u");
  for(TH1* h : {cuts, cuts_u}){
    if(!h) continue;
    TAxis* x = h->GetXaxis();
    x->SetBinLabel(1, "All");
    x->SetBinLabel(2, "ZZ || ZW");
    x->SetBinLabel(3, "2l2j || 2l1J");
    x->SetBinLabel(4, "#gamma loose");
    x->SetBinLabel(5, "analyzed");
    x->SetBinLabel(6, "#gamma medium");
  }
  TH1* kinPhotons_ID = theHistograms->get("kinPhotons_ID");
  TH1* kinPhotons_ID_u = theHistograms->get("kinPhotons_ID_u");
  for(TH1* h : {kinPhotons_ID, kinPhotons_ID_u}){
    if(!h) continue;
    TAxis* x = h->GetXaxis();
    x->SetBinLabel(1, "Total");
    x->SetBinLabel(2, "Loose");
    x->SetBinLabel(3, "Medium");
    x->SetBinLabel(4, "Tight");
  }
  TH1* ph_eScale_count = theHistograms->get("ph_eScale_count");
  if(ph_eScale_count){
    TAxis* x = ph_eScale_count->GetXaxis();
    x->SetBinLabel(1, "Kin");
    x->SetBinLabel(2, "Kin_Edown");
    x->SetBinLabel(3, "Kin_Eup");
    x->SetBinLabel(4, "Good");
    x->SetBinLabel(5, "Good_Edown");
    x->SetBinLabel(6, "Good_Eup");
  }
}


void VVGammaAnalyzer::finish(){
  // cout<<Form("\nTotal events: %lu (weighted: %.2f)", evtN_, totEvtW_);
  // cout<<Form("\nPassing cut:  %lu (weighted: %.2f)", analyzedN_, analyzedW_);
  // cout<<Form("\nFraction:     %.1f %% (weighted: %.1f %%)", 100.*analyzedN_/evtN_, 100.*analyzedW_/totEvtW_)<<'\n';	
	
	
  float elapsedSec = (float)(clock()-startTime_)/CLOCKS_PER_SEC;
  int elapsedSecInt = (int)elapsedSec;
  cout<<"\nElapsed Time: "<<elapsedSec<<" s\t\t("<<elapsedSecInt/60<<"\' "<<elapsedSecInt%60<<"\")\n";
  for(char i=0; i<25; ++i) cout<<'-';
  cout<<" End of VZZAnalyzer ";
  for(char i=0; i<25; ++i) cout<<'-';
  cout<<"\n\n";
}


void VVGammaAnalyzer::initCherryPick(){
  for(auto R : regions_){
    FILE* cherryFile = fopen(Form("data/2016D_%s.txt", regionType(R).c_str()), "r");
    if(!cherryFile){
      cout << "Warning: no cherry pick file for region " << regionType(R) << '\n';
      continue;
    }
    unsigned long r, l, e;
    while( fscanf(cherryFile, "%lu:%lu:%lu", &r, &l, &e) == 3 )
      cherryEvents[R][r][l].insert(e);
    fclose(cherryFile);
  }
}


void VVGammaAnalyzer::genEventSetup(){
  genQuarks_->clear();
  genChLeptons_->clear();
  genNeutrinos_->clear();
  genPhotons_->clear();
	
  genZlepCandidates_->clear();
  genWlepCandidates_->clear();
  genZhadCandidates_->clear();
  genWhadCandidates_->clear();
	
  genZZ_ = DiBoson<Particle, Particle>();
  genWZ_ = DiBoson<Particle, Particle>();
	
  // Sort gen particles
  for(auto p : *genParticles){
    if(abs(p.id()) < 9)
      genQuarks_->push_back(p);
    else if(abs(p.id()) == 11 || abs(p.id()) == 13){
      genChLeptons_->push_back(p);
    }
    else if(abs(p.id()) == 12 || abs(p.id()) == 14)
      genNeutrinos_->push_back(p);
    else if( p.id() == 22 && p.pt() > 18.)
      genPhotons_->push_back(p);
  }
	
  // Gen W --> l nu
  if(genNeutrinos_->size() > 0 && genChLeptons_->size() > 0){
    for(auto l : *genChLeptons_){
      for(auto v : *genNeutrinos_){
	if( abs(l.id() + v.id()) == 1 ){
	  Boson<Particle> Wcand(l,v);
	  if(GenWBosonDefinition(Wcand))
	    genWlepCandidates_->push_back(Wcand);
	}
      }
    }
  }
	
  // Gen W --> q q'bar
  if(genQuarks_->size() >= 2){
    for(size_t i = 0  ; i < genQuarks_->size(); ++i){
      Particle& q1 = genQuarks_->at(i);
      if(q1.id() > 5) continue;
      for(size_t j = i+1; j < genQuarks_->size(); ++j){
	Particle& q2 = genQuarks_->at(j);
	if(q2.id() > 5) continue;
				
	if( (q1.id() * q2.id() < 0) && ( abs(q1.id()+q2.id()) % 2 ==1 ) ){
	  Boson<Particle> Wcand(q1,q2);
	  if(GenWBosonDefinition(Wcand))
	    genWhadCandidates_->push_back(Wcand);
	}
      }
    }
  }
	
  // Gen Z --> q qbar
  if(genQuarks_->size() >= 2){
    for(size_t i = 0  ; i < genQuarks_->size(); ++i){
      Particle& q1 = genQuarks_->at(i);
      if(q1.id() > 5) continue;
      for(size_t j = i+1; j < genQuarks_->size(); ++j){
	Particle& q2 = genQuarks_->at(j);
	if(q2	.id() > 5) continue;
				
	if( q1.id() + q2.id() == 0 ){
	  Boson<Particle> Zcand(q1,q2);
	  if(ZBosonDefinition(Zcand))
	    genZhadCandidates_->push_back(Zcand);
	}
      }
    }
  }
	
  // Gen Z --> l lbar
  if(genChLeptons_->size() >= 2){
    for(size_t i = 0 ; i < genChLeptons_->size(); ++i){
      Particle& l1 = genChLeptons_->at(i);
      for(size_t j = i+1; j < genChLeptons_->size(); ++j){
	Particle& l2 = genChLeptons_->at(j);
				
	if( l1.id() + l2.id() == 0 ){
	  Boson<Particle> Zcand(l1,l2);
	  if(ZBosonDefinition(Zcand))
	    genZhadCandidates_->push_back(Zcand);
	}
      }
    }
  }
	
  // genZZ --> 4l
  if(genChLeptons_->size() >= 4 && genZlepCandidates_->size() >= 2){
    std::sort(genZlepCandidates_->begin(), genZlepCandidates_->end(), MassComparator(phys::ZMASS));
    Boson<Particle>& Z0 = genZlepCandidates_->front();
		
    // Vector containing the rest of the Zll candidates
    vector<Boson<Particle>> Zll(genZlepCandidates_->begin()+1, genZlepCandidates_->end());
    std::sort(Zll.begin(), Zll.end(), ScalarSumPtComparator());
    Boson<Particle>* pZ1 = nullptr;
    for(size_t i = 0; i < Zll.size(); ++i){
      if(! haveCommonDaughter(Z0, Zll.at(i))){
	pZ1 = &(Zll.at(i));
	break;
      }
    }
    if(pZ1)
      genZZ_ = DiBoson<Particle, Particle>(Z0, *pZ1);
  }
	
  // genZW --> 3l nu
  if(genChLeptons_->size() >= 3 && genZlepCandidates_->size() >= 1 && genWlepCandidates_->size() >= 1){	
    std::sort(genZlepCandidates_->begin(), genZlepCandidates_->end(), MassComparator(phys::ZMASS));
    Boson<Particle>& Z0 = genZlepCandidates_->front();
		
    std::sort(genWlepCandidates_->begin(), genWlepCandidates_->end(), MassComparator(phys::WMASS));
    Boson<Particle>& W0 = genWlepCandidates_->front();
		
    genWZ_ = DiBoson<Particle, Particle>(Z0, W0);
  }
	
}


void VVGammaAnalyzer::genEventHistos(){
  theHistograms->fill("GEN_genZlepCandidates", "# genZlepCandidates_", 4,-0.5,3.5, genZlepCandidates_->size());
  theHistograms->fill("GEN_genWlepCandidates", "# genWlepCandidates_", 4,-0.5,3.5, genWlepCandidates_->size());
  theHistograms->fill("GEN_genZhadCandidates", "# genZhadCandidates_", 4,-0.5,3.5, genZhadCandidates_->size());
  theHistograms->fill("GEN_genWhadCandidates", "# genWhadCandidates_", 4,-0.5,3.5, genWhadCandidates_->size());
	
  for(auto v : *genZlepCandidates_)
    theHistograms->fill("GEN_genZlepCandidates_mass", "mass genZlepCandidates;[GeV/c^{2}]", 35.,50.,120., v.mass());
  for(auto v : *genWlepCandidates_)
    theHistograms->fill("GEN_genWlepCandidates_mass", "mass genWlepCandidates;[GeV/c^{2}]", 35.,50.,120., v.mass());
  for(auto v : *genZhadCandidates_)
    theHistograms->fill("GEN_genZhadCandidates_mass", "mass genZhadCandidates;[GeV/c^{2}]", 35.,50.,120., v.mass());
  for(auto v : *genWhadCandidates_)
    theHistograms->fill("GEN_genWhadCandidates_mass", "mass genWhadCandidates;[GeV/c^{2}]", 35.,50.,120., v.mass());
}


void VVGammaAnalyzer::baseHistos_cut(){
  theHistograms->fill("theWeight_cut", "Event weights in cut()", 120, -0.1, 1.1, theWeight);
  // Photons
  for(const Photon& ph: *kinPhotons_){
    theHistograms->fill("kinPhotons_pt"  , "p_{t} of kinPhotons" , 50, 0., 250.,      ph.pt()  , theWeight);
    theHistograms->fill("kinPhotons_aeta", "|#eta| of kinPhotons", 50, 0., 2.5 , fabs(ph.eta()), theWeight);

                              theHistograms->fill("kinPhotons_ID_u", "Cut Based ID of kinPhotons", 4,-0.5,3.5, 0, 1.);
    if(ph.cutBasedIDLoose())  theHistograms->fill("kinPhotons_ID_u", "Cut Based ID of kinPhotons", 4,-0.5,3.5, 1, 1.);
    if(ph.cutBasedIDMedium()) theHistograms->fill("kinPhotons_ID_u", "Cut Based ID of kinPhotons", 4,-0.5,3.5, 2, 1.);
    if(ph.cutBasedIDTight())  theHistograms->fill("kinPhotons_ID_u", "Cut Based ID of kinPhotons", 4,-0.5,3.5, 3, 1.);
    
                              theHistograms->fill("kinPhotons_ID", "Cut Based ID of kinPhotons", 4,-0.5,3.5, 0, theWeight);
    if(ph.cutBasedIDLoose())  theHistograms->fill("kinPhotons_ID", "Cut Based ID of kinPhotons", 4,-0.5,3.5, 1, theWeight);
    if(ph.cutBasedIDMedium()) theHistograms->fill("kinPhotons_ID", "Cut Based ID of kinPhotons", 4,-0.5,3.5, 2, theWeight);
    if(ph.cutBasedIDTight())  theHistograms->fill("kinPhotons_ID", "Cut Based ID of kinPhotons", 4,-0.5,3.5, 3, theWeight);
  }
  
  // Leptons
  auto lead_ele = std::max_element(electrons->begin(), electrons->end(), [](const Lepton& a, const Lepton& b){ return a.pt() < b.pt(); } );
  auto lead_muo = std::max_element(muons->begin()    , muons->end()    , [](const Lepton& a, const Lepton& b){ return a.pt() < b.pt(); } );
  if(lead_ele != electrons->end())
    theHistograms->fill("lead_ele_pt", "Leading electron p_{t}:p_{t} [GeV/c]", 50, 0., 250., lead_ele->pt(), theWeight);
  if(lead_muo != muons->end())
    theHistograms->fill("lead_muo_pt", "Leading muoon p_{t}:p_{t} [GeV/c]"   , 50, 0., 250., lead_muo->pt(), theWeight);
  

  // JetsAK8
  vector<std::pair<const Particle*, const Jet*>> JetsAK8genrec = matchGenRec(*genJetsAK8, *jetsAK8);
  for(auto pair: JetsAK8genrec){
    const Particle& gen = * pair.first;
    theHistograms->fill("AK8_gen_pt"  , "AK8 generated;p_{t} [GeV/c]"             , 40,100,500, gen.pt(), theWeight);
    theHistograms->fill("AK8_gen_pt_u", "AK8 generated (unweighted);p_{t} [GeV/c]", 40,100,500, gen.pt());
    
    if(! pair.second)
      continue;
    const Jet& rec = * pair.second;
    theHistograms->fill("AK8_rec_gen_pt"  , "AK8 reconstructed, gen p_{t};p_{t} [GeV/c]"             , 40,100,500, gen.pt(), theWeight);
    theHistograms->fill("AK8_rec_gen_pt_u", "AK8 reconstructed, gen p_{t} (unweighted);p_{t} [GeV/c]", 40,100,500, gen.pt());
    theHistograms->fill("AK8_resolution_dR"        ,"AK8: #DeltaR(reco,gen)"                   ,40,  0.,0.2,physmath::deltaR(gen, rec)    );
    theHistograms->fill("AK8_resolution_pt"        ,"AK8: pt - genpt;[GeV/c]"                  ,41,-41.,41.,rec.pt()           -gen.pt()  );
    theHistograms->fill("AK8_resolution_corrPruned","AK8: corrPrunedMass - genMass;[GeV/c^{2}]",41,-41.,41.,rec.corrPrunedMass()-gen.mass());
    theHistograms->fill("AK8_resolution_pruned"    ,"AK8: prunedMass - genMass;[GeV/c^{2}]"    ,41,-41.,41.,rec.prunedMass()   -gen.mass());
    theHistograms->fill("AK8_resolution_softDrop"  ,"AK8: softDropMass - genMass;[GeV/c^{2}]"  ,41,-41.,41.,rec.softDropMass() -gen.mass());
    theHistograms->fill("AK8_resolution_mass"      ,"AK8: mass - genMass;[GeV/c^{2}]"          ,41,-41.,41.,rec.mass()         -gen.mass());
  }
  
  // Jets AK4
  vector<std::pair<const Particle*, const Jet*>> JetsAK4genrec = matchGenRec(*genJets, *jets);
  for(auto pair: JetsAK4genrec){
    const Particle& gen = * pair.first;
    theHistograms->fill("AK4_gen_pt"  , "AK4 generated;p_{t} [GeV/c]"             , 40,100,500, gen.pt(), theWeight);
    theHistograms->fill("AK4_gen_pt_u", "AK4 generated (unweighted);p_{t} [GeV/c]", 40,100,500, gen.pt());
    
    if(! pair.second)
      continue;
    const Jet& rec = * pair.second;
    theHistograms->fill("AK4_rec_gen_pt"  , "AK4 reconstructed, gen p_{t};p_{t} [GeV/c]"             , 47,30,500, gen.pt(), theWeight);
    theHistograms->fill("AK4_rec_gen_pt_u", "AK4 reconstructed, gen p_{t} (unweighted);p_{t} [GeV/c]", 47,30,500, gen.pt());
    theHistograms->fill("AK4_resolution_dR"        ,"AK4: #DeltaR(reco,gen)"                   ,40,  0.,0.2,physmath::deltaR(gen, rec)    );
    theHistograms->fill("AK4_resolution_pt"        ,"AK4: pt - genpt;[GeV/c]"                  ,41,-20.5,20.5,rec.pt()           -gen.pt()  );
    theHistograms->fill("AK4_resolution_mass"      ,"AK4: mass - genMass;[GeV/c^{2}]"          ,41,-20.5,20.5,rec.mass()         -gen.mass());
  }
}


void VVGammaAnalyzer::PKU_comparison(){
  std::string channel_gen("??");

  if(theSampleInfo.isMC()){
    genEventSetup();
  
    vector<Particle> genEle, genMuo;
    for(const Particle& l : *genChLeptons_){
      if     (abs(l.id()) == 11) genEle.push_back(l);
      else if(abs(l.id()) == 13) genMuo.push_back(l);
      else cout<<">>>genLept: "<<l.id()<<'\n';
    }
    
    if     (genMuo.size() == 2) channel_gen = "mm";
    else if(genEle.size() == 2) channel_gen = "ee";
    else channel_gen = Form("%lue%lum", genEle.size(), genMuo.size());
  }
  
  std::string title_N( Form("Z #rightarrow %s: Number of %s", channel_gen.c_str(), "%s") );
  theHistograms->fill(Form("PKU_%s_kinPh_N" , channel_gen.c_str()), Form(title_N.c_str(), "photons passing kinematics")      ,5,-0.5,4.5, kinPhotons_->size() , 1);
  theHistograms->fill(Form("PKU_%s_goodPh_N", channel_gen.c_str()), Form(title_N.c_str(), "photons passing cutBasedIDLoose") ,5,-0.5,4.5, goodPhotons_->size(), 1);
  theHistograms->fill(Form("PKU_%s_POGele_N", channel_gen.c_str()), Form(title_N.c_str(), "electrons")                       ,5,-0.5,4.5, electrons->size()   , 1);
  theHistograms->fill(Form("PKU_%s_POGmuo_N", channel_gen.c_str()), Form(title_N.c_str(), "muon")                            ,5,-0.5,4.5, muons->size()       , 1);
  
  std::string title_pt( Form("Z #rightarrow %s: p_{t} of %s;p_{t} [GeV/c]", channel_gen.c_str(), "%s") );
  for(const Photon& ph : *kinPhotons_ )
    theHistograms->fill(Form("PKU_%s_kinPh_pt" , channel_gen.c_str()), Form(title_pt.c_str(), "photons passing kinematics")     , 20,0.,200., ph.pt(), 1);
  for(const Photon& ph : *goodPhotons_)
    theHistograms->fill(Form("PKU_%s_goodPh_pt", channel_gen.c_str()), Form(title_pt.c_str(), "photons passing cutBasedIDLoose"), 20,0.,200., ph.pt(), 1);
  
  for(const Lepton& el : *electrons )
    theHistograms->fill(Form("PKU_%s_POGele_pt", channel_gen.c_str()), Form(title_pt.c_str(), "electrons")                      , 20,0.,200., el.pt(), 1);
  for(const Lepton& mu : *muons)
    theHistograms->fill(Form("PKU_%s_POGmuo_pt", channel_gen.c_str()), Form(title_pt.c_str(), "muons")                          , 20,0.,200., mu.pt(), 1);
}


void VVGammaAnalyzer::LeptonFakeRate(){
  if(region_ != phys::CRLFR || !ZL || ZL->first.pt() * ZL->second.pt() <1.)
    return;  // We must be in the LFR control region
  if(met->pt() > 30)
    return;  // Exclude WZ phase space
  
  const Lepton& lep = ZL->second;

  std::string flavour;
  switch( abs(lep.id()) ){
  case 11:
    flavour = "electron";
    break;
  case 13:
    flavour = "muon";
    break;
  default:
    return;
  }
  std::string status = lep.isGood() ? "PASS" : "FAIL";
  std::string eta_range = fabs(lep.eta()) < 1.5 ? "EB" : "EE";
  
  
  theHistograms->fill(Form("LFR_%s_%s_%s", flavour.c_str(), eta_range.c_str(), status.c_str()), 
		      Form("Fake %s in %s: %s;p_{t} [GeV/c]", flavour.c_str(), eta_range.c_str(), status.c_str()),
		      pt_bins_LFR, lep.pt(), theWeight);
  
  theHistograms->fill(Form("LFR_%s_%s", flavour.c_str(), status.c_str()), 
		      Form("Fake %s: %s;p_{t} [GeV/c]", flavour.c_str(), status.c_str()),
		      pt_bins_LFR, aeta_bins, lep.pt(), fabs(lep.eta()), theWeight);
  return;
}


void VVGammaAnalyzer::PhotonFakeRate(){
  if( !(region_ == phys::CRLFR || region_ == phys::SR2P || region_ == phys::SR2P_1L) )
    return;
  if(kinPhotons_->size() < 1)
    return;
  
  if( (region_ == phys::SR2P || region_ == phys::SR2P_1L) ){
    if(jets->size() >= 2){
      //Try to reconstruct V hadronic
      vector<Boson<Jet>> VtojjCands;  VtojjCands.reserve( jets->size()*(jets->size()-1)/2 );
      for(auto it_i = jets->begin(); it_i <= jets->end() ; ++it_i){
	for(auto it_j = it_i+1; it_j <= jets->end() ; ++it_j){
	  Boson<Jet> cand(*it_i, *it_j);
	  if(physmath::minDM(cand.mass()) < 30.)
	    VtojjCands.push_back(std::move(cand));
	}
      }
    
      if(VtojjCands.size() >= 1)
	return;
    }
  
    if(  std::any_of(jetsAK8->begin(), jetsAK8->end(), [](const Jet& cand){ return physmath::minDM(cand.mass()) < 30.; } )  )
      return;
  }
  
  // We're not in the signal region
  for(const Photon& ph : *kinPhotons_){
    std::string status = ph.cutBasedIDLoose() ? "PASS" : "FAIL";
    std::string eta_range = fabs(ph.eta()) > 1.5 ? "EE" : "EB";
    double reg_pt = ph.pt() > pt_bins.back() ? pt_bins.back() : ph.pt();
    
    theHistograms->fill(Form("PhFR_%s_%s", eta_range.c_str(), status.c_str()), 
			Form("Photons that pass kinematic cuts in %s and %s loose ID;p_{t} [GeV/c]", eta_range.c_str(), status.c_str()),
			pt_bins, ph.pt(), theWeight);

    theHistograms->fill(Form("PhFR_%s", status.c_str()),
			Form("Photons that pass kinematic cuts and %s loose ID;p_{t} [GeV/c];|#eta|", status.c_str()),
			pt_bins, aeta_bins, ph.pt(), fabs(ph.eta()), theWeight);
  }
}


void VVGammaAnalyzer::effPhotons(){
  // Photon gen-reco matching
  vector<Photon> goodPhotonsC(*goodPhotons_);
	
  for(auto gPh : *genPhotons_){
    float ph_aeta = fabs(gPh.eta());
    //if(ph_aeta > 1.4442 && ph_aeta < 1.566) continue;
		
    theHistograms->fill("eff: den G pt", "p_{t} gen #gamma", 25,0.,250., gPh.pt(), theWeight);
    theHistograms->fill("eff: den G eta","#eta gen #gamma", eta_bins, gPh.eta(),theWeight);
		
    if(goodPhotonsC.size() == 0 ) continue;
		
    std::sort(goodPhotonsC.begin(), goodPhotonsC.end(), DeltaRComparator(gPh));
    Photon& rPh = goodPhotonsC.front();
    if(deltaR(rPh, gPh) < 0.1){
      theHistograms->fill("res: dR gen rec","#DeltaR(#gamma_{GEN}, #gamma_{REC})", 20,0.,0.1, deltaR(rPh,gPh), theWeight);
      theHistograms->fill("eff: num G pt","p_{t} gen #gamma", 25,0.,250., gPh.pt(), theWeight);
      theHistograms->fill("eff: num G eta","#eta gen #gamma",eta_bins, gPh.eta(),theWeight);
    }
  }
	
	
  unsigned int tightPh = 0;
  for(auto ph: *goodPhotons_){
    theHistograms->fill("A: goodPh pt","p_{t} loose #gamma", 25,0.,250., ph.pt(), theWeight);
    if(ph.cutBasedIDTight()){
      tightPh++;
      theHistograms->fill("A: goodPh tight pt","p_{t} tight #gamma", 25,0.,250., ph.pt(), theWeight);
      theHistograms->fill("A: goodPh tight eta","#eta tight #gamma",25,-2.5,2.5, ph.eta(),theWeight);
    }
  }
	
  theHistograms->fill("A: tightG num", "# tight #gamma", 5,-0.5,4.5, tightPh, theWeight);
	
  if(tightPh > 0){
    theHistograms->fill("AAA_cuts"  , "Cuts weighted"  , CUT_LAYOUT, 6, theWeight);
    theHistograms->fill("AAA_cuts_u", "Cuts unweighted", CUT_LAYOUT, 6);
  }
}


bool VVGammaAnalyzer::cherrypickEvt() const {
  auto R = cherryEvents.find(region_);
  if(R == cherryEvents.end()) return false;

  auto r = R->second.find(run);
  if(r == R->second.end())    return false;

  auto l = r->second.find(lumiBlock);
  if(l == r->second.end())    return false;

  auto e = l->second.find(event);
  if(e == l->second.end())    return false;

  return true;
}


const vector<double> VVGammaAnalyzer::pt_bins(
  {20., 35., 50., 100., 200., 500.}
					      );

const vector<double> VVGammaAnalyzer::pt_bins_LFR(
  {5., 7., 10., 20., 30., 40., 50., 80}
					      );

const vector<double> VVGammaAnalyzer::eta_bins(
  {-2.4, -2., -1.566, -1.4442, -0.8, 0., 0.8, 1.4442, 1.566, 2., 2.4}
  //{-2.4, -2.1, -1.8, -1.566, -1.4442, -1.2, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.2, 1.4442, 1.566, 1.8, 2.1, 2.4}
					       );

const vector<double> VVGammaAnalyzer::aeta_bins(
  {0., 0.8, 1.4442, 1.566, 2., 2.4}
  //{0., 0.3, 0.6, 0.9, 1.2, 1.4442, 1.566, 1.8, 2.1, 2.4}
						);
