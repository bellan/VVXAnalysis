#include "VVXAnalysis/TreeAnalysis/interface/VVGammaAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include "VVXAnalysis/Commons/interface/GenVBHelper.h"

#include "TTree.h"
#include "Math/Vector2D.h"
#include "TRandom.h"

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

// #define DEBUG

#define BINS_CUTFLOW 8,-0.5,7.5
//{"All", "ZZ || ZW", "2l2j || 2l1J", "#gamma kin", "#gamma good", "Analyzed", "#gamma medium"}
#define BINS_PHCUTFLOW 6,-0.5,5.5
#define BINS_KINPHID 4,-0.5,3.5

std::pair<TLorentzVector, TLorentzVector> solveNuPz(const Boson<Lepton>& W, int& error);

void VVGammaAnalyzer::begin(){
  cout<<'\n';
  for(char i=0; i<25; ++i) cout<<'-';
  cout<<" Start of VVGammaAnalyzer ";
  for(char i=0; i<25; ++i) cout<<'-';
  cout<<'\n';
  
  const char* file_LOCAL = "data/LtoT_FR_data-ZGToLLG.root";
  TFile fileFR( file_LOCAL, "READ" );
  hPhotonFR_.reset( std::move((TH2F*) fileFR.Get("PhFR")) );
  hPhotonFR_->SetDirectory(nullptr);  // prevent ROOT from deleting it
  cout << "INFO: retrieved Photon FR histogram from \""<<fileFR.GetName()<<"\"\n";
  fileFR.Close();
  
  // initCherryPick();
  for(const char* sys : {"central", "EScale_Up", "EScale_Down", "ESigma_Up", "ESigma_Down"}){
    kinPhotons_ [sys]  = std::make_unique<vector<Photon>>();
    loosePhotons_[sys] = std::make_unique<vector<Photon>>();
    goodPhotons_[sys]  = std::make_unique<vector<Photon>>();
  }
  
  #ifndef DEBUG
  size_t digits = std::to_string(tree()->GetEntries()).length();
  std::string spaces( digits, ' ' );
  cout<<"Analyzed:\t"<<spaces<<'/'<<tree()->GetEntries()<<std::flush;
  #else
  cout<<'\n';
  #endif

  return;
}


void VVGammaAnalyzer::initEvent(){
  // Cleanup
  leptons_    ->clear();
  for(auto const& it: kinPhotons_ )  it.second->clear();
  for(auto const& it: loosePhotons_) it.second->clear();
  for(auto const& it: goodPhotons_)  it.second->clear();
  candVTojj_ = Boson<Jet>();
  fakeVTojj_ = Boson<Jet>();
  candVToJ_ = Jet();
  fakeVToJ_ = Jet();

  // Contruct vector with leptons from dibosons
  if     (is4Lregion(region_))  //(ZZ && ZZ->pt() > 1.){
    leptons_->insert(leptons_->end(), {
    	  ZZ->first().daughter(0), 
    	  ZZ->first().daughter(1), 
    	  ZZ->second().daughter(0), 
    	  ZZ->second().daughter(1)
        });
  else if(is3Lregion(region_))  //(ZW && ZW->pt() > 1.){
    leptons_->insert(leptons_->end(), {
  	  ZW->first().daughter(0), 
  	  ZW->first().daughter(1), 
  	  ZW->second().daughter(0)
	});
  else if(is2Lregion(region_)){
    leptons_->insert(leptons_->cend(), electrons->begin(), electrons->end());
    leptons_->insert(leptons_->cend(), muons    ->begin(), muons    ->end());
  }
  else if(region_ == CRLFR && ZL && ZL->first.pt() > 1.)
    leptons_->insert(leptons_->end(), {
  	  ZL->first.daughter(0),
  	  ZL->first.daughter(1),
  	  ZL->second
	});

  
  // Photon selection
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
    for(const Lepton lep : *leptons_){
      if(deltaR(ph,lep) < 0.3){
	match = true;
	break;
      }
    }
    if(match) continue;
    
    Photon::IDwp wp = Photon::IDwp::Loose;
    TLorentzVector p4_EScale_Up = ph.p4() * (ph.energyScaleUp()  /ph.e());
    TLorentzVector p4_EScale_Dn = ph.p4() * (ph.energyScaleDown()/ph.e());
    TLorentzVector p4_ESigma_Up = ph.p4() * (ph.energySigmaUp()  /ph.e());
    TLorentzVector p4_ESigma_Dn = ph.p4() * (ph.energySigmaDown()/ph.e());
    
    if(p4_EScale_Up.Pt() > 20){
      Photon copy(ph);
      copy.setP4(p4_EScale_Up);
      kinPhotons_["EScale_Up"]->push_back(copy);
      if(copy.passHoverE(wp) && copy.passPhotonIsolation(wp) && copy.passNeutralIsolation(wp))
	loosePhotons_["EScale_Up"]->push_back(copy);
      if(copy.cutBasedIDLoose())
	goodPhotons_["EScale_Up"]->push_back(std::move(copy));
    }
    if(p4_EScale_Dn.Pt() > 20){
      Photon copy(ph);
      copy.setP4(p4_EScale_Dn);
      kinPhotons_["EScale_Down"]->push_back(copy);
      if(copy.passHoverE(wp) && copy.passPhotonIsolation(wp) && copy.passNeutralIsolation(wp))
    	loosePhotons_["EScale_Down"]->push_back(copy);
      if(copy.cutBasedIDLoose())
    	goodPhotons_["EScale_Down"]->push_back(std::move(copy));
    }
    if(p4_ESigma_Up.Pt() > 20){
      Photon copy(ph);
      copy.setP4(p4_ESigma_Up);
      kinPhotons_["ESigma_Up"]->push_back(copy);
      if(copy.passHoverE(wp) && copy.passPhotonIsolation(wp) && copy.passNeutralIsolation(wp))
        loosePhotons_["ESigma_Up"]->push_back(copy);
      if(copy.cutBasedIDLoose())
    	goodPhotons_["ESigma_Up"]->push_back(std::move(copy));
    }
    if(p4_ESigma_Dn.Pt() > 20){
      Photon copy(ph);
      copy.setP4(p4_ESigma_Dn);
      kinPhotons_["ESigma_Down"]->push_back(copy);
      if(copy.passHoverE(wp) && copy.passPhotonIsolation(wp) && copy.passNeutralIsolation(wp))
    	loosePhotons_["ESigma_Down"]->push_back(copy);
      if(copy.cutBasedIDLoose())
    	goodPhotons_["ESigma_Down"]->push_back(std::move(copy));
    }
    if(ph.pt() > 20){
      kinPhotons_["central"]->push_back(ph);
      if(ph.passHoverE(wp) && ph.passPhotonIsolation(wp) && ph.passNeutralIsolation(wp))
	loosePhotons_["central"]->push_back(ph);
      if(ph.cutBasedIDLoose())
	goodPhotons_["central"]->push_back(ph);
    }
    
  }
  
}


Int_t VVGammaAnalyzer::cut() {
  evtN_++; evtNInReg_[region_]++; evtWInReg_[region_] += theWeight;
  #ifndef DEBUG
  cout<<"\r\t\t"<<evtN_;
  #endif
  //cout << regionType(region_) << ':' << run << ':' << lumiBlock << ':' << event << '\n';
  // if(cherrypickEvt()){
  //   theHistograms->fill("cherry_ZZ_mass"    , "Events in {UL#backslashLegacy}: m_{ZZ};m_{ZZ} [GeV/c^{2}]", 25,0.,500., ZZ->mass()                      , theWeight);
  //   theHistograms->fill("cherry_ZZ_goodLept", "Events in {UL#backslashLegacy}: # good leptons"           , 5,-0.5,4.5, ZZ->numberOfGoodGrandDaughters(), theWeight);
  //   theHistograms->fill("cherry_ZZ_badLept" , "Events in {UL#backslashLegacy}: # bad leptons"            , 5,-0.5,4.5, ZZ->numberOfBadGrandDaughters() , theWeight);
  // }
  
  theHistograms->fill("AAA_cuts"  , "cuts weighted"  , BINS_CUTFLOW, 0, theWeight);
  theHistograms->fill("AAA_cuts_u", "cuts unweighted", BINS_CUTFLOW, 0, 1);

  bool haveZVlep = false;
  bool have2l2j = false;
  bool haveGoodPhoton = false;
  
  initEvent();
  
  theHistograms->fill("POG_leptons", "n leptons", 7,-0.5,6.5, electrons->size()+muons->size(), theWeight);  
  
  baseHistos_cut();
  photonHistos();
  jetHistos();
  // PKU_comparison();
  photonIsolation(*    photons            , "all" );
  photonIsolation(* kinPhotons_["central"], "kin" );
  photonIsolation(*goodPhotons_["central"], "good");
  photonEfficiency(*    photons            ,"photons");
  photonEfficiency(* kinPhotons_["central"],  "kinPh");
  photonEfficiency(*goodPhotons_["central"], "goodPh");
  
  
  // ----- BASELINE SELECTION -----
  // ----- Cut1: require at least a ZZ or a WZ candidate
  haveZVlep = (ZZ && ZZ->pt() > 1.) || (ZW && ZW->pt() > 1.);	
  if(haveZVlep){
    theHistograms->fill("AAA_cuts"  , "cuts weighted"  , BINS_CUTFLOW, 1, theWeight);
    theHistograms->fill("AAA_cuts_u", "cuts unweighted", BINS_CUTFLOW, 1, 1);
  }

  have2l2j = (muons->size()+electrons->size()==2) && (jets->size()==2 || jetsAK8->size()>=1);
  if(have2l2j){
    theHistograms->fill("AAA_cuts"  , "cuts weighted"  , BINS_CUTFLOW, 2, theWeight);
    theHistograms->fill("AAA_cuts_u", "cuts unweighted", BINS_CUTFLOW, 2, 1);
  }
  
  if(kinPhotons_["central"]->size() >= 1){
    theHistograms->fill("AAA_cuts"  , "cuts weighted"  , BINS_CUTFLOW, 3, theWeight);
    theHistograms->fill("AAA_cuts_u", "cuts unweighted", BINS_CUTFLOW, 3, 1);
  }
  if(loosePhotons_["central"]->size() >= 1){
    theHistograms->fill("AAA_cuts"  , "cuts weighted"  , BINS_CUTFLOW, 4, theWeight);
    theHistograms->fill("AAA_cuts_u", "cuts unweighted", BINS_CUTFLOW, 4, 1);
  }
  // ----- Cut2: Require at least 1 loose photon with pt > 20 GeV
  // haveGoodPhoton = std::any_of(goodPhotons_->begin(), goodPhotons->end(), 
  // 			       [](const std::pair<const char*, std::unique_ptr<vector<Photon>>>& p){ return p.second->size() >= 1; } 
  // 			       );
  haveGoodPhoton = goodPhotons_["central"]->size() >= 1;
  if(haveGoodPhoton){
    theHistograms->fill("AAA_cuts"  , "cuts weighted"  , BINS_CUTFLOW, 5, theWeight);
    theHistograms->fill("AAA_cuts_u", "cuts unweighted", BINS_CUTFLOW, 5, 1);
  }
  
  
  return 1 ; //TEMP
  
  if(!haveZVlep) 
    return -1;
  else if(!haveGoodPhoton)
    return -1;
  else
    return 1;
}


void VVGammaAnalyzer::analyze(){
  analyzedNInReg_[region_]++; analyzedWInReg_[region_] += theWeight;
  theHistograms->fill("AAA_cuts"  , "cuts weighted"  , BINS_CUTFLOW, 6, theWeight);
  theHistograms->fill("AAA_cuts_u", "cuts unweighted", BINS_CUTFLOW, 6, 1);
  
  if( std::any_of(goodPhotons_["central"]->cbegin(),
		  goodPhotons_["central"]->cend(), 
		  [](const Photon& ph){return ph.cutBasedIDMedium();}) )
  {
    theHistograms->fill("AAA_cuts"  , "cuts weighted"  , BINS_CUTFLOW, 7, theWeight);
    theHistograms->fill("AAA_cuts_u", "cuts unweighted", BINS_CUTFLOW, 7, 1);
  }
  
  genEventSetup();
  
  bool four_lep  = is4Lregion(region_);  // && ZZ && ZZ->pt() > 1.;
  bool three_lep = is3Lregion(region_);  // && ZW && ZW->pt() > 1.;
  bool two_lep   = is2Lregion(region_);
  bool LFR_lep   = region_ == CRLFR;     // && ZL && ZL->first.pt() > 1.;
  bool b_WZpaperSel = false;
  
  if(two_lep)
    hadronicObjectsReconstruction();  

  if(three_lep)
    debug3Lregion();
  
  // Decide channel name depending on leptons
  makeChannelReco();
  
  // leptonFakeRate();
  photonFakeRate();
  systematicsStudy();
  
  // Basic histograms on leptonic side
  if     (four_lep){
    // const Boson<Lepton>& Z0 = ZZ->first();
    // const Boson<Lepton>& Z1 = ZZ->second();
    // cout << "region: " << regionType(region_) << '\n';
    // cout << ">>> [eleEffSFUnc]      ZZ: " << ZZ->eleEffSFUnc() << "   Z0: " << Z0.eleEffSFUnc() << "   Z1: " << Z1.eleEffSFUnc() << '\n';
    // cout << ">>> [muEffSFUnc]       ZZ: " << ZZ->muEffSFUnc()  << "   Z0: " << Z0.muEffSFUnc()  << "   Z1: " << Z1.muEffSFUnc()  << '\n';
    // cout << "\tZ0_l0: "<<Z0.daughter(0).efficiencySFUnc()<<"   Z0_l1: "<<Z0.daughter(1).efficiencySFUnc()<<"   Z1_l0: "<<Z1.daughter(0).efficiencySFUnc()<<"   Z1_l1: "<<Z1.daughter(1).efficiencySFUnc() << '\n';
    // cout << "--------------------------------------------------------------------------------\n";
    // cout << ">>> [fakeRateSF]       ZZ: " << ZZ->fakeRateSF()       << "   Z0: "<<Z0.fakeRateSF()        << "   Z1: " << Z1.fakeRateSF() << '\n';
    // cout << "\tZ0_l0: "<<Z0.daughter(0).fakeRateSF()   <<"   Z0_l1: "<<Z0.daughter(1).fakeRateSF()   <<"   Z1_l0: "<<Z1.daughter(0).fakeRateSF()   <<"   Z1_l1: "<<Z1.daughter(1).fakeRateSF()    << '\n';
    // cout << "\tZ0_l0: "<<Z0.daughter(0).fakeRateSFUnc()<<"   Z0_l1: "<<Z0.daughter(1).fakeRateSFUnc()<<"   Z1_l0: "<<Z1.daughter(0).fakeRateSFUnc()<<"   Z1_l1: "<<Z1.daughter(1).fakeRateSFUnc() << '\n';
    // cout << ">>> [eleFakeRateSFUnc] ZZ: " << ZZ->eleFakeRateSFUnc() << "   Z1: "<<Z1.eleFakeRateSFUnc() << '\n';
    // cout << ">>> [muFakeRateSFUnc]  ZZ: " << ZZ->muoFakeRateSFUnc() << "   Z1: "<<Z1.muFakeRateSFUnc()  << '\n';
    // cout << "################################################################################\n";
    
    theHistograms->fill("ZZ_mass"               , "m_{4l};GeV/c^{2}", mVV_bins   , ZZ->mass()                   , theWeight);
    theHistograms->fill("Z0_mass"               , "m_{Z0};GeV/c^{2}", 35,55.,125., ZZ->first().mass()           , theWeight);
    theHistograms->fill("Z1_mass"               , "m_{Z1};GeV/c^{2}", 35,55.,125., ZZ->second().mass()          , theWeight);
    theHistograms->fill("ZZ_pt"                 , "p_{t,ZZ};GeV/c"  , 20,0.,400. , ZZ->pt()                     , theWeight);
    theHistograms->fill("Z0_l0_pt"              , "p_{t,l00};GeV/c" , 20,0.,400. , ZZ->first().daughter(0).pt() , theWeight);
    theHistograms->fill("Z0_l1_pt"              , "p_{t,l01};GeV/c" , 20,0.,400. , ZZ->first().daughter(1).pt() , theWeight);
    theHistograms->fill("Z1_l0_pt"              , "p_{t,l10};GeV/c" , 20,0.,400. , ZZ->second().daughter(0).pt(), theWeight);
    theHistograms->fill("Z1_l1_pt"              , "p_{t,l11};GeV/c" , 20,0.,400. , ZZ->second().daughter(1).pt(), theWeight);

    theHistograms->fill("ZZ_mass_" +channelReco_, "m_{4l};GeV/c^{2}", mVV_bins   , ZZ->mass()                   , theWeight);
    theHistograms->fill("Z0_mass"  +channelReco_, "m_{Z0};GeV/c^{2}", 35,55.,125., ZZ->first().mass()           , theWeight);
    theHistograms->fill("Z1_mass_" +channelReco_, "m_{Z1};GeV/c^{2}", 35,55.,125., ZZ->second().mass()          , theWeight);
    theHistograms->fill("ZZ_pt_"   +channelReco_, "p_{t,ZZ};GeV/c"  , 20,0.,400. , ZZ->pt()                     , theWeight);
    theHistograms->fill("Z0_l0_pt_"+channelReco_, "p_{t,l00};GeV/c" , 20,0.,400. , ZZ->first().daughter(0).pt() , theWeight);
    theHistograms->fill("Z0_l1_pt_"+channelReco_, "p_{t,l01};GeV/c" , 20,0.,400. , ZZ->first().daughter(1).pt() , theWeight);
    theHistograms->fill("Z1_l0_pt_"+channelReco_, "p_{t,l10};GeV/c" , 20,0.,400. , ZZ->second().daughter(0).pt(), theWeight);
    theHistograms->fill("Z1_l1_pt_"+channelReco_, "p_{t,l11};GeV/c" , 20,0.,400. , ZZ->second().daughter(1).pt(), theWeight);
  }

  else if(three_lep){
    double m_lll = (ZW->first().p4() + ZW->second().daughter(0).p4()).M();
    theHistograms->fill("ZW_massT"              , "m_{T,3l};GeV/c^{2}", mVV_bins   , ZW->p4().Mt()                , theWeight);
    theHistograms->fill("Z_mass"                , "m_{Z};GeV/c^{2}"   , 35,55.,125., ZW->first().mass()           , theWeight);
    theHistograms->fill("W_massT"               , "m_{T,W};GeV/c^{2}" , 35,55.,125., ZW->second().p4().Mt()       , theWeight);
    theHistograms->fill("ZW_pt"                 , "p_{t,ZW};GeV/c"    , 20,0.,400. , ZW->pt()                     , theWeight);
    theHistograms->fill("Z_l0_pt"               , "p_{t,l00};GeV/c"   , 20,0.,400. , ZW->first().daughter(0).pt() , theWeight);
    theHistograms->fill("Z_l1_pt"               , "p_{t,l01};GeV/c"   , 20,0.,400. , ZW->first().daughter(1).pt() , theWeight);
    theHistograms->fill("W_l_pt"                , "p_{t,l10};GeV/c"   , 20,0.,400. , ZW->second().daughter(0).pt(), theWeight);
    theHistograms->fill("W_MET_pt"              , "p_{t,MET};GeV/c"   , 20,0.,400. , ZW->second().daughter(1).pt(), theWeight);
    theHistograms->fill("lll_mass"              , "m_{lll};GeV/c^{2}" , 20,0.,400. , m_lll                        , theWeight);

    theHistograms->fill("ZW_massT_"+channelReco_, "m_{T,3l};GeV/c^{2}", mVV_bins   , ZW->p4().Mt()                , theWeight);
    theHistograms->fill("Z_mass_"  +channelReco_, "m_{Z};GeV/c^{2}"   , 35,55.,125., ZW->first().mass()           , theWeight);
    theHistograms->fill("W_massT_" +channelReco_, "m_{T,W};GeV/c^{2}" , 35,55.,125., ZW->second().p4().Mt()       , theWeight);
    theHistograms->fill("ZW_pt_"   +channelReco_, "p_{t,ZW};GeV/c"    , 20,0.,400. , ZW->pt()                     , theWeight);
    theHistograms->fill("Z_l0_pt_" +channelReco_, "p_{t,l00};GeV/c"   , 20,0.,400. , ZW->first().daughter(0).pt() , theWeight);
    theHistograms->fill("Z_l1_pt_" +channelReco_, "p_{t,l01};GeV/c"   , 20,0.,400. , ZW->first().daughter(1).pt() , theWeight);
    theHistograms->fill("W_l_pt_"  +channelReco_, "p_{t,l10};GeV/c"   , 20,0.,400. , ZW->second().daughter(0).pt(), theWeight);
    theHistograms->fill("W_MET_pt_"+channelReco_, "p_{t,MET};GeV/c"   , 20,0.,400. , ZW->second().daughter(1).pt(), theWeight);
    theHistograms->fill("lll_mass_"+channelReco_, "m_{lll};GeV/c^{2}" , 20,0.,400. , m_lll                        , theWeight);

    // Neutrino pz reconstruction
    // int error;
    // const Boson<Lepton>& Wboson = ZW->second();
    // std::pair<TLorentzVector, TLorentzVector> solutions = solveNuPz(Wboson, error);
    // if(!error){
    //   const TLorentzVector& pl = Wboson.daughter(0).p4();
    //   const TLorentzVector& pv = (gRandom->Integer(2) ? solutions.first : solutions.second);
    //   theHistograms->fill("ZW_mass_reconstructed", "m_{3l}, p_{z}^{#nu} from m_{W} contraint; GeV/c^2", mVV_bins, (pl + pv + ZW->first().p4()).M(), theWeight);
    // }
    
    // replicate the selection from CMS-SMP-20-014
    theHistograms->fill("WZ_cutflow", "Cutflow for WZ", 7,-0.5,6.5, 0, theWeight);
    const Lepton& lZ1 = ZW->first ().daughter(0);
    const Lepton& lZ2 = ZW->first ().daughter(1);
    const Lepton& lW  = ZW->second().daughter(0);
    bool b_lepton_pt = (lZ1.pt() > 25 &&
			lZ2.pt() > 10 &&
			lW .pt() > 25);
    bool b_Z_mass = fabs(ZW->first().mass() - phys::ZMASS) < 15;
    bool b_MET = met->pt() > 30;
    bool b_trilep_mass = (lZ1.p4() + lZ2.p4() + lW.p4()).M() > 100;
    bool b_bveto = ! std::any_of(jets->begin(), jets->end(), [](const Jet& j){ return j.csvtagger() > 0.7665; }); // TODO: use DeepCSV
    bool b_bvetoAK8 = ! std::any_of(jetsAK8->begin(), jetsAK8->end(), [](const Jet& j){ return j.csvtagger() > 0.7665; }); // TODO: use DeepCSV
    
    if(b_lepton_pt  ) theHistograms->fill("WZ_cutflow", "Cutflow for WZ", 7,-0.5,6.5, 1, theWeight);
    if(b_Z_mass     ) theHistograms->fill("WZ_cutflow", "Cutflow for WZ", 7,-0.5,6.5, 2, theWeight);
    if(b_MET        ) theHistograms->fill("WZ_cutflow", "Cutflow for WZ", 7,-0.5,6.5, 3, theWeight);
    if(b_trilep_mass) theHistograms->fill("WZ_cutflow", "Cutflow for WZ", 7,-0.5,6.5, 4, theWeight);
    if(b_bveto      ) theHistograms->fill("WZ_cutflow", "Cutflow for WZ", 7,-0.5,6.5, 5, theWeight);
    if(b_lepton_pt && b_Z_mass && b_MET && b_trilep_mass && b_bveto)
      theHistograms->fill("WZ_cutflow", "Cutflow for WZ", 7,-0.5,6.5, 6, theWeight);

    // Is the third lepton the problem?
    const char* leptFlav = ( abs(lW.id()) == 11 ? "e" : "m"  );
    theHistograms->fill(Form("l3_%s_pt" , leptFlav), "3rd lepton p_{t};p_{t} [GeV/c]"                , 50,0.,200., lW.pt()                 , theWeight);
    theHistograms->fill(Form("l3_%s_Iso", leptFlav), "3rd lepton combRelIsoFSRCorr;combRelIsoFSRCorr", 50,0.,0.1 , lW.pfCombRelIsoFSRCorr(), theWeight);
    if(b_MET){
      theHistograms->fill(Form("l3_%s_pt_MET" , leptFlav),"3rd lepton p_{t} - MET > 30;p_{t} [GeV/c]"                , 50,0.,200.,lW.pt()                 ,theWeight);
      theHistograms->fill(Form("l3_%s_Iso_MET", leptFlav),"3rd lepton combRelIsoFSRCorr - MET > 30;combRelIsoFSRCorr", 50,0.,0.1 ,lW.pfCombRelIsoFSRCorr(),theWeight);
    }
    
    b_WZpaperSel = b_Z_mass && b_MET && b_trilep_mass;
    if(b_WZpaperSel){
      theHistograms->fill("paperSel_ZW_massT"              , "m_{T,3l};GeV/c^{2}", mVV_bins   , ZW->p4().Mt()                , theWeight);
      theHistograms->fill("paperSel_Z_mass"                , "m_{Z};GeV/c^{2}"   , 35,55.,125., ZW->first().mass()           , theWeight);
      theHistograms->fill("paperSel_W_massT"               , "m_{T,W};GeV/c^{2}" , 35,55.,125., ZW->second().p4().Mt()       , theWeight);
      theHistograms->fill("paperSel_ZW_pt"                 , "p_{t,ZW};GeV/c"    , 20,0.,400. , ZW->pt()                     , theWeight);
      theHistograms->fill("paperSel_Z_l0_pt"               , "p_{t,l00};GeV/c"   , 20,0.,400. , ZW->first().daughter(0).pt() , theWeight);
      theHistograms->fill("paperSel_Z_l1_pt"               , "p_{t,l01};GeV/c"   , 20,0.,400. , ZW->first().daughter(1).pt() , theWeight);
      theHistograms->fill("paperSel_W_l_pt"                , "p_{t,l10};GeV/c"   , 20,0.,400. , ZW->second().daughter(0).pt(), theWeight);
      theHistograms->fill("paperSel_W_MET_pt"              , "p_{t,MET};GeV/c"   , 20,0.,400. , ZW->second().daughter(1).pt(), theWeight);
      theHistograms->fill("paperSel_lll_mass"              , "m_{lll};GeV/c^{2}" , 20,0.,400. , m_lll                        , theWeight);

      theHistograms->fill("paperSel_ZW_massT_"+channelReco_, "m_{T,3l};GeV/c^{2}", mVV_bins   , ZW->p4().Mt()                , theWeight);
    }
  }
  
  else if(two_lep){
    // TEMP: reconstruct here the Z candidates from POG leptons
    vector<Boson<Lepton>> Zcandidates;
    for(auto leps : {electrons, muons}){
      if(leps->size() >= 2){
	for(auto i = leps->begin() ; i != leps->end() ; ++ i)
	  for(auto j = i+1 ; j != leps->end() ; ++ j)
	    Zcandidates.push_back(Boson<Lepton>(*i, *j));
      }
    }
    
    if(Zcandidates.size() > 0){
      std::sort(Zcandidates.begin(), Zcandidates.end(), MassComparator(phys::ZMASS));
      Boson<Lepton>& Zcand = Zcandidates.front();
      switch( abs(Zcand.daughter(0).id()) ){
      case 11:
	channelReco_ = "2e"; break;
      case 13:
	channelReco_ = "2m"; break;
      }
      theHistograms->fill("Z_pt", "Z p_{T}", 30, 0.,300., Zcand.pt(), theWeight);
      
      // jet histos
      theHistograms->fill("AK4_N" , "# jets AK4", 6, -0.5, 5.5, jets->size(), theWeight);
      for(const Jet& jet : *jets)
	theHistograms->fill("AK4_pt", "p_{T} jets AK4;GeV/c", 50, 0., 500., jet.pt(), theWeight);
      theHistograms->fill("AK8_N" , "# jets AK8", 6, -0.5, 5.5, jetsAK8->size(), theWeight);
      for(const Jet& jet : *jetsAK8)
	theHistograms->fill("AK8_pt", "p_{T} jets AK8;GeV/c", 50, 0., 500., jet.pt(), theWeight);
      theHistograms->fill("Z_mass_"+channelReco_, "m_{ll};GeV/c^{2}", 30, 60., 120., Zcand.mass(), theWeight);
      
      // Vhad histograms
      if(candVTojj_.isValid()){
	theHistograms->fill("VTojj_mass", "V#rightarrowjj mass;GeV/c^{2}", 30, 50, 125, candVTojj_.mass(), theWeight);
	theHistograms->fill("VTojj_pt"  , "V#rightarrowjj p_{t};GeV/c"   , 25,  0, 500, candVTojj_.pt()  , theWeight);
	theHistograms->fill("VTojj_eta" , "V#rightarrowjj #eta;#eta"     , 25, -5,   5, candVTojj_.eta() , theWeight);
      }
      else if(candVToJ_.isValid()){
	theHistograms->fill("VToJ_mass" , "V#rightarrowJ mass;GeV/c^{2}" , 30, 50, 125, candVToJ_.mass(), theWeight);
	theHistograms->fill("VToJ_pt"   , "V#rightarrowJ p_{t};GeV/c"    , 25,  0, 500, candVToJ_.pt()  , theWeight);
	theHistograms->fill("VToJ_eta"  , "V#rightarrowJ #eta;#eta"      , 25, -5,   5, candVToJ_.eta() , theWeight);
	theHistograms->fill("VToJ_PNet_TvsQCD"     , "V#rightarrowJ ParticleNet: T vs QCD;T vs QCD", 40,0,1, candVToJ_.particleNet().TvsQCD, theWeight);
	theHistograms->fill("VToJ_PNet_WvsQCD"     , "V#rightarrowJ ParticleNet: W vs QCD;W vs QCD", 40,0,1, candVToJ_.particleNet().WvsQCD, theWeight);
	theHistograms->fill("VToJ_PNet_ZvsQCD"     , "V#rightarrowJ ParticleNet: Z vs QCD;Z vs QCD", 40,0,1, candVToJ_.particleNet().ZvsQCD, theWeight);
	theHistograms->fill("VToJ_deepAK8_TvsQCD"  , "V#rightarrowJ deepAK8: T vs QCD;T vs QCD"    , 40,0,1, candVToJ_.deepAK8()    .TvsQCD, theWeight);
	theHistograms->fill("VToJ_deepAK8_WvsQCD"  , "V#rightarrowJ deepAK8: W vs QCD;W vs QCD"    , 40,0,1, candVToJ_.deepAK8()    .WvsQCD, theWeight);
	theHistograms->fill("VToJ_deepAK8_ZvsQCD"  , "V#rightarrowJ deepAK8: Z vs QCD;Z vs QCD"    , 40,0,1, candVToJ_.deepAK8()    .ZvsQCD, theWeight);
	theHistograms->fill("VToJ_deepAK8MD_TvsQCD", "V#rightarrowJ deepAK8 MD: T vs QCD;T vs QCD" , 40,0,1, candVToJ_.deepAK8_MD() .TvsQCD, theWeight);
	theHistograms->fill("VToJ_deepAK8MD_WvsQCD", "V#rightarrowJ deepAK8 MD: W vs QCD;W vs QCD" , 40,0,1, candVToJ_.deepAK8_MD() .WvsQCD, theWeight);
	theHistograms->fill("VToJ_deepAK8MD_ZvsQCD", "V#rightarrowJ deepAK8 MD: Z vs QCD;Z vs QCD" , 40,0,1, candVToJ_.deepAK8_MD() .ZvsQCD, theWeight);
      }
      else{
	if(fakeVTojj_.isValid()){
	  theHistograms->fill("VTojjFake_mass", "V#rightarrowjj mass;GeV/c^{2}", 30, 0, 180, fakeVTojj_.mass(), theWeight);
	  theHistograms->fill("VTojjFake_pt"  , "V#rightarrowjj p_{T};GeV/c"   , 25, 0, 500, fakeVTojj_.pt()  , theWeight);
	  theHistograms->fill("VTojjFake_eta" , "V#rightarrowjj #eta;#eta"     , 25, -5,  5, fakeVTojj_.eta() , theWeight);
	}
	else if(fakeVToJ_.isValid()){
	  theHistograms->fill("VToJFake_mass" , "V#rightarrowJ mass;GeV/c^{2}" , 30, 0, 180, fakeVToJ_.mass(), theWeight);
	  theHistograms->fill("VToJFake_pt"   , "V#rightarrowJ p_{T};GeV/c"    , 25, 0, 500, fakeVToJ_.pt()  , theWeight);
	  theHistograms->fill("VToJFake_eta"  , "V#rightarrowJ #eta;#eta"      , 25, -5,  5, fakeVToJ_.eta() , theWeight);
	  theHistograms->fill("VToJFake_PNet_TvsQCD"     , "V#rightarrowJ ParticleNet: T vs QCD;T vs QCD", 40,0,1, candVToJ_.particleNet().TvsQCD, theWeight);
	  theHistograms->fill("VToJFake_PNet_WvsQCD"     , "V#rightarrowJ ParticleNet: W vs QCD;W vs QCD", 40,0,1, candVToJ_.particleNet().WvsQCD, theWeight);
	  theHistograms->fill("VToJFake_PNet_ZvsQCD"     , "V#rightarrowJ ParticleNet: Z vs QCD;Z vs QCD", 40,0,1, candVToJ_.particleNet().ZvsQCD, theWeight);
	  theHistograms->fill("VToJFake_deepAK8_TvsQCD"  , "V#rightarrowJ deepAK8: T vs QCD;T vs QCD"    , 40,0,1, candVToJ_.deepAK8()    .TvsQCD, theWeight);
	  theHistograms->fill("VToJFake_deepAK8_WvsQCD"  , "V#rightarrowJ deepAK8: W vs QCD;W vs QCD"    , 40,0,1, candVToJ_.deepAK8()    .WvsQCD, theWeight);
	  theHistograms->fill("VToJFake_deepAK8_ZvsQCD"  , "V#rightarrowJ deepAK8: Z vs QCD;Z vs QCD"    , 40,0,1, candVToJ_.deepAK8()    .ZvsQCD, theWeight);
	  theHistograms->fill("VToJFake_deepAK8MD_TvsQCD", "V#rightarrowJ deepAK8 MD: T vs QCD;T vs QCD" , 40,0,1, candVToJ_.deepAK8_MD() .TvsQCD, theWeight);
	  theHistograms->fill("VToJFake_deepAK8MD_WvsQCD", "V#rightarrowJ deepAK8 MD: W vs QCD;W vs QCD" , 40,0,1, candVToJ_.deepAK8_MD() .WvsQCD, theWeight);
	  theHistograms->fill("VToJFake_deepAK8MD_ZvsQCD", "V#rightarrowJ deepAK8 MD: Z vs QCD;Z vs QCD" , 40,0,1, candVToJ_.deepAK8_MD() .ZvsQCD, theWeight);
	}
      }
      
      // Photon region plots with Z->ll
      if( kinPhotons_["central"]->size() > 0){
	theHistograms->fill(  "Z_mass_kinG"   , "Z mass with kin #gamma"   , 30, 60, 120 , Zcand.mass(),theWeight);
	theHistograms->fill(  "ZG_mass_kinG"  , "ZG mass with kin #gamma"  , mVVG_bins, (Zcand.p4() + kinPhotons_["central"]->front().p4()).M(), theWeight);
      }
      else
	theHistograms->fill(  "Z_mass_noG"    , "Z mass without #gamma"    , 30, 60, 120 , Zcand.mass(),theWeight);
      if(loosePhotons_["central"]->size() > 0){
	theHistograms->fill(  "Z_mass_looseG" , "Z mass with loose #gamma" , 30, 60, 120 , Zcand.mass(),theWeight);
	theHistograms->fill(  "ZG_mass_looseG", "ZG mass with loose #gamma", mVVG_bins, (Zcand.p4() + loosePhotons_["central"]->front().p4()).M(), theWeight);
	if(goodPhotons_["central"]->size() > 0){
	  theHistograms->fill("Z_mass_tightG" , "Z mass with tight #gamma" , 30, 60, 120 , Zcand.mass(),theWeight);
	  theHistograms->fill("ZG_mass_tightG", "ZG mass with tight #gamma", mVVG_bins, (Zcand.p4() + goodPhotons_["central"]->front().p4()).M(), theWeight);
	}
	else{
	  theHistograms->fill("Z_mass_failG"  , "Z mass with loose && !tight #gamma"  , 30, 60, 120 , Zcand.mass(),theWeight);
	  theHistograms->fill("ZG_mass_failG" , "ZG mass with loose && !tight #gamma" , mVVG_bins, (Zcand.p4() + kinPhotons_["central"]->front().p4()).M(), theWeight);
	}
      }
    }
  }
  else if(LFR_lep){
    theHistograms->fill("ZL_mass"               , "m_{3l};GeV/c^{2}", 25,0.,500. , (ZL->first.p4()+ZL->second.p4()).M(), theWeight);
    theHistograms->fill("Z_mass"                , "m_{Z};GeV/c^{2}" , 35,55.,125., ZL->first.mass()                    , theWeight);
    theHistograms->fill("Z_l0_pt"               , "p_{t,l00};GeV/c" , 20,0.,400. , ZL->first.daughter(0).pt()          , theWeight);
    theHistograms->fill("Z_l1_pt"               , "p_{t,l00};GeV/c" , 20,0.,400. , ZL->first.daughter(1).pt()          , theWeight);
    theHistograms->fill("L_pt"                  , "p_{t,l3};GeV/c"  , 20,0.,400. , ZL->second.pt()                     , theWeight);

    theHistograms->fill("ZL_mass_" +channelReco_, "m_{3l};GeV/c^{2}", 25,0.,500. , (ZL->first.p4()+ZL->second.p4()).M(), theWeight);
    theHistograms->fill("Z_mass_"  +channelReco_, "m_{Z};GeV/c^{2}" , 35,55.,125., ZL->first.mass()                    , theWeight);
    theHistograms->fill("Z_l0_pt_" +channelReco_, "p_{t,l00};GeV/c" , 20,0.,400. , ZL->first.daughter(0).pt()          , theWeight);
    theHistograms->fill("Z_l1_pt_" +channelReco_, "p_{t,l00};GeV/c" , 20,0.,400. , ZL->first.daughter(1).pt()          , theWeight);
    theHistograms->fill("L_pt_"    +channelReco_, "p_{t,l3};GeV/c"  , 20,0.,400. , ZL->second.pt()                     , theWeight);
  }
  
  /*
    theHistograms->fill("nZtoChLep"    , "Number of Z->ll per event" , 7,0,7, genVBHelper_.ZtoChLep().size());
    theHistograms->fill("nZtoNeutrinos", "Number of Z->nn per event" , 7,0,7, genVBHelper_.ZtoNeutrinos().size());
    theHistograms->fill("nWtoLep"      , "Number of W->lnu per event", 7,0,7, genVBHelper_.WtoLep().size());
    theHistograms->fill("nZtoQ"        , "Number of Z->qq per event" , 7,0,7, genVBHelper_.ZtoQ().size());
    theHistograms->fill("nWtoQ"        , "Number of W->qq' per event", 7,0,7, genVBHelper_.WtoQ().size());
  */
  
  // Photon region plots + V leptonic
  if(kinPhotons_["central"]->size() > 0){
    if     (four_lep){
      theHistograms->fill("ZZ_mass_kinG"  , "ZZ mass with kin #gamma"  , mVV_bins , ZZ->mass(), theWeight);
      theHistograms->fill("ZZG_mass_kinG" , "ZZG mass with kin #gamma" , mVVG_bins, (ZZ->p4()+kinPhotons_["central"]->front().p4()).M(), theWeight);
    }
    else if(three_lep){
      theHistograms->fill("ZW_massT_kinG" , "ZW massT with kin #gamma" , mVV_bins , ZW->p4().Mt(),theWeight);
      theHistograms->fill("ZWG_massT_kinG", "ZWG massT with kin #gamma", mVVG_bins, (ZW->p4()+kinPhotons_["central"]->front().p4()).Mt(), theWeight);
      if(b_WZpaperSel){
	theHistograms->fill("paperSel_ZW_massT_kinG" , "ZW massT with kin #gamma" , mVV_bins , ZW->p4().Mt(),theWeight);
	theHistograms->fill("paperSel_ZWG_massT_kinG", "ZWG massT with kin #gamma", mVVG_bins, (ZW->p4()+kinPhotons_["central"]->front().p4()).Mt(), theWeight);
      }
    }
  }
  // No photon
  else {
    if     (four_lep){
      theHistograms->fill("ZZ_mass_noG"  , "ZZ mass without #gamma"  , mVV_bins, ZZ->mass(), theWeight);
    }
    else if(three_lep){
      theHistograms->fill("ZW_massT_noG" , "ZW massT without #gamma" , mVV_bins, ZW->p4().Mt(),theWeight);
      if(b_WZpaperSel)
	theHistograms->fill("paperSel_ZW_massT_noG" , "ZW massT without #gamma" , mVV_bins, ZW->p4().Mt(),theWeight);
    }
  }
  // Loose photons
  if(loosePhotons_["central"]->size() > 0){
    if     (four_lep){
      theHistograms->fill("ZZ_mass_looseG"  , "ZZ mass with loose #gamma"  , mVV_bins , ZZ->mass(), theWeight);
      theHistograms->fill("ZZG_mass_looseG" , "ZZG mass with loose #gamma" , mVVG_bins, (ZZ->p4()+loosePhotons_["central"]->front().p4()).M(), theWeight);
    }
    else if(three_lep){
      theHistograms->fill("ZW_massT_looseG" , "ZW massT with loose #gamma" , mVV_bins , ZW->p4().Mt(),theWeight);
      theHistograms->fill("ZWG_massT_looseG", "ZWG massT with loose #gamma", mVVG_bins, (ZW->p4()+loosePhotons_["central"]->front().p4()).Mt(), theWeight);
      if(b_WZpaperSel){
        theHistograms->fill("paperSel_ZW_massT_looseG" , "ZW massT with loose #gamma" , mVV_bins , ZW->p4().Mt(),theWeight);
	theHistograms->fill("paperSel_ZWG_massT_looseG", "ZWG massT with loose #gamma", mVVG_bins, (ZW->p4()+loosePhotons_["central"]->front().p4()).Mt(), theWeight);
      }
    }
    // Tight photons
    if(goodPhotons_["central"]->size() > 0){
      if     (four_lep){
	theHistograms->fill("ZZ_mass_tightG"  , "ZZ mass with tight #gamma"  , mVV_bins , ZZ->mass(), theWeight);
	theHistograms->fill("ZZG_mass_tightG" , "ZZG mass with tight #gamma" , mVVG_bins, (ZZ->p4()+goodPhotons_["central"]->front().p4()).M(), theWeight);
      }
      else if(three_lep){
	theHistograms->fill("ZW_massT_looseG" , "ZW massT with tight #gamma" , mVV_bins , ZW->p4().Mt(),theWeight);
	theHistograms->fill("ZWG_massT_looseG", "ZWG massT with tight #gamma", mVVG_bins, (ZW->p4()+goodPhotons_["central"]->front().p4()).Mt(), theWeight);
	if(b_WZpaperSel){
	  theHistograms->fill("paperSel_ZW_massT_looseG" , "ZW massT with tight #gamma" , mVV_bins , ZW->p4().Mt(),theWeight);
	  theHistograms->fill("paperSel_ZWG_massT_looseG", "ZWG massT with tight #gamma", mVVG_bins, (ZW->p4()+goodPhotons_["central"]->front().p4()).Mt(), theWeight);
	}
      }
    }
    // Fail photon (loose && !tight)
    else{
      if     (four_lep){
	theHistograms->fill("ZZ_mass_failG"  , "ZZ mass with loose && !tight #gamma"  , mVV_bins , ZZ->mass(), theWeight);
	theHistograms->fill("ZZG_mass_failG" , "ZZG mass with loose && !tight #gamma" , mVVG_bins, (ZZ->p4()+loosePhotons_["central"]->front().p4()).M(), theWeight);
      }
      else if(three_lep){
	theHistograms->fill("ZW_massT_failG" , "ZW massT with loose && !tight #gamma" , mVV_bins , ZW->p4().Mt(),theWeight);
	theHistograms->fill("ZWG_massT_failG", "ZWG massT with loose && !tight #gamma", mVVG_bins, (ZW->p4()+loosePhotons_["central"]->front().p4()).Mt(), theWeight);
	if(b_WZpaperSel){
	  theHistograms->fill("paperSel_ZW_massT_failG" , "ZW massT with loose && !tight #gamma" , mVV_bins , ZW->p4().Mt(),theWeight);
	  theHistograms->fill("paperSel_ZWG_massT_failG", "ZWG massT with loose && !tight #gamma", mVVG_bins, (ZW->p4()+loosePhotons_["central"]->front().p4()).Mt(), theWeight);
	}
      }
      theHistograms->fill("failG_pt", "Loose && !tight #gamma p_{T}", 15, 0., 300, kinPhotons_["central"]->front().pt(), theWeight);
    }
  }
  

  //int nVBs = genVBHelper_.ZtoChLep().size() + genVBHelper_.ZtoNeutrinos().size() + genVBHelper_.WtoLep().size() + genVBHelper_.ZtoQ().size() + genVBHelper_.WtoQ().size();
  //theHistograms->fill("nVBs", "Number of VB per event", 7,0,7, nVBs);

  // Check Lepton fake rate
  // if(is4Lregion(region_)){
  //   for(const Lepton& lep : *leptons_)
  //     theHistograms->fill("leptonsFromZZ_fakeRateSFUnc", ";fakeRateSFUnc", 40, 0., 2., lep.fakeRateSFUnc());
  //   for(const Boson<Lepton>& bos : {ZZ->first(), ZZ->second()}){
  //     theHistograms->fill("bosonsFromZZ_fakeRateSFUnc", ";fakeRateSFUnc", 40, 0., 2., bos.fakeRateSFUnc());
  //   }
  //   theHistograms->fill("ZZ_fakeRateSFUnc", ";fakeRateSFUnc", 40, 0., 2., ZZ->fakeRateSFUnc());
  // }
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
  for(const char* name : {"AAA_cuts", "AAA_cuts_u", "kinPhotons_ID", "ph_eScale_N", "kinPhotons_cutflow"}){
    TH1* h = theHistograms->get(name);
    if(h) h->LabelsDeflate("X");
  }
  
  TH1* AAA_cuts   = theHistograms->get("AAA_cuts"  );
  TH1* AAA_cuts_u = theHistograms->get("AAA_cuts_u");
  for(TH1* h : {AAA_cuts, AAA_cuts_u}){
    if(!h) continue;
    TAxis* axis = h->GetXaxis();
    axis->SetBinLabel(1, "All");
    axis->SetBinLabel(2, "ZZ || ZW");
    axis->SetBinLabel(3, "2l2j || 2l1J");
    axis->SetBinLabel(4, "#gamma kin");
    axis->SetBinLabel(5, "#gamma loose");
    axis->SetBinLabel(6, "#gamma good");
    axis->SetBinLabel(7, "Analyzed");
    axis->SetBinLabel(8, "#gamma medium");
  }
  
  TH1* kinPhotons_cutflow = theHistograms->get("kinPhotons_cutflow");
  TH1* kinPhotons_Nm1     = theHistograms->get("kinPhotons_Nm1"    );
  if(kinPhotons_cutflow) kinPhotons_cutflow->GetXaxis()->SetBinLabel(1, "No #gamma" );
  if(kinPhotons_Nm1    ) kinPhotons_Nm1    ->GetXaxis()->SetBinLabel(1, "All #gamma");
  for(TH1* h : {kinPhotons_cutflow, kinPhotons_Nm1}){
    if(!h) continue;
    TAxis* axis = h->GetXaxis();
    axis->SetBinLabel(2, "#sigma_{i#etai#eta}");
    axis->SetBinLabel(3, "HoverE");
    axis->SetBinLabel(4, "IsoCH");
    axis->SetBinLabel(5, "IsoNE");
    axis->SetBinLabel(6, "IsoPh");
  }
  
  TH1* kinPhotons_ID = theHistograms->get("kinPhotons_ID");
  if(kinPhotons_ID){
    TAxis* axis = kinPhotons_ID->GetXaxis();
    axis->SetBinLabel(1, "Total");
    axis->SetBinLabel(2, "Loose");
    axis->SetBinLabel(3, "Medium");
    axis->SetBinLabel(4, "Tight");
  }
  
  TH1* WZ_cutflow = theHistograms->get("WZ_cutflow");
  if(WZ_cutflow){
    TAxis* xaxis = WZ_cutflow->GetXaxis();
    xaxis->SetBinLabel(1, "Total"                );
    xaxis->SetBinLabel(2, "lepton p_{t}"         );
    xaxis->SetBinLabel(3, "|m_{ll} - m_{Z}| < 15");
    xaxis->SetBinLabel(4, "MET > 30"             );
    xaxis->SetBinLabel(5, "m_{lll} > 100"        );
    xaxis->SetBinLabel(6, "bveto"                );
    xaxis->SetBinLabel(7, "All cuts"             );
  }
  
  TH1* channel_lep = theHistograms->get("channel_lep");
  if(channel_lep){
    TAxis* axis = channel_lep->GetXaxis();
    if     (is4Lregion(region_)){
      axis->SetBinLabel(1, "4e"  );
      axis->SetBinLabel(2, "2e2m");
      axis->SetBinLabel(3, "4m"  );
    }
    else if(is3Lregion(region_)){
      axis->SetBinLabel(1, "3e"  );
      axis->SetBinLabel(2, "2e1m");
      axis->SetBinLabel(3, "2m1e");
      axis->SetBinLabel(4, "3m"  );
    }
    else if(is2Lregion(region_)){
      axis->SetBinLabel(1, "2e");
      axis->SetBinLabel(2, "2m");
    }
    else if(region_ == CRLFR){
      axis->SetBinLabel(1, "2e+e");
      axis->SetBinLabel(2, "2m+e");
      axis->SetBinLabel(3, "2e+m");
      axis->SetBinLabel(4, "2m+m");
    }
  }
}


void VVGammaAnalyzer::finish(){
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
    FILE* cherryFile = fopen(Form("data/2016D_%s_manual.txt", regionType(R).c_str()), "r");
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


template <class T, class UnaryPredicate>
vector<Boson<T>> VVGammaAnalyzer::makeBosons(const std::vector<T>& vec, UnaryPredicate selection){
  if(vec.size() < 2)
    return vector<Boson<T>>();

  vector<Boson<T>> out;
  out.reserve(vec.size() * (vec.size()-1) / 2);;
  for  (auto i = vec.cbegin(); i != vec.end(); ++i)
    for(auto j = i+1         ; j != vec.end(); ++j){
      Boson<T> b(*i, *j, 1);  // last parameter is ID; if not set to something != 0 will cause isValid() to return false and me to waste an afternoon debugging
      if(selection(b))
	out.push_back(b);
    }
  return out;
}


vector<Boson<Jet>> VVGammaAnalyzer::candidatesVTojj(const std::vector<phys::Jet>& theJets){
  vector<Boson<Jet>> pairs = makeBosons(theJets, 
					[](const Boson<Jet> b){ return minDM(b.mass()) < 30 ;}
					);
  return pairs;
}


void VVGammaAnalyzer::hadronicObjectsReconstruction(){
  vector<Boson<Jet>> candsjj = candidatesVTojj(*jets);
  
  if(candsjj.size() > 0){
    std::sort(candsjj.begin(), candsjj.end(), Mass2Comparator(phys::ZMASS, phys::WMASS));
    candVTojj_ = Boson<Jet>(std::move(candsjj.front()));
  }
  else if(jetsAK8->size() > 0){
    auto it = std::min_element(jetsAK8->begin(), jetsAK8->end(), Mass2Comparator(phys::ZMASS, phys::WMASS));
    if(minDM(it->mass()) < 30)
      candVToJ_ = *it;
    else
      fakeVToJ_ = *it;
  }
  else{  // try to find pairs for CRs
    if(jets->size() >= 2){
      vector<Boson<Jet>> allPairs = makeBosons(*jets, [](const Boson<Jet>&){ return true; } );
      std::sort(allPairs.begin(), allPairs.end(), Mass2Comparator(phys::ZMASS, phys::WMASS));
      fakeVTojj_ = Boson<Jet>(std::move(allPairs.front()));
    }
  }
}


void VVGammaAnalyzer::makeChannelReco(){
  // sets channelReco_ for each event
  unsigned int nEl(0), nMu(0);  
  for(const Lepton l : *leptons_){
    switch( abs(l.id()) ){
    case 11:
      nEl++;
      break;
    case 13:
      nMu++;
      break;
    default:
      cout<<"Error: In region: "<<phys::regionType(region_)<<" found lept from ZZ/ZW/ZL with ID: "<<l.id()<<" \taddr: "<<l<<'\n';
      if(ZZ) cout<<'\t'<<ZZ->first().daughter(0).id()<<" ; "<<ZZ->first().daughter(1).id()<<" ; "<<ZZ->second().daughter(0).id()<<" ; "<<ZZ->second().daughter(1).id()<<'\n';
      if(ZW) cout<<ZW->first().daughter(0).id()<<" ; "<<ZW->first().daughter(1).id()<<" ; "<<ZW->second().daughter(0).id()<<'\n';
    }
  }
  
  channelReco_ = "??"; // = Form("%ue%um", nEl, nMu);
  
  if     (is4Lregion(region_)){
    if     (nEl == 4)             { channelReco_ = "4e"  ; theHistograms->fill("channel_lep", "Lepton channel", 3,-0.5,2.5, 0, theWeight); }
    else if(nEl == 2 && nMu == 2) { channelReco_ = "2e2m"; theHistograms->fill("channel_lep", "Lepton channel", 3,-0.5,2.5, 1, theWeight); }
    else if(nMu == 4)             { channelReco_ = "4m"  ; theHistograms->fill("channel_lep", "Lepton channel", 3,-0.5,2.5, 2, theWeight); }
  }
  else if(is3Lregion(region_)){
    if     (nEl == 3)             { channelReco_ = "3e"  ; theHistograms->fill("channel_lep", "Lepton channel", 4,-0.5,3.5, 0, theWeight); }
    else if(nEl == 2 && nMu == 1) { channelReco_ = "2e1m"; theHistograms->fill("channel_lep", "Lepton channel", 4,-0.5,3.5, 1, theWeight); }
    else if(nEl == 1 && nMu == 2) { channelReco_ = "2m1e"; theHistograms->fill("channel_lep", "Lepton channel", 4,-0.5,3.5, 2, theWeight); }
    else if(nMu == 3)             { channelReco_ = "3m"  ; theHistograms->fill("channel_lep", "Lepton channel", 4,-0.5,3.5, 3, theWeight); }
  }
  else if(is2Lregion(region_)){
    if     (nEl == 2)             { channelReco_ = "2e"  ; theHistograms->fill("channel_lep", "Lepton channel", 2,-0.5,1.5, 0, theWeight); }
    else if(nMu == 2)             { channelReco_ = "2m"  ; theHistograms->fill("channel_lep", "Lepton channel", 2,-0.5,1.5, 1, theWeight); }
  }
  else if(region_ == CRLFR){
    if     (nEl == 3)             { channelReco_ = "2e+e"; theHistograms->fill("channel_lep", "Lepton channel", 4,-0.5,3.5, 0, theWeight); }
    else if(nEl == 1 && nMu == 2) { channelReco_ = "2m+e"; theHistograms->fill("channel_lep", "Lepton channel", 4,-0.5,3.5, 1, theWeight); }
    else if(nEl == 2 && nMu == 1) { channelReco_ = "2e+m"; theHistograms->fill("channel_lep", "Lepton channel", 4,-0.5,3.5, 2, theWeight); }
    else if(nMu == 3)             { channelReco_ = "2m+m"; theHistograms->fill("channel_lep", "Lepton channel", 4,-0.5,3.5, 3, theWeight); }
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


// const phys::Jet* candAK8(const std::vector<phys::Jet>* collection) const{
//   return &(*std::min_element(collection->begin(), collection->end(),
// 			     [](const phys::Jet& a, const){
// 			       return minDM();
// 				 })
// 	   );
// }


void VVGammaAnalyzer::baseHistos_cut(){
  // Leptons
  auto lead_ele = std::max_element(electrons->begin(), electrons->end(), [](const Lepton& a, const Lepton& b){ return a.pt() < b.pt(); } );
  auto lead_muo = std::max_element(muons->begin()    , muons->end()    , [](const Lepton& a, const Lepton& b){ return a.pt() < b.pt(); } );
  if(lead_ele != electrons->end())
    theHistograms->fill("lead_ele_pt", "Leading electron p_{t}:p_{t} [GeV/c]", 50, 0., 250., lead_ele->pt(), theWeight);
  if(lead_muo != muons->end())
    theHistograms->fill("lead_muo_pt", "Leading muoon p_{t}:p_{t} [GeV/c]"   , 50, 0., 250., lead_muo->pt(), theWeight);
  
  // MET
  theHistograms->fill("MET", "MET;#slash{E}_{T} [GeV/c]", 50, 0., 250., met->pt(), theWeight);
}


void VVGammaAnalyzer::photonHistos(){
  // No photons passing kinematic cuts
  if(theSampleInfo.isMC() && kinPhotons_["central"]->size() == 0 && genPhotons_->size() > 0){
    theHistograms->fill("noKinPh_all_genPh_N" , "Number of #gamma_{GEN} when no #gamma_{REC};# #gamma_{GEN}", 5,-0.5,4.5, genPhotons_->size(), theWeight);
    for(const Particle& genPh : *genPhotons_){
      theHistograms->fill("noKinPh_all_genPh_pt" , "#gamma_{GEN} p_{T} when no #gamma_{REC};#gamma_{GEN} p_{T}", 50,0. ,200., genPh.pt() , theWeight);
      theHistograms->fill("noKinPh_all_genPh_eta", "#gamma_{GEN} #eta when no #gamma_{REC};#gamma_{GEN} #eta"  , 50,-5.,  5., genPh.eta(), theWeight);
      auto it_ph = std::min_element(photons->begin(), photons->end(), DeltaRComparator(genPh));
      if(it_ph != photons->end() && deltaR(*it_ph, genPh) < 0.2){
	theHistograms->fill("noKinPh_rec_genPh_pt" , "#gamma_{GEN} p_{T} when no #gamma_{REC};#gamma_{GEN} p_{T}", 50,0. ,200., genPh.pt() , theWeight);
	theHistograms->fill("noKinPh_rec_genPh_eta", "#gamma_{GEN} #eta when no #gamma_{REC};#gamma_{GEN} #eta"  , 50,-5.,  5., genPh.eta(), theWeight);
      }
    }
  }
  
  // Photons passing kinematic cuts; pt, eta, ID variables
  for(const Photon& ph : *kinPhotons_["central"]){
    theHistograms->fill("kinPh_pt"  , "p_{t} of kinPhotons" , 50, 0., 250.,      ph.pt()  , theWeight);
    theHistograms->fill("kinPh_aeta", "|#eta| of kinPhotons", 50, 0., 2.5 , fabs(ph.eta()), theWeight);

    if(ph.isBarrel()){
      vector<double> sieie_bins(30);
      for(size_t i = 0; i < sieie_bins.size(); i++) sieie_bins[i] = 0.004 + 0.001*i;
      vector<double> chIso_bins({0., 0.1, 0.25, 0.65, 0.8, 1.141, 1.3, 1.694, 2, 3, 4});
      theHistograms->fill("kinPh_sieie_chIso_EB", "kinPhotons in Barrel;#sigma_{i#etai#eta};chIso",
    			  sieie_bins,
    			  chIso_bins,
    			  ph.sigmaIetaIeta(), ph.chargedIsolation(), theWeight);
      theHistograms->fill("kinPh_sieie_EB" , "kinPhotons in Barrel;#sigma_{i#etai#eta}", sieie_bins, ph.sigmaIetaIeta()         , theWeight);
      theHistograms->fill("kinPh_chIso_EB" , "kinPhotons in Barrel;chIso"              , chIso_bins, ph.chargedIsolation()      , theWeight);
      theHistograms->fill("kinPh_HoverE_EB", "kinPhotons in Barrel;HoverE"             , 75,0.,0.15, ph.HoverE()                , theWeight);
      theHistograms->fill("kinPh_neIso_EB" , "kinPhotons in Barrel;neIso"              , 180,0., 30, ph.neutralHadronIsolation(), theWeight);
      theHistograms->fill("kinPh_phIso_EB" , "kinPhotons in Barrel;phIso"              , 200,0.,100, ph.photonIsolation()       , theWeight);
    }
    else{
      vector<double> sieie_bins(28);
      for(size_t i = 0; i < sieie_bins.size(); i++) sieie_bins[i] = 0.01 + 0.0025*i;
      vector<double> chIso_bins({0., 0.1, 0.25, 0.517, 0.8, 1.051, 1.3, 1.6, 2.089, 3, 4});
      theHistograms->fill("kinPh_sieie_chIso_EE", "kinPhotons in Endcap;#sigma_{i#etai#eta};chIso",
			  sieie_bins,
			  chIso_bins,
			  ph.sigmaIetaIeta(), ph.chargedIsolation(), theWeight);
      theHistograms->fill("kinPh_sieie_EE" , "kinPhotons in Endcap;#sigma_{i#etai#eta}", sieie_bins, ph.sigmaIetaIeta()         , theWeight);
      theHistograms->fill("kinPh_chIso_EE" , "kinPhotons in Endcap;chIso"              , chIso_bins, ph.chargedIsolation()      , theWeight);
      theHistograms->fill("kinPh_HoverE_EE", "kinPhotons in Endcap;HoverE"             , 75,0.,0.15, ph.HoverE()                , theWeight);
      theHistograms->fill("kinPh_neIso_EE" , "kinPhotons in Endcap;neIso"              , 200,0., 40, ph.neutralHadronIsolation(), theWeight);
      theHistograms->fill("kinPh_phIso_EE" , "kinPhotons in Endcap;phIso"              , 200,0.,100, ph.photonIsolation()       , theWeight);
    }
    
    theHistograms->fill("kinPhotons_ID", "Cut Based ID", BINS_KINPHID, 0, theWeight);
    if(ph.cutBasedIDLoose())  theHistograms->fill("kinPhotons_ID", "Cut Based ID", BINS_KINPHID, 1, theWeight);
    if(ph.cutBasedIDMedium()) theHistograms->fill("kinPhotons_ID", "Cut Based ID", BINS_KINPHID, 2, theWeight);
    if(ph.cutBasedIDTight())  theHistograms->fill("kinPhotons_ID", "Cut Based ID", BINS_KINPHID, 3, theWeight);
  }
  
  // Systematics
  for(auto & [syst, phVect] : kinPhotons_)
    if(phVect->size() > 0)
      // Creating void histograms, then filling alphanumeric labels --> new ones are created as they are encountered
      theHistograms->book<TH1F>("ph_eScale_N", "Number of #gamma passing selection", 1,0,0)->Fill(Form("kin_%s" , syst), theWeight);
  for(auto & [syst, phVect] : goodPhotons_)
    if(phVect->size() > 0)
      theHistograms->book<TH1F>("ph_eScale_N", "Number of #gamma passing selection", 1,0,0)->Fill(Form("good_%s", syst), theWeight);
  
  // How important are subleading photons?
  for(const auto& [syst, phVect] : goodPhotons_)
    theHistograms->fill(Form("tightPhotons_%s_N", syst), Form("Number of tight photons (%s)", syst), 5,-0.5,4.5, phVect->size(), theWeight);
  if(goodPhotons_["central"]->size() > 0)
    theHistograms->fill("TighPh_lead_pt", "Leading tight #gamma;p_{T} [GeV/c]", 20,0.,100., goodPhotons_["central"]->at(0).pt());
  if(goodPhotons_["central"]->size() > 1)
    theHistograms->fill("subleadTighPh_pt", "Subleading tight #gamma;p_{T} [GeV/c]", 20,0.,100., goodPhotons_["central"]->at(1).pt());
  
  for(const auto& [syst, phVect] : loosePhotons_)
    theHistograms->fill(Form("loosePhotons_%s_N", syst), Form("Number of loose photons (%s)", syst), 5,-0.5,4.5, phVect->size(), theWeight);
  if(loosePhotons_["central"]->size() > 0)
    theHistograms->fill("LoosePh_lead_pt", "Leading loose #gamma;p_{T} [GeV/c]", 20,0.,100., loosePhotons_["central"]->at(0).pt());
  if(loosePhotons_["central"]->size() > 1)
    theHistograms->fill("subleadLoosePh_pt", "Subleading loose #gamma;p_{T} [GeV/c]", 20,0.,100., loosePhotons_["central"]->at(1).pt());

  for(const auto& [syst, phVect] : kinPhotons_)
    theHistograms->fill(Form(  "kinPhotons_%s_N", syst), Form("Number of kin photons (%s)"  , syst), 5,-0.5,4.5, phVect->size(), theWeight);
  if( kinPhotons_["central"]->size() > 0)
    theHistograms->fill("KinPh_lead_pt", "Leading kin #gamma;p_{T} [GeV/c]", 20,0.,100., kinPhotons_["central"]->at(0).pt());
  if( kinPhotons_["central"]->size() > 1)
    theHistograms->fill(  "subleadKinPh_pt", "Subleading kin #gamma;p_{T} [GeV/c]"  , 20,0.,100.,  kinPhotons_["central"]->at(1).pt());
  
  // Test data/MC for each cut separately
  if(kinPhotons_["central"]->size() == 0){
    theHistograms->fill("kinPhotons_cutflow", "Cutflow of photons;;Events", BINS_PHCUTFLOW, 0, theWeight);  // No photons
    return;
  }
  else{
    // Select photon that passes more cuts
    Photon::IDwp wp = Photon::IDwp::Loose;
    auto bestG2 =std::max_element(kinPhotons_["central"]->begin(), kinPhotons_["central"]->end(),
				 [wp](const Photon& a, const Photon& b){
				   int na = a.passSigmaiEtaiEta(wp)+a.passHoverE(wp)+a.passChargedIsolation(wp)+a.passNeutralIsolation(wp)+a.passPhotonIsolation(wp);
				   int nb = b.passSigmaiEtaiEta(wp)+b.passHoverE(wp)+b.passChargedIsolation(wp)+b.passNeutralIsolation(wp)+b.passPhotonIsolation(wp);
				   return na < nb;
				 });  // max_element returns the first among those with max value --> preserve pt ordering
    const Photon* bestG = nullptr;
    int best_ncuts = -1;
    for(const Photon& ph : *kinPhotons_["central"]){
      int ncuts = ph.passSigmaiEtaiEta(wp)+ph.passHoverE(wp)+ph.passChargedIsolation(wp)+ph.passNeutralIsolation(wp)+ph.passPhotonIsolation(wp);
      if(ncuts > best_ncuts){
	bestG = &ph;
	best_ncuts = ncuts;
      }
      if(ncuts == 5)
	break;
    }
    if(bestG != &*bestG2)
      cout << ">>> max_element and loop returned different photons: " << bestG << ' ' << &*bestG2 << '\n';
    bool b_sieie  = bestG->passSigmaiEtaiEta(wp);
    bool b_HoverE = bestG->passHoverE(wp);
    bool b_chIso  = bestG->passChargedIsolation(wp);
    bool b_neIso  = bestG->passNeutralIsolation(wp);
    bool b_phIso  = bestG->passPhotonIsolation(wp);
    
    // Single cut efficiency
    if(b_sieie)
      theHistograms->fill("kinPhotons_cutflow", "Cutflow of photons;;Events", BINS_PHCUTFLOW, 1, theWeight);
    if(b_HoverE)
      theHistograms->fill("kinPhotons_cutflow", "Cutflow of photons;;Events", BINS_PHCUTFLOW, 2, theWeight);
    if(b_chIso)
      theHistograms->fill("kinPhotons_cutflow", "Cutflow of photons;;Events", BINS_PHCUTFLOW, 3, theWeight);
    if(b_neIso)
      theHistograms->fill("kinPhotons_cutflow", "Cutflow of photons;;Events", BINS_PHCUTFLOW, 4, theWeight);
    if(b_phIso)
      theHistograms->fill("kinPhotons_cutflow", "Cutflow of photons;;Events", BINS_PHCUTFLOW, 5, theWeight);
    
    // N-1 efficiency of the cuts
    theHistograms->fill("kinPhotons_Nm1", "N-1 cutflow of photons;;Events", BINS_PHCUTFLOW, 0, theWeight);  // All photons
    if(!b_sieie &&  b_HoverE &&  b_chIso &&  b_neIso &&  b_phIso)
      theHistograms->fill("kinPhotons_Nm1", "N-1 cutflow of photons;;Events", BINS_PHCUTFLOW, 1, theWeight);
    if( b_sieie && !b_HoverE &&  b_chIso &&  b_neIso &&  b_phIso)
      theHistograms->fill("kinPhotons_Nm1", "N-1 cutflow of photons;;Events", BINS_PHCUTFLOW, 2, theWeight);
    if( b_sieie &&  b_HoverE && !b_chIso &&  b_neIso &&  b_phIso)
      theHistograms->fill("kinPhotons_Nm1", "N-1 cutflow of photons;;Events", BINS_PHCUTFLOW, 3, theWeight);
    if( b_sieie &&  b_HoverE &&  b_chIso && !b_neIso &&  b_phIso)
      theHistograms->fill("kinPhotons_Nm1", "N-1 cutflow of photons;;Events", BINS_PHCUTFLOW, 4, theWeight);
    if( b_sieie &&  b_HoverE &&  b_chIso &&  b_neIso && !b_phIso)
      theHistograms->fill("kinPhotons_Nm1", "N-1 cutflow of photons;;Events", BINS_PHCUTFLOW, 5, theWeight);
    
    theHistograms->fill("kinPhotons_MVA", "kinematic #gamma MVA", 41,-1.025,1.025, bestG->MVAvalue(), theWeight);
  }  // END if(kinPhoton["central"]->size() == 0)
  
  // Systematics histos
  // for(const Photon& ph : *photons){
  //   theHistograms->fill("ph_E"                , "E;[GeV]"                             , 50,0.,500., ph.e()                           , theWeight);
  //   theHistograms->fill("ph_E_m_ScaleUp"    , "energyScaleUp - E;#DeltaE[GeV]"    , 60,-1.,5.,  ph.energyScaleUp() - ph.e()      , theWeight);
  //   theHistograms->fill("ph_E_m_ScaleDn"    , "E - energyScaleDn;#DeltaE[GeV]"    , 60,-1.,5.,  ph.e() - ph.energyScaleDown()    , theWeight);
  //   theHistograms->fill("ph_E_m_ScaleStatUp", "energyScaleStatUp - E;#DeltaE[GeV]", 60,-1.,5.,  ph.energyScaleStatUp() - ph.e()  , theWeight);
  //   theHistograms->fill("ph_E_m_ScaleStatDn", "E - energyScaleStatDn;#DeltaE[GeV]", 60,-1.,5.,  ph.e() - ph.energyScaleStatDown(), theWeight);
  //   theHistograms->fill("ph_E_m_ScaleSystUp", "energyScaleSystUp - E;#DeltaE[GeV]", 60,-1.,5.,  ph.energyScaleSystUp() - ph.e()  , theWeight);
  //   theHistograms->fill("ph_E_m_ScaleSystDn", "E - energyScaleSystDn;#DeltaE[GeV]", 60,-1.,5.,  ph.e() - ph.energyScaleSystDown(), theWeight);
  //   theHistograms->fill("ph_E_m_ScaleGainUp", "energyScaleGainUp - E#;DeltaE[GeV]", 60,-1.,5.,  ph.energyScaleGainUp() - ph.e()  , theWeight);
  //   theHistograms->fill("ph_E_m_ScaleGainDn", "E - energyScaleGainDn;#DeltaE[GeV]", 60,-1.,5.,  ph.e() - ph.energyScaleGainDown(), theWeight);
  //   theHistograms->fill("ph_E_m_ScaleEtUp"  , "energyScaleEtUp - E;#DeltaE[GeV]"  , 60,-1.,5.,  ph.energyScaleEtUp() - ph.e()    , theWeight);
  //   theHistograms->fill("ph_E_m_ScaleEtDn"  , "E - energyScaleEtDn;#DeltaE[GeV]"  , 60,-1.,5.,  ph.e() - ph.energyScaleEtDown()  , theWeight);
  // }
}


void VVGammaAnalyzer::jetHistos(){
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
    theHistograms->fill("AK4_gen_pt"  , "AK4 generated;p_{t} [GeV/c]"             , 47, 30, 500, gen.pt(), theWeight);
    theHistograms->fill("AK4_gen_pt_u", "AK4 generated (unweighted);p_{t} [GeV/c]", 47, 30, 500, gen.pt());
    
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
  theHistograms->fill(Form("PKU_%s_kinPh_N" , channel_gen.c_str()), Form(title_N.c_str(), "photons passing kinematics")      ,5,-0.5,4.5, kinPhotons_["central"]->size() , 1);
  theHistograms->fill(Form("PKU_%s_goodPh_N", channel_gen.c_str()), Form(title_N.c_str(), "photons passing cutBasedIDLoose") ,5,-0.5,4.5, goodPhotons_["central"]->size(), 1);
  theHistograms->fill(Form("PKU_%s_POGele_N", channel_gen.c_str()), Form(title_N.c_str(), "electrons")                       ,5,-0.5,4.5, electrons->size()   , 1);
  theHistograms->fill(Form("PKU_%s_POGmuo_N", channel_gen.c_str()), Form(title_N.c_str(), "muon")                            ,5,-0.5,4.5, muons->size()       , 1);
  
  std::string title_pt( Form("Z #rightarrow %s: p_{t} of %s;p_{t} [GeV/c]", channel_gen.c_str(), "%s") );
  for(const Photon& ph : *kinPhotons_["central"] )
    theHistograms->fill(Form("PKU_%s_kinPh_pt" , channel_gen.c_str()), Form(title_pt.c_str(), "photons passing kinematics")     , 20,0.,200., ph.pt(), 1);
  for(const Photon& ph : *goodPhotons_["central"])
    theHistograms->fill(Form("PKU_%s_goodPh_pt", channel_gen.c_str()), Form(title_pt.c_str(), "photons passing cutBasedIDLoose"), 20,0.,200., ph.pt(), 1);
  
  for(const Lepton& el : *electrons )
    theHistograms->fill(Form("PKU_%s_POGele_pt", channel_gen.c_str()), Form(title_pt.c_str(), "electrons")                      , 20,0.,200., el.pt(), 1);
  for(const Lepton& mu : *muons)
    theHistograms->fill(Form("PKU_%s_POGmuo_pt", channel_gen.c_str()), Form(title_pt.c_str(), "muons")                          , 20,0.,200., mu.pt(), 1);
}


void VVGammaAnalyzer::leptonFakeRate(){
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


void VVGammaAnalyzer::photonFakeRate(){
  const vector<Photon>& thePhVect = *loosePhotons_["central"];  // *kinPhotons_["central"];
  if(thePhVect.size() == 0)
    return;
  
  Photon::IDwp wp = Photon::IDwp::Loose;
  auto bestG = std::max_element(thePhVect.begin(), thePhVect.end(),
				[wp](const Photon& a, const Photon& b){
				  return a.passSigmaiEtaiEta(wp) + a.passChargedIsolation(wp) < b.passSigmaiEtaiEta(wp) + b.passChargedIsolation(wp);
				});  // max_element returns the first among those with max value --> preserve pt ordering
  bool isPrompt = false;
  if(theSampleInfo.isMC())
    isPrompt = std::any_of(genPhotons_->begin(), genPhotons_->end(), 
			   [bestG](const Particle& gen){ return physmath::deltaR(*bestG, gen) < 0.4; }
			   );
  
  // Actual fake rate histograms
  double theAeta = fabs(bestG->eta());
  double thePt   = bestG->pt();
  if(thePt > ph_pt_bins.back())
    thePt = ph_pt_bins.back() - 0.1;
  
  // Fake rate with ABCD
  if( !bestG->passHoverE(wp) || !bestG->passNeutralIsolation(wp) || !bestG->passPhotonIsolation(wp))
    throw std::logic_error("Inconsistent photon selection");
  
  char ABCD = phABCD(*bestG, wp);
  //TEMP check
  if(ABCD == 'A' && !bestG->cutBasedIDLoose()){
  //   cout << "\tHoverE: " << bestG->passHoverE(wp) << " - sigmaiEtaiEta: " << bestG->passSigmaiEtaiEta(wp) << " - chargedIsolation:" << bestG->passChargedIsolation(wp) << " - neutralIsolation: " << bestG->passNeutralIsolation(wp) << " - photonIsolation: " << bestG->passPhotonIsolation(wp) << '\n';
    throw std::logic_error("Inconsistency in Photon ID!");
  }
  const char* name = Form("PhFR_%c", ABCD);
  if(theSampleInfo.isMC())
    name = Form("%s_%s", name, (isPrompt ? "nonprompt" : "prompt"));
  
  theHistograms->fill(name, Form("Photons: %c;p_{T} [GeV/c];#eta", ABCD), ph_pt_bins, ph_aeta_bins, thePt, theAeta, theWeight);
  
  // Fake rate with Loose-to-Tight (loose = 3/5 cutBasedIDLoose, tight = 5/5 cutBasedIDLoose)
  int nPass = 3 + bestG->passSigmaiEtaiEta(wp) + bestG->passChargedIsolation(wp);
  name = Form("PhFR_LToT_%d", nPass);
  if(theSampleInfo.isMC())
    name = Form("%s_%s", name, (isPrompt ? "nonprompt" : "prompt"));
  
  theHistograms->fill(name, "Photon fake rate with loose to tight;p_{T} [GeV/c];#eta", ph_pt_bins, ph_aeta_bins, thePt, theAeta, theWeight);

  // Closure test
  const char* varName;
  double varValue;
  if     (is4Lregion(region_)){
    varName = "mZZG";
    varValue = ( ZZ->p4() + bestG->p4() ).M();
  }
  else if(is3Lregion(region_)){
    varName = "mWZG";
    varValue = ( ZW->p4() + bestG->p4() ).M();
  }
  else
    return;
  
  double f = getPhotonFR(*bestG); //photonScaleFactors_.fakeRate(*bestG);
  double weight = theWeight * f/(1-f);
  // cout << region_ << ' ' << varName << "=" << varValue << "  theWeight=" << theWeight;
  // cout << "  f=" << f << "  weight=" << weight << "  pt=" << bestG->pt() << "  eta=" << bestG->eta() << '\n';
  
  if     (nPass == 5)   // TIGHT           //bestG->cutBasedIDLoose()
    theHistograms->fill(Form("PhFRClosure_PASS_%s"      , varName), "Closure test: PASS"    , mVVG_bins, varValue, theWeight);
  else if(nPass == 3 || nPass == 4){  // LOOSE && !TIGHT
    theHistograms->fill(Form("PhFRClosure_FAIL_%s"      , varName), "Closure test: FAIL"    , mVVG_bins, varValue, theWeight);
    theHistograms->fill(Form("PhFRClosure_reweighted_%s", varName), "Closure test: (3+4)*SF", mVVG_bins, varValue, weight   );
  }
}


void studyJetsChoice(std::ofstream& fout){
  // // Study how to choose the correct jet AK8 in case it exists. The problem of rejecting 
  // // events in which it does not exist is part of the event selection
  
  // // Assuming genEventSetup has been called for this event
  // if( !theSampleInfo.isMC() )
  //   return;
  
  // // Find the gen quarks
  // if(genQuarks_->size() != 2)
  //   return;
  // Boson<Particle> diquark(genQuarks_->at(0), genQuarks_->at(1));
  
  // // Check that the corresponding gen jet exists
  // vector<Particle>::const_iterator genAK8 = std::min_element(genJetsAK8->begin(), genJetsAK8->end(), DeltaRComparator(diquark));
  // if(genAK8 == genJetsAK8->end() || deltaR(*genAK8, diquark) < 0.4)
  //   return;
  
  // // Check that at least one reco jet exists
  // vector<Particle>::const_iterator recAK8 = std::min_element(jetsAK8->begin(), jetsAK8->end(), DeltaRComparator(*genAK8));
  // if(recAK8 == jetsAK8->end() || deltaR(*recAK8, *genAK8) < 0.4)
  //   return;
  
  // // There should be a file
  // for(vector<Jet>::const_iterator rec_it = jetsAK8->cbegin(); jetsAK8 != jetsAK8.cend(); ++rec_it){
  //   fout << (rec_it == recAK8)           << ',' 
  // 	 << rec_it->mass()               << ','
  // 	 << rec_it->particleNet().WvsQCD << ','
  // 	 << rec_it->particleNet().ZvsQCD << ','
  // 	 << rec_it->particleNet().TvsQCD << ','
  // 	 << rec_it->deepAK8_MD().WvsQCD  << ','
  // 	 << rec_it->deepAK8_MD().ZvsQCD  << ','
  // 	 << rec_it->deepAK8_MD().TvsQCD  << '\n';
  // }
}


void VVGammaAnalyzer::photonEfficiency(const std::vector<phys::Photon>& vPh, const char* label){
  // Reconstruction efficiency and resolution
  if(theSampleInfo.isMC()){
    vector<std::pair<const Particle*, const Photon*>> genRecPh = matchGenRec(*genPhotons_, vPh, 0.4);  // intrinsic threshold of deltaR = 0.2
    for(auto & [gen, rec] : genRecPh){
      theHistograms->fill(Form("%sEff_den_pt" , label), "DEN photons p_{T};p_{T} [GeV/c]", 25,0.,250., gen->pt() , theWeight);
      theHistograms->fill(Form("%sEff_den_E"  , label), "DEN photons energy;E [GeV]"     , 25,0.,500., gen->e()  , theWeight);
      theHistograms->fill(Form("%sEff_den_eta", label), "DEN photons eta;#eta"           , eta_bins  , gen->eta(), theWeight);
      
      if(rec == nullptr)
  	continue;
      theHistograms->fill(Form("%sEff_num_pt" , label), "NUM photons p_{T};p_{T} [GeV/c]", 25,0.,250., gen->pt() , theWeight);
      theHistograms->fill(Form("%sEff_num_E"  , label), "NUM photons energy;E [GeV]"     , 25,0.,500., gen->e()  , theWeight);
      theHistograms->fill(Form("%sEff_num_eta", label), "NUM photons eta;#eta"           , eta_bins  , gen->eta(), theWeight);
      
      double deltaR = physmath::deltaR(*rec, *gen);
      double deltaEoverE = (rec->e() - gen->e()) / gen->e();
      double deltapToverpT = (rec->pt() - gen->pt()) / gen->pt();
      theHistograms->fill(Form("%sRes_dR"  , label), "photon resolution #DeltaR;#DeltaR"        , 20, 0.,0.4, deltaR       , theWeight);
      theHistograms->fill(Form("%sRes_E"   , label), "photon resolution Energy;#DeltaE/E"       , 20,-1.,1. , deltaEoverE  , theWeight);
      theHistograms->fill(Form("%sRes_pt"  , label), "photon resolution p_{T};#Deltap_{T}/p_{T}", 20,-1.,1. , deltapToverpT, theWeight);
      theHistograms->fill(Form("%sRes_EvsE", label), "photon resolution;E;#DeltaE/E", 20,0.,400., 16,-0.4,0.4, gen->e(), deltaEoverE , theWeight);
    }
  }
  
  // Old method for check
  // vector<Photon> photonsC(*kinPhotons_["central"]);
	
  // for(auto gPh : *genPhotons_){
  //   float ph_aeta = fabs(gPh.eta());
  //   if(ph_aeta > 1.4442 && ph_aeta < 1.566) continue;
		
  //   theHistograms->fill("effG_den_pt" , "p_{t} gen #gamma", 25,0.,250., gPh.pt() , theWeight);
  //   theHistograms->fill("effG_den_E"  , "E gen #gamma"    , 25,0.,250., gPh.e()  , theWeight);
  //   theHistograms->fill("effG_den_eta", "#eta gen #gamma" , eta_bins  , gPh.eta(), theWeight);
		
  //   if(photonsC.size() == 0 ) continue;
		
  //   std::sort   (photonsC.begin(), photonsC.end(), DeltaRComparator(gPh));
  //   std::reverse(photonsC.begin(), photonsC.end());
    
  //   Photon& rPh = photonsC.back();
  //   if(deltaR(rPh, gPh) < 0.2){
  //     theHistograms->fill("resG_dR", "#DeltaR(#gamma_{GEN}, #gamma_{REC})", 50,0.,0.1, deltaR(rPh,gPh), theWeight);
  //     theHistograms->fill("effG_num_pt" , "p_{t} gen #gamma", 25,0.,250., gPh.pt(),  theWeight);
  //     theHistograms->fill("effG_num_E"  , "E gen #gamma"    , 25,0.,250., gPh.e()  , theWeight);
  //     theHistograms->fill("effG_num_eta", "#eta gen #gamma" , eta_bins  , gPh.eta(), theWeight);
  //     photonsC.pop_back();  // Same rec photon cannot be matched to more than one gen ph
  //   }
  // }
}


void VVGammaAnalyzer::photonIsolation(const vector<Photon>& vPh, const char* label){
  // DeltaR(lep, ph): How many events would we loose if we excluded photons with DR(g, any l) < DR0 ?
  double maxG_minL_DR = -1;

  for(const Photon& ph : vPh){
    double minL_DR = 10;
    for(const Lepton& lep: *leptons_){
      double DR = physmath::deltaR(ph, lep);
      minL_DR = DR < minL_DR ? DR : minL_DR;
    }
    // const std::vector<Lepton>::iterator closestLep = std::min_element(leptons_->begin(), leptons_->end(), [](const Lepton& a, const Lepton& b){ return physmath::deltaR(ph, a) < physmath::deltaR(ph, b); } );
    theHistograms->fill(Form("minL_DR_%s", label), Form("min_{l}(#DeltaR(#gamma, l)) for each #gamma %s;#DeltaR;# of #gamma", label), 50, 0.,1., minL_DR, theWeight);
    maxG_minL_DR = minL_DR > maxG_minL_DR ? minL_DR : maxG_minL_DR;
  }
  if(maxG_minL_DR > 0)
    theHistograms->fill(Form("maxG_minL_DR_%s", label), Form("max_{#gamma}(min_{l}(#DeltaR(#gamma, l))) %s;#DeltaR;Events",label), 50, 0.,1., maxG_minL_DR, theWeight);
}


std::tuple<double, double, double> _nuEquationCoefficients(const TLorentzVector& pl, const TLorentzVector& pv){
  double El  = pl.E();
  ROOT::Math::XYVector plT(pl.Px(), pl.Py());  // Polar2DVector(pl.Pt(), )
  double plz = pl.Pz();
  double Ev  = pv.E();
  ROOT::Math::XYVector pvT(pv.Px(), pv.Py());
  //double pvz = pv.Pz();
  double pTdot = plT.Dot(pvT);

  double a = plz*plz - El*El;
  double b = phys::WMASS*phys::WMASS*plz + 2*plz*pTdot;
  double c = phys::WMASS*phys::WMASS*phys::WMASS*phys::WMASS/4 + pTdot*pTdot + phys::WMASS*phys::WMASS*pTdot - El*El*pvT.Mag2();
  
  return std::tuple<double, double, double>(a, b, c);
}


std::pair<TLorentzVector, TLorentzVector> solveNuPz(const Boson<Lepton>& W, int& error){
  std::tuple<double, double, double> coeff = _nuEquationCoefficients(W.daughter(0).p4(), W.daughter(1).p4());
  double a(std::get<0>(coeff)), b(std::get<1>(coeff)), c(std::get<2>(coeff));
  
  double delta = b*b - 4*a*c;
  if(delta < 0){
    error = 1;
    delta = 0;
  }
  double pz1 = (-b + sqrt(delta)) / (2*a);
  double pz2 = (-b - sqrt(delta)) / (2*a);
  
  TLorentzVector sol1(W.daughter(1).p4()), sol2(W.daughter(1).p4());
  sol1.SetPz(pz1);
  sol2.SetPz(pz2);
  sol1.SetE(sol1.Vect().Mag());
  sol2.SetE(sol2.Vect().Mag());
  
  return std::make_pair( sol1, sol2 );
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


bool isVeryLooseSiEiE(const Photon& ph, const double& barrel_thr=0.012, const double& endcap_thr=0.034){  // Additional separation between ph passing loose and fakes
  return ph.sigmaIetaIeta() > ( ph.isBarrel() ? barrel_thr : endcap_thr );
}


char VVGammaAnalyzer::phABCD(const Photon& ph, const Photon::IDwp wp){
  char ABCD = '0';  // Defaults to something recognizable
  if( ph.passChargedIsolation(wp) ){    // Signal region: either A or B
    if     (ph.passSigmaiEtaiEta(wp)) ABCD = 'A';
    else if(isVeryLooseSiEiE(ph)    ) ABCD = 'B';
  }
  else{                                 // Measurement region
    if     (ph.passSigmaiEtaiEta(wp)) ABCD = 'C';
    else if(isVeryLooseSiEiE(ph)    ) ABCD = 'D';
  }
  return ABCD;
}


char phABCD_study(const phys::Photon&, const double& barrel_thr, const double& endcap_thr){
  return '0';
}


void VVGammaAnalyzer::doPlots(const char* syst, const double& weight, const Photon* ph){
  if(is4Lregion(region_)){
    if(ph){
      theHistograms->fill(Form("SYS_mZZG_%s", syst), Form("m_{ZZ#gamma} %s", syst), mVVG_bins, (ZZ->p4() + ph->p4()).M(), weight);
      if(ph->cutBasedIDMedium())
      	theHistograms->fill(Form("SYS_mZZGmedium_%s", syst), Form("m_{ZZ#gamma} %s", syst), mVVG_bins, (ZZ->p4() + ph->p4()).M(), weight);
      if(ph->cutBasedIDTight())
      	theHistograms->fill(Form("SYS_mZZGtight_%s" , syst), Form("m_{ZZ#gamma} %s", syst), mVVG_bins, (ZZ->p4() + ph->p4()).M(), weight);
    }
    theHistograms->fill(  Form("SYS_mZZ_%s" , syst), Form("m_{ZZ} %s"      , syst), mVV_bins , ZZ->mass()               , weight);
  }
  else if(is3Lregion(region_)){
    if(ph){
      theHistograms->fill(Form("SYS_mWZG_%s", syst), Form("m_{WZ#gamma} %s", syst), mVVG_bins, (ZW->p4() + ph->p4()).M(), weight);
      if(ph->cutBasedIDMedium())
      	theHistograms->fill(Form("SYS_mWZGmedium_%s", syst), Form("m_{ZZ#gamma} %s", syst), mVVG_bins, (ZW->p4() + ph->p4()).M(), weight);
      if(ph->cutBasedIDTight())
      	theHistograms->fill(Form("SYS_mWZGtight_%s" , syst), Form("m_{ZZ#gamma} %s", syst), mVVG_bins, (ZW->p4() + ph->p4()).M(), weight);
    }
    theHistograms->fill(  Form("SYS_mWZ_%s" , syst), Form("m_{WZ} %s"      , syst), mVV_bins , ZW->mass()               , weight);
  }
}

void VVGammaAnalyzer::systematicsStudy(){
  double base_w = theWeight; // * ph.efficiencySF_
  
  const Photon* ph = nullptr;
  if(goodPhotons_["central"]->size() >= 1)
    ph = & (goodPhotons_["central"]->front());

  // Photon::IDwp wp = Photon::IDwp::Loose;
  // if( ph && !(ph->passSigmaiEtaiEta(wp) * ph->passHoverE(wp) * ph->passChargedIsolation(wp) * ph->passNeutralIsolation(wp) * ph->passPhotonIsolation(wp))){
  //   double mVVG = 0.;
  //   if     (is4Lregion(region_)) mVVG = (ZZ->p4() + ph->p4()).M();
  //   else if(is3Lregion(region_)) mVVG = (ZW->p4() + ph->p4()).M();
  //   cout << ">>> "<< evtN_ << " pt:"<< ph->pt() << " eta:" << ph->eta() << " IDloose?"<< ph->cutBasedIDLoose() << "   cuts: ";
  //   cout << ph->passSigmaiEtaiEta(wp)<<' '<<ph->passHoverE(wp)<<' '<<ph->passChargedIsolation(wp)<<' '<<ph->passNeutralIsolation(wp)<<' '<<ph->passPhotonIsolation(wp)<<'\n';
  //   cout << "\tchIso:" << ph->chargedIsolation() << " threshold:1.694" << "    phIso:" << ph->photonIsolation() << " threshold:"<<2.876 + 0.004017 *ph->pt() << '\n';
  //   cout << "\tnGoodPh:" << goodPhotons_["central"]->size() << "   mVVG:" << mVVG << '\n';
  //   cout << "--------------------------------------------------------------------------------\n";
  // }
  
  // central
  doPlots("central", base_w, ph);
  
  // Photons energy scale and resolution
  for(const auto& [syst, phVect] : goodPhotons_){
    if(strcmp(syst, "central") == 0) continue;
    if(phVect->size() == 0) continue;
    doPlots(Form("Ph%s", syst), base_w, & phVect->front());
  }
  
  bool isMC = theSampleInfo.isMC();
  // puWeightUnc
  doPlots("puWeight_Up"  , base_w * ( isMC ? theSampleInfo.puWeightUncUp() / theSampleInfo.puWeight() : 1.), ph);
  doPlots("puWeight_Down", base_w * ( isMC ? theSampleInfo.puWeightUncDn() / theSampleInfo.puWeight() : 1.), ph);
  
  // L1PrefiringWeight
  doPlots("L1Prefiring_Up"  , base_w * ( isMC ? theSampleInfo.L1PrefiringWeightUp() / theSampleInfo.L1PrefiringWeight() : 1.), ph);
  doPlots("L1Prefiring_Down", base_w * ( isMC ? theSampleInfo.L1PrefiringWeightDn() / theSampleInfo.L1PrefiringWeight() : 1.), ph);
  
  // QCD scale
  doPlots("QCDscaleF_Up"    , base_w * ( isMC ? theSampleInfo.QCDscale_muR1F2()   : 1.), ph);  // "QCDscale_muR1F2"
  doPlots("QCDscaleF_Down"  , base_w * ( isMC ? theSampleInfo.QCDscale_muR1F0p5() : 1.), ph);  // "QCDscale_muR1F0p5"
  doPlots("QCDscalemuR_Up"  , base_w * ( isMC ? theSampleInfo.QCDscale_muR2F1()   : 1.), ph);  // "QCDscale_muR2F1"
  doPlots("QCDscalemuR_Down", base_w * ( isMC ? theSampleInfo.QCDscale_muR0p5F1() : 1.), ph);  // "QCDscale_muR0p5F1"
  
  // PDF var
  doPlots("PDFVar_Up"  , base_w * ( isMC ? theSampleInfo.PDFVar_Up()   : 1.), ph);
  doPlots("PDFVar_Down", base_w * ( isMC ? theSampleInfo.PDFVar_Down() : 1.), ph);
  
  // alphas MZ
  doPlots("alphas_Up"  , base_w * ( isMC ? theSampleInfo.alphas_MZ_Up()   : 1.), ph);
  doPlots("alphas_Down", base_w * ( isMC ? theSampleInfo.alphas_MZ_Down() : 1.), ph);
  
  // lepton efficiency SF
  doPlots("eleEffSF_Up"  , base_w * (1 + ZZ->eleEffSFUnc()/ZZ->efficiencySF()), ph);
  doPlots("eleEffSF_Down", base_w * (1 - ZZ->eleEffSFUnc()/ZZ->efficiencySF()), ph);
  doPlots("muoEffSF_Up"  , base_w * (1 + ZZ->muEffSFUnc() /ZZ->efficiencySF()), ph);
  doPlots("muoEffSF_Down", base_w * (1 - ZZ->muEffSFUnc() /ZZ->efficiencySF()), ph);

  // lepton fake rate SF
  doPlots("eleFakeRateSF_Up"  , base_w * (1 + ZZ->eleFakeRateSFUnc()/ZZ->fakeRateSF()), ph);
  doPlots("eleFakeRateSF_Down", base_w * (1 - ZZ->eleFakeRateSFUnc()/ZZ->fakeRateSF()), ph);
  doPlots("muoFakeRateSF_Up"  , base_w * (1 + ZZ->muoFakeRateSFUnc()/ZZ->fakeRateSF()), ph);
  doPlots("muoFakeRateSF_Down", base_w * (1 - ZZ->muoFakeRateSFUnc()/ZZ->fakeRateSF()), ph);
  
  if(ph){
    // Photons ID efficiency
    doPlots("phEffSF_Up"  , base_w * (1 + ph->efficiencySFUnc()/ph->efficiencySF()), ph);
    doPlots("phEffSF_Down", base_w * (1 - ph->efficiencySFUnc()/ph->efficiencySF()), ph);
    // TODO: photon FR uncertaintiy
  }
  
  // doPlots("_Up"  , base_w * theSampleInfo.() / theSampleInfo.(), ph);
  // doPlots("_Down", base_w * theSampleInfo.() / theSampleInfo.(), ph);
}


void VVGammaAnalyzer::debug3Lregion(){
  const Lepton& l1 = ZW->first ().daughter(0);
  const Lepton& l2 = ZW->first ().daughter(1);
  const Lepton& l3 = ZW->second().daughter(0);
  theHistograms->fill("debug3L_l1_FRSF", "fake rate scale factor 1" , 101,-1.01,1.01, l1.fakeRateSF(), 1.);
  theHistograms->fill("debug3L_l2_FRSF", "fake rate scale factor 2" , 101,-1.01,1.01, l2.fakeRateSF(), 1.);
  theHistograms->fill("debug3L_l3_FRSF", "fake rate scale factor 3" , 101,-1.01,1.01, l3.fakeRateSF(), 1.);
  theHistograms->fill("debug3L_ZW_FRSF", "fake rate scale factor ZW", 101,-1.01,1.01,ZW->fakeRateSF(), 1.);

  double base_w = theSampleInfo.weight() * ZW->efficiencySF();
  double mWZ = ZW->mass();
  theHistograms->fill("debug3L_ZW_mass"    , "ZW mass, base weight;m_{ZW} [GeV/c^{2}]"                , mVV_bins, mWZ, base_w * l1.fakeRateSF());
  
  theHistograms->fill("debug3L_ZW_mass_w1" , "ZW mass, weighted only by FR_{l1};m_{ZW} [GeV/c^{2}]"   , mVV_bins, mWZ, base_w * l1.fakeRateSF());
  theHistograms->fill("debug3L_ZW_mass_w2" , "ZW mass, weighted only by FR_{l2};m_{ZW} [GeV/c^{2}]"   , mVV_bins, mWZ, base_w * l2.fakeRateSF());
  theHistograms->fill("debug3L_ZW_mass_w3" , "ZW mass, weighted only by FR_{l3};m_{ZW} [GeV/c^{2}]"   , mVV_bins, mWZ, base_w * l3.fakeRateSF());

  theHistograms->fill("debug3L_ZW_mass_w12", "ZW mass, weighted by FR_{l1}*FR_{l2};m_{ZW} [GeV/c^{2}]", mVV_bins, mWZ, base_w * l1.fakeRateSF()*l2.fakeRateSF());
  theHistograms->fill("debug3L_ZW_mass_w13", "ZW mass, weighted by FR_{l1}*FR_{l3};m_{ZW} [GeV/c^{2}]", mVV_bins, mWZ, base_w * l1.fakeRateSF()*l3.fakeRateSF());
  theHistograms->fill("debug3L_ZW_mass_w23", "ZW mass, weighted by FR_{l2}*FR_{l3};m_{ZW} [GeV/c^{2}]", mVV_bins, mWZ, base_w * l2.fakeRateSF()*l3.fakeRateSF());

  theHistograms->fill("debug3L_ZW_mass_w123", "ZW mass, weighted by FR_{l1}*FR_{l2}*FR_{l3};m_{ZW} [GeV/c^{2}]", mVV_bins, mWZ, base_w * l1.fakeRateSF()*l2.fakeRateSF()*l3.fakeRateSF());
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

const vector<double> VVGammaAnalyzer::ph_aeta_bins(
						   {0., 0.783, 1.4442, 1.566, 2, 2.5}
						   // {0., 0.435, 0.783, 1.13, 1.4442, 1.566, 1.8, 2.1, 2.5}
						   );

const vector<double> VVGammaAnalyzer::ph_pt_bins(
						 {20., 45., 70., 120}
						 // {20., 30., 45., 70., 120}
						 // {20., 25., 30, 45., 70., 100}
						 );


const vector<double> VVGammaAnalyzer::mVV_bins(
					       {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000}
					       );

const vector<double> VVGammaAnalyzer::mVVG_bins(
						{0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000}
						//{150, 250, 350, 450, 850}
						);
double VVGammaAnalyzer::getPhotonFR   (const phys::Photon& ph) const{
  return hPhotonFR_->GetBinContent(hPhotonFR_->FindBin(ph.pt(), abs(ph.eta())));
}
double VVGammaAnalyzer::getPhotonFRUnc(const phys::Photon& ph) const{
  double FR      = hPhotonFR_->GetBinContent(hPhotonFR_->FindBin(ph.pt(), abs(ph.eta())));
  double FRError = hPhotonFR_->GetBinError  (hPhotonFR_->FindBin(ph.pt(), abs(ph.eta())));
  return FR == 0 ? 0. : FRError/FR ;
}
