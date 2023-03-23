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

#define BIN_GENCATEGORY 10,-0.5,9.5
#define BINS_CUTFLOW 9,-0.5,8.5
//{"All", "ZZ || ZW", "2l2j || 2l1J", "#gamma kin", "#gamma good", "Analyzed", "#gamma medium"}
#define BINS_PHCUTFLOW 7,-0.5,6.5
#define BINS_PHCUTNM1 6,-0.5,5.5
#define BINS_KINPHID 5,-0.5,4.5

std::pair<TLorentzVector, TLorentzVector> solveNuPz(const Boson<Lepton>& W, int& error);

void VVGammaAnalyzer::begin(){
  cout<<'\n';
  for(char i=0; i<25; ++i) cout<<'-';
  cout<<" Start of VVGammaAnalyzer ";
  for(char i=0; i<25; ++i) cout<<'-';
  cout<<'\n';
  
  // Photon FR
  TFile fileFR(Form("data/LtoT_FR_data-ZGToLLG_%d.root", year), "READ" );
  if(fileFR.IsOpen()){
    hPhotonFR_.reset( std::move((TH2F*) fileFR.Get("PhFR")) );
    hPhotonFR_->SetDirectory(nullptr);  // prevent ROOT from deleting it
    cout << "INFO: retrieved Photon FR histogram from \""<<fileFR.GetName()<<"\"\n";
    fileFR.Close();
  }
  else{
    cout << colour::Red("WARN") << ": Photon FR file not found in \""<<fileFR.GetName()<<"\"\n";
    hPhotonFR_.reset( new TH2F("PhFR", "", 1,0.,1., 1,0.,1.) );
  }
  
  // FR extended
  TFile fileFR_KtoVL(Form("data/KtoVL_FR_data-ZGToLLG_%d.root", year), "READ" );
  if(fileFR_KtoVL.IsOpen()){
    hPhotonFR_KtoVL_.reset( std::move((TH2F*) fileFR_KtoVL.Get("PhFR")) );
    hPhotonFR_KtoVL_->SetDirectory(nullptr);  // prevent ROOT from deleting it
    cout << "INFO: retrieved Photon FR KtoVL histogram from \""<<fileFR_KtoVL.GetName()<<"\"\n";
    fileFR_KtoVL.Close();
  }
  else{
    cout << colour::Red("WARN") << ": Photon FR KtoVL file not found in \""<<fileFR_KtoVL.GetName()<<"\"\n";
    hPhotonFR_KtoVL_.reset( new TH2F("PhFR", "", 1,0.,1., 1,0.,1.) );
  }
  
  TFile fileFR_KtoVLexcl(Form("data/KtoVLexcl_FR_data-ZGToLLG_%d.root", year), "READ" );
  if(fileFR_KtoVLexcl.IsOpen()){
    hPhotonFR_KtoVLexcl_.reset( std::move((TH2F*) fileFR_KtoVLexcl.Get("PhFR")) );
    hPhotonFR_KtoVLexcl_->SetDirectory(nullptr);  // prevent ROOT from deleting it
    cout << "INFO: retrieved Photon FR KtoVLexcl histogram from \""<<fileFR_KtoVLexcl.GetName()<<"\"\n";
    fileFR_KtoVLexcl.Close();
  }
  else{
    cout << colour::Red("WARN") << ": Photon FR KtoVLexcl file not found in \""<<fileFR_KtoVLexcl.GetName()<<"\"\n";
    hPhotonFR_KtoVLexcl_.reset( new TH2F("PhFR", "", 1,0.,1., 1,0.,1.) );
  }
  
  
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
  cout<<"\tSizes: Particle="<<sizeof(Particle)<<" , Lepton="<<sizeof(Lepton)<<" , Photon="<<sizeof(Photon)<<" Jet="<<sizeof(Jet)<<" , Boson<Particle>="<<sizeof(Boson<Particle>)<<" , DiBoson<Particle,Particle>="<<sizeof(DiBoson<Particle,Particle>)<<'\n';
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
  candVToJ_ = Jet();

  // Contruct vector with leptons from dibosons
  if     (is4Lregion(region_) || (region_ == MC && ZZ && ZZ->pt() > 1.))
    leptons_->insert(leptons_->end(), {
    	  ZZ->first().daughter(0), 
    	  ZZ->first().daughter(1), 
    	  ZZ->second().daughter(0), 
    	  ZZ->second().daughter(1)
        });
  else if(is3Lregion(region_) || (region_ == MC && ZW && ZW->pt() > 1.))
    leptons_->insert(leptons_->end(), {
  	  ZW->first().daughter(0), 
  	  ZW->first().daughter(1), 
  	  ZW->second().daughter(0)
	});
  else if((region_ == CRLFR || region_ == MC) && ZL && ZL->first.pt() > 1.)
    leptons_->insert(leptons_->end(), {
  	  ZL->first.daughter(0),
  	  ZL->first.daughter(1),
  	  ZL->second
	});
  else if(is2Lregion(region_) || region_ == MC){
    leptons_->insert(leptons_->cend(), electrons->begin(), electrons->end());
    leptons_->insert(leptons_->cend(),     muons->begin(),     muons->end());
  }

  
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
      if(deltaR(ph,lep) < 0.07){  // 0.3
	match = true;
	break;
      }
    }
    if(match) continue;
    
    // Photon::IdWp wp = Photon::IdWp::Loose;
    TLorentzVector p4_EScale_Up = ph.p4() * (ph.energyScaleUp()  /ph.e());
    TLorentzVector p4_EScale_Dn = ph.p4() * (ph.energyScaleDown()/ph.e());
    TLorentzVector p4_ESigma_Up = ph.p4() * (ph.energySigmaUp()  /ph.e());
    TLorentzVector p4_ESigma_Dn = ph.p4() * (ph.energySigmaDown()/ph.e());
    
    if(p4_EScale_Up.Pt() > 20){
      Photon copy(ph);
      copy.setP4(p4_EScale_Up);
      kinPhotons_["EScale_Up"]->push_back(copy);
      if(passVeryLoose(ph))
	loosePhotons_["EScale_Up"]->push_back(copy);
      if(copy.cutBasedIDLoose())
	goodPhotons_["EScale_Up"]->push_back(std::move(copy));
    }
    if(p4_EScale_Dn.Pt() > 20){
      Photon copy(ph);
      copy.setP4(p4_EScale_Dn);
      kinPhotons_["EScale_Down"]->push_back(copy);
      if(passVeryLoose(ph))
    	loosePhotons_["EScale_Down"]->push_back(copy);
      if(copy.cutBasedIDLoose())
    	goodPhotons_["EScale_Down"]->push_back(std::move(copy));
    }
    if(p4_ESigma_Up.Pt() > 20){
      Photon copy(ph);
      copy.setP4(p4_ESigma_Up);
      kinPhotons_["ESigma_Up"]->push_back(copy);
      if(passVeryLoose(ph))
        loosePhotons_["ESigma_Up"]->push_back(copy);
      if(copy.cutBasedIDLoose())
    	goodPhotons_["ESigma_Up"]->push_back(std::move(copy));
    }
    if(p4_ESigma_Dn.Pt() > 20){
      Photon copy(ph);
      copy.setP4(p4_ESigma_Dn);
      kinPhotons_["ESigma_Down"]->push_back(copy);
      if(passVeryLoose(ph))
    	loosePhotons_["ESigma_Down"]->push_back(copy);
      if(copy.cutBasedIDLoose())
    	goodPhotons_["ESigma_Down"]->push_back(std::move(copy));
    }
    if(ph.pt() > 20){
      kinPhotons_["central"]->push_back(ph);
      if(passVeryLoose(ph))
	loosePhotons_["central"]->push_back(ph);
      if(ph.cutBasedIDLoose())
	goodPhotons_["central"]->push_back(ph);
    }
    
  } // end loop on *photons
  
  if(kinPhotons_["central"]->size() > 0)
    bestKinPh_ = &*std::max_element(kinPhotons_["central"]->begin(), kinPhotons_["central"]->end(),
				    [](const Photon& a, const Photon& b){ return a.nCutsPass(Photon::IdWp::Loose) < b.nCutsPass(Photon::IdWp::Loose); }
				    );  // max_element returns the first among those with max value --> preserve pt ordering
  else
    bestKinPh_ = nullptr;
  
}


Int_t VVGammaAnalyzer::cut() {
  evtN_++; evtNInReg_[region_]++; evtWInReg_[region_] += theWeight;
  #ifndef DEBUG
  // cout<<"\r\t\t"<<evtN_;
  #else
  cout << regionType(region_) << ':' << run << ':' << lumiBlock << ':' << event << '\t'
       << theWeight << '\n';
  // for(const Particle& p : *genParticles)
  //   cout << "\tid: " << p.id() << " \tgenStatusFlag = " << p.genStatusFlags() << " --> " << p.genStatusFlags().to_ulong() << std::endl;
  #endif
       // << ZZ->mass() << ' ' << ZW->p4().Mt() << ' ' << kinPhotons_["central"]->size() << ' ' << goodPhotons_["central"]->size() << '\n';
  // if(cherrypickEvt()){
    // theHistograms->fill("cherry_ZZ_mass"    , "Events in {UL#backslashLegacy}: m_{ZZ};m_{ZZ} [GeV/c^{2}]", 25,0.,500., ZZ->mass()                      , theWeight);
    // theHistograms->fill("cherry_ZZ_goodLept", "Events in {UL#backslashLegacy}: # good leptons"           , 5,-0.5,4.5, ZZ->numberOfGoodGrandDaughters(), theWeight);
    // theHistograms->fill("cherry_ZZ_badLept" , "Events in {UL#backslashLegacy}: # bad leptons"            , 5,-0.5,4.5, ZZ->numberOfBadGrandDaughters() , theWeight);
  // }
  
  std::bitset<32> genCategory_bits(genCategory);
  theHistograms->fill("AAA_genCategory"  , "gen category weighted"  , BIN_GENCATEGORY, 0, theWeight);
  theHistograms->fill("AAA_genCategory_u", "gen category unweighted", BIN_GENCATEGORY, 0);
  for(int i = 0; i < 8; ++i){
    if(genCategory_bits.test(i)){
      theHistograms->fill("AAA_genCategory"  , "gen category weighted"  , BIN_GENCATEGORY, i+1, theWeight);
      theHistograms->fill("AAA_genCategory_u", "gen category unweighted", BIN_GENCATEGORY, i+1);
    }
  }
  if((genCategory_bits.test(0) || genCategory_bits.test(1)) && genCategory_bits.test(4) && genCategory_bits.test(5) && genCategory_bits.test(6) && genCategory_bits.test(7)){
     theHistograms->fill("AAA_genCategory"  , "gen category weighted"  , BIN_GENCATEGORY, 9, theWeight);
     theHistograms->fill("AAA_genCategory_u", "gen category unweighted", BIN_GENCATEGORY, 9);
  }

  theHistograms->fill("AAA_cuts"  , "cuts weighted"  , BINS_CUTFLOW, 0, theWeight);
  theHistograms->fill("AAA_cuts_u", "cuts unweighted", BINS_CUTFLOW, 0, 1);

  bool haveZVlep = false;
  bool have2l2j = false;
  bool haveGoodPhoton = false;
  
  if(theSampleInfo.isMC()){
    genEventSetup();
    genEventHistos();
  }
  initEvent();
  
  theHistograms->fill("POG_leptons", "n leptons", 7,-0.5,6.5, electrons->size()+muons->size(), theWeight);  
  
  baseHistos_cut();
  photonHistos();
  jetHistos();
  // PKU_comparison();
  photonIsolation(*    photons            , "all" );
  photonIsolation(* kinPhotons_["central"], "kin" );
  photonIsolation(*goodPhotons_["central"], "good");

  efficiency(*    photons            , *genPhotons_, "photons"    , "gen", 0.2);
  efficiency(* kinPhotons_["central"], *genPhotons_, "kinPhotons" , "gen", 0.2);
  efficiency(*goodPhotons_["central"], *genPhotons_, "goodPhotons", "gen", 0.2);
  
  efficiency(*jets, *genJets, "AK4", "genJets", 0.4);
  
  
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
  
  if(is3Lregion(region_)){
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
    bool b_lll_mass = (lZ1.p4() + lZ2.p4() + lW.p4()).M() > 100;
    bool b_bveto = ! std::any_of(jets->begin(), jets->end(), [](const Jet& j){ auto dF = j.deepFlavour(); return dF.probb + dF.probbb + dF.problepb > 0.2770; });    // Note: this is the medium WP for Legacy samples (102X)
    
    bool b_WZpaperSel = b_Z_mass && b_MET && b_lll_mass && b_lepton_pt && b_bveto;
    
    if(b_lepton_pt) theHistograms->fill("WZ_cutflow", "Cutflow for WZ", 7,-0.5,6.5, 1, theWeight);
    if(b_Z_mass   ) theHistograms->fill("WZ_cutflow", "Cutflow for WZ", 7,-0.5,6.5, 2, theWeight);
    if(b_MET      ) theHistograms->fill("WZ_cutflow", "Cutflow for WZ", 7,-0.5,6.5, 3, theWeight);
    if(b_lll_mass ) theHistograms->fill("WZ_cutflow", "Cutflow for WZ", 7,-0.5,6.5, 4, theWeight);
    if(b_bveto    ) theHistograms->fill("WZ_cutflow", "Cutflow for WZ", 7,-0.5,6.5, 5, theWeight);
    if(b_lepton_pt && b_Z_mass && b_MET && b_lll_mass && b_bveto)
      theHistograms->fill("WZ_cutflow", "Cutflow for WZ", 7,-0.5,6.5, 6, theWeight);
    
    if(!b_WZpaperSel)
      return -1;
    else{
      theHistograms->fill("AAA_cuts"  , "cuts weighted"  , BINS_CUTFLOW, 6, theWeight);
      theHistograms->fill("AAA_cuts_u", "cuts unweighted", BINS_CUTFLOW, 6, 1);
    }
  }
  
  
  if(!haveZVlep) 
    return -1;
  return 1;  //TEMP
  
  if(!haveGoodPhoton)
    return -1;
  else
    return 1;
}


void VVGammaAnalyzer::analyze(){
  analyzedNInReg_[region_]++; analyzedWInReg_[region_] += theWeight;
  theHistograms->fill("AAA_cuts"  , "cuts weighted"  , BINS_CUTFLOW, 7, theWeight);
  theHistograms->fill("AAA_cuts_u", "cuts unweighted", BINS_CUTFLOW, 7, 1);
  
  if( std::any_of(goodPhotons_["central"]->cbegin(),
		  goodPhotons_["central"]->cend(), 
		  [](const Photon& ph){return ph.cutBasedIDMedium();}) )
  {
    theHistograms->fill("AAA_cuts"  , "cuts weighted"  , BINS_CUTFLOW, 8, theWeight);
    theHistograms->fill("AAA_cuts_u", "cuts unweighted", BINS_CUTFLOW, 8, 1);
  }
  
  bool four_lep  = is4Lregion(region_);  // && ZZ && ZZ->pt() > 1.;
  bool three_lep = is3Lregion(region_);  // && ZW && ZW->pt() > 1.;
  bool two_lep   = is2Lregion(region_);
  bool LFR_lep   = region_ == CRLFR;     // && ZL && ZL->first.pt() > 1.;
  
  if(two_lep){
    hadronicObjectsReconstruction();
    studyJetsChoice();
  }

  if(three_lep)
    debug3Lregion();
  
  // Decide channel name depending on leptons
  makeChannelReco();
  
  // leptonFakeRate();
  photonGenStudy();
  photonFakeRate();
  photonFakeRate_extended();
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
    
    plotsVVGstatus("ZZ", "ZZ", ZZ->p4(), "mass");
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
    
    plotsVVGstatus("ZW", "ZW", ZW->p4(), "massT");
    
    // Check triggers in signal region
    if(goodPhotons_["central"]->size() > 0){
      std::bitset<16> triggerBits(triggerWord);
      // bool passOneTrigger = triggerBits.test(0);  // triggerWord << 0 & 0x1
      bool passDiMu       = triggerBits.test(1);  // triggerWord << 1 & 0x1
      bool passDiEle      = triggerBits.test(2);  // triggerWord << 2 & 0x1
      bool passMuEle      = triggerBits.test(3);  // triggerWord << 3 & 0x1
      bool passTriEle     = triggerBits.test(4);  // triggerWord << 4 & 0x1
      bool passTriMu      = triggerBits.test(5);  // triggerWord << 5 & 0x1
      bool passSingleEle  = triggerBits.test(6);  // triggerWord << 6 & 0x1
      bool passSingleMu   = triggerBits.test(7);  // triggerWord << 7 & 0x1
      
      // triggers in common: 0000 0000 1100 1111 = 0xCF
      // extra triggers    : 0000 0000 0011 0001 = 0x31
      bool commonTriggers = passSingleMu || passSingleEle || passDiMu || passMuEle || passDiEle;  // triggerWord & 0xCF
      bool onlyExtraTriggers = (passTriEle || passTriMu) && (! commonTriggers);
      if(true)
	theHistograms->fill("trig_Cutflow"  , "Trigger study"  , 4,0,4, 0, theWeight);
      if(commonTriggers)
	theHistograms->fill("trig_Cutflow"  , "Trigger study"  , 4,0,4, 1, theWeight);
      if(passTriEle || passTriMu)
	theHistograms->fill("trig_Cutflow"  , "Trigger study"  , 4,0,4, 2, theWeight);
      if(onlyExtraTriggers)
	theHistograms->fill("trig_Cutflow"  , "Trigger study"  , 4,0,4, 3, theWeight);  
    }
    
    // Neutrino pz reconstruction
    // int error;
    // const Boson<Lepton>& Wboson = ZW->second();
    // std::pair<TLorentzVector, TLorentzVector> solutions = solveNuPz(Wboson, error);
    // if(!error){
    //   const TLorentzVector& pl = Wboson.daughter(0).p4();
    //   const TLorentzVector& pv = (gRandom->Integer(2) ? solutions.first : solutions.second);
    //   theHistograms->fill("ZW_mass_reconstructed", "m_{3l}, p_{z}^{#nu} from m_{W} contraint; GeV/c^2", mVV_bins, (pl + pv + ZW->first().p4()).M(), theWeight);
    // }
    
    for(size_t i = 0; i < leptons_->size(); ++i){
      for(size_t j = i+1; j < leptons_->size(); ++j){
	double mij = (leptons_->at(i).p4() + leptons_->at(j).p4()).M();
	if(mij < 4){
	  theHistograms->fill("ERROR_smartCut", "Smart cut not applied;m_{ll} [GeV/c^2];# pairs", 10,0.,10, mij);
	}
      }
    }

    // Is the third lepton the problem?
    const Lepton& lZ1 = ZW->first ().daughter(0);
    const Lepton& lZ2 = ZW->first ().daughter(1);
    const Lepton& lW  = ZW->second().daughter(0);
    
    const char* leptFlav = ( abs(lW.id()) == 11 ? "e" : "m"  );
    theHistograms->fill(Form("l3_%s_pt" , leptFlav), "3rd lepton p_{t};p_{t} [GeV/c]"                , 50,0.,200., lW.pt()                 , theWeight);
    theHistograms->fill(Form("l3_%s_Iso", leptFlav), "3rd lepton combRelIsoFSRCorr;combRelIsoFSRCorr", 50,0.,0.1 , lW.pfCombRelIsoFSRCorr(), theWeight);
    // if(b_MET){
    //   theHistograms->fill(Form("l3_%s_pt_MET" ,leptFlav),"3rd lepton p_{t} - MET > 30;p_{t} [GeV/c]"                ,50,0.,200.,lW.pt()                 ,theWeight);
    //   theHistograms->fill(Form("l3_%s_Iso_MET",leptFlav),"3rd lepton combRelIsoFSRCorr - MET > 30;combRelIsoFSRCorr",50,0.,0.1 ,lW.pfCombRelIsoFSRCorr(),theWeight);
    // }
    
    // if(b_WZpaperSel){
    //   theHistograms->fill("paperSel_ZW_massT"              , "m_{T,3l};GeV/c^{2}", mVV_bins   , ZW->p4().Mt()                , theWeight);
    //   theHistograms->fill("paperSel_Z_mass"                , "m_{Z};GeV/c^{2}"   , 35,55.,125., ZW->first().mass()           , theWeight);
    //   theHistograms->fill("paperSel_W_massT"               , "m_{T,W};GeV/c^{2}" , 35,55.,125., ZW->second().p4().Mt()       , theWeight);
    //   theHistograms->fill("paperSel_ZW_pt"                 , "p_{t,ZW};GeV/c"    , 20,0.,400. , ZW->pt()                     , theWeight);
    //   theHistograms->fill("paperSel_Z_l0_pt"               , "p_{t,l00};GeV/c"   , 20,0.,400. , ZW->first().daughter(0).pt() , theWeight);
    //   theHistograms->fill("paperSel_Z_l1_pt"               , "p_{t,l01};GeV/c"   , 20,0.,400. , ZW->first().daughter(1).pt() , theWeight);
    //   theHistograms->fill("paperSel_W_l_pt"                , "p_{t,l10};GeV/c"   , 20,0.,400. , ZW->second().daughter(0).pt(), theWeight);
    //   theHistograms->fill("paperSel_W_MET_pt"              , "p_{t,MET};GeV/c"   , 20,0.,400. , ZW->second().daughter(1).pt(), theWeight);
    //   theHistograms->fill("paperSel_lll_mass"              , "m_{lll};GeV/c^{2}" , 20,0.,400. , m_lll                        , theWeight);

    //   theHistograms->fill("paperSel_ZW_massT_"+channelReco_, "m_{T,3l};GeV/c^{2}", mVV_bins   , ZW->p4().Mt()                , theWeight);
    //   plotsVVGstatus("paperSel_ZW", "paper selection ZW", ZW->p4(), "massT");
    // }
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
      
      theHistograms->fill("Z_mass_"+channelReco_, "m_{ll};GeV/c^{2}", 30, 60., 120., Zcand.mass(), theWeight);
      
      // Vhad histograms
      if(candVTojj_.isValid()){
	theHistograms->fill("VTojj_mass", "V#rightarrowjj mass;GeV/c^{2}", 30, 50, 125, candVTojj_.mass(), theWeight);
	theHistograms->fill("VTojj_pt"  , "V#rightarrowjj p_{t};GeV/c"   , 25,  0, 500, candVTojj_.pt()  , theWeight);
	theHistograms->fill("VTojj_eta" , "V#rightarrowjj #eta;#eta"     , 25, -5,   5, candVTojj_.eta() , theWeight);
	
	TLorentzVector pZjj = Zcand.p4() + candVTojj_.p4();
	plotsVVGstatus("Zjj", "Zjj", pZjj, "mass");
      }
      if(candVToJ_.isValid()){
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
	
	TLorentzVector pZJ = Zcand.p4() + candVToJ_.p4();
	plotsVVGstatus("ZJ", "ZJ", pZJ, "mass");
      }
      
      plotsVVGstatus("Z", "Z", Zcand.p4(), "mass");
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
  for(const char* name : {"AAA_cuts", "AAA_cuts_u", "kinPhotons_ID", "kinPh_eScale_N", "kinPhotons_cuts"}){
    TH1* h = theHistograms->get(name);
    if(h) h->LabelsDeflate("X");
  }
  
  TH1* AAA_genCategory   = theHistograms->get("AAA_genCategory"  );
  TH1* AAA_genCategory_u = theHistograms->get("AAA_genCategory_u");
  for(TH1* h : {AAA_genCategory, AAA_genCategory_u}){
    if(!h) continue;
    TAxis* axis = h->GetXaxis();
    axis->SetBinLabel(1, "All");
    axis->SetBinLabel(2, "ZZ");
    axis->SetBinLabel(3, "ZW");
    axis->SetBinLabel(4, "unused");
    axis->SetBinLabel(5, "Z");
    axis->SetBinLabel(6, "Fid accept");
    axis->SetBinLabel(7, "Detect accept");
    axis->SetBinLabel(8, "Trig plateau");
    axis->SetBinLabel(9, "mVV > 100");
    axis->SetBinLabel(10, "0|1 & 4 & 5 & 6 & 7");
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
    axis->SetBinLabel(5, "#gamma VeryLoose");
    axis->SetBinLabel(6, "#gamma Loose (Good)");
    axis->SetBinLabel(7, "WZ_paperSel");
    axis->SetBinLabel(8, "Analyzed");
    axis->SetBinLabel(9, "#gamma Medium");
  }
  
  TH1* trig_Cutflow = theHistograms->get("trig_Cutflow");
  if(trig_Cutflow){
    TAxis* axis = trig_Cutflow->GetXaxis();
    axis->SetBinLabel(1, "all");
    axis->SetBinLabel(2, "commonTrig");
    axis->SetBinLabel(3, "triEle|triMu");
    axis->SetBinLabel(4, "extraTrig");
  }
  
  TH1* kinPhotons_cuts = theHistograms->get("kinPhotons_cuts");
  if(kinPhotons_cuts){
    TAxis* axis = kinPhotons_cuts->GetXaxis();
    axis->SetBinLabel(1, "No #gamma");
    axis->SetBinLabel(2, "All #gamma");
    axis->SetBinLabel(3, "#sigma_{i#etai#eta}");
    axis->SetBinLabel(4, "HoverE");
    axis->SetBinLabel(5, "IsoCH");
    axis->SetBinLabel(6, "IsoNE");
    axis->SetBinLabel(7, "IsoPh");
  }
  
  TH1* kinPhotons_Nm1  = theHistograms->get("kinPhotons_Nm1" );
  if(kinPhotons_Nm1){
    TAxis* axis = kinPhotons_Nm1->GetXaxis();
    axis->SetBinLabel(1, ">= 4");
    axis->SetBinLabel(2, "#sigma_{i#etai#eta}");
    axis->SetBinLabel(3, "HoverE");
    axis->SetBinLabel(4, "IsoCH");
    axis->SetBinLabel(5, "IsoNE");
    axis->SetBinLabel(6, "IsoPh");
  }
  
  TH1* kinPhotons_ID = theHistograms->get("kinPhotons_ID");
  if(kinPhotons_ID){
    TAxis* axis = kinPhotons_ID->GetXaxis();
    axis->SetBinLabel(1, "Kinematic");
    axis->SetBinLabel(2, "VeryLoose");
    axis->SetBinLabel(3, "Loose");
    axis->SetBinLabel(4, "Medium");
    axis->SetBinLabel(5, "Tight");
  }
  
  if(is3Lregion(region_)){
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
  
  TH1* GEN_photons_genStatusFlags = theHistograms->get("GEN_photons_genStatusFlags");
  if(GEN_photons_genStatusFlags){
    TAxis* axis = GEN_photons_genStatusFlags->GetXaxis();
    axis->SetBinLabel( 1, "isPrompt");
    axis->SetBinLabel( 2, "isDecayedLeptonHadron");
    axis->SetBinLabel( 3, "isTauDecayProduct");
    axis->SetBinLabel( 4, "isPromptTauDecayProduct");
    axis->SetBinLabel( 5, "isDirectTauDecayProduct");
    axis->SetBinLabel( 6, "isDirectPromptTauDecayProduct");
    axis->SetBinLabel( 7, "isDirectHadronDecayProduct");
    axis->SetBinLabel( 8, "isHardProcess");
    axis->SetBinLabel( 9, "fromHardProcess");
    axis->SetBinLabel(10, "isHardProcessTauDecayProduct");
    axis->SetBinLabel(11, "isDirectHardProcessTauDecayProduct");
    axis->SetBinLabel(12, "fromHardProcessBeforeFSR");
    axis->SetBinLabel(13, "isFirstCopy");
    axis->SetBinLabel(14, "isLastCopy");
    axis->SetBinLabel(15, "isLastCopyBeforeFSR");
    axis->SetBinLabel(16, "All");
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
      cout << colour::Red("WARN") << ": no cherry pick file for region " << regionType(R) << '\n';
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
  for  (auto i = vec.cbegin(); i != vec.end()-1; ++i)
    for(auto j = i+1         ; j != vec.end()  ; ++j){
      Boson<T> b(*i, *j, -23);  // last parameter is ID; if not set to something != 0 will cause isValid() to return false and me to waste an afternoon debugging
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
  vector<Boson<Jet>> candsjj = candidatesVTojj(*jets);  // all pairs with mW - 30 < m(jj) < mZ + 30
  
  if(candsjj.size() > 0){
    std::sort(candsjj.begin(), candsjj.end(), Mass2Comparator(phys::ZMASS, phys::WMASS));
    candVTojj_ = Boson<Jet>(std::move(candsjj.front()));
  }
  if(jetsAK8->size() > 0){
    auto it = std::min_element(jetsAK8->begin(), jetsAK8->end(), Mass2Comparator(phys::ZMASS, phys::WMASS));
    if(minDM(it->mass()) < 30)
      candVToJ_ = *it;
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
    unsigned int aPID = abs(p.id());
    if(aPID < 9)
      genQuarks_->push_back(p);
    else if(aPID == 11 || aPID == 13){
      genChLeptons_->push_back(p);
    }
    else if(aPID == 12 || aPID == 14)
      genNeutrinos_->push_back(p);
    else if(p.id() == 22){
      genPhotons_->push_back(p);
    }
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
  
  // Gen Z --> l lbar
  if(genChLeptons_->size() >= 2){
    for(size_t i = 0 ; i < genChLeptons_->size(); ++i){
      Particle& l1 = genChLeptons_->at(i);
      for(size_t j = i+1; j < genChLeptons_->size(); ++j){
	Particle& l2 = genChLeptons_->at(j);
	
	if( l1.id() + l2.id() == 0 ){
	  Boson<Particle> Zcand(l1,l2);
	  if(ZBosonDefinition(Zcand))
	    genZlepCandidates_->push_back(Zcand);
	}
      }
    }
  }
  
  if(genQuarks_->size() >= 2){
    for(size_t i = 0  ; i < genQuarks_->size(); ++i){
      Particle& q1 = genQuarks_->at(i);
      if(q1.id() > 5) continue;
      for(size_t j = i+1; j < genQuarks_->size(); ++j){
	Particle& q2 = genQuarks_->at(j);
	if(q2.id() > 5) continue;

	// Gen W --> q q'bar
	if( (q1.id() * q2.id() < 0) && ( abs(q1.id()+q2.id()) % 2 ==1 ) ){
	  Boson<Particle> Wcand(q1,q2);
	  if(GenWBosonDefinition(Wcand))
	    genWhadCandidates_->push_back(Wcand);
	}

	// Gen Z --> q qbar
	if( q1.id() + q2.id() == 0 ){
	  Boson<Particle> Zcand(q1,q2);
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
  theHistograms->fill("GEN_ZlepCandidates", "# genZlepCandidates_", 4,-0.5,3.5, genZlepCandidates_->size());
  theHistograms->fill("GEN_WlepCandidates", "# genWlepCandidates_", 4,-0.5,3.5, genWlepCandidates_->size());
  theHistograms->fill("GEN_ZhadCandidates", "# genZhadCandidates_", 4,-0.5,3.5, genZhadCandidates_->size());
  theHistograms->fill("GEN_WhadCandidates", "# genWhadCandidates_", 4,-0.5,3.5, genWhadCandidates_->size());
  
  theHistograms->fill("GEN_quarks"   , "# genQuarks"   , 10,-0.5,9.5, genQuarks_   ->size());
  theHistograms->fill("GEN_chLeptons", "# genChLeptons", 10,-0.5,9.5, genChLeptons_->size());
  theHistograms->fill("GEN_neutrinos", "# genNeutrinos", 10,-0.5,9.5, genNeutrinos_->size());
  theHistograms->fill("GEN_photons"  , "# genPhotons"  , 10,-0.5,9.5, genPhotons_  ->size());
  
  for(auto v : *genZlepCandidates_)
    theHistograms->fill("GEN_genZlepCandidates_mass", "mass genZlepCandidates;[GeV/c^{2}]", 35.,50.,120., v.mass());
  for(auto v : *genWlepCandidates_)
    theHistograms->fill("GEN_genWlepCandidates_mass", "mass genWlepCandidates;[GeV/c^{2}]", 35.,50.,120., v.mass());
  for(auto v : *genZhadCandidates_)
    theHistograms->fill("GEN_genZhadCandidates_mass", "mass genZhadCandidates;[GeV/c^{2}]", 35.,50.,120., v.mass());
  for(auto v : *genWhadCandidates_)
    theHistograms->fill("GEN_genWhadCandidates_mass", "mass genWhadCandidates;[GeV/c^{2}]", 35.,50.,120., v.mass());
  
  theHistograms->fill("GEN_n_jets"   , "Number of genJets;# genJets"   , 6,-0.5,5.5, genJets->size()   );
  theHistograms->fill("GEN_n_jetsAK8", "Number of genJetsAK8;# genJets", 6,-0.5,5.5, genJetsAK8->size());
  
  for(const Particle& p : *genPhotons_){
    theHistograms->fill("GEN_photons_genStatusFlags", "#gamma_{GEN} pass flag", 16, -0.5, 15.5, 15, theWeight);
    std::bitset<15> flags = p.genStatusFlags();
    theHistograms->book<TH1F>("GEN_photons_genStatus_ulong", "#gamma_{GEN} flag ulong", vector<double>({0.,1.}))->Fill(Form("%zu", flags.to_ulong()), theWeight);
    for(size_t i = 0; i < 15; ++i)
      if(flags.test(i))
	theHistograms->fill("GEN_photons_genStatusFlags", "#gamma_{GEN} pass flag", 16, -0.5, 15.5, i, theWeight);
  }
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
  
  // pt and eta distribution of photons
  if(kinPhotons_["central"]->size() > 0){
    const Photon& ph = kinPhotons_["central"]->at(0);
    theHistograms->fill("kinPh_lead_pt"  , "p_{t} of Kinematic photons;p_{T} [GeV/c];Events", ph_pt_bins  ,      ph.pt()  , theWeight);
    theHistograms->fill("kinPh_lead_aeta", "|#eta| of Kinematic photons;#eta;Events"        , ph_aeta_bins, fabs(ph.eta()), theWeight);
    theHistograms->fill("kinPh_lead_pt_aeta", "Kinematic photons;p_{t};|#eta|;Events", ph_pt_bins, ph_aeta_bins, ph.pt(), fabs(ph.eta()), theWeight);
    
    if(kinPhotons_["central"]->size() > 1)
      theHistograms->fill("kinPh_sublead_pt", "Subleading Kin #gamma;p_{T} [GeV/c]", 20,0.,100., kinPhotons_["central"]->at(1).pt());
  }
  
  if(loosePhotons_["central"]->size() > 0){
    const Photon& ph = loosePhotons_["central"]->at(0);
    theHistograms->fill("veryLoosePh_lead_pt"  , "p_{t} of VeryLoose photons;p_{T} [GeV/c];Events", ph_pt_bins  ,      ph.pt()  , theWeight);
    theHistograms->fill("veryLoosePh_lead_aeta", "|#eta| of VeryLoose photons;#eta;Events"        , ph_aeta_bins, fabs(ph.eta()), theWeight);
    theHistograms->fill("veryLoosePh_lead_pt_aeta", "VeryLoose photons;p_{t};|#eta|;Events", ph_pt_bins, ph_aeta_bins, ph.pt(), fabs(ph.eta()), theWeight);
    
    if(ph.cutBasedIDLoose()){
      theHistograms->fill("loosePh_lead_pt"  , "p_{t} of LooseID photons;p_{T} [GeV/c];Events", ph_pt_bins  ,      ph.pt()  , theWeight);
      theHistograms->fill("loosePh_lead_aeta", "|#eta| of LooseID photons;#eta;Events"        , ph_aeta_bins, fabs(ph.eta()), theWeight);
      theHistograms->fill("loosePh_lead_pt_aeta", "LooseID photons;p_{t};|#eta|;Events", ph_pt_bins, ph_aeta_bins, ph.pt(), fabs(ph.eta()), theWeight);
    }
    else{
      theHistograms->fill("failPh_lead_pt"  , "p_{t} of Fail photons;p_{T} [GeV/c];Events", ph_pt_bins  ,      ph.pt()  , theWeight);
      theHistograms->fill("failPh_lead_aeta", "|#eta| of Fail photons;#eta;Events"        , ph_aeta_bins, fabs(ph.eta()), theWeight);
      theHistograms->fill("failPh_lead_pt_aeta", "Fail photons;p_{t};|#eta|;Events", ph_pt_bins, ph_aeta_bins, ph.pt(), fabs(ph.eta()), theWeight);
    }
    
    if(loosePhotons_["central"]->size() > 1){
      const Photon& ph_sublead = loosePhotons_["central"]->at(1);
      theHistograms->fill("veryLoosePh_sublead_pt", "Subleading VeryLoose #gamma;p_{T} [GeV/c]", ph_pt_bins, ph_sublead.pt());
      if(ph_sublead.cutBasedIDLoose())
	theHistograms->fill("loosePh_sublead_pt", "Subleading Loose #gamma;p_{T} [GeV/c]", ph_pt_bins, ph_sublead.pt());
      else
	theHistograms->fill("failPh_sublead_pt" , "Subleading Fail #gamma;p_{T} [GeV/c]" , ph_pt_bins, ph_sublead.pt());
    }
  }

  
  // Photons passing kinematic cuts; pt, eta, ID variables
  for(const Photon& ph : *kinPhotons_["central"]){
    if(ph.isBarrel()){
      vector<double> sieie_bins(30);
      for(size_t i = 0; i < sieie_bins.size(); i++) sieie_bins[i] = 0.004 + 0.001*i;
      vector<double> chIso_bins({0., 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4}); //({0., 0.1, 0.25, 0.65, 0.8, 1.141, 1.3, 1.694, 2, 3, 4});
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
      vector<double> chIso_bins({0., 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4}); //({0., 0.1, 0.25, 0.517, 0.8, 1.051, 1.3, 1.6, 2.089, 3, 4});
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
    if(passVeryLoose(ph))     theHistograms->fill("kinPhotons_ID", "Cut Based ID", BINS_KINPHID, 1, theWeight);
    if(ph.cutBasedIDLoose())  theHistograms->fill("kinPhotons_ID", "Cut Based ID", BINS_KINPHID, 2, theWeight);
    if(ph.cutBasedIDMedium()) theHistograms->fill("kinPhotons_ID", "Cut Based ID", BINS_KINPHID, 3, theWeight);
    if(ph.cutBasedIDTight())  theHistograms->fill("kinPhotons_ID", "Cut Based ID", BINS_KINPHID, 4, theWeight);
  }
  
  // Systematics
  for(auto & [syst, phVect] : kinPhotons_)
    if(phVect->size() > 0)
      // Creating void histograms, then filling alphanumeric labels --> new ones are created as they are encountered
      theHistograms->book<TH1F>("kinPh_eScale_N", "Number of #gamma passing selection", 1,0,0)->Fill(Form("kin_%s" , syst), theWeight);
  
  for(auto & [syst, phVect] : goodPhotons_)
    if(phVect->size() > 0)
      theHistograms->book<TH1F>("kinPh_eScale_N", "Number of #gamma passing selection", 1,0,0)->Fill(Form("good_%s", syst), theWeight);
  
  // How many photons are there in each event?
  for(const auto& [syst, phVect] : goodPhotons_)
    theHistograms->fill(Form("loosePh_%s_N", syst), Form("Number of Loose photons (%s)", syst), 5,-0.5,4.5, phVect->size(), theWeight);
  
  for(const auto& [syst, phVect] : loosePhotons_)
    theHistograms->fill(Form("veryLoosePh_%s_N", syst), Form("Number of VeryLoose photons (%s)", syst), 5,-0.5,4.5, phVect->size(), theWeight);

  for(const auto& [syst, phVect] : kinPhotons_)
    theHistograms->fill(Form(  "kinPh_%s_N", syst), Form("Number of Kinematic photons (%s)"  , syst), 5,-0.5,4.5, phVect->size(), theWeight);
  
  // Test data/MC for each cut separately
  if(kinPhotons_["central"]->size() == 0){
    theHistograms->fill("kinPhotons_cuts", "Single cut efficiency;;Events", BINS_PHCUTFLOW, 0, theWeight);  // No photons
    return;
  }
  else{
    
    Photon::IdWp wp = Photon::IdWp::Loose;
    
    bool b_sieie  = bestKinPh_->cutBasedID(wp, Photon::IDcut::sieie );
    bool b_HoverE = bestKinPh_->cutBasedID(wp, Photon::IDcut::HoverE);
    bool b_chIso  = bestKinPh_->cutBasedID(wp, Photon::IDcut::chIso );
    bool b_neIso  = bestKinPh_->cutBasedID(wp, Photon::IDcut::neIso );
    bool b_phIso  = bestKinPh_->cutBasedID(wp, Photon::IDcut::phIso );
    UInt_t nCutsPass = bestKinPh_->nCutsPass(wp);
    
    // Single cut efficiency
    theHistograms->fill("kinPhotons_cuts", "Single cut efficiency;;Events", BINS_PHCUTFLOW, 1, theWeight);  // All photons
    if(b_sieie)
      theHistograms->fill("kinPhotons_cuts", "Single cut efficiency;;Events", BINS_PHCUTFLOW, 2, theWeight);
    if(b_HoverE)
      theHistograms->fill("kinPhotons_cuts", "Single cut efficiency;;Events", BINS_PHCUTFLOW, 3, theWeight);
    if(b_chIso)
      theHistograms->fill("kinPhotons_cuts", "Single cut efficiency;;Events", BINS_PHCUTFLOW, 4, theWeight);
    if(b_neIso)
      theHistograms->fill("kinPhotons_cuts", "Single cut efficiency;;Events", BINS_PHCUTFLOW, 5, theWeight);
    if(b_phIso)
      theHistograms->fill("kinPhotons_cuts", "Single cut efficiency;;Events", BINS_PHCUTFLOW, 6, theWeight);
    
    // N-1 efficiency of the cuts
    if(nCutsPass >= 4)
      theHistograms->fill("kinPhotons_Nm1", "N-1 cut efficiency;;Events", BINS_PHCUTNM1, 0, theWeight);
    if(           b_HoverE && b_chIso && b_neIso && b_phIso)
      theHistograms->fill("kinPhotons_Nm1", "N-1 cut efficiency;;Events", BINS_PHCUTNM1, 1, theWeight);
    if(b_sieie             && b_chIso && b_neIso && b_phIso)
      theHistograms->fill("kinPhotons_Nm1", "N-1 cut efficiency;;Events", BINS_PHCUTNM1, 2, theWeight);
    if(b_sieie && b_HoverE            && b_neIso && b_phIso)
      theHistograms->fill("kinPhotons_Nm1", "N-1 cut efficiency;;Events", BINS_PHCUTNM1, 3, theWeight);
    if(b_sieie && b_HoverE && b_chIso            && b_phIso)
      theHistograms->fill("kinPhotons_Nm1", "N-1 cut efficiency;;Events", BINS_PHCUTNM1, 4, theWeight);
    if(b_sieie && b_HoverE && b_chIso && b_neIso           )
      theHistograms->fill("kinPhotons_Nm1", "N-1 cut efficiency;;Events", BINS_PHCUTNM1, 5, theWeight);
    
    theHistograms->fill("kinPhotons_MVA", "kinematic #gamma MVA", 41,-1.025,1.025, bestKinPh_->MVAvalue(), theWeight);
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
  theHistograms->fill("AK8_N" , "# jets AK8", 6, -0.5, 5.5, jetsAK8->size(), theWeight);
  for(const Jet& jet : *jetsAK8)
    theHistograms->fill("AK8_pt", "p_{T} jets AK8;GeV/c", 50, 0., 500., jet.pt(), theWeight);
  
  vector<std::pair<const Particle*, const Jet*>> JetsAK8genrec = matchDeltaR(*genJetsAK8, *jetsAK8);
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
  theHistograms->fill("AK4_N" , "# jets AK4", 6, -0.5, 5.5, jets->size(), theWeight);
  for(const Jet& jet : *jets)
    theHistograms->fill("AK4_pt", "p_{T} jets AK4;GeV/c", 50, 0., 500., jet.pt(), theWeight);

  vector<std::pair<const Particle*, const Jet*>> JetsAK4genrec = matchDeltaR(*genJets, *jets);
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


void VVGammaAnalyzer::plotsVVGstatus(const char* name, const char* title, const TLorentzVector& p4, const char* mType){
  // return either mass or transverse mass, depending on the request
  auto mValue = (!strcmp(mType, "massT") ? [](const TLorentzVector& v){ return v.Mt() ; } : [](const TLorentzVector& v){ return v.M() ; });
  
  bool isSingleBoson = !strcmp(name, "Z");
  const vector<double>& binsVV  = isSingleBoson ? mZ_bins  : mVV_bins;
  const vector<double>& binsVVG = isSingleBoson ? mZG_bins : mVVG_bins;
  
  // Kin photons
  if(kinPhotons_["central"]->size() > 0){
    const TLorentzVector& ph_p4 = kinPhotons_["central"]->front().p4();
    theHistograms->fill(Form("%s_%s_kinPh" , name, mType), Form("%s %s with Kin #gamma" , title, mType), binsVV , mValue(p4      ), theWeight);
    theHistograms->fill(Form("%sG_%s_kinPh", name, mType), Form("%sG %s with Kin #gamma", title, mType), binsVVG, mValue(p4+ph_p4), theWeight);
  }
  // No photon
  else
    theHistograms->fill(Form("%s_%s_noPh", name, mType), Form("%s mass with No #gamma", title), binsVV, mValue(p4), theWeight);

  // Loose photons
  if(loosePhotons_["central"]->size() > 0){
    const TLorentzVector& ph_p4 = loosePhotons_["central"]->front().p4();
    theHistograms->fill(Form("%s_%s_veryLoosePh" , name, mType), Form("%s %s with VeryLoose #gamma" , title, mType), binsVV , mValue(p4      ), theWeight);
    theHistograms->fill(Form("%sG_%s_veryLoosePh", name, mType), Form("%sG %s with VeryLoose #gamma", title, mType), binsVVG, mValue(p4+ph_p4), theWeight);
    // Tight photons
    if(goodPhotons_["central"]->size() > 0){
      const TLorentzVector& ph_p4 = goodPhotons_["central"]->front().p4();
      theHistograms->fill(Form("%s_%s_loosePh" , name, mType), Form("%s %s with LooseID #gamma" , title, mType), binsVV , mValue(p4      ), theWeight);
      theHistograms->fill(Form("%sG_%s_loosePh", name, mType), Form("%sG %s with LooseID #gamma", title, mType), binsVVG, mValue(p4+ph_p4), theWeight);
    }
    // Fail photon (loose && !tight)
    else{
      theHistograms->fill(Form("%s_%s_failPh" , name, mType), Form("%s %s with Fail #gamma" , title, mType), binsVV , mValue(p4      ), theWeight);
      theHistograms->fill(Form("%sG_%s_failPh", name, mType), Form("%sG %s with Fail #gamma", title, mType), binsVVG, mValue(p4+ph_p4), theWeight);
    }
  }
}


void debugPhotonID(const Photon& ph){
  Photon::IdWp wp;
  if(ph.cutBasedIDLoose()){
    wp = Photon::IdWp::Loose;
    bool HoverE           = ph.passHoverE(wp);
    bool sigmaiEtaiEta    = ph.passSigmaiEtaiEta(wp);
    bool chargedIsolation = ph.passChargedIsolation(wp);
    bool neutralIsolation = ph.passNeutralIsolation(wp);
    bool photonIsolation  = ph.passPhotonIsolation(wp);
    if( HoverE+sigmaiEtaiEta+chargedIsolation+neutralIsolation+photonIsolation < 5 )
      cout << ">>>Loose\t" <<
	"HoverE: "           << HoverE           << " - " <<
	"sigmaiEtaiEta: "    << sigmaiEtaiEta    << " - " <<
	"chargedIsolation: " << chargedIsolation << " - " <<
	"neutralIsolation: " << neutralIsolation << " - " <<
	"photonIsolation: "  << photonIsolation  << '\n';
  }
  else
    return;
  
  if(ph.cutBasedIDMedium()){
    wp = Photon::IdWp::Medium;
    bool HoverE           = ph.passHoverE(wp);
    bool sigmaiEtaiEta    = ph.passSigmaiEtaiEta(wp);
    bool chargedIsolation = ph.passChargedIsolation(wp);
    bool neutralIsolation = ph.passNeutralIsolation(wp);
    bool photonIsolation  = ph.passPhotonIsolation(wp);
    if( HoverE+sigmaiEtaiEta+chargedIsolation+neutralIsolation+photonIsolation < 5 )
      cout << ">>>Medium\t" <<
	"HoverE: "           << HoverE           << " - " <<
	"sigmaiEtaiEta: "    << sigmaiEtaiEta    << " - " <<
	"chargedIsolation: " << chargedIsolation << " - " <<
	"neutralIsolation: " << neutralIsolation << " - " <<
	"photonIsolation: "  << photonIsolation  << '\n';
  }
  else
    return;
  
  if(ph.cutBasedIDTight()){
    wp = Photon::IdWp::Tight;
    bool HoverE           = ph.passHoverE(wp);
    bool sigmaiEtaiEta    = ph.passSigmaiEtaiEta(wp);
    bool chargedIsolation = ph.passChargedIsolation(wp);
    bool neutralIsolation = ph.passNeutralIsolation(wp);
    bool photonIsolation  = ph.passPhotonIsolation(wp);
    if( HoverE+sigmaiEtaiEta+chargedIsolation+neutralIsolation+photonIsolation < 5 )
      cout << ">>>Tight\t" <<
	"HoverE: "           << HoverE           << " - " <<
	"sigmaiEtaiEta: "    << sigmaiEtaiEta    << " - " <<
	"chargedIsolation: " << chargedIsolation << " - " <<
	"neutralIsolation: " << neutralIsolation << " - " <<
	"photonIsolation: "  << photonIsolation  << '\n';
  }
}


void VVGammaAnalyzer::photonFakeRate(){
  const vector<Photon>& thePhVect = *loosePhotons_["central"];  // *kinPhotons_["central"];
  if(thePhVect.size() == 0)
    return;
  
  Photon::IdWp wp = Photon::IdWp::Loose;
  auto bestG = std::max_element(thePhVect.begin(), thePhVect.end(),
				[wp](const Photon& a, const Photon& b){ return a.nCutsPass(wp) < b.nCutsPass(wp);
				});  // max_element returns the first among those with max value --> preserve pt ordering
  
  if( !bestG->cutBasedID(wp,Photon::IDcut::HoverE) || !bestG->cutBasedID(wp,Photon::IDcut::neIso) || !bestG->cutBasedID(wp,Photon::IDcut::phIso) )
    throw std::logic_error("Best of loose photon does not pass one of {HoverE, neutralIso, photonIso}");
  
  bool isPrompt = false;
  if(theSampleInfo.isMC())
    isPrompt = std::any_of(genPhotons_->begin(), genPhotons_->end(), 
			   [bestG](const Particle& gen){ return gen.genStatusFlags().test(phys::isPrompt) && physmath::deltaR(*bestG, gen) < 0.2; }
			   );
  
  // Debug photon ID
  // debugPhotonID(*bestG);
  
  // Actual fake rate histograms
  double theAeta = fabs(bestG->eta());
  double thePt   = bestG->pt();
  if(thePt > ph_pt_bins.back())
    thePt = ph_pt_bins.back() - 0.1;
  
  // Fake rate with ABCD
  char ABCD = phABCD(*bestG, wp);

  const char* name = Form("PhFR_%c", ABCD);
  if(theSampleInfo.isMC())
    name = Form("%s_%s", name, (isPrompt ? "nonprompt" : "prompt"));
  
  theHistograms->fill(name, Form("Photons: %c;p_{T} [GeV/c];#eta", ABCD), ph_pt_bins, ph_aeta_bins, thePt, theAeta, theWeight);
  
  // Fake rate with Loose-to-Tight (loose = 3/5 cutBasedIDLoose, tight = 5/5 cutBasedIDLoose)
  int nPass = bestG->nCutsPass(wp);

  name = Form("PhFR_LtoT_%d", nPass);
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
    theHistograms->fill(Form("PhFRClosure_LtoT_PASS_%s"      , varName), "Closure test: PASS"    , mVVG_bins, varValue, theWeight);
  else if(nPass == 3 || nPass == 4){  // LOOSE && !TIGHT
    theHistograms->fill(Form("PhFRClosure_LtoT_FAIL_%s"      , varName), "Closure test: FAIL"    , mVVG_bins, varValue, theWeight);
    theHistograms->fill(Form("PhFRClosure_LtoT_reweighted_%s", varName), "Closure test: (3+4)*SF", mVVG_bins, varValue, weight   );
  }
}


void VVGammaAnalyzer::photonFakeRate_extended(){
  if(!bestKinPh_)
    return;

  double theAeta = fabs(bestKinPh_->eta());
  double thePt   = bestKinPh_->pt();
  if(thePt > ph_pt_bins.back())
    thePt = ph_pt_bins.back() - 0.1;

  bool isPassVL = bestKinPh_->cutBasedID(Photon::IdWp::VeryLoose);
  bool isPassLoose = bestKinPh_->cutBasedIDLoose();
  const char* strPassVL = isPassVL ? "PASS" : "FAIL";
  const char* strPrompt = "";
  if(theSampleInfo.isMC()){
    bool isPrompt = std::any_of(genPhotons_->begin(), genPhotons_->end(), 
  				[this](const Particle& gen){ return gen.genStatusFlags().test(phys::isPrompt) && physmath::deltaR(*bestKinPh_, gen) < 0.2; }
  				);
    strPrompt = isPrompt ? "_nonprompt" : "_prompt";
  }

  // PhFR_KToL_[FAIL/PASS]
  const char* name_KtoVL = Form("PhFR_KtoVL%s_%s", strPrompt, strPassVL);
  theHistograms->fill(name_KtoVL, "Photon fake rate Kin to VeryLoose;p_{T} [GeV/c];#eta", ph_pt_bins, ph_aeta_bins, thePt, theAeta, theWeight);

  // PhFR_KToLexcl_[FAIL/PASS]
  if(!isPassLoose){
    const char* name_KtoVLexcl = Form("PhFR_KtoVLexcl%s_%s", strPrompt, strPassVL);
    theHistograms->fill(name_KtoVLexcl, "Photon fake rate Kin to VeryLoose&&!Loose;p_{T} [GeV/c];#eta", ph_pt_bins, ph_aeta_bins, thePt, theAeta, theWeight);
  }

  // Closure test extended
  const char* varName;
  double varValue;
  if     (is4Lregion(region_)){
    varName = "mZZG";
    varValue = ( ZZ->p4() + bestKinPh_->p4() ).M();
  }
  else if(is3Lregion(region_)){
    varName = "mWZG";
    varValue = ( ZW->p4() + bestKinPh_->p4() ).M();
  }
  else
    return;
  
  // Kin --> VeryLoose
  double f_KToVL = getPhotonFR_KtoVL(*bestKinPh_);
  double weight_KtoVL = theWeight * f_KToVL/(1-f_KToVL);
  
  if(isPassVL)
    theHistograms->fill(Form("PhFRClosure_KtoVL_PASS_%s"      , varName), "Closure test KtoVL: PASS"   , mVVG_bins, varValue, theWeight);
  else{
    theHistograms->fill(Form("PhFRClosure_KtoVL_FAIL_%s"      , varName), "Closure test KtoVL: FAIL"   , mVVG_bins, varValue, theWeight);
    theHistograms->fill(Form("PhFRClosure_KtoVL_reweighted_%s", varName), "Closure test KtoVL: FAIL*SF", mVVG_bins, varValue, weight_KtoVL);
  }
  
  // Kin && !Loose --> VeryLoose && !Loose
  if(!isPassLoose){
    double f_KtoVLexcl = getPhotonFR_KtoVLexcl(*bestKinPh_);
    double weight_KtoVLexcl = theWeight * f_KtoVLexcl/(1-f_KtoVLexcl);
  
    if     (isPassVL)
      theHistograms->fill(Form("PhFRClosure_KtoVLexcl_PASS_%s"      , varName), "Closure test KtoVLexcl: PASS"   , mVVG_bins, varValue, theWeight);
    else{
      theHistograms->fill(Form("PhFRClosure_KtoVLexcl_FAIL_%s"      , varName), "Closure test KtoVLexcl: FAIL"   , mVVG_bins, varValue, theWeight);
      theHistograms->fill(Form("PhFRClosure_KtoVLexcl_reweighted_%s", varName), "Closure test KtoVLexcl: FAIL*SF", mVVG_bins, varValue, weight_KtoVLexcl);
    }
  }
}

void VVGammaAnalyzer::studyJetsChoice(){
  // Study how to choose the correct jet AK8 in case it exists. The problem of rejecting 
  // events in which it does not exist is part of the event selection
  if( !theSampleInfo.isMC() )
    return;
  Boson<Particle>* qq = nullptr;
  if     (genZhadCandidates_->size() > 0)
    qq = &(genZhadCandidates_->at(0));
  else if(genWhadCandidates_->size() > 0)
    qq = &(genWhadCandidates_->at(0));
  else
    return;
  
  // int statusAK8;
  int statusAK4 = studyAK4Choice(fAK4_, *qq, 0.4);
  if(statusAK4 != 0)
    // statusAK8 = studyAK8Choice(fAK8_, *qq, 0.4);
    studyAK8Choice(fAK8_, *qq, 0.4);
  
}


int VVGammaAnalyzer::studyAK4Choice(std::ofstream& fout, const phys::Boson<phys::Particle>& diquark, const double& tolerance){
  vector<Particle> vQuarks({*diquark.daughterPtr(0), *diquark.daughterPtr(1)});
  efficiency(*jets, vQuarks, "AK4"   , "quarks", tolerance);
  efficiency(*genJets, vQuarks, "genAK4", "quarks", tolerance);

  if(jets->size() <= 2)
    return 1;  // Not enough jets
  
  // vector<std::pair<const Particle*, const Jet*>> vGenRec = matchDeltaR(vQuarks, *jets, tolerance);
  
  return 0;
}


int VVGammaAnalyzer::studyAK8Choice(std::ofstream& fout, const phys::Boson<phys::Particle>& diquark, const double& tolerance){
  vector<Jet>::const_iterator rec = std::min_element(jetsAK8->cbegin(), jetsAK8->cend(), DeltaRComparator(diquark));
  if(rec == jets->cend() || deltaR(*rec, diquark) > tolerance)
    return 1;  // either jetsAK8->size() == 0 or dR > tol

  // for(vector<Jet>::const_iterator rec_it = jetsAK8->cbegin(); jetsAK8 != jetsAK8.cend(); ++rec_it){
  //   fout << (rec_it == recAK8)           << ','
  // 	 << theWeight                    << ','
  // 	 << rec_it->mass()               << ','
  // 	 << rec_it->particleNet().WvsQCD << ','
  // 	 << rec_it->particleNet().ZvsQCD << ','
  // 	 << rec_it->particleNet().TvsQCD << ','
  // 	 << rec_it->deepAK8_MD().WvsQCD  << ','
  // 	 << rec_it->deepAK8_MD().ZvsQCD  << ','
  // 	 << rec_it->deepAK8_MD().TvsQCD  << '\n';
  
  return 0;
}


template <class PAR>
void VVGammaAnalyzer::efficiency(const vector<PAR>& vRec, const vector<Particle>& vGen, const char* recLabel, const char* genLabel, double tolerance){
  // Reconstruction efficiency and resolution
  // cout << "\tefficiency " << recLabel << ' ' << genLabel << '\n';

  TString tNameEff = TString::Format("Eff_%s_%s_%s_%s", recLabel, genLabel, "%s", "%s");  // e.g. Eff_photos_loose_den_pt
  TString tNameRes = TString::Format("Res_%s_%s_%s"   , recLabel, genLabel, "%s");         // e.g. Res_AK8_all_eta
  const char* nameEff = tNameEff.Data();  // The Form() in the global namespace uses a circulary buffer that can be overwritten
  const char* nameRes = tNameRes.Data();  // storing and using the returned pointer is unsafe

  vector<std::pair<const Particle*, const PAR*>> vGenRec = matchDeltaR(vGen, vRec, tolerance);  // intrinsic threshold of deltaR = 0.2
  size_t nRec = vRec.size();
  
  for(auto & [gen, rec] : vGenRec){
    theHistograms->fill(Form(nameEff, "DEN", "pt" ), "DEN p_{T};p_{T} GEN", 25,0.,250., gen->pt() , theWeight);
    theHistograms->fill(Form(nameEff, "DEN", "E"  ), "DEN energy;E GEN"   , 25,0.,500., gen->e()  , theWeight);
    theHistograms->fill(Form(nameEff, "DEN", "eta"), "DEN eta;#eta GEN"   , eta_bins  , gen->eta(), theWeight);
    theHistograms->fill(Form(nameEff, "DEN", "N"  ), "DEN # rec;# rec"    , 6,-0.5,5.5, nRec      , theWeight);
    
    if(rec == nullptr)
      continue;
    theHistograms->fill(Form(nameEff, "NUM", "pt" ), "NUM p_{T};p_{T} GEN", 25,0.,250., gen->pt() , theWeight);
    theHistograms->fill(Form(nameEff, "NUM", "E"  ), "NUM energy;E GEN"   , 25,0.,500., gen->e()  , theWeight);
    theHistograms->fill(Form(nameEff, "NUM", "eta"), "NUM eta;#eta GEN"   , eta_bins  , gen->eta(), theWeight);
    theHistograms->fill(Form(nameEff, "NUM", "N"  ), "NUM # rec;# rec"    , 6,-0.5,5.5, nRec      , theWeight);
    
    double deltaR = physmath::deltaR(*rec, *gen);
    double deltaEoverE   = rec->e()  / gen->e()  - 1;  // (rec->e() - gen->e()) / gen->e();
    double deltapToverpT = rec->pt() / gen->pt() - 1;  // (rec->pt() - gen->pt()) / gen->pt();
    theHistograms->fill(Form(nameRes, "dR"  ), "Resolution #DeltaR;#DeltaR"        , 20, 0.,0.4, deltaR       , theWeight);
    theHistograms->fill(Form(nameRes, "E"   ), "Resolution Energy;#DeltaE/E"       , 20,-1.,1. , deltaEoverE  , theWeight);
    theHistograms->fill(Form(nameRes, "pt"  ), "Resolution p_{T};#Deltap_{T}/p_{T}", 20,-1.,1. , deltapToverpT, theWeight);
    theHistograms->fill(Form(nameRes, "EvsE"), "Resolution;E;#DeltaE/E", 20,0.,400., 16,-0.4,0.4, gen->e(), deltaEoverE , theWeight);
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
  // double Ev  = pv.E();
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


char VVGammaAnalyzer::phABCD(const Photon& ph, const Photon::IdWp wp){
  char ABCD = '0';  // Defaults to something recognizable
  if( ph.cutBasedID(wp, Photon::IDcut::chIso) ){    // Signal region: either A or B
    if     (ph.cutBasedID(wp, Photon::IDcut::sieie)) ABCD = 'A';
    else if(isVeryLooseSiEiE(ph)                   ) ABCD = 'B';
  }
  else{                                 // Measurement region
    if     (ph.cutBasedID(wp, Photon::IDcut::sieie)) ABCD = 'C';
    else if(isVeryLooseSiEiE(ph)                   ) ABCD = 'D';
  }
  return ABCD;
}


char phABCD_study(const phys::Photon&, const double& barrel_thr, const double& endcap_thr){
  return '0';
}


void VVGammaAnalyzer::SYSplots(const char* syst, const double& weight, const Photon* ph){
  if(is4Lregion(region_)){
    theHistograms->fill(  Form("SYS_mZZ_%s" , syst), Form("m_{ZZ} %s"      , syst), mVV_bins , ZZ->mass()               , weight);
    
    if(ph){
      double w = weight * ph->efficiencySF();
      theHistograms->fill(Form("SYS_mZZGloose_%s", syst), Form("m_{ZZ#gamma} %s", syst), mVVG_bins, (ZZ->p4() + ph->p4()).M(), w);
      if(ph->cutBasedIDMedium())
      	theHistograms->fill(Form("SYS_mZZGmedium_%s", syst), Form("m_{ZZ#gamma} %s", syst), mVVG_bins, (ZZ->p4() + ph->p4()).M(), w);
      if(ph->cutBasedIDTight())
      	theHistograms->fill(Form("SYS_mZZGtight_%s" , syst), Form("m_{ZZ#gamma} %s", syst), mVVG_bins, (ZZ->p4() + ph->p4()).M(), w);
    }
  }
  else if(is3Lregion(region_)){
    theHistograms->fill(  Form("SYS_mWZ_%s" , syst), Form("m_{WZ} %s"      , syst), mVV_bins , ZW->mass()               , weight);
    
    if(ph){
      double w = weight * ph->efficiencySF();
      theHistograms->fill(Form("SYS_mWZGloose_%s", syst), Form("m_{WZ#gamma} %s", syst), mVVG_bins, (ZW->p4() + ph->p4()).M(), w);
      if(ph->cutBasedIDMedium())
      	theHistograms->fill(Form("SYS_mWZGmedium_%s", syst), Form("m_{ZZ#gamma} %s", syst), mVVG_bins, (ZW->p4() + ph->p4()).M(), w);
      if(ph->cutBasedIDTight())
      	theHistograms->fill(Form("SYS_mWZGtight_%s" , syst), Form("m_{ZZ#gamma} %s", syst), mVVG_bins, (ZW->p4() + ph->p4()).M(), w);
    }
    
  }
}

void VVGammaAnalyzer::systematicsStudy(){
  double base_w = theWeight; // * ph.efficiencySF_
  
  const Photon* ph = nullptr;
  if(goodPhotons_["central"]->size() >= 1)
    ph = & (goodPhotons_["central"]->front());
  
  // central
  SYSplots("central", base_w, ph);
  
  // Photons energy scale and resolution
  for(const auto& [syst, phVect] : goodPhotons_){
    if(strcmp(syst, "central") == 0) continue;
    if(phVect->size() == 0) continue;
    SYSplots(Form("ph%s", syst), base_w, & phVect->front());
  }
  
  bool isMC = theSampleInfo.isMC();
  // puWeightUnc
  SYSplots("puWeight_Up"  , base_w * ( isMC ? theSampleInfo.puWeightUncUp() / theSampleInfo.puWeight() : 1.), ph);
  SYSplots("puWeight_Down", base_w * ( isMC ? theSampleInfo.puWeightUncDn() / theSampleInfo.puWeight() : 1.), ph);
  
  // L1PrefiringWeight
  SYSplots("L1Prefiring_Up"  , base_w * ( isMC ? theSampleInfo.L1PrefiringWeightUp() / theSampleInfo.L1PrefiringWeight() : 1.), ph);
  SYSplots("L1Prefiring_Down", base_w * ( isMC ? theSampleInfo.L1PrefiringWeightDn() / theSampleInfo.L1PrefiringWeight() : 1.), ph);
  
  // QCD scale
  SYSplots("QCDscaleF_Up"    , base_w * ( isMC ? theSampleInfo.QCDscale_muR1F2()   : 1.), ph);  // "QCDscale_muR1F2"
  SYSplots("QCDscaleF_Down"  , base_w * ( isMC ? theSampleInfo.QCDscale_muR1F0p5() : 1.), ph);  // "QCDscale_muR1F0p5"
  SYSplots("QCDscalemuR_Up"  , base_w * ( isMC ? theSampleInfo.QCDscale_muR2F1()   : 1.), ph);  // "QCDscale_muR2F1"
  SYSplots("QCDscalemuR_Down", base_w * ( isMC ? theSampleInfo.QCDscale_muR0p5F1() : 1.), ph);  // "QCDscale_muR0p5F1"
  
  // PDF var
  SYSplots("PDFVar_Up"  , base_w * ( isMC ? theSampleInfo.PDFVar_Up()   : 1.), ph);
  SYSplots("PDFVar_Down", base_w * ( isMC ? theSampleInfo.PDFVar_Down() : 1.), ph);
  
  // alphas MZ
  SYSplots("alphas_Up"  , base_w * ( isMC ? theSampleInfo.alphas_MZ_Up()   : 1.), ph);
  SYSplots("alphas_Down", base_w * ( isMC ? theSampleInfo.alphas_MZ_Down() : 1.), ph);
  

  double eleEff_w=0., muoEff_w=0., eleFake_w=0., muoFake_w=0.;
  if     (is4Lregion(region_)){
    eleEff_w  = ZZ->eleEffSFUnc()/ZZ->efficiencySF();
    muoEff_w  = ZZ->muoEffSFUnc()/ZZ->efficiencySF();
    eleFake_w = ZZ->eleFakeRateSFUnc()/ZZ->fakeRateSF();
    muoFake_w = ZZ->muoFakeRateSFUnc()/ZZ->fakeRateSF();
  }
  else if(is3Lregion(region_)){
    eleEff_w  = ZW->eleEffSFUnc()/ZW->efficiencySF();
    muoEff_w  = ZW->muoEffSFUnc()/ZW->efficiencySF();
    eleFake_w = ZW->eleFakeRateSFUnc()/ZW->fakeRateSF();
    muoFake_w = ZW->muoFakeRateSFUnc()/ZW->fakeRateSF();
  }
  else if(is2Lregion(region_)){
    eleEff_w  = Z->eleEffSFUnc()/Z->efficiencySF();
    muoEff_w  = Z->muoEffSFUnc()/Z->efficiencySF();
    eleFake_w = Z->eleFakeRateSFUnc()/Z->fakeRateSF();
    muoFake_w = Z->muoFakeRateSFUnc()/Z->fakeRateSF();
  }
  
  // lepton efficiency SF
  SYSplots("eleEffSF_Up"  , base_w * (1 + eleEff_w), ph);
  SYSplots("eleEffSF_Down", base_w * (1 - eleEff_w), ph);
  SYSplots("muoEffSF_Up"  , base_w * (1 + muoEff_w), ph);
  SYSplots("muoEffSF_Down", base_w * (1 - muoEff_w), ph);
  
  // lepton fake rate SF
  SYSplots("eleFakeRateSF_Up"  , base_w * (1 + eleFake_w), ph);
  SYSplots("eleFakeRateSF_Down", base_w * (1 - eleFake_w), ph);
  SYSplots("muoFakeRateSF_Up"  , base_w * (1 + muoFake_w), ph);
  SYSplots("muoFakeRateSF_Down", base_w * (1 - muoFake_w), ph);
  
  if(ph){
    // Photons ID efficiency
    double phEff_w = ph->efficiencySF()/ph->efficiencySFUnc();
    SYSplots("phEffSF_Up"  , base_w * (1 + phEff_w), ph);
    SYSplots("phEffSF_Down", base_w * (1 - phEff_w), ph);
    // Photon FR uncertaintiy  WARN: for this to have a meaning, the photon FR SF should be applied
    double phFake_w = getPhotonFRUnc(*ph)/getPhotonFR(*ph);
    SYSplots("phFakeRate_Up"  , base_w * (1 + phFake_w), ph);
    SYSplots("phFakeRate_Down", base_w * (1 - phFake_w), ph);
  }
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


void VVGammaAnalyzer::photonGenStudy(){
  // Study the efficiency/background rejection of a cut in deltaR(photon, lepton)
  // Make collection of genPhotons from hard process
  vector<Particle> genPhotonsPrompt;
  genPhotonsPrompt.reserve(genPhotons_->size());
  copy_if(genPhotons_->begin(), genPhotons_->end(), std::back_inserter(genPhotonsPrompt),
	  [](const Particle& p){
	    std::bitset<15> flags = p.genStatusFlags();
	    return /*status == 1 &&*/ (flags.test(isHardProcess) || flags.test(fromHardProcess)) /*&& fromHardProcessFinalState*/; }
	  );
  
  // Since the base selection already contains this cut, it is necessary to do it again
  std::map<const char*, vector<Photon>> mapWPtoPhotons;
  // mapWPtoPhotons["kin"]       = vector<Photon>();
  // mapWPtoPhotons["veryLoose"] = vector<Photon>();
  mapWPtoPhotons["loose"]     = vector<Photon>();
    
  // Select reco photons except cut on DR(ph, lep)
  for(const Photon& ph: *photons){
    //Pixel seed and electron veto
    if(ph.hasPixelSeed() || !ph.passElectronVeto()) continue;
		
    //Kinematic selection
    if(ph.pt() < 20) continue;
    float ph_aeta = fabs(ph.eta());
    if(ph_aeta > 2.4 || (ph_aeta > 1.4442 && ph_aeta < 1.566)) continue;
    
    // mapWPtoPhotons["kin"].push_back(ph)  // kin.push_back(ph);
    // if(passVeryLoose(ph))
    //   mapWPtoPhotons["veryLoose"].push_back(ph)  // veryLoose.push_back(ph);
    if(ph.cutBasedIDLoose())
      mapWPtoPhotons["loose"].push_back(ph);  // loose.push_back(ph);
  }
  
  // The rec-gen matching is done separately for kin, veryLoose and loose, since we don't want e.g. a kin reco photon 
  // that fails the LooseID to be selected just because it comes before in the photon vector
  for(auto & [wp, phVect] : mapWPtoPhotons){
    if(phVect.size() == 0)
      continue;
    
    const Photon* best = & phVect.at(0);
    const char* matched = "nomatch";
    if(genPhotonsPrompt.size() > 0){
      for(const Photon& rec : phVect){
        auto closestGen = closestDeltaR(rec, *genPhotons_);  // std::vector<phys::Particle>::const_iterator
	if(deltaR(rec, *closestGen) < 0.1){  //closestGen != genPhotons_->cend() &&
	  best = &rec;
	  matched = "promptm";  // Matched to prompt photon
	  break;
	}
      }
    }
    else if(genPhotons_->size() > 0){
      for(const Photon& rec : phVect){
	auto closestGen = closestDeltaR(rec, *genPhotons_);
	if(deltaR(rec, *closestGen) < 0.1){
	  best = &rec;
	  matched = "matched";  // Matched to gen photon, not from hard scattering
	  break;
	}
      }
    }
    
    const char* hname = Form("PhGenStudy_%s_%s_%s", "%s", wp, matched);
    
    // Find the closest lep and get the DR
    auto closestLep = closestDeltaR(*best, *leptons_);
    float dR = closestLep != leptons_->cend() ? deltaR(*best, *closestLep) : 10;
    theHistograms->fill(Form(hname, "DRLep"), Form("%s;min #DeltaR(#gamma, l);Events", wp), 50, 0., 1., dR   , theWeight);
    
    // Closest jet
    auto closestJet = closestDeltaR(*best, *jets);
    float dRJet = closestJet != jets->cend() ? deltaR(*best, *closestJet) : 10;
    theHistograms->fill(Form(hname, "DRJet"), Form("%s;#DeltaR(#gamma, j);Events"    , wp), 60, 0., 3 , dRJet, theWeight);
    
    // 2D plot
    theHistograms->fill(Form(hname, "DRLepJet"), Form("%s;#DeltaR(#gamma, l);#DeltaR(#gamma, j);Events", wp), 
			50,0.,1.,
			30,0.,3.,
			dR, dRJet, theWeight);
    
    // Distance to ZZ, Z0, Z1
    theHistograms->fill(Form(hname, "DRZZ"), Form("%s;#DeltaR(#gamma, ZZ);Events", wp), 60, 0., 3., deltaR(*best, *ZZ         ), theWeight);
    theHistograms->fill(Form(hname, "DRZ0"), Form("%s;#DeltaR(#gamma, Z0);Events", wp), 60, 0., 3., deltaR(*best, ZZ->first() ), theWeight);
    theHistograms->fill(Form(hname, "DRZ1"), Form("%s;#DeltaR(#gamma, Z1);Events", wp), 60, 0., 3., deltaR(*best, ZZ->second()), theWeight);
    
    // Pt
    double pt = best->pt() > ph_pt_bins.back() ? ph_pt_bins.back() : best->pt();
    theHistograms->fill(Form(hname, "ptPh") , Form("%s;#gamma p_{T};Events"          , wp), 20,20.,120., pt   , theWeight);
    
    // mVVG
    if(is4Lregion(region_))
      theHistograms->fill(Form(hname, "m4lG"), Form("%s;m_{4l#gamma};Events"         , wp), mVVG_bins, (ZZ->p4() + best->p4()).M() , theWeight);
    if(is3Lregion(region_))
      theHistograms->fill(Form(hname, "mT3lG"), Form("%s;mT_{3l#gamma};Events"       , wp), mVVG_bins, (ZW->p4() + best->p4()).Mt(), theWeight);
    
    // # of Jets
    theHistograms->fill(Form(hname, "nJets"), Form("%s;nJets;Events"                 , wp), 5,0,5, jets->size(), theWeight);
  } // end loop on mapWPtoPhotons
}


// Utilities
double VVGammaAnalyzer::getPhotonFR   (const phys::Photon& ph) const{
  return hPhotonFR_->GetBinContent(hPhotonFR_->FindBin(ph.pt(), abs(ph.eta())));
}

double VVGammaAnalyzer::getPhotonFRUnc(const phys::Photon& ph) const{
  double FRError = hPhotonFR_->GetBinError  (hPhotonFR_->FindBin(ph.pt(), abs(ph.eta())));
  return FRError;
}

double VVGammaAnalyzer::getPhotonFR_KtoVL(const phys::Photon& ph) const{
  return hPhotonFR_KtoVL_->GetBinContent(hPhotonFR_KtoVL_->FindBin(ph.pt(), abs(ph.eta())));
}

double VVGammaAnalyzer::getPhotonFR_KtoVLexcl(const phys::Photon& ph) const{
  return hPhotonFR_KtoVLexcl_->GetBinContent(hPhotonFR_KtoVLexcl_->FindBin(ph.pt(), abs(ph.eta())));
}

bool VVGammaAnalyzer::passVeryLoose(const Photon& ph){
  Photon::IdWp wp = Photon::IdWp::Loose;
  return ph.cutBasedID(wp, Photon::IDcut::HoverE) && ph.cutBasedID(wp, Photon::IDcut::phIso) && ph.cutBasedID(wp, Photon::IDcut::neIso);
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
						   {0., 0.8, 1.4442, 1.566, 2, 2.5}
						   // {0., 0.435, 0.783, 1.13, 1.4442, 1.566, 1.8, 2.1, 2.5}
						   );

const vector<double> VVGammaAnalyzer::ph_pt_bins(
						 {20., 35., 50., 80., 120}
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
const vector<double> VVGammaAnalyzer::mZ_bins(
					      {60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120}
					      );

const vector<double> VVGammaAnalyzer::mZG_bins(
					       {0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250}
					       );
