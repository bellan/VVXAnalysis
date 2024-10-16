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
#define BINS_PHCUTFLOW 7,-0.5,6.5
#define BINS_KINPHID 5,-0.5,4.5

namespace {
  // Anonymous namespace that holds constants specific to this analyzer
  enum class FSRcutType { MLL_MIN, MLL_IMPROVES, MLLG_MIN };
  constexpr FSRcutType FSR_CUT_TYPE = FSRcutType::MLL_MIN;
  constexpr float CUT_MLL_MIN = 81.;
  constexpr float CUT_PTG_MIN = 20.;
  constexpr float CUT_G_AETA_MAX     = 2.4;
  constexpr float CUT_G_AETA_GAP_MIN = 1.4442;
  constexpr float CUT_G_AETA_GAP_MAX = 1.566;
  constexpr float CUT_MLLG_MIN=100.;
  constexpr bool APPLY_FSR_CUT       = false;
  constexpr bool APPLY_PIXELSEED_CUT = false;
  constexpr bool PHFR_SPLIT          = true;
}

std::pair<TLorentzVector, TLorentzVector> solveNuPz(const Boson<Lepton>& W, int& error);
bool inPhotonEtaAcceptance(double eta);

void VVGammaAnalyzer::begin(){
  cout<<'\n';
  for(char i=0; i<25; ++i) cout<<'-';
  cout<<" Start of VVGammaAnalyzer ";
  for(char i=0; i<25; ++i) cout<<'-';
  cout<<'\n';

  const char* CMSSW_BASE = getenv("CMSSW_BASE");
  std::string VVXAnalysis_dir = CMSSW_BASE == nullptr ? ".." : Form("%s/src/VVXAnalysis", CMSSW_BASE);

  std::string year_str;
  if(PHFR_SPLIT){
    year_str = std::to_string(year);
    if(year == 2016)
      year_str += subEra_;
  }
  else{
    year_str = "Run2";
  }
  
  // Photon FR
  hPhotonFR_VLtoL_data_   = getHistfromFile(Form("data/FR_VLtoL_pt-aeta_data_%s.root"        , year_str.c_str()), "PhFR", " VLtoL (data)"   );

  // FR extended
  hPhotonFR_VLtoL_dataZG_ = getHistfromFile(Form("data/FR_VLtoL_pt-aeta_data-ZGToLLG_%s.root", year_str.c_str()), "PhFR", " VLtoL (data-ZG)");
  hPhotonFR_KtoVLexcl_    = getHistfromFile(Form("data/FR_KtoVLexcl_pt-aeta_data-ZGToLLG_%s.root", year_str.c_str()), "PhFR", " KtoVLexcl (data-ZG)");

  // FR SF
  hPhotonFRSF_VLtoL_      = getHistfromFile(Form("data/ratio_VLtoL_pt-aeta_data_over_ZZ_%s.root", year_str.c_str()), "PhFRSF");

  // Photon efficiency SF for cut-based ID (temporary)
  std::string pathPhotonEffSF(Form("%s/Commons/data/egammaEffi.txt_EGM2D_Pho_Loose_UL%d%s.root", VVXAnalysis_dir.c_str(), year%100, (subEra_.size() > 0 ? ('_'+subEra_).c_str() : "")));
  hPhotonEffSF_           = getHistfromFile(pathPhotonEffSF.c_str(), "EGamma_SF2D");
  hPhotonEffSF_maxPt_     = hPhotonEffSF_->GetYaxis()->GetBinUpEdge(hPhotonEffSF_->GetNbinsY());

  // Photon MVA SF
  mapPhotonMVASF_[Photon::MVAwp::wp80] = getHistfromFile(Form("%s/Commons/data/%d_PhotonsMVAwp80.root", VVXAnalysis_dir.c_str(), year), "EGamma_SF2D");
  mapPhotonMVASF_[Photon::MVAwp::wp90] = getHistfromFile(Form("%s/Commons/data/%d_PhotonsMVAwp90.root", VVXAnalysis_dir.c_str(), year), "EGamma_SF2D");
  for(auto& it: mapPhotonMVASF_)  // std::pair<const Photon::MVAwp, std::unique_ptr<TH2F>>
    mapPhotonMVASF_maxPt_[it.first] = it.second->GetYaxis()->GetBinUpEdge(it.second->GetNbinsY());

  // initCherryPick();
  for(const char* sys : photonSystKeys_){
    kinPhotons_ [sys]  = std::make_unique<vector<Photon>>();
    loosePhotons_[sys] = std::make_unique<vector<Photon>>();
    goodPhotons_[sys]  = std::make_unique<vector<Photon>>();
  }

#ifndef DEBUG
  size_t digits = std::to_string(tree()->GetEntries()).length();
  std::string spaces( digits, ' ' );
  cout<<"Analyzed:\t"<<spaces<<'/'<<tree()->GetEntries()<<std::flush;
#else
  cout<<"\tSizes: Particle="<<sizeof(Particle)<<" , Lepton="<<sizeof(Lepton)<<" , Photon="<<sizeof(Photon)<<" , Jet="<<sizeof(Jet)<<" , Boson<Particle>="<<sizeof(Boson<Particle>)<<" , DiBoson<Particle,Particle>="<<sizeof(DiBoson<Particle,Particle>)<<'\n';
#endif

  return;
}


void VVGammaAnalyzer::initEvent(){
  // Cleanup
  leptons_    ->clear();
  jets_noph_  ->clear();
  fsrPhotons_ ->clear();
  for(auto const& it: kinPhotons_ )  it.second->clear();
  for(auto const& it: loosePhotons_) it.second->clear();
  for(auto const& it: goodPhotons_)  it.second->clear();
  candVTojj_ = Boson<Jet>();
  candVToJ_ = Jet();

  // Contruct vector with leptons from dibosons and FSR photons
  if     (is4Lregion(region_) || (region_ == MC && ZZ && ZZ->pt() > 0.001)){
    leptons_->insert(leptons_->end(), {
    	  ZZ->first().daughter(0), 
    	  ZZ->first().daughter(1), 
    	  ZZ->second().daughter(0), 
    	  ZZ->second().daughter(1)
        });

    for(auto bos : {ZZ->first(), ZZ->second()}){
      std::bitset<2> fsrIndex = std::bitset<2>(bos.daughtersWithFSR());
      if(fsrIndex.test(0)) fsrPhotons_->push_back(bos.fsrPhoton(0));
      if(fsrIndex.test(1)) fsrPhotons_->push_back(bos.fsrPhoton(1));
    }

    theHistograms->fill("Z0_dRll", "leptons of Z_{0};#DeltaR(l_{0},l_{1});Events", 40,0.,1., deltaR(ZZ->first() .daughter(0), ZZ->first() .daughter(1)), theWeight);
    theHistograms->fill("Z1_dRll", "leptons of Z_{1};#DeltaR(l_{0},l_{1});Events", 40,0.,1., deltaR(ZZ->second().daughter(0), ZZ->second().daughter(1)), theWeight);
  }
  else if(is3Lregion(region_) || (region_ == MC && ZW && ZW->pt() > 0.001)){
    leptons_->insert(leptons_->end(), {
  	  ZW->first().daughter(0), 
  	  ZW->first().daughter(1), 
  	  ZW->second().daughter(0)
	});

    for(auto bos : {ZW->first(), ZW->second()}){
      std::bitset<2> fsrIndex = std::bitset<2>(bos.daughtersWithFSR());
      if(fsrIndex.test(0)) fsrPhotons_->push_back(bos.fsrPhoton(0));
      if(fsrIndex.test(1)) fsrPhotons_->push_back(bos.fsrPhoton(1));
    }

    theHistograms->fill("Z_dRll", "leptons of Z;#DeltaR(l_{0},l_{1});Events", 40,0.,1., deltaR(ZW->first().daughter(0), ZW->first().daughter(1)), theWeight);
  }
  else if(region_ == CRLFR || (region_ == MC && ZL && ZL->first.pt() > 0.001)){
    leptons_->insert(leptons_->end(), {
  	  ZL->first.daughter(0),
  	  ZL->first.daughter(1),
  	  ZL->second
	});

    const Boson<Lepton>& bos = ZL->first;
    std::bitset<2> fsrIndex = std::bitset<2>(bos.daughtersWithFSR());
    if(fsrIndex.test(0)) fsrPhotons_->push_back(bos.fsrPhoton(0));
    if(fsrIndex.test(1)) fsrPhotons_->push_back(bos.fsrPhoton(1));

    for(auto ph : *photons){
      if(canBeFSR(ph, std::vector<Lepton> {ZL->second}))
	fsrPhotons_->push_back(ph);
    }

    theHistograms->fill("Z_dRll", "leptons of Z;#DeltaR(l_{0},l_{1});Events", 40,0.,1., deltaR(ZL->first.daughter(0), ZL->first.daughter(1)), theWeight);
  }
  else if(is2Lregion(region_) || region_ == MC){
    leptons_->insert(leptons_->cend(), electrons->begin(), electrons->end());
    leptons_->insert(leptons_->cend(),     muons->begin(),     muons->end());

    if(leptons_->size() > 2){
      throw std::runtime_error(Form("%lu leptons in 2L region %s", leptons_->size(), regionType(region_).c_str()));
    }

    if(!(Z && Z->pt() > 0.001) && leptons_->size() == 2)
      *Z = Boson<Lepton>(leptons_->at(0), leptons_->at(1));

    // Correct method for when we'll have ntuples with ZCand filled correctly
    // const Boson<Lepton>& bos = *Z;
    // std::bitset<2> fsrIndex = std::bitset<2>(bos.daughtersWithFSR());
    // if(fsrIndex.test(0)) fsrPhotons_->push_back(bos.fsrPhoton(0));
    // if(fsrIndex.test(1)) fsrPhotons_->push_back(bos.fsrPhoton(1));

    // By hand for now
    for(auto ph : *photons){
      if(canBeFSR(ph, *leptons_))
	fsrPhotons_->push_back(ph);
    }
  }

  
  // Photon selection
  unsigned int nKinPh_0p07 (0), nVLPh_0p07 (0), nFailPh_0p07 (0), nLoosePh_0p07 (0);
  unsigned int nKinPh_0p3  (0), nVLPh_0p3  (0), nFailPh_0p3  (0), nLoosePh_0p3  (0);
  unsigned int nKinPh_0p5  (0), nVLPh_0p5  (0), nFailPh_0p5  (0), nLoosePh_0p5  (0);
  unsigned int nKinPh_0p7  (0), nVLPh_0p7  (0), nFailPh_0p7  (0), nLoosePh_0p7  (0);
  unsigned int nKinPh_noFSR(0), nVLPh_noFSR(0), nFailPh_noFSR(0), nLoosePh_noFSR(0);

  // DEBUG
  vector<Photon> fsrMatched; fsrMatched.reserve(fsrPhotons_->size());
  vector<Photon> kinPhotonsWithFSR;  // do not veto FSR photons

  for(auto ph : *photons){
    //Pixel seed and electron veto
    if(!ph.passElectronVeto()) continue;
    if(APPLY_PIXELSEED_CUT && ph.hasPixelSeed()) continue;

    //Kinematic selection
    // if(ph.pt() < 20) continue;
    float ph_aeta = fabs(ph.eta());
    if(!inPhotonEtaAcceptance(ph_aeta)) continue;

    // Check ID
    bool isPassVL    = ph.cutBasedID(Photon::IdWp::VeryLoose);
    bool isPassLoose = ph.cutBasedIDLoose();

    auto closestLep = closestDeltaR(ph, *leptons_);
    float minDR_lep = closestLep != leptons_->cend() ? deltaR(ph, *closestLep) : 10.;
    if(minDR_lep > 0.07 && ph.pt() > CUT_PTG_MIN){
      if(true)        ++nKinPh_0p07;
      if(isPassVL)    ++nVLPh_0p07;
      if(isPassVL && !isPassLoose) ++nFailPh_0p07;
      if(isPassLoose) ++nLoosePh_0p07;
    }
    if(minDR_lep > 0.3  && ph.pt() > CUT_PTG_MIN){
      if(true)        ++nKinPh_0p3;
      if(isPassVL)    ++nVLPh_0p3;
      if(isPassVL && !isPassLoose) ++nFailPh_0p3;
      if(isPassLoose) ++nLoosePh_0p3;
    }
    if(minDR_lep > 0.5  && ph.pt() > CUT_PTG_MIN){
      if(true)        ++nKinPh_0p5;
      if(isPassVL)    ++nVLPh_0p5;
      if(isPassVL && !isPassLoose) ++nFailPh_0p5;
      if(isPassLoose) ++nLoosePh_0p5;
    }
    if(minDR_lep > 0.7  && ph.pt() > CUT_PTG_MIN){
      if(true)        ++nKinPh_0p7;
      if(isPassVL)    ++nVLPh_0p7;
      if(isPassVL && !isPassLoose) ++nFailPh_0p7;
      if(isPassLoose) ++nLoosePh_0p7;
    }

    if(ph.pt() > CUT_PTG_MIN)
      kinPhotonsWithFSR.push_back(ph);
    // Remove photons that were used for FSR
    auto closestFSR = closestDeltaR(ph, *fsrPhotons_);
    if(closestFSR != fsrPhotons_->cend() && deltaR(ph, *closestFSR) < 1e-3){
      if(ph.pt() > CUT_PTG_MIN)
	fsrMatched.push_back(ph);
      continue;
    }

    if(ph.pt() > CUT_PTG_MIN){
      if(true)        ++nKinPh_noFSR;
      if(isPassVL)    ++nVLPh_noFSR;
      if(isPassLoose) ++nLoosePh_noFSR;
      if(isPassVL && !isPassLoose) ++nFailPh_noFSR;
    }

    // CUT: dRl > 0.5
    if(minDR_lep < 0.5)
      continue;

    TLorentzVector p4_EScale_Up = ph.p4() * (ph.energyScaleUp()  /ph.e());
    TLorentzVector p4_EScale_Dn = ph.p4() * (ph.energyScaleDown()/ph.e());
    TLorentzVector p4_ESigma_Up = ph.p4() * (ph.energySigmaUp()  /ph.e());
    TLorentzVector p4_ESigma_Dn = ph.p4() * (ph.energySigmaDown()/ph.e());
    
    if(p4_EScale_Up.Pt() > CUT_PTG_MIN){
      Photon copy(ph);
      copy.setP4(p4_EScale_Up);
      kinPhotons_["EScale_Up"]->push_back(copy);
      if(isPassVL)
	loosePhotons_["EScale_Up"]->push_back(copy);
      if(isPassLoose)
	goodPhotons_["EScale_Up"]->push_back(std::move(copy));
    }
    if(p4_EScale_Dn.Pt() > CUT_PTG_MIN){
      Photon copy(ph);
      copy.setP4(p4_EScale_Dn);
      kinPhotons_["EScale_Down"]->push_back(copy);
      if(isPassVL)
    	loosePhotons_["EScale_Down"]->push_back(copy);
      if(isPassLoose)
    	goodPhotons_["EScale_Down"]->push_back(std::move(copy));
    }
    if(p4_ESigma_Up.Pt() > CUT_PTG_MIN){
      Photon copy(ph);
      copy.setP4(p4_ESigma_Up);
      kinPhotons_["ESigma_Up"]->push_back(copy);
      if(isPassVL)
        loosePhotons_["ESigma_Up"]->push_back(copy);
      if(isPassLoose)
    	goodPhotons_["ESigma_Up"]->push_back(std::move(copy));
    }
    if(p4_ESigma_Dn.Pt() > CUT_PTG_MIN){
      Photon copy(ph);
      copy.setP4(p4_ESigma_Dn);
      kinPhotons_["ESigma_Down"]->push_back(copy);
      if(isPassVL)
    	loosePhotons_["ESigma_Down"]->push_back(copy);
      if(isPassLoose)
    	goodPhotons_["ESigma_Down"]->push_back(std::move(copy));
    }
    if(ph.pt() > CUT_PTG_MIN){
      kinPhotons_["central"]->push_back(ph);
      if(isPassVL)
	loosePhotons_["central"]->push_back(ph);
      if(isPassLoose)
	goodPhotons_["central"]->push_back(ph);
    }
    
  } // end loop on *photons

  // DEBUG
  theHistograms->fill("kinPhotonsWithFSR_N"  , "Kin photons including FSR;N #gamma;Events", 5,-0.5,4.5, kinPhotonsWithFSR.size(), theWeight);
  if(kinPhotonsWithFSR.size() > 0){
    const Photon& ph = kinPhotonsWithFSR.at(0);
    fillPhotonPlots(ph, "lead_kinPhotonsWithFSR", "Leading Kin including FSR");
  }
  theHistograms->fill("fsrMatched_N"         , "FSR matched photons;N #gamma;Events", 5,-0.5,4.5, fsrMatched.size(), theWeight);
  if(fsrMatched.size() > 0){
    const Photon& ph = fsrMatched.at(0);
    fillPhotonPlots(ph, "lead_fsrMatched", "Leading FSR matched");
  }

  // plots on surviving photons
  if(true)           theHistograms->fill("cuts_kinPhIso"      , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, "All"  , theWeight);
  if(nKinPh_0p07)    theHistograms->fill("cuts_kinPhIso"      , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, ">0.07", theWeight);
  if(nKinPh_noFSR)   theHistograms->fill("cuts_kinPhIso"      , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, "noFSR", theWeight);
  if(nKinPh_0p3)     theHistograms->fill("cuts_kinPhIso"      , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, ">0.3" , theWeight);
  if(nKinPh_0p5)     theHistograms->fill("cuts_kinPhIso"      , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, ">0.5" , theWeight);
  if(nKinPh_0p7)     theHistograms->fill("cuts_kinPhIso"      , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, ">0.7" , theWeight);

  if(true)           theHistograms->fill("cuts_veryLoosePhIso", ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, "All"  , theWeight);
  if(nVLPh_0p07)     theHistograms->fill("cuts_veryLoosePhIso", ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, ">0.07", theWeight);
  if(nVLPh_noFSR)    theHistograms->fill("cuts_veryLoosePhIso", ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, "noFSR", theWeight);
  if(nVLPh_0p3)      theHistograms->fill("cuts_veryLoosePhIso", ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, ">0.3" , theWeight);
  if(nVLPh_0p5)      theHistograms->fill("cuts_veryLoosePhIso", ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, ">0.5" , theWeight);
  if(nVLPh_0p7)      theHistograms->fill("cuts_veryLoosePhIso", ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, ">0.7" , theWeight);

  if(true)           theHistograms->fill("cuts_failPhIso"     , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, "All"  , theWeight);
  if(nFailPh_0p07)   theHistograms->fill("cuts_failPhIso"     , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, ">0.07", theWeight);
  if(nFailPh_noFSR)  theHistograms->fill("cuts_failPhIso"     , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, "noFSR", theWeight);
  if(nFailPh_0p3)    theHistograms->fill("cuts_failPhIso"     , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, ">0.3" , theWeight);
  if(nFailPh_0p5)    theHistograms->fill("cuts_failPhIso"     , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, ">0.5" , theWeight);
  if(nFailPh_0p7)    theHistograms->fill("cuts_failPhIso"     , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, ">0.7" , theWeight);

  if(true)           theHistograms->fill("cuts_loosePhIso"    , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, "All"  , theWeight);
  if(nLoosePh_0p07)  theHistograms->fill("cuts_loosePhIso"    , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, ">0.07", theWeight);
  if(nLoosePh_noFSR) theHistograms->fill("cuts_loosePhIso"    , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, "noFSR", theWeight);
  if(nLoosePh_0p3)   theHistograms->fill("cuts_loosePhIso"    , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, ">0.3" , theWeight);
  if(nLoosePh_0p5)   theHistograms->fill("cuts_loosePhIso"    , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, ">0.5" , theWeight);
  if(nLoosePh_0p7)   theHistograms->fill("cuts_loosePhIso"    , ";cut;Events", {"All", ">0.07", "noFSR", ">0.3", ">0.5", ">0.7"}, ">0.7" , theWeight);

  studyFSRregion(fsrMatched);

  if(kinPhotons_["central"]->size() > 0){
    bestKinPh_ = &*std::max_element(kinPhotons_["central"]->begin(), kinPhotons_["central"]->end(),
				    [](const Photon& a, const Photon& b){ return a.nCutsPass(Photon::IdWp::Loose) < b.nCutsPass(Photon::IdWp::Loose); }
				    );  // max_element returns the first among those with max value --> preserve pt ordering
    bestMVAPh_ = &*std::max_element(kinPhotons_["central"]->begin(), kinPhotons_["central"]->end(),
                                    [](const Photon& a, const Photon& b){
                                        int score_a = 2*(a.passMVA(Photon::MVAwp::wp80)) + (a.passMVA(Photon::MVAwp::wp90));
                                        int score_b = 2*(b.passMVA(Photon::MVAwp::wp80)) + (b.passMVA(Photon::MVAwp::wp90));
                                        if(score_a == score_b)
                                          return a.MVAvalue() < b.MVAvalue();
                                        else
                                         return score_a < score_b;
                                      }
				    );
  }
  else{
    bestKinPh_ = nullptr;
    bestMVAPh_ = nullptr;
  }

  // Decide channel name depending on leptons
  makeChannelReco();

  // Initialize the collection of jets disambiguated from photons with deltaR = 0.2
  vector<Photon>& phvect = *kinPhotons_["central"];
  if(phvect.size() == 0)
    std::copy(   jets->begin(), jets->end(), std::back_inserter(*jets_noph_));
  else
    std::copy_if(jets->begin(), jets->end(), std::back_inserter(*jets_noph_),
		 [&phvect](const Jet& j){
		   return deltaR(j, *closestDeltaR(j, phvect)) > 0.2;
		 });
}


Int_t VVGammaAnalyzer::cut() {
  evtN_++; evtNInReg_[region_]++; evtWInReg_[region_] += theWeight;
  #ifndef DEBUG
  // cout<<"\r\t\t"<<evtN_;
  #else
  // cout << regionType(region_) << ':' << run << ':' << lumiBlock << ':' << event << '\t' << theWeight << '\n';
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

  theHistograms->fill("AAA_cuts"  , "cuts weighted"  , {"All"}, "All", theWeight);
  theHistograms->fill("AAA_cuts_u", "cuts unweighted", {"All"}, "All", 1);
  
  if(theSampleInfo.isMC()){
    genEventSetup();
    genEventHistos();
    // ----- SIGNAL DEFINITION -----
    sigdefHelper.eval();
      if(sigdefHelper.pass()) theHistograms->fill("AAA_cuts_sigdef",";;Events", {"All"}, "All", theWeight);
  }
  initEvent();
  
  theHistograms->fill("POG_leptons", "n leptons", 7,-0.5,6.5, electrons->size()+muons->size(), theWeight);  
  
  baseHistos_cut();
  photonHistos();
  jetHistos();
  // PKU_comparison();

  std::vector<Particle> genPhotonsKin;
  std::vector<Particle> genPhotonsKinPrompt;
  std::vector<Particle> genPhotonsDRl;
  std::vector<Particle> genPhotonsKinDRl;
  for(const Particle& p : *genPhotons_){
    double aeta = fabs(p.eta());
    bool isPrompt = p.genStatusFlags().test(phys::isPrompt);
    bool passKin  = p.pt() > CUT_PTG_MIN && inPhotonEtaAcceptance(aeta);
    bool passDRl  = physmath::deltaR(p, *std::min_element(genChLeptons_->begin(), genChLeptons_->end(), DeltaRComparator(p))) > 0.5;
    if(passKin            ) genPhotonsKin      .push_back(p);
    if(passKin && isPrompt) genPhotonsKinPrompt.push_back(p);
    if(passDRl            ) genPhotonsDRl      .push_back(p);
    if(passKin && passDRl ) genPhotonsKinDRl   .push_back(p);
  }

  photonIsolation(*    photons            , "all" );
  photonIsolation(* kinPhotons_["central"], "kin" );
  photonIsolation(*goodPhotons_["central"], "good");

  efficiency(*    photons            , *genPhotons_, "photons"    , "gen", 0.2);
  efficiency(* kinPhotons_["central"], *genPhotons_, "kinPhotons" , "gen", 0.2);
  efficiency(*goodPhotons_["central"], *genPhotons_, "goodPhotons", "gen", 0.2);
  efficiency(*goodPhotons_["central"],  genPhotonsKin      , "goodPhotons", "genKin"      , 0.2);
  efficiency(*goodPhotons_["central"],  genPhotonsDRl      , "goodPhotons", "genDRl"      , 0.2);
  efficiency(*goodPhotons_["central"], *genPhotonsPrompt_  , "goodPhotons", "genPrompt"   , 0.2);
  efficiency(*goodPhotons_["central"],  genPhotonsKinPrompt, "goodPhotons", "genKinPrompt", 0.2);
  efficiency(*goodPhotons_["central"],  genPhotonsKinDRl   , "goodPhotons", "genKinDRl"   , 0.2);

  efficiency(*jets, *genJets, "AK4", "genJets", 0.4);


  // ----- BASELINE SELECTION -----
  vector<Photon>& phVect_CUT_mllimprov = *kinPhotons_["central"];  // Photon vector used for the "Improves mll" CUT
  bool passFSRcut = false;

  // -----  4L   -----
  if     (is4Lregion(region_)){
    // Cut 4L.ZZ: require a ZZ candidate
    bool haveZZlep = (ZZ && ZZ->pt() > 0.001);
    if(haveZZlep){
      theHistograms->fill("AAA_cuts"  , "cuts weighted"  , {}, "ZZ", theWeight);
      theHistograms->fill("AAA_cuts_u", "cuts unweighted", {}, "ZZ", 1);
      if(sigdefHelper.pass()) theHistograms->fill("AAA_cuts_sigdef", "", {}, "ZZ", theWeight);
    }
    else return -1;

    // Cut 4L.FSR
    if(FSR_CUT_TYPE == FSRcutType::MLL_MIN){
      // Cut 4L.FSR.mll_min: Require that m_{ll} is greater than threshold
      passFSRcut = ZZ->first().mass() > CUT_MLL_MIN && ZZ->second().mass() > CUT_MLL_MIN;
    }
    else{
      if(phVect_CUT_mllimprov.size() > 0){
	double Zll_mass(0.), ZllG_mass(0.);
	std::tie(Zll_mass, ZllG_mass) = getZllAndZllgMasses(phVect_CUT_mllimprov);

	if     (FSR_CUT_TYPE == FSRcutType::MLLG_MIN)
	  // Cut 4L.FSR.mllg_min: Require that m_{llg} is greater than threshold
	  passFSRcut = ZllG_mass > CUT_MLLG_MIN;
	else if(FSR_CUT_TYPE == FSRcutType::MLL_IMPROVES)
	  // Cut 4L.FSR.mll_noimprov: Require that neither of the Z masses improves with the addition of the momentum of a photon
	  passFSRcut = (Zll_mass > 0) && (Zll_mass < CUT_MLL_MIN) && ( fabs(ZllG_mass - phys::ZMASS) < fabs(Zll_mass - phys::ZMASS) );
      }
    }
  }

  // -----  3L   -----
  else if(is3Lregion(region_)){
    // Cut 3L.WZ: require a WZ candidate
    bool haveWZlep = (ZW && ZW->pt() > 0.001);
    if(haveWZlep){
      theHistograms->fill("AAA_cuts"  , "cuts weighted"  , {}, "WZ", theWeight);
      theHistograms->fill("AAA_cuts_u", "cuts unweighted", {}, "WZ", 1);
      if(sigdefHelper.pass()) theHistograms->fill("AAA_cuts_sigdef", "", {}, "WZ", theWeight);
    }
    else return -1;

    // Cut 3L.paperSel: replicate the selection from CMS-SMP-20-014
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
    
    if(b_WZpaperSel){
      theHistograms->fill("AAA_cuts"  , "cuts weighted"  , {}, "paperSel", theWeight);
      theHistograms->fill("AAA_cuts_u", "cuts unweighted", {}, "paperSel", 1);
      if(sigdefHelper.pass()) theHistograms->fill("AAA_cuts_sigdef", "", {}, "paperSel", theWeight);
    }
    else return -1;

    // Cut 3L.FSR
    if(FSR_CUT_TYPE == FSRcutType::MLL_MIN){
      // Cut 3L.FSR.mll_min: Require that m_{ll} is greater than threshold
      passFSRcut = ZW->first().mass() > CUT_MLL_MIN;
    }
    else{
      auto closestPhoLep = closestPairDeltaR(phVect_CUT_mllimprov, *leptons_);
      auto best_lep_index = std::distance(leptons_->cbegin(), closestPhoLep.second);
      if(best_lep_index <= 2){  // The lepton belongs to the Z
	double Zll_mass  = ZW->first().mass();
	double ZllG_mass = (ZW->first().p4() + closestPhoLep.first->p4()).M();

	if     (FSR_CUT_TYPE == FSRcutType::MLLG_MIN)
	  // Cut 3L.FSR.mllg_min: Require that m_{llg} is greater than threshold
	  passFSRcut = ZllG_mass > CUT_MLLG_MIN;
	else if(FSR_CUT_TYPE == FSRcutType::MLL_IMPROVES)
	  // Cut 3L.FSR.mll_noimprov: Require that neither of the Z masses improves with the addition of the momentum of a photon
	  passFSRcut = (Zll_mass < CUT_MLL_MIN) && fabs(ZllG_mass - phys::ZMASS) < fabs(Zll_mass - phys::ZMASS);
      }
    }
  }

  // -----  2L   -----
  else if(is2Lregion(region_)){
    // Cut 2L.Zlep: require a Z lep candidate
    bool haveZlep = (Z && Z->pt() > 0.001);
    if(haveZlep){
      theHistograms->fill("AAA_cuts"  , "cuts weighted"  , {}, "Z lep", theWeight);
      theHistograms->fill("AAA_cuts_u", "cuts unweighted", {}, "Z lep", 1);
      if(sigdefHelper.pass()) theHistograms->fill("AAA_cuts_sigdef", "", {}, "Z lep", theWeight);
    }
    else return -1;

    // Cut 2L.Jets: require two AK4 or an AK8
    bool haveJets = (jets->size() >= 2 || jetsAK8->size() >= 1);
    if(haveJets){
      theHistograms->fill("AAA_cuts"  , "cuts weighted"  , {}, "2l2j || 2l1J", theWeight);
      theHistograms->fill("AAA_cuts_u", "cuts unweighted", {}, "2l2j || 2l1J", 1);
      if(sigdefHelper.pass()) theHistograms->fill("AAA_cuts_sigdef", "", {}, "2l2j || 2l1J", theWeight);
    }
    else return -1;

    // Cut 2L.FSR
    if(FSR_CUT_TYPE == FSRcutType::MLL_MIN){
      // Cut 2L.FSR.mll_min: Require that m_{ll} is greater than threshold
      passFSRcut = Z->mass() > CUT_MLL_MIN;
    }
    else{
      if(phVect_CUT_mllimprov.size() > 0){
	const Lepton& l0 = Z->daughter(0);
	const Lepton& l1 = Z->daughter(1);
	auto ph0 = closestDeltaR(l0, phVect_CUT_mllimprov);
	auto ph1 = closestDeltaR(l1, phVect_CUT_mllimprov);
	auto thePh = deltaR(l0, *ph0) < deltaR(l1, *ph1) ? ph0 : ph1;
	double Zll_mass = Z->mass();
	double ZllG_mass = (Z->p4() + thePh->p4()).M();

	if     (FSR_CUT_TYPE == FSRcutType::MLLG_MIN)
	  // Cut 2L.FSR.mllg_min: Require that m_{llg} is greater than threshold
	  passFSRcut = ZllG_mass > CUT_MLLG_MIN;
	else if(FSR_CUT_TYPE == FSRcutType::MLL_IMPROVES)
	  // Cut 2L.FSR.mll_noimprov: Require that neither of the Z masses improves with the addition of the momentum of a photon
	  passFSRcut = (Zll_mass < CUT_MLL_MIN) && fabs(ZllG_mass - phys::ZMASS) < fabs(Zll_mass - phys::ZMASS);
      }
    }
  }

  // ----- CRLFR -----
  else if(region_ == CRLFR){
    // Cut CRLFR.ZL: require a ZL candidate
    bool haveZL = (ZL && ZL->first.pt()+ZL->second.pt() > 0.001);
    if(haveZL){
      theHistograms->fill("AAA_cuts"  , "cuts weighted"  , {}, "ZL", theWeight);
      theHistograms->fill("AAA_cuts_u", "cuts unweighted", {}, "ZL", 1);
      if(sigdefHelper.pass()) theHistograms->fill("AAA_cuts_sigdef", "", {}, "ZL", theWeight);
    }
    else return -1;

    // Cut CRLFR.FSR
    if(FSR_CUT_TYPE == FSRcutType::MLL_MIN){
      // Cut CRLFR.FSR.mll_min: Require that m_{ll} is greater than threshold
      passFSRcut = ZL->first.mass() > CUT_MLL_MIN;
    }
    else{
      if(phVect_CUT_mllimprov.size() > 0){
	const Lepton& l0 = ZL->first.daughter(0);
	const Lepton& l1 = ZL->first.daughter(1);
	auto ph0 = closestDeltaR(l0, phVect_CUT_mllimprov);
	auto ph1 = closestDeltaR(l1, phVect_CUT_mllimprov);
	auto thePh = deltaR(l0, *ph0) < deltaR(l1, *ph1) ? ph0 : ph1;
	double Zll_mass = ZL->first.mass();
	double ZllG_mass = (ZL->first.p4() + thePh->p4()).M();

	if     (FSR_CUT_TYPE == FSRcutType::MLLG_MIN)
	  // Cut CRLFR.FSR.mllg_min: Require that m_{llg} is greater than threshold
	  passFSRcut = ZllG_mass > CUT_MLLG_MIN;
	else if(FSR_CUT_TYPE == FSRcutType::MLL_IMPROVES)
	  // Cut CRLFR.FSR.mll_noimprov: Require that neither of the Z masses improves with the addition of the momentum of a photon
	  passFSRcut = (Zll_mass < CUT_MLL_MIN) && fabs(ZllG_mass - phys::ZMASS) < fabs(Zll_mass - phys::ZMASS);
      }
    }
  }

  // ----- Common cuts -----
  // Cut: common.FSR
  if(passFSRcut){
    theHistograms->fill("AAA_cuts"  , "cuts weighted"  , {}, "FSR cut", theWeight);
    theHistograms->fill("AAA_cuts_u", "cuts unweighted", {}, "FSR cut", 1);
    if(sigdefHelper.pass()) theHistograms->fill("AAA_cuts_sigdef", "", {}, "FSR cut", theWeight);
  }
  else if(APPLY_FSR_CUT) return -1;

  return 1;
}


void VVGammaAnalyzer::analyze(){
  analyzedNInReg_[region_]++; analyzedWInReg_[region_] += theWeight;
  theHistograms->fill("AAA_cuts"  , "cuts weighted"  , {}, "Analyzed", theWeight);
  theHistograms->fill("AAA_cuts_u", "cuts unweighted", {}, "Analyzed", 1);
  if(sigdefHelper.pass()) theHistograms->fill("AAA_cuts_sigdef", "", {}, "Analyzed", theWeight);

  // Disabled cuts
  // disabled.gamma_kin: Require at least 1 kin photon with pt > 20 GeV
  bool haveKinPhoton = kinPhotons_["central"]->size() >= 1;
  if(haveKinPhoton){
    theHistograms->fill("AAA_cuts"  , "cuts weighted"  , {}, "#gamma kin", theWeight);
    theHistograms->fill("AAA_cuts_u", "cuts unweighted", {}, "#gamma kin", 1);
    if(sigdefHelper.pass()) theHistograms->fill("AAA_cuts_sigdef", "", {}, "#gamma kin", theWeight);
  }

  // disabled.gamma_veryloose: Require at least 1 veryLoose photon with pt > 20 GeV
  bool haveVeryLoosePhoton = loosePhotons_["central"]->size() >= 1;
  if(haveVeryLoosePhoton){
    theHistograms->fill("AAA_cuts"  , "cuts weighted"  , {}, "#gamma veryLoose", theWeight);
    theHistograms->fill("AAA_cuts_u", "cuts unweighted", {}, "#gamma veryLoose", 1);
    if(sigdefHelper.pass()) theHistograms->fill("AAA_cuts_sigdef", "", {}, "#gamma VeryLoose", theWeight);
  }

  // disabled.gamma_loose: Require at least 1 loose photon with pt > 20 GeV
  bool haveLoosePhoton = goodPhotons_["central"]->size() >= 1;
  if(haveLoosePhoton){
    theHistograms->fill("AAA_cuts"  , "cuts weighted"  , {}, "#gamma loose", theWeight);
    theHistograms->fill("AAA_cuts_u", "cuts unweighted", {}, "#gamma loose", 1);
    if(sigdefHelper.pass()) theHistograms->fill("AAA_cuts_sigdef", "", {}, "#gamma loose", theWeight);
  }

  if( std::any_of(goodPhotons_["central"]->cbegin(),
		  goodPhotons_["central"]->cend(), 
		  [](const Photon& ph){return ph.cutBasedIDMedium();}) )
  {
    theHistograms->fill("AAA_cuts"  , "cuts weighted"  , {}, "#gamma medium", theWeight);
    theHistograms->fill("AAA_cuts_u", "cuts unweighted", {}, "#gamma medium", 1);
    if(sigdefHelper.pass()) theHistograms->fill("AAA_cuts_sigdef", "", {}, "#gamma medium", theWeight);
  }

  bool havewp90Photon = bestMVAPh_ && bestMVAPh_->passMVA(Photon::MVAwp::wp90);
  if(havewp90Photon){
    theHistograms->fill("AAA_cuts_extra"  , "cuts weighted"  , {"#gamma wp90"}, "#gamma wp90", theWeight);
    theHistograms->fill("AAA_cuts_u_extra", "cuts unweighted", {"#gamma wp90"}, "#gamma wp90", 1);
    if(sigdefHelper.pass()) theHistograms->fill("AAA_cuts_sigdef_extra", "", {"#gamma wp90"}, "#gamma wp90", theWeight);
  }

  bool havewp80Photon = bestMVAPh_ && bestMVAPh_->passMVA(Photon::MVAwp::wp80);
  if(havewp80Photon){
    theHistograms->fill("AAA_cuts_extra"  , "cuts weighted"  , {}, "#gamma wp80", theWeight);
    theHistograms->fill("AAA_cuts_u_extra", "cuts unweighted", {}, "#gamma wp80", 1);
    if(sigdefHelper.pass()) theHistograms->fill("AAA_cuts_sigdef_extra", "", {}, "#gamma wp80", theWeight);
  }

  fillCutFlow("AAA_cuts_genreco_cutID", ";cut;Events", {
      {"#gamma_{GEN}"      , sigdefHelper.pass_photon()},
      {"full sig def"      , sigdefHelper.pass()},
      {"#gamma_{REC} kin"  , haveKinPhoton},
      {"#gamma_{REC} VL"   , haveVeryLoosePhoton},
      {"#gamma_{REC} Loose", haveLoosePhoton}
    }, theWeight);
  fillCutFlow("AAA_cuts_genreco_MVAID", ";cut;Events", {
      {"#gamma_{GEN}"      , sigdefHelper.pass_photon()},
      {"full sig def"      , sigdefHelper.pass()},
      {"#gamma_{REC} kin"  , haveKinPhoton},
      {"#gamma_{REC} wp90" , havewp90Photon},
      {"#gamma_{REC} wp80" , havewp80Photon}
    }, theWeight);

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
  
  // leptonFakeRate();
  photonGenStudy();
  photonIsolation_bestKin();
  photonFakeRate_ABCD();

  // Photon fake rate loose to tight with various working points and closure tests
  if(bestKinPh_){
    bool isPassVL    = bestKinPh_->cutBasedID(Photon::IdWp::VeryLoose);
    bool isPassLoose = bestKinPh_->cutBasedIDLoose();
    double cutEffSF = isPassLoose ? getPhotonEffSF(*bestKinPh_) : 1.;

    if(isPassVL){
      photonFakeRate_LtoT("VLtoL", *bestKinPh_, isPassLoose, cutEffSF);

      double f_VLtoL_data = getPhotonFR_VLtoL_data(*bestKinPh_);
      photonFRClosure("VLtoL_pt-aeta_data"  , *bestKinPh_, isPassLoose, f_VLtoL_data  );

      double f_VLtoL_dataZG = getPhotonFR_VLtoL_dataZG(*bestKinPh_);
      photonFRClosure("VLtoL_pt-aeta_dataZG", *bestKinPh_, isPassLoose, f_VLtoL_dataZG);
    }

    if(!isPassLoose){
      photonFakeRate_LtoT("KtoVLexcl", *bestKinPh_, isPassVL, 1.);

      double f_KtoVLexcl = getPhotonFR_KtoVLexcl(*bestKinPh_);
      photonFRClosure("KtoVLexcl_pt-aeta" /*_data*/ , *bestKinPh_, isPassVL, f_KtoVLexcl  );
    }
    // photonFakeRate_LtoT("KtoVL", *bestKinPh_, isPassVL, 1.);
  }

  if(bestMVAPh_){
    bool pass80 = bestMVAPh_->passMVA(Photon::MVAwp::wp80);
    bool pass90 = bestMVAPh_->passMVA(Photon::MVAwp::wp90);
    double mvaEffSF = 1.;

    if(pass90){
      mvaEffSF = pass80 ? getPhotonEffSF_MVA(*bestMVAPh_, Photon::MVAwp::wp80) : getPhotonEffSF_MVA(*bestMVAPh_, Photon::MVAwp::wp90);
      photonFakeRate_LtoT("90to80", *bestMVAPh_, pass80, mvaEffSF);

      // double f_90to80_data = getPhotonFR_90to80_data(thePh);
      // photonFRClosure("90to80_pt-aeta_data", *bestMVA, pass80, f_90to80_data);
    }
  }

  if(bestKinPh_ && bestMVAPh_){
    if(bestMVAPh_ != bestKinPh_){
      theHistograms->fill("best_CutVsMVA_nKinPh_disagree", ";# #gamma_{KIN};Events", 5,-0.5,4.5, kinPhotons_["central"]->size(), theWeight);
      theHistograms->fill("best_CutVsMVA_nCuts_disagree" , ";# cuts passed;Events" , 6,-0.5,5.5, bestKinPh_->nCutsPass(Photon::IdWp::Loose), theWeight);
      theHistograms->fill("best_CutVsMVA_MVA_disagree"   , ";# MVA score;Events"   , 10,-1,1   , bestMVAPh_->MVAvalue(), theWeight);
    }
    else{
      theHistograms->fill("best_CutVsMVA_nKinPh_agree"   , ";# #gamma_{KIN};Events", 5,-0.5,4.5, kinPhotons_["central"]->size(), theWeight);
      theHistograms->fill("best_CutVsMVA_nCuts_agree"    , ";# cuts passed;Events" , 6,-0.5,5.5, bestKinPh_->nCutsPass(Photon::IdWp::Loose), theWeight);
      theHistograms->fill("best_CutVsMVA_MVA_agree"      , ";# MVA score;Events"   , 10,-1,1   , bestMVAPh_->MVAvalue(), theWeight);
    }
  }

  orphanPhotonStudy();
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
    theHistograms->fill("Z0_vs_Z1_mass"         , "m_{Z0} [GeV/c^{2}];m_{Z1} [GeV/c^{2}]", 30,60.,90., 30,60.,90., ZZ->first().mass(), ZZ->second().mass(), theWeight);

    theHistograms->fill("ZZ_mass_" +channelReco_, "m_{4l};GeV/c^{2}", mVV_bins   , ZZ->mass()                   , theWeight);
    theHistograms->fill("Z0_mass"  +channelReco_, "m_{Z0};GeV/c^{2}", 35,55.,125., ZZ->first().mass()           , theWeight);
    theHistograms->fill("Z1_mass_" +channelReco_, "m_{Z1};GeV/c^{2}", 35,55.,125., ZZ->second().mass()          , theWeight);
    theHistograms->fill("ZZ_pt_"   +channelReco_, "p_{t,ZZ};GeV/c"  , 20,0.,400. , ZZ->pt()                     , theWeight);
    theHistograms->fill("Z0_l0_pt_"+channelReco_, "p_{t,l00};GeV/c" , 20,0.,400. , ZZ->first().daughter(0).pt() , theWeight);
    theHistograms->fill("Z0_l1_pt_"+channelReco_, "p_{t,l01};GeV/c" , 20,0.,400. , ZZ->first().daughter(1).pt() , theWeight);
    theHistograms->fill("Z1_l0_pt_"+channelReco_, "p_{t,l10};GeV/c" , 20,0.,400. , ZZ->second().daughter(0).pt(), theWeight);
    theHistograms->fill("Z1_l1_pt_"+channelReco_, "p_{t,l11};GeV/c" , 20,0.,400. , ZZ->second().daughter(1).pt(), theWeight);

    const char* ph_cutID = "noph";
    if     (goodPhotons_["central"]->size() > 0) ph_cutID = "loose";
    else if( kinPhotons_["central"]->size() > 0) ph_cutID = "kinVetoL";

    const char* jets_str = (jets_noph_->size() > 0 ? "gt1j" : "0j");

    theHistograms->fill(Form("Z0_mass_%s", ph_cutID), "m_{Z0};GeV/c^{2}", 35,55.,125., ZZ->first().mass()           , theWeight);
    theHistograms->fill(Form("Z1_mass_%s", ph_cutID), "m_{Z1};GeV/c^{2}", 35,55.,125., ZZ->second().mass()          , theWeight);

    theHistograms->fill(Form("Z0_mass_%s", jets_str), "m_{Z0};GeV/c^{2}", 35,55.,125., ZZ->first().mass()           , theWeight);
    theHistograms->fill(Form("Z1_mass_%s", jets_str), "m_{Z1};GeV/c^{2}", 35,55.,125., ZZ->second().mass()          , theWeight);

    if(theSampleInfo.isMC()){
	bool signaldef = sigdefHelper.pass_photon();  // Test if the event passes the GEN signal definition
	const char* sigdef_str = signaldef ? "prompt" : "nonpro";

	theHistograms->fill(Form("Z0_mass_%s"             , sigdef_str), "m_{Z0};GeV/c^{2}", 35,55.,125., ZZ->first().mass()           , theWeight);
	theHistograms->fill(Form("Z1_mass_%s"             , sigdef_str), "m_{Z1};GeV/c^{2}", 35,55.,125., ZZ->second().mass()          , theWeight);

	theHistograms->fill(Form("Z0_mass_%s_%s", ph_cutID, sigdef_str), "m_{Z0};GeV/c^{2}", 35,55.,125., ZZ->first().mass()           , theWeight);
	theHistograms->fill(Form("Z1_mass_%s_%s", ph_cutID, sigdef_str), "m_{Z1};GeV/c^{2}", 35,55.,125., ZZ->second().mass()          , theWeight);

	theHistograms->fill(Form("Z0_mass_%s_%s", jets_str, sigdef_str), "m_{Z0};GeV/c^{2}", 35,55.,125., ZZ->first().mass()           , theWeight);
	theHistograms->fill(Form("Z1_mass_%s_%s", jets_str, sigdef_str), "m_{Z1};GeV/c^{2}", 35,55.,125., ZZ->second().mass()          , theWeight);
    }

    ZllVsZllGstudy(*kinPhotons_["central"] , "kin"  );
    ZllVsZllGstudy(*goodPhotons_["central"], "loose");

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
    theHistograms->fill("Z_l0_pt"               , "p_{t,lZ0};GeV/c" , 20,0.,400. , ZL->first.daughter(0).pt()          , theWeight);
    theHistograms->fill("Z_l1_pt"               , "p_{t,lZ1};GeV/c" , 20,0.,400. , ZL->first.daughter(1).pt()          , theWeight);
    theHistograms->fill("L_pt"                  , "p_{t,l3};GeV/c"  , 20,0.,400. , ZL->second.pt()                     , theWeight);

    theHistograms->fill("ZL_mass_" +channelReco_, "m_{3l};GeV/c^{2}", 25,0.,500. , (ZL->first.p4()+ZL->second.p4()).M(), theWeight);
    theHistograms->fill("Z_mass_"  +channelReco_, "m_{Z};GeV/c^{2}" , 35,55.,125., ZL->first.mass()                    , theWeight);
    theHistograms->fill("Z_l0_pt_" +channelReco_, "p_{t,lZ0};GeV/c" , 20,0.,400. , ZL->first.daughter(0).pt()          , theWeight);
    theHistograms->fill("Z_l1_pt_" +channelReco_, "p_{t,lZ1};GeV/c" , 20,0.,400. , ZL->first.daughter(1).pt()          , theWeight);
    theHistograms->fill("L_pt_"    +channelReco_, "p_{t,l3};GeV/c"  , 20,0.,400. , ZL->second.pt()                     , theWeight);

    if(theSampleInfo.isMC()){
      bool isPrompt = sigdefHelper.pass_photon();
      std::string strPrompt = isPrompt ? "_prompt" : "_nonpro";
      theHistograms->fill("ZL_mass"               +strPrompt, ";m_{3l} [GeV/c^{2}]" , 25,0.,500. , (ZL->first.p4()+ZL->second.p4()).M(), theWeight);
      theHistograms->fill("Z_mass"                +strPrompt, ";m_{Z} [GeV/c^{2}]"  , 35,55.,125., ZL->first.mass()                    , theWeight);
      theHistograms->fill("Z_l0_pt"               +strPrompt, ";p_{T}^{lZ0} [GeV/c]", 20,0.,400. , ZL->first.daughter(0).pt()          , theWeight);
      theHistograms->fill("Z_l1_pt"               +strPrompt, ";p_{T}^{lZ1} [GeV/c]", 20,0.,400. , ZL->first.daughter(1).pt()          , theWeight);
      theHistograms->fill("L_pt"                  +strPrompt, ";p_{T}^{l3} [GeV/c]" , 20,0.,400. , ZL->second.pt()                     , theWeight);

      theHistograms->fill("ZL_mass_" +channelReco_+strPrompt, ";m_{3l} [GeV/c^{2}]" , 25,0.,500. , (ZL->first.p4()+ZL->second.p4()).M(), theWeight);
      theHistograms->fill("Z_mass_"  +channelReco_+strPrompt, ";m_{Z} [GeV/c^{2}]"  , 35,55.,125., ZL->first.mass()                    , theWeight);
      theHistograms->fill("Z_l0_pt_" +channelReco_+strPrompt, ";p_{T}^{lZ0} [GeV/c]", 20,0.,400. , ZL->first.daughter(0).pt()          , theWeight);
      theHistograms->fill("Z_l1_pt_" +channelReco_+strPrompt, ";p_{T}^{lZ1} [GeV/c]", 20,0.,400. , ZL->first.daughter(1).pt()          , theWeight);
      theHistograms->fill("L_pt_"    +channelReco_+strPrompt, ";p_{T}^{l3} [GeV/c]" , 20,0.,400. , ZL->second.pt()                     , theWeight);
    }
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
  
  cout<<"\n\t----- "<<phys::regionType(region_)<<" -----\n";
  cout<<Form("Total events: %lu (weighted: %.3g)\n", evtN, evtW);
  cout<<Form("Passing cut:  %lu (weighted: %.3g)\n", analyzedN, analyzedW);
  cout<<Form("Fraction:     %.1f %% (weighted: %.1f %%)\n", 100.*analyzedN/evtN, 100.*analyzedW/evtW);
  
  // Label names
  endNameHistos();

  for(const char* name_main : {"AAA_cuts", "AAA_cuts_u", "AAA_cuts_sigdef"}){
    const char* name_extra = Form("%s_extra", name_main);
    TH1* h_main  = theHistograms->get(name_main );
    TH1* h_extra = theHistograms->get(name_extra);
    if(h_main && h_extra){
      for(int bx = 1; bx <= h_extra->GetNbinsX(); ++bx){
	auto val = h_extra->GetBinContent(bx);
	auto err = h_extra->GetBinError  (bx);
	const char* new_label = h_extra->GetXaxis()->GetBinLabel(bx);
	int bnew = h_main->GetXaxis()->FindBin(new_label);
	h_main->SetBinContent(bnew, val);
	h_main->SetBinError  (bnew, err);
      }
      if(theHistograms->erase(name_extra))
	cout << "Error: problem deleting "<< name_extra << '\n';
    }
  }
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


bool VVGammaAnalyzer::SignalDefinitionHelper::pass() const{
  if     (is4Lregion(analyzer_->region_)){
    return photon_ && ZZ4L_;
  }
  else if(is3Lregion(analyzer_->region_)){
    // 3L signal definition
    return photon_ && WZ3L_;
  }
  else if(is2Lregion(analyzer_->region_)){
    // 2L 2j signal definition
    return photon_ && ZV2L2j_;
  }
  else if(analyzer_->region_ == CRLFR){
    // TODO
    return photon_;
  }
  else{
    std::cerr << "SignalDefinitionHelper::pass(): Unknown region: " << regionType(analyzer_->region_);
    return false;
  }
}

void VVGammaAnalyzer::SignalDefinitionHelper::eval(/*const VVGammaAnalyzer* analyzer*/){
  VVGammaAnalyzer* analyzer = analyzer_;
  Histogrammer* theHistograms = analyzer_->theHistograms;
  double theWeight = analyzer_->theWeight;
  phys::RegionTypes region_ = analyzer->region_;

  // General cuts
  bool passllLowMass_ = true;
  if(analyzer->genChLeptons_->size() >= 2){
    for(auto l1 = analyzer->genChLeptons_->begin(); l1 != analyzer->genChLeptons_->end() && passllLowMass_; ++l1)
      for(auto l2 = l1 + 1; l2 != analyzer->genChLeptons_->end(); ++l2)
	if( (l1->p4() + l2->p4()).M() ){
	  passllLowMass_ = false;
	  break;
	}
  }

  // Photon
  photon_ = eval_photon();

  bool regional = false;
  if     (is4Lregion(region_)) // 4L signal definition
    regional = eval_ZZ4L();
  else if(is3Lregion(region_)) // 3L signal definition
    regional = eval_WZ3L();
  else if(is2Lregion(region_)) // 2L 2j signal definition
    regional = eval_ZV2L2j();
  else if(region_ == CRLFR)    // CRLFR
    regional = true; // Nothing so far

  if(photon_) theHistograms->fill("DEBUG_sigdef", "Signal definition", {}, "photon", theWeight);

  if(photon_ && regional)
    theHistograms->fill("DEBUG_sigdef", "Signal definition", {}, "passVVG", theWeight);

}

bool VVGammaAnalyzer::SignalDefinitionHelper::eval_photon(){
  return std::any_of(analyzer_->genPhotonsPrompt_->begin(), analyzer_->genPhotonsPrompt_->end(),
		     [](const Particle& ph){
		       float aeta = fabs(ph.eta());
		       return ph.pt() > 20 && aeta < CUT_G_AETA_MAX;
		     });
}

bool VVGammaAnalyzer::SignalDefinitionHelper::eval_ZZ4L(){
  double theWeight = analyzer_->theWeight;

  // TODO analyzer->genChLeptons_->size() == 4
  bool topology = analyzer_->topology.test(TopologyBit::ZZto4L);
  bool l0_pt    = analyzer_->genChLeptons_->size() > 0 && analyzer_->genChLeptons_->at(0).pt() > 20;
  bool l1_pt    = analyzer_->genChLeptons_->size() > 1 && analyzer_->genChLeptons_->at(1).pt() > 10;
  bool fiducial = analyzer_->topology.test(TopologyBit::Fiducial);

  bool Z0def = analyzer_->GenZtoLLDefinition(analyzer_->genZZ_.first ());
  bool Z1def = analyzer_->GenZtoLLDefinition(analyzer_->genZZ_.second());

  ZZ4L_ = topology && l0_pt && l1_pt && fiducial && Z0def && Z1def;

  analyzer_->theHistograms->fill("DEBUG_sigdef", "", {"All","photon","topology","l0_pt","l1_pt","fiducial","Z0def","Z1def","passVV","passVVG"}, "All", theWeight);
  if(topology) analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "topology", theWeight);
  if(l0_pt)    analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "l0_pt"   , theWeight);
  if(l1_pt)    analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "l1_pt"   , theWeight);
  if(fiducial) analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "fiducial", theWeight);
  if(Z0def)    analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "Z0def"   , theWeight);
  if(Z1def)    analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "Z1def"   , theWeight);
  if(ZZ4L_)    analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "passVV"  , theWeight);

  analyzer_->fillCutsNm1("DEBUG_sigdef_ZZ_Nm1"    , "N-1 (ZZ);cut excluded;Events"       , {
      {"topology"      , topology},
      {"l0_pt && l1_pt", l0_pt && l1_pt},
      {"fiducial"      , fiducial},
      {"Z0def && Z1def", Z0def && Z1def}
    }, theWeight);
  analyzer_->fillCutsNm1("DEBUG_sigdef_Nm1"       , "N-1 (ZZ+#gamma);cut excluded;Events", {
      {"ZZ", ZZ4L_}, {"photon"  , photon_ }
    }, theWeight);
  analyzer_->fillCutFlow("DEBUG_sigdef_cutflow"   , "GEN cutflow (ZZ+#gamma);;Events", {
      {"topology", topology},
      {"l0_pt"   , l0_pt   },
      {"l1_pt"   , l1_pt   },
      {"fiducial", fiducial},
      {"Z0def"   , Z0def   },
      {"Z1def"   , Z1def   },
      {"photon"  , photon_ }
    }, theWeight);


  return ZZ4L_;
}

bool VVGammaAnalyzer::SignalDefinitionHelper::eval_WZ3L(){
  double theWeight = analyzer_->theWeight;

  // TODO analyzer->genChLeptons_->size() == 3
  bool topology = analyzer_->topology.test(TopologyBit::WZto3LNu);
  bool lZ0_pt   = analyzer_->genZW_.first() .daughter(0).pt() > 20;
  bool lZ1_pt   = analyzer_->genZW_.first() .daughter(1).pt() > 10;
  bool lW_pt    = analyzer_->genZW_.second().daughter(0).pt() > 20;
  bool fiducial = analyzer_->topology.test(TopologyBit::Fiducial);

  bool Zdef = analyzer_->GenZtoLLDefinition (analyzer_->genZW_.first ());
  bool Wdef = analyzer_->GenWtoLNuDefinition(analyzer_->genZW_.second());

  WZ3L_ = topology && lZ0_pt && lZ1_pt && lW_pt && Zdef && fiducial;

  analyzer_->theHistograms->fill("DEBUG_sigdef", "", {"All","photon","topology","lZ0_pt","lZ1_pt","lW_pt","fiducial","Zdef","Wdef","passVV","passVVG"}, "All", theWeight);
  if(topology) analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "topology", theWeight);
  if(lZ0_pt)   analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "lZ0_pt"  , theWeight);
  if(lZ1_pt)   analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "lZ1_pt"  , theWeight);
  if(lW_pt)    analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "lW_pt"   , theWeight);
  if(fiducial) analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "fiducial", theWeight);
  if(Zdef)     analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "Zdef"    , theWeight);
  if(Wdef)     analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "Wdef"    , theWeight);
  if(WZ3L_)    analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "passVV"  , theWeight);

  return WZ3L_;
}

bool VVGammaAnalyzer::SignalDefinitionHelper::eval_ZV2L2j(){
  double theWeight = analyzer_->theWeight;
  // TODO analyzer_->genChLeptons_->size() == 2

  // Z --> ll
  Zll_ = analyzer_->topology.test(TopologyBit::ZtoLL);//(analyzer_->genZlepCandidates_->size() == 1) && analyzer->GenZtoLLDefinition(analyzer_->genZlepCandidates_->at(0));

  // V --> qq
  bool Zhad = analyzer_->genZhadCandidates_->size() >= 1;
  bool Whad = analyzer_->genWhadCandidates_->size() >= 1;
  // TODO: fix quark pt and aeta requirements
  bool qq_accept = false;
  if(Zhad)
    qq_accept = std::any_of(analyzer_->genZhadCandidates_->begin(), analyzer_->genZhadCandidates_->end(),
				 [this](const Boson<Particle>& VB){ return analyzer_->GenVtoQQDefinition(VB); }
				 );

  if(Whad && !qq_accept)
    qq_accept = std::any_of(analyzer_->genWhadCandidates_->begin(), analyzer_->genWhadCandidates_->end(),
				 [this](const Boson<Particle>& VB){ return analyzer_->GenVtoQQDefinition(VB); }
				 );

  Vhad_ = (Zhad || Whad) && qq_accept;
  ZV2L2j_ = Zll_ && Vhad_;

  analyzer_->theHistograms->fill("DEBUG_sigdef", "", {"All","photon","Zll_topo","Vhad","qq_accept","passVV","passVVG"}, "All", theWeight);
  if(Zll_)     analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "Zll_topo" , theWeight);
  if(Vhad_)    analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "Vhad"     , theWeight);
  if(qq_accept)analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "qq_accept",theWeight);
  if(ZV2L2j_)  analyzer_->theHistograms->fill("DEBUG_sigdef", "", {}, "passVV"   , theWeight);

  return ZV2L2j_;
}


std::unique_ptr<TH2F> VVGammaAnalyzer::getHistfromFile(const char* fname, const char* hname, const char* info){
  TFile fileFR(Form(fname, "READ" ));
  if(!fileFR.IsOpen()){
    cout << colour::Red("WARN") << ": file"<<info<<" not found in \""<<fname<<"\"\n";
    return std::unique_ptr<TH2F> (new TH2F(hname, "DEFAULT", 1,0.,1., 1,0.,1.) );
  }

  TObject* retrievedObj = fileFR.Get(hname);
  if(!retrievedObj){
    cout << colour::Red("WARN") << ": histogram \""<<hname<<"\""<<info<<" not found in \""<<fname<<"\"\n"
	 << "\tavailable keys:" ;
    for(auto key : *fileFR.GetListOfKeys()) cout << ' ' << key->GetName() << ',';
    cout << std::endl;
    fileFR.Close();
    return std::unique_ptr<TH2F> (new TH2F(hname, "DEFAULT", 1,0.,1., 1,0.,1.) );
  }

  std::unique_ptr<TH2F> result((TH2F*) retrievedObj);
  result->SetDirectory(nullptr);  // prevent ROOT from deleting it
  cout << "INFO: retrieved histogram \""<<hname<<"\""<<info<<" from \""<<fname<<"\"\n";
  fileFR.Close();
  return result;
}


bool VVGammaAnalyzer::canBeFSR(const Photon& ph, const vector<Lepton>& leptons) const{
  // See: https://github.com/CJLST/ZZAnalysis/blob/Run2UL_22/AnalysisStep/plugins/LeptonPhotonMatcher.cc#L187-L250
  float pt = ph.pt();

  for(const Lepton& lep : leptons){
    // See:
    //   ZZAnalysis/AnalysisStep/test/MasterPy/ZZ4lAnalysis.py:321        :SIP =  "abs(dB('PV3D')/edB('PV3D')) < 4"
    //   ZZAnalysis/AnalysisStep/test/MasterPy/ZZ4lAnalysis.py:396,481,509:        isSIP = cms.string(SIP),
    if(! (lep.sip() < 4))
      continue;

    float dR = deltaR(ph, lep);
    if(dR < 0.5  && dR/(pt*pt) < 0.012){
      // TODO? Check the relative isolation
      // Trying to replicate https://github.com/CJLST/ZZAnalysis/blob/Run2UL_22/AnalysisStep/src/LeptonIsoHelper.cc#L210-245
      return true;
    }
  }

  return false;
}


void VVGammaAnalyzer::initCherryPick(){
  FILE* cherryFile = fopen(Form("data/cherry-pick.txt"), "r");
  if(!cherryFile){
    cout << colour::Red("WARN") << ": no cherry pick file\n";
    return;
  }
  unsigned long r, l, e;
  while( fscanf(cherryFile, "%lu:%lu:%lu", &r, &l, &e) == 3 )
    cherryEvents[r][l].insert(e);
  fclose(cherryFile);
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
  genPhotonsPrompt_->clear();
	
  genZlepCandidates_->clear();
  genWlepCandidates_->clear();
  genZhadCandidates_->clear();
  genWhadCandidates_->clear();
	
  genZZ_ = DiBoson<Particle, Particle>();
  genZW_ = DiBoson<Particle, Particle>();
	
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
      if(p.genStatusFlags().test(phys::isPrompt))
	genPhotonsPrompt_->push_back(p);
    }
  }
	
  // Gen W --> l nu
  if(genNeutrinos_->size() > 0 && genChLeptons_->size() > 0){
    for(auto l : *genChLeptons_){
      for(auto v : *genNeutrinos_){
	if( abs(l.id() + v.id()) == 1 ){
	  Boson<Particle> Wcand(l,v);
	  if(GenWtoLNuDefinition(Wcand))
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
	  if(GenVtoQQDefinition(Wcand))
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
    auto it_Z1 = std::find_if(genZlepCandidates_->begin(), genZlepCandidates_->end(),
			      [&Z0](auto Z){ return !haveCommonDaughter(Z0, Z); }
			      );
    if(it_Z1 != genZlepCandidates_->end()){
      genZZ_ = DiBoson<Particle, Particle>(Z0, *it_Z1);
    }
  }

  // genZW --> 3l nu
  if(genChLeptons_->size() >= 3 && genZlepCandidates_->size() >= 1 && genWlepCandidates_->size() >= 1){	
    std::sort(genZlepCandidates_->begin(), genZlepCandidates_->end(), MassComparator(phys::ZMASS));
    Boson<Particle>& Z0 = genZlepCandidates_->front();

    std::sort(genWlepCandidates_->begin(), genWlepCandidates_->end(), MassComparator(phys::WMASS));
    auto it_W = std::find_if(genWlepCandidates_->begin(), genWlepCandidates_->end(),
			     [&Z0](auto W){ return !haveCommonDaughter(Z0, W); }
			     );
    if(it_W != genWlepCandidates_->end()){
      genZW_ = DiBoson<Particle, Particle>(Z0, *it_W);
    }
  }

  std::sort(genQuarks_       ->begin(), genQuarks_       ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genChLeptons_    ->begin(), genChLeptons_    ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genNeutrinos_    ->begin(), genNeutrinos_    ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genPhotons_      ->begin(), genPhotons_      ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genPhotonsPrompt_->begin(), genPhotonsPrompt_->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
}


void VVGammaAnalyzer::genEventHistos(){
  theHistograms->fill("GEN_ZlepCandidates", "# genZlepCandidates_", 4,-0.5,3.5, genZlepCandidates_->size(), theWeight);
  theHistograms->fill("GEN_WlepCandidates", "# genWlepCandidates_", 4,-0.5,3.5, genWlepCandidates_->size(), theWeight);
  theHistograms->fill("GEN_ZhadCandidates", "# genZhadCandidates_", 4,-0.5,3.5, genZhadCandidates_->size(), theWeight);
  theHistograms->fill("GEN_WhadCandidates", "# genWhadCandidates_", 4,-0.5,3.5, genWhadCandidates_->size(), theWeight);
  
  theHistograms->fill("GEN_quarks"   , "# genQuarks"   , 10,-0.5,9.5, genQuarks_   ->size(), theWeight);
  theHistograms->fill("GEN_chLeptons", "# genChLeptons", 10,-0.5,9.5, genChLeptons_->size(), theWeight);
  theHistograms->fill("GEN_neutrinos", "# genNeutrinos", 10,-0.5,9.5, genNeutrinos_->size(), theWeight);
  theHistograms->fill("GEN_photons"  , "# genPhotons"  , 10,-0.5,9.5, genPhotons_  ->size(), theWeight);
  theHistograms->fill("GEN_photonsPrompt", "# genPhotonsPrompt", 10,-0.5,9.5, genPhotonsPrompt_->size(), theWeight);
  
  for(auto v : *genZlepCandidates_)
    theHistograms->fill("GEN_genZlepCandidates_mass", "mass genZlepCandidates;[GeV/c^{2}]", 35.,50.,120., v.mass(), theWeight);
  for(auto v : *genWlepCandidates_)
    theHistograms->fill("GEN_genWlepCandidates_mass", "mass genWlepCandidates;[GeV/c^{2}]", 35.,50.,120., v.mass(), theWeight);
  for(auto v : *genZhadCandidates_){
    theHistograms->fill("GEN_genZhadCandidates_mass", "mass genZhadCandidates;[GeV/c^{2}]", 35.,50.,120., v.mass(), theWeight);
    theHistograms->fill("GEN_genZhadCandidates_pt"  , "pt   genZhadCandidates;[GeV/c]    ", 60., 0.,300., v.pt()  , theWeight);
    theHistograms->fill("GEN_genZhadCandidates_dRqq", "dRqq genZhadCandidates;#DeltaR(q,q)",60., 0.,  6., deltaR(v.daughter(0), v.daughter(1)), theWeight);
  }
  for(auto v : *genWhadCandidates_){
    theHistograms->fill("GEN_genWhadCandidates_mass", "mass genWhadCandidates;[GeV/c^{2}]", 35.,50.,120., v.mass(), theWeight);
    theHistograms->fill("GEN_genWhadCandidates_pt"  , "pt   genWhadCandidates;[GeV/c]    ", 60.,0. ,300., v.pt()  , theWeight);
    theHistograms->fill("GEN_genWhadCandidates_dRqq", "dRqq genWhadCandidates;#DeltaR(q,q)",60., 0.,  6., deltaR(v.daughter(0), v.daughter(1)), theWeight);
  }

  theHistograms->fill("GEN_n_jets"   , "Number of genJets;# genJets"   , 6,-0.5,5.5, genJets->size()   , theWeight);
  theHistograms->fill("GEN_n_jetsAK8", "Number of genJetsAK8;# genJets", 6,-0.5,5.5, genJetsAK8->size(), theWeight);
  
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
    theHistograms->fill("lead_muo_pt", "Leading muon p_{t}:p_{t} [GeV/c]"   , 50, 0., 250., lead_muo->pt(), theWeight);
  
  // MET
  theHistograms->fill("MET", "MET;#slash{E}_{T} [GeV/c]", 50, 0., 250., met->pt(), theWeight);
}


void VVGammaAnalyzer::fillPhotonPlots(const Photon& ph, const char* name, const char* title){
    float pt    = ph.pt();
    float aeta  = fabs(ph.eta());
    float MVAv  = ph.MVAvalue();
    auto closestLep = closestDeltaR(ph, *leptons_);
    float dRl = closestLep != leptons_->cend() ? deltaR(ph, *closestLep) : 10;
    float chIso = ph.chargedIsolation();
    float sieie = ph.sigmaIetaIeta();

    theHistograms->fill(Form("%s_pt_fine", name), Form("%s photon;p_{T} [GeV/c];Events"     , title), ph_ptExtended_bins  , pt   , theWeight);
    theHistograms->fill(Form("%s_pt"     , name), Form("%s photon;p_{T} [GeV/c];Events"     , title), ph_pt_bins          , pt   , theWeight);
    theHistograms->fill(Form("%s_aeta"   , name), Form("%s photon;#eta;Events"              , title), ph_aeta_bins        , aeta , theWeight);
    theHistograms->fill(Form("%s_aeta_fine",name),Form("%s photon;#eta;Events"              , title), ph_aetaExtended_bins, aeta , theWeight);
    theHistograms->fill(Form("%s_dRl"    , name), Form("%s photon;#DeltaR(#gamma, l);Events", title), 40,0.,1.            , dRl  , theWeight);
    theHistograms->fill(Form("%s_dRl_fine",name),Form("%s photon;#DeltaR(#gamma, l);Events", title), 100,0.,1.           , dRl  , theWeight);
    theHistograms->fill(Form("%s_MVA"    , name), Form("%s photon;MVA"                      , title), 40,-1.,1.           , MVAv , theWeight);
    theHistograms->fill(Form("%s_chIso"  , name), Form("%s photon;chIso (uncorrected)"      , title), 40, 0., 10          , chIso, theWeight);
    theHistograms->fill(Form("%s_sieie"  , name), Form("%s photon;#sigma_{i#etai#eta}"      , title), 40, 0., .08         , sieie, theWeight);

    if(theSampleInfo.isMC()){
      const char* genStatus = sigdefHelper.pass_photon() ? "prompt" : "nonpro" ;
      theHistograms->fill(Form("%s_pt_fine_%s"  , name, genStatus), Form("%s photon;p_{T} [GeV/c];Events"      , title), ph_ptExtended_bins  , pt   , theWeight);
      theHistograms->fill(Form("%s_pt_%s"       , name, genStatus), Form("%s photon;p_{T} [GeV/c];Events"      , title), ph_pt_bins          , pt   , theWeight);
      theHistograms->fill(Form("%s_aeta_%s"     , name, genStatus), Form("%s photon;#eta;Events"               , title), ph_aeta_bins        , aeta , theWeight);
      theHistograms->fill(Form("%s_aeta_fine_%s", name, genStatus), Form("%s photon;#eta;Events"               , title), ph_aetaExtended_bins, aeta , theWeight);
      theHistograms->fill(Form("%s_dRl_%s"      , name, genStatus), Form("%s photon;#DeltaR(#gamma, l);Events" , title),  40, 0., 1.         , dRl  , theWeight);
      theHistograms->fill(Form("%s_dRl_fine_%s" , name, genStatus), Form("%s photon;#DeltaR(#gamma, l);Events" , title), 100, 0., 1.         , dRl  , theWeight);
      theHistograms->fill(Form("%s_MVA_%s"      , name, genStatus), Form("%s photon;MVA;Events"                , title),  40,-1., 1.         , MVAv , theWeight);
      theHistograms->fill(Form("%s_chIso_%s"    , name, genStatus), Form("%s photon;chIso (uncorrected);Events", title),  40, 0., 10         , chIso, theWeight);
      theHistograms->fill(Form("%s_sieie_%s"    , name, genStatus), Form("%s photon;#sigma_{i#etai#eta};Events", title),  40, 0., .08        , sieie, theWeight);
    }
    else if(strcmp(name, "lead_fail") == 0){  // is data and photon is fail (VL && !L) --> reweight
      double f_VLtoL_data   = getPhotonFR_VLtoL_data(  ph);
      double f_VLtoL_dataZG = getPhotonFR_VLtoL_dataZG(ph);
      double weight_VLtoL_data   = theWeight * f_VLtoL_data  /(1-f_VLtoL_data  );
      double weight_VLtoL_dataZG = theWeight * f_VLtoL_dataZG/(1-f_VLtoL_dataZG);

      theHistograms->fill(Form("%s_pt_fine_%s"  , name, "reweight_data"), Form("%s photon;p_{T} [GeV/c];Events"      , title), ph_ptExtended_bins  , pt   , weight_VLtoL_data);
      theHistograms->fill(Form("%s_pt_%s"       , name, "reweight_data"), Form("%s photon;p_{T} [GeV/c];Events"      , title), ph_pt_bins          , pt   , weight_VLtoL_data);
      theHistograms->fill(Form("%s_aeta_%s"     , name, "reweight_data"), Form("%s photon;#eta;Events"               , title), ph_aeta_bins        , aeta , weight_VLtoL_data);
      theHistograms->fill(Form("%s_aeta_fine_%s", name, "reweight_data"), Form("%s photon;#eta;Events"               , title), ph_aetaExtended_bins, aeta , weight_VLtoL_data);
      theHistograms->fill(Form("%s_dRl_%s"      , name, "reweight_data"), Form("%s photon;#DeltaR(#gamma, l);Events" , title),  40, 0., 1.         , dRl  , weight_VLtoL_data);
      theHistograms->fill(Form("%s_dRl_fine_%s" , name, "reweight_data"), Form("%s photon;#DeltaR(#gamma, l);Events" , title), 100, 0., 1.         , dRl  , weight_VLtoL_data);
      theHistograms->fill(Form("%s_MVA_%s"      , name, "reweight_data"), Form("%s photon;MVA;Events"                , title),  40,-1., 1.         , MVAv , weight_VLtoL_data);
      theHistograms->fill(Form("%s_chIso_%s"    , name, "reweight_data"), Form("%s photon;chIso (uncorrected);Events", title),  40, 0., 10         , chIso, weight_VLtoL_data);
      theHistograms->fill(Form("%s_sieie_%s"    , name, "reweight_data"), Form("%s photon;#sigma_{i#etai#eta};Events", title),  40, 0., .08        , sieie, weight_VLtoL_data);

      theHistograms->fill(Form("%s_pt_fine_%s"  , name, "reweight_dataZG"), Form("%s photon;p_{T} [GeV/c];Events"      , title), ph_ptExtended_bins  , pt   , weight_VLtoL_dataZG);
      theHistograms->fill(Form("%s_pt_%s"       , name, "reweight_dataZG"), Form("%s photon;p_{T} [GeV/c];Events"      , title), ph_pt_bins          , pt   , weight_VLtoL_dataZG);
      theHistograms->fill(Form("%s_aeta_%s"     , name, "reweight_dataZG"), Form("%s photon;#eta;Events"               , title), ph_aeta_bins        , aeta , weight_VLtoL_dataZG);
      theHistograms->fill(Form("%s_aeta_fine_%s", name, "reweight_dataZG"), Form("%s photon;#eta;Events"               , title), ph_aetaExtended_bins, aeta , weight_VLtoL_dataZG);
      theHistograms->fill(Form("%s_dRl_%s"      , name, "reweight_dataZG"), Form("%s photon;#DeltaR(#gamma, l);Events" , title),  40, 0., 1.         , dRl  , weight_VLtoL_dataZG);
      theHistograms->fill(Form("%s_dRl_fine_%s" , name, "reweight_dataZG"), Form("%s photon;#DeltaR(#gamma, l);Events" , title), 100, 0., 1.         , dRl  , weight_VLtoL_dataZG);
      theHistograms->fill(Form("%s_MVA_%s"      , name, "reweight_dataZG"), Form("%s photon;MVA;Events"                , title),  40,-1., 1.         , MVAv , weight_VLtoL_dataZG);
      theHistograms->fill(Form("%s_chIso_%s"    , name, "reweight_dataZG"), Form("%s photon;chIso (uncorrected);Events", title),  40, 0., 10         , chIso, weight_VLtoL_dataZG);
      theHistograms->fill(Form("%s_sieie_%s"    , name, "reweight_dataZG"), Form("%s photon;#sigma_{i#etai#eta};Events", title),  40, 0., .08        , sieie, weight_VLtoL_dataZG);
    }

    static vector<double> dRl_bins {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
    float ptl        = closestLep != leptons_->cend() ? closestLep->pt()                  : -1.;
    float RelIso     = closestLep != leptons_->cend() ? closestLep->pfCombRelIso()        : -1.;
    float RelIsoCorr = closestLep != leptons_->cend() ? closestLep->pfCombRelIsoFSRCorr() : -1.;
    theHistograms->fill(Form("%s_pt_aeta", name), Form("%s photon;p_{T};|#eta|;Events"                   , title), ph_pt_bins, ph_aeta_bins, pt , aeta      , theWeight);
    theHistograms->fill(Form("%s_dRl_pt" , name), Form("%s photon;#DeltaR(#gamma, l);p_{T};Events"       , title), dRl_bins  , ph_pt_bins  , dRl, pt        , theWeight);
    theHistograms->fill(Form("%s_dRl_ptl", name), Form("%s photon;#DeltaR(#gamma, l);p_{T}^{l};Events"   , title), dRl_bins  , ph_pt_bins  , dRl, ptl       , theWeight);
    theHistograms->fill(Form("%s_dRl_pt" , name), Form("%s photon;#DeltaR(#gamma, l);l RelIso;Events"    , title), 10,0.,1.  , 16,0.,0.16  , dRl, RelIso    , theWeight);
    theHistograms->fill(Form("%s_dRl_pt" , name), Form("%s photon;#DeltaR(#gamma, l);l RelIsoCorr;Events", title), 10,0.,1.  , 16,0.,0.16  , dRl, RelIsoCorr, theWeight);
}


void VVGammaAnalyzer::photonHistos(){
  std::string strPrompt = "";
  if(theSampleInfo.isMC())
    strPrompt = sigdefHelper.pass_photon() ? "_prompt" : "_nonpro";

  theHistograms->fill(  "fsrPhotons_N"          , "# #gamma_{FSR};;Events", 5, -0.5, 4.5, fsrPhotons_->size(), theWeight);
  if(theSampleInfo.isMC())
    theHistograms->fill("fsrPhotons_N"+strPrompt, "# #gamma_{FSR};;Events", 5, -0.5, 4.5, fsrPhotons_->size(), theWeight);

  if(fsrPhotons_->size() >= 1){
    const Particle& ph = fsrPhotons_->at(0);
    auto closestLep = closestDeltaR(ph, *leptons_);
    float dRl = closestLep != leptons_->cend() ? deltaR(ph, *closestLep) : 10;
    float pt = ph.pt();

    theHistograms->fill(  "lead_fsrPhotons_pt"             , ";p_{T} [GeV];Events"       , 60, 0., 180, pt < 180 ? pt : 179, theWeight);
    theHistograms->fill(  "lead_fsrPhotons_aeta"           , ";|#eta|;Events"            , 60, 0., CUT_G_AETA_MAX, fabs(ph.eta()), theWeight);
    theHistograms->fill(  "lead_fsrPhotons_dRl"            , ";#DeltaR(#gamma, l);Events", 60, 0., 0.6, dRl < 1  ? dRl: .99, theWeight);

    if(theSampleInfo.isMC()){
      theHistograms->fill("lead_fsrPhotons_pt"  + strPrompt, ";p_{T} [GeV];Events"       , 60, 0., 180, pt < 180 ? pt : 179, theWeight);
      theHistograms->fill("lead_fsrPhotons_aeta"+ strPrompt, ";|#eta|;Events"            , 60, 0., CUT_G_AETA_MAX, fabs(ph.eta()), theWeight);
      theHistograms->fill("lead_fsrPhotons_dRl" + strPrompt, ";#DeltaR(#gamma, l);Events", 60, 0., 0.6, dRl < 1  ? dRl: .99, theWeight);
    }
  }
  if(fsrPhotons_->size() >= 2){
    const Particle& ph = fsrPhotons_->at(1);
    auto closestLep = closestDeltaR(ph, *leptons_);
    float dRl = closestLep != leptons_->cend() ? deltaR(ph, *closestLep) : 10;
    float pt = ph.pt();

    theHistograms->fill(  "sublead_fsrPhotons_pt"            , ";p_{T};Events"              , 60, 0., 180, pt < 180 ? pt : 179., theWeight);
    theHistograms->fill(  "sublead_fsrPhotons_aeta"          , ";#eta;Events"               , 60, 0., CUT_G_AETA_MAX, fabs(ph.eta()), theWeight);
    theHistograms->fill(  "sublead_fsrPhotons_dRl"           , ";#DeltaR(#gamma, l);Events" , 60, 0., 0.6, dRl < 1  ? dRl: .99 , theWeight);
    if(theSampleInfo.isMC()){
      theHistograms->fill("sublead_fsrPhotons_pt"  +strPrompt, ";p_{T};Events"              , 60, 0., 180, pt < 180 ? pt : 179., theWeight);
      theHistograms->fill("sublead_fsrPhotons_aeta"+strPrompt, ";#eta;Events"               , 60, 0., CUT_G_AETA_MAX, fabs(ph.eta()), theWeight);
      theHistograms->fill("sublead_fsrPhotons_dRl" +strPrompt, ";#DeltaR(#gamma, l);Events" , 60, 0., 0.6, dRl < 1  ? dRl: .99 , theWeight);
    }
  }

  // No photons passing kinematic cuts
  if(theSampleInfo.isMC() && kinPhotons_["central"]->size() == 0){
    size_t nGenPh = 0;
    for(const Particle& genPh : *genPhotons_){
      if(genPh.pt() < 15)
	continue;
      ++nGenPh;
      theHistograms->fill("noKinPh_all_genPh_pt" , "#gamma_{GEN} (p_{T}>15) when no #gamma_{REC};#gamma_{GEN} p_{T}", 60,15.,165., genPh.pt() , theWeight);
      theHistograms->fill("noKinPh_all_genPh_eta", "#gamma_{GEN} (p_{T}>15) when no #gamma_{REC};#gamma_{GEN} #eta" , 40,-4.,  4., genPh.eta(), theWeight);
      auto it_ph = std::min_element(photons->begin(), photons->end(), DeltaRComparator(genPh));
      if(it_ph != photons->end() && deltaR(*it_ph, genPh) < 0.2){
	theHistograms->fill("noKinPh_rec_genPh_pt" , "#gamma_{GEN} (p_{T}>15) when no #gamma_{REC};#gamma_{GEN} p_{T}", 60,15.,165., genPh.pt() , theWeight);
	theHistograms->fill("noKinPh_rec_genPh_eta", "#gamma_{GEN} (p_{T}>15) when no #gamma_{REC};#gamma_{GEN} #eta" , 40,-4.,  4., genPh.eta(), theWeight);
      }
    }
    theHistograms->fill("noKinPh_all_genPh_N" , "#gamma_{GEN} (p_{T}>15) when no #gamma_{REC};# #gamma_{GEN}", 5,-0.5,4.5, nGenPh, theWeight);
  }
  
  // pt and eta distribution of photons
  theHistograms->fill(  "kinPhotons_N"                , "Kin photons;N #gamma;Events"      , 5,-0.5,4.5,   kinPhotons_["central"]->size(), theWeight);
  theHistograms->fill(  "veryLoosePhotons_N"          , "VeryLoose photons;N #gamma;Events", 5,-0.5,4.5, loosePhotons_["central"]->size(), theWeight);
  theHistograms->fill(  "loosePhotons_N"              , "Loose photons;N #gamma;Events"    , 5,-0.5,4.5,  goodPhotons_["central"]->size(), theWeight);
  if(theSampleInfo.isMC()){
    theHistograms->fill("kinPhotons_N"      +strPrompt, "Kin photons;N #gamma;Events"      , 5,-0.5,4.5,   kinPhotons_["central"]->size(), theWeight);
    theHistograms->fill("veryLoosePhotons_N"+strPrompt, "VeryLoose photons;N #gamma;Events", 5,-0.5,4.5, loosePhotons_["central"]->size(), theWeight);
    theHistograms->fill("loosePhotons_N"    +strPrompt, "Loose photons;N #gamma;Events"    , 5,-0.5,4.5,  goodPhotons_["central"]->size(), theWeight);
  }

  if  ( goodPhotons_["central"]->size() > 0){
    fillPhotonPlots(goodPhotons_["central"]->at(0), "lead_loose", "Leading Loose");
    if( goodPhotons_["central"]->size() > 1)
	theHistograms->fill("sublead_loose_pt", "Subleading Loose #gamma;p_{T} [GeV/c]", ph_pt_bins, goodPhotons_["central"]->at(1).pt());
  }
  else{
    if(loosePhotons_["central"]->size() > 0){
      const Photon& ph = loosePhotons_["central"]->at(0);
      fillPhotonPlots(ph, "lead_fail", "Leading Fail");
      if(loosePhotons_["central"]->size() > 1)
	theHistograms->fill("sublead_fail_pt", "Subleading Fail #gamma;p_{T} [GeV/c]", ph_pt_bins, loosePhotons_["central"]->at(1).pt());
      bool passChIso = ph.cutBasedID(Photon::IdWp::Loose, Photon::IDcut::chIso);
      bool passSieie = ph.cutBasedID(Photon::IdWp::Loose, Photon::IDcut::sieie);

      if( passChIso && !passSieie)  // 4a = passChIso
	  fillPhotonPlots(ph, "lead_fail4a", "Leading 4a: pass chIso");
      if(!passChIso &&  passSieie)  // 4b = passSieie
	  fillPhotonPlots(ph, "lead_fail4b", "Leading 4a: pass sieie");
      if(!passChIso && !passSieie)
	fillPhotonPlots(ph, "lead_fail3", "Leading fail3");
    }

    if(  kinPhotons_["central"]->size() > 0){
      fillPhotonPlots(kinPhotons_["central"]->at(0), "lead_kinVetoL", "Leading Kin-L");
      if(kinPhotons_["central"]->size() > 1)
	theHistograms->fill("sublead_kinVetoL_pt", "Subleading Kin-L #gamma;p_{T} [GeV/c]", ph_pt_bins, kinPhotons_["central"]->at(1).pt());
    }
  }
  if(bestMVAPh_){
    if(bestMVAPh_->passMVA(Photon::MVAwp::wp90)){
      fillPhotonPlots(*bestMVAPh_, "lead_wp90", "Leading wp90");
      if(bestMVAPh_->passMVA(Photon::MVAwp::wp80))
	fillPhotonPlots(*bestMVAPh_, "lead_wp80", "Leading wp80");
      else
	fillPhotonPlots(*bestMVAPh_, "lead_90not80", "Leading 90not80");
    }
  }

  // SCEta effect on VeryLoose wp
  auto n_SCEta = std::count_if(loosePhotons_["central"]->begin(), loosePhotons_["central"]->end(),
			       [](const Photon& ph) { return ph.cutBasedID(Photon::IdWp::VeryLoose, Photon::IDcut::SCEta); }
			       );
  theHistograms->fill("veryLoosePh_SCEta_N", "Number of VeryLoose (+SCEta) photons", 5,-0.5,4.5, n_SCEta, theWeight);

  // Furthest from leptons
  {
    double dRl;
    std::vector<Photon>::const_iterator itPh;
    if(kinPhotons_["central"]->size() > 0){
      std::tie(itPh, dRl) = furthestFromAny(*kinPhotons_["central"], *leptons_);
      theHistograms->fill("furthestKinPh" , "Kin photon furthest from any lepton;#DeltaR(#gamma, l);Events"      , 40,0.,1., dRl, theWeight);
    }
    if(loosePhotons_["central"]->size() > 0){
      std::tie(itPh, dRl) = furthestFromAny(*loosePhotons_["central"], *leptons_);
      theHistograms->fill("furthestVLPh"  , "VeryLoose photon furthest from any lepton;#DeltaR(#gamma, l);Events", 40,0.,1., dRl, theWeight);
      if(!itPh->cutBasedIDLoose())
	theHistograms->fill("furthestFailPh", "Fail photon furthest from any lepton;#DeltaR(#gamma, l);Events"   , 40,0.,1., dRl, theWeight);
    }
    if(goodPhotons_["central"]->size() > 0){
      std::tie(itPh, dRl) = furthestFromAny(*goodPhotons_["central"], *leptons_);
      theHistograms->fill("furthestLoosePh", "Loose photon furthest from any lepton;#DeltaR(#gamma, l);Events"   , 40,0.,1., dRl, theWeight);
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

      theHistograms->fill(  "kinPh_sieie_EB"           , "kinPhotons in Barrel;#sigma_{i#etai#eta}", sieie_bins, ph.sigmaIetaIeta()        , theWeight);
      theHistograms->fill(  "kinPh_chIso_EB"           , "kinPhotons in Barrel;chIso"              , 60,0.,60., ph.chargedIsolation()      , theWeight);
      theHistograms->fill(  "kinPh_HoverE_EB"          , "kinPhotons in Barrel;HoverE"             , 60,0.,.15, ph.HoverE()                , theWeight);
      theHistograms->fill(  "kinPh_neIso_EB"           , "kinPhotons in Barrel;neIso"              , 60,0.,60., ph.neutralHadronIsolation(), theWeight);
      theHistograms->fill(  "kinPh_phIso_EB"           , "kinPhotons in Barrel;phIso"              , 60,0.,60., ph.photonIsolation()       , theWeight);
      if(theSampleInfo.isMC()){
	theHistograms->fill("kinPh_sieie_EB" +strPrompt, "kinPhotons in Barrel;#sigma_{i#etai#eta}", sieie_bins, ph.sigmaIetaIeta()        , theWeight);
	theHistograms->fill("kinPh_chIso_EB" +strPrompt, "kinPhotons in Barrel;chIso"              , 60,0.,60., ph.chargedIsolation()      , theWeight);
	theHistograms->fill("kinPh_HoverE_EB"+strPrompt, "kinPhotons in Barrel;HoverE"             , 60,0.,.15, ph.HoverE()                , theWeight);
	theHistograms->fill("kinPh_neIso_EB" +strPrompt, "kinPhotons in Barrel;neIso"              , 60,0.,60., ph.neutralHadronIsolation(), theWeight);
	theHistograms->fill("kinPh_phIso_EB" +strPrompt, "kinPhotons in Barrel;phIso"              , 60,0.,60., ph.photonIsolation()       , theWeight);
      }
    }
    else{
      vector<double> sieie_bins(28);
      for(size_t i = 0; i < sieie_bins.size(); i++) sieie_bins[i] = 0.01 + 0.0025*i;
      vector<double> chIso_bins({0., 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4}); //({0., 0.1, 0.25, 0.517, 0.8, 1.051, 1.3, 1.6, 2.089, 3, 4});
      theHistograms->fill("kinPh_sieie_chIso_EE", "kinPhotons in Endcap;#sigma_{i#etai#eta};chIso",
			  sieie_bins,
			  chIso_bins,
			  ph.sigmaIetaIeta(), ph.chargedIsolation(), theWeight);

      theHistograms->fill(  "kinPh_sieie_EE"           , "kinPhotons in Endcap;#sigma_{i#etai#eta}", sieie_bins, ph.sigmaIetaIeta()        , theWeight);
      theHistograms->fill(  "kinPh_chIso_EE"           , "kinPhotons in Endcap;chIso"              , 60,0.,60., ph.chargedIsolation()      , theWeight);
      theHistograms->fill(  "kinPh_HoverE_EE"          , "kinPhotons in Endcap;HoverE"             , 60,0.,.15, ph.HoverE()                , theWeight);
      theHistograms->fill(  "kinPh_neIso_EE"           , "kinPhotons in Endcap;neIso"              , 60,0.,60., ph.neutralHadronIsolation(), theWeight);
      theHistograms->fill(  "kinPh_phIso_EE"           , "kinPhotons in Endcap;phIso"              , 60,0.,60., ph.photonIsolation()       , theWeight);

      if(theSampleInfo.isMC()){
	theHistograms->fill("kinPh_sieie_EE" +strPrompt, "kinPhotons in Endcap;#sigma_{i#etai#eta}", sieie_bins, ph.sigmaIetaIeta()        , theWeight);
	theHistograms->fill("kinPh_chIso_EE" +strPrompt, "kinPhotons in Endcap;chIso"              , 60,0.,60., ph.chargedIsolation()      , theWeight);
	theHistograms->fill("kinPh_HoverE_EE"+strPrompt, "kinPhotons in Endcap;HoverE"             , 60,0.,.15, ph.HoverE()                , theWeight);
	theHistograms->fill("kinPh_neIso_EE" +strPrompt, "kinPhotons in Endcap;neIso"              , 60,0.,60., ph.neutralHadronIsolation(), theWeight);
	theHistograms->fill("kinPh_phIso_EE" +strPrompt, "kinPhotons in Endcap;phIso"              , 60,0.,60., ph.photonIsolation()       , theWeight);
      }
    }

    theHistograms->fill(  "kinPhotons_ID"          , "Cut Based ID", BINS_KINPHID, 0, theWeight);
    if(passVeryLoose(ph))     theHistograms->fill(  "kinPhotons_ID"          , "Cut Based ID", BINS_KINPHID, 1, theWeight);
    if(ph.cutBasedIDLoose())  theHistograms->fill(  "kinPhotons_ID"          , "Cut Based ID", BINS_KINPHID, 2, theWeight);
    if(ph.cutBasedIDMedium()) theHistograms->fill(  "kinPhotons_ID"          , "Cut Based ID", BINS_KINPHID, 3, theWeight);
    if(ph.cutBasedIDTight())  theHistograms->fill(  "kinPhotons_ID"          , "Cut Based ID", BINS_KINPHID, 4, theWeight);
    if(theSampleInfo.isMC()){
      theHistograms->fill("kinPhotons_ID"+strPrompt, "Cut Based ID", BINS_KINPHID, 0, theWeight);
      if(passVeryLoose(ph))     theHistograms->fill("kinPhotons_ID"+strPrompt, "Cut Based ID", BINS_KINPHID, 1, theWeight);
      if(ph.cutBasedIDLoose())  theHistograms->fill("kinPhotons_ID"+strPrompt, "Cut Based ID", BINS_KINPHID, 2, theWeight);
      if(ph.cutBasedIDMedium()) theHistograms->fill("kinPhotons_ID"+strPrompt, "Cut Based ID", BINS_KINPHID, 3, theWeight);
      if(ph.cutBasedIDTight())  theHistograms->fill("kinPhotons_ID"+strPrompt, "Cut Based ID", BINS_KINPHID, 4, theWeight);
    }
  }
  
  // Systematics
  for(auto & [syst, phVect] : kinPhotons_)
    if(phVect->size() > 0)
      // Creating void histograms, then filling alphanumeric labels --> new ones are created as they are encountered
      theHistograms->book<TH1F>("kinPh_eScale_N", "Number of #gamma passing selection", 1,0,0)->Fill(Form("kin_%s" , syst.c_str()), theWeight);
  
  for(auto & [syst, phVect] : goodPhotons_)
    if(phVect->size() > 0)
      theHistograms->book<TH1F>("kinPh_eScale_N", "Number of #gamma passing selection", 1,0,0)->Fill(Form("good_%s", syst.c_str()), theWeight);
  
  // How many photons are there in each event?
  for(const auto& [syst, phVect] : goodPhotons_)
    theHistograms->fill(Form("loosePh_%s_N", syst.c_str()), Form("Number of Loose photons (%s)", syst.c_str()), 5,-0.5,4.5, phVect->size(), theWeight);
  
  for(const auto& [syst, phVect] : loosePhotons_)
    theHistograms->fill(Form("veryLoosePh_%s_N", syst.c_str()), Form("Number of VeryLoose photons (%s)", syst.c_str()), 5,-0.5,4.5, phVect->size(), theWeight);

  for(const auto& [syst, phVect] : kinPhotons_)
    theHistograms->fill(Form(  "kinPh_%s_N", syst.c_str()), Form("Number of Kinematic photons (%s)"  , syst.c_str()), 5,-0.5,4.5, phVect->size(), theWeight);
  
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
    fillCutsNm1("kinPhotons_Nm1", "N-1 cut efficiency;;Events", {
      {"sieie" , b_sieie },
      {"HoverE", b_HoverE},
      {"chIso" , b_chIso },
      {"neIso" , b_neIso },
      {"phIso" , b_phIso }
      }, theWeight);
    
    theHistograms->fill(  "kinPhoton_MVA"                    , "Kin #gamma MVA"                    , 40,-1.,1., bestKinPh_->MVAvalue(), theWeight);
    theHistograms->fill(  "kinPhoton_MVA_"      +channelReco_, "Kin #gamma MVA in "   +channelReco_, 40,-1.,1., bestKinPh_->MVAvalue(), theWeight);
    if(bestKinPh_->cutBasedID(Photon::IdWp::VeryLoose)){
      theHistograms->fill("veryLoosePhoton_MVA"              , "VeryLoose #gamma MVA"              , 40,-1.,1., bestKinPh_->MVAvalue(), theWeight);
      theHistograms->fill("veryLoosePhoton_MVA_"+channelReco_, "VeryLoose #gamma MVA "+channelReco_, 40,-1.,1., bestKinPh_->MVAvalue(), theWeight);
    }
    if(bestKinPh_->cutBasedID(Photon::IdWp::Loose    )){
      theHistograms->fill("loosePhoton_MVA"                  , "Loose #gamma MVA"                  , 40,-1.,1., bestKinPh_->MVAvalue(), theWeight);
      theHistograms->fill("loosePhoton_MVA_"    +channelReco_, "Loose #gamma MVA "    +channelReco_, 40,-1.,1., bestKinPh_->MVAvalue(), theWeight);
    }
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
  bool haveKinPh   = bestKinPh_ != nullptr;
  bool haveVLPh    = bestKinPh_ && bestKinPh_->cutBasedID(Photon::IdWp::VeryLoose);
  bool haveLoosePh = bestKinPh_ && bestKinPh_->cutBasedID(Photon::IdWp::Loose    );

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
    theHistograms->fill("AK8_resolution_dR"        ,"AK8: #DeltaR(reco,gen);#DeltaR"                 , 40,  0.,0.2,physmath::deltaR(gen, rec)    );
    theHistograms->fill("AK8_resolution_pt"        ,"AK8;p_{T}^{REC}-p_{T}^{GEN} [GeV/c]"            , 80,-40.,40.,rec.pt()           -gen.pt()  );
    theHistograms->fill("AK8_resolution_corrPruned","AK8: corrPrunedMass;m_{REC}-m_{GEN} [GeV/c^{2}]", 80,-40.,40.,rec.corrPrunedMass()-gen.mass());
    theHistograms->fill("AK8_resolution_pruned"    ,"AK8: prunedMass;m_{REC}-m_{GEN} [GeV/c^{2}]"    , 80,-40.,40.,rec.prunedMass()   -gen.mass());
    theHistograms->fill("AK8_resolution_softDrop"  ,"AK8: softDropMass;m_{REC}-m_{GEN} [GeV/c^{2}]"  , 80,-40.,40.,rec.softDropMass() -gen.mass());
    theHistograms->fill("AK8_resolution_mass"      ,"AK8: mass;m_{REC}-m_{GEN} [GeV/c^{2}]"          , 80,-40.,40.,rec.mass()         -gen.mass());
  }
  
  // Jets AK4
  theHistograms->fill("AK4_N" , "# jets AK4", 6, -0.5, 5.5, jets->size(), theWeight);
  theHistograms->fill("AK4_noph_N" , "# jets AK4 (disambiguated)", 6, -0.5, 5.5, jets_noph_->size(), theWeight);
  if(haveKinPh  ) theHistograms->fill("AK4_N_KinPh"   , "# jets AK4 when there is a Kin photon"  , 6, -0.5, 5.5, jets->size(), theWeight);
  if(haveVLPh   ) theHistograms->fill("AK4_N_VLPh"    , "# jets AK4 when there is a VL photon"   , 6, -0.5, 5.5, jets->size(), theWeight);
  if(haveLoosePh) theHistograms->fill("AK4_N_LoosePh" , "# jets AK4 when there is a Loose photon", 6, -0.5, 5.5, jets->size(), theWeight);

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
    theHistograms->fill("AK4_resolution_dR"        ,"AK4: #DeltaR(reco,gen);#DeltaR"                 , 40,0.,0.2, physmath::deltaR(gen, rec)   );
    theHistograms->fill("AK4_resolution_pt"        ,"AK4;p_{T}^{REC}-p_{T}^{GEN} [GeV/c]"            , 60,-30,30, rec.pt()          -gen.pt()  );
    theHistograms->fill("AK4_resolution_mass"      ,"AK4;m_{REC}-m_{GEN} [GeV/c^{2}]"                , 60,-30,30, rec.mass()        -gen.mass());
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
    const Photon& ph = kinPhotons_["central"]->front();
    const TLorentzVector& ph_p4 = ph.p4();
    theHistograms->fill(Form("%s_%s_kinPh" , name, mType), Form("%s %s with Kin #gamma" , title, mType), binsVV , mValue(p4      ), theWeight);
    theHistograms->fill(Form("%sG_%s_kinPh", name, mType), Form("%sG %s with Kin #gamma", title, mType), binsVVG, mValue(p4+ph_p4), theWeight);

    const char* genStatus = sigdefHelper.pass_photon() ? "prompt" : "nonpro" ;
    theHistograms->fill(Form("%s_%s_kinPh_%s" , name, mType, genStatus), Form("%s %s with Kin #gamma" , title, mType), binsVV , mValue(p4      ), theWeight);
    theHistograms->fill(Form("%sG_%s_kinPh_%s", name, mType, genStatus), Form("%sG %s with Kin #gamma", title, mType), binsVVG, mValue(p4+ph_p4), theWeight);

    if(goodPhotons_["central"]->size() == 0){
      theHistograms->fill(Form("%s_%s_kinVetoL" , name, mType), Form("%s %s with Kin #gamma" , title, mType), binsVV , mValue(p4      ), theWeight);
      theHistograms->fill(Form("%sG_%s_kinVetoL", name, mType), Form("%sG %s with Kin #gamma", title, mType), binsVVG, mValue(p4+ph_p4), theWeight);
      theHistograms->fill(Form("%s_%s_kinVetoL_%s" , name, mType, genStatus), Form("%s %s with Kin #gamma" , title, mType), binsVV , mValue(p4      ), theWeight);
      theHistograms->fill(Form("%sG_%s_kinVetoL_%s", name, mType, genStatus), Form("%sG %s with Kin #gamma", title, mType), binsVVG, mValue(p4+ph_p4), theWeight);
    }
  }
  // No photon
  else
    theHistograms->fill(Form("%s_%s_noPh", name, mType), Form("%s mass with No #gamma", title), binsVV, mValue(p4), theWeight);

  // Loose photons
  if(loosePhotons_["central"]->size() > 0){
    const Photon& ph = loosePhotons_["central"]->front();
    const TLorentzVector& ph_p4 = ph.p4();
    theHistograms->fill(Form("%s_%s_veryLoosePh" , name, mType), Form("%s %s with VeryLoose #gamma" , title, mType), binsVV , mValue(p4      ), theWeight);
    theHistograms->fill(Form("%sG_%s_veryLoosePh", name, mType), Form("%sG %s with VeryLoose #gamma", title, mType), binsVVG, mValue(p4+ph_p4), theWeight);

    const char* genStatus = sigdefHelper.pass_photon() ? "prompt" : "nonpro" ;
    theHistograms->fill(Form("%s_%s_veryLoosePh_%s" , name, mType, genStatus), Form("%s %s with VeryLoose #gamma" , title, mType), binsVV , mValue(p4      ), theWeight);
    theHistograms->fill(Form("%sG_%s_veryLoosePh_%s", name, mType, genStatus), Form("%sG %s with VeryLoose #gamma", title, mType), binsVVG, mValue(p4+ph_p4), theWeight);

    // Tight photons
    if(goodPhotons_["central"]->size() > 0){
      const Photon& ph = goodPhotons_["central"]->front();
      const TLorentzVector& ph_p4 = ph.p4();
      theHistograms->fill(Form("%s_%s_loosePh" , name, mType), Form("%s %s with LooseID #gamma" , title, mType), binsVV , mValue(p4      ), theWeight);
      theHistograms->fill(Form("%sG_%s_loosePh", name, mType), Form("%sG %s with LooseID #gamma", title, mType), binsVVG, mValue(p4+ph_p4), theWeight);

      theHistograms->fill(Form("%s_%s_loosePh_%s" , name, mType, genStatus), Form("%s %s with LooseID #gamma" , title, mType), binsVV , mValue(p4      ), theWeight);
      theHistograms->fill(Form("%sG_%s_loosePh_%s", name, mType, genStatus), Form("%sG %s with LooseID #gamma", title, mType), binsVVG, mValue(p4+ph_p4), theWeight);
    }
    // Fail photon (loose && !tight)
    else{
      float f_VLtoL = getPhotonFR_VLtoL(ph);
      float w_VLtoL = f_VLtoL / (1 - f_VLtoL);

      theHistograms->fill(Form("%s_%s_failPh" , name, mType), Form("%s %s with Fail #gamma" , title, mType), binsVV , mValue(p4      ), theWeight);
      theHistograms->fill(Form("%sG_%s_failPh", name, mType), Form("%sG %s with Fail #gamma", title, mType), binsVVG, mValue(p4+ph_p4), theWeight);
      theHistograms->fill(Form("%s_%s_reweightPh" , name, mType), Form("%s %s with Reweighted #gamma" , title, mType), binsVV , mValue(p4      ), theWeight * w_VLtoL);
      theHistograms->fill(Form("%sG_%s_reweightPh", name, mType), Form("%sG %s with Reweighted #gamma", title, mType), binsVVG, mValue(p4+ph_p4), theWeight * w_VLtoL);

      const char* genStatus = sigdefHelper.pass_photon() ? "prompt" : "nonpro" ;
      theHistograms->fill(Form("%s_%s_failPh_%s" , name, mType, genStatus), Form("%s %s with Fail #gamma" , title, mType), binsVV , mValue(p4      ), theWeight);
      theHistograms->fill(Form("%sG_%s_failPh_%s", name, mType, genStatus), Form("%sG %s with Fail #gamma", title, mType), binsVVG, mValue(p4+ph_p4), theWeight);
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


void VVGammaAnalyzer::photonFakeRate_ABCD(){
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
    isPrompt = sigdefHelper.pass_photon();
  // std::any_of(genPhotonsPrompt_->begin(), genPhotonsPrompt_->end(), [bestG](const Particle& gen){ return physmath::deltaR(*bestG, gen) < 0.2; });
  
  // Debug photon ID
  // debugPhotonID(*bestG);
  
  // Actual fake rate histograms
  double theAeta = fabs(bestG->eta());
  double thePt   = bestG->pt();
  if(thePt > ph_pt_bins.back())
    thePt = ph_pt_bins.back() - 0.1;
  
  // Fake rate with ABCD
  char ABCD = phABCD(*bestG, wp);

  const char* name;
  if(theSampleInfo.isMC())
    name = Form("PhFR_%c_%s", ABCD, (isPrompt ? "prompt" : "nonprompt"));
  else
    name = Form("PhFR_%c_%s", ABCD, "data");
  
  theHistograms->fill(name, Form("Photons: %c;p_{T} [GeV/c];#eta", ABCD), ph_pt_bins, ph_aeta_bins, thePt, theAeta, theWeight);
}


void VVGammaAnalyzer::photonFakeRate_LtoT(const char* method, const Photon& thePh, bool isPass, double effSF){
  double theAeta = fabs(thePh.eta());
  double thePt   = thePh.pt();
  if(thePt > ph_pt_bins.back())
    thePt = ph_pt_bins.back() - 0.1;

  char phFSRch = canBeFSR(thePh, *leptons_) ? 'Y' : 'N';

  char channel[8]; int e;
  if(region_ == CRLFR){
    bool lepPass = ZL->second.passFullSel();
    e = snprintf(channel, 8, "%s%c-%c", channelReco_.c_str(), lepPass ? 'P' : 'F', phFSRch);
  }
  else
    e = snprintf(channel, 8, "%s-%c", channelReco_.c_str(), phFSRch);
  if(!(e >= 0 && e < 8)) std::cerr << "WARN: problem encoding channel string\n";

  const char* strPass = isPass ? "PASS" : "FAIL";

  const char* strPrompt = "data";
  if(theSampleInfo.isMC()){
    bool isPrompt = sigdefHelper.pass_photon();
    strPrompt = isPrompt ? "prompt" : "nonprompt";
  }

  char phEtaRegion[4]; sprintf(phEtaRegion, thePh.isBarrel() ? "EB" : "EE");

  // Closest lep
  std::vector<Lepton>::const_iterator closestLep = closestDeltaR(thePh, *leptons_);
  float dR_l = physmath::deltaR(*closestLep, thePh);

  // char lepFlavour = '?';
  // unsigned int aLepID = abs(closestLep->id());
  // if     (aLepID == 11) lepFlavour = 'e';
  // else if(aLepID == 13) lepFlavour = 'm';

  // char lepStatus = closestLep->passFullSel() ? 'P' : 'F';

  // char all_str[4]; sprintf(all_str, "all");
  // char all_char = 'a';
  static vector<double> edges_dR {0., 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00};
  if(dR_l > edges_dR.back()) dR_l = edges_dR.back() - 0.001;

  unsigned int njets = jets_noph_->size();
  static vector<double> edges_njets {-0.5, 0.5, 1.5, 2.5, 3.5};
  float dR_j = edges_dR.back() - 0.001;
  if(njets > 0)
    dR_j = physmath::deltaR(*closestDeltaR(thePh, *jets_noph_), thePh);

  // Fill photon FR plots
  const char* name_aeta_inclusive = Form("PhFR_%s_pt-aeta_%s_%s"   , method,          strPrompt, strPass);
  theHistograms->fill(name_aeta_inclusive, "Photon fake rate VeryLoose to Loose;p_{T} [GeV/c];#eta;Events"             , ph_pt_bins, ph_aeta_bins, thePt, theAeta, theWeight*effSF);
  const char* name_aeta_channel   = Form("PhFR_%s_pt-aeta_%s_%s_%s", method, channel, strPrompt, strPass);
  theHistograms->fill(name_aeta_channel  , "Photon fake rate VeryLoose to Loose;p_{T} [GeV/c];#eta;Events"             , ph_pt_bins, ph_aeta_bins, thePt, theAeta, theWeight*effSF);

  const char* name_dRl_inclusive  = Form("PhFR_%s_pt-dRl_%s_%s"   , method,          strPrompt, strPass);
  theHistograms->fill(name_dRl_inclusive ,"Photon fake rate VeryLoose to Loose;p_{T} [GeV/c];#DeltaR(#gamma, l);Events", ph_pt_bins, edges_dR    , thePt, dR_l   , theWeight*effSF);
  const char* name_dRl_channel    = Form("PhFR_%s_pt-dRl_%s_%s_%s", method, channel, strPrompt, strPass);
  theHistograms->fill(name_dRl_channel   ,"Photon fake rate VeryLoose to Loose;p_{T} [GeV/c];#DeltaR(#gamma, l);Events", ph_pt_bins, edges_dR    , thePt, dR_l   , theWeight*effSF);

  const char* name_njets_inclusive= Form("PhFR_%s_%s_pt-njets_%s_%s"   , method, phEtaRegion,          strPrompt, strPass);
  theHistograms->fill(name_njets_inclusive,"Photon fake rate VeryLoose to Loose;p_{T} [GeV/c];# jets);Events"          , ph_pt_bins, edges_njets , thePt, njets  , theWeight*effSF);
  const char* name_njets_channel  = Form("PhFR_%s_%s_pt-njets_%s_%s_%s", method, phEtaRegion, channel, strPrompt, strPass);
  theHistograms->fill(name_njets_channel ,"Photon fake rate VeryLoose to Loose;p_{T} [GeV/c];# jets);Events"           , ph_pt_bins, edges_njets , thePt, njets  , theWeight*effSF);

  const char* name_dRj_inclusive  = Form("PhFR_%s_%s_pt-dRj_%s_%s"   , method, phEtaRegion,          strPrompt, strPass);
  theHistograms->fill(name_dRj_inclusive ,"Photon fake rate VeryLoose to Loose;p_{T} [GeV/c];#DeltaR(#gamma, j);Events", ph_pt_bins, edges_dR    , thePt, dR_j   , theWeight*effSF);
  const char* name_dRj_channel    = Form("PhFR_%s_%s_pt-dRj_%s_%s_%s", method, phEtaRegion, channel, strPrompt, strPass);
  theHistograms->fill(name_dRj_channel   ,"Photon fake rate VeryLoose to Loose;p_{T} [GeV/c];#DeltaR(#gamma, j);Events", ph_pt_bins, edges_dR    , thePt, dR_j   , theWeight*effSF);

  // for(char lepSt : {all_char, lepStatus}){
  //   for(char lepFl : {all_char, lepFlavour}){
  //     for(char* phEta : {all_str, phEtaRegion}){
  // 	for(char phFSR : {all_char, phFSRch}){
  // 	  const char* name_byChannel = Form("PhFR_%s_pt-dRl_%c-%c-%s-%c_%s_%s" ,
  // 					    method,
  // 					    lepSt,
  // 					    lepFl,
  // 					    phEta,
  // 					    phFSR,
  // 					    strPrompt,
  // 					    strPass
  // 					    );
  // 	  theHistograms->fill(name_byChannel, Form("FR #gamma %s;p_{T}^{#gamma};#DeltaR(#gamma, l);Events", method), ph_pt_bins, edges_dR, thePt, dR_l, theWeight*effSF);
  // 	}
  //     }
  //   }
  // }
}


void VVGammaAnalyzer::photonFRClosure(const char* method, const Photon& thePh, bool isPass, double f_FR){
  // TEMP: this has to be moved to a separate function or computed once per event, without using the photon
  char phFSRch = canBeFSR(thePh, *leptons_) ? 'Y' : 'N';

  char channel[8]; int e;
  if(region_ == CRLFR){
    bool lepPass = ZL->second.passFullSel();
    e = snprintf(channel, 8, "%s%c-%c", channelReco_.c_str(), lepPass ? 'P' : 'F', phFSRch);
  }
  else
    e = snprintf(channel, 8, "%s-%c", channelReco_.c_str(), phFSRch);
  if(!(e >= 0 && e < 8)) std::cerr << "WARN: problem encoding channel string\n";

  // ##### Closure tests #####
  const char* varName;
  double varValue;
  if     (is4Lregion(region_)){
    varName = "mZZG";
    varValue = ( ZZ->p4() + thePh.p4() ).M();
  }
  else if(is3Lregion(region_)){
    varName = "mWZG";
    varValue = ( ZW->p4() + thePh.p4() ).Mt();
  }
  else if(is2Lregion(region_)){
    varName = "mZG";
    varValue = ( Z->p4() + thePh.p4() ).M();
  }
  else
    return;

  // e.g. method = VLtoL_pt-aeta_data
  double weight_FR = theWeight * f_FR/(1-f_FR);
  if(isPass){
    theHistograms->fill(Form("PhFRClosure_%s_PASS_%s"         , method ,          varName), Form("Closure test %s: PASS"   , method), mVVG_bins, varValue, theWeight);
    theHistograms->fill(Form("PhFRClosure_%s_%s_PASS_%s"      , method , channel, varName), Form("Closure test %s: PASS"   , method), mVVG_bins, varValue, theWeight);
  }
  else{
    theHistograms->fill(Form("PhFRClosure_%s_FAIL_%s"         , method ,          varName), Form("Closure test %s: FAIL"   , method), mVVG_bins, varValue, theWeight);
    theHistograms->fill(Form("PhFRClosure_%s_reweighted_%s"   , method ,          varName), Form("Closure test %s: FAIL*TF", method), mVVG_bins, varValue, weight_FR);
    theHistograms->fill(Form("PhFRClosure_%s_%s_FAIL_%s"      , method , channel, varName), Form("Closure test %s: FAIL"   , method), mVVG_bins, varValue, theWeight);
    theHistograms->fill(Form("PhFRClosure_%s_%s_reweighted_%s", method , channel, varName), Form("Closure test %s: FAIL*TF", method), mVVG_bins, varValue, weight_FR);
  }
}


void VVGammaAnalyzer::studyJetsChoice(){
  // Study how to choose the correct jet AK8 in case it exists. The problem of rejecting 
  // events in which it does not exist is part of the event selection
  if( !theSampleInfo.isMC() )
    return;
  // Boson<Particle>* qq = nullptr;
  // if     (genZhadCandidates_->size() > 0)
  //   qq = &(genZhadCandidates_->at(0));
  // else if(genWhadCandidates_->size() > 0)
  //   qq = &(genWhadCandidates_->at(0));
  // else
  //   return;
  
  // int statusAK8;
  // int statusAK4 = studyAK4Choice(fAK4_, *qq, 0.4);
  // if(statusAK4 != 0)
  //   // statusAK8 = studyAK8Choice(fAK8_, *qq, 0.4);
  //   studyAK8Choice(fAK8_, *qq, 0.4);
  
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

void VVGammaAnalyzer::studyFSRregion(const vector<Photon>& fsrMatched){
  /*
    Study the overalp between ZZTo4l and ZZGTo4LG in the FSR region (dRl < 0.5)
   */
  if(fsrMatched.size() == 0)
    return;

  // Kin
  if(kinPhotons_["central"]->size() == 0){ // Exclude events which may be in the signal region
    const Photon& ph = fsrMatched[0];
    fillPhotonPlots(ph, "lead_FSRkin", "FSR (kin) veto others");
  }

  if(loosePhotons_["central"]->size() == 0){ // Exclude events which may be in the signal region
    const Photon& ph = fsrMatched[0];
    fillPhotonPlots(ph, "lead_FSRloose", "FSR (loose) veto others");
  }
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
  vector<Photon> isolatedPh;

  for(const Photon& ph : vPh){
    const std::vector<Lepton>::iterator closestLep = std::min_element(leptons_->begin(), leptons_->end(),
								      [ph](const Lepton& a, const Lepton& b){
									return physmath::deltaR(ph, a) < physmath::deltaR(ph, b);
								      } );
    double minL_DR = closestLep != leptons_->end() ? deltaR(ph, *closestLep) : 10;

    theHistograms->fill(Form("minL_DR_%s", label), Form("min_{l}(#DeltaR(#gamma, l)) for each #gamma %s;#DeltaR;# of #gamma", label), 50, 0.,1., minL_DR, theWeight);
    maxG_minL_DR = std::max(maxG_minL_DR, minL_DR);
  }
  if(maxG_minL_DR > 0)
    theHistograms->fill(Form("maxG_minL_DR_%s", label), Form("max_{#gamma}(min_{l}(#DeltaR(#gamma, l))) %s;#DeltaR;Events",label), 50, 0.,1., maxG_minL_DR, theWeight);
}

void VVGammaAnalyzer::photonIsolation_bestKin(){
  if(!bestKinPh_)
    return;

  // Photon ID: Kin, VeryLooseL, Loose
  bool isPassVL = bestKinPh_->cutBasedID(Photon::IdWp::VeryLoose);
  bool isPassLoose = bestKinPh_->cutBasedIDLoose();
  char phStatus[8] = "Kin";
  if     (isPassLoose) sprintf(phStatus, "Pass");
  else if(isPassVL   ) sprintf(phStatus, "Fail");

  // Photon eta and pt
  char phEtaRegion[4]; sprintf(phEtaRegion, bestKinPh_->isBarrel() ? "EB" : "EE");
  float phPtValue = bestKinPh_->pt();
  char phPtRegion[8] = "20-35";
  if     (phPtValue > 50) sprintf(phPtRegion, "50-"  );
  else if(phPtValue > 35) sprintf(phPtRegion, "35-50");

  // Closest lep
  std::vector<Lepton>::const_iterator closestLep = closestDeltaR(*bestKinPh_, *leptons_);
  float dR_l = physmath::deltaR(*closestLep, *bestKinPh_);

  // Closest Z
  const Boson<Lepton> *closestZ = nullptr;
  if(is4Lregion(region_)){
    vector<Boson<Lepton>*> Zs {ZZ->firstPtr(), ZZ->secondPtr()};
    closestZ = *closestDeltaR_p(*bestKinPh_, Zs);
  }
  else if(is3Lregion(region_)){
    closestZ = ZW->firstPtr();
  }
  else if(is2Lregion(region_)){
    closestZ = Z;
  }
  else if(region_ == CRLFR){
    closestZ = &ZL->first;
  }
  float dM_Z = closestZ->mass() - phys::ZMASS;
  float dR_Z = physmath::deltaR(*closestZ, *bestKinPh_);

  // Fill the plots
  char all_str[4]; sprintf(all_str, "all");
  vector<double> edges_dR {0., 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00};
  for(auto phEta : {all_str, phEtaRegion}){
    for(auto phSt : {all_str, phStatus}){
      for(auto phPt : {all_str, phPtRegion}){
	theHistograms->fill(Form("photon_dRl_dMZ_%s_%s_%s", phSt, phPt, phEta), ";#DeltaR(#gamma, l);m_{Z}-m_{Z0};Events", 20,0.,1., 30,-20.,10., dR_l, dM_Z, theWeight);
	theHistograms->fill(Form("photon_dRZ_dMZ_%s_%s_%s", phSt, phPt, phEta), ";#DeltaR(#gamma, Z);m_{Z}-m_{Z0};Events", 20,0.,1., 30,-20.,10., dR_Z, dM_Z, theWeight);
      }
      theHistograms->fill(Form("photon_dRl_pt_%s_%s", phSt, phEta), ";#DeltaR(#gamma, l);p_{T}^{#gamma};Events", edges_dR, ph_pt_bins, dR_l, phPtValue, theWeight);
      theHistograms->fill(Form("photon_dRZ_pt_%s_%s", phSt, phEta), ";#DeltaR(#gamma, Z);p_{T}^{#gamma};Events", edges_dR, ph_pt_bins, dR_Z, phPtValue, theWeight);
    }
  }
}


void VVGammaAnalyzer::orphanPhotonStudy(){
  if(!theSampleInfo.isMC() || !bestKinPh_)
    return;

  // One photon per event
  bool isPassVL    = bestKinPh_->cutBasedID(Photon::IdWp::VeryLoose);
  bool isPassLoose = bestKinPh_->cutBasedID(Photon::IdWp::Loose    );

  vector<Particle>::const_iterator gen_it;
  gen_it = closestDeltaR(*bestKinPh_, *genPhotons_);

  float dR;
  if(gen_it != genPhotons_->cend() && (dR = deltaR(*gen_it, *bestKinPh_)) < 0.1){
    float pr = bestKinPh_->p() / gen_it->p();  // p-ratio
    theHistograms->fill("orphanKin_22_dR", "bestKinPh matched to genPhoton;#DeltaR(gen, #gamma);Events", 40,0.,0.4, dR, theWeight);
    theHistograms->fill("orphanKin_22_pr", "bestKinPh matched to genPhoton;p^{#gamma}/p^{GEN};Events"  , 40,0.,2. , pr, theWeight);
    if(isPassVL){
      theHistograms->fill("orphanVL_22_dR", "bestVLPh matched to genPhoton;#DeltaR(gen, #gamma);Events", 40,0.,0.4, dR, theWeight);
      theHistograms->fill("orphanVL_22_dR", "bestVLPh matched to genPhoton;p^{#gamma}/p^{GEN};Events"  , 40,0.,2. , pr, theWeight);
    }
    if(isPassLoose){
      theHistograms->fill("orphanLoose_22_dR", "bestLoosePh matched to genPhoton;#DeltaR(gen, #gamma);Events", 40,0.,0.4, dR, theWeight);
      theHistograms->fill("orphanLoose_22_pr", "bestLoosePh matched to genPhoton;p^{#gamma}/p^{GEN};Events"  , 40,0.,2. , pr, theWeight);
    }
    return;
  }

  // Any gen particle
  gen_it = closestDeltaR(*bestKinPh_, *genParticles);
  if(gen_it != genParticles->cend() && (dR = deltaR(*gen_it, *bestKinPh_)) < 0.4){
    float pr = bestKinPh_->p() / gen_it->p();  // p-ratio
    theHistograms->fill(Form("orphanKin_%d_dR", abs(gen_it->id())), "bestKinPh matched with genParticle;#DeltaR(gen, #gamma);Events", 40,0.,0.4, dR, theWeight);
    theHistograms->fill(Form("orphanKin_%d_pr", abs(gen_it->id())), "bestKinPh matched with genParticle;p^{#gamma}/p^{GEN};Events"  , 40,0.,2. , pr, theWeight);
    if(isPassVL){
      theHistograms->fill(Form("orphanVL_%d_dR", abs(gen_it->id())), "bestVLPh matched with genParticle;#DeltaR(gen, #gamma);Events", 40,0.,0.4, dR, theWeight);
      theHistograms->fill(Form("orphanVL_%d_pr", abs(gen_it->id())), "bestVLPh matched with genParticle;p^{#gamma}/p^{GEN};Events"  , 40,0.,2. , pr, theWeight);
    }
    if(isPassLoose){
      theHistograms->fill(Form("orphanLoose_%d_dR", abs(gen_it->id())), "bestLoosePh matched with genParticle;#DeltaR(gen, #gamma);Events", 40,0.,0.4, dR, theWeight);
      theHistograms->fill(Form("orphanLoose_%d_pr", abs(gen_it->id())), "bestLoosePh matched with genParticle;p^{#gamma}/p^{GEN};Events"  , 40,0.,2. , pr, theWeight);
    }
    return;
  }

  // GenJets
  gen_it = closestDeltaR(*bestKinPh_, *genJets);
  if(gen_it != genJets->cend() && (dR = deltaR(*gen_it, *bestKinPh_)) < 0.4){
    float pr = bestKinPh_->p() / gen_it->p();  // p-ratio
    theHistograms->fill(Form("orphanKin_%d_dR", abs(gen_it->id())), "bestKinPh matched with genJet;#DeltaR(gen, #gamma);Events", 40,0.,0.4, dR, theWeight);
    theHistograms->fill(Form("orphanKin_%d_pr", abs(gen_it->id())), "bestKinPh matched with genJet;p^{#gamma}/p^{GEN};Events"  , 40,0.,2.0, pr, theWeight);
    if(isPassVL){
      theHistograms->fill(Form("orphanVL_%d_dR", abs(gen_it->id())), "bestVLPh matched with genJet;#DeltaR(gen, #gamma);Events", 40,0.,0.4, dR, theWeight);
      theHistograms->fill(Form("orphanVL_%d_pr", abs(gen_it->id())), "bestVLPh matched with genJet;p^{#gamma}/p^{GEN};Events"  , 40,0.,2.0, pr, theWeight);
    }
    if(isPassLoose){
      theHistograms->fill(Form("orphanLoose_%d_dR", abs(gen_it->id())), "bestLoosePh matched with genJet;#DeltaR(gen, #gamma);Events", 40,0.,0.4, dR, theWeight);
      theHistograms->fill(Form("orphanLoose_%d_pr", abs(gen_it->id())), "bestLoosePh matched with genJet;p^{#gamma}/p^{GEN};Events"  , 40,0.,2.0, pr, theWeight);
    }
    return;
  }

  // If we're still here, these photons are not associated with anything; are they real?
  auto closestLep = closestDeltaR(*bestKinPh_, *leptons_);
  float dRlep = deltaR(*closestLep, *bestKinPh_);
  theHistograms->fill("orphanKin_unmatched_pt"   , "bestKinPh unmatched;p_{T}^#gamma;Events"       , 40,20.,180., bestKinPh_->pt() , theWeight);
  theHistograms->fill("orphanKin_unmatched_eta"  , "bestKinPh unmatched;#eta^#gamma;Events"        , 40,-CUT_G_AETA_MAX,CUT_G_AETA_MAX, bestKinPh_->eta(), theWeight);
  theHistograms->fill("orphanKin_unmatched_dRlep", "bestKinPh unmatched;#DeltaR(#gamma, l);Events" , 40,0.,1.   , dRlep            , theWeight);
  if(isPassVL){
      theHistograms->fill("orphanVL_unmatched_pt"   , "bestVLPh unmatched;p_{T}^#gamma;Events"       , 40,20.,180., bestKinPh_->pt() , theWeight);
      theHistograms->fill("orphanVL_unmatched_eta"  , "bestVLPh unmatched;#eta^#gamma;Events"        , 40,-CUT_G_AETA_MAX,CUT_G_AETA_MAX, bestKinPh_->eta(), theWeight);
      theHistograms->fill("orphanVL_unmatched_dRlep", "bestVLPh unmatched;#DeltaR(#gamma, l);Events" , 40,0.,1.   , dRlep            , theWeight);
  }
  if(isPassLoose){
      theHistograms->fill("orphanLoose_unmatched_pt"   , "bestLoosePh unmatched;p_{T}^#gamma;Events"       , 40,20.,180., bestKinPh_->pt() , theWeight);
      theHistograms->fill("orphanLoose_unmatched_eta"  , "bestLoosePh unmatched;#eta^#gamma;Events"        , 40,-CUT_G_AETA_MAX,CUT_G_AETA_MAX, bestKinPh_->eta(), theWeight);
      theHistograms->fill("orphanLoose_unmatched_dRlep", "bestLoosePh unmatched;#DeltaR(#gamma, l);Events" , 40,0.,1.   , dRlep            , theWeight);
  }
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
  auto r = cherryEvents.find(run);
  if(r == cherryEvents.end()) return false;

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


std::pair<double, double> VVGammaAnalyzer::getZllAndZllgMasses(const phys::Photon& ph){
  vector<Boson<Lepton>*> tmp {ZZ->firstPtr(), ZZ->secondPtr()};
  const Boson<Lepton>* closestZ = *std::min_element(tmp.begin(), tmp.end(),
                                             [&ph](const Boson<Lepton>* pZ1, const Boson<Lepton>* pZ2){
						double minZ1 = std::min(deltaR(pZ1->daughter(0), ph), deltaR(pZ1->daughter(1), ph));
						double minZ2 = std::min(deltaR(pZ2->daughter(0), ph), deltaR(pZ2->daughter(1), ph));
						return minZ1 < minZ2;
					      });
  // const Boson<Lepton> *furthestZ = (closestZ == ZZ->firstPtr() ? ZZ->secondPtr() : ZZ->firstPtr());
  double Zll_mass = closestZ->mass();
  double ZllG_mass = (closestZ->p4() + ph.p4()).M();
  return std::make_pair(Zll_mass, ZllG_mass);
}


std::pair<double, double> VVGammaAnalyzer::getZllAndZllgMasses(const std::vector<phys::Photon>& phVect){
  if(phVect.size() == 0)
    return std::make_pair(-1., -1.);

  const Photon* closestPhoton = nullptr;
  const Boson<Lepton>* closestZ = nullptr;
  float dRmin(10.);
  for(auto ZbosonP : {ZZ->firstPtr(), ZZ->secondPtr()}){
    for(auto lepton : {ZbosonP->daughter(0) , ZbosonP->daughter(1) }){
      auto closestPh = std::min_element(phVect.begin(), phVect.end(), DeltaRComparator(lepton));
      float dR = deltaR(*closestPh, lepton);
      if(dR < dRmin){
	dRmin = dR;
	closestPhoton = &*closestPh;
	closestZ = ZbosonP;
      }
    }
  }

  if(closestPhoton){
    // const Boson<Lepton> *furthestZ = (closestZ == ZZ->firstPtr() ? ZZ->secondPtr() : ZZ->firstPtr());
    double Zll_mass = closestZ->mass();
    double ZllG_mass = (closestZ->p4() + closestPhoton->p4()).M();
    return std::make_pair(Zll_mass, ZllG_mass);
  }
  else
    return std::make_pair(-2., -2.);
}


void VVGammaAnalyzer::SYSplots_inclusive(const char* syst, double weight){
  const char* strPrompt = "";
  if(theSampleInfo.isMC())
    strPrompt = sigdefHelper.pass_photon() ? "prompt" : "nonpro" ;

  if(is4Lregion(region_)){
    theHistograms->fill(  Form("SYS_mZZ_%s" , syst), Form("m_{ZZ} %s"      , syst), mVV_bins , ZZ->mass()               , weight);
    if(theSampleInfo.isMC())
      theHistograms->fill(Form("SYS_mZZ_%s_%s", syst, strPrompt), Form("m_{ZZ} %s", syst), mVV_bins , ZZ->mass()        , weight);
  }
  else if(is3Lregion(region_)){
    theHistograms->fill(  Form("SYS_mWZ_%s" , syst), Form("m_{WZ} %s"      , syst), mVV_bins , ZW->mass()               , weight);
    if(theSampleInfo.isMC())
      theHistograms->fill(Form("SYS_mWZ_%s_%s", syst, strPrompt), Form("m_{WZ} %s", syst), mVV_bins , ZW->mass()        , weight);
  }
  else if(region_ == CRLFR){
    Boson<Lepton>& theZ = ZL->first;
    Lepton&        theL = ZL->second;
    theHistograms->fill(  Form("SYS_mZ_%s"  , syst), Form("m_{Z} %s"       , syst), mZ_bins  , theZ.mass()               , weight);
    theHistograms->fill(  Form("SYS_mZL_%s" , syst), Form("m_{ZL} %s"      , syst), mZG_bins , (theZ.p4()+theL.p4()).M() , weight);
    if(theSampleInfo.isMC()){
      theHistograms->fill(Form("SYS_mZ_%s_%s"  , syst, strPrompt), Form("m_{Z} %s" , syst), mZ_bins , theZ.mass()              , weight);
      theHistograms->fill(Form("SYS_mZL_%s_%s" , syst, strPrompt), Form("m_{ZL} %s", syst), mZG_bins, (theZ.p4()+theL.p4()).M(), weight);
    }
  }
}


void VVGammaAnalyzer::SYSplots_photon(const char* syst, double weight, const Photon& ph, const char* ph_selection){
  const char* phGenStatus = "";
  if(theSampleInfo.isMC())
    phGenStatus = sigdefHelper.pass_photon() ? "prompt" : "nonpro" ;  // it is actually the event status (pass/fail the signal definition)

  theHistograms->fill(  Form("SYS_%sMVA_%s"   , ph_selection             , syst), Form("MVA %s %s"   , ph_selection             , syst), 40,-1,1   , ph.MVAvalue(),weight);
  theHistograms->fill(  Form("SYS_%spt_%s"    , ph_selection             , syst), Form("pt %s %s"    , ph_selection             , syst), ph_pt_bins, ph.pt()      ,weight);
  if(theSampleInfo.isMC()){
    theHistograms->fill(Form("SYS_%sMVA-%s_%s", ph_selection, phGenStatus, syst), Form("MVA %s %s %s", ph_selection, phGenStatus, syst), 40,-1,1   , ph.MVAvalue(),weight);
    theHistograms->fill(Form("SYS_%spt-%s_%s" , ph_selection, phGenStatus, syst), Form("pt %s %s %s" , ph_selection, phGenStatus, syst), ph_pt_bins, ph.pt()      ,weight);
  }

  // We do not care about the distrbution of mV(V)G with these photon selections
  if(strcmp(ph_selection, "kin") == 0 ||
     strcmp(ph_selection, "veryLoose") == 0)
    return;

  if(is4Lregion(region_)){
    double mZZG = (ZZ->p4() + ph.p4()).M();
    theHistograms->fill(  Form("SYS_mZZG%s_%s"   , ph_selection             , syst), Form("m_{ZZ#gamma %s} %s", ph_selection, syst), mVVG_bins, mZZG, weight);
    if(theSampleInfo.isMC())
      theHistograms->fill(Form("SYS_mZZG%s-%s_%s", ph_selection, phGenStatus, syst), Form("m_{ZZ#gamma %s} %s", ph_selection, syst), mVVG_bins, mZZG, weight);

    double Zll_mass(0.), ZllG_mass(0.);
    std::tie(Zll_mass, ZllG_mass) = getZllAndZllgMasses(ph);
    theHistograms->fill(  Form("SYS_mZllplusZllG%s_%s"   , ph_selection             , syst), Form("m_{ZZ} + m_{ZZ#gamma %s} %s", ph_selection, syst), 10,120.,320., Zll_mass+ZllG_mass, weight);
    if(theSampleInfo.isMC())
      theHistograms->fill(Form("SYS_mZllplusZllG%s-%s_%s", ph_selection, phGenStatus, syst), Form("m_{ZZ} + m_{ZZ#gamma %s} %s", ph_selection, syst), 10,120.,320., Zll_mass+ZllG_mass, weight);
  }

  else if(is3Lregion(region_)){
    double mtWZG = (ZW->p4() + ph.p4()).Mt();
    theHistograms->fill(  Form("SYS_mWZG%s_%s"   , ph_selection             , syst), Form("m_{T}^{WZ#gamma %s}; %s", ph_selection, syst), mVVG_bins, mtWZG, weight);
    if(theSampleInfo.isMC())
      theHistograms->fill(Form("SYS_mWZG%s-%s_%s", ph_selection, phGenStatus, syst), Form("m_{T}^{WZ#gamma %s}; %s", ph_selection, syst), mVVG_bins, mtWZG, weight);
  }

  else if(region_ == CRLFR){
    Boson<Lepton>& theZ = ZL->first;
    Lepton&        theL = ZL->second;
    double mZG  = (theZ.p4()             + ph.p4()).M();
    double mZLG = (theZ.p4() + theL.p4() + ph.p4()).M();
    theHistograms->fill(Form("SYS_mZG%s_%s" , ph_selection, syst), Form("m_{Z#gamma %s} %s" , ph_selection, syst), mZG_bins, mZG , weight);
    theHistograms->fill(Form("SYS_mZLG%s_%s", ph_selection, syst), Form("m_{ZL#gamma %s} %s", ph_selection, syst), mZG_bins, mZLG, weight);
    if(theSampleInfo.isMC()){
      theHistograms->fill(Form("SYS_mZG%s-%s_%s" , ph_selection, phGenStatus, syst), Form("m_{Z#gamma %s} %s %s" , ph_selection, phGenStatus, syst), mZG_bins, mZG , weight);
      theHistograms->fill(Form("SYS_mZLG%s-%s_%s", ph_selection, phGenStatus, syst), Form("m_{ZL#gamma %s} %s %s", ph_selection, phGenStatus, syst), mZG_bins, mZLG, weight);
    }
  }
}


void VVGammaAnalyzer::SYSplots_phCut(const char* syst, double weight, const Photon& phCut){
  SYSplots_photon(syst, weight, phCut, "kin");

  if(phCut.cutBasedID(Photon::IdWp::VeryLoose)){
    SYSplots_photon(syst, weight, phCut, "veryLoose");
    if(phCut.cutBasedIDLoose())
      SYSplots_photon(syst, weight, phCut, "loose");
    else{
      SYSplots_photon(syst, weight, phCut, "fail");

      double f_VLtoL = getPhotonFR_VLtoL(phCut);
      double w_VLtoL = f_VLtoL / (1 - f_VLtoL);
      SYSplots_photon(syst, weight * w_VLtoL, phCut, "failReweight");
    }
  }
}


void VVGammaAnalyzer::SYSplots_phMVA(const char* syst, double weight, const Photon& phMVA){
  // Note: unlinke the plots for phCut, there's no category for kin photons to avoid duplication

  double effSF = 1.;
  std::string MVACut_s("none");

  if  (phMVA.passMVA(Photon::MVAwp::wp90)){
    MVACut_s = "wp90";
    if(theSampleInfo.isMC())
      effSF = getPhotonEffSF_MVA(phMVA, Photon::MVAwp::wp90);

    SYSplots_photon(syst, weight * effSF, phMVA, "wp90");

    if  (phMVA.passMVA(Photon::MVAwp::wp80)){
      MVACut_s = "wp80";
      if(theSampleInfo.isMC())
        effSF = getPhotonEffSF_MVA(phMVA, Photon::MVAwp::wp80);

      SYSplots_photon(syst, weight * effSF, phMVA, "wp80");
    }
    else{
      SYSplots_photon(syst, weight * effSF, phMVA, "90not80");

      // Here the reweight when we will have the fake rate transfer factors
    }
  }

  theHistograms->fill(  Form("SYS_MVAcut_%s"                , syst), Form("MVAcut %s"               ,syst), {"none","wp90","wp80"}, MVACut_s.c_str(), weight*effSF);
  if(theSampleInfo.isMC()){
    const char* phGenStatus;
    phGenStatus = sigdefHelper.pass_photon() ? "prompt" : "nonpro" ;  // it is actually the event status (pass/fail the signal definition)
    theHistograms->fill(Form("SYS_MVAcut-%s_%s", phGenStatus, syst), Form("MVAcut %s %s",phGenStatus,syst), {"none","wp90","wp80"}, MVACut_s.c_str(), weight*effSF);
  }
}


void VVGammaAnalyzer::SYSplots(const char* syst, double weight, const Photon* phCut, const Photon* phMVA){
  SYSplots_inclusive(syst, weight);

  if(phCut)
    SYSplots_phCut(syst, weight, *phCut);

  if(phMVA)
    SYSplots_phMVA(syst, weight, *phMVA);
}


void VVGammaAnalyzer::systematicsStudy(){
  double base_w = theWeight;

  const Photon* ph = nullptr;
  if(goodPhotons_["central"]->size() >= 1)
    ph = & (goodPhotons_["central"]->front());
  else if(loosePhotons_["central"]->size() >= 1)
    ph = & (loosePhotons_["central"]->front());
  else if(kinPhotons_["central"]->size() >= 1)
    ph = & (kinPhotons_["central"]->front());

  // central
  SYSplots("central", base_w, ph, bestMVAPh_);

  // Photons energy scale and resolution
  for(const char* syst : photonSystKeys_){
    if(strcmp(syst, "central") == 0) continue;

    std::string syst_name = Form("ph%s", syst);

    SYSplots_inclusive( syst_name.c_str(), base_w);

    if(kinPhotons_[syst]->size() > 0)
      SYSplots_photon(  syst_name.c_str(), base_w, kinPhotons_[syst]->front()  , "kin");

    if(loosePhotons_[syst]->size() > 0){
      SYSplots_photon(  syst_name.c_str(), base_w, loosePhotons_[syst]->front(), "veryLoose");

      if(goodPhotons_[syst]->size() > 0)
	SYSplots_photon(syst_name.c_str(), base_w, goodPhotons_[syst]->front() , "loose");
      else{
	const Photon& failPh = loosePhotons_[syst]->front();
	SYSplots_photon(syst_name.c_str(), base_w          , failPh, "fail"        );
	double f_VLtoL = getPhotonFR_VLtoL(failPh);
	double w_VLtoL = f_VLtoL / (1 - f_VLtoL);
	SYSplots_photon(syst_name.c_str(), base_w * w_VLtoL, failPh, "failReweight");
      }
    }

    auto bestMVA_syst = std::max_element(kinPhotons_[syst]->begin(), kinPhotons_[syst]->end(),
					[](const Photon& a, const Photon& b){ return a.MVAvalue() < b.MVAvalue(); }
					);
    if(bestMVA_syst != kinPhotons_[syst]->end())
      SYSplots_phMVA(syst_name.c_str(), base_w, *bestMVA_syst);
  }

  bool isMC = theSampleInfo.isMC();
  // puWeightUnc
  SYSplots("puWeight_Up"  , base_w * ( isMC ? theSampleInfo.puWeightUncUp() / theSampleInfo.puWeight() : 1.), ph, bestMVAPh_);
  SYSplots("puWeight_Down", base_w * ( isMC ? theSampleInfo.puWeightUncDn() / theSampleInfo.puWeight() : 1.), ph, bestMVAPh_);
  
  // L1PrefiringWeight
  SYSplots("L1Prefiring_Up"  , base_w * ( isMC ? theSampleInfo.L1PrefiringWeightUp() / theSampleInfo.L1PrefiringWeight() : 1.), ph, bestMVAPh_);
  SYSplots("L1Prefiring_Down", base_w * ( isMC ? theSampleInfo.L1PrefiringWeightDn() / theSampleInfo.L1PrefiringWeight() : 1.), ph, bestMVAPh_);
  
  // QCD scale
  // envelope: consider the six variations: {Do, Central, Up} x {Dn, Central, Up} - (central, central) - (Dn, Dn) - (Up, Up) and use the max and min
  float QCDscale_Up(1.), QCDscale_Dn(1.);
  if(isMC){
    std::vector<float> envelope {
      theSampleInfo.QCDscale_muR0p5F1(),
      theSampleInfo.QCDscale_muR0p5F2(),
      theSampleInfo.QCDscale_muR1F0p5(),
      theSampleInfo.QCDscale_muR1F2(),
      theSampleInfo.QCDscale_muR2F0p5(),
      theSampleInfo.QCDscale_muR2F1()
    };
    QCDscale_Up = *max_element(envelope.begin(), envelope.end());
    QCDscale_Dn = *min_element(envelope.begin(), envelope.end());
  }
  SYSplots("QCDscale_Up"  , base_w * QCDscale_Up, ph, bestMVAPh_);
  SYSplots("QCDscale_Down", base_w * QCDscale_Dn, ph, bestMVAPh_);
  
  // PDF var
  SYSplots("PDFVar_Up"  , base_w * ( isMC ? theSampleInfo.PDFVar_Up()   : 1.), ph, bestMVAPh_);
  SYSplots("PDFVar_Down", base_w * ( isMC ? theSampleInfo.PDFVar_Down() : 1.), ph, bestMVAPh_);
  
  // alphas MZ
  SYSplots("alphas_Up"  , base_w * ( isMC ? theSampleInfo.alphas_MZ_Up()   : 1.), ph, bestMVAPh_);
  SYSplots("alphas_Down", base_w * ( isMC ? theSampleInfo.alphas_MZ_Down() : 1.), ph, bestMVAPh_);
  

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
  SYSplots("eleEffSF_Up"  , base_w * (1 + eleEff_w), ph, bestMVAPh_);
  SYSplots("eleEffSF_Down", base_w * (1 - eleEff_w), ph, bestMVAPh_);
  SYSplots("muoEffSF_Up"  , base_w * (1 + muoEff_w), ph, bestMVAPh_);
  SYSplots("muoEffSF_Down", base_w * (1 - muoEff_w), ph, bestMVAPh_);
  
  // lepton fake rate SF
  SYSplots("eleFakeRateSF_Up"  , base_w * (1 + eleFake_w), ph, bestMVAPh_);
  SYSplots("eleFakeRateSF_Down", base_w * (1 - eleFake_w), ph, bestMVAPh_);
  SYSplots("muoFakeRateSF_Up"  , base_w * (1 + muoFake_w), ph, bestMVAPh_);
  SYSplots("muoFakeRateSF_Down", base_w * (1 - muoFake_w), ph, bestMVAPh_);

  // Photon cut-based ID efficiency
  SYSplots_inclusive("phEffSF_Up"  , base_w);
  SYSplots_inclusive("phEffSF_Down", base_w);

  if(ph){
    double phEff_dw = 0.;
    if(ph->cutBasedID(Photon::IdWp::Loose) && getPhotonEffSF(*ph) != 0)
      phEff_dw = getPhotonEffSFUnc(*ph)/getPhotonEffSF(*ph);

    SYSplots_phCut("phEffSF_Up"  , base_w * (1 + phEff_dw), *ph);
    SYSplots_phCut("phEffSF_Down", base_w * (1 - phEff_dw), *ph);
  }
  if(bestMVAPh_){  // Should be true if the previous ph != nullptr is true
    SYSplots_phMVA("phEffSF_Up"  , base_w, *bestMVAPh_);
    SYSplots_phMVA("phEffSF_Down", base_w, *bestMVAPh_);
  }

  // Photon MVA ID efficiency
  // - inclusive
  SYSplots_inclusive("phEffMVASF_Up"  , base_w);
  SYSplots_inclusive("phEffMVASF_Down", base_w);
  // - cut based IDs
  if(ph){
    SYSplots_phCut("phEffMVASF_Up"  , base_w, *ph);
    SYSplots_phCut("phEffMVASF_Down", base_w, *ph);
  }
  // - MVA based ID --> fill plots "manually"
  if(bestMVAPh_){
    bool pass80 = bestMVAPh_->passMVA(Photon::MVAwp::wp80);
    bool pass90 = bestMVAPh_->passMVA(Photon::MVAwp::wp90);
    double effSF_wp80 = getPhotonEffSF_MVA(*bestMVAPh_, Photon::MVAwp::wp80);
    double effSF_wp90 = getPhotonEffSF_MVA(*bestMVAPh_, Photon::MVAwp::wp90);
    std::string MVAcut_s = "none";
    double effSF(1.), effdw(0.);

    if(pass90 && effSF_wp90 != 0){
      effSF = effSF_wp90;
      MVAcut_s = "wp90";
      getPhotonEffSFUnc_MVA(*bestMVAPh_, Photon::MVAwp::wp90);
      double effdw90 = effdw = getPhotonEffSFUnc(*ph)/effSF_wp90;
      SYSplots_photon("phEffMVASF_Up"  , base_w * effSF_wp90 * (1 + effdw90), *bestMVAPh_, "wp90");
      SYSplots_photon("phEffMVASF_Down", base_w * effSF_wp90 * (1 - effdw90), *bestMVAPh_, "wp90");

      if(pass80 && effSF_wp80 != 0){
	effSF = effSF_wp80;
	MVAcut_s = "wp80";
	getPhotonEffSFUnc_MVA(*bestMVAPh_, Photon::MVAwp::wp80);
	double effdw80 = effdw = getPhotonEffSFUnc(*ph)/effSF_wp80;
	SYSplots_photon("phEffMVASF_Up"  , base_w * effSF_wp80 * (1 + effdw80), *bestMVAPh_, "wp80");
	SYSplots_photon("phEffMVASF_Down", base_w * effSF_wp80 * (1 - effdw80), *bestMVAPh_, "wp80");
      }
      else{
	SYSplots_photon("phEffMVASF_Up"  , base_w * effSF_wp90 * (1 + effdw90), *bestMVAPh_, "90not80");
	SYSplots_photon("phEffMVASF_Down", base_w * effSF_wp90 * (1 - effdw90), *bestMVAPh_, "90not80");
      }
    }

    theHistograms->fill("SYS_MVAcut_phEffMVASF_Up"  , "MVAcut phEffMVASF_Up"  , {"none","wp90","wp80"}, MVAcut_s.c_str(), base_w * effSF * (1 + effdw));
    theHistograms->fill("SYS_MVAcut_phEffMVASF_Down", "MVAcut phEffMVASF_Down", {"none","wp90","wp80"}, MVAcut_s.c_str(), base_w * effSF * (1 - effdw));
    if(theSampleInfo.isMC()){
      const char* phGenStatus;
      phGenStatus = sigdefHelper.pass_photon() ? "prompt" : "nonpro" ;
      theHistograms->fill(Form("SYS_MVAcut-%s_phEffMVASF_Up"  , phGenStatus), Form("MVAcut %s phEffMVASF_Up"  , phGenStatus), {"none","wp90","wp80"}, MVAcut_s.c_str(), base_w * effSF * (1 + effdw));
      theHistograms->fill(Form("SYS_MVAcut-%s_phEffMVASF_Down", phGenStatus), Form("MVAcut %s phEffMVASF_Down", phGenStatus), {"none","wp90","wp80"}, MVAcut_s.c_str(), base_w * effSF * (1 - effdw));
    }
  }

  // Photon FR uncertaintiy  WARN: for this to have a meaning, the photon FR SF should be applied
  SYSplots_inclusive("phFakeRate_Up"  , base_w);
  SYSplots_inclusive("phFakeRate_Down", base_w);
  if(ph){
    double w_up = 1.;
    double w_dn = 1.;
    if(ph->cutBasedID(Photon::IdWp::VeryLoose) && ! ph->cutBasedID(Photon::IdWp::Loose)){
      double func = getPhotonFRUnc_VLtoL(*ph);
      double f_ce = getPhotonFR_VLtoL(*ph);
      double f_up = f_ce + func;
      double f_dn = std::max(f_ce - func, 0.);
      double sf_ce = f_ce/(1 - f_ce);
      double sf_up = f_up/(1 - f_up);
      double sf_dn = f_dn/(1 - f_dn);
      w_up = sf_up/sf_ce;
      w_dn = sf_dn/sf_ce;
    }
    SYSplots_phCut("phFakeRate_Up"  , base_w * w_up, *ph);
    SYSplots_phCut("phFakeRate_Down", base_w * w_dn, *ph);
  }
  if(bestMVAPh_){
    SYSplots_phMVA("phFakeRate_Up"  , base_w, *bestMVAPh_);
    SYSplots_phMVA("phFakeRate_Down", base_w, *bestMVAPh_);
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
  // vector<Particle> genPhotonsHard;
  // genPhotonsHard.reserve(genPhotonsPrompt_->size());
  // copy_if(genPhotonsPrompt_->begin(), genPhotonsPrompt_->end(), std::back_inserter(genPhotonsHard),
  //         [](const Particle& p){
  //           std::bitset<15> flags = p.genStatusFlags();
  //           return (flags.test(isHardProcess) || flags.test(fromHardProcess)); }
  //         );
  
  std::map<const char*, const vector<Photon>*> mapWPtoPhotons;
  mapWPtoPhotons["kin"]       = kinPhotons_["central"].get();
  mapWPtoPhotons["veryLoose"] = loosePhotons_["central"].get();
  mapWPtoPhotons["loose"]     = goodPhotons_["central"].get();
  
  // The rec-gen matching is done separately for kin, veryLoose and loose, since we don't want e.g. a kin reco photon 
  // that fails the LooseID to be selected just because it comes before in the photon vector
  for(auto & [wp, pPhVect] : mapWPtoPhotons){
    const vector<Photon>& phVect = *pPhVect;
    if(phVect.size() == 0)
      continue;
    
    auto best = phVect.cbegin();  // vector<Photon>::const_iterator
    const char* matched = "nomatch";
    bool isMatched = false;

    // First try matching with prompt gen photons
    for(auto gen : *genPhotonsPrompt_){
      auto closestRec = closestDeltaR(gen, phVect);
      if(closestRec != phVect.cend() && deltaR(gen, *closestRec) < 0.1){
	best = closestRec;
	matched = "promptm";
	isMatched = true;
	break;
      }
    }
    if(!isMatched){       // Try again with all gen photons
      for(auto gen : *genPhotons_){
	auto closestRec = closestDeltaR(gen, phVect);
	if(closestRec != phVect.cend() && deltaR(gen, *closestRec) < 0.1){
	  best = closestRec;
	  matched = "matched";
	  break;
	}
      }
    }
    
    theHistograms->fill(Form("PhGenStudy_status_%s", wp), Form("%s;;Events", wp), {"promptm", "matched", "nomatch"}, matched, theWeight);

    const char* hname = Form("PhGenStudy_%s_%s_%s", "%s", wp, matched);
    
    // Find the closest lep and get the DR
    auto closestLep = closestDeltaR(*best, *leptons_);
    float dR = closestLep != leptons_->cend() ? deltaR(*best, *closestLep) : 10;
    theHistograms->fill(Form(hname, "DRLep"), Form("%s;min #DeltaR(#gamma, l);Events", wp), 50, 0., 1., dR   , theWeight);
    
    // Closest jet
    auto closestJet = closestDeltaR(*best, *jets);
    float dRJet = closestJet != jets->cend() ? deltaR(*best, *closestJet) : 10;
    theHistograms->fill(Form(hname, "DRJet"), Form("%s;#DeltaR(#gamma, j);Events"    , wp), 50, 0., 1 , dRJet, theWeight);
    if(closestJet != jets->cend()){
      float dPJet = (closestJet->p4()     - best->p4()    ).P() / closestJet->p4().P();
      float fPJet = best->p4().P() / closestJet->p4().P();
      theHistograms->fill(Form(hname, "vsJet"), Form("%s;P^{#gamma}/P^{j};#DeltaR(#gamma,j);Events"    , wp), 10,0., 1., 10,0.,1., fPJet, dRJet, theWeight);
      if(dRJet < 0.2){
	theHistograms->fill(Form(hname, "DPJet"), Form("%s;(#vec{P^{#gamma}}-#vec{P^{j}})/P^{j};Events", wp), 40,  0., 1 , dPJet, theWeight);
	theHistograms->fill(Form(hname, "FPJet"), Form("%s;P^{#gamma}/P^{j};Events"                    , wp), 40,  0., 1., fPJet, theWeight);
      }
    }
    
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


void VVGammaAnalyzer::ZllVsZllGstudy(const std::vector<phys::Photon>& phvect, const char* label){
  double Zll_mass(0.), ZllG_mass(0.);
  std::tie(Zll_mass, ZllG_mass) = getZllAndZllgMasses(phvect);
  if(Zll_mass <= 0)
    return;

  bool improves = ( fabs(ZllG_mass - phys::ZMASS) < fabs(Zll_mass - phys::ZMASS) ); // The addition of the photon momentum improved the mass
  const char* improv_str = improves ? "improves" : "noimprov";
  theHistograms->fill(Form("Zll_mass_%s_%s" , label, improv_str), ";m_{ll} [GeV/c^{2}]"      , 30,60,120, Zll_mass , theWeight);
  theHistograms->fill(Form("ZllG_mass_%s_%s", label, improv_str), ";m_{ll#gamma} [GeV/c^{2}]", 30,60,210, ZllG_mass, theWeight);
  if(theSampleInfo.isMC()){
    bool signaldef = sigdefHelper.pass_photon();  // Test if the event passes the GEN signal definition
    const char* sigdef_str = signaldef ? "prompt" : "nonpro";
    theHistograms->fill(Form("Zll_mass_%s_%s_%s" , label, improv_str, sigdef_str), ";m_{ll} [GeV/c^{2}]"      , 30,60,120, Zll_mass , theWeight);
    theHistograms->fill(Form("ZllG_mass_%s_%s_%s", label, improv_str, sigdef_str), ";m_{ll#gamma} [GeV/c^{2}]", 30,60,210, ZllG_mass, theWeight);
  }

  theHistograms->fill(Form("ZllG_mass_vs_Zll_mass_%s"  , label), ";m_{ll#gamma};m_{ll}", 60,60,150., 30,60.,120., ZllG_mass, Zll_mass, theWeight);
  theHistograms->fill(Form("Zll_mass_%s"               , label), ";m_{ll}"             , 30,60.,120. , Zll_mass, theWeight);
  theHistograms->fill(Form("ZllG_mass_%s"              , label), ";m_{ll#gamma}"       , 40,60.,160. , ZllG_mass, theWeight);
  theHistograms->fill(Form("Zll_mass_plus_ZllG_mass_%s", label), ";m_{ll}+m_{ll#gamma}", 40,120.,320., Zll_mass+ZllG_mass, theWeight);
}


void VVGammaAnalyzer::fillCutsNm1(const std::string& name, const std::string& title, const std::vector<std::pair<std::string, bool>>& cuts, const double& weight = 1){
  std::vector<std::string> cut_names;
  cut_names.reserve(cuts.size()+2);
  // cut_names.push_back("all debug");
  cut_names.push_back("Nm1");
  for(auto& name_val : cuts) cut_names.push_back(name_val.first);

  unsigned int ncuts = std::count_if(cuts.begin(), cuts.end(), [](auto& p){ return p.second; });
  // theHistograms->fill(name, title, cut_names, "all debug", weight);

  if(ncuts < cuts.size() - 1)
    return;

  theHistograms->fill(name, title, cut_names, "Nm1", weight);
  bool passAll = (ncuts == cuts.size());
  for(auto& name_val : cuts){
    if(passAll || !name_val.second)
      theHistograms->fill(name, title, cut_names, name_val.first, weight);
  }

  // if(passAll)
  //   theHistograms->fill(name, title, cut_names, "N", weight);
}


void VVGammaAnalyzer::fillCutFlow(const std::string& name, const std::string& title, const std::vector<std::pair<std::string, bool>>& cuts, const double& weight){
  std::vector<std::string> cut_names;
  cut_names.reserve(cuts.size()+1);
  cut_names.push_back("All");
  for(auto& name_val : cuts) cut_names.push_back(name_val.first);

  theHistograms->fill(name, title, cut_names, "All", weight);
  for(auto& name_val : cuts){
    if(name_val.second)
      theHistograms->fill(name, title, cut_names, name_val.first, weight);
    else
      break;
  }
}


// Utilities
double VVGammaAnalyzer::getPhotonFR_VLtoL   (const phys::Photon& ph) const{
  const TH2F* h = hPhotonFR_VLtoL_data_.get();  // TODO: move it to dataZG
  return h->GetBinContent(h->FindFixBin(
					ph.pt() < 120 ? ph.pt() : 119.9,
					abs(ph.eta())
					));
}

double VVGammaAnalyzer::getPhotonFRUnc_VLtoL(const phys::Photon& ph) const{
  const TH2F* h = hPhotonFR_VLtoL_data_.get();  // TODO: move it to dataZG
  return h->GetBinError(h->FindFixBin(
				      ph.pt() < 120 ? ph.pt() : 119.9,
				      abs(ph.eta())
				      ));
}


double VVGammaAnalyzer::getPhotonFR_VLtoL_data(const phys::Photon& ph) const{
  const TH2F* h = hPhotonFR_VLtoL_data_.get();  // TODO: move it to dataZG
  return h->GetBinContent(h->FindFixBin(
					ph.pt() < 120 ? ph.pt() : 119.9,
					abs(ph.eta())
					));
}


double VVGammaAnalyzer::getPhotonFR_VLtoL_dataZG(const phys::Photon& ph) const{
  const TH2F* h = hPhotonFR_VLtoL_dataZG_.get();
  return h->GetBinContent(h->FindFixBin(
					ph.pt() < 120 ? ph.pt() : 119.9,
					abs(ph.eta())
					));
}

double VVGammaAnalyzer::getPhotonFR_KtoVLexcl(const phys::Photon& ph) const{
  const TH2F* h = hPhotonFR_KtoVLexcl_.get();
  return h->GetBinContent(h->FindFixBin(
					ph.pt() < 120 ? ph.pt() : 119.9,
					abs(ph.eta())
					));
}

double VVGammaAnalyzer::getPhotonFRSF_VLtoL(const phys::Photon& ph) const{
  const TH2F* h = hPhotonFRSF_VLtoL_.get();
  return h->GetBinContent(h->FindFixBin(
					ph.pt() < 120 ? ph.pt() : 119.9,
					abs(ph.eta())
					));
}

double VVGammaAnalyzer::getPhotonEffSF_MVA(const phys::Photon& ph, Photon::MVAwp wp) const{
  const TH2F* hSF = mapPhotonMVASF_.at(wp).get();
  float maxPt = mapPhotonMVASF_maxPt_.at(wp);
  Int_t bin = hSF->FindFixBin(
			      ph.eta(),
			      ph.pt() < maxPt ? ph.pt() : maxPt-0.1
			      );
  return hSF->GetBinContent(bin);
}

double VVGammaAnalyzer::getPhotonEffSFUnc_MVA(const phys::Photon& ph, Photon::MVAwp wp) const{
  const TH2F* hSF = mapPhotonMVASF_.at(wp).get();
  float maxPt = mapPhotonMVASF_maxPt_.at(wp);
  Int_t bin = hSF->FindFixBin(
			      ph.eta(),
			      ph.pt() < maxPt ? ph.pt() : maxPt-0.1
			      );
  return hSF->GetBinError(bin);
}

bool VVGammaAnalyzer::passVeryLoose(const Photon& ph){
  Photon::IdWp wp = Photon::IdWp::Loose;
  return ph.cutBasedID(wp, Photon::IDcut::HoverE) && ph.cutBasedID(wp, Photon::IDcut::phIso) && ph.cutBasedID(wp, Photon::IDcut::neIso);
}

bool inPhotonEtaAcceptance(double eta){
  double aeta = fabs(eta);
  if(aeta > CUT_G_AETA_MAX || (aeta > CUT_G_AETA_GAP_MIN && aeta < CUT_G_AETA_GAP_MAX))
    return false;
  return true;
}


const vector<double> VVGammaAnalyzer::pt_bins(
  {20., 35., 50., 100., 200., 500.}
					      );

const vector<double> VVGammaAnalyzer::pt_bins_LFR(
  {5., 7., 10., 20., 30., 40., 50., 80}
					      );

const vector<double> VVGammaAnalyzer::eta_bins(
  {-CUT_G_AETA_MAX, -2., -CUT_G_AETA_GAP_MAX, -CUT_G_AETA_GAP_MIN, -0.8, 0., 0.8, CUT_G_AETA_GAP_MIN, CUT_G_AETA_GAP_MAX, 2., CUT_G_AETA_MAX}
  //{-2.4, -2.1, -1.8, -1.566, -1.4442, -1.2, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.2, 1.4442, 1.566, 1.8, 2.1, 2.4}
					       );

const vector<double> VVGammaAnalyzer::aeta_bins(
  {0., 0.8, CUT_G_AETA_GAP_MIN, CUT_G_AETA_GAP_MAX, 2., CUT_G_AETA_MAX}
  //{0., 0.3, 0.6, 0.9, 1.2, 1.4442, 1.566, 1.8, 2.1, 2.4}
						);

const vector<double> VVGammaAnalyzer::ph_aeta_bins(
						   {0., 0.8, CUT_G_AETA_GAP_MIN, CUT_G_AETA_GAP_MAX, 2, CUT_G_AETA_MAX}
						   // {0., 0.435, 0.783, 1.13, 1.4442, 1.566, 1.8, 2.1, 2.5}
						   );

const vector<double> VVGammaAnalyzer::ph_aetaExtended_bins {
  0.0, 0.072, 0.144, 0.216, 0.288, 0.36, 0.432, 0.504, 0.576, 0.648, 0.72, 0.792, 0.864, 0.936, 1.008, 1.08, 1.152, 1.224, 1.296, 1.368, 1.4442, // 20
    1.566, 1.642, 1.72, 1.798, 1.876, 1.954, 2.032, 2.11, 2.188, 2.266, 2.344, 2.422, 2.5  // 13
							     };

const vector<double> VVGammaAnalyzer::ph_pt_bins(
						 {20., 25., 35., 50., 80., 120}
						 );

const vector<double> VVGammaAnalyzer::ph_ptExtended_bins {
  20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200
    };


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
