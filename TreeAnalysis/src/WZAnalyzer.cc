#include "VVXAnalysis/TreeAnalysis/interface/WZAnalyzer.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include <TSpline.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 

using namespace boost::assign;
using namespace colour;
using namespace physmath;
using namespace phys;
using namespace std;

// !!!! to analyze data, comment line 489 !!!! still valid in 2021?

void WZAnalyzer::begin() {
  
  cout << "\n--------------------------------------------------------------------------" << endl;
  cout << "\n Starting WZAnalysis on sample in " << fileName << endl;

  // events counters
  eventGen = 0;
  eventReco = 0;
  eventGenReco = 0;
  eventSample = 0;

  genAfterCut = 0; 
  recoZlempty = 0;
  recoJetless2 = 0;
  recoAfterCut = 0;

  gen3e = 0;
  gen3m = 0;
  gen2e1m = 0;
  gen2m1e = 0;
  reco3e = 0;
  reco3m = 0;
  reco2e1m = 0;
  reco2m1e = 0;
  
  // free counters
  counter1 = 0; // gen firstly Z, secondly W -> Z's daughters of different flavour (em or me)
  counter2 = 0; // gen firstly W, secondly Z -> Z's daughters of different flavour (em or me)
  counter3 = 0; // firstly Z, secondly W -> choosed Z's daughter of different flavour (even if not final Z)
  counter4 = 0; // firstly W, secondly Z -> choosed Z's daughter of different flavour (even if not final Z)
  counter5 = 0; // gen Z and reco W have different IDs
  counter6 = 0;
  
  // time begins
  begintime = ((float)clock())/CLOCKS_PER_SEC;
}


Int_t WZAnalyzer::cut() {
  return 1;
}


void WZAnalyzer::GenAnalysis(DiBosonParticle &WZ, Particle &Jet0, Particle &Jet1){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~ Begin of gen Analysis ~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  counter1++;
  string eventkind = "WZ";
  int cut = 0;
  theHistograms.fill("GenCuts", "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);

  tuple<bool, BosonParticle, BosonParticle> WZtuple;
  

  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Z & W ~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  // ----- filter on leptons, total mass >100GeV, smart cut on
  zz::SignalTopology topology = zz::getSignalTopology(*genParticles, *genJets, *genJetsAK8);

  //VVjjHelper::test(2);
  
  BosonParticle Zet = get<1>(topology);
  BosonParticle Weh = get<5>(topology);
  bitset<16> top = (bitset<16>)get<0>(topology);

  theHistograms.fill("GenTopology", "Topology di WZ", 10, -0.5, 9.5, top.test(0));
  theHistograms.fill("GenTopology", "Topology di WZ", 10, -0.5, 9.5, top.test(1)+2);
  theHistograms.fill("GenTopology", "Topology di WZ", 10, -0.5, 9.5, top.test(5)+4);
  theHistograms.fill("GenTopology", "Topology di WZ", 10, -0.5, 9.5, top.test(6)+6);
  theHistograms.fill("GenTopology", "Topology di WZ", 10, -0.5, 9.5, top.test(7)+8);
  
  if(top.test(0) == 1 || top.test(1) == 0 || top.test(5) == 0 || top.test(6) == 0 || top.test(7) == 0){
    return;
  }
  
  cut++;
  theHistograms.fill("GenCuts", "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);

  //VVjjHelper::test(3);
  
  vector<pairBosonParticle> Zls; Zls.push_back(pairBosonParticle(Zet, Weh.daughter(0)));

  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Jets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ----- filter on jets' number and eta  
  helper_->FindLeadingJets(genJets, Jet0, Jet1, genParticles);

  if(Jet0.p4().Mt() == 0){
    wrongjet++;
    return;
  }

  cut++;
  theHistograms.fill("GenCuts", "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);

  //VVjjHelper::test(1);
  
  // leading jets
  TLorentzVector JJp4 = Jet0.p4() + Jet1.p4();
  double JJdeltaEta = Jet0.eta() - Jet1.eta();
  
  eventGen++;
  theHistograms.fill("WZ_Events", "Weighted counters", 10, -0.5, 10.5, 2, theWeight);

  
  // ----------------------- W & Z histograms ----------------------  
  helper_->PlotBoson(Zet, "GenZ", theWeight, "BC");
  helper_->PlotBoson(Weh, "GenW", theWeight, "BC");
  helper_->PlotDiBoson(DiBosonParticle(Weh, Zet), "GenWZ", theWeight, "BC");
  helper_->PlotJets(Jet0, Jet1, "Gen", theWeight, "BC");
  
  // Zl
  int ZlsID = abs(Zls[0].first.daughter(0).id()) + abs(Zls[0].first.daughter(1).id()) + abs(Zls[0].second.id());

  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Cut ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  float zeppenllJ0 = Jet0.eta() - (Zet.daughter(0).eta() + Zet.daughter(1).eta())/2;

  
  // --------------- Cuts implemented in reco analysis -------------
  bool filtromassaZ = abs(ZMASS - Zet.mass()) > 15.;
  bool filtroMET = Weh.daughter(1).pt() < 30.;
  bool filtrotrmassaWmin = Weh.p4().Mt() < 30.;
  bool filtrotrmassaWmax = Weh.p4().Mt() > 500.;
  bool filtroJJdeltaEta = abs(JJdeltaEta) < 2.4;
  bool filtroJJpt = JJp4.M() < 285.5;
  bool filtrozeppenfeldllJ0 = abs(zeppenllJ0) < 0.6;
  bool filtroJ1pt = Jet1.pt() < 50.;
  bool filtroWleppt = Weh.daughter(0).pt() < 30.;

  /*
  if(filtromassaZ || filtroMET || filtrotrmassaWmin || filtrotrmassaWmax || filtroJJdeltaEta || filtroJJpt || filtrozeppenfeldllJ0 || filtroJ1pt || filtroWleppt){
    return;
  }
  */

  if(filtroMET) return;
  cut++;
  theHistograms.fill("GenCuts", "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);

  if(filtrotrmassaWmin || filtrotrmassaWmax) return;
  cut++;
  theHistograms.fill("GenCuts", "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);
  
  if(filtroJJdeltaEta) return;
  cut++;
  theHistograms.fill("GenCuts", "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);
  
  if(filtroJJpt) return;
  cut++;
  theHistograms.fill("GenCuts", "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);
  
  if(filtrozeppenfeldllJ0) return;
  cut++;
  theHistograms.fill("GenCuts", "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);
  
  if(filtroJ1pt) return;
  cut++;
  theHistograms.fill("GenCuts", "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);
  
  if(filtroWleppt) return;
  cut++;
  theHistograms.fill("GenCuts", "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);
  
  if(filtromassaZ) return;
  cut++;
  theHistograms.fill("GenCuts", "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);


  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ All ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TLorentzVector WZjjp4 = Weh.p4() + Zet.p4() + Jet0.p4() + Jet1.p4();
  
  // W&Z diboson
  genAfterCut++;
  theHistograms.fill("WZ_Events", "Weighted counters", 10, -0.5, 10.5, 3, theWeight);
  WZ = DiBosonParticle(Weh, Zet);
  
  // ----------------------- W & Z histograms ----------------------
  
  helper_->PlotBoson(Zet, "GenZ", theWeight, "AC");
  helper_->PlotBoson(Weh, "GenW", theWeight, "AC");
  helper_->PlotDiBoson(WZ, "GenWZ", theWeight, "AC");
  helper_->PlotJets(Jet0, Jet1, "Gen", theWeight, "AC");

  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~ End of gen Analysis ~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}


void WZAnalyzer::RecoAnalysis(DiBosonLepton &WZ, Particle &Jet0, Particle &Jet1){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~ Begin of reco Analysis ~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int cut = 0;
  int cut2 = 0;
  theHistograms.fill("RecoCuts", "Reco events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("RecoCuts_wei", "Reco events after cuts", 14, -0.5, 13.5, cut, theWeight);

  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Filters ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ZLCompositeCandidates recoZls;
  ZLCompositeCandidates tempZls;
  ZLCompositeCandidates temp2Zls;
  
  // ----- filter on ZLCand's size > 0
  theHistograms.fill("recoZl_size_ZLCand", 12, -0.5, 11.5, ZLCand->size());

  if(ZLCand->size() == 0){
    recoZlempty++;
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCuts_wei", "Events after cuts", 13, -0.5, 12.5, cut, theWeight);

  // ----- 3 leptons
  if(electrons->size() + muons->size() != 3)
    return;
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCuts_wei", "Events after cuts", 13, -0.5, 12.5, cut, theWeight);
    
  // ----- filter on 3rd lepton full selection
  foreach(const ZLCompositeCandidate Zl, *ZLCand){
    if(Zl.second.passFullSel()){
      tempZls.push_back(Zl);
    }
  }

  if(tempZls.size() == 0){
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCuts_wei", "Events after cuts", 13, -0.5, 12.5, cut, theWeight);
  
  // ----- filter on 3rd lepton pt > 30
    foreach(const ZLCompositeCandidate Zl, tempZls){
    if(Zl.second.pt() > 30.){
      temp2Zls.push_back(Zl);
    }
  }

  if(temp2Zls.size() == 0){
    recoZlempty++;
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCuts_wei", "Events after cuts", 13, -0.5, 12.5, cut, theWeight);
  
  // ----- filter on MET's pt
  if(met->pt() < 30.){
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCuts_wei", "Events after cuts", 13, -0.5, 12.5, cut, theWeight);

  // ----- filter on Z's mass 60 < m < 120
  foreach(const ZLCompositeCandidate Zl, temp2Zls){
    if(Zl.first.mass() > 60. && Zl.first.mass() < 120.){
      recoZls.push_back(Zl);
    }
  }

  if(recoZls.size() == 0){
    recoZlempty++;
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCuts_wei", "Events after cuts", 13, -0.5, 12.5, cut, theWeight);
  
  // ----- filter on jet's number >= 2
  if(jets->size() < 2){
    recoJetless2++;
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCuts_wei", "Events after cuts", 13, -0.5, 12.5, cut, theWeight);

  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Z & W ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BosonLepton recoW;
  BosonLepton recoWtemp;
  BosonLepton recoZ;
  Lepton MET = Lepton(met->p4());
  vector<DiBosonLepton> recoWZs;
  
  // ------------------------ W construction -----------------------
  
  if(recoZls.size() == 1){    
    recoZ = recoZls[0].first;
    recoW = BosonLepton(recoZls[0].second, MET, copysign(24, recoZls[0].second.charge()));
  
    // ----- filter on W's trmass 30 < trm < 500
    if(recoW.p4().Mt() >= 30. && recoW.p4().Mt() <= 500.){
      recoWZs.push_back(DiBosonLepton(recoW, recoZ));
    }    
    else{ 
      return;
    }    
  }
  
  else if(recoZls.size() > 1){    
    // Z is the one with closer mass to Z
    sort(recoZls.begin(), recoZls.end(), pairMassComparator(0, ZMASS));
    
    // W is the one with the lepton with higher pt
    bool choosingW = kTRUE;
    
    for(int i = 0; choosingW; i++){
      recoWtemp = BosonLepton(recoZls[i].second, MET, copysign(24, recoZls[i].second.charge()));

      // ----- filter on W's trmass 30 < trm < 500
      if(recoWtemp.p4().Mt() >= 30. && recoWtemp.p4().Mt() <= 500.){
	recoWZs.push_back(DiBosonLepton(recoWtemp, recoZls[i].first));
      }
      
      if(i == (int)recoZls.size() - 1 || recoZls[i].first.mass() - recoZls[i + 1].first.mass() > 0.1){
	choosingW = kFALSE;
      }
    }

    if(recoWZs.size() == 0){
      return;
    }
    
    sort(recoWZs.begin(), recoWZs.end(), PtComparator());
    
    recoW = recoWZs[0].first();
    recoZ = recoWZs[0].second();
  }
  
  cut++;
  theHistograms.fill("RecoCuts",      "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCuts_pres", "Events after cuts", 13, -0.5, 12.5, cut2);
  theHistograms.fill("RecoCuts_wei",  "Events after cuts", 13, -0.5, 12.5, cut, theWeight);
  
  eventReco++;
  theHistograms.fill("WZ_Events", "Weighted counters", 10, -0.5, 10.5, 4, theWeight);
  

  
  // ~~~~~~ Refine event weight with efficiency scale factors ~~~~~~  
  //theWeight *= WZ.efficiencySF();


  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ------------------------ W & Z variables ----------------------
  // Z
  float recoZpzlp; //lp = positive lepton
  float recoZpzlm; //lm = negative lepton
  if(recoZ.daughter(0).id() > 0){
    recoZpzlp = recoZ.daughter(0).p4().Pz();
    recoZpzlm = recoZ.daughter(1).p4().Pz();
  }
  else{
    recoZpzlm = recoZ.daughter(0).p4().Pz();
    recoZpzlp = recoZ.daughter(1).p4().Pz();
  }

  
  // Zl
  int recoZlID = abs(recoZ.daughter(0).id()) + abs(recoZ.daughter(1).id()) + abs(recoW.daughter(0).id());
  TLorentzVector recoZlp4 = recoZ.daughter(0).p4() + recoZ.daughter(1).p4() + recoW.daughter(0).p4();
  
  if(recoZlID == 33)      reco3e++;
  else if(recoZlID == 35) reco2e1m++;
  else if(recoZlID == 37) reco2m1e++;
  else if(recoZlID == 39) reco3m++;

  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Jets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  helper_->FindLeadingJets(jets, Jet0, Jet1);
  
  // leading jets
  TLorentzVector recoJJptot = Jet0.p4() + Jet1.p4();
  double recoJJdeltaEta = Jet0.eta() - Jet1.eta();
  double recoJJdeltaPhi = physmath::deltaPhi(Jet0.phi(), Jet1.phi());
  double recoJJdeltaR = abs(physmath::deltaR(Jet0, Jet1));


  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ All ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TLorentzVector recoPtot = recoW.p4() + recoZ.p4() + Jet0.p4() + Jet1.p4();

  double WZjptDeltaPhi = physmath::deltaPhi(recoWZs[0].phi(), Jet0.phi());
  
  float zeppenfeldJJZ = recoZ.eta() - (Jet0.eta() + Jet1.eta())/2;
  float zeppenfeldJJl = recoW.daughter(0).eta() - (Jet0.eta() + Jet1.eta())/2;
  float zeppenfeldllJ0 = Jet0.eta() - (recoZ.daughter(0).eta() + recoZ.daughter(1).eta())/2;
  float zeppenfeldllJ1 = Jet1.eta() - (recoZ.daughter(0).eta() + recoZ.daughter(1).eta())/2;
  float zeppenfeldlll = recoW.daughter(0).eta() - (recoZ.daughter(0).eta() + recoZ.daughter(1).eta())/2;
  float zeppenfeldArticle = recoZlp4.Eta() - (Jet0.eta() + Jet1.eta())/2;


  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Cuts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // interesting histos before cuts
  theHistograms.fill("recoAll_trmass_BC",         "Transverse mass recoW,Z,J",          34, 220, 7095   , recoPtot.Mt()         , theWeight);
  theHistograms.fill("recoWZ_trmass_BC",          "recoWZ trmass",                      48,   0, 2400   , recoWZs[0].p4().Mt()  , theWeight);
  theHistograms.fill("recoJJ_mass_BC",            "Jets' mass",                          9,  85, 4585   , recoJJptot.M()        , theWeight);
  theHistograms.fill("recoJJ_mass_BC_2",          "Jets' mass",                         45,  85, 4585   , recoJJptot.M()        , theWeight);
  theHistograms.fill("recoJJ_deltaEta_BC",        "Jets' #Delta#eta",                   50,  -9,    9   , recoJJdeltaEta        , theWeight); 
  theHistograms.fill("recoJJ_deltaEta_abs_BC",    "Jets' #Delta#eta",                   25,   0,    9   , abs(recoJJdeltaEta)   , theWeight);
  theHistograms.fill("recoZ_mass_BC",             "recoZ's mass",                       35,  55,  125   , recoZ.mass()          , theWeight);
  theHistograms.fill("Zeppenfeld_llJ0_BC",        "Zeppenfeld variable for Z and Jet0", 50,  -6,    6   , zeppenfeldllJ0        , theWeight);
  theHistograms.fill("Zeppenfeld_llJ0_abs_BC",    "Zeppenfeld variable for Z and Jet0", 26,   0,    6.24, abs(zeppenfeldllJ0)   , theWeight);
  theHistograms.fill("Zeppenfeld_article_BC",     "Zeppenfeld variable from article",   50,  -6,    6   , zeppenfeldArticle     , theWeight);
  theHistograms.fill("Zeppenfeld_article_abs_BC", "Zeppenfeld variable from article",   26,   0,    6.24, abs(zeppenfeldArticle), theWeight);
  theHistograms.fill("recoAll_Zj0_deltaEta_BC",   "#Delta#eta between recoZ and Jet0",  50,  -6,    6   , recoZ.eta()-Jet0.eta(), theWeight);
  theHistograms.fill("recoJ0_pt_BC",              "Jet0's p_{t}",                       40,   0,  400   , Jet0.pt()             , theWeight);
  theHistograms.fill("recoJ1_pt_BC",              "Jet1's p_{t}",                       40,   0,  400   , Jet1.pt()             , theWeight);
  theHistograms.fill("recoZl_3rd_pt_BC",          "W's lepton p_{t}",                   50,   0,  400   , recoW.daughter(0).pt(), theWeight);
  theHistograms.fill("recoZl_size",   "Reco Zl's size",           15, -0.5,   14.5, recoZls.size()        , theWeight);
  theHistograms.fill("recoWZ_size",     "Reco WZ's size",        15,  -0.5,   14.5, recoWZs.size()      , theWeight); 

  
  ///*
  if(abs(recoJJdeltaEta) < 2.4){ //2.592; filter on JJ's deltaEta 2.106
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCuts_wei", "Events after cuts", 13, -0.5, 12.5, cut, theWeight);
  cut2++;
  theHistograms.fill("RecoCuts_pres", "Events after cuts", 13, -0.5, 12.5, cut2);
  
  theHistograms.fill("recoJJ_mass_AEta",         "Jets' mass",                         45, 85, 4585   , recoJJptot.M()        , theWeight);
  theHistograms.fill("recoZ_mass_AEta",          "recoZ's mass",                       35, 55,  125   , recoZ.mass()          , theWeight);
  theHistograms.fill("Zeppenfeld_llJ0_abs_AEta", "Zeppenfeld variable for Z and Jet0", 26,  0,    6.24, abs(zeppenfeldllJ0)   , theWeight);
  theHistograms.fill("recoJ0_pt_BC",             "Jet0's p_{t}",                       40,  0,  400   , Jet0.pt()             , theWeight);
  theHistograms.fill("recoJ1_pt_BC",             "Jet1's p_{t}",                       40,  0,  400   , Jet1.pt()             , theWeight);
  theHistograms.fill("recoZl_3rd_pt_BC",         "W's lepton p_{t}",                   50,  0,  400   , recoW.daughter(0).pt(), theWeight);
  //*/
  
  ///*
  if(recoJJptot.M() < 285.5){ //397.773; filter on JJ's mass
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCuts_wei", "Events after cuts", 13, -0.5, 12.5, cut, theWeight);
  cut2++;
  theHistograms.fill("RecoCuts_pres", "Events after cuts", 13, -0.5, 12.5, cut2);
  //*/

  
  theHistograms.fill("ZLCand_Zeppenfeld", "ZLCand_Zeppenfeld", 50, -6, 6, zeppenfeldllJ0);
  theHistograms.fill("ZLCand_EtaJ0", "ZLCand_EtaJ0", 50, -9, 9, Jet0.eta());
  theHistograms.fill("ZLCand_EtaL1", "ZLCand_EtaL1", 50, -9, 9, recoZ.daughter(0).eta());
  theHistograms.fill("ZLCand_EtaL2", "ZLCand_EtaL2", 50, -9, 9, recoZ.daughter(1).eta());

  ///*
  if(abs(zeppenfeldllJ0) < 0.6){ //0.72; filter on Zeppenfeld variable
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCuts_wei", "Events after cuts", 13, -0.5, 12.5, cut, theWeight);
  cut2++;
  theHistograms.fill("RecoCuts_pres", "Events after cuts", 13, -0.5, 12.5, cut2);
  //*/

  theHistograms.fill("recoJJ_mass_BJ1pt",     "Jets' mass",        70, 280, 4510, recoJJptot.M(), theWeight);
  theHistograms.fill("recoJJ_deltaEta_BJ1pt", "Jets' #Delta#eta",  50,  -9,    9, recoJJdeltaEta, theWeight);
  
  ///*
  if(Jet1.pt() < 50){ //0 for QCD; filter on second jet's pt
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCuts_wei", "Events after cuts", 13, -0.5, 12.5, cut, theWeight);
  cut2++;
  theHistograms.fill("RecoCuts_pres", "Events after cuts", 13, -0.5, 12.5, cut2);
  //*/
  /*
  theHistograms.fill("recoJJ_mass_B3lpt",     "Jets' mass",        70, 280, 4510, recoJJptot.M(), theWeight);
  theHistograms.fill("recoJJ_deltaEta_B3lpt", "Jets' #Delta#eta",  50,  -9,    9, recoJJdeltaEta, theWeight);

  
  if(recoW.daughter(0).pt() < 30){ //10; filter on 3rd lepton pt
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCuts_wei", "Events after cuts", 13, -0.5, 12.5, cut, theWeight);
  cut2++;
  theHistograms.fill("RecoCuts_pres", "Events after cuts", 13, -0.5, 12.5, cut2);
  */

  theHistograms.fill("recoJJ_mass_BZm",     "Jets' mass",        70, 280, 4510, recoJJptot.M(), theWeight);
  theHistograms.fill("recoJJ_deltaEta_BZm", "Jets' #Delta#eta",  50,  -9,    9, recoJJdeltaEta, theWeight);

  ///*
  if(abs(ZMASS - recoZ.mass()) > 15){ //25; filter on Zmass
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCuts_wei", "Events after cuts", 13, -0.5, 12.5, cut, theWeight);
  cut2++;
  theHistograms.fill("RecoCuts_pres", "Events after cuts", 13, -0.5, 12.5, cut2);
  // */

  recoAfterCut++;
  theHistograms.fill("WZ_Events", "Weighted counters", 10, -0.5, 10.5, 5, theWeight);


  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~ Histograms ~~~~~~~~~~~~~~~~~~~~~~~~~
  // ------------------------------ MET ----------------------------
  theHistograms.fill("recoMET_trmass", "MET's trmass",  100, 25  , 425  , met->p4().Mt(), theWeight);
  theHistograms.fill("recoMET_pt",     "MET's p_{t}",   140,  0  , 700  , MET.pt()      , theWeight);
  theHistograms.fill("recoMET_ID",     "MET's p_{t}",   140,  0  , 700  , MET.id()      , theWeight);
  //theHistograms.fill("recoMET_phi",    "MET's #phi",     50, -3.5,   3.5, MET.phi()     , theWeight);
  
  // ------------------------------- W -----------------------------
  helper_->PlotBoson(recoW, "RecoW", theWeight, "AC");
  
  // ------------------------------- Z -----------------------------
  helper_->PlotBoson(recoZ, "RecoZ", theWeight, "AC");  
  theHistograms.fill("recoZ_pzp_pzm",  "recoZ's leptons' pz", 320, -400, 400, 320, -400, 400, recoZpzlp, recoZpzlm, theWeight);
  
  // ------------------------------ Zl -----------------------------
  theHistograms.fill("recoZl_mass",   "3 leptons mass",          400,  0  , 1200  , recoZlp4.M()          , theWeight);
  theHistograms.fill("recoZl_1st_pt", "Z's 1^{st} lepton p_{t}", 200,  0  ,  400  , recoZ.daughter(0).pt(), theWeight);
  theHistograms.fill("recoZl_2nd_pt", "Z's 2^{nd} lepton p_{t}", 200,  0  ,  400  , recoZ.daughter(1).pt(), theWeight);
  theHistograms.fill("recoZl_3rd_pt", "W's lepton p_{t}",         50,  0  ,  400  , recoW.daughter(0).pt(), theWeight);
  
  // ------------------------------ WZ -----------------------------
  helper_->PlotDiBoson(recoWZs[0], "RecoWZ", theWeight, "AC");
  
  // ----------------------------- Jets ----------------------------
  helper_->PlotJets(Jet0, Jet1, "Reco", theWeight, "AC");
  /*
  foreach(const Particle jet, recoJets){
    theHistograms.fill("recoJet_charge", "Jets' charge",  19, -8.5,  10.5, jet.charge()  , theWeight);
    theHistograms.fill("recoJet_pt",     "Jets' p_{t}",  160,  0  , 400  , jet.pt()      , theWeight);
    theHistograms.fill("recoJet_eta",    "Jets' #eta",    70, -5  ,   5  , jet.eta()     , theWeight);
    theHistograms.fill("recoJet_Y",      "Jets' Y",       80, -5  ,   5  , jet.rapidity(), theWeight);
    theHistograms.fill("recoJet_phi",    "Jets' #phi",    50, -3.5,   3.5, jet.phi()     , theWeight);
  }
  
  theHistograms.fill("recoJet_number", "Jets' number", 10, -0.5, 9.5, recoJets.size(), theWeight);
  */

  // third jet
  /*
  if(recoJets.size() > 2){
    for(int i = 2; i < (int)recoJets.size(); i++){
      theHistograms.fill("recoJet_pt_M3", "Not-leading jets' p_{t}", 64, 30, 350, recoJets[2].pt(), theWeight);
    }
  }
  */
  /*
  // leading jets  
  theHistograms.fill("recoJJ_pt",              "Jets' p_{t}",               200,   0  ,  1200  , recoJJptot.Pt()                        , theWeight);
  theHistograms.fill("recoJJ_ptpt",            "Jets' p_{t}*p_{t}}",        800, 500  , 80500  , Jet0.pt()*Jet1.pt()                    , theWeight);
  theHistograms.fill("recoJJ_pt0_pt1_sum",     "Jets' p_{t}+p_{t}",         320,   0  ,   800  , Jet0.pt()+Jet1.pt()                    , theWeight);
  theHistograms.fill("recoJJ_pt0_pt1_sumquad", "Jets' p_{t}^{2}+p_{t}^{2}", 800,   0  , 80000  , Jet0.pt()*Jet0.pt()+Jet1.pt()*Jet1.pt(), theWeight);
  theHistograms.fill("recoJJ_mass",            "Jets' mass",                  9,  85  ,  4585  , recoJJptot.M()                         , theWeight);
  theHistograms.fill("recoJJ_mass_2",          "Jets' mass",                 45,  85  ,  4585  , recoJJptot.M()                         , theWeight);
  theHistograms.fill("recoJJ_massquad",        "Jets' mass^{2}",            600,   0  , 10000  , recoJJptot.M()*recoJJptot.M()          , theWeight);
  theHistograms.fill("recoJJ_trmass",          "Jets' transverse mass",     300,   0  ,  4500  , recoJJptot.Mt()                        , theWeight);
  theHistograms.fill("recoJJ_deltaEta",        "Jets' #Delta#eta",           50,  -9  ,     9  , recoJJdeltaEta                         , theWeight);
  theHistograms.fill("recoJJ_deltaEta_abs",    "Jets' |#Delta#eta|",         25,   0  ,     8.1, abs(recoJJdeltaEta)                    , theWeight);
  theHistograms.fill("recoJJ_deltaEtaquad",    "Jets' #Delta#eta^{2}",      400, -81  ,    81  , recoJJdeltaEta*recoJJdeltaEta          , theWeight);
  theHistograms.fill("recoJJ_deltaR",          "Jets' #DeltaR",              25,  -0.5,     9  , recoJJdeltaR                           , theWeight);
  theHistograms.fill("recoJJ_deltaPhi",        "Jets' #Delta#phi",           50,  -3.5,     3.5, recoJJdeltaPhi                         , theWeight);
  theHistograms.fill("recoJJ_deltapt",         "Jets' #Deltap_{t}",         600,   0  ,   600  , Jet0.pt()-Jet1.pt()                    , theWeight);
  
  theHistograms.fill("recoJJ_pt0_pt1",           "Jet0's p_{t} (x) vs Jet1's p_{t} (y)",        160,    0,   400, 160,    0  ,   400  , Jet0.pt()          , Jet1.pt()          , theWeight);
  theHistograms.fill("recoJJ_pz0_pz1",           "Jet0's p_{z} (x) vs Jet1's p_{z} (y)",        320, -400,   400, 320, -400  ,   400  , Jet0.p4().Pz()     , Jet1.p4().Pz()     , theWeight);
  theHistograms.fill("recoJJ_J0pt_ptpt",         "Jet0's p_{t} (x) vs Jets' p_{t}^{2} (y)",     160,   20,   420, 800,  500  , 80500  , Jet0.pt()          , Jet0.pt()*Jet1.pt(), theWeight);
  theHistograms.fill("recoJJ_J1pt_ptpt",         "Jet1's p_{t} (x) vs Jets' p_{t}^{2} (y)",     112,   20,   300, 800,  500  , 80500  , Jet1.pt()          , Jet0.pt()*Jet1.pt(), theWeight);
  theHistograms.fill("recoJJ_pt_ptpt",           "Jets' p_{t} (x) vs Jets' p_{t}^{2} (y)",      300,    0,  1000, 800,  500  , 80500  , recoJJptot.Pt()    , Jet0.pt()*Jet1.pt(), theWeight);
  theHistograms.fill("recoJJ_ptpt_deltaEta",     "Jets' p_{t}^{2} (x) vs Jets' #Delta#eta (y)", 800,  500, 80500, 100,   -9  ,     9  , Jet0.pt()*Jet1.pt(), recoJJdeltaEta     , theWeight);
  theHistograms.fill("recoJJ_mass_deltaEta",     "Jets' #Delta#eta (y) vs m_{t} (x)",           300,    0,  1000, 100,   -9  ,     9  , recoJJptot.M()     , recoJJdeltaEta     , theWeight);
  theHistograms.fill("recoJJ_pt_deltaEta",       "Jets' #Delta#eta (y) vs p_{t} (x)",           300,    0,  1000, 100,   -9  ,     9  , recoJJptot.Pt()    , recoJJdeltaEta     , theWeight);
  theHistograms.fill("recoJJ_J0pt_deltaEta",     "Jets' #Delta#eta (y) vs Jet0's p_{t} (x)",    300,    0,  1000, 100,   -9  ,     9  , Jet0.pt()          , recoJJdeltaEta     , theWeight);
  theHistograms.fill("recoJJ_J1pt_deltaEta",     "Jets' #Delta#eta (y) vs Jet1's p_{t} (x)",    240,    0,   800, 100,   -9  ,     9  , Jet1.pt()          , recoJJdeltaEta     , theWeight);
  theHistograms.fill("recoJJ_deltaEta_deltaPhi", "Jets' #Delta#eta (x) vs #Delta#phi (y)",      100,   -9,     9,  50,   -3.5,     3.5, recoJJdeltaEta     , recoJJdeltaPhi     , theWeight);
  */
  // deltaphi when the jet is closer by deltaR
  /*
  if(abs(WZj0DeltaR) < abs(WZj1DeltaR)){
    theHistograms.fill("recoAll_WZj0_vicinoDR_DeltaPhi", "#Delta#phi between recoWZ and Jet0 closer than rJ1", 50, -3.5, 3.5, WZj0DeltaPhi, theWeight);
    theHistograms.fill("recoAll_WZj_vicinoDR_DeltaPhi",  "#Delta#phi between recoWZ and recoJets' jet closer",        50, -3.5, 3.5, WZj0DeltaPhi, theWeight);

    if(recoJets.size() == 2){
      theHistograms.fill("recoAll_WZj0_vicinoDR_DeltaPhi_2jet", "#Delta#phi between recoWZ and Jet0 closer than rJ1", 50, -3.5, 3.5, WZj0DeltaPhi, theWeight);
      theHistograms.fill("recoAll_WZj_vicinoDR_DeltaPhi_2jet",  "#Delta#phi between recoWZ and Jet0 closer than rJ1", 50, -3.5, 3.5, WZj0DeltaPhi, theWeight);
    }
  }
  else{
    theHistograms.fill("recoAll_WZj1_vicinoDR_DeltaPhi", "#Delta#phi between recoWZ and Jet1 closer than rJ0", 50, -3.5, 3.5, WZj1DeltaPhi, theWeight);
    theHistograms.fill("recoAll_WZj_vicinoDR_DeltaPhi",  "#Delta#phi between recoWZ and recoJets' jet closer",        50, -3.5, 3.5, WZj1DeltaPhi, theWeight);   

    if(recoJets.size() == 2){
      theHistograms.fill("recoAll_WZj1_vicinoDR_DeltaPhi_2jet", "#Delta#phi between recoWZ and Jet0 closer than rJ1", 50, -3.5, 3.5, WZj1DeltaPhi, theWeight);
      theHistograms.fill("recoAll_WZj_vicinoDR_DeltaPhi_2jet",  "#Delta#phi between recoWZ and Jet0 closer than rJ1", 50, -3.5, 3.5, WZj1DeltaPhi, theWeight);
    } 
  }
  */

  // deltaphi when the jet is closer by deltaPhi
  /*
  if(abs(WZj0DeltaPhi) < abs(WZj1DeltaPhi)){
    theHistograms.fill("recoAll_WZj0_vicinoDp_DeltaPhi", "#Delta#phi between recoWZ and Jet0 closer than rJ1", 50, -3.5, 3.5, WZj0DeltaPhi, theWeight);
    theHistograms.fill("recoAll_WZj_vicinoDp_DeltaPhi",  "#Delta#phi between recoWZ and recoJets' jet closer",        50, -3.5, 3.5, WZj0DeltaPhi, theWeight);

    if(recoJets.size() == 2){
      theHistograms.fill("recoAll_WZj0_vicinoDp_DeltaPhi_2jet", "#Delta#phi between recoWZ and Jet0 closer than rJ1", 50, -3.5, 3.5, WZj0DeltaPhi, theWeight);
      theHistograms.fill("recoAll_WZj_vicinoDp_DeltaPhi_2jet",  "#Delta#phi between recoWZ and Jet0 closer than rJ1", 50, -3.5, 3.5, WZj0DeltaPhi, theWeight);
    }
  }
  else{
    theHistograms.fill("recoAll_WZj1_vicinoDp_DeltaPhi", "#Delta#phi between recoWZ and Jet1 closer than rJ0", 50, -3.5, 3.5, WZj1DeltaPhi, theWeight);
    theHistograms.fill("recoAll_WZj_vicinoDp_DeltaPhi",  "#Delta#phi between recoWZ and recoJets' jet closer",        50, -3.5, 3.5, WZj1DeltaPhi, theWeight);   

    if(recoJets.size() == 2){
      theHistograms.fill("recoAll_WZj1_vicinoDp_DeltaPhi_2jet", "#Delta#phi between recoWZ and Jet0 closer than rJ1", 50, -3.5, 3.5, WZj1DeltaPhi, theWeight);
      theHistograms.fill("recoAll_WZj_vicinoDp_DeltaPhi_2jet",  "#Delta#phi between recoWZ and Jet0 closer than rJ1", 50, -3.5, 3.5, WZj1DeltaPhi, theWeight);
    } 
  }
  */
  /*
  // ---------------------- WZ & leading jets ----------------------
  //theHistograms.fill("recoAll_mass",           "Mass recoW,Z,J,J",                   150, 220  , 7720  , recoPtot.M()                               , theWeight);
  theHistograms.fill("recoAll_trmass",         "Transverse mass recoW,Z,J",           34, 220  , 7095  , recoPtot.Mt()                              , theWeight);
  theHistograms.fill("recoAll_WZjpt_deltaPhi", "#Delta#phi between recoWZ and Jet0",  50,  -3.5,    3.5, WZjptDeltaPhi                              , theWeight);
  theHistograms.fill("recoAll_Zj0_deltaEta",   "#Delta#eta between recoZ and Jet0" ,  50,  -6  ,    6  , recoZ.eta()-Jet0.eta()                     , theWeight);
  theHistograms.fill("recoAll_Zj1_deltaEta",   "#Delta#eta between recoZ and Jet1" , 100,  -7  ,    7  , recoZ.eta()-Jet1.eta()                     , theWeight);
  theHistograms.fill("recoAll_Zj0_deltaPhi",   "#Delta#phi between recoZ and Jet0" ,  50,  -3.5,    3.5, physmath::deltaPhi(recoZ.phi(), Jet0.phi()), theWeight);
  theHistograms.fill("recoAll_Zj1_deltaPhi",   "#Delta#phi between recoZ and Jet1" ,  50,  -3.5,    3.5, physmath::deltaPhi(recoZ.phi(), Jet1.phi()), theWeight);
  theHistograms.fill("recoAll_Wj0_deltaPhi",   "#Delta#phi between recoW and Jet0" ,  50,  -3.5,    3.5, physmath::deltaPhi(recoW.phi(), Jet0.phi()), theWeight);
  theHistograms.fill("recoAll_Wj1_deltaPhi",   "#Delta#phi between recoW and Jet1" ,  50,  -3.5,    3.5, physmath::deltaPhi(recoW.phi(), Jet1.phi()), theWeight);
  /*
  if(recoJets.size() == 2){
    theHistograms.fill("recoAll_WZjpt_DeltaPhi_2jet", "#Delta#phi between recoWZ and Jet0", 50, -3.5, 3.5, WZjptDeltaPhi, theWeight);
  }
  *//*
  theHistograms.fill("recoAll_trmass_deltaEta",   200, 220  , 8220  , 100, -9  ,    9  , recoPtot.Mt()                              , recoJJdeltaEta                             , theWeight);
  theHistograms.fill("recoAll_trmass_Zltrmass",   200, 220  , 8220  , 400,  0  , 1200  , recoPtot.Mt()                              , recoZlp4.Mt()                              , theWeight);
  theHistograms.fill("recoAll_trmass_JJtrmass",   200, 220  , 8220  , 600,  0  , 4500  , recoPtot.Mt()                              , recoJJptot.Mt()                            , theWeight);
  theHistograms.fill("recoAll_trmass_METtrmass",  200, 220  , 8220  , 400,  0  ,  400  , recoPtot.Mt()                              , MET.p4().Mt()                              , theWeight);
  theHistograms.fill("recoAll_METpt_JJdeltaEta",   70,   0  ,  350  , 100, -9  ,    9  , MET.pt()                                   , recoJJdeltaEta                             , theWeight);
  theHistograms.fill("recoAll_WZpt_JJdeltaEta",   200, 200  , 1000  , 100, -9  ,    9  , recoWZs[0].pt()                            , recoJJdeltaEta                             , theWeight);
  theHistograms.fill("recoAll_Zlmass_JJdeltaEta", 400,   0  , 1200  , 100, -9  ,    9  , recoZlp4.M()                               , recoJJdeltaEta                             , theWeight);
  theHistograms.fill("recoAll_WJ0vsZJ0_deltaPhi",  50,  -3.5,    3.5,  50, -3.5,    3.5, physmath::deltaPhi(recoW.phi(), Jet0.phi()), physmath::deltaPhi(recoZ.phi(), Jet0.phi()), theWeight);
  theHistograms.fill("recoAll_WJ1vsZJ1_deltaPhi",  50,  -3.5,    3.5,  50, -3.5,    3.5, physmath::deltaPhi(recoW.phi(), Jet1.phi()), physmath::deltaPhi(recoZ.phi(), Jet1.phi()), theWeight);
  theHistograms.fill("recoAll_WJ0vsWJ1_deltaPhi",  50,  -3.5,    3.5,  50, -3.5,    3.5, physmath::deltaPhi(recoW.phi(), Jet0.phi()), physmath::deltaPhi(recoW.phi(), Jet1.phi()), theWeight);
  theHistograms.fill("recoAll_ZJ0vsZJ1_deltaPhi",  50,  -3.5,    3.5,  50, -3.5,    3.5, physmath::deltaPhi(recoZ.phi(), Jet0.phi()), physmath::deltaPhi(recoZ.phi(), Jet1.phi()), theWeight);
  theHistograms.fill("recoAll_JJvsZJ0_deltaPhi",   50,  -3.5,    3.5,  50, -3.5,    3.5, recoJJdeltaPhi                             , physmath::deltaPhi(recoZ.phi(), Jet0.phi()), theWeight);
  theHistograms.fill("recoAll_JJvsZJ1_deltaPhi",   50,  -3.5,    3.5,  50, -3.5,    3.5, recoJJdeltaPhi                             , physmath::deltaPhi(recoZ.phi(), Jet1.phi()), theWeight);
  theHistograms.fill("recoAll_JJvsWJ0_deltaPhi",   50,  -3.5,    3.5,  50, -3.5,    3.5, recoJJdeltaPhi                             , physmath::deltaPhi(recoW.phi(), Jet0.phi()), theWeight);
  theHistograms.fill("recoAll_JJvsWJ1_deltaPhi",   50,  -3.5,    3.5,  50, -3.5,    3.5, recoJJdeltaPhi                             , physmath::deltaPhi(recoW.phi(), Jet1.phi()), theWeight);

  // Zeppenfeld
  theHistograms.fill("Zeppenfeld_JJZ",         25, -6, 6   , zeppenfeldJJZ      , theWeight);
  //theHistograms.fill("Zeppenfeld_JJZ_abs",     13,  0, 6.24, abs(zeppenfeldJJZ) , theWeight);
  theHistograms.fill("Zeppenfeld_JJl",         25, -6, 6   , zeppenfeldJJl      , theWeight);
  //theHistograms.fill("Zeppenfeld_JJl_abs",     13,  0, 6.24, abs(zeppenfeldJJl) , theWeight);
  theHistograms.fill("Zeppenfeld_llJ0",        50, -6, 6   , zeppenfeldllJ0     , theWeight);
  theHistograms.fill("Zeppenfeld_llJ0_abs",    26,  0, 6.24, abs(zeppenfeldllJ0), theWeight);
  theHistograms.fill("Zeppenfeld_llJ1",        25, -6, 6   , zeppenfeldllJ1     , theWeight);
  //theHistograms.fill("Zeppenfeld_llJ1_abs",    13,  0, 6.24, abs(zeppenfeldllJ1), theWeight);
  theHistograms.fill("Zeppenfeld_lll",         25, -6, 6   , zeppenfeldlll      , theWeight);
  //theHistograms.fill("Zeppenfeld_lll_abs",     13,  0, 6.24, abs(zeppenfeldlll) , theWeight);
  theHistograms.fill("Zeppenfeld_article",     50, -6, 6   , zeppenfeldArticle  , theWeight);
  theHistograms.fill("Zeppenfeld_article_abs", 26,  0, 6.24, zeppenfeldArticle  , theWeight);
    */

  
  WZ = recoWZs[0];
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~ End of reco Analysis ~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}


void WZAnalyzer::RecoZWCand(DiBosonLepton &recoZW){
  int cut = 0;
  theHistograms.fill("RecoCuts_ZWCand", "Reco events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("RecoCuts_ZWCand_wei", "Reco events after cuts", 14, -0.5, 13.5, cut, theWeight);
  
  DiBosonLepton tempZW;
  
  // ZWCand size > 0
  if(ZW->mass() == 0)
    return;

  cut++;
  theHistograms.fill("RecoCuts_ZWCand", "Reco events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("RecoCuts_ZWCand_wei", "Reco events after cuts", 14, -0.5, 13.5, cut, theWeight);

  tempZW = DiBosonLepton(*ZW);

  // 3 leptons
  if(electrons->size() + muons->size() != 3) //&& ZLCand.size() == 1
    return;
  cut++;
  theHistograms.fill("RecoCuts_ZWCand", "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCuts_ZWCand_wei", "Events after cuts", 13, -0.5, 12.5, cut, theWeight);
  
  // 3rd lepton full selection
  if(!tempZW.second().daughter(0).passFullSel())
    return;

  cut++;
  theHistograms.fill("RecoCuts_ZWCand", "Reco events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("RecoCuts_ZWCand_wei", "Reco events after cuts", 14, -0.5, 13.5, cut, theWeight);

  // 3rd lepton pt > 30
  if(tempZW.second().daughter(0).pt() <30.)
    return;

  cut++;
  theHistograms.fill("RecoCuts_ZWCand", "Reco events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("RecoCuts_ZWCand_wei", "Reco events after cuts", 14, -0.5, 13.5, cut, theWeight);

  // MET pt > 30
  if(tempZW.second().daughter(1).pt() <30.)
    return;

  cut++;
  theHistograms.fill("RecoCuts_ZWCand", "Reco events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("RecoCuts_ZWCand_wei", "Reco events after cuts", 14, -0.5, 13.5, cut, theWeight);

  //Z's mass 60 < m < 120
  if(tempZW.first().mass() < 60. || tempZW.first().mass() > 120.)
    return;

  cut++;
  theHistograms.fill("RecoCuts_ZWCand", "Reco events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("RecoCuts_ZWCand_wei", "Reco events after cuts", 14, -0.5, 13.5, cut, theWeight);

  //jet.size >= 2
  if(jets->size() < 2)
    return;

  cut++;
  theHistograms.fill("RecoCuts_ZWCand", "Reco events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("RecoCuts_ZWCand_wei", "Reco events after cuts", 14, -0.5, 13.5, cut, theWeight);

  // filter on W's trmass 30 < trm < 500
  if(tempZW.second().p4().Mt() < 30. || tempZW.second().p4().Mt() > 500.)
    return;

  cut++;
  theHistograms.fill("RecoCuts_ZWCand", "Reco events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("RecoCuts_ZWCand_wei", "Reco events after cuts", 14, -0.5, 13.5, cut, theWeight);

  // Jets
  Particle Jet0;
  Particle Jet1;

  helper_->FindLeadingJets(jets, Jet0, Jet1);
  
  TLorentzVector recoPtot = tempZW.p4() + Jet0.p4() + Jet1.p4();
  theHistograms.fill("recoAll_ZWCand_trmass_BC", "Transverse mass recoW,Z,J", 34, 220, 7095, recoPtot.Mt(), theWeight);
  
  double recoJJdeltaEta = Jet0.eta() - Jet1.eta();
  TLorentzVector recoJJptot = Jet0.p4() + Jet1.p4();
  float zeppenfeldllJ0 = Jet0.eta() - (tempZW.first().daughter(0).eta() + tempZW.first().daughter(1).eta())/2;
  
  // filter deltaEta JJ > 2.4
  if(abs(recoJJdeltaEta) < 2.4){
    return;
  }

  cut++;
  theHistograms.fill("RecoCuts_ZWCand", "Reco events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("RecoCuts_ZWCand_wei", "Reco events after cuts", 14, -0.5, 13.5, cut, theWeight);

  // filter JJmass > 285.5
  if(recoJJptot.M() < 285.5)
    return;

  cut++;
  theHistograms.fill("RecoCuts_ZWCand", "Reco events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("RecoCuts_ZWCand_wei", "Reco events after cuts", 14, -0.5, 13.5, cut, theWeight);

  
  theHistograms.fill("ZWCand_Zeppenfeld", "ZWCand_Zeppenfeld", 50, -6, 6, zeppenfeldllJ0);
  theHistograms.fill("ZWCand_EtaJ0", "ZWCand_EtaJ0", 50, -9, 9, Jet0.eta());
  theHistograms.fill("ZWCand_EtaL1", "ZWCand_EtaL1", 50, -9, 9, tempZW.first().daughter(0).eta());
  theHistograms.fill("ZWCand_EtaL2", "ZWCand_EtaL2", 50, -9, 9, tempZW.first().daughter(1).eta());

  // filter Zeppenfeld > 0.6
  if(abs(zeppenfeldllJ0) < 0.6)
    return;

  cut++;
  theHistograms.fill("RecoCuts_ZWCand", "Reco events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("RecoCuts_ZWCand_wei", "Reco events after cuts", 14, -0.5, 13.5, cut, theWeight);

  // filter Jet1 pt > 50
  if(Jet1.pt() < 50.)
    return;

  cut++;
  theHistograms.fill("RecoCuts_ZWCand", "Reco events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("RecoCuts_ZWCand_wei", "Reco events after cuts", 14, -0.5, 13.5, cut, theWeight);

  // filter Zmass < 15
  if(abs(ZMASS - tempZW.first().mass() > 15))
    return;

  cut++;
  theHistograms.fill("RecoCuts_ZWCand", "Reco events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("RecoCuts_ZWCand_wei", "Reco events after cuts", 14, -0.5, 13.5, cut, theWeight);
  
  recoZW = tempZW;
  
  theHistograms.fill("recoAll_ZWCand_trmass_AC", "Transverse mass recoW,Z,J", 34, 220, 7095, recoPtot.Mt(), theWeight);
}


void WZAnalyzer::GenRecoAnalysis(const DiBosonParticle genWZ, const Particle genJet0, const Particle genJet1, const DiBosonLepton recoWZ, const Particle recoJet0, const Particle recoJet1){
  eventGenReco++;
  theHistograms.fill("WZ_Events", "Weighted counters", 10, -0.5, 10.5, 6, theWeight);
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~ Begin of Reco vs Gen ~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /*
  bool choosingW = kTRUE;
  BosonLepton recoWtemp;
  TLorentzVector genWZjjp4 = genWZ.first().p4() + genWZ.second().p4() + genJet0.p4() + genJet1.p4();
  TLorentzVector recoWZjjp4 = recoWZ.first().p4() + recoWZ.second().p4() + recoJet0.p4() + recoJet1.p4();
  vector<DiBosonLepton> recoWZs;
  ZLCompositeCandidates recoZls;
 
  // test if recoW is genW by trmass
  foreach(const ZLCompositeCandidate Z, *ZLCand){
    if(Z.first.mass() > 60. && Z.second.mass() < 120. && Z.second.passFullSel()){
      recoZls.push_back(Z);
    }
  }

  if(recoZls.size() > 0){
    sort(recoZls.begin(), recoZls.end(), pairMassComparator(0, ZMASS));
    
    for(int i = 0; choosingW; i++){
      recoWtemp = BosonLepton(recoZls[i].second, Lepton(met->p4()), copysign(24, recoZls[i].second.charge()));

      if(recoWtemp.p4().Mt() > 30 && recoWtemp.p4().Mt() < 500){
	recoWZs.push_back(DiBosonLepton(recoWtemp, recoZls[i].first));
      }
      
      if(i + 1 == (int)recoZls.size() || recoZls[i].first.mass() - recoZls[i + 1].first.mass() <= 0.1){
	choosingW = kFALSE;
      }
    }
  }
  else{
    return;
  }

  if(recoWZs.size() == 0){
    return;
  }
  
  sort(recoWZs.begin(), recoWZs.end(), pairTrmassComparator(0, genWZ.first().p4().Mt()));
  */
  
  // check if gen and reco IDs are the same
  int genZID = abs(genWZ.second().daughter(0).id()) + abs(genWZ.second().daughter(1).id());
  int recoZID = abs(recoWZ.second().daughter(0).id()) + abs(recoWZ.second().daughter(1).id());
  
  // check if gen and reco IDs are the same
  int genWZID = abs(genWZ.first().daughter(0).id()) + abs(genWZ.second().daughter(0).id());
  int recoWZID = abs(recoWZ.first().daughter(0).id()) + abs(recoWZ.second().daughter(0).id());

  if(recoZID != genZID){
    counter5++;
  }


  bool WdifferentID = (abs(genWZ.first().daughter(0).id()) != abs(recoWZ.first().daughter(0).id()));
  bool ZdifferentID = (abs(genWZ.second().daughter(0).id()) != abs(recoWZ.second().daughter(0).id()));
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~ Histograms ~~~~~~~~~~~~~~~~~~~~~~~~~
  /*
  // genN vs MET
  theHistograms.fill("GR_MET_Deltatrmass", "#Deltam_{t} between neutrino and MET", 600, -300, 300, genWZ.first().daughter(1).p4().Mt() - recoWZ.first().daughter(1).p4().Mt(), theWeight);
  
  // genW vs recoW
  theHistograms.fill("GR_W_Deltatrmass",  "W's #Deltam_{t} in Gen&Reco events",  500, -250, 250, genWZ.first().p4().Mt() - recoWZ.first().p4().Mt(), theWeight);
  
  //theHistograms.fill("GR_recoWtm_trmass", "recoW by tr mass",                600,  0  , 600  , recoWZs[0].first().p4().Mt(), theWeight);
  //theHistograms.fill("GR_recoWtm_WvsW",   "recoW and recoWtm are the same?",   4, -1.5,   2.5, isTheSame(recoWZ.first()    , recoWZs[0].first()));

  // genZ vs recoZ
  theHistograms.fill("GR_Z_Deltatrmass",  "Z's #Deltam_{t} in Gen&Reco events",  200, -100, 100, genWZ.second().p4().Mt() - recoWZ.second().p4().Mt(), theWeight);

  // genWZ vs recoWZ
  theHistograms.fill("GR_WZ_Deltatrmass", "WZ's #Deltam_{t} in Gen&Reco events", 600, -300  , 300  , genWZ.p4().Mt() - recoWZ.p4().Mt(), theWeight);
  //theHistograms.fill("GR_WZ_sameID",      "Are genWZ and recoWZ the same?",         4,   -1.5,   2.5, genZlID == recoZlID               , theWeight);

  // genWZJJ vs recoWZJJ
  theHistograms.fill("GR_WZJJ_Deltatrmass", "WZ and jets' #Deltam_{t} in Gen&Reco events", 600, -600, 600, genWZjjp4.Mt() - recoWZjjp4.Mt(), theWeight);

  */
  // genZ vs recoZ ID
  theHistograms.fill("GR_ID_genZ_vs_recoZ", "GenZ's and RecoZ's ID", 5, 19, 29, 5, 19, 29, genZID, recoZID);  
  theHistograms.fill("GR_ID_genWZ_vs_recoWZ", "GenWZ and RecoWZ's daughter's ID", 5, 19, 29, 5, 19, 29, genWZID, recoWZID);  
  if(genWZID != recoWZID)
    theHistograms.fill("GR_ID_genWZ_vs_recoWZ_2", "GenWZ and RecoWZ's daughter's ID", 5, 19, 29, 5, 19, 29, genWZID, recoWZID);  
  
  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~ End of Reco vs Gen ~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}


void WZAnalyzer::CheckBuildWZ(){
  vector<Particle> electron;
  vector<Particle> lepton;
  vector<Particle> muon;
  vector<Particle> neutrino;
  
  foreach(const Particle &gen, *genParticles){
    if((abs(gen.id()) != 11 && abs(gen.id()) != 13 && abs(gen.id()) != 12 && abs(gen.id()) != 14) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) 
      continue;
    
    if(abs(gen.id()) == 11) electron.push_back(gen);
    else if(abs(gen.id()) == 13) muon.push_back(gen);
    else if(abs(gen.id()) == 12 || abs(gen.id()) == 14) neutrino.push_back(gen);
  }

  // ----- filter on leptons' number
  if(electron.size() + muon.size() + neutrino.size() != 4 && neutrino.size() != 1){
    return;
  }
  
  // ----- filters on leptons' pt and eta
  foreach(const Particle ele, electron){
    if(ele.pt() < 7 || abs(ele.eta()) > 2.5){
      return;
    }
    lepton.push_back(ele);
  }
  
  foreach(const Particle mu, muon){
    if(mu.pt() < 5 || abs(mu.eta()) > 2.4){
      return;
    }
    lepton.push_back(mu);
  }

  if(lepton.size() != 3){
    return;
  }
  
  sort(electron.begin(), electron.end(), PtComparator());
  sort(lepton.begin(), lepton.end(), PtComparator());
  sort(muon.begin(), muon.end(), PtComparator());
  
  if(lepton[0].pt() < 20){
    return;
  }
  
  if(abs(lepton[1].id()) == 11 && lepton[1].pt() < 12){
    return;
  }

  if(abs(lepton[1].id()) == 13 && lepton[1].pt() < 10){
    return;
  }
  
  // ----- filter on total mass: Z and W must be at least on shell  
  TLorentzVector Ptot = neutrino[0].p4();
  
  foreach(const Particle lep, lepton)
    Ptot += lep.p4();
  
  theHistograms.fill("AllGenlllnu_mass",   "m 3 leptons and #nu",     150, 0, 1500, Ptot.M() , theWeight);
  theHistograms.fill("AllGenlllnu_trmass", "m_{T} 3 leptons and #nu", 150, 0, 1500, Ptot.Mt(), theWeight);
  
  if(Ptot.M() < 150.){
    return;
  }// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Z & W ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BosonParticle Weh;
  BosonParticle Zet;
  vector<pairBosonParticle> Zls;
  vector<DiBosonParticle> WZs;
  bool mixed = false;
  
  // ---------------------- W & Z construction ---------------------
  // case 2e 1m
  if(electron.size()==2 && muon.size()==1 && electron[0].charge() != electron[1].charge()){
    gen2e1m++;
    mixed = true;/*
    Zet = BosonParticle(electron[0], electron[1], 23);
    Zls.push_back(pairBosonParticle(Zet, muon[0])); // useless if not cut on Zls size
    Weh = BosonParticle(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()) );*/
  }
  
  // case 1e 2m
  if(electron.size()==1 && muon.size()==2 && muon[0].charge() != muon[1].charge()){
    gen2m1e++;
    mixed = true;/*
    Zet = BosonParticle(muon[0], muon[1], 23);
    Zls.push_back(pairBosonParticle(Zet, electron[0])); // useless if not cut on Zls size
    Weh = BosonParticle(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()) );*/
  }

  
  // ~~~~~ Mixed cases, test: see if better identify W from Z or Z from W
  if(mixed){
    double differenceZ = 0;
    double differenceW = 0;

    for(int i = 0; i < (int)lepton.size() -1; i++){
      for(int j = i + 1; j < (int)lepton.size(); j++){
	for(int k = 0; k < (int)lepton.size(); k++){
	  if(k != i && k != j){
	    if(lepton[i].charge() != lepton[j].charge()){
	      Zls.push_back(pairBosonParticle(BosonParticle(lepton[i], lepton[j], 23), lepton[k]));
	    }
	  }
	}
      }
    }
    
    if(Zls.size() < 1){   
      return;
    }
    
    // Z is made up of the couple which gives a better Zmass 
    sort(Zls.begin(), Zls.end(), pairMassComparator(0, ZMASS));
    differenceZ = fabs(ZMASS - Zls[0].first.p4().M());

    // check if choosen Z has different daughters
    if(abs(Zls[0].first.daughter(0).id()) != abs(Zls[0].first.daughter(1).id())){
	choosedZwrongID++;
    }

    for(int i = 0; i < (int)Zls.size(); i++){
      WZs.push_back(DiBosonParticle(BosonParticle(neutrino[0], Zls[i].second, copysign(24, Zls[i].second.charge())), Zls[i].first));
    }

    // check if choosen Z&W have different mass
    if(abs(WMASS - WZs[0].first().mass()) > 30. || abs(ZMASS - WZs[0].second().mass()) > 30.){
      choosedZoutsiderange++;
    }
    
    // W is made up of the couple which gives a better Wmass 
    sort(WZs.begin(), WZs.end(), pairMassComparator(0, WMASS));
    differenceW = fabs(WMASS - WZs[0].first().p4().M());
    if(abs(WZs[0].second().daughter(0).id()) != abs(WZs[0].second().daughter(1).id())){
	choosedWwrongID++;
    }
    
    if(abs(WMASS - WZs[0].first().mass()) > 30. || abs(ZMASS - WZs[0].second().mass()) > 30.){
      choosedWoutsiderange++;
    }
    
    // Best couple has less difference in mass from main boson
    if(differenceZ < differenceW){ // Z chosen first
      Zet = Zls[0].first;
      Weh = BosonParticle(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()));

      choosedZfirst++;

      if(abs(Zet.daughter(0).id()) != abs(Zet.daughter(1).id())){
	counter1++;
      }
      
      if(Weh.mass() < 50.00 || Weh.mass() > 110.00 || Zet.mass() < 60.00 || Zet.mass() > 120.00){
	counter3++;
	//return;
      }
      
    } else{ // W chosen first
      Weh = WZs[0].first();
      Zet = WZs[0].second();

      choosedWfirst++;

      if(abs(Zet.daughter(0).id()) != abs(Zet.daughter(1).id())){
	counter2++;
      }
    
      if(Weh.mass() < 50.00 || Weh.mass() > 110.00 || Zet.mass() < 60.00 || Zet.mass() > 120.00){
	counter3++;
	//return;
      }
    }

    if(isTheSame(Zls[0].first, WZs[0].second())){
      counter6++;
      if(abs(Zet.daughter(0).id()) != abs(Zet.daughter(1).id()))
	sicheso5++;
    }
    else {
      if(differenceZ < differenceW){
	sicheso1++;
	if(abs(Zet.daughter(0).id()) != abs(Zet.daughter(1).id()))
	  sicheso3++;
      }
      else{
	sicheso2++;
	if(abs(Zet.daughter(0).id()) != abs(Zet.daughter(1).id()))
	  sicheso4++;
      }
    }    
    
  }
  
  
  // ----- choosing Z for non mixed cases
  // case 3e
  if(electron.size()==3){
    gen3e++;/*
    for(int i = 0; i < (int)electron.size() -1; i++){
      for(int j = i + 1; j < (int)electron.size(); j++){
	for(int k = 0; k < (int)electron.size(); k++){
	  if(k != i && k != j){
	    if(electron[i].charge() != electron[j].charge()){
	      Zls.push_back(pairBosonParticle(BosonParticle(electron[i], electron[j], 23), electron[k]));
	    }
	  }
	}
      }
    }
    
    // Z is made up of the couple which gives a better Zmass 
    sort(Zls.begin(), Zls.end(), pairMassComparator(0, ZMASS));
    Zet = Zls[0].first;
    
    // W is made up of the remaining lepton and the neutrino
    Weh = BosonParticle(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()));*/
  }
  
  // case 3m
  else if(muon.size()==3){
    gen3m++;/*
    for(int i = 0; i < (int)muon.size() -1; i++){
      for(int j = i + 1; j < (int)muon.size(); j++){
	for(int k = 0; k < (int)muon.size(); k++){
	  if(k != i && k != j){
	    if(muon[i].charge() != muon[j].charge()){
	      Zls.push_back(pairBosonParticle(BosonParticle(muon[i], muon[j], 23), muon[k]));
	    }
	  }
	}
      }
    }
    
    sort(Zls.begin(), Zls.end(), pairMassComparator(0, ZMASS));
    Zet = Zls[0].first;
    Weh = BosonParticle(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()));*/
  }
  /*
  if(!mixed){
    double differenceZ = 0;
    double differenceW = 0;

    for(int i = 0; i < (int)lepton.size() -1; i++){
      for(int j = i + 1; j < (int)lepton.size(); j++){
	for(int k = 0; k < (int)lepton.size(); k++){
	  if(k != i && k != j){
	    if(lepton[i].charge() != lepton[j].charge()){
	      Zls.push_back(pairBosonParticle(BosonParticle(lepton[i], lepton[j], 23), lepton[k]));
	    }
	  }
	}
      }
    }
    
    // Z is made up of the couple which gives a better Zmass 
    sort(Zls.begin(), Zls.end(), pairMassComparator(0, ZMASS));
    differenceZ = fabs(ZMASS - Zls[0].first.p4().M());

    for(int i = 0; i < (int)Zls.size(); i++){
      WZs.push_back(DiBosonParticle(BosonParticle(Zls[i].second, neutrino[0], copysign(24, Zls[i].second.charge())), Zls[i].first));
    }
    
    // W is made up of the couple which gives a better Wmass 
    sort(WZs.begin(), WZs.end(), pairMassComparator(0, WMASS));
    differenceW = fabs(WMASS - WZs[0].first().p4().M());
    
    // Best couple has less difference in mass from main boson
    if(differenceZ < differenceW){ // Z chosen first
      Zet = Zls[0].first;
      Weh = BosonParticle(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()));
      counter1++;      
    } else{ // W chosen first
      Weh = WZs[0].first();
      Zet = WZs[0].second();
      counter2++;
    }
    
    if(Weh.mass() < 50.00 || Weh.mass() > 110.00 || Zet.mass() < 60.00 || Zet.mass() > 120.00){
      counter3++;
      //return;
    }

    if(isTheSame(Zls[0].first, WZs[0].second())){
      counter6++;
    }
  }
  */
}


void WZAnalyzer::analyze(){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~ Main analysis ~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  eventSample++;
  theHistograms.fill("WZ_Events", "Weighted counters", 10, -0.5, 9.5, 1, theWeight);
  
  //gen variables	
  DiBosonParticle genWZ;
  Particle genJet0;
  Particle genJet1;
  
  //reco variables
  DiBosonLepton recoWZ;
  DiBosonLepton recoZW;
  Jet recoJet0;
  Jet recoJet1;

  // Provvisoria
  WZAnalyzer::CheckBuildWZ();
  
  //Gen analysis
  WZAnalyzer::GenAnalysis(genWZ, genJet0, genJet1);
  
  //Reco analysis
  //if(genWZ.pt() != 0){
    WZAnalyzer::RecoAnalysis(recoWZ, recoJet0, recoJet1);
  //}

    WZAnalyzer::RecoZWCand(recoZW);

  ///*
  if(genWZ.pt() != 0. && recoWZ.pt() != 0.){
    //Reco vs Gen analysis
    WZAnalyzer::GenRecoAnalysis(genWZ, genJet0, genJet1, recoWZ, recoJet0, recoJet1);
  }
  //*/

  if(genWZ.pt() != 0. && recoWZ.pt() == 0.){
    //Gen event not reconstructed for reasons
    eventGenNOReco++;
    theHistograms.fill("WZ_Events", "Weighted counters", 10, -0.5, 9.5, 7, theWeight);
  }

  if(genWZ.pt() == 0. && recoWZ.pt() != 0.){
    //Reco event not generated for reasons
    eventRecoNOGen++;
    theHistograms.fill("WZ_Events", "Weighted counters", 10, -0.5, 9.5, 8, theWeight);
  }
}


void WZAnalyzer::end(TFile &){    
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~ End banner ~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout << "\n--------------------------------------------------------------------------" << endl;
  
  cout << "\nEvents of the sample analyzed:                       " << setw(9) << eventSample << endl;
  cout << "Gen events analyzed:                                 " << setw(9) << eventGen << "\t" << eventGen*100./eventSample << "%" << "\t" << eventGen*theWeight*1.0 << endl;
  cout << "Gen events after all cuts:                           " << setw(9) << genAfterCut << "\t" << genAfterCut*100./eventGen << "%" << "\t" << genAfterCut*theWeight*1.0 << endl;
  cout << "Reco events analyzed:                                " << setw(9) << eventReco << "\t" << eventReco*100./eventSample << "%" << "\t" << eventReco*theWeight*1.0 << endl;
  /*
  cout << "RecoZls empty:                                       " << setw(9) << recoZlempty << "\t" << recoZlempty*100./eventSample << "%" << endl;
  cout << "RecoJets.size < 2:                                   " << setw(9) << recoJetless2 << "\t" << recoJetless2*100./recoZlempty << "%" << endl;
  */
  cout << "Reco events after all cuts:                          " << setw(9) << recoAfterCut << "\t" << recoAfterCut*100./eventReco << "%" << "\t" << recoAfterCut*theWeight*1.0 << endl;
  cout << "Gen&Reco events analyzed:                            " << setw(9) << eventGenReco << "\t" << eventGenReco*100./eventSample << "%" << "\t" << eventGenReco*theWeight*1.0 << endl;
  cout << "Gen events not reconstructed:                        " << setw(9) << eventGenNOReco << "\t" << eventGenNOReco*100./eventGen << "%" << "\t" << eventGenNOReco*theWeight*1.0 << endl;
  cout << "Reco events that weren't gen events:                 " << setw(9) << eventRecoNOGen << "\t" << eventRecoNOGen*100./eventReco << "%" << "\t" << eventRecoNOGen*theWeight*1.0 << endl;

  ///*
  cout << "\nNumber of couples reconstructed per type: " << endl;
  cout << "- Generated:     3 electrons        " << setw(9) << gen3e    << "\n                 3 muons            " << setw(9) << gen3m << endl;
  cout << "                 2 electrons 1 muon " << setw(9) << gen2e1m  << "\n                 2 muons 1 electron " << setw(9) << gen2m1e << endl;
  cout << "- Reconstructed: 3 electrons        " << setw(9) << reco3e   << "\n                 3 muons            " << setw(9) << reco3m << endl;
  cout << "                 2 electrons 1 muon " << setw(9) << reco2e1m << "\n                 2 muons 1 electron " << setw(9) << reco2m1e << endl;
  cout << "Number of W and/or Z with mass out of range " << setw(9) << counter4 << endl;
  //*/
  

  int total = gen2e1m+gen2m1e; //gen3e+gen3m

  if(total != 0){
    cout << "\n When event has mixed flavour leptons: " << endl;
    cout << "GEN PARTICLES: (total events == mixed events for numbers and percentage below)";
    cout << "\nTotal number of mixed events                                  " << setw(9) << total << endl;
    cout << "Z identified first                                            " << setw(9) << choosedZfirst << " -> " << choosedZfirst*100./total << "%" << endl;
    cout << "  with Z's wrong daughters                                    " << setw(9) << counter1 << " -> " << counter1*100./(total) << "%" << endl;
    cout << "W identified first                                            " << setw(9) << choosedWfirst << " -> " << choosedWfirst*100./total << "%" << endl;
    cout << "  with Z's wrong daughters                                    " << setw(9) << counter2 << " -> " << counter2*100./(total) << "%" << endl;
    cout << "Total number of Z whoose daughters are with different flavour " << setw(9) << counter1+counter2 << " -> " << (counter2+counter1)*100./(total) << "%" << endl;
    cout << "Mixed events with W's and/or Z's mass out of range            " << setw(9) << counter3 << " -> " << counter3*100./(total) << "% (cases with different flavour daughters are included)" << endl;
    cout << "Events in which Z or W identification are the same            " << setw(9) << counter6 << " -> " << counter6*100./(total) << "%" << endl;
    cout << "  Z has wrong daughters                                       " << setw(9) << sicheso5  << " -> " << (sicheso5)*100./(counter6) << "%" << endl;
    cout << "Events in which Z or W identification aren't the same         " << setw(9) << sicheso1 + sicheso2  << " -> " << (sicheso1 + sicheso2)*100./(total) << "%" << endl;
    cout << "  Z chosen first                                              " << setw(9) << sicheso1  << " -> " << (sicheso1)*100./(sicheso1 + sicheso2) << "%" << endl;
    cout << "    Z has wrong daughters                                     " << setw(9) << sicheso3  << " -> " << (sicheso3)*100./(sicheso1) << "%" << endl;
    cout << "  W chosen first                                              " << setw(9) << sicheso2  << " -> " << (sicheso2)*100./(sicheso1 + sicheso2) << "%" << endl;
    cout << "    Z has wrong daughters                                     " << setw(9) << sicheso4  << " -> " << (sicheso4)*100./(sicheso2) << "%" << endl;

    cout << "\n If Z were always choosen first:" << endl;
    cout << "- Z with wrong daughters       " << setw(9) << choosedZwrongID << " -> " << choosedZwrongID*100./total << "%" << endl;
    cout << "- Z and/or W out of mass range " << setw(9) << choosedZoutsiderange << " -> " << choosedZoutsiderange*100./total << "%"  << endl;
    cout << "\n If W were always choosen first:" << endl;
    cout << "- Z with wrong daughters       " << setw(9) << choosedWwrongID << " -> " << choosedWwrongID*100./total << "%" << endl;
    cout << "- Z and/or W out of mass range " << setw(9) << choosedWoutsiderange << " -> " << choosedWoutsiderange*100./total << "%"  << endl;

    cout << "\nRECO vs GEN PARTICLES: ";
    cout << "\nNumber of Zgen and Zreco with different ID: " << setw(9) << counter5 << " -> "<< counter5*100./eventGenReco << "%" << endl;
  }  
  
  /*
  cout << "\nNumber of                                           " << setw(9) << counter1 << endl;
  cout << "Number of                                            " << setw(9) << counter2 << endl;
  cout << "Number of                                            " << setw(9) << counter3 << endl;
  cout << "Number of                                            " << setw(9) << counter4 << endl;
  cout << "Number of                                            " << setw(9) << counter5 << endl;
  */
    
  // execution time
  endtime = ((float)clock())/CLOCKS_PER_SEC;
  helper_->printTime(begintime, endtime);
  cout << "\n--------------------------------------------------------------------------" << endl;
}
