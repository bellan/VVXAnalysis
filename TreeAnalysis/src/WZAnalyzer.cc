#include "VVXAnalysis/TreeAnalysis/interface/WZAnalyzer.h"
#include "VVXAnalysis/Commons/interface/AriEle.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 

using namespace boost::assign;
using namespace colour;
using namespace physmath;
using namespace phys;
using namespace std;

void WZAnalyzer::begin() {

  //events counters
  eventGen = 0;
  eventReco = 0;
  eventGenReco = 0;
  eventSample = 0;
  
  recoZlempty = 0;
  recoJetless2 = 0;

  //free counters
  counter1 = 0;
  counter2 = 0;
  counter3 = 0;
  counter4 = 0;

  begintime = ((float)clock())/CLOCKS_PER_SEC;
}

Int_t WZAnalyzer::cut() {
  return 1;
}

void WZAnalyzer::GenAnalysis(ZZtype &WZ, Particle &Jet0, Particle &Jet1){
  // ~~~~~~~~~~~~~~~~~~~~ Begin of gen Analysis ~~~~~~~~~~~~~~~~~~~~
  
  vector<Particle> electron;
  vector<Particle> genjets;
  vector<Particle> lepton;
  vector<Particle> muon;
  vector<Particle> neutrino;
  
  foreach(const Particle &gen, *genParticles){
    if((abs(gen.id()) != 11 && abs(gen.id()) != 13 && abs(gen.id()) != 12 && abs(gen.id()) != 14) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) 
      continue;
    
    if(abs(gen.id()) == 11)      electron.push_back(gen);
    else if(abs(gen.id()) == 13) muon.push_back(gen);
    else if(abs(gen.id()) == 12 || abs(gen.id()) == 14)  neutrino.push_back(gen);
  }
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Filters ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // filter on jets' number
  foreach(const Particle jet, *genJets){ //in genJets are included leptons
    bool leptonMatch = false; 
    foreach(const phys::Particle &gen, *genParticles){
      if(physmath::deltaR(gen,jet) < 0.4 && (abs(gen.id()) == 11 || abs(gen.id()) == 13)) leptonMatch = true;
    }
    if(!leptonMatch){
      if(fabs(jet.eta()) < 4.7) genjets.push_back(jet);
    }
  }
  
  if(genjets.size() < 2){
    return;
  }
  
  // filter on leptons' number
  if(electron.size() + muon.size() + neutrino.size() != 4){
    return;
  }
  
  // filters on leptons' pt and eta
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
  
  // Z and W must be at least on shell
  TLorentzVector Ptot = neutrino[0].p4();
  foreach(const Particle lep, lepton)
    Ptot += lep.p4();
  
  theHistograms.fill("AllGenlllnu_mass",   "m 3 leptons and #nu",     300, 0, 1500, Ptot.M() , theWeight);
  theHistograms.fill("AllGenlllnu_trmass", "m_{T} 3 leptons and #nu", 300, 0, 1500, Ptot.Mt(), theWeight);
  
  if(Ptot.M() < 165){
    return;
  }
  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Z & W ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  eventGen++;
  
  Vtype Weh;
  Vtype Zet;
  vector<Zltype> Zls;
  
  // ---------------------- W & Z construction ---------------------
  // case 2e 1m
  if(electron.size()==2 && muon.size()==1 && electron[0].charge() != electron[1].charge()){
    Zet = Vtype(electron[0], electron[1], 23);
    Zls.push_back(Zltype(Zet, muon[0]));
    Weh = Vtype(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()) );
  }
  
  // case 1e 2m
  if(electron.size()==1 && muon.size()==2 && muon[0].charge() != muon[1].charge()){
    Zet = Vtype(muon[0], muon[1], 23);
    Zls.push_back(Zltype(Zet, electron[0]));
    Weh = Vtype(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()) );
  }
  
  // case 3e
  if(electron.size()==3){
    for(int i = 0; i < (int)electron.size() -1; i++){
      for(int j = i + 1; j < (int)electron.size(); j++){
	for(int k = 0; k < (int)electron.size(); k++){
	  if(k != i && k != j){
	    if(electron[i].charge() != electron[j].charge()){
	      Zls.push_back(Zltype(Vtype(electron[i], electron[j], 23), electron[k]));
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
    Zet = Zls[0].first;
    
    // W is made up of the remaining lepton and the neutrino
    Weh = Vtype(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()));
  }
  
  // case 3m
  else if(muon.size()==3){
    for(int i = 0; i < (int)muon.size() -1; i++){
      for(int j = i + 1; j < (int)muon.size(); j++){
	for(int k = 0; k < (int)muon.size(); k++){
	  if(k != i && k != j){
	    if(muon[i].charge() != muon[j].charge()){
	      Zls.push_back(Zltype(Vtype(muon[i], muon[j], 23), muon[k]));
	    }
	  }
	}
      }
    }
    
    if(Zls.size() < 1){
      return;
    }
    
    sort(Zls.begin(), Zls.end(), pairMassComparator(0, ZMASS));
    Zet = Zls[0].first;
    Weh = Vtype(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()));
  }
  
  // W&Z diboson
  WZ = ZZtype(Weh, Zet);
  
  // ------------------------ W & Z variables ----------------------
  // Z
  //double ZdeltaEta = Zet.daughter(0).eta() - Zet.daughter(1).eta();
  //double ZdeltaPhi = physmath::deltaPhi(Zet.daughter(0).phi(), Zet.daughter(1).phi());
  //double ZdeltaR = abs(physmath::deltaR(Zet.daughter(0), Zet.daughter(1)));
  
  // W
  //double WdeltaEta = Weh.daughter(0).eta() - Weh.daughter(1).eta();
  //double WdeltaPhi = physmath::deltaPhi(Weh.daughter(0).phi(), Weh.daughter(1).phi());
  //double WdeltaR = abs(physmath::deltaR(Weh.daughter(0), Weh.daughter(1)));
  
  // Zl
  //int ZlsID = abs(Zls[0].first.daughter(0).id()) + abs(Zls[0].first.daughter(1).id()) + abs(Zls[0].second.id());
  
  // WZ
  //double WZdeltaEta = Zet.eta() - Weh.eta();
  //double WZdeltaPhi = physmath::deltaPhi(Zet.phi(), Weh.phi());
  //double WZdeltaR = abs(physmath::deltaR(Zet, Weh));
  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Jets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // leading jets are those with higher pt
  sort(genjets.begin(), genjets.end(), PtComparator());
  Jet0 = genjets[0];
  Jet1 = genjets[1];
  
  // leading jets
  TLorentzVector JJp4 = genjets[0].p4() + genjets[1].p4();
  double JJdeltaEta = genjets[0].eta() - genjets[1].eta();
  double JJdeltaPhi = physmath::deltaPhi(genjets[0].phi(), genjets[1].phi());
  double JJdeltaR = abs(physmath::deltaR(genjets[0], genjets[1]));
  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ All ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TLorentzVector WZjjp4 = WZ.p4() + genjets[0].p4() + genjets[1].p4();
  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~ Histograms ~~~~~~~~~~~~~~~~~~~~~~~~~
  /*
  // Nu
  theHistograms.fill("GenN_charge", "#nu's charge",  5, -2.5,   2.5, neutrino[0].charge()  , theWeight); //just to be sure...
  theHistograms.fill("GenN_pt",     "#nu's p_{t}", 140,  0  , 700  , neutrino[0].pt()      , theWeight);
  theHistograms.fill("GenN_Y",      "#nu's Y",      90, -6.5,   6.5, neutrino[0].rapidity(), theWeight);
  theHistograms.fill("GenN_eta",    "#nu's #eta",   90, -6.5,   6.5, neutrino[0].eta()     , theWeight);
  theHistograms.fill("GenN_phi",    "#nu's #phi",   50, -4  ,   4  , neutrino[0].phi()     , theWeight);
  */
  
  // W
  //theHistograms.fill("GenW_charge",       "W's charge",       5, -2.5,   2.5, Weh.charge()  , theWeight);
  theHistograms.fill("GenW_mass",         "W's mass",       400,  0  , 400  , Weh.mass()    , theWeight);
  theHistograms.fill("GenW_trmass",       "W's trmass",     400,  0  , 400  , Weh.p4().Mt() , theWeight);
  theHistograms.fill("GenW_pt",           "W's p_{t}",      130,  0  , 650  , Weh.pt()      , theWeight);
  //theHistograms.fill("GenW_Y",            "W's Y",           50, -5  ,   5  , Weh.rapidity(), theWeight);
  //theHistograms.fill("GenW_eta",          "W's #eta",        50, -9  ,   9  , Weh.eta()     , theWeight);
  //theHistograms.fill("GenW_phi",          "W's #phi",        50, -4  ,   4  , Weh.phi()     , theWeight);
  //theHistograms.fill("GenW_deltaEta",     "W's #Delta#eta",  50, -6.5,   6.5, WdeltaEta     , theWeight);
  //theHistograms.fill("GenW_deltaR",       "W's #DeltaR",     50, -0.5,   6.5, WdeltaR       , theWeight);
  //theHistograms.fill("GenW_deltaPhi",     "W's #Delta#phi",  50, -4  ,   4  , WdeltaPhi     , theWeight);
  
  theHistograms.fill("GenW_massvstrmass", "W's mass(x) vs trmass(y)", 400, 0, 400, 400, 0, 400, Weh.mass(), Weh.p4().Mt(), theWeight);
  
  // Z
  //theHistograms.fill("GenZ_charge",   "Z's charge",       5, -2.5,   2.5, Zet.charge()  , theWeight); //just to be sure...
  theHistograms.fill("GenZ_mass",     "Z's mass",        30, 60  , 120  , Zet.mass()    , theWeight);
  theHistograms.fill("GenZ_trmass",   "Z's trmass",      30, 60  , 120  , Zet.p4().Mt() , theWeight);
  theHistograms.fill("GenZ_pt",       "Z's p_{t}",      130,  0  , 650  , Zet.pt()      , theWeight);
  //theHistograms.fill("GenZ_Y",        "Z's Y",           45, -3  ,   3  , Zet.rapidity(), theWeight);
  //theHistograms.fill("GenZ_eta",      "Z's #eta",       100, -7  ,   7  , Zet.eta()     , theWeight);
  //theHistograms.fill("GenZ_phi",      "Z's #phi",        50, -4  ,   4  , Zet.phi()     , theWeight);
  //theHistograms.fill("GenZ_deltaEta", "Z's #Delta#eta",  30, -5  ,   5  , ZdeltaEta     , theWeight);
  //theHistograms.fill("GenZ_deltaR",   "Z's #DeltaR",     20, -0.5,   5.5, ZdeltaR       , theWeight);
  //theHistograms.fill("GenZ_deltaPhi", "Z's #Delta#phi",  50, -4  ,   4  , ZdeltaPhi     , theWeight);
  
  theHistograms.fill("GenZ_massvstrmass", "Z's mass(x) vs trmass(y)", 400, 0, 400, 400, 0, 400, Zet.mass(), Zet.p4().Mt(), theWeight);
  
  // Zl
  //theHistograms.fill("GenZL_ID",   "Zls' ID",   9, 31.5, 40.5, ZlsID     , theWeight);
  //theHistograms.fill("GenZL_size", "Zl's size", 4, -0.5,  3.5, Zls.size(), theWeight);
  
  // WZ
  theHistograms.fill("GenWZ_mass",         "W and Z mass",        268, 160  , 1500  , WZ.mass()   , theWeight);
  theHistograms.fill("GenWZ_trmass",       "W and Z trmass",      268, 160  , 1500  , WZ.p4().Mt(), theWeight);
  theHistograms.fill("GenWZ_pt",           "W and Z p_{t}",       500,   0  , 1000  , WZ.pt()     , theWeight);
  //theHistograms.fill("GenWZ_deltaEta",     "W and Z #Delta#eta",   50,  -9  ,    9  , WZdeltaEta  , theWeight);
  //theHistograms.fill("GenWZ_deltaR",       "W and Z #Delta R",     25,  -0.5,    9  , WZdeltaR    , theWeight);
  //theHistograms.fill("GenWZ_deltaPhi",     "W and Z #Delta#phi",   50,  -4  ,    4  , WZdeltaPhi  , theWeight);
  
  theHistograms.fill("GenWZ_massvstrmass", "WZ's mass(x) vs trmass(y)", 268, 160, 1500, 268, 160, 1500, WZ.mass(), WZ.p4().Mt(), theWeight);
  
  // Jets
  foreach(const Particle jet, genjets){  
    //theHistograms.fill("GenJet_charge", "charge jets",   5, -2.5,   2.5, jet.charge()  , theWeight);  
    //theHistograms.fill("GenJet_mass",   "mass jets",   120,  0  , 120  , jet.mass()    , theWeight);
    //theHistograms.fill("GenJet_pt",     "p_{t} jets",   80,  0  , 400  , jet.pt()      , theWeight);
    //theHistograms.fill("GenJet_Y",      "Y jets",       70, -5  ,   5  , jet.rapidity(), theWeight);
    //theHistograms.fill("GenJet_eta",    "#eta jets",    70, -5  ,   5  , jet.eta()     , theWeight);
    //theHistograms.fill("GenJet_phi",    "#phi jets",    50, -4  ,   4  , jet.phi()     , theWeight);
  }
  
  theHistograms.fill("GenJets_number", "number of all gen jets",  10,  -0.5,  9.5, genjets.size(), theWeight);
  
  theHistograms.fill("GenJJ_mass",     "Leading Jets' mass",       200,  0  , 500, JJp4.M()  , theWeight);
  theHistograms.fill("GenJJ_trmass",   "Leading Jets' trmass",     200,  0  , 500, JJp4.Mt() , theWeight);
  theHistograms.fill("GenJJ_deltaEta", "Leading Jets' #Delta#eta", 100, -9  ,   9, JJdeltaEta, theWeight); 
  theHistograms.fill("GenJJ_deltaR",   "Leading Jets' #DeltaR",     25, -0.5,   9, JJdeltaR  , theWeight);
  theHistograms.fill("GenJJ_deltaPhi", "Leading Jets' #Delta#phi",  50, -4  ,   4, JJdeltaPhi, theWeight);
  
  if(genjets.size() > 2){
    for(int i = 2; i < (int)genjets.size(); i++){
      theHistograms.fill("genJet_pt_M3", "Not-leading jets' p_{t}", 336, 15, 351, genjets[2].pt(), theWeight);
    }
  }
  
  // WZ & leading jets
  theHistograms.fill("GenAll_mass",         "Mass W,Z,J,J",            200, 220, 8220, WZjjp4.M() , theWeight);
  theHistograms.fill("GenAll_trmass",       "Transverse mass W,Z,J,J", 200, 220, 8220, WZjjp4.Mt(), theWeight);
  
  theHistograms.fill("GenAll_massvstrmass", "WZjj's mass(x) vs trmass(y)", 200, 220, 1200, 200, 220, 1200, WZjjp4.M(), WZjjp4.Mt(), theWeight);
  
  
  // ~~~~~~~~~~~~~~~~~~~~~ End of gen Analysis ~~~~~~~~~~~~~~~~~~~~~
}

void WZAnalyzer::RecoAnalysis(DiBosonLepton &WZ, Particle &Jet0, Particle &Jet1){
  // ~~~~~~~~~~~~~~~~~~~~ Begin of reco Analysis ~~~~~~~~~~~~~~~~~~~
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Filters ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // filter on ZLCand's size
  if(ZLCand->size() == 0){
    recoZlempty++;
    return;
  }
  
  // filter on jet's number
  if(jets->size() < 2){
    recoJetless2++;
    return;
  }
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Z & W ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  eventReco++;
  
  BosonLepton recoW;
  BosonLepton recoZ;
  Lepton MET = Lepton(met->p4());
  vector<DiBosonLepton> recoWZs;
  ZLCompositeCandidates recoZls;
  
  // ------------------------ W construction -----------------------
  foreach(const ZLCompositeCandidate Z, *ZLCand){
    recoZls.push_back(Z);
  }
  
  if(recoZls.size() == 1){    
    recoZ = recoZls[0].first;
    recoW = BosonLepton(recoZls[0].second, MET, copysign(24, recoZls[0].second.charge()));
    recoWZs.push_back(DiBosonLepton(recoW, recoZ));
  }
  
  else if(recoZls.size() > 1){
    
    // Z is the one with closer mass to Z
    sort(recoZls.begin(), recoZls.end(), pairMassComparator(0, ZMASS));
    
    // W is the one with the lepton with higher pt
    bool choosingW = kTRUE;
    
    for(int i = 0; choosingW; i++){
      recoWZs.push_back(DiBosonLepton(BosonLepton(recoZls[i].second, MET, copysign(24, recoZls[i].second.charge())), recoZls[i].first));
      
      if(i == (int)recoZls.size() || recoZls[i].first.mass() != recoZls[i + 1].first.mass()){
	choosingW = kFALSE;
      }
    }
    
    sort(recoWZs.begin(), recoWZs.end(), WZPtComparator());
    
    recoW = recoWZs[0].first();
    recoZ = recoWZs[0].second();
  }
  
  WZ = recoWZs[0];
  
  // ------------------------ W & Z variables ----------------------
  // Z
  //double recoZdeltaEta = recoZ.daughter(0).eta() - recoZ.daughter(1).eta();
  //double recoZdeltaPhi = physmath::deltaPhi(recoZ.daughter(0).phi(), recoZ.daughter(1).phi());
  //double recoZdeltaR = abs(physmath::deltaR(recoZ.daughter(0), recoZ.daughter(1)));
  
  // W
  //double recoWdeltaEta = recoW.daughter(0).eta() - recoW.daughter(1).eta();
  //double recoWdeltaPhi = physmath::deltaPhi(recoW.daughter(0).phi(), recoW.daughter(1).phi());
  //double recoWdeltaR = abs(physmath::deltaR(recoW.daughter(0), recoW.daughter(1)));
  
  // WZ
  double recoWZdeltaPhi = physmath::deltaPhi(recoW.phi(), recoZ.phi());
  double recoWZdeltaEta = recoW.eta() - recoZ.eta();
  double recoWZdeltaR = abs(physmath::deltaR(recoW, recoZ));
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Jets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  vector<Particle> recoJets;
  
  foreach(const Particle jet, *jets){
    recoJets.push_back(jet);
  }
  
  // leading jets are those with higher pt
  sort(recoJets.begin(), recoJets.end(), PtComparator());
  Jet0 = recoJets[0];
  Jet1 = recoJets[1];
  
  // leading jets
  TLorentzVector recoJJptot = recoJets[0].p4() + recoJets[1].p4();
  double recoJJdeltaEta = recoJets[0].eta() - recoJets[1].eta();
  double recoJJdeltaPhi = physmath::deltaPhi(recoJets[0].phi(), recoJets[1].phi());
  double recoJJdeltaR = abs(physmath::deltaR(recoJets[0], recoJets[1]));
  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ All ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TLorentzVector recoPtot = recoW.p4() + recoZ.p4() + recoJets[0].p4() + recoJets[1].p4();

  double WZjptDeltaPhi = physmath::deltaPhi(recoWZs[0].phi(), recoJets[0].phi());
  
  double WZj0DeltaR = physmath::deltaR(recoWZs[0], recoJets[0]);
  double WZj1DeltaR = physmath::deltaR(recoWZs[0], recoJets[1]);
  
  if(WZj0DeltaR < WZj1DeltaR){
    double WZj0DeltaPhi = physmath::deltaPhi(recoWZs[0].phi(), recoJets[0].phi());
    theHistograms.fill("recoAll_WZj0vDeltaPhi", "#Delta#phi between recoWZ and recoJets[0] closer than rJ1", 50, -4, 4, WZj0DeltaPhi, theWeight);
    theHistograms.fill("recoAll_WZjvDeltaPhi", "#Delta#phi between recoWZ and recoJets' jet closer",  50, -4, 4, WZj0DeltaPhi, theWeight);
  }
  else{
    double WZj1DeltaPhi = physmath::deltaPhi(recoWZs[0].phi(), recoJets[1].phi());
    theHistograms.fill("recoAll_WZj1vDeltaPhi", "#Delta#phi between recoWZ and recoJets[1] closer than rJ0", 50, -4, 4, WZj1DeltaPhi, theWeight);
    theHistograms.fill("recoAll_WZjvDeltaPhi", "#Delta#phi between recoWZ and recoJets' jet closer",  50, -4, 4, WZj1DeltaPhi, theWeight);    
  }
  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~ Histograms ~~~~~~~~~~~~~~~~~~~~~~~~~
  // MET
  theHistograms.fill("recoMET_pt",  "MET's p_{t}", 325,  0, 650, met->pt() , theWeight);
  //theHistograms.fill("recoMET_phi", "MET's #phi",   50, -4,   4, met->phi(), theWeight);
  
  // W
 // theHistograms.fill("recoW_charge",   "recoW's charge",       5, -2.5,   2.5, recoW.charge()  , theWeight);
  theHistograms.fill("recoW_trmass",   "recoW's trmass",     400,  0  , 400  , recoW.p4().Mt() , theWeight);
  theHistograms.fill("recoW_pt",       "recoW's p_{t}",      260,  0  , 650  , recoW.pt()      , theWeight);
  //theHistograms.fill("recoW_Y",        "recoW's Y",           45, -3  ,   3  , recoW.rapidity(), theWeight);
  //theHistograms.fill("recoW_eta",      "recoW's #eta",       100, -7  ,   7  , recoW.eta()     , theWeight);
  //theHistograms.fill("recoW_phi",      "recoW's #phi",        50, -4  ,   4  , recoW.phi()     , theWeight);
  //theHistograms.fill("recoW_deltaEta", "recoW's #Delta#eta",  55, -5  ,   5  , recoWdeltaEta   , theWeight);
  //theHistograms.fill("recoW_deltaR",   "recoW's #DeltaR",     20, -0.5,   5.5, recoWdeltaR     , theWeight);
  //theHistograms.fill("recoW_deltaPhi", "recoW's #Delta#phi",  50, -4  ,   4  , recoWdeltaPhi   , theWeight);
  
  // Z
  //theHistograms.fill("recoZ_charge",   "recoZ's charge",       5, -2.5,   2.5, recoZ.charge()  , theWeight); //just to be sure...
  theHistograms.fill("recoZ_mass",     "recoZ's mass",       400,  0  , 400  , recoZ.mass()    , theWeight);
  theHistograms.fill("recoZ_trmass",   "recoZ's trmass",     400,  0  , 400  , recoZ.p4().Mt() , theWeight);
  //theHistograms.fill("recoZ_pt",       "recoZ's p_{t}",      260,  0  , 650  , recoZ.pt()      , theWeight);
  //theHistograms.fill("recoZ_Y",        "recoZ's Y",           45, -3  ,   3  , recoZ.rapidity(), theWeight);
  //theHistograms.fill("recoZ_eta",      "recoZ's #eta",       100, -7  ,   7  , recoZ.eta()     , theWeight);
  //theHistograms.fill("recoZ_phi",      "recoZ's #phi",        50, -4  ,   4  , recoZ.phi()     , theWeight);
  //theHistograms.fill("recoZ_deltaEta", "recoZ's #Delta#eta",  55, -5  ,   5  , recoZdeltaEta   , theWeight);
  //theHistograms.fill("recoZ_deltaR",   "recoZ's #DeltaR",     20, -0.5,   5.5, recoZdeltaR     , theWeight);
  //theHistograms.fill("recoZ_deltaPhi", "recoZ's #Delta#phi",  50, -4  ,   4  , recoZdeltaPhi   , theWeight);
  
  // Zl
  theHistograms.fill("recoZl_size", "Reco Zl's size", 15, -0.5, 14.5, recoZls.size(), theWeight);  
  
  // WZ
  theHistograms.fill("recoWZ_trmass",   "recoWZ trmass",        300,   0  ,  900  , recoWZs[0].p4().Mt(), theWeight);
  theHistograms.fill("recoWZ_pt",       "recoWZ p_{T}",         200, 200  , 3000  , recoWZs[0].pt()     , theWeight);
  theHistograms.fill("recoWZ_deltaPhi", "recoW & Z #Delta#phi",  50,  -4  ,    4  , recoWZdeltaPhi      , theWeight);
  theHistograms.fill("recoWZ_deltaEta", "recoW & Z #Delta#eta",  55,  -5  ,    5  , recoWZdeltaEta      , theWeight);
  theHistograms.fill("recoWZ_deltaR",   "recoW & Z #DeltaR",     20,  -0.5,    5.5, recoWZdeltaR        , theWeight);
  
  // Jets
  foreach(const Particle jet, recoJets){
    theHistograms.fill("recoJet_charge", "Jets' charge",  19, -8.5,  10.5, jet.charge()  , theWeight);
    theHistograms.fill("recoJet_pt",     "Jets' p_{t}",  160,  0  , 400  , jet.pt()      , theWeight);
    theHistograms.fill("recoJet_eta",    "Jets' #eta",    70, -5  ,   5  , jet.eta()     , theWeight);
    theHistograms.fill("recoJet_Y",      "Jets' Y",       80, -5  ,   5  , jet.rapidity(), theWeight);
    theHistograms.fill("recoJet_phi",    "Jets' #phi",    50, -4  ,   4  , jet.phi()     , theWeight);
  }
  
  theHistograms.fill("recoJet_number", "Jets' number", 10, -0.5, 9.5, recoJets.size(), theWeight);  
  
  theHistograms.fill("recoJJ_pt",       "Jets p_{tot}",         600,  0  , 2000, recoJJptot.Pt(), theWeight);
  theHistograms.fill("recoJJ_mass",     "Jets mass",            600,  0  , 4500, recoJJptot.M() , theWeight);
  theHistograms.fill("recoJJ_trmass",   "Jets transverse mass", 600,  0  , 4500, recoJJptot.Mt(), theWeight);
  theHistograms.fill("recoJJ_deltaEta", "Jets #Delta#eta",      100, -9  ,    9, recoJJdeltaEta , theWeight); 
  theHistograms.fill("recoJJ_deltaR",   "Jets #DeltaR",          25, -0.5,    9, recoJJdeltaR   , theWeight);
  theHistograms.fill("recoJJ_deltaPhi", "Jets #Delta#phi",       50, -4  ,    4, recoJJdeltaPhi , theWeight);
  
  if(recoJets.size() > 2){
    for(int i = 2; i < (int)recoJets.size(); i++){
      theHistograms.fill("recoJet_pt_M3", "Not-leading jets' p_{t}", 336, 15, 351, recoJets[2].pt(), theWeight);
    }
  }
  
  // WZ & leading jets
  theHistograms.fill("recoAll_trmass",        "Transverse mass recoW,Z,J,J",               200, 220, 8220, recoPtot.Mt(), theWeight);
  theHistograms.fill("recoAll_WZjptDeltaPhi", "#Delta#phi between recoWZ and recoJets[0]",  50,  -4,    4, WZjptDeltaPhi, theWeight);
  
  
  // ~~~~~~~~~~~~~~~~~~~~~ End of reco Analysis ~~~~~~~~~~~~~~~~~~~~
}

void WZAnalyzer::GenRecoAnalysis(const ZZtype genWZ, const Particle genJet0, const Particle genJet1, const DiBosonLepton recoWZ, const Particle recoJet0, const Particle recoJet1){
  eventGenReco++;
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~ Reco vs Gen ~~~~~~~~~~~~~~~~~~~~~~~~~
  bool choosingW = kTRUE;
  vector<DiBosonLepton> recoWZs;
  ZLCompositeCandidates recoZls;
 
  // test if recoW is genW by trmass
  foreach(const ZLCompositeCandidate Z, *ZLCand){
    recoZls.push_back(Z);
  }
  
  sort(recoZls.begin(), recoZls.end(), pairMassComparator(0, ZMASS));
    
  for(int i = 0; choosingW; i++){
    recoWZs.push_back(DiBosonLepton(BosonLepton(recoZls[i].second, Lepton(met->p4()), copysign(24, recoZls[i].second.charge())), recoZls[i].first));
    
    if(i == (int)recoZls.size() || recoZls[i].first.mass() != recoZls[i + 1].first.mass()){
      choosingW = kFALSE;
    }
  }
  
  sort(recoWZs.begin(), recoWZs.end(), pairTrmassComparator(0, genWZ.first().p4().Mt()));

  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~ Histograms ~~~~~~~~~~~~~~~~~~~~~~~~~
  // genW vs recoW
  theHistograms.fill("GR_W_Deltatrmass",  "W's #Deltam_{t} in Genn&Reco events",  600, -300, 300, genWZ.first().p4().Mt() - recoWZ.first().p4().Mt()  , theWeight);
  
  theHistograms.fill("GR_recoWtm_trmass", "recoW by tr mass", 600, 0, 600, recoWZs[0].first().p4().Mt(), theWeight);
  theHistograms.fill("GR_recoWtm_WvsW",  "recoW and recoWtm are the same?", 4, -1.5, 2.5, isTheSame(recoWZ.first(), recoWZs[0].first()));

  // genZ vs recoZ
  theHistograms.fill("GR_Z_Deltatrmass",  "Z's #Deltam_{t} in Genn&Reco events",  600, -300, 300, genWZ.second().p4().Mt() - recoWZ.second().p4().Mt(), theWeight);

  // genWZ vs recoWZ
  theHistograms.fill("GR_WZ_Deltatrmass", "WZ's #Deltam_{t} in Genn&Reco events", 600, -300, 300, genWZ.p4().Mt() - recoWZ.p4().Mt()                  , theWeight);
  
}

void WZAnalyzer::analyze(){
  
  eventSample++;
  
  //gen variables	
  ZZtype genWZ;
  Particle genJet0;
  Particle genJet1;
  
  //reco variables
  DiBosonLepton recoWZ;
  Particle recoJet0;
  Particle recoJet1;
  
  //Gen analysis
  WZAnalyzer::GenAnalysis(genWZ, genJet0, genJet1);
  
  //Reco analysis
  WZAnalyzer::RecoAnalysis(recoWZ, recoJet0, recoJet1);
  
  if(genWZ.pt() != 0. && recoWZ.pt() != 0.){
    //Reco vs Gen analysis
    WZAnalyzer::GenRecoAnalysis(genWZ, genJet0, genJet1, recoWZ, recoJet0, recoJet1);
  }
}

void WZAnalyzer::end(TFile &){
  
  cout << "\n--------------------------------------------------------------------------" << endl;
  
  cout << "\nEvents of the sample analyzed:                       " << setw(9) << eventSample << endl;
  cout << "Gen events analyzed:                                 " << setw(9) << eventGen << endl;
  cout << "Reco events analyzed:                                " << setw(9) << eventReco << endl;
  cout << "Gen&Reco events analyzed:                            " << setw(9) << eventGenReco << endl;
  cout << "RecoZls empty:                                       " << setw(9) << recoZlempty << endl;
  cout << "RecoJets.size < 2:                                   " << setw(9) << recoJetless2 << endl;
  
  // /*
  cout << "\nNumber of                                            " << setw(9) << counter1 << endl;
  cout << "Number of                                            " << setw(9) << counter2 << endl;
  cout << "Number of                                            " << setw(9) << counter3 << endl;
  cout << "Number of                                            " << setw(9) << counter4 << endl;
  // */
  
  // execution time
  endtime = ((float)clock())/CLOCKS_PER_SEC;
  WZAnalyzer::printTime(begintime, endtime);
  cout << "\n--------------------------------------------------------------------------" << endl;
}
