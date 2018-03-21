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
  recoAfterCut = 0;

  gen3e = 0;
  gen3m = 0;
  gen2e1m = 0;
  gen2m1e = 0;
  reco3e = 0;
  reco3m = 0;
  reco2e1m = 0;
  reco2m1e = 0;
  
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
  if(electron.size() + muon.size() + neutrino.size() != 4 && neutrino.size() == 1){
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
  
  // Z and W must be at least on shell  
  TLorentzVector Ptot = neutrino[0].p4();
  
  foreach(const Particle lep, lepton)
    Ptot += lep.p4();
  
  theHistograms.fill("AllGenlllnu_mass",   "m 3 leptons and #nu",     150, 0, 1500, Ptot.M() , theWeight);
  theHistograms.fill("AllGenlllnu_trmass", "m_{T} 3 leptons and #nu", 150, 0, 1500, Ptot.Mt(), theWeight);
  
  if(Ptot.M() < 165){
    return;
  }

  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Z & W ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Vtype Weh;
  Vtype Zet;
  vector<Zltype> Zls;
  
  // ---------------------- W & Z construction ---------------------
  // case 2e 1m
  if(electron.size()==2 && muon.size()==1 && electron[0].charge() != electron[1].charge()){
    gen2e1m++;
    Zet = Vtype(electron[0], electron[1], 23);
    Zls.push_back(Zltype(Zet, muon[0]));
    Weh = Vtype(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()) );
  }
  
  // case 1e 2m
  if(electron.size()==1 && muon.size()==2 && muon[0].charge() != muon[1].charge()){
    gen2m1e++;
    Zet = Vtype(muon[0], muon[1], 23);
    Zls.push_back(Zltype(Zet, electron[0]));
    Weh = Vtype(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()) );
  }
  
  // case 3e
  if(electron.size()==3){
    gen3e++;
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
    
    // Z is made up of the couple which gives a better Zmass 
    sort(Zls.begin(), Zls.end(), pairMassComparator(0, ZMASS));
    Zet = Zls[0].first;
    
    // W is made up of the remaining lepton and the neutrino
    Weh = Vtype(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()));
  }
  
  // case 3m
  else if(muon.size()==3){
    gen3m++;
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
    
    sort(Zls.begin(), Zls.end(), pairMassComparator(0, ZMASS));
    Zet = Zls[0].first;
    Weh = Vtype(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()));
  }
  
  if(Zls.size() < 1){
    return;
  }
  
  // W&Z diboson
  WZ = ZZtype(Weh, Zet);
  eventGen++;
  
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
  //double JJdeltaEta = genjets[0].eta() - genjets[1].eta();
  //double JJdeltaPhi = physmath::deltaPhi(genjets[0].phi(), genjets[1].phi());
  //double JJdeltaR = abs(physmath::deltaR(genjets[0], genjets[1]));
  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ All ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TLorentzVector WZjjp4 = WZ.p4() + genjets[0].p4() + genjets[1].p4();

  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~ Histograms ~~~~~~~~~~~~~~~~~~~~~~~~~
  /*
  // Nu
  theHistograms.fill("GenN_charge", "#nu's charge",  5, -2.5,   2.5, neutrino[0].charge()  , theWeight); //just to be sure...
  theHistograms.fill("GenN_trmass", "W's trmass",  400,  0  , 400  , neutrino[0].p4().Mt() , theWeight);
  theHistograms.fill("GenN_pt",     "#nu's p_{t}", 140,  0  , 700  , neutrino[0].pt()      , theWeight);
  theHistograms.fill("GenN_Y",      "#nu's Y",      90, -6.5,   6.5, neutrino[0].rapidity(), theWeight);
  theHistograms.fill("GenN_eta",    "#nu's #eta",   90, -6.5,   6.5, neutrino[0].eta()     , theWeight);
  theHistograms.fill("GenN_phi",    "#nu's #phi",   50, -3.5,   3.5, neutrino[0].phi()     , theWeight);
  
  // W
  theHistograms.fill("GenW_charge",       "W's charge",       5, -2.5,   2.5, Weh.charge()  , theWeight);
  theHistograms.fill("GenW_mass",         "W's mass",       400,  0  , 400  , Weh.mass()    , theWeight);
  theHistograms.fill("GenW_trmass",       "W's trmass",     400,  0  , 400  , Weh.p4().Mt() , theWeight);
  theHistograms.fill("GenW_pt",           "W's p_{t}",      130,  0  , 650  , Weh.pt()      , theWeight);
  theHistograms.fill("GenW_Y",            "W's Y",           50, -5  ,   5  , Weh.rapidity(), theWeight);
  theHistograms.fill("GenW_eta",          "W's #eta",        50, -9  ,   9  , Weh.eta()     , theWeight);
  theHistograms.fill("GenW_phi",          "W's #phi",        50, -3.5,   3.5, Weh.phi()     , theWeight);
  theHistograms.fill("GenW_deltaEta",     "W's #Delta#eta",  50, -6.5,   6.5, WdeltaEta     , theWeight);
  theHistograms.fill("GenW_deltaR",       "W's #DeltaR",     50, -0.5,   6.5, WdeltaR       , theWeight);
  theHistograms.fill("GenW_deltaPhi",     "W's #Delta#phi",  50, -3.5,   3.5, WdeltaPhi     , theWeight);
  
  theHistograms.fill("GenW_massvstrmass", "W's mass(x) vs trmass(y)", 400, 0, 400, 400, 0, 400, Weh.mass(), Weh.p4().Mt(), theWeight);
  
  // Z
  theHistograms.fill("GenZ_charge",   "Z's charge",       5, -2.5,   2.5, Zet.charge()  , theWeight); //just to be sure...
  theHistograms.fill("GenZ_mass",     "Z's mass",       400,  0  , 400  , Zet.mass()    , theWeight);
  theHistograms.fill("GenZ_trmass",   "Z's trmass",     400,  0  , 400  , Zet.p4().Mt() , theWeight);
  theHistograms.fill("GenZ_pt",       "Z's p_{t}",      130,  0  , 650  , Zet.pt()      , theWeight);
  theHistograms.fill("GenZ_Y",        "Z's Y",           45, -3  ,   3  , Zet.rapidity(), theWeight);
  theHistograms.fill("GenZ_eta",      "Z's #eta",       100, -7  ,   7  , Zet.eta()     , theWeight);
  theHistograms.fill("GenZ_phi",      "Z's #phi",        50, -3.5,   3.5, Zet.phi()     , theWeight);
  theHistograms.fill("GenZ_deltaEta", "Z's #Delta#eta",  30, -5  ,   5  , ZdeltaEta     , theWeight);
  theHistograms.fill("GenZ_deltaR",   "Z's #DeltaR",     20, -0.5,   5.5, ZdeltaR       , theWeight);
  theHistograms.fill("GenZ_deltaPhi", "Z's #Delta#phi",  50, -3.5,   3.5, ZdeltaPhi     , theWeight);
  
  theHistograms.fill("GenZ_massvstrmass", "Z's mass(x) vs trmass(y)", 400, 0, 400, 400, 0, 400, Zet.mass(), Zet.p4().Mt(), theWeight);
  
  // Zl
  theHistograms.fill("GenZL_ID",   "Zls' ID",   9, 31.5, 40.5, ZlsID     , theWeight);
  theHistograms.fill("GenZL_size", "Zl's size", 4, -0.5,  3.5, Zls.size(), theWeight);
  
  // WZ
  theHistograms.fill("GenWZ_mass",         "W and Z mass",        268, 160  , 1500  , WZ.mass()   , theWeight);
  theHistograms.fill("GenWZ_trmass",       "W and Z trmass",      268, 160  , 1500  , WZ.p4().Mt(), theWeight);
  theHistograms.fill("GenWZ_pt",           "W and Z p_{t}",       500,   0  , 1000  , WZ.pt()     , theWeight);
  theHistograms.fill("GenWZ_deltaEta",     "W and Z #Delta#eta",   50,  -9  ,    9  , WZdeltaEta  , theWeight);
  theHistograms.fill("GenWZ_deltaR",       "W and Z #Delta R",     25,  -0.5,    9  , WZdeltaR    , theWeight);
  theHistograms.fill("GenWZ_deltaPhi",     "W and Z #Delta#phi",   50,  -3.5,    3.5, WZdeltaPhi  , theWeight);
  
  theHistograms.fill("GenWZ_massvstrmass", "WZ's mass(x) vs trmass(y)", 268, 160, 1500, 268, 160, 1500, WZ.mass(), WZ.p4().Mt(), theWeight);
  */
  // Jets
  /*
  foreach(const Particle jet, genjets){  
    theHistograms.fill("GenJet_charge", "charge jets",   5, -2.5,   2.5, jet.charge()  , theWeight);  
    theHistograms.fill("GenJet_mass",   "mass jets",   120,  0  , 120  , jet.mass()    , theWeight);
    theHistograms.fill("GenJet_pt",     "p_{t} jets",   80,  0  , 400  , jet.pt()      , theWeight);
    theHistograms.fill("GenJet_Y",      "Y jets",       70, -5  ,   5  , jet.rapidity(), theWeight);
    theHistograms.fill("GenJet_eta",    "#eta jets",    70, -5  ,   5  , jet.eta()     , theWeight);
    theHistograms.fill("GenJet_phi",    "#phi jets",    50, -3.5,   3.5, jet.phi()     , theWeight);
  }
  
  theHistograms.fill("GenJets_number", "number of all gen jets",  10,  -0.5,  9.5, genjets.size(), theWeight);
  */
  /*
  theHistograms.fill("GenJJ_mass",     "Leading Jets' mass",       600,  0  , 4500  , JJp4.M()  , theWeight);
  theHistograms.fill("GenJJ_trmass",   "Leading Jets' trmass",     600,  0  , 4500  , JJp4.Mt() , theWeight);
  theHistograms.fill("GenJJ_deltaEta", "Leading Jets' #Delta#eta", 100, -9  ,    9  , JJdeltaEta, theWeight); 
  theHistograms.fill("GenJJ_deltaR",   "Leading Jets' #DeltaR",     25, -0.5,    9  , JJdeltaR  , theWeight);
  theHistograms.fill("GenJJ_deltaPhi", "Leading Jets' #Delta#phi",  50, -3.5,    3.5, JJdeltaPhi, theWeight);
  */
  /*
  if(genjets.size() > 2){
    for(int i = 2; i < (int)genjets.size(); i++){
      theHistograms.fill("GenJet_pt_M3", "Not-leading jets' p_{t}", 64, 30, 350, genjets[2].pt(), theWeight);
    }
  }
  */
  
  // WZ & leading jets
  /*
  theHistograms.fill("GenAll_mass",         "Mass W,Z,J,J",            200, 220, 8220, WZjjp4.M() , theWeight);
  theHistograms.fill("GenAll_trmass",       "Transverse mass W,Z,J,J", 200, 220, 8220, WZjjp4.Mt(), theWeight);
  
  theHistograms.fill("GenAll_massvstrmass", "WZjj's mass(x) vs trmass(y)", 200, 220, 1200, 200, 220, 1200, WZjjp4.M(), WZjjp4.Mt(), theWeight);
  */

  
  // ~~~~~~~~~~~~~~~~~~~~~ End of gen Analysis ~~~~~~~~~~~~~~~~~~~~~
}

void WZAnalyzer::RecoAnalysis(DiBosonLepton &WZ, Particle &Jet0, Particle &Jet1){
  // ~~~~~~~~~~~~~~~~~~~~ Begin of reco Analysis ~~~~~~~~~~~~~~~~~~~
  int cut = 0;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);

  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Filters ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ZLCompositeCandidates recoZls;
  ZLCompositeCandidates tempZls;
  
  // filter on ZLCand's size > 0
  theHistograms.fill("recoZl_size_ZLCand", 12, -0.5, 11.5, ZLCand->size());

  if(ZLCand->size() == 0){
    recoZlempty++;
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
    
  // filter on 3rd lepton full selection
  foreach(const ZLCompositeCandidate Zl, *ZLCand){
    if(Zl.second.passFullSel()){
      tempZls.push_back(Zl);
    }
  }

  if(tempZls.size() == 0){
    recoZlempty++;
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);

  // filter on Z's mass 60 < m < 120
  foreach(const ZLCompositeCandidate Zl, tempZls){
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
  
  // filter on jet's number > 2
  if(jets->size() < 2){
    recoJetless2++;
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);

  // filter on MET's pt
  if(met->pt() < 30.){
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);

  
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

    // filter on W's trmass 30 < trm < 500
    if(recoW.p4().Mt() > 30. && recoW.p4().Mt() < 500.){
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

      // filter on W's trmass 30 < trm < 500
      if(recoWtemp.p4().Mt() > 30. && recoWtemp.p4().Mt() < 500.){
	recoWZs.push_back(DiBosonLepton(recoWtemp, recoZls[i].first));
      }
      
      if(i == (int)recoZls.size() - 1 || recoZls[i].first.mass() - recoZls[i + 1].first.mass() < 0.1){
	choosingW = kFALSE;
      }
    }

    if(recoWZs.size() == 0){
      return;
    }
    
    sort(recoWZs.begin(), recoWZs.end(), WZPtComparator());
    
    recoW = recoWZs[0].first();
    recoZ = recoWZs[0].second();
  }
  
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  
  WZ = recoWZs[0];
  eventReco++;
  
  
  // ------------------------ W & Z variables ----------------------
  // Z
  double recoZdeltaEta = recoZ.daughter(0).eta() - recoZ.daughter(1).eta();
  double recoZdeltaPhi = physmath::deltaPhi(recoZ.daughter(0).phi(), recoZ.daughter(1).phi());
  double recoZdeltaR = abs(physmath::deltaR(recoZ.daughter(0), recoZ.daughter(1)));
  
  float recoZpzlp;
  float recoZpzlm;
  if(recoZ.daughter(0).id() > 0){
    recoZpzlp = recoZ.daughter(0).p4().Pz();
    recoZpzlm = recoZ.daughter(1).p4().Pz();
  }
  else{
    recoZpzlm = recoZ.daughter(0).p4().Pz();
    recoZpzlp = recoZ.daughter(1).p4().Pz();
  }
  
  // W
  //double recoWdeltaEta = recoW.daughter(0).eta() - recoW.daughter(1).eta();
  double recoWdeltaPhi = physmath::deltaPhi(recoW.daughter(0).phi(), recoW.daughter(1).phi());
  //double recoWdeltaR = abs(physmath::deltaR(recoW.daughter(0), recoW.daughter(1)));

  // Zl
  int recoZlID = abs(recoZ.daughter(0).id()) + abs(recoZ.daughter(1).id()) + abs(recoW.daughter(0).id());
  TLorentzVector recoZlp4 = recoZ.daughter(0).p4() + recoZ.daughter(1).p4() + recoW.daughter(0).p4();
  
  if(recoZlID == 33)      reco3e++;
  else if(recoZlID == 35) reco2e1m++;
  else if(recoZlID == 37) reco2m1e++;
  else if(recoZlID == 39) reco3m++;
    
    
  // WZ
  double recoWZdeltaPhi = physmath::deltaPhi(recoW.phi(), recoZ.phi());
  //double recoWZdeltaEta = recoW.eta() - recoZ.eta();
  //double recoWZdeltaR = abs(physmath::deltaR(recoW, recoZ));
  
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

  //double WZj0DeltaR = physmath::deltaR(recoWZs[0], recoJets[0]);
  //double WZj1DeltaR = physmath::deltaR(recoWZs[0], recoJets[1]);
  //double WZj0DeltaPhi = physmath::deltaPhi(recoWZs[0].phi(), recoJets[0].phi());
  //double WZj1DeltaPhi = physmath::deltaPhi(recoWZs[0].phi(), recoJets[1].phi());
  
  float zeppenfeldJJZ = recoZ.eta() - (recoJets[0].eta() + recoJets[1].eta())/2;
  float zeppenfeldJJl = recoW.daughter(0).eta() - (recoJets[0].eta() + recoJets[1].eta())/2;
  float zeppenfeldllJ0 = recoJets[0].eta() - (recoZ.daughter(0).eta() + recoZ.daughter(1).eta())/2;
  float zeppenfeldllJ1 = recoJets[1].eta() - (recoZ.daughter(0).eta() + recoZ.daughter(1).eta())/2;
  float zeppenfeldlll = recoW.daughter(0).eta() - (recoZ.daughter(0).eta() + recoZ.daughter(1).eta())/2;

  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Cuts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //histos before cuts
  theHistograms.fill("recoJJ_mass_BC",         "Jets' mass",                         300, 10, 4510, recoJJptot.M()     , theWeight);
  theHistograms.fill("recoJJ_deltaEta_BC",     "Jets' #Delta#eta",                    50, -9,    9, recoJJdeltaEta     , theWeight); 
  theHistograms.fill("recoJJ_deltaEta_abs_BC", "Jets' #Delta#eta",                    25,  0,    9, abs(recoJJdeltaEta), theWeight);
  theHistograms.fill("Zeppenfeld_llJ0_BC",     "Zeppenfeld variable for Z and Jet0", 100, -6,    6, zeppenfeldllJ0     , theWeight);
  theHistograms.fill("Zeppenfeld_llJ0_abs_BC", "Zeppenfeld variable for Z and Jet0", 100,  0,    6, abs(zeppenfeldllJ0), theWeight);

  ///*
  if(recoJJptot.M() < 285.5){ //397.773; filter on JJ's mass
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  //*/
  
  ///*
  if(abs(recoJJdeltaEta) < 2.6){ //2.592; filter on JJ's deltaEta 2.106
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  //*/

  ///*
  if(abs(zeppenfeldllJ0) < 0.6){ //0.72; filter on Zeppenfeld variable
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  //*/

  theHistograms.fill("recoJJ_mass_BJ1pt",     "Jets' mass",        70, 280, 4510, recoJJptot.M(), theWeight);
  theHistograms.fill("recoJJ_deltaEta_BJ1pt", "Jets' #Delta#eta",  50,  -9,    9, recoJJdeltaEta, theWeight);
  
  ///*
  if(recoJets[1].pt() < 50){ //0 for QCD; filter on second jet's pt
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  //*/

  theHistograms.fill("recoJJ_mass_B3lpt",     "Jets' mass",        70, 280, 4510, recoJJptot.M(), theWeight);
  theHistograms.fill("recoJJ_deltaEta_B3lpt", "Jets' #Delta#eta",  50,  -9,    9, recoJJdeltaEta, theWeight);

  ///*
  if(recoW.daughter(0).pt() < 30){ //10; filter on 3rd lepton pt
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  //*/

  theHistograms.fill("recoJJ_mass_BZm",     "Jets' mass",        70, 280, 4510, recoJJptot.M(), theWeight);
  theHistograms.fill("recoJJ_deltaEta_BZm", "Jets' #Delta#eta",  50,  -9,    9, recoJJdeltaEta, theWeight);

  ///*
  if(abs(ZMASS - recoZ.mass()) > 15){ //25; filter on Zmass
    return;
  }
  cut++;
  theHistograms.fill("RecoCuts", "Events after cuts", 13, -0.5, 12.5, cut);
  //*/

  recoAfterCut++;
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~ Histograms ~~~~~~~~~~~~~~~~~~~~~~~~~
  // ------------------------------ MET ----------------------------
  theHistograms.fill("recoMET_trmass", "MET's trmass",  100, 25  , 425  , met->p4().Mt(), theWeight);
  theHistograms.fill("recoMET_pt",     "MET's p_{t}",   140,  0  , 700  , MET.pt()      , theWeight);
  //theHistograms.fill("recoMET_phi",    "MET's #phi",     50, -3.5,   3.5, MET.phi()     , theWeight);
  
  // ------------------------------- W -----------------------------
  //theHistograms.fill("recoW_charge",   "recoW's charge",       5, -2.5,   2.5, recoW.charge()  , theWeight);
  theHistograms.fill("recoW_trmass",   "recoW's trmass",      80, 25  , 505  , recoW.p4().Mt() , theWeight);
  theHistograms.fill("recoW_pt",       "recoW's p_{t}",      260,  0  , 650  , recoW.pt()      , theWeight);
  //theHistograms.fill("recoW_Y",        "recoW's Y",           45, -3  ,   3  , recoW.rapidity(), theWeight);
  //theHistograms.fill("recoW_eta",      "recoW's #eta",       100, -7  ,   7  , recoW.eta()     , theWeight);
  //theHistograms.fill("recoW_phi",      "recoW's #phi",        50, -3.5,   3.5, recoW.phi()     , theWeight);
  //theHistograms.fill("recoW_deltaEta", "recoW's #Delta#eta",  55, -5  ,   5  , recoWdeltaEta   , theWeight);
  //theHistograms.fill("recoW_deltaR",   "recoW's #DeltaR",     20, -0.5,   5.5, recoWdeltaR     , theWeight);
  theHistograms.fill("recoW_deltaPhi", "recoW's #Delta#phi",  50, -3.5,   3.5, recoWdeltaPhi   , theWeight);
  
  // ------------------------------- Z -----------------------------
  //theHistograms.fill("recoZ_charge",   "recoZ's charge",       5, -2.5,   2.5, recoZ.charge()  , theWeight); //just to be sure...
  theHistograms.fill("recoZ_mass",     "recoZ's mass",        35, 55  , 125  , recoZ.mass()    , theWeight);
  theHistograms.fill("recoZ_trmass",   "recoZ's trmass",     175, 55  , 405  , recoZ.p4().Mt() , theWeight);
  theHistograms.fill("recoZ_pt",       "recoZ's p_{t}",      130,  0  , 650  , recoZ.pt()      , theWeight);
  //theHistograms.fill("recoZ_Y",        "recoZ's Y",           45, -3  ,   3  , recoZ.rapidity(), theWeight);
  //theHistograms.fill("recoZ_eta",      "recoZ's #eta",       100, -7  ,   7  , recoZ.eta()     , theWeight);
  //theHistograms.fill("recoZ_phi",      "recoZ's #phi",        50, -3.5,   3.5, recoZ.phi()     , theWeight);
  theHistograms.fill("recoZ_deltaEta", "recoZ's #Delta#eta",  55, -5  ,   5  , recoZdeltaEta   , theWeight);
  theHistograms.fill("recoZ_deltaR",   "recoZ's #DeltaR",     20, -0.5,   5.5, recoZdeltaR     , theWeight);
  theHistograms.fill("recoZ_deltaPhi", "recoZ's #Delta#phi",  50, -3.5,   3.5, recoZdeltaPhi   , theWeight);
  
  theHistograms.fill("recoZ_pzp_pzm",  "recoZ's leptons' pz", 320, -400, 400, 320, -400, 400, recoZpzlp, recoZpzlm, theWeight);
  
  // ------------------------------ Zl -----------------------------
  theHistograms.fill("recoZl_size",   "Reco Zl's size",           15, -0.5,   14.5, recoZls.size()        , theWeight);
  theHistograms.fill("recoZl_mass",   "3 leptons mass",          400,  0  , 1200  , recoZlp4.M()          , theWeight);
  theHistograms.fill("recoZl_1st_pt", "Z's 1^{st} lepton p_{t}", 200,  0  ,  400  , recoZ.daughter(0).pt(), theWeight);
  theHistograms.fill("recoZl_2nd_pt", "Z's 2^{nd} lepton p_{t}", 200,  0  ,  400  , recoZ.daughter(1).pt(), theWeight);
  theHistograms.fill("recoZl_3rd_pt", "W's lepton p_{t}",        200,  0  ,  400  , recoW.daughter(0).pt(), theWeight);
  
  // ------------------------------ WZ -----------------------------
  theHistograms.fill("recoWZ_trmass",   "recoWZ trmass",         50,   0  , 2400  , recoWZs[0].p4().Mt(), theWeight);
  theHistograms.fill("recoWZ_pt",       "recoWZ p_{T}",         100, 200  , 2400  , recoWZs[0].pt()     , theWeight);
  theHistograms.fill("recoWZ_deltaPhi", "recoW & Z #Delta#phi",  50,  -3.5,    3.5, recoWZdeltaPhi      , theWeight);
  //theHistograms.fill("recoWZ_deltaEta", "recoW & Z #Delta#eta",  55,  -5  ,    5  , recoWZdeltaEta      , theWeight);
  //theHistograms.fill("recoWZ_deltaR",   "recoW & Z #DeltaR",     20,  -0.5,    5.5, recoWZdeltaR        , theWeight);
  
  theHistograms.fill("recoWZ_ptWvsptZ", "W's p_{t} (x) vs Z's p_{t} (y)", 260,  0,  650, 260,  0,  650, recoW.pt(), recoZ.pt(), theWeight);
  
  // ----------------------------- Jets ----------------------------
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

  // first jet
  theHistograms.fill("recoJ0_pt",  "Jet0's p_{t}", 160,    0, 400, recoJets[0].pt()     , theWeight);
  theHistograms.fill("recoJ0_pz",  "Jet0's p_{z}", 400, -800, 800, recoJets[0].p4().Pz(), theWeight);
  theHistograms.fill("recoJ0_Eta", "Jet0's #eta",  100,   -7,   7, recoJets[0].eta()    , theWeight);

  // second jet
  theHistograms.fill("recoJ1_pt",  "Jet1's p_{t}", 160,    0, 400, recoJets[1].pt()     , theWeight);
  theHistograms.fill("recoJ1_pz",  "Jet1's p_{z}", 400, -800, 800, recoJets[1].p4().Pz(), theWeight);
  theHistograms.fill("recoJ1_Eta", "Jet1's #eta",  100,   -7,   7, recoJets[1].eta()    , theWeight);

  // third jet
  /*
  if(recoJets.size() > 2){
    for(int i = 2; i < (int)recoJets.size(); i++){
      theHistograms.fill("recoJet_pt_M3", "Not-leading jets' p_{t}", 64, 30, 350, recoJets[2].pt(), theWeight);
    }
  }
  */

  // leading jets  
  theHistograms.fill("recoJJ_pt",              "Jets' p_{t}",               200,   0  ,  1200  , recoJJptot.Pt()                                                    , theWeight);
  theHistograms.fill("recoJJ_ptpt",            "Jets' p_{t}*p_{t}}",        800, 500  , 80500  , recoJets[0].pt()*recoJets[1].pt()                                  , theWeight);
  theHistograms.fill("recoJJ_pt0_pt1_sum",     "Jets' p_{t}+p_{t}",         320,   0  ,   800  , recoJets[0].pt()+recoJets[1].pt()                                  , theWeight);
  theHistograms.fill("recoJJ_pt0_pt1_sumquad", "Jets' p_{t}^{2}+p_{t}^{2}", 800,   0  , 80000  , recoJets[0].pt()*recoJets[0].pt()+recoJets[1].pt()*recoJets[1].pt(), theWeight);
  theHistograms.fill("recoJJ_mass",            "Jets' mass",                 42, 280  ,  4480  , recoJJptot.M()                                                     , theWeight);
  theHistograms.fill("recoJJ_massquad",        "Jets' mass^{2}",            600,   0  , 10000  , recoJJptot.M()*recoJJptot.M()                                      , theWeight);
  theHistograms.fill("recoJJ_trmass",          "Jets' transverse mass",     300,   0  ,  4500  , recoJJptot.Mt()                                                    , theWeight);
  theHistograms.fill("recoJJ_deltaEta",        "Jets' #Delta#eta",           50,  -9  ,     9  , recoJJdeltaEta                                                     , theWeight);
  theHistograms.fill("recoJJ_deltaEta_abs",    "Jets' |#Delta#eta|",         25,   0  ,     8.1, abs(recoJJdeltaEta)                                                , theWeight);
  theHistograms.fill("recoJJ_deltaEtaquad",    "Jets' #Delta#eta^{2}",      400, -81  ,    81  , recoJJdeltaEta*recoJJdeltaEta                                      , theWeight);
  theHistograms.fill("recoJJ_deltaR",          "Jets' #DeltaR",              25,  -0.5,     9  , recoJJdeltaR                                                       , theWeight);
  theHistograms.fill("recoJJ_deltaPhi",        "Jets' #Delta#phi",           50,  -3.5,     3.5, recoJJdeltaPhi                                                     , theWeight);
  theHistograms.fill("recoJJ_deltapt",         "Jets' #Deltap_{t}",         600,   0  ,   600  , recoJets[0].pt()-recoJets[1].pt()                                  , theWeight);
  
  theHistograms.fill("recoJJ_pt0_pt1",           "Jet0's p_{t} (x) vs Jet1's p_{t} (y)",        160,    0,   400, 160,    0  ,   400  , recoJets[0].pt()                 , recoJets[1].pt()                 , theWeight);
  theHistograms.fill("recoJJ_pz0_pz1",           "Jet0's p_{z} (x) vs Jet1's p_{z} (y)",        320, -400,   400, 320, -400  ,   400  , recoJets[0].p4().Pz()            , recoJets[1].p4().Pz()            , theWeight);
  theHistograms.fill("recoJJ_J0pt_ptpt",         "Jet0's p_{t} (x) vs Jets' p_{t}^{2} (y)",     160,   20,   420, 800,  500  , 80500  , recoJets[0].pt()                 , recoJets[0].pt()*recoJets[1].pt(), theWeight);
  theHistograms.fill("recoJJ_J1pt_ptpt",         "Jet1's p_{t} (x) vs Jets' p_{t}^{2} (y)",     112,   20,   300, 800,  500  , 80500  , recoJets[1].pt()                 , recoJets[0].pt()*recoJets[1].pt(), theWeight);
  theHistograms.fill("recoJJ_pt_ptpt",           "Jets' p_{t} (x) vs Jets' p_{t}^{2} (y)",      300,    0,  1000, 800,  500  , 80500  , recoJJptot.Pt()                  , recoJets[0].pt()*recoJets[1].pt(), theWeight);
  theHistograms.fill("recoJJ_ptpt_deltaEta",     "Jets' p_{t}^{2} (x) vs Jets' #Delta#eta (y)", 800,  500, 80500, 100,   -9  ,     9  , recoJets[0].pt()*recoJets[1].pt(), recoJJdeltaEta                   , theWeight);
  theHistograms.fill("recoJJ_mass_deltaEta",     "Jets' #Delta#eta (y) vs m_{t} (x)",           300,    0,  1000, 100,   -9  ,     9  , recoJJptot.M()                   , recoJJdeltaEta                   , theWeight);
  theHistograms.fill("recoJJ_pt_deltaEta",       "Jets' #Delta#eta (y) vs p_{t} (x)",           300,    0,  1000, 100,   -9  ,     9  , recoJJptot.Pt()                  , recoJJdeltaEta                   , theWeight);
  theHistograms.fill("recoJJ_J0pt_deltaEta",     "Jets' #Delta#eta (y) vs Jet0's p_{t} (x)",    300,    0,  1000, 100,   -9  ,     9  , recoJets[0].pt()                 , recoJJdeltaEta                   , theWeight);
  theHistograms.fill("recoJJ_J1pt_deltaEta",     "Jets' #Delta#eta (y) vs Jet1's p_{t} (x)",    240,    0,   800, 100,   -9  ,     9  , recoJets[1].pt()                 , recoJJdeltaEta                   , theWeight);
  theHistograms.fill("recoJJ_deltaEta_deltaPhi", "Jets' #Delta#eta (x) vs #Delta#phi (y)",      100,   -9,     9,  50,   -3.5,     3.5, recoJJdeltaEta                   , recoJJdeltaPhi                   , theWeight);

  // deltaphi when the jet is closer by deltaR
  /*
  if(abs(WZj0DeltaR) < abs(WZj1DeltaR)){
    theHistograms.fill("recoAll_WZj0_vicinoDR_DeltaPhi", "#Delta#phi between recoWZ and recoJets[0] closer than rJ1", 50, -3.5, 3.5, WZj0DeltaPhi, theWeight);
    theHistograms.fill("recoAll_WZj_vicinoDR_DeltaPhi",  "#Delta#phi between recoWZ and recoJets' jet closer",        50, -3.5, 3.5, WZj0DeltaPhi, theWeight);

    if(recoJets.size() == 2){
      theHistograms.fill("recoAll_WZj0_vicinoDR_DeltaPhi_2jet", "#Delta#phi between recoWZ and recoJets[0] closer than rJ1", 50, -3.5, 3.5, WZj0DeltaPhi, theWeight);
      theHistograms.fill("recoAll_WZj_vicinoDR_DeltaPhi_2jet",  "#Delta#phi between recoWZ and recoJets[0] closer than rJ1", 50, -3.5, 3.5, WZj0DeltaPhi, theWeight);
    }
  }
  else{
    theHistograms.fill("recoAll_WZj1_vicinoDR_DeltaPhi", "#Delta#phi between recoWZ and recoJets[1] closer than rJ0", 50, -3.5, 3.5, WZj1DeltaPhi, theWeight);
    theHistograms.fill("recoAll_WZj_vicinoDR_DeltaPhi",  "#Delta#phi between recoWZ and recoJets' jet closer",        50, -3.5, 3.5, WZj1DeltaPhi, theWeight);   

    if(recoJets.size() == 2){
      theHistograms.fill("recoAll_WZj1_vicinoDR_DeltaPhi_2jet", "#Delta#phi between recoWZ and recoJets[0] closer than rJ1", 50, -3.5, 3.5, WZj1DeltaPhi, theWeight);
      theHistograms.fill("recoAll_WZj_vicinoDR_DeltaPhi_2jet",  "#Delta#phi between recoWZ and recoJets[0] closer than rJ1", 50, -3.5, 3.5, WZj1DeltaPhi, theWeight);
    } 
  }
  */

  // deltaphi when the jet is closer by deltaPhi
  /*
  if(abs(WZj0DeltaPhi) < abs(WZj1DeltaPhi)){
    theHistograms.fill("recoAll_WZj0_vicinoDp_DeltaPhi", "#Delta#phi between recoWZ and recoJets[0] closer than rJ1", 50, -3.5, 3.5, WZj0DeltaPhi, theWeight);
    theHistograms.fill("recoAll_WZj_vicinoDp_DeltaPhi",  "#Delta#phi between recoWZ and recoJets' jet closer",        50, -3.5, 3.5, WZj0DeltaPhi, theWeight);

    if(recoJets.size() == 2){
      theHistograms.fill("recoAll_WZj0_vicinoDp_DeltaPhi_2jet", "#Delta#phi between recoWZ and recoJets[0] closer than rJ1", 50, -3.5, 3.5, WZj0DeltaPhi, theWeight);
      theHistograms.fill("recoAll_WZj_vicinoDp_DeltaPhi_2jet",  "#Delta#phi between recoWZ and recoJets[0] closer than rJ1", 50, -3.5, 3.5, WZj0DeltaPhi, theWeight);
    }
  }
  else{
    theHistograms.fill("recoAll_WZj1_vicinoDp_DeltaPhi", "#Delta#phi between recoWZ and recoJets[1] closer than rJ0", 50, -3.5, 3.5, WZj1DeltaPhi, theWeight);
    theHistograms.fill("recoAll_WZj_vicinoDp_DeltaPhi",  "#Delta#phi between recoWZ and recoJets' jet closer",        50, -3.5, 3.5, WZj1DeltaPhi, theWeight);   

    if(recoJets.size() == 2){
      theHistograms.fill("recoAll_WZj1_vicinoDp_DeltaPhi_2jet", "#Delta#phi between recoWZ and recoJets[0] closer than rJ1", 50, -3.5, 3.5, WZj1DeltaPhi, theWeight);
      theHistograms.fill("recoAll_WZj_vicinoDp_DeltaPhi_2jet",  "#Delta#phi between recoWZ and recoJets[0] closer than rJ1", 50, -3.5, 3.5, WZj1DeltaPhi, theWeight);
    } 
  }
  */

  // ---------------------- WZ & leading jets ----------------------
  //theHistograms.fill("recoAll_mass",           "Mass recoW,Z,J,J",                          150, 220  , 7720  , recoPtot.M()                                      , theWeight);
  theHistograms.fill("recoAll_trmass",         "Transverse mass recoW,Z,J",                  34, 220  , 7095  , recoPtot.Mt()                                     , theWeight);
  theHistograms.fill("recoAll_WZjpt_deltaPhi", "#Delta#phi between recoWZ and recoJets[0]",  50,  -3.5,    3.5, WZjptDeltaPhi                                     , theWeight);
  theHistograms.fill("recoAll_Zj0_deltaEta",   "#Delta#eta between reco Z and recoJets[0]", 100,  -7  ,    7  , recoZ.eta()-recoJets[0].eta()                     , theWeight);
  theHistograms.fill("recoAll_Zj1_deltaEta",   "#Delta#eta between reco Z and recoJets[1]", 100,  -7  ,    7  , recoZ.eta()-recoJets[1].eta()                     , theWeight);
  theHistograms.fill("recoAll_Zj0_deltaPhi",   "#Delta#phi between reco Z and recoJets[0]",  50,  -3.5,    3.5, physmath::deltaPhi(recoZ.phi(), recoJets[0].phi()), theWeight);
  theHistograms.fill("recoAll_Zj1_deltaPhi",   "#Delta#phi between reco Z and recoJets[1]",  50,  -3.5,    3.5, physmath::deltaPhi(recoZ.phi(), recoJets[1].phi()), theWeight);
  theHistograms.fill("recoAll_Wj0_deltaPhi",   "#Delta#phi between reco Z and recoJets[0]",  50,  -3.5,    3.5, physmath::deltaPhi(recoW.phi(), recoJets[0].phi()), theWeight);
  theHistograms.fill("recoAll_Wj1_deltaPhi",   "#Delta#phi between reco Z and recoJets[1]",  50,  -3.5,    3.5, physmath::deltaPhi(recoW.phi(), recoJets[1].phi()), theWeight);
  /*
  if(recoJets.size() == 2){
    theHistograms.fill("recoAll_WZjpt_DeltaPhi_2jet", "#Delta#phi between recoWZ and recoJets[0]", 50, -3.5, 3.5, WZjptDeltaPhi, theWeight);
  }
  */
  theHistograms.fill("recoAll_trmass_deltaEta",   200, 220  , 8220  , 100, -9  ,    9  , recoPtot.Mt()                                     , recoJJdeltaEta                                    , theWeight);
  theHistograms.fill("recoAll_trmass_Zltrmass",   200, 220  , 8220  , 400,  0  , 1200  , recoPtot.Mt()                                     , recoZlp4.Mt()                                     , theWeight);
  theHistograms.fill("recoAll_trmass_JJtrmass",   200, 220  , 8220  , 600,  0  , 4500  , recoPtot.Mt()                                     , recoJJptot.Mt()                                   , theWeight);
  theHistograms.fill("recoAll_trmass_METtrmass",  200, 220  , 8220  , 400,  0  ,  400  , recoPtot.Mt()                                     , MET.p4().Mt()                                     , theWeight);
  theHistograms.fill("recoAll_METpt_JJdeltaEta",   70,   0  ,  350  , 100, -9  ,    9  , MET.pt()                                          , recoJJdeltaEta                                    , theWeight);
  theHistograms.fill("recoAll_WZpt_JJdeltaEta",   200, 200  , 1000  , 100, -9  ,    9  , recoWZs[0].pt()                                   , recoJJdeltaEta                                    , theWeight);
  theHistograms.fill("recoAll_Zlmass_JJdeltaEta", 400,   0  , 1200  , 100, -9  ,    9  , recoZlp4.M()                                      , recoJJdeltaEta                                    , theWeight);
  theHistograms.fill("recoAll_WJ0vsZJ0_deltaPhi",  50,  -3.5,    3.5,  50, -3.5,    3.5, physmath::deltaPhi(recoW.phi(), recoJets[0].phi()), physmath::deltaPhi(recoZ.phi(), recoJets[0].phi()), theWeight);
  theHistograms.fill("recoAll_WJ1vsZJ1_deltaPhi",  50,  -3.5,    3.5,  50, -3.5,    3.5, physmath::deltaPhi(recoW.phi(), recoJets[1].phi()), physmath::deltaPhi(recoZ.phi(), recoJets[1].phi()), theWeight);
  theHistograms.fill("recoAll_WJ0vsWJ1_deltaPhi",  50,  -3.5,    3.5,  50, -3.5,    3.5, physmath::deltaPhi(recoW.phi(), recoJets[0].phi()), physmath::deltaPhi(recoW.phi(), recoJets[1].phi()), theWeight);
  theHistograms.fill("recoAll_ZJ0vsZJ1_deltaPhi",  50,  -3.5,    3.5,  50, -3.5,    3.5, physmath::deltaPhi(recoZ.phi(), recoJets[0].phi()), physmath::deltaPhi(recoZ.phi(), recoJets[1].phi()), theWeight);
  theHistograms.fill("recoAll_JJvsZJ0_deltaPhi",   50,  -3.5,    3.5,  50, -3.5,    3.5, recoJJdeltaPhi                                    , physmath::deltaPhi(recoZ.phi(), recoJets[0].phi()), theWeight);
  theHistograms.fill("recoAll_JJvsZJ1_deltaPhi",   50,  -3.5,    3.5,  50, -3.5,    3.5, recoJJdeltaPhi                                    , physmath::deltaPhi(recoZ.phi(), recoJets[1].phi()), theWeight);
  theHistograms.fill("recoAll_JJvsWJ0_deltaPhi",   50,  -3.5,    3.5,  50, -3.5,    3.5, recoJJdeltaPhi                                    , physmath::deltaPhi(recoW.phi(), recoJets[0].phi()), theWeight);
  theHistograms.fill("recoAll_JJvsWJ1_deltaPhi",   50,  -3.5,    3.5,  50, -3.5,    3.5, recoJJdeltaPhi                                    , physmath::deltaPhi(recoW.phi(), recoJets[1].phi()), theWeight);

  // Zeppenfeld
  theHistograms.fill("Zeppenfeld_JJZ",      100, -6, 6, zeppenfeldJJZ      , theWeight);
  //theHistograms.fill("Zeppenfeld_JJZ_abs",   50,  0, 6, abs(zeppenfeldJJZ) , theWeight);
  theHistograms.fill("Zeppenfeld_JJl",      100, -6, 6, zeppenfeldJJl      , theWeight);
  //theHistograms.fill("Zeppenfeld_JJl_abs",   50,  0, 6, abs(zeppenfeldJJl) , theWeight);
  theHistograms.fill("Zeppenfeld_llJ0",     100, -6, 6, zeppenfeldllJ0     , theWeight);
  theHistograms.fill("Zeppenfeld_llJ0_abs",  50,  0, 6, abs(zeppenfeldllJ0), theWeight);
  theHistograms.fill("Zeppenfeld_llJ1",     100, -6, 6, zeppenfeldllJ1     , theWeight);
  //theHistograms.fill("Zeppenfeld_llJ1_abs",  50,  0, 6, abs(zeppenfeldllJ1), theWeight);
  theHistograms.fill("Zeppenfeld_lll",      100, -6, 6, zeppenfeldlll      , theWeight);
  //theHistograms.fill("Zeppenfeld_lll_abs",   50,  0, 6, abs(zeppenfeldlll) , theWeight);
  
  
  // ~~~~~~~~~~~~~~~~~~~~~ End of reco Analysis ~~~~~~~~~~~~~~~~~~~~
}

void WZAnalyzer::GenRecoAnalysis(const ZZtype genWZ, const Particle genJet0, const Particle genJet1, const DiBosonLepton recoWZ, const Particle recoJet0, const Particle recoJet1){
  eventGenReco++;
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~ Reco vs Gen ~~~~~~~~~~~~~~~~~~~~~~~~~
  bool choosingW = kTRUE;
  BosonLepton recoWtemp;
  TLorentzVector genWZjjp4 = genWZ.first().p4() + genWZ.second().p4() + genJet0.p4() + genJet1.p4();
  TLorentzVector recoWZjjp4 = recoWZ.first().p4() + recoWZ.second().p4() + recoJet0.p4() + recoJet1.p4();
  vector<DiBosonLepton> recoWZs;
  ZLCompositeCandidates recoZls;
 
  // test if recoW is genW by trmass
  foreach(const ZLCompositeCandidate Z, *ZLCand){
    //if(Z.first.mass() > 60. && Z.second.mass() < 120. && Z.second.passFullSel()){
      recoZls.push_back(Z);
      //}
  }

  if(recoZls.size() > 0){
    sort(recoZls.begin(), recoZls.end(), pairMassComparator(0, ZMASS));
    
    for(int i = 0; choosingW; i++){
      recoWtemp = BosonLepton(recoZls[i].second, Lepton(met->p4()), copysign(24, recoZls[i].second.charge()));

      //if(recoWtemp.p4().Mt() > 30 && recoWtemp.p4().Mt() < 500){
	recoWZs.push_back(DiBosonLepton(recoWtemp, recoZls[i].first));
	//}
      
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
  
  
  // check if gen and reco IDs are the same
  //int genZlID = abs(genWZ.second().daughter(0).id()) + abs(genWZ.second().daughter(1).id()) + abs(genWZ.first().daughter(0).id());
  //int recoZlID = abs(recoWZ.second().daughter(0).id()) + abs(recoWZ.second().daughter(1).id()) + abs(recoWZ.first().daughter(0).id());

  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~ Histograms ~~~~~~~~~~~~~~~~~~~~~~~~~
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

  ///*
  if(genWZ.pt() != 0. && recoWZ.pt() != 0.){
    //Reco vs Gen analysis
    WZAnalyzer::GenRecoAnalysis(genWZ, genJet0, genJet1, recoWZ, recoJet0, recoJet1);
  }
  //*/
}

void WZAnalyzer::end(TFile &){
  
  cout << "\n--------------------------------------------------------------------------" << endl;
  
  cout << "\nEvents of the sample analyzed:                       " << setw(9) << eventSample << endl;
  cout << "Gen events analyzed:                                 " << setw(9) << eventGen << endl;
  cout << "Reco events analyzed:                                " << setw(9) << eventReco << endl;
  cout << "Gen&Reco events analyzed:                            " << setw(9) << eventGenReco << endl;
  cout << "RecoZls empty:                                       " << setw(9) << recoZlempty << endl;
  cout << "RecoJets.size < 2:                                   " << setw(9) << recoJetless2 << endl;
  cout << "Reco events after all cuts:                          " << setw(9) << recoAfterCut << "\t" << recoAfterCut*100./eventReco << endl;
  
  /*
  cout << "\nNumber of                                            " << setw(9) << counter1 << endl;
  cout << "Number of                                            " << setw(9) << counter2 << endl;
  cout << "Number of                                            " << setw(9) << counter3 << endl;
  cout << "Number of                                            " << setw(9) << counter4 << endl;
  */
  
  // execution time
  endtime = ((float)clock())/CLOCKS_PER_SEC;
  WZAnalyzer::printTime(begintime, endtime);
  cout << "\n--------------------------------------------------------------------------" << endl;
}
