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
  nunumber = 0;
  totalevent = 0;
  WZevent = 0;
  zahl = 0;
  recoZlempty = 0;
  wrongrecognition = 0;

  //free counters
  threemuonsplus = 0;
  threemuonsminus = 0;
  threeelesplus = 0;
  threeelesminus = 0;

  begintime = ((float)clock())/CLOCKS_PER_SEC;
}

Int_t WZAnalyzer::cut() {
  return 1;
}


void WZAnalyzer::analyze(){
  zahl++;
  
  vector<Particle> electron;
  vector<Particle> genjets;
  vector<Particle> lepton;
  vector<Particle> muon;
  vector<Particle> neutrino;

  cout << "\n--------------------------------------------------------------------------"<< endl;
  cout << "Run: " << run << " event: " << event << endl;

  foreach(const Particle &gen, *genParticles){
    if((abs(gen.id()) != 11 && abs(gen.id()) != 13 && abs(gen.id()) != 12 && abs(gen.id()) != 14) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) continue;
    //cout << "id: " << gen.id() << " pt: " << gen.pt() << "\t eta: " << gen.eta() << endl;
    /*
    theHistograms.fill("AllGenParticle_id",  "ids all particles",     4, 10.5,  14.5, abs(gen.id()) , theWeight);
    theHistograms.fill("AllGenParticle_pt",  "p_{t} all particles", 300,  0  , 600  , gen.pt()      , theWeight);
    theHistograms.fill("AllGenParticle_Y",   "Y all particles",     100,-10  ,  10  , gen.rapidity(), theWeight);
    theHistograms.fill("AllGenParticle_eta", "#eta all particles",  100,-10  ,  10  , gen.eta()     , theWeight);
    */
    
    if(abs(gen.id()) == 11)      electron.push_back(gen);
    else if(abs(gen.id()) == 13) muon.push_back(gen);
    else if(abs(gen.id()) == 12 || abs(gen.id()) == 14)  neutrino.push_back(gen);
  }
  

  // ~~~~~~ Gen Analysis ~~~~~~
  
  // /*
  // filter on leptons number
  if(electron.size() + muon.size() + neutrino.size() != 4){
    return;
  }
  
  // filter on jets number
  foreach(const Particle jet, *genJets){
    bool leptonMatch = false; //in genJets are included leptons
    foreach(const phys::Particle &gen, *genParticles)
      if(physmath::deltaR(gen,jet) < 0.4 && (abs(gen.id()) == 11 || abs(gen.id()) == 13)) leptonMatch = true;
    
    if(!leptonMatch){
      if(fabs(jet.eta()) < 4.7) genjets.push_back(jet);
    }
  }

  if(genjets.size() < 2){
    cout << "Not enough jets" << endl;
    return;
  }
  // */

  nunumber++;
  
  // ------ filters on leptons ------
  // first filter on pt and eta
  foreach(const Particle ele, electron){
    if(ele.pt() < 7 || abs(ele.eta()) > 2.5){
      cout << Violet("\nElectrons: pt less than 7 GeV or eta's absolute value greater than 2.5") << endl;
      return;
    }
    lepton.push_back(ele);
  }
  
  foreach(const Particle mu, muon){
    if(mu.pt() < 5 || abs(mu.eta()) > 2.4){
      cout << Violet("\nMuons: pt less than 5 GeV or eta's absolute value greater than 2.4") << endl;
      return;
    }
    lepton.push_back(mu);
  }
  
  sort(electron.begin(), electron.end(), PtComparator());
  sort(lepton.begin(), lepton.end(), PtComparator());
  sort(muon.begin(), muon.end(), PtComparator());

  // second filter on pt and eta
  // /*
  if(lepton[0].pt() < 20){
    cout << Violet("\nFirst lepton pt less than 20 GeV") << endl;
    return;
  }
  if(abs(lepton[1].id()) == 11 && lepton[1].pt() < 12){
    cout << Violet("\nSecond lepton is an electron and has pt less than 12 GeV") << endl;
    return;
  }
  if(abs(lepton[1].id()) == 13 && lepton[1].pt() < 10){
    cout << Violet("\nSecond lepton is a muon and has pt less than 10 GeV") << endl;
    return;
  }
  // */
  
  // Z and W must be on shell
  TLorentzVector Ptot = neutrino[0].p4();
  foreach(const Particle lep, lepton)
    Ptot += lep.p4();

  /*
  theHistograms.fill("AllGenlllnu_mass",   "m 3 leptons and #nu",     300, 0, 1500, Ptot.M()  , theWeight);
  theHistograms.fill("AllGenlllnu_trmass", "m_{T} 3 leptons and #nu", 300, 0, 1500, Ptot.Mt(), theWeight);
  */

  if(Ptot.M() < 165){
    cout << Yellow("\nTotal mass of the products insufficient for the WZ analysis.") << endl;
    return;
  }
  
  foreach(const Particle &gen, *genParticles){
    if((abs(gen.id()) != 11 && abs(gen.id()) != 13 && abs(gen.id()) != 12 && abs(gen.id()) != 14) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) continue;
    cout << "id: " << gen.id() << " pt: " << setw(5) << gen.pt() << "\t eta: " << gen.eta() << endl;
  }

  /*
  foreach(const Particle lep, lepton){
    theHistograms.fill("GenL_charge", "leptons charge", 3, -1.5, 1.5, lep.charge(), theWeight);
  }
  */
    
  // ------ Z & W ------
  WZevent++;

  Vtype Weh;
  Vtype Zet;
  vector<Vtype> possibleZ;
  vector<Zltype> Zls;  
  ZZtype WZ;

  // Construction of the two possible Zs and W
  
  if(electron.size()==2 && muon.size()==1 && electron[0].charge() != electron[1].charge()){
    // return;
    // /*
    Zet = Vtype(electron[0], electron[1], 23);
    Zls.push_back(Zltype(Zet, muon[0]));
    Weh = Vtype(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()) );
    // */
  }
  
  if(electron.size()==1 && muon.size()==2 && muon[0].charge() != muon[1].charge()){
    // return;
    // /*
    Zet = Vtype(muon[0], muon[1], 23);
    Zls.push_back(Zltype(Zet, electron[0]));
    Weh = Vtype(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()) );
    // */
  }
  
  if(electron.size()==3){
    //return;
    // /*
    for(int i = 0; i < (int)electron.size() -1; i++){
      for(int j = i + 1; j < (int)electron.size(); j++){
	for(int k = 0; k < (int)electron.size(); k++){
	  if(k != i && k != j){
	    if(electron[i].charge() != electron[j].charge()){
	      possibleZ.push_back(Vtype(electron[i], electron[j], 23));
	      Zls.push_back(Zltype (possibleZ.back(), electron[k]));
	    }
	  }
	}
      }
    }
    // */
    
    // Z is made up of the couple which gives a better Zmass 
    if(Zls.size() < 1){
      cout << Red("No Z formed.") << endl;
      return;
    }
    
    sort(Zls.begin(), Zls.end(), ZlMassComparator(ZMASS));

    Zet = Zls[0].first;

    // W is made up of the remaining lepton and the neutrino
    Weh = Vtype(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()));
    // */
  }
  
  else if(muon.size()==3){
    //return;
    // /*
    for(int i = 0; i < (int)muon.size() -1; i++){
      for(int j = i + 1; j < (int)muon.size(); j++){
	for(int k = 0; k < (int)muon.size(); k++){
	  if(k != i && k != j){
	    if(muon[i].charge() != muon[j].charge()){
	      possibleZ.push_back(Vtype(muon[i], muon[j], 23));
	      Zls.push_back(Zltype (possibleZ.back(), muon[k]));
	    }
	  }
	}
      }
    }
    // */
    
    if(Zls.size() < 1){
      cout << Red("No Z formed.") << endl;
      return;
    }
    
    sort(Zls.begin(), Zls.end(), ZlMassComparator(ZMASS));

    Zet = Zls[0].first;
    Weh = Vtype(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()));

    // */
  }

  //filter on Z and W mass
  /*
  if(Zet.mass() < 60 || Zet.mass() > 120 || Weh.mass() < 50 || Weh.mass() > 110){
    return;
  }
  */

  // W&Z
  WZ = ZZtype(Weh, Zet);
  
  int ZlsID = abs(Zls[0].first.daughter(0).id()) + abs(Zls[0].first.daughter(1).id()) + abs(Zls[0].second.id());
  /*
  double ZdeltaEta = Zet.daughter(0).eta() - Zet.daughter(1).eta();
  //double ZdeltaPhi = physmath::deltaPhi(Zet.daughter(0).phi(), Zet.daughter(1).phi());
  double ZdeltaR = abs(physmath::deltaR(Zet.daughter(0), Zet.daughter(1)));
  
  double WdeltaEta = Weh.daughter(0).eta() - Weh.daughter(1).eta();
  double WdeltaPhi = physmath::deltaPhi(Weh.daughter(0).phi(), Weh.daughter(1).phi());
  double WdeltaR = abs(physmath::deltaR(Weh.daughter(0), Weh.daughter(1)));
  
  double WZdeltaEta = Zet.eta() - Weh.eta();
  //double WZdeltaPhi = physmath::deltaPhi(Zet.phi(), Weh.phi());
  double WZdeltaR = abs(physmath::deltaR(Zet, Weh));
  */
  
  // Histograms and printouts
  
  cout << "\n~~~~~~~~~~~~~~~~~ Gen ~~~~~~~~~~~~~~~~~" << endl;
  cout << "---------------- Z + l ----------------\n" << endl;
  cout << "\nZl candidates are: " << Zls.size() << endl;
  //theHistograms.fill("GenZL_size", "Zl's size", 4, -0.5, 3.5, Zls.size(), theWeight);
  
  foreach(const Zltype zl, Zls){
    cout << "   Z " << zl.first << "\n   l " << zl.second << endl << endl;
  }
  
  cout << "Z is: " << Zet << endl;
  cout << "W is: " << Weh << "\n  her lepton daughter is: " << Weh.daughter(0) << endl;

  /*
  //neutrino histograms
  theHistograms.fill("GenN_pt",     "#nu's p_{t}", 140,  0  , 700  , neutrino[0].pt()    , theWeight);
  theHistograms.fill("GenN_eta",    "#nu's #eta",   90, -6.5,   6.5, neutrino[0].eta()   , theWeight);
  //theHistograms.fill("GenN_charge", "#nu's charge",  5, -2.5,   2.5, neutrino[0].charge(), theWeight); //just to be sure...
  //theHistograms.fill("GenN_phi",    "#nu's #phi",   50, -4  ,   4  , neutrino[0].phi()   , theWeight);
  //theHistograms.fill("GenN_Y",      "#nu's Y",      90, -6.5,   6.5, neutrino[0].rapidity(), theWeight);
  
  //W histograms
  theHistograms.fill("GenW_charge",   "W's charge",       5, -2.5,   2.5, Weh.charge()  , theWeight);
  theHistograms.fill("GenW_mass",     "W's mass",        30, 50  , 110  , Weh.mass()    , theWeight);
  theHistograms.fill("GenW_pt",       "W's p_{t}",      130,  0  , 650  , Weh.pt()      , theWeight);
  theHistograms.fill("GenW_Y",        "W's Y",           50, -5  ,   5  , Weh.rapidity(), theWeight);
  theHistograms.fill("GenW_deltaEta", "W's #Delta#eta",  50, -6.5,   6.5, WdeltaEta     , theWeight);
  theHistograms.fill("GenW_deltaR",   "W's #DeltaR",     50, -0.5,   6.5, WdeltaR       , theWeight);
  //theHistograms.fill("GenW_eta",      "W's #eta",       50, -9, 9, Weh.eta(), theWeight);
  //theHistograms.fill("GenW_phi",      "W's #phi",       50, -4, 4, Weh.phi(), theWeight);
  //theHistograms.fill("GenW_deltaPhi", "W's #Delta#phi", 50, -4, 4, WdeltaPhi, theWeight);
  
  //Z histograms
  theHistograms.fill("GenZ_mass",     "Z's mass",        30, 60  , 120  , Zet.mass()    , theWeight);
  theHistograms.fill("GenZ_pt",       "Z's p_{t}",      130,  0  , 650  , Zet.pt()      , theWeight);
  theHistograms.fill("GenZ_Y",        "Z's Y",           45, -3  ,   3  , Zet.rapidity(), theWeight);
  theHistograms.fill("GenZ_deltaEta", "Z's #Delta#eta",  30, -5  ,   5  , ZdeltaEta     , theWeight);
  theHistograms.fill("GenZ_deltaR",   "Z's #DeltaR",     20, -0.5,   5.5, ZdeltaR       , theWeight);
  //theHistograms.fill("GenZ_charge",   "Z's charge",       5, -2.5,   2.5, Zet.charge()  , theWeight); //just to be sure...
  //theHistograms.fill("GenZ_eta",      "Z's #eta",       100, -7, 7, Zet.eta(), theWeight);
  //theHistograms.fill("GenZ_phi",      "Z's #phi",        50, -4, 4, Zet.phi(), theWeight);
  //theHistograms.fill("GenZ_deltaPhi", "Z's #Delta#phi",  50, -4, 4, ZdeltaPhi, theWeight);
  
  //W&Z histograms
  theHistograms.fill("AllGenWZ_mass", "m 3 leptons and #nu", 268, 160  , 1500  , masslllnu , theWeight);
  theHistograms.fill("WZ_mass",       "W and Z mass",        268, 160  , 1500  , WZ.mass() , theWeight);
  theHistograms.fill("WZ_deltaEta",   "W and Z #Delta#eta",   50,  -9  ,    9  , WZdeltaEta, theWeight);
  theHistograms.fill("WZ_deltaR",     "W and Z #Delta R",     25,  -0.5,    9  , WZdeltaR  , theWeight);
  theHistograms.fill("WZ_ID",         "Zls' ID",               9,  31.5,   40.5, ZlsID     , theWeight);
  //theHistograms.fill("WZ_deltaPhi", "W and Z #Delta#phi", 50, -4, 4, WZdeltaPhi, theWeight);
  */
  
  // ------- Jets -------

  //theHistograms.fill("GenJetsp_number", "number of jets",  13,  -0.5,  12.5, pgenJets->size());

  //jets, pT > 30 GeV and |eta| < 4.7
  
  cout << "\n---------------- Jets -----------------\n" << endl;
  foreach(const Particle jet, genjets){
    
    cout << "ID: " << jet.id() << " pt: " << jet.pt() << "\t eta: " << jet.eta() << endl;

    /*
    theHistograms.fill("GenJet_pt",     "p_{t} jets",   80,  0  , 400  , jet.pt()    , theWeight);
    theHistograms.fill("GenJet_eta",    "#eta jets",    70, -5  ,   5  , jet.eta()   , theWeight);    
    theHistograms.fill("GenJet_charge", "charge jets",   5, -2.5,   2.5, jet.charge(), theWeight);
    //theHistograms.fill("GenJet_mass",   "mass jets",   120,  0  , 120  , jet.mass()  , theWeight);
    //theHistograms.fill("GenJet_Y",   "Y jets",    70, -5, 5, jet.rapidity(), theWeight);
    //theHistograms.fill("GenJet_phi", "#phi jets", 50, -4, 4, jet.phi()     , theWeight);
    */
  }
  
  //theHistograms.fill("GenJets_number", "number of all gen jets",  10,  -0.5,  9.5, genjets.size(), theWeight);

  sort(genjets.begin(), genjets.end(), PtComparator());

  /*
  double JJdeltaEta = genjets[0].eta() - genjets[1].eta();
  double JJdeltaPhi = physmath::deltaPhi(genjets[0].phi(), genjets[1].phi());
  double JJdeltaR = abs(physmath::deltaR(genjets[0], genjets[1]));
  double JJmT = mT(genjets[0], genjets[1]);

  theHistograms.fill("GenJet_deltaEta", "Jets' #Delta#eta", 100, -9  ,   9, JJdeltaEta, theWeight); 
  theHistograms.fill("GenJet_deltaR",   "Jets' #DeltaR",     25, -0.5,   9, JJdeltaR  , theWeight);
  theHistograms.fill("GenJet_deltaPhi", "Jets' #Delta#phi",  50, -4  ,   4, JJdeltaPhi, theWeight);
  theHistograms.fill("GenJet_mT",       "Jets' m_{T}",      200,  0  , 500, JJmT      , theWeight);
  */

  
  // ~~~~~~ Reco Analysis ~~~~~~

  BosonLepton recoZ;
  BosonLepton recoW;
  Lepton MET = Lepton(met->p4());
  vector<Particle> recoJets;
  vector<BosonLepton> possibleW;
  vector<DiBosonLepton> recoWZs;
  ZLCompositeCandidates recoZls;
  
  cout << "\n~~~~~~~~~~~~~~~~~ Reco ~~~~~~~~~~~~~~~~" << endl;  
  cout << "----------------- MET -----------------\n" << endl;
  cout << "Met -> pt = " << MET.pt() << "\t phi: " << MET.phi() << endl;

  //MET histograms
  theHistograms.fill("recoMET_pt",     "MET's p_{t}", 325,  0  , 650  , met->pt()      , theWeight);
  //theHistograms.fill("recoMET_phi",    "MET's #phi",   50, -4  ,   4  , met->phi()     , theWeight);
  //theHistograms.fill("recoMET_charge", "MET's charge",  5, -2.5,   2.5, met->charge()  , theWeight);
  //theHistograms.fill("recoMET_Y",      "MET's Y",      50, -4  ,   4  , met->rapidity(), theWeight);
  //theHistograms.fill("recoMET_eta",    "MET's #eta",  100, -7  ,   7  , met->eta()     , theWeight);

  
  //cout << "---------------- Z + l ----------------\n" << endl;
  foreach(const ZLCompositeCandidate Z, *ZLCand){
    recoZls.push_back(Z);
  }

  theHistograms.fill("recoZl_size", "Reco Zl's size", 15, -0.5, 14.5, recoZls.size(), theWeight);  

  //recoZ and recoW reconstruction
  if(recoZls.size() == 0){
    cout << Yellow("No pair reconstructed") << endl;
    recoZlempty++;
    return;
  }

  else if(recoZls.size() == 1){
    
    int recoZlsID = abs(recoZls[0].first.daughter(0).id()) + abs(recoZls[0].first.daughter(1).id()) + abs(recoZls[0].second.id());
    
    if(recoZlsID != ZlsID || Zls[0].second.id() != recoZls[0].second.id()){
      wrongrecognition++;
      //return;
    }
    
    recoZ = recoZls[0].first;
    recoW = BosonLepton(recoZls[0].second, MET, copysign(24, recoZls[0].second.charge()));
    recoWZs.push_back(DiBosonLepton(recoW, recoZ));
  }
    
  else if(recoZls.size() > 1){

    sort(recoZls.begin(), recoZls.end(), ZlMassComparator(ZMASS));
    int recoZlsID;
    bool choosingZ = kTRUE;
    bool choosingW = kTRUE;
    bool isthesame;

    // /*
    for(int i = 0; choosingZ; i++){
      recoZlsID = abs(recoZls[i].first.daughter(0).id()) + abs(recoZls[i].first.daughter(1).id()) + abs(recoZls[i].second.id());
      isthesame = recoZls[i].second.id() == Zls[0].second.id() && (abs(recoZls[i].second.pt() - Zls[0].second.pt()) < 2.5);

      if(recoZlsID == ZlsID && isthesame){
	choosingZ = kFALSE;
      }

      if(i == (int)recoZls.size() && choosingZ){
	choosingZ = kFALSE;
	wrongrecognition++;
	//return;
      }
    }
    // */
    
    for(int i = 0; choosingW; i++){
      possibleW.push_back(BosonLepton(recoZls[i].second, MET, copysign(24, recoZls[i].second.charge())));
      recoWZs.push_back(DiBosonLepton(possibleW[i], recoZls[i].first));

      if(i == (int)recoZls.size() || recoZls[i].first.mass() != recoZls[i + 1].first.mass()){
	choosingW = kFALSE;
      }
    }
    
    sort(recoWZs.begin(), recoWZs.end(), ZWMassComparator(WMASS));
    
    recoW = recoWZs[0].first();
    recoZ = recoWZs[0].second();
  }
  
  double recoZdeltaEta = recoZ.daughter(0).eta() - recoZ.daughter(1).eta();
  //double recoZdeltaPhi = physmath::deltaPhi(recoZ.daughter(0).phi(), recoZ.daughter(1).phi());
  double recoZdeltaR = abs(physmath::deltaR(recoZ.daughter(0), recoZ.daughter(1)));
  
  double recoWdeltaEta = recoW.daughter(0).eta() - recoW.daughter(1).eta();
  //double recoWdeltaPhi = physmath::deltaPhi(recoW.daughter(0).phi(), recoW.daughter(1).phi());
  double recoWdeltaR = abs(physmath::deltaR(recoW.daughter(0), recoW.daughter(1)));

  double recoWZdeltaPhi = physmath::deltaPhi(recoW.phi(), recoZ.phi());
  double recoWZdeltaEta = recoW.eta() - recoZ.eta();
  double recoWZdeltaR = abs(physmath::deltaR(recoW, recoZ));
  
  //recoZ histograms
  //theHistograms.fill("recoZ_charge",   "recoZ's charge",       5, -2.5,   2.5, recoZ.charge()  , theWeight); //just to be sure...
  theHistograms.fill("recoZ_mass",     "recoZ's mass",       400,  0  , 400  , recoZ.mass()    , theWeight);
  theHistograms.fill("recoZ_pt",       "recoZ's p_{t}",      260,  0  , 650  , recoZ.pt()      , theWeight);
  theHistograms.fill("recoZ_Y",        "recoZ's Y",           45, -3  ,   3  , recoZ.rapidity(), theWeight);
  theHistograms.fill("recoZ_deltaEta", "recoZ's #Delta#eta",  55, -5  ,   5  , recoZdeltaEta   , theWeight);
  theHistograms.fill("recoZ_deltaR",   "recoZ's #DeltaR",     20, -0.5,   5.5, recoZdeltaR     , theWeight);
  //theHistograms.fill("recoZ_eta",      "recoZ's #eta",       100, -7  ,   7  , recoZ.eta()     , theWeight);
  //theHistograms.fill("recoZ_phi",      "recoZ's #phi",        50, -4  ,   4  , recoZ.phi()     , theWeight);
  //theHistograms.fill("recoZ_deltaPhi", "recoZ's #Delta#phi",  50, -4  ,   4  , recoZdeltaPhi   , theWeight);
  
  //recoW histograms
  //theHistograms.fill("recoW_charge",   "recoW's charge",       5, -2.5,   2.5, recoW.charge()  , theWeight);
  theHistograms.fill("recoW_mass",     "recoW's mass",       400,  0  , 400  , recoW.mass()    , theWeight);
  theHistograms.fill("recoW_pt",       "recoW's p_{t}",      260,  0  , 650  , recoW.pt()      , theWeight);
  theHistograms.fill("recoW_Y",        "recoW's Y",           45, -3  ,   3  , recoW.rapidity(), theWeight);
  theHistograms.fill("recoW_deltaEta", "recoW's #Delta#eta",  55, -5  ,   5  , recoWdeltaEta   , theWeight);
  theHistograms.fill("recoW_deltaR",   "recoW's #DeltaR",     20, -0.5,   5.5, recoWdeltaR     , theWeight);
  //theHistograms.fill("recoW_eta",      "recoW's #eta",       100, -7  ,   7  , recoW.eta()     , theWeight);
  //theHistograms.fill("recoW_phi",      "recoW's #phi",        50, -4  ,   4  , recoW.phi()     , theWeight);
  //theHistograms.fill("recoW_deltaPhi", "recoW's #Delta#phi",  50, -4  ,   4  , recoWdeltaPhi   , theWeight);

  //recoW&Z histograms
  theHistograms.fill("recoWZ_deltaPhi", "recoW & Z #Delta#phi", 50, -4, 4, recoWZdeltaPhi, theWeight);
  theHistograms.fill("recoWZ_deltaEta", "recoW & Z #Delta#eta", 55, -5, 5, recoWZdeltaEta, theWeight);
  theHistograms.fill("recoWZ_deltaR",   "recoW & Z #DeltaR",    20, -0.5, 5.5, recoWZdeltaR, theWeight);
  theHistograms.fill("recoWZ_ptot",     "recoW & Z p_{T}",     200, 200, 3000, recoWZs[0].pt(), theWeight);
  theHistograms.fill("recoWZ_mass",     "recoW & Z mass",      300, 0, 900, recoWZs[0].mass(), theWeight);
  

  // ------- Jets -------
  //cout << "----------------- JET -----------------\n" << endl;
  foreach(const Particle jet, *jets){
    recoJets.push_back(jet);
    
    theHistograms.fill("recoJet_charge", "Jets' charge",  19, -8.5,  10.5, jet.charge(), theWeight);
    theHistograms.fill("recoJet_pt",     "Jets' p_{t}",  160,  0  , 400  , jet.pt()    , theWeight);
    theHistograms.fill("recoJet_eta",    "Jets' #eta",    70, -5  ,   5  , jet.eta()   , theWeight);
    //theHistograms.fill("recoJet_Y",      "Jets' Y",       80,  -5  ,   5  , jet.rapidity(), theWeight);
    //theHistograms.fill("recoJet_phi",    "Jets' #phi",    50,  -4  ,   4  , jet.phi()     , theWeight);
  }

  theHistograms.fill("recoJet_number", "Jets' number", 10, -0.5, 9.5, recoJets.size(), theWeight);
  
  if(recoJets.size() < 2){
    totalevent++;
    return;
  }

  sort(recoJets.begin(), recoJets.end(), PtComparator());
  double recoJJdeltaEta = recoJets[0].eta() - recoJets[1].eta();
  double recoJJdeltaPhi = physmath::deltaPhi(recoJets[0].phi(), recoJets[1].phi());
  double recoJJdeltaR = abs(physmath::deltaR(recoJets[0], recoJets[1]));
  TLorentzVector recoJJptot = recoJets[0].p4() + recoJets[1].p4();
  
  theHistograms.fill("recoJets_deltaEta", "Jets #Delta#eta",      100, -9  , 9, recoJJdeltaEta, theWeight); 
  theHistograms.fill("recoJets_deltaR",   "Jets #DeltaR",          25, -0.5, 9, recoJJdeltaR  , theWeight);
  theHistograms.fill("recoJets_deltaPhi", "Jets #Delta#phi",       50, -4  , 4, recoJJdeltaPhi, theWeight);
  theHistograms.fill("recoJets_mass",     "Jets mass",            600,  0  , 4500, recoJJptot.M(), theWeight);
  theHistograms.fill("recoJets_trmass",   "Jets transverse mass", 600,  0  , 4500, recoJJptot.Mt(), theWeight);
  theHistograms.fill("recoJets_pt",       "Jets p_{tot}",         600,  0  , 2000, recoJJptot.Pt(), theWeight);

  /*
  bool softjets01 = recoJets[0].pt() < 50 && recoJets[1].pt() < 50;
  bool softjets2;
  bool softjets3;
  bool softjets4;
  bool softjets5;
  
  theHistograms.fill("AAA_softjets_01", "test 1", 4, -1.5, 2.5, softjets01);
  
  if(recoJets.size() > 2){
    softjets2 = recoJets[2].pt() < 50;
    theHistograms.fill("AAA_softjets_2", "test 2", 4, -1.5, 2.5, softjets2);
    
    if(recoJets.size() > 3){
      softjets3 = recoJets[3].pt() < 50;
      theHistograms.fill("AAA_softjets_3", "test 3", 4, -1.5, 2.5, softjets3);
      
      if(recoJets.size() > 4){
	softjets4 = recoJets[4].pt() < 50; 
	theHistograms.fill("AAA_softjets_4", "test 4", 4, -1.5, 2.5, softjets4);
      
	if(recoJets.size() > 5){
	  softjets5 = recoJets[5].pt() < 50; 
	  theHistograms.fill("AAA_softjets_5", "test 4", 4, -1.5, 2.5, softjets5);
	}
      }
    }
  }
  */
  
  // ----- All reco -----
  // /*
  TLorentzVector recoPtot = recoW.p4() + recoZ.p4() + recoJets[0].p4() + recoJets[1].p4();

  theHistograms.fill("recoAll_mass",   "Mass recoW,Z,J,J",            200, 220, 8220, recoPtot.M() , theWeight);
  theHistograms.fill("recoAll_trmass", "Transverse mass recoW,Z,J,J", 200, 220, 8220, recoPtot.Mt(), theWeight);
  // */
}
  
void WZAnalyzer::end(TFile &){

  cout << "\n--------------------------------------------------------------------------" << endl;
  
  cout << "\nNumber of events of the sample:                                 " << setw(7) << zahl << endl;
  cout << "Number of events ending with 3 leptons, 1 neutrino and >2 jets: " << setw(7) << nunumber << endl;
  cout << "Number of events useful for WZ analysis:                        " << setw(7) << WZevent << endl;
  cout << "Number of events with recoZls empty:                            " << setw(7) << recoZlempty << endl;
  cout << "Number of recoZls' ids != Zls' ids:                             " << setw(7) << wrongrecognition << endl;
  cout << "Number of recoJets.size < 2:                                    " << setw(7) << totalevent << endl;
  
  /*
  cout << "\nNumber of                         " << setw(7) << threeelesplus << endl;
  cout << "Number of                         " << setw(7) << threeelesminus << endl;
  cout << "Number of                         " << setw(7) << threemuonsplus << endl;
  cout << "Number of                         " << setw(7) << threemuonsminus << endl;
  */
  
  // execution time
  endtime = ((float)clock())/CLOCKS_PER_SEC;
  WZAnalyzer::printTime(begintime, endtime);
  cout << "\n--------------------------------------------------------------------------" << endl;
}
