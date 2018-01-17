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

  nunumber = 0;
  totalevent = 0;
  WZevent = 0;
  zahl = 0;

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
    theHistograms.fill("AllGenParticle_id",  "ids all particles",     4, 10.5,  14.5, abs(gen.id()) , theWeight);
    theHistograms.fill("AllGenParticle_pt",  "p_{t} all particles", 600,  0  , 600  , gen.pt()      , theWeight);
    theHistograms.fill("AllGenParticle_Y",   "Y all particles",     100,-10  ,  10  , gen.rapidity(), theWeight);
    theHistograms.fill("AllGenParticle_eta", "#eta all particles",  100,-10  ,  10  , gen.eta()     , theWeight);
    
    if(abs(gen.id()) == 11)      electron.push_back(gen);
    else if(abs(gen.id()) == 13) muon.push_back(gen);
    else if(abs(gen.id()) == 12 || abs(gen.id()) == 14)  neutrino.push_back(gen);
  }
  

  // ~~~~~~ WZ Analysis ~~~~~~
  
  // /*
  if(electron.size() + muon.size() + neutrino.size() != 4){
    return;
  }
  
  theHistograms.fill("GenJets_number", "number of all gen jets",  13,  -0.5,  12.5, genJets->size(), theWeight);
  if(genJets->size() != 2){
    cout << "Wrong number of jets" << endl;
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
  
  stable_sort(electron.begin(), electron.end(), PtComparator());
  stable_sort(lepton.begin(), lepton.end(), PtComparator());
  stable_sort(muon.begin(), muon.end(), PtComparator());

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
  
  totalevent++;
  
  // Z and W must be on shell
  TLorentzVector Ptot = neutrino[0].p4();
  foreach(const Particle lep, lepton)
    Ptot += lep.p4();

  double masslllnu = Ptot.M();
  double trmasslllnu = Ptot.Mt();

  theHistograms.fill("AllGenlllnu_mass",   "m 3 leptons and #nu",     500, 0, 1500, masslllnu  , theWeight); //what happens between 90 and 160 GeV?
  theHistograms.fill("AllGenlllnu_trmass", "m_{T} 3 leptons and #nu", 500, 0, 1500, trmasslllnu, theWeight);

  if(masslllnu < 165){
    cout << Yellow("\nTotal mass of the products insufficient for the WZ analysis.") << endl;
    return;
  }
  
  foreach(const Particle &gen, *genParticles){
    if((abs(gen.id()) != 11 && abs(gen.id()) != 13 && abs(gen.id()) != 12 && abs(gen.id()) != 14) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) continue;
    cout << "id: " << gen.id() << " pt: " << setw(5) << gen.pt() << "\t eta: " << gen.eta() << endl;
  }
  
  foreach(const Particle lep, lepton){
    theHistograms.fill("GenL_charge", "leptons charge", 3, -1.5, 1.5, lep.charge(), theWeight);
  }
  
    
  // ------ Z & W ------
  WZevent++;
  
  Vtype Weh;
  Vtype Zet;
  vector<Vtype> possibleZ;
  vector<Zltype> Zls;  

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
    
    if(Zls.size() > 1){
      stable_sort(Zls.begin(), Zls.end(), ZlMassComparator(ZMASS));
    }

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
    
    if(Zls.size() > 1){
      stable_sort(Zls.begin(), Zls.end(), ZlMassComparator(ZMASS));
    }

    Zet = Zls[0].first;
    Weh = Vtype(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()));

    // */
  }

  double ZdeltaEta = Zet.daughter(0).eta() - Zet.daughter(1).eta();
  //double ZdeltaPhi = physmath::deltaPhi(Zet.daughter(0).phi(), Zet.daughter(1).phi());
  double ZdeltaR = physmath::deltaR(Zet.daughter(0), Zet.daughter(1));
  
  double WdeltaEta = Weh.daughter(0).eta() - Weh.daughter(1).eta();
  //double WdeltaPhi = physmath::deltaPhi(Weh.daughter(0).phi(), Weh.daughter(1).phi());
  double WdeltaR = physmath::deltaR(Weh.daughter(0), Weh.daughter(1));
  
  double WZdeltaEta = Zet.eta() - Weh.eta();
  //double WZdeltaPhi = physmath::deltaPhi(Zet.phi(), Weh.phi());
  double WZdeltaR = physmath::deltaR(Zet, Weh);

  
  // Histograms and printouts
  
  cout << "\n~~~~~~~~~~~~~~~~~ Gen ~~~~~~~~~~~~~~~~~" << endl;
  cout << "---------------- Z + l ----------------\n" << endl;
  cout << "\nZl candidates are: " << Zls.size() << endl;
  theHistograms.fill("GenZL_size", "Zl's size", 4, -0.5, 3.5, Zls.size(), theWeight);
  
  foreach(const Zltype zl, Zls){
    cout << "   Z " << zl.first << "\n   l " << zl.second << endl << endl;
  }
  
  cout << "Z is: " << Zet << endl;
  cout << "W is: " << Weh << "\n  her lepton daughter is: " << Weh.daughter(0) << endl;

  //neutrino histograms
  theHistograms.fill("GenN_charge", "#nu's charge",  5, -2.5,   2.5, neutrino[0].charge(), theWeight); //just to be sure...
  theHistograms.fill("GenN_pt",     "#nu's p_{t}", 350,  0  , 700  , neutrino[0].pt()    , theWeight);
  theHistograms.fill("GenN_eta",    "#nu's #eta",   90, -6.5,   6.5, neutrino[0].eta()   , theWeight);
  //theHistograms.fill("GenN_Y",      "#nu's Y",      90, -6.5,   6.5, neutrino[0].rapidity(), theWeight);
  //theHistograms.fill("GenN_phi",    "#nu's #phi",   50, -4  ,   4  , neutrino[0].phi()     , theWeight);
  
  //W histograms
  theHistograms.fill("GenW_charge",   "W's charge",       5, -2.5,   2.5, Weh.charge()  , theWeight);
  theHistograms.fill("GenW_mass",     "W's mass",       400,  0  , 400  , Weh.mass()    , theWeight);
  theHistograms.fill("GenW_pt",       "W's p_{t}",      325,  0  , 650  , Weh.pt()      , theWeight);
  theHistograms.fill("GenW_Y",        "W's Y",           50, -5  ,   5  , Weh.rapidity(), theWeight);
  theHistograms.fill("GenW_deltaEta", "W's #Delta#eta",  50, -6.5,   6.5, WdeltaEta     , theWeight);
  theHistograms.fill("GenW_deltaR",   "W's #DeltaR",     50, -0.5,   6.5, WdeltaR       , theWeight);
  //theHistograms.fill("GenW_eta",      "W's #eta",       50, -9, 9, Weh.eta(), theWeight);
  //theHistograms.fill("GenW_phi",      "W's #phi",       50, -4, 4, Weh.phi(), theWeight);
  //theHistograms.fill("GenW_deltaPhi", "W's #Delta#phi", 50, -4, 4, WdeltaPhi, theWeight);
  
  //Z histograms
  theHistograms.fill("GenZ_charge",   "Z's charge",       5, -2.5,   2.5, Zet.charge()  , theWeight); //just to be sure...
  theHistograms.fill("GenZ_mass",     "Z's mass",       400,  0  , 400  , Zet.mass()    , theWeight);
  theHistograms.fill("GenZ_pt",       "Z's p_{t}",      325,  0  , 650  , Zet.pt()      , theWeight);
  theHistograms.fill("GenZ_Y",        "Z's Y",           45, -3  ,   3  , Zet.rapidity(), theWeight);
  theHistograms.fill("GenZ_deltaEta", "Z's #Delta#eta",  30, -5  ,   5  , ZdeltaEta     , theWeight);
  theHistograms.fill("GenZ_deltaR",   "Z's #DeltaR",     20, -0.5,   5.5, ZdeltaR       , theWeight);
  //theHistograms.fill("GenZ_eta",      "Z's #eta",       100, -7, 7, Zet.eta(), theWeight);
  //theHistograms.fill("GenZ_phi",      "Z's #phi",        50, -4, 4, Zet.phi(), theWeight);
  //theHistograms.fill("GenZ_deltaPhi", "Z's #Delta#phi",  50, -4, 4, ZdeltaPhi, theWeight);
  
  //W&Z histograms
  theHistograms.fill("AllGenWZ_mass", "m 3 leptons and #nu", 450, 160  , 1500, masslllnu , theWeight);
  theHistograms.fill("WZ_deltaEta",   "W and Z #Delta#eta",   50,  -9  ,    9, WZdeltaEta, theWeight);
  theHistograms.fill("WZ_deltaR",     "W and Z #Delta R",     25,  -0.5,    9, WZdeltaR  , theWeight);
  //theHistograms.fill("WZ_deltaPhi", "W and Z #Delta#phi", 50, -4, 4, WZdeltaPhi, theWeight);
  
  
  // ------- Jets -------

  //theHistograms.fill("GenJetsp_number", "number of jets",  13,  -0.5,  12.5, pgenJets->size());

  //jets, pT > 30 GeV and |eta| < 4.7
  
  cout << "\n---------------- Jets -----------------\n" << endl;
  foreach(const Particle jet, *genJets){
    genjets.push_back(jet);

    cout << "ID: " << jet.id() << " pt: " << jet.pt() << "\t eta: " << jet.eta() << endl;

    theHistograms.fill("GenJet_pt",     "p_{t} jets",  350,  0  , 400  , jet.pt()    , theWeight);
    theHistograms.fill("GenJet_eta",    "#eta jets",    70, -5  ,   5  , jet.eta()   , theWeight);    
    theHistograms.fill("GenJet_charge", "charge jets",   5, -2.5,   2.5, jet.charge(), theWeight);
    theHistograms.fill("GenJet_mass",   "mass jets",   120,  0  , 120  , jet.mass()  , theWeight);
    //theHistograms.fill("GenJet_Y",   "Y jets",    70, -5, 5, jet.rapidity(), theWeight);
    //theHistograms.fill("GenJet_phi", "#phi jets", 50, -4, 4, jet.phi()     , theWeight);
  }

  if(genjets.size() == 2){
    double JJdeltaEta = genjets[0].eta() - genjets[1].eta();
    //double JJdeltaPhi = physmath::deltaPhi(genjets[0].phi(), genjets[1].phi());
    double JJdeltaR = physmath::deltaR(genjets[0], genjets[1]);

    theHistograms.fill("GenJets_deltaEta", "Jets #Delta#eta", 50, -9  , 9, JJdeltaEta, theWeight); 
    theHistograms.fill("GenJets_deltaR",   "Jets #DeltaR",    25, -0.5, 9, JJdeltaR  , theWeight);
    //theHistograms.fill("GenJets_deltaPhi", "Jets #Delta#phi", 50, -4, 4, JJdeltaPhi, theWeight);  
  }

  // ------- Reco -------

  vector<Particle> recoJets;
  
  Vtype recoW; //to be modified
  BosonLepton recoZ;
  ZLCompositeCandidates recoZls;
  
  cout << "\n~~~~~~~~~~~~~~~~~ Reco ~~~~~~~~~~~~~~~~" << endl;
  //cout << "---------------- Z + l ----------------\n" << endl;
  foreach(const ZLCompositeCandidate Z, *ZLCand){
    recoZls.push_back(Z);
    
    theHistograms.fill("recoZl_Z_charge", "Z's charge",   5, -2.5,   2.5, Z.first.charge()  , theWeight); //just to be sure...
    theHistograms.fill("recoZl_Z_mass",   "Z's mass",   400,  0  , 400  , Z.first.mass()    , theWeight);
    theHistograms.fill("recoZl_Z_pt",     "Z's p_{t}",  325,  0  , 650  , Z.first.pt()      , theWeight);
    theHistograms.fill("recoZl_Z_Y",      "Z's Y",       45, -3  ,   3  , Z.first.rapidity(), theWeight);
    //theHistograms.fill("recoZl_Z_eta",    "Z's #eta",   100, -7  ,   7  , Z.first.eta(), theWeight);
    //theHistograms.fill("recoZl_Z_phi",    "Z's #phi",    50, -4  ,   4  , Z.first.phi(), theWeight);
    
    theHistograms.fill("recoZl_l_charge", "l's charge",  5, -2.5,   2.5, Z.second.charge() , theWeight);
    theHistograms.fill("recoZl_l_id",     "l's id",      5,  9.5,  14.5, abs(Z.second.id()), theWeight);
    theHistograms.fill("recoZl_l_pt",     "l's p_{t}", 325,  0  , 650  , Z.second.pt()     , theWeight);
    theHistograms.fill("recoZl_l_eta",    "l's #eta",   45, -3  ,   3  , Z.second.eta()    , theWeight);
    //theHistograms.fill("recoZl_l_Y",      "l's Y",      45, -3  ,   3  , Z.second.rapidity(), theWeight);
    //theHistograms.fill("recoZl_l_phi",    "l's #phi",   50, -4  ,   4  , Z.second.phi()     , theWeight);
  }

  theHistograms.fill("recoZl_size", "Reco Zl's size", 9, -0.5, 8.5, recoZls.size());
  
  //cout << "----------------- MET -----------------\n" << endl;
  theHistograms.fill("recoMET_charge", "MET's charge",  5, -2.5,   2.5, met->charge()  , theWeight);
  theHistograms.fill("recoMET_pt",     "MET's p_{t}", 325,  0  , 650  , met->pt()      , theWeight);
  theHistograms.fill("recoMET_Y",      "MET's Y",      50, -4  ,   4  , met->rapidity(), theWeight);
  theHistograms.fill("recoMET_eta",    "MET's #eta",  100, -7  ,   7  , met->eta()     , theWeight);

  //recoZ and recoW reconstruction
  if(recoZls.size() == 0){
    cout << Yellow("No pair reconstructed") << endl;
    if(Zls.size() == 1)
      threemuonsplus++;
    if(Zls.size() == 2)
      threemuonsminus++;
  }
  
  if(recoZls.size() == 1){
    
    //Z vs recoZ
    recoZ = recoZls[0].first;
    double ZrecoZdeltaPt = Zet.pt() - recoZ.pt();

    if(Zls.size() == 1){
      //threeelesplus++;
      theHistograms.fill("ZvsrecoZ_21_deltapt", "#Deltap_{t} between Z and reco Z", 75, -11, 11, ZrecoZdeltaPt, theWeight);
    }
    /* //no recoZls.size()=1 matches Zls.size()=2
    if(Zls.size() == 2)
      threeelesminus++;
    */

    //W vs recoW
    //recoW = Boson(recoZls[0].second, met, copysign(24, recoZls[0].second.charge())); //doesn't work: recoZls.second is Lepton and met is Particle
  }

  if(recoZls.size() == 4){
    if(Zls.size() == 1)
      threeelesplus++;
    if(Zls.size() == 4)
      threeelesminus++;
  }
  
  //cout << "----------------- JET -----------------\n" << endl;
  foreach(const Particle jet, *jets){
    recoJets.push_back(jet);

    theHistograms.fill("recoJet_charge", "Jets' charge",  19, -8.5,  10.5, jet.charge(), theWeight);
    theHistograms.fill("recoJet_pt",     "Jets' p_{t}",  350,  0  , 400  , jet.pt()    , theWeight);
    theHistograms.fill("recoJet_eta",    "Jets' #eta",    80, -5  ,   5  , jet.eta()   , theWeight);
    //theHistograms.fill("recoJet_Y",      "Jets' Y",       80,  -5  ,   5  , jet.rapidity(), theWeight);
    //theHistograms.fill("recoJet_phi",    "Jets' #phi",    50,  -4  ,   4  , jet.phi()     , theWeight);
  }

  theHistograms.fill("recoJet_number", "Jets' number", 4, -0.5, 3.5, recoJets.size(), theWeight);    
}
  
void WZAnalyzer::end(TFile &){

  cout << "\n--------------------------------------------------------------------------" << endl;
  
  cout << "\nNumber of events of the sample:                         " << setw(7) << zahl << endl;
  cout << "Number of events ending with 3 leptons and 1 neutrino:  " << setw(7) << nunumber << endl;
  cout << "Number of events useful for Wlllnu and WZ analysis:     " << setw(7) << totalevent << endl;
  cout << "Number of events useful for WZ analysis:                " << setw(7) << WZevent << endl;

  
  cout << "\nNumber of recoZl > 1, Zl = 1 " << setw(7) << threeelesplus << endl;
  cout << "Number of recoZl > 1, Zl = 2 " << setw(7) << threeelesminus << endl;
  cout << "Number of recoZl = 0, Zl = 1 " << setw(7) << threemuonsplus << endl;
  cout << "Number of recoZl = 0, Zl = 2 " << setw(7) << threemuonsminus << endl;
  
  
  // execution time
  endtime = ((float)clock())/CLOCKS_PER_SEC;
  WZAnalyzer::printTime(begintime, endtime);
  cout << "\n--------------------------------------------------------------------------" << endl;
}
