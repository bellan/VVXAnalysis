#include "VVXAnalysis/TreeAnalysis/interface/WZAnalyzer.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include <chrono>
#include <ctime> 

using namespace physmath;
using namespace phys;
using namespace std;


void WZAnalyzer::GenAnalysis(DiBosonParticle &WZ, Particle &Jet0, Particle &Jet1){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~ Begin of gen Analysis ~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int cut = 0;
  theHistograms.fill("GenCuts",     "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);
  

  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Z & W ~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  // DiBoson variables
  zz::SignalTopology topology = zz::getSignalTopology(*genParticles, *genJets, *genJetsAK8);

  BosonParticle Zet = get<1>(topology);
  BosonParticle Weh = get<5>(topology);
  bitset<16> top = (bitset<16>)get<0>(topology);
  tuple<bool, BosonParticle, BosonParticle> WZtuple;  

  theHistograms.fill("GenTopology", "Topology di WZ", 10, -0.5, 9.5, top.test(0));
  theHistograms.fill("GenTopology", "Topology di WZ", 10, -0.5, 9.5, top.test(1)+2);
  theHistograms.fill("GenTopology", "Topology di WZ", 10, -0.5, 9.5, top.test(5)+4);
  theHistograms.fill("GenTopology", "Topology di WZ", 10, -0.5, 9.5, top.test(6)+6);
  theHistograms.fill("GenTopology", "Topology di WZ", 10, -0.5, 9.5, top.test(7)+8);


  
  // ----- Filter on leptons, total mass >100GeV, smart cut on
  if(top.test(0) == 1 || top.test(1) == 0 || top.test(5) == 0 || top.test(6) == 0 || top.test(7) == 0)
    return;
  
  cut++;
  theHistograms.fill("GenCuts",     "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);

  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Jets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ----- filter on jets' number and eta  
  helper_->FindLeadingJets(genJets, Jet0, Jet1, genParticles);

  if(Jet0.p4().Mt() == 0)
    return;

  cut++;
  theHistograms.fill("GenCuts",     "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);
  
  // Leading jets variables
  TLorentzVector JJp4 = Jet0.p4() + Jet1.p4();
  double JJdeltaEta = Jet0.eta() - Jet1.eta();



  // ~~~~~~~~~~~~~~~~~~~ Histograms before cuts ~~~~~~~~~~~~~~~~~~~~  
  eventGen++;
  weightGen += theWeight;
  theHistograms.fill("WZ_Events", "Weighted counters", 10, -0.5, 10.5, 2, theWeight);
  
  helper_->PlotBoson(Zet, "GenZ", theWeight, "BC");
  helper_->PlotBoson(Weh, "GenW", theWeight, "BC");
  helper_->PlotDiBoson(DiBosonParticle(Weh, Zet), "GenWZ", theWeight, "BC");
  helper_->PlotJets(Jet0, Jet1, "Gen", theWeight, "BC");

  
  
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

  
  if(filtromassaZ || filtroMET || filtrotrmassaWmin || filtrotrmassaWmax || filtroJJdeltaEta || filtroJJpt || filtrozeppenfeldllJ0 || filtroJ1pt || filtroWleppt)
    return;
  
  /*
  if(filtroMET) return;
  cut++;
  theHistograms.fill("GenCuts",     "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);

  if(filtrotrmassaWmin || filtrotrmassaWmax) return;
  cut++;
  theHistograms.fill("GenCuts",     "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);
  
  if(filtroJJdeltaEta) return;
  cut++;
  theHistograms.fill("GenCuts",     "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);
  
  if(filtroJJpt) return;
  cut++;
  theHistograms.fill("GenCuts",     "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);
  
  if(filtrozeppenfeldllJ0) return;
  cut++;
  theHistograms.fill("GenCuts",     "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);
  
  if(filtroJ1pt) return;
  cut++;
  theHistograms.fill("GenCuts",     "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);
  
  if(filtroWleppt) return;
  cut++;
  theHistograms.fill("GenCuts",     "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);
  
  if(filtromassaZ) return;*/
  cut++;
  theHistograms.fill("GenCuts",     "Gen events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("GenCuts_wei", "Gen events after cuts", 14, -0.5, 13.5, cut, theWeight);



  // ~~~~~~~~~~~~~~~~~~~~ Histograms after cuts ~~~~~~~~~~~~~~~~~~~~
  eventGenaftercut++;
  weightGenaftercut += theWeight;
  theHistograms.fill("WZ_Events", "Weighted counters", 10, -0.5, 10.5, 3, theWeight);
  WZ = DiBosonParticle(Weh, Zet);
  
  helper_->PlotBoson(Zet, "GenZ", theWeight, "AC");
  helper_->PlotBoson(Weh, "GenW", theWeight, "AC");
  helper_->PlotDiBoson(WZ, "GenWZ", theWeight, "AC");
  helper_->PlotJets(Jet0, Jet1, "Gen", theWeight, "AC"); 

  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~ End of gen Analysis ~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}


void WZAnalyzer::RecoAnalysis(DiBosonLepton &recoWZ, Particle &Jet0, Particle &Jet1){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~ Begin of reco Analysis ~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int cut = 0;
  theHistograms.fill("RecoCuts",     "Reco events after cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("RecoCuts_wei", "Reco events after cuts", 20, -0.5, 19.5, cut, theWeight);
  

  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Z & W ~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  // DiBoson variables
  DiBosonLepton tempZW;


  
  // ----- Filter on WZCand in SR
  // ZWCand size > 0
  if(ZW->mass() == 0)
    return;

  cut++;
  theHistograms.fill("RecoCuts",     "Reco events after cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("RecoCuts_wei", "Reco events after cuts", 20, -0.5, 19.5, cut, theWeight);

  // 3 leptons
  if(electrons->size() + muons->size() != 3)
    return;
  
  cut++;
  theHistograms.fill("RecoCuts",     "Reco events after cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("RecoCuts_wei", "Reco events after cuts", 20, -0.5, 19.5, cut, theWeight);

  // 3rd lepton full selection
  tempZW = DiBosonLepton(*ZW);
  if(!tempZW.second().daughter(0).passFullSel())
    return;
  
  cut++;
  theHistograms.fill("RecoCuts",     "Reco events after cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("RecoCuts_wei", "Reco events after cuts", 20, -0.5, 19.5, cut, theWeight);

  // 3rd lepton pt > 30
  if(tempZW.second().daughter(0).pt() <30.)
    return;
  
  cut++;
  theHistograms.fill("RecoCuts",     "Reco events after cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("RecoCuts_wei", "Reco events after cuts", 20, -0.5, 19.5, cut, theWeight);

  // MET pt > 30
  if(tempZW.second().daughter(1).pt() <30.)
    return;
  
  cut++;
  theHistograms.fill("RecoCuts",     "Reco events after cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("RecoCuts_wei", "Reco events after cuts", 20, -0.5, 19.5, cut, theWeight);

  // Z's mass 60 < m < 120
  if(tempZW.first().mass() < 60. || tempZW.first().mass() > 120.)
    return;
  
  cut++;
  theHistograms.fill("RecoCuts",     "Reco events after cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("RecoCuts_wei", "Reco events after cuts", 20, -0.5, 19.5, cut, theWeight);

  // W's trmass 30 < trm < 500
  if(tempZW.second().p4().Mt() < 30. || tempZW.second().p4().Mt() > 500.)
    return;
  
  cut++;
  theHistograms.fill("RecoCuts",     "Reco events after cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("RecoCuts_wei", "Reco events after cuts", 20, -0.5, 19.5, cut, theWeight);

  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Jets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ----- filter on jets' number and eta  
  if(jets->size() < 2)
    return;
  
  cut++;
  theHistograms.fill("RecoCuts",     "Reco events after cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("RecoCuts_wei", "Reco events after cuts", 20, -0.5, 19.5, cut, theWeight);

  helper_->FindLeadingJets(jets, Jet0, Jet1);
  
  cut++;
  theHistograms.fill("RecoCuts",     "Reco events after cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("RecoCuts_wei", "Reco events after cuts", 20, -0.5, 19.5, cut, theWeight);


  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ All ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  eventReco++;
  weightReco++;
  
  // Variables
  TLorentzVector recoPtot = tempZW.p4() + Jet0.p4() + Jet1.p4();
  TLorentzVector recoJJptot = Jet0.p4() + Jet1.p4();
  
  float zeppenfeldllJ0 = Jet0.eta() - (tempZW.first().daughter(0).eta() + tempZW.first().daughter(1).eta())/2;
  float recoJJdeltaEta = Jet0.eta() - Jet1.eta();

  
  
  // ------------------- Histograms before cuts --------------------
  helper_->PlotJets(Jet0, Jet1, "Reco", theWeight, "BC");
  
  theHistograms.fill("WZ_Events", "Weighted counters", 10, -0.5, 10.5, 4, theWeight);
  
  theHistograms.fill("RecoAll_trmass_BC",     "Transverse mass recoW,Z,J", 34, 220, 7095, recoPtot.Mt() , theWeight);
  theHistograms.fill("RecoAll_Zeppenfeld_BC", "Zeppenfeld variable",       50,  -6,    6, zeppenfeldllJ0, theWeight);

  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Cut ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  // DeltaEta JJ > 2.4
  if(abs(recoJJdeltaEta) < 2.4)
    return;
  
  cut++;
  theHistograms.fill("RecoCuts",     "Reco events after cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("RecoCuts_wei", "Reco events after cuts", 20, -0.5, 19.5, cut, theWeight);

  // JJmass > 280
  if(recoJJptot.M() < 280)
    return;

  cut++;
  theHistograms.fill("RecoCuts",     "Reco events after cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("RecoCuts_wei", "Reco events after cuts", 20, -0.5, 19.5, cut, theWeight);

  // filter Zeppenfeld > 0.6
  if(abs(zeppenfeldllJ0) < 0.6)
    return;

  cut++;
  theHistograms.fill("RecoCuts",     "Reco events after cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("RecoCuts_wei", "Reco events after cuts", 20, -0.5, 19.5, cut, theWeight);

  // filter Jet1 pt > 50
  if(Jet1.pt() < 50.)
    return;

  cut++;
  theHistograms.fill("RecoCuts",     "Reco events after cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("RecoCuts_wei", "Reco events after cuts", 20, -0.5, 19.5, cut, theWeight);

  // filter Zmass < 15
  if(abs(ZMASS - tempZW.first().mass() > 15))
    return;

  cut++;
  theHistograms.fill("RecoCuts",     "Reco events after cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("RecoCuts_wei", "Reco events after cuts", 20, -0.5, 19.5, cut, theWeight);


  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ All ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  eventRecoaftercut++;
  weightRecoaftercut += theWeight;

  // Variables
  recoWZ = DiBosonLepton(tempZW.second(), tempZW.first());
  
  TLorentzVector recoZlp4 = recoWZ.first().daughter(0).p4() + recoWZ.second().daughter(0).p4() + recoWZ.second().daughter(1).p4();
  
  int recoZlID = abs(recoWZ.first().daughter(0).id()) + abs(recoWZ.first().daughter(1).id()) + abs(recoWZ.second().daughter(0).id());


  
  // -------------------- Histograms after cuts --------------------  
  theHistograms.fill("WZ_Events", "Weighted counters", 10, -0.5, 10.5, 5, theWeight);
  
  theHistograms.fill("RecoAll_trmass_AC", "Transverse mass recoW,Z,J", 34, 220, 7095, recoPtot.Mt(), theWeight);


  
  // ----- Particles
  helper_->PlotParticle(*met, "RecoMET", theWeight, "AC");
  
  helper_->PlotBoson(recoWZ.first(),  "RecoW", theWeight, "AC");
  helper_->PlotBoson(recoWZ.second(), "RecoZ", theWeight, "AC");
  
  helper_->PlotDiBoson(recoWZ, "RecoWZ", theWeight, "AC");
  
  helper_->PlotJets(Jet0, Jet1, "Reco", theWeight, "AC");


  
  // ----- Zl 
  theHistograms.fill("RecoZl_mass",   "3 leptons mass",           400,  0  , 1200  , recoZlp4.M()                    , theWeight);
  theHistograms.fill("RecoZl_1st_pt", "Z's 1^{st} lepton p_{t}",  200,  0  ,  400  , recoWZ.second().daughter(0).pt(), theWeight);
  theHistograms.fill("RecoZl_2nd_pt", "Z's 2^{nd} lepton p_{t}",  200,  0  ,  400  , recoWZ.second().daughter(1).pt(), theWeight);
  theHistograms.fill("RecoZl_3rd_pt", "W's lepton p_{t}",          50,  0  ,  400  , recoWZ.first().daughter(0).pt() , theWeight);
  theHistograms.fill("RecoZl_ID",     "ID of WZ leptons daughters", 6, 30.5,   41.5, recoZlID                        , theWeight);

  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~ End of reco Analysis ~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}


void WZAnalyzer::GenRecoAnalysis(const DiBosonParticle genWZ, const Particle genJet0, const Particle genJet1, const DiBosonLepton recoWZ, const Particle recoJet0, const Particle recoJet1){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~ Begin of Reco vs Gen ~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  eventGenReco++;
  weightGenReco += theWeight;
  theHistograms.fill("WZ_Events", "Weighted counters", 10, -0.5, 10.5, 6, theWeight);
  
  // check if genZ and recoZ IDs are the same
  int genZID = abs(genWZ.second().daughter(0).id()) + abs(genWZ.second().daughter(1).id());
  int recoZID = abs(recoWZ.second().daughter(0).id()) + abs(recoWZ.second().daughter(1).id());
  
  // check if genWZ and recoWZ IDs are the same
  int genWZID = abs(genWZ.first().daughter(0).id()) + abs(genWZ.second().daughter(0).id());
  int recoWZID = abs(recoWZ.first().daughter(0).id()) + abs(recoWZ.second().daughter(0).id());
  
  theHistograms.fill("GR_ID_genZ_vs_recoZ", "GenZ's and RecoZ's ID", 5, 19, 29, 5, 19, 29, genZID, recoZID);  
  theHistograms.fill("GR_ID_genWZ_vs_recoWZ", "GenWZ and RecoWZ's daughter's ID", 5, 19, 29, 5, 19, 29, genWZID, recoWZID);  
  if(genWZID != recoWZID)
    theHistograms.fill("GR_ID_genWZ_vs_recoWZ_2", "GenWZ and RecoWZ's daughter's ID", 5, 19, 29, 5, 19, 29, genWZID, recoWZID);


  bool WdifferentID = (abs(genWZ.first().daughter(0).id()) != abs(recoWZ.first().daughter(0).id()));
  bool ZdifferentID = (abs(genWZ.second().daughter(0).id()) != abs(recoWZ.second().daughter(0).id()));
  bool WdifferentGR = genWZ.first().daughter(0).id() != recoWZ.first().daughter(0).id();
  bool ZdifferentGR = genZID != recoZID;
  bool WZdifferentGR = genWZID != recoWZID;

  
  theHistograms.fill("GR_ID_differences", "Events with WZ with different IDs", 10, -0.5, 9.5, WdifferentID);  
  theHistograms.fill("GR_ID_differences", "Events with WZ with different IDs", 10, -0.5, 9.5, ZdifferentID +2);
  theHistograms.fill("GR_ID_differences", "Events with WZ with different IDs", 10, -0.5, 9.5, WdifferentGR +4);
  theHistograms.fill("GR_ID_differences", "Events with WZ with different IDs", 10, -0.5, 9.5, ZdifferentGR +6);
  theHistograms.fill("GR_ID_differences", "Events with WZ with different IDs", 10, -0.5, 9.5, WZdifferentGR +8);
  
  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~ End of Reco vs Gen ~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}


void WZAnalyzer::BuildingWZ(){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~ Begin of Build WZ algorithm ~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  int cut = 0;
  theHistograms.fill("GenBuilding_Cuts", "Cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("GenBuilding_Cuts_wei", "Cuts", 20, -0.5, 19.5, cut, theWeight);

  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Jets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ----- filter on jets' number and eta
  if(genJets->size() < 2)
    return;

  Particle Jet0;
  Particle Jet1;
  
  helper_->FindLeadingJets(genJets, Jet0, Jet1, genParticles);

  if(Jet0.p4().Mt() == 0)
    return;

  cut++;
  theHistograms.fill("GenBuilding_Cuts", "Cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("GenBuilding_Cuts_wei", "Cuts", 20, -0.5, 19.5, cut, theWeight);


  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Leptons ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Leptons vectors
  vector<Particle> electron;
  vector<Particle> lepton;
  vector<Particle> muon;
  vector<Particle> neutrino;

  // Leptons variables
  int leptonsID = 0;
  TLorentzVector Ptot;

  for(const Particle &gen: *genParticles){
    if((abs(gen.id()) != 11 && abs(gen.id()) != 13 && abs(gen.id()) != 12 && abs(gen.id()) != 14) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) 
      continue;

    // Leptons accepted are:
    // - electrons with pt >= 7. and |eta| <= 2.5
    if(abs(gen.id()) == 11 && gen.pt() >= 7. && abs(gen.eta()) <= 2.5){
      lepton.push_back(gen);
      leptonsID += abs(gen.id());
      electron.push_back(gen);
    }
    // - muons with pt >= 5. and |eta| <= 2.4    
    else if(abs(gen.id()) == 13 && gen.pt() >= 5. && abs(gen.eta()) <= 2.4){
      lepton.push_back(gen);
      leptonsID += abs(gen.id());
      muon.push_back(gen);
    }
    // - electronic or muonic neutrinos
    else if(abs(gen.id()) == 12 || abs(gen.id()) == 14){
      neutrino.push_back(gen);
    }
  }
  
  // There must be 4 leptons, 3l + 1nu
  if(lepton.size() != 3 || neutrino.size() != 1){
    return;
  }

  cut++;
  theHistograms.fill("GenBuilding_Cuts", "Cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("GenBuilding_Cuts_wei", "Cuts", 20, -0.5, 19.5, cut, theWeight);

  // Additional cuts on leptons' pt
  sort(electron.begin(), electron.end(), PtComparator());
  sort(lepton.begin(),   lepton.end(),   PtComparator());
  sort(muon.begin(),     muon.end(),     PtComparator());

  if(lepton[0].pt() < 20.){
    return;
  }
  
  if(abs(lepton[1].id()) == 11 && lepton[1].pt() < 12){
    return;
  }

  if(abs(lepton[1].id()) == 13 && lepton[1].pt() < 10){
    return;
  }

  cut++;
  theHistograms.fill("GenBuilding_Cuts", "Cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("GenBuilding_Cuts_wei", "Cuts", 20, -0.5, 19.5, cut, theWeight);

  // W and Z must be on shell
  Ptot = neutrino[0].p4();
  for(const Particle lep: lepton){
    Ptot += lep.p4();
  }
  
  theHistograms.fill("GenAll_lllnu_mass",   "m 3 leptons and #nu",     150, 0, 1500, Ptot.M(), theWeight);
  theHistograms.fill("GenAll_lllnu_trmass", "m_{T} 3 leptons and #nu", 150, 0, 1500, Ptot.Mt(), theWeight);

  if(Ptot.M() < 150.){
    return;
  }

  cut++;
  theHistograms.fill("GenBuilding_Cuts", "Cuts", 20, -0.5, 19.5, cut);
  theHistograms.fill("GenBuilding_Cuts_wei", "Cuts", 20, -0.5, 19.5, cut, theWeight);


  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Z & W ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Bosons
  BosonParticle Weh;
  BosonParticle Zet;
  vector<pairBosonParticle> Zls;
  vector<DiBosonParticle> WZs;

  // Boson variables
  bool mixed = false;
  double differenceZ = 999.;
  double differenceW = 999.;

  if(electron.size() == 0 && muon.size() == 3){
    theHistograms.fill("GenBuilding_Topology", "Leptons type", 8, 31.5, 40.5, leptonsID);
    theHistograms.fill("GenBuilding_Topology_wei", "Leptons type", 8, 31.5, 40.5, leptonsID, theWeight);
  }
  else if(electron.size() == 1 && muon.size() == 2){
    theHistograms.fill("GenBuilding_Topology", "Leptons type", 8, 31.5, 40.5, leptonsID);
    theHistograms.fill("GenBuilding_Topology_wei", "Leptons type", 8, 31.5, 40.5, leptonsID, theWeight);
    mixed = true;
  }
  else if(electron.size() == 2 && muon.size() == 1){
    theHistograms.fill("GenBuilding_Topology", "Leptons type", 8, 31.5, 40.5, leptonsID);
    theHistograms.fill("GenBuilding_Topology_wei", "Leptons type", 8, 31.5, 40.5, leptonsID, theWeight);
    mixed = true;
  }
  else if(electron.size() == 3 && muon.size() == 0){
    theHistograms.fill("GenBuilding_Topology", "Leptons type", 8, 31.5, 40.5, leptonsID);
    theHistograms.fill("GenBuilding_Topology_wei", "Leptons type", 8, 31.5, 40.5, leptonsID, theWeight);
  }

  if(mixed){ // ~~~~~ Mixed cases, W and Z identification without caring about ID
    int problems = 0;
    theHistograms.fill("GenBuilding_problems_mixed", "Problems in creating WZ in mixed cases", 30, -0.5, 29.5, problems);
    
    // Creating a Z+l couple
    for(int i = 0; i < (int)lepton.size() -1; i++){
      for(int j = i + 1; j < (int)lepton.size(); j++){
	for(int k = 0; k < (int)lepton.size(); k++){
	  if(k != i && k != j){
	    if(lepton[i].charge() != lepton[j].charge()){
	      Zls.push_back(pairBosonParticle(BosonParticle(lepton[i], lepton[j], 23), lepton[k]));
	    } }	} } }
    
    problems++;
    if(Zls.size() < 1){   
      theHistograms.fill("GenBuilding_problems_mixed", "Problems in creating WZ in mixed cases", 30, -0.5, 29.5, problems);
      return;
    }
    
    // Z is made up of the couple which gives a better Zmass 
    sort(Zls.begin(), Zls.end(), pairMassComparator(0, ZMASS));
    differenceZ = fabs(ZMASS - Zls[0].first.p4().M());

    // check if Z has mixed daughters when Z is always choosen first
    problems++;
    if(abs(Zls[0].first.daughter(0).id()) != abs(Zls[0].first.daughter(1).id())){
      theHistograms.fill("GenBuilding_problems_mixed", "Problems in creating WZ in mixed cases", 30, -0.5, 29.5, problems);
    }

    for(int i = 0; i < (int)Zls.size(); i++){
      WZs.push_back(DiBosonParticle(BosonParticle(neutrino[0], Zls[i].second, copysign(24, Zls[i].second.charge())), Zls[i].first));
    }

    // check if choosen Z&W have mass outside range
    problems++;
    if(abs(WMASS - WZs[0].first().mass()) > 30.){      
      theHistograms.fill("GenBuilding_problems_mixed", "Problems in creating WZ in mixed cases", 30, -0.5, 29.5, problems);
    }
    
    problems++;
    if(abs(ZMASS - WZs[0].second().mass()) > 30.){      
      theHistograms.fill("GenBuilding_problems_mixed", "Problems in creating WZ in mixed cases", 30, -0.5, 29.5, problems);
    }
    
    // W is made up of the couple which gives a better Wmass 
    sort(WZs.begin(), WZs.end(), pairMassComparator(0, WMASS));
    differenceW = fabs(WMASS - WZs[0].first().p4().M());

    problems++;
    if(abs(WZs[0].second().daughter(0).id()) != abs(WZs[0].second().daughter(1).id())){
      theHistograms.fill("GenBuilding_problems_mixed", "Problems in creating WZ in mixed cases", 30, -0.5, 29.5, problems);
    }

    // check if choosen Z&W have mass outside range
    problems++;
    if(abs(WMASS - WZs[0].first().mass()) > 30.){      
      theHistograms.fill("GenBuilding_problems_mixed", "Problems in creating WZ in mixed cases", 30, -0.5, 29.5, problems);
    }
    
    problems++;
    if(abs(ZMASS - WZs[0].second().mass()) > 30.){      
      theHistograms.fill("GenBuilding_problems_mixed", "Problems in creating WZ in mixed cases", 30, -0.5, 29.5, problems);
    }
    
    // Best couple has less difference in mass from main boson
    if(differenceZ < differenceW){ // Z chosen first
      Zet = Zls[0].first;
      Weh = BosonParticle(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()));
      
      problems++;
      theHistograms.fill("GenBuilding_problems_mixed", "Problems in creating WZ in mixed cases", 30, -0.5, 29.5, problems);
      
      problems++;
      if(abs(Zet.daughter(0).id()) != abs(Zet.daughter(1).id())){
	theHistograms.fill("GenBuilding_problems_mixed", "Problems in creating WZ in mixed cases", 30, -0.5, 29.5, problems);
      }
      
      problems++;
      if(abs(WMASS - Weh.mass()) > 30.){      
	theHistograms.fill("GenBuilding_problems_mixed", "Problems in creating WZ in mixed cases", 30, -0.5, 29.5, problems);
      }
      
      problems++;
      if(abs(ZMASS - Zet.mass()) > 30.){      
	theHistograms.fill("GenBuilding_problems_mixed", "Problems in creating WZ in mixed cases", 30, -0.5, 29.5, problems);
      }
      
    } else{ // W chosen first
      problems +=5;
      
      Weh = WZs[0].first();
      Zet = WZs[0].second();

      theHistograms.fill("GenBuilding_problems_mixed", "Problems in creating WZ in mixed cases", 30, -0.5, 29.5, problems);
      
      problems++;
      if(abs(Zet.daughter(0).id()) != abs(Zet.daughter(1).id())){
	theHistograms.fill("GenBuilding_problems_mixed", "Problems in creating WZ in mixed cases", 30, -0.5, 29.5, problems);
      }
      
      problems++;
      if(abs(WMASS - Weh.mass()) > 30.){      
	theHistograms.fill("GenBuilding_problems_mixed", "Problems in creating WZ in mixed cases", 30, -0.5, 29.5, problems);
      }
      
      problems++;
      if(abs(ZMASS - Zet.mass()) > 30.){      
	theHistograms.fill("GenBuilding_problems_mixed", "Problems in creating WZ in mixed cases", 30, -0.5, 29.5, problems);
      }
    }

    if(isTheSame(Zls[0].first, WZs[0].second())){
      theHistograms.fill("GenBuilding_problems_mixed", "Problems in creating WZ in mixed cases", 30, -0.5, 29.5, problems);
    }
  }
  else{ // ~~~~~ Homogeneous cases, W and Z identification
    int problems = 0;
    theHistograms.fill("GenBuilding_problems_homogeneous", "Problems in creating WZ in homogeneous cases", 30, -0.5, 29.5, problems);

    for(int i = 0; i < (int)lepton.size() -1; i++){
      for(int j = i + 1; j < (int)lepton.size(); j++){
	for(int k = 0; k < (int)lepton.size(); k++){
	  if(k != i && k != j){
	    if(lepton[i].charge() != lepton[j].charge()){
	      Zls.push_back(pairBosonParticle(BosonParticle(lepton[i], lepton[j], 23), lepton[k]));
	    } } } } }
    
    problems++;
    if(Zls.size() < 1){   
      theHistograms.fill("GenBuilding_problems_homogeneous", "Problems in creating WZ in homogeneous cases", 30, -0.5, 29.5, problems);
      return;
    }
    
    // Z is made up of the couple which gives a better Zmass 
    sort(Zls.begin(), Zls.end(), pairMassComparator(0, ZMASS));
    differenceZ = fabs(ZMASS - Zls[0].first.p4().M());

    // check if Z has mixed daughters when Z is always choosen first
    problems++;
    if(abs(Zls[0].first.daughter(0).id()) != abs(Zls[0].first.daughter(1).id())){
      theHistograms.fill("GenBuilding_problems_homogeneous", "Problems in creating WZ in homogeneous cases", 30, -0.5, 29.5, problems);
    }

    for(int i = 0; i < (int)Zls.size(); i++){
      WZs.push_back(DiBosonParticle(BosonParticle(Zls[i].second, neutrino[0], copysign(24, Zls[i].second.charge())), Zls[i].first));
    }

    // check if choosen Z&W have mass outside range
    problems++;
    if(abs(WMASS - WZs[0].first().mass()) > 30.){      
      theHistograms.fill("GenBuilding_problems_homogeneous", "Problems in creating WZ in homogeneous cases", 30, -0.5, 29.5, problems);
    }
    
    problems++;
    if(abs(ZMASS - WZs[0].second().mass()) > 30.){      
      theHistograms.fill("GenBuilding_problems_homogeneous", "Problems in creating WZ in homogeneous cases", 30, -0.5, 29.5, problems);
    }
    
    // W is made up of the couple which gives a better Wmass 
    sort(WZs.begin(), WZs.end(), pairMassComparator(0, WMASS));
    differenceW = fabs(WMASS - WZs[0].first().p4().M());

    problems++;
    if(abs(WZs[0].second().daughter(0).id()) != abs(WZs[0].second().daughter(1).id())){
      theHistograms.fill("GenBuilding_problems_homogeneous", "Problems in creating WZ in homogeneous cases", 30, -0.5, 29.5, problems);
    }

    // check if choosen Z&W have mass outside range
    problems++;
    if(abs(WMASS - WZs[0].first().mass()) > 30.){      
      theHistograms.fill("GenBuilding_problems_homogeneous", "Problems in creating WZ in homogeneous cases", 30, -0.5, 29.5, problems);
    }
    
    problems++;
    if(abs(ZMASS - WZs[0].second().mass()) > 30.){      
      theHistograms.fill("GenBuilding_problems_homogeneous", "Problems in creating WZ in homogeneous cases", 30, -0.5, 29.5, problems);
    }
    
    // Best couple has less difference in mass from main boson
    if(differenceZ < differenceW){ // Z chosen first
      Zet = Zls[0].first;
      Weh = BosonParticle(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()));
      
      problems++;
      theHistograms.fill("GenBuilding_problems_homogeneous", "Problems in creating WZ in homogeneous cases", 30, -0.5, 29.5, problems);
      
      problems++;
      if(abs(Zet.daughter(0).id()) != abs(Zet.daughter(1).id())){
	theHistograms.fill("GenBuilding_problems_homogeneous", "Problems in creating WZ in homogeneous cases", 30, -0.5, 29.5, problems);
      }
      
      problems++;
      if(abs(WMASS - Weh.mass()) > 30.){      
	theHistograms.fill("GenBuilding_problems_homogeneous", "Problems in creating WZ in homogeneous cases", 30, -0.5, 29.5, problems);
      }
      
      problems++;
      if(abs(ZMASS - Zet.mass()) > 30.){      
	theHistograms.fill("GenBuilding_problems_homogeneous", "Problems in creating WZ in homogeneous cases", 30, -0.5, 29.5, problems);
      }
      
    } else{ // W chosen first
      problems +=5;
      
      Weh = WZs[0].first();
      Zet = WZs[0].second();

      theHistograms.fill("GenBuilding_problems_homogeneous", "Problems in creating WZ in homogeneous cases", 30, -0.5, 29.5, problems);
      
      problems++;
      if(abs(Zet.daughter(0).id()) != abs(Zet.daughter(1).id())){
	theHistograms.fill("GenBuilding_problems_homogeneous", "Problems in creating WZ in homogeneous cases", 30, -0.5, 29.5, problems);
      }
      
      problems++;
      if(abs(WMASS - Weh.mass()) > 30.){      
	theHistograms.fill("GenBuilding_problems_homogeneous", "Problems in creating WZ in homogeneous cases", 30, -0.5, 29.5, problems);
      }
      
      problems++;
      if(abs(ZMASS - Zet.mass()) > 30.){      
	theHistograms.fill("GenBuilding_problems_homogeneous", "Problems in creating WZ in homogeneous cases", 30, -0.5, 29.5, problems);
      }
    }

    if(isTheSame(Zls[0].first, WZs[0].second())){
      theHistograms.fill("GenBuilding_problems_homogeneous", "Problems in creating WZ in homogeneous cases", 30, -0.5, 29.5, problems);
    }
  }

  
    
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~ End of Build WZ algorithm ~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}


void WZAnalyzer::begin(){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~ Begin function ~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout << "\n--------------------------------------------------------------------------" << endl;
  cout << "\n Starting WZAnalysis on sample in " << fileName << endl;
  
  auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  cout << " Beginnin analysis at: " << ctime(&now) << endl;

  // events counters
  eventSample = 0;
  eventGen = 0;
  eventReco = 0;
  eventGenReco = 0;
  eventGenaftercut = 0;
  eventRecoaftercut = 0;
  eventGenNOReco = 0;
  eventRecoNOGen = 0;

  // weight counters
  weightGen = 0;
  weightReco = 0;
  weightGenReco = 0;
  weightGenaftercut = 0;
  weightRecoaftercut = 0;
  weightGenNOReco = 0;
  weightRecoNOGen = 0;

  // free counters
  
  // time begins
  begintime = ((float)clock())/CLOCKS_PER_SEC;
}


void WZAnalyzer::analyze(){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~ Main analysis ~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  eventSample++;
  theHistograms.fill("WZ_Events", "Weighted counters", 10, -0.5, 9.5, 1, theWeight);

  if(eventSample % 10000 == 0)
    cout << "Event: " << eventSample << endl;

  //gen variables	
  DiBosonParticle genWZ;
  Particle genJet0;
  Particle genJet1;
  
  //reco variables
  DiBosonLepton recoWZ;
  Jet recoJet0;
  Jet recoJet1;

  //Check WZ construction
  WZAnalyzer::BuildingWZ();
  
  //Gen analysis
  WZAnalyzer::GenAnalysis(genWZ, genJet0, genJet1);
  
  //Reco analysis
  WZAnalyzer::RecoAnalysis(recoWZ, recoJet0, recoJet1);

  ///*
  if(genWZ.pt() != 0. && recoWZ.pt() != 0.){
    //Reco vs Gen analysis
    WZAnalyzer::GenRecoAnalysis(genWZ, genJet0, genJet1, recoWZ, recoJet0, recoJet1);
  }

  if(genWZ.pt() != 0. && recoWZ.pt() == 0.){
    //Gen event not reconstructed
    eventGenNOReco++;
    theHistograms.fill("WZ_Events", "Weighted counters", 10, -0.5, 9.5, 7, theWeight);
  }

  if(genWZ.pt() == 0. && recoWZ.pt() != 0.){
    //Reco event not generated
    eventRecoNOGen++;
    theHistograms.fill("WZ_Events", "Weighted counters", 10, -0.5, 9.5, 8, theWeight);
  }  
  //*/
}


Int_t WZAnalyzer::cut() {
  return 1;
}


void WZAnalyzer::end(TFile &){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~ End banner ~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout << "\n--------------------------------------------------------------------------" << endl;
  
  cout << "\nEvents of the sample analyzed:                       " << setw(9) << eventSample << endl;
  cout << "Gen events analyzed:                                 " << setw(9) << eventGen << "\t" << eventGen*100./eventSample << "%" << "\t" << eventGen*weightGen*1.0 << endl;
  cout << "Gen events after all cuts:                           " << setw(9) << eventGenaftercut << "\t" << eventGenaftercut*100./eventGen << "%" << "\t" << eventGenaftercut*weightGenaftercut*1.0 << endl;
  cout << "Reco events analyzed:                                " << setw(9) << eventReco << "\t" << eventReco*100./eventSample << "%" << "\t" << eventReco*weightReco*1.0 << endl;
  cout << "Reco events after all cuts:                          " << setw(9) << eventRecoaftercut << "\t" << eventRecoaftercut*100./eventReco << "%" << "\t" << eventRecoaftercut*weightRecoaftercut*1.0 << endl;
  cout << "Gen&Reco events analyzed:                            " << setw(9) << eventGenReco << "\t" << eventGenReco*100./eventSample << "%" << "\t" << eventGenReco*weightGenReco*1.0 << endl;
  cout << "Gen events not reconstructed:                        " << setw(9) << eventGenNOReco << "\t" << eventGenNOReco*100./eventGen << "%" << "\t" << eventGenNOReco*weightGenNOReco*1.0 << endl;
  cout << "Reco events that weren't gen events:                 " << setw(9) << eventRecoNOGen << "\t" << eventRecoNOGen*100./eventReco << "%"<< "\t" << eventRecoNOGen*weightRecoNOGen*1.0  << endl;
    
  // execution time  
  auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  cout << "\n Finished analysis at: " << ctime(&now) << endl;
  
  endtime = ((float)clock())/CLOCKS_PER_SEC;
  helper_->printTime(begintime, endtime);
  cout << "\n--------------------------------------------------------------------------" << endl;
}
