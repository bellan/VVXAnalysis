#include "VVXAnalysis/TreeAnalysis/interface/VZGAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"

#include "VVXAnalysis/Commons/interface/GenVBHelper.h"


#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// #include <TString.h>

#include <boost/assign/std/vector.hpp>
using namespace boost::assign;

using std::cout;
using std::endl;

using namespace phys;


double etacut=4.7;

double ptcut=30;

bool KinematicsOK(phys::Particle p, float pt,float eta)
{
       if (fabs(p.eta()) < eta && fabs(p.pt()) > pt)
       {
              return true;
       }
       else
              return false;
}


void VZGAnalyzer::PlotJets(const phys::Particle &Jet0, const phys::Particle &Jet1, std::string prename, const float weight, std::string suffix)
{
       std::string name = "J0";
       theHistograms->fill(prename + name + "_charge_" + suffix, prename + name + "'s charge", 5, -2.5, 2.5, Jet0.charge(), weight);
       theHistograms->fill(prename + name + "_mass_" + suffix, prename + name + "'s mass", 63, 0, 252, Jet0.mass(), weight);
       theHistograms->fill(prename + name + "_trmass_" + suffix, prename + name + "'s trmass", 50, 0, 400, Jet0.p4().Mt(), weight);
       theHistograms->fill(prename + name + "_pt_" + suffix, prename + name + "'s p_{t}", 50, ptcut, 600, Jet0.pt(), weight);
       theHistograms->fill(prename + name + "_Y_" + suffix, prename + name + "'s Y", 50, -5, 5, Jet0.rapidity(), weight);
       theHistograms->fill(prename + name + "_eta_" + suffix, prename + name + "'s #eta", 50, -9, 9, Jet0.eta(), weight);
       theHistograms->fill(prename + name + "_phi_" + suffix, prename + name + "'s #phi", 50, -3.5, 3.5, Jet0.phi(), weight);

       name = "J1";
       theHistograms->fill(prename + name + "_charge_" + suffix, prename + name + "'s charge", 5, -2.5, 2.5, Jet1.charge(), weight);
       theHistograms->fill(prename + name + "_mass_" + suffix, prename + name + "'s mass", 63, 0, 252, Jet1.mass(), weight);
       theHistograms->fill(prename + name + "_trmass_" + suffix, prename + name + "'s trmass", 50, 0, 400, Jet1.p4().Mt(), weight);
       theHistograms->fill(prename + name + "_pt_" + suffix, prename + name + "'s p_{t}", 50, 30, 600, Jet1.pt(), weight);
       theHistograms->fill(prename + name + "_Y_" + suffix, prename + name + "'s Y", 50, -5, 5, Jet1.rapidity(), weight);
       theHistograms->fill(prename + name + "_eta_" + suffix, prename + name + "'s #eta", 50, -9, 9, Jet1.eta(), weight);
       theHistograms->fill(prename + name + "_phi_" + suffix, prename + name + "'s #phi", 50, -3.5, 3.5, Jet1.phi(), weight);

       name = "JJ";
       TLorentzVector JJp4 = Jet0.p4() + Jet1.p4();
       double JJdeltaEta = Jet0.eta() - Jet1.eta();
       double JJdeltaPhi = physmath::deltaPhi(Jet0.phi(), Jet1.phi());
       double JJdeltaR = abs(physmath::deltaR(Jet0, Jet1));

       theHistograms->fill(prename + name + "_mass_" + suffix, " Jets' mass", 10, 50, 120, JJp4.M(), weight);
       theHistograms->fill(prename + name + "_trmass_" + suffix, " Jets' trmass", 50, 0, 200, JJp4.Mt(), weight);
       theHistograms->fill(prename + name + "_pt_" + suffix, " Jets' p_{t}", 50, 0, 600, JJp4.Pt(), weight);
       theHistograms->fill(prename + name + "_deltaEta_" + suffix, " Jets' #Delta#eta", 50, -9, 9, JJdeltaEta, weight);
       theHistograms->fill(prename + name + "_deltaEtaabs_" + suffix, " Jets' |#Delta#eta|", 25, 0, 9, abs(JJdeltaEta), weight);
       theHistograms->fill(prename + name + "_deltaR_" + suffix, " Jets' #DeltaR", 25, -0.5, 9, JJdeltaR, weight);
       theHistograms->fill(prename + name + "_deltaPhi_" + suffix, " Jets' #Delta#phi", 50, -3.5, 3.5, JJdeltaPhi, weight);

       theHistograms->fill(prename + name + "_massvsdeltaEta_" + suffix, prename + name + "'s mass(x) vs #Delta#eta(y)", 12, 160, 1780, 10, -6.5, 6.5, JJp4.M(), JJdeltaEta, weight);
       theHistograms->fill(prename + name + "_massvsdeltaEtaabs_" + suffix, prename + name + "'s mass(x) vs |#Delta#eta|(y)", 12, 160, 1780, 10, -6.5, 6.5, JJp4.M(), abs(JJdeltaEta), weight);
}

void VZGAnalyzer::PlotJet(const phys::Particle &Jet, std::string prename, const float weight, std::string suffix)
{
       std::string name = " ";
       theHistograms->fill(prename + name + "_charge_" + suffix, prename + name + "'s charge", 5, -2.5, 2.5, Jet.charge(), weight);
       theHistograms->fill(prename + name + "_mass_" + suffix, prename + name + "'s mass", 63, 50, 120, Jet.mass(), weight);
       theHistograms->fill(prename + name + "_trmass_" + suffix, prename + name + "'s trmass", 50, 0, 200, Jet.p4().Mt(), weight);
       theHistograms->fill(prename + name + "_pt_" + suffix, prename + name + "'s p_{t}", 50, ptcut, 600, Jet.pt(), weight);
       theHistograms->fill(prename + name + "_Y_" + suffix, prename + name + "'s Y", 50, -5, 5, Jet.rapidity(), weight);
       theHistograms->fill(prename + name + "_eta_" + suffix, prename + name + "'s #eta", 50, -9, 9, Jet.eta(), weight);
       theHistograms->fill(prename + name + "_phi_" + suffix, prename + name + "'s #phi", 50, -3.5, 3.5, Jet.phi(), weight);
}

bool VZGAnalyzer::LeptonicSignalCostraint()
{
  //int LeptonicZdecays=genVBHelper_.ZtoChLep().size();//+genVBHelper_.ZtoNeutrinos().size();
       if (genVBHelper_.ZtoChLep().size()==1 && KinematicsOK(genVBHelper_.ZtoChLep()[0].daughter(0), 5.,2.5) && KinematicsOK(genVBHelper_.ZtoChLep()[0].daughter(1), 5.,2.5) && genVBHelper_.ZtoChLep()[0].mass()>60 && genVBHelper_.ZtoChLep()[0].mass()<120)
       {
              return true;
       }
       else
       {
              return false;
       }
}
bool VZGAnalyzer::HadronicSignalCostraint()
{
       if ((genVBHelper_.ZtoQ().size()==1 && KinematicsOK(genVBHelper_.ZtoQ()[0].daughter(0), ptcut,etacut) && KinematicsOK(genVBHelper_.ZtoQ()[0].daughter(1), ptcut,etacut)&& genVBHelper_.ZtoQ()[0].mass()>50 && genVBHelper_.ZtoQ()[0].mass()<120)||(genVBHelper_.WtoQ().size()==1 && KinematicsOK(genVBHelper_.WtoQ()[0].daughter(0), ptcut,etacut)&& KinematicsOK(genVBHelper_.WtoQ()[0].daughter(1), ptcut,etacut)&& genVBHelper_.WtoQ()[0].mass()>50 && genVBHelper_.WtoQ()[0].mass()<120))
       {
              return true;
       }
       else
       {
              return false;
       }
}

bool VZGAnalyzer::PhotonSignalCostraint()
{
       std::vector<phys::Particle> selectedphotons;
       for (auto p : *genParticles)
       {
              if (p.id() == 22 && KinematicsOK(p, 20, 2.4) && p.genStatusFlags().test(phys::isPrompt) &&  p.genStatusFlags().test(phys::fromHardProcess))
              {
                  selectedphotons.push_back(p);
              }
       }
       std::cout<< "Number of selected gen photons = "<<selectedphotons.size()<<std::endl;
       if (selectedphotons.size()>=1)
       {
              return true;
       }
       else
       {
              return false;
       }
}

bool VZGAnalyzer::RECOsignalCostraint()
{
       //----------------------------------------Building jj pairs ----------------------------------------//

       std::vector<phys::Jet> selectedRECOjets;
       std::vector<phys::Boson<phys::Jet>> DiJets;
       foreach (const phys::Jet &jet, *jets)
       {
              if (KinematicsOK(jet, ptcut, etacut)) // KinematicsOK(jet)
              {
                     selectedRECOjets.push_back(jet);
              }
       }
       for (int i = 0; i < selectedRECOjets.size(); i++)
       {
              phys::Jet jetA = selectedRECOjets[i];
              for (int j = i + 1; j < selectedRECOjets.size(); j++)
              {
                     phys::Jet jetB = selectedRECOjets[j];
                     float mjj = (jetA.p4() + jetB.p4()).M();
                     std::cout << "mjj= " << mjj << std::endl;
                     if (mjj > 50 && mjj < 120) 
                     {
                       DiJets.push_back(phys::Boson<phys::Jet>(jetA, jetB));
                     }
              }
       }

      std::cout << "DiJets size: " << DiJets.size() << std::endl;

       std::vector<phys::Photon> selectedphotons;

       for (auto p : *photons)
       {
              // Pixel seed and electron veto
              // if (ph.hasPixelSeed() || !ph.passElectronVeto())
              //        continue;

              if (p.id() == 22 && KinematicsOK(p, 20, 2.4) && !p.hasPixelSeed() && p.passElectronVeto())
              {
                     selectedphotons.push_back(p);
              }
       }
       std::cout<< "Number of selected RECO photons = "<<selectedphotons.size()<<std::endl;

       bool goodZ= (Z->mass() > 60 && Z->mass() < 120 && KinematicsOK(Z->daughter(0), 5, 2.5) && KinematicsOK(Z->daughter(1), 5, 2.5)); // KinematicsOK(jet)

       
       if (DiJets.size() > 0 && selectedphotons.size() >0 && goodZ )
       {
              return true;
       }
       else
       {
              return false;
       }
}


Bool_t VZGAnalyzer::cut(Int_t n, phys::Boson<phys::Jet> recoV, std::vector<phys::Photon> selectedphotons)
{ // returns false if the event has to be cut

  switch (n)
  {

  case 1:
    if ((jets->size() > 1 && recoV.mass() > 50 && recoV.mass() < 120  && KinematicsOK(recoV.daughter(0), ptcut,etacut) && KinematicsOK(recoV.daughter(1), ptcut,etacut) ) && (Z->mass() > 50 && Z->mass() < 120  && KinematicsOK(Z->daughter(0), 5.,etacut) && KinematicsOK(Z->daughter(1), 5.,etacut) )&& selectedphotons.size()>0)
      return true;
    break;

  // case 2:
  //   if (recoV.mass() > 65)
  //     return true;

  //   break;

  // case 3:
  //   if (recoV.mass() < 105)
  //     return true;

  //   break;

  default:
    return true;
  }
  return false;
}

void VZGAnalyzer::analyze()
{ // It's the only member function running each event.

  cout << "----------------------------------------------------------------" << endl;
  cout << "Run: " << run << " event: " << event << endl;

  if (LeptonicSignalCostraint()  && HadronicSignalCostraint() && PhotonSignalCostraint() )
  {
    theHistograms->fill("Signal_fraction_QLG", "Signal_fraction_QLG", 2, 0, 2, 1., theWeight);
  }
  else
  {
    theHistograms->fill("Signal_fraction_QLG", "Signal_fraction_QLG", 2, 0, 2, 0., theWeight);
  }

  if (HadronicSignalCostraint() )
  {
    theHistograms->fill("Signal_fraction_Q", "Signal_fraction_Q", 2, 0, 2, 1., theWeight);
  }
  else
  {
    theHistograms->fill("Signal_fraction_Q", "Signal_fraction_Q", 2, 0, 2, 0., theWeight);
  }

  if (LeptonicSignalCostraint() )
  {
    theHistograms->fill("Signal_fraction_L", "Signal_fraction_L", 2, 0, 2, 1., theWeight);
  }
  else
  {
    theHistograms->fill("Signal_fraction_L", "Signal_fraction_L", 2, 0, 2, 0., theWeight);
  }

  if (PhotonSignalCostraint() )
  {
    theHistograms->fill("Signal_fraction_G", "Signal_fraction_G", 2, 0, 2, 1., theWeight);
  }
  else
  {
    theHistograms->fill("Signal_fraction_G", "Signal_fraction_G", 2, 0, 2, 0., theWeight);
  }
  //___________________________________________________________________________________

  bool GENsignal= (HadronicSignalCostraint() && LeptonicSignalCostraint());
  bool RECOsignal=RECOsignalCostraint();
  if (RECOsignal)
  {
    theHistograms->fill("Signal_fraction_RECO", "Signal_fraction_RECO", 2, 0, 2, 1., theWeight);
  }
  else
  {
    theHistograms->fill("Signal_fraction_RECO", "Signal_fraction_RECO", 2, 0, 2, 0., theWeight);
  }
  //_____________________________________________________________________________________________//

  if (GENsignal && RECOsignal)
  {
    theHistograms->fill("GENRECO_11", "GENRECO_11", 2, 0, 2, 1., theWeight);
  }
  else
  {
    theHistograms->fill("GENRECO_11", "GENRECO_11", 2, 0, 2, 0., theWeight);
  }

  if (!GENsignal && RECOsignal)
  {
    theHistograms->fill("GENRECO_01", "GENRECO_01", 2, 0, 2, 1., theWeight);
  }
  else
  {
    theHistograms->fill("GENRECO_01", "GENRECO_01", 2, 0, 2, 0., theWeight);
  }

  if (GENsignal && !RECOsignal)
  {
    theHistograms->fill("GENRECO_10", "GENRECO_10", 2, 0, 2, 1., theWeight);
  }
  else
  {
    theHistograms->fill("GENRECO_10", "GENRECO_10", 2, 0, 2, 0., theWeight);
  }

  if (!GENsignal && !RECOsignal)
  {
    theHistograms->fill("GENRECO_00", "GENRECO_00", 2, 0, 2, 1., theWeight);
  }
  else
  {
    theHistograms->fill("GENRECO_00", "GENRECO_00", 2, 0, 2, 0., theWeight);
  }


//___________________________________________________________________________________
  // genVBAnalyzer();
  //genAnalyze();
    PhotonvsJet();

  phys::Boson<phys::Jet> recoV;
  Reconstruct(&recoV);
  std::vector<phys::Photon> selectedphotons;
  PhotonSelection(&selectedphotons);

  // std::vector<phys::Jet> selectedRECOjets;
  // foreach (const phys::Jet &jet, *jets)
  // {
  //   if (KinematicsOK(jet, ptcut, etacut)) // KinematicsOK(jet)
  //   {
  //     selectedRECOjets.push_back(jet);
  //   }
  // }

  // phys::Photon mostEnergeticPhoton;
  //  if (selectedphotons.size() > 0)
  // {
  //   std::stable_sort(selectedphotons.begin(), selectedphotons.end(), phys::EComparator());
  //   mostEnergeticPhoton = selectedphotons[0];
  // }

  // theHistograms->fill("pt_mostenergeticphoton", "pt_mostenergeticphoton", 30, 0, 300, mostEnergeticPhoton.pt(), (theWeight));
  // std::vector<std::pair<phys::Photon, phys::Jet>> nearestRECOjetstoPhoton;

  // if (selectedRECOjets.size() > 0)
  // {
  //   std::stable_sort(selectedRECOjets.begin(), selectedRECOjets.end(), phys::DeltaRComparator(mostEnergeticPhoton));
  //   phys::Jet nearestRECOjet = selectedRECOjets.at(0);
  //   nearestRECOjetstoPhoton.push_back({mostEnergeticPhoton, nearestRECOjet});
  // }
  // for (auto pair : nearestRECOjetstoPhoton)
  // {
  //   ResolutionPlots(pair.first, pair.second, "Photon_vs_jet_", theWeight);
  //   theHistograms->fill("DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet", "DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet; #DeltaR", 50, 0, 5, abs(physmath::deltaR(pair.first, pair.second)), theWeight);
  //   theHistograms->fill("DeltaR_vs_Deltapt", "DeltaR_vs_Deltapt;#Delta pt [GeV/c] ; #DeltaR", 20, -100, 100, 20, 0, 2, pair.first.pt()-pair.second.pt(),abs(physmath::deltaR(pair.first, pair.second)), theWeight);
  // }


  std::vector<phys::Jet> selectedRECOjets;
  foreach (const phys::Jet &jet, *jets)
  {
    if (KinematicsOK(jet, ptcut, etacut)) // KinematicsOK(jet)
    {
      selectedRECOjets.push_back(jet);
    }
  }

  phys::Photon mostEnergeticPhoton;
   if (selectedphotons.size() > 0)
  {
    std::stable_sort(selectedphotons.begin(), selectedphotons.end(), phys::EComparator());
    mostEnergeticPhoton = selectedphotons[0];
  }

  theHistograms->fill("pt_mostenergeticphoton", "pt_mostenergeticphoton", 30, 0, 300, mostEnergeticPhoton.pt(), (theWeight));
  std::vector<std::pair<phys::Photon, phys::Jet>> nearestRECOjetstoPhoton;

  if (selectedRECOjets.size() > 0)
  {
    std::stable_sort(selectedRECOjets.begin(), selectedRECOjets.end(), phys::DeltaRComparator(mostEnergeticPhoton));
    phys::Jet nearestRECOjet = selectedRECOjets.at(0);
    nearestRECOjetstoPhoton.push_back({mostEnergeticPhoton, nearestRECOjet});
  }
  for (auto pair : nearestRECOjetstoPhoton)
  {
    //ResolutionPlots(pair.first, pair.second, "Photon_vs_jet_", theWeight);
    theHistograms->fill("DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet", "DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet; #DeltaR", 50, 0, 5, abs(physmath::deltaR(pair.first, pair.second)), theWeight);
    theHistograms->fill("DeltaR_vs_Deltapt", "DeltaR_vs_Deltapt;#Delta pt [GeV/c] ; #DeltaR", 20, -100, 100, 20, 0, 2, pair.first.pt()-pair.second.pt(),abs(physmath::deltaR(pair.first, pair.second)), theWeight);
  }



  //---------------ALL-THE-EVENTS---------------//

 
  printHistos(0, "all", recoV,selectedphotons);

  //---------------SIGNAL-EVENTS---------------//

  if (HadronicSignalCostraint() == 1 && LeptonicSignalCostraint() == 1 && PhotonSignalCostraint() == 1)
  {
    // std::vector phys::Photon selectedphotons;
    // for (auto p : *photons)
    // {
    //   if (p.id() == 22 && KinematicsOK(p, 20, 2.4) && !p.hasPixelSeed() && p.passElectronVeto())
    //   {
    //     selectedphotons->push_back(p);
    //   }
    // }



    theHistograms->fill("Signal_V_vs_Z_pt", "Signal_V_vs_Z_pt;V p_{t} [GeV/c]; #Z p_{t} [GeV/c]", 30, 0, 300, 30, 0, 300, recoV.pt(), Z->pt(), theWeight);
    printHistos(0, "sign", recoV, selectedphotons);
  } // not signalCostraint anymore

  //---------------BACKGROUND-EVENTS---------------//

  else
  {
    printHistos(0, "bckg", recoV,selectedphotons);
  }
}


void VZGAnalyzer::genPhotonsAnalyzer()
{}


void VZGAnalyzer::genVBAnalyzer()

{
  theHistograms->fill("nZtoChLep"    ,"Number of Z->ll per event" , 7,0,7, genVBHelper_.ZtoChLep().size());
  theHistograms->fill("nZtoNeutrinos","Number of Z->nn per event" , 7,0,7, genVBHelper_.ZtoNeutrinos().size());
  theHistograms->fill("nWtoLep"      ,"Number of W->lnu per event", 7,0,7, genVBHelper_.WtoLep().size());
  theHistograms->fill("nZtoQ"        ,"Number of Z->qq per event" , 7,0,7, genVBHelper_.ZtoQ().size());
  theHistograms->fill("nWtoQ"        ,"Number of W->qq' per ev""ent", 7,0,7, genVBHelper_.WtoQ().size());

  int nVBs = genVBHelper_.ZtoChLep().size() + genVBHelper_.ZtoNeutrinos().size() + genVBHelper_.WtoLep().size() + genVBHelper_.ZtoQ().size() + genVBHelper_.WtoQ().size();
  theHistograms->fill("nVBs", "Number of VB per event", 7,0,7, nVBs);
  for (auto VB : genVBHelper_.ZtoQ())
  {
    theHistograms->fill("ZtoQ mass", "ZtoQ mass", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass", "VB mass", 10, 50, 120, VB.mass());
  }
  for (auto VB : genVBHelper_.WtoQ())
  {
    theHistograms->fill("WtoQ mass", "WtoQ mass", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass", "VB mass", 10, 50, 120, VB.mass());
  }
  for (auto VB : genVBHelper_.ZtoChLep())
  {
    theHistograms->fill("ZtoChLep mass", "ZtoChLep mass", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass", "VB mass", 10, 50, 120, VB.mass());
  }
  for (auto VB : genVBHelper_.WtoLep())
  {
    theHistograms->fill("WtoLep mass", "WtoLep mass", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass", "VB mass", 10, 50, 120, VB.mass());
  }
  for (auto VB : genVBHelper_.ZtoNeutrinos())
  {
    theHistograms->fill("ZtoNeutrinos mass", "ZtoNeutrinos mass", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass", "VB mass", 10, 50, 120, VB.mass());
  }


  theHistograms->fill("nZtoChLep Alternative"    ,"Number of Z->ll per event Alternative" , 7,0,7, genZlepCandidates_->size());
  theHistograms->fill("nWtoLep Alternative"      ,"Number of W->lnu per event Alternative", 7,0,7, genWlepCandidates_->size());
  theHistograms->fill("nZtoQ Alternative"        ,"Number of Z->qq per event Alternative" , 7,0,7, genZhadCandidates_->size());
  theHistograms->fill("nWtoQ Alternative"        ,"Number of W->qq' per event Alternative", 7,0,7, genWlepCandidates_->size());

  int nVBs_alternative = genZlepCandidates_->size() + genWlepCandidates_->size() + genZhadCandidates_->size() + genWhadCandidates_->size();
  theHistograms->fill("nVBs Alternative", "Number of VB per event Alternative", 7,0,7, nVBs_alternative);
  foreach (auto & VB , *genZhadCandidates_)
  {
    theHistograms->fill("ZtoQ mass Alternative", "ZtoQ mass Alternative", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass Alternative", "VB mass Alternative", 10, 50, 120, VB.mass());
  }
  foreach (auto & VB , *genWhadCandidates_)
  {
    theHistograms->fill("WtoQ mass Alternative", "WtoQ mass Alternative", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass Alternative", "VB mass Alternative", 10, 50, 120, VB.mass());
  }
  foreach (auto & VB , *genZlepCandidates_)
  {
    theHistograms->fill("ZtoChLep mass Alternative", "ZtoChLep mass Alternative", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass Alternative", "VB mass Alternative", 10, 50, 120, VB.mass());
  }
  foreach (auto & VB , *genWlepCandidates_)
  {
    theHistograms->fill("WtoLep mass Alternative", "WtoLep mass Alternative", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass Alternative", "VB mass Alternative", 10, 50, 120, VB.mass());
  }
  // for (auto VB : genVBHelper_.ZtoNeutrinos())
  // {
  //   theHistograms->fill("ZtoNeutrinos mass Alternative", "ZtoNeutrinos mass Alternative", 10, 50, 120, VB.mass());
  //   theHistograms->fill("VB mass Alternative", "VB mass Alternative", 10, 50, 120, VB.mass());
  // }
}


void VZGAnalyzer::Reconstruct(phys::Boson<phys::Jet> *V_JJCandidate)
{
  // Building of every jets pairs combination
  std::vector<phys::Boson<phys::Jet>> DiJets;

  for (uint i = 0; i < jets->size(); i++) // Warning: size can be 0
    for (uint j = i+1; j < jets->size(); j++)
        DiJets.push_back(phys::Boson<phys::Jet>(jets->at(i), jets->at(j)));

  if (jets->size() > 1)
  {
    std::stable_sort(DiJets.begin(), DiJets.end(), phys::Mass2Comparator(phys::ZMASS, phys::WMASS));
    *V_JJCandidate = DiJets.at(0);
  }
}

void VZGAnalyzer::PhotonSelection(std::vector<phys::Photon> *phot)
{
  for (auto p : *photons)
  {
    // Pixel seed and electron veto
    // if (ph.hasPixelSeed() || !ph.passElectronVeto())
    //        continue;

    if (p.id() == 22 && KinematicsOK(p, 20, 2.4) && !p.hasPixelSeed() && p.passElectronVeto())
    {
        phot->push_back(p);
    }
  }
  std::cout << "Number of selected RECO photons = " << phot->size() << std::endl;
}


void VZGAnalyzer::CompatibilityTest(phys::Boson<phys::Jet> bestCandidate, phys::Boson<phys::Particle> genVB, std::string sample, std::string algorithm)
{
}

void VZGAnalyzer::printHistos(uint i, std::string histoType, phys::Boson<phys::Jet> recoV, std::vector<phys::Photon> selectedphotons)
{

  std::vector<std::string> cuts = {"0", "1"};//, "2", "3"};

  if (i < cuts.size() && cut(i, recoV,selectedphotons))
  {

    theHistograms->fill("recoVMass_" + histoType + cuts.at(i), "mass of recoV", 10, 0, 200, recoV.mass(), (theWeight));
    //theHistograms->fill("recoVTot_" + histoType + cuts.at(i), "mass of recoV", 1, 50, 120, recoV.mass(), (theWeight));

    theHistograms->fill("recoVDaughter0Pt_" + histoType + cuts.at(i), "pt of recoVDaughter0", 50, 0, 600, recoV.daughter(0).pt(), (theWeight));
    theHistograms->fill("recoVDaughter1Pt_" + histoType + cuts.at(i), "pt of recoVDaughter1", 50, 0, 600, recoV.daughter(1).pt(), (theWeight));
    //PlotJets(recoV.daughter(0),recoV.daughter(1), "", theWeight, histoType + cuts.at(i));

    theHistograms->fill("recoZMass_" + histoType + cuts.at(i), "mass of recoZ", 10, 0, 200, Z->mass(), (theWeight));

    theHistograms->fill("recoVPt_" + histoType + cuts.at(i), "pt of recoV", 50, 0, 300, recoV.pt(), (theWeight));
    theHistograms->fill("recoVEta_" + histoType + cuts.at(i), "eta of recoV", 30, 0, 3.5, fabs(recoV.eta()), (theWeight));
    theHistograms->fill("recoVPhi_" + histoType + cuts.at(i), "phi of recoV", 30, -3.2, 3.2, recoV.phi(), (theWeight));
    theHistograms->fill("recoVEnergy_" + histoType + cuts.at(i), "energy of  recoV", 60, 0, 2400, fabs(recoV.e()), (theWeight));
    theHistograms->fill("recoVDaughtersDeltaPhi_" + histoType + cuts.at(i), "dPhi of recoVDaughters", 30, 0, 3.2, fabs(physmath::deltaPhi(recoV.daughter(0).phi(), recoV.daughter(1).phi())), (theWeight));

    theHistograms->fill("recoZPt_" + histoType + cuts.at(i), "pt of recoZ", 50, 0, 600, Z->pt(), (theWeight));
    theHistograms->fill("recoZEta_" + histoType + cuts.at(i), "eta of recoZ", 70, 0, 3.5, fabs(Z->eta()), (theWeight));
    theHistograms->fill("recoZEnergy_" + histoType + cuts.at(i), "energy of  recoZ", 120, 0, 400, fabs(Z->e()), (theWeight));
    theHistograms->fill("recoZDeltaPhi_" + histoType + cuts.at(i), "dPhi of recoZ", 30, 0, 3.2, fabs(physmath::deltaPhi(Z->daughter(0).phi(), Z->daughter(0).phi())), (theWeight));




    
  // std::vector<phys::Jet> selectedRECOjets;
  // foreach (const phys::Jet &jet, *jets)
  // {
  //   if (KinematicsOK(jet, ptcut, etacut)) // KinematicsOK(jet)
  //   {
  //     selectedRECOjets.push_back(jet);
  //   }
  // }

  // phys::Photon mostEnergeticPhoton;
  //  if (selectedphotons.size() > 0)
  // {
  //   std::stable_sort(selectedphotons.begin(), selectedphotons.end(), phys::EComparator());
  //   mostEnergeticPhoton = selectedphotons[0];
  // }

  // theHistograms->fill("pt_mostenergeticphoton"+histoType + cuts.at(i), "pt_mostenergeticphoton"+histoType + cuts.at(i), 30, 0, 300, mostEnergeticPhoton.pt(), (theWeight));
  // std::vector<std::pair<phys::Photon, phys::Jet>> nearestRECOjetstoPhoton;

  // if (selectedRECOjets.size() > 0)
  // {
  //   std::stable_sort(selectedRECOjets.begin(), selectedRECOjets.end(), phys::DeltaRComparator(mostEnergeticPhoton));
  //   phys::Jet nearestRECOjet = selectedRECOjets.at(0);
  //   nearestRECOjetstoPhoton.push_back({mostEnergeticPhoton, nearestRECOjet});
  // }
  // for (auto pair : nearestRECOjetstoPhoton)
  // {
  //   ResolutionPlots(pair.first, pair.second, "Photon_vs_jet_", theWeight,"_"+ histoType + cuts.at(i));
  //   theHistograms->fill("DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet_"+histoType + cuts.at(i), "DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet"+histoType + cuts.at(i)+"; #DeltaR", 50, 0, 5, abs(physmath::deltaR(pair.first, pair.second)), theWeight);
  //   theHistograms->fill("DeltaR_vs_Deltapt_"+histoType + cuts.at(i), "DeltaR_vs_Deltapt"+histoType + cuts.at(i)+";#Delta pt [GeV/c] ; #DeltaR", 20, -100, 100, 20, 0, 2, pair.first.pt()-pair.second.pt(),abs(physmath::deltaR(pair.first, pair.second)), theWeight);
  // }
  








    printHistos(++i, histoType, recoV, selectedphotons); 
  }
  return;
}



void VZGAnalyzer::genEventSetup(){
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
    std::vector<Boson<Particle>> Zll(genZlepCandidates_->begin()+1, genZlepCandidates_->end());
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

























void VZGAnalyzer::genAnalyze()
{
  std::vector<phys::Boson<phys::Particle>> genV;
  //std::vector<phys::Boson<phys::Particle>> genV(genVBHelper_.ZtoQ().size()+genVBHelper_.WtoQ().size());
  if(genVBHelper_.ZtoQ().size()>0)
  {
    genV=genVBHelper_.ZtoQ();
    genV.insert(genV.end(), genVBHelper_.WtoQ().begin(), genVBHelper_.WtoQ().end());
  }
  else if(genVBHelper_.WtoQ().size()>0)
  {
         genV=genVBHelper_.WtoQ();
    genV.insert(genV.end(), genVBHelper_.ZtoQ().begin(), genVBHelper_.ZtoQ().end());
  }
    //---------------------------------------- Single q analysis ----------------------------------------//
  std::vector<phys::Particle> genQuarksfromV;
  for (auto VB : genV)
  {
    theHistograms->fill("quarkfromV charge", "quarkfromV charge", 7, -7. / 6., 7. / 6., VB.daughter(0).charge(), theWeight);
    theHistograms->fill("quarkfromV charge", "quarkfromV charge", 7, -7. / 6., 7. / 6., VB.daughter(1).charge(), theWeight);

    theHistograms->fill("quarkfromV pt", "quarkfromV pt", 50, 0, 600, VB.daughter(0).pt(), theWeight);
    theHistograms->fill("quarkfromV pt", "quarkfromV pt", 50, 0, 600, VB.daughter(1).pt(), theWeight);

    genQuarksfromV.push_back(VB.daughter(0));
    genQuarksfromV.push_back(VB.daughter(1));
  }
  theHistograms->fill("0size_GENQuarksfromV_beforecuts", "0size_GENQuarksfromV_beforecuts", 10, -0.5, 9.5, genQuarksfromV.size(), theWeight);
  genQuarksfromV.erase(std::remove_if(genQuarksfromV.begin(), genQuarksfromV.end(), [](phys::Particle p)
                                      { return !KinematicsOK(p, ptcut, etacut); }),
                       genQuarksfromV.end());
  theHistograms->fill("0size_GENQuarksfromV_aftercuts", "0size_GENQuarksfromV_aftercuts", 10, -0.5, 9.5, genQuarksfromV.size(), theWeight);
 
  //----------------------------------------Kinematic Cuts on VB(qq)----------------------------------------//
  std::vector<phys::Boson<phys::Particle>> DiQuarks=genV;

  DiQuarks.erase(std::remove_if(DiQuarks.begin(), DiQuarks.end(), [](phys::Boson<phys::Particle> VB)
                                      { return !(KinematicsOK(VB.daughter(0), ptcut, etacut)&&KinematicsOK(VB.daughter(1), ptcut, etacut)); }),
                       DiQuarks.end());

 //----------------------------------------Kinematic Cuts GEN Jets AK4 & RECO Jets AK4----------------------------------------//
  std::vector<phys::Particle> selectedGENjets;
  foreach (const phys::Particle &jet, *genJets)
  {
    if (KinematicsOK(jet,ptcut,etacut)) // KinematicsOK(jet)
    {
      selectedGENjets.push_back(jet);
    }
  }
  std::vector<phys::Jet> selectedRECOjets;
  foreach (const phys::Jet &jet, *jets)
  {
    if (KinematicsOK(jet,ptcut,etacut)) // KinematicsOK(jet)
    {
      selectedRECOjets.push_back(jet);
    }
  }

  //----------------------------------------Matching efficiency ______ SINGLE QUARK/SINGLE GENJET--------------//
  std::vector<phys::Particle> jetsfromquarks;
  std::vector<std::pair<phys::Particle, phys::Particle>> nearestjetstoquark;

  for (auto quark : genQuarksfromV)
  {
    phys::Particle nearestjet;
    bool makesjet = false;
    theHistograms->fill("Pt_quark_den", " Pt_quark_den; GeV/c", 10, ptcut, 300, quark.pt(), theWeight);

    if (selectedGENjets.size() > 0)
    {
      std::stable_sort(selectedGENjets.begin(), selectedGENjets.end(), phys::DeltaRComparator(quark));
      nearestjet = selectedGENjets.at(0);
      nearestjetstoquark.push_back({quark, nearestjet});
      if (abs(physmath::deltaR(quark, nearestjet)) < 0.4)
      {
        jetsfromquarks.push_back(nearestjet);
        makesjet = true;
      }
    }
    if (makesjet && selectedGENjets.size() > 0)
    {
      theHistograms->fill("#QUARK=>GEN", "#QUARK=>GEN", 2, 0, 2, 1., theWeight);
      theHistograms->fill("Pt_quark_num", " Pt_quark_num; GeV/c", 10, ptcut, 300, quark.pt(), theWeight);
    }
    else
    {
      theHistograms->fill("#QUARK=>GEN", "#QUARK=>GEN", 2, 0, 2, 0., theWeight);
    }
  }
  for (auto pair : nearestjetstoquark)
  {
    ResolutionPlots(pair.first,pair.second,"Hadronization_",theWeight,"");
    theHistograms->fill("DeltaR_quark_vs_BestMatchedGENJet", "DeltaR_quark_vs_BestMatchedGENJet; #DeltaR", 20, 0, 0.5, abs(physmath::deltaR(pair.first, pair.second)), theWeight);
    theHistograms->fill("DeltaR_quark_jet_vs_pt", "DeltaR vs pt;pt [GeV/c] ; #DeltaR", 10, ptcut, 300, 20, 0, 0.2, pair.first.pt(),abs(physmath::deltaR(pair.first, pair.second)), theWeight);

  }
  theHistograms->fill("1size_GENjetsfromquarks", "size_GENjetsfromquarks", 10, -0.5, 9.5, jetsfromquarks.size(), theWeight);

  //----------------------------------------Matching efficiency ______ SINGLE GENJET/SINGLE RECOJET--------------//
  std::vector<phys::Particle> RECOjetsfromGENjets;
  std::vector<std::pair<phys::Particle, phys::Particle>> nearestRECOjetstoGENjets;

  for (auto genJet : selectedGENjets)
  {
    phys::Particle nearestRECOjet;
    bool isreconstructed = false;
    theHistograms->fill("Pt_genJet_den", " Pt_genJet_den; GeV/c", 10, ptcut, 300, genJet.pt(), theWeight);

    if (selectedRECOjets.size() > 0)
    {
      std::stable_sort(selectedRECOjets.begin(), selectedRECOjets.end(), phys::DeltaRComparator(genJet));
      nearestRECOjet = selectedRECOjets.at(0);
      nearestRECOjetstoGENjets.push_back({genJet, nearestRECOjet});
      if (abs(physmath::deltaR(genJet, nearestRECOjet)) < 0.4)
      {
        RECOjetsfromGENjets.push_back(nearestRECOjet);
        isreconstructed = true;
      }
    }
    if (isreconstructed && selectedRECOjets.size() > 0)
    {
      theHistograms->fill("#GEN=>RECO", "#GEN=>RECO", 2, 0, 2, 1., theWeight);
      theHistograms->fill("Pt_genJet_num", " Pt_genJet_num; GeV/c", 10, ptcut, 300, genJet.pt(), theWeight);

    }
    else
    {
      theHistograms->fill("#GEN=>RECO", "#GEN=>RECO", 2, 0, 2, 0., theWeight);
    }
  }
  for (auto pair : nearestRECOjetstoGENjets)
  {
    ResolutionPlots(pair.first,pair.second,"SingleJetsReconstruction_",theWeight,"");
    theHistograms->fill("DeltaR_GENjet_vs_BestMatchedRECOJet", "DeltaR_GENjet_vs_BestMatchedRECOJet; #DeltaR", 20, 0, 0.5, abs(physmath::deltaR(pair.first, pair.second)), theWeight);
    theHistograms->fill("DeltaR_jets_vs_pt", "DeltaR jets vs pt;pt [GeV/c] ; #DeltaR", 10, ptcut, 300, 20, 0, 0.2, pair.first.pt(),abs(physmath::deltaR(pair.first, pair.second)), theWeight);
  }
  theHistograms->fill("2size_GENjetsRECONSTRUCTED", "2size_GENjetsRECONSTRUCTED", 10, -0.5, 9.5, RECOjetsfromGENjets.size(), theWeight);

  //----------------------------------------Matching efficiency ______ QUARKS PAIR------------------------//

  std::vector<std::pair<phys::Particle, phys::Particle>> DijetsmatchedtoDiquark;
  std::vector<phys::Boson<phys::Particle>> DiJetsGEN;
  std::cout << ".................GEN to QUARKS MATCHING..............." << std::endl;

  std::cout << "DiQuarks size: " << DiQuarks.size() << std::endl;
  for (auto Diquark : DiQuarks)
  {

    bool firstmatches = false;
    phys::Particle jetmatchedtoFIRSTquark;

    bool secondmatches = false;
    phys::Particle jetmatchedtoSECONDquark;

    bool atleastonematches = false;
    bool bothmatch = false;
    std::cout << "selectedGENjets size: " << selectedGENjets.size() << std::endl;
    for (auto genJet : selectedGENjets)
    {
      double deltaR1 = abs(physmath::deltaR(Diquark.daughter(0), genJet));
      std::cout << "deltaR1= " << deltaR1 << std::endl;
      double deltaR2 = abs(physmath::deltaR(Diquark.daughter(1), genJet));
      std::cout << "deltaR2= " << deltaR2 << std::endl;

      if (deltaR1 < 0.4)
      {
        std::cout << "first matched" << std::endl;
        jetmatchedtoFIRSTquark = genJet;
        firstmatches = true;
      }
      if (deltaR2 < 0.4)
      {
        jetmatchedtoSECONDquark = genJet;
        std::cout << "second matched" << std::endl;
        secondmatches = true;
      }
    }

    bothmatch = (firstmatches && secondmatches);
    atleastonematches = (firstmatches || secondmatches);

    std::cout << "bothmatch= " << bothmatch << std::endl;
    std::cout << "atleastonematches= " << atleastonematches << std::endl;

    if (bothmatch)
    {
      std::cout << "both matched" << std::endl;
      std::cout << "reconstructing a boson from dijets matched to diquark" << std::endl;
      DiJetsGEN.push_back(phys::Boson<phys::Particle>(jetmatchedtoFIRSTquark, jetmatchedtoSECONDquark));
      DijetsmatchedtoDiquark.push_back({Diquark, phys::Boson<phys::Particle>(jetmatchedtoFIRSTquark, jetmatchedtoSECONDquark)});
      theHistograms->fill("#Bothmatched", "#Bothmatched", 2, 0, 2, 1., theWeight);
    }
    if (!bothmatch)
    {
      theHistograms->fill("#Bothmatched", "#Bothmatched", 2, 0, 2, 0., theWeight);
    }

    if (atleastonematches)
    {
      theHistograms->fill("#AtLeastONEmatches", "#AtLeastONEmatches", 2, 0, 2, 1., theWeight);
    }
    if (!atleastonematches)
    {
      theHistograms->fill("#AtLeastONEmatches", "#AtLeastONEmatches", 2, 0, 2, 0., theWeight);
    }
    if (atleastonematches && bothmatch)
    {
      theHistograms->fill("#Bothmatched|AtLeastONEmatches", "#Bothmatched|AtLeastONEmatches", 2, 0, 2, 1., theWeight);
    }
    if (atleastonematches && !bothmatch)
    {
      theHistograms->fill("#Bothmatched|AtLeastONEmatches", "#Bothmatched|AtLeastONEmatches", 2, 0, 2, 0., theWeight);
    }
  }
  theHistograms->fill("1.1size_GEN_DiJets", "1size_GEN_DiJets", 10, -0.5, 9.5, DiJetsGEN.size(), theWeight);
  for (auto DiJet : DiJetsGEN)
  {
    float mjj = (DiJet.daughter(0).p4() + DiJet.daughter(1).p4()).M();
    theHistograms->fill("mjj_GEN", "mjj_GEN", 10, 50, 120, mjj, theWeight);
  }

  //----------------------------------------Matching efficiency ______ RECO to GEN  PAIR------------------------//
  std::cout << ".................RECO to GEN MATCHING..............." << std::endl;

  // std::vector<std::pair<phys::Particle,phys::Particle>> DiRECOjetsmatchedtoDiGENjets;
  // std::vector<phys::Boson<phys::Particle>> DiJetsRECO;
  std::vector<phys::Boson<phys::Particle>> DiJetsGENreconstructed;
  std::cout << "GEN Dijets size: " << DiJetsGEN.size() << std::endl;
  for (auto DiJet : DiJetsGEN)
  {

    bool firstmatches = false;
    // phys::Particle jetmatchedtoFIRSTgen;

    bool secondmatches = false;
    // phys::Particle jetmatchedtoSECONDgen;

    bool atleastonematches = false;
    bool bothmatch = false;
    std::cout << "selectedRECOjets size: " << selectedRECOjets.size() << std::endl;
    for (auto recoJet : selectedRECOjets)
    {
      double deltaR1 = abs(physmath::deltaR(DiJet.daughter(0), recoJet));
      std::cout << "deltaR1= " << deltaR1 << std::endl;
      double deltaR2 = abs(physmath::deltaR(DiJet.daughter(1), recoJet));
      std::cout << "deltaR2= " << deltaR2 << std::endl;

      if (deltaR1 < 0.4)
      {
        std::cout << "first matched" << std::endl;
        // jetmatchedtoFIRSTgen=recoJet;
        firstmatches = true;
      }
      if (deltaR2 < 0.4)
      {
        std::cout << "second matched" << std::endl;
        // jetmatchedtoSECONDgen=recoJet;
        secondmatches = true;
      }
    }

    bothmatch = (firstmatches && secondmatches);
    atleastonematches = (firstmatches || secondmatches);

    std::cout << "bothmatch= " << bothmatch << std::endl;
    std::cout << "atleastonematches= " << atleastonematches << std::endl;

    if (bothmatch)
    {
      std::cout << "both matched" << std::endl;
      // std::cout << "reconstructing a boson from diRECOjets matched to diGENjets" << std::endl;
      // DiJetsRECO.push_back(phys::Boson<phys::Particle>(jetmatchedtoFIRSTgen, jetmatchedtoSECONDgen));
      DiJetsGENreconstructed.push_back(DiJet);
      // DiRECOjetsmatchedtoDiGENjets.push_back({DiJet,phys::Boson<phys::Particle>(jetmatchedtoFIRSTgen, jetmatchedtoSECONDgen)});

      theHistograms->fill("#RECOGEN_Bothmatched", "#RECOGEN_Bothmatched", 2, 0, 2, 1., theWeight);
    }
    if (!bothmatch)
    {
      theHistograms->fill("#RECOGEN_Bothmatched", "#RECOGEN_Bothmatched", 2, 0, 2, 0., theWeight);
    }

    if (atleastonematches)
    {
      theHistograms->fill("#RECOGEN_AtLeastONEmatches", "#RECOGEN_AtLeastONEmatches", 2, 0, 2, 1., theWeight);
    }
    if (!atleastonematches)
    {
      theHistograms->fill("#RECOGEN_AtLeastONEmatches", "#RECOGEN_AtLeastONEmatches", 2, 0, 2, 0., theWeight);
    }
    if (atleastonematches && bothmatch)
    {
      theHistograms->fill("#RECOGEN_Bothmatched|AtLeastONEmatches", "#RECOGEN_Bothmatched|AtLeastONEmatches", 2, 0, 2, 1., theWeight);
    }
    if (atleastonematches && !bothmatch)
    {
      theHistograms->fill("#RECOGEN_Bothmatched|AtLeastONEmatches", "#RECOGEN_Bothmatched|AtLeastONEmatches", 2, 0, 2, 0., theWeight);
    }
  }
  // theHistograms->fill("2size_RECO_DiJets", "1size_RECO_DiJets", 10, -0.5, 9.5, DiJetsRECO.size(), theWeight);
  //  for (auto DiJet:DiJetsRECO)
  //  {
  //      float mjj = (DiJet.daughter(0).p4() + DiJet.daughter(1).p4()).M();
  //      theHistograms->fill("mjj_RECO", "mjj_RECO", 10, 50, 120, mjj, theWeight);
  //  }

  //----------------------------------------Matching efficiency ______ ALGORITHM------------------------//
  std::cout << ".................Algorithm efficiency..............." << std::endl;

  std::vector<std::pair<phys::Particle, phys::Particle>> DiRECOjetsmatchedtoDiGENjets;
  std::vector<phys::Boson<phys::Particle>> DiJetsRECO;


// reconstructing pairs


    std::map<std::string, Boson<phys::Jet>> Candidates;

    std::vector<phys::Boson<phys::Jet>> JetPairs;

    for (int i = 0; i < selectedRECOjets.size(); i++)
    {

      for (int j = i + 1; j < selectedRECOjets.size(); j++)
      {
        JetPairs.push_back(phys::Boson<phys::Jet>(selectedRECOjets.at(i), selectedRECOjets.at(j)));
      }
    }

    theHistograms->fill("2size_Jetspairs_RECO", "size_Jetsparis_RECO", 10, -0.5, 9.5, JetPairs.size(), theWeight);
    std::cout << "#reco jets pairs : " << JetPairs.size() << std::endl;

    if (JetPairs.size() > 0)
    {

      // 1st reconstruction model: comparison with WMass
      std::stable_sort(JetPairs.begin(), JetPairs.end(), phys::MassComparator(phys::WMASS));
      Candidates["mW"] = JetPairs.at(0);

      //2nd reconstruction model: comparison with ZMass
      std::stable_sort(JetPairs.begin(), JetPairs.end(), phys::MassComparator(phys::ZMASS));
      Candidates["mZ"] = JetPairs.at(0);

      // 3rd reconstruction model: maximization of candidate Pt
      std::stable_sort(JetPairs.begin(), JetPairs.end(), phys::ScalarSumPtComparator());
      Candidates["maxVPt"] = JetPairs.at(0);

      // // // 4th reconstruction model: minimization of total Pt of Zjj system
      // std::vector<phys::Particle> ZZjj;
      // phys::Particle ZZjjCandidate;
      // for (uint i = 0; i < JetPairs.size(); i++)
      // {
      //   phys::Particle totState(ZZ->p4() + (JetPairs.at(i)).p4());
      //   ZZjj.push_back(totState.p4());
      // }
      // std::stable_sort(ZZjj.begin(), ZZjj.end(), phys::PtComparator());
      // ZZjjCandidate = ZZjj.back();
      // for (uint i = 0; i < JetPairs.size(); i++)
      //   if ((JetPairs.at(i)).p4() == (ZZjjCandidate.p4() - ZZ->p4()))
      //     Candidates["minTotPt"] = JetPairs.at(i);

      // 4th reconstruction model: comparison with a mean value between ZMass and WMass
      std::stable_sort(JetPairs.begin(), JetPairs.end(), phys::MassComparator(0.2 * phys::ZMASS + 0.8 * phys::WMASS));
      Candidates["m8W2Z"] = JetPairs.at(0);
      // 5th reconstruction model: comparison with a mean value between ZMass and WMass

      std::stable_sort(JetPairs.begin(), JetPairs.end(), phys::Mass2Comparator(phys::ZMASS, phys::WMASS));
      Candidates["mWZ"] = JetPairs.at(0);

      for (auto Candidate : Candidates)
      {
        theHistograms->fill("mjj_" + Candidate.first + "_Candidate", "mjj_" + Candidate.first + "_Candidate", 10, 50, 120, Candidate.second.mass(), theWeight);
      }
    }
  
  theHistograms->fill("2size_RECO_JetPairs", "2size_RECO_JetPairs", 10, -0.5, 9.5, JetPairs.size(), theWeight);


  std::cout << "#true gen jets pairs reconstructed: " << DiJetsGENreconstructed.size() << std::endl;

  for (auto DiJet : DiJetsGENreconstructed)
  {
    theHistograms->fill("mjj_den" , " mjj_den; GeV/c^{2}" , 10, 50, 120, DiJet.mass(), theWeight);
    theHistograms->fill("Pt_den" , " Pt_den; GeV/c" , 10, 0, 300, DiJet.pt(), theWeight);
    std::cout<<""<<std::endl;
    int k=0;


    //----------------------------------------MATCHED Jets total mass----------------------------------//
    for (auto Candidate : Candidates)
    {

        std::cout<<""<<std::endl;

        bool truepair = false;
        std::cout << "Algorithm: " << Candidate.first << std::endl;
        phys::Jet jetRECOA = Candidate.second.daughter(0);
        phys::Jet jetRECOB = Candidate.second.daughter(1);
        phys::Particle jetGENA = DiJet.daughter(0);
        phys::Particle jetGENB = DiJet.daughter(1);
        double deltaRAA = abs(physmath::deltaR(jetGENA, jetRECOA));
        double deltaRAB = abs(physmath::deltaR(jetGENA, jetRECOB));
        double deltaRBA = abs(physmath::deltaR(jetGENB, jetRECOA));
        double deltaRBB = abs(physmath::deltaR(jetGENB, jetRECOB));
        if ((deltaRAA < 0.4 && deltaRBB < 0.4) || (deltaRAB < 0.4 && deltaRBA < 0.4))
        {
          truepair = true;
        }
        if (truepair)
        {
          std::cout << "the algorithm selected a reco pair matched to a gen pair" << std::endl;
          theHistograms->fill("#Algorithm_" + Candidate.first, "#Algorithm_" + Candidate.first, 2, 0, 2, 1., theWeight);
          theHistograms->fill("PASSED mjj_" + Candidate.first + "_Candidate", " PASSED mjj_" + Candidate.first + "_Candidate", 10, 50, 120, DiJet.mass(), theWeight);
          theHistograms->fill("mjj_" + Candidate.first + "_num", " mjj_" + Candidate.first + "_num; GeV/c^{2}", 10, 50, 120, DiJet.mass(), theWeight);
          theHistograms->fill("Pt_" + Candidate.first + "_num", " Pt_" + Candidate.first + "_num; GeV/c", 10, 0, 300, DiJet.pt(), theWeight);

        }
        else
        {
          theHistograms->fill("FAILED mjj_" + Candidate.first + "_Candidate", " FAILED mjj_" + Candidate.first + "_Candidate", 10, 50, 120, DiJet.mass(), theWeight);
          std::cout << "the algorithm selected a reco pair NOT matched to a gen pair" << std::endl;
          theHistograms->fill("#Algorithm_" + Candidate.first, "#Algorithm_" + Candidate.first, 2, 0, 2, 0., theWeight);
        }
        if (truepair && Candidate.first=="mWZ")
        {
          ResolutionPlots(DiJet,Candidate.second,"VectorBosonReconstruction_",theWeight,"");
        }
    }
  }

}

void VZGAnalyzer::ResolutionPlots(const phys::Particle &gen, const phys::Particle &reco, std::string prename, const float weight, std::string suffix)
{
  std::string where;
  if (fabs(gen.eta()) < 2.4)
  {
    where = "Barrel";
  }
  else if (fabs(gen.eta()) > 2.4 && fabs(gen.eta()) < 4.7)
  {
    where = "Endcap";
  }
  double delta_charge = (reco.charge() - gen.charge());
  double delta_mass = (reco.mass() - gen.mass());
  double delta_trmass = (reco.p4().Mt() - gen.p4().Mt());
  double delta_pt = (reco.pt() - gen.pt());
  double delta_eta = (reco.eta() - gen.eta());
  double delta_phi = (reco.phi() - gen.phi());
  // double JJdeltaR = abs(physmath::deltaR(gen, reco));

  double res_mass = (reco.mass() - gen.mass()) / gen.mass();
  double res_trmass = (reco.p4().Mt() - gen.p4().Mt()) / gen.p4().Mt();
  double res_pt = (reco.pt() - gen.pt()) / gen.pt();

  if (where == "Barrel" || where == "Endcap")
  {
    std::string name = "Delta";
    theHistograms->fill(prename + name + "_charge_" + where + suffix, prename + "#" + name + "_charge_" + where + ";#Delta charge", 9, -4.5, 4.5, delta_charge, weight);
    theHistograms->fill(prename + name + "_mass_" + where + suffix, prename + "#" + name + "_mass_" + where + ";#Delta mass [GeV/c^2]", 10, -20, 20, delta_mass, weight);
    theHistograms->fill(prename + name + "_trmass_" + where + suffix, prename + "#" + name + "_trmass_" + where + ";#Delta Tr mass [GeV/c^2]", 10, -20, 20, delta_trmass, weight);
    theHistograms->fill(prename + name + "_pt_" + where + suffix, prename + "#" + name + "_pt_" + where + ";#Delta pt [GeV/c]", 10, -20, 20, delta_pt, weight);

    theHistograms->fill(prename + name + "_eta_" + where + suffix, prename + "#" + name + "_eta_" + where + ";#Delta eta", 10, -1, 1, delta_eta, weight);
    theHistograms->fill(prename + name + "_phi_" + where + suffix, prename + "#" + name + "_#phi_" + where + ";#Delta phi", 10, -1, 1, delta_phi, weight);

    theHistograms->fill(prename + name + "_mass_vs_mass_" + where + suffix, prename + "#" + name + "_mass_vs_mass_" + where + ";mass[GeV/c^2]; #Delta mass [GeV/c^2]", 20, 0, 200, 10, -20, 20, gen.mass(), delta_mass, weight);
    theHistograms->fill(prename + name + "_trmass_vs_trmass_" + where + suffix, prename + "#" + name + "_trmass_vs_trmass_" + where + ";Trmass[GeV/c^2]; #Delta Trmass [GeV/c^2]", 20, 0, 200, 10, -20, 20, gen.p4().Mt(), delta_trmass, weight);
    theHistograms->fill(prename + name + "_pt_vs_pt_" + where + suffix, prename + "#" + name + "_pt_vs_p_{t}_" + where + ";p_{t} [GeV/c^2]; #Delta p_{t} [GeV/c^2]", 30, 0, 300, 10, -20, 20, gen.pt(), delta_pt, weight);

    name = "Res";
    //  theHistograms->fill(prename + name + "_charge_" + suffix, prename + name + "'s charge", 5, -2.5, 2.5, Jet1.charge(), weight);
    theHistograms->fill(prename + name + "_mass_" + where + suffix, prename + name + "_mass" + where, 10, -1, 1, res_mass, weight);
    theHistograms->fill(prename + name + "_trmass_" + where + suffix, prename + name + "_trmass_" + where, 10, -1, 1, res_trmass, weight);
    theHistograms->fill(prename + name + "_pt_" + where + suffix, prename + name + "_p_{t}_" + where, 10, -1, 1, res_pt, weight);
    //  theHistograms->fill(prename + name + "_Y_" + suffix, prename + name + "'s Y", 50, -5, 5, Jet1.rapidity(), weight);
    //  theHistograms->fill(prename + name + "_eta_" + suffix, prename + name + "'s #eta", 50, -9, 9, Jet1.eta(), weight);
    //  theHistograms->fill(prename + name + "_phi_" + suffix, prename + name + "'s #phi", 50, -3.5, 3.5, Jet1.phi(), weight);
  }
}




void VZGAnalyzer::PhotonvsJet()
{
    std::vector<phys::Particle> genQuarks;
       foreach (const Particle &p, *genParticles)
       {
              if  (abs(p.id()) < 10) // Is it a quark? 
              {
                     theHistograms->fill("quark charge", "quark charge", 7, -7. / 6., 7. / 6., p.charge(), theWeight);
                     theHistograms->fill("quark pt", "quark pt", 50, 0, 600, p.pt(), theWeight);
                     genQuarks.push_back(Particle(p));
                     

              }
       }
    //---------------------------------------- Single q analysis and cuts ----------------------------------------//
  // std::vector<phys::Particle> genQuarksfromV;
  // for (auto VB : genV)
  // {
  //   theHistograms->fill("quarkfromV charge", "quarkfromV charge", 7, -7. / 6., 7. / 6., VB.daughter(0).charge(), theWeight);
  //   theHistograms->fill("quarkfromV charge", "quarkfromV charge", 7, -7. / 6., 7. / 6., VB.daughter(1).charge(), theWeight);

  //   theHistograms->fill("quarkfromV pt", "quarkfromV pt", 50, 0, 600, VB.daughter(0).pt(), theWeight);
  //   theHistograms->fill("quarkfromV pt", "quarkfromV pt", 50, 0, 600, VB.daughter(1).pt(), theWeight);

  //   genQuarksfromV.push_back(VB.daughter(0));
  //   genQuarksfromV.push_back(VB.daughter(1));
  // }
  // theHistograms->fill("0size_GENQuarksfromV_beforecuts", "0size_GENQuarksfromV_beforecuts", 10, -0.5, 9.5, genQuarksfromV.size(), theWeight);
  // genQuarksfromV.erase(std::remove_if(genQuarksfromV.begin(), genQuarksfromV.end(), [](phys::Particle p)
  //                                     { return !KinematicsOK(p, ptcut, etacut); }),
  //                      genQuarksfromV.end());
  // theHistograms->fill("0size_GENQuarksfromV_aftercuts", "0size_GENQuarksfromV_aftercuts", 10, -0.5, 9.5, genQuarksfromV.size(), theWeight);
 

 //----------------------------------------Kinematic Cuts GEN Jets AK4 & RECO Jets AK4----------------------------------------//
  std::vector<phys::Particle> selectedGENjets;
  foreach (const phys::Particle &jet, *genJets)
  {
    if (true) // KinematicsOK(jet,ptcut,etacut)
    {
      selectedGENjets.push_back(jet);
    }
  }
  std::vector<phys::Jet> selectedRECOjets;
  foreach (const phys::Jet &jet, *jets)
  {
    if (true) // KinematicsOK(jet,ptcut,etacut)
    {
      selectedRECOjets.push_back(jet);
    }
  }

  //----------------------------------------Matching efficiency ______ SINGLE QUARK/SINGLE GENJET--------------//
  std::vector<phys::Particle> jetsfromquarks;
  std::vector<std::pair<phys::Particle, phys::Particle>> nearestjetstoquark;

  for (auto quark : genQuarks)
  {
    phys::Particle nearestjet;
    bool makesjet = false;
    theHistograms->fill("Pt_quark_den", " Pt_quark_den; GeV/c", 10, 0, 300, quark.pt(), theWeight);

    if (selectedGENjets.size() > 0)
    {
      std::stable_sort(selectedGENjets.begin(), selectedGENjets.end(), phys::DeltaRComparator(quark));
      nearestjet = selectedGENjets.at(0);
      nearestjetstoquark.push_back({quark, nearestjet});
      if (abs(physmath::deltaR(quark, nearestjet)) < 0.4)
      {
        jetsfromquarks.push_back(nearestjet);
        makesjet = true;
      }
    }
    if (makesjet && selectedGENjets.size() > 0)
    {
      theHistograms->fill("#QUARK=>GEN", "#QUARK=>GEN", 2, 0, 2, 1., theWeight);
      theHistograms->fill("Pt_quark_num", " Pt_quark_num; GeV/c", 10, 0, 300, quark.pt(), theWeight);
    }
    else
    {
      theHistograms->fill("#QUARK=>GEN", "#QUARK=>GEN", 2, 0, 2, 0., theWeight);
    }
  }
  for (auto pair : nearestjetstoquark)
  {
    ResolutionPlots(pair.first,pair.second,"Hadronization_",theWeight,"");
    theHistograms->fill("DeltaR_quark_vs_BestMatchedGENJet", "DeltaR_quark_vs_BestMatchedGENJet; #DeltaR", 20, 0, 0.5, abs(physmath::deltaR(pair.first, pair.second)), theWeight);
    theHistograms->fill("DeltaR_quark_jet_vs_pt", "DeltaR vs pt;pt [GeV/c] ; #DeltaR", 10, 0, 300, 20, 0, 0.2, pair.first.pt(),abs(physmath::deltaR(pair.first, pair.second)), theWeight);

  }
  theHistograms->fill("1size_GENjetsfromquarks", "size_GENjetsfromquarks", 10, -0.5, 9.5, jetsfromquarks.size(), theWeight);


  //----------------------------------------Matching efficiency ______ SINGLE GENJET/SINGLE RECOJET--------------//
  std::vector<phys::Particle> RECOjetsfromGENjets;
  std::vector<std::pair<phys::Particle, phys::Particle>> nearestRECOjetstoGENjets;

  for (auto genJet : jetsfromquarks)
  {
    phys::Particle nearestRECOjet;
    bool isreconstructed = false;
    theHistograms->fill("Pt_genJet_den", " Pt_genJet_den; GeV/c", 10, 0, 300, genJet.pt(), theWeight);

    if (selectedRECOjets.size() > 0)
    {
      std::stable_sort(selectedRECOjets.begin(), selectedRECOjets.end(), phys::DeltaRComparator(genJet));
      nearestRECOjet = selectedRECOjets.at(0);
      nearestRECOjetstoGENjets.push_back({genJet, nearestRECOjet});
      if (abs(physmath::deltaR(genJet, nearestRECOjet)) < 0.4)
      {
        RECOjetsfromGENjets.push_back(nearestRECOjet);
        isreconstructed = true;
      }
    }
    if (isreconstructed && selectedRECOjets.size() > 0)
    {
      theHistograms->fill("#GEN=>RECO", "#GEN=>RECO", 2, 0, 2, 1., theWeight);
      theHistograms->fill("Pt_genJet_num", " Pt_genJet_num; GeV/c", 10, 0, 300, genJet.pt(), theWeight);

    }
    else
    {
      theHistograms->fill("#GEN=>RECO", "#GEN=>RECO", 2, 0, 2, 0., theWeight);
    }
  }
  for (auto pair : nearestRECOjetstoGENjets)
  {
    ResolutionPlots(pair.first,pair.second,"SingleJetsReconstruction_",theWeight,"");
    theHistograms->fill("DeltaR_GENjet_vs_BestMatchedRECOJet", "DeltaR_GENjet_vs_BestMatchedRECOJet; #DeltaR", 20, 0, 0.5, abs(physmath::deltaR(pair.first, pair.second)), theWeight);
    theHistograms->fill("DeltaR_jets_vs_pt", "DeltaR jets vs pt;pt [GeV/c] ; #DeltaR", 10, 0, 300, 20, 0, 0.2, pair.first.pt(),abs(physmath::deltaR(pair.first, pair.second)), theWeight);
  }
  theHistograms->fill("2size_GENjetsRECONSTRUCTED", "2size_GENjetsRECONSTRUCTED", 10, -0.5, 9.5, RECOjetsfromGENjets.size(), theWeight);
}
