#include "VVXAnalysis/TreeAnalysis/interface/VVjjHelper.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/DataFormats/interface/TypeDefs.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using namespace std;
using namespace phys;
using namespace physmath;


void VVjjHelper::test(int number = 0) {
  cout << "pippo!! test number " << number << endl;
}



void VVjjHelper::LeptonSearch(const vector<Particle> &genparticles, string eventkind, vector<Particle> &lepm, vector<Particle> &lepp, vector<Particle> &neutrino){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~ Function to find Leptons ~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  // Reset Data Member containers
  leptons_.clear();
  
  // Get leptons from genParticles
  foreach(const Particle &gen, genparticles){
    if((abs(gen.id()) != 11 && abs(gen.id()) != 13 && abs(gen.id()) != 12 && abs(gen.id()) != 14) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) 
      continue;

    // Leptons accepted are:
    // - electrons with pt >= 7. and |eta| <= 2.5
    if(abs(gen.id()) == 11 && gen.pt() >= 7. && abs(gen.eta()) <= 2.5){
      leptons_.push_back(gen);
    }
    // - muons with pt >= 5. and |eta| <= 2.4    
    else if(abs(gen.id()) == 13 && gen.pt() >= 5. && abs(gen.eta()) <= 2.4){
      leptons_.push_back(gen);
    }
    // - electronic or muonic neutrinos
    else if(abs(gen.id()) == 12 || abs(gen.id()) == 14){
      neutrino.push_back(gen);
    }
  }

  // sorting leptons by pt
  sort(leptons_.begin(), leptons_.end(), PtComparator());

  // different selection for differents event types: 1 for WZ and 0 for ZZ
  int switchcase = -1;
  if(eventkind == "WZ" && leptons_.size() == 3){
    switchcase = 1;
  }
  else if(eventkind == "ZZ" && leptons_.size() == 4){
    switchcase = 0;
  }

  // first leptons selection is the same
  switch(switchcase){
  case 0:{// ZZevent
    // things for ZZ event
  }
    break;
   
  case 1:{// WZ event
    // first lepton must have at least 20GeV pt
    if(leptons_[0].pt() < 20.){
      leptons_.clear();
      //cout << "Leptons_ size " << leptons_.size() << endl;
    }

    // second lepton must have at least 12GeV pt if electron or 10GeV pt if muon
    switch(abs(leptons_[1].id())){
    case 11:{
      if(leptons_[1].pt() < 12.){
	leptons_.clear();
	//cout << "Leptons_ size " << leptons_.size() << endl;
      }
    }
      break;

    case 13:{
      if(leptons_[1].pt() < 10.){
	leptons_.clear();
	//cout << "Leptons_ size " << leptons_.size() << endl;
      }
    }
      break;

    default:{
      cout << "ERROR: second lepton's ID is: " << leptons_[1].id() << endl << endl;
      return;
    }
    }
  }
    break;
     
  default:{
    //cout << "ERROR: lepton's number is: " << leptons_.size() << endl << endl;
    return;
  }
  }

  TLorentzVector Ptot;

  foreach(const Particle neu, neutrino){
    Ptot += neu.p4();
    neutrinos_.push_back(neu);
  }

  if(leptons_.size() != 0){
    foreach(const Particle lep, leptons_){
      Ptot += lep.p4();
      
      if(lep.charge() < 0)
	lepm.push_back(lep);
      else
	lepp.push_back(lep);
    }
  }
  
  histo_->fill("AllGenlllnu_mass",   "m 3 leptons and #nu",     150, 0, 1500, Ptot.M());
  histo_->fill("AllGenlllnu_trmass", "m_{T} 3 leptons and #nu", 150, 0, 1500, Ptot.Mt());
  
  if(Ptot.M() < 100.){ //before was 165.
    leptons_.clear();
    lepm.clear();
    lepp.clear();
  }
}



void VVjjHelper::FindLeadingJets(const vector<Particle> *jetcollection, Particle &Jet0, Particle &Jet1, const vector<Particle> *particlecollection){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~ Function to find leading jets ~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  vector<Particle> jets;
  bool leptonMatch = false; 
  
  foreach(const Particle &jet, *jetcollection){ // in genJets are included genleptons
    leptonMatch = false;
    
    foreach(const Particle &gen, *particlecollection){ // get jets only
      if(deltaR(gen,jet) < 0.4 && (abs(gen.id()) == 11 || abs(gen.id()) == 13))
	leptonMatch = true;
    }
    
    if(!leptonMatch){ // jets must have rapidity less than 4.7
      if(fabs(jet.eta()) < 4.7 && jet.pt() > 30.)
	jets.push_back(jet);
    }
  }
  
  if(jets.size() < 2){
    return;
  }
  
  // leading jets are those with higher pt
  sort(jets.begin(), jets.end(), PtComparator());
  Jet0 = jets[0];
  Jet1 = jets[1];  
}



void VVjjHelper::FindLeadingJets(const vector<Jet> *jetcollection, Particle &Jet0, Particle &Jet1){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~ Function to find leading jets ~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  vector<Particle> jets;

  foreach(const Particle &jet, *jetcollection){
    if(fabs(jet.eta()) < 4.7 && jet.pt() > 30.)
      jets.push_back(jet);
  }
  
  if(jets.size() < 2){
    return;
  }
  
  // leading jets are those with higher pt
  sort(jets.begin(), jets.end(), PtComparator());
  Jet0 = jets[0];
  Jet1 = jets[1];
  
}



  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~ Histograms functions ~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void VVjjHelper::PlotParticle(const Particle &particle, string name, const float weight, string suffix){
  histo_->fill(name + "_charge_" + suffix, name + "'s charge",   5, -2.5,   2.5, particle.charge()  , weight);
  histo_->fill(name + "_pt_" + suffix,     name + "'s p_{t}",  140,  0  , 700  , particle.pt()      , weight);
  histo_->fill(name + "_Y_" + suffix,      name + "'s Y",       50, -5  ,   5  , particle.rapidity(), weight);
  histo_->fill(name + "_eta_" + suffix,    name + "'s #eta",    50, -9  ,   9  , particle.eta()     , weight);
  histo_->fill(name + "_phi_" + suffix,    name + "'s #phi",    50, -3.5,   3.5, particle.phi()     , weight);

  if(particle.charge() != 0){
  histo_->fill(name + "_mass_" + suffix,   name + "'s mass",   200,  0  , 400  , particle.mass()    , weight);
  histo_->fill(name + "_trmass_" + suffix, name + "'s trmass", 200,  0  , 400  , particle.p4().Mt() , weight);
  }
}



void VVjjHelper::PlotJets(const Particle &Jet0, const Particle &Jet1, string prename, const float weight, string suffix){
  string name = "J0";  
  histo_->fill(prename + name + "_charge_" + suffix, prename + name + "'s charge",   5, -2.5,   2.5, Jet0.charge()  , weight);
  histo_->fill(prename + name + "_mass_" + suffix,   prename + name + "'s mass",    25,  0  , 250  , Jet0.mass()    , weight);
  histo_->fill(prename + name + "_trmass_" + suffix, prename + name + "'s trmass", 200,  0  , 400  , Jet0.p4().Mt() , weight);
  histo_->fill(prename + name + "_pt_" + suffix,     prename + name + "'s p_{t}",  140,  0  , 700  , Jet0.pt()      , weight);
  histo_->fill(prename + name + "_Y_" + suffix,      prename + name + "'s Y",       50, -5  ,   5  , Jet0.rapidity(), weight);
  histo_->fill(prename + name + "_eta_" + suffix,    prename + name + "'s #eta",    50, -9  ,   9  , Jet0.eta()     , weight);
  histo_->fill(prename + name + "_phi_" + suffix,    prename + name + "'s #phi",    50, -3.5,   3.5, Jet0.phi()     , weight);
  
  name = "J1";  
  histo_->fill(prename + name + "_charge_" + suffix, prename + name + "'s charge",   5, -2.5,   2.5, Jet1.charge()  , weight);
  histo_->fill(prename + name + "_mass_" + suffix,   prename + name + "'s mass",    25,  0  , 250  , Jet1.mass()    , weight);
  histo_->fill(prename + name + "_trmass_" + suffix, prename + name + "'s trmass", 200,  0  , 400  , Jet1.p4().Mt() , weight);
  histo_->fill(prename + name + "_pt_" + suffix,     prename + name + "'s p_{t}",  140,  0  , 700  , Jet1.pt()      , weight);
  histo_->fill(prename + name + "_Y_" + suffix,      prename + name + "'s Y",       50, -5  ,   5  , Jet1.rapidity(), weight);
  histo_->fill(prename + name + "_eta_" + suffix,    prename + name + "'s #eta",    50, -9  ,   9  , Jet1.eta()     , weight);
  histo_->fill(prename + name + "_phi_" + suffix,    prename + name + "'s #phi",    50, -3.5,   3.5, Jet1.phi()     , weight);

  name = "JJ";
  TLorentzVector JJp4 = Jet0.p4() + Jet1.p4();
  double JJdeltaEta = Jet0.eta() - Jet1.eta();
  double JJdeltaPhi = physmath::deltaPhi(Jet0.phi(), Jet1.phi());
  double JJdeltaR = abs(physmath::deltaR(Jet0, Jet1));
  
  histo_->fill(prename + name + "_mass_" + suffix,     "Leading Jets' mass",       200,  0  , 4000  , JJp4.M()  , weight);
  histo_->fill(prename + name + "_trmass_" + suffix,   "Leading Jets' trmass",     300,  0  , 4500  , JJp4.Mt() , weight);
  histo_->fill(prename + name + "_deltaEta_" + suffix, "Leading Jets' #Delta#eta", 100, -9  ,    9  , JJdeltaEta, weight); 
  histo_->fill(prename + name + "_deltaR_" + suffix,   "Leading Jets' #DeltaR",     25, -0.5,    9  , JJdeltaR  , weight);
  histo_->fill(prename + name + "_deltaPhi_" + suffix, "Leading Jets' #Delta#phi",  50, -3.5,    3.5, JJdeltaPhi, weight);  
}



template <class BOS>
void VVjjHelper::PlotBoson(const BOS &particle, string name, const float weight, string suffix){
  histo_->fill(name + "_charge_" + suffix, name + "'s charge",   5, -2.5,   2.5, particle.charge()  , weight);
  histo_->fill(name + "_trmass_" + suffix, name + "'s trmass", 200, 50  , 450  , particle.p4().Mt() , weight);
  histo_->fill(name + "_pt_" + suffix,     name + "'s p_{t}",  140,  0  , 700  , particle.pt()      , weight);
  histo_->fill(name + "_Y_" + suffix,      name + "'s Y",       50, -5  ,   5  , particle.rapidity(), weight);
  histo_->fill(name + "_eta_" + suffix,    name + "'s #eta",    50, -9  ,   9  , particle.eta()     , weight);
  histo_->fill(name + "_phi_" + suffix,    name + "'s #phi",    50, -3.5,   3.5, particle.phi()     , weight);
  histo_->fill(name + "_massvstrmass_" + suffix, name + "'s mass(x) vs trmass(y)", 400, 0, 400, 400, 0, 400, particle.mass(), particle.p4().Mt(), weight);

  double dEta = particle.daughter(0).eta() - particle.daughter(1).eta();
  double dPhi = deltaPhi(particle.daughter(0).phi(), particle.daughter(1).phi());
  double dR = abs(deltaR(particle.daughter(0), particle.daughter(1)));
  
  histo_->fill(name + "_deltaEta_" + suffix,     name + "'s #Delta#eta", 50, -6.5, 6.5, dEta, weight);
  histo_->fill(name + "_deltaR_" + suffix,       name + "'s #DeltaR",    50, -0.5, 6.5, dR  , weight);
  histo_->fill(name + "_deltaPhi_" + suffix,     name + "'s #Delta#phi", 50, -3.5, 3.5, dPhi, weight);

  if(particle.charge() == 0){
  histo_->fill(name + "_mass_" + suffix,   name + "'s mass",   100, 50  , 130  , particle.mass()    , weight);
  }
}

// Explicit template instantiation to avoid errors in linking process while compiling
template void VVjjHelper::PlotBoson<BosonParticle>(const BosonParticle& particle, string name, const float weight, string suffix);
template void VVjjHelper::PlotBoson<BosonLepton>(const BosonLepton& particle, string name, const float weight, string suffix);



template <class DiBOS>
void VVjjHelper::PlotDiBoson(const DiBOS& particle, string name, const float weight, string suffix){
  histo_->fill(name + "_trmass_" + suffix, name + "'s trmass", 100, 150, 1650, particle.p4().Mt(), weight);
  histo_->fill(name + "_pt_" + suffix,     name + "'s p_{t}",  140,   0,  700, particle.pt()     , weight);
  histo_->fill(name + "_massvstrmass_" + suffix, name + "'s mass(x) vs trmass(y)",            268, 160, 1500, 268, 160, 1500, particle.mass(),            particle.p4().Mt(),          weight);
  histo_->fill(name + "_pt1vspt2_" + suffix,     "first's p_{t} (x) vs second's p_{t} (y)",   260,   0,  650, 260,   0,  650, particle.first().pt(),      particle.second().pt(),      weight);
  histo_->fill(name + "_mass1vsmass2_" + suffix, "first's trmass (x) vs second's trmass (y)", 260,   0,  650, 260,   0,  650, particle.first().p4().Mt(), particle.second().p4().Mt(), weight);
  
  double dEta = particle.first().eta() - particle.second().eta();
  double dPhi = deltaPhi(particle.first().phi(), particle.second().phi());
  double dR = abs(deltaR(particle.first(), particle.second()));
  
  histo_->fill(name + "_deltaEta_" + suffix, name + "'s #Delta#eta", 50, -6.5, 6.5, dEta, weight);
  histo_->fill(name + "_deltaR_" + suffix,   name + "'s #DeltaR",    50, -0.5, 6.5, dR  , weight);
  histo_->fill(name + "_deltaPhi_" + suffix, name + "'s #Delta#phi", 50, -3.5, 3.5, dPhi, weight);
  
  int firstID = abs(particle.first().daughter(0).id()) + abs(particle.first().daughter(1).id());
  int secondID = abs(particle.second().daughter(0).id()) + abs(particle.second().daughter(1).id());
  
  histo_->fill(name + "_IDs_" + suffix, name + "'s IDs", 5, 19, 29, 5, 20, 30, firstID, secondID, weight);

  if(particle.charge() == 0){
  histo_->fill(name + "_mass_" + suffix,   name + "'s mass",   100, 100, 600, particle.mass()   , weight);
  }
}

// Explicit template instantiation to avoid errors in linking process while compiling
template void VVjjHelper::PlotDiBoson<DiBosonParticle>(const DiBosonParticle& particle, string name, const float weight, string suffix);
template void VVjjHelper::PlotDiBoson<DiBosonLepton>(const DiBosonLepton& particle, string name, const float weight, string suffix);



  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~ Getter functions ~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
unsigned int VVjjHelper::GetLeptonsNumber(){
  return leptons_.size();
}

unsigned int VVjjHelper::GetNeutrinosNumber(){
  return neutrinos_.size();
}
