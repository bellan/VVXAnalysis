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


void VVjjHelper::test(int number = 0){
  cout << "Test number " << number << endl;
}



void VVjjHelper::printTime(float btime, float etime){
  cout << "\nExecution time: " << (int)((etime - btime)/3600) << " h " << (((int)(etime - btime)%3600)/60) << " m " << etime - btime - (int)((etime - btime)/3600)*3600 - (((int)(etime - btime)%3600)/60)*60 << " s." << endl;
}



float VVjjHelper::getSpline(double ZZmass) const {
  return MELASpline.Eval(ZZmass);
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
  histo_->fill(name + "_mass_" + suffix,   name + "'s mass",   100,  0  , 400  , particle.mass()    , weight);
  histo_->fill(name + "_trmass_" + suffix, name + "'s trmass", 100,  0  , 400  , particle.p4().Mt() , weight);
  }
}



void VVjjHelper::PlotJets(const Particle &Jet0, const Particle &Jet1, string prename, const float weight, string suffix){
  string name = "J0";  
  histo_->fill(prename + name + "_charge_" + suffix, prename + name + "'s charge",   5, -2.5,   2.5, Jet0.charge()  , weight);
  histo_->fill(prename + name + "_mass_" + suffix,   prename + name + "'s mass",    63,  0  , 252  , Jet0.mass()    , weight);
  histo_->fill(prename + name + "_trmass_" + suffix, prename + name + "'s trmass", 100,  0  , 400  , Jet0.p4().Mt() , weight);
  histo_->fill(prename + name + "_pt_" + suffix,     prename + name + "'s p_{t}",  140,  0  , 700  , Jet0.pt()      , weight);
  histo_->fill(prename + name + "_Y_" + suffix,      prename + name + "'s Y",       50, -5  ,   5  , Jet0.rapidity(), weight);
  histo_->fill(prename + name + "_eta_" + suffix,    prename + name + "'s #eta",    50, -9  ,   9  , Jet0.eta()     , weight);
  histo_->fill(prename + name + "_phi_" + suffix,    prename + name + "'s #phi",    50, -3.5,   3.5, Jet0.phi()     , weight);
  
  name = "J1";  
  histo_->fill(prename + name + "_charge_" + suffix, prename + name + "'s charge",   5, -2.5,   2.5, Jet1.charge()  , weight);
  histo_->fill(prename + name + "_mass_" + suffix,   prename + name + "'s mass",    63,  0  , 252  , Jet1.mass()    , weight);
  histo_->fill(prename + name + "_trmass_" + suffix, prename + name + "'s trmass", 100,  0  , 400  , Jet1.p4().Mt() , weight);
  histo_->fill(prename + name + "_pt_" + suffix,     prename + name + "'s p_{t}",  140,  0  , 700  , Jet1.pt()      , weight);
  histo_->fill(prename + name + "_Y_" + suffix,      prename + name + "'s Y",       50, -5  ,   5  , Jet1.rapidity(), weight);
  histo_->fill(prename + name + "_eta_" + suffix,    prename + name + "'s #eta",    50, -9  ,   9  , Jet1.eta()     , weight);
  histo_->fill(prename + name + "_phi_" + suffix,    prename + name + "'s #phi",    50, -3.5,   3.5, Jet1.phi()     , weight);

  name = "JJ";
  TLorentzVector JJp4 = Jet0.p4() + Jet1.p4();
  double JJdeltaEta = Jet0.eta() - Jet1.eta();
  double JJdeltaPhi = physmath::deltaPhi(Jet0.phi(), Jet1.phi());
  double JJdeltaR = abs(physmath::deltaR(Jet0, Jet1));
  
  histo_->fill(prename + name + "_mass_" + suffix,     "Leading Jets' mass",       100,  0  , 4000  , JJp4.M()  , weight);
  histo_->fill(prename + name + "_trmass_" + suffix,   "Leading Jets' trmass",     100,  0  , 4500  , JJp4.Mt() , weight);
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
