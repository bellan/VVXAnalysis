#include "VVXAnalysis/TreeAnalysis/interface/VZGammaAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Colours.h"

#include "VVXAnalysis/Commons/interface/GenVBHelper.h"


#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;
using namespace colour;

using namespace phys;

#define CUT_LAYOUT 6,0.5,6.5

Int_t VZGammaAnalyzer::cut() {
  
  return 1;
}

void VZGammaAnalyzer::analyze(){

  theHistograms->fill("nZtoChLep"    ,"Number of Z->ll per event" , 7,0,7, genVBHelper_.ZtoChLep().size());
  theHistograms->fill("nZtoNeutrinos","Number of Z->nn per event" , 7,0,7, genVBHelper_.ZtoNeutrinos().size());
  theHistograms->fill("nWtoLep"      ,"Number of W->lnu per event", 7,0,7, genVBHelper_.WtoLep().size());
  theHistograms->fill("nZtoQ"        ,"Number of Z->qq per event" , 7,0,7, genVBHelper_.ZtoQ().size());
  theHistograms->fill("nWtoQ"        ,"Number of W->qq' per event", 7,0,7, genVBHelper_.WtoQ().size());


  int nVBs = genVBHelper_.ZtoChLep().size() + genVBHelper_.ZtoNeutrinos().size() + genVBHelper_.WtoLep().size() + genVBHelper_.ZtoQ().size() + genVBHelper_.WtoQ().size();
  theHistograms->fill("nVBs", "Number of VB per event", 7,0,7, nVBs);


  
  
  std::vector<phys::Lepton> leptons;
  if     (region_ >= SR4P && region_ <= CR4P_1F)  //(ZZ && ZZ->pt() > 1.){
    leptons.insert(leptons.end(), {
    	  ZZ->first().daughter(0), 
    	  ZZ->first().daughter(1), 
    	  ZZ->second().daughter(0), 
    	  ZZ->second().daughter(1)
        });
  else if(ZW && ZW->pt() > 1.)
    leptons.insert(leptons.end(), {
  	  ZW->first().daughter(0), 
  	  ZW->first().daughter(1), 
  	  ZW->second().daughter(0)
	});
  else if(ZL && ZL->first.pt() > 1.)
    leptons.insert(leptons.end(), {
  	  ZL->first.daughter(0),
  	  ZL->first.daughter(1),
  	  ZL->second
	});
  
  
  std::vector<phys::Photon> kinPhotons ;
  std::vector<phys::Photon> goodPhotons;
  for(auto ph : *photons){
    //Pixel seed and electron veto
    if(ph.hasPixelSeed() || !ph.passElectronVeto()) continue;
		
    //KinemaAtic selection
    if(ph.pt() < 20) continue;
    float ph_aeta = fabs(ph.eta());
    if(ph_aeta > 2.4) continue;
    if(ph_aeta > 1.4442 && ph_aeta < 1.566) continue;
		
    //Electrons and muons maAtching
    bool match = false;
    for(const Lepton& lep : leptons){
      if(physmath::deltaR(ph,lep) < 0.3){
	match = true;
	break;
      }
    }
    if(match) continue;

    kinPhotons.push_back(ph);
    if(ph.cutBasedIDLoose())
      goodPhotons.push_back(ph);
  }
  

}











