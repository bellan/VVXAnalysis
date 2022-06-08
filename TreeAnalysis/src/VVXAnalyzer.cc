#include "VVXAnalysis/TreeAnalysis/interface/VVXAnalyzer.h"
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

Int_t VVXAnalyzer::cut() {
  
  return 1;
}

void VVXAnalyzer::analyze(){

  //theHistograms->fill("nZtoChLep"    ,"Number of Z->ll per event" , 7,0,7, genVBHelper_.ZtoChLep().size());
  //theHistograms->fill("nZtoNeutrinos","Number of Z->nn per event" , 7,0,7, genVBHelper_.ZtoNeutrinos().size());
  //theHistograms->fill("nWtoLep"      ,"Number of W->lnu per event", 7,0,7, genVBHelper_.WtoLep().size());
  //theHistograms->fill("nZtoQ"        ,"Number of Z->qq per event" , 7,0,7, genVBHelper_.ZtoQ().size());
  //theHistograms->fill("nWtoQ"        ,"Number of W->qq' per event", 7,0,7, genVBHelper_.WtoQ().size());

  //int nVBs = genVBHelper_.ZtoChLep().size() + genVBHelper_.ZtoNeutrinos().size() + genVBHelper_.WtoLep().size() + genVBHelper_.ZtoQ().size() + genVBHelper_.WtoQ().size();
  //theHistograms->fill("nVBs", "Number of VB per event", 7,0,7, nVBs);

  theHistograms->fill("ZZ4l_mass"      , "ZZ4l mass"  , 15 , 0, 1500, ZZ->mass(), theWeight);
  theHistograms->fill("ZZ4l_mass_u"    , "ZZ4l mass"  , 15 , 0, 1500, ZZ->mass());
  int ZZsumid = abs(ZZ->first().daughter(0).id())+abs(ZZ->second().daughter(0).id());
  if     (ZZsumid == 22){
    theHistograms->fill("ZZ4e_mass"    , "ZZ4e mass"  , 15 , 0, 1500, ZZ->mass(), theWeight);
    theHistograms->fill("ZZ4e_mass_u"  , "ZZ4e mass"  , 15 , 0, 1500, ZZ->mass());
  }
  else if(ZZsumid == 24){
    theHistograms->fill("ZZ2e2m_mass"  , "ZZ2e2m mass", 15 , 0, 1500, ZZ->mass(), theWeight);
    theHistograms->fill("ZZ2e2m_mass_u", "ZZ2e2m mass", 15 , 0, 1500, ZZ->mass());
  }
  else if(ZZsumid == 26){
    theHistograms->fill("ZZ4m_mass"    , "ZZ4m mass"  , 15 , 0, 1500, ZZ->mass(), theWeight);
    theHistograms->fill("ZZ4m_mass_u"  , "ZZ4m mass"  , 15 , 0, 1500, ZZ->mass());
  }

  theHistograms->fill("ZW3l_tmass"  , "ZW tmass", 15 , 0, 1500, ZW->p4().Mt(), theWeight);
  theHistograms->fill("ZW3l_tmass_u", "ZW tmass", 15 , 0, 1500, ZW->p4().Mt());
  
  int ZWsumid = abs(ZW->first().daughter(0).id())+abs(ZW->second().daughter(0).id());
  if     (ZWsumid == 22){
    theHistograms->fill("ZW3e_tmass"      , "ZW3e tmass"  , 15 , 0, 1500, ZW->p4().Mt(), theWeight);
    theHistograms->fill("ZW3e_tmass"      , "ZW3e tmass"  , 15 , 0, 1500, ZW->p4().Mt());
  }
  else if(ZWsumid == 24){
    if     (abs(ZW->first().daughter(0).id()) == 11){
      theHistograms->fill("ZZ2e1m_tmass"  , "ZW2e1m tmass", 15 , 0, 1500, ZW->p4().Mt(), theWeight);
      theHistograms->fill("ZZ2e1m_tmass_u", "ZW2e1m tmass", 15 , 0, 1500, ZW->p4().Mt());
    }
    else if(abs(ZW->first().daughter(0).id()) == 13){
      theHistograms->fill("ZZ2m1e_tmass"  , "ZW2m1e tmass", 15 , 0, 1500, ZW->p4().Mt(), theWeight);
      theHistograms->fill("ZZ2m1e_tmass_u", "ZW2m1e tmass", 15 , 0, 1500, ZW->p4().Mt());
    }
  }
  else if(ZWsumid == 26){
    theHistograms->fill("ZW3m_tmass"      , "ZW4m tmass"  , 15 , 0, 1500, ZW->p4().Mt(), theWeight);
    theHistograms->fill("ZW3m_tmass_u"    , "ZW4m tmass"  , 15 , 0, 1500, ZW->p4().Mt());
  }
  
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
		
    //Kinematic selection
    if(ph.pt() < 20) continue;
    float ph_aeta = fabs(ph.eta());
    if(ph_aeta > 2.4) continue;
    if(ph_aeta > 1.4442 && ph_aeta < 1.566) continue;
		
    //Electrons and muons matching
    bool match = false;
    for(const Lepton lep : leptons){
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
  
  
  theHistograms->fill("AAA_cuts"  , "Cuts weighted"  , CUT_LAYOUT, 1, theWeight);
  theHistograms->fill("AAA_cuts_u", "Cuts unweighted", CUT_LAYOUT, 1);
  bool haveZVlep = (ZZ && ZZ->pt() > 1.) || (ZW && ZW->pt() > 1.);
  if(haveZVlep){
    theHistograms->fill("AAA_cuts"  , "Cuts weighted"  , CUT_LAYOUT, 2, theWeight);
    theHistograms->fill("AAA_cuts_u", "Cuts unweighted", CUT_LAYOUT, 2);
  }
  bool have2l2j = (muons->size()+electrons->size()==2) && (jets->size()==2 || jetsAK8->size()>=1);
  if(have2l2j){
    theHistograms->fill("AAA_cuts"  , "Cuts weighted"  , CUT_LAYOUT, 3, theWeight);
    theHistograms->fill("AAA_cuts_u", "Cuts unweighted", CUT_LAYOUT, 3);
  }
  if(kinPhotons.size() >= 1){
    theHistograms->fill("AAA_cuts"  , "Cuts weighted"  , CUT_LAYOUT, 4, theWeight);
    theHistograms->fill("AAA_cuts_u", "Cuts unweighted", CUT_LAYOUT, 4);
  }
  bool haveGoodPhoton = goodPhotons.size() >= 1;
  if(haveGoodPhoton){
    theHistograms->fill("AAA_cuts"  , "Cuts weighted"  , CUT_LAYOUT, 5, theWeight);
    theHistograms->fill("AAA_cuts_u", "Cuts unweighted", CUT_LAYOUT, 5);
  }
  if(haveZVlep && haveGoodPhoton){
    theHistograms->fill("AAA_cuts"  , "Cuts weighted"  , CUT_LAYOUT, 6, theWeight);
    theHistograms->fill("AAA_cuts_u", "Cuts unweighted", CUT_LAYOUT, 6);
  }
}

void VVXAnalyzer::end(TFile& fout){
  TH1* hCuts   = theHistograms->get("AAA_cuts");
  TH1* hCuts_u = theHistograms->get("AAA_cuts_u");
  for(TH1* h : {hCuts, hCuts_u}){
    TAxis* x = h->GetXaxis();
    x->SetBinLabel(1, "All");
    x->SetBinLabel(2, "ZZ || ZW");
    x->SetBinLabel(3, "2l2j || 2l1J");
    x->SetBinLabel(4, "#gamma kin");
    x->SetBinLabel(5, "#gamma IDloose");
    x->SetBinLabel(6, "ZZ/ZW && #gamma IDloose");
  }
}










