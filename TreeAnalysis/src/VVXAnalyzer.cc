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

  theHistograms->fill("ZZ4l_mass", "ZZ mass", 15 , 0, 1500, ZZ->mass(), theWeight);
  int ZZsumid = abs(ZZ->first().daughter(0).id())+abs(ZZ->second().daughter(0).id());
  if(ZZsumid == 22)      theHistograms->fill("ZZ4e_mass"  , "ZZ mass", 15 , 0, 1500, ZZ->mass(),theWeight);
  else if(ZZsumid == 23) theHistograms->fill("ZZ2e2m_mass", "ZZ mass", 15 , 0, 1500, ZZ->mass(),theWeight);
  else if(ZZsumid == 26) theHistograms->fill("ZZ4m_mass"  , "ZZ mass", 15 , 0, 1500, ZZ->mass(),theWeight);


  theHistograms->fill("ZW3l_tmass", "ZW tmass", 15 , 0, 1500, ZW->p4().Mt(),theWeight);
  
  int ZWsumid = abs(ZW->first().daughter(0).id())+abs(ZW->second().daughter(0).id());
  if(ZWsumid == 22)                                   theHistograms->fill("ZW3e_tmass"  , "ZW tmass", 15 , 0, 1500, ZW->p4().Mt(),theWeight);
  else if(ZWsumid == 23){
    if(abs(ZW->first().daughter(0).id()) == 11)       theHistograms->fill("ZZ2e1m_tmass", "ZW tmass", 15 , 0, 1500, ZW->p4().Mt(),theWeight);
    else if (abs(ZW->first().daughter(0).id()) == 13) theHistograms->fill("ZZ2m1e_tmass", "ZW tmass", 15 , 0, 1500, ZW->p4().Mt(),theWeight);
  }
  else if(ZWsumid == 26)                              theHistograms->fill("ZW3m_tmass"  , "ZW tmass", 15 , 0, 1500, ZW->p4().Mt(),theWeight);
    


  
}












