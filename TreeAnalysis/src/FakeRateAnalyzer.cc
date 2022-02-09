#include "VVXAnalysis/TreeAnalysis/interface/FakeRateAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include <boost/format.hpp>
#include "VVXAnalysis/Commons/interface/StringTools.h"
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;

using namespace phys;


Int_t FakeRateAnalyzer::cut() {
  return 1;
}

void FakeRateAnalyzer::addOptions(){
  //if(!ZZ->isValid()) return;
  
  //cout << "SF before: " <<  ZZ1.fakeRateSF() << std::endl;
  //cout << ZZ1.firstPtr() ->daughterPtr(0)->fakeRateSF() << " " << lepSF.fakeRateScaleFactor(ZZ1.first() .daughter(0)).first << " " 
  //     << ZZ1.firstPtr() ->daughterPtr(1)->fakeRateSF() << " " << lepSF.fakeRateScaleFactor(ZZ1.first() .daughter(1)).first << " " 
  //     << ZZ1.secondPtr()->daughterPtr(0)->fakeRateSF() << " " << lepSF.fakeRateScaleFactor(ZZ1.second().daughter(0)).first << " " 
  //     << ZZ1.secondPtr()->daughterPtr(1)->fakeRateSF() << " " << lepSF.fakeRateScaleFactor(ZZ1.second().daughter(1)).first << endl;

  //ZZ->firstPtr()->daughterPtr(0)->setFakeRateSF (lepSF.fakeRateScaleFactor(ZZ->first().daughter(0)));
  //ZZ->firstPtr()->daughterPtr(1)->setFakeRateSF (lepSF.fakeRateScaleFactor(ZZ->first().daughter(1)));
  //ZZ->secondPtr()->daughterPtr(0)->setFakeRateSF(lepSF.fakeRateScaleFactor(ZZ->second().daughter(0)));
  //ZZ->secondPtr()->daughterPtr(1)->setFakeRateSF(lepSF.fakeRateScaleFactor(ZZ->second().daughter(1)));

  //cout << " after: " <<  ZZ->fakeRateSF() << endl;
}

void FakeRateAnalyzer::ZZplots(int id){

  if(ZZ->id() != id && id != -1) return; // -1 here means generic 4l final state

  std::string decay  = "4l";
  std::string decay1 = "l";
  std::string decay2 = "l";
  if      (id == 52) {decay = "4m"  ; decay1 = "m"; decay1 = "m";}
  else if (id == 48) {decay = "2e2m"; decay1 = "e"; decay2 = "m";}
  else if (id == 44) {decay = "4e"  ; decay1 = "e"; decay2 = "e";}


  theHistograms->fill(std::string("ZZTo")+decay+std::string("_mZ1To2")+decay1, std::string("Invariant mass of Z_{1}#rightarrow 2")+decay1,  15, 60,  120, ZZ->first().mass() , theWeight); 
  theHistograms->fill(std::string("ZZTo")+decay+std::string("_mZ2To2")+decay2, std::string("Invariant mass of Z_{2}#rightarrow 2")+decay2,  15, 60,  120, ZZ->second().mass(), theWeight); 
  theHistograms->fill(std::string("ZZTo")+decay+std::string("_mZZTo") +decay , std::string("Invariant mass of ZZ#rightarrow ")    +decay ,  40,  0, 1000, ZZ->mass()         , theWeight); 

  theHistograms->fill(std::string("ZZTo")+decay+"_nJets"       , "Number of jets (|#eta|<4.7 and p_T > 30 GeV)"        , 10, 0, 10, jets->size(), theWeight); 
  theHistograms->fill(std::string("ZZTo")+decay+"_nCentralJets", "Number of central jets (|#eta|<2.5 and p_T > 30 GeV)", 10, 0, 10, centralJets->size(), theWeight); 

  if(jets->size() >= 2)
    theHistograms->fill(std::string("ZZTo")+decay+"_deltaEtaJJ", "#Delta #eta(j,j) between the two most energetic jets",  10, 0, 8, fabs(jets->at(0).eta() - jets->at(1).eta()), theWeight); 
  
  
  if(centralJets->size() >= 2){
    theHistograms->fill(std::string("ZZTo")+decay+"_deltaEtaJJcentral", "#Delta #eta(j,j) between the two most energetyc central jets",  10, 0, 8, fabs(centralJets->at(0).eta() - centralJets->at(1).eta()), theWeight); 
    theHistograms->fill(std::string("ZZTo")+decay+"_mJJ", "m_{jj}",  20, 0, 1000, (centralJets->at(0).p4() + centralJets->at(1).p4()).M(), theWeight); 
  }
  theHistograms->fill(std::string("ZZTo")+decay+"_nExtraMuons"    , "Number of extra muons in the event"    , 10, 0, 10, muons->size(), theWeight); 
  theHistograms->fill(std::string("ZZTo")+decay+"_nExtraElectrons", "Number of extra electrons in the event", 10, 0, 10, electrons->size(), theWeight); 
  theHistograms->fill(std::string("ZZTo")+decay+"_nExtraLeptons"  , "Number of extra leptons in the event"  , 10, 0, 10, muons->size()+electrons->size(), theWeight); 
}


void FakeRateAnalyzer::analyze(){
  
  theHistograms->fill("fakeRateWeight", "FakeRate",  100, -2, 2, ZZ->fakeRateSF() , 1); 
  if(ZZ->numberOfGoodGrandDaughters() == 3) theHistograms->fill("fakeRateWeight_3p1F", "FakeRate 3P1F",  100, -2, 2, ZZ->fakeRateSF() , 1); 
  if(ZZ->numberOfGoodGrandDaughters() == 2) theHistograms->fill("fakeRateWeight_2p2F", "FakeRate 2P2F",  100, -2, 2, ZZ->fakeRateSF() , 1); 
  
  // Some basic plots on ZZ
  ZZplots();   // ZZ --> 4l
  ZZplots(52); // ZZ --> 4m
  ZZplots(48); // ZZ --> 2e2m
  ZZplots(44); // ZZ --> 4e
  
    
  if(ZL->first.isValid()){
    
    Float_t mt = sqrt(2*met->pt()*(ZL->second.pt())*(1-cos(physmath::deltaPhi(ZL->second.phi(),met->phi()))));

    if( abs(ZL->second.id()) == 13)   theHistograms->fill("MT_mu","Transverse mass",200, 0, 200, mt, theWeight);
    else if( abs(ZL->second.id()) == 11)   theHistograms->fill("MT_ele","Transverse mass",200, 0, 200, mt, theWeight);
    else std::cout<<"ERROR check lepton id"<<std::endl;
    

    std::string eventstr=std::to_string(run)+":"+std::to_string(lumiBlock)+":"+std::to_string(event);
    
    theWeight=theMCInfo.weight(); //TheWeight have problems


        if((met->pt()<25) && (mt < 30)) { //ZZ
	  //   if(met->pt()<25) {  //HZZ4L

      
    

      std::bitset<16> trigger(triggerWord);
      if(ZL->first.isValid() == 1 && trigger.test(9)){
 
	//	std::vector<double> xbins;xbins += 0,10,20,30,40,50,60,70,80,100,200;
	// std::vector<double> xbins;xbins += 5,7,10,20,30,40,50,80; //HZZ?
	// std::vector<double> xbins;xbins += 5,10,30,60,200; //UW binning
	// std::vector<double> xbins;xbins += 5,7,10,20,30,80;
  

	std::string eventstr=std::to_string(run)+":"+std::to_string(lumiBlock)+":"+std::to_string(event);
	
	
	Float_t	m3L = sqrt(((ZL->first.p4())+(ZL->second.p4()))*((ZL->first.p4())+(ZL->second.p4())));
	
	std::string  eventStrInf = strtool::sRound(m3L)+":"+strtool::sRound(ZL->first.mass())+":"+strtool::sRound(ZL->second.pt());

	Float_t pt = ZL->second.pt();
	   
	if(abs(ZL->second.id()) == 13){
	  

	  if( pt > 200 ) pt = 199;
	  theHistograms->fill("ZL_Muon_pt","Pt of the loose lepton",150,0,300,ZL->second.pt(),theWeight);	  
	  theHistograms->fill("ZL_Muon_Iso","Iso of the loose lepton",50,0,0.5,ZL->second.pfCombRelIso(),theWeight);
	  theHistograms->fill("ZL_Muon_sip","Sip of the loose lepton",50,0,4.5,ZL->second.sip(),theWeight);	  
	  
	  if(abs(ZL->first.daughterPtr(0)->id()) == 13) { eventsD_mmm.push_back(eventstr); eventsStr.push_back(eventstr+":mmm:"+eventStrInf+":"+std::to_string(ZL->second.passFullSel()));} //eventcouting
	  else { eventsD_eem.push_back(eventstr); eventsStr.push_back(eventstr+":eem:"+eventStrInf+":"+std::to_string(ZL->second.passFullSel()));} 
	  
	  if(fabs(ZL->second.eta()) < 1.2)   theHistograms->fill("FakeRate_denom_muons_barrel_pt","Total number of soft leptons in the barrel",ybins_mu,pt,theWeight);
	  else                        theHistograms->fill("FakeRate_denom_muons_endcap_pt","Total number of soft leptons in the endcaps",ybins_mu,pt,theWeight);


	  theHistograms->fill("FakeRate_denom_muons","Total number of soft leptons in the barrel",xbins_mu,ybins_mu,ZL->second.eta(),pt,theWeight);
	  
	  if(ZL->second.passFullSel()){
	    
	    if(abs(ZL->first.daughterPtr(0)->id()) == 13)  eventsN_mmm.push_back(eventstr); //eventcounting
	    else  eventsN_eem.push_back(eventstr);
	    
	    if(fabs(ZL->second.eta()) < 1.2) theHistograms->fill("FakeRate_num_muons_barrel_pt","Number of tight leptons in the barrel",ybins_mu,pt,theWeight);
	    else theHistograms->fill("FakeRate_num_muons_endcap_pt","Number of tight leptons in the endcaps",ybins_mu,pt,theWeight);

	    theHistograms->fill("FakeRate_num_muons","Total number of soft leptons in the barrel",xbins_mu,ybins_mu,ZL->second.eta(),pt,theWeight);
   	  }
	}
	
	if(abs(ZL->second.id()) == 11){
	  
	  if( pt > 80 ) pt = 79;

	  theHistograms->fill("ZL_Ele_pt","Pt of the loose lepton",150,0,300,pt,theWeight);	  
	  theHistograms->fill("ZL_Ele_Iso","Iso of the loose lepton",50,0,0.5,ZL->second.pfCombRelIso(),theWeight);
	  theHistograms->fill("ZL_Ele_sip","Sip of the loose lepton",50,0,4.5,ZL->second.sip(),theWeight);
	   
	  if(abs(ZL->first.daughterPtr(0)->id()) == 13) { eventsD_mme.push_back(eventstr); eventsStr.push_back(eventstr+":mme:"+eventStrInf+":"+std::to_string(ZL->second.passFullSel()));} //eventcounting
	  else  {eventsD_eee.push_back(eventstr); eventsStr.push_back(eventstr+":mmm:"+eventStrInf+":"+std::to_string(ZL->second.passFullSel()));}
	  
	  if(fabs(ZL->second.eta()) < 1.45)   theHistograms->fill("FakeRate_denom_electrons_barrel_pt","Total number of soft leptons in the barrel",ybins_ele,pt,theWeight);
	  else                         theHistograms->fill("FakeRate_denom_electrons_endcap_pt","Total number of soft leptons in the endcaps",ybins_ele,pt,theWeight);
	  
	  theHistograms->fill("FakeRate_denom_electrons","Total number of soft leptons in the barrel",xbins_ele,ybins_ele,ZL->second.eta(),pt,theWeight);

	  if(ZL->second.passFullSel()){
	    
	    if(abs(ZL->first.daughterPtr(0)->id()) == 13)  eventsN_mme.push_back(eventstr); //eventcounting
	    else  eventsN_eee.push_back(eventstr);
	    
	    if(fabs(ZL->second.eta()) < 1.45) theHistograms->fill("FakeRate_num_electrons_barrel_pt","Number of tight leptons in the barrel",ybins_ele,pt,theWeight);
	    else                       theHistograms->fill("FakeRate_num_electrons_endcap_pt","Number of tight leptons in the endcaps",ybins_ele,pt,theWeight);	
    
	    theHistograms->fill("FakeRate_num_electrons","Total number of soft leptons in the barrel",xbins_ele,ybins_ele,ZL->second.eta(),pt,theWeight);
	  }
	}
      } 
    }
  }
}

void FakeRateAnalyzer::begin() {


  ybins_ele +=  7, 8, 9, 10, 11, 12, 14, 16, 18, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80;
  xbins_ele +=  0,0.8,1.2,2.,2.5;
  
  ybins_mu += 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0, 30.0, 70.0, 200.0;
  xbins_mu  += 0,1.2,2.4;
  
}

void FakeRateAnalyzer::end( TFile &){
  // theHistograms->clone("FakeRate_num_muons_barrel_pt","FakeRate_muons_barrel_pt");
  // theHistograms->get("FakeRate_muons_barrel_pt")->Divide(theHistograms->get("FakeRate_denom_muons_barrel_pt"));
  
  // theHistograms->clone("FakeRate_num_muons_endcap_pt","FakeRate_muons_endcap_pt");
  // theHistograms->get("FakeRate_muons_endcap_pt")->Divide(theHistograms->get("FakeRate_denom_muons_endcap_pt"));
  
  // theHistograms->clone("FakeRate_num_electrons_barrel_pt","FakeRate_electrons_barrel_pt");
  // theHistograms->get("FakeRate_electrons_barrel_pt")->Divide(theHistograms->get("FakeRate_denom_electrons_barrel_pt"));
  
  // theHistograms->clone("FakeRate_num_electrons_endcap_pt","FakeRate_electrons_endcap_pt");
  // theHistograms->get("FakeRate_electrons_endcap_pt")->Divide(theHistograms->get("FakeRate_denom_electrons_endcap_pt"));
  
  std::sort (eventsStr.begin(), eventsStr.end(), strtool::sortEvents);
  
  std::cout<<"eee:"<<std::endl;  
  std::cout<<"Num: "<< eventsN_eee.size()<<std::endl;
  std::cout<<"Den: "<<eventsD_eee.size()<<std::endl;
  std::cout<<"eem:"<<std::endl;
  std::cout<<"Num: "<<eventsN_eem.size()<<std::endl;
  std::cout<<"Den: "<<eventsD_eem.size()<<std::endl;
  std::cout<<"mme:"<<std::endl;
  std::cout<<"Num: "<<eventsN_mme.size()<<std::endl;
  std::cout<<"Den: "<<eventsD_mme.size()<<std::endl;
  std::cout<<"mmm:"<<std::endl;
  std::cout<<"Num: "<<eventsN_mmm.size()<<std::endl;
  std::cout<<"Den: "<<eventsD_mmm.size()<<std::endl;
  std::cout<<"Sum "<<eventsD_eee.size()+eventsD_eem.size()+eventsD_mme.size()+eventsD_mmm.size()<<std::endl;
  //  for (std::vector<std::string>::iterator it = eventsStr.begin() ; it != eventsStr.end(); ++it) std::cout<<*it<<std::endl;

}
