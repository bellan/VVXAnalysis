#include "VVXAnalysis/TreeAnalysis/interface/FakeRateAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"


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




  theWeight = theMCInfo.weight(); 

  theHistograms.fill(std::string("ZZTo")+decay+std::string("_mZ1To2")+decay1, std::string("Invariant mass of Z_{1}#rightarrow 2")+decay1,  15, 60,  120, ZZ->first().mass() , theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+std::string("_mZ2To2")+decay2, std::string("Invariant mass of Z_{2}#rightarrow 2")+decay2,  15, 60,  120, ZZ->second().mass(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+std::string("_mZZTo") +decay , std::string("Invariant mass of ZZ#rightarrow ")    +decay ,  40,  0, 1000, ZZ->mass()         , theWeight); 

  theHistograms.fill(std::string("ZZTo")+decay+"_nJets"       , "Number of jets (|#eta|<4.7 and p_T > 30 GeV)"        , 10, 0, 10, jets->size(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_nCentralJets", "Number of central jets (|#eta|<2.5 and p_T > 30 GeV)", 10, 0, 10, centralJets->size(), theWeight); 

  if(jets->size() >= 2)
    theHistograms.fill(std::string("ZZTo")+decay+"_deltaEtaJJ", "#Delta #eta(j,j) between the two most energetic jets",  10, 0, 8, fabs(jets->at(0).eta() - jets->at(1).eta()), theWeight); 
  
  
  if(centralJets->size() >= 2){
    theHistograms.fill(std::string("ZZTo")+decay+"_deltaEtaJJcentral", "#Delta #eta(j,j) between the two most energetyc central jets",  10, 0, 8, fabs(centralJets->at(0).eta() - centralJets->at(1).eta()), theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_mJJ", "m_{jj}",  20, 0, 1000, (centralJets->at(0).p4() + centralJets->at(1).p4()).M(), theWeight); 
  }

  theHistograms.fill(std::string("ZZTo")+decay+"_nExtraMuons"    , "Number of extra muons in the event"    , 10, 0, 10, muons->size(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_nExtraElectrons", "Number of extra electrons in the event", 10, 0, 10, electrons->size(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_nExtraLeptons"  , "Number of extra leptons in the event"  , 10, 0, 10, muons->size()+electrons->size(), theWeight); 

}


void FakeRateAnalyzer::analyze(){

  theHistograms.fill("fakeRateWeight", "FakeRate",  100, -2, 2, ZZ->fakeRateSF() , 1); 
  if(ZZ->numberOfGoodGrandDaughters() == 3) theHistograms.fill("fakeRateWeight_3p1F", "FakeRate 3P1F",  100, -2, 2, ZZ->fakeRateSF() , 1); 
  if(ZZ->numberOfGoodGrandDaughters() == 2) theHistograms.fill("fakeRateWeight_2p2F", "FakeRate 2P2F",  100, -2, 2, ZZ->fakeRateSF() , 1); 



  
  // Some basic plots on ZZ
  ZZplots();   // ZZ --> 4l
  ZZplots(52); // ZZ --> 4m
  ZZplots(48); // ZZ --> 2e2m
  ZZplots(44); // ZZ --> 4e
  

  
  theHistograms.fill("NumberOfZLCandidates","Number of ZL candidates",10, 0, 10, ZLCand->size(), theWeight);
  theHistograms.fill("NumberOfZLCandidatesVsMET","Number of ZL candidates vs MET",200,0,800,10, 0, 10, met->pt(), ZLCand->size(), theWeight);
  theHistograms.fill("NumberOfZLCandidatesVsRegion","Number of ZL candidates vs region type",4,0,4,10, 0, 10,
		     regionWord.test(3) ? 0 : (regionWord.test(22) ? 1 : (regionWord.test(23) ? 2 : 3 ))
		     , ZLCand->size(), theWeight);
  

  theHistograms.fill("NumberOfZL","Number of selected ZL",10, 0, 10, ZL->size(), theWeight);
  theHistograms.fill("NumberOfZLVsMET","Number of selected ZL vs MET",200,0,800,10, 0, 10, met->pt(), ZL->size(), theWeight);
  theHistograms.fill("NumberOfZLVsRegion","Number of selected ZL vs region type",4,0,4,10, 0, 10,
		     regionWord.test(3) ? 0 : (regionWord.test(22) ? 1 : (regionWord.test(23) ? 2 : 3 ))
		     , ZL->size(), theWeight);
  
  

  if(ZL->size()){

  Float_t mt = sqrt(2*met->pt()*ZL->front().second.pt()*(1-cos(physmath::deltaPhi(ZL->front().second.phi(),met->phi()))));
  

  if( abs(ZL->front().second.id()) == 13)   theHistograms.fill("MT_mu","Transverse mass",200, 0, 200, mt, theWeight);
  else if( abs(ZL->front().second.id()) == 11)   theHistograms.fill("MT_ele","Transverse mass",200, 0, 200, mt, theWeight);
  else std::cout<<"ERROR check elp id"<<std::endl;


  std::string eventstr=std::to_string(run)+":"+std::to_string(lumiBlock)+":"+std::to_string(event);
  


  if((met->pt()<25) && (ZL->front().second.pt()>10)) {
  
    //    std::vector<double> xbins;xbins += 5,7,10,20,30,40,50,80;
    //    std::vector<double> xbins;xbins += 5,7,10,20,30,80;  

    std::vector<double> xbins;xbins += 0,30,60,200;
	if(abs(ZL->front().second.id()) == 13){
	  
	   
	  if(fabs(ZL->front().second.eta()) < 1.2)   theHistograms.fill("FakeRate_denom_muons_NoMt_barrel_pt","Total number of soft leptons in the barrel",xbins,ZL->front().second.pt(),theWeight);
	  else                        theHistograms.fill("FakeRate_denom_muons_NoMt_endcap_pt","Total number of soft leptons in the endcaps",xbins,ZL->front().second.pt(),theWeight);
	  
	  if(ZL->front().second.passFullSel()){
	     
	    
	    if(fabs(ZL->front().second.eta()) < 1.2) theHistograms.fill("FakeRate_num_muons_NoMt_barrel_pt","Number of tight leptons in the barrel",xbins,ZL->front().second.pt(),theWeight);
	    else theHistograms.fill("FakeRate_num_muons_NoMt_endcap_pt","Number of tight leptons in the endcaps",xbins,ZL->front().second.pt(),theWeight);
	    
	  }
	}
  }


  if((met->pt()<25) && (mt < 30) && (ZL->front().second.pt()>10)) {
    // if(met->pt()<25) {  //HZZ4L

        if ((ZL->front().first.daughterPtr(0)->pt() <10) || (ZL->front().first.daughterPtr(1)->pt() <10)) std::cout<<ZL->front().first.daughterPtr(0)->pt()<<" "<<(ZL->front().first.daughterPtr(1)->pt())<<std::endl;
     
 
    if(ZL->front().second.sip()>4) std::cout<<"bad sip "<<ZL->front().second.sip()<<std::endl;


      theHistograms.fill("MET25_NumberOfZLCandidates","Number of ZL candidates MET<25",10, 0, 10, ZLCand->size(), theWeight);
      theHistograms.fill("MET25_NumberOfZLCandidatesVsMET","Number of ZL candidates vs MET MET<25",200,0,800,10, 0, 10, met->pt(), ZLCand->size(), theWeight);
      theHistograms.fill("MET25_NumberOfZLCandidatesVsRegion","Number of ZL candidates vs region type MET<25",4,0,4,10, 0, 10,
			 regionWord.test(3) ? 0 : (regionWord.test(22) ? 1 : (regionWord.test(23) ? 2 : 3 ))
			 , ZLCand->size(), theWeight);
      
      theHistograms.fill("MET25_NumberOfZL","Number of selected ZL MET<25",10, 0, 10, ZL->size(), theWeight);
      theHistograms.fill("MET25_NumberOfZLVsMET","Number of selected ZL vs MET MET<25",200,0,800,10, 0, 10, met->pt(), ZL->size(), theWeight);
      theHistograms.fill("MET25_NumberOfZLVsRegion","Number of selected ZL vs region type MET<25",4,0,4,10, 0, 10,
			 regionWord.test(3) ? 0 : (regionWord.test(22) ? 1 : (regionWord.test(23) ? 2 : 3 ))
			 , ZL->size(), theWeight);
      
      std::bitset<16> trigger(triggerWord);
      if(ZL->size() == 1 && trigger.test(9)){
	//std::vector<double> xbins;xbins += 5,7,10,20,30,40,50,80;
	std::vector<double> xbins;xbins += 0,30,60,200;
	std::string eventstr=std::to_string(run)+":"+std::to_string(lumiBlock)+":"+std::to_string(event);
	
	


  	if(abs(ZL->front().second.id()) == 13){

	  theHistograms.fill("ZL_Muon_pt","Pt of the loose lepton",150,0,300,ZL->front().second.pt(),theWeight);	  
	  theHistograms.fill("ZL_Muon_Iso","Iso of the loose lepton",50,0,0.5,ZL->front().second.pfCombRelIso(),theWeight);
	  theHistograms.fill("ZL_Muon_sip","Sip of the loose lepton",50,0,4.5,ZL->front().second.sip(),theWeight);	  

	  if(abs(ZL->front().first.daughterPtr(0)->id()) == 13)  eventsD_mmm.push_back(eventstr);  //eventcouting
	  else  eventsD_eem.push_back(eventstr);
	  
	  if(fabs(ZL->front().second.eta()) < 1.2)   theHistograms.fill("FakeRate_denom_muons_barrel_pt","Total number of soft leptons in the barrel",xbins,ZL->front().second.pt(),theWeight);
	  else                        theHistograms.fill("FakeRate_denom_muons_endcap_pt","Total number of soft leptons in the endcaps",xbins,ZL->front().second.pt(),theWeight);
	  
	  //FIXME
	  //theHistograms.fill("FakeRate_denom_muons_barrel_pt","Total number of soft leptons in the barrel",xbins,ZL->front().second.pt(),theWeight);
	  //theHistograms.fill("FakeRate_denom_muons_endcap_pt","Total number of soft leptons in the endcaps",xbins,ZL->front().second.pt(),theWeight);


	  if(ZL->front().second.passFullSel()){
	    
	    if(abs(ZL->front().first.daughterPtr(0)->id()) == 13)  eventsN_mmm.push_back(eventstr); //eventcounting
	    else  eventsN_eem.push_back(eventstr);
	    
	    
	     if(fabs(ZL->front().second.eta()) < 1.2) theHistograms.fill("FakeRate_num_muons_barrel_pt","Number of tight leptons in the barrel",xbins,ZL->front().second.pt(),theWeight);
	     else theHistograms.fill("FakeRate_num_muons_endcap_pt","Number of tight leptons in the endcaps",xbins,ZL->front().second.pt(),theWeight);
	    //FIXME
	    //theHistograms.fill("FakeRate_num_muons_barrel_pt","Number of tight leptons in the barrel",xbins,ZL->front().second.pt(),theWeight);
	    //theHistograms.fill("FakeRate_num_muons_endcap_pt","Number of tight leptons in the endcaps",xbins,ZL->front().second.pt(),theWeight);

	    
	  }
	}
	
	if(abs(ZL->front().second.id()) == 11){
	  

	  theHistograms.fill("ZL_Ele_pt","Pt of the loose lepton",150,0,300,ZL->front().second.pt(),theWeight);	  
	  theHistograms.fill("ZL_Ele_Iso","Iso of the loose lepton",50,0,0.5,ZL->front().second.pfCombRelIso(),theWeight);
	  theHistograms.fill("ZL_Ele_sip","Sip of the loose lepton",50,0,4.5,ZL->front().second.sip(),theWeight);


	  if(abs(ZL->front().first.daughterPtr(0)->id()) == 13)  eventsD_mme.push_back(eventstr); //eventcounting
	  else  eventsD_eee.push_back(eventstr);
	  
	  if(fabs(ZL->front().second.eta()) < 1.45)   theHistograms.fill("FakeRate_denom_electrons_barrel_pt","Total number of soft leptons in the barrel",xbins,ZL->front().second.pt(),theWeight);
	  else                         theHistograms.fill("FakeRate_denom_electrons_endcap_pt","Total number of soft leptons in the endcaps",xbins,ZL->front().second.pt(),theWeight);

	  //theHistograms.fill("FakeRate_denom_electrons_barrel_pt","Total number of soft leptons in the barrel",xbins,ZL->front().second.pt(),theWeight);
	  //theHistograms.fill("FakeRate_denom_electrons_endcap_pt","Total number of soft leptons in the endcaps",xbins,ZL->front().second.pt(),theWeight);
	  if(ZL->front().second.passFullSel()){
	    
	    if(abs(ZL->front().first.daughterPtr(0)->id()) == 13)  eventsN_mme.push_back(eventstr); //eventcounting
	    else  eventsN_eee.push_back(eventstr);
	    
	    if(fabs(ZL->front().second.eta()) < 1.45) theHistograms.fill("FakeRate_num_electrons_barrel_pt","Number of tight leptons in the barrel",xbins,ZL->front().second.pt(),theWeight);
	    else                       theHistograms.fill("FakeRate_num_electrons_endcap_pt","Number of tight leptons in the endcaps",xbins,ZL->front().second.pt(),theWeight);
	    
	    //theHistograms.fill("FakeRate_num_electrons_barrel_pt","Number of tight leptons in the barrel",xbins,ZL->front().second.pt(),theWeight);
	    //theHistograms.fill("FakeRate_num_electrons_endcap_pt","Number of tight leptons in the endcaps",xbins,ZL->front().second.pt(),theWeight);
	  }
	}
      } 
  }
}
}

    void FakeRateAnalyzer::end( TFile &){
      // theHistograms.clone("FakeRate_num_muons_barrel_pt","FakeRate_muons_barrel_pt");
      // theHistograms.get("FakeRate_muons_barrel_pt")->Divide(theHistograms.get("FakeRate_denom_muons_barrel_pt"));
      
      // theHistograms.clone("FakeRate_num_muons_endcap_pt","FakeRate_muons_endcap_pt");
      // theHistograms.get("FakeRate_muons_endcap_pt")->Divide(theHistograms.get("FakeRate_denom_muons_endcap_pt"));
      
      // theHistograms.clone("FakeRate_num_electrons_barrel_pt","FakeRate_electrons_barrel_pt");
      // theHistograms.get("FakeRate_electrons_barrel_pt")->Divide(theHistograms.get("FakeRate_denom_electrons_barrel_pt"));
      
      // theHistograms.clone("FakeRate_num_electrons_endcap_pt","FakeRate_electrons_endcap_pt");
      // theHistograms.get("FakeRate_electrons_endcap_pt")->Divide(theHistograms.get("FakeRate_denom_electrons_endcap_pt"));

    
    std::cout<<"eee:"<<std::endl;

    std::cout<<"Num: "<< eventsN_eee.size()<<std::endl;
    //    for (std::vector<std::string>::iterator it = eventsN_eee.begin() ; it != eventsN_eee.end(); ++it) std::cout<<*it<<std::endl;
    std::cout<<"Den: "<<eventsD_eee.size()<<std::endl;

    //for (std::vector<std::string>::iterator it = eventsD_eee.begin() ; it != eventsD_eee.end(); ++it) std::cout<<*it<<std::endl;

    std::cout<<"eem:"<<std::endl;

    std::cout<<"Num: "<<eventsN_eem.size()<<std::endl;
    //for (std::vector<std::string>::iterator it = eventsN_eem.begin() ; it != eventsN_eem.end(); ++it) std::cout<<*it<<std::endl;
    std::cout<<"Den: "<<eventsD_eem.size()<<std::endl;
    //for (std::vector<std::string>::iterator it = eventsD_eem.begin() ; it != eventsD_eem.end(); ++it) std::cout<<*it<<std::endl;


    std::cout<<"mme:"<<std::endl;

    std::cout<<"Num: "<<eventsN_mme.size()<<std::endl;
    //for (std::vector<std::string>::iterator it = eventsN_mme.begin() ; it != eventsN_mme.end(); ++it) std::cout<<*it<<std::endl;
    std::cout<<"Den: "<<eventsD_mme.size()<<std::endl;
    //for (std::vector<std::string>::iterator it = eventsD_mme.begin() ; it != eventsD_mme.end(); ++it) std::cout<<*it<<std::endl;



    std::cout<<"mmm:"<<std::endl;

    std::cout<<"Num: "<<eventsN_mmm.size()<<std::endl;
    //for (std::vector<std::string>::iterator it = eventsN_mmm.begin() ; it != eventsN_mmm.end(); ++it) std::cout<<*it<<std::endl;
    std::cout<<"Den: "<<eventsD_mmm.size()<<std::endl;
    //for (std::vector<std::string>::iterator it = eventsD_mmm.begin() ; it != eventsD_mmm.end(); ++it) std::cout<<*it<<std::endl;
    


}
