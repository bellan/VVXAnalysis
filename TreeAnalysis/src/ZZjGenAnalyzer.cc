#include "VVXAnalysis/TreeAnalysis/interface/ZZjGenAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <string> 
#include "VVXAnalysis/Commons/interface/StringTools.h"                                                                

#define foreach BOOST_FOREACH

using std::cout;
using std::endl;

using namespace phys;
using namespace physmath;


void ZZjGenAnalyzer::ZZplots(std::string decay){

  //if(ZZ->id() != id && id != -1) return; // -1 here means generic 4l final state


  std::string eventstr=std::to_string(run)+":"+std::to_string(lumiBlock)+":"+std::to_string(event);
  std::string channel  = "";
  
  if (decay  == "4m"){
    channel  = "mmmm";
    events4mu.push_back(eventstr);
  }
  
  else if ( decay  == "2e2m"){
    channel  = "eemm";
    events2e2mu.push_back(eventstr);
  }
  else if ( decay == "4e"){
    channel  = "eeee";
    events4e.push_back(eventstr);
  }
  else if ( decay == "4l") {
    events4l.push_back(eventstr);
  }
  else {std::cout<<"wrong decay "<<decay<<std::endl; abort();}


  //  std::cout<<"Z0Mass "<<genVBParticles->at(0).mass()<<std::endl;
  //  std::cout<<"Z1Mass "<<genVBParticles->at(1).mass()<<std::endl;  
  m4L_gen = sqrt((genVBParticles->at(0).p4()+genVBParticles->at(1).p4())*(genVBParticles->at(0).p4()+genVBParticles->at(1).p4()));

  mZ1_gen = genVBParticles->at(0).mass();
  mZ2_gen = genVBParticles->at(1).mass();

  std::string eventStr = std::to_string(run)+":"+std::to_string(lumiBlock)+":"+std::to_string(event)+":"+channel+":"+strtool::sRound(m4L_gen)+":"+strtool::sRound(mZ1_gen)+":"+strtool::sRound(mZ2_gen)+":"+std::to_string(genJets->size());

  if(genJets->size()>0){ 
    eventStr+=":";
    eventStr+=strtool::sRound(genJets->at(0).pt());
  }
  else  eventStr+=":0.00";


  theHistograms->fill(std::string("ZZTo")+decay+"_Z0Mass"       , std::string("Z0  mass  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,mZ1_gen,theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z1Mass"       , std::string("Z1  mass  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,mZ2_gen,theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z0pt"         , std::string("Z0  p_{T}  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,genVBParticles->at(0).pt(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z1pt"         , std::string("Z1  p_{T}  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,genVBParticles->at(1).pt(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Mass"         , std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay            ,  200, 50,  1000,  m4L_gen, theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z0lep0_pt"    , std::string("pt of  Z0 lep0 of ZZ_{1}#rightarrow ")+decay            ,  300, 0,  300,genVBParticles->at(0).daughterPtr(0)->pt(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z0lep1_pt"    , std::string("pt of  Z0 lep1 of ZZ_{1}#rightarrow ")+decay            ,  300, 0,  300,genVBParticles->at(0).daughterPtr(1)->pt(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z1lep0_pt"    , std::string("pt of  Z1 lep0 of ZZ_{1}#rightarrow ")+decay            ,  300, 0,  300,genVBParticles->at(1).daughterPtr(0)->pt(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z1lep1_pt"    , std::string("pt of  Z1 lep1 of ZZ_{1}#rightarrow ")+decay            ,  300, 0,  300,genVBParticles->at(1).daughterPtr(1)->pt(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Met"          , std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay            ,  200, 0,  300, met->pt(),theWeight);

  
    if((decay != "4l") && ((region_ == phys::MC && topology.test(3)) || (region_ == phys::MC_HZZ && topology.test(1)))) eventsFull.push_back(eventStr);

    theHistograms->fill(std::string("ZZTo")+decay+"_Z0Mass_Sig"   , std::string("Z0  mass  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,mZ1_gen,theWeight);
    theHistograms->fill(std::string("ZZTo")+decay+"_Z1Mass_Sig"   , std::string("Z1  mass  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,mZ2_gen,theWeight);
    theHistograms->fill(std::string("ZZTo")+decay+"_Z0pt_Sig"     , std::string("Z0  p_{T}  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,genVBParticles->at(0).pt(),theWeight);
    theHistograms->fill(std::string("ZZTo")+decay+"_Z1pt_Sig"     , std::string("Z1  p_{T}  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,genVBParticles->at(1).pt(),theWeight);
    theHistograms->fill(std::string("ZZTo")+decay+"_Mass_Sig"     , std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay            ,  200, 50,  1000,  m4L_gen, theWeight);
    theHistograms->fill(std::string("ZZTo")+decay+"_Z0lep0_pt_Sig", std::string("pt of  Z0 lep0 of ZZ_{1}#rightarrow ")+decay            ,  300, 0,  300,genVBParticles->at(0).daughterPtr(0)->pt(),theWeight);
    theHistograms->fill(std::string("ZZTo")+decay+"_Z0lep1_pt_Sig", std::string("pt of  Z0 lep1 of ZZ_{1}#rightarrow ")+decay            ,  300, 0,  300,genVBParticles->at(0).daughterPtr(1)->pt(),theWeight);
    theHistograms->fill(std::string("ZZTo")+decay+"_Z1lep0_pt_Sig", std::string("pt of  Z1 lep0 of ZZ_{1}#rightarrow ")+decay            ,  300, 0,  300,genVBParticles->at(1).daughterPtr(0)->pt(),theWeight);
    theHistograms->fill(std::string("ZZTo")+decay+"_Z1lep1_pt_Sig", std::string("pt of  Z1 lep1 of ZZ_{1}#rightarrow ")+decay            ,  300, 0,  300,genVBParticles->at(1).daughterPtr(1)->pt(),theWeight);
    theHistograms->fill(std::string("ZZTo")+decay+"_Met_Sig"      , std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay            ,  200, 0,  300, met->pt(),theWeight);
  

 
  Float_t   nJets = 0;
  nJets = genJets->size();

  Float_t  ptzz_gen =  (genVBParticles->at(0).p4()+genVBParticles->at(1).p4()).Pt();

  if(nJets>3) nJets=3;

  theHistograms->fill(std::string("ZZTo")+decay+"_nJets"      , "Number of jets (|#eta|<4.7 and p_T > 30 GeV)"            , 4, 0, 4, nJets, theWeight);   

  if(genJets->size()>0){
    theHistograms->fill(std::string("ZZTo")+decay+"_ptJet"      , "#Delta Y(j,j) between the two most energetyc central jets"            ,  100, 20, 400,genJets->at(0).pt(), theWeight);     
    theHistograms->fill(std::string("ZZTo")+decay+"_etaJet"      , "#Delta Y(j,j) between the two most energetyc central jets"            ,  100, -5, 5,genJets->at(0).eta(), theWeight);     
  }
  
  if(genJets->size() >= 2) {
    
    Float_t  mjj =  (genJets->at(0).p4() + genJets->at(1).p4()).M();
    eventStr+=":";
    eventStr+=strtool::sRound(genJets->at(1).pt());
    eventStr+=":"+strtool::sRound(mjj); 
    
    
    if(genJets->at(0).pt() < genJets->at(1).pt()) std::cout<<"jet1 pt "<<genJets->at(0).pt() <<" jet2 pt "<<genJets->at(1).pt()<<std::endl;  
    
    Float_t zZ1 = genVBParticles->at(0).eta()-(genJets->at(0).eta() + genJets->at(1).eta())/2;
    Float_t zZ2 = genVBParticles->at(1).eta()-(genJets->at(0).eta() + genJets->at(1).eta())/2;
    Float_t zZZ = (genVBParticles->at(0).p4()+genVBParticles->at(1).p4()).Eta()-(genJets->at(0).eta() + genJets->at(1).eta())/2;
    
    
    
    //    Float_t SumVecPt = (genVBParticles->at(0).p4()+genVBParticles->at(1).p4()+genJets->at(0).p4()+genJets->at(1).p4()).Pt();
		   //Float_t SumPt = (genVBParticles->at(0).pt()+genVBParticles->at(1).pt()+genJets->at(0).pt()+genJets->at(1).pt());
    
    Float_t PtRatio = ((genVBParticles->at(0).p4()+genVBParticles->at(1).p4()+genJets->at(0).p4()+genJets->at(1).p4()).Pt())/(genVBParticles->at(0).pt()+genVBParticles->at(1).pt()+genJets->at(0).pt()+genJets->at(1).pt());
    Float_t PtJRatio = ((genJets->at(0).p4()+genJets->at(1).p4()).Pt())/(genJets->at(0).pt()+genJets->at(1).pt());

    Float_t Dphi = physmath::deltaPhi(genJets->at(0).phi(),genJets->at(1).phi());    


    theHistograms->fill(std::string("ZZTo")+decay+"_PtZZ"     , "Transverse momentum of ZZ",  100,0,300,ptzz_gen, theWeight); 
    theHistograms->fill(std::string("ZZTo")+decay+"_Dphi"     , "Delta phi of the two leading jets",  100, 0, 3.20, Dphi, theWeight); 
    theHistograms->fill(std::string("ZZTo")+decay+"_DphiVsPtZZ"     , "Delta phi of the two leading jets",  100, 0, 3.20, 100,0,300,Dphi, ptzz_gen, theWeight); 

    theHistograms->fill(std::string("ZZTo")+decay+"_Mjj"      , "Invariant mass of the leading jets",  100, 0, 900, mjj, theWeight); 
    theHistograms->fill(std::string("ZZTo")+decay+"_ptRatio"  , "The ratio of transverse momentum of the vector sum of Z1, Z2, tj1, tj2 to the sum of pTs"        ,  50, 0, 0.8, PtRatio, theWeight); 
    theHistograms->fill(std::string("ZZTo")+decay+"_ptJRatio" , "The ratio of transverse momentum of the vector sum of tj1, tj2 to the sum of pTs"        ,  50, 0, 1., PtJRatio, theWeight); 
    
    theHistograms->fill(std::string("ZZTo")+decay+"_Z1z"      , "Zeppenfeld Variable fo Z1 wrt the two leading jets"            ,  50, -6, 6, zZ1, theWeight); 
    theHistograms->fill(std::string("ZZTo")+decay+"_Z2z"      , "Zeppenfeld Variable fo Z2 wrt the two leading jets"            ,  49, -6, 6, zZ2, theWeight); 
    theHistograms->fill(std::string("ZZTo")+decay+"_ZZz"      , "Zeppenfeld Variable for ZZ wrt the two leading jets"            ,  49, -6, 6, zZZ, theWeight);     
    
    if((region_ == phys::MC && topology.test(2)) || (region_ == phys::MC_HZZ && topology.test(0))){


      theHistograms->fill(std::string("ZZTo")+decay+"_PtZZ_Sig"     , "Transverse momentum of ZZ",  100,0,300,ptzz_gen, theWeight); 
      theHistograms->fill(std::string("ZZTo")+decay+"_Dphi_Sig"     , "Delta phi of the two leading jets",  100, 0, 7, Dphi, theWeight); 
      theHistograms->fill(std::string("ZZTo")+decay+"_DphiVsPtZZ_Sig"     , "Delta phi of the two leading jets",  100, 0, 7, 100,0,300,Dphi, ptzz_gen, theWeight); 
                  
      theHistograms->fill(std::string("ZZTo")+decay+"_Mjj_Sig"      , "Invariant mass of the leading jets",  100, 0, 900, mjj, theWeight); 
      theHistograms->fill(std::string("ZZTo")+decay+"_ptRatio_Sig"      , "The ratio of transverse momentum of the vector sum of Z1, Z2, tj1, tj2 to the sum of pTs"        ,  50, 0, 0.8, PtRatio, theWeight); 
      theHistograms->fill(std::string("ZZTo")+decay+"_ptJRatio_Sig"      , "The ratio of transverse momentum of the vector sum of tj1, tj2 to the sum of pTs"        ,  50, 0, 1., PtJRatio, theWeight); 
      
      theHistograms->fill(std::string("ZZTo")+decay+"_Z1z_Sig"      , "Zeppenfeld Variable for Z1 wrt the two leading jets"            ,  49, -6, 6, zZ1, theWeight); 
      theHistograms->fill(std::string("ZZTo")+decay+"_Z2z_Sig"      , "Zeppenfeld Variable for Z2 wrt the two leading jets"            ,  49, -6, 6, zZ2, theWeight); 
      theHistograms->fill(std::string("ZZTo")+decay+"_ZZz_Sig"      , "Zeppenfeld Variable for ZZ wrt the two leading jets"            ,  49, -6, 6, zZZ, theWeight); 
    }
    
    
    theHistograms->fill(std::string("ZZTo")+decay+"_deltaEtaJJ"      , "#Delta #eta(j,j) between the two most energetic jets"            ,  10, 0, 8, fabs(genJets->at(0).eta() - genJets->at(1).eta()), theWeight); 
    
  }
  else  eventStr+=":0.00:0.00";

  if(nJets >= 3){
    //    Float_t z =  (genJets->at(2).eta()-(genJets->at(0).eta() + genJets->at(1).eta()))/fabs(genJets->at(0).eta() - genJets ->at(1).eta());
    Float_t z =  genJets->at(2).eta()-(genJets->at(0).eta() + genJets->at(1).eta())/2;
    //      std::cout<<<z1<<std::endl;
    theHistograms->fill(std::string("ZZTo")+decay+"_z"      , "z between the two most energetic jets"            ,  21, -10, 10, z, theWeight); 
  }

  theHistograms->fill(std::string("ZZTo")+decay+"_nExtraMuons"    , "Number of extra muons in the event"    , 10, 0, 10, muons->size(), theWeight); 
  theHistograms->fill(std::string("ZZTo")+decay+"_nExtraElectrons", "Number of extra electrons in the event", 10, 0, 10, electrons->size(), theWeight); 
  theHistograms->fill(std::string("ZZTo")+decay+"_nExtraLeptons"  , "Number of extra leptons in the event"  , 10, 0, 10, muons->size()+electrons->size(), theWeight);  




  if((region_ == phys::MC && topology.test(2)) || (region_ == phys::MC_HZZ && topology.test(0))){

 if(nJets>3) nJets=3;

  theHistograms->fill(std::string("ZZTo")+decay+"_nJets_Sig"      , "Number of jets (|#eta|<4.7 and p_T > 30 GeV)"            , 4, 0, 4, nJets, theWeight);   

  if(genJets->size()>0){
    theHistograms->fill(std::string("ZZTo")+decay+"_ptJet_Sig"      , "#Delta Y(j,j) between the two most energetyc central jets"            ,  100, 20, 400,genJets->at(0).pt(), theWeight);     
    theHistograms->fill(std::string("ZZTo")+decay+"_etaJet_Sig"      , "#Delta Y(j,j) between the two most energetyc central jets"            ,  100, -5, 5,genJets->at(0).eta(), theWeight);     
  }
  
  if(genJets->size() >= 2) {
  if(genJets->at(0).pt() < genJets->at(1).pt()) std::cout<<"jet1 pt "<<genJets->at(0).pt() <<" jet2 pt "<<genJets->at(1).pt()<<std::endl;  

    theHistograms->fill(std::string("ZZTo")+decay+"_deltaEtaJJ_Sig"      , "#Delta #eta(j,j) between the two most energetic jets"            ,  10, 0, 8, fabs(genJets->at(0).eta() - genJets->at(1).eta()), theWeight); 
    
}


  if(nJets >= 3){
    //    Float_t z =  (genJets->at(2).eta()-(genJets->at(0).eta() + genJets->at(1).eta()))/fabs(genJets->at(0).eta() - genJets ->at(1).eta());
    Float_t z =  genJets->at(2).eta()-(genJets->at(0).eta() + genJets->at(1).eta())/2.;
    //      std::cout<<<z1<<std::endl;
    theHistograms->fill(std::string("ZZTo")+decay+"_z_Sig"      , "z between the two most energetic jets"            ,  21, -10, 10, z, theWeight); 
  }

  theHistograms->fill(std::string("ZZTo")+decay+"_nExtraMuons_Sig"    , "Number of extra muons in the event"    , 10, 0, 10, muons->size(), theWeight); 
  theHistograms->fill(std::string("ZZTo")+decay+"_nExtraElectrons_Sig", "Number of extra electrons in the event", 10, 0, 10, electrons->size(), theWeight); 
  theHistograms->fill(std::string("ZZTo")+decay+"_nExtraLeptons_Sig"  , "Number of extra leptons in the event"  , 10, 0, 10, muons->size()+electrons->size(), theWeight);  
  }


}

void ZZjGenAnalyzer::analyze(){


  if((region_ == phys::MC && topology.test(2)) || (region_ == phys::MC_HZZ && topology.test(0))){  
  std::string eventstr=std::to_string(run)+":"+std::to_string(lumiBlock)+":"+std::to_string(event);
 
  Float_t w_kf = 1.;
  
  //  if((theMCInfo.fileName()=="ggZZ2e2mu") || (theMCInfo.fileName()=="ggZZ4e") || (theMCInfo.fileName()=="ggZZ4mu") || (theMCInfo.fileName()=="ggH125") )   w_kf = theMCInfo.kF_ggZZ() ; 
  //  else if((theMCInfo.fileName()=="ZZTo4l") || (theMCInfo.fileName()=="ZZTo4lamcatnlo")) w_kf = theMCInfo.kF_qqZZM() * theMCInfo.kF_EWKqqZZ() ; 
  
 if((theMCInfo.fileName()=="ggZZ2e2mu") || (theMCInfo.fileName()=="ggZZ4e") || (theMCInfo.fileName()=="ggZZ4mu") || (theMCInfo.fileName()=="ggTo2e2mu_Contin_MCFM701") || (theMCInfo.fileName()=="ggTo4e_Contin_MCFM701") || (theMCInfo.fileName()=="ggTo4mu_Contin_MCFM701"))   w_kf = 1.7 ; 
 else if(theMCInfo.fileName()=="ZZTo4l") w_kf = 1.1; 

 theWeight*=w_kf;

  int Ele  = 0;
  int Muon = 0; 
  int lep = 0;
  
  foreach(const phys::Particle &gen, *genParticles){
    if((!gen.genStatusFlags().test(GenStatusBit::isPrompt)) || (!gen.genStatusFlags().test(GenStatusBit::fromHardProcess)))  continue;
    if(abs(gen.id())==13){	
      Muon += 1; 
      lep +=1;
    }
    else if(abs(gen.id())==11){
      Ele += 1;
      lep+=1; 
    }
  }

  bool doAnalyse =kTRUE;  
  std::string decay="None";
  if(Ele==Muon && Ele!=0)  {decay = "2e2m";} 
  else if(Ele<Muon) {decay = "4m";}    
  else if(Ele>Muon) {decay = "4e";}   
  else if(Ele==Muon && Ele == 0) std::cout<<"NO FINALSTATE"<<std::endl;
  if(lep<4){doAnalyse = kFALSE;}
  
  if(doAnalyse){
    ZZplots("4l");
    ZZplots(decay);   
  } 
}

}
void ZZjGenAnalyzer::end( TFile &) {
  
  //   std::cout<<"eeee:"<<std::endl;
  //   for (std::vector<std::string>::iterator it = events4e.begin() ; it != events4e.end(); ++it) std::cout<<"    "<<*it<<std::endl;
  //   std::cout<<"eemm:"<<std::endl;
  //   for (std::vector<std::string>::iterator it = events2e2mu.begin() ; it != events2e2mu.end(); ++it) std::cout<<"    "<<*it<<std::endl;
  //   std::cout<<"mmmm:"<<std::endl;
  //   for (std::vector<std::string>::iterator it = events4mu.begin() ; it != events4mu.end(); ++it) std::cout<<"    "<<*it<<std::endl;

  std::sort (eventsFull.begin(), eventsFull.end(), strtool::sortEvents);
  for (std::vector<std::string>::iterator it = eventsFull.begin() ; it != eventsFull.end(); ++it) std::cout<<*it<<std::endl;

  //std::cout<<"Final \neeee "<<events4e.size()<<std::endl<<"eemm "<<events2e2mu.size()<<std::endl<<"mmmm "<<events4mu.size()<<"\nsum "<<events4e.size()+events2e2mu.size()+events4mu.size()<<" total "<<events4l.size()<<std::endl;
    
}

  
