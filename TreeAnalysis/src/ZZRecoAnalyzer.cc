#include "VVXAnalysis/TreeAnalysis/interface/ZZRecoAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "TRandom.h"
#include "TTree.h"
#include <boost/assign/std/vector.hpp> 
#include <boost/assert.hpp> 
using namespace std;
using namespace boost::assign; // bring 'operator+=()' into scope

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using std::cout;
using std::endl;

using namespace phys;


void ZZRecoAnalyzer::analyze(){

  UpJER_jets->clear();
  DownJER_jets->clear();
  UpJER_centraljets->clear();
  DownJER_centraljets->clear();
  UpJES_jets->clear();
  DownJES_jets->clear();
  UpJES_centraljets->clear();
  DownJES_centraljets->clear();
  UpJESData_jets->clear();
  DownJESData_jets->clear();
  UpJESData_centraljets->clear();
  DownJESData_centraljets->clear();
  
  e++;
  
  std::string decay  = "";
  std::string sample = "";
  std::string region = "";

  Int_t id = ZZ->id();

  if      (id == 52) {decay = "4m";}
  else if (id == 48) {decay = "2e2m";}
  else if (id == 44) {decay = "4e";}

  if(e < nentries/2){sample = "0";}
  else {sample = "1";}

  m4L = ZZ->mass();  
  //  if(m4L > 800) m4L = 799;

  drzz =physmath::deltaR(ZZ->first(),ZZ->second());
  //  if(ptzz>=300) ptzz = 299;

  ptzz = ZZ->pt();
  //  if(ptzz>=300) ptzz = 299;

  Float_t w_kf = 1.;

  //  if((theMCInfo.fileName()=="ggZZ2e2mu") || (theMCInfo.fileName()=="ggZZ4e") || (theMCInfo.fileName()=="ggZZ4mu"))   w_kf = kFactor_ggZZ ; 
  //  else if((theMCInfo.fileName()=="ZZTo4l") || (theMCInfo.fileName()=="ZZTo4lamcatnlo")) w_kf = kFactor_qqZZM * kFactor_EWKqqZZ ; 

  if((theMCInfo.fileName()=="ggZZ2e2mu") || (theMCInfo.fileName()=="ggZZ4e") || (theMCInfo.fileName()=="ggZZ4mu"))  w_kf = 1.7 ;
  else if(theMCInfo.fileName()=="ZZTo4l") w_kf = 1.1;

  theWeight*=w_kf; 
  
  Float_t scaleFacErrSq    = ZZ->efficiencySFUnc(); 
  Float_t scaleFacMuErrSq  = ZZ->muEffSFUnc(); 
  Float_t scaleFacEleErrSq = ZZ->eleEffSFUnc(); 

  //Data 
  if(!theMCInfo.isMC()){    
     
    foreach(const phys::Jet &dataJet, *pjets){
      //  if(!dataJet.fullPuId()) continue; //to activate Pu id
      double dataJetPt = 0;
      dataJetPt = dataJet.pt();
      
      //  JES correction: Up and down velues used to assess systematic uncertainty on jet energy resolution
      double newJetPtJESData_up   =0;  
      double newJetPtJESData_down =0;

      newJetPtJESData_up   = dataJetPt*(1 + dataJet.jecUncertainty());
      newJetPtJESData_down = dataJetPt*(1 - dataJet.jecUncertainty());

      if(newJetPtJESData_up > 30)     UpJESData_jets->push_back(dataJet);       
      if(newJetPtJESData_down > 30)   DownJESData_jets->push_back(dataJet);

      if(newJetPtJESData_up > 30 && fabs(dataJet.eta())<2.4)   UpJESData_centraljets->push_back(dataJet); 
      if(newJetPtJESData_down > 30 && fabs(dataJet.eta())<2.4) DownJESData_centraljets->push_back(dataJet);      
   
      stable_sort(UpJESData_jets->begin(), UpJESData_jets->end(), PtComparator());
      stable_sort(DownJESData_jets->begin(), DownJESData_jets->end(), PtComparator());

    }

    FillHistosJets(decay,theWeight,UpJESData_jets,"JESDataUp_01");
    FillHistosJets(decay,theWeight,DownJESData_jets,"JESDataDn_01");
    FillHistosJets(decay,theWeight,UpJESData_centraljets,"Central_JESDataUp_01");
    FillHistosJets(decay,theWeight,DownJESData_centraljets,"Central_JESDataDn_01");
 }
  
  else{

    foreach(const phys::Jet &jet, *pjets){
      // if(!jet.fullPuId()) continue;

      //      if(jet.ptJerDn() < jet.pt()) cout<<"jet.ptJerDn() "<<jet.ptJerDn() <<"jet.ptJerUp() "<<jet.ptJerUp() <<" jet.pt() "<<jet.pt()<<endl;
      //       phys::Jet *newjet =  new  phys::Jet();
      // //      phys

      //       *newjet = jet;
      // 	    TLorentzVector newvec = jet->p4();
      // 	    newvec.SetPt
      // 	    newjet->
      // cout<<newjet->pt()<<endl;

      // jet
      if(jet.ptJerUp() > 30.) 	UpJER_jets->push_back(jet); 
      if(jet.ptJerDn() > 30.)	DownJER_jets->push_back(jet);
      
      if(jet.ptJerDn()  > 30. && fabs(jet.eta())<2.4)   UpJER_centraljets->push_back(jet); 
      if(jet.ptJerDn()  > 30. && fabs(jet.eta())<2.4) DownJER_centraljets->push_back(jet);

     

      //JES correction: Up and down velues used to assess systematic uncertainty on jet energy resolution
      double newJetPtJES_up   = 0;  
      double newJetPtJES_down = 0;
           
      newJetPtJES_up   = jet.pt()*(1+jet.jecUncertainty());
      newJetPtJES_down = jet.pt()*(1-jet.jecUncertainty());
      
      if(newJetPtJES_up > 30)  	  UpJES_jets->push_back(jet);
      if(newJetPtJES_down > 30)   DownJES_jets->push_back(jet);
      if(newJetPtJES_up > 30 && fabs(jet.eta())<2.4 ) UpJES_centraljets->push_back(jet); 
      if(newJetPtJES_down > 30 &&  fabs(jet.eta())<2.4) DownJES_centraljets->push_back(jet);
      }


    FillHistosBase(decay,theWeight,sample);
    FillHistosJets(decay,theWeight,jets,sample);
    FillHistosJets(decay,theWeight,centralJets,"Central_"+sample);   

    FillHistosJets(decay,theWeight,UpJER_jets,"JERUp_01");
    FillHistosJets(decay,theWeight,DownJER_jets,"JERDn_01");

    FillHistosJets(decay,theWeight,UpJER_centraljets,"Central_JERUp_01");
    FillHistosJets(decay,theWeight,DownJER_centraljets,"Central_JERDn_01");

    FillHistosJets(decay,theWeight,UpJES_jets,"JESUp_01");
    FillHistosJets(decay,theWeight,DownJES_jets,"JESDn_01");

    FillHistosJets(decay,theWeight,UpJES_centraljets,"Central_JESUp_01");
    FillHistosJets(decay,theWeight,DownJES_centraljets,"Central_JESDn_01");

    FillHistosBase(decay,theWeight*(1-scaleFacErrSq),"SFSqDn_01");
    FillHistosBase(decay,theWeight*(1+scaleFacErrSq),"SFSqUp_01");

    FillHistosJets(decay,theWeight*(1-scaleFacErrSq),jets,"SFSqDn_01");
    FillHistosJets(decay,theWeight*(1+scaleFacErrSq),jets,"SFSqUp_01");

    FillHistosJets(decay,theWeight*(1-scaleFacErrSq),centralJets,"Central_SFSqDn_01");
    FillHistosJets(decay,theWeight*(1+scaleFacErrSq),centralJets,"Central_SFSqUp_01");

    FillHistosBase(decay,theWeight*(1-scaleFacMuErrSq),"MuSFSqDn_01");
    FillHistosBase(decay,theWeight*(1+scaleFacMuErrSq),"MuSFSqUp_01");

    FillHistosJets(decay,theWeight*(1-scaleFacMuErrSq),jets,"MuSFSqDn_01");
    FillHistosJets(decay,theWeight*(1+scaleFacMuErrSq),jets,"MuSFSqUp_01");

    FillHistosJets(decay,theWeight*(1-scaleFacMuErrSq),centralJets,"Central_MuSFSqDn_01");
    FillHistosJets(decay,theWeight*(1+scaleFacMuErrSq),centralJets,"Central_MuSFSqUp_01");

    FillHistosBase(decay,theWeight*(1-scaleFacEleErrSq),"EleSFSqDn_01");
    FillHistosBase(decay,theWeight*(1+scaleFacEleErrSq),"EleSFSqUp_01");

    FillHistosJets(decay,theWeight*(1-scaleFacEleErrSq),jets,"EleSFSqDn_01");
    FillHistosJets(decay,theWeight*(1+scaleFacEleErrSq),jets,"EleSFSqUp_01");

    FillHistosJets(decay,theWeight*(1-scaleFacEleErrSq),centralJets,"Central_EleSFSqDn_01");
    FillHistosJets(decay,theWeight*(1+scaleFacEleErrSq),centralJets,"Central_EleSFSqUp_01");
    
    FillHistosBase(decay,theWeight*theMCInfo.puWeightUncDn(),"PuDn_01");
    FillHistosBase(decay,theWeight*theMCInfo.puWeightUncUp(),"PuUp_01");
    
    FillHistosJets(decay,theWeight*theMCInfo.puWeightUncDn(),jets,"PuDn_01");
    FillHistosJets(decay,theWeight*theMCInfo.puWeightUncUp(),jets,"PuUp_01");
    
    FillHistosJets(decay,theWeight*theMCInfo.puWeightUncDn(),centralJets,"Central_PuDn_01");
    FillHistosJets(decay,theWeight*theMCInfo.puWeightUncUp(),centralJets,"Central_PuUp_01");

    //pdf alpha s
    FillHistosBase(decay,theWeight*LHEweight_PDFVariation_Dn,"PDFDn_01");
    FillHistosBase(decay,theWeight*LHEweight_PDFVariation_Up,"PDFUp_01");

    FillHistosJets(decay,theWeight*LHEweight_PDFVariation_Dn,jets,"PDFDn_01");
    FillHistosJets(decay,theWeight*LHEweight_PDFVariation_Up,jets,"PDFUp_01");

    FillHistosJets(decay,theWeight*LHEweight_PDFVariation_Dn,centralJets,"Central_PDFDn_01");
    FillHistosJets(decay,theWeight*LHEweight_PDFVariation_Up,centralJets,"Central_PDFUp_01");

    FillHistosBase(decay,theWeight*LHEweight_AsMZ_Dn,"AsDn_01");
    FillHistosBase(decay,theWeight*LHEweight_AsMZ_Up,"AsUp_01");

    FillHistosJets(decay,theWeight*LHEweight_AsMZ_Dn,jets,"AsDn_01");
    FillHistosJets(decay,theWeight*LHEweight_AsMZ_Up,jets,"AsUp_01");

    FillHistosJets(decay,theWeight*LHEweight_AsMZ_Dn,centralJets,"Central_AsDn_01");
    FillHistosJets(decay,theWeight*LHEweight_AsMZ_Up,centralJets,"Central_AsUp_01");


    //Is Signal
   if((region_ == phys::SR && topology.test(2)) || (region_ == phys::SR_HZZ && topology.test(0) ) ){       

     m4L_gen = sqrt((genVBParticles->at(0).p4()+genVBParticles->at(1).p4())*(genVBParticles->at(0).p4()+genVBParticles->at(1).p4()));
     //     if (m4L_gen>=800) m4L_gen = 799;
     
     drzz_gen =physmath::deltaR(genVBParticles->at(0),genVBParticles->at(1));
     
     ptzz_gen =  (genVBParticles->at(0).p4()+genVBParticles->at(1).p4()).Pt();
     //  if (ptzz_gen>300) ptzz_gen=299;
     dphizz_gen = fabs(physmath::deltaPhi(genVBParticles->at(0).phi(),genVBParticles->at(1).phi())); 
     
     
    FillMatrixHistosBase(decay,theWeight,sample);
    FillMatrixHistosJets(decay,theWeight, jets, genJets,sample);
    FillMatrixHistosJets(decay,theWeight, centralJets, centralGenJets,"Central_"+sample);

    FillMatrixHistosBase(decay,theWeight,"01");
    FillMatrixHistosJets(decay,theWeight,jets,genJets,"01");
    FillMatrixHistosJets(decay,theWeight,centralJets,centralGenJets,"Central_01");

    FillMatrixHistosJets(decay,theWeight,UpJER_jets,genJets,"JERUp_01");
    FillMatrixHistosJets(decay,theWeight,DownJER_jets,genJets,"JERDn_01");

    FillMatrixHistosJets(decay,theWeight,UpJER_centraljets,centralGenJets,"Central_JERUp_01");
    FillMatrixHistosJets(decay,theWeight,DownJER_centraljets,centralGenJets,"Central_JERDn_01");

    FillMatrixHistosJets(decay,theWeight,UpJES_jets,genJets,"JESUp_01");
    FillMatrixHistosJets(decay,theWeight,DownJES_jets,genJets,"JESDn_01");

    FillMatrixHistosJets(decay,theWeight,UpJES_centraljets,centralGenJets,"Central_JESUp_01");
    FillMatrixHistosJets(decay,theWeight,DownJES_centraljets,centralGenJets,"Central_JESDn_01");

    FillMatrixHistosBase(decay,theWeight*(1-scaleFacErrSq),"SFSqDn_01");
    FillMatrixHistosBase(decay,theWeight*(1+scaleFacErrSq),"SFSqUp_01");

    FillMatrixHistosJets(decay,theWeight*(1-scaleFacErrSq),jets,genJets,"SFSqDn_01");
    FillMatrixHistosJets(decay,theWeight*(1+scaleFacErrSq),jets,genJets,"SFSqUp_01");

    FillMatrixHistosJets(decay,theWeight*(1-scaleFacErrSq),centralJets,centralGenJets,"Central_SFSqDn_01");
    FillMatrixHistosJets(decay,theWeight*(1+scaleFacErrSq),centralJets,centralGenJets,"Central_SFSqUp_01");

    FillMatrixHistosBase(decay,theWeight*(1-scaleFacMuErrSq),"MuSFSqDn_01");
    FillMatrixHistosBase(decay,theWeight*(1+scaleFacMuErrSq),"MuSFSqUp_01");

    FillMatrixHistosJets(decay,theWeight*(1-scaleFacMuErrSq),jets,genJets,"MuSFSqDn_01");
    FillMatrixHistosJets(decay,theWeight*(1+scaleFacMuErrSq),jets,genJets,"MuSFSqUp_01");

    FillMatrixHistosJets(decay,theWeight*(1-scaleFacMuErrSq),centralJets,centralGenJets,"Central_MuSFSqDn_01");
    FillMatrixHistosJets(decay,theWeight*(1+scaleFacMuErrSq),centralJets,centralGenJets,"Central_MuSFSqUp_01");

    FillMatrixHistosBase(decay,theWeight*(1-scaleFacEleErrSq),"EleSFSqDn_01");
    FillMatrixHistosBase(decay,theWeight*(1+scaleFacEleErrSq),"EleSFSqUp_01");

    FillMatrixHistosJets(decay,theWeight*(1-scaleFacEleErrSq),jets,genJets,"EleSFSqDn_01");
    FillMatrixHistosJets(decay,theWeight*(1+scaleFacEleErrSq),jets,genJets,"EleSFSqUp_01");

    FillMatrixHistosJets(decay,theWeight*(1-scaleFacEleErrSq),centralJets,centralGenJets,"Central_EleSFSqDn_01");
    FillMatrixHistosJets(decay,theWeight*(1+scaleFacEleErrSq),centralJets,centralGenJets,"Central_EleSFSqUp_01");
 
    FillMatrixHistosBase(decay,theWeight*theMCInfo.puWeightUncDn(),"PuDn_01");
    FillMatrixHistosBase(decay,theWeight*theMCInfo.puWeightUncUp(),"PuUp_01");

    FillMatrixHistosJets(decay,theWeight*theMCInfo.puWeightUncDn(),jets,genJets,"PuDn_01");
    FillMatrixHistosJets(decay,theWeight*theMCInfo.puWeightUncUp(),jets,genJets,"PuUp_01");

    FillMatrixHistosJets(decay,theWeight*theMCInfo.puWeightUncDn(),centralJets,centralGenJets,"Central_PuDn_01");
    FillMatrixHistosJets(decay,theWeight*theMCInfo.puWeightUncUp(),centralJets,centralGenJets,"Central_PuUp_01");

    //pdf alpha s
    FillMatrixHistosBase(decay,theWeight*LHEweight_PDFVariation_Dn,"PDFDn_01");
    FillMatrixHistosBase(decay,theWeight*LHEweight_PDFVariation_Up,"PDFUp_01");

    FillMatrixHistosJets(decay,theWeight*LHEweight_PDFVariation_Dn,jets,genJets,"PDFDn_01");
    FillMatrixHistosJets(decay,theWeight*LHEweight_PDFVariation_Up,jets,genJets,"PDFUp_01");

    FillMatrixHistosJets(decay,theWeight*LHEweight_PDFVariation_Dn,centralJets,centralGenJets,"Central_PDFDn_01");
    FillMatrixHistosJets(decay,theWeight*LHEweight_PDFVariation_Up,centralJets,centralGenJets,"Central_PDFUp_01");

    FillMatrixHistosBase(decay,theWeight*LHEweight_AsMZ_Dn,"AsDn_01");
    FillMatrixHistosBase(decay,theWeight*LHEweight_AsMZ_Up,"AsUp_01");

    FillMatrixHistosJets(decay,theWeight*LHEweight_AsMZ_Dn,jets,genJets,"AsDn_01");
    FillMatrixHistosJets(decay,theWeight*LHEweight_AsMZ_Up,jets,genJets,"AsUp_01");

    FillMatrixHistosJets(decay,theWeight*LHEweight_AsMZ_Dn,centralJets,centralGenJets,"Central_AsDn_01");
    FillMatrixHistosJets(decay,theWeight*LHEweight_AsMZ_Up,centralJets,centralGenJets,"Central_AsUp_01");

    }


    if((region_ == phys::SR && topology.test(3)) || (region_ == phys::SR_HZZ && topology.test(1))){

       inFiducialRegion ++;
       
       FillMatrixHistosBase(decay,theWeight,"01_fr");     
       FillMatrixHistosBase(decay,theWeight,sample+"_fr");

       FillMatrixHistosJets(decay,theWeight,jets,genJets,"01_fr");
       FillMatrixHistosJets(decay,theWeight,jets,genJets,sample+"_fr");

       FillMatrixHistosJets(decay,theWeight,centralJets,centralGenJets,"Central_01_fr");
       FillMatrixHistosJets(decay,theWeight,centralJets,centralGenJets,"Central_"+sample+"_fr");

       FillMatrixHistosJets(decay,theWeight,UpJER_jets,genJets,"JERUp_01_fr");
       FillMatrixHistosJets(decay,theWeight,DownJER_jets,genJets,"JERDn_01_fr");
       
       FillMatrixHistosJets(decay,theWeight,UpJER_centraljets,centralGenJets,"Central_JERUp_01_fr");
       FillMatrixHistosJets(decay,theWeight,DownJER_centraljets,centralGenJets,"Central_JERDn_01_fr");
       
       FillMatrixHistosJets(decay,theWeight,UpJES_jets,genJets,"JESUp_01_fr");
       FillMatrixHistosJets(decay,theWeight,DownJES_jets,genJets,"JESDn_01_fr");
       
       FillMatrixHistosJets(decay,theWeight,UpJES_centraljets,centralGenJets,"Central_JESUp_01_fr");
       FillMatrixHistosJets(decay,theWeight,DownJES_centraljets,centralGenJets,"Central_JESDn_01_fr");
       
       FillMatrixHistosBase(decay,theWeight*(1-scaleFacErrSq),"SFSqDn_01_fr");
       FillMatrixHistosBase(decay,theWeight*(1+scaleFacErrSq),"SFSqUp_01_fr");
       
       FillMatrixHistosJets(decay,theWeight*(1-scaleFacErrSq),jets,genJets,"SFSqDn_01_fr");
       FillMatrixHistosJets(decay,theWeight*(1+scaleFacErrSq),jets,genJets,"SFSqUp_01_fr");
       
       FillMatrixHistosJets(decay,theWeight*(1-scaleFacErrSq),centralJets,centralGenJets,"Central_SFSqDn_01_fr");
       FillMatrixHistosJets(decay,theWeight*(1+scaleFacErrSq),centralJets,centralGenJets,"Central_SFSqUp_01_fr");

       FillMatrixHistosBase(decay,theWeight*theMCInfo.puWeightUncDn(),"PuDn_01_fr");
       FillMatrixHistosBase(decay,theWeight*theMCInfo.puWeightUncUp(),"PuUp_01_fr");
       
       FillMatrixHistosJets(decay,theWeight*theMCInfo.puWeightUncDn(),jets,genJets,"PuDn_01_fr");
       FillMatrixHistosJets(decay,theWeight*theMCInfo.puWeightUncUp(),jets,genJets,"PuUp_01_fr");
       
       FillMatrixHistosJets(decay,theWeight*theMCInfo.puWeightUncDn(),centralJets,centralGenJets,"Central_PuDn_01_fr");
       FillMatrixHistosJets(decay,theWeight*theMCInfo.puWeightUncUp(),centralJets,centralGenJets,"Central_PuUp_01_fr");       

       FillMatrixHistosBase(decay,theWeight*(1-scaleFacMuErrSq),"MuSFSqDn_01_fr");
       FillMatrixHistosBase(decay,theWeight*(1+scaleFacMuErrSq),"MuSFSqUp_01_fr");
       
       FillMatrixHistosJets(decay,theWeight*(1-scaleFacMuErrSq),jets,genJets,"MuSFSqDn_01_fr");
       FillMatrixHistosJets(decay,theWeight*(1+scaleFacMuErrSq),jets,genJets,"MuSFSqUp_01_fr");
       
       FillMatrixHistosJets(decay,theWeight*(1-scaleFacMuErrSq),centralJets,centralGenJets,"Central_MuSFSqDn_01_fr");
       FillMatrixHistosJets(decay,theWeight*(1+scaleFacMuErrSq),centralJets,centralGenJets,"Central_MuSFSqUp_01_fr");
       
       FillMatrixHistosBase(decay,theWeight*(1-scaleFacEleErrSq),"EleSFSqDn_01_fr");
       FillMatrixHistosBase(decay,theWeight*(1+scaleFacEleErrSq),"EleSFSqUp_01_fr");
       
       FillMatrixHistosJets(decay,theWeight*(1-scaleFacEleErrSq),jets,genJets,"EleSFSqDn_01_fr");
       FillMatrixHistosJets(decay,theWeight*(1+scaleFacEleErrSq),jets,genJets,"EleSFSqUp_01_fr");
       
       FillMatrixHistosJets(decay,theWeight*(1-scaleFacEleErrSq),centralJets,centralGenJets,"Central_EleSFSqDn_01_fr");
       FillMatrixHistosJets(decay,theWeight*(1+scaleFacEleErrSq),centralJets,centralGenJets,"Central_EleSFSqUp_01_fr");

       //pdf alpha s       
       
       FillMatrixHistosBase(decay,theWeight*LHEweight_PDFVariation_Dn,"PDFDn_01_fr");
       FillMatrixHistosBase(decay,theWeight*LHEweight_PDFVariation_Up,"PDFUp_01_fr");
       
       FillMatrixHistosJets(decay,theWeight*LHEweight_PDFVariation_Dn,jets,genJets,"PDFDn_01_fr");
       FillMatrixHistosJets(decay,theWeight*LHEweight_PDFVariation_Up,jets,genJets,"PDFUp_01_fr");
       
       FillMatrixHistosJets(decay,theWeight*LHEweight_PDFVariation_Dn,centralJets,centralGenJets,"Central_PDFDn_01_fr");
       FillMatrixHistosJets(decay,theWeight*LHEweight_PDFVariation_Up,centralJets,centralGenJets,"Central_PDFUp_01_fr");
       
       FillMatrixHistosBase(decay,theWeight*LHEweight_AsMZ_Dn,"AsDn_01_fr");
       FillMatrixHistosBase(decay,theWeight*LHEweight_AsMZ_Up,"AsUp_01_fr");
       
       FillMatrixHistosJets(decay,theWeight*LHEweight_AsMZ_Dn,jets,genJets,"AsDn_01_fr");
       FillMatrixHistosJets(decay,theWeight*LHEweight_AsMZ_Up,jets,genJets,"AsUp_01_fr");
       
       FillMatrixHistosJets(decay,theWeight*LHEweight_AsMZ_Dn,centralJets,centralGenJets,"Central_AsDn_01_fr");
       FillMatrixHistosJets(decay,theWeight*LHEweight_AsMZ_Up,centralJets,centralGenJets,"Central_AsUp_01_fr");
       
       //1D

       FillHistosBase(decay,theWeight,"01_fr");
       FillHistosBase(decay,theWeight*(1-scaleFacMuErrSq),"MuSFSqDn_01_fr");
       FillHistosBase(decay,theWeight*(1+scaleFacMuErrSq),"MuSFSqUp_01_fr");
            
       FillHistosBase(decay,theWeight*(1-scaleFacEleErrSq),"EleSFSqDn_01_fr");
       FillHistosBase(decay,theWeight*(1+scaleFacEleErrSq),"EleSFSqUp_01_fr");
      
       FillHistosBase(decay,theWeight*theMCInfo.puWeightUncDn(),"PuDn_01_fr");
       FillHistosBase(decay,theWeight*theMCInfo.puWeightUncUp(),"PuUp_01_fr");
      
       //pdf alpha s
       
       FillHistosBase(decay,theWeight*LHEweight_PDFVariation_Dn,"PDFDn_01_fr");
       FillHistosBase(decay,theWeight*LHEweight_PDFVariation_Up,"PDFUp_01_fr");
                  
       FillHistosBase(decay,theWeight*LHEweight_AsMZ_Dn,"AsDn_01_fr");
       FillHistosBase(decay,theWeight*LHEweight_AsMZ_Up,"AsUp_01_fr");
            
    } // end fiducial region
    else{
      FillHistosBase(decay,theWeight,"01_nofr");

       FillHistosBase(decay,theWeight*(1-scaleFacMuErrSq),"MuSFSqDn_01_nofr");
       FillHistosBase(decay,theWeight*(1+scaleFacMuErrSq),"MuSFSqUp_01_nofr");
       
       FillHistosBase(decay,theWeight*(1-scaleFacEleErrSq),"EleSFSqDn_01_nofr");
       FillHistosBase(decay,theWeight*(1+scaleFacEleErrSq),"EleSFSqUp_01_nofr");
       
       FillHistosBase(decay,theWeight*theMCInfo.puWeightUncDn(),"PuDn_01_nofr");
       FillHistosBase(decay,theWeight*theMCInfo.puWeightUncUp(),"PuUp_01_nofr");

       //pdf alpha s
       
       FillHistosBase(decay,theWeight*LHEweight_PDFVariation_Dn,"PDFDn_01_nofr");
       FillHistosBase(decay,theWeight*LHEweight_PDFVariation_Up,"PDFUp_01_nofr");
       
       FillHistosBase(decay,theWeight*LHEweight_AsMZ_Dn,"AsDn_01_nofr");
       FillHistosBase(decay,theWeight*LHEweight_AsMZ_Up,"AsUp_01_nofr");
    }
  } //end is MC
  
  
  FillHistosBase(decay,theWeight,"01");
  FillHistosJets(decay,theWeight,centralJets,"Central_01");
  FillHistosJets(decay,theWeight,jets,"01");

  if(region_ == phys::CR3P1F || region_ == phys::CR2P2F) {
    FillHistosBase(decay,ZZ->fakeRateSFVar(),"FRVar");
    FillHistosJets(decay,ZZ->fakeRateSFVar(),jets,"FRVar");
    FillHistosJets(decay,ZZ->fakeRateSFVar(),centralJets,"Central_FRVar");
  }
}


 void ZZRecoAnalyzer::FillHistosBase(std::string decay, float Wh,std::string type ){

  theHistograms.fill("ZZTo"+decay+"_Mass_"+type,"", Xbins, m4L,Wh);
  theHistograms.fill("ZZTo"+decay+"_PtZZ_"+type, Xbins_ptzz, ptzz,Wh);
  theHistograms.fill("ZZTo"+decay+"_DphiZZ_"+type, Xbins_dphizz, dphizz,Wh);
  theHistograms.fill("ZZTo"+decay+"_dRZZ_"+type,"", Xbins_drzz, drzz,Wh);
  
}

void ZZRecoAnalyzer::FillHistosJets(std::string decay,float Wh,std::vector<phys::Jet> *jetsVec,std::string type){
  
  njets =  jetsVec->size(); 
  
  if (njets>3) njets=3;
  theHistograms.fill("ZZTo"+decay+"_nJets_"+type, "", Xbins_nJets,njets, Wh);
  
  if(njets>0){  
  
    stable_sort(jetsVec->begin(), jetsVec->end(), PtComparator());

     if (type.find("JESUp") != std::string::npos)     ptJet1  =  jetsVec->at(0).ptJerUp();
     else if(type.find("JESDn") != std::string::npos) ptJet1  =  jetsVec->at(0).ptJerDn();
     else if(type.find("JERUp") != std::string::npos) ptJet1  =  jetsVec->at(0).pt()*(1+jetsVec->at(0).jecUncertainty());
     else if(type.find("JERDn") != std::string::npos) ptJet1  =  jetsVec->at(0).pt()*(1-jetsVec->at(0).jecUncertainty());
     else{
       ptJet1 = jetsVec->at(0).pt(); 
     }



    //    if (ptJet1>=500) ptJet1 = 499;    

    theHistograms.fill("ZZTo"+decay+"_PtJet1_"+type," ", Xbins_ptJet1, ptJet1, Wh); 
    
    etaJet1 = fabs(jetsVec->at(0).eta());
    //if (etaJet1>=4.7) etaJet1 = 4.6;
    theHistograms.fill("ZZTo"+decay+"_EtaJet1_"+type, "", Xbins_etaJet1, etaJet1, Wh);
  }
   
   if(njets>1){  
     deta = fabs(jetsVec->at(0).eta() - jetsVec->at(1).eta());
     
     if (deta>=4.7) deta = 4.6;
     mjj =  (jetsVec->at(0).p4() + jetsVec->at(1).p4()).M();

     if (type.find("JESUp") != std::string::npos)     ptJet2  =  jetsVec->at(1).ptJerUp();
     else if(type.find("JESDn") != std::string::npos) ptJet2  =  jetsVec->at(1).ptJerDn();
     else if(type.find("JERUp") != std::string::npos) ptJet2  =  jetsVec->at(1).pt()*(1+jetsVec->at(1).jecUncertainty());
     else if(type.find("JERDn") != std::string::npos) ptJet2  =  jetsVec->at(1).pt()*(1-jetsVec->at(1).jecUncertainty());
     else{
       ptJet2 = jetsVec->at(1).pt(); 
     }


     //if (ptJet2>=300) ptJet2 = 299;

     dphi = physmath::deltaPhi(jetsVec->at(0).phi(),jetsVec->at(1).phi());
     //if (dphi>=6) dphi = 5;
     
     theHistograms.fill("ZZTo"+decay+"_PtJet2_"+type, "", Xbins_ptJet2, ptJet2, Wh);
     
     etaJet2 = fabs(jetsVec->at(1).eta());
     //if (etaJet2>=4.7) etaJet2 = 4.6;
     
     theHistograms.fill("ZZTo"+decay+"_EtaJet2_"+type, "", Xbins_etaJet2, etaJet2, Wh);
     theHistograms.fill("ZZTo"+decay+"_Mjj_"+type,"",Xbins_mjj,mjj,Wh); 
     theHistograms.fill("ZZTo"+decay+"_Deta_"+type,"",Xbins_deta,deta,Wh); 
     theHistograms.fill("ZZTo"+decay+"_Phi_"+type,"",Xbins_dphi,dphi,Wh); 
   }
   if(njets>2){  
     
     ptJet3 = jetsVec->at(2).pt();
     //if (ptJet3>=300) ptJet3 = 299;
     theHistograms.fill("ZZTo"+decay+"_PtJet3_"+type, "", Xbins_ptJet3, ptJet3, Wh);
   }
}


void ZZRecoAnalyzer::FillMatrixHistosBase(std::string decay, float Wh,std::string type ){
  theHistograms.fill("ResMat_ZZTo"+decay+"_Mass_"+type,"", Xbins,Xbins, m4L,m4L_gen,Wh);
  theHistograms.fill("ResMat_ZZTo"+decay+"_PtZZ_"+type, Xbins_ptzz, Xbins_ptzz, ptzz,ptzz_gen,Wh);
  theHistograms.fill("ResMat_ZZTo"+decay+"_DphiZZ_"+type, Xbins_dphizz,Xbins_dphizz, dphizz,dphizz_gen,Wh);
  theHistograms.fill("ResMat_ZZTo"+decay+"_dRZZ_"+type,"", Xbins_drzz, Xbins_drzz,  drzz, drzz_gen,Wh);  
}

 void ZZRecoAnalyzer::FillMatrixHistosJets(std::string decay,float Wh,std::vector<phys::Jet> *jetsVec,std::vector<phys::Particle> *jetsGenVec,std::string type){


   njets_gen = jetsGenVec->size();    
   njets     = jetsVec->size(); 
   if (njets_gen>3) njets_gen=3;
   if (njets>3)     njets=3;

   theHistograms.fill("ResMat_ZZTo"+decay+"_nJets_"+type, "", Xbins_nJets,Xbins_nJets,njets,njets_gen, Wh);      
   if(njets>0) {
     stable_sort(jetsVec->begin(), jetsVec->end(), PtComparator());
     etaJet1 = fabs(jetsVec->at(0).eta());
     
     if (type.find("JESUp") != std::string::npos)     ptJet1  =  jetsVec->at(0).ptJerUp();
     else if(type.find("JESDn") != std::string::npos) ptJet1  =  jetsVec->at(0).ptJerDn();
     else if(type.find("JERUp") != std::string::npos) ptJet1  =  jetsVec->at(0).pt()*(1+jetsVec->at(0).jecUncertainty());
     else if(type.find("JERDn") != std::string::npos) ptJet1  =  jetsVec->at(0).pt()*(1-jetsVec->at(0).jecUncertainty());
     else{
       ptJet1 = jetsVec->at(0).pt(); 
     }

   }
   
   if(njets_gen>0){       
     
     stable_sort(jetsGenVec->begin(), jetsGenVec->end(), PtComparator());
     
     ptJet1_gen = genJets->at(0).pt(); 
     etaJet1_gen =  fabs(genJets->at(0).eta());
     
     //     if (ptJet1_gen >=500)  ptJet1_gen =499;
     //if (etaJet1_gen >=4.7) etaJet1_gen =4.6;
     
     theHistograms.fill("ResMat_ZZTo"+decay+"_PtJet1_"+type," ", Xbins_ptJet1,Xbins_ptJet1, ptJet1,ptJet1_gen, Wh); 
     theHistograms.fill("ResMat_ZZTo"+decay+"_EtaJet1_"+type, "", Xbins_etaJet1,Xbins_etaJet1, etaJet1,etaJet1_gen, Wh);
   }
   
   if(njets>1){ 

     deta = fabs(jetsVec->at(0).eta() - jetsVec->at(1).eta());
     mjj =  (jetsVec->at(0).p4() + jetsVec->at(1).p4()).M();
     etaJet2 = fabs(jetsVec->at(1).eta());
     dphi = physmath::deltaPhi(jetsVec->at(0).phi(),jetsVec->at(1).phi());

     if (type.find("JESUp") != std::string::npos)     ptJet2  =  jetsVec->at(1).ptJerUp();
     else if(type.find("JESDn") != std::string::npos) ptJet2  =  jetsVec->at(1).ptJerDn();
     else if(type.find("JERUp") != std::string::npos) ptJet2  =  jetsVec->at(1).pt()*(1+jetsVec->at(1).jecUncertainty());
     else if(type.find("JERDn") != std::string::npos) ptJet2  =  jetsVec->at(1).pt()*(1-jetsVec->at(1).jecUncertainty());
     else{
       ptJet2 = jetsVec->at(1).pt(); 
     }
     
     // if (deta>=4.7)    deta    = 4.6;
     // if (mjj>=800)     mjj     = 799;
     // if (ptJet2>=500)  ptJet2  = 499;
     // if (dphi>=6)      dphi    = 5;
     // if (etaJet2>=4.7) etaJet2 = 4.6;
   }

   if(njets_gen>1){  
     
     deta_gen    = fabs(genJets->at(0).eta() - genJets->at(1).eta());
     mjj_gen     = (genJets->at(0).p4() + genJets->at(1).p4()).M();
     ptJet2_gen  = genJets->at(1).pt(); 
     etaJet2_gen = fabs(genJets->at(1).eta());
     dphi        = physmath::deltaPhi(jetsGenVec->at(0).phi(),jetsGenVec->at(1).phi());     

     // if (ptJet2_gen >=500)  ptJet2_gen  = 499;
     // if (etaJet2_gen >=4.7) etaJet2_gen = 4.6;
     // if (deta_gen>=4.7)     deta_gen    = 4.6;
     // if (mjj_gen>=800)      mjj_gen     = 799;
     
     theHistograms.fill("ResMat_ZZTo"+decay+"_PtJet2_"+type, "", Xbins_ptJet2,Xbins_ptJet2, ptJet2,ptJet2_gen, Wh);     
     theHistograms.fill("ResMat_ZZTo"+decay+"_EtaJet2_"+type, "", Xbins_etaJet2, Xbins_etaJet2, etaJet2, etaJet2_gen, Wh);
     theHistograms.fill("ResMat_ZZTo"+decay+"_Mjj_" +type,"",Xbins_mjj,Xbins_mjj,mjj,mjj_gen,Wh); 
     theHistograms.fill("ResMat_ZZTo"+decay+"_Deta_"+type,"",Xbins_deta,Xbins_deta, deta,deta_gen,Wh); 
     theHistograms.fill("ResMat_ZZTo"+decay+"_Phi_" +type,"",Xbins_dphi,Xbins_dphi,dphi,dphi_gen,Wh); 
   }
 }

 
void ZZRecoAnalyzer::begin() {

  nentries = tree()->GetEntries("ZZCand.regionWord_ & (1<<26)"); 

  UpJESData_jets           = new std::vector<phys::Jet>();
  DownJESData_jets         = new std::vector<phys::Jet>();
  UpJESData_centraljets    = new std::vector<phys::Jet>();
  DownJESData_centraljets  = new std::vector<phys::Jet>();
  
  UpJER_jets               = new std::vector<phys::Jet>();
  DownJER_jets             = new std::vector<phys::Jet>();
  UpJER_centraljets        = new std::vector<phys::Jet>();
  DownJER_centraljets      = new std::vector<phys::Jet>();

  UpJES_jets               = new std::vector<phys::Jet>();
  DownJES_jets             = new std::vector<phys::Jet>();
  UpJES_centraljets        = new std::vector<phys::Jet>();
  DownJES_centraljets      = new std::vector<phys::Jet>();


  UpJESData_jets->reserve(2);           
  DownJESData_jets->reserve(2);         
  UpJESData_centraljets->reserve(2);    
  DownJESData_centraljets->reserve(2);  
  
  UpJER_jets->reserve(2);               
  DownJER_jets->reserve(2);             
  UpJER_centraljets->reserve(2);        
  DownJER_centraljets->reserve(2);      

  UpJES_jets->reserve(2);               
  DownJES_jets->reserve(2);             
  UpJES_centraljets->reserve(2);        
  DownJES_centraljets->reserve(2);      


  Xbins += 100,200,250,300,350,400,500,600,800;
  Xbins_nJets += 0,1,2,3,4;
  Xbins_ptJet1 += 30,50,100,200,300,500;
  Xbins_ptJet2 += 30,50,100,170,300;
  Xbins_ptJet3 += 30,50,100,170,300;

  Xbins_etaJet1 += 0,1.5,2.4,3.2,4.7;
  Xbins_etaJet2 += 0,1.5,3,4.7;
  Xbins_mjj     += 0.,200,400,600,1000;

  Xbins_deta += 0,1.2,2.4,3.6,4.7;
  Xbins_dphizz += 0,1.5,2.,2.25,2.5,2.75,3,3.25;
  Xbins_dphi +=  0,2,3,4; 
  Xbins_drzz += 0,1,2,3,4,5,6;
  Xbins_ptzz += 0,25,50,75,100,150,200,300;


  m4L = 0; 
  drzz= 0;
  m4L_gen = 0; 
  drzz_gen= 0;
  njets_gen = 0;
  ngencentraljets =0; 
  mjj_gen = 0;
  deta_gen = 0; 
  mjj_gen_cj = 0;
  deta_gen_cj = 0; 
  ptJet1_gen = 0;
  ptJet2_gen = 0;
  etaJet1_gen = 0;
  etaJet2_gen = 0; 
  ptzz = 0;
  ptzz_gen = 0;
  dphizz = 0; 
  dphizz_gen = 0;
  dphi = 0;
  dphi_gen = 0;
  
  njets   =0;
  ptJet1  =0;
  ptJet2  =0;
  ptJet3  =0;
  etaJet1 =0;
  etaJet2 =0;
  deta    =0;
  mjj     =0;
}


void ZZRecoAnalyzer::end( TFile &) {
  cout<<theMCInfo.analyzedEvents()<< " " << e <<endl;
  cout << "Events in the Fiducial Region " << inFiducialRegion << endl;
  
  if(region_ == phys::CR3P1F || region_ == phys::CR2P2F) {
    vector<std::string>  FinalState = {"4m","4e","2e2m"};
    vector<std::string>  Variable = {"Mass","nJets","nJets_Central","PtJet1","EtaJet1","PtJet2","EtaJet2","Deta","Mjj","Deta_Central","Mjj_Central","dRZZ", "PtZZ","DphiZZ"};
    
    for (std::vector<std::string>::iterator var = Variable.begin() ; var != Variable.end(); ++var){
      for (std::vector<std::string>::iterator it = FinalState.begin() ; it != FinalState.end(); ++it){

	TH1 *hvar =  theHistograms.get(("ZZTo"+*it+"_"+*var+"_FRVar").c_str());	
  	TH1 *h = theHistograms.get(("ZZTo"+*it+"_"+*var+"_01").c_str());
  	if(!h || !hvar) continue;
  	for(int i = 1; i<=h->GetNbinsX();i++){  
  	  Float_t Err = h->GetBinError(i);
	  //std::cout<<"Err "<<Err<<" var "<<Err*Err+hvar->GetBinContent(i)<<" tot "<<Err*Err+hvar->GetBinContent(i)<<std::endl;
  	  h->SetBinError(i,sqrt(Err*Err+hvar->GetBinContent(i)));
  	}
      }
    }  
  }
}
