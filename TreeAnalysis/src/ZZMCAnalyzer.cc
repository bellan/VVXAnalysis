#include "VVXAnalysis/TreeAnalysis/interface/ZZMCAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include <boost/assign/std/vector.hpp> 
#include <boost/assert.hpp> 
using namespace std;
using namespace boost::assign; // bring 'operator+=()' into scope

#include <boost/foreach.hpp>
#include <sstream> 
#include <string> 
#define foreach BOOST_FOREACH

using std::cout;
using std::endl;


using namespace phys;


void ZZMCAnalyzer::analyze(){
  
  PreCounter+=1;
  if((region_ == phys::MC && topology.test(2)) || (region_ == phys::MC_HZZ && topology.test(0) ) ){       
    
    int z1 = abs(genVBParticles->at(0).daughter(0).id()); 
    int z2 = abs(genVBParticles->at(1).daughter(0).id()); 
    
    std::string decay="None";                                                                                                                                                                                
    
    if((z1==11 && z2==13) || (z1==13 && z2==11)) {decay = "2e2m";} 
    else if(z1==13 && z2==13) {decay = "4m";} 
    else if(z1==11 && z2==11) {decay = "4e";}
    else {cout<<"Wrong decay, check z doughters: Z0 l0"<<genVBParticles->at(0).daughter(0).id()<<" Z0 l1 "<<genVBParticles->at(0).daughter(1).id()<<" Z1 l0 "<<genVBParticles->at(1).daughter(0).id()<<" Z1 l1 "<<genVBParticles->at(1).daughter(1).id()<<endl; abort();} 


    isTightFr =kFALSE;
    
    nEvent ++;
    
    if(decay == "None"){
      std::cout<<"Check decay channel"<<std::endl;
      return;
    }
    
    sample = "01";
    if(PreCounter < nentries/2) {sample = "0";} 
    else {sample = "1";}
    
    
    if (m4L_gen>=800) m4L_gen = 799;
    
    drzz_gen =physmath::deltaR(genVBParticles->at(0),genVBParticles->at(1));
    //if(drzz_gen>6) drzz_gen = 5.9; //overflow bin
    
    njets = genJets->size();
    if (njets>3) njets=3;
    
    ncentraljets = centralGenJets->size();
    if (ncentraljets>3) ncentraljets=3;
    
    ptzz_gen =  (genVBParticles->at(0).p4()+genVBParticles->at(1).p4()).Pt();
    if (ptzz_gen>300) ptzz_gen=299;
    
    dphizz_gen = 0;
    dphizz_gen = fabs(physmath::deltaPhi(genVBParticles->at(0).phi(),genVBParticles->at(1).phi())); 
    
    
    if(njets >=1){
      ptjet1_gen = genJets->at(0).pt();
      if(ptjet1_gen>=500) ptjet1_gen=499;
      etajet1_gen = fabs(genJets->at(0).eta());
      if(etajet1_gen>=4.7) etajet1_gen=4.6;
      
      if(njets>=2){  
	
	deta_gen = fabs(genJets->at(0).eta() - genJets->at(1).eta());
	mjj_gen =  (genJets->at(0).p4() + genJets->at(1).p4()).M();
	ptjet2_gen = genJets->at(1).pt();
	etajet2_gen = fabs(genJets->at(1).eta());
	
	if (deta_gen>=4.7) deta_gen = 4.6;
	if (mjj_gen>=800) mjj_gen = 799;
	if(ptjet2_gen>=500) ptjet2_gen=499; 
	if(etajet2_gen>=4.7) etajet2_gen=4.6;
      }
    }
    
    if(ncentraljets>=2){  
      
      deta_gen_cj = fabs(centralGenJets->at(0).eta() - centralGenJets->at(1).eta());
      mjj_gen_cj =  (centralGenJets->at(0).p4() + centralGenJets->at(1).p4()).M();
      if (deta_gen_cj>=4.7) deta_gen_cj = 4.6;
      if (mjj_gen_cj>=800) mjj_gen_cj = 799;
    }
    scaleFacErrSq = ZZ->efficiencySFUnc();
    // scaleFacMuErrSq = ZZ->muEffSFUnc(); 
    // scaleFacEleErrSq = ZZ->eleEffSFUnc(); 
   
    ScalVarVal = {LHEweight_QCDscale_muR1_muF1,LHEweight_QCDscale_muR1_muF2,LHEweight_QCDscale_muR1_muF0p5,LHEweight_QCDscale_muR2_muF1,LHEweight_QCDscale_muR2_muF2,LHEweight_QCDscale_muR0p5_muF1,LHEweight_QCDscale_muR0p5_muF0p5};
    
    
    w_kf = 1; 
    // if((theMCInfo.fileName()=="ggZZ2e2mu") || (theMCInfo.fileName()=="ggZZ4e") || (theMCInfo.fileName()=="ggZZ4mu"))   w_kf = kFactor_ggZZ ; 
    // else if((theMCInfo.fileName()=="ZZTo4l") || (theMCInfo.fileName()=="ZZTo4lamcatnlo")) w_kf = kFactor_qqZZM * kFactor_EWKqqZZ ; 
    if((theMCInfo.fileName()=="ggZZ2e2mu") || (theMCInfo.fileName()=="ggZZ4e") || (theMCInfo.fileName()=="ggZZ4mu"))  w_kf = 1.7 ; 
    else if(theMCInfo.fileName()=="ZZTo4l") w_kf = 1.1; 
    
    m4L_gen  = sqrt((genVBParticles->at(0).p4()+genVBParticles->at(1).p4())*(genVBParticles->at(0).p4()+genVBParticles->at(1).p4()));
    if (m4L_gen>=800) m4L_gen = 799;

  
  // theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_pdfUp_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,Wh*LHEweight_PDFVariation_Up);
  // theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_pdfDn_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,Wh*LHEweight_PDFVariation_Dn);
  
  // theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_AsMzUp_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,Wh*LHEweight_AsMZ_Up);
  // theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_AsMzDn_01", std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,Wh*LHEweight_AsMZ_Dn);

    
    FillHistosBase(decay,theMCInfo.sampleWeight()*w_kf,"Gen_01");
    FillHistosBase(decay,theMCInfo.sampleWeight()*w_kf,"Gen_"+sample);

    FillHistosBase(decay,theMCInfo.sampleWeight()*w_kf*LHEweight_PDFVariation_Up,"Gen_01_pdfUp");
    FillHistosBase(decay,theMCInfo.sampleWeight()*w_kf*LHEweight_PDFVariation_Dn,"Gen_01_pdfDn");
    FillHistosBase(decay,theMCInfo.sampleWeight()*w_kf*LHEweight_AsMZ_Up,"Gen_01_asMZUp");
    FillHistosBase(decay,theMCInfo.sampleWeight()*w_kf*LHEweight_AsMZ_Dn,"Gen_01_asMZDn");


    if((region_ == phys::MC && regionWord.test(26)) || ((region_ == phys::MC_HZZ) && regionWord.test(3))){
      FillHistosBase(decay,theWeight*w_kf,"GenReco_"+sample);
      FillHistosBase(decay,theWeight*w_kf,"GenReco_01");
      FillHistosBase(decay,theWeight*w_kf*(1-scaleFacErrSq),"GenRecoSFErrSqMinus_01");
      FillHistosBase(decay,theWeight*w_kf*(1+scaleFacErrSq),"GenRecoSFErrSqPlus_01");
    }
    if((region_ == phys::MC && topology.test(3)) || ( region_ == phys::MC_HZZ && topology.test(1)) ){
      FillHistosBase(decay,theMCInfo.sampleWeight()*w_kf,"Gen_01_fr");
      FillHistosBase(decay,theMCInfo.sampleWeight()*w_kf,"Gen_"+sample+"_fr");

      FillHistosBase(decay,theMCInfo.sampleWeight()*w_kf*LHEweight_PDFVariation_Up,"Gen_01_fr_pdfUp");
      FillHistosBase(decay,theMCInfo.sampleWeight()*w_kf*LHEweight_PDFVariation_Dn,"Gen_01_fr_pdfDn");
      FillHistosBase(decay,theMCInfo.sampleWeight()*w_kf*LHEweight_AsMZ_Up,"Gen_01_fr_asMZUp");
      FillHistosBase(decay,theMCInfo.sampleWeight()*w_kf*LHEweight_AsMZ_Dn,"Gen_01_fr_asMZDn");

      if((region_ == phys::MC && regionWord.test(26)) || ((region_ == phys::MC_HZZ) && regionWord.test(3))){
	FillHistosBase(decay,theWeight*w_kf,"GenReco_01_fr");
	FillHistosBase(decay,theWeight*w_kf,"GenReco_"+sample+"_fr");
	FillHistosBase(decay,theWeight*w_kf*(1-scaleFacErrSq),"GenRecoSFErrSqMinus_01_fr");
	FillHistosBase(decay,theWeight*w_kf*(1+scaleFacErrSq),"GenRecoSFErrSqPlus_01_fr");
      }
    }
  }
}  

void ZZMCAnalyzer::FillHistosBase(std::string decay, float Wh,std::string type ){

  // theHistograms.fill(std::string("ZZTo")+decay+"_MassGen_"+sample, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,theMCInfo.sampleWeight()*w_kf);

  theHistograms.fill("ZZTo"+decay+"_Mass"+type,"", Xbins, m4L_gen,Wh);
  theHistograms.fill("ZZTo"+decay+"_PtZZ"+type, Xbins_ptzz, ptzz_gen,Wh);
  theHistograms.fill("ZZTo"+decay+"_DphiZZ"+type, Xbins_dphizz, dphizz_gen,Wh);
  theHistograms.fill("ZZTo"+decay+"_dRZZ"+type,"", Xbins_drzz, drzz_gen,Wh);
  theHistograms.fill(std::string("ZZTo")+decay+"_nJets"+type, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,Xbins_nJets,njets,Wh);    
  theHistograms.fill(std::string("ZZTo")+decay+"_nJets_Central"+type, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,Xbins_nJets,ncentraljets,Wh);    
  
  //  cout<<endl;
  if(type=="Gen_01" || type=="Gen_01_fr"){ 
    for(unsigned int i =0; i<ScalVar.size(); i++){    
      // cout<<ScalVar.at(i)<<" "<<ScalVarVal.at(i)<<endl;
      theHistograms.fill(std::string("ZZTo"+decay+"_Mass"+type+ScalVar.at(i)), std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, m4L_gen,Wh*ScalVarVal.at(i));
      theHistograms.fill(std::string("ZZTo"+decay+"_PtZZ"+type+ScalVar.at(i)), std::string("Generated  PtZZ #rightarrow ")+decay , Xbins_ptzz, ptzz_gen,Wh*ScalVarVal.at(i));
      theHistograms.fill(std::string("ZZTo"+decay+"_DphiZZ"+type+ScalVar.at(i)), std::string("Generated #Delta #Phi ZZ #rightarrow")+decay , Xbins_dphizz, dphizz_gen,Wh*ScalVarVal.at(i));
      theHistograms.fill(std::string("ZZTo"+decay+"_dRZZ"+type+ScalVar.at(i)), std::string("Generated  #Delta R  ZZ #rightarrow")+decay , Xbins_dphizz, drzz_gen,Wh*ScalVarVal.at(i));      
      theHistograms.fill(std::string("ZZTo")+decay+"_nJets"+type+ScalVar.at(i), std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,Xbins_nJets,njets,Wh*ScalVarVal.at(i));    
      theHistograms.fill(std::string("ZZTo")+decay+"_nJets_Central"+type+ScalVar.at(i), std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,Xbins_nJets,ncentraljets,Wh*ScalVarVal.at(i));    

    }
  }
  
  if(njets >=1){
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1"+type,"",Xbins_ptjet1,ptjet1_gen,Wh);  
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1"+type,"",Xbins_etajet1,etajet1_gen,Wh);  
        
    if(type=="Gen_01" || type=="Gen_01_fr"){ 
      for(unsigned int i =0; i<ScalVar.size(); i++){    
	theHistograms.fill(std::string("ZZTo")+decay+"_PtJet1"+type+ScalVar.at(i),"",Xbins_ptjet1,ptjet1_gen,Wh*ScalVarVal.at(i));  
	theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet1"+type+ScalVar.at(i),"",Xbins_etajet1,etajet1_gen,Wh*ScalVarVal.at(i));  
      }
    }
  } 
  
  if(njets>=2){  
    
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj"+type, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,Wh);  
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta"+type, std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,Wh);  
    theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2"+type,"",Xbins_ptjet2,ptjet2_gen,Wh);  
    theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2"+type,"",Xbins_etajet2,etajet2_gen,Wh);  
    
    
    if(type=="Gen_01" || type=="Gen_01_fr"){ 
      for(unsigned int i =0; i<ScalVar.size(); i++){    
	
	theHistograms.fill(std::string("ZZTo")+decay+"_Mjj"+type+ScalVar.at(i), std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen,Wh*ScalVarVal.at(i));  
	theHistograms.fill(std::string("ZZTo")+decay+"_Deta"+type+ScalVar.at(i), std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen,Wh*ScalVarVal.at(i)); 
	theHistograms.fill(std::string("ZZTo")+decay+"_PtJet2"+type+ScalVar.at(i),"",Xbins_ptjet2,ptjet2_gen,Wh*ScalVarVal.at(i));  
	theHistograms.fill(std::string("ZZTo")+decay+"_EtaJet2"+type+ScalVar.at(i),"",Xbins_etajet2,etajet2_gen,Wh*ScalVarVal.at(i));  
      }
    }
  }
  if(ncentraljets>=2){  
    
    theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_Central"+type, std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,Wh);     
    theHistograms.fill(std::string("ZZTo")+decay+"_Deta_Central"+type, std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,Wh);  
    
    if(type=="Gen_01" || type=="Gen_01_fr"){ 
      for(unsigned int i =0; i<ScalVar.size(); i++){    
	
	theHistograms.fill(std::string("ZZTo")+decay+"_Mjj_Central"+type+ScalVar.at(i), std::string("m_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_mjj,mjj_gen_cj,Wh*ScalVarVal.at(i));
	theHistograms.fill(std::string("ZZTo")+decay+"_Deta_Central"+type+ScalVar.at(i), std::string("#Delta#eta_{jj} of ZZ_{1}#rightarrow ")+decay,Xbins_deta,deta_gen_cj,Wh*ScalVarVal.at(i));  
	
      }
    }
  }
}



void ZZMCAnalyzer::begin() {
  nentries =  tree()->GetEntries();
  PreCounter = 0;
  nEvent = 0;
  inFiducialRegion=0;

  Xbins_nJets += 0,1,2,3,4;
  Xbins += 100,200,250,300,350,400,500,600,800; 
  Xbins_deta += 0,2.4,4.7;
  Xbins_mjj += 0.,200,800;
  Xbins_ptjet1 += 30,50,100,200,300,500;
  Xbins_ptjet2 += 30,100,200,500;
  Xbins_etajet1 += 0,1.5,3,4.7;
  Xbins_etajet2 += 0,1.5,3,4.7; 
  Xbins_drzz += 0,1,2,3,4,5,6;
  Xbins_ptzz += 0,25,50,75,100,150,200,300;
  Xbins_dphizz += 0,1.5,2.,2.25,2.5,2.75,3,3.25;
  ScalVarVal = {};

  ScalVar = {"_mf1mr1","_mf1mr2","_mf1mr0p5","_mf2mr1","_mf2mr2","_mf0p5mr1","_mf0p5mr0p5"};
  
  m4L_gen = 0;
  njets = 0;
  mjj_gen = 0;
  deta_gen = 0;
  ncentraljets = 0;
  mjj_gen_cj = 0;
  deta_gen_cj = 0;
  ptjet1_gen = 0;
  ptjet2_gen = 0; 
  etajet1_gen = 0;
  etajet2_gen = 0; 
  drzz_gen =0;
  ptzz_gen =0;
  dphizz_gen =0;

  isTightFr =kFALSE;

  region = "";
  sample = "01";
}

void ZZMCAnalyzer::end( TFile &) {
  cout <<"Tree Entries"<<nentries<< endl;
  cout <<"events in the fiducial region"<<inFiducialRegion<< endl;

  std::vector<double> binsVec;
  
  vector<TH1F>  ScalHistos = {};
  vector<std::string>  FinalState = {"4m","4e","2e2m"};
  vector<std::string>  Variable = {"Mass","nJets","nJets_Central","PtJet1","EtaJet1","PtJet2","EtaJet2","Deta","Mjj","Deta_Central","Mjj_Central","dRZZ", "PtZZ","DphiZZ"};
  vector<std::string>  Regions = {"","_fr"};
  
  for (std::vector<std::string>::iterator var = Variable.begin() ; var != Variable.end(); ++var){

    if(*var=="Mass")                binsVec=Xbins;
    else if(*var=="nJets")          binsVec=Xbins_nJets;
    else if(*var=="nJets_Central")  binsVec=Xbins_nJets;
    else if(*var=="PtJet1")         binsVec=Xbins_ptjet1;
    else if(*var=="EtaJet1")        binsVec=Xbins_etajet1;
    else if(*var=="PtJet2")         binsVec=Xbins_ptjet2;
    else if(*var=="EtaJet2")        binsVec=Xbins_etajet2;
    else if(*var=="Deta")           binsVec=Xbins_deta;
    else if(*var=="Mjj")            binsVec=Xbins_mjj;
    else if(*var=="Deta_Central")   binsVec=Xbins_deta;
    else if(*var=="Mjj_Central")    binsVec=Xbins_mjj;
    else if(*var=="dRZZ")           binsVec=Xbins_drzz;
    else if(*var=="PtZZ")           binsVec=Xbins_ptzz;
    else if(*var=="DphiZZ")         binsVec=Xbins_dphizz;
    else {cout<<"Error: wron variable in endjob"<<endl; abort();}
   

    for (std::vector<std::string>::iterator it = FinalState.begin() ; it != FinalState.end(); ++it){
      for (std::vector<std::string>::iterator rg = Regions.begin() ; rg != Regions.end(); ++rg){	
	//	cout<<"var "<<*var<<" fn "<<*it<<" rg "<<*rg<<endl;
	for (std::vector<std::string>::iterator scal = ScalVar.begin() ; scal != ScalVar.end(); ++scal){

	  TH1 *h = theHistograms.get(("ZZTo"+*it+"_"+*var+"Gen_01"+*rg+*scal).c_str());
	  if(!h){
	    continue;
	  }
	  ScalHistos.push_back(  *dynamic_cast<TH1F *>(h) );
	}

	theHistograms.fill(std::string("ZZTo")+*it+"_"+*var+"Gen_01_scaleUp"+*rg, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+*it , binsVec, 0,0);
	theHistograms.fill(std::string("ZZTo")+*it+"_"+*var+"Gen_01_scaleDn"+*rg, std::string("Generated invariant mass of ZZ_{1}#rightarrow ")+*it , binsVec, 0,0);
	TH1 * hvarUp = theHistograms.get(("ZZTo"+*it+"_"+*var+"Gen_01_scaleUp"+*rg).c_str());
	TH1 * hvarDn = theHistograms.get(("ZZTo"+*it+"_"+*var+"Gen_01_scaleDn"+*rg).c_str());

	for(int b = 1; b<=hvarUp->GetNbinsX();b++){  
	  bool isFirst = true;
	  for (std::vector<TH1F>::iterator hvar = ScalHistos.begin() ; hvar != ScalHistos.end(); ++hvar){  
	    if(isFirst){
	      hvarUp->SetBinContent(b,hvar->GetBinContent(b));
	      hvarDn->SetBinContent(b,hvar->GetBinContent(b));
	      isFirst = false;
	    }
	    else{
	      //	      cout<<"var "<<hvar->GetBinContent(b)<<" vardn "<<hvarDn->GetBinContent(b)<<" "<<endl;
	      if(hvar->GetBinContent(b) > hvarUp->GetBinContent(b)) hvarUp->SetBinContent(b,hvar->GetBinContent(b));
	      if(hvar->GetBinContent(b) < hvarDn->GetBinContent(b)) hvarDn->SetBinContent(b,hvar->GetBinContent(b));
	    }
	  }
	}
      ScalHistos.clear();
      }
    }
  }
}

 
