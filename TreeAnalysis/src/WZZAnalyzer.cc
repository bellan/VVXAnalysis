#include "VVXAnalysis/TreeAnalysis/interface/WZZAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

//#include <TString.h>

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;


using namespace phys;

Int_t WZZAnalyzer::signalCostraint() {
  int nGenHadBoson=0;
  int nGenLepBoson=0;
  if(topology.test(0)) {
    foreach(const phys::Boson<phys::Particle>& gen, *genVBParticles){

      if(abs(gen.daughter(0).id()) < 10 && nGenHadBoson<1 && ((gen.id()==23 && gen.mass()>60 && gen.mass()<120) || (abs(gen.id())==24 && gen.mass()>50 && gen.mass()<110)) && fabs(gen.daughter(0).pt())>30 &&  fabs(gen.daughter(1).pt())>30 && fabs(gen.daughter(0).eta())<2.5 && fabs(gen.daughter(1).eta())<2.5)
	nGenHadBoson++;

      if((abs(gen.daughter(0).id())==11 || abs(gen.daughter(0).id())==13) &&
	 gen.id()==23 && gen.mass()>60 && gen.mass()<120 &&
	 fabs(gen.daughter(0).eta())<2.5 && fabs(gen.daughter(1).eta())<2.5 && gen.daughter(0).pt()>5 && gen.daughter(1).pt()>5 &&
	 nGenLepBoson<2)
	nGenLepBoson++;
    }
  }
  if(nGenHadBoson!=1 || nGenLepBoson!=2)
    return 0;

  return 1;

}

Bool_t WZZAnalyzer::cut(Int_t n,  phys::Boson<phys::Jet> recoV) { //returns false if the event has to be cut
  
  switch(n){
    
  case 1:
    if(jets->size()>1 && recoV.mass()>50 && recoV.mass()<120 &&
       fabs(recoV.daughter(0).pt())>30 && fabs(recoV.daughter(1).pt())>30 && fabs(recoV.daughter(0).eta())<2.5 && fabs(recoV.daughter(1).eta())<2.5 &&
       ZZ->first().daughter(0).pt()>5 && ZZ->second().daughter(0).pt()>5 && ZZ->first().daughter(1).pt()>5 && ZZ->second().daughter(1).pt()>5 &&
       ZZ->first().daughter(0).eta()<2.5 && ZZ->second().daughter(0).eta()<2.5 && ZZ->first().daughter(1).eta()<2.5 && ZZ->second().daughter(1).eta()<2.5)
       return true;
    break;


  case 2: 
    if(recoV.mass()>65)
      return true;

    break;

  case 3: 
    if(recoV.mass()<105)
      return true;

    break;

  default:
    return true;

  }
  return false;
}





void WZZAnalyzer::analyze(){ //It's the only member function running each event.

  //cout << "----------------------------------------------------------------"<<endl;
  //cout << "Run: " << run << " event: " << event << endl;

  //genAnalyze();
  
  phys::Boson<phys::Jet> recoV;
  Reconstruct(&recoV);
  

  //---------------ALL-THE-EVENTS---------------//


  
  foreach(const phys::Boson<phys::Particle>& gen, *genVBParticles){
    if(abs(gen.id())==24 && abs(gen.daughter(0).id()) < 10) {
      theHistograms.fill("genVMass_WZZ_den","mass of gen W",60,0.1,150,gen.mass());
      theHistograms.fill("genVTot_WZZ_den","mass of gen W",1,0.1,150,gen.mass());
    }
    if(abs(gen.id())==23 && abs(gen.daughter(0).id()) < 10) {
      theHistograms.fill("genVMass_ZZZ_den","mass of gen Z",60,0.1,150,gen.mass());
      theHistograms.fill("genVTot_ZZZ_den","mass of gen Z",1,0.1,150,gen.mass());
    }
  }
  //printHistos(0,"all", recoV);



  //---------------SIGNAL-EVENTS---------------//
  
  if(signalCostraint()==1){
    foreach(const phys::Boson<phys::Particle>& gen, *genVBParticles){
      if(abs(gen.id())==24 && abs(gen.daughter(0).id()) < 10) {
	theHistograms.fill("genVMass_WZZ_num_sign","mass of gen W",60,0.1,150,gen.mass(), theWeight);
	theHistograms.fill("genVTot_WZZ_num_sign","mass of gen W",1,0.1,150,gen.mass(), theWeight);
      }
      if(abs(gen.id())==23 && abs(gen.daughter(0).id()) < 10) {
	theHistograms.fill("genVMass_ZZZ_num_sign","mass of gen Z",60,0.1,150,gen.mass(), theWeight);
	theHistograms.fill("genVTot_ZZZ_num_sign","mass of gen Z",1,0.1,150,gen.mass(), theWeight);
      }

    }

    printHistos(0,"sign", recoV);
    genAnalyze();
  }//not signalCostraint anymore 
    
  
  //---------------BACKGROUND-EVENTS---------------//
  
  else{
    foreach(const phys::Boson<phys::Particle>& gen, *genVBParticles){
      if(abs(gen.id())==24 && abs(gen.daughter(0).id()) < 10) {
	theHistograms.fill("genVMass_WZZ_num_bckg","mass of gen W",60,0.1,150,gen.mass(), theWeight);
	theHistograms.fill("genVTot_WZZ_num_bckg","mass of gen W",1,0.1,150,gen.mass(), theWeight);
      }
      if(abs(gen.id())==23 && abs(gen.daughter(0).id()) < 10) {
	theHistograms.fill("genVMass_ZZZ_num_bckg","mass of gen Z",60,0.1,150,gen.mass(), theWeight);
	theHistograms.fill("genVTot_ZZZ_num_bckg","mass of gen Z",1,0.1,150,gen.mass(), theWeight);
      }

    }
  
    printHistos(0,"bckg", recoV);

  }


}


void WZZAnalyzer::genAnalyze(){

  //cout << "-----------------------------------------------------------------"<<endl;
  //cout << "Run: " << run << " event: " << event << endl;





  //----------------------------------------GenParticles----------------------------------------//
  
  int ngenElectrons=0; 
  int ngenMuons=0;
  
  foreach(const phys::Particle& gen, *genParticles){
    if(abs(gen.id()) == 13 ){//for muons: mc pdg id = 13 by notation
      theHistograms.fill("genMuonPt","pt of gen muons",20,0,200,gen.pt(),theWeight);
      theHistograms.fill("genMuonEta","eta of gen muons",30,0,2.5,fabs(gen.eta()),theWeight);
      theHistograms.fill("genMuonPhi","phi of gen muons",30,-3.2,3.2,gen.phi(),theWeight);
      theHistograms.fill("genMuonRapidity","rapidity of gen muons",30,0,2.5,fabs(gen.rapidity()),theWeight);
      theHistograms.fill("genMuonEnergy","energy of gen muons",120,0,400,fabs(gen.e()),theWeight);

      ngenMuons++;
    }


    
    if(abs(gen.id()) == 11 ){//for electrons: mc pdg id = 11 by notation
      theHistograms.fill("genElectronPt","pt of gen electrons",20,0,200,gen.pt(),theWeight);
      theHistograms.fill("genElectronEta","eta of gen electrons",30,0,2.5,fabs(gen.eta()),theWeight);
      theHistograms.fill("genElectronPhi","phi of gen electrons",30,-3.2,3.2,gen.phi(),theWeight);
      theHistograms.fill("genElectronRapidity","rapidity of gen electrons",30,0,2.5,fabs(gen.rapidity()),theWeight);
      theHistograms.fill("genElectronEnergy","energy of gen electrons",120,0,400,fabs(gen.e()),theWeight);
      
      ngenElectrons++;
    }
   
  }
  theHistograms.fill("ngenElectrons","Number of genElectrons",  10, 0, 10 , ngenElectrons, theWeight);
  theHistograms.fill("ngenMuons","Number of genMuons",  10, 0, 10 , ngenMuons, theWeight);
  





  
  //----------------------------------------MUONS-------------------------------------------//


  theHistograms.fill("nmuons","Number of muons",  10, 0, 10 , muons->size(), theWeight);

  
  foreach(const phys::Lepton& mu, *muons){
    theHistograms.fill("muonsPt","pt of muons",20,0,200,mu.pt());
    theHistograms.fill("muonsPtEta","eta vs pt of muons",60,0,200,30,0,2.5, mu.pt(),fabs(mu.eta()));
    theHistograms.fill("muonsPtRapidity","rapidity vs pt of muons",60,0,200,30,0,2.5, mu.pt(),fabs(mu.rapidity()));
    theHistograms.fill("muonsPhiEta","eta vs phi of muons",60,-3.2,3.2,30,0,2.5, mu.phi(),fabs(mu.eta()));
    theHistograms.fill("muonsPhiPt","pt vs phi of muons",60,-3.2,3.2,60,0,200, mu.phi(),mu.pt());

    
  }


  //----------------------------------------ELECTRONS----------------------------------------//

  
  theHistograms.fill("nelectrons","Number of electrons",  10, 0, 10 , electrons->size(), theWeight);

  
  foreach(const phys::Lepton& e, *electrons){
    theHistograms.fill("electronsPt","pt of electrons",20,0,200,e.pt());
    theHistograms.fill("electronsPtEta","eta vs pt of electrons",60,0,200,30,0,2.5, e.pt(),fabs(e.eta()));
    theHistograms.fill("electronsPtRapidity","rapidity vs pt of electrons",60,0,200,30,0,2.5, e.pt(),fabs(e.rapidity()));
    theHistograms.fill("electronsPhiEta","eta vs phi of electrons",60,-3.2,3.2,30,0,2.5, e.phi(),fabs(e.eta()));
    theHistograms.fill("electronsPhiPt","pt vs phi of electrons",60,-3.2,3.2,60,0,200, e.phi(),e.pt());

    
  }

  theHistograms.fill("nelectronsmuons","Number of electrons and muons",  10, 0, 10 , 10, 0 , 10 , electrons->size(), muons->size(), theWeight);





  //----------------------------------------JETS----------------------------------------//


  foreach(const phys::Particle& gen, *genJets){
    theHistograms.fill("genJetPt","pt of gen jets",20,0,200,gen.pt(),theWeight);
    theHistograms.fill("genJetEta","eta of gen jets",30,0,2.5,fabs(gen.eta()),theWeight);
    theHistograms.fill("genJetPhi","phi of gen jets",30,-3.2,3.2,gen.phi(),theWeight);
    theHistograms.fill("genJetRapidity","rapidity of gen jets",30,0,2.5,fabs(gen.rapidity()),theWeight);


      
  }

  theHistograms.fill("njets","Number of jets",  10, 0, 10 , jets->size(), theWeight);

  
  foreach(const phys::Jet& j, *jets){
    theHistograms.fill("jetsPt","pt of jets",20,30,200,j.pt());
    theHistograms.fill("jetsEta","eta of jets",30,0,4.7,fabs(j.eta()));
    theHistograms.fill("jetsPtEta","eta vs pt of jets",20,30,200,30,0,4.7, j.pt(),fabs(j.eta()));
    theHistograms.fill("jetsRapidity","rapidity of jets",30,0,5,fabs(j.rapidity()));
    theHistograms.fill("jetsPtRapidity","rapidity vs pt of jets",20,30,200,30,0,5, j.pt(),fabs(j.rapidity()));
    theHistograms.fill("jetsPhi","phi of jets",60,-3.2,3.2,j.phi());
    theHistograms.fill("jetsPtPhi","phi vs pt of jets",20,30,200,60,-3.2,3.2, j.pt(),j.phi());
    theHistograms.fill("jetsEtaPhi","phi vs eta of jets",30,0,4.7,60,-3.2,3.2,fabs(j.eta()),j.phi());




  }


  //----------------------------------------JETS(AK8)----------------------------------------//

  foreach(const phys::Particle& gen, *genJetsAK8){
    theHistograms.fill("genJetAK8Pt","pt of gen jetsAK8",20,0,200,gen.pt(),theWeight);
    theHistograms.fill("genJetAK8Eta","eta of gen jetsAK8",30,0,2.5,fabs(gen.eta()),theWeight);
    theHistograms.fill("genJetAK8Phi","phi of gen jetsAK8",30,-3.2,3.2,gen.phi(),theWeight);
    theHistograms.fill("genJetAK8Rapidity","rapidity of gen jetsAK8",30,0,2.5,fabs(gen.rapidity()),theWeight);

      
  }


  foreach(const phys::Jet& j, *jetsAK8){
    theHistograms.fill("jetsAK8Pt","pt of jetsAK8",20,30,200,j.pt());
    theHistograms.fill("jetsAK8Eta","eta of jetsAK8",30,0,4.7,fabs(j.eta()));
    theHistograms.fill("jetsAK8PtEta","eta vs pt of jetsAK8",20,30,200,30,0,4.7, j.pt(),fabs(j.eta()));
    theHistograms.fill("jetsAK8Rapidity","rapidity of jetsAK8",30,0,5,fabs(j.rapidity()));
    theHistograms.fill("jetsAK8PtRapidity","rapidity vs pt of jetsAK8",20,30,200,30,0,5, j.pt(),fabs(j.rapidity()));
    theHistograms.fill("jetsAK8Phi","phi of jetsAK8",60,-3.2,3.2,j.phi());
    theHistograms.fill("jetsAK8PtPhi","phi vs pt of jetsAK8",20,30,200,60,-3.2,3.2, j.pt(),j.phi());
    theHistograms.fill("jetsAK8EtaPhi","phi vs eta of jetsAK8",30,0,4.7,60,-3.2,3.2,fabs(j.eta()),j.phi());

  }



  //----------------------------------------MET_ANALYSIS----------------------------------------//

 
  theHistograms.fill("met","MET",  100, 0, 500 , met->pt(), theWeight);

  TLorentzVector pTot = TLorentzVector(0., 0., 0., 0.);
  double genMet=0.0;
 
  foreach(const phys::Particle& gen, *genParticles)//neutrinos
    if(abs(gen.id()) == 12 || abs(gen.id()) == 14 || abs(gen.id()) == 16)
      pTot+=gen.p4();

  genMet = sqrt(pTot.Px()*pTot.Px()+pTot.Py()*pTot.Py());
  theHistograms.fill("genMet","genMET",  100, 0, 300 , genMet , theWeight);

  theHistograms.fill("metDeviation","deviation of MET",  100, -150, 150, genMet-met->pt(), theWeight);
 




  //----------------------------------------RESOLUTION&EFFICIENCY----------------------------------------//



  bool electronMatched[electrons->size()];
  for(uint k=0; k<electrons->size(); k++)
    electronMatched[k]=false;

  bool muonMatched[muons->size()];
  for(uint k=0; k<electrons->size(); k++)
    muonMatched[k]=false;

  bool genParticleMatched;

  foreach(const phys::Particle& gen, *genParticles){

    double dptInv=0;
    double dpt=0;
    double deta=0;
    double dE=0;
  
    double dRminEl=99.9;
    double dRminMu=99.9;

    double dRmax=0.1;

    double dREl[ngenElectrons][electrons->size()];
    double dRMu[ngenMuons][muons->size()];

    int j=0;
    int i=0;//only for electrons
    int m=0;//only for muons
    int electronIndex=0;
    int muonIndex=0;

    
    genParticleMatched=false;



    //-----electrons-----//


    if(abs(gen.id()) == 11){//electrons mc pdg id
          
      theHistograms.fill("genElectronsEta_den","eta of gen electrons",30,0,2.5,fabs(gen.eta()));//for efficiency
      theHistograms.fill("genElectronsPt_den","pt of gen electrons",20,0,200,gen.pt());

      theHistograms.fill("genElectronsPhi_den","phi of gen electrons",30,-3.2,3.2,gen.phi());
      theHistograms.fill("genElectronsE_den","energy of gen electrons",100,0,400,gen.e());


      theHistograms.fill("genElectronsPtTot_den","pt of gen electrons",1,0,200,gen.pt());
      theHistograms.fill("genElectronsETot_den","energy of gen electrons",1,0,400,gen.e());
     

     
      foreach(const phys::Lepton& e, *electrons){
	dREl[i][j]=physmath::deltaR(gen, e);
	if(dREl[i][j]<dRminEl && !electronMatched[j]){
	  dRminEl=dREl[i][j];
	  electronIndex=j;
	  genParticleMatched=true;
	  dptInv=1/gen.pt()-1/e.pt();
	  dpt=gen.pt()-e.pt();
	  deta=gen.eta()-e.eta();
	  dE=gen.e()-e.e();
	}
	j++;
      }
      electronMatched[electronIndex]=true;
      i++;
      if(genParticleMatched){
	theHistograms.fill("electronsDeltaR","dR of e",  50, 0, 2 , dRminEl, theWeight);
	if(dRminEl<=dRmax){
	  theHistograms.fill("electronsCutDeltaR","cut dR of e",  50, 0, 0.05 , dRminEl, theWeight);
	  theHistograms.fill("genElectronsEta_num","eta of gen electrons",30,0,2.5,fabs(gen.eta()));
	  theHistograms.fill("genElectronsPt_num","pt of gen electrons",20,0,200,gen.pt());
	  theHistograms.fill("genElectronsPhi_num","phi of gen electrons",30,-3.2,3.2,gen.phi());
	  theHistograms.fill("genElectronsE_num","energy of gen electrons",100,0,400,gen.e());



	  theHistograms.fill("genElectronsPtTot_num","pt of gen electrons",1,0,200,gen.pt());
	  theHistograms.fill("genElectronsETot_num","energy of gen electrons",1,0,400,gen.e());

	}
      }


      if(fabs(dpt)>0){
	theHistograms.fill("electronsDeltaPt","dpt of e",  50, -10, 10 , dpt, theWeight);
	theHistograms.fill("electronsAbsResolution(pt)","absResolution(pt) of e",  50, -0.2, 0.2, dpt/gen.pt() , theWeight);
	theHistograms.fill("electronsAbsResolution(pt)vsPt","absResolution(pt) of e vs pt",  50, -0.2, 0.2 ,  10, 0, 20, dpt/gen.pt(), gen.pt(), theWeight);
      }
      if(fabs(dptInv)>0){
	theHistograms.fill("electronsDeltaPtInverse","dptInv of e",  50, -0.2, 0.2 , dptInv, theWeight);
	theHistograms.fill("electronsAbsResolution(ptInv)","absResolution(ptInv) of e",  50, -0.2, 0.2 , dptInv/(1/gen.pt()), theWeight);
      }

      if(fabs(deta)>0){
	theHistograms.fill("electronsDeltaEta","deta of e",  50, -10, 10, deta, theWeight);
	theHistograms.fill("electronsAbsResolution(eta)","absResolution(eta) of e",  50, -0.2, 0.2 , deta/gen.eta(), theWeight);
      }
      if(fabs(dE)>0){
	theHistograms.fill("electronsDeltaE","dE of e",  50, -30, 30 , dE, theWeight);
	theHistograms.fill("electronsAbsResolution(E)","absResolution(E) of e",  50, -0.2, 0.2 , dE/gen.e(), theWeight);
      }
    }

    
      //-----muons-----//
     
    else if(abs(gen.id()) == 13){//muons mc pdg id
      theHistograms.fill("genMuonsEta_den","eta of gen muons",30,0,2.5,fabs(gen.eta()));
      theHistograms.fill("genMuonsPt_den","pt of gen muons",20,0,200,gen.pt());
      theHistograms.fill("genMuonsPhi_den","phi of gen muons",30,-3.2,3.2,gen.phi());
      theHistograms.fill("genMuonsE_den","energy of gen muons",120,0,400,gen.e());


      theHistograms.fill("genMuonsPtTot_den","pt of gen muons",1,0,200,gen.pt());
      theHistograms.fill("genMuonsETot_den","energy of gen muons",1,0,400,gen.e());


      
      foreach(const phys::Lepton& mu, *muons){
	dRMu[m][j]=physmath::deltaR(gen, mu);
	if(dRMu[m][j]<dRminMu && !muonMatched[j]){
	  dRminMu=dRMu[m][j];
	  muonIndex=j;
	  genParticleMatched=true;
	  dptInv=1/gen.pt()-1/mu.pt();
	  dpt=gen.pt()-mu.pt();
	  deta=gen.eta()-mu.eta();
	  dE=gen.e()-mu.e();
	}
	j++;
      }
      muonMatched[muonIndex]=true;
      m++;
      if(genParticleMatched){
	theHistograms.fill("muonsDeltaR","dR of mu",  50, 0, 2 , dRminMu, theWeight);
	if(dRminMu<=dRmax){
	  theHistograms.fill("muonsCutDeltaR","cut dR of mu",  50, 0, 0.05 , dRminMu, theWeight);
	  theHistograms.fill("genMuonsEta_num","eta of gen muons",30,0,2.5,fabs(gen.eta()));
	  theHistograms.fill("genMuonsPt_num","pt of gen muons",20,0,200,gen.pt());
	  theHistograms.fill("genMuonsPhi_num","phi of gen muons",30,-3.2,3.2,gen.phi());
	  theHistograms.fill("genMuonsE_num","energy of gen muons",120,0,400,gen.e());



	  theHistograms.fill("genMuonsPtTot_num","pt of gen muons",1,0,200,gen.pt());
	  theHistograms.fill("genMuonsETot_num","energy of gen muons",1,0,400,gen.e());

	}
      }
     
      if(fabs(dpt)>0){
	theHistograms.fill("muonsDeltaPt","dpt of mu",  50, -10, 10 , dpt, theWeight);
	theHistograms.fill("muonsAbsResolution(pt)","absResolution(pt) of mu",  50, -0.2, 0.2, dpt/gen.pt() , theWeight);
	theHistograms.fill("muonsAbsResolution(pt)vsPt","absResolution(pt) of mu vs pt",  50, -0.2, 0.2 ,  10, 0, 20, dpt/gen.pt(), gen.pt(), theWeight);
      }
      if(fabs(dptInv)>0){
	theHistograms.fill("muonsDeltaPtInverse","dptInv of mu",  50, -0.2, 0.2 , dptInv, theWeight);
	theHistograms.fill("muonsAbsResolution(ptInv)","absResolution(ptInv) of mu",  50, -0.2, 0.2 , dptInv/(1/gen.pt()), theWeight);
      }

      if(fabs(deta)>0){
	theHistograms.fill("muonsDeltaEta","deta of mu",  50, -10, 10, deta, theWeight);
	theHistograms.fill("muonsAbsResolution(eta)","absResolution(eta) of mu",  50, -0.2, 0.2 , deta/gen.eta(), theWeight);
      }
      if(fabs(dE)>0){
	theHistograms.fill("muonsDeltaE","dE of mu",  50, -30, 30 , dE, theWeight);
	theHistograms.fill("muonsAbsResolution(E)","absResolution(E) of mu",  50, -0.2, 0.2 , dE/gen.e(), theWeight);
      }
    }

      

  }



  //----------------------------------------VECTOR_BOSONS_RECONSTRUCTION----------------------------------------//

  TLorentzVector genP4Tot = TLorentzVector(0., 0. , 0., 0.);
  
  int nGenHadBoson=0;
  int nGenLepBoson=0;
  if(topology.test(0)){
    foreach(const phys::Boson<phys::Particle>& gen, *genVBParticles){
      if(abs(gen.daughter(0).id()) < 10 && nGenHadBoson<1 && ((gen.id()==23 && gen.mass()>60 && gen.mass()<120) || (gen.id()==24 && gen.mass()>50 && gen.mass()<110))  && (gen.daughter(0)).pt()>30 &&  (gen.daughter(1)).pt()>30){
	genP4Tot+=gen.p4();
	nGenHadBoson++;
      }
      

      if((abs(gen.daughter(0).id())==11 || abs(gen.daughter(0).id())==13) && gen.id()==23 && nGenLepBoson<2){
	genP4Tot+=gen.p4();
	nGenLepBoson++;
      }

    }
    phys::Particle genZZjj(genP4Tot);    
    theHistograms.fill("genZZjjPt","Pt of gen system ZZjj",40,0,200,genZZjj.pt());
    
  }
  

  
  //Declaration of different algorithm's candidate
  phys::Boson<phys::Jet> mWCandidate;
  phys::Boson<phys::Jet> mZCandidate;
  phys::Boson<phys::Jet> maxVPtCandidate;
  phys::Boson<phys::Jet> minTotPtCandidate;
  phys::Boson<phys::Jet> mWZCandidate;

  
  //Building of every jets pairs combination
  std::vector<phys::Boson<phys::Jet> > DiJets;

  for(uint i=0; i</*centralJets*/jets->size(); i++)//Warning: size can be 0
    for(uint j=i; j</*centralJets*/jets->size(); j++)  
      if(j!=i)
	DiJets.push_back(phys::Boson<phys::Jet>(/*centralJets*/jets->at(i), /*centralJets*/jets->at(j)));

  if(topology.test(0)){  //For the signal definition, the categories with the first bit on (i.e. bit0=true), are the only ones that have to be considered

    if(/*centralJets*/jets->size()>1){	
  
      //1st reconstruction model: comparison with WMass
      std::stable_sort(DiJets.begin(), DiJets.end(), phys::MassComparator(phys::WMASS));
      mWCandidate = DiJets.at(0);

      //2nd reconstruction model: comparison with ZMass
      std::stable_sort(DiJets.begin(), DiJets.end(), phys::MassComparator(phys::ZMASS));
      mZCandidate = DiJets.at(0);
    
      //3rd reconstruction model: maximization of candidate Pt
      std::stable_sort(DiJets.begin(), DiJets.end(), phys::ScalarSumPtComparator());
      maxVPtCandidate=DiJets.at(0);

      //4th reconstruction model: minimization of total Pt of ZZjj system
      std::vector<phys::Particle> ZZjj;
      phys::Particle ZZjjCandidate;
      for(uint i=0; i<DiJets.size(); i++){
	phys:: Particle totState(ZZ->p4()+(DiJets.at(i)).p4());
	ZZjj.push_back(totState.p4());
      }
      std::stable_sort(ZZjj.begin(), ZZjj.end(), phys::PtComparator());
      ZZjjCandidate=ZZjj.back();
      for (uint i=0; i<DiJets.size(); i++)
	if((DiJets.at(i)).p4()==(ZZjjCandidate.p4()-ZZ->p4()))
	  minTotPtCandidate = DiJets.at(i);

      //5th reconstruction model: comparison with a mean value between ZMass and WMass
      std::stable_sort(DiJets.begin(), DiJets.end(), phys::MassComparator(0.2*phys::ZMASS+0.8*phys::WMASS));
      mWZCandidate = DiJets.at(0);
      
      foreach(const phys::Boson<phys::Particle>& gen, *genVBParticles){
	if(abs(gen.daughter(0).id()) < 10){
	  if(abs(gen.id())==23  && gen.mass()>60 && gen.mass()<120 && (gen.daughter(0)).pt()>30 &&  (gen.daughter(1)).pt()>30){//for Z bosons: mc pdg id = 23 by notation
	    theHistograms.fill("genZPt","pt of gen Z",20,30,200,gen.pt(),theWeight);
	    theHistograms.fill("genZEta","eta of gen Z",30,0,2.5,fabs(gen.eta()),theWeight);
	    theHistograms.fill("genZPhi","phi of gen Z",30,-3.2,3.2,gen.phi(),theWeight);
	    theHistograms.fill("genZRapidity","rapidity of gen Z",30,0,2.5,fabs(gen.rapidity()),theWeight);
	    theHistograms.fill("genZEnergy","energy of gen Z",120,0,400,fabs(gen.e()),theWeight);
	    theHistograms.fill("genZmass","mass of gen Z",40,40,130,fabs(gen.mass()),theWeight);
	    theHistograms.fill("genJetDeltaPhi","dPhi of gen jets",30,0,3.2,fabs(physmath::deltaPhi(gen.daughter(0).phi(), gen.daughter(1).phi())),theWeight);
	    theHistograms.fill("genZDaughtersMass","mass of gen Z daughters",100,0,100, (gen.daughter(0)).mass(),theWeight);
	    theHistograms.fill("genZDaughtersMass","mass of gen Z daughters",100,0,100, (gen.daughter(1)).mass(),theWeight);


	    CompatibilityTest(mWCandidate, gen, "ZZZ", "mW");
	    CompatibilityTest(mZCandidate, gen, "ZZZ", "mZ");
	    CompatibilityTest(maxVPtCandidate, gen, "ZZZ", "maxVPt");
	    CompatibilityTest(minTotPtCandidate, gen, "ZZZ", "minTotPt");
	    CompatibilityTest(mWZCandidate, gen, "ZZZ", "mWZ");

	  }
	  else if(abs(gen.id())==24 && gen.mass()>50 && gen.mass()<110 && (gen.daughter(0)).pt()>30 &&  (gen.daughter(1)).pt()>30){
	    theHistograms.fill("genWPt","pt of gen W",20,30,200,gen.pt(),theWeight);
	    theHistograms.fill("genWEta","eta of gen W",30,0,2.5,fabs(gen.eta()),theWeight);
	    theHistograms.fill("genWPhi","phi of gen W",30,-3.2,3.2,gen.phi(),theWeight);
	    theHistograms.fill("genWRapidity","rapidity of gen W",30,0,2.5,fabs(gen.rapidity()),theWeight);
	    theHistograms.fill("genWEnergy","energy of gen W",120,0,400,fabs(gen.e()),theWeight);
	    theHistograms.fill("genWmass","mass of gen W",40,40,130,fabs(gen.mass()),theWeight);
  	    theHistograms.fill("genJetDeltaPhi","dPhi of gen jets",30,0,3.2,fabs(physmath::deltaPhi(gen.daughter(0).phi(), gen.daughter(1).phi())),theWeight);
	    theHistograms.fill("genWDaughtersMass","mass of gen W daughters",100,0,100, (gen.daughter(0)).mass(),theWeight);
	    theHistograms.fill("genWDaughtersMass","mass of gen W daughters",100,0,100, (gen.daughter(1)).mass(),theWeight);


	    CompatibilityTest(mWCandidate, gen, "WZZ", "mW");
	    CompatibilityTest(mZCandidate, gen, "WZZ", "mZ");
	    CompatibilityTest(maxVPtCandidate, gen, "WZZ", "maxVPt");
	    CompatibilityTest(minTotPtCandidate, gen, "WZZ", "minTotPt");
	    CompatibilityTest(mWZCandidate, gen, "WZZ", "mWZ");

	  }

	}
      }

    
    }



  }//topology.test(0) closed
}




void WZZAnalyzer::Reconstruct(phys::Boson<phys::Jet>* mWCandidate){
  //Building of every jets pairs combination
  std::vector<phys::Boson<phys::Jet> > DiJets;

  for(uint i=0; i<jets->size(); i++)//Warning: size can be 0
    for(uint j=i; j<jets->size(); j++)  
      if(j!=i)
	DiJets.push_back(phys::Boson<phys::Jet>(jets->at(i), jets->at(j)));
  /*
  if(topology.test(0)){  //For the signal definition, the categories with the first bit on (i.e. bit0=true), are the only ones that have to be considered
  */
  //Deleted because it may be used after the signalCostraint() running 

  if(jets->size()>1){	
      std::stable_sort(DiJets.begin(), DiJets.end(), phys::MassComparator(phys::WMASS));
      *mWCandidate = DiJets.at(0);
  }
}

void WZZAnalyzer::CompatibilityTest(phys::Boson<phys::Jet> bestCandidate, phys::Boson<phys::Particle> genVB, std::string sample, std::string algorithm){
  double dRMax=0.4;
  std::string genVBId;

  if(genVB.id()==24) // && genVB.mass()>50 && genVB.mass()<110 && (genVB.daughter(0)).mass()>30 &&  (genVB.daughter(1)).mass()>30)
    genVBId="W";
  else if(genVB.id()==23) // && genVB.mass()>60 && genVB.mass()<120 && (genVB.daughter(0)).mass()>30 &&  (genVB.daughter(1)).mass()>30)
    genVBId="Z";
  else return;

  if(genVBId+"ZZ"==sample){
    
    theHistograms.fill("gen"+genVBId+"Eta_"+sample+"_"+algorithm+"_den","eta of gen "+genVBId,20,0,2.5,fabs(genVB.eta()));
    theHistograms.fill("gen"+genVBId+"Pt_"+sample+"_"+algorithm+"_den","pt of gen "+genVBId,20,0,200,genVB.pt());
    theHistograms.fill("gen"+genVBId+"Phi_"+sample+"_"+algorithm+"_den","phi of gen "+genVBId,30,-3.2,3.2,genVB.phi());
    theHistograms.fill("gen"+genVBId+"E_"+sample+"_"+algorithm+"_den","energy of gen "+genVBId,20,0,400,genVB.e());
    theHistograms.fill("gen"+genVBId+"Mass_"+sample+"_"+algorithm+"_den","mass of gen "+genVBId,10,50,120,genVB.mass());
    theHistograms.fill("gen"+genVBId+"Tot_"+sample+"_"+algorithm+"_den","tot gen "+genVBId,1,0,200,genVB.mass());

    
    double dRJ0d0=physmath::deltaR(bestCandidate.daughter(0),genVB.daughter(0));
    double dRJ0d1=physmath::deltaR(bestCandidate.daughter(0),genVB.daughter(1));
    double dRJ1d0=physmath::deltaR(bestCandidate.daughter(1),genVB.daughter(0));
    double dRJ1d1=physmath::deltaR(bestCandidate.daughter(1),genVB.daughter(1));

    if(dRJ1d0+dRJ0d1<dRJ0d0+dRJ1d1){
      theHistograms.fill(genVBId+"jets_deltaR","deltaR of "+genVBId+" reco", 50, 0, 1, dRJ1d0);
      theHistograms.fill(genVBId+"jets_deltaR","deltaR of "+genVBId+" reco", 50, 0, 1, dRJ0d1);
      if(dRJ1d0<dRMax && dRJ0d1<dRMax){
	theHistograms.fill("gen"+genVBId+"Eta_"+sample+"_"+algorithm+"_num","eta of gen "+genVBId,20,0,2.5,fabs(genVB.eta()));
	theHistograms.fill("gen"+genVBId+"Pt_"+sample+"_"+algorithm+"_num","pt of gen "+genVBId,20,0,200,genVB.pt());
	theHistograms.fill("gen"+genVBId+"Phi_"+sample+"_"+algorithm+"_num","phi of gen "+genVBId,30,-3.2,3.2,genVB.phi());
	theHistograms.fill("gen"+genVBId+"E_"+sample+"_"+algorithm+"_num","energy of gen "+genVBId,20,0,400,genVB.e());
	theHistograms.fill("gen"+genVBId+"Mass_"+sample+"_"+algorithm+"_num","mass of gen "+genVBId,10,50,120,genVB.mass());
	theHistograms.fill("gen"+genVBId+"Tot_"+sample+"_"+algorithm+"_num","tot gen "+genVBId,1,0,200,genVB.mass());

      }
    }
    else{
      theHistograms.fill(genVBId+"jets_deltaR","deltaR of "+genVBId+" reco", 50, 0, 1, dRJ0d0);
      theHistograms.fill(genVBId+"jets_deltaR","deltaR of "+genVBId+" reco", 50, 0, 1, dRJ1d1);
      if(dRJ1d1<dRMax && dRJ1d1<dRMax){
	theHistograms.fill("gen"+genVBId+"Eta_"+sample+"_"+algorithm+"_num","eta of gen "+genVBId,20,0,2.5,fabs(genVB.eta()));
	theHistograms.fill("gen"+genVBId+"Pt_"+sample+"_"+algorithm+"_num","pt of gen "+genVBId,20,0,200,genVB.pt());
	theHistograms.fill("gen"+genVBId+"Phi_"+sample+"_"+algorithm+"_num","phi of gen "+genVBId,30,-3.2,3.2,genVB.phi());
	theHistograms.fill("gen"+genVBId+"E_"+sample+"_"+algorithm+"_num","energy of gen "+genVBId,20,0,400,genVB.e());
	theHistograms.fill("gen"+genVBId+"Mass_"+sample+"_"+algorithm+"_num","mass of gen "+genVBId,10,50,120,genVB.mass());
	theHistograms.fill("gen"+genVBId+"Tot_"+sample+"_"+algorithm+"_num","tot gen "+genVBId,1,0,200,genVB.mass());

      }
    }
  }
  return;
}

void WZZAnalyzer::printHistos(uint i, std::string histoType,  phys::Boson<phys::Jet> recoV){
  
  std::vector<std::string> cuts ={"0","1","2","3"};
  
  if(i<cuts.size() && cut(i, recoV)){
    theHistograms.fill("recoVMass_"+histoType+cuts.at(i),"mass of recoV",30,0.1,200, recoV.mass(), (theWeight));
    theHistograms.fill("recoVTot_"+histoType+cuts.at(i),"mass of recoV",1,0.1,200, recoV.mass(), (theWeight));

    theHistograms.fill("recoVDaughter0Pt_"+histoType+cuts.at(i),"pt of recoVDaughter0",50,0,600, recoV.daughter(0).pt(), (theWeight));
    theHistograms.fill("recoVDaughter1Pt_"+histoType+cuts.at(i),"pt of recoVDaughter1",50,0,600, recoV.daughter(1).pt(), (theWeight));


    
    theHistograms.fill("recoZZMass_"+histoType+cuts.at(i),"mass of recoZZ",25,0.1,400, ZZ->mass(), (theWeight));

    theHistograms.fill("recoVPt_"+histoType+cuts.at(i),"pt of recoV",50,0,300,recoV.pt(), (theWeight));
    theHistograms.fill("recoVEta_"+histoType+cuts.at(i),"eta of recoV",30,0,3.5,fabs(recoV.eta()), (theWeight));
    theHistograms.fill("recoVPhi_"+histoType+cuts.at(i),"phi of recoV",30,-3.2,3.2,recoV.phi(), (theWeight));
    theHistograms.fill("recoVEnergy_"+histoType+cuts.at(i),"energy of  recoV",60,0,2400,fabs(recoV.e()), (theWeight));
    theHistograms.fill("recoVDaughtersDeltaPhi_"+histoType+cuts.at(i),"dPhi of recoVDaughters",30,0,3.2,fabs(physmath::deltaPhi(recoV.daughter(0).phi(), recoV.daughter(1).phi())), (theWeight));

    theHistograms.fill("recoZZPt_"+histoType+cuts.at(i),"pt of recoZZ",50,0,600,ZZ->pt(), (theWeight));
    theHistograms.fill("recoZZEta_"+histoType+cuts.at(i),"eta of recoZZ",70,0,3.5,fabs(ZZ->eta()), (theWeight));
    theHistograms.fill("recoZZEnergy_"+histoType+cuts.at(i),"energy of  recoZZ",120,0,400,fabs(ZZ->e()), (theWeight));
    theHistograms.fill("recoZZDeltaPhi_"+histoType+cuts.at(i),"dPhi of recoZZ",30,0,3.2,fabs(physmath::deltaPhi(ZZ->first().phi(), ZZ->second().phi())), (theWeight));

    theHistograms.fill("recoZZDeltaPhi_vs_recoV_"+histoType+cuts.at(i),"phi of recoZZ vs recoV",100,0.0,3.14,fabs(physmath::deltaPhi(ZZ->first().phi(),recoV.daughter(0).phi())),100,0.0,3.14, fabs(physmath::deltaPhi(ZZ->second().phi(), recoV.phi())) ,theWeight);

    
    theHistograms.fill("recoFirstZVDeltaPhi_"+histoType+cuts.at(i),"phi of recoZV",100,0.0,3.14,fabs(physmath::deltaPhi(ZZ->first().phi(),(recoV.daughter(0).phi()+recoV.phi()))), theWeight);

        
    theHistograms.fill("recoSecondZVDeltaPhi_"+histoType+cuts.at(i),"phi of recoZV",100,0.0,3.14,fabs(physmath::deltaPhi(ZZ->second().phi(),(recoV.daughter(0).phi()+recoV.phi()))), theWeight);
 
    printHistos(++i, histoType, recoV);
      
 
  }
  

  return;
}





