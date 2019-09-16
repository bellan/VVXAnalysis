#include "VVXAnalysis/TreeAnalysis/interface/WZZAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;


using namespace phys;

Int_t WZZAnalyzer::cut() {
  
  return 1;
}



void WZZAnalyzer::analyze(){

  cout << "-----------------------------------------------------------------"<<endl;
  cout << "Run: " << run << " event: " << event << endl;





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
  //double genMetX=0.0;
  //double genMetY=0.0;
 
  foreach(const phys::Particle& gen, *genParticles)//neutrinos
    if(abs(gen.id()) == 12 || abs(gen.id()) == 14 || abs(gen.id()) == 16)
      pTot+=gen.p4();
 
  //genMet = sqrt(genMetX*genMetX+genMetY*genMetY);
  genMet = sqrt(pTot.Px()*pTot.Px()+pTot.Py()*pTot.Py());
  theHistograms.fill("genMet","genMET",  100, 0, 300 , genMet , theWeight);

  theHistograms.fill("metDeviation","deviation of MET",  100, -150, 150, genMet-met->pt(), theWeight);
 




  //----------------------------------------RESOLUTION&EFFICIENCY----------------------------------------//



  bool electronMatched[electrons->size()];
  for(int k=0; k<electrons->size(); k++)
    electronMatched[k]=false;

  bool muonMatched[muons->size()];
  for(int k=0; k<electrons->size(); k++)
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
      theHistograms.fill("genElectronsPt_den","pt of gen electrons",20,0,200,gen.pt(),theWeight);
      theHistograms.fill("genElectronsPhi_den","phi of gen electrons",30,-3.2,3.2,gen.phi(),theWeight);
      theHistograms.fill("genElectronsE_den","energy of gen electrons",100,0,400,gen.e(),theWeight);
     

     
      foreach(const phys::Lepton& e, *electrons){
	dREl[i][j]=deltaR(gen.eta(), gen.phi(), e.eta(), e.phi());
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
	  theHistograms.fill("genElectronsPt_num","pt of gen electrons",20,0,200,gen.pt(),theWeight);
	  theHistograms.fill("genElectronsPhi_num","phi of gen electrons",30,-3.2,3.2,gen.phi(),theWeight);
	  theHistograms.fill("genElectronsE_num","energy of gen electrons",100,0,400,gen.e(),theWeight);
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
      theHistograms.fill("genMuonsPt_den","pt of gen muons",20,0,200,gen.pt(),theWeight);
      theHistograms.fill("genMuonsPhi_den","phi of gen muons",30,-3.2,3.2,gen.phi(),theWeight);
      theHistograms.fill("genMuonsE_den","energy of gen muons",90,0,400,gen.e(),theWeight);

      foreach(const phys::Lepton& mu, *muons){
	dRMu[m][j]=deltaR(gen.eta(), gen.phi(), mu.eta(), mu.phi());
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
	  theHistograms.fill("genMuonsPt_num","pt of gen muons",20,0,200,gen.pt(),theWeight);
	  theHistograms.fill("genMuonsPhi_num","phi of gen muons",30,-3.2,3.2,gen.phi(),theWeight);
	  theHistograms.fill("genMuonsE_num","energy of gen muons",90,0,400,gen.e(),theWeight);
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
      

}
