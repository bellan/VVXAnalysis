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

  
  //This part of the analysis can be used both for WZZ and for ZZZ samples. 
  //For the signal definition, the categories with the first bit on (i.e. bit0=true), are the only ones that have to be considered

  //vector <phys::Boson<phys::Jet>> ZDecay;
  //vector <phys::Boson<phys::Jet>> WDecay;

  
  std::vector <bool>jetsCoupleMatchedW((jets->size()*(jets->size()-1))/2,false);//parallel vectors
  std::vector <bool>jetsCoupleMatchedZ(jetsCoupleMatchedW);

  std::vector <std::vector <int>> jetsCouplesW((jets->size()*(jets->size()-1))/2, std::vector <int>(2) );
  /*int k=0;
  for(int i=0; i<jets->size()-1; i++)
    for(int j=i+1; j<jets->size(); j++){
      jetsCouplesW[k][0]=i;
      jetsCouplesW[k][1]=j;
      k++;
    }
  */
  std::vector <std::vector <int>> jetsCouplesZ(jetsCouplesW);

  /*  
  int n=0;
  for(int i=0; i<jets->size()-1; i++)
    for(int j=i+1; j<jets->size(); j++){
      jetsCouplesZ[n][0]=i;
      jetsCouplesZ[n][1]=j;
      n++;
    }
  */
  std::vector <std::vector <double>> WMatching;//matrix of compatibility parameters
  std::vector <std::vector <double>> ZMatching;
  
  double etaTot, mTot,  dEta ,dm, dEtaMax=30, dmMax=20, etaWeight=0.3, mWeight=0.7/*,  compatibility */ ;
  
  double WBosonIndex=0, WJetsCoupleIndex=0;
  double ZBosonIndex=0, ZJetsCoupleIndex=0;
  
  if(topology.test(0)){
       
    foreach(const phys::Boson<phys::Particle>& gen, *genVBParticles){

      int Zj0Counter=0, Zj1Counter=0;

      //-------------Vector Boson Z-------------//
      
      if(abs(gen.id()) == 23 && abs(gen.daughter(0).id()) < 10){//for Z bosons: mc pdg id = 23 by notation
	theHistograms.fill("genZPt","pt of gen Z",20,0,200,gen.pt(),theWeight);
	theHistograms.fill("genZEta","eta of gen Z",30,0,2.5,fabs(gen.eta()),theWeight);
	theHistograms.fill("genZPhi","phi of gen Z",30,-3.2,3.2,gen.phi(),theWeight);
	theHistograms.fill("genZRapidity","rapidity of gen Z",30,0,2.5,fabs(gen.rapidity()),theWeight);
	theHistograms.fill("genZEnergy","energy of gen Z",120,0,400,fabs(gen.e()),theWeight);
	theHistograms.fill("genZmass","mass of gen Z",40,40,130,fabs(gen.mass()),theWeight);

	foreach(const phys::Jet& j0, *jets){
	  foreach(const phys::Jet& j1, *jets){

	    if(Zj1Counter>Zj0Counter){

	      etaTot=j0.eta()+j1.eta();
	      dEta=fabs((gen.eta()-etaTot)/gen.eta());

	      mTot=j0.mass()+j1.mass();
	      dm=fabs((mTot-phys::ZMASS)/phys::ZMASS);

	      ZMatching[ZJetsCoupleIndex][ZBosonIndex]=1-(etaWeight*dEta+mWeight*dm);

	      if(ZMatching[ZJetsCoupleIndex][ZBosonIndex]<0)
		ZMatching[ZJetsCoupleIndex][ZBosonIndex]=0;
	      
	      jetsCouplesZ[ZJetsCoupleIndex][0]=Zj0Counter;
	      jetsCouplesZ[ZJetsCoupleIndex][1]=Zj1Counter;
	      ZJetsCoupleIndex++;
	      
	    }

	    Zj1Counter++;

	  }	  

	  Zj0Counter++;

	}

	ZBosonIndex++;
    
      }

      


      


      

      //-------------Vector Boson W-------------//

      
      if(abs(gen.id()) == 24 && abs(gen.daughter(0).id()) < 10){//for W bosons: mc pdg id = 24 by notation
	theHistograms.fill("genWPt","pt of gen W",20,0,200,gen.pt(),theWeight);
	theHistograms.fill("genWEta","eta of gen W",30,0,2.5,fabs(gen.eta()),theWeight);
	theHistograms.fill("genWPhi","phi of gen W",30,-3.2,3.2,gen.phi(),theWeight);
	theHistograms.fill("genWRapidity","rapidity of gen W",30,0,2.5,fabs(gen.rapidity()),theWeight);
	theHistograms.fill("genWEnergy","energy of gen W",120,0,400,fabs(gen.e()),theWeight);
	theHistograms.fill("genWmass","mass of gen W",40,40,130,fabs(gen.mass()),theWeight);

      }
      
    }


    //-------------Vector Boson Z best jets couple-------------//

    //Best matching searching

    int  ZMatchedIndex, jetsCoupleMatchedIndex;
    std::vector<std::vector<int>> bestZJetsCoupleIndex(ZBosonIndex, std::vector<int>(2));//each boson there is the best jets pair
    std::vector<bool>ZMatched(ZBosonIndex);
    //std::vector<std::vector<int>>bestZJetsCoupleIndex(ZBosonIndex,<std::vector<int>(2));

    double bestMatching;
    for(int count=0; count<ZBosonIndex; count++){//i.e. while there is at least one ZBoson not matched yet
      bestMatching=-1;
      for(int i=0; i<ZBosonIndex; i++)
	for(int j=0; j<ZJetsCoupleIndex; j++)
	  if(!jetsCoupleMatchedZ[j] && !ZMatched[i] && ZMatching[j][i]>bestMatching){
	    bestMatching=ZMatching[j][i];
	    ZMatchedIndex=i;
	    jetsCoupleMatchedIndex=j;
	  }
      
      bestZJetsCoupleIndex[ZMatchedIndex][0]=jetsCouplesZ[jetsCoupleMatchedIndex][0];//setting matrix of Bosons-Jets (index) 
      bestZJetsCoupleIndex[ZMatchedIndex][1]=jetsCouplesZ[jetsCoupleMatchedIndex][1];

      ZMatched[ZMatchedIndex]=true;

      jetsCoupleMatchedZ[jetsCoupleMatchedIndex]=true;
      for(int i=0; i<ZJetsCoupleIndex; i++)
	for(int j=0; j<2; j++)
	  if(jetsCouplesZ[i][j]==jetsCouplesZ[jetsCoupleMatchedIndex][0] || jetsCouplesZ[i][j]==jetsCouplesZ[jetsCoupleMatchedIndex][1])
	    jetsCoupleMatchedZ[i]=true;
    }
    


    
  }

  

}
