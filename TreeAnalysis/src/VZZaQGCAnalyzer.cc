#include "VVXAnalysis/TreeAnalysis/interface/VZZaQGCAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using namespace std;
using namespace phys;

Int_t VZZaQGCAnalyzer::cut() {
  return 1;
}
void VZZaQGCAnalyzer::analyze(){
  
  //assignment and properties of generated bosons
  
  TLorentzVector a;
  phys::Boson<phys::Particle> z1,z2;
  foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
    if((abs(genVBParticle.daughter(0).id())==11||abs(genVBParticle.daughter(0).id())==13)&&genVBParticle.id()==23){
      theHistograms.fill("theta of generated bosons","Generated boson theta",75 ,0,3.5,genVBParticle.p4().Theta());
      theHistograms.fill("eta of generated bosons","Generated boson eta",60,-6,6,genVBParticle.eta());
      double angle1=genVBParticle.daughter(0).p4().Angle(genVBParticle.daughter(1).p4().Vect());
      theHistograms.fill("angle of generated leptons","Generated lepton angle",75,0,3.5,angle1);
      theHistograms.fill("energy of major lepton","Generated major lepton energy",200,0,2000,std::max(genVBParticle.daughter(0).e(),genVBParticle.daughter(1).e()));
      theHistograms.fill("energy of minor lepton","Generated minor lepton energy",200,0,800,std::min(genVBParticle.daughter(0).e(),genVBParticle.daughter(1).e()));
      a+=genVBParticle.p4();
     }}
  if(topology.test(0)){
    double dR=9999;
    double dR2=9999;
    double mass1=sqrt((a.E()*a.E())-(a.Px()*a.Px())-(a.Py()*a.Py())-(a.Pz()*a.Pz()));
    theHistograms.fill("mass of generated dibosons","Generated ZZ mass",78,80,600,mass1);

    foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
      if((abs(genVBParticle.daughter(0).id())==11||abs(genVBParticle.daughter(0).id())==13)&&genVBParticle.id()==23){
	if(z1.mass()==0){
	  z1=genVBParticle;}
	else{if(abs(genVBParticle.mass()-phys::ZMASS)<abs(z1.mass()-phys::ZMASS)){
	    z2=z1;
	    z1=genVBParticle;}
	  else{z2=genVBParticle;}}}}
    
    theHistograms.fill("mass of generated Z1","Generated Z1 mass",180,50,130,z1.mass());
    theHistograms.fill("mass of generated Z2","Generated Z2 mass",90,50,130,z2.mass());
    theHistograms.fill("mass of well generated Z1","Generated Z1 mass",10,80,100,z1.mass());
    theHistograms.fill("mass of well generated Z2","Generated Z2 mass",10,70,110,z2.mass());
    
    theHistograms.fill("energy of generated Z1","Generated Z1 energy",135,0,2700,z1.e());
    theHistograms.fill("energy of generated Z2","Generated Z2 energy",135,0,2700,z2.e());
    theHistograms.fill("energy of well generated Z1","Generated Z1 energy",10,0,1400,z1.e());
    theHistograms.fill("energy of well generated Z2","Generated Z2 energy",10,0,1400,z2.e());
    
    theHistograms.fill("pt of generated Z1","Generated Z1 pt",150,0,900,z1.pt());
    theHistograms.fill("pt of generated Z2","Generated Z2 pt",150,0,900,z2.pt());
    theHistograms.fill("pt of well generated Z1","Generated Z1 pt",10,0,350,z1.pt());
    theHistograms.fill("pt of well generated Z2","Generated Z2 pt",10,0,350,z2.pt());
    
    theHistograms.fill("eta of well generated Z1","Generated Z1 eta",10,-4.5,4.5,z1.eta());
    theHistograms.fill("eta of well generated Z2","Generated Z2 eta",10,-4.5,4.5,z2.eta());
    
    theHistograms.fill("energy of good Z1 major lepton","Generated Z1 major lepton's energy",10,0,1000,std::max(z1.daughter(0).e(),z1.daughter(1).e()));
    theHistograms.fill("energy of good Z1 minor lepton","Generated Z1 minor lepton's energy",10,0,300,std::min(z1.daughter(0).e(),z1.daughter(1).e()));
    theHistograms.fill("energy of good Z2 major lepton","Generated Z2 major lepton's energy",10,0,1000,std::max(z2.daughter(0).e(),z2.daughter(1).e()));
    theHistograms.fill("energy of good Z2 minor lepton","Generated Z2 minor lepton's energy",10,0,300,std::min(z2.daughter(0).e(),z2.daughter(1).e()));
    
    double ZZangle1=z1.p4().Angle(z2.p4().Vect());
    theHistograms.fill("angolo bosoni generati","Angolo bosoni generati",35,0,3.5,ZZangle1);

    //studio proprietà bosoni ricostruiti
    
    if(ZZ->first().mass()!=0&&ZZ->passFullSelection()){
      theHistograms.fill("mass of reconstructed Z1","Reconstructed Z1 mass",180,50,130,ZZ->first().mass());
      theHistograms.fill("mass of reconstructed Z2","Reconstructed Z2 mass",180,50,130,ZZ->second().mass());
      
      theHistograms.fill("pt of reconstructed Z1","Reconstructed Z1 pt",150,0,900,ZZ->first().pt());
      theHistograms.fill("pt of reconstructed Z2","Reconstructed Z2 pt",150,0,900,ZZ->second().pt());
      
      theHistograms.fill("theta of reconstructed bosons","Reconstructed boson theta",75 ,0,3.5,ZZ->first().p4().Theta());
      theHistograms.fill("theta of reconstructed bosons","Reconstructed boson theta",75 ,0,3.5,ZZ->second().p4().Theta());
      
      theHistograms.fill("energy of reconstructed Z1","Reconstructed Z1 energy",135,0,2700,ZZ->first().e());
      theHistograms.fill("energy of reconstructed Z2","Reconstructed Z2 energy",135,0,2700,ZZ->second().e());
      
      theHistograms.fill("eta of reconstructed bosons","Reconstructed boson eta",60,-6,6,ZZ->first().eta());
      theHistograms.fill("eta of reconstructed bosons","Reconstructed boson eta",60,-6,6,ZZ->second().eta());
      
      double angle2=ZZ->first().daughter(0).p4().Angle(ZZ->first().daughter(1).p4().Vect());
      theHistograms.fill("angle of reconstructed leptons","Reconstructed lepton angle",75,0,3.5,angle2);
      double angle3=ZZ->second().daughter(0).p4().Angle(ZZ->second().daughter(1).p4().Vect());
      theHistograms.fill("angle of reconstructed leptons","Reconstructed lepton angle",75,0,3.5,angle3);
      double ZZangle2=ZZ->first().p4().Angle(ZZ->second().p4().Vect());
      theHistograms.fill("angle of reconstructed bosons","Reconstructed boson angle",35,0,3.5,ZZangle2);

      TLorentzVector b=ZZ->first().p4()+ZZ->second().p4();
      double mass2=sqrt((b.E()*b.E())-(b.Px()*b.Px())-(b.Py()*b.Py())-(b.Pz()*b.Pz()));
      theHistograms.fill("mass of reconstructed dibosons","Reconstructed ZZ mass",78,80,600,mass2);

      theHistograms.fill("energy of major reconstructed leptons","Major reconstructed lepton energy",200,0,2000,std::max(ZZ->first().daughter(0).e(),ZZ->first().daughter(1).e()));
      theHistograms.fill("energy of minor reconstructed leptons","Minor reconstructd lepton energy",200,0,800,std::min(ZZ->first().daughter(0).e(),ZZ->first().daughter(1).e()));
      theHistograms.fill("energy of major reconstructed leptons","Major reconstructed lepton energy",200,0,2000,std::max(ZZ->second().daughter(0).e(),ZZ->second().daughter(1).e()));
      theHistograms.fill("energy of minor reconstructed leptons","Minor reconstructed lepton energy",200,0,800,std::min(ZZ->second().daughter(0).e(),ZZ->second().daughter(1).e()));

      //coupling generated and reconstructed bosons
      
      double enlep1,enlep2,enlep3,enlep4;
      phys::Boson<phys::Particle> z3,z4;
      foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
	if((abs(genVBParticle.daughter(0).id())==11||abs(genVBParticle.daughter(0).id())==13)&&genVBParticle.id()==23){
	  if(physmath::deltaR(genVBParticle,ZZ->first())<dR){
	    dR=physmath::deltaR(genVBParticle,ZZ->first());
	    z3=genVBParticle;
	    enlep1=std::max(genVBParticle.daughter(0).e(),genVBParticle.daughter(1).e());
	    enlep2=std::min(genVBParticle.daughter(0).e(),genVBParticle.daughter(1).e());}}}
      theHistograms.fill("dR Z1","dR Z1",100,0,0.5,dR);
      
      foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
	if((abs(genVBParticle.daughter(0).id())==11||abs(genVBParticle.daughter(0).id())==13)&&genVBParticle.id()==23){
	  if(physmath::deltaR(genVBParticle,ZZ->second())<dR2){
	    dR2=physmath::deltaR(genVBParticle,ZZ->second());
	    z4=genVBParticle;
	    enlep3=std::max(genVBParticle.daughter(0).e(),genVBParticle.daughter(1).e());
	    enlep4=std::min(genVBParticle.daughter(0).e(),genVBParticle.daughter(1).e());}}}
      theHistograms.fill("dR Z2","dR Z2",100,0,0.5,dR2);

      if(dR<0.05&&dR2<0.1&&z3!=z4){
	
	theHistograms.fill("comparison of Z1 mass","Difference between generated/reconstructed Z1 mass",40,-8,8,ZZ->first().mass()-z3.mass());
        theHistograms.fill("comparison of Z1 energy","Difference between generated/reconstructed Z1 energy",200,-40,40,ZZ->first().e()-z3.e());
        theHistograms.fill("comparison of Z1 pt","Difference between generated/reconstructed Z1 pt",150,-30,30,ZZ->first().pt()-z3.pt());
        theHistograms.fill("comparison of Z1 eta","Difference between generated/reconstructed Z1 eta",50,-0.1,0.1,ZZ->first().eta()-z3.eta());
	
        theHistograms.fill("mass of well reconstructed Z1","Reconstructed Z1 mass",10,80,100,z3.mass());
        theHistograms.fill("pt of well reconstructed Z1","Reconstructed Z1 pt",10,0,350,z3.pt());
        theHistograms.fill("energy of well reconstructed Z1","Reconstructed Z1 energy",10,0,1400,z3.e());
        theHistograms.fill("eta of well reconstructed Z1","Reconstructed Z1 eta",10,-4.5,4.5,z3.eta());
	
        theHistograms.fill("comparison of Z1 major lepton energy","Difference between generated/reconstructed Z1 major lepton energy",100,-20,20,std::max(ZZ->first().daughter(0).e(),ZZ->first().daughter(1).e())-enlep1);
        theHistograms.fill("comparison of Z1 minor lepton energy","Difference between generated/reconstructed Z1 minor lepton energy",00,-20,20,std::min(ZZ->first().daughter(0).e(),ZZ->first().daughter(1).e())-enlep2);
        theHistograms.fill("energy of well reconstructed Z1 major lepton","Reconstructed major Z1 lepton energy",10,0,1000,enlep1);
        theHistograms.fill("energy of well reconstructed Z1 minor lepton","Reconstructed minor Z1 lepton energy",10,0,300,enlep2);
       
        theHistograms.fill("comparison of Z2 mass","Difference between generated/reconstructed Z2 mass",40,-8,8,ZZ->second().mass()-z4.mass());
        theHistograms.fill("comparison of Z2 energy","Difference between generated/reconstructed Z2 energy",200,-40,40,ZZ->second().e()-z4.e());
        theHistograms.fill("comparison of Z2 pt","Difference between generated/reconstructed Z2 pt",150,-30,30,ZZ->second().pt()-z4.pt());
        theHistograms.fill("comparison of Z2 eta","Difference between generated/reconstructed Z2 eta",50,-0.1,0.1,ZZ->second().eta()-z4.eta());
	
        theHistograms.fill("mass of well reconstructed Z2","Reconstructed Z2 mass",10,70,110,z4.mass());
        theHistograms.fill("pt of well reconstructed Z2","Reconstructed Z2 pt",10,0,350,z4.pt());
        theHistograms.fill("energy of well reconstructed Z2","Reconstructed Z2 energy",10,0,1400,z4.e());
        theHistograms.fill("eta of well reconstructed Z2","Reconstructed Z2 eta",10,-4.5,4.5,z4.eta());
        theHistograms.fill("comparison of Z2 major lepton energy","Difference between generated/reconstructed Z2 major lepton energy",100,-20,20,std::max(ZZ->second().daughter(0).e(),ZZ->second().daughter(1).e())-enlep3);
        theHistograms.fill("comparison of Z2 minor lepton energy","Difference between generated/reconstructed Z2 minor lepton energy",100,-20,20,std::min(ZZ->second().daughter(0).e(),ZZ->second().daughter(1).e())-enlep4);
        theHistograms.fill("energy of well reconstructed Z2 major lepton","Reconstructed Z2 major lepton energy",10,0,1000,enlep3);
        theHistograms.fill("energy of well reconstructed Z2 minor lepton","Reconstructed Z2 minor lepton energy",10,0,300,enlep4); }
      }}
  
  //analysis of hadronic bosons

  int toomany=0;
  double dRq,dRq1,dRq2;
  TLorentzVector pquark;
  phys::Particle jet1,jet2;
  foreach(const phys::Particle genParticle,*genParticles){
    dRq=9999;
    if(abs(genParticle.id())>0&&abs(genParticle.id())<7){
      phys::Particle jet;
      foreach(const phys::Particle genJet,*genJets){
	if(physmath::deltaR(genJet,genParticle)<dRq){
	  dRq=physmath::deltaR(genJet,genParticle);
	  jet=genJet;}}
      theHistograms.fill("dR quark-jet","dR quark-jet",100,0,0.5,dRq);
      if(dRq<0.3){
	if(dRq1==0){
	  dRq1=dRq;
	  pquark+=genParticle.p4();
	  jet1=jet;}
	else if(dRq2==0&&jet.e()!=jet1.e()){
	  dRq2=dRq;
	  pquark+=genParticle.p4();
	  jet2=jet;}
	else if(dRq1!=0&&dRq2!=0){toomany+=1;}}
     }}
  
  if(dRq1!=0&&dRq2!=0&&toomany==0){
    
    phys::Boson<phys::Particle> jetboson=phys::Boson<phys::Particle>(jet1,jet2,0);
    
    double massqq=sqrt((pquark.E()*pquark.E())-(pquark.Px()*pquark.Px())-(pquark.Py()*pquark.Py())-(pquark.Pz()*pquark.Pz()));
    theHistograms.fill("mass of qq","qq mass",120,50,130,massqq);
    theHistograms.fill("mass of hadronic bosons","Hadronic boson mass",120,50,130,jetboson.mass());
    theHistograms.fill("deltam qq-jetjet","Difference between qq/jetjet mass",50,-25,25,jetboson.mass()-massqq);
    
    theHistograms.fill("energy of hadronic bosons","Hadronic boson energy",135,0,2700,jetboson.e());
    theHistograms.fill("eta of hadronic bosons","Hadronic boson eta",60,-6,6,jetboson.eta());
    theHistograms.fill("pt of hadronic bosons","Hadronic boson pt",150,0,900,jetboson.pt());
    
    if(topology.test(0)){
      
      TLorentzVector WZ1=z1.p4()+jetboson.p4();
      TLorentzVector WZ2=z2.p4()+jetboson.p4();
      TLorentzVector Z1Z2=z1.p4()+jetboson.p4();
      TLorentzVector WZZ=z1.p4()+z2.p4()+jetboson.p4();
      double massWZ1=sqrt((WZ1.E()*WZ1.E())-(WZ1.Px()*WZ1.Px())-(WZ1.Py()*WZ1.Py())-(WZ1.Pz()*WZ1.Pz()));
      double massWZ2=sqrt((WZ2.E()*WZ2.E())-(WZ2.Px()*WZ2.Px())-(WZ2.Py()*WZ2.Py())-(WZ2.Pz()*WZ2.Pz()));
      double massZ1Z2=sqrt((Z1Z2.E()*Z1Z2.E())-(Z1Z2.Px()*Z1Z2.Px())-(Z1Z2.Py()*Z1Z2.Py())-(Z1Z2.Pz()*Z1Z2.Pz()));
      double massWZZ=sqrt((WZZ.E()*WZZ.E())-(WZZ.Px()*WZZ.Px())-(WZZ.Py()*WZZ.Py())-(WZZ.Pz()*WZZ.Pz()));
      theHistograms.fill("mass of ZZ","ZZ mass",100,150,1000,massZ1Z2);
      theHistograms.fill("mass of WZ1","WZ1 mass",50,150,1000,massWZ1);
      theHistograms.fill("mass of WZ2","WZ2 mass",50,150,1000,massWZ2);
      theHistograms.fill("mass of WZ","WZ mass",100,150,1000,massWZ1);
      theHistograms.fill("mass of WZ","WZ mass",100,150,1000,massWZ2);
      theHistograms.fill("mass of tribosons","Triboson mass",50,200,1200,massWZZ); 
  
      theHistograms.fill("energy of all bosons","Total boson energy",50,0,3000,WZZ.E());
      theHistograms.fill("energy of ZZ","ZZ energy",50,0,2000,Z1Z2.E());
      theHistograms.fill("energy of WZ1","WZ1 energy",50,0,2000,WZ1.E());
      theHistograms.fill("energy of WZ2","WZ2 energy",50,0,2000,WZ2.E());
      theHistograms.fill("energy of Z1","Z1 energy",50,0,1500,z1.e());
      theHistograms.fill("energy of Z2","Z2 energy",50,0,1500,z2.e());
      theHistograms.fill("energy of W","W energy",50,0,1500,jetboson.e());
  
      double angleWZ1=z1.p4().Angle(jetboson.p4().Vect());
      double angleWZ2=z2.p4().Angle(jetboson.p4().Vect());
      double angleZZ=z1.p4().Angle(z2.p4().Vect());
      theHistograms.fill("boson relative angle","Boson relative angle",50,0,3.5,angleWZ1);
      theHistograms.fill("boson relative angle","Boson relative angle",50,0,3.5,angleWZ2);
      theHistograms.fill("boson relative angle","Boson relative angle",50,0,3.5,angleZZ);
  
      theHistograms.fill("total pt scalar sum","Scalar pt sum",100,0,1200,jetboson.pt()+z1.pt()+z2.pt());
      theHistograms.fill("total pt vector sum","Vector pt sum",100,0,300,sqrt(WZZ.Px()*WZZ.Px()+WZZ.Py()*WZZ.Py()));
  
      theHistograms.fill("pt of Z1","Z1 pt",100,0,400,z1.pt());
      theHistograms.fill("pt of Z2","Z2 pt",100,0,400,z2.pt());
      theHistograms.fill("pt of W","W pt",100,0,400,jetboson.pt());
  
      double elepZ1min= min(z1.daughter(0).e(),z1.daughter(1).e());
      double elepZ1max= max(z1.daughter(0).e(),z1.daughter(1).e());
      theHistograms.fill("energy of major Z1 leptons","Major Z1 lepton energy",50,0,2500,elepZ1max);
      theHistograms.fill("energy of minor Z1 leptons","Minor Z1 lepton energy",50,0,1000,elepZ1min);
      double elepZ2min= min(z2.daughter(0).e(),z2.daughter(1).e());
      double elepZ2max= max(z2.daughter(0).e(),z2.daughter(1).e());
      theHistograms.fill("energy of major Z2 leptons","Major Z2 lepton energy",50,0,2500,elepZ2max);
      theHistograms.fill("energy of minor Z2 leptons","Minor Z2 lepton energy",50,0,1000,elepZ2min);}
  }
}
