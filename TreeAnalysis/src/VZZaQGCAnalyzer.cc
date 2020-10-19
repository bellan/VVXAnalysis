#include "VVXAnalysis/TreeAnalysis/interface/VZZaQGCAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;

using namespace phys;

Int_t VZZaQGCAnalyzer::cut() {
  return 1;
}
void VZZaQGCAnalyzer::analyze(){
  //assegnazione e proprietà bosoni generati
   TLorentzVector a;
   phys::Boson<phys::Particle> z1,z2;
   foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
     if((abs(genVBParticle.daughter(0).id())==11||abs(genVBParticle.daughter(0).id())==13)&&genVBParticle.id()==23){
      theHistograms.fill("theta bosoni generati","Theta bosoni generati",75 ,0,3.5,genVBParticle.p4().Theta());
      theHistograms.fill("eta bosoni generati","Eta bosoni generati",60,-6,6,genVBParticle.eta());
      double angolo1=genVBParticle.daughter(0).p4().Angle(genVBParticle.daughter(1).p4().Vect());
      theHistograms.fill("angolo leptoni generati","Angolo leptoni generati",75,0,3.5,angolo1);
      theHistograms.fill("E leptone maggiore","Energia leptone piu' energetico",200,0,2000,std::max(genVBParticle.daughter(0).e(),genVBParticle.daughter(1).e()));
      theHistograms.fill("E leptone minore","Energia leptone meno energetico",200,0,800,std::min(genVBParticle.daughter(0).e(),genVBParticle.daughter(1).e()));
      a+=genVBParticle.p4();
     }}
   if(topology.test(0)){
      double dR=9999;
      double dR2=9999;
      double mass1=sqrt((a.E()*a.E())-(a.Px()*a.Px())-(a.Py()*a.Py())-(a.Pz()*a.Pz()));
       theHistograms.fill("massa dibosoni generati","Massa ZZ generati",78,80,600,mass1);
       
	 foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
	   if((abs(genVBParticle.daughter(0).id())==11||abs(genVBParticle.daughter(0).id())==13)&&genVBParticle.id()==23){
	   if(z1.mass()==0){
	     z1=genVBParticle;}
	   else{if(abs(genVBParticle.mass()-phys::ZMASS)<abs(z1.mass()-phys::ZMASS)){
	       z2=z1;
	       z1=genVBParticle;}
	     else{z2=genVBParticle;}}}}
	    theHistograms.fill("massa Z1 generati","Massa Z1 generati",180,50,130,z1.mass());
	    theHistograms.fill("massa Z2 generati","Massa Z2 generati",90,50,130,z2.mass());
	    theHistograms.fill("massa Z1 generati bene","Massa Z1 generati bene",10,80,100,z1.mass());
	    theHistograms.fill("massa Z2 generati bene","Massa Z2 generati bene",10,70,110,z2.mass());
	    theHistograms.fill("energia Z1 generati","Energia Z1 generati",135,0,2700,z1.e());
	    theHistograms.fill("energia Z2 generati","Energia Z2 generati",135,0,2700,z2.e());
	    theHistograms.fill("energia Z1 generati bene","Energia Z1 generati bene",10,0,1400,z1.e());
	    theHistograms.fill("energia Z2 generati bene","Energia Z2 generati bene",10,0,1400,z2.e());
	    theHistograms.fill("pt Z1 generati","Pt Z1 generati",150,0,900,z1.pt());
	    theHistograms.fill("pt Z2 generati","Pt Z2 generati",150,0,900,z2.pt());
	    theHistograms.fill("pt Z1 generati bene","Pt Z1 generati bene",10,0,350,z1.pt());
	    theHistograms.fill("pt Z2 generati bene","Pt Z2 generati bene",10,0,350,z2.pt());
	    theHistograms.fill("eta Z1 generati bene","Eta Z1 generati bene",10,0,4.5,z1.eta());
	    theHistograms.fill("eta Z2 generati bene","Eta Z2 generati bene",10,0,4.5,z2.eta());
	    theHistograms.fill("E leptone maggiore Z1 buono","Energia leptone piu' energetico Z1",10,0,1000,std::max(z1.daughter(0).e(),z1.daughter(1).e()));
	    theHistograms.fill("E leptone minore Z1 buono","Energia leptone meno energetico Z1",10,0,300,std::min(z1.daughter(0).e(),z1.daughter(1).e()));
	    theHistograms.fill("E leptone maggiore Z2 buono","Energia leptone piu' energetico Z2",10,0,1000,std::max(z2.daughter(0).e(),z2.daughter(1).e()));
	    theHistograms.fill("E leptone minore Z2 buono","Energia leptone meno energetico Z2",10,0,300,std::min(z2.daughter(0).e(),z2.daughter(1).e()));
	    double angoloZZ1=z1.p4().Angle(z2.p4().Vect());
	    theHistograms.fill("angolo bosoni generati","Angolo bosoni generati",35,0,3.5,angoloZZ1);

	    //studio proprietà bosoni ricostruiti	    
	    if(ZZ->first().mass()!=0&&ZZ->passFullSelection()){
   theHistograms.fill("massa Z1 ricostruiti","Massa Z1 ricostruiti",180,50,130,ZZ->first().mass());
   theHistograms.fill("massa Z2 ricostruiti","Massa Z2 ricostruiti",180,50,130,ZZ->second().mass());
   theHistograms.fill("pt Z1 ricostruiti","Pt Z1 ricostruiti",150,0,900,ZZ->first().pt());
   theHistograms.fill("pt Z2 ricostruiti","Pt Z2 ricostruiti",150,0,900,ZZ->second().pt());
   theHistograms.fill("theta bosoni ricostruiti","Theta bosoni ricostruiti",75 ,0,3.5,ZZ->first().p4().Theta());
   theHistograms.fill("theta bosoni ricostruiti","Theta bosoni ricostruiti",75 ,0,3.5,ZZ->second().p4().Theta());
   theHistograms.fill("energia Z1 ricostruiti","Energia Z1 ricostruiti",135,0,2700,ZZ->first().e());
   theHistograms.fill("energia Z2 ricostruiti","Energia Z2 ricostruiti",135,0,2700,ZZ->second().e());
   theHistograms.fill("eta bosoni ricostruiti","Eta bosoni ricostruiti",60,-6,6,ZZ->first().eta());
   theHistograms.fill("eta bosoni ricostruiti","Eta bosoni ricostruiti",60,-6,6,ZZ->second().eta()); 

   double angolo2=ZZ->first().daughter(0).p4().Angle(ZZ->first().daughter(1).p4().Vect());
   theHistograms.fill("angolo leptoni ricostruiti","Angolo leptoni ricostruiti",75,0,3.5,angolo2);
   double angolo3=ZZ->second().daughter(0).p4().Angle(ZZ->second().daughter(1).p4().Vect());
   theHistograms.fill("angolo leptoni ricostruiti","Angolo leptoni ricostruiti",75,0,3.5,angolo3);
   double angoloZZ2=ZZ->first().p4().Angle(ZZ->second().p4().Vect());
   theHistograms.fill("angolo bosoni ricostruiti","Angolo bosoni ricostruiti",35,0,3.5,angoloZZ2);
   
   TLorentzVector b=ZZ->first().p4()+ZZ->second().p4();
   double mass2=sqrt((b.E()*b.E())-(b.Px()*b.Px())-(b.Py()*b.Py())-(b.Pz()*b.Pz()));
   theHistograms.fill("massa dibosoni ricostruiti","Massa ZZ ricostruiti",78,80,600,mass2);
	
   theHistograms.fill("E leptone maggiore ricostruito","Energia leptone piu' energetico ricostruito",200,0,2000,std::max(ZZ->first().daughter(0).e(),ZZ->first().daughter(1).e()));
   theHistograms.fill("E leptone minore ricostruito","Energia leptone meno energetico ricostruito",200,0,800,std::min(ZZ->first().daughter(0).e(),ZZ->first().daughter(1).e()));
   theHistograms.fill("E leptone maggiore ricostruito","Energia leptone piu' energetico ricostruito",200,0,2000,std::max(ZZ->second().daughter(0).e(),ZZ->second().daughter(1).e()));
   theHistograms.fill("E leptone minore ricostruito","Energia leptone meno energetico ricostruito",200,0,800,std::min(ZZ->second().daughter(0).e(),ZZ->second().daughter(1).e()));
   //accoppiamento bosoni generati/ricostruiti
   if(topology.test(0)){
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
       theHistograms.fill("confronto massa Z1","Differenza massa generata/ricostruita Z1",40,-8,8,ZZ->first().mass()-z3.mass());
       theHistograms.fill("confronto energia Z1","Differenza energia generata/ricostruita Z1",200,-40,40,ZZ->first().e()-z3.e());
       theHistograms.fill("confronto pt Z1","Differenza pt generata/ricostruita Z1",150,-30,30,ZZ->first().pt()-z3.pt());
       theHistograms.fill("confronto eta Z1","Differenza eta generata/ricostruita Z1",50,-0.1,0.1,ZZ->first().eta()-z3.eta());
       theHistograms.fill("massa Z1 ricostruiti bene","Massa Z1 ricostruiti bene",10,80,100,z3.mass());
       theHistograms.fill("pt Z1 ricostruiti bene","Pt Z1 ricostruiti bene",10,0,350,z3.pt());
       theHistograms.fill("energia Z1 ricostruiti bene","Energia Z1 ricostruiti bene",10,0,1400,z3.e());
       theHistograms.fill("eta Z1 ricostruiti bene","Eta Z1 ricostruiti bene",10,-4.5,4.5,z3.eta());
       theHistograms.fill("confronto leptone maggiore Z1","DeltaE leptone maggiore Z1",100,-20,20,std::max(ZZ->first().daughter(0).e(),ZZ->first().daughter(1).e())-enlep1);
       theHistograms.fill("confronto leptone minore Z1","DeltaE leptone minore Z1",100,-20,20,std::min(ZZ->first().daughter(0).e(),ZZ->first().daughter(1).e())-enlep2);
       theHistograms.fill("E leptone maggiore Z1 buono ricostruito","Energia leptone piu' energetico Z1",10,0,1000,enlep1);
       theHistograms.fill("E leptone minore Z1 buono ricostruito","Energia leptone meno energetico Z1",10,0,300,enlep2);
       
       theHistograms.fill("confronto massa Z2","Differenza massa generata/ricostruita Z2",40,-8,8,ZZ->second().mass()-z4.mass());
       theHistograms.fill("confronto energia Z2","Differenza energia generata/ricostruita Z2",200,-40,40,ZZ->second().e()-z4.e());
       theHistograms.fill("confronto pt Z2","Differenza pt generata/ricostruita Z2",150,-30,30,ZZ->second().pt()-z4.pt());
       theHistograms.fill("confronto eta Z2","Differenza eta generata/ricostruita Z2",50,-0.1,0.1,ZZ->second().eta()-z4.eta());
       theHistograms.fill("massa Z2 ricostruiti bene","Massa Z2 ricostruiti bene",10,70,110,z4.mass());
       theHistograms.fill("pt Z2 ricostruiti bene","Pt Z2 ricostruiti bene",10,0,350,z4.pt());
       theHistograms.fill("energia Z2 ricostruiti bene","Energia Z2 ricostruiti bene",10,0,1400,z4.e());
       theHistograms.fill("eta Z2 ricostruiti bene","Eta Z2 ricostruiti bene",10,-4.5,4.5,z4.eta());
       theHistograms.fill("confronto leptone maggiore Z2","DeltaE leptone maggiore Z2",100,-20,20,std::max(ZZ->second().daughter(0).e(),ZZ->second().daughter(1).e())-enlep3);
       theHistograms.fill("confronto leptone minore Z2","DeltaE leptone minore Z2",100,-20,20,std::min(ZZ->second().daughter(0).e(),ZZ->second().daughter(1).e())-enlep4);
       theHistograms.fill("E leptone maggiore Z2 buono ricostruito","Energia leptone piu' energetico Z2",10,0,1000,enlep3);
       theHistograms.fill("E leptone minore Z2 buono ricostruito","Energia leptone meno energetico Z2",10,0,300,enlep4); }
   }}}
   //studio bosoni adronici
   int troppi;
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
	 else if(dRq1!=0&&dRq2!=0){troppi+=1;}}
     }}
     if(dRq1!=0&&dRq2!=0&&troppi==0){
       phys::Boson<phys::Particle> jetboson=phys::Boson<phys::Particle>(jet1,jet2,0);
       double massaqq=sqrt((pquark.E()*pquark.E())-(pquark.Px()*pquark.Px())-(pquark.Py()*pquark.Py())-(pquark.Pz()*pquark.Pz()));
       theHistograms.fill("massa qq","Massa qq",120,50,130,massaqq);
       theHistograms.fill("massa bosone adronico","Massa bosone adronico",120,50,130,jetboson.mass());
       theHistograms.fill("deltam qq-jetjet","Differenza di massa qq-jetjet",50,-25,25,jetboson.mass()-massaqq);
       theHistograms.fill("energia bosone adronico","Energia bosone adronico",135,0,2700,jetboson.e());
       theHistograms.fill("eta bosone adronico","Eta bosone adronico",60,-6,6,jetboson.eta());
       theHistograms.fill("pt bosone adronico","Pt bosone adronico",150,0,900,jetboson.pt());
       if(topology.test(0)){
	 theHistograms.fill("massa W","Massa W",40,60,100,massaqq);
	 TLorentzVector q=jetboson.p4()+z1.p4()+z2.p4();
	 TLorentzVector wz1=jetboson.p4()+z1.p4();
	 TLorentzVector wz2=jetboson.p4()+z2.p4();
	 double masswz1=sqrt((wz1.E()*wz1.E())-(wz1.Px()*wz1.Px())-(wz1.Py()*wz1.Py())-(wz1.Pz()*wz1.Pz()));
	 double masswz2=sqrt((wz2.E()*wz2.E())-(wz2.Px()*wz2.Px())-(wz2.Py()*wz2.Py())-(wz2.Pz()*wz2.Pz()));
	 double masstot=sqrt((q.E()*q.E())-(q.Px()*q.Px())-(q.Py()*q.Py())-(q.Pz()*q.Pz()));
	 theHistograms.fill("massa tribosoni","Massa tribosoni",10,200,1000,masstot);
	 theHistograms.fill("massa WZ1","Massa WZ1",20,0,1000,masswz1);
	 theHistograms.fill("massa WZ2","Massa WZ2",20,0,1000,masswz2);
         double angolo4=z1.p4().Angle(jetboson.p4().Vect());
	 double angolo5=z2.p4().Angle(jetboson.p4().Vect());
	 double angolo6=z1.p4().Angle(z2.p4().Vect());
         theHistograms.fill("angolo bosoni WZ1","Angolo bosoni WZ1",10,0,3.5,angolo4);
	 theHistograms.fill("angolo bosoni WZ2","Angolo bosoni WZ2",10,0,3.5,angolo5);
	 theHistograms.fill("angolo bosoni WZ","Angolo bosoni WZ",10,0,3.5,angolo4);
	 theHistograms.fill("angolo bosoni WZ","Angolo bosoni WZ",10,0,3.5,angolo5);
	 theHistograms.fill("angolo bosoni ZZ","Angolo bosoni ZZ",10,0,3.5,angolo6);
	 double Emax=std::max(jetboson.e(),std::max(z1.e(),z2.e()));
         double Emin=std::min(jetboson.e(),std::min(z1.e(),z2.e()));
         theHistograms.fill("energia bosone piu' energetico","Energia bosone piu' energetico",100,100,1800,Emax);
         theHistograms.fill("energia bosone meno energetico","Energia bosone meno energetico",100,0,1000,Emin);
         theHistograms.fill("deltae bosoni piu'/meno","Differenza di energia bosoni piu'/meno energetico",100,0,1800,Emax-Emin);}}
}
