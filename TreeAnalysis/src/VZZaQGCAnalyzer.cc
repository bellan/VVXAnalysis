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

double theta1,theta2,theta3,theta4,theta5,theta6,phi1,phi2,phi3,phi4,phi5,phi6,angolo1,angolo2,angolo,dR,dR2,dRq,dRq1,dRq2,enlep1,enlep2,enlep3,enlep4,mass1,mass2,massaqq,massajetjet;
TLorentzVector a,b,pquark,pjet,zero;
int nquark,troppi;
phys::Particle jet;
phys::Boson<phys::Particle> z1,z2,z3,z4,bosonzero;

void VZZaQGCAnalyzer::analyze(){
   a=zero;
   foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
     if((abs(genVBParticle.daughter(0).id())==11||abs(genVBParticle.daughter(0).id())==13)&&genVBParticle.id()==23){
      theHistograms.fill("theta bosoni generati","Theta bosoni generati",75 ,0,3.5,genVBParticle.p4().Theta());
      theHistograms.fill("eta bosoni generati","Eta bosoni generati",60,0,6,genVBParticle.eta());
      theta1=genVBParticle.daughter(0).p4().Theta();
      theta2=genVBParticle.daughter(1).p4().Theta();
      phi1=genVBParticle.daughter(0).phi();
      phi2=genVBParticle.daughter(1).phi();
      angolo=acos(sin(theta1)*sin(theta2)*cos(physmath::deltaPhi(phi1,phi2))+cos(theta1)*cos(theta2));
      theHistograms.fill("angolo leptoni generati","Angolo leptoni generati",75,0,3.5,angolo);
      if(genVBParticle.daughter(1).e()>genVBParticle.daughter(0).e()){
	theHistograms.fill("E leptone maggiore","Energia leptone piu' energetico",200,0,2000,genVBParticle.daughter(1).e());
      	theHistograms.fill("E leptone minore","Energia leptone meno energetico",200,0,800,genVBParticle.daughter(0).e());}
      else{theHistograms.fill("E leptone maggiore","Energia leptone piu' energetico",200,0,2000,genVBParticle.daughter(0).e());
	theHistograms.fill("E leptone minore","Energia leptone meno energetico",200,0,800,genVBParticle.daughter(1).e());}
      a+=genVBParticle.p4();
     }}
   z1=z2=bosonzero;
   if(topology.test(0)){
     mass1=sqrt((a.E()*a.E())-(a.Px()*a.Px())-(a.Py()*a.Py())-(a.Pz()*a.Pz()));
       theHistograms.fill("massa dibosoni generati","Massa ZZ generati",78,80,600,mass1);
       
	 foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
	   if(z1.mass()==0){
	     z1=genVBParticle;}
	   else{if(abs(genVBParticle.mass()-phys::ZMASS)<abs(z1.mass()-phys::ZMASS)){
	       z2=z1;}
		else{z2=genVBParticle;}}}
	 
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
	    if(z1.daughter(0).e()>z1.daughter(1).e()){
	      theHistograms.fill("E leptone maggiore Z1 buono","Energia leptone piu' energetico Z1",10,0,1000,z1.daughter(0).e());
	      theHistograms.fill("E leptone minore Z1 buono","Energia leptone meno energetico Z1",10,0,300,z1.daughter(1).e());}
	    else{theHistograms.fill("E leptone maggiore Z1 buono","Energia leptone piu' energetico Z1",10,0,1000,z1.daughter(1).e());
	      theHistograms.fill("E leptone minore Z1 buono","Energia leptone meno energetico Z1",10,0,300,z1.daughter(0).e());}
            if(z2.daughter(0).e()>z2.daughter(1).e()){
	      theHistograms.fill("E leptone maggiore Z2 buono","Energia leptone piu' energetico Z2",10,0,1000,z2.daughter(0).e());
	      theHistograms.fill("E leptone minore Z2 buono","Energia leptone meno energetico Z2",10,0,300,z2.daughter(1).e());}
	    else{theHistograms.fill("E leptone maggiore Z2 buono","Energia leptone piu' energetico Z2",10,0,1000,z2.daughter(1).e());
	      theHistograms.fill("E leptone minore Z2 buono","Energia leptone meno energetico Z2",10,0,300,z2.daughter(0).e());}}
     
   if(ZZ->first().mass()!=0){
   theHistograms.fill("massa Z1 ricostruiti","Massa Z1 ricostruiti",180,50,130,ZZ->first().mass());
   theHistograms.fill("massa Z2 ricostruiti","Massa Z2 ricostruiti",180,50,130,ZZ->second().mass());
   theHistograms.fill("pt Z1 ricostruiti","Pt Z1 ricostruiti",150,0,900,ZZ->first().pt());
   theHistograms.fill("pt Z2 ricostruiti","Pt Z2 ricostruiti",150,0,900,ZZ->second().pt());
   theHistograms.fill("theta bosoni ricostruiti","Theta bosoni ricostruiti",75 ,0,3.5,ZZ->first().p4().Theta());
   theHistograms.fill("theta bosoni ricostruiti","Theta bosoni ricostruiti",75 ,0,3.5,ZZ->second().p4().Theta());
   theHistograms.fill("energia Z1 ricostruiti","Energia Z1 ricostruiti",135,0,2700,ZZ->first().e());
   theHistograms.fill("energia Z2 ricostruiti","Energia Z2 ricostruiti",135,0,2700,ZZ->second().e());
   theHistograms.fill("eta bosoni ricostruiti","Eta bosoni ricostruiti",60,0,6,ZZ->first().eta());
   theHistograms.fill("eta bosoni ricostruiti","Eta bosoni ricostruiti",60,0,6,ZZ->second().eta()); 

   theta3=ZZ->first().daughter(0).p4().Theta();
   theta4=ZZ->first().daughter(1).p4().Theta();
   phi3=ZZ->first().daughter(0).phi();
   phi4=ZZ->first().daughter(1).phi();
   angolo1=acos(sin(theta3)*sin(theta4)*cos(physmath::deltaPhi(phi3,phi4))+cos(theta3)*cos(theta4));
   theHistograms.fill("angolo leptoni ricostruiti","Angolo leptoni ricostruiti",75,0,3.5,angolo1);
   
   theta5=ZZ->second().daughter(0).p4().Theta();
   theta6=ZZ->second().daughter(1).p4().Theta();
   phi5=ZZ->second().daughter(0).phi();
   phi6=ZZ->second().daughter(1).phi();
   angolo2=acos(sin(theta5)*sin(theta6)*cos(physmath::deltaPhi(phi5,phi6))+cos(theta5)*cos(theta6));
   theHistograms.fill("angolo leptoni ricostruiti","Angolo leptoni ricostruiti",75,0,3.5,angolo2);
   
   b=ZZ->first().p4()+ZZ->second().p4();
   mass2=sqrt((b.E()*b.E())-(b.Px()*b.Px())-(b.Py()*b.Py())-(b.Pz()*b.Pz()));
   theHistograms.fill("massa dibosoni ricostruiti","Massa ZZ ricostruiti",78,80,600,mass2);
	
   if(ZZ->first().daughter(1).e()>ZZ->first().daughter(0).e()){
     theHistograms.fill("E leptone maggiore ricostruito","Energia leptone piu' energetico ricostruito",200,0,2000,ZZ->first().daughter(1).e());}
   else{theHistograms.fill("E leptone maggiore ricostruito","Energia leptone piu' energetico ricostruito",200,0,2000,ZZ->first().daughter(0).e());};
   
   if(ZZ->second().daughter(1).e()>ZZ->second().daughter(0).e()){
     theHistograms.fill("E leptone maggiore ricostruito","Energia leptone piu' energetico ricostruito",200,0,2000,ZZ->second().daughter(1).e());}
   else{theHistograms.fill("E leptone maggiore ricostruito","Energia leptone piu' energetico ricostruito",200,0,2000,ZZ->second().daughter(0).e());};
   
   if(ZZ->first().daughter(1).e()>ZZ->first().daughter(0).e()){
     theHistograms.fill("E leptone minore ricostruito","Energia leptone meno energetico ricostruito",200,0,800,ZZ->first().daughter(0).e());}
   else{theHistograms.fill("E leptone minore ricostruito","Energia leptone meno energetico ricostruito",200,0,800,ZZ->first().daughter(1).e());};
   
   if(ZZ->second().daughter(1).e()>ZZ->second().daughter(0).e()){
     theHistograms.fill("E leptone minore ricostruito","Energia leptone meno energetico ricostruito",200,0,800,ZZ->second().daughter(0).e());}
   else{theHistograms.fill("E leptone minore ricostruito","Energia leptone meno energetico ricostruito",200,0,800,ZZ->second().daughter(1).e());};
   
   dR=dR2=9999;
   enlep1=enlep2=enlep3=enlep4=0;
   if(topology.test(0)){
    foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
      if((abs(genVBParticle.daughter(0).id())==11||abs(genVBParticle.daughter(0).id())==13)&&genVBParticle.id()==23){
       if(physmath::deltaR(genVBParticle,ZZ->first())<dR){
	   dR=physmath::deltaR(genVBParticle,ZZ->first());
	   z3=genVBParticle;
	   if(genVBParticle.daughter(0).e()>genVBParticle.daughter(1).e()){
	     enlep1=genVBParticle.daughter(0).e();
	     enlep2=genVBParticle.daughter(1).e();}
	   else{enlep3=genVBParticle.daughter(1).e();
	     enlep4=genVBParticle.daughter(0).e();}}}}}
   theHistograms.fill("dR Z1","dR Z1",100,0,0.5,dR);
       foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
	 if((abs(genVBParticle.daughter(0).id())==11||abs(genVBParticle.daughter(0).id())==13)&&genVBParticle.id()==23){
           if(physmath::deltaR(genVBParticle,ZZ->second())<dR2){
	   dR2=physmath::deltaR(genVBParticle,ZZ->second());
	   z4=genVBParticle;
	   if(genVBParticle.daughter(0).e()>genVBParticle.daughter(1).e()){
	   enlep3=genVBParticle.daughter(0).e();
	   enlep4=genVBParticle.daughter(1).e();}
         else{enlep3=genVBParticle.daughter(1).e();
	   enlep4=genVBParticle.daughter(0).e();}}}}
       theHistograms.fill("dR Z2","dR Z2",100,0,0.5,dR2);

     if(dR<0.05&&dR2<0.1&&z3!=z4){
       theHistograms.fill("confronto massa Z1","Differenza massa generata/ricostruita Z1",40,-8,8,ZZ->first().mass()-z3.mass());
       theHistograms.fill("confronto energia Z1","Differenza energia generata/ricostruita Z1",200,-40,40,ZZ->first().e()-z3.e());
       theHistograms.fill("confronto pt Z1","Differenza pt generata/ricostruita Z1",150,-30,30,ZZ->first().pt()-z3.pt());
       theHistograms.fill("confronto eta Z1","Differenza eta generata/ricostruita Z1",50,-0.1,0.1,ZZ->first().eta()-z3.eta());
       theHistograms.fill("massa Z1 ricostruiti bene","Massa Z1 ricostruiti bene",10,80,100,z3.mass());
       theHistograms.fill("pt Z1 ricostruiti bene","Pt Z1 ricostruiti bene",10,0,350,z3.pt());
       theHistograms.fill("energia Z1 ricostruiti bene","Energia Z1 ricostruiti bene",10,0,1400,z3.e());
       theHistograms.fill("eta Z1 ricostruiti bene","Eta Z1 ricostruiti bene",10,0,4.5,z3.eta());
     if(ZZ->first().daughter(0).e()>ZZ->first().daughter(1).e()){
       theHistograms.fill("confronto leptone maggiore Z1","DeltaE leptone maggiore Z1",100,-20,20,ZZ->first().daughter(0).e()-enlep1);
       theHistograms.fill("confronto leptone minore Z1","DeltaE leptone minore Z1",100,-20,20,ZZ->first().daughter(1).e()-enlep2);
       theHistograms.fill("E leptone maggiore Z1 buono ricostruito","Energia leptone piu' energetico Z1",10,0,1000,enlep1);
       theHistograms.fill("E leptone minore Z1 buono ricostruito","Energia leptone meno energetico Z1",10,0,300,enlep2);}
     else{theHistograms.fill("confronto leptone maggiore Z1","DeltaE leptone maggiore Z1",100,-20,20,ZZ->first().daughter(1).e()-enlep1);
       theHistograms.fill("confronto leptone minore Z1","DeltaE leptone minore Z1",100,-20,20,ZZ->first().daughter(0).e()-enlep2);
       theHistograms.fill("E leptone maggiore Z1 buono ricostruito","Energia leptone piu' energetico Z1",10,0,1000,enlep1);
       theHistograms.fill("E leptone minore Z1 buono ricostruito","Energia leptone meno energetico Z1",10,0,300,enlep2);}
     theHistograms.fill("confronto massa Z2","Differenza massa generata/ricostruita Z2",40,-8,8,ZZ->second().mass()-z4.mass());
     theHistograms.fill("confronto energia Z2","Differenza energia generata/ricostruita Z2",200,-40,40,ZZ->second().e()-z4.e());
     theHistograms.fill("confronto pt Z2","Differenza pt generata/ricostruita Z2",150,-30,30,ZZ->second().pt()-z4.pt());
     theHistograms.fill("confronto eta Z2","Differenza eta generata/ricostruita Z2",50,-0.1,0.1,ZZ->second().eta()-z4.eta());
     theHistograms.fill("massa Z2 ricostruiti bene","Massa Z2 ricostruiti bene",10,70,110,z4.mass());
     theHistograms.fill("pt Z2 ricostruiti bene","Pt Z2 ricostruiti bene",10,0,350,z4.pt());
     theHistograms.fill("energia Z2 ricostruiti bene","Energia Z2 ricostruiti bene",10,0,1400,z4.e());
     theHistograms.fill("eta Z2 ricostruiti bene","Eta Z2 ricostruiti bene",10,0,4.5,z4.eta());
       if(ZZ->second().daughter(0).e()>ZZ->second().daughter(1).e()){
       theHistograms.fill("confronto leptone maggiore Z2","DeltaE leptone maggiore Z2",100,-20,20,ZZ->second().daughter(0).e()-enlep3);
       theHistograms.fill("confronto leptone minore Z2","DeltaE leptone minore Z2",100,-20,20,ZZ->second().daughter(1).e()-enlep4);
       theHistograms.fill("E leptone maggiore Z2 buono ricostruito","Energia leptone piu' energetico Z2",10,0,1000,enlep3);
       theHistograms.fill("E leptone minore Z2 buono ricostruito","Energia leptone meno energetico Z2",10,0,300,enlep4);}
     else{theHistograms.fill("confronto leptone maggiore Z2","DeltaE leptone maggiore Z2",100,-20,20,ZZ->second().daughter(1).e()-enlep3);
       theHistograms.fill("confronto leptone minore Z2","DeltaE leptone minore Z2",100,-20,20,ZZ->second().daughter(0).e()-enlep4);
       theHistograms.fill("E leptone maggiore Z2 buono ricostruito","Energia leptone piu' energetico Z2",10,0,1000,ZZ->second().daughter(1).e());
       theHistograms.fill("E leptone minore Z2 buono ricostruito","Energia leptone meno energetico Z2",10,0,300,ZZ->second().daughter(0).e());} }
   }
   troppi=0;
   pquark=pjet=zero;
   dRq1=dRq2=0;
   foreach(const phys::Particle genParticle,*genParticles){
     dRq=9999;
     if(abs(genParticle.id())>0&&abs(genParticle.id())<7){
       foreach(const phys::Particle genJet,*genJets){
	 if(physmath::deltaR(genJet,genParticle)<dRq){
	   dRq=physmath::deltaR(genJet,genParticle);
	   jet=genJet;}}
       theHistograms.fill("dR quark-jet","dR quark-jet",100,0,0.5,dRq);
       if(dRq<0.1){
	 theHistograms.fill("deltam quark-jet","Differenza di massa quark-jet",100,-5,30,jet.mass()-genParticle.mass());
	 if(dRq1==0){
	   dRq1=dRq;
	   pquark+=genParticle.p4();
	   pjet+=jet.p4();}
	 else if(dRq2==0&&jet.e()!=pjet.E()){
	   dRq2=dRq;
	   pquark+=genParticle.p4();
	   pjet+=jet.p4();}
	 else if(dRq1!=0&&dRq2!=0){troppi+=1;}}
     }}
     if(dRq1!=0&&dRq2!=0&&troppi==0){
       massaqq=sqrt((pquark.E()*pquark.E())-(pquark.Px()*pquark.Px())-(pquark.Py()*pquark.Py())-(pquark.Pz()*pquark.Pz()));
       massajetjet=sqrt((pjet.E()*pjet.E())-(pjet.Px()*pjet.Px())-(pjet.Py()*pjet.Py())-(pjet.Pz()*pjet.Pz()));
       theHistograms.fill("massa qq","Massa qq",180,50,130,massaqq);
       theHistograms.fill("massa jet-jet","Massa jet-jet",180,50,130,massajetjet);
       theHistograms.fill("deltam qq-jetjet","Differenza di massa qq-jetjet",50,-25,25,massajetjet-massaqq);}
}
