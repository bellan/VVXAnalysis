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

double a1,b1,c1,d1,mass1,a2,b2,c2,d2,a3,b3,c3,d3,mass2,theta1,theta2,theta3,theta4,theta5,theta6,phi1,phi2,phi3,phi4,phi5,phi6,angolo1,angolo2,angolo,dR;
int good;

void VZZaQGCAnalyzer::analyze(){
   foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
     if((genVBParticle.daughter(0).id()==11&&genVBParticle.daughter(1).id()==-11)||(genVBParticle.daughter(0).id()==13&&genVBParticle.daughter(1).id()==-13)){
      theHistograms.fill("massa bosoni generati","Massa bosoni generati",180,50,130,genVBParticle.mass());
      theHistograms.fill("pt bosoni generati","Pt bosoni generati",150,0,900,genVBParticle.pt());
      theHistograms.fill("theta bosoni generati","Theta bosoni generati",75 ,0,3.5,genVBParticle.p4().Theta());
      theHistograms.fill("energia bosoni generati","Energia bosoni generati",135,0,2700,genVBParticle.e());
      theHistograms.fill("eta bosoni generati","Eta bosoni generati",60,0,6,genVBParticle.eta());
      theta1=genVBParticle.daughter(0).p4().Theta();
      theta2=genVBParticle.daughter(1).p4().Theta();
      phi1=genVBParticle.daughter(0).phi();
      phi2=genVBParticle.daughter(1).phi();
      angolo=acos(sin(theta1)*sin(theta2)*cos(phi1-phi2)+cos(theta1)*cos(theta2));
      theHistograms.fill("angolo leptoni generati","Angolo leptoni generati",75,0,3.5,angolo);
      if(genVBParticle.daughter(1).e()>genVBParticle.daughter(0).e()){
	theHistograms.fill("E leptone maggiore","Energia leptone piu' energetico",200,0,2000,genVBParticle.daughter(1).e());}
      else{theHistograms.fill("E leptone maggiore","Energia leptone piu' energetico",200,0,2000,genVBParticle.daughter(0).e());};
      if(genVBParticle.daughter(1).e()>genVBParticle.daughter(0).e()){
	theHistograms.fill("E leptone minore","Energia leptone meno energetico",200,0,800,genVBParticle.daughter(0).e());}
      else{theHistograms.fill("E leptone minore","Energia leptone meno energetico",200,0,800,genVBParticle.daughter(1).e());};
     }
   }  
     good=0;
     a1=b1=c1=d1=mass1=0;
     foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
	  if((genVBParticle.daughter(0).id()==11&&genVBParticle.daughter(1).id()==-11)||(genVBParticle.daughter(0).id()==13&&genVBParticle.daughter(1).id()==-13)){
	a1+=genVBParticle.p4().Px();
	b1+=genVBParticle.p4().Py();
        c1+=genVBParticle.p4().Pz();
	d1+=genVBParticle.p4().E();
	good+=1;
	  }}
	  mass1=sqrt((d1*d1)-(a1*a1)-(b1*b1)-(c1*c1));
	  if(good==2){
	    theHistograms.fill("massa dibosoni generati","Massa dibosoni generati",78,80,600,mass1);}
   if(ZZ->first().mass()!=0){
   theHistograms.fill("massa bosoni ricostruiti","Massa bosoni ricostruiti",180,50,130,ZZ->first().mass());
   theHistograms.fill("massa bosoni ricostruiti","Massa bosoni ricostruiti",180,50,130,ZZ->second().mass());
   theHistograms.fill("pt bosoni ricostruiti","Pt bosoni ricostruiti",150,0,900,ZZ->first().pt());
   theHistograms.fill("pt bosoni ricostruiti","Pt bosoni ricostruiti",150,0,900,ZZ->second().pt());
   theHistograms.fill("theta bosoni ricostruiti","Theta bosoni ricostruiti",75 ,0,3.5,ZZ->first().p4().Theta());
   theHistograms.fill("theta bosoni ricostruiti","Theta bosoni ricostruiti",75 ,0,3.5,ZZ->second().p4().Theta());
   theHistograms.fill("energia bosoni ricostruiti","Energia bosoni ricostruiti",135,0,2700,ZZ->first().e());
   theHistograms.fill("energia bosoni ricostruiti","Energia bosoni ricostruiti",135,0,2700,ZZ->second().e());
   theHistograms.fill("eta bosoni ricostruiti","Eta bosoni ricostruiti",60,0,6,ZZ->first().eta());
   theHistograms.fill("eta bosoni ricostruiti","Eta bosoni ricostruiti",60,0,6,ZZ->second().eta()); 

   theta3=ZZ->first().daughter(0).p4().Theta();
   theta4=ZZ->first().daughter(1).p4().Theta();
   phi3=ZZ->first().daughter(0).phi();
   phi4=ZZ->first().daughter(1).phi();
   angolo1=acos(sin(theta3)*sin(theta4)*cos(phi3-phi4)+cos(theta3)*cos(theta4));
   theHistograms.fill("angolo leptoni ricostruiti","Angolo leptoni ricostruiti",75,0,3.5,angolo1);
   
   theta5=ZZ->second().daughter(0).p4().Theta();
   theta6=ZZ->second().daughter(1).p4().Theta();
   phi5=ZZ->second().daughter(0).phi();
   phi6=ZZ->second().daughter(1).phi();
   angolo2=acos(sin(theta5)*sin(theta6)*cos(phi5-phi6)+cos(theta5)*cos(theta6));
   theHistograms.fill("angolo leptoni ricostruiti","Angolo leptoni ricostruiti",75,0,3.5,angolo2);
   
      a2=b2=c2=d2=a3=b3=c3=d3=mass2=0;
      a2=ZZ->first().p4().Px();
      b2=ZZ->first().p4().Py();
      c2=ZZ->first().p4().Pz();
      d2=ZZ->first().p4().E();
      a3=ZZ->second().p4().Px();
      b3=ZZ->second().p4().Py();
      c3=ZZ->second().p4().Pz();
      d3=ZZ->second().p4().E();
      mass2=sqrt((d2+d3)*(d2+d3)-(a2+a3)*(a2+a3)-(b2+b3)*(b2+b3)-(c2+c3)*(c2+c3));
      theHistograms.fill("massa dibosoni ricostruiti","Massa dibosoni ricostruiti",78,80,600,mass2);
	
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
   
   foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
     if((genVBParticle.daughter(0).id()==11&&genVBParticle.daughter(1).id()==-11)||(genVBParticle.daughter(0).id()==13&&genVBParticle.daughter(1).id()==-13)){
       dR=physmath::deltaR(genVBParticle,ZZ->first());
       if(abs(dR)<0.05){
	 theHistograms.fill("delta M","Delta m bosoni generati/ricostruiti",40,-10,10,ZZ->first().mass()-genVBParticle.mass());}
       dR=physmath::deltaR(genVBParticle,ZZ->second());
       if(abs(dR)<0.05){
	 theHistograms.fill("delta M","Delta m bosoni generati/ricostruiti",40,-10,10,ZZ->second().mass()-genVBParticle.mass());}}}
	 }
} 
