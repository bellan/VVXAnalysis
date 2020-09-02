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

double a,b,c,d,mass,a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4,d1,d2,d3,d4,eta1,eta2,eta3,eta4,phi1,phi2,phi3,phi4,mass1,mass2;
int good;

void VZZaQGCAnalyzer::analyze(){
   foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
     if((genVBParticle.daughter(0).id()==11&&genVBParticle.daughter(1).id()==-11)||(genVBParticle.daughter(0).id()==13&&genVBParticle.daughter(1).id()==-13)){
      theHistograms.fill("massa bosoni generati","Massa bosoni generati",180,50,130,genVBParticle.mass());
      theHistograms.fill("pt bosoni generati","Pt bosoni generati",150,0,900,genVBParticle.pt());
      theHistograms.fill("theta bosoni generati","Theta bosoni generati",75 ,0,3.5,genVBParticle.p4().Theta());
      theHistograms.fill("energia bosoni generati","Energia bosoni generati",135,0,2700,genVBParticle.e());
      if(genVBParticle.daughter(1).e()>genVBParticle.daughter(0).e()){
	theHistograms.fill("E leptone maggiore","Energia leptone piu' energetico",200,0,2000,genVBParticle.daughter(1).e());}
      else{theHistograms.fill("E leptone maggiore","Energia leptone piu' energetico",200,0,2000,genVBParticle.daughter(0).e());};
      if(genVBParticle.daughter(1).e()>genVBParticle.daughter(0).e()){
	theHistograms.fill("E leptone minore","Energia leptone meno energetico",200,0,800,genVBParticle.daughter(0).e());}
      else{theHistograms.fill("E leptone minore","Energia leptone meno energetico",200,0,800,genVBParticle.daughter(1).e());};
     }
   }
   theHistograms.fill("numero bosoni","Numero bosoni",5,0,5,genVBParticles->size());
   
   if(genVBParticles->size()==2){
     good=0;
     foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
	  if((genVBParticle.daughter(0).id()==11&&genVBParticle.daughter(1).id()==-11)||(genVBParticle.daughter(0).id()==13&&genVBParticle.daughter(1).id()==-13)){
	    good+=1;}
	  if(good==2){
     a=b=c=d=mass=0;
      foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
	a+=genVBParticle.p4().Px();
	b+=genVBParticle.p4().Py();
        c+=genVBParticle.p4().Pz();
	d+=genVBParticle.p4().E();}
      mass=sqrt((d*d)-(a*a)-(b*b)-(c*c));
      theHistograms.fill("massa dibosoni generati","Massa dibosoni generati",55,80,520,mass);}
     }}
   
   if(genVBParticles->size()==4){
      foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
	  if((genVBParticle.daughter(0).id()==11&&genVBParticle.daughter(1).id()==-11)||(genVBParticle.daughter(0).id()==13&&genVBParticle.daughter(1).id()==-13)){
	    good+=1;}
	  if(good==4){
     a1=a2=a3=a4=b1=b2=b3=b4=c1=c2=c3=c4=d1=d2=d3=d4=eta1=eta2=eta3=eta4=phi1=phi2=phi3=phi4=mass1=mass2=0;
     foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
       if(a1==0){
	 a1=genVBParticle.p4().Px();
	 b1=genVBParticle.p4().Py();
         c1=genVBParticle.p4().Pz();
	 d1=genVBParticle.p4().E();
	 eta1=genVBParticle.p4().Eta();
	 phi1=genVBParticle.p4().Phi();}
       else{if(a2==0){
	 a2=genVBParticle.p4().Px();
	 b2=genVBParticle.p4().Py();
         c2=genVBParticle.p4().Pz();
	 d2=genVBParticle.p4().E();
	 eta2=genVBParticle.p4().Eta();
	 phi2=genVBParticle.p4().Phi();
	 }
	 else{if(a3==0){
	   if(physmath::deltaPhi(genVBParticle.p4().Phi(),phi1)<physmath::deltaPhi(phi2,phi1)){
	     a3=a2;
	     b3=b2;
	     c3=c2;
	     d3=d2;
	     eta3=eta2;
	     phi3=phi2;
	     a2=genVBParticle.p4().Px();
	     b2=genVBParticle.p4().Py();
             c2=genVBParticle.p4().Pz();
	     d2=genVBParticle.p4().E();
	     eta2=genVBParticle.p4().Eta();
	     phi2=genVBParticle.p4().Phi();}
	   else{a3=genVBParticle.p4().Px();
	        b3=genVBParticle.p4().Py();
                c3=genVBParticle.p4().Pz();
	        d3=genVBParticle.p4().E();
	        eta3=genVBParticle.p4().Eta();
	        phi3=genVBParticle.p4().Phi();}
	 }
	   else{if(physmath::deltaPhi(genVBParticle.p4().Phi(),phi1)<physmath::deltaPhi(phi2,phi1)){
	       a4=a2;
	       b4=b2;
	       c4=c2;
	       d4=d2;
	       eta4=eta2;
	       phi4=phi2;
	       a2=genVBParticle.p4().Px();
	       b2=genVBParticle.p4().Py();
               c2=genVBParticle.p4().Pz();
	       d2=genVBParticle.p4().E();
	       eta2=genVBParticle.p4().Eta();
	       phi2=genVBParticle.p4().Phi();}
	     else{a4=genVBParticle.p4().Px();
	       b4=genVBParticle.p4().Py();
               c4=genVBParticle.p4().Pz();
	       d4=genVBParticle.p4().E();
	       eta4=genVBParticle.p4().Eta();
	       phi4=genVBParticle.p4().Phi();}
	   }
	 }
       }
     }
       mass1=sqrt((d1+d2)*(d1+d2)-(a1+a2)*(a1+a2)-(b1+b2)*(b1+b2)-(c1+c2)*(c1+c2));
       mass2=sqrt((d3+d4)*(d3+d4)-(a3+a4)*(a3+a4)-(b3+b4)*(b3+b4)-(c3+c4)*(c3+c4));
       theHistograms.fill("massa dibosoni generati","Massa dibosoni generati",55,80,520,mass1);
       theHistograms.fill("massa dibosoni generati","Massa dibosoni generati",55,80,520,mass2);}}}
   theHistograms.fill("massa bosoni ricostruiti","Massa bosoni ricostruiti",180,50,130,ZZ->first().mass());
   theHistograms.fill("massa bosoni ricostruiti","Massa bosoni ricostruiti",180,50,130,ZZ->second().mass());
   theHistograms.fill("pt bosoni ricostruiti","Pt bosoni ricostruiti",150,0,900,ZZ->first().pt());
   theHistograms.fill("pt bosoni ricostruiti","Pt bosoni ricostruiti",150,0,900,ZZ->second().pt());
   theHistograms.fill("theta bosoni ricostruiti","Theta bosoni ricostruiti",75 ,0,3.5,ZZ->first().eta());
   theHistograms.fill("theta bosoni ricostruiti","Theta bosoni ricostruiti",75 ,0,3.5,ZZ->second().eta());
   theHistograms.fill("energia bosoni ricostruiti","Energia bosoni ricostruiti",135,0,2700,ZZ->first().e());
   theHistograms.fill("energia bosoni ricostruiti","Energia bosoni ricostruiti",135,0,2700,ZZ->second().e());
   if(ZZ->first().daughter(1).e()>ZZ->first().daughter(0).e()){
     theHistograms.fill("E leptone maggiore ricostruito","Energia leptone piu' energetico ricostruito",200,0,2000,ZZ->first().daughter(1).e());}
   else{theHistograms.fill("E leptone maggiore ricostruito","Energia leptone piu' energetico ricostruito",200,0,2000,ZZ->first().daughter(0).e());};
   
   if(ZZ->second().daughter(1).e()>ZZ->second().daughter(0).e()){
     theHistograms.fill("E leptone maggiore ricostruito","Energia leptone piu' energetico ricostruito",200,0,2000,ZZ->second().daughter(1).e());}
   else{theHistograms.fill("E leptone maggiore ricostruito","Energia leptone piu' energetico ricostruito",200,0,2000,ZZ->second().daughter(0).e());};

   
   if(ZZ->first().daughter(1).e()>ZZ->first().daughter(0).e()){
     theHistograms.fill("E leptone minore ricostruito","Energia leptone meno energetico ricostruito",200,0,2000,ZZ->first().daughter(0).e());}
   else{theHistograms.fill("E leptone minore ricostruito","Energia leptone meno energetico ricostruito",200,0,2000,ZZ->first().daughter(1).e());};
   if(ZZ->second().daughter(1).e()>ZZ->second().daughter(0).e()){
     theHistograms.fill("E leptone minore ricostruito","Energia leptone meno energetico ricostruito",200,0,2000,ZZ->second().daughter(0).e());}
   else{theHistograms.fill("E leptone minore ricostruito","Energia leptone meno energetico ricostruito",200,0,2000,ZZ->second().daughter(1).e());};
   
   
   
   } 
