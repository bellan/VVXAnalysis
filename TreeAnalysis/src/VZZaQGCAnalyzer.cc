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
     a=0;
     b=0;
     c=0;
     d=0;
     mass=0;
      foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
	a+=genVBParticle.p4().Px();
	b+=genVBParticle.p4().Py();
        c+=genVBParticle.p4().Pz();
	d+=genVBParticle.p4().E();}
      mass=sqrt((d*d)-(a*a)-(b*b)-(c*c));
      theHistograms.fill("massa dibosoni generati","Massa dibosoni generati",90,80,440,mass);
   }
  
   }

     

 
