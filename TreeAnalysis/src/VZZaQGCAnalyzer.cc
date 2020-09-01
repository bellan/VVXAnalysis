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
}
 
