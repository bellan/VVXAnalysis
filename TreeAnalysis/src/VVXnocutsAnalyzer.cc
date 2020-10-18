#include "VVXAnalysis/TreeAnalysis/interface/VVXnocutsAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Colours.h"

#include "VVXAnalysis/Commons/interface/GenVBHelper.h"


#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;
using std::min;
using std::max;
using namespace colour;

using namespace phys;


Int_t VVXnocutsAnalyzer::cut() {
  
  return 1;
}


void VVXnocutsAnalyzer::analyze(){
  
  std::vector<phys::Boson<phys::Particle>>  bosoni=genVBHelper_.ZtoChLep();
  foreach(const phys::Boson<phys::Particle> bos,genVBHelper_.ZtoNeutrinos()){
    bosoni.push_back(bos);}
  foreach(const phys::Boson<phys::Particle> bos,genVBHelper_.ZtoQ()){
    bosoni.push_back(bos);}
  foreach(const phys::Boson<phys::Particle> bos,genVBHelper_.WtoLep()){
    bosoni.push_back(bos);}
  foreach(const phys::Boson<phys::Particle> bos,genVBHelper_.WtoQ()){
    bosoni.push_back(bos);}

  phys::Boson<phys::Particle> Z1,Z2,WZ,nullboson;
  if(bosoni.size()==3){
    foreach(const phys::Boson<phys::Particle> bosone,bosoni){
      if(WZ==nullboson&&abs(bosone.id())==24){
	WZ=bosone;}
      else if(Z1==nullboson&&abs(bosone.id())==23){
        Z1=bosone;}
      else if(Z2==nullboson&&abs(bosone.id())==23){
	Z2=bosone;}
      else{WZ=bosone;}}
    
    TLorentzVector WZ1=Z1.p4()+WZ.p4();
    TLorentzVector WZ2=Z2.p4()+WZ.p4();
    TLorentzVector Z1Z2=Z1.p4()+Z2.p4();
    TLorentzVector WZZ=Z1.p4()+Z2.p4()+WZ.p4();
    double massaWZ1=sqrt((WZ1.E()*WZ1.E())-(WZ1.Px()*WZ1.Px())-(WZ1.Py()*WZ1.Py())-(WZ1.Pz()*WZ1.Pz()));
    double massaWZ2=sqrt((WZ2.E()*WZ2.E())-(WZ2.Px()*WZ2.Px())-(WZ2.Py()*WZ2.Py())-(WZ2.Pz()*WZ2.Pz()));
    double massaZ1Z2=sqrt((Z1Z2.E()*Z1Z2.E())-(Z1Z2.Px()*Z1Z2.Px())-(Z1Z2.Py()*Z1Z2.Py())-(Z1Z2.Pz()*Z1Z2.Pz()));
    double massaWZZ=sqrt((WZZ.E()*WZZ.E())-(WZZ.Px()*WZZ.Px())-(WZZ.Py()*WZZ.Py())-(WZZ.Pz()*WZZ.Pz()));
    theHistograms.fill("massa ZZ","Massa ZZ",100,150,1000,massaZ1Z2);
    if(WZ.id()==23){
      theHistograms.fill("massa ZZ","Massa ZZ",100,150,1000,massaWZ1);
      theHistograms.fill("massa ZZ","Massa ZZ",100,150,1000,massaWZ2);}
    if(WZ.id()==24){
      theHistograms.fill("massa WZ","Massa WZ",100,150,1000,massaWZ1);
      theHistograms.fill("massa WZ","Massa WZ",100,150,1000,massaWZ2);}
    theHistograms.fill("massa tribosoni","Massa tribosoni",100,200,2500,massaWZZ);
    theHistograms.fill("massa bosoni a coppie","Massa bosoni a coppie",100,150,1000,massaWZ1);
    theHistograms.fill("massa bosoni a coppie","Massa bosoni a coppie",100,150,1000,massaWZ2);
    theHistograms.fill("massa bosoni a coppie","Massa bosoni a coppie",100,150,1000,massaZ1Z2);
    theHistograms.fill("energia totale bosoni","Energia totale bosoni",100,0,5000,WZZ.E());
    theHistograms.fill("energia bosoni a coppie","Energia bosoni a coppie",100,0,3000,Z1.e()+Z2.e());
    theHistograms.fill("energia bosoni a coppie","Energia bosoni a coppie",100,0,3000,Z1.e()+WZ.e());
    theHistograms.fill("energia bosoni a coppie","Energia bosoni a coppie",100,0,3000,WZ.e()+Z2.e());
    double angoloWZ1=Z1.p4().Angle(WZ.p4().Vect());
    double angoloWZ2=Z2.p4().Angle(WZ.p4().Vect());
    double angoloZZ=Z1.p4().Angle(Z2.p4().Vect());
    double angolomax=max(angoloWZ1,max(angoloWZ2,angoloZZ));
    double angolomin=min(angoloWZ1,min(angoloWZ2,angoloZZ));
    theHistograms.fill("angolo relativo bosoni","Angolo relativo bosoni",100,0,3.5,angoloWZ1);
    theHistograms.fill("angolo relativo bosoni","Angolo relativo bosoni",100,0,3.5,angoloWZ2);
    theHistograms.fill("angolo relativo bosoni","Angolo relativo bosoni",100,0,3.5,angoloZZ);
    theHistograms.fill("angolo relativo massimo bosoni","Angolo relativo massimo bosoni",100,0,3.5,angolomax);
    theHistograms.fill("angolo relativo minimo bosoni","Angolo relativo minimo bosoni",100,0,3.5,angolomin);
    theHistograms.fill("differenza angolo max/min bosoni","Differenza angolo max/min bosoni",100,0,3.5,angolomax-angolomin);
    theHistograms.fill("rapporto angolo max/min bosoni","Rapporto angolo max/min bosoni",200,0,30,angolomax/angolomin);
    theHistograms.fill("somma angoli relativi","Somma angoli relativi",300,0,10.5,angoloWZ1+angoloWZ2+angoloZZ);
    double emax=max(WZ.e(),max(Z1.e(),Z2.e()));
    double emin=min(WZ.e(),min(Z1.e(),Z2.e()));
    theHistograms.fill("energia bosone piu' energetico","Energia bosone piu' energetico",100,100,1800,emax);
    theHistograms.fill("energia bosone meno energetico","Energia bosone meno energetico",100,0,1000,emin);
    theHistograms.fill("deltae bosoni max/min","Differenza di energia bosoni piu'/meno energetico",100,0,1800,emax-emin);
    theHistograms.fill("rapporto energia bosoni max/min","Rapporto tra energia bosoni piu'/meno energetico",150,0,30,emax/emin);
    double ptmax=max(WZ.pt(),max(Z1.pt(),Z2.pt()));
    double ptmin=min(WZ.pt(),min(Z1.pt(),Z2.pt()));
    theHistograms.fill("pt bosone maggiore","Pt bosone maggiore",100,100,1000,ptmax);
    theHistograms.fill("pt bosone minore","Pt bosone minore",100,0,400,ptmin);
    theHistograms.fill("deltapt bosoni max/min","Differenza di pt bosoni massima/minima",100,0,800,ptmax-ptmin);
    theHistograms.fill("rapporto pt bosoni max/min","Rapporto tra pt bosoni massima/minima",150,0,30,ptmax/ptmin);
    theHistograms.fill("somma totale pt(in modulo)","Somma totale pt (in modulo)",100,0,2000,WZ.pt()+Z1.pt()+Z2.pt());
    theHistograms.fill("rapporto e max/pt max","Rapporto tra e massima e pt minima",190,1,20,emax/ptmax);
    theHistograms.fill("rapporto e min/pt min","Rapporto tra e massima e pt minima",190,1,20,emin/ptmin);
    theHistograms.fill("rapporto e max/pt min","Rapporto tra e massima e pt minima",200,0,80,emax/ptmin);
    theHistograms.fill("rapporto e min/pt max","Rapporto tra e minima e pt massima",100,0,7,emin/ptmax);
    double etamax=max(abs(WZ.eta()),max(abs(Z1.eta()),abs(Z2.eta())));
    double etamin=min(abs(WZ.eta()),min(abs(Z1.eta()),abs(Z2.eta())));
    theHistograms.fill("eta massima","Eta massima (in modulo) bosoni",100,0,7,etamax);
    theHistograms.fill("eta minima","Eta minima (in modulo) bosoni",100,0,4,etamin);
    theHistograms.fill("deltaeta max/min","Differenza di eta bosoni max/min (in modulo)",100,0,7,etamax-etamin);
    theHistograms.fill("rapporto eta max/min","Rapporto eta max/min (in modulo)",100,1,20,etamax/etamin);
  }}
