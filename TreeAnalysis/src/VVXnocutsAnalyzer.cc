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

  gErrorIgnoreLevel=kFatal; //to avoid errors due to nonexistent tree branches
  
  std::vector<phys::Boson<phys::Particle>>  bosons=genVBHelper_.ZtoChLep();
  foreach(const phys::Boson<phys::Particle> bos,genVBHelper_.ZtoNeutrinos()){
    bosons.push_back(bos);}
  foreach(const phys::Boson<phys::Particle> bos,genVBHelper_.ZtoQ()){
    bosons.push_back(bos);}
  foreach(const phys::Boson<phys::Particle> bos,genVBHelper_.WtoLep()){
    bosons.push_back(bos);}
  foreach(const phys::Boson<phys::Particle> bos,genVBHelper_.WtoQ()){
    bosons.push_back(bos);}

  phys::Boson<phys::Particle> Z1,Z2,WZ,nullboson;
  if(bosons.size()==3){
    foreach(const phys::Boson<phys::Particle> bosone,bosons){
      if(WZ==nullboson&&abs(bosone.id())==24){
	WZ=bosone;}
      else if(Z1==nullboson&&abs(bosone.id())==23){
        Z1=bosone;}
      else if(Z2==nullboson&&abs(bosone.id())==23){
	Z2=bosone;}
      else{WZ=bosone;}}


    // analysis of three generic bosons
    
    TLorentzVector WZ1=Z1.p4()+WZ.p4();
    TLorentzVector WZ2=Z2.p4()+WZ.p4();
    TLorentzVector Z1Z2=Z1.p4()+Z2.p4();
    TLorentzVector WZZ=Z1.p4()+Z2.p4()+WZ.p4();
    double massWZ1=sqrt((WZ1.E()*WZ1.E())-(WZ1.Px()*WZ1.Px())-(WZ1.Py()*WZ1.Py())-(WZ1.Pz()*WZ1.Pz()));
    double massWZ2=sqrt((WZ2.E()*WZ2.E())-(WZ2.Px()*WZ2.Px())-(WZ2.Py()*WZ2.Py())-(WZ2.Pz()*WZ2.Pz()));
    double massZ1Z2=sqrt((Z1Z2.E()*Z1Z2.E())-(Z1Z2.Px()*Z1Z2.Px())-(Z1Z2.Py()*Z1Z2.Py())-(Z1Z2.Pz()*Z1Z2.Pz()));
    double massWZZ=sqrt((WZZ.E()*WZZ.E())-(WZZ.Px()*WZZ.Px())-(WZZ.Py()*WZZ.Py())-(WZZ.Pz()*WZZ.Pz()));
    theHistograms.fill("mass of ZZ","ZZ mass",100,150,1000,massZ1Z2);
    if(WZ.id()==23){
      theHistograms.fill("mass of ZZ","ZZ mass",100,150,1000,massWZ1);
      theHistograms.fill("mass of ZZ","ZZ mass",100,150,1000,massWZ2);}
    if(WZ.id()==24){
      theHistograms.fill("mass of WZ","WZ mass",100,150,1000,massWZ1);
      theHistograms.fill("mass of WZ","WZ mass",100,150,1000,massWZ2);}
    theHistograms.fill("mass of tribosons","Triboson mass",50,200,2500,massWZZ);
    theHistograms.fill("mass of coupled bosons","Coupled boson mass",50,150,1000,massWZ1);
    theHistograms.fill("mass of coupled bosons","Coupled boson mass",50,150,1000,massWZ2);
    theHistograms.fill("mass of coupled bosons","Coupled boson mass",50,150,1000,massZ1Z2);
    theHistograms.fill("energy of all bosons","Total boson energy",50,0,5000,WZZ.E());
    theHistograms.fill("energy of coupled bosons","Coupled boson energy",50,0,3000,Z1.e()+Z2.e());
    theHistograms.fill("energy of coupled bosons","Coupled boson energy",50,0,3000,Z1.e()+WZ.e());
    theHistograms.fill("energy of coupled bosons","Coupled boson energy",50,0,3000,WZ.e()+Z2.e());
    double angleWZ1=Z1.p4().Angle(WZ.p4().Vect());
    double angleWZ2=Z2.p4().Angle(WZ.p4().Vect());
    double angleZZ=Z1.p4().Angle(Z2.p4().Vect());
    double anglemax=max(angleWZ1,max(angleWZ2,angleZZ));
    double anglemin=min(angleWZ1,min(angleWZ2,angleZZ));
    theHistograms.fill("boson relative angle","Boson relative angle",50,0,3.5,angleWZ1);
    theHistograms.fill("boson relative angle","Boson relative angle",50,0,3.5,angleWZ2);
    theHistograms.fill("boson relative angle","Boson relative angle",50,0,3.5,angleZZ);
    theHistograms.fill("maximum boson relative angle","Maximum boson relative angle",100,0,3.5,anglemax);
    theHistograms.fill("minimum boson relative angle","Minimum boson relative angle",100,0,3.5,anglemin);
    theHistograms.fill("difference between max/min relative angle","Difference between max/min relative boson angle",100,0,3.5,anglemax-anglemin);
    theHistograms.fill("ratio between max/min relative angle","Ratio between max/min relative boson angle",200,0,30,anglemax/anglemin);
    theHistograms.fill("sum of relative angles","Sum of boson relative angles",300,0,10.5,angleWZ1+angleWZ2+angleZZ);
    double emax=max(WZ.e(),max(Z1.e(),Z2.e()));
    double emin=min(WZ.e(),min(Z1.e(),Z2.e()));
    theHistograms.fill("energy of major bosons","Major boson energy",50,100,1800,emax);
    theHistograms.fill("energy of minor bosons","Minor boson energy",100,0,1000,emin);
    theHistograms.fill("deltae between major/minor bosons","Energy difference between major/minor bosons",100,0,1800,emax-emin);
    theHistograms.fill("energy ratio between major major/minor bosons","Ratio between major/minor boson energy",150,0,30,emax/emin);
    double ptmax=max(WZ.pt(),max(Z1.pt(),Z2.pt()));
    double ptmin=min(WZ.pt(),min(Z1.pt(),Z2.pt()));
    theHistograms.fill("pt of major bosons","Major boson pt",50,100,1000,ptmax);
    theHistograms.fill("pt of minor bosons","Minor boson pt",100,0,400,ptmin);
    theHistograms.fill("deltapt between major/minor bosons","Pt difference between major/minor bosons",100,0,800,ptmax-ptmin);
    theHistograms.fill("pt ratio between major/minor bosons","Pt ratio between major/minor bosons",150,0,30,ptmax/ptmin);
    theHistograms.fill("total pt(abs value) sum","Sum of total pt (in abs value)",100,0,2000,WZ.pt()+Z1.pt()+Z2.pt());
    theHistograms.fill("ratio e max/pt max","Ratio between maximum energy and maximum pt",190,1,20,emax/ptmax);
    theHistograms.fill("ratio e min/pt min","Ratio between minimum energy and minimum pt",190,1,20,emin/ptmin);
    theHistograms.fill("rapporto e max/pt min","Ratio between maximum energy and minimum pt",200,0,80,emax/ptmin);
    theHistograms.fill("rapporto e min/pt max","Ratio between minimum energy and minimum pt",100,0,7,emin/ptmax);
    double etamax=max(abs(WZ.eta()),max(abs(Z1.eta()),abs(Z2.eta())));
    double etamin=min(abs(WZ.eta()),min(abs(Z1.eta()),abs(Z2.eta())));
    theHistograms.fill("eta max","Maximum eta (in abs value) of bosons",100,0,7,etamax);
    theHistograms.fill("eta min","Minimum eta (in abs value) of bosons",100,0,4,etamin);
    theHistograms.fill("deltaeta max/min","Difference between max/min boson eta(in abs value)",100,0,7,etamax-etamin);
    theHistograms.fill("ratio eta max/min","Ratio eta max/min (in abs value)",100,1,20,etamax/etamin);}

  //analisys of leptonic Z bosons
  
  foreach(const phys::Boson<phys::Particle> lepbos,genVBHelper_.ZtoChLep()){
    double elepmin= min(lepbos.daughter(0).e(),lepbos.daughter(1).e());
    double elepmax= max(lepbos.daughter(0).e(),lepbos.daughter(1).e());
    theHistograms.fill("energy of major leptons","Major lepton energy",25,0,1000,elepmax);
    theHistograms.fill("energy of minor leptons","Minor lepton energy",25,0,600,elepmin);
    theHistograms.fill("mass of leptonic bosons","Leptonic boson mass",25,70,110,lepbos.mass());
    theHistograms.fill("energy of leptonic bosons","Leptonic boson energy",25,0,1000,lepbos.e());
    theHistograms.fill("pt of leptonic bosons","Leptonic boson pt",25,0,500,lepbos.pt());
  }
  
    
  }
