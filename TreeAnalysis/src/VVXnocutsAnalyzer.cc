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

using namespace std;
using namespace colour;
using namespace phys;


Int_t VVXnocutsAnalyzer::cut() {
  
  return 1;
}


void VVXnocutsAnalyzer::analyze(){

  gErrorIgnoreLevel=kFatal; //to avoid errors due to nonexistent tree branches
  
  phys::Boson<phys::Particle> Za,Zb,Z1,Z2,WZ;

  //first choice: analysis for WZZ (uncomment if needed)

  /*
  std::vector<phys::Boson<phys::Particle>>  bosons=genVBHelper_.ZtoChLep();
  foreach(const phys::Boson<phys::Particle> bos,genVBHelper_.ZtoNeutrinos()){
    bosons.push_back(bos);}
  foreach(const phys::Boson<phys::Particle> bos,genVBHelper_.ZtoQ()){
    bosons.push_back(bos);}
  foreach(const phys::Boson<phys::Particle> bos,genVBHelper_.WtoLep()){
    bosons.push_back(bos);}
  foreach(const phys::Boson<phys::Particle> bos,genVBHelper_.WtoQ()){
    bosons.push_back(bos);}
  
  if(bosons.size()==3){
    foreach(const phys::Boson<phys::Particle> bosone,bosons){
      if(WZ==nullboson&&abs(bosone.id())==24){
	WZ=bosone;}
      else if(Z1==nullboson&&abs(bosone.id())==23){
        Z1=bosone;}
      else if(Z2==nullboson&&abs(bosone.id())==23){
	Z2=bosone;}
      else{WZ=bosone;}}
    */

  //second choice: analysis for eemumujj (uncomment if needed)
  
  phys::Particle eminus,eplus,muminus,muplus;
  phys::Particle j1,j2,nullj;
  foreach(const phys::Particle par,*genParticles){
    if(par.id()==11){
      eminus=par;}
    if(par.id()==-11){
      eplus=par;}
    if(par.id()==13){
      muminus=par;}
    if(par.id()==-13){
      muplus=par;}
    if(abs(par.id())<7){
      if(j1==nullj){
	j1=par;}
      else{j2=par;}}
  }
  
  const double zmass=91.1876;
  Za=Boson<phys::Particle>(eminus,eplus);
  Zb=Boson<phys::Particle>(muminus,muplus);
  if(abs(Za.mass()-zmass)<abs(Zb.mass()-zmass)){
    Z1=Za;
    Z2=Zb;}
  else{
    Z1=Zb;
    Z2=Za;}
  WZ=Boson<phys::Particle>(j1,j2);

  //distribution analysis
  
  TLorentzVector WZ1=Z1.p4()+WZ.p4();
  TLorentzVector WZ2=Z2.p4()+WZ.p4();
  TLorentzVector Z1Z2=Z1.p4()+Z2.p4();
  TLorentzVector WZZ=Z1.p4()+Z2.p4()+WZ.p4();
  double massWZ1=sqrt((WZ1.E()*WZ1.E())-(WZ1.Px()*WZ1.Px())-(WZ1.Py()*WZ1.Py())-(WZ1.Pz()*WZ1.Pz()));
  double massWZ2=sqrt((WZ2.E()*WZ2.E())-(WZ2.Px()*WZ2.Px())-(WZ2.Py()*WZ2.Py())-(WZ2.Pz()*WZ2.Pz()));
  double massZ1Z2=sqrt((Z1Z2.E()*Z1Z2.E())-(Z1Z2.Px()*Z1Z2.Px())-(Z1Z2.Py()*Z1Z2.Py())-(Z1Z2.Pz()*Z1Z2.Pz()));
  double massWZZ=sqrt((WZZ.E()*WZZ.E())-(WZZ.Px()*WZZ.Px())-(WZZ.Py()*WZZ.Py())-(WZZ.Pz()*WZZ.Pz()));
  theHistograms.fill("mass of ZZ","ZZ mass",100,150,1000,massZ1Z2);
  theHistograms.fill("mass of WZ1","WZ1 mass",50,150,1000,massWZ1);
  theHistograms.fill("mass of WZ2","WZ2 mass",50,150,1000,massWZ2);
  theHistograms.fill("mass of WZ","WZ mass",100,150,1000,massWZ1);
  theHistograms.fill("mass of WZ","WZ mass",100,150,1000,massWZ2);
  theHistograms.fill("mass of tribosons","Triboson mass",50,200,1200,massWZZ,theWeight); 
  
  theHistograms.fill("energy of all bosons","Total boson energy",50,0,3000,WZZ.E());
  theHistograms.fill("energy of ZZ","ZZ energy",50,0,2000,Z1.e()+Z2.e());
  theHistograms.fill("energy of WZ1","WZ1 energy",50,0,2000,Z1.e()+WZ.e());
  theHistograms.fill("energy of WZ2","WZ2 energy",50,0,2000,WZ.e()+Z2.e());
  theHistograms.fill("energy of Z1","Z1 energy",50,0,1500,Z1.e());
  theHistograms.fill("energy of Z2","Z2 energy",50,0,1500,Z2.e());
  theHistograms.fill("energy of W","W energy",50,0,1500,WZ.e());
  
  double angleWZ1=Z1.p4().Angle(WZ.p4().Vect());
  double angleWZ2=Z2.p4().Angle(WZ.p4().Vect());
  double angleZZ=Z1.p4().Angle(Z2.p4().Vect());
  theHistograms.fill("WZ1 angle","WZ1 angle",50,0,3.5,angleWZ1);
  theHistograms.fill("WZ2 angle","WZ2 angle",50,0,3.5,angleWZ2);
  theHistograms.fill("ZZ angle","ZZ angle",50,0,3.5,angleZZ);
  
  theHistograms.fill("total pt scalar sum","Scalar pt sum",100,0,1200,Z1.pt()+Z2.pt()+WZ.pt());
  theHistograms.fill("total pt vector sum","Vector pt sum",100,0,300,sqrt(WZZ.Px()*WZZ.Px()+WZZ.Py()*WZZ.Py()));
  
  theHistograms.fill("pt of Z1","Z1 pt",100,0,400,Z1.pt());
  theHistograms.fill("pt of Z2","Z2 pt",100,0,400,Z2.pt());
  theHistograms.fill("pt of W","W pt",100,0,400,WZ.pt());
  
  double elepZ1min= min(Z1.daughter(0).e(),Z1.daughter(1).e());
  double elepZ1max= max(Z1.daughter(0).e(),Z1.daughter(1).e());
  theHistograms.fill("energy of major Z1 leptons","Major Z1 lepton energy",50,0,2500,elepZ1max);
  theHistograms.fill("energy of minor Z1 leptons","Minor Z1 lepton energy",50,0,1000,elepZ1min);
  double elepZ2min= min(Z2.daughter(0).e(),Z2.daughter(1).e());
  double elepZ2max= max(Z2.daughter(0).e(),Z2.daughter(1).e());
  theHistograms.fill("energy of major Z2 leptons","Major Z2 lepton energy",50,0,2500,elepZ2max);
  theHistograms.fill("energy of minor Z2 leptons","Minor Z2 lepton energy",50,0,1000,elepZ2min);
}
