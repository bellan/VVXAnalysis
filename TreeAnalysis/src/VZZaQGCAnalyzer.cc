/*
In order to use this analysis, one needs to follow these steps:
1. Analyze the weighting samples of generated data (SM,linear and quadratic) with VVXnocutsAnalyzer. The python script multipleAnalysis.py could be useful to save time.
2. Use the ROOT macro VZZDistributionComparison to obtain the weighting histogram (the variable used for the weighting may be changed). 
3. Analyze the sample with this analysis.
4. The result of this analysis, as far as aGCs go, can and should be extracted with the ROOT macro aGCRecosntructedPlot.

Naming conventions: 
histograms that are meant to represent aGCs data (and therefore have theWeight) have the same name as their SM counterpart except with "weighted" in front. This is useful to have one less parameter in aGCReconstructedPlot.
histograms that include "good","well generated" or "well reconstructed" in their names are the ones used to study efficiency. In order to do that, it is advisable to use the macro VZZaQGCEfficiency. 
*/


#include "VVXAnalysis/TreeAnalysis/interface/VZZaQGCAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include <TFile.h>
#include <TRandom3.h>

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

  //to check whether the event is a signal one for generated bosons
  
  int weightedevent=0;

  //to reweight aGCs histograms
  
  double weight=0.;
  
  //assignment and properties of generated bosons
  
  TLorentzVector a;
  phys::Boson<phys::Particle> z1,z2;
  foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
    if((abs(genVBParticle.daughter(0).id())==11||abs(genVBParticle.daughter(0).id())==13)&&genVBParticle.id()==23){
      theHistograms->fill("theta of generated bosons","Generated boson theta",75 ,0,3.5,genVBParticle.p4().Theta(),theWeight);
      theHistograms->fill("eta of generated bosons","Generated boson eta",60,-6,6,genVBParticle.eta(),theWeight);
      double angle1=genVBParticle.daughter(0).p4().Angle(genVBParticle.daughter(1).p4().Vect());
      theHistograms->fill("angle of generated leptons","Generated lepton angle",75,0,3.5,angle1,theWeight);
      theHistograms->fill("energy of major lepton","Generated major lepton energy",200,0,2000,std::max(genVBParticle.daughter(0).e(),genVBParticle.daughter(1).e()),theWeight);
      theHistograms->fill("energy of minor lepton","Generated minor lepton energy",200,0,800,std::min(genVBParticle.daughter(0).e(),genVBParticle.daughter(1).e()),theWeight);
      a+=genVBParticle.p4();
     }}

  //study of generated ZZ pair
  
  if(topology.test(0)){
    
    double massZZ=sqrt((a.E()*a.E())-(a.Px()*a.Px())-(a.Py()*a.Py())-(a.Pz()*a.Pz()));
    theHistograms->fill("mass of generated dibosons","Generated ZZ mass",78,80,600,massZZ,theWeight);

    foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
      if((abs(genVBParticle.daughter(0).id())==11||abs(genVBParticle.daughter(0).id())==13)&&genVBParticle.id()==23){
	if(z1.mass()<0.01){
	  z1=genVBParticle;}
	else{if(abs(genVBParticle.mass()-phys::ZMASS)<abs(z1.mass()-phys::ZMASS)){
	    z2=z1;
	    z1=genVBParticle;}
	  else{z2=genVBParticle;}}}}
    
    theHistograms->fill("mass of generated Z1","Generated Z1 mass",180,50,130,z1.mass(),theWeight);
    theHistograms->fill("mass of generated Z2","Generated Z2 mass",90,50,130,z2.mass(),theWeight);
    theHistograms->fill("mass of well generated Z1","Generated Z1 mass",10,80,100,z1.mass(),theWeight);
    theHistograms->fill("mass of well generated Z2","Generated Z2 mass",10,70,110,z2.mass(),theWeight);
    
    theHistograms->fill("energy of generated Z1","Generated Z1 energy",135,0,2700,z1.e(),theWeight);
    theHistograms->fill("energy of generated Z2","Generated Z2 energy",135,0,2700,z2.e(),theWeight);
    theHistograms->fill("energy of well generated Z1","Generated Z1 energy",10,0,1400,z1.e(),theWeight);
    theHistograms->fill("energy of well generated Z2","Generated Z2 energy",10,0,1400,z2.e(),theWeight);
    
    theHistograms->fill("pt of generated Z1","Generated Z1 pt",150,0,900,z1.pt(),theWeight);
    theHistograms->fill("pt of generated Z2","Generated Z2 pt",150,0,900,z2.pt(),theWeight);
    theHistograms->fill("pt of well generated Z1","Generated Z1 pt",10,0,350,z1.pt(),theWeight);
    theHistograms->fill("pt of well generated Z2","Generated Z2 pt",10,0,350,z2.pt(),theWeight);
    
    theHistograms->fill("eta of well generated Z1","Generated Z1 eta",10,-4.5,4.5,z1.eta(),theWeight);
    theHistograms->fill("eta of well generated Z2","Generated Z2 eta",10,-4.5,4.5,z2.eta(),theWeight);
    
    theHistograms->fill("energy of good Z1 major lepton","Generated Z1 major lepton's energy",10,0,1000,std::max(z1.daughter(0).e(),z1.daughter(1).e()),theWeight);
    theHistograms->fill("energy of good Z1 minor lepton","Generated Z1 minor lepton's energy",10,0,300,std::min(z1.daughter(0).e(),z1.daughter(1).e()),theWeight);
    theHistograms->fill("energy of good Z2 major lepton","Generated Z2 major lepton's energy",10,0,1000,std::max(z2.daughter(0).e(),z2.daughter(1).e()),theWeight);
    theHistograms->fill("energy of good Z2 minor lepton","Generated Z2 minor lepton's energy",10,0,300,std::min(z2.daughter(0).e(),z2.daughter(1).e()),theWeight);
    
    double ZZangle1=z1.p4().Angle(z2.p4().Vect());
    theHistograms->fill("angolo bosoni generati","Angolo bosoni generati",35,0,3.5,ZZangle1,theWeight);}

  //study of reconstructed boson properties
    
  if(ZZ->passFullSelection()){
    
    theHistograms->fill("mass of reconstructed Z1","Reconstructed Z1 mass",180,50,130,ZZ->first().mass(),theWeight);
    theHistograms->fill("mass of reconstructed Z2","Reconstructed Z2 mass",180,50,130,ZZ->second().mass(),theWeight);
      
    theHistograms->fill("pt of reconstructed Z1","Reconstructed Z1 pt",150,0,900,ZZ->first().pt(),theWeight);
    theHistograms->fill("pt of reconstructed Z2","Reconstructed Z2 pt",150,0,900,ZZ->second().pt(),theWeight);
      
    theHistograms->fill("theta of reconstructed bosons","Reconstructed boson theta",75 ,0,3.5,ZZ->first().p4().Theta(),theWeight);
    theHistograms->fill("theta of reconstructed bosons","Reconstructed boson theta",75 ,0,3.5,ZZ->second().p4().Theta(),theWeight);
      
    theHistograms->fill("energy of reconstructed Z1","Reconstructed Z1 energy",135,0,2700,ZZ->first().e(),theWeight);
    theHistograms->fill("energy of reconstructed Z2","Reconstructed Z2 energy",135,0,2700,ZZ->second().e(),theWeight);
      
    theHistograms->fill("eta of reconstructed bosons","Reconstructed boson eta",60,-6,6,ZZ->first().eta(),theWeight);
    theHistograms->fill("eta of reconstructed bosons","Reconstructed boson eta",60,-6,6,ZZ->second().eta(),theWeight);
      
    double angle2=ZZ->first().daughter(0).p4().Angle(ZZ->first().daughter(1).p4().Vect());
    theHistograms->fill("angle of reconstructed leptons","Reconstructed lepton angle",75,0,3.5,angle2,theWeight);
    double angle3=ZZ->second().daughter(0).p4().Angle(ZZ->second().daughter(1).p4().Vect());
    theHistograms->fill("angle of reconstructed leptons","Reconstructed lepton angle",75,0,3.5,angle3,theWeight);
    double ZZangle2=ZZ->first().p4().Angle(ZZ->second().p4().Vect());
    theHistograms->fill("angle of reconstructed bosons","Reconstructed boson angle",35,0,3.5,ZZangle2,theWeight);

    TLorentzVector b=ZZ->first().p4()+ZZ->second().p4();
    double mass2=sqrt((b.E()*b.E())-(b.Px()*b.Px())-(b.Py()*b.Py())-(b.Pz()*b.Pz()));
    theHistograms->fill("mass of reconstructed dibosons","Reconstructed ZZ mass",78,80,600,mass2,theWeight);

    theHistograms->fill("energy of major reconstructed leptons","Major reconstructed lepton energy",200,0,2000,std::max(ZZ->first().daughter(0).e(),ZZ->first().daughter(1).e()),theWeight);
    theHistograms->fill("energy of minor reconstructed leptons","Minor reconstructd lepton energy",200,0,800,std::min(ZZ->first().daughter(0).e(),ZZ->first().daughter(1).e()),theWeight);
    theHistograms->fill("energy of major reconstructed leptons","Major reconstructed lepton energy",200,0,2000,std::max(ZZ->second().daughter(0).e(),ZZ->second().daughter(1).e()),theWeight);
    theHistograms->fill("energy of minor reconstructed leptons","Minor reconstructed lepton energy",200,0,800,std::min(ZZ->second().daughter(0).e(),ZZ->second().daughter(1).e()),theWeight);

    //coupling generated and reconstructed leptonic bosons

    double dR=9999;
    double dR2=9999;
    double enlep1,enlep2,enlep3,enlep4;
    phys::Boson<phys::Particle> z3,z4;
    foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
      if((abs(genVBParticle.daughter(0).id())==11||abs(genVBParticle.daughter(0).id())==13)&&genVBParticle.id()==23){
	if(physmath::deltaR(genVBParticle,ZZ->first())<dR){
	  dR=physmath::deltaR(genVBParticle,ZZ->first());
	  z3=genVBParticle;
	  enlep1=std::max(genVBParticle.daughter(0).e(),genVBParticle.daughter(1).e());
	  enlep2=std::min(genVBParticle.daughter(0).e(),genVBParticle.daughter(1).e());}}}
    theHistograms->fill("dR Z1","dR Z1",100,0,0.5,dR);
      
    foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
      if((abs(genVBParticle.daughter(0).id())==11||abs(genVBParticle.daughter(0).id())==13)&&genVBParticle.id()==23){
	if(physmath::deltaR(genVBParticle,ZZ->second())<dR2){
	  dR2=physmath::deltaR(genVBParticle,ZZ->second());
	  z4=genVBParticle;
	  enlep3=std::max(genVBParticle.daughter(0).e(),genVBParticle.daughter(1).e());
	  enlep4=std::min(genVBParticle.daughter(0).e(),genVBParticle.daughter(1).e());}}}
    theHistograms->fill("dR Z2","dR Z2",100,0,0.5,dR2);

    if(dR<0.05&&dR2<0.1&&z3!=z4){
	
      theHistograms->fill("comparison of Z1 mass","Difference between generated/reconstructed Z1 mass",40,-8,8,ZZ->first().mass()-z3.mass(),theWeight);
      theHistograms->fill("comparison of Z1 energy","Difference between generated/reconstructed Z1 energy",200,-40,40,ZZ->first().e()-z3.e(),theWeight);
      theHistograms->fill("comparison of Z1 pt","Difference between generated/reconstructed Z1 pt",150,-30,30,ZZ->first().pt()-z3.pt());
      theHistograms->fill("comparison of Z1 eta","Difference between generated/reconstructed Z1 eta",50,-0.1,0.1,ZZ->first().eta()-z3.eta(),theWeight);
	
      theHistograms->fill("mass of well reconstructed Z1","Reconstructed Z1 mass",10,80,100,z3.mass(),theWeight);
      theHistograms->fill("pt of well reconstructed Z1","Reconstructed Z1 pt",10,0,350,z3.pt(),theWeight);
      theHistograms->fill("energy of well reconstructed Z1","Reconstructed Z1 energy",10,0,1400,z3.e(),theWeight);
      theHistograms->fill("eta of well reconstructed Z1","Reconstructed Z1 eta",10,-4.5,4.5,z3.eta(),theWeight);
	
      theHistograms->fill("comparison of Z1 major lepton energy","Difference between generated/reconstructed Z1 major lepton energy",100,-20,20,std::max(ZZ->first().daughter(0).e(),ZZ->first().daughter(1).e())-enlep1,theWeight);
      theHistograms->fill("comparison of Z1 minor lepton energy","Difference between generated/reconstructed Z1 minor lepton energy",00,-20,20,std::min(ZZ->first().daughter(0).e(),ZZ->first().daughter(1).e())-enlep2,theWeight);
      theHistograms->fill("energy of well reconstructed Z1 major lepton","Reconstructed major Z1 lepton energy",10,0,1000,enlep1,theWeight);
      theHistograms->fill("energy of well reconstructed Z1 minor lepton","Reconstructed minor Z1 lepton energy",10,0,300,enlep2,theWeight);
       
      theHistograms->fill("comparison of Z2 mass","Difference between generated/reconstructed Z2 mass",40,-8,8,ZZ->second().mass()-z4.mass(),theWeight);
      theHistograms->fill("comparison of Z2 energy","Difference between generated/reconstructed Z2 energy",200,-40,40,ZZ->second().e()-z4.e(),theWeight);
      theHistograms->fill("comparison of Z2 pt","Difference between generated/reconstructed Z2 pt",150,-30,30,ZZ->second().pt()-z4.pt(),theWeight);
      theHistograms->fill("comparison of Z2 eta","Difference between generated/reconstructed Z2 eta",50,-0.1,0.1,ZZ->second().eta()-z4.eta(),theWeight);

      theHistograms->fill("mass of well reconstructed Z2","Reconstructed Z2 mass",10,70,110,z4.mass(),theWeight);
      theHistograms->fill("pt of well reconstructed Z2","Reconstructed Z2 pt",10,0,350,z4.pt(),theWeight);
      theHistograms->fill("energy of well reconstructed Z2","Reconstructed Z2 energy",10,0,1400,z4.e(),theWeight);
      theHistograms->fill("eta of well reconstructed Z2","Reconstructed Z2 eta",10,-4.5,4.5,z4.eta(),theWeight);
      theHistograms->fill("comparison of Z2 major lepton energy","Difference between generated/reconstructed Z2 major lepton energy",100,-20,20,std::max(ZZ->second().daughter(0).e(),ZZ->second().daughter(1).e())-enlep3,theWeight);
      theHistograms->fill("comparison of Z2 minor lepton energy","Difference between generated/reconstructed Z2 minor lepton energy",100,-20,20,std::min(ZZ->second().daughter(0).e(),ZZ->second().daughter(1).e())-enlep4,theWeight);
      theHistograms->fill("energy of well reconstructed Z2 major lepton","Reconstructed Z2 major lepton energy",10,0,1000,enlep3,theWeight);
      theHistograms->fill("energy of well reconstructed Z2 minor lepton","Reconstructed Z2 minor lepton energy",10,0,300,enlep4,theWeight); }
      }
  
  //analysis of generated hadronic bosons
  
  int nquark=0;
  phys::Particle quark1,quark2;

  foreach(const phys::Particle genParticle,*genParticles){
    if(abs(genParticle.id())>0&&abs(genParticle.id())<7){
      nquark++;
      if(quark1.e()<0.01){
	quark1=genParticle;}
      else if(quark2.e()<0.01){
	quark2=genParticle;}}}

  double totalenergy;

  if(nquark==2){
    phys::Boson<phys::Particle> jetboson=phys::Boson<phys::Particle>(quark1,quark2);

    theHistograms->fill("generated hadronic boson mass","generated mass of boson with hadronic decay",50,70,110,jetboson.mass(),theWeight);

    //event weighting
    
    if(topology.test(0)){

      weightedevent++;

      totalenergy=z1.e()+z2.e()+jetboson.e();

      TFile *weightfile = TFile::Open("./aGCweighthist.root");
      TH1F *weighthist= (TH1F*)weightfile->Get("energy of all bosons");
      double weightarray[50];
      for(int i=0;i<50;i++){
        weightarray[i]=weighthist->GetBinContent(i+1);}
      for(int j=0;j<50;j++){
        if(totalenergy>60*j&&totalenergy<60*(j+1)){
	  weight=weightarray[j];}}
      weightfile->Close();

      theHistograms->fill("generated energy of all bosons","Sum of generated boson energies",10,0,3000,totalenergy,theWeight);
      theHistograms->fill("good generated energy of all bosons","Sum of generated boson energies",5,0,3000,totalenergy,theWeight);
      theHistograms->fill("weighted generated energy of all bosons","Sum of generated boson energies",10,0,3000,totalenergy,theWeight*weight);
    }
  } 
  
  //hadronic boson reconstruction

  Particle* hadVB = nullptr;
	
  // First look for a pair of AK4
  
  if(jets->size() > 2){
    size_t j4_size = jets->size();
    vector<Boson<Jet>> AK4pairs;
    AK4pairs.reserve( j4_size*(j4_size - 1)/2 );
    for(size_t i = 0; i < j4_size; ++i){
      for(size_t j = i+1; j < j4_size; ++j){
	AK4pairs.emplace_back(Boson<Jet>(jets->at(i), jets->at(j))); }}
	       
    std::sort(AK4pairs.begin(), AK4pairs.end(), Mass2Comparator(phys::WMASS, phys::ZMASS));
		
    if(AK4pairs.front().mass() > 60 && AK4pairs.front().mass() < 120){
      hadVB = new Boson<Jet>(AK4pairs.front());}
	}
	
  // If no pair of AK4, look for single AK8
  
  if(jetsAK8->size() > 1 && !hadVB){
    std::sort(jetsAK8->begin(), jetsAK8->end(), Mass2Comparator(phys::WMASS, phys::ZMASS));
    if(jetsAK8->front().chosenAlgoMass() > 60 && jetsAK8->front().chosenAlgoMass() < 120){
      hadVB = new Jet(jetsAK8->front());}
	}
  
  Particle WZ;
  if(hadVB){
    WZ=*hadVB;
    theHistograms->fill("reconstructed hadronic boson mass","reconstructed mass of boson with hadronic decay",15,60,120,WZ.mass(),theWeight);}

  //coupling of generated and reconstructed hadronic bosons

  double dRhad=9999;
  phys::Boson<phys::Particle> Wgen;
  foreach(const phys::Boson<phys::Particle> genVBParticle,*genVBParticles){
    if(abs(genVBParticle.daughter(0).id())<7&&abs(genVBParticle.daughter(1).id())<7){
      if(physmath::deltaR(genVBParticle,WZ)<dRhad){
	dRhad=physmath::deltaR(genVBParticle,Wgen);
	Wgen=genVBParticle;
      }}}
  theHistograms->fill("dR hadronic boson","dR hadronic boson",20,0,5,dRhad,theWeight);
  if(dRhad<5){
    theHistograms->fill("comparison of w mass","Difference between generated/reconstructed mass of boson with hadronic decay",40,-40,40,WZ.mass()-Wgen.mass(),theWeight);}

  //study of the final state with three reconstructed bosons 

  if(hadVB&&(weightedevent==1&&ZZ->first().e()>0.01)){

    phys::Boson<phys::Lepton> Z1=ZZ->first();
    phys::Boson<phys::Lepton> Z2=ZZ->second();
      
    TLorentzVector WZ1=Z1.p4()+WZ.p4();
    TLorentzVector WZ2=Z2.p4()+WZ.p4();
    TLorentzVector Z1Z2=Z1.p4()+Z2.p4();
    TLorentzVector WZZ=Z1.p4()+Z2.p4()+WZ.p4();
    double massWZ1=sqrt((WZ1.E()*WZ1.E())-(WZ1.Px()*WZ1.Px())-(WZ1.Py()*WZ1.Py())-(WZ1.Pz()*WZ1.Pz()));
    double massWZ2=sqrt((WZ2.E()*WZ2.E())-(WZ2.Px()*WZ2.Px())-(WZ2.Py()*WZ2.Py())-(WZ2.Pz()*WZ2.Pz()));
    double massZ1Z2=sqrt((Z1Z2.E()*Z1Z2.E())-(Z1Z2.Px()*Z1Z2.Px())-(Z1Z2.Py()*Z1Z2.Py())-(Z1Z2.Pz()*Z1Z2.Pz()));
    double massWZZ=sqrt((WZZ.E()*WZZ.E())-(WZZ.Px()*WZZ.Px())-(WZZ.Py()*WZZ.Py())-(WZZ.Pz()*WZZ.Pz()));
    theHistograms->fill("mass of ZZ","ZZ mass",10,150,1000,massZ1Z2,theWeight);
    theHistograms->fill("mass of WZ1","WZ1 mass",10,150,1000,massWZ1,theWeight);
    theHistograms->fill("mass of WZ2","WZ2 mass",10,150,1000,massWZ2,theWeight);
    theHistograms->fill("mass of WZ","WZ mass",10,150,1000,massWZ1,theWeight);
    theHistograms->fill("mass of WZ","WZ mass",10,150,1000,massWZ2,theWeight);
    theHistograms->fill("mass of tribosons","Triboson mass",5,200,1200,massWZZ,theWeight);
    theHistograms->fill("weighted mass of tribosons","Triboson mass",5,200,1200,massWZZ,theWeight*weight);

    theHistograms->fill("energy of all bosons","Total boson energy",5,0,2200,WZZ.E(),theWeight);
    theHistograms->fill("weighted energy of all bosons","Weighted total boson energy",5,0,2200,WZZ.E(),theWeight*weight);
    theHistograms->fill("resolution of total energy","Total energy resolution",1000,-500,500,WZZ.E()-totalenergy,theWeight);
    
    theHistograms->fill("energy of ZZ","ZZ energy",10,0,2000,Z1Z2.E(),theWeight);
    theHistograms->fill("energy of WZ1","WZ1 energy",10,0,2000,WZ1.E(),theWeight);
    theHistograms->fill("energy of WZ2","WZ2 energy",10,0,2000,WZ2.E(),theWeight);
    theHistograms->fill("energy of Z1","Z1 energy",10,0,1500,Z1.e(),theWeight);
    theHistograms->fill("weighted energy of Z1","Z1 energy",10,0,1500,Z1.e(),theWeight*weight);
    theHistograms->fill("energy of Z2","Z2 energy",10,0,1500,Z2.e(),theWeight);
    theHistograms->fill("energy of W","W energy",10,0,1500,WZ.e(),theWeight);
  
    double angleWZ1=Z1.p4().Angle(WZ.p4().Vect());
    double angleWZ2=Z2.p4().Angle(WZ.p4().Vect());
    double angleZZ=Z1.p4().Angle(Z2.p4().Vect());
    theHistograms->fill("WZ1 relative angle","WZ1 relative angle",5,0,3.5,angleWZ1,theWeight);
    theHistograms->fill("WZ2 relative angle","WZ2 relative angle",5,0,3.5,angleWZ2,theWeight);
    theHistograms->fill("ZZ relative angle","ZZ relative angle",5,0,3.5,angleZZ,theWeight);
    theHistograms->fill("weighted WZ1 relative angle","WZ1 relative angle",5,0,3.5,angleWZ1,theWeight*weight);
    theHistograms->fill("weighted WZ2 relative angle","WZ2 relative angle",5,0,3.5,angleWZ2,theWeight*weight);
    theHistograms->fill("weighted ZZ relative angle","ZZ relative angle",5,0,3.5,angleZZ,theWeight*weight);
  
    theHistograms->fill("total pt scalar sum","Scalar pt sum",10,0,1200,Z1.pt()+Z2.pt()+WZ.pt(),theWeight);
    theHistograms->fill("weighted total pt scalar sum","Scalar pt sum",10,0,1200,Z1.pt()+Z2.pt()+WZ.pt(),theWeight*weight);
    theHistograms->fill("total pt vector sum","Vector pt sum",10,0,300,sqrt(WZZ.Px()*WZZ.Px()+WZZ.Py()*WZZ.Py()),theWeight);
  
    theHistograms->fill("pt of Z1","Z1 pt",10,0,400,Z1.pt(),theWeight);
    theHistograms->fill("pt of Z2","Z2 pt",10,0,400,Z2.pt(),theWeight);
    theHistograms->fill("pt of W","W pt",10,0,400,WZ.pt(),theWeight);
  
    double elepZ1min= min(Z1.daughter(0).e(),Z1.daughter(1).e());
    double elepZ1max= max(Z1.daughter(0).e(),Z1.daughter(1).e());
    theHistograms->fill("energy of major Z1 leptons","Major Z1 lepton energy",10,0,800,elepZ1max,theWeight);
    theHistograms->fill("weighted energy of major Z1 leptons","Major Z1 lepton energy",10,0,800,elepZ1max,theWeight*weight);
    theHistograms->fill("energy of minor Z1 leptons","Minor Z1 lepton energy",10,0,1000,elepZ1min,theWeight);
    double elepZ2min= min(Z2.daughter(0).e(),Z2.daughter(1).e());
    double elepZ2max= max(Z2.daughter(0).e(),Z2.daughter(1).e());
    theHistograms->fill("energy of major Z2 leptons","Major Z2 lepton energy",10,0,2500,elepZ2max,theWeight);
    theHistograms->fill("energy of minor Z2 leptons","Minor Z2 lepton energy",10,0,1000,elepZ2min,theWeight);}
}
