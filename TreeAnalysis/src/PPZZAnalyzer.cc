#include "VVXAnalysis/TreeAnalysis/interface/PPZZAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <TMath.h>
#include <TF1.h>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/assign/std/vector.hpp>
using namespace boost::assign;

using std::cout;
using std::endl;


using namespace phys;

Double_t acoplanarity(phys::DiBoson <phys::Lepton,phys::Lepton> *ZZ);
bool regionselection(double x, double y);
double expfunc(double *x, double *par);

Int_t PPZZAnalyzer::cut() {
  
  return 1;
}

void PPZZAnalyzer::analyze(){

  std::vector<phys::Proton> eventprotons;
  //forw: sector 45 (ALICE), back: sector 56 (LHCb)
  std::vector<phys::Proton> forweventprotons;
  std::vector<phys::Proton> backeventprotons;
  
  foreach(const phys::Proton mproton, *multiRPprotons){
    if(mproton.appropriatexi()){
      eventprotons.push_back(mproton);
      if(mproton.LHCSector()){
	forweventprotons.push_back(mproton);
        theHistograms->fill("xiforw","xi of forward protons",30,0,0.22,mproton.xi(),theWeight);
      }
      if(!mproton.LHCSector()){
	backeventprotons.push_back(mproton);
        theHistograms->fill("xiback","xi of backward protons",30,0,0.22,mproton.xi(),theWeight);
      }
    }
  }

  theHistograms->fill("protonmult","proton multiplicity",10,-0.5,9.5,eventprotons.size(),theWeight);
  theHistograms->fill("forwprotonmult","forward proton multiplicity",10,-0.5,9.5,forweventprotons.size(),theWeight);
  theHistograms->fill("backprotonmult","backward proton multiplicity",10,-0.5,9.5,backeventprotons.size(),theWeight);

  
  if(ZZ->passFullSelection()){
    theHistograms->fill("mZZ","mass of ZZ pair",50,300,10000,ZZ->mass(),theWeight);
    theHistograms->fill("yZZ","rapidity of ZZ pair",50,-1,1,ZZ->rapidity(),theWeight);
    if(eventprotons.size()>0){
     theHistograms->fill("mZZoneproton","mass of ZZ pair with one reconstructed proton",50,300,10000,ZZ->mass(),theWeight);
    }
    if(eventprotons.size()>1){
     theHistograms->fill("mZZtwoprotons","mass of ZZ pair with two reconstructed protons",50,300,10000,ZZ->mass(),theWeight);
    }
  }
  

  phys::Proton p1, p2;
  phys::ProtonPair pp;
  
  if(forweventprotons.size()&&backeventprotons.size()){
    int nmatches=0;
    bool matched=false;
    bool matcheddelta=false;
    double bestmassdif=9999.;
    double bestydif=9999.;
    for(unsigned int i=0; i<eventprotons.size();i++){
      p1=eventprotons[i];
      if(p1.LHCSector()) continue;
      for(unsigned int j=0; j<eventprotons.size();j++){
        p2=eventprotons[j];
        if(!p2.LHCSector()) continue;
	pp=phys::ProtonPair(p1,p2);
        theHistograms->fill("mpp","Expected mass of ZZ pair",20,450,2400,pp.mpp(),theWeight);
        theHistograms->fill("ypp","Expected rapidity of ZZ pair",15,0,0.75,pp.ypp(),theWeight);

        if(ZZ->passFullSelection()){
	  double a = acoplanarity(ZZ);
	  theHistograms->fill("acoplanarity","Acoplanarity distribution",10,0,0.002,a,theWeight);
	  double massdif = 1-ZZ->mass()/pp.mpp();
	  double ydif = pp.ypp()-ZZ->rapidity();
	  nmatches++;
	  theHistograms->fill("th2","2D matching distribution",60,-2.5,0.5,60,-1.5,1.5,massdif,ydif,theWeight);
	  bool region = regionselection(massdif,ydif);
	  if(massdif*massdif+ydif*ydif<bestmassdif*bestmassdif+bestydif*bestydif && matcheddelta==false){
	      bestmassdif=massdif;
	      bestydif=ydif;}
	  if(region){
	    if(region) theHistograms->fill("deltaPhi","deltaPhi distribution",10,0,3.14,physmath::deltaPhi(pp.ppp4(),ZZ->p4()));
	    if(matched==false){
	      matched=true;}		
	    if(massdif>-0.05){
	      matcheddelta=true;}
	}
      }}}
    theHistograms->fill("nmatches","Number of possible matches for event",10,-0.5,9.5,nmatches,theWeight);
    theHistograms->fill("th2good","2D matching distribution",60,-2.5,0.5,60,-1.5,1.5,bestmassdif,bestydif,theWeight);
    if(matched==true){
      theHistograms->fill("counteromicron","counteromicron",1,0,1,0.5,theWeight);
      theHistograms->fill("th2matched","2D matching distribution",60,-2.5,0.5,60,-1.5,1.5,bestmassdif,bestydif,theWeight);
      theHistograms->fill("goodmZZ","Mass of matched ZZ pairs",10,650,4500,ZZ->mass(),theWeight);
      theHistograms->fill("goodyZZ","Rapidity of matched ZZ pairs",15,-1.3,1.3,ZZ->rapidity(),theWeight);
    }
    if(matcheddelta==true) theHistograms->fill("counterdelta","counterdelta",1,0,1,0.5,theWeight);
}

  //test: use the highest xi proton every time
  
  if(forweventprotons.size()&&backeventprotons.size()){
    double ximaxforw=0;
    double ximaxback=0;
    phys::Proton forwproton, backproton;
    foreach(const phys::Proton forwp, forweventprotons){
      if(forwp.xi()>ximaxforw){
	ximaxforw=forwp.xi();
	forwproton=forwp;
      }
    }
    foreach(const phys::Proton backp, backeventprotons){
      if(backp.xi()>ximaxback){
	ximaxback=backp.xi();
	backproton=backp;
      }
    }
    phys::ProtonPair ppxi=phys::ProtonPair(backproton,forwproton);
    double massdifxi = 1-ZZ->mass()/ppxi.mpp();
    double ydifxi = ppxi.ypp()-ZZ->rapidity();
    theHistograms->fill("th2xi","2D matching distribution",60,-2.5,0.5,60,-1.5,1.5,massdifxi,ydifxi,theWeight);
    bool regionxi = regionselection(massdifxi,ydifxi);
    if(regionxi){
      theHistograms->fill("counterxi","counterxi",1,0,1,0.5,theWeight);
      if(massdifxi>-0.05) theHistograms->fill("counterdeltaxi","counterdeltaxi",1,0,1,0.5,theWeight);} 
  }

  
  eventprotons.clear();
  forweventprotons.clear();
  backeventprotons.clear();
}

void PPZZAnalyzer::finish(){ 
  TH1 *counteromicron = theHistograms->get("counteromicron");
  TH1 *counterdelta = theHistograms->get("counterdelta");
  cout<<"Algorithm A (min distance):"<<endl;
  cout<<"Total # of events in the signal region: "<<counteromicron->GetEntries()*theWeight<<endl;
  cout<<"# of events in the delta signal region: "<<counterdelta->GetEntries()*theWeight<<endl;
  cout<<"# of events in the omicron signal region: "<<(counteromicron->GetEntries()-counterdelta->GetEntries())*theWeight<<endl<<endl;
  TH1 *counteromicronxi = theHistograms->get("counterxi");
  TH1 *counterdeltaxi = theHistograms->get("counterdeltaxi");
  cout<<"Algorithm B (max xi):"<<endl;
  cout<<"Total # of events in the signal region: "<<counteromicronxi->GetEntries()*theWeight<<endl;
  cout<<"# of events in the delta signal region: "<<counterdeltaxi->GetEntries()*theWeight<<endl;
  cout<<"# of events in the omicron signal region: "<<(counteromicronxi->GetEntries()-counterdeltaxi->GetEntries())*theWeight<<endl<<endl;
}

Double_t acoplanarity(phys::DiBoson <phys::Lepton,phys::Lepton> *ZZ){
  double_t a = abs(1-physmath::deltaPhi(ZZ->first().phi(),ZZ->second().phi())/TMath::Pi());
  return a;
}

//out=0 -> not signal
//out=1 -> signal

bool regionselection(double x, double y){
  int out=0;
  TF1 *f1 = new TF1("func1",expfunc,-2,-0.05,1);
  f1->SetParameter(0,0.95);
  TF1 *f2 = new TF1("func2",expfunc,-2,-0.05,1);
  f2->SetParameter(0,1.15);
  if(x<-2||x>0.1) {out=false;}
  else if(x<=-0.05) {out=abs(y) >= f1->Eval(x) && abs(y) <= f2->Eval(x);}
  else if(x<=0) {out=abs(y)<=x+0.1;}
  else if(x<=0.05) {out=abs(y)<=0.1;}
  else {out=abs(y)<=-x+0.15;}
  return out;
}

double expfunc(double *x, double *par){
  double xx=x[0];
  return log(par[0]*(1-xx));
}
