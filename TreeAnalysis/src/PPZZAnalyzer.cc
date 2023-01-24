#include "VVXAnalysis/TreeAnalysis/interface/PPZZAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <TMath.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TRandom.h>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/assign/std/vector.hpp>
using namespace boost::assign;

using std::cout;
using std::endl;


using namespace phys;

Double_t acoplanarity(phys::DiBoson <phys::Lepton,phys::Lepton> *ZZ);
bool regionselection(double x, double y, double kmin, double kmax);
double logfunc(double *x, double *par);

Int_t PPZZAnalyzer::cut() {
  
  return 1;
}

void PPZZAnalyzer::analyze(){
  Bool_t effsyst=false;
  Bool_t xisyst=false;
  Bool_t isPUbg=false;
  double kmin=0.90;
  double kmax=1.15;
  double acomax=0.005;
  double xilimit=0.05;
  
  //forw: sector 45 (towards ALICE), back: sector 56 (towards LHCb)
  std::vector<phys::Proton> forweventprotons;
  std::vector<phys::Proton> backeventprotons;

  
  if(!isPUbg){
    foreach(const phys::Proton mproton, *multiRPprotons){
      if(mproton.xi()>=xilimit){
        if(mproton.LHCSector()){
          if(!xisyst){
	    if(!effsyst){forweventprotons.push_back(mproton);
		}
	    else if(gRandom->Rndm()>0.0206) forweventprotons.push_back(mproton);
	  }
	  else {
	    double rand=gRandom->Rndm();
	    std::string era;
	    TFile *filesyst=TFile::Open("data/reco_charactersitics_version1.root","r");
	    if(rand<0.3321) era="2018_preTS1";
	    else if (rand<0.8133) era="2018_TS1_TS2";
	    else if (rand>0.813) era="2018_postTS2";
	    std::string histname=era+"/multi rp-0/"+"xi/g_systematics_vs_xi";
	    TGraphErrors* hs= (TGraphErrors*)filesyst->Get(histname.c_str());
	    double sigma=hs->Eval(mproton.xi());
	    double shiftedxiforw=mproton.xi()+gRandom->Gaus(0.,sigma);
	    const phys::Proton shiftedproton= phys::Proton(shiftedxiforw,true);
	    forweventprotons.push_back(shiftedproton);
	    filesyst->Close();
	    delete filesyst;
	    delete hs;
	  }
          theHistograms->fill("xiforw","xi of forward protons",25,0,0.22,mproton.xi(),theWeight);
        }
        if(!mproton.LHCSector()){
	  if(!xisyst) {
	    if(!effsyst) {backeventprotons.push_back(mproton);
	    }
	    else if(gRandom->Rndm()>0.0221) backeventprotons.push_back(mproton);
	  }
	  else {
	    double rand=gRandom->Rndm();
	    std::string era;
	    TFile *filesyst=TFile::Open("data/reco_charactersitics_version1.root","r");
	    if(rand<0.3321) era="2018_preTS1";
	    else if (rand<0.8133) era="2018_TS1_TS2";
	    else if (rand>0.813) era="2018_postTS2";
	    std::string histname=era+"/multi rp-1/"+"xi/g_systematics_vs_xi";
	    TGraphErrors* hs= (TGraphErrors*)filesyst->Get(histname.c_str());
	    double sigma=hs->Eval(mproton.xi());
	    double shiftedxiback=mproton.xi()+gRandom->Gaus(0.,sigma);
	    const phys::Proton shiftedproton= phys::Proton(shiftedxiback,false);
	    backeventprotons.push_back(shiftedproton);
	    filesyst->Close();
	    delete filesyst;
	    delete hs;
	  }
	  double shiftedxiback=mproton.xi();
          theHistograms->fill("xiback","xi of backward protons",25,0,0.22,mproton.xi(),theWeight);
        }
      }
    }
  }
  else{
    
   TFile *f= new TFile;
   f=TFile::Open("data/NewMixingDistributions2018.root","r");
   TH1D *hn45 = (TH1D *)f->Get("mtpl_multi_S45");
   TH1D *hn56 = (TH1D *)f->Get("mtpl_multi_S56");
   TH1D *hxi45 = (TH1D *)f->Get("xi_multi_S45");
   TH1D *hxi56 = (TH1D *)f->Get("xi_multi_S56");
   f->Close();
   delete f;
   
   double xi45,xi56;
   int n45 = (int)(hn45->GetRandom()+0.5);
   if(n45>0) for(int i=0;i<n45;i++){
       double xi45=hxi45->GetRandom();
       phys::Proton PUP45 = phys::Proton(xi45,true);
       forweventprotons.push_back(PUP45);
     }
   int n56 = (int)(hn56->GetRandom()+0.5);
   if(n56>0) for(int i=0;i<n56;i++){
       double xi56=hxi56->GetRandom();
       phys::Proton PUP56 = phys::Proton(xi56,false);
       backeventprotons.push_back(PUP56);
     } 
   delete hn45;
   delete hn56;
   delete hxi45;
   delete hxi56;
   
  }

  theHistograms->fill("forwprotonmult","forward proton multiplicity",6,-0.5,5.5,forweventprotons.size(),theWeight);
  theHistograms->fill("backprotonmult","backward proton multiplicity",6,-0.5,5.5,backeventprotons.size(),theWeight);
  
  if(ZZ->passFullSelection()){
    theHistograms->fill("mZZ","mass of ZZ pair",20,400,8000,ZZ->mass(),theWeight);
    theHistograms->fill("yZZ","rapidity of ZZ pair",100,-5,5,ZZ->rapidity(),theWeight);
    if(forweventprotons.size()&&backeventprotons.size()) theHistograms->fill("mZZtwoprotons","mass of ZZ pair",20,400,8000,ZZ->mass(),theWeight);
    if(forweventprotons.size()||backeventprotons.size()) theHistograms->fill("mZZoneproton","mass of ZZ pair",20,400,8000,ZZ->mass(),theWeight);
    if(forweventprotons.size()&&backeventprotons.size()) theHistograms->fill("yZZtwoprotons","ZZ rapidity",100,-5,5,ZZ->rapidity(),theWeight);
    if(forweventprotons.size()||backeventprotons.size()) theHistograms->fill("yZZoneproton","ZZ rapidity",100,-5,5,ZZ->rapidity(),theWeight);
  }

  bool zetamarker=false;
  if(forweventprotons.size()&& !(backeventprotons.size())){
    zetamarker=true;
    TFile *f= new TFile;
    f=TFile::Open("data/NewMixingDistributions2018.root","r");
    TH1D *hxi56 = (TH1D *)f->Get("xi_multi_S56");
    f->Close();
    delete f;  
    double xi56=hxi56->GetRandom();
    phys::Proton PUP56 = phys::Proton(xi56,false);
    backeventprotons.push_back(PUP56);
    delete hxi56;
  }

  if(backeventprotons.size()&& !(forweventprotons.size())){
    zetamarker=true;
    TFile *f= new TFile;
    f=TFile::Open("data/NewMixingDistributions2018.root","r");
    TH1D *hxi45 = (TH1D *)f->Get("xi_multi_S45");
    f->Close();
    delete f;  
    double xi45=hxi45->GetRandom();
    phys::Proton PUP45 = phys::Proton(xi45,true);
    forweventprotons.push_back(PUP45);
    delete hxi45;
  }

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
  
  //if(forweventprotons.size()>2) forweventprotons.clear();
  //if(backeventprotons.size()>2) backeventprotons.clear();
  
  phys::ProtonPair pp;
  phys::Proton p1,p2;
  
  int nmatches=0;
  bool matched=false;
  bool matcheddelta=false;
  double bestmassdif=9999.;
  double bestydif=9999.;
  
  if(forweventprotons.size()&&backeventprotons.size()){
    for(unsigned int i=0; i<backeventprotons.size();i++){
      p1=backeventprotons[i];
      //if(p1.LHCSector()) continue;
      for(unsigned int j=0; j<forweventprotons.size();j++){
	p2=forweventprotons[j];
        //if(!p2.LHCSector()) continue;
	pp=phys::ProtonPair(p1,p2);
        theHistograms->fill("mpp","Expected mass of ZZ pair",50,400,2400,pp.mpp(),theWeight);
        theHistograms->fill("ypp","Expected rapidity of ZZ pair",15,0,0.75,pp.ypp(),theWeight);

        if(ZZ->passFullSelection()){
	  double a = acoplanarity(ZZ);
	  theHistograms->fill("acoplanarity","Acoplanarity distribution",25,0,1,a,theWeight);
	  double massdif = 1-ZZ->mass()/pp.mpp();
	  double ydif = pp.ypp()-ZZ->rapidity();
	  nmatches++;
	  theHistograms->fill("th2","2D matching distribution",60,-2.5,0.5,60,-1.5,1.5,massdif,ydif,theWeight);
	  bool region = regionselection(massdif,ydif,kmin,kmax);
	  if(massdif*massdif+ydif*ydif<bestmassdif*bestmassdif+bestydif*bestydif && !matcheddelta){
	    bestmassdif=massdif;
	    bestydif=ydif;}
	  if(region){
	    if(!matcheddelta && massdif>-0.05){
	      bestmassdif=massdif;
	      bestydif=ydif;
	      matcheddelta=true;
	      //theHistograms->fill("xitest","xitest",100,-0.1,0.1,p1.xi()-ZZ->expectedxi2(),theWeight);
	      //theHistograms->fill("xitest","xitest",100,-0.1,0.1,p2.xi()-ZZ->expectedxi1(),theWeight);
	    }
	    if(!matched){
	      matched=true;}		
	  }
	}
      }}
    theHistograms->fill("nmatches","Number of possible matches for event",10,-0.5,9.5,nmatches,theWeight);
    theHistograms->fill("th2good","2D matching distribution",60,-2.5,0.5,60,-1.5,1.5,bestmassdif,bestydif,theWeight);
    if(regionselection(bestmassdif,bestydif,kmin,kmax)){
      theHistograms->fill("th1good","1D matching distribution",20,-2.5,0.5,bestmassdif,theWeight);
      theHistograms->fill("counteromicron2","counteromicron",1,0,1,0.5,theWeight);
      theHistograms->fill("th2matched","2D matching distribution",60,-2.5,0.5,60,-1.5,1.5,bestmassdif,bestydif,theWeight);
      theHistograms->fill("goodmZZ","Mass of matched ZZ pairs",10,650,4500,ZZ->mass(),theWeight);
      theHistograms->fill("goodyZZ","Rapidity of matched ZZ pairs",15,-1.3,1.3,ZZ->rapidity(),theWeight);
    }
    if(matcheddelta&&acoplanarity(ZZ)<=acomax) {    
      if(!zetamarker) theHistograms->fill("counterdelta2","counterdelta",1,0,1,0.5,theWeight);
      else theHistograms->fill("counterzeta","counterzeta",1,0,1,0.5,theWeight);
      if(abs(ZZ->first().daughter(0).id())==13&&abs(ZZ->second().daughter(0).id())==13) theHistograms->fill("counterdelta4mu","counterdelta4mu",1,0,1,0.5,theWeight);
      if(abs(ZZ->first().daughter(0).id())==11&&abs(ZZ->second().daughter(0).id())==11) theHistograms->fill("counterdelta4e","counterdelta4e",1,0,1,0.5,theWeight);
      if(abs(ZZ->first().daughter(0).id())==13&&abs(ZZ->second().daughter(0).id())==11) theHistograms->fill("counterdelta2e2mu","counterdelta2e2mu",1,0,1,0.5,theWeight);
      if(abs(ZZ->first().daughter(0).id())==11&&abs(ZZ->second().daughter(0).id())==13) theHistograms->fill("counterdelta2e2mu","counterdelta2e2mu",1,0,1,0.5,theWeight);
      if(!zetamarker) theHistograms->fill("th2goodC","2D matching distribution",70,-2.5,1,60,-1.5,1.5,bestmassdif,bestydif,theWeight);
      if(zetamarker) theHistograms->fill("th2goodCZeta","2D matching distribution",70,-2.5,1,60,-1.5,1.5,bestmassdif,bestydif,theWeight);
    }
  }
  
  //Algorithm B and C: use the highest xi proton every time
  
  if(forweventprotons.size()&&backeventprotons.size()){
    phys::ProtonPair ppxi= phys::ProtonPair(backproton,forwproton);
    double massdifxi = 1-ZZ->mass()/ppxi.mpp();
    double ydifxi = ppxi.ypp()-ZZ->rapidity();
    //theHistograms->fill("th2goodB","2D matching distribution",60,-2.5,0.5,60,-1.5,1.5,massdifxi,ydifxi,theWeight);
    bool regionxi = regionselection(massdifxi,ydifxi,kmin,kmax);
    if(!matcheddelta&&acoplanarity(ZZ)<=acomax){ if(!zetamarker) theHistograms->fill("th2goodC","2D matching distribution",70,-2.5,1,60,-1.5,1.5,massdifxi,ydifxi,theWeight);
      else theHistograms->fill("th2goodCZeta","2D matching distribution",70,-2.5,1,60,-1.5,1.5,bestmassdif,bestydif,theWeight);}
    if(regionxi&&ZZ->passFullSelection()){
      theHistograms->fill("counteromicronxi","counteromicronxi",1,0,1,0.5,theWeight);
      if(massdifxi>-0.05) {
	theHistograms->fill("counterdeltaxi","counterdeltaxi",1,0,1,0.5,theWeight);
      }
      if(!matcheddelta&&acoplanarity(ZZ)<=acomax){
	if(!zetamarker) theHistograms->fill("counteromicronC","counteromicronC",1,0,1,0.5,theWeight);
	if(abs(ZZ->first().daughter(0).id())==13&&abs(ZZ->second().daughter(0).id())==13) theHistograms->fill("counteromicron4mu","counteromicron4mu",1,0,1,0.5,theWeight);
        if(abs(ZZ->first().daughter(0).id())==11&&abs(ZZ->second().daughter(0).id())==11) theHistograms->fill("counteromicron4e","counteromicron4e",1,0,1,0.5,theWeight);
        if(abs(ZZ->first().daughter(0).id())==13&&abs(ZZ->second().daughter(0).id())==11) theHistograms->fill("counteromicron2e2mu","counteromicron2e2mu",1,0,1,0.5,theWeight);
        if(abs(ZZ->first().daughter(0).id())==11&&abs(ZZ->second().daughter(0).id())==13) theHistograms->fill("counteromicron2e2mu","counteromicron2e2mu",1,0,1,0.5,theWeight);
      }
      if(zetamarker&&acoplanarity(ZZ)<=acomax) theHistograms->fill("counterzeta","counterzeta",1,0,1,0.5,theWeight);
    }
  }

  phys::Proton genproton;
  foreach(const phys::Particle particle, *genParticles){
    if(particle.id()==2212){
      phys::Proton genproton= phys::Proton(particle);
      theHistograms->fill("genprotonxi","gen proton xi",50,0,1,genproton.xi(),theWeight);
    }
  }

  
}



Double_t acoplanarity(phys::DiBoson <phys::Lepton,phys::Lepton> *ZZ){
  double_t aco = abs(1-abs(physmath::deltaPhi(ZZ->first().phi(),ZZ->second().phi()))/TMath::Pi());
  return aco;
}

//out=0 -> not signal
//out=1 -> signal

bool regionselection(double x, double y, double kmin, double kmax){
  
  bool out=false;
  TF1 *f1 = new TF1("func1",logfunc,-2,-0.05,1);
  f1->SetParameter(0,kmin);
  TF1 *f2 = new TF1("func2",logfunc,-2,-0.05,1);
  f2->SetParameter(0,kmax);
  if(x<-2||x>0.1) {return false;}
  else if(x<=-0.05) {
    return abs(y) >= f1->Eval(x) && abs(y) <= f2->Eval(x);
  }
  else if(x<=0&&x>-0.05) {return abs(y)<=x+0.1;}
  else if(x<=0.05&&x>0) {return abs(y)<=0.1;}
  else if(x>0.05&&x<=0.1) {return abs(y)<=-x+0.15;}
  return out;
}


double logfunc(double *x, double *par){
  double xx=x[0];
  return log(par[0]*(1-xx));
}
