#include "VVXAnalysis/TreeAnalysis/interface/WlllnuAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"

#include <TCanvas.h>
#include <TH1F.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp>
using namespace boost::assign;

using std::cout;
using std::endl;


using namespace phys;

void WlllnuAnalyzer::begin(){
  nevents = 0;
}

Int_t WlllnuAnalyzer::cut() {
  nevents++;
  return 1;
}


void WlllnuAnalyzer::analyze(){
  //  if(event != 7701979 && event != 7701973 && event != 7701836  && event != 7701827  && event != 7701817  && event !=  7701753  ) return; 
  cout << "------------------------------------------------------------------"<<endl;
  cout << "Run: " << run << " event: " << event << "  #: " << nevents << endl;
  
  //leptons
  
  std::vector<phys::Particle > leptons;
  std::vector<phys::Particle > neutrinos;
  std::vector<Particle> electrons; // Z daughters
  std::vector<Particle> muons; // Z daughters
  int finalid = 0;

   //find leptons
  foreach(const phys::Particle &gen, *genParticles){    
    if( (abs(gen.id()) != 11 && abs(gen.id()) != 13 && abs(gen.id()) != 12 && abs(gen.id()) != 14) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) continue;
    //if(abs(gen.id()) != 11 && abs(gen.id()) != 13) continue;
    finalid += abs(gen.id());
    //    cout << " genLepton: " << gen << endl;
    cout << "id: " << gen.id() << " pt: " << gen.pt() << " mass: " << gen.mass() << " eta: " << gen.eta() << endl;
    theHistograms.fill("ptAllGenParticle",   "pt ",   100, 0,   200,  gen.pt());
    theHistograms.fill("etaAllGenParticle",  "eta ",  100, -10, 10,   gen.eta());
    //theHistograms.fill("YAllGenParticle",    "Y ",    100, -10, 10,  gen.rapidity());
    bool isLepton = abs(gen.id()) == 11 || abs(gen.id()) == 13;
    isLepton ? leptons.push_back(gen) : neutrinos.push_back(gen);
    if (gen.id() == 11){
      electrons.insert(electrons.begin(),gen); //first all e-, then all e+
      //theHistograms.fill("ptElectrons","pt e", 100, 0, 200, gen.pt());
    }
    else if (gen.id() == -11){
      electrons.push_back(gen);
      //theHistograms.fill("ptElectrons","pt e", 100, 0, 200, gen.pt());
    }
    else if (gen.id() == 13){
      muons.insert(muons.begin(),gen); //first all mu-, then all mu+
      //theHistograms.fill("ptMuons","pt mu", 100, 0, 200, gen.pt());
    }
    else if (gen.id() == -13){
      muons.push_back(gen);
      //theHistograms.fill("ptMuons","pt mu", 100, 0, 200, gen.pt());
    }
  }
  theHistograms.fill("leptonsNumber",  "number of leptons ",  5, -0.5, 4.5, leptons.size());
  
  //ZZ
  if( (electrons.size() + muons.size() == 4) && (finalid == 44 || finalid == 48 || finalid == 52) ){
    cout << "\nZZ analysis" << endl; 
    //build z0, z1
    std::vector<Boson<Particle> > Zcandidates;
    phys::Boson<phys::Particle> z0;
    phys::Boson<phys::Particle> z1;
    
    if (finalid == 44){ //4e
      Zcandidates.push_back(phys::Boson<Particle>(electrons[0],electrons[2], 23));
      Zcandidates.push_back(phys::Boson<Particle>(electrons[1],electrons[3], 23));
      Zcandidates.push_back(phys::Boson<Particle>(electrons[0],electrons[3], 23));
      Zcandidates.push_back(phys::Boson<Particle>(electrons[1],electrons[2], 23));
    }
    else if (finalid == 52){ //4mu
      Zcandidates.push_back(phys::Boson<Particle>(muons[0],muons[2], 23));
      Zcandidates.push_back(phys::Boson<Particle>(muons[1],muons[3], 23));
      Zcandidates.push_back(phys::Boson<Particle>(muons[0],muons[3], 23));
      Zcandidates.push_back(phys::Boson<Particle>(muons[1],muons[2], 23));
    }
    else if (finalid == 48){ //2e2mu
      Zcandidates.push_back(phys::Boson<Particle>(electrons[0],electrons[1], 23));
      Zcandidates.push_back(phys::Boson<Particle>(muons[0],muons[1], 23));
    }

    //making z0, z1 with mass criteria (mass closer to ZMASS)
    std::stable_sort(Zcandidates.begin(), Zcandidates.end(), phys::MassComparator(ZMASS));
    z0 = Zcandidates[0];
    foreach(const phys::Boson<phys::Particle> cand, Zcandidates)
      if ( cand.daughter(0).pt() != z0.daughter(0).pt() &&  cand.daughter(1).pt() != z0.daughter(1).pt() ) z1 = cand;
    //!!!need to find a way to check if 2 particles are the same!!!!
    theHistograms.fill("massGenZ0", "z0 mass ", 1000, 50, 150, z0.mass());
    theHistograms.fill("massGenZ1", "z1 mass ", 1000, 50, 150, z1.mass());
    theHistograms.fill("massGenZParticles", "Z mass ", 1000, 50, 150, z0.mass());
    theHistograms.fill("massGenZParticles", "Z mass ", 1000, 50, 150, z1.mass());

    cout << "\nz0 " "\tmass: " << z0.mass() << "\tmTdaughters-mZ0: "<< mT(z0.daughter(0), z0.daughter(1)) - z0.mass() << "\n\tz0.daughter(0) \tid(): " << z0.daughter(0).id() <<  " pt(): " << z0.daughter(0).pt() << "\n\tz0.daughter(1)\tid(): " << z0.daughter(1).id() << " pt(): " << z0.daughter(1).pt() <<  endl;
    cout << "\nz1 " "\tmass: " << z1.mass() << "\tmTdaughters-mZ1: "<< mT(z1.daughter(0), z1.daughter(1)) - z1.mass() << "\n\tz1.daughter(0) \tid(): " << z1.daughter(0).id() <<  " pt(): " << z1.daughter(0).pt() << "\n\tz1.daughter(1)\tid(): " << z1.daughter(1).id() << " pt(): " << z1.daughter(1).pt() <<  endl;
    //comparing mZ with mT daughters
    theHistograms.fill("mZmT", "mZ vs mT daughters", 600, 0, 150, 600, 0, 150, mT(z0.daughter(0),z0.daughter(1)), z0.mass());
    theHistograms.fill("mZmT", "mZ vs mT daughters", 600, 0, 150, 600, 0, 150, mT(z1.daughter(0),z1.daughter(1)), z1.mass());
    theHistograms.fill("deltamZmT", "mZ - mT daughters", 100, -100, 100, z0.mass() - mT(z0.daughter(0),z0.daughter(1)));
    theHistograms.fill("deltamZmT", "mZ - mT daughters", 100, -100, 100, z1.mass() - mT(z1.daughter(0),z1.daughter(1)));
    
    
    //making z0,z1 with mT criteria (mT daughters closer to ZMASS)
    std::stable_sort(Zcandidates.begin(), Zcandidates.end(), mTComparator(ZMASS));
    phys::Boson<phys::Particle> z0mT = Zcandidates[0];
    phys::Boson<phys::Particle> z1mT;
    foreach(const phys::Boson<phys::Particle> cand, Zcandidates)
      if ( cand.daughter(0).pt() != z0mT.daughter(0).pt() &&  cand.daughter(1).pt() != z0mT.daughter(1).pt() ) z1mT = cand;
    theHistograms.fill("massGenZmT", "Z (selected by mT) mass ", 1000, 50, 150, z0mT.mass());
    theHistograms.fill("massGenZmT", "Z (selected by mT) mass ", 1000, 50, 150, z1mT.mass());
    
    //comparing z0mT, z1mT to z0, z1
    int sameZandZmT = ((isTheSame(z0, z0mT) && isTheSame(z1,z1mT)) || (isTheSame(z0, z1mT) && isTheSame(z1,z0mT))) ? 1 : 0;
    theHistograms.fill("sameGenZ", "Z and ZmT are the same? ", 2 , -0.5, 1.5, sameZandZmT);
    int sameZ0Z1 = (isTheSame(z0, z0mT) && isTheSame(z1,z1mT)) ? 1 : 0;
    theHistograms.fill("sameGenZ0z1", "Z and ZmT are the same? (also the order) ", 2 , -0.5, 1.5, sameZ0Z1);
        
    /*   //deltaR, deltaEta, deltaPhi
    double deltaEta0 = z0.daughter(0).eta() - z1.eta(); //between z1 and z0's first daughter
    double deltaEta1 = z0.daughter(1).eta() - z1.eta(); //between z1 and z0's first daughter
    // some print out  
    cout << "\n\ndeltaR daughter 0: " << deltaR(z1.p4(), z0.daughter(0).p4()) << endl;
    cout << "deltaR daughter 1: "     << deltaR(z1.p4(), z0.daughter(1).p4()) << endl; 
    cout << "deltaPhi daughter 0: "   << deltaPhi(z0.daughter(0).p4(), z1.p4()) << endl;
    cout << "deltaPhi daughter 1: "   << deltaPhi(z0.daughter(1).p4(), z1.p4()) << endl;
    cout << "deltaEta daughter 0: "   << deltaEta0 << endl;
    cout << "deltaEta daughter 1: "   << deltaEta1 << endl;
  */
    //how do we know z0 is the first in the diagram?? I know
    
    // ZZ
    ZZtype ZZ(z0,z1);
    cout << "\n\n ZZ: " << ZZ.id() << " pt: " << ZZ.pt() << " mass: " << ZZ.mass() << "\n daughters: " << ZZ.daughter<Particle>(0).id() << ", " << ZZ.second().id() << " Y: " << ZZ.rapidity() << endl;
    
    theHistograms.fill("ptGenZZ",   "pt ",   100,  0,   100, ZZ.pt());
    //    theHistograms.fill("etaGenZZ",  "eta ",  100,  -10, 10,  ZZ.eta());
    theHistograms.fill("massGenZZ", "mass ", 100,  0,   400, ZZ.mass());
    theHistograms.fill("YGenZZ",    "Y ",    100,  -10, 10,  ZZ.rapidity());
    
    return;
  }
    
  //  cout << " \n\t\t  -----------------------------------\n " << endl;
 
  else if(leptons.size()<3){
    cout << "\nLess than 3 leptons" << endl;
    return;
  }
  else if(neutrinos.size()>1){
    cout << "\nMore than 1 nu" << endl;
    return;
  }
  
  //ZL
  else if(leptons.size()==3 && neutrinos.size() == 0){// 3l(finalid == 33 || finalid == 35 || finalid == 37 || finalid == 39)){ //3l

    //----RECO particles----//
        
    cout << "\nZL analysis" << endl; 
    cout << "\n # of ZL candidates: "  << (*ZL).size() << endl; //reco ZL
    foreach(const ZLCompositeCandidate &zl, *ZL){
      cout << "\nZlcand: \t" << std::get<0>(zl) << "\n\t\t" << std::get<1>(zl) << endl;
      theHistograms.fill("massBosonRecoZl","mass Z", 100, 0, 500, (std::get<0>(zl)).mass());
    }
    if((*ZL).size() > 0) theHistograms.fill("idParticlesRecoZL"," leptons & daughters id in ZL", 10 , 30.5, 40.5, finalid);
        
    //----my particles----//
    cout << "\nmyZl analysis" << endl; 
    //forms Zl every time there are 3 leptons (1 or 2 Zl depending on presence of both e and mu) if leptons pt and eta is good enough
    
    std::vector<Zltype > Zl; //change Particle to Lepton    
    std::vector<Particle> eptSort; // electrons sorted by pt
    std::vector<Particle> muptSort;
    std::vector<Particle> lepptSort;
    
    //check pt(e)>7Gev, pt(mu)>5Gev, abs(eta(e))<2.5, abs(eta(mu))<2.4
    foreach(const Particle e, electrons){
      if(e.pt()<7 || abs(e.eta())>2.5 ){
	cout<<"pt leptons not sufficient (e 7Gev) or eta over limit (|eta| 2.5)"<< endl;
	return;
      }
      eptSort.push_back(e);
      lepptSort.push_back(e);
    }
    foreach(const Particle mu, muons){
      if(mu.pt()<5 || abs(mu.eta())>2.4 ){
	cout<<"pt leptons not sufficient (mu 5Gev) or eta over limit (|eta| 2.5)"<< endl;
	return;
      }
      muptSort.push_back(mu);
      lepptSort.push_back(mu);
    }
    std::stable_sort(eptSort.begin(), eptSort.end(), PtComparator());
    std::stable_sort(muptSort.begin(), muptSort.end(), PtComparator());
    std::stable_sort(lepptSort.begin(), lepptSort.end(), PtComparator());
        
    //check 1st lepton has pt>20 and 2nd pt>10(mu)/12(e)
    if(lepptSort[0].pt()<20 || (abs(lepptSort[1].id())==11 && lepptSort[1].pt()<12) || (abs(lepptSort[1].id())==13 && lepptSort[1].pt()<10)){
      cout<<"pt leptons not sufficient (first 20Gev, second 10Gev(mu)/12Gev(e))"<< endl;
      return;
    }
    
    if(finalid == 33){//3e     
      Zl.push_back(Zltype(phys::Boson<phys::Particle>(electrons[0],electrons[2], 23), electrons[1]));
      electrons[1].id()<0 ? Zl.push_back(Zltype(phys::Boson<phys::Particle>(electrons[0],electrons[1], 23), electrons[2])) : Zl.push_back(Zltype(phys::Boson<phys::Particle>(electrons[1],electrons[2], 23), electrons[0]));
    }
    else if(finalid == 39){//3mu
      Zl.push_back(Zltype(phys::Boson<phys::Particle>(muons[0],muons[2], 23), muons[1]));
      muons[1].id()<0 ? Zl.push_back(Zltype(phys::Boson<phys::Particle>(muons[0],muons[1], 23), muons[2])) : Zl.push_back(Zltype(phys::Boson<phys::Particle>(muons[1],muons[2], 23), muons[0]));	
    }
    else if(finalid == 35){//2e1mu
      Zl.push_back(Zltype(phys::Boson<phys::Particle>(electrons[0],electrons[1], 23), muons[0]));
    }
    else if(finalid == 37){//2mu1e
      Zl.push_back(Zltype(phys::Boson<phys::Particle>(muons[0],muons[1], 23), electrons[0]));
    }
    
    cout << "\n # of my Zl candidates: "  << Zl.size() << endl;
    foreach(const Zltype zl, Zl){
      cout << "\nZlcand: \t" << std::get<0>(zl) << "\n\t\t" << std::get<1>(zl) << endl;
      theHistograms.fill("massBosonGenZl","mass Z", 100, 0, 500, (std::get<0>(zl)).mass());
         
      //cout << "\nZl good cand (right deltaEta and pt)\t" << " leptons forming Z: " << abs(zl.first.daughter(0).id()) << " other lepton: " << abs(zl.second.id()) << endl;
    }
    if (Zl.size() > 0) theHistograms.fill("idParticlesGenZl"," leptons & daughters id in Zl", 10 , 30.5, 40.5, finalid);

    //--------------------------------------------------//
    //comparing mZ with mT daughters
    foreach(const Zltype zl, Zl){
      if(Zl.size() == 0) continue;
      theHistograms.fill("mZlmT", "mZ vs mT daughters", 150, 0, 150, 150, 0, 150, mT((std::get<0>(zl)).daughter(0),(std::get<0>(zl)).daughter(1)), (std::get<0>(zl)).mass());
      theHistograms.fill("deltamZlmT", "mZ - mT daughters", 100, -100, 100, (std::get<0>(zl)).mass() - mT((std::get<0>(zl)).daughter(0),(std::get<0>(zl)).daughter(1)));
    }
    
    return;
  }

  //Wlllnu
  else if(leptons.size() == 3 && neutrinos.size() == 1){//3l1nu (finalid == 45 || finalid == 49 || finalid == 53)){
    //  Particle nu = neutrinos[0];
    TLorentzVector Ptot = neutrinos[0].p4();
    foreach(const Particle lep, leptons)
      Ptot += lep.p4();
    double masslllnu = Ptot.M();
    theHistograms.fill("massGenlllnu","mass lllnu", 100, 0, 600, masslllnu);
    //cout << "\n masslllnu: " << masslllnu << endl;
    
    return;
  }


}
