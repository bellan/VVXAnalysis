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

Int_t WlllnuAnalyzer::cut() {
  return 1;
}


template<typename T> double mT(const T& p1, const T& p2){

  return sqrt( 2*p1.pt()*p2.pt()*(1-TMath::Cos(physmath::deltaPhi(p1.phi(), p2.phi()))) );

}


void WlllnuAnalyzer::analyze(){
  
  cout << "------------------------------------------------------------------"<<endl;
  cout << "Run: " << run << " event: " << event << endl;
  
  //leptons
  
  std::vector<phys::Particle>  leptons;
  std::vector<Particle> electrons; // Z daughters
  std::vector<Particle> muons; // Z daughters  
  int finalid = 0;

   //find leptons
  foreach(const phys::Particle &gen, *genParticles){    
    if( (abs(gen.id()) != 11 && abs(gen.id()) != 13) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) continue;
    //if(abs(gen.id()) != 11 && abs(gen.id()) != 13) continue;
    finalid += abs(gen.id());
    //cout << " genLepton: " << gen << endl;
    cout << "id: " << gen.id() << " pt: " << gen.pt() << " mass: " << gen.mass() << " eta: " << gen.eta() << endl;
    theHistograms.fill("ptAllGenParticle",   "pt ",   100, 0,   200,  gen.pt());
    theHistograms.fill("etaAllGenParticle",  "eta ",  100, -10, 10,   gen.eta());
    theHistograms.fill("massAllGenParticle", "mass ", 100, 0,   0.12, gen.mass());
    theHistograms.fill("YAllGenParticle",    "Y ",    100, 0,   100,  gen.rapidity());
    leptons.push_back(gen);
    if (gen.id() == 11){
      electrons.insert(electrons.begin(),gen); //fist all e-, then all e+
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
    std::stable_sort(Zcandidates.begin(), Zcandidates.end(), phys::MassComparator(ZMASS));
    z0 = Zcandidates[0];
    //foreach(const phys::Boson<phys::Particle> cand, Zcandidates)
    //if ( abs(cand.mass() - ZMASS) < abs(z0.mass() - ZMASS)) z0 = cand;
    foreach(const phys::Boson<phys::Particle> cand, Zcandidates)
      if ( cand.daughter(0).pt() != z0.daughter(0).pt() &&  cand.daughter(1).pt() != z0.daughter(1).pt() ) z1 = cand;
    //!!!need to find a way to check if 2 particles are the same!!!!
    
    theHistograms.fill("massZParticles", "Z mass ", 1000, 50, 150, z0.mass());
    theHistograms.fill("massZparticles", "Z mass ", 1000, 50, 150, z1.mass());

    //some printout   
    //cout << "\n z0: " << z0 << endl;
    cout << "\nz0: " << z0.id() << " mass: " << z0.mass() << "\tz0.daughter(0).id(): " << z0.daughter(0).id() << " z0.daughter(1).id(): " << z0.daughter(1).id() << endl;
    //cout << " z1: " << z1 << endl;
    cout << "\nz1: " << z1.id() << " mass: " << z1.mass() << "\tz1.daughter(0).id(): " << z1.daughter(0).id() << " z1.daughter(1).id(): " << z1.daughter(1).id() << endl;
      
    //comparing mT leptons with z0, z1 (a bit messy!!)
      
    cout << "\nComparing mT leptons and mass z0/z1" << endl;
    for(unsigned int i=0; i<leptons.size(); i++){
      for(unsigned int j=i+1; j<leptons.size(); j++){
	if(leptons[j].id() == -leptons[i].id()){
	  if(abs(leptons[i].id()) == abs(z0.daughter(0).id()) )
	    cout << "mT(leptons " << i << ", " << j << ") - mT(z0): " << abs(mT(leptons[i], leptons[j]) - z0.mass()) << endl;
	  if(abs(leptons[i].id()) == abs(z1.daughter(0).id()) )
	    cout << "mT(leptons " << i << ", " << j << ") - mT(z1): " << abs(mT(leptons[i], leptons[j]) - z1.mass()) << endl;
	}
	else cout << "leptons " << i << ", " << j << ": no Z daughters candididates" << endl;
      }
    }
    
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
    
    // ZZ
    DiBoson<phys::Particle,phys::Particle>  ZZ(z0,z1);
    cout << "\n\n ZZ: " << ZZ.id() << " pt: " << ZZ.pt() << " mass: " << ZZ.mass() << "\n daughters: " << ZZ.daughter<Particle>(0).id() << ", " << ZZ.second().id() << " Y: " << ZZ.rapidity() << endl; //why daughter(0) instead of first() is not working?? It's templated!! Try daughter<Particle>(0)
    
    theHistograms.fill("ptZZ",   "pt ",   100,  0,   100, ZZ.pt());
    theHistograms.fill("etaZZ",  "eta ",  100,  -10, 10,  ZZ.eta());
    theHistograms.fill("massZZ", "mass ", 1000, 0,   400, ZZ.pt());
    theHistograms.fill("YZZ",    "Y ",    100,  0,   100, ZZ.rapidity());
    //<phys::Boson<phys::Particle>,phys::Boson<phys::Particle> >

    return;
  }
    
  //  cout << " \n\t\t  -----------------------------------\n " << endl;
  
  else if(leptons.size()<3){
    cout << "\nLess than 3 leptons" << endl;
    return;
  }
  
  //ZL
  else if(leptons.size()==3){

    //----RECO particles----//
    // std::pair<Boson<phys::Lepton> ,phys::Lepton>
    
    cout << "\nZL analysis" << endl; 
    cout << "\n # of ZL candidates: "  << (*ZL).size() << endl;
    foreach(const ZLCompositeCandidate &zl, *ZL){
      cout << "\nZlcand: \t" << std::get<0>(zl) << "\n\t\t" << std::get<1>(zl) << endl;
      //theHistograms.fill("massBosonZl","mass Z", 100, 0, 500, (std::get<0>(zl)).mass());
      //theHistograms.fill("massLeptonZl","mass l", 100, 0, 0.12, (std::get<1>(zl)).mass());
      theHistograms.fill("idDaughterZl"," Z daughters id", 5 , 9.5, 14.5, abs((zl.first).daughter(0).id()));
      theHistograms.fill("idLeptonZl"," leptons id", 5 , 9.5, 14.5, abs((zl.second).id()));

    }
    if ((*ZL).size() > 0) theHistograms.fill("idAllParticlesZL"," leptons & daughters id in ZL", 10 , 30.5, 40.5, finalid);
    //  if(ZL.size()>0) theHistograms.fill("daughtersZl"," Z daughters id", 100, 0, 500, abs((ZL[0].first).daughter(0).id()));
    // if(ZL.size()>1) theHistograms.fill("daughtersZl"," Z daughters id", 100, 0, 500, abs((ZL[1].first).daughter(0).id()));
    
    //----my particles----//
    cout << "\nmyZl analysis" << endl; 
    //forms Zl every time there are 3 leptons (1 or 2 Zl depending on presence of both e and mu) if leptons pt is good enough
    
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
        
    //check 1st lepton has pt>20 and 2nd pt>10
    if(lepptSort[0].pt()<20 || lepptSort[1].pt()<10){
      cout<<"pt leptons not sufficient (first 20Gev, second 10Gev)"<< endl;
      return;
    }

    if(finalid == 33){//3e     
      Zl.push_back(Zltype(phys::Boson<phys::Particle>(electrons[0],electrons[2], 23), electrons[1]));
      if(electrons[1].id()<0)
	Zl.push_back(Zltype(phys::Boson<phys::Particle>(electrons[0],electrons[1], 23), electrons[2]));
      else
	Zl.push_back(Zltype(phys::Boson<phys::Particle>(electrons[1],electrons[2], 23), electrons[0]));
    }
    else if(finalid == 39){//3mu
      Zl.push_back(Zltype(phys::Boson<phys::Particle>(muons[0],muons[2], 23), muons[1]));
      if(muons[1].id()<0)
	Zl.push_back(Zltype(phys::Boson<phys::Particle>(muons[0],muons[1], 23), muons[2]));
      else
	Zl.push_back(Zltype(phys::Boson<phys::Particle>(muons[1],muons[2], 23), muons[0]));	
    }
    else if(finalid == 35){//2e1mu
      Zl.push_back(Zltype(phys::Boson<phys::Particle>(electrons[0],electrons[1], 23), muons[0]));
    }
    else if(finalid == 37){//2mu1e
      Zl.push_back(Zltype(phys::Boson<phys::Particle>(muons[0],muons[1], 23), electrons[0]));
    }
    
    cout << "\n # of my Zl candidates: "  << Zl.size() << endl;
    foreach(const Zltype zl, Zl){
      //cout << "\nZlcand: \t" << std::get<0>(zl) << "\n\t\t" << std::get<1>(zl) << endl;
      //double deltaEtaZl =zl.first.eta()-zl.second.eta();
      //theHistograms.fill("massBosonMyZl","mass Z", 100, 0, 500, (std::get<0>(zl)).mass());
      //theHistograms.fill("massLeptonMyZl","mass l", 100, 0, 0.12, (std::get<1>(zl)).mass());
      theHistograms.fill("idDaughterMyZl", " Z daughters id", 5,   9.5, 14.5, abs((zl.first).daughter(0).id()));
      theHistograms.fill("idLeptonMyZl",   " leptons id",     5,   9.5, 14.5, abs((zl.second).id()));
      //theHistograms.fill("deltaEtaMyZl",   " delta Eta ",     100, -10, 10,   deltaEtaZl );
      
      cout << "\nZl good cand (right deltaEta and pt)\t" << " leptons forming Z: " << abs(zl.first.daughter(0).id()) << " other lepton: " << abs(zl.second.id()) << endl;
    }
    if (Zl.size() > 0) theHistograms.fill("idAllParticlesMyZl"," leptons & daughters id in Zl", 10 , 30.5, 40.5, finalid);

    //some histograms
    if(eptSort.size()  > 0) theHistograms.fill("pte0",  "pt e0",  100, 0, 200, eptSort[0].pt());
    if(eptSort.size()  > 1) theHistograms.fill("pte1",  "pt e1",  100, 0, 200, eptSort[1].pt());
    if(eptSort.size()  > 2) theHistograms.fill("pte2",  "pt e2",  100, 0, 200, eptSort[2].pt());
    if(muptSort.size() > 0) theHistograms.fill("ptmu0", "pt mu0", 100, 0, 200, muptSort[0].pt());
    if(muptSort.size() > 1) theHistograms.fill("ptmu1", "pt mu1", 100, 0, 200, muptSort[1].pt());
    if(muptSort.size() > 2) theHistograms.fill("ptmu2", "pt mu2", 100, 0, 200, muptSort[2].pt());

    
    return;
 
  }
  
}

//there are still some events where I form Zl but it has no ZL!! (the ones with 3e)
