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
    if( (abs(gen.id()) != 11 && abs(gen.id()) != 13 && (abs(gen.id()) != 12) && (abs(gen.id()) != 14)) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) continue;
    //if(abs(gen.id()) != 11 && abs(gen.id()) != 13) continue;
    finalid += abs(gen.id());
    //    cout << " genLepton: " << gen << endl;
    cout << "id: " << gen.id() << " pt: " << gen.pt() << " mass: " << gen.mass() << " eta: " << gen.eta() << endl;
    theHistograms.fill("ptAllGenParticle",   "pt ",   100, 0,   200,  gen.pt());
    theHistograms.fill("etaAllGenParticle",  "eta ",  100, -10, 10,   gen.eta());
    theHistograms.fill("massAllGenParticle", "mass ", 100, 0,   0.12, gen.mass());
    theHistograms.fill("YAllGenParticle",    "Y ",    100, 0,   100,  gen.rapidity());
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
    //foreach(const phys::Boson<phys::Particle> cand, Zcandidates)
    //if (abs(cand.mass() - ZMASS) < abs(z0.mass() - ZMASS)) z0 = cand;
    foreach(const phys::Boson<phys::Particle> cand, Zcandidates)
      if ( cand.daughter(0).pt() != z0.daughter(0).pt() &&  cand.daughter(1).pt() != z0.daughter(1).pt() ) z1 = cand;
    //!!!need to find a way to check if 2 particles are the same!!!!
    theHistograms.fill("massGenZ0", "z0 mass ", 1000, 50, 150, z0.mass());
    theHistograms.fill("massGenZ1", "z1 mass ", 1000, 50, 150, z1.mass());
    theHistograms.fill("massGenZParticles", "Z mass ", 1000, 50, 150, z0.mass());
    theHistograms.fill("massGenZParticles", "Z mass ", 1000, 50, 150, z1.mass());

    //cout << "\n z0: " << z0 << endl;
    cout << "\nz0: " << z0.id() << " mass: " << z0.mass() << "\tz0.daughter(0).id(): " << z0.daughter(0).id() << " z0.daughter(1).id(): " << z0.daughter(1).id() << "\tz0.daughter(0).pt(): " << z0.daughter(0).pt() << " z0.daughter(1).pt(): " << z0.daughter(1).pt() << "\n  mTdaughters-mZ0: "<< mT(z0.daughter(0), z0.daughter(1)) - z0.mass() << endl;
    //cout << " z1: " << z1 << endl;
    cout << "\nz1: " << z1.id() << " mass: " << z1.mass() << "\tz1.daughter(0).id(): " << z1.daughter(0).id() << " z1.daughter(1).id(): " << z1.daughter(1).id() << "\tz1.daughter(0).pt(): " << z1.daughter(0).pt() << " z1.daughter(1).pt(): " << z1.daughter(1).pt() << "\n  mTdaughters-mZ1: "<< mT(z1.daughter(0), z1.daughter(1)) - z1.mass() << endl;

    //comparing mZ with mT daughters
    theHistograms.fill("mZmT", "mZ vs mT daughters", 300, 0, 150, 300, 0, 150, mT(z0.daughter(0),z0.daughter(1)), z0.mass());
    theHistograms.fill("mZmT", "mZ vs mT daughters", 300, 0, 150, 300, 0, 150, mT(z1.daughter(0),z1.daughter(1)), z1.mass());
    
    /*
    //making z0,z1 with mT criteria (mT daughters closer to ZMASS)
    std::stable_sort(Zcandidates.begin(), Zcandidates.end(), mTComparator(ZMASS));
    phys::Boson<phys::Particle> z0mT = Zcandidates[0];
    phys::Boson<phys::Particle> z1mT;
    foreach(const phys::Boson<phys::Particle> cand, Zcandidates)
      if ( cand.daughter(0).pt() != z0mT.daughter(0).pt() &&  cand.daughter(1).pt() != z0mT.daughter(1).pt() ) z1mT = cand;
    theHistograms.fill("massGenZmT", "Z (selected by mT) mass ", 1000, 50, 150, z0mT.mass());
    theHistograms.fill("massGenZmT", "Z (selected by mT) mass ", 1000, 50, 150, z1mT.mass());
    
    cout << "\nz0mT: " << z0mT.id() << " mass: " << z0mT.mass() << "\tz0mT.daughter(0).id(): " << z0mT.daughter(0).id() << " z0mT.daughter(1).id(): " << z0mT.daughter(1).id() << "\tz0mT.daughter(0).pt(): " << z0mT.daughter(0).pt() << " z0mT.daughter(1).pt(): " << z0mT.daughter(1).pt() << endl;
    cout << "\nz1mT: " << z1mT.id() << " mass: " << z1mT.mass() << "\tz1mT.daughter(0).id(): " << z1mT.daughter(0).id() << " z1mT.daughter(1).id(): " << z1mT.daughter(1).id() << "\tz1mT.daughter(0).pt(): " << z1mT.daughter(0).pt() << " z1mT.daughter(1).pt(): " << z1mT.daughter(1).pt() << endl;

    //comparing z0mT, z1mT to z0, z1
    int sameZandZmT = ((isTheSame(z0, z0mT) && isTheSame(z1,z1mT)) || (isTheSame(z0, z1mT) && isTheSame(z1,z0mT))) ? 1 : 0;
    theHistograms.fill("sameGenZ", "Z and ZmT are the same? ", 2 , -0.5, 1.5, sameZandZmT);
       */
    /*
    //comparing mT leptons with z0, z1 (a bit messy!!)      
    cout << "\nComparing mT leptons and mass z0/z1" << endl;
    for(unsigned int i=0; i<leptons.size(); i++){
      for(unsigned int j=i+1; j<leptons.size(); j++){
	if(!(leptons[j].id() == -leptons[i].id())){//no opposite id
	  //cout << "leptons " << i << ", " << j << ": no opposite id!" << endl;
	  continue;
	}
	if(abs(leptons[i].id()) == abs(z0.daughter(0).id()) )
	  cout << "mT(leptons " << i << ", " << j << ") - m(z0): " << abs(mT(leptons[i], leptons[j]) - z0.mass()) << /* "\tid: " << leptons[i].id() << " " << leptons[j].id() << endl;
	if(abs(leptons[i].id()) == abs(z1.daughter(0).id()) )
	  cout << "mT(leptons " << i << ", " << j << ") - m(z1): " << abs(mT(leptons[i], leptons[j]) - z1.mass()) <</* "\tid: " << leptons[i].id() << " " << leptons[j].id() << endl;
      }
    }
    */
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
    //how do we know z0 is the first in the diagram?? I assumed it, but I'm not sure!
    
    // ZZ
    //DiBoson<phys::Particle,phys::Particle>  ZZ(z0,z1);
    ZZtype ZZ(z0,z1);
    cout << "\n\n ZZ: " << ZZ.id() << " pt: " << ZZ.pt() << " mass: " << ZZ.mass() << "\n daughters: " << ZZ.daughter<Particle>(0).id() << ", " << ZZ.second().id() << " Y: " << ZZ.rapidity() << endl;
    
    theHistograms.fill("ptGenZZ",   "pt ",   100,  0,   100, ZZ.pt());
    theHistograms.fill("etaGenZZ",  "eta ",  100,  -10, 10,  ZZ.eta());
    theHistograms.fill("massGenZZ", "mass ", 1000, 0,   400, ZZ.mass());
    theHistograms.fill("YGenZZ",    "Y ",    100,  0,   100, ZZ.rapidity());
    
    return;
  }
    
  //  cout << " \n\t\t  -----------------------------------\n " << endl;
 
  else if(leptons.size()<3){
    cout << "\nLess than 3 leptons" << endl;
    return;
  }
  else if(neutrinos.size()<1){
    cout << "\nMore than 1 nu" << endl;
    return;
  }
  
  //ZL
  else if(leptons.size()==3 && neutrinos.size() == 0){// 3l(finalid == 33 || finalid == 35 || finalid == 37 || finalid == 39)){ //3l

    //----RECO particles----//
        
    cout << "\nZL analysis" << endl; 
    cout << "\n # of ZL candidates: "  << (*ZL).size() << endl; //reco ZL
    //there're some negative or "strange" electrons' mass values: why?? eg  event: 3553456  #: 555,  event: 1542753  #: 502,  event: 7702004  #: 1973 -> depends on p4, in any case so small (mass) that doesn't matter
    foreach(const ZLCompositeCandidate &zl, *ZL){
      cout << "\nZlcand: \t" << std::get<0>(zl) << "\n\t\t" << std::get<1>(zl) << endl;
      theHistograms.fill("massBosonRecoZl","mass Z", 100, 0, 500, (std::get<0>(zl)).mass());
      //theHistograms.fill("massLeptonZl","mass l", 100, 0, 0.12, (std::get<1>(zl)).mass());
      //theHistograms.fill("idDaughterRecoZl"," Z daughters id", 5 , 9.5, 14.5, abs((zl.first).daughter(0).id()));
      //theHistograms.fill("idLeptonRecoZl"," leptons id", 5 , 9.5, 14.5, abs((zl.second).id()));

    }
    if((*ZL).size() > 0) theHistograms.fill("idParticlesRecoZL"," leptons & daughters id in ZL", 10 , 30.5, 40.5, finalid);
    // if(ZL.size()>0) theHistograms.fill("daughtersZl"," Z daughters id", 100, 0, 500, abs((ZL[0].first).daughter(0).id()));
    // if(ZL.size()>1) theHistograms.fill("daughtersZl"," Z daughters id", 100, 0, 500, abs((ZL[1].first).daughter(0).id()));
    
    //----my particles----//
    cout << "\nmyZl analysis" << endl; 
    //forms Zl every time there are 3 leptons (1 or 2 Zl depending on presence of both e and mu) if leptons pt and eta is good enough
    
    // std::pair<Boson<phys::Lepton> ,phys::Lepton>
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
      //theHistograms.fill("massLeptonMyZl","mass l", 100, 0, 0.12, (std::get<1>(zl)).mass());
      //theHistograms.fill("idDaughterMyZl", " Z daughters id", 5,   9.5, 14.5, abs((zl.first).daughter(0).id()));
      //theHistograms.fill("idLeptonMyZl",   " leptons id",     5,   9.5, 14.5, abs((zl.second).id()));
      
      //cout << "\nZl good cand (right deltaEta and pt)\t" << " leptons forming Z: " << abs(zl.first.daughter(0).id()) << " other lepton: " << abs(zl.second.id()) << endl;
    }
    if (Zl.size() > 0) theHistograms.fill("idParticlesGenZl"," leptons & daughters id in Zl", 10 , 30.5, 40.5, finalid);

    //some histograms
    //if(eptSort.size()  > 0) theHistograms.fill("pte0",  "pt e0",  100, 0, 200, eptSort[0].pt());
    //if(eptSort.size()  > 1) theHistograms.fill("pte1",  "pt e1",  100, 0, 200, eptSort[1].pt());
    //if(eptSort.size()  > 2) theHistograms.fill("pte2",  "pt e2",  100, 0, 200, eptSort[2].pt());
    //if(muptSort.size() > 0) theHistograms.fill("ptmu0", "pt mu0", 100, 0, 200, muptSort[0].pt());
    //if(muptSort.size() > 1) theHistograms.fill("ptmu1", "pt mu1", 100, 0, 200, muptSort[1].pt());
    //if(muptSort.size() > 2) theHistograms.fill("ptmu2", "pt mu2", 100, 0, 200, muptSort[2].pt());

    return;
  }

  //Wlllnu
  else if(leptons.size() == 3 && neutrinos.size() == 1){//3l1nu (finalid == 45 || finalid == 49 || finalid == 53)){
    //  Particle nu = neutrinos[0];
    TLorentzVector Ptot = neutrinos[0].p4();
    foreach(const Particle lep, leptons)
      Ptot += lep.p4();
    double masslllnu = Ptot.M();
    theHistograms.fill("massGenlllnu","mass lllnu", 75, 0, 600, masslllnu);
    cout << "\n masslllnu: " << masslllnu << endl;

    return;
  }


}

//there are still some events where I form Zl but it has no ZL!! (all the ones with 3e)
// event: 7702027  #: 1988 -> 3mu, 2nd mu has pt<10 BUT there're 2 ZL reco
// event: 7701979  #: 1957 -> 3e, all good (eta and pt) BUT no ZL reco
// event: 7701973  #: 1953 -> 3e, IDEM
// event: 7701934  #: 1928 -> 2e1mu, IDEM
// event: 7701871  #: 1894 -> 3e, IDEM
// event: 7701864  #: 1891 -> 2mu1e, IDEM
// event: 7701836  #: 1873 -> 3e, IDEM
// event: 7701827  #: 1867 -> 3e, IDEM
// event: 7701826  #: 1866 -> 3e, IDEM
// event: 7701817  #: 1862 -> 3e, IDEM
// event: 7701815  #: 1860 -> 2mu1e, IDEM
// event: 7701765  #: 1828 -> 2mu1e
// event: 7701753  #: 1823 -> 3e, IDEM
// event: 7701730  #: 1808 -> 2mu1e
// event: 6879026  #: 1741 -> 2mu1e
// event: 6879024  #: 1739 -> 3e
// event: 6878960  #: 1704 -> 3e
// event: 6878948  #: 1699 -> 3mu, 2nd mu has has pt<10 BUT there're 2 ZL reco
// event: 6878934  #: 1691 -> 3e, all good (eta and pt) BUT no ZL reco
// event: 6878932  #: 1690 -> 3mu (!)
// event: 6878916  #: 1680 -> 3e
// event: 6878907  #: 1673 -> 3mu
// event: 6878894  #: 1664 -> 3e
// event: 6878873  #: 1655 -> 2mu1e
// event: 6878834  #: 1633 -> 3e
// event: 6878751  #: 1588 -> 3mu
// event: 6878746  #: 1586 -> 2mu1e
// event: 6878742  #: 1583 -> 3e
// event: 4903971  #: 1541 -> 2e1mu
// event: 4903893  #: 1491 -> 2e1mu
// event: 4903857  #: 1472 -> 3e
// event: 4903788  #: 1440 -> 2e1mu
// event: 4903762  #: 1424 -> 2e1mu
// event: 4903750  #: 1418 -> 3mu
// event: 4903747  #: 1417 -> 2e1mu
// event: 4903741  #: 1416 -> 3e
// event: 4903724  #: 1410 -> 2mu1e
// event: 4903717  #: 1406 -> 3e
// event: 4903708  #: 1402 -> 3mu
// event: 4903668  #: 1375 -> 3e
// event: 4903653  #: 1369 -> 3e
// event: 4903566  #: 1322 -> 3e
// event: 4899910  #: 1287 -> 3e
// event: 4899893  #: 1276 -> 2e1mu
// event: 4899865  #: 1264 -> 3mu
// event: 4899789  #: 1215 -> 2e1mu
// event: 4899779  #: 1210 -> 3mu
// event: 4899772  #: 1204 -> 3e
// event: 4899732  #: 1184 -> 2e1mu
// event: 4899722  #: 1179 -> 3e
// event: 4899711  #: 1174 -> 3mu
// event: 4899689  #: 1165 -> 2e1mu
// event: 4899688  #: 1164 -> 3e
// event: 4899677  #: 1158 -> 2e1mu
// event: 4899662  #: 1150 -> 2e1mu
// event: 4899626  #: 1127 -> 3e
// event: 4899552  #: 1078 -> 2e1mu
// event: 4899496  #: 1042 -> 2mu1e
// event: 4899492  #: 1038 -> 3e
// event: 4877486  #: 1034 -> 3e
// event: 4877485  #: 1033 -> 3e
// event: 4877466  #: 1022 -> 3e
// event: 4877435  #: 1008 -> 3e
// event: 4877369  #: 968 -> 3e
// event: 4877311  #: 932 -> 3mu, 2nd mu has has pt<10 BUT there're 2 ZL reco
// event: 4877304  #: 928 -> 3e
// event: 4877293  #: 921 -> 3e
// event: 4877226  #: 886 -> 2e1mu
// event: 4877210  #: 876 -> 3e
// event: 4877186  #: 864 -> 2mu1e
// event: 4877118  #: 825 -> 3mu
// event: 4877106  #: 816 -> 3e
// event: 4877065  #: 795 -> 3e
// event: 4877045  #: 780 -> 3e
// event: 3553805  #: 757 -> 3mu, 2nd mu has has pt<10 BUT there're 2 ZL reco
// event: 3553753  #: 725 -> 3mu
// event: 3553701  #: 694 -> 3mu
// event: 3553640  #: 662 -> 2mu1e
// event: 3553623  #: 656 -> 3e
// event: 3553572  #: 629 -> 3mu
// event: 3553553  #: 618 -> 2e1mu
// event: 3553543  #: 611 -> 3mu, one mu has eta>2.4 BUT there's 1 ZL reco WITH e AS OTHER LEPTON!!!
// event: 3553469  #: 565 -> 3e
// event: 3553388  #: 509 -> 3e
// event: 1542758  #: 505 -> 3mu
// event: 1542728  #: 486 -> 2mu1e
// event: 1542692  #: 465 -> 3e
// event: 1542684  #: 460 -> 3e
// event: 1542660  #: 441 -> 3e
// event: 1542657  #: 439 -> 2mu1e
// event: 1542575  #: 390 -> 2e1mu
// event: 1542574  #: 389 -> 3e
// event: 1542424  #: 314 -> 2mu1e
// event: 1542419  #: 310 -> 3e
// event: 1542408  #: 302 -> 3e
// event: 1542397  #: 294 -> 3e
// event: 1542387  #: 288 -> 3e
// event: 1542381  #: 284 -> 2e1mu
// event: 1542378  #: 282 -> 2e1mu
// event: 1542370  #: 279 -> 3e
// event: 1542361  #: 272 -> 3e
// event: 1542356  #: 268 -> 3e
// event: 1542332  #: 256 -> 3e
// event: 367724  #: 242 -> 3e
// event: 367620  #: 184 -> 3mu
// event: 367570  #: 162 -> 2mu1e
// event: 367520  #: 135 -> 3e
// event: 367506  #: 127 -> 3mu
// event: 367470  #: 112 -> 3e
// event: 367412  #: 80 -> 3e
// event: 367409  #: 77 -> 3e
// event: 367348  #: 41 -> 3e
// event: 367338  #: 32 -> 3e
// event: 367326  #: 25 -> 2e1mu
// event: 367318  #: 21 -> 3mu
// event: 367313  #: 17 -> 3e
// event: 367283  #: 1 -> 3e
// TOT anomalous events (WZ sample): 117/2000
