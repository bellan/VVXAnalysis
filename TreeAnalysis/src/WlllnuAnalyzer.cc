#include "VVXAnalysis/TreeAnalysis/interface/WlllnuAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include "VVXAnalysis/DataFormats/interface/TypeDefs.h"
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
using namespace colour;

using std::cout;
using std::endl;


using namespace phys;
using namespace physmath;

void WlllnuAnalyzer::begin(){
  nevents = 0;
  mass80Counter = 0;
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
    theHistograms->fill("ptAllGenParticle",   "pt ",   100, 0,   200,  gen.pt());
    //theHistograms->fill("etaAllGenParticle",  "eta ",  100, -10, 10,   gen.eta());
    //theHistograms->fill("YAllGenParticle",    "Y ",    100, -10, 10,  gen.rapidity());
    bool isLepton = abs(gen.id()) == 11 || abs(gen.id()) == 13;
    isLepton ? leptons.push_back(gen) : neutrinos.push_back(gen);
    if (gen.id() == 11){
      electrons.insert(electrons.begin(),gen); //first all e-, then all e+
      //theHistograms->fill("ptElectrons","pt e", 100, 0, 200, gen.pt());
    }
    else if (gen.id() == -11){
      electrons.push_back(gen);
      //theHistograms->fill("ptElectrons","pt e", 100, 0, 200, gen.pt());
    }
    else if (gen.id() == 13){
      muons.insert(muons.begin(),gen); //first all mu-, then all mu+
      //theHistograms->fill("ptMuons","pt mu", 100, 0, 200, gen.pt());
    }
    else if (gen.id() == -13){
      muons.push_back(gen);
      //theHistograms->fill("ptMuons","pt mu", 100, 0, 200, gen.pt());
    }
  }
  //theHistograms->fill("leptonsNumber",  "number of leptons ",  5, -0.5, 4.5, leptons.size());
  //theHistograms->fill("idAllGenParticles"," finalid ", 60 , -0.5, 59.5, finalid);

  //std::vector<pairBosonParticle > Zl; //change Particle to Lepton    
    
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
    theHistograms->fill("massGenZ0", "z0 mass ", 1000, 50, 150, z0.mass());
    theHistograms->fill("massGenZ1", "z1 mass ", 1000, 50, 150, z1.mass());
    theHistograms->fill("massGenZParticles", "Z mass ", 1000, 50, 150, z0.mass());
    theHistograms->fill("massGenZParticles", "Z mass ", 1000, 50, 150, z1.mass());

    cout << "\nz0 " "\tmass: " << z0.mass() << "\tmTdaughters-mZ0: "<< mT(z0.daughter(0), z0.daughter(1)) - z0.mass()
	 << "\n\tz0.daughter(0) \tid: " << z0.daughter(0).id() << "     pt: " << z0.daughter(0).pt()
	 << "\n\tz0.daughter(1) \tid: " << z0.daughter(1).id() << "    pt: " << z0.daughter(1).pt() <<  endl;
    cout << "\nz1 " "\tmass: " << z1.mass() << "\tmTdaughters-mZ1: "<< mT(z1.daughter(0), z1.daughter(1)) - z1.mass()
	 << "\n\tz1.daughter(0) \tid: " << z1.daughter(0).id() << "     pt: " << z1.daughter(0).pt()
	 << "\n\tz1.daughter(1) \tid: " << z1.daughter(1).id() << "    pt: " << z1.daughter(1).pt() <<  endl;
    //comparing mZ with mT daughters
    theHistograms->fill("mZmT", "mZ vs mT daughters", 600, 0, 150, 600, 0, 150, mT(z0.daughter(0),z0.daughter(1)), z0.mass());
    theHistograms->fill("mZmT", "mZ vs mT daughters", 600, 0, 150, 600, 0, 150, mT(z1.daughter(0),z1.daughter(1)), z1.mass());
    theHistograms->fill("deltamZmT", "mZ - mT daughters", 100, -100, 100, z0.mass() - mT(z0.daughter(0),z0.daughter(1)));
    theHistograms->fill("deltamZmT", "mZ - mT daughters", 100, -100, 100, z1.mass() - mT(z1.daughter(0),z1.daughter(1)));
    
    
    //making z0,z1 with mT criteria (mT daughters closer to ZMASS)
    std::stable_sort(Zcandidates.begin(), Zcandidates.end(), mTComparator(ZMASS));
    phys::Boson<phys::Particle> z0mT = Zcandidates[0];
    phys::Boson<phys::Particle> z1mT;
    foreach(const phys::Boson<phys::Particle> cand, Zcandidates)
      if ( cand.daughter(0).pt() != z0mT.daughter(0).pt() &&  cand.daughter(1).pt() != z0mT.daughter(1).pt() ) z1mT = cand;
    theHistograms->fill("massGenZmT", "Z (selected by mT) mass ", 1000, 50, 150, z0mT.mass());
    theHistograms->fill("massGenZmT", "Z (selected by mT) mass ", 1000, 50, 150, z1mT.mass());
    
    //comparing z0mT, z1mT to z0, z1
    int sameZandZmT = ((isTheSame(z0, z0mT) && isTheSame(z1,z1mT)) || (isTheSame(z0, z1mT) && isTheSame(z1,z0mT))) ? 1 : 0;
    theHistograms->fill("sameGenZ", "Z and ZmT are the same? ", 2 , -0.5, 1.5, sameZandZmT);
    int sameZ0Z1 = (isTheSame(z0, z0mT) && isTheSame(z1,z1mT)) ? 1 : 0;
    theHistograms->fill("sameGenZ0z1", "Z and ZmT are the same? (also the order) ", 2 , -0.5, 1.5, sameZ0Z1);
        
    /*  //deltaR, deltaEta, deltaPhi
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
    DiBosonParticle ZZ(z0,z1);
    cout << "\n\n ZZ: " << ZZ.id() << " pt: " << ZZ.pt() << " mass: " << ZZ.mass() << "\n daughters: " << ZZ.daughter<Particle>(0).id() << ", " << ZZ.second().id() << " Y: " << ZZ.rapidity() << endl;
    
    theHistograms->fill("ptGenZZ",   "pt ",   100,  0,   100, ZZ.pt());
    //    theHistograms->fill("etaGenZZ",  "eta ",  100,  -10, 10,  ZZ.eta());
    theHistograms->fill("massGenZZ", "mass ", 100,  0,   400, ZZ.mass());
    theHistograms->fill("YGenZZ",    "Y ",    100,  -10, 10,  ZZ.rapidity());
    
    return;
  }
    
  //  cout << " \n\t\t  -----------------------------------\n " << endl;
 
  else if(leptons.size()<3){
    cout << colour::Red("\nLess than 3 leptons") << endl;
    return;
  }
  else if(neutrinos.size()>1){
    cout << colour::Red("\nMore than 1 nu") << endl;
    return;
  }
  
  //ZL
  else if(leptons.size()==3 && neutrinos.size() == 0){// 3l(finalid == 33 || finalid == 35 || finalid == 37 || finalid == 39)){ //3l

    //----RECO particles----//
        
    cout << Violet("\nZL analysis") << endl; 
    cout << "\nZlcand: \t" << std::get<0>(*ZL) << "\n\t\t" << std::get<1>(*ZL) << endl;
      theHistograms->fill("massBosonRecoZl","mass Z", 100, 0, 500, (std::get<0>(*ZL)).mass());
    
    if(ZL->first.isValid()) theHistograms->fill("idParticlesRecoZL"," leptons & daughters id in ZL", 10 , 30.5, 40.5, finalid);
        
    //----my particles----//
    cout << "\nmyZl analysis" << endl; 
    //forms Zl every time there are 3 leptons (1 or 2 Zl depending on presence of both e and mu) if leptons pt and eta is good enough
    
    std::vector<pairBosonParticle > Zl; //change Particle to Lepton    
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
    /* 
    //check 1st lepton has pt>20 and 2nd pt>10(mu)/12(e)
    if(lepptSort[0].pt()<20 || (abs(lepptSort[1].id())==11 && lepptSort[1].pt()<12) || (abs(lepptSort[1].id())==13 && lepptSort[1].pt()<10)){
      cout<<"pt leptons not sufficient (first 20Gev, second 10Gev(mu)/12Gev(e))"<< endl;
      return;
    }*/
    
    if(finalid == 33){//3e     
      Zl.push_back(pairBosonParticle(phys::Boson<phys::Particle>(electrons[0],electrons[2], 23), electrons[1]));
      electrons[1].id()<0 ? Zl.push_back(pairBosonParticle(phys::Boson<phys::Particle>(electrons[0],electrons[1], 23), electrons[2])) : Zl.push_back(pairBosonParticle(phys::Boson<phys::Particle>(electrons[1],electrons[2], 23), electrons[0]));
    }
    else if(finalid == 39){//3mu
      Zl.push_back(pairBosonParticle(phys::Boson<phys::Particle>(muons[0],muons[2], 23), muons[1]));
      muons[1].id()<0 ? Zl.push_back(pairBosonParticle(phys::Boson<phys::Particle>(muons[0],muons[1], 23), muons[2])) : Zl.push_back(pairBosonParticle(phys::Boson<phys::Particle>(muons[1],muons[2], 23), muons[0]));	
    }
    else if(finalid == 35){//2e1mu
      Zl.push_back(pairBosonParticle(phys::Boson<phys::Particle>(electrons[0],electrons[1], 23), muons[0]));
    }
    else if(finalid == 37){//2mu1e
      Zl.push_back(pairBosonParticle(phys::Boson<phys::Particle>(muons[0],muons[1], 23), electrons[0]));
    }
    
    cout << "\n # of my Zl candidates: "  << Zl.size() << endl;
    foreach(const pairBosonParticle zl, Zl){
      cout << "\nZlcand: \t" << std::get<0>(zl) << "\n\t\t" << std::get<1>(zl) << endl;
      theHistograms->fill("massBosonGenZl","mass Z", 100, 0, 500, (std::get<0>(zl)).mass());
         
      //cout << "\nZl good cand (right deltaEta and pt)\t" << " leptons forming Z: " << abs(zl.first.daughter(0).id()) << " other lepton: " << abs(zl.second.id()) << endl;
    }
    if (Zl.size() > 0) theHistograms->fill("idParticlesGenZl"," leptons & daughters id in Zl", 10 , 30.5, 40.5, finalid);

    //--------------------------------------------------//
    //comparing mZ with mT daughters
    foreach(const pairBosonParticle zl, Zl){
      if(Zl.size() == 0) continue;
      theHistograms->fill("mZlmT", "mZ vs mT daughters", 150, 0, 150, 150, 0, 150, mT((std::get<0>(zl)).daughter(0),(std::get<0>(zl)).daughter(1)), (std::get<0>(zl)).mass());
      theHistograms->fill("deltamZlmT", "mZ - mT daughters", 100, -100, 100, (std::get<0>(zl)).mass() - mT((std::get<0>(zl)).daughter(0),(std::get<0>(zl)).daughter(1)));
    }
    
    //return;
  }

  //Wlllnu
  if(finalid == 45 || finalid == 49 || finalid == 53){///*(leptons.size() == 3 && neutrinos.size() == 1){3l1nu*/ 
    cout << Blue("\nWlllnu analysis") << endl; 
    //theHistograms->fill("idlllnuNOcut"," total id lllnu ", 13 , 43.5, 56.5, finalid);
    
    //----------------------------------------
    //check pt and eta
    std::vector<Particle> lepptSort;
    foreach(const Particle e, electrons){
      if(e.pt()<7 || abs(e.eta())>2.5 ){
	cout<<"pt leptons not sufficient (e 7Gev) or eta over limit (|eta| 2.5)"<< endl;
	return;
      }
      lepptSort.push_back(e);
    }
    foreach(const Particle mu, muons){
      if(mu.pt()<5 || abs(mu.eta())>2.4 ){
	cout<<"pt leptons not sufficient (mu 5Gev) or eta over limit (|eta| 2.5)"<< endl;
	return;
      }
      lepptSort.push_back(mu);
    }
    std::stable_sort(lepptSort.begin(), lepptSort.end(), PtComparator());

    
    /*
    //check 1st lepton has pt>20 and 2nd pt>10(mu)/12(e)
    if(lepptSort[0].pt()<20 || (abs(lepptSort[1].id())==11 && lepptSort[1].pt()<12) || (abs(lepptSort[1].id())==13 && lepptSort[1].pt()<10)){
      cout<<"pt leptons not sufficient (first 20Gev, second 10Gev(mu)/12Gev(e))"<< endl;
      return;
    }*/


    //-----------------------------------------------------
    TLorentzVector Ptot = neutrinos[0].p4();
    foreach(const Particle lep, leptons)
      Ptot += lep.p4();
    double masslllnu = Ptot.M();
    if (masslllnu > 70 && masslllnu < 90) mass80Counter++;
    theHistograms->fill("massGenlllnu","mass lllnu", 1040, 0, 1040, masslllnu);
    //theHistograms->fill("mTGenlllnu","mT lllnu", 400, 40, 440, Ptot.Mt());
    //theHistograms->fill("massMtlllnu", "mass vs mT lllnu", 400, 40, 440, 400, 40, 440, Ptot.Mt(),masslllnu) ;
    //    cout <<"\n events with mass < 165 Gev: " << mass80Counter << endl;
    //cout << "\n masslllnu: " << masslllnu << endl;
    //theHistograms->fill("idlllnu"," total id lllnu ", 13 , 43.5, 56.5, finalid);
    cout << "mass lllnu: " << masslllnu << "\n" << endl;
    if (masslllnu <= 70 || masslllnu >= 90){
      cout << Yellow("Mass lllnu off shell (70-90 Gev)") << endl;
      return;
    }
    Particle nu = neutrinos[0];
    theHistograms->fill("ptnu",   "pt nu",   100, 0,   200,  nu.pt());
    std::vector<boolInt> nuCut;
    nuCut.push_back(boolInt(true, 0));//no nu cut
    nuCut.push_back(boolInt(nu.pt()>5 && nu.eta()<2.4, 5));//"democratic" cut
    nuCut.push_back(boolInt(nu.pt()<8, 8)); //"undemocratic" cut
    nuCut.push_back(boolInt(nu.pt()<10, 10)); //"undemocratic" cut
    nuCut.push_back(boolInt(nu.pt()<15, 15)); //"undemocratic" cut
    nuCut.push_back(boolInt(nu.pt()>54, 54)); //"undemocratic" cut
    // nuCut.push_back(boolInt(nu.pt()<10 && nu.eta()<2.4, 10)); //"undemocratic" cut
    
    foreach(const boolInt theCut, nuCut){
      cout << Violet("\nnuCut: ") << Violet(theCut.second) << endl;
      if(!(theCut.first)){
	cout <<  "nu not valid\n-------------------" << endl;
	continue;
      }
          
      std::vector<pairParticle> pCombos; //all the possible combos with 3l
      pairParticle theCombo; //the best one
      std::vector <pairBosonParticle> Zl; //Z and lepton
      BosonParticle W;
      bool isNuAlone = NULL; //in order to mark W -> nu + lll or W -> l + llnu
      //int diagramId = 0;
      
      if(finalid == 49 && electrons.size() == 2){//2e 1mu 1nu_mu
	Zl.push_back(pairBosonParticle(BosonParticle(electrons[0], electrons[1], 23),muons[0]));
	pCombos.push_back(pairParticle(nu,Particle(electrons[0].p4()+electrons[1].p4()+muons[0].p4())));
	pCombos.push_back(pairParticle(muons[0],Particle(electrons[0].p4()+electrons[1].p4()+nu.p4())));
	std::stable_sort(pCombos.begin(), pCombos.end(), mTComparator(masslllnu));
	theCombo = pCombos[0];
	isNuAlone = isTheSame(nu, theCombo.first);
      }
      else if(finalid == 49 && muons.size() == 2){//2mu 1e 1nu_e
	Zl.push_back(pairBosonParticle(BosonParticle(muons[0], muons[1], 23), electrons[0]));
	pCombos.push_back(pairParticle(nu,Particle(muons[0].p4()+muons[1].p4()+electrons[0].p4())));
	pCombos.push_back(pairParticle(electrons[0],Particle(muons[0].p4()+muons[1].p4()+nu.p4())));
	std::stable_sort(pCombos.begin(), pCombos.end(), mTComparator(masslllnu));
	theCombo = pCombos[0];
	isNuAlone = isTheSame(nu, theCombo.first);
      }
      
      else if(finalid == 45){//3e 1nu_e
	pCombos.push_back(pairParticle(nu,Particle(electrons[0].p4()+electrons[1].p4()+electrons[2].p4())));
	pCombos.push_back(pairParticle(electrons[1],Particle(nu.p4()+electrons[0].p4()+electrons[2].p4())));
	electrons[1].id()>0 ?
	  pCombos.push_back(pairParticle(electrons[0],Particle(nu.p4()+electrons[1].p4()+electrons[2].p4()))) :
	  pCombos.push_back(pairParticle(electrons[2],Particle(nu.p4()+electrons[0].p4()+electrons[1].p4())));
	std::stable_sort(pCombos.begin(), pCombos.end(), mTComparator(masslllnu));
	theCombo = pCombos[0];
	// foreach(const pairParticle p, pCombos)
	// if(abs(mT(p.first, p.second)-masslllnu) < abs(mT(theCombo.first, theCombo.second)-masslllnu))
	//   theCombo = p; //the best Combo has closest mT to lllnu (W)
	isNuAlone = isTheSame(nu, theCombo.first);
	if(isNuAlone){
	  Zl.push_back(pairBosonParticle(BosonParticle(electrons[0], electrons[2], 23), electrons[1]));
	  electrons[1].id()>0 ?
	    Zl.push_back(pairBosonParticle(BosonParticle(electrons[1], electrons[2], 23), electrons[0])) :
	    Zl.push_back(pairBosonParticle(BosonParticle(electrons[0], electrons[1], 23), electrons[2]));
	}
	else if(isTheSame(electrons[0],theCombo.first))
	  Zl.push_back(pairBosonParticle(BosonParticle(electrons[1], electrons[2], 23), electrons[0]));
	else if(isTheSame(electrons[1],theCombo.first))
	  Zl.push_back(pairBosonParticle(BosonParticle(electrons[0], electrons[2], 23), electrons[1]));
	else if(isTheSame(electrons[2],theCombo.first))
	  Zl.push_back(pairBosonParticle(BosonParticle(electrons[0], electrons[1], 23), electrons[2]));
      }
      
      else if(finalid == 53){//3mu 1nu_mu
	pCombos.push_back(pairParticle(nu,Particle(muons[0].p4()+muons[1].p4()+muons[2].p4())));
	pCombos.push_back(pairParticle(muons[1],Particle(nu.p4()+muons[0].p4()+muons[2].p4())));
	muons[1].id()>0 ?
	  pCombos.push_back(pairParticle(muons[0],Particle(nu.p4()+muons[1].p4()+muons[2].p4()))) :
	  pCombos.push_back(pairParticle(muons[2],Particle(nu.p4()+muons[0].p4()+muons[1].p4())));
	std::stable_sort(pCombos.begin(), pCombos.end(), mTComparator(masslllnu));
	theCombo = pCombos[0];
	// foreach(const pairParticle p, pCombos)
	//   if(abs(mT(p.first, p.second) - masslllnu) < abs(mT(theCombo.first, theCombo.second) - masslllnu))
	//   theCombo = p;
	isNuAlone = isTheSame(nu, theCombo.first);
	if(isNuAlone){
	  Zl.push_back(pairBosonParticle(BosonParticle(muons[0], muons[2], 23), muons[1]));
	  muons[1].id()>0 ?
	    Zl.push_back(pairBosonParticle(BosonParticle(muons[1], muons[2], 23), muons[0])) :
	    Zl.push_back(pairBosonParticle(BosonParticle(muons[0], muons[1], 23), muons[2]));
	}
	else if(isTheSame(muons[0],theCombo.first))
	  Zl.push_back(pairBosonParticle(BosonParticle(muons[1], muons[2], 23), muons[0]));
	else if(isTheSame(muons[1],theCombo.first))
	  Zl.push_back(pairBosonParticle(BosonParticle(muons[0], muons[2], 23), muons[1]));
	else if(isTheSame(muons[2],theCombo.first))
	  Zl.push_back(pairBosonParticle(BosonParticle(muons[0], muons[1], 23), muons[2]));
      }
      
      
      else cout << Red("invalid total id") << endl;
      cout << Green("\n # of valid Zl: ")  << Green(Zl.size()) << endl;
      
      theHistograms->fill("isNuAlone","type of diagram", 56, -0.5, 55.5, isNuAlone+theCut.second);
      
      if(theCut.second != 0) continue;
      
      foreach(const pairParticle pcombos, pCombos){
	theHistograms->fill("mTAllCombos","mT all combo", 100, 0, 100, abs(mT(pcombos.first, pcombos.second)-masslllnu));
      }
      
      theHistograms->fill("mTcombo","mT combo", 100, 0, 100, mT(theCombo.first, theCombo.second));
      theHistograms->fill("mTcombo1","mT second combo", 100, 0, 100, mT(pCombos[1].first, pCombos[1].second));
      theHistograms->fill("mTcombo2","mT third combo", 100, 0, 100, mT(pCombos[2].first, pCombos[2].second));
      theHistograms->fill("mTcomboMasslllnu", "mT theCombo vs mass lllnu", 100, 50, 100, 200, 0, 100, masslllnu,mT(theCombo.first, theCombo.second));
      //diagramId = abs(nu.id()) + abs(Zl[0].second.id()) + abs(Zl[0].first.daughter(0).id());
      //isNuAlone ? theHistograms->fill("diagramId","type of diagram", 21, -10.5, 10.5, diagramId-30) : theHistograms->fill("diagramId","type of diagram", 21, -10.5, 10.5, -(diagramId-30));
      //theHistograms->fill("absDiagId","type of diagram", 7, 3.5, 10.5, diagramId-30);
      isNuAlone ? theHistograms->fill("ptnuAlone",   "pt nu when is alone",   100, 0,   200,  nu.pt()) : theHistograms->fill("ptnuNoAlone",   "pt nu when is forming Z",   100, 0,   200,  nu.pt());
    
      if(Zl.size() > 1){
	//std::stable_sort(Zl.begin(), Zl.end(), deltaRComparator());
	//theHistograms->fill("deltaRZ0","deltaR best Zcand", 100, 0, 10, physmath::deltaR(Zl[0].first.daughter(0), Zl[0].first.daughter(1)));
	//theHistograms->fill("deltaRZ1","deltaR other Zcand", 100, 0, 10, physmath::deltaR(Zl[1].first.daughter(0), Zl[1].first.daughter(1)));
	//cout << "\n Z0: " << Zl[0].first << endl;
	//cout << "\n l0: " << Zl[0].second << endl;
	//cout << "\n deltaRZ0:\t " << physmath::deltaR(Zl[0].first.daughter(0), Zl[0].first.daughter(1)) << endl;
	//cout << "\n\n Z1: " << Zl[1].first << endl;
	//cout << "\n l1: " << Zl[1].second << endl;
	//cout << "\n deltaRZ1:\t " << physmath::deltaR(Zl[1].first.daughter(0), Zl[1].first.daughter(1)) << endl;
	continue;
      }
      
      //theHistograms->fill("deltaRl","deltaR leptons couples", 100, 0, 10, deltaR(leptons[0].p4(), leptons[1].p4()));
      //theHistograms->fill("deltaRl","deltaR leptons couples", 100, 0, 10, deltaR(leptons[0].p4(), leptons[2].p4()));
      //theHistograms->fill("deltaRl","deltaR leptons couples", 100, 0, 10, deltaR(leptons[1].p4(), leptons[2].p4()));
      
      foreach(const pairBosonParticle zl, Zl){
	
	W = BosonParticle(theCombo.first, theCombo.second, copysign(24, zl.second.charge()));
	cout << "is nu Alone:\t";
	isNuAlone ? cout << "yes" : cout << "no";
	cout << endl;
	cout << "\n Z: " << zl.first << endl;
	cout << "\n l: " << zl.second << endl;
	cout << "\n nu: " << nu << endl;
	cout << "\n W: " << W << endl;
	
	if(finalid != 49) continue;
	isNuAlone ? theHistograms->fill("deltaRZtriplet","deltaR Z and related leptons ", 100, 0, 10, physmath::deltaR(zl.first, zl.second)) : theHistograms->fill("deltaRZtriplet","deltaR Z and related leptons ", 100, 0, 10, physmath::deltaR(zl.first, nu)); //are Z and its "related" lepton collinear?
	!(isNuAlone) ? theHistograms->fill("deltaRZsinglet","deltaR Z and NOT related leptons ", 100, 0, 10, physmath::deltaR(zl.first, zl.second)) : theHistograms->fill("deltaRZsinglet","deltaR Z and NOT related leptons ", 100, 0, 10, physmath::deltaR(zl.first, nu)); //are Z and its "related" lepton collinear?
	theHistograms->fill("deltaRZdaughters","deltaR Z daughters", 100, 0, 10, physmath::deltaR(zl.first.daughter(0), zl.first.daughter(1)));
	theHistograms->fill("massZ","mass Z", 300, 0, 150, zl.first.mass());
	//theHistograms->fill("mTZ","direct mT Z", 300, 0, 150, zl.first.p4().Mt());
	theHistograms->fill("mTZdaughters","mT Z from daughters", 300, 0, 150, mT(zl.first.daughter(0), zl.first.daughter(1)));
	theHistograms->fill("massZmTZdaughters", "mass Z vs mT Z daughters", 100, 0, 100, 100, 0, 100, mT(zl.first.daughter(0), zl.first.daughter(1)), zl.first.mass());
	theHistograms->fill("ptZ",   "pt Z",   100, 0,   200,  zl.first.pt());
	theHistograms->fill("ptl",   "pt l",   125, 0,   125,  zl.second.pt());
	isNuAlone ? theHistograms->fill("ptlAlone",   "pt l when is alone",   125, 0,   125,  zl.second.pt()) : theHistograms->fill("ptlNoAlone",   "pt l when is forming Z",   125, 0,   125,  zl.second.pt());
	theHistograms->fill("massW","mass W", 300, 0, 150, W.mass());
	//	theHistograms->fill("chargeW","charge W", 3, -1.5, 1.5, copysign(1, W.id()));
      }
    }
    
  }
    return;
    
}


void WlllnuAnalyzer::end(TFile &){
  if (mass80Counter > 0)
    cout << Yellow("\n-----------------------------------------------------------\n")
	 << " Events with lllnu mass between 70 and 90 Gev: " << mass80Counter
	 << Yellow("\n-----------------------------------------------------------\n") << endl;
  return;
}
