#include "VVXAnalysis/TreeAnalysis/interface/VZZAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/DataFormats/interface/Jet.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/TypeDefs.h"

//#define PY_SSIZE_T_CLEAN
#ifdef USE_PYTHON
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#endif

#include <iostream>
#include <fstream>			// open(), close(), <<
#include <string>				// find_last_of()
#include <time.h>				// clock_t, clock()
#include <utility>			// std::pair, std::make_pair()
#include <tuple>        // std::tuple, std::make_tuple()
#include <list>         // std::list
#include <algorithm>    // std::min_element(), std::max()
#include <new>          // ::operator new
#include "TSystem.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "THStack.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TString.h"

#include "boost/foreach.hpp"
#define foreach BOOST_FOREACH
#include "boost/assign/std/vector.hpp" 

#define PT_SIZE 40,0.,400.
#define E_SIZE 50,0.,1500.
#define ETA_SIZE 51,-5.05,5.05
#define M_SIZE 35,50.,120.
#define CAND_M_SIZE 40,50.,150.
//Dimension of 2D plots
#define MASS_2D_SIZE 12,60.,120.
#define P_2D_SIZE 25,0.,1000.
#define PT_2D_SIZE 18,0.,540.
#define ZZMASS_2D_SIZE 30,100.,700.
#define S_2D_SIZE 30,100.,1600.
#define BINARY_SIZE 2,-1.5,1.5
#define CUTS_SIZE 3,-0.5,2.5
#define CUT_AN_SIZE 4,-0.5,3.5

using namespace boost::assign;

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;
using std::pair;
using std::make_pair;
using std::tuple;
using std::make_tuple;

using namespace phys;


void VZZAnalyzer::begin(){
	cout<<'\n';
	for(char i=0; i<25; ++i) cout<<'-';
	cout<<" Start of VZZAnalyzer ";
	for(char i=0; i<25; ++i) cout<<'-';
	cout<<'\n';
	
	startTime_ = clock();
	
	//Memory allocation
	AK4pairs_         = new vector<Boson<Jet>>;
	genHadVBs_        = new vector<Boson<Particle>>;
	AllGenVBjj_       = new vector<Boson<Particle>>;
	AK4GenRec_        = new vector<pair<Boson<Particle>, Boson<Jet>>>;
	AK8WithPred_      = new vector<pair<Jet, double>>;
	AK4pairsWithPred_ = new vector<pair<Boson<Jet>, double>>;
	genZZ_ = (DiBoson<Particle, Particle>*) ::operator new (sizeof(DiBoson<Particle, Particle>));  // Allocates memory without calling constructor
	//qq_    = (Boson<Particle>*) ::operator new (sizeof(Boson<Particle>));
	sigVB_ = (Boson<Particle>*) ::operator new (sizeof(Boson<Particle>));
	
	
	//Python: scikit predictors
	#ifdef USE_PYTHON
	int pyStatus = initPy();
	if(pyStatus){
		std::cerr<<"Failed to initialize py predictors"<<std::endl;
		exit(pyStatus);
	}
	#endif
	
	std::string spaces( ceil(log10(tree()->GetEntries())), ' ' );
	cout<<"Analyzed:\t"<<spaces<<'/'<<tree()->GetEntries()<<std::flush;
	return;
}


Int_t VZZAnalyzer::cut(){
	++evtN_;
	cout<<"\r\t\t"<<evtN_;
	totEvtW_ += theWeight;
	
	theHistograms->fill("TOT_weight", "Total weight", 1,0.,1., 0.5, theWeight);
	
	//Cleanup of previous event
	if(genHadVBs_)     genHadVBs_->clear(); //Destroys objects but keeps memory allocated
	if(AK4pairs_)       AK4pairs_->clear();
	if(AK4GenRec_)     AK4GenRec_->clear();
	if(AllGenVBjj_)   AllGenVBjj_->clear();
	if(AK8WithPred_) AK8WithPred_->clear();
	if(AK4pairsWithPred_)AK4pairsWithPred_->clear(); 
	if(qq_.p() > 1.)      qq_    = Boson<Particle>();  // Making invalid
	if(sigVB_->p() > 1.) *sigVB_ = Particle();
	
	//test topology(0) --> is signal region ZZ
	//AND
	//4 --> Fiducial acceptance 
	// per avere segnale serve^ (bisogna girare su MC, altrimenti non ha senso)
	//if( !topology.test(0) || !topology.test(4) ) 
	
	//Preparation for this event
	fillGenHadVBs();
	fillRecHadVBs();
	//fillAK4GenRec(false);
	fillGenVBtoAK4();
	makeGenZZ();
	makeQQ(true);
	
	
	baseHistos();
	//ZZGraphs();
	ZZRecoEfficiency();
	jetRecoEfficiency();
	
	
	// GEN SIGNAL
	sigType_ = isSignal(false/*Uses genJets(AK4) and genJetsAK8*/);
	theHistograms->fill("Signal Type","Signal Type;type;events",6,-1.5,4.5, -sigType_,theWeight);
	
	theHistograms->fill("Cuts gen_f_", "Cuts on gen variables;;weighted evts", CUTS_SIZE, 0., theWeight);
	if(topology.test(0) && topology.test(4)){  // equivalent to genZZ_->isValid()
		theHistograms->fill("Cuts gen_f_", "Cuts on gen variables;;weighted evts", CUTS_SIZE, 1., theWeight);
		if(sigType_ > 0)
			theHistograms->fill("Cuts gen_f_", "Cuts on gen variables;;weighted evts", CUTS_SIZE, 2., theWeight);
	 }
	
	// REC BASELINE
	theHistograms->fill("Cuts rec_f_", "Cuts on rec variables;;weighted evts", CUTS_SIZE, 0., theWeight);
	if(!ZZ || ZZ->p() < 1.) return -1;  //Maybe move to analyze()?
	theHistograms->fill("Cuts rec_f_", "Cuts on rec variables;;weighted evts", CUTS_SIZE, 1., theWeight);
	if(jets->size() < 2. && jetsAK8->size() == 0) return -2;
	theHistograms->fill("Cuts rec_f_", "Cuts on rec variables;;weighted evts", CUTS_SIZE, 2., theWeight);
	
	return 1;
}


void VZZAnalyzer::analyze(){
	++analyzedN_; analyzedW_ += theWeight;
	
	#ifdef USE_PYTHON	
	// AK4pairsWithPred_
	if(!AK4_classifier_) {
		cout<<"!AK4_classifier_\n";
		return;
	}
	foreach(const Boson<Jet>& jj, *AK4pairs_)
		if(jj.mass() > 50. && jj.mass() < 120.){
			vector<double>* feat = getAK4features(jj);
			//double pred = getPyPrediction(*feat, AK4_classifier_);
			//AK4pairsWithPred_->emplace_back(pair<Boson<Jet>, double>(jj, pred));
			delete feat;
		}
	vector<double>* predAK4 = getPyPredAll(*AK4pairs_, AK4_classifier_);
	delete predAK4;
	
	
	return; // TODO TEMP
	
	
	if(AK4pairsWithPred_->size() > 0){
		//auto it4 = max_element(AK4pairsWithPred_->begin(), AK4pairsWithPred_->end(), [](pair<Boson<Jet>, double>& a, pair<Boson<Jet>, double>& b) { return a.second < b.second;});  //lambda expressions are useful but cumbersome
		auto it4 = max_element(AK4pairsWithPred_->begin(), AK4pairsWithPred_->end(), PairComparator());
		theHistograms->fill("BestAK4 score", "BestAK4 score", 10,0.,1., it4->second, theWeight);
	}
	
	// AK8WithPred_
	/*
	foreach(const Jet& j, *jetsAK8)
		if(j.mass() > 50. && j.mass() < 120.){
			vector<double>* feat = getAK8features(j);
			double pred = getPyPrediction(*feat, AK8_classifier_);
			AK8WithPred_->emplace_back(pair<Jet, double>(j, pred));
			delete feat;
		}
	if(AK8WithPred_->size() > 0){
		auto it8 = max_element(AK8WithPred_->begin(), AK8WithPred_->end(), PairComparator());
		theHistograms->fill("BestAK8 score", "BestAK8 score", 10,0.,1., it8->second, theWeight);
	}*/
	#endif
	
	//genSignalGraphs();
	//recSignalGraphs();
	//recSignalAnalysis();
	
	//genHadVBCategorization();  //test topology bit 8 and 9
	//genTauAnalisys();
	//AK8nearGenHadVB();  //Are there any gen/rec AK8 near a genHadVB (pair of gen AK4)?
	
	//reconstructionAK();
	
	//AK8MassAlgorithms();
	//ptCutMVA();
	//closestJetAnalisys();
	//furthestJetMVA();
	//minPtJetMVA();
	
	//bestZMassJetMVA();
	resolutionZmass();  // same as bestZMassJetMVA(), but without dR(jet, ZZ) > 2. Also, gen part is done with respect to genZZ_ and not ZZ
	
	simpleGraphs();
	return;
}

void VZZAnalyzer::end(TFile& fout){
	cout<<'\n';
	
	//Final cleanup
	if(genHadVBs_)   delete genHadVBs_; //Deallocates memory
	if(AK4pairs_)    delete AK4pairs_;
	if(AllGenVBjj_)  delete AllGenVBjj_;
	if(AK4GenRec_)   delete AK4GenRec_;
	if(AK8WithPred_) delete AK8WithPred_;
	if(AK4pairsWithPred_)delete AK4pairsWithPred_;
	
	if(win4_ || win8_)  // This counter was incremented if closestJetAnalisys() was called
		endClosestJetAn();
	
	
	endGenHadVbCateg();
	endNameHistos();
	//endResolutionAnalisys(fout);
	//endSignalEff(fout);
	
	cout<<Form("\nTotal events: %lu (weighted: %.2f)", evtN_, totEvtW_);
	cout<<Form("\nPassing cut:  %lu (weighted: %.2f)", analyzedN_, analyzedW_)<<'\n';
	endReconstructionAK();  // Writes to cout only if reconstructionAK4() has been called
	
	
	float elapsedSec = (float)(clock()-startTime_)/CLOCKS_PER_SEC;
	int elapsedSecInt = (int)elapsedSec;
	cout<<"\nElapsed Time: "<<elapsedSec<<" s\t\t("<<elapsedSecInt/60<<"\' "<<elapsedSecInt%60<<"\")\n";
	for(char i=0; i<25; ++i) cout<<'-';
	cout<<" End of VZZAnalyzer ";
	for(char i=0; i<25; ++i) cout<<'-';
	cout<<"\n\n";
}


int VZZAnalyzer::isSignal(bool doGraphs){ //isHadSignal
	if(qq_.p() < 1.) return -1;      // -1: no qq_
	
	if(genJets->size() < 1) return -2;  // -2: no 2 genAK4
	
	vector<Particle> cpGenJets(*genJets);
	size_t i4_0 = 99;
	const Particle* gen4_0 = closestSing(&cpGenJets, qq_.daughter(0), i4_0);
	size_t i4_1 = 99;
	const Particle* gen4_1 = closestSing(&cpGenJets, qq_.daughter(1), i4_1);
	
	if(gen4_0 && gen4_1){
		float dR0 = physmath::deltaR(*gen4_0, qq_.daughter(0));
		float dR1 = physmath::deltaR(*gen4_1, qq_.daughter(1));
		
		if(physmath::deltaR(*gen4_0, *gen4_1) < 0.01){ // If they're matched to the same jet
			if(dR0 < dR1){ // Search for a second, different jet
				cpGenJets.erase(cpGenJets.begin() + i4_0);
				gen4_1 = closestSing(&cpGenJets, qq_.daughter(1), i4_1);
			} else if(dR1 < dR0){
				cpGenJets.erase(cpGenJets.begin() + i4_1);
				gen4_0 = closestSing(&cpGenJets, qq_.daughter(1), i4_0);
			}
		}
	}
	
	if(!gen4_0 || !gen4_1) return -3;  // -3: no 2 distinct genAK4
	
	float mass = (gen4_0->p4() + gen4_1->p4()).M();
	if( phys::WMASS-30 < mass && mass < phys::ZMASS+30 ){
		new (sigVB_) Boson<Particle>(*gen4_0, *gen4_1);
		return 1;
	}
	else
		return -4;                       // -4: genAK4 with mass out of range
	
	return 0;  // Should never reach this point
}


void VZZAnalyzer::genSignalGraphs(){
	theHistograms->fill("Sig Gen: ZZ eta _f_", "Sig Gen: |ZZ eta|;|eta|", 20,0.,6., fabs(ZZ->eta()), theWeight);
	theHistograms->fill("Sig Gen: ZZ pt _f_", "Sig Gen: ZZ pt;pt [GeV/c]", 20,0.,400., fabs(ZZ->pt()), theWeight);
	
	theHistograms->fill("Sig Gen: eta _r_","Sig Gen: |eta|;|eta|", 20,0.,6.,fabs(sigVB_->eta()), theWeight);
	theHistograms->fill("Sig Gen: pt _f_", "Sig Gen: pt;pt [GeV/c]", 20,0.,400., sigVB_->pt(), theWeight);
	theHistograms->fill("Sig Gen: E _f_", "Sig Gen: E;E [GeV]", 20,0.,600., sigVB_->e(), theWeight);
	
	theHistograms->fill("Sig Gen: pt V+ZZ _r_", "Sig Gen: pt V+ZZ;pt [GeV/c]", 20,0.,400., (sigVB_->p4() + ZZ->p4()).Pt(), theWeight);
	
	theHistograms->fill("Sig Gen: min dM _r_", "Sig Gen: min #DeltaM;#DeltaM [GeV/c^{2}]", 20,0.,30, std::min( fabs(sigVB_->mass() - phys::WMASS), fabs(sigVB_->mass() - phys::WMASS) ), theWeight);
	
	//All the other jets
	vector<Particle> allOtherJets;
	switch(sigType_){
		case 1:{  // {}: see https://stackoverflow.com/questions/5685471/error-jump-to-case-label
			Boson<Particle>* BosSigVB = static_cast<Boson<Particle>*>(sigVB_); //safe: dynamic_cast
			foreach(const Particle& j, *genJets)
				if(  std::min( physmath::deltaR(BosSigVB->daughter(0), j), physmath::deltaR(BosSigVB->daughter(0), j) ) > 0.1  )
					allOtherJets.push_back(Particle(j));
			break;
		}
		case 2:
			foreach(const Particle& J, *genJetsAK8){
				if(physmath::deltaR(*sigVB_, J) > 0.1)
					allOtherJets.push_back(Particle(J));
			}
			break;
	}
	
	if(sigType_ && allOtherJets.size() > 0){
		stable_sort(allOtherJets.begin(), allOtherJets.end(), EComparator());
		theHistograms->fill("Sig Gen: highest E (other) _f_", "Sig Gen: highest E (other jets);E [GeV]", 20,0.,600., allOtherJets.front().e(), theWeight);
		if(allOtherJets.size() > 1)
			theHistograms->fill("Sig Gen: 2nd highest E (other) _f_", "Sig Gen: 2^{nd} highest E (other jets);E [GeV]", 20,0.,600., allOtherJets.at(1).e(), theWeight);
		
		stable_sort(allOtherJets.begin(), allOtherJets.end(), PtComparator());
		theHistograms->fill("Sig Gen: highest pt (other) _f_", "Sig Gen: highest pt (other jets);pt [GeV/c]", 20,0.,300., allOtherJets.front().pt(), theWeight);
		if(allOtherJets.size() > 1)
			theHistograms->fill("Sig Gen: 2nd highest pt (other) _f_", "Sig Gen: 2^{nd} highest pt (other jets);pt [GeV/c]", 20,0.,300., allOtherJets.at(1).pt(), theWeight);
	}
}


void VZZAnalyzer::recSignalGraphs(){
	Particle* candClosest = nullptr;
	switch(sigType_){
		case 1:
			if(AK4pairs_->size() == 0) break;  // && jetsAK8->size() > 0) goto ak8_label;
			stable_sort(AK4pairs_->begin(), AK4pairs_->end(), DeltaRComparator(*sigVB_) );
			if(physmath::deltaR(AK4pairs_->front(), *sigVB_) < 0.4)
				candClosest = &(AK4pairs_->front());
			break;
		case 2:
			if(jetsAK8->size() == 0) break;  // && AK4pairs_->size() > 0) goto ak4_label;
			stable_sort(jetsAK8->begin(), jetsAK8->end(), DeltaRComparator(*sigVB_) );
			if(physmath::deltaR(jetsAK8->front(), *sigVB_) < 0.4)
				candClosest = &(jetsAK8->front());
			break;
		default:
			cout<<"Error, sigType_ has an invalid value in recSignalGraphs()\n";
			exit(1);
	}
	
	if(!candClosest) return;  // Either there's no rec VB or it's too far in dR
	
	
	TLorentzVector candp4(candClosest->p4());
	TLorentzVector Z1p4(ZZ->first().p4());
	TLorentzVector Z2p4(ZZ->second().p4());
	theHistograms->fill("Rec: pt Had_f_", "Rec: Had pt;pt [GeV/c]", 25,0.,500.,candp4.Pt(), theWeight);
	float refinedM = getRefinedMass(*candClosest);
	theHistograms->fill("Rec: minDM Had_r_", "Rec: minDM Had;#DeltaM [GeV/c^{2}]", 30,0.,30., physmath::minDM(refinedM), theWeight);
	
	// Angles
	float ang0 = Z1p4.Angle(Z2p4.Vect());
	float ang1 = Z1p4.Angle(candp4.Vect());
	float ang2 = candp4.Angle(Z2p4.Vect());
	theHistograms->fill("Rec: ang0 _f_", "Rec: ang0;", 20,0.,2., ang0, theWeight);
	theHistograms->fill("Rec: ang1 _f_", "Rec: ang1;", 20,0.,2., ang1, theWeight);
	theHistograms->fill("Rec: ang2 _f_", "Rec: ang2;", 20,0.,2., ang2, theWeight);
	/*
	// phi
	float dPhi0 = physmath::deltaPhi(ZZ->first(), ZZ->second());
	float dPhi1 = physmath::deltaPhi(ZZ->first(), *candClosest);
	float dPhi2 = physmath::deltaPhi(*candClosest, ZZ->second());
	
	// eta
	float dEta0 = ZZ->first().eta() - ZZ->second().eta();
	float dEta1 = ZZ->first().eta() - candClosest->eta();
	float dEta2 = candClosest->eta() - ZZ->second().eta();
	
	// R
	float dR0 = physmath::deltaPhi(ZZ->first(), ZZ->second());
	float dR1 = physmath::deltaPhi(ZZ->first(), *candClosest);
	float dR2 = physmath::deltaPhi(*candClosest, ZZ->second());
	*/
	// projection - relative pt
	float pt0 = candp4.P() * sin( candp4.Angle(ZZ->p4().Vect()) );
	float pt1 = Z1p4.P() * sin( Z1p4.Angle((candp4 + Z2p4).Vect()) );
	float pt2 = Z2p4.P() * sin( Z2p4.Angle((candp4 + Z1p4).Vect()) );
	theHistograms->fill("Rec: pt0 _r_", "Rec: pt0;pt [GeV/c]", 20,0.,200., pt0, theWeight);
	theHistograms->fill("Rec: pt1 _r_", "Rec: pt1;pt [GeV/c]", 20,0.,200., pt1, theWeight);
	theHistograms->fill("Rec: pt2 _r_", "Rec: pt2;pt [GeV/c]", 20,0.,200., pt2, theWeight);
}


void VZZAnalyzer::recSignalAnalysis(){
	Particle* candClosest;
	switch(sigType_){
		case 1:
			candClosest = &( *std::min_element(AK4pairs_->begin(), AK4pairs_->end(), DeltaRComparator(*sigVB_)) );
			break;
		case 2:
			candClosest = &( *std::min_element(jetsAK8->begin(), jetsAK8->end(), DeltaRComparator(*sigVB_)) );
			break;
		default:
			cout<<"Error, sigType_ has an invalid value in recSignalAnalysis()\n";
			exit(1);
	}
	
	const char* name_dR = Form("Sig Rec closest AK%d: #DeltaR", sigType_*4);
	float dR_closest = physmath::deltaR(*sigVB_, *candClosest);
	theHistograms->fill(name_dR, name_dR, 20,0.,2., dR_closest, 1.);
	
	//Efficiency
	if(dR_closest < 0.5){
		const char* name_eta = Form("Sig Rec AK%d: #eta of gen", sigType_*4);  // 4*1 = 4, 4*2 = 8
		theHistograms->fill(name_eta, name_eta, 20,-5.,5., sigVB_->eta(), 1.);
		const char* name_pt = Form("Sig Rec AK%d: pt of gen", sigType_*4);
		theHistograms->fill(name_pt, name_pt, 20,0.,500., sigVB_->pt(), 1.);
		
		//Resolution
		const char* name_resE = Form("Sig Resolution AK%d: E", sigType_*4);
		theHistograms->fill(name_resE, name_resE, 20,-100.,100., candClosest->e()-sigVB_->e(), 1.);
		const char* name_resM = Form("Sig Resolution AK%d: M", sigType_*4);
		theHistograms->fill(name_resM, name_resM, 20,-50.,50., candClosest->mass()-sigVB_->mass());
	}
}


void VZZAnalyzer::endSignalEff(TFile& f_out){
	TH1* gen_eta = theHistograms->get("Sig Gen AK4: #eta");
	TH1* gen_pt  = theHistograms->get("Sig Gen AK4: pt");
	TH1* rec_eta = theHistograms->get("Sig Rec AK4: #eta of gen");
	TH1* rec_pt  = theHistograms->get("Sig Rec AK4: pt of gen");
	
	f_out.cd();
	if(gen_eta != nullptr && rec_eta != nullptr){
		TGraphAsymmErrors* etaEff  = new TGraphAsymmErrors(rec_eta,  gen_eta,  "cp");
		etaEff->GetYaxis()->SetRangeUser(0.,1.01);
		etaEff->SetTitle("Sig efficiency vs #eta");
		etaEff->SetName("Sig efficiency vs #eta");
		etaEff->Write();
	}
	if(gen_pt != nullptr && rec_pt != nullptr){
		TGraphAsymmErrors* ptEff  = new TGraphAsymmErrors(rec_pt,  gen_pt,  "cp");
		ptEff->GetYaxis()->SetRangeUser(0.,1.01);
		ptEff->SetTitle("Sig efficiency vs pt");
		ptEff->SetName("Sig efficiency vs pt");
		ptEff->Write();
	}
}


void VZZAnalyzer::baseHistos(){
	//Topology
	if(topology.test(0))
		theHistograms->fill("topology_u", "topology (unweighted)", 3,-0.5,2.5, 2);
	if(topology.test(4))
		theHistograms->fill("topology_u", "topology (unweighted)", 3,-0.5,2.5, 2);
	if(topology.test(6))
		theHistograms->fill("topology_u", "topology (unweighted)", 3,-0.5,2.5, 2);
	
	
	Mass2Comparator comp = Mass2Comparator(phys::WMASS, phys::ZMASS);
	auto it4 = std::min_element(AK4pairs_->begin(), AK4pairs_->end(), comp);
	auto it8 = std::min_element(jetsAK8->begin()  , jetsAK8->end()  , comp);
	//VCandType temp;
	Particle*       cand4 = &(*it4); //findBestVFromPair(jets,    temp);
	const Particle* cand8 = &(*it8); //findBestVFromSing(jetsAK8, temp);
	
	bool genHAD = sigType_ > 0 ? true : false;
	bool recHAD = (cand4 || cand8);
	if(genHAD)
		theHistograms->fill("Signal_HAD", "Signal (HAD)", 1,0.,1., 0.5, theWeight);
	if(recHAD)
		theHistograms->fill("Baseline_HAD", "Baseline (HAD)", 1,0.,1., 0.5, theWeight);
	if(genHAD && recHAD)
		theHistograms->fill("Sig_Base_HAD", "Signal + Baseline (HAD)", 1,0.,1., 0.5, theWeight);
	
	bool genLEP = topology.test(0) && topology.test(4);  //genZZ_ && genZZ_->isValid();
	bool recLEP = ZZ && ZZ->p() > 1. && ZZ->passFullSelection();
	if(genLEP)
		theHistograms->fill("Signal_LEP", "Signal (LEP)", 1,0.,1., 0.5, theWeight);
	if(recLEP)
		theHistograms->fill("Baseline_LEP", "Baseline (LEP)", 1,0.,1., 0.5, theWeight);
	if(genLEP && recLEP)
		theHistograms->fill("Sig_Base_LEP", "Signal + Baseline (LEP)", 1,0.,1., 0.5, theWeight);
	
	if(genLEP && genHAD)
		theHistograms->fill("Signal_BOTH", "Signal (LEP + HAD)", 1,0.,1., 0.5, theWeight);
	if(recLEP && recHAD)
		theHistograms->fill("Baseline_BOTH", "Baseline (LEP + HAD)", 1,0.,1., 0.5, theWeight);
	if(genLEP && recLEP && genHAD && recHAD)
		theHistograms->fill("Sig_Base_BOTH", "Signal + Baseline (LEP + HAD)", 1,0.,1., 0.5, theWeight);
	
	/*if( qq_.p() > 1. ){// qq_.daughter(0).pt() > 30 && qq_.daughter(1).pt() > 30){
		theHistograms->fill<TH2F>("DeltaR vs eta qq", "#DeltaR(qq) vs |#eta(qq)|;#DeltaR;|#eta|", 14,0.,3.5, 10,0.,5., physmath::deltaR(qq_.daughter(0), qq_.daughter(1)), fabs(qq_.eta()), theWeight);
		genQuarksAnalisys();
	}
	if(qq_.p() > 1. && genZZ_->p() > 1.){
		theHistograms->fill<TH2F>("genZZ pt vs qq pt", ";qq p_{T} [GeV/c];ZZ_{GEN} p_{T} [GeV/c]", 25,0.,250., 25,0.,250., qq_.pt(), genZZ_->pt(), theWeight);
	}*/
}


void VZZAnalyzer::ZZRecoEfficiency(){
	// ----- Efficiency of ZZ -----
	if(genZZ_ && genZZ_->p() > 1.){
		theHistograms->fill("ZZgen: genZZpt","ZZgen: genZZpt", 25,0.,250., genZZ_->pt(), theWeight);
		float max_pt_lep = std::max({genZZ_->first().daughter(0).pt(), genZZ_->first().daughter(1).pt(), genZZ_->second().daughter(0).pt(), genZZ_->second().daughter(0).pt()});
		theHistograms->fill("ZZgen: max_pt_lep","ZZgen: max_pt_lep", 25,0.,250., max_pt_lep, theWeight);
		float min_pt_lep = std::min({genZZ_->first().daughter(0).pt(), genZZ_->first().daughter(1).pt(), genZZ_->second().daughter(0).pt(), genZZ_->second().daughter(0).pt()});
		theHistograms->fill("ZZgen: min_pt_lep","ZZgen: min_pt_lep", 25,0.,80., min_pt_lep, theWeight);
		
		if(ZZ && ZZ->p() > 1.){
			theHistograms->fill("ZZrec: genZZpt","ZZrec: genZZpt", 25,0.,250., genZZ_->pt(), theWeight);
			theHistograms->fill("ZZrec: max_pt_lep","ZZrec: max_pt_lep", 25,0.,250., max_pt_lep, theWeight);
			theHistograms->fill("ZZrec: min_pt_lep","ZZrec: min_pt_lep", 25,0.,80., min_pt_lep, theWeight);
			
			theHistograms->fill("ZZ mass res","ZZ mass res", 40,-15.,15., genZZ_->mass() - ZZ->mass(), theWeight);
			theHistograms->fill("Z mass res","Z mass res", 40,-15.,15., genZZ_->first().mass() - ZZ->first().mass(), theWeight);
			theHistograms->fill("Z mass res","Z mass res", 40,-15.,15., genZZ_->first().mass() - ZZ->first().mass(), theWeight);
		}
	}
}


void VZZAnalyzer::jetRecoEfficiency(){
	// ----- Efficiency of AK4 jets -----
	std::list<Jet>     rAK4(jets->begin(), jets->end());
	
	foreach(const Particle& gen, *genJets){
		if(gen.pt() < 30) continue;
		theHistograms->fill("effJ genAK4: eta", "genAK4: eta;#eta;# jets", ETA_SIZE, gen.eta(), theWeight);
		theHistograms->fill("effJ genAK4: pt", "genAK4: pt;p_{T} [Gev/c];# jets", PT_SIZE, gen.pt(), theWeight);
		
		for(auto rec = rAK4.begin(); rec != rAK4.end(); ++rec){
			if(physmath::deltaR(gen, *rec) > 0.4)
				continue;
			
			theHistograms->fill("effJ recAK4: eta", "recAK4: gen eta;#eta;# jets", ETA_SIZE, gen.eta(), theWeight);
			theHistograms->fill("effJ recAK4: pt", "recAK4: gen pt;p_{T} [Gev/c];# jets", PT_SIZE, gen.pt(), theWeight);
			rAK4.erase(rec);
			break;
		}
	}
	
	// ----- Efficiency of AK8 jets -----
	std::list<Jet>     rAK8(jetsAK8->begin(), jetsAK8->end());
	
	foreach(const Particle& gen, *genJetsAK8){
		if(gen.pt() < 150) continue;
		theHistograms->fill("effJ genAK8: eta", "genAK8: eta;#eta;# jets", ETA_SIZE, gen.eta(), theWeight);
		theHistograms->fill("effJ genAK8: pt", "genAK8: pt;p_{T} [Gev/c];# jets", PT_SIZE, gen.pt(), theWeight);
		
		for(auto rec = rAK8.begin(); rec != rAK8.end(); ++rec){
			if(physmath::deltaR(gen, *rec) > 0.4)
				continue;
			
			theHistograms->fill("effJ recAK8: eta", "genAK8: gen eta;#eta;# jets", ETA_SIZE, gen.eta(), theWeight);
			theHistograms->fill("effJ recAK8: pt", "genAK8: gen pt;p_{T} [Gev/c];# jets", PT_SIZE, gen.pt(), theWeight);
			rAK8.erase(rec);
			break;
		}
	}
}


void VZZAnalyzer::endNameHistos(){
	TH1* massAlgorithms = theHistograms->get("Resolution AK8: winner");
	if(massAlgorithms)
		for(unsigned int i = 0; i < 6; ++i)
			massAlgorithms->GetXaxis()->SetBinLabel(i+1, massAlgsNames_[i]);
	
	TH1* sigType = theHistograms->get("Signal Type");
	if(sigType){
		sigType->GetXaxis()->SetBinLabel(1, "OK");
		sigType->GetXaxis()->SetBinLabel(2, "Problem");
		sigType->GetXaxis()->SetBinLabel(3, "No qq");
		sigType->GetXaxis()->SetBinLabel(4, "genAK4 < 2");
		sigType->GetXaxis()->SetBinLabel(5, "no distinct AK4");
		sigType->GetXaxis()->SetBinLabel(6, "mass range");
	}
	
	TH1* top_u = theHistograms->get("topology_u");
	if(top_u){
		top_u->GetXaxis()->SetBinLabel(1, "0");
		top_u->GetXaxis()->SetBinLabel(2, "4");
		top_u->GetXaxis()->SetBinLabel(3, "6");
	}
	
	TH1* gCutsG = theHistograms->get("Cuts gen_f_");
	if(gCutsG){
		gCutsG->GetXaxis()->SetBinLabel(1, "No cut");
		gCutsG->GetXaxis()->SetBinLabel(2, "Lep signal");
		gCutsG->GetXaxis()->SetBinLabel(3, "Had signal");
	}
	
	TH1* gCutsR = theHistograms->get("Cuts rec_f_");
	if(gCutsR){
		gCutsR->GetXaxis()->SetBinLabel(1, "No cut");
		gCutsR->GetXaxis()->SetBinLabel(2, "Lep signal");
		gCutsR->GetXaxis()->SetBinLabel(3, "Had signal");
	}
	
	TH1* gCutsAnalysis = theHistograms->get("Cut analysis_f_");
	if(gCutsAnalysis){
		gCutsAnalysis->GetXaxis()->SetBinLabel(1, "Baseline");
		gCutsAnalysis->GetXaxis()->SetBinLabel(2, "Any minDM < 30");
		gCutsAnalysis->GetXaxis()->SetBinLabel(3, "AK4 minDM < 30");
		gCutsAnalysis->GetXaxis()->SetBinLabel(4, "Any minDM < 13");
	}
}


void VZZAnalyzer::genHadVBCategorization(){
	if(!topology.test(0)) return;
	const char name[] = "GenHad VBs: category";
	//None
	if(!topology.test(8) && !topology.test(9))
		theHistograms->fill(name, name, 4,-0.5,3.5, 0., theWeight);
	//Only Wjj
	else if(topology.test(8) && !topology.test(9))
		theHistograms->fill(name, name, 4,-0.5,3.5, 1., theWeight);
	//Only Zjj
	else if(!topology.test(8) && topology.test(9))
		theHistograms->fill(name, name, 4,-0.5,3.5, 2., theWeight);
	//Both
	else if(topology.test(8) && topology.test(9)){
		theHistograms->fill(name, name, 4,-0.5,3.5, 3., theWeight);
		
		//Are they made from the same jets?
		const Boson<Particle>& V1 = genHadVBs_->at(0);
		const Boson<Particle>& V2 = genHadVBs_->at(1);
		float dP00 = (V1.daughter(0).p4() - V2.daughter(0).p4()).P();
		float dP01 = (V1.daughter(0).p4() - V2.daughter(1).p4()).P();
		float dP10 = (V1.daughter(1).p4() - V2.daughter(0).p4()).P();
		float dP11 = (V1.daughter(1).p4() - V2.daughter(1).p4()).P();
		float dP0 = std::min(dP00, dP01);  //minimum distance of V1(0) from any daugt. of V2
		float dP1 = std::min(dP10, dP11);  //                    V1(1)  
		const char name2[] = "GenHadVBs(2): # jets in common";
		
		if(dP0 < 1. && dP1 < 1.){       // 2 matches
			theHistograms->fill(name2, name2, 3,-0.5,2.5, 2., theWeight);
			float mass = genHadVBs_->front().mass();
			char val = ( fabs(mass - phys::WMASS) < fabs(mass - phys::ZMASS) ? 0 : 1 );
			theHistograms->fill<TH2F>("GenHadVBs(2): mass", "GenHadVBs(2): mass;m [GeV/c^{2}]", 18,50.,122., 2,-0.5,1.5, mass, val, theWeight);
		}
		else if(dP0 < 1. || dP1 < 1.)  // 1 match
			theHistograms->fill(name2, name2, 3,-0.5,2.5, 1., theWeight);
		else                           // 0 matches
			theHistograms->fill(name2, name2, 3,-0.5,2.5, 0., theWeight);
	}
}


void VZZAnalyzer::endGenHadVbCateg(){
	TH1* hist = theHistograms->get("GenHad VBs: category");
	if(hist){
		hist->GetXaxis()->SetBinLabel(1, "No Vjj");
		hist->GetXaxis()->SetBinLabel(2, "Wjj");
		hist->GetXaxis()->SetBinLabel(3, "Zjj");
		hist->GetXaxis()->SetBinLabel(4, "Both");
	}
	
	TH1* hist2 = theHistograms->get("GenHadVBs(2): mass");
	if(hist2){
		hist2->GetYaxis()->SetBinLabel(1, "W");
		hist2->GetYaxis()->SetBinLabel(2, "Z");
	}
}


void VZZAnalyzer::AK8recover(){
	++NnoAK4_;  // There wasn't a V->jj
	if(genHadVBs_->size() == 0) return;
	++N8genVB_;  // Number of events in which there's a VB->had, not found in V->jj
	if(jetsAK8->size() == 0)    return;
	
	vector<Jet> jetsAK8pTau35;
	foreach(const Jet& j, *jetsAK8)
		if(j.puppiTau2()/j.puppiTau1() < 0.35) 
			jetsAK8pTau35.push_back(j);
	
	// Search for a V-->J
	foreach(const Boson<Particle>& genVB, *genHadVBs_){  // maybe AllGenVBjj_?
		auto closestAK8 = std::min_element(jetsAK8->begin(), jetsAK8->end(), DeltaRComparator(genVB));
		float dR = physmath::deltaR(*closestAK8, genVB);
		if(dR < 0.8)
			theHistograms->fill("Extra AK8: #DeltaR","Extra AK8: #DeltaR", 40,0.,0.8, dR, theWeight);
		if(dR < 0.2){
			++N8recover_;
			theHistograms->fill("Extra AK8: genVB pt","Extra AK8: genVB pt", 40,0.,400, genVB.pt(), theWeight);
			float dM = closestAK8->corrPrunedMass() - genVB.mass();
			theHistograms->fill("Extra AK8: #DeltaMcorrPruned","Extra AK8: #DeltaMcorrPruned", 26,-40.,25., dM, theWeight);
			if(ZZ)
				theHistograms->fill("Extra AK8: ZZ pt","Extra AK8: ZZ pt", 40,0.,400, ZZ->pt(), theWeight);
		}
		
		// PuppiTau21 < 0.35
		if(jetsAK8pTau35.size() == 0) continue;
		auto closestAK8pTau35 = std::min_element(jetsAK8pTau35.begin(), jetsAK8pTau35.end(), DeltaRComparator(genVB));
		float dRpt = physmath::deltaR(*closestAK8pTau35, genVB);
		if(dRpt < 0.8)
			theHistograms->fill("Extra AK8_{pTau<.35}: #DeltaR", "Extra AK8_{pTau<.35}: #DeltaR", 40,0.,0.8, dRpt, theWeight);
		if(dRpt < 0.2){
			theHistograms->fill("Extra AK8_{pTau<.35}: genVB pt", "Extra AK8_{pTau<.35}: genVB pt", 40,0.,400, genVB.pt(), theWeight);
			float dM = closestAK8pTau35->corrPrunedMass() - genVB.mass();
			theHistograms->fill("Extra AK8_{pTau<.35}: #DeltaMcorrPruned", "Extra AK8_{pTau<.35}: #DeltaMcorrPruned", 26,-40.,25., dM, theWeight);
			if(ZZ)
				theHistograms->fill("Extra AK8_{pTau<.35}: ZZ pt", "Extra AK8_{pTau<.35}: ZZ pt", 40,0.,400, ZZ->pt(), theWeight);
		}
	}
}


void VZZAnalyzer::genTauAnalisys(){
	vector<Jet> rec4fromTau;
	if(genTaus->size() == 0){
		theHistograms->fill("Tau daughters", "Number of #tau daughters of a genHadVB", 3,-0.5,2.5, 0., theWeight * genHadVBs_->size()*2); //None of the 2 daughters of the VB is a tau
		if(topology.test(0))
			theHistograms->fill("Tau daughters topology0", "Number of #tau daughters of a genHadVB with topology(0)==true", 3,-0.5,2.5, 0., theWeight * genHadVBs_->size()*2);
		return;  // Total number of genAK4 and recAK4 can be read from other graphs
	}
	
	foreach(const Particle& tau, *genTaus)
		theHistograms->fill("genTaus #eta", "genTaus #eta", 25,-5.,5., tau.eta(), 1.);
	
	foreach(const Boson<Particle>& hadVB, *genHadVBs_){
		unsigned char tau_daughters = 0;
		for(int d = 0; d <= 1; ++d){ //Loop on VB daughters
			//theHistograms->fill("Closest #tau: AK4_{VB} #eta", "Closest #tau: AK4_{VB} #eta", 25,-5.,5., hadVB.daughter(d).eta(), 1.);
			std::sort(genTaus->begin(), genTaus->end(), DeltaRComparator(hadVB.daughter(d)));
			const Particle& tau = genTaus->front();
			float dR_gen = physmath::deltaR(hadVB.daughter(d), tau);
			if(dR_gen < 0.2){
				++tau_daughters;
				theHistograms->fill("Closest #tau: #DeltaR(AK4_{VB})", "Closest #tau: #DeltaR(AK4_{VB})", 20,0.,0.2, dR_gen, 1.);
				theHistograms->fill("Closest #tau: match4_{VB} #tau #eta", "Closest #tau: match4_{VB} #tau #eta", 25,-5.,5., tau.eta(), 1.);
				//float dM = hadVB.daughter(d).mass() - tau.mass();
				//theHistograms->fill("Closest #tau: #DeltaM(AK4_{VB})", "Closest #tau: #DeltaM(AK4_{VB})", 20,-25.,25., dM, 1.);
			
				// Reconstructed AK4
				if(jets->size() == 0) continue;
				std::sort(jets->begin(), jets->end(), DeltaRComparator(hadVB.daughter(d)));
				float dR_genrec = physmath::deltaR(hadVB.daughter(d), jets->front());
				if(dR_genrec < 0.2){
					//theHistograms->fill("Closest #tau: #DeltaR(AK4_{VB}, AK4_{rec})", "Closest #tau: #DeltaR(AK4_{VB}, AK4_{rec})", 20,0.,0.4, dR_genrec, 1.);
					//theHistograms->fill("Closest #tau: AK4_{rec} #eta", "Closest #tau: AK4_{rec} #eta", 25,-5.,5., hadVB.daughter(0).eta(), 1.);
					
					float dR_rec = physmath::deltaR(tau, jets->front());
					if(dR_rec < 0.2){
						//The same tau is close to the jet that reconstructs the genJet from the VB
						theHistograms->fill("Closest #tau: #DeltaR(AK4_{rec})", "Closest #tau: #DeltaR(AK4_{rec})", 20,0.,0.2, dR_rec, 1.);
						theHistograms->fill("Closest #tau: match4_{rec} #tau #eta", "Closest #tau: match4_{rec} #tau #eta", 25,-5.,5., tau.eta(), 1.);
					}
				}
			}
		}// End Loop on VB daughters
		theHistograms->fill("Tau daughters", "Number of #tau daughters of a genHadVB", 3,-0.5,2.5, tau_daughters, theWeight);
		if(topology.test(0))
			theHistograms->fill("Tau daughters topology0", "Number of #tau daughters of a genHadVB with topology(0)==true", 3,-0.5,2.5, tau_daughters, theWeight);
	} 
	
	
	//Search for a reconstructed AK4 near the tau, indipendently of the genVB
	foreach(const Particle& tau, *genTaus){  //foreach(const Jet& jet, *jets){
		//theHistograms->fill("Closest #tau: All AK4_{rec} #eta", "Closest #tau: All AK4_{rec} #eta", 25,-5.,5., jet.eta(), 1.);
		stable_sort(jets->begin(), jets->end(), DeltaRComparator(tau));
		const Jet& jet = jets->front(); //const Particle& tau = genTaus->front();
		if(physmath::deltaR(jet, tau) < 0.2){
			theHistograms->fill("Closest #tau: #DeltaR(All AK4_{rec})", "Closest #tau: #DeltaR(All AK4_{rec})", 25,0.,0.5, physmath::deltaR(jet, tau), 1.);
			theHistograms->fill("Closest #tau: match4_{All rec} #tau #eta", "Closest #tau: match4_{All rec} #tau #eta", 25,-5.,5., tau.eta(), 1.);
			float dMr = jet.mass() - tau.mass();
			theHistograms->fill("Closest #tau: #DeltaM(All AK4_{rec})", "Closest #tau: #DeltaM(All AK4_{rec})", 20,-25.,25., dMr, 1.);
		}
	}
	
	
	// AK8
	/*vector<Particle> gen8fromTaus;
	vector<Jet> rec8fromTaus;
	for(size_t i = 0; i < genTaus->size(); ++i){
		for(size_t j = i+1; j < genTaus->size(); ++j){
			TLorentzVector p4Taus = genTaus->at(i).p4() + genTaus->at(j).p4();
			//Generated
			foreach(const Particle& jet8, *genJetsAK8)
				if(physmath::deltaR(jet8.p4(), p4Taus) < 0.5)
					gen8fromTaus.push_back( jet8 );
			
			//Reconstructed
			foreach(const Jet& jet8, *jetsAK8)
				if(physmath::deltaR(jet8.p4(), p4Taus) < 0.5)
					rec8fromTaus.push_back( jet8 );
		}
	}*/
	
	//Combinatory Zs
	//vector<Boson<Particle>> genZtoTau;
	for(size_t i = 0; i < genTaus->size(); ++i){
		for(size_t j = i+1; j < genTaus->size(); ++j){
			TLorentzVector p4tau = genTaus->at(i).p4() + genTaus->at(j).p4();
			float mass = p4tau.M();
			if(60 < mass && mass < 120){
				//genZtoTau.push_back(Boson<Particle>(genTaus->at(i), genTaus->at(j)));
				theHistograms->fill("Combinatory Z (#tau): mass", "Combinatory Z (#tau): mass;[GeV/c^{2}]", 24,60.,120., mass, 1.);
				theHistograms->fill("Combinatory Z (#tau): pt", "Combinatory Z (#tau): pt;[GeV/c]", 20,0.,400., p4tau.Pt(), 1.);
			}
		}
	}
	
	//vector<Boson<Jet>> genZtoTauTo4;
	for(size_t i = 0; i < rec4fromTau.size(); ++i){
		for(size_t j = i+1; j < rec4fromTau.size(); ++j){
			TLorentzVector p4tau = rec4fromTau.at(i).p4() + rec4fromTau.at(j).p4();
			float mass = p4tau.M();
			if(60 < mass && mass < 120){
				//genZtoTauTo4.push_back(Boson<Jet>(rec4fromTau.at(i), rec4fromTau.at(j)));
				theHistograms->fill("Combinatory Z (#tau-->AK4): mass", "Combinatory Z (#tau-->AK4): mass;[GeV/c^{2}]", 24,60.,120., mass, 1.);
				theHistograms->fill("Combinatory Z (#tau-->AK4): pt", "Combinatory Z (#tau-->AK4): pt;[GeV/c]", 20,0.,400., p4tau.Pt(), 1.);
			}
		}
	}
	
}


void VZZAnalyzer::AK8nearGenHadVB(){
	if(genHadVBs_->size() == 0)
		return;
	
	if(genJetsAK8->size() > 0){
		std::sort(genJetsAK8->begin(), genJetsAK8->end(), DeltaRComparator(genHadVBs_->front()));
		float dRg = physmath::deltaR(genJetsAK8->front(), genHadVBs_->front());
		if(dRg < 1.){
			theHistograms->fill("genHadVB0 closest genAK8: deltaR", "#DeltaR between genHadVB(0) and the closest gen AK8", 25,0.,1., dRg, theWeight);
			float dMg = genJetsAK8->front().mass() - genHadVBs_->front().mass();
			theHistograms->fill("genHadVB0 closest genAK8: deltaM", "mass (closest gen AK8) - mass(genHadVB(0))", 25,-50.,50., dMg, theWeight);
		}
	}
	
	if(jetsAK8->size() > 0){
		std::sort( jetsAK8->begin(), jetsAK8->end(), DeltaRComparator(genHadVBs_->front()) );
		float dRr = physmath::deltaR(jetsAK8->front(), genHadVBs_->front());
		if(dRr < 1.){
			theHistograms->fill("genHadVB0 closest recAK8: deltaR", "#DeltaR between genHadVB(0) and the closest rec AK8", 25,0.,1., dRr, theWeight);
			float dMr = jetsAK8->front().mass() - genHadVBs_->front().mass();
			theHistograms->fill("genHadVB0 closest recAK8: deltaM", "mass (closest rec AK8) - mass(genHadVB(0))", 25,-50.,50., dMr, theWeight);
		}
	}
}


void VZZAnalyzer::AK8MassAlgorithms(){
	if(jetsAK8->size() == 0 || genJetsAK8->size() == 0) return;
	
	
	// Define a more precise signal region
	vector<Jet> jetsAK8_puppiTauCut;
	foreach(const Jet& jet, *jetsAK8)
		if(jet.puppiTau2()/jet.puppiTau1() < 0.35) 
			jetsAK8_puppiTauCut.push_back(jet);
	if(jetsAK8_puppiTauCut.size() == 0) return;
	
	
	foreach(const Particle& genj, *genJetsAK8){
		if(genj.mass() < 60. || genj.mass() > 120.) continue;
		stable_sort(jetsAK8_puppiTauCut.begin(), jetsAK8_puppiTauCut.end(), phys::DeltaRComparator(genj)); //DeltaRComparator is defined in Commons/interface/Utils.h
		float dR = physmath::deltaR(jetsAK8_puppiTauCut.front(), genj);
		if(dR < 0.5 /*or 0.8?*/){
			theHistograms->fill("Resolution AK8: #DeltaR", "Resolution AK8: #DeltaR", 50,0.,0.5, dR, theWeight);
			//We have a match! This AK8 should have the same mass of the gen one
			float mass           = jetsAK8_puppiTauCut.front().mass();
			float secvtxMass     = jetsAK8_puppiTauCut.front().secvtxMass();
			float corrPrunedMass = jetsAK8_puppiTauCut.front().corrPrunedMass();
			float prunedMass     = jetsAK8_puppiTauCut.front().prunedMass();
			float softDropMass   = jetsAK8_puppiTauCut.front().softDropMass();
			float puppiMass      = jetsAK8_puppiTauCut.front().puppiMass();
			float massesVal[6] = {mass, secvtxMass, corrPrunedMass, prunedMass, softDropMass, puppiMass};
			float absDiffs[6];
			for(unsigned int i = 0; i<6; ++i){
				char* name = Form("Resolution AK8: #DeltaM with %s", massAlgsNames_[i]);
				char* title = Form("%s;#DeltaM [GeV/c^{2}]", name);
				theHistograms->fill(name, title, 40,-50.,50., massesVal[i] - genj.mass(), theWeight);
				absDiffs[i] = fabs( massesVal[i] - genj.mass() );
			}
			//find the minimum difference
			int minpos = std::min_element(absDiffs, absDiffs + 6) - absDiffs; //pointer arithmetic
			theHistograms->fill("Resolution AK8: winner", "Resolution AK8: winner;Algorithm", 6,-0.5,5.5, minpos, theWeight);
		}
	}
}


pair<const Particle*, Jet*> VZZAnalyzer::reconstructionAK8(/*bool doGraphs*/){
	if(jetsAK8->size() == 0 || genJetsAK8->size() == 0) 
		return make_pair(nullptr, nullptr);
	
	// Define a more precise signal region (60<m_gen<120 && puppiTau21_rec<0.35)
	vector<Jet> jetsAK8_puppiTauCut;
	foreach(const Jet& jet, *jetsAK8)
		if(jet.puppiTau2()/jet.puppiTau1() < 0.35) 
			jetsAK8_puppiTauCut.push_back(jet);
	
	// Reconstruction algorithm
	const Jet* recAK8cand = findBestVFromSing(&jetsAK8_puppiTauCut);
	if(recAK8cand) ++Nall_rec8_pTau35_;  // rec by alg (60<m_gen<120 && puppiTau21_rec<0.35)
	
	const Particle* theGen = nullptr;
	foreach(const Particle& genj, *genJetsAK8){
		if(genj.mass() < 60. || genj.mass() > 120.)
			continue;
		N_gen8_++;
		//theHistograms->fill("AK8_gen in 60-120: pt_gen", "AK8_gen in 60-120: pt", 25,0.,500., genj.pt(), theWeight);
		
		// Detector efficiency in the loose signal region (60<m_gen<120)
		stable_sort(jetsAK8->begin(), jetsAK8->end(), phys::DeltaRComparator(genj));
		if(physmath::deltaR(jetsAK8->front(), genj) < 0.5){
			++Ncms_rec8_;
			//theHistograms->fill("AK8_rec in 60-120: pt_gen", "AK8_rec in 60-120: pt_gen", 25,0.,500., genj.pt(), theWeight);
		}
	
		// Algorithm in tight signal region (60<m_gen<120 && puppiTau21_rec<0.35)
		if(jetsAK8_puppiTauCut.size() == 0) continue;
		stable_sort(jetsAK8_puppiTauCut.begin(), jetsAK8_puppiTauCut.end(), phys::DeltaRComparator(genj)); //DeltaRComparator is defined in Commons/interface/Utils.h
		float dR = physmath::deltaR(jetsAK8_puppiTauCut.front(), genj);
		if(dR < 0.5 /*or 0.8?*/){
			++Ncms_rec8_pTau35_;  // Reconstructed by the detector in this particular region
			//theHistograms->fill("AK8_rec in 60-120, ptau21<.35: pt_gen", "AK8_rec in 60-120, ptau21<.35: pt_gen", 25,0.,500., genj.pt(), theWeight);
			if(recAK8cand && physmath::deltaR(jetsAK8_puppiTauCut.front(), *recAK8cand) < 0.1){
				++Ngood_rec8_pTau35_;
				theGen = &genj;
			}
		}
	}
	
	if(jetsAK8_puppiTauCut.size() == 0) 
		return make_pair(nullptr, nullptr);
	
	Jet* candCopy = nullptr;
	if(recAK8cand)
		candCopy = new Jet(*recAK8cand); //recAK8cand is part of a temporary vector, it would be lost when the program returns from this function
	return make_pair(theGen, candCopy);
}


void VZZAnalyzer::makeQQ(bool doGraphs){
	vector<Particle> genQuarks;
	vector<Particle> genLeptons;
	foreach(const Particle& p, *genParticles){
		unsigned int aID = abs(p.id());
		if(doGraphs)
			theHistograms->fill("particle ID", "particle ID", 61,-30.5,30.5, p.id(), theWeight);
		
		if(aID < 10){       // Is it a quark?
			genQuarks.push_back(Particle(p));
			//theHistograms->fill("quark ID", "quark ID", 21,-10.5,10.5, p.id(), theWeight);
			if(doGraphs){
				theHistograms->fill("quark charge", "quark charge", 7,-7./6.,7./6., p.charge(), theWeight);
				theHistograms->fill("quark mother ID", "quark mother ID", 61,-30.5,30.5, p.motherId(), theWeight);
			}
			//if( p.motherId() == 23 || abs(p.motherId()) == 24 )  // Is a daugther of a VB?
		}
		else if(aID < 20){  // Is it a lepton?
			genLeptons.push_back(Particle(p));
			//theHistograms->fill("lepton ID", "lepton ID", 41,-20.5,20.5, p.id(), theWeight);
		}
		//else
			//theHistograms->fill("Extra ID", "Extra ID", 201,-100.5,100.5, p.id(), theWeight);
	}
	
	if(doGraphs){
		theHistograms->fill("genLeptons and genQuarks", "Number of genLeptons and genQuarks;genLeptons;genQuarks", 7,-0.5,6.5, 7,-0.5,6.5, genLeptons.size(), genQuarks.size(), theWeight);
		theHistograms->fill("Type of event", "Type of event (all, 4l, 2q, 4l+2q)", 4,-0.5,3.5, 0., theWeight);
		if(genLeptons.size() == 4)
			theHistograms->fill("Type of event", "Type of event (all, 4l, 2q, 4l+2q)", 4,-0.5,3.5, 1., theWeight);
		if(genQuarks.size() == 2)
			theHistograms->fill("Type of event", "Type of event (all, 4l, 2q, 4l+2q)", 4,-0.5,3.5, 2., theWeight);
		if(genLeptons.size() == 4 && genQuarks.size() == 2)
			theHistograms->fill("Type of event", "Type of event (all, 4l, 2q, 4l+2q)", 4,-0.5,3.5, 3., theWeight);
	}
	
	// 4l 2q: my final state
	if(genLeptons.size() == 4 && genQuarks.size() == 2){
		float mqq = ( genQuarks.at(0).p4() + genQuarks.at(1).p4() ).M();
		float chqq = genQuarks.at(0).charge() + genQuarks.at(1).charge();
		
		if(doGraphs){
			bool sfos = genQuarks.at(0).id() + genQuarks.at(1).id() == 0;
			bool ch1 = fabs(chqq + 1) < 0.3 || fabs(chqq - 1) < 0.3;
			
			if( sfos )  //from Z?
				theHistograms->fill("qq SFOS mass","SFOS;GeV/c^{2}", 25,60.,110., mqq,theWeight);
			if( ch1 )  //from W?
				theHistograms->fill("qq ch1 mass", "ch = #pm1;GeV/c^{2}", 25,60.,110., mqq, theWeight);
			if(!sfos && !ch1){
				theHistograms->fill("qq extra mass", "extra", 25,60.,110., mqq, theWeight);
				theHistograms->fill("qq extra charge", "extra", 9,-1.5,1.5, chqq, theWeight);
			}
		}
		//return new Boson<Particle>(genQuarks.at(0), genQuarks.at(1));
		qq_ = Boson<Particle>(genQuarks.at(0), genQuarks.at(1));
	}
	//else
		//new (qq_) Boson<Particle>(); //return nullptr;	
}


void VZZAnalyzer::genQuarksAnalisys(){
	// ----- Efficiency, resolution, ecc. -----
	float dRqq = physmath::deltaR(qq_.daughter(0), qq_.daughter(1));
	theHistograms->fill("deltaR quark", "#DeltaR quark;#DeltaR", 25,-0.,5., dRqq, theWeight);  // Denominator of "efficiency AK4/8 vs dR(q,q)"
	theHistograms->fill("pt quarks", "pt quarks;pt[GeV/c]", 25,0.,500., qq_.pt(), theWeight);  // Denominator of "efficiency AK4/8 vs pt(qq)"
	
	// Best mass pair of gen AK4
	phys::Boson<phys::Particle>* candGen4 = findBestVFromPair(genJets);
	if(candGen4)
		theHistograms->fill("pt quarks: AK4s algo", "pt quarks (genAK4s found);pt[GeV/c]", 25,0.,500., qq_.pt(), theWeight);  // Denominator of "Purity genAK4s vs pt(qq)"
	
	// GenJets AK4  IMPORTANT: eff/res ONLY for the 4l 2q channel (no 4q/6q)!
	Boson<Particle>* existingGen4 = nullptr;
	theHistograms->fill("Eta quark", "#eta quark;#eta", 25,-5.,5., qq_.daughter(0).eta(), theWeight);  // Denominator of "efficiency AK4 vs eta(q)"
	theHistograms->fill("Eta quark", "#eta quark;#eta", 25,-5.,5., qq_.daughter(1).eta(), theWeight);  // Denominator of "efficiency AK4 vs eta(q)"
	
	if(genJets->size() > 1){
		vector<Particle> cpGenJets(*genJets);
		size_t i4_0 = 99;
		const Particle* gen4_0 = closestSing(genJets, qq_.daughter(0), i4_0);
		size_t i4_1 = 99;
		const Particle* gen4_1 = closestSing(genJets, qq_.daughter(1), i4_1);
		
		if(gen4_0 && gen4_1){
			float dR0 = physmath::deltaR(*gen4_0, qq_.daughter(0));
			float dR1 = physmath::deltaR(*gen4_1, qq_.daughter(1));
			
			if(physmath::deltaR(*gen4_0, *gen4_1) < 0.01){ // If they're matched to the same jet
				if(dR0 < dR1){ // Search for a second, different jet
					cpGenJets.erase(cpGenJets.begin()+i4_0);
					gen4_1 = closestSing(&cpGenJets, qq_.daughter(1), i4_1);
				} else if(dR1 < dR0){
					cpGenJets.erase(cpGenJets.begin()+i4_1);
					gen4_0 = closestSing(&cpGenJets, qq_.daughter(1), i4_0);
				}
			}
		}
		
		if(gen4_0){
			theHistograms->fill("Resolution q-genAK4 dR", "Resolution q-genAK4: #DeltaR;#DeltaR", 20,0.,0.4, physmath::deltaR(*gen4_0, qq_.daughter(0)), theWeight);
			theHistograms->fill("Eta quark: AK4", "#eta quark (genAK4 exists);#eta", 25,-5.,5., qq_.daughter(0).eta(), theWeight);  // Numerator of "efficiency AK4 vs eta(q)"
		}
		if(gen4_1){
			theHistograms->fill("Resolution q-genAK4 dR", "Resolution q-genAK4: #DeltaR;#DeltaR", 20,0.,0.4, physmath::deltaR(*gen4_1, qq_.daughter(1)), theWeight);
			theHistograms->fill("Eta quark: AK4", "#eta quark (genAK4 exists);#eta", 25,-5.,5., qq_.daughter(1).eta(), theWeight);  // Numerator of "efficiency AK4 vs eta(q)"
		}
		
		if(gen4_0 && gen4_1){ // Both valid AFTER checking they are not the same jet
			existingGen4 = new Boson<Particle>(*gen4_0, *gen4_1);
			float dR = physmath::deltaR(gen4_0->p4() + gen4_1->p4(), qq_.p4());
			float dM = (gen4_0->p4() + gen4_1->p4()).M() - qq_.mass();
			
			theHistograms->fill<TH2F>("Resolution qq-genAK4s", "Resolution genAK4s-qq: #DeltaR and #Deltam;#DeltaR;#Deltam [GeV/c^2]", 20,0.,0.4, 10,-50.,50., dR, dM, theWeight);
			theHistograms->fill("Resolution qq-genAK4s mass", "Resolution: #Deltam genAK4s-qq;#Deltam [GeV/c^2]", 25,-50.,50., dM, theWeight);
			theHistograms->fill("deltaR quark: AK4s", "#DeltaR quark (genAK4s exist);#DeltaR", 25,-0.,5., dRqq, theWeight);  // Numerator of "efficiency AK4 vs dR(q,q)"
			theHistograms->fill("pt quarks: AK4s", "pt quarks (genAK4s exist);pt[GeV/c]", 25,0.,500., qq_.pt(), theWeight);  // Numerator of "efficiency AK4 vs pt(qq)"
		}
		
		if(existingGen4 && candGen4){
			float dR00 = physmath::deltaR(*gen4_0, candGen4->daughter(0));
			float dR01 = physmath::deltaR(*gen4_0, candGen4->daughter(1));
			float dR10 = physmath::deltaR(*gen4_1, candGen4->daughter(0));
			float dR11 = physmath::deltaR(*gen4_1, candGen4->daughter(1));
			theHistograms->fill("pt quarks: AK4s found", "pt quarks (genAK4s exist and found);pt[GeV/c]", 25,0.,500., qq_.pt(), theWeight); // Numerator in "eff choice (loose)"
			
			if(std::min(dR00, dR01) < 0.01 && std::min(dR10, dR11) < 0.01){ // chosen the same jets?
				theHistograms->fill("pt quarks: AK4s found good", "pt quarks (genAK4s exist and chosen correctly);pt[GeV/c]", 25,0.,500., qq_.pt(), theWeight); // Numerator in "eff choice (strict)"
				
				theHistograms->fill("Res dR genAK4s-qq found good", "Resolution: #DeltaR genAK4s-qq (found good)", 20,0.,0.4, physmath::deltaR(qq_, *candGen4), theWeight);
				theHistograms->fill("Res mass genAK4s-qq found good", "Resolution: #Deltam genAK4s-qq (found good);#Deltam [GeV/c^2]", 25,-50.,50., qq_.mass() - candGen4->mass(), theWeight);
			}
			//else
				//theHistograms->fill("pt quarks: AK4s chosen bad", "pt quarks (genAK4s found, but chosen others);pt[GeV/c]", 25,0.,500., qq_.pt(), theWeight);
		}
	}
	else if(genJets->size() == 1){
		float dR0 = physmath::deltaR(genJets->front(), qq_.daughter(0));
		float dR1 = physmath::deltaR(genJets->front(), qq_.daughter(1));
		float dRmin = (dR0 < dR1 ? dR0 : dR1);
		if(dRmin < 0.4){
			theHistograms->fill("Resolution q-genAK4 dR", "Resolution q-genAK4: #DeltaR;#DeltaR", 20,0.,0.4, dRmin, theWeight);
			theHistograms->fill("Eta quark: AK4", "#eta quark (genAK4 exists);#eta", 25,-5.,5., qq_.daughter(dR0 < dR1 ? 0 : 1).eta(), theWeight);  // Numerator of "efficiency AK4 vs eta(q)"
		}
	}
	
	
	// Best mass gen AK8
	const Particle* candGen8 = findBestVFromSing(genJetsAK8);
	if(candGen8)
		theHistograms->fill("pt quarks: AK8 algo", "pt quarks (genAK8 found);pt[GeV/c]", 25,0.,500., qq_.pt(), theWeight);  // Denominator of "Purity genAK4s vs pt(qq)"
	
	// GenJets AK8  IMPORTANT: eff/res ONLY for the 4l 2q channel (no 4q/6q)!
	const Particle* existingGen8 = nullptr;
	theHistograms->fill("Eta quarks", "#eta quarks;#eta", 25,-5.,5., qq_.eta(), theWeight);  // Denominator of "efficiency AK8 vs eta(qq)"
	
	if(genJetsAK8->size() > 0){
		sort(genJetsAK8->begin(), genJetsAK8->end(), phys::DeltaRComparator(qq_));
		float dR = physmath::deltaR(genJetsAK8->front(), qq_);
		float dM = genJetsAK8->front().mass() - qq_.mass();
		
		if(dR < 0.8){
			if(physmath::deltaR(genJetsAK8->front(), qq_))
			existingGen8 = &(genJetsAK8->front());
			
			theHistograms->fill("Resolution qq-genAK8 dR", "Resolution qq-genAK8: #DeltaR;#DeltaR", 20,0.,0.8, dR, theWeight);
			theHistograms->fill<TH2F>("Resolution qq-genAK8", "Resolution genAK8-qq: #DeltaR and #Deltam;#DeltaR;#Deltam [GeV/c^2]", 10,0.,0.8, 10,-50.,50., dR, dM, theWeight);
			theHistograms->fill("Resolution qq-genAK8 mass", "Resolution: #Deltam genAK8-qq;#Deltam [GeV/c^2]", 25,-50.,50., dM, theWeight);
			
			theHistograms->fill("Eta quarks: AK8", "#eta quarks (genAK8 exist);#eta", 25,-5.,5., qq_.eta(), theWeight);  // Numerator of "efficiency AK8 vs eta(qq)"
			theHistograms->fill("deltaR quark: AK8", "#DeltaR quark (genAK8 exist);#DeltaR", 25,-0.,5., dRqq, theWeight);  // Numerator of "efficiency AK8 vs dR(q,q)"
			theHistograms->fill("pt quarks: AK8", "pt quarks (genAK8 exist);pt [GeV/c]", 25,0.,500., qq_.pt(), theWeight);  // Numerator of "efficiency AK8 vs pt(qq)"
		}
	}
	if(existingGen8 && candGen8){
		theHistograms->fill("pt quarks: AK8 found", "pt quarks (genAK8 exists and found);pt[GeV/c]", 25,0.,500., qq_.pt(), theWeight); // Numerator in "eff choice (loose)"
		
		if(physmath::deltaR(*existingGen8, *candGen8) < 0.01){ // chosen the same jets?
			theHistograms->fill("pt quarks: AK8 found good", "pt quarks (genAK8 exist and chosen correctly);pt[GeV/c]", 25,0.,500., qq_.pt(), theWeight); // Numerator in "eff choice (strict)"
		}
	}
	
	
	// Cleanup
	if(candGen4)     delete candGen4;
	if(existingGen4) delete existingGen4;
}


void VZZAnalyzer::reconstructionAK(){
	pair<const Boson<Particle>*, Boson<Jet>*> matchAK4 = reconstructionAK4();
	pair<const Particle*, Jet*> matchAK8 = reconstructionAK8();
	
	if(matchAK8.second != nullptr && matchAK8.first != nullptr && matchAK4.second != nullptr && matchAK4.first != nullptr)
		AKrace(matchAK4, matchAK8); //Compare the resolution of AK4 and 8 when both have been rec.
	
	if(matchAK4.second == nullptr || matchAK4.first == nullptr)
		AK8recover();
	
	if(matchAK4.second) delete matchAK4.second;  //Cleanup
	if(matchAK8.second) delete matchAK8.second;
}


pair<const Boson<Particle>*, Boson<Jet>*> VZZAnalyzer::reconstructionAK4(){
	// Reconstruction algorithm
	Boson<Jet>* recAK4cand = findBestVFromPair(jets);
	if(recAK4cand){
		++Nall_recVB_60_120_;
		theHistograms->fill("Purity AK4 pairs (all alg-rec): rec mass", "Purity AK4 pairs (all alg-rec): rec mass", 30,60.,120., recAK4cand->mass(), theWeight);
		theHistograms->fill("Purity AK4 pairs (all alg-rec): rec pt", "Purity AK4 pairs (all alg-rec): rec pt", 25,0.,500., recAK4cand->pt(), theWeight);
		theHistograms->fill("Purity AK4 pairs (all alg-rec): # of rec AK4", "Purity AK4 pairs (all alg-rec): # of rec AK4", 11,-0.5,10.5, jets->size(), theWeight);
		if(phys::WMASS - 20. < recAK4cand->mass() && recAK4cand->mass() < phys::ZMASS + 20.)
			++Nall_recVB_m20_;
	}
	
	
	// Gen vs rec matching
	if(genHadVBs_->size() == 0){
		//if(recAK4cand) delete recAK4cand;
		return make_pair(nullptr, recAK4cand);
	}
	const Boson<Particle>* genMatch = nullptr;
	//++Nevt_with_genVB;
	
	foreach(const Boson<Particle>& genVB, *genHadVBs_){
		++Nhad_genVB_;
		theHistograms->fill("GenHadVB: mass", "VB_{gen} --> jj: mass;mass [GeV/c^{2}]", 30,60.,120., genVB.mass(), theWeight);
		theHistograms->fill("GenHadVB: pt", "VB_{gen} --> jj: pt;pt [GeV/c]", 25,0.,500., genVB.pt(), theWeight);
		theHistograms->fill("GenHadVB: # of rec AK4", "VB_{gen} --> jj: # of rec AK4", 11,-0.5,10.5, jets->size(), theWeight);
		
		std::stable_sort( jets->begin(), jets->end(), DeltaRComparator(genVB.daughter(0)) );
		Jet rec0 = jets->front();  // Copy for future use
		float dR0 = physmath::deltaR(rec0, genVB.daughter(0));
		if(dR0 < 0.2)
			theHistograms->fill("Resolution AK4_VB: #DeltaR", "Resolution AK4_VB: #DeltaR", 20,0.,0.2, dR0, theWeight);
		
		std::stable_sort( jets->begin(), jets->end(), DeltaRComparator(genVB.daughter(1)) );
		Jet rec1 = jets->front();  // Copy for future use
		float dR1 = physmath::deltaR(rec1, genVB.daughter(1));
		if(dR1 < 0.2)
			theHistograms->fill("Resolution AK4_VB: #DeltaR", "Resolution AK4_VB: #DeltaR", 20,0.,0.2, dR1, theWeight);
		
		// Gen VB can be correctly reconstructed
		if(dR0 < 0.2 && dR1 < 0.2){
			++Ncms_recVB_;
			Boson<Jet> bestPossibleAK4(rec0, rec1);
			//float dR = physmath::deltaR(bestPossibleAK4, genVB);
			// reconstructible --> reconstructible using detector-reconstructed particles
			// therefore "CMS-rec" --> best possible reconstructible by algorithm: the ideal algorithm would always choose this as the VB-candidate, as its constituent jet are the reconstructed counterparts to the genJets of the genVB
			
			
			// How often do we match the gen VB with rec AK4, without knowing MC truth?
			if(recAK4cand){
				float dRrec00 =physmath::deltaR(recAK4cand->daughter(0), bestPossibleAK4.daughter(0));
				float dRrec01 =physmath::deltaR(recAK4cand->daughter(0), bestPossibleAK4.daughter(1));
				float dRrec10 =physmath::deltaR(recAK4cand->daughter(1), bestPossibleAK4.daughter(0));
				float dRrec11 =physmath::deltaR(recAK4cand->daughter(1), bestPossibleAK4.daughter(1));
				
				// they should be the same jet, so the tolerance can be very low
				if((dRrec00 < 0.01 && dRrec11 < 0.01) || (dRrec01 < 0.01 && dRrec10 < 0.01)){
					// We correctly reconstructed the best possible reconstructible Boson for this particular genBoson with hadronic daughters
					if(!genMatch) genMatch = &genVB; //points to the VB correctly reconstructed by the algorithm.
					theHistograms->fill("Purity AK4 pairs (good): genVB mass", "Purity AK4 pairs (good): gen VB mass", 30,60.,120., genVB.mass(), theWeight);
					theHistograms->fill("Purity AK4 pairs (good): rec mass", "Purity AK4 pairs (good): rec mass", 30,60.,120., recAK4cand->mass(), theWeight);
					theHistograms->fill("Purity AK4 pairs (good): rec pt", "Purity AK4 pairs (good): rec pt", 25,0.,500., recAK4cand->pt(), theWeight);
					theHistograms->fill("Purity AK4 pairs (good): # of rec AK4", "Purity AK4 pairs (good): # of rec AK4", 11,-0.5,10.5, jets->size(), theWeight);
					theHistograms->fill("Purity AK4 pairs: #DeltaM (good rec-gen)", "Purity AK4 pairs: #DeltaM (good rec-gen)", 25,-25.,25., recAK4cand->mass() - genVB.mass(), theWeight);
				}
			}
		}
	}
	
	if(genMatch && recAK4cand){
		++Ngood_recVB_60_120_; // "good" --> the algorithm has picked the correct rec jets
		if(phys::WMASS - 20. < recAK4cand->mass() && recAK4cand->mass() < phys::ZMASS + 20.)
			++Ngood_recVB_m20_;
	}
	return make_pair(genMatch, recAK4cand);
	// NOTE: recAK4cand must be deleted outside
}


void VZZAnalyzer::AKrace(pair<const Boson<Particle>*, Boson<Jet>*> match4, pair<const Particle*, Jet*> match8){
	float dM4 = match4.second->mass() - match4.first->mass();
	float dM8 = match8.second->mass() - match8.first->mass(); //maybe use another mass alg.
	float dM8corrPruned = match8.second->corrPrunedMass() - match8.first->mass();
	float genDM = match4.first->mass() - match8.first->mass();
	float recDM = match4.second->mass() - match8.second->corrPrunedMass();
	
	theHistograms->fill("AKrace: mass res. 4", "AKrace: mass res. 4", 25,-25.,25.,dM4,theWeight);
	theHistograms->fill("AKrace: mass res. 8", "AKrace: mass res. 8", 25,-25.,25.,dM8,theWeight);
	theHistograms->fill("AKrace: corrPrunedMass res. 8", "AKrace: corrPrunedMass res. 8", 25,-25.,25., dM8corrPruned, theWeight);
	theHistograms->fill("AKrace: gen(mass(4) -mass(8))", "AKrace: gen(mass(4) -mass(8))", 25,-25.,25., genDM, theWeight);
	theHistograms->fill("AKrace: rec(mass(4) -mass(8))", "AKrace: rec(mass(4) -mass(8))", 25,-25.,25., recDM, theWeight);
	
	theHistograms->fill<TH2F>("AKrace: #DeltaM vs gen_pt 4", "AKrace: #DeltaM vs gen_pt 4", PT_2D_SIZE, 25,-25.,25., match4.first->pt(), dM4,theWeight);
	theHistograms->fill<TH2F>("AKrace: #DeltaM vs gen_pt 8", "AKrace: #DeltaM vs gen_pt 4", PT_2D_SIZE, 25,-25.,25., match8.first->pt(), dM8,theWeight);
	
	if(ZZ){
		//theHistograms->fill("AKrace: ZZ pt 4", "AKrace: ZZ pt 4", PT_2D_SIZE, ZZ->pt(), )
		theHistograms->fill<TH2F>("AKrace: #DeltaM vs ZZ pt 4", "AKrace: #DeltaM vs gen_pt 4", PT_2D_SIZE, 16,-20.,20., ZZ->pt(), dM4,theWeight);
		theHistograms->fill<TH2F>("AKrace: #DeltaM vs gen_pt 8", "AKrace: #DeltaM vs gen_pt 4", PT_2D_SIZE, 16,-20.,20., ZZ->pt(), dM8,theWeight);
	}
}


void VZZAnalyzer::endReconstructionAK(){
	if(!Nhad_genVB_) return;
	cout<<"\t----- AK4 -----\n";
	cout<<"Generated VB->jj: "<<Nhad_genVB_<<"    Reconstructed by detector: "<<Ncms_recVB_<< Form("    Detector eff.: %.1f %%", 100.*Ncms_recVB_/Nhad_genVB_)<<'\n';
	cout<<"(bestV, 60-120)    Total rec. by alg.: "<<Nall_recVB_60_120_<<"    Matching: "<<Ngood_recVB_60_120_<<Form("    Alg. eff.: %.1f %%", 100.*Ngood_recVB_60_120_/Ncms_recVB_) <<Form("    Alg. purity: %.1f %%", 100.*Ngood_recVB_60_120_/Nall_recVB_60_120_)<<'\n';
	cout<<"(bestV, |mV-20|)   Total rec. by alg.: "<<Nall_recVB_m20_<<"    Matching: "<<Ngood_recVB_m20_<<Form("    Alg. eff.: %.1f %%", 100.*Ngood_recVB_m20_/Ncms_recVB_)<<Form("    Alg. purity: %.1f %%", 100.*Ngood_recVB_m20_/Nall_recVB_m20_)<<'\n';
	cout<<"\t----- AK8 -----\n";
	cout<<"Generated AK8 (60<m<120): "<<N_gen8_<<"    Reconstructed by detector: "<<Ncms_rec8_<< Form("    Detector eff.: %.1f %%", 100.*Ncms_rec8_/N_gen8_)<<'\n';
	cout<<"(bestV, 60-120, pTau<.35)    Total rec. by alg.: "<<Nall_rec8_pTau35_<<"    Matching with a gen: "<<Ngood_rec8_pTau35_<<Form("    Alg. purity: %.1f %%", 100.*Ngood_rec8_pTau35_/Nall_rec8_pTau35_)<<'\n';
	
	if(!NnoAK4_) return;
	cout<<"\t----- EXTRA AK8 -----\n";
	cout<<"No AK4 found: "<<NnoAK4_<<"    There's a genVB->jj_gen: "<<N8genVB_<<"    	Found an AK8_rec with DR<0.2: "<<N8recover_<<'\n';
}


void VZZAnalyzer::makeGenZZ(){
	if(!topology.test(0) || genVBParticles->size() < 2){
		new (genZZ_) DiBoson<Particle, Particle>(); //does not allocate memory but constructs an object at genZZ_ -- see https://www.cplusplus.com/reference/new/operator%20new/
		return;
	}
	
	vector<Boson<Particle>> Zll;
	foreach(const Boson<Particle>& genVB, *genVBParticles){
		int id0 = genVB.daughter(0).id();
		int id1 = genVB.daughter(1).id();
		if(genVB.id() == 23 && 10<abs(id0) && abs(id0)<20 && 10<abs(id1) && abs(id1)<20)
			Zll.push_back(genVB);  // Some conditions are redundant
	}
	
	if(Zll.size() < 2)
		new (genZZ_) DiBoson<Particle, Particle>();
	else
		new (genZZ_) DiBoson<Particle, Particle>(Zll.at(0), Zll.at(1));
}


void VZZAnalyzer::fillGenVBtoAK4(){
	using std::get;
	typedef std::tuple<size_t,size_t,double> myTuple;
	
	if(genJets->size() < 2) return;
	
	vector<Particle> cleanedJets;
	foreach(const Particle& jet, *genJets)
		if(jet.pt() > 30) cleanedJets.push_back(jet);
	
	if(cleanedJets.size() < 2) return;
	
	vector< myTuple > genP_indices;  // Mass in range
	vector< myTuple > goodP_indices; // Mass in range AND no daug. in common
	for(size_t i = 0; i < cleanedJets.size(); ++i)
		for(size_t j = i+1; j < cleanedJets.size(); ++j){
			float mass = (cleanedJets.at(i).p4() + cleanedJets.at(j).p4()).M();
			if(phys::WMASS-30. < mass && mass < phys::ZMASS+30.){ //same range of SignalDefinitions
				genP_indices.push_back( make_tuple(i, j, mass) );
				//theHistograms->fill("Pairs_mass", "Pairs_mass", 35,70.,105., mass, theWeight);
 			}
		}
	
	
 	if(genP_indices.size() > 0){ //genPairs.size() > 0){
 		sort(genP_indices.begin(), genP_indices.end(), Mass2Comparator(phys::ZMASS, phys::WMASS));
		//stable_sort(genPairs.begin(), genPairs.end(), ScalarSumPtComparator());
 		
 		foreach(const myTuple& cand_i, genP_indices){
 			bool match = false;
 			foreach(const myTuple& good_i, goodP_indices){
 				if(get<0>(cand_i) == get<0>(good_i) || get<0>(cand_i) == get<1>(good_i) || get<1>(cand_i) == get<1>(good_i) || get<1>(cand_i) == get<0>(good_i)){
 					match = true;
 					//theHistograms->fill("Conflicts: # of genJets", "Conflicts: # of genJets", 8,2.5,10.5, genJets->size(), theWeight);
					//theHistograms->fill("Conflicts: # of genJets norm", "Conflicts: # of genJets norm", 8,2.5,10.5, genJets->size(), theWeight/genJets->size());
 					break;
 				}
 			}
 			if(!match){
				Boson<Particle> cand_VB(genJets->at(get<0>(cand_i)), genJets->at(get<1>(cand_i)));
				AllGenVBjj_->push_back(cand_VB);
 			}
		}
 	}
 	//theHistograms->fill("Tot_genPairs","Tot_genPairs", 12,-0.5,11.5, genP_indices.size(), theWeight);
 	//theHistograms->fill("Tot unique genPairs", "Tot unique genPairs", 8,-0.5,7.5, AllGenVBjj_->size(), theWeight);
 	//theHistograms->fill<TH2F>("Unique genPairs (|m-30|) vs cleanedJets", "Unique genPairs (|m-30|) vs cleanedJets", 8,-0.5,7.5, 8,-0.5,7.5, cleanedJets.size(), AllGenVBjj_->size(), theWeight);
 	
 	//foreach(const Boson<Particle>& uVB, *AllGenVBjj_)
 		//theHistograms->fill("UniqueVBjj_mass", "Unique VB->jj: mass", 35,70.,105., uVB.mass(), theWeight);
}


void VZZAnalyzer::resolutionZmass(){  // same as bestZMassJetMVA() but without dR(jet, ZZ)>2.
	MassComparator ZmassComp = MassComparator(phys::ZMASS);
	/*
	//------------------	GEN AK8	------------------
	if(genJetsAK8->size() > 0){
		stable_sort( genJetsAK8->begin(), genJetsAK8->end(), ZmassComp );
		float mass_8g = genJetsAK8->front().mass();
		if( fabs(mass_8g - phys::ZMASS) < 30.){
			theHistograms->fill<TH2F>("BestZ_AK8_gen", "BestZ AK8_{gen} mass vs ZZ_{gen} pt;pt ZZ [GeV/c];AK8_{gen} mass [GeV/c]", PT_2D_SIZE, MASS_2D_SIZE, genZZ_->pt(), mass_8g, theWeight);
		}
	}
	//------------------	GEN AK4	------------------
	if(genHadVBs_->size() > 0){
		stable_sort( genHadVBs_->begin(), genHadVBs_->end(), ZmassComp );
		float mass_4g = genHadVBs_->front().mass();
		if( fabs(mass_4g - phys::ZMASS) < 30.){
			theHistograms->fill<TH2F>("BestZ_AK4s_gen", "BestZ AK4s_{gen} mass vs ZZ_{gen} pt;pt ZZ [GeV/c];AK4s_{gen} mass [GeV/c]", PT_2D_SIZE, MASS_2D_SIZE, genZZ_->pt(), mass_4g, theWeight);
		}
	}
	*/
	
	// ------------------	REC AK8	------------------
	/*if(jetsAK8->size() > 0){
		stable_sort( jetsAK8->begin(), jetsAK8->end(), ZmassComp );
		float mass_8r = jetsAK8->front().mass();
		if( fabs(mass_8r - phys::ZMASS) < 30.){
			theHistograms->fill<TH2F>("BestZ_AK8_rec", "BestZ AK8_{rec} mass vs ZZ_{rec} pt;pt ZZ [GeV/c];AK8_{gen} mass [GeV/c]", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), mass_8r, theWeight);
		}
	}*/
	//------------------	REC AK4	------------------
	if(AK4pairs_->size() > 0){
		stable_sort( AK4pairs_->begin(), AK4pairs_->end(), ZmassComp );
		float mass_4r = AK4pairs_->front().mass();
		theHistograms->fill("BestZ r4 mass", "BestZ AK4s rec mass;mass [GeV/c^{2}];# evts", 30,60.,120., mass_4r, theWeight);
		theHistograms->fill("BestZ r4 dM _r_", "BestZ AK4s rec dM;mass [GeV/c^{2}];# evts", 32,0.,32., fabs(mass_4r - phys::ZMASS), theWeight);
		//if( fabs(mass_4r - phys::ZMASS) < 30.){
			//theHistograms->fill<TH2F>("BestZ_AK4s_rec","BestZ AK4s_{rec} mass vs ZZ_{rec} pt;pt ZZ [GeV/c];AK4s_{rec} mass [GeV/c]", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), mass_4r, theWeight);
		//}
	}
}


void VZZAnalyzer::specialPeakAnalisys(const Particle& theGenAK8){
	stable_sort(AK4pairs_->begin(), AK4pairs_->end(), DeltaRComparator(theGenAK8) );
	
	if(physmath::deltaR(AK4pairs_->front(), theGenAK8) < 0.4){
		theHistograms->fill("Special Peak: AK4 #DeltaR", "Special Peak: AK4 #DeltaR", 10,0.,2., physmath::deltaR( AK4pairs_->front().daughter(0), AK4pairs_->front().daughter(1) ) , 1.);
		theHistograms->fill("Special Peak: #DeltaR(AK8, Boson)", "Special Peak: #DeltaR(AK8, Boson)", 8,0.,0., physmath::deltaR( AK4pairs_->front(), theGenAK8 ) , 1.);
		theHistograms->fill("Special Peak: #DeltaM(AK8, Boson)", "Special Peak: #DeltaM(AK8, Boson)", 10,-20.,20., AK4pairs_->front().mass() - theGenAK8.mass(), 1.);
		theHistograms->fill("Special Peak: #DeltaP(AK8, Boson)", "Special Peak: #DeltaP(AK8, Boson)", 10,0.,50., (AK4pairs_->front().p4() - theGenAK8.p4()).P() , 1.);
	}
	else //Counting events without AK4 pairs
		theHistograms->fill("Special Peak: AK4 #DeltaR", "Special Peak: AK4 #DeltaR", 10,0.,2., -100., 1.);
		
	return;
}


void VZZAnalyzer::bestZMassJetMVA(){
	if(!ZZ || ZZ->p() < 1.) return;
	
	//-----------------------------	GEN PARTICLES	-----------------------------
	//Mass2Comparator VBmassComp(phys::ZMASS, phys::WMASS);
	MassComparator ZmassComp = MassComparator(phys::ZMASS);
	//------------------	GEN AK8	------------------
	auto it_8g = genJetsAK8->begin();
	bool found_8g = false;
	if(genJetsAK8->size() > 0){
		stable_sort( genJetsAK8->begin(), genJetsAK8->end(), ZmassComp );
		while(it_8g != genJetsAK8->end()){
			if(physmath::deltaR(*it_8g, *ZZ) > 2. && ZBosonDefinition(*it_8g)){
				found_8g = true;
				break;
			}
			else ++it_8g;
		}
		if(found_8g){
			theHistograms->fill("bestZ AK8_{gen} mass", "bestZ AK8_{gen} mass", 30,0.,150, it_8g->mass(), theWeight);
			theHistograms->fill<TH2F>("BestZ AK8_{gen} mass vs ZZ pt", "BestZ AK8_{gen} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), it_8g->mass(), theWeight);
			//theHistograms->fill<TH2F>("BestZ AK8_{gen} mass vs ZZ mass", "BestZ AK8_{gen} mass vs ZZ mass", ZZMASS_2D_SIZE, MASS_2D_SIZE, ZZ->mass(), it_8g->mass(), theWeight);
			//theHistograms->fill<TH2F>("BestZ AK8_{gen} mass vs s", "BestZ AK8_{gen} mass vs s", S_2D_SIZE, MASS_2D_SIZE, sAK8g_, it_8g->mass(), theWeight);
			
			//Special region
			if( 80. < ZZ->pt() && ZZ->pt() < 120. && 90. < it_8g->mass() && it_8g->mass() < 110.)
				specialPeakAnalisys(*it_8g);
		}
	}
	
	//------------------	GEN AK4	------------------
	//fillGenHadVBs();
	auto it_4g = genHadVBs_->begin();
	bool found_4g = false;
	if(genHadVBs_->size() > 0){
		stable_sort( genHadVBs_->begin(), genHadVBs_->end(), ZmassComp );
		while(it_4g != genHadVBs_->end()){
			if(physmath::deltaR(*it_4g, *ZZ) > 2. && ZBosonDefinition(*it_4g)){
				found_4g = true;
				break;
			}
			else ++it_4g;
		}
		if(found_4g){
			theHistograms->fill("bestZ AK4_{gen} mass", "bestZ AK4_{gen} mass", 30,0.,150, it_4g->mass(), theWeight);
			theHistograms->fill<TH2F>("BestZ AK4_{gen} mass vs ZZ pt", "BestZ AK4_{gen} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), it_4g->mass(), theWeight);
			//theHistograms->fill<TH2F>("BestZ AK4_{gen} mass vs ZZ mass", "BestZ AK4_{gen} mass vs ZZ mass", ZZMASS_2D_SIZE, MASS_2D_SIZE, ZZ->mass(), it_4g->mass(), theWeight);
			//theHistograms->fill<TH2F>("BestZ AK4_{gen} mass vs s", "BestZ AK4_{gen} mass vs s", S_2D_SIZE, MASS_2D_SIZE, sAK4g_, it_4g->mass(), theWeight);
		}
	}
	
	return;  //TEMP
	
	//------------------	WINNER GEN	------------------
	if(found_4g || found_8g){
		float type = 0.;  // -1 --> AK4 pair,   1 --> AK8
		if(found_4g && found_8g){
			type = (fabs(it_8g->mass() -phys:: ZMASS) < fabs(it_4g->mass() -phys:: ZMASS) ? 1.:-1.);
			theHistograms->fill<TH2F>("Winner Z mass", "Winner Z mass", MASS_2D_SIZE, BINARY_SIZE, (type > 0.5 ? *it_8g : *it_4g).mass(), type, theWeight);
			theHistograms->fill<TH2F>("Loser Z mass", "Loser Z mass", MASS_2D_SIZE, BINARY_SIZE, (type < 0.5 ? *it_8g : *it_4g).mass(), type, theWeight);
		}
		else if(found_8g)
			type = 1.;
		else if(found_4g)
			type = -1.;
		
		double sHat_g = ( ZZ->p4() + (type > 0.5 ? *it_8g : *it_4g).p4() ).M();
		theHistograms->fill<TH2F>("Best Z gen vs #hat{s}", "Best Z gen vs #hat{s}", S_2D_SIZE, BINARY_SIZE, sHat_g, type, theWeight);
		theHistograms->fill<TH2F>("Best Z gen vs s", "Best Z gen vs s", S_2D_SIZE, BINARY_SIZE, (type > 0.5 ? sAK8g_ : sAK4g_), type, theWeight);
		theHistograms->fill<TH2F>("Best Z gen vs ZZ pt", "Best Z gen vs ZZ pt", PT_2D_SIZE, BINARY_SIZE, ZZ->pt(), type, theWeight);
		theHistograms->fill<TH2F>("Best Z gen vs #Delta#phi(ZZ, Jet)_div_pi", "Best Z gen vs #Delta#phi(ZZ, Jet)_div_#pi", 40,-1.,1., BINARY_SIZE, physmath::deltaPhi( *ZZ, (type > 0.5 ? *it_8g : *it_4g) )/M_PI, type, theWeight); //32,-3.2,3.2
		theHistograms->fill<TH2F>("Best Z gen vs #Delta#eta(ZZ, Jet)", "Best Z gen vs #Delta#eta(ZZ, Jet)", 20,-4.,4., BINARY_SIZE, ZZ->eta()-(type > 0.5 ? *it_8g : *it_4g).eta(), type, theWeight);
		theHistograms->fill<TH2F>("Best Z gen vs #DeltaR(ZZ, Jet)", "Best Z gen vs #DeltaR(ZZ, Jet)", 40,2.,6., BINARY_SIZE, physmath::deltaR( *ZZ, (type > 0.5 ? *it_8g : *it_4g) ), type, theWeight);
		
		theHistograms->fill<TH2F>("Best Z vs deltaR(4_{gen})", "Best Z vs deltaR(4_{gen})", 50,0.,5., BINARY_SIZE, (found_4g ? physmath::deltaR( it_4g->daughter(0), it_4g->daughter(1) ) : -1.), type /*filling underflow bin --> bestZ found with AK8*/, theWeight);
	theHistograms->fill<TH2F>("Best Z vs deltaPhi(4_{gen})","Best Z vs deltaPhi(4_{gen})", 32,-3.2,3.2, BINARY_SIZE, (found_4g ? physmath::deltaPhi(it_4g->daughter(0), it_4g->daughter(1)) : -9.), type, theWeight);
	}
	
	if(found_8g && !found_4g){	//Special_174 AK8 that win because they're alone
		theHistograms->fill("Special_174 mass", "Special_174 mass", 16,50.,130., it_8g->mass(), theWeight);
		theHistograms->fill("Special_174 pt",   "Special_174 pt",  16,200.,1000., it_8g->pt(), theWeight);
		theHistograms->fill("Special_174 #Delta#phi(ZZ)_div_pi", "Special_174 #Delta#phi(ZZ)_div_pi",  20,-1.,1., physmath::deltaPhi(*ZZ, *it_8g)/M_PI, theWeight);
	}
	
	
	
	//-----------------------------	REC PARTICLES	-----------------------------
	//------------------	REC AK8	------------------
	auto it_8r = jetsAK8->begin();
	bool found_8r = false;
	if(jetsAK8->size() > 0){
		stable_sort( jetsAK8->begin(), jetsAK8->end(), MassComparator(phys::ZMASS) );
		while(it_8r != jetsAK8->end()){
			if(physmath::deltaR(*it_8r, *ZZ) > 2. && ZBosonDefinition(*it_8r) /*&& tauCut(it_8r)*/){
				found_8r = true;
				break;
			}
			else ++it_8r;
		}
		if(found_8r){
			theHistograms->fill("bestZ AK8_{rec} mass", "bestZ AK8_{rec} mass", 30,0.,150, it_8r->mass(), theWeight);
			theHistograms->fill<TH2F>("BestZ AK8_{rec} mass vs ZZ pt", "BestZ AK8_{rec} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), it_8r->corrPrunedMass(), theWeight);
			theHistograms->fill<TH2F>("BestZ AK8_{rec} mass vs ZZ mass", "BestZ AK8_{rec} mass vs ZZ mass", PT_2D_SIZE, MASS_2D_SIZE, ZZ->mass(), it_8r->corrPrunedMass(), theWeight);
			
			theHistograms->fill("Best Z (8_{rec}) tau21", "Best Z (8_{rec}) tau21", 20,0.,1., it_8r->tau2(), theWeight);
		}
	}
	
	//------------------	REC AK4	------------------
	auto it_4r = AK4pairs_->begin();
	bool found_4r = false;
	if(AK4pairs_->size() > 0){
		stable_sort( AK4pairs_->begin(), AK4pairs_->end(), MassComparator(phys::ZMASS) );
		while(it_4r != AK4pairs_->end()){
			if(physmath::deltaR(*it_4r, *ZZ) > 2. && ZBosonDefinition(*it_4r)){
				found_4r = true;
				break;
			}
			else ++it_4r;
		}
		if(found_4r){
			theHistograms->fill("bestZ AK4_{rec} mass", "bestZ AK4_{rec} mass", 30,0.,150, it_4r->mass(), theWeight);
			theHistograms->fill<TH2F>("BestZ AK4_{rec} mass vs ZZ pt", "BestZ AK4_{rec} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), it_4r->mass(), theWeight);
			theHistograms->fill<TH2F>("BestZ AK4_{rec} mass vs ZZ mass", "BestZ AK4_{rec} mass vs ZZ mass", PT_2D_SIZE, MASS_2D_SIZE, ZZ->mass(), it_4r->mass(), theWeight);
		}
	}
	
	//------------------	WINNER REC	------------------
	if(found_4r || found_8r){
		float type = 0.;  // -1 --> AK4 pair,   1 --> AK8
		if(found_4r && found_8r)
			type = (fabs(it_8r->corrPrunedMass() - phys:: ZMASS) < fabs(it_4r->mass() - phys:: ZMASS) ?1.:-1.);
		else if(found_8r)
			type = 1.;
		else if(found_4r)
			type = -1.;
		
		theHistograms->fill<TH2F>("Best Z rec vs ZZ pt", "Best Z rec vs ZZ pt", PT_2D_SIZE, BINARY_SIZE, ZZ->pt(), type, theWeight);
		theHistograms->fill<TH2F>("Best Z rec vs s", "Best Z rec vs s", S_2D_SIZE, BINARY_SIZE, (type > 0.5 ? sAK8r_ : sAK4r_), type, theWeight);
		
		theHistograms->fill<TH2F>("Best Z vs deltaR(4_{rec})", "Best Z vs deltaR(4_{rec})", 50,0.,5., BINARY_SIZE, (found_4r ? physmath::deltaR( it_4r->daughter(0), it_4r->daughter(1) ) : -1.), type /*filling underflow bin --> bestZ found with AK8*/, theWeight);
		theHistograms->fill<TH2F>("Best Z vs deltaPhi(4_{rec})","Best Z vs deltaPhi(4_{rec})", 32,-3.2,3.2, BINARY_SIZE, (found_4r ? physmath::deltaPhi(it_4r->daughter(0), it_4r->daughter(1)) : -9.), type, theWeight);
	
	}
}


void VZZAnalyzer::minPtJetMVA(){
	if(!ZZ || ZZ->p() < 1.) return;
	
	//------------------	GEN AK8	------------------
	if(genJetsAK8->size() > 0){
		stable_sort( genJetsAK8->begin(), genJetsAK8->end(), PtTotRefComparator(*ZZ) );
		auto it = genJetsAK8->begin();
		bool found = false;
		while(it != genJetsAK8->end()){
			if(physmath::deltaR(*it, *ZZ) > 2. && it->mass() > 50.){
				found = true;
				break;
			}
			else ++it;
		}
		if(found){
			theHistograms->fill("#DeltaR(ZZ, minPt AK8_{gen})", "#DeltaR(ZZ, minPt AK8_{gen})", 35,0.,7., physmath::deltaR(*it, *ZZ), 1.);
			theHistograms->fill<TH2F>("Min PtTot AK8_{gen} mass vs ZZ pt", "Min PtTot AK8_{gen} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), it->mass(), 1.);
			
			//Special region
			/*
			if( 80. < ZZ->pt() && ZZ->pt() < 120. && 90. < it->mass() && it->mass() < 110.)
				specialPeakAnalisys(*it);*/
		}
	}
	//------------------	GEN AK4	------------------
	//fillGenHadVBs();
	if(genHadVBs_->size() > 0){
		stable_sort( genHadVBs_->begin(), genHadVBs_->end(), PtTotRefComparator(*ZZ) );
		auto it = genHadVBs_->begin();
		bool found = false;
		while(it != genHadVBs_->end()){
			if(physmath::deltaR(*it, *ZZ) > 2. && it->mass() > 50.){
				found = true;
				break;
			}
			else ++it;
		}
		if(found){
			theHistograms->fill("#DeltaR(ZZ, minPt AK4_{gen})", "#DeltaR(ZZ, minPt AK4_{gen})", 35,0.,7., physmath::deltaR(*it, *ZZ), 1.);
			theHistograms->fill<TH2F>("Min PtTot AK4_{gen} mass vs ZZ pt", "Min PtTot AK4_{gen} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), it->mass(), 1.);
		}
	}
	
	
	//------------------	REC AK8	------------------
	if(jetsAK8->size() > 0){
		stable_sort( jetsAK8->begin(), jetsAK8->end(), PtTotRefComparator(*ZZ) );
		auto it = jetsAK8->begin();
		bool found = false;
		while(it != jetsAK8->end()){
			if(physmath::deltaR(*it, *ZZ) > 2. && it->corrPrunedMass() > 50.){
				found = true;
				break;
			}
			else ++it;
		}
		if(found){
			theHistograms->fill("#DeltaR(ZZ, minPt AK8_{rec})", "#DeltaR(ZZ, minPt AK8_{rec})", 35,0.,7., physmath::deltaR(*it, *ZZ), 1.);
			theHistograms->fill<TH2F>("Min PtTot AK8_{rec} mass vs ZZ pt", "Min PtTot AK8_{rec} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), it->corrPrunedMass(), 1.);
		}
	}
	//------------------	REC AK4	------------------
	if(AK4pairs_->size() > 0){
		stable_sort( AK4pairs_->begin(), AK4pairs_->end(), PtTotRefComparator(*ZZ) );
		auto it = AK4pairs_->begin();
		bool found = false;
		while(it != AK4pairs_->end()){
			if(physmath::deltaR(*it, *ZZ) > 2. && it->mass() > 50.){
				found = true;
				break;
			}
			else ++it;
		}
		if(found){
			theHistograms->fill("#DeltaR(ZZ, minPt AK4_{rec})", "#DeltaR(ZZ, minPt AK4_{rec})", 35,0.,7., physmath::deltaR(*it, *ZZ), 1.);
			theHistograms->fill<TH2F>("Min PtTot AK4_{rec} mass vs ZZ pt", "Min PtTot AK4_{rec} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), it->mass(), 1.);
		}
	}
}


void VZZAnalyzer::furthestJetMVA(){
	if(!ZZ || ZZ->p() < 1.) return;
	
	//------------------	GEN PARTICLES	------------------
	//For the AK4 part: genVB are constructed from pairs of AK4 (FOR NOW!) See SignalDefinitions
	//fillGenHadVBs();
	
	if(genHadVBs_->size() > 0){ //Looks like sometimes ZZ is a void DiBoson
		if(genHadVBs_->size() > 1)
			std::stable_sort(genHadVBs_->begin(), genHadVBs_->end(), DeltaRComparator(*ZZ));
		
		//Now the furthest from ZZ is the LAST of genHadVBs
		theHistograms->fill<TH2F>("Furthest AK4_{gen} mass vs ZZ pt", "Furthest AK4_{gen} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), genHadVBs_->back().mass(), 1.);
		theHistograms->fill<TH2F>("Furthest AK4_{gen} mass vs ZZ P",  "Furthest AK4_{gen} mass vs ZZ P",  P_2D_SIZE,  MASS_2D_SIZE, ZZ->p(),  genHadVBs_->back().mass(), 1.);
		theHistograms->fill("Max #DeltaR(ZZ, AK4_{gen})", "Max #DeltaR(ZZ, AK4_{gen})", 30,2.,8., physmath::deltaR(genJets->back(), *ZZ), 1.);
	}
	
	// AK8
	Particle* furthGenAK8 = furthestSing(genJetsAK8, *ZZ, 2., make_pair(50.,130.));
	if(furthGenAK8 != nullptr){
		theHistograms->fill<TH2F>("Furthest AK8_{gen} mass vs ZZ pt", "Furthest AK8_{gen} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), furthGenAK8->mass(), 1.);
		theHistograms->fill<TH2F>("Furthest AK8_{gen} mass vs ZZ P",  "Furthest AK8_{gen} mass vs ZZ P",  P_2D_SIZE,  MASS_2D_SIZE, ZZ->p(),  furthGenAK8->mass(), 1.);
		theHistograms->fill("Max #DeltaR(ZZ, AK8_{gen})", "Max #DeltaR(ZZ, AK8_{gen})", 30,2.,8., physmath::deltaR(*furthGenAK8, *ZZ), 1.);
		delete furthGenAK8;
	}
	
	//------------------	REC PARTICLES	------------------
	// AK8
	Jet* furthestAK8 = furthestSing(jetsAK8, *ZZ, 2., make_pair(50.,130.));
	if(furthestAK8 != nullptr){
		if(furthestAK8->corrPrunedMass() < 50.)
			cout<<"  Outside "<<furthestAK8->corrPrunedMass()<<'\n';
		theHistograms->fill<TH2F>("Furthest AK8_{rec} mass vs ZZ pt", "Furthest AK8_{rec} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), furthestAK8->corrPrunedMass(), 1.);
		theHistograms->fill<TH2F>("Furthest AK8_{rec} mass vs ZZ P",  "Furthest AK8_{rec} mass vs ZZ P",  P_2D_SIZE,  MASS_2D_SIZE, ZZ->p(),  furthestAK8->corrPrunedMass(), 1.);
		theHistograms->fill("Max #DeltaR(ZZ, AK8_{rec})", "Max #DeltaR(ZZ, AK8_{rec})", 30,2.,8., physmath::deltaR(*furthestAK8, *ZZ), 1.);
		delete furthestAK8;
	}
	
	// AK4
	Boson<Jet>* furthestAK4 = furthestPair(jets, *ZZ, 2./*min deltaR*/, make_pair(5.,200.)); 
	if(furthestAK4 != nullptr){
		theHistograms->fill<TH2F>("Furthest AK4_{rec} mass vs ZZ pt", "Furthest AK4_{rec} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), furthestAK4->mass(), 1.);
		theHistograms->fill<TH2F>("Furthest AK4_{rec} mass vs ZZ P",  "Furthest AK4_{rec} mass vs ZZ P", P_2D_SIZE,  MASS_2D_SIZE, ZZ->p(),  furthestAK4->mass(), 1.);
		theHistograms->fill("Max #DeltaR(ZZ, AK4_{rec})", "Max #DeltaR(ZZ, AK8_{rec})", 30,2.,8., physmath::deltaR(*furthestAK4, *ZZ), 1.);
		delete furthestAK4;
	}
}


void VZZAnalyzer::closestJetAnalisys(){
	//1: find the VBs with hadronic deacay:
	//fillGenHadVBs();
	if(genHadVBs_->size() == 0) return; //No analysis can be done without a V-->JJ
	
	
	foreach(const Boson<Particle>& hadVB, *genHadVBs_){
		//2.a: find the closest genAK8
		const Particle* closestGenAK8 = closestSing(genJetsAK8, hadVB);
		if(closestGenAK8){
			float dR = physmath::deltaR(*closestGenAK8, hadVB);
			float dM = closestGenAK8->mass() - hadVB.mass(); //genParticles --> mass()
			theHistograms->fill("#DeltaM (AK8_{gen}, V)", "#DeltaM (AK8_{gen}, V)", 40,-20,20, dM,1.);
			theHistograms->fill("#DeltaR (AK8_{gen}, V)", "#DeltaR (AK8_{gen}, V)", 40,0.,0.4, dR,1.);
		}
	
		//2.b: find the closest pair of genAK4
		Boson<Particle>* closestGenAK4s = closestPair(genJets, hadVB);
		if(closestGenAK4s){
			float dR = physmath::deltaR(*closestGenAK4s, hadVB);
			float dM = closestGenAK4s->mass() - hadVB.mass(); //genParticles --> mass()
			theHistograms->fill("#DeltaM (AK4_{gen} pair, V)", "#DeltaM (AK4_{gen} pair, V)", 40,-20,20, dM, 1.);
			theHistograms->fill("#DeltaR (AK4_{gen} pair, V)", "#DeltaR (AK4_{gen} pair, V)", 40,0.,0.4, dR, 1.);
		}
		
		//3.a: find the closest recAK8
		const Jet* closestAK8 = closestSing(jetsAK8, hadVB);
		if(closestAK8){
			float dR = physmath::deltaR(*closestAK8, hadVB);
			float dM = closestAK8->corrPrunedMass() - hadVB.mass(); //Jet -->corrPrunedmass()
			theHistograms->fill("#DeltaM_{corr-prun} (AK8_{rec}, V)", "#DeltaM_{corr-prun} (AK8_{rec}, V)", 40,-20,20, dM, 1.);
			theHistograms->fill("#DeltaR (AK8_{rec}, V)", "#DeltaR (AK8_{rec}, V)", 40,0.,0.4, dR,1.);
		}
		
		//3.b: find the closest pair of recAK4
		Boson<Jet>* closestAK4s = closestPair(jets, hadVB);
		if(closestAK4s != nullptr){
			float dR = physmath::deltaR(*closestAK4s, hadVB);
			float dM = closestAK4s->mass() - hadVB.mass(); //Boson is Particle --> mass()
			theHistograms->fill("#DeltaM (AK4_{rec} pair, V)", "#DeltaM (AK4_{rec} pair, V)", 40,-20,20, dM, 1.);
			theHistograms->fill("#DeltaR (AK4_{rec} pair, V)", "#DeltaR (AK4_{rec} pair, V)", 40,0.,0.4, dR, 1.);
		}
		
		//4: How often AK8 are better? Are there some variables that discriminate?
		if(closestAK4s != nullptr && closestAK8 != nullptr){
			float dM8 = closestAK8->corrPrunedMass() - hadVB.mass();
			float dM4 = closestAK4s->mass() - hadVB.mass();
			if(fabs(dM8) < fabs(dM4)){
				++win8_;
				theHistograms->fill("V_pt AK8 wins",  "V_pt AK8 wins",  21,0.,630,  hadVB.pt(),  1.);
				theHistograms->fill("V_P AK8 wins",   "V_P AK8 wins",   20,0.,2000, hadVB.p(),   1.);
				theHistograms->fill("V_eta AK8 wins", "V_eta AK8 wins",21,-3.15,3.15,hadVB.eta(),1.);
			}
			else{
				++win4_;
				theHistograms->fill("V_pt AK4s win",  "V_pt AK4s win",  21,0.,630,  hadVB.pt(),  1.);
				theHistograms->fill("V_P AK4s win",   "V_P AK4s win",   20,0.,2000, hadVB.p(),   1.);
				theHistograms->fill("V_eta AK4s win", "V_eta AK4s win",21,-3.15,3.15,hadVB.eta(),1.);
			}
		}
		
		//Cleanup of allocated particles
		if(closestGenAK4s) delete closestGenAK4s;
		if(closestAK4s)    delete closestAK4s;
	}
}


template <class P, class R = Boson<Particle>> // P = Jet or Particle
const P* VZZAnalyzer::closestSing(vector<P>* cands, const R& reference, size_t& k){
	if(cands->size() == 0) return nullptr;
	if(cands->size() == 1) { k = 0; return &(cands->front());}
	
	auto it_best = std::min_element(cands->begin(), cands->end(), Mass2Comparator(phys::WMASS, phys::ZMASS));
	k = it_best - cands->begin();
	
	if(k < 99 && physmath::deltaR(reference, cands->at(k)) < 0.4) //cands->front()) < 0.4 )
		return &(*it_best);
	else return nullptr;
}


template <class P, class R = Boson<Particle>> // P = Jet or Particle
Boson<P>* VZZAnalyzer::closestPair(vector<P>* cands, const R& reference){
	if(cands->size() < 2) return nullptr;
	if(cands->size() == 2){
		Boson<P>* res = new Boson<P>( cands->at(0), cands->at(1) );
		if(physmath::deltaR( *res, reference ) < 0.4 )
			return res;
		else{
			delete res;
			return nullptr;
		}
	}
	
	//Find the pair with the closest p4
	pair<size_t, size_t> indices(0,0);
	float minDR = 0.4; //starting value = the threshold we use
	for(size_t i = 0; i < cands->size(); ++i)
		for(size_t j = i+1; j < cands->size(); ++j){
			TLorentzVector p4Cand = cands->at(i).p4() + cands->at(j).p4();
			float DR = physmath::deltaR( p4Cand, reference.p4() );
			if(DR < minDR){
				minDR = DR;
				indices = std::make_pair(i,j);
			}
		}
	if(indices.second != 0) //then we've found a good pair
		return new Boson<P>( cands->at(indices.first), cands->at(indices.second) );
	else return nullptr;
}

template <class P, class R> //P = Particle or Jet
P* VZZAnalyzer::furthestSing(vector<P>* cands, const R& reference, const float& minDR, const pair<float,float>& massLimits){
	if(cands->size() < 1)
		return nullptr;
	
	vector<P> goodCands;
	for(size_t i = 0; i < cands->size(); ++i){
		float thisMass = getRefinedMass(cands->at(i));
		if(thisMass > massLimits.first && thisMass < massLimits.second){
			goodCands.push_back(cands->at(i));
		}
	}
	if(goodCands.size() > 0){	
		if(goodCands.size() > 1)
			std::sort( goodCands.begin(), goodCands.end(), phys::DeltaRComparator(reference) );
		if(physmath::deltaR(reference, goodCands.back()) > minDR ){
			P* result = new P(goodCands.back());
			//double resMass = getRefinedMass(*result);
			return result;
		}
		else return nullptr;
	}
	else{
		return nullptr;
	}
}


template <class P, class R> // P = Jet or Particle
Boson<P>* VZZAnalyzer::furthestPair(vector<P>* cands, const R& reference, const float& minDR, const pair<float,float>& massLimits){
	if(cands->size() < 2)
		return nullptr;
	
	//Find the pair with the furthest p4 that has minMass < mass < maxMass
	vector<Boson<P>> pairs; //goodPairs
	for(size_t i = 0; i < cands->size(); ++i)
		for(size_t j = i+1; j < cands->size(); ++j){
			Boson<P> thisPair(cands->at(i), cands->at(j));
			if(thisPair.mass() > massLimits.first && thisPair.mass() < massLimits.second){
				pairs.push_back(thisPair);
			}
		}
	if(pairs.size() > 0){
		if(pairs.size() > 1)
			std::stable_sort(pairs.begin(), pairs.end(), phys::DeltaRComparator(reference));
		return new Boson<P>( pairs.back() ); //Copy local object
	}
	else return nullptr;
	//pairs.clear();  // Delete other local Boson objects
}


void VZZAnalyzer::endClosestJetAn(){
	cout<<"Reconstruction of a VB --> AK8: "<<win8_<<" \tAK4 pairs: "<<win4_;
	const char* hNames[] = {"V_pt AK8 wins", "V_P AK8 wins", "V_eta AK8 wins", "V_pt AK4s win", "V_P AK4s win", "V_eta AK4s win"};
	foreach(const char* name, hNames){
		TH1* h = theHistograms->get(name);
		Double_t scale = 1./(h->Integral());
		h->Scale(scale);
	}
}


void VZZAnalyzer::ptCutMVA(){
	if(genHadVBs_->size() == 0) return;
	vector<float> vector_cuts;
	for(int i = 0; i < 10; ++i){
		vector_cuts.push_back((float)(30.+30*i)); //there's no variation before 170 GeV
	}
	
	// ALL GENERATED
	foreach(const Particle& jet, *genJetsAK8){
		float genMass = jet.mass();
		foreach(const float& cut, vector_cuts){
			if(genJetsAK8->front().pt() < cut)
				break;
			const char* name = Form("AK8 gen Mass (all): pt > %d", (int)(cut));
			theHistograms->fill(name, Form("%s;mass [GeV/c^{2}]", name), 30,0.,120., genMass, theWeight);
		}
	}
		// genZZ_ pt
	if(genZZ_ && genZZ_->p() > 1.){
		foreach(const float& cut, vector_cuts){
			if(genZZ_->pt() < cut)
				break;
			const char* name = Form("AK8 gen Mass (all): genZZ pt > %d", (int)(cut));
			const char* title = Form("%s;mass [GeV/c^{2}]", name);
			foreach(const Particle& jet, *genJetsAK8)
				theHistograms->fill(name, title, 30,0.,120., getRefinedMass(jet), theWeight);
		}
	}
	
	
	Boson<Particle>& genVB = genHadVBs_->front();
	// GENERATED CLOSE
	stable_sort(genJetsAK8->begin(), genJetsAK8->end(), phys::DeltaRComparator(genVB));
	if(physmath::deltaR(genJetsAK8->front(), genVB) < 0.2){
		//It matches, should have the mass of a VB
		float genMass = genJetsAK8->front().mass();
		foreach(const float& cut, vector_cuts){
			if(genJetsAK8->front().pt() < cut)
				break;
			const char* name = Form("AK8 gen Mass: pt > %d", (int)(cut));
			theHistograms->fill(name, Form("%s;mass [GeV/c^{2}]", name), 30,0.,120., genMass, theWeight);
		}
		/*
		// RECONSTRUCTED
		//if(jetsAK8->size() == 0)
			//return; //Give up: can't analyze mass functions for AK8s if there are none
		stable_sort(jetsAK8->begin(), jetsAK8->end(), phys::DeltaRComparator(genVB));
		//DeltaRComparator is defined in Commons/interface/Utils.h
		if(physmath::deltaR(jetsAK8->front(), genVB) < 0.5){
			//We have a match! This AK8 should have the mass of a VB (80-90 GeV)
			//Compute the mass once and for all
			float mass           = jetsAK8->front().mass();
			float secvtxMass     = jetsAK8->front().secvtxMass();
			float corrPrunedMass = jetsAK8->front().corrPrunedMass();
			float prunedMass     = jetsAK8->front().prunedMass();
			float softDropMass   = jetsAK8->front().softDropMass();
			float puppiMass      = jetsAK8->front().puppiMass();
			float massesVal[6] = {mass, secvtxMass, corrPrunedMass, prunedMass, softDropMass, puppiMass}; //[1] = {corrPrunedMass};
			
			foreach(const float& cut, vector_cuts){
				if(jetsAK8->front().pt() < cut)
					break; //If it doesn't pass this pt selection, it won't pass the others
				for(int i = 0; i < 6 ; ++i){  // massesVal's size
					const char* name = Form("%s: pt > %d", massAlgsNames_[i], (int)(cut));
					theHistograms->fill(name, Form("%s;mass [GeV/c^{2}]", name), 30,0.,120., massesVal[i], theWeight);
				}
			}
		}*/
	}
}


void VZZAnalyzer::endResolutionAnalisys(TFile& f_out){
	TH2* genAK8_ZZpt = dynamic_cast<TH2*>(theHistograms->get("BestZ AK8_{gen} mass vs ZZ pt"));
	if(genAK8_ZZpt == nullptr) return;
	TH2* genAK4_ZZpt = dynamic_cast<TH2*>(theHistograms->get("BestZ AK4_{gen} mass vs ZZ pt"));
	
	f_out.cd();
	// Resolution of GEN 8
	genAK8_ZZpt->RebinX(2);
	vector<double> std_8g;
	vector<double> std_err_8g;
	vector<double> bins_8g;
	vector<double> bins_w_8g;
	for(int x = 1; x <= genAK8_ZZpt->GetNbinsX(); ++x){
		TH1D* px = genAK8_ZZpt->ProjectionY(Form("px_%d",x), x,x+1, "e");
		double bcx = ((TAxis*)genAK8_ZZpt->GetXaxis())->GetBinCenter(x);
		double bwx = ((TAxis*)genAK8_ZZpt->GetXaxis())->GetBinWidth(x);
		bins_8g.push_back(bcx);
		bins_w_8g.push_back(bwx/sqrt(12));
		std_8g.push_back(px->GetStdDev());
		std_err_8g.push_back(px->GetStdDevError());
	}
	TGraphErrors stdDev8g(bins_8g.size(), bins_8g.data(), std_8g.data(),  bins_w_8g.data(), std_err_8g.data());
	stdDev8g.SetName("Resolution AK8 with error");
	stdDev8g.SetTitle("Resolution AK8 vs ZZ pt;ZZ pt [GeV/c];Resolution [GeV/c^{2}]");
	stdDev8g.SetMinimum(0.);
	stdDev8g.Write();
	
	genAK8_ZZpt->ProfileX("Precision AK8 vs ZZ pt")->Write();
	
	// Resolution of GEN 4
	vector<double> std_4g;
	vector<double> std_err_4g;
	vector<double> bins_4g;
	vector<double> bins_w_4g;
	for(int x = 1; x <= genAK4_ZZpt->GetNbinsX(); ++x){
		TH1D* px = genAK4_ZZpt->ProjectionY(Form("px_%d",x), x,x+1, "e");
		double bcx = ((TAxis*)genAK4_ZZpt->GetXaxis())->GetBinCenter(x);
		double bwx = ((TAxis*)genAK4_ZZpt->GetXaxis())->GetBinWidth(x);
		bins_4g.push_back(bcx);
		bins_w_4g.push_back(bwx/sqrt(12));
		std_4g.push_back(px->GetStdDev());
		std_err_4g.push_back(px->GetStdDevError());
	}
	TGraphErrors stdDev4g(bins_4g.size(), bins_4g.data(), std_4g.data(),  bins_w_4g.data(), std_err_4g.data());
	stdDev4g.SetName("Resolution AK4 with error");
	stdDev4g.SetTitle("Resolution AK4 vs ZZ pt;ZZ pt [GeV/c];Resolution [GeV/c^{2}]");
	stdDev4g.SetMinimum(0.);
	stdDev4g.Write();
	
	TProfile* prof_4g = genAK4_ZZpt->ProfileX("Precision AK4 vs ZZ pt");
	prof_4g->SetTitle("Precision AK4 vs ZZ pt;ZZ pt [GeV/c];Precision [GeV/c^{2}]");
	prof_4g->Write();
	
	
	// ----- Ratio of AK8 to AK4-resolved gen bosons vs s
	TH2* genType_s = dynamic_cast<TH2*>(theHistograms->get("Best Z gen vs s"));
	if(!genType_s) return; //TGraphErrors or TGraphAsymmErrors?
	TH1D* g4_s = genType_s->ProjectionX("AK4_proj", 1,2);
	TH1D* g8_s = genType_s->ProjectionX("AK8_proj", 2,3);
	g4_s->Rebin(2);
	g8_s->Rebin(2);
	TH1D* sum_s = (TH1D*)g4_s->Clone("total");
	sum_s->Add(g8_s);
	
	TGraphAsymmErrors ratio_g8_s(g8_s, sum_s);
	ratio_g8_s.SetName("Ratio of best Z AK8 vs s"); 
	ratio_g8_s.SetTitle("Ratio of best Z AK8;s [GeV/c^{2}]");
	ratio_g8_s.SetMinimum(0.);
	ratio_g8_s.Write();
	TGraphAsymmErrors ratio_g4_s(g4_s, sum_s);
	ratio_g4_s.SetName("Ratio of best Z AK4 vs s"); 
	ratio_g4_s.SetTitle("Ratio of best Z AK4;s [GeV/c^{2}]");
	ratio_g4_s.SetMinimum(0.);
	ratio_g4_s.Write();
	
	// ----- Ratio of AK8 to AK4-resolved gen bosons vs ZZ pt
	TH2* genType_pt = dynamic_cast<TH2*>(theHistograms->get("Best Z gen vs ZZ pt"));
	if(!genType_pt) return; //TGraphErrors or TGraphAsymmErrors?
	TH1D* g4_pt = genType_pt->ProjectionX("AK4_proj_pt", 1,2);
	TH1D* g8_pt = genType_pt->ProjectionX("AK8_proj_pt", 2,3);
	g4_pt->Rebin(2);
	g8_pt->Rebin(2);
	TH1D* sum_pt = (TH1D*)g4_pt->Clone("total_pt");
	sum_pt->Add(g8_pt);
	
	TGraphAsymmErrors ratio_g8_pt(g8_pt, sum_pt);
	ratio_g8_pt.SetName("Ratio of best Z AK8 vs ZZ pt");
	ratio_g8_pt.SetTitle("Ratio of best Z AK8;ZZ pt [GeV/c]");
	ratio_g8_pt.SetMinimum(0.);
	ratio_g8_pt.Write();
	TGraphAsymmErrors ratio_g4_pt(g4_pt, sum_pt);
	ratio_g4_pt.SetName("Ratio of best Z AK4 vs ZZ pt"); 
	ratio_g4_pt.SetTitle("Ratio of best Z AK4;ZZ pt [GeV/c]");
	ratio_g4_pt.SetMinimum(0.);
	ratio_g4_pt.Write();
}


void VZZAnalyzer::simpleGraphs(){
	foreach(const Jet& j, *jetsAK8){
		theHistograms->fill("All8_t21_f_", ";#tau_{2}/#tau_{1};# jets", 25,0.,1., j.tau2()/j.tau1(), theWeight);
		theHistograms->fill("All8_t32_f_", ";#tau_{3}/#tau_{2};# jets", 25,0.,1., j.tau3()/j.tau2(), theWeight);
		theHistograms->fill("All8_PUPPIt21_f_", ";PUPPI #tau_{2}/#tau_{1};# jets", 25,0.,1., j.puppiTau2()/j.puppiTau1(), theWeight);
	}
	
	foreach(const Boson<Particle>& b, *genHadVBs_)
		theHistograms->fill("All AK4_{gen} mass", "All AK4_{gen} mass", 50,0.,200., b.mass(), 1.);
	
	/*
	theHistograms->fill("genVBParticles size","genVBParticles size", 5,-0.5,4.5, genVBParticles->size(), theWeight);
	theHistograms->fill("genHadVBs_ size", "genHadVBs_ size", 5,-0.5,4.5, genHadVBs_->size(), theWeight);
	if(genHadVBs_->size() > 0){
		theHistograms->fill("genHadVB0_daugh deltaR", "genHadVB #DeltaR between daughters", 40,0.,4., physmath::deltaR(genHadVBs_->front().daughter(0), genHadVBs_->front().daughter(1)), theWeight);
		theHistograms->fill("genHadVB0_daugh deltaPhi", "genHadVB #Delta#phi between daughters", 40,0.,4., physmath::deltaPhi(genHadVBs_->front().daughter(0), genHadVBs_->front().daughter(1)), theWeight);
		theHistograms->fill("genHadVB0_daugh deltaEta", "genHadVB #Delta#eta between daughters", 40,0.,4., fabs(genHadVBs_->front().daughter(0).eta() - genHadVBs_->front().daughter(1).eta()), theWeight);
	}
	unsigned int nW = 0;
	unsigned int nZ = 0;
	foreach(const phys::Boson<phys::Particle>& gen, *genVBParticles){
    if(abs(gen.id()) == 23 && abs(gen.daughter(0).id()) < 10) {
      theHistograms->fill("genMass_Z", "mass of gen Z", 35,50.,120., gen.mass(), theWeight);
      ++nZ;
    }
    else if(abs(gen.id()) == 24 && abs(gen.daughter(0).id()) < 10) {
      theHistograms->fill("genMass_W", "mass of gen W", 35,50.,120., gen.mass(), theWeight);
      ++nW;
    }
    theHistograms->fill("number of gen W", "number of gen W", 4,-0.5,3.5, nW, 1.);
    theHistograms->fill("number of gen Z", "number of gen Z", 4,-0.5,3.5, nZ, 1.);    
  }
	
	theHistograms->fill("Size AK4_{rec}","Size AK4_{rec}", 8,-0.5,7.5, jets->size(), theWeight);
	theHistograms->fill("Size AK8_{rec}","Size AK8_{rec}", 8,-0.5,7.5, jetsAK8->size(), theWeight);
	theHistograms->fill("Size AK4_{gen}","Size AK4_{gen}", 8,-0.5,7.5, genJets->size(), theWeight);
	theHistograms->fill("Size AK8_{gen}","Size AK8_{gen}", 8,-0.5,7.5, genJetsAK8->size(), theWeight);
	*/
	foreach(const Particle& jet, *genJetsAK8)
		theHistograms->fill("genAK8 pt", "AK8_{gen} pt;pt [GeV]", 50,20.,120.,jet.pt(), theWeight);
	foreach(const Jet& jet, *jetsAK8)
		theHistograms->fill("recAK8 pt", "AK8_{rec} pt;pt [GeV]",50,100.,200.,jet.pt(), theWeight);
	foreach(const Particle& jet, *genJets)
		theHistograms->fill("genAK4 pt", "AK4_{gen} pt;pt [GeV]", 50,20.,70., jet.pt(), theWeight);
	foreach(const Jet& jet, *jets)
		theHistograms->fill("recAK4 pt", "AK4_{rec} pt;pt [GeV]", 50,20.,70., jet.pt(), theWeight);
	
	//foreach(const Jet& jet, *jets)
		//theHistograms->fill("All AK4_{rec} mass", "All AK4_{rec} mass", 24,-0.,120., jet.mass(), theWeight);
	
	//foreach(const Jet& jet, *jetsAK8){
		//theHistograms->fill("All AK8_{rec} chosenAlgoMass", "All AK8_{rec} chosenAlgoMass;mass [GeV/c^2]", 24,-0.,120., jet.chosenAlgoMass(), theWeight);
		/*if(ZZ != nullptr)
			if(ZZ->p() > 1.)
				theHistograms->fill<TH2F>("All AK8_{rec} mass vs ZZ pt", "All AK8_{rec} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), jet.mass(), 1.);*/
	//}
	
	//foreach(const Particle& jet, *genJets)
		//theHistograms->fill("All AK4_{gen} mass", "All AK4_{gen} mass", 24,-0.,120., jet.mass(), theWeight);
	
	//foreach(const Particle& jet, *genJetsAK8){
		//theHistograms->fill("All AK8_{gen} mass", "All AK8_{gen} mass", 24,-0.,120., jet.mass(), theWeight);
		/*if(ZZ != nullptr)
			if(ZZ->p() > 1.)
				theHistograms->fill<TH2F>("All AK8_{gen} mass vs ZZ pt", "All AK8_{gen} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), jet.mass(), 1.);*/
	//}
}


void VZZAnalyzer::ZZGraphs(){
	if(!ZZ || ZZ->p() < 1.)
		return;
	
	theHistograms->fill("ZZ_rec mass", "ZZ_{rec} mass (m4l)", 30,60.,360., ZZ->mass(), theWeight);
	theHistograms->fill("ZZ_rec pt", "ZZ_{rec} pt", 30,0.,240., ZZ->pt(), theWeight);
	
	theHistograms->fill("Z0_rec mass", "Z0_{rec} mass", 30,60.,120., ZZ->first().mass(), theWeight);
	theHistograms->fill("Z1_rec mass", "Z1_{rec} mass", 30,60.,120., ZZ->second().mass(), theWeight);
	
	theHistograms->fill("l00_rec pt", "l00_{rec} mass", 30,0.,240., ZZ->first().daughter(0).pt(), theWeight);
	theHistograms->fill("l01_rec pt", "l01_{rec} mass", 30,0.,240., ZZ->first().daughter(1).pt(), theWeight);
	theHistograms->fill("l10_rec pt", "l10_{rec} mass", 30,0.,240., ZZ->second().daughter(0).pt(), theWeight);
	theHistograms->fill("l11_rec pt", "l11_{rec} mass", 30,0.,240., ZZ->second().daughter(1).pt(), theWeight);
}


template <class J = phys::Jet>
phys::Boson<J>* VZZAnalyzer::findBestVFromPair(const std::vector<J>* js	){
	if(js->size() < 2)
		return nullptr;
		
	pair<size_t, size_t> indicesZ(0,0);
	pair<size_t, size_t> indicesW(0,0);
	float minDifZ = 50.;
	float minDifW = 50.;
	float tmpMass = 0.;
	for(size_t i = 0; i < js->size(); ++i){
		for(size_t j = i+1; j < js->size(); ++j){
			tmpMass = (js->at(i).p4() + js->at(j).p4()).M();
			float diffZa = fabs(tmpMass - phys::ZMASS);
			float diffWa = fabs(tmpMass - phys::WMASS);
			if(diffZa < minDifZ){
				minDifZ = diffZa;
				indicesZ = std::make_pair(i,j);
			}
			if(diffWa < minDifW){
				minDifW = diffWa;
				indicesW = std::make_pair(i,j);
			}
		}
	}
	
	phys::Boson<J>* thisCandidate = nullptr;
	if(minDifZ < minDifW){
		thisCandidate = new Boson<J>(js->at(indicesZ.first), js->at(indicesZ.second));
		if(!ZBosonDefinition(*thisCandidate))
			return nullptr;
	}
	else{
		thisCandidate = new Boson<J>(js->at(indicesW.first), js->at(indicesW.second));
		if(!WBosonDefinition(*thisCandidate))
			return nullptr;
	}
	return thisCandidate;
}


template<class P = phys::Jet>
const P* VZZAnalyzer::findBestVFromSing(std::vector<P>* js){
	if(js->size() < 1)
		return nullptr;
	
	const P* thisCandidate = nullptr;
	auto it_best = std::min_element(js->begin(), js->end(), Mass2Comparator(phys::WMASS, phys::ZMASS));
		
	if(physmath::minDM(getRefinedMass(*it_best)) < 30.){
		thisCandidate = &(*it_best);
	}
	
	return thisCandidate;
}


template <class P = phys::Particle>
const P* VZZAnalyzer::findBestVPoint(std::vector<const P*>& js){
	if(js.size() < 1)
		return nullptr;
	size_t indexZ = 0;
	size_t indexW = 0;
	float minDifZ = 50.;
	float minDifW = 50.;
	float tmpMass = 0.;
	for(size_t i = 0; i < js.size(); ++i){
		tmpMass = getRefinedMass(js.at(i));
		float diffZa = fabs(tmpMass - phys::ZMASS);
		float diffWa = fabs(tmpMass - phys::WMASS);
		if(diffZa < minDifZ){
			minDifZ = diffZa;
			indexZ = i;
		}
		if(diffWa < minDifW){
			minDifW = diffWa;
			indexW = i;
		}
	}
	
	const P* thisCandidate = nullptr;
	if(minDifZ < minDifW){
		thisCandidate = js.at(indexZ);
		if(!ZBosonDefinition(*thisCandidate))
			return nullptr;
	}
	else{
		thisCandidate = js.at(indexW);
		if(!WBosonDefinition(*thisCandidate))
			return nullptr;
	}
	return thisCandidate;
}

#ifdef USE_PYTHON
int VZZAnalyzer::initPy(){
	Py_Initialize();
	if(! Py_IsInitialized() ){
		std::cerr<<"initPy: Py interpreter not initialized"<<std::endl;
		return 1;
	}
	
	PyRun_SimpleString(
	 "print('Importing sys...')\n"
   "import sys\n"
   "print('Python sys.version =', sys.version)\n"
   "sys.path.append('./python')\n"
	);
	//cout<<"Importing sklearn...\n";
	//PyObject * skModule = PyImport_ImportModule("sklearn");
	//cout<<"Succesul import\n";
	
	//Py_Finalize(); return(0);
	
	helper_module_ = PyImport_ImportModule("VZZhelper"); // import module
	if (!helper_module_){
		std::cerr<<"Error: could not load \"VZZhelper\"."<<std::endl;
		Py_Finalize();
		return(2);
	}
	
	AK4_classifier_ = PyObject_CallMethod(helper_module_, (char*)"load_object", (char*)"s", (char*)"predictors/VZZ_AK4_tree.pkl");
	if(!AK4_classifier_ || AK4_classifier_ == Py_None){
		std::cerr<<"Error: could not load AK4_classifier_."<<std::endl;
		Py_DECREF(helper_module_);
		Py_Finalize();
		return(3);
	}
	
	AK8_classifier_ = PyObject_CallMethod(helper_module_, (char*)"load_object", (char*)"s", (char*)"predictors/VZZ_AK8_tree.pkl");
	if(!AK8_classifier_ || AK8_classifier_ == Py_None){
		std::cerr<<"Error: could not load AK8_classifier_."<<std::endl;
		Py_DECREF(AK4_classifier_);
		Py_DECREF(helper_module_);
		Py_Finalize();
		return(3);
	}
	
	EVT_classifier_ = PyObject_CallMethod(helper_module_, (char*)"load_object", (char*)"s", (char*)"predictors/VZZ_Evt_tree.pkl");
	if(!EVT_classifier_ || EVT_classifier_ == Py_None){
		std::cerr<<"Error: could not load VZZ_Evt_tree."<<std::endl;
		Py_DECREF(AK4_classifier_);
		Py_DECREF(AK8_classifier_);
		Py_DECREF(helper_module_);
		Py_Finalize();
		return(3);
	}
	
	_import_array();
	
	return 0;
}


template<class P>
vector<double>* VZZAnalyzer::getPyPredAll(const vector<P>& vect, PyObject* predictor) const {
	
	size_t vsize = vect.size();
	if(vsize == 0)
		return new vector<double>();
	
	auto temp = getFeatures(vect.front());
	size_t nfeat = temp->size();
	delete temp;
	
	double* data = new double[vsize*nfeat];
	for(size_t i = 0; i < vsize; ++i){
		vector<double>* feat = getFeatures(vect.at(i));
		std::copy( feat->begin(), feat->end(), &(data[i*nfeat]) );
		delete feat;
	}
	
	vector<double>* preds = predAll(data, vsize, nfeat, predictor);
	
	delete[] data;
	
	return preds;
}


vector<double>* VZZAnalyzer::predAll(double* data, size_t vsize, size_t nfeat, PyObject* predictor) const{
	npy_intp dims[2] = {(long) vsize, (long)nfeat};
	PyObject* ndarray = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, data);
	
	PyObject* raw_preds = PyObject_CallMethod(predictor, (char*)"predict_proba", (char*)"O", ndarray);
	double *p = (double*)PyArray_DATA((PyArrayObject*) raw_preds);
	
	vector<double>* preds = new vector<double>(vsize, -2.);
	for(size_t i = 0; i < vsize; ++i){
		preds->at(i) = p[2*i+1];
	}
	
	Py_DECREF(raw_preds);
	Py_DECREF(ndarray);
	
	return preds;
}


double VZZAnalyzer::getPyPrediction(const vector<double>& vect, PyObject* predictor){
	if(!predictor){
		std::cerr<<"Warning (getPyPrediction): predictor is NULL"<<std::endl;
		return -2.;
	}
	
	size_t nfeat = vect.size();
	double* data = new double[nfeat];
	std::copy( vect.begin(), vect.end(), data );
	
	npy_intp dims[2] = {(long)1, (long)nfeat};
	PyObject* ndarray = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, data);
	PyObject* raw_preds = PyObject_CallMethod(predictor, (char*)"predict_proba", (char*)"O", ndarray);
	double *p = (double*)PyArray_DATA((PyArrayObject*) raw_preds);
	
	Py_DECREF(raw_preds);
	Py_DECREF(ndarray);
	delete[] data;
	return p[1];
}
#endif

vector<double>* VZZAnalyzer::getAK4features(const Boson<Jet>& jj){
	vector<double>* buffer = new vector<double>;
	buffer->reserve(19);
	const Jet& d0 = jj.daughter(0);
	const Jet& d1 = jj.daughter(1);
	
	buffer->push_back(physmath::deltaR(d0, d1));
	buffer->push_back(jj.pt());
	buffer->push_back(d0.pt() + d1.pt());
	buffer->push_back(jj.mass());
	buffer->push_back(physmath::minDM(jj.mass()));     // 5
	
	//Still not included
	buffer->push_back(d0.mass());
	buffer->push_back(d1.mass());            // 7
	buffer->push_back(d0.qgLikelihood());
	buffer->push_back(d1.qgLikelihood());    // 9
	buffer->push_back(d0.csvtagger());
	buffer->push_back(d1.csvtagger());       // 11
	buffer->push_back(d0.girth());
	buffer->push_back(d1.girth());           // 13
	buffer->push_back(d0.girth_charged());
	buffer->push_back(d1.girth_charged());   // 15
	buffer->push_back(d0.jetArea());
	buffer->push_back(d1.jetArea());         // 17
	buffer->push_back(d0.ptd());
	buffer->push_back(d1.ptd());             //19 (+14)
	
	// Useless variables (values are always default) for AK4:
	// tau(1,2,3), puppiTau(1,2,3), mass functions except mass(), ptd, secvtxMass
	return buffer;
}


vector<double>* VZZAnalyzer::getAK8features(const Jet& j){
	vector<double>* buffer = new vector<double>();
	buffer->reserve(15);
	
	buffer->push_back(j.pt());
	buffer->push_back(j.chosenAlgoMass());  //softDropMass_
	buffer->push_back(physmath::minDM(j.chosenAlgoMass()));
	buffer->push_back(j.corrPrunedMass());
	buffer->push_back(j.prunedMass());      // 5
	
	//float tau21 = j.tau2()/j.tau1();  // This should be high
	//float tau32 = j.tau3()/j.tau2();  // This should be low
	
	buffer->push_back(j.tau1());
	buffer->push_back(j.tau2());
	buffer->push_back(j.tau3());            //8
	buffer->push_back(j.puppiTau1());
	buffer->push_back(j.puppiTau2());
	buffer->push_back(j.puppiTau3());       //11
	buffer->push_back(j.girth());
	buffer->push_back(j.girth_charged());   //13
	buffer->push_back(j.jetArea());
	buffer->push_back(j.ptd());             //15
	
	// Useless variables (values are always default) for AK8:
	// qgLikelihood, csvtagger
	return buffer;
}


void VZZAnalyzer::fillAK4GenRec(bool doGraphs){
	if(genHadVBs_->size() == 0)
		return;
	
	foreach(const Boson<Particle>& genVB, *genHadVBs_){
		//++Nhad_genVB_;
		std::stable_sort( jets->begin(), jets->end(), DeltaRComparator(genVB.daughter(0)) );
		Jet rec0 = jets->front();  // Copy for future use
		float dR0 = physmath::deltaR(rec0, genVB.daughter(0));
		if(doGraphs && dR0 < 0.4){
			theHistograms->fill("Resolution CMS-rec AK4: #DeltaR", "Resolution CMS-rec AK4: #DeltaR", 20,0.,0.4, dR0);
			float dM0 = rec0.mass() - genVB.daughter(0).mass();
			theHistograms->fill("Resolution CMS-rec AK4: #DeltaM", "Resolution of CMS-rec: #DeltaM", 31,-15.5,15.5, dM0);
			theHistograms->fill<TH2F>("Resolution CMS-rec AK4: #DeltaM vs pt", "Resolution CMS-rec AK4: #DeltaM vs pt", 10,20.,220., 13,-16.25,16.25, genVB.daughter(0).pt(), dM0, 1.);
		}
		
		std::stable_sort( jets->begin(), jets->end(), DeltaRComparator(genVB.daughter(1)) );
		Jet rec1 = jets->front();  // Copy for future use
		float dR1 = physmath::deltaR(rec1, genVB.daughter(1));
		if(doGraphs && dR1 < 0.4){
			theHistograms->fill("Resolution CMS-rec AK4: #DeltaR", "Resolution CMS-rec AK4: #DeltaR", 20,0.,0.4, dR1);
			float dM1 = rec1.mass() - genVB.daughter(1).mass();
			theHistograms->fill("Resolution CMS-rec AK4: #DeltaM", "Resolution CMS-rec AK4: #DeltaM", 31,-15.5,15.5, dM1);
			theHistograms->fill<TH2F>("Resolution CMS-rec AK4: #DeltaM vs pt", "Resolution CMS-rec AK4: #DeltaM vs pt", 10,20.,220., 13,-16.25,16.25, genVB.daughter(1).pt(), dM1, 1.);
		}
		
		// Gen VB can be correctly reconstructed
		if(dR0 < 0.2 && dR1 < 0.2){
			//++Ncms_recVB_;
			Boson<Jet> bestPossibleAK4(rec0, rec1);
			//float dR = physmath::deltaR(bestPossibleAK4, genVB);
			// reconstructible --> reconstructible using detector-reconstructed particles
			// therefore "CMS-rec" --> best possible reconstructible by algorithm: the ideal algorithm would always choose this as the VB-candidate, as its constituent jet are the reconstructed counterparts to the genJets of the genVB
			if(doGraphs){
				theHistograms->fill("Resolution CMS-rec VB: #DeltaR", "Resolution CMS-rec VB: #DeltaR", 20,0.,0.4, physmath::deltaR(bestPossibleAK4, genVB));
				float dM = bestPossibleAK4.mass() - genVB.mass();
				theHistograms->fill("Resolution CMS-rec VB: #DeltaM", "Resolution CMS-rec VB: #DeltaM", 31,-31.,31., dM);
				theHistograms->fill<TH2F>("Resolution CMS-rec VB: #DeltaM vs pt", "Resolution CMS-rec VB: #DeltaM vs pt", PT_2D_SIZE, 13,-32.5,32.5, genVB.pt(), dM, 1.);
			}
			AK4GenRec_->push_back(make_pair(genVB, bestPossibleAK4));
		}
	}
}


void VZZAnalyzer::fillGenHadVBs(){
	foreach(const Boson<Particle>& genVB, *genVBParticles)
		if( abs(genVB.daughter(0).id()) < 10 && abs(genVB.daughter(1).id()) < 10 )
			genHadVBs_->push_back(genVB); //The second condition is redundant
}


void VZZAnalyzer::fillRecHadVBs(){
	for(size_t i = 0; i < jets->size(); ++i)
		for(size_t j = i+1; j < jets->size(); ++j)
			AK4pairs_->push_back( Boson<Jet>(jets->at(i), jets->at(j)) );
}


