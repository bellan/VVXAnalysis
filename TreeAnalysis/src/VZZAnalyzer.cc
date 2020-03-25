#include "VVXAnalysis/TreeAnalysis/interface/VZZAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/DataFormats/interface/Jet.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/TypeDefs.h"

#include <iostream>
#include <fstream>			//open(), close(), <<
#include <string>				//find_last_of()
#include <time.h>				//clock_t, clock()
#include <utility>			//std::pair(), std::make_pair()
#include "TSystem.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1.h"
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
#define MASS_2D_SIZE 16,50.,130.
#define P_2D_SIZE 25,0.,1000.
#define PT_2D_SIZE 18,0.,540.
#define ZZMASS_2D_SIZE 30,100.,700.
#define S_2D_SIZE 30,100.,1600.
#define BINARY_SIZE 2,-1.5,1.5

using namespace boost::assign;

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;
using std::pair;
using std::make_pair;

using namespace phys;


void VZZAnalyzer::begin(){
	cout<<"\n";
	for(char i=0; i<25; ++i) cout<<"-";
	cout<<" \tStart of VZZAnalyzer\t";
	for(char i=0; i<25; ++i) cout<<"-";
	
	startTime_ = clock();
	
	//Memory allocation
	if(AK4pairs_ == nullptr)   AK4pairs_ = new vector<Boson<Jet>>;
	if(genHadVBs_ == nullptr) genHadVBs_ = new vector<Boson<Particle>>;
	
	cout<<"\nAnalyzed:\t      /"<<tree()->GetEntries()<<std::flush;
	return;
}


Int_t VZZAnalyzer::cut(){
	++evtN_;
	//cout<<"\r\t\t"<<evtN_;
	
	//Cleanup of previous event
	if(genHadVBs_) genHadVBs_->clear(); //Destroys objects but keeps memory allocated
	if(AK4pairs_)   AK4pairs_->clear();
	
	
	if(!ZZ)
		return -1;
	else if(ZZ->pt() < 3.)
		return -1;
	
	//Preparation for this event
	fillGenHadVBs();
	fillRecHadVBs();
	
	int sigType = isSignal(/*Uses genJets(AK4) and genJetsAK8*/);
	theHistograms.fill("Signal Type", "Signal Type", 3,-0.5,2.5, sigType);
	if(sigType) return 1;
	else return -1;
}


void VZZAnalyzer::analyze(){
	++analyzedN_;
	
	calcS();
	//bestCandidateAnalysis();
	//ptCutMVA();
	//closestJetAnalisys();
	//furthestJetMVA();
	//minPtJetMVA();
	
	bestZMassJetMVA();
	simpleGraphs();
	
	
	foreach(const Boson<Particle>& b, *genHadVBs_)
		theHistograms.fill("All AK4_{gen} mass", "All AK4_{gen} mass", 50,0.,200., b.mass(), 1./*theWeight*/);
	/*
	theHistograms.fill("s AK4 gen", "s AK4 gen", 100,0.,2000., sAK4g_, theWeight);
	theHistograms.fill("s AK8 gen", "s AK8 gen", 100,0.,2000., sAK8g_, theWeight);
	theHistograms.fill("s AK4 rec", "s AK4 rec", 100,0.,2000., sAK4r_, theWeight);
	theHistograms.fill("s AK8 rec", "s AK8 rec", 100,0.,2000., sAK8r_, theWeight);
	*/
	
	return;
}

void VZZAnalyzer::end(TFile & fout){
	//Final cleanup
	if(genHadVBs_) delete genHadVBs_; //Deallocate memory
	if(AK4pairs_)  delete AK4pairs_;
	
	cout<<'\n';
	if(withGenVB_) //This counter was incremented if bestCandidateAnalysis() was called
		endBestCandAnalysis(fout); //Dividing histograms to get efficiency and stuff
	
	if(win4_ || win8_) //This counter was incremented if closestJetAnalisys() was called
		endClosestJetAn();
	
	
	if(withGenVB_ /*if the candJetAnalysis has been done*/){
		cout<<"Candidate W: JetPair = "<<pairWFromJets_<<" \tJetSing = "<<singWFromJets_<<" \tAK8Pair = "<<pairWFromJetsAK8_<<" \tAK8Sing = "<<singWFromJetsAK8_;
		cout<<"\nCandidate Z: JetPair = "<<pairZFromJets_<<" \tJetSing = "<<singZFromJets_<<" \tAK8Pair = "<<pairZFromJetsAK8_<<" \tAK8Sing = "<<singZFromJetsAK8_;
		cout<<"\nEvents with a reconstructed VB: "<<recVBtot_<<" -> "<<recVBtot_*100./evtN_<<'%';
		cout<<"\nReconstructed near (deltaR < 0.4) a genVB: "<<goodRec_<<" \tevents with a genVB: "<<withGenVB_<<" \t--> "<<goodRec_*100./withGenVB_<<'%';
	}
	if(win8_ || win4_) //if closestJetAnalisys has been done
		cout<<"Reconstruction of a VB --> AK8: "<<win8_<<" \tAK4 pairs: "<<win4_;
	
	endResolutionAnalisys(fout);
	
	cout<<"\nAnalized: "<<analyzedN_<<'\n';
	
	float elapsedSec = (float)(clock()-startTime_)/CLOCKS_PER_SEC;
	int elapsedSecInt = (int)elapsedSec;
	cout<<"\nElapsed Time: "<<elapsedSec<<" s\t\t("<<elapsedSecInt/60<<"\' "<<elapsedSecInt%60<<"\")\n";
	for(char i=0; i<25; ++i) cout<<"-";
	cout<<" \tEnd of VZZAnalyzer\t";
	for(char i=0; i<25; ++i) cout<<"-";
	cout<<"\n\n";
}


int VZZAnalyzer::isSignal() const{
	bool twoAK4 = ( genHadVBs_->size() >= 2 );
	bool oneAK8 = ( genJetsAK8->size() >= 1 );
	if(!twoAK4 || !oneAK8)
		return 0;
	
	std::stable_sort(genHadVBs_->begin(), genHadVBs_->end(), MassComparator(phys::ZMASS));
	float mass4 = genHadVBs_->front().mass();
	if(60. < mass4 && mass4 < 120.)
		return 1;
	
	std::stable_sort(genJetsAK8->begin(), genJetsAK8->end(), MassComparator(phys::ZMASS));
	float mass8 = genJetsAK8->front().mass();
	if(60. < mass4 && mass4 < 120.)
		return 2;
		
	return 0;
}


void VZZAnalyzer::specialPeakAnalisys(const Particle& theGenAK8){
	//fillGenHadVBs();
	static int n_calls = 1;
	stable_sort(AK4pairs_->begin(), AK4pairs_->end(), DeltaRComparator(theGenAK8) );
	/*
	cout<<"\n---------------------------------------- "<<n_calls<<'\n';
	cout<<"\t#DeltaR(AK8, B) =   "<<physmath::deltaR( AK4pairs_->front(), theGenAK8 )<<'\n';
	cout<<"\t#DeltaPhi(AK8, B) = "<<physmath::deltaPhi( AK4pairs_->front(), theGenAK8 )*180. /M_PI <<"°\n";
	cout<<"\tAK4 #DeltaR =       "<<physmath::deltaR( AK4pairs_->front().daughter(0), AK4pairs_->front().daughter(1) )<<'\n';
	cout<<"\t#DeltaM(AK8, B) =   "<<AK4pairs_->front().mass() - theGenAK8.mass()<<'\n';
	cout<<"\t#DeltaP(AK8, B) =   "<<(AK4pairs_->front().p4() - theGenAK8.p4()).P()<<'\n';
	*/
	
	if(physmath::deltaR(AK4pairs_->front(), theGenAK8) < 0.4){
		theHistograms.fill("Special Peak: AK4 #DeltaR", "Special Peak: AK4 #DeltaR", 10,0.,2., physmath::deltaR( AK4pairs_->front().daughter(0), AK4pairs_->front().daughter(1) ) , 1.);
		theHistograms.fill("Special Peak: #DeltaR(AK8, Boson)", "Special Peak: #DeltaR(AK8, Boson)", 8,0.,0., physmath::deltaR( AK4pairs_->front(), theGenAK8 ) , 1.);
		theHistograms.fill("Special Peak: #DeltaM(AK8, Boson)", "Special Peak: #DeltaM(AK8, Boson)", 10,-20.,20., AK4pairs_->front().mass() - theGenAK8.mass(), 1.);
		theHistograms.fill("Special Peak: #DeltaP(AK8, Boson)", "Special Peak: #DeltaP(AK8, Boson)", 10,0.,50., (AK4pairs_->front().p4() - theGenAK8.p4()).P() , 1.);
	}
	else //Counting events without AK4 pairs
		theHistograms.fill("Special Peak: AK4 #DeltaR", "Special Peak: AK4 #DeltaR", 10,0.,2., -100., 1.);
		
	++n_calls;
	return;
}


void VZZAnalyzer::bestZMassJetMVA(){
	if(!ZZ) return;
	if(ZZ->pt() < 1.) return;
	
	//-----------------------------	GEN PARTICLES	-----------------------------
	//------------------	GEN AK8	------------------
	auto it_8g = genJetsAK8->begin();
	bool found_8g = false;
	if(genJetsAK8->size() > 0){
		stable_sort( genJetsAK8->begin(), genJetsAK8->end(), MassComparator(phys::ZMASS) );
		while(it_8g != genJetsAK8->end()){
			if(physmath::deltaR(*it_8g, *ZZ) > 2. && ZBosonDefinition(*it_8g)){
				found_8g = true;
				break;
			}
			else ++it_8g;
		}
		if(found_8g){
			theHistograms.fill("bestZ AK8_{gen} mass", "bestZ AK8_{gen} mass", 30,0.,150, it_8g->mass(), 1.);
			theHistograms.fill<TH2F>("BestZ AK8_{gen} mass vs ZZ pt", "BestZ AK8_{gen} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), it_8g->mass(), 1.);
			theHistograms.fill<TH2F>("BestZ AK8_{gen} mass vs ZZ mass", "BestZ AK8_{gen} mass vs ZZ mass", ZZMASS_2D_SIZE, MASS_2D_SIZE, ZZ->mass(), it_8g->mass(), 1.);
			theHistograms.fill<TH2F>("BestZ AK8_{gen} mass vs s", "BestZ AK8_{gen} mass vs s", S_2D_SIZE, MASS_2D_SIZE, sAK8g_, it_8g->mass());
			
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
		stable_sort( genHadVBs_->begin(), genHadVBs_->end(), MassComparator(phys::ZMASS) );
		while(it_4g != genHadVBs_->end()){
			if(physmath::deltaR(*it_4g, *ZZ) > 2. && ZBosonDefinition(*it_4g)){
				found_4g = true;
				break;
			}
			else ++it_4g;
		}
		if(found_4g){
			theHistograms.fill("bestZ AK4_{gen} mass", "bestZ AK4_{gen} mass", 30,0.,150, it_4g->mass(), 1.);
			theHistograms.fill<TH2F>("BestZ AK4_{gen} mass vs ZZ pt", "BestZ AK4_{gen} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), it_4g->mass(), 1.);
			theHistograms.fill<TH2F>("BestZ AK4_{gen} mass vs ZZ mass", "BestZ AK4_{gen} mass vs ZZ mass", ZZMASS_2D_SIZE, MASS_2D_SIZE, ZZ->mass(), it_4g->mass(), 1.);
			theHistograms.fill<TH2F>("BestZ AK4_{gen} mass vs s", "BestZ AK4_{gen} mass vs s", S_2D_SIZE, MASS_2D_SIZE, sAK4g_, it_4g->mass());
		}
	}
	
	
	//------------------	WINNER GEN	------------------
	if(found_4g || found_8g){
		float type = 0.;  // -1 --> AK4 pair,   1 --> AK8
		if(found_4g && found_8g){
			type = (fabs(it_8g->mass() -phys:: ZMASS) < fabs(it_4g->mass() -phys:: ZMASS) ? 1.:-1.);
			theHistograms.fill<TH2F>("Winner Z mass", "Winner Z mass", MASS_2D_SIZE, BINARY_SIZE, (type > 0.5 ? *it_8g : *it_4g).mass(), type);
			theHistograms.fill<TH2F>("Loser Z mass", "Loser Z mass", MASS_2D_SIZE, BINARY_SIZE, (type < 0.5 ? *it_8g : *it_4g).mass(), type);
		}
		else if(found_8g)
			type = 1.;
		else if(found_4g)
			type = -1.;
		
		double sHat_g = ( ZZ->p4() + (type > 0.5 ? *it_8g : *it_4g).p4() ).M();
		theHistograms.fill<TH2F>("Best Z gen vs #hat{s}", "Best Z gen vs #hat{s}", S_2D_SIZE, BINARY_SIZE, sHat_g, type);
		theHistograms.fill<TH2F>("Best Z gen vs s", "Best Z gen vs s", S_2D_SIZE, BINARY_SIZE, (type > 0.5 ? sAK8g_ : sAK4g_), type);
		theHistograms.fill<TH2F>("Best Z gen vs ZZ pt", "Best Z gen vs ZZ pt", PT_2D_SIZE, BINARY_SIZE, ZZ->pt(), type);
		theHistograms.fill<TH2F>("Best Z gen vs #Delta#phi(ZZ, Jet)_div_pi", "Best Z gen vs #Delta#phi(ZZ, Jet)_div_#pi", 40,-1.,1., BINARY_SIZE, physmath::deltaPhi( *ZZ, (type > 0.5 ? *it_8g : *it_4g) )/M_PI, type); //32,-3.2,3.2
		theHistograms.fill<TH2F>("Best Z gen vs #Delta#eta(ZZ, Jet)", "Best Z gen vs #Delta#eta(ZZ, Jet)", 20,-4.,4., BINARY_SIZE, ZZ->eta()-(type > 0.5 ? *it_8g : *it_4g).eta(), type);
		theHistograms.fill<TH2F>("Best Z gen vs #DeltaR(ZZ, Jet)", "Best Z gen vs #DeltaR(ZZ, Jet)", 40,2.,6., BINARY_SIZE, physmath::deltaR( *ZZ, (type > 0.5 ? *it_8g : *it_4g) ), type);
		
		theHistograms.fill<TH2F>("Best Z vs deltaR(4_{gen})", "Best Z vs deltaR(4_{gen})", 50,0.,5., BINARY_SIZE, (found_4g ? physmath::deltaR( it_4g->daughter(0), it_4g->daughter(1) ) : -1.), type /*filling underflow bin --> bestZ found with AK8*/ );
	theHistograms.fill<TH2F>("Best Z vs deltaPhi(4_{gen})","Best Z vs deltaPhi(4_{gen})", 32,-3.2,3.2, BINARY_SIZE, (found_4g ? physmath::deltaPhi(it_4g->daughter(0), it_4g->daughter(1)) : -9.), type);
	}
	
	if(found_8g && !found_4g){	//Special_174 AK8 that win because they're alone
		theHistograms.fill("Special_174 mass", "Special_174 mass", 16,50.,130., it_8g->mass());
		theHistograms.fill("Special_174 pt",   "Special_174 pt",  16,200.,1000., it_8g->pt());
		theHistograms.fill("Special_174 #Delta#phi(ZZ)_div_pi", "Special_174 #Delta#phi(ZZ)_div_pi",  20,-1.,1., physmath::deltaPhi(*ZZ, *it_8g)/M_PI );
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
			theHistograms.fill("bestZ AK8_{rec} mass", "bestZ AK8_{rec} mass", 30,0.,150, it_8r->mass(), 1.);
			theHistograms.fill<TH2F>("BestZ AK8_{rec} mass vs ZZ pt", "BestZ AK8_{rec} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), it_8r->corrPrunedMass(), 1.);
			theHistograms.fill<TH2F>("BestZ AK8_{rec} mass vs ZZ mass", "BestZ AK8_{rec} mass vs ZZ mass", PT_2D_SIZE, MASS_2D_SIZE, ZZ->mass(), it_8r->corrPrunedMass(), 1.);
			
			theHistograms.fill("Best Z (8_{rec}) tau21", "Best Z (8_{rec}) tau21", 20,0.,1., it_8r->tau2());
		}
	}
	
	//------------------	REC AK4	------------------
	//fillRecHadVBs();
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
			theHistograms.fill("bestZ AK4_{rec} mass", "bestZ AK4_{rec} mass", 30,0.,150, it_4r->mass(), 1.);
			theHistograms.fill<TH2F>("BestZ AK4_{rec} mass vs ZZ pt", "BestZ AK4_{rec} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), it_4r->mass(), 1.);
			theHistograms.fill<TH2F>("BestZ AK4_{rec} mass vs ZZ mass", "BestZ AK4_{rec} mass vs ZZ mass", PT_2D_SIZE, MASS_2D_SIZE, ZZ->mass(), it_4r->mass(), 1.);
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
		
		theHistograms.fill<TH2F>("Best Z rec vs ZZ pt", "Best Z rec vs ZZ pt", PT_2D_SIZE, BINARY_SIZE, ZZ->pt(), type);
		theHistograms.fill<TH2F>("Best Z rec vs s", "Best Z rec vs s", S_2D_SIZE, BINARY_SIZE, (type > 0.5 ? sAK8r_ : sAK4r_), type);
		
		theHistograms.fill<TH2F>("Best Z vs deltaR(4_{rec})", "Best Z vs deltaR(4_{rec})", 50,0.,5., BINARY_SIZE, (found_4r ? physmath::deltaR( it_4r->daughter(0), it_4r->daughter(1) ) : -1.), type /*filling underflow bin --> bestZ found with AK8*/ );
		theHistograms.fill<TH2F>("Best Z vs deltaPhi(4_{rec})","Best Z vs deltaPhi(4_{rec})", 32,-3.2,3.2, BINARY_SIZE, (found_4r ? physmath::deltaPhi(it_4r->daughter(0), it_4r->daughter(1)) : -9.), type);
	}
}


void VZZAnalyzer::minPtJetMVA(){
	if(!ZZ) return;
	if(ZZ->pt() < 1.) return;
	
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
			theHistograms.fill("#DeltaR(ZZ, minPt AK8_{gen})", "#DeltaR(ZZ, minPt AK8_{gen})", 35,0.,7., physmath::deltaR(*it, *ZZ), 1.);
			theHistograms.fill<TH2F>("Min PtTot AK8_{gen} mass vs ZZ pt", "Min PtTot AK8_{gen} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), it->mass(), 1.);
			
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
			theHistograms.fill("#DeltaR(ZZ, minPt AK4_{gen})", "#DeltaR(ZZ, minPt AK4_{gen})", 35,0.,7., physmath::deltaR(*it, *ZZ), 1.);
			theHistograms.fill<TH2F>("Min PtTot AK4_{gen} mass vs ZZ pt", "Min PtTot AK4_{gen} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), it->mass(), 1.);
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
			theHistograms.fill("#DeltaR(ZZ, minPt AK8_{rec})", "#DeltaR(ZZ, minPt AK8_{rec})", 35,0.,7., physmath::deltaR(*it, *ZZ), 1.);
			theHistograms.fill<TH2F>("Min PtTot AK8_{rec} mass vs ZZ pt", "Min PtTot AK8_{rec} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), it->corrPrunedMass(), 1.);
		}
	}
	//------------------	REC AK4	------------------
	//fillRecHadVBs();
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
			theHistograms.fill("#DeltaR(ZZ, minPt AK4_{rec})", "#DeltaR(ZZ, minPt AK4_{rec})", 35,0.,7., physmath::deltaR(*it, *ZZ), 1.);
			theHistograms.fill<TH2F>("Min PtTot AK4_{rec} mass vs ZZ pt", "Min PtTot AK4_{rec} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), it->mass(), 1.);
		}
	}
}


void VZZAnalyzer::furthestJetMVA(){
	if(!ZZ) return;
	if(ZZ->pt() < 1.) return;  // similar to Particle::IsValid()
	
	//------------------	GEN PARTICLES	------------------
	//For the AK4 part: genVB are constructed from pairs of AK4 (FOR NOW!) See SignalDefinitions
	//fillGenHadVBs();
	
	if(genHadVBs_->size() > 0){ //Looks like sometimes ZZ is a void DiBoson
		if(genHadVBs_->size() > 1)
			std::stable_sort(genHadVBs_->begin(), genHadVBs_->end(), DeltaRComparator(*ZZ));
		
		//Now the furthest from ZZ is the LAST of genHadVBs
		theHistograms.fill<TH2F>("Furthest AK4_{gen} mass vs ZZ pt", "Furthest AK4_{gen} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), genHadVBs_->back().mass(), 1.);
		theHistograms.fill<TH2F>("Furthest AK4_{gen} mass vs ZZ P",  "Furthest AK4_{gen} mass vs ZZ P",  P_2D_SIZE,  MASS_2D_SIZE, ZZ->p(),  genHadVBs_->back().mass(), 1.);
		theHistograms.fill("Max #DeltaR(ZZ, AK4_{gen})", "Max #DeltaR(ZZ, AK4_{gen})", 30,2.,8., physmath::deltaR(genJets->back(), *ZZ), 1.);
	}
	
	// AK8
	Particle* furthGenAK8 = furthestSing(genJetsAK8, *ZZ, 2., make_pair(50.,130.));
	if(furthGenAK8 != nullptr){
		theHistograms.fill<TH2F>("Furthest AK8_{gen} mass vs ZZ pt", "Furthest AK8_{gen} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), furthGenAK8->mass(), 1.);
		theHistograms.fill<TH2F>("Furthest AK8_{gen} mass vs ZZ P",  "Furthest AK8_{gen} mass vs ZZ P",  P_2D_SIZE,  MASS_2D_SIZE, ZZ->p(),  furthGenAK8->mass(), 1.);
		theHistograms.fill("Max #DeltaR(ZZ, AK8_{gen})", "Max #DeltaR(ZZ, AK8_{gen})", 30,2.,8., physmath::deltaR(*furthGenAK8, *ZZ), 1.);
		delete furthGenAK8;
	}
	
	//------------------	REC PARTICLES	------------------
	// AK8
	Jet* furthestAK8 = furthestSing(jetsAK8, *ZZ, 2., make_pair(50.,130.));
	if(furthestAK8 != nullptr){
		if(furthestAK8->corrPrunedMass() < 50.)
			cout<<"  Outside "<<furthestAK8->corrPrunedMass()<<'\n';
		theHistograms.fill<TH2F>("Furthest AK8_{rec} mass vs ZZ pt", "Furthest AK8_{rec} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), furthestAK8->corrPrunedMass(), 1.);
		theHistograms.fill<TH2F>("Furthest AK8_{rec} mass vs ZZ P",  "Furthest AK8_{rec} mass vs ZZ P",  P_2D_SIZE,  MASS_2D_SIZE, ZZ->p(),  furthestAK8->corrPrunedMass(), 1.);
		theHistograms.fill("Max #DeltaR(ZZ, AK8_{rec})", "Max #DeltaR(ZZ, AK8_{rec})", 30,2.,8., physmath::deltaR(*furthestAK8, *ZZ), 1.);
		delete furthestAK8;
	}
	
	// AK4
	Boson<Jet>* furthestAK4 = furthestPair(jets, *ZZ, 2./*min deltaR*/, make_pair(5.,200.)); 
	if(furthestAK4 != nullptr){
		theHistograms.fill<TH2F>("Furthest AK4_{rec} mass vs ZZ pt", "Furthest AK4_{rec} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), furthestAK4->mass(), 1.);
		theHistograms.fill<TH2F>("Furthest AK4_{rec} mass vs ZZ P",  "Furthest AK4_{rec} mass vs ZZ P", P_2D_SIZE,  MASS_2D_SIZE, ZZ->p(),  furthestAK4->mass(), 1.);
		theHistograms.fill("Max #DeltaR(ZZ, AK4_{rec})", "Max #DeltaR(ZZ, AK8_{rec})", 30,2.,8., physmath::deltaR(*furthestAK4, *ZZ), 1.);
		delete furthestAK4;
	}
}


void VZZAnalyzer::closestJetAnalisys(){
	//1: find the VBs with hadronic deacay:
	//fillGenHadVBs();
	theHistograms.fill("# of V --> jj", "# of V --> jj", 4,0,4, genHadVBs_->size());
	if(genHadVBs_->size() == 0) return; //No analysis can be done without a V-->JJ
	
	
	foreach(const Boson<Particle>& hadVB, *genHadVBs_){
		//2.a: find the closest genAK8
		const Particle* closestGenAK8 = closestSing(genJetsAK8, hadVB);
		if(closestGenAK8){
			float dR = physmath::deltaR(*closestGenAK8, hadVB);
			float dM = closestGenAK8->mass() - hadVB.mass(); //genParticles --> mass()
			theHistograms.fill("#DeltaM (AK8_{gen}, V)", "#DeltaM (AK8_{gen}, V)", 40,-20,20, dM,1.);
			theHistograms.fill("#DeltaR (AK8_{gen}, V)", "#DeltaR (AK8_{gen}, V)", 40,0.,0.4, dR,1.);
		}
	
		//2.b: find the closest pair of genAK4
		Boson<Particle>* closestGenAK4s = closestPair(genJets, hadVB);
		if(closestGenAK4s){
			float dR = physmath::deltaR(*closestGenAK4s, hadVB);
			float dM = closestGenAK4s->mass() - hadVB.mass(); //genParticles --> mass()
			theHistograms.fill("#DeltaM (AK4_{gen} pair, V)", "#DeltaM (AK4_{gen} pair, V)", 40,-20,20, dM, 1.);
			theHistograms.fill("#DeltaR (AK4_{gen} pair, V)", "#DeltaR (AK4_{gen} pair, V)", 40,0.,0.4, dR, 1.);
		}
		
		//3.a: find the closest recAK8
		const Jet* closestAK8 = closestSing(jetsAK8, hadVB);
		if(closestAK8){
			float dR = physmath::deltaR(*closestAK8, hadVB);
			float dM = closestAK8->corrPrunedMass() - hadVB.mass(); //Jet -->corrPrunedmass()
			theHistograms.fill("#DeltaM_{corr-prun} (AK8_{rec}, V)", "#DeltaM_{corr-prun} (AK8_{rec}, V)", 40,-20,20, dM, 1.);
			theHistograms.fill("#DeltaR (AK8_{rec}, V)", "#DeltaR (AK8_{rec}, V)", 40,0.,0.4, dR,1.);
		}
		
		//3.b: find the closest pair of recAK4
		Boson<Jet>* closestAK4s = closestPair(jets, hadVB);
		if(closestAK4s != nullptr){
			float dR = physmath::deltaR(*closestAK4s, hadVB);
			float dM = closestAK4s->mass() - hadVB.mass(); //Boson is Particle --> mass()
			theHistograms.fill("#DeltaM (AK4_{rec} pair, V)", "#DeltaM (AK4_{rec} pair, V)", 40,-20,20, dM, 1.);
			theHistograms.fill("#DeltaR (AK4_{rec} pair, V)", "#DeltaR (AK4_{rec} pair, V)", 40,0.,0.4, dR, 1.);
		}
		
		//4: How often AK8 are better? Are there some variables that discriminate?
		if(closestAK4s != nullptr && closestAK8 != nullptr){
			float dM8 = closestAK8->corrPrunedMass() - hadVB.mass();
			float dM4 = closestAK4s->mass() - hadVB.mass();
			if(fabs(dM8) < fabs(dM4)){
				++win8_;
				theHistograms.fill("V_pt AK8 wins",  "V_pt AK8 wins",  21,0.,630,  hadVB.pt(),  1.);
				theHistograms.fill("V_P AK8 wins",   "V_P AK8 wins",   20,0.,2000, hadVB.p(),   1.);
				theHistograms.fill("V_eta AK8 wins", "V_eta AK8 wins",21,-3.15,3.15,hadVB.eta(),1.);
			}
			else{
				++win4_;
				theHistograms.fill("V_pt AK4s win",  "V_pt AK4s win",  21,0.,630,  hadVB.pt(),  1.);
				theHistograms.fill("V_P AK4s win",   "V_P AK4s win",   20,0.,2000, hadVB.p(),   1.);
				theHistograms.fill("V_eta AK4s win", "V_eta AK4s win",21,-3.15,3.15,hadVB.eta(),1.);
			}
		}
		
		//Cleanup of allocated particles
		if(closestGenAK4s) delete closestGenAK4s;
		if(closestAK4s)    delete closestAK4s;
	}
}


template <class P, class R = Boson<Particle>> // P = Jet or Particle
P* VZZAnalyzer::closestSing(vector<P>* cands, const R& reference){
	if(cands->size() < 1) return nullptr;
	else{
		if(cands->size() > 1)
			std::sort( cands->begin(), cands->end(), phys::DeltaRComparator(reference) );
		if(physmath::deltaR(reference, cands->front()) < 0.4 )
			return &(cands->front());
		else return nullptr;
	}
}


template <class P, class R = Boson<Particle>> // P = Jet or Particle
Boson<P>* VZZAnalyzer::closestPair(vector<P>* cands, const R& reference){
	if(cands->size() < 2) return nullptr;
	else{
		if(cands->size() == 2){
			Boson<P>* res = new Boson<P>( cands->at(0), cands->at(1) );
			if(physmath::deltaR( *res, reference ) < 0.4 )
				return res;
			else return nullptr;
		}
		else{
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
	}
}


template <class P, class R = phys::DiBoson<phys::Lepton, phys::Lepton>> //P = Particle or Jet
P* VZZAnalyzer::furthestSing(vector<P>* cands, const R& reference, const float& minDR, const pair<float,float>& massLimits){
	if(cands->size() < 1)
		return nullptr;
	
	vector<P> goodCands;
	//cout<<"\n\nLooping... ";
	for(size_t i = 0; i < cands->size(); ++i){
		float thisMass = getRefinedMass(cands->at(i));
		//cout<<thisMass<<"   ";
		if(thisMass > massLimits.first && thisMass < massLimits.second){
			goodCands.push_back(cands->at(i));
		}
		//else cout<<" -discarding "<<thisMass<<' '<<(thisMass > massLimits.first)<<' '<<(thisMass < massLimits.second)<<"- ";
	}
	if(goodCands.size() > 0){	
		if(goodCands.size() > 1)
			std::sort( goodCands.begin(), goodCands.end(), phys::DeltaRComparator(reference) );
		if(physmath::deltaR(reference, goodCands.back()) > minDR ){
			P* result = new P(goodCands.back());
			//cout<<"\n\tResult: ";
			double resMass = getRefinedMass(*result);
			//cout<<resMass;
			return result;
		}
		else return nullptr;
		//cout<<"\n\tNo Result: deltaRMin = "<<physmath::deltaR(reference, goodCands.back());
	}
	else{
		//cout<<"\n\tNo Result: no good cands";
		return nullptr;
	}
}


template <class P, class R = phys::DiBoson<phys::Lepton, phys::Lepton>> // P = Jet or Particle
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
	const char* hNames[] = {"V_pt AK8 wins", "V_P AK8 wins", "V_eta AK8 wins", "V_pt AK4s win", "V_P AK4s win", "V_eta AK4s win"};
	foreach(const char* name, hNames){
		TH1* h = theHistograms.get(name);
		Double_t scale = 1./(h->Integral());
		h->Scale(scale);
	}
}


void VZZAnalyzer::ptCutMVA(){
	if(jetsAK8->size() == 0)
		return; //Give up: can't analyze mass functions for AK8s if there are none
	vector<float> vector_cuts;
	for(int i = 0; i < 10; ++i){
		vector_cuts.push_back((float)(170.+10*i)); //there's no variation before 170 GeV
	}
	const char* massesNames[6] = {"simpleMass", "secvtxMass", "corrPrunedMass", "prunedMass", "softDropMass", "puppiMass"}; //[1] = {"corrPrunedMass"};
	//vector<float (*)(const Jet&)> VectMassFunction(5);
	
	foreach(const Boson<Particle>& genVB, *genVBParticles){
		stable_sort(jetsAK8->begin(), jetsAK8->end(), phys::DeltaRComparator(genVB));
		//DeltaRComparator is defined in Commons/interface/Utils.h
		if(physmath::deltaR(jetsAK8->at(0), genVB) < 0.4 /*or 0.8?*/){
			//We have a match! This AK8 should have the mass of a VB (80-90 GeV)
			//Compute the mass once and for all
			float simpleMass     = jetsAK8->at(0).mass();
			float secvtxMass     = jetsAK8->at(0).secvtxMass();
			float corrPrunedMass = jetsAK8->at(0).corrPrunedMass();
			float prunedMass     = jetsAK8->at(0).prunedMass();
			float softDropMass   = jetsAK8->at(0).softDropMass();
			float puppiMass      = jetsAK8->at(0).puppiMass();
			float massesVal[6] = {simpleMass, secvtxMass, corrPrunedMass, prunedMass, softDropMass, puppiMass}; //[1] = {corrPrunedMass};
			
			foreach(const float& cut, vector_cuts){
				if(jetsAK8->at(0).pt() > cut){
					for(int i = 0; i < 6 /*1*/ /*massesVal's size*/; ++i){
						char name[32];
						sprintf(name, "%s: pt > %d", massesNames[i], (int)(cut));
						theHistograms.fill(name, name, 30,0.,120., massesVal[i], 1./*theWeight*/);
					}
				}
				//else break; //If it doesn't pass this pt selection, it won't pass the others
			}
		}
	}
}


void VZZAnalyzer::bestCandidateAnalysis(){
	//Find the best candidate in Jets and JetsAK8 as single fat Jet or as a pair (Boson<Jet>)
	phys::Boson<phys::Jet>* candPairJets = nullptr; //initialization
	VCandType typePairJets = VCandType::None; //initialization
	bool pairFromJets = findBestVFromPair(jets, candPairJets, typePairJets); //modifies VCandidate and candType
	if(pairFromJets){
		( typePairJets == VCandType::W ?  ++pairWFromJets_  : ++pairZFromJets_  );
		theHistograms.fill("pairFromJets_M","pairFromJets_M",CAND_M_SIZE, candPairJets->mass(), 1./*theWeight*/);
	}
	
	const phys::Jet* candSingJets = nullptr;
	VCandType typeSingJets = VCandType::None;
	bool singFromJets = findBestVFromSing(jets, candSingJets, typeSingJets);
	if(singFromJets){
		( typeSingJets == VCandType::W ?  ++singWFromJets_  : ++singZFromJets_  );
		theHistograms.fill("singFromJets_M","singFromJets_M",CAND_M_SIZE, candSingJets->mass(), 1./*theWeight*/);
	}
	
	phys::Boson<phys::Jet>* candPairJetsAK8 = nullptr;
	VCandType typePairJetsAK8 = VCandType::None;
	bool pairFromJetsAK8 = findBestVFromPair(jetsAK8, candPairJetsAK8, typePairJetsAK8);
	if(pairFromJetsAK8){
		(typePairJetsAK8==VCandType::W ? ++pairWFromJetsAK8_ : ++pairZFromJetsAK8_);
		theHistograms.fill("pairFromAK8_M","pairFromAK8_M",CAND_M_SIZE, candPairJetsAK8->mass(), 1./*theWeight*/);
	}
	
	const phys::Jet* candSingJetsAK8 = nullptr;
	VCandType typeSingJetsAK8 = VCandType::None;
	bool singFromJetsAK8 = findBestVFromSing(jetsAK8, candSingJetsAK8, typeSingJetsAK8);
	if(singFromJetsAK8){
		(typeSingJetsAK8==VCandType::W ? ++singWFromJetsAK8_ : ++singZFromJetsAK8_);
		theHistograms.fill("singFromAK8_M","singFromAK8_M",CAND_M_SIZE, candSingJetsAK8->mass(), 1./*theWeight*/);
	}
	
	//Choose best candidate
	vector<const phys::Particle*>* candidates = new vector<const Particle*>(); // Both Jet and Boson inherit from Particle
	if(pairFromJets)    candidates->push_back(candPairJets);   //0
	if(singFromJets)    candidates->push_back(candSingJets);   //1
	if(pairFromJetsAK8) candidates->push_back(candPairJetsAK8);//2
	if(singFromJetsAK8) candidates->push_back(candSingJetsAK8);//3
	
	const phys::Particle* bestV = nullptr;
	VCandType bestVtype = VCandType::None;
	bool b = findBestVPoint(candidates, bestV, bestVtype);
	
	// ----- ----- PLOTS! ----- -----
	if(b && ZZ != nullptr){
		++recVBtot_; //Total VB reconstructed
		//if(ZZ->pt() > 1.){
		TLorentzVector p4Cand = bestV->p4();
		TLorentzVector p4Tot = p4Cand + ZZ->p4();
		float ptTot = p4Tot.Pt();
		float mTot = p4Tot.M(); //Doesn't make much sense, but let's try anyway
		float deltaEtaTot = fabs(p4Cand.Eta() - ZZ->eta());
		float deltaRTot = ZZ->p4().DeltaR(p4Cand); //TLorentzVector::DeltaR()
		theHistograms.fill("bestV_Mass", "bestV_Mass;[GeV]", CAND_M_SIZE,bestV->mass(), 1./*theWeight*/);
		theHistograms.fill("ptTot_r_", "ptTot;[GeV]", 200,0.,400., ptTot, 1./*theWeight*/);
		theHistograms.fill("mTot", "mTot;[GeV]", 200,0.,2000., mTot, 1./*theWeight*/);
		theHistograms.fill("deltaEtaTot_r_", "#DeltaEtaTot", 100,0.,5., deltaEtaTot, 1./*theWeight*/);
		theHistograms.fill("deltaRTot_r_", "#DeltaRTot;[GeV]", 140,0.,7., deltaRTot, 1./*theWeight*/);

		float mWJJ_norm = fabs(p4Cand.M() - phys::WMASS)/phys::WMASS;
		theHistograms.fill("mWJJ_norm_r_","M_W-M_JJ_norm;[m_W]", 200,0.,.8, mWJJ_norm,1./*theWeight*/);
		theHistograms.fill("mWJJ_norm_sq_r_","(M_W-M_JJ_norm)^2;[m_W]", 300,0.,.6, mWJJ_norm*mWJJ_norm, 1./*theWeight*/);
		//}
	}
	
	//Let's find the best genV    TODO: define signal region(?)
	Boson<Particle>* genVB = nullptr;
	if(b && genVBParticles->size() >= 1){
		if(genVBParticles->size() >= 2){
			sort( genVBParticles->begin(), genVBParticles->end(), phys::DeltaRComparator(bestV->p4()) );
		}
		genVB = &(genVBParticles->at(0));
		++withGenVB_;
		theHistograms.fill("_N_genVB_pt",   "genVB_pt",   PT_SIZE,  genVB->pt(),  1.);
		theHistograms.fill("_N_genVB_E",    "genVB_E",    E_SIZE,   genVB->e(),   1.);
		theHistograms.fill("_N_genVB_#eta", "genVB_#eta", ETA_SIZE, genVB->eta(), 1.);
		theHistograms.fill("_N_genVB_M",    "genVB_M",    M_SIZE,   genVB->mass(),1.);
		
		//Match of reconstructed VB with the closest generated VB
		float deltaRmin = physmath::deltaR(bestV->p4(), genVB->p4());
		float deltaMclosest = fabs(bestV->mass() - genVB->mass());
		theHistograms.fill("#DeltaR_gen_rec_VB_n_", "#DeltaR_gen_rec_VB_n_", 120,0.,6., deltaRmin);
		theHistograms.fill("#DeltaM_gen_rec_VB_n_", "#DeltaM_gen_rec_VB_n_", 100,0.,50., deltaMclosest);
		if(deltaRmin < 0.4){//Filling with the corresponding generated VB properties, for efficiency analysis
			++goodRec_;
			theHistograms.fill("_N_recVB_pt",   "recVB_pt",   PT_SIZE,  genVB->pt(),  1.);
			theHistograms.fill("_N_recVB_E",    "recVB_E",    E_SIZE,   genVB->e(),   1.);
			theHistograms.fill("_N_recVB_#eta", "recVB_#eta", ETA_SIZE, genVB->eta(), 1.);
			theHistograms.fill("_N_recVB_M",    "recVB_M",    M_SIZE,   genVB->mass(),1.);
		}
	}
	
	
	if(candPairJets)    delete candPairJets;   //Cleaning the Bosons created in this event
	//if(candSingJets)    delete candSingJets; //It's an element of jet, no need to delete
	if(candPairJetsAK8) delete candPairJetsAK8;
	//if(candSingJetsAK8) delete candSingJetsAK8;
	if(candidates) delete candidates;
}


void VZZAnalyzer::endResolutionAnalisys(TFile& fout){
	TH2* genAK8_ZZpt = dynamic_cast<TH2*>(theHistograms.get("BestZ AK8_{gen} mass vs ZZ pt"));
	if(genAK8_ZZpt == nullptr) return;
	TH2* genAK4_ZZpt = dynamic_cast<TH2*>(theHistograms.get("BestZ AK4_{gen} mass vs ZZ pt"));
	
	fout.cd();
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
	TH2* genType_s = dynamic_cast<TH2*>(theHistograms.get("Best Z gen vs s"));
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
	TH2* genType_pt = dynamic_cast<TH2*>(theHistograms.get("Best Z gen vs ZZ pt"));
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


void VZZAnalyzer::endBestCandAnalysis(TFile & fout){
	TH1* genVB_pt  = theHistograms.get("_N_genVB_pt");
	TH1* genVB_E   = theHistograms.get("_N_genVB_E");
	TH1* genVB_eta = theHistograms.get("_N_genVB_#eta");
	TH1* genVB_M   = theHistograms.get("_N_genVB_M");
	TH1* recVB_pt  = theHistograms.get("_N_recVB_pt");
	TH1* recVB_E   = theHistograms.get("_N_recVB_E");
	TH1* recVB_eta = theHistograms.get("_N_recVB_#eta");
	TH1* recVB_M   = theHistograms.get("_N_recVB_M");
	
	fout.cd();
	if(genVB_pt != nullptr && recVB_pt != nullptr){
		int nbins = genVB_pt->GetNbinsX();
		genVB_pt->Fill(genVB_pt->GetBinCenter(nbins), genVB_pt->GetBinContent(nbins+1));
		genVB_pt->SetBinError(nbins, sqrt(genVB_pt->GetBinContent(nbins)));
		recVB_pt->Fill(recVB_pt->GetBinCenter(nbins), recVB_pt->GetBinContent(nbins+1));
		recVB_pt->SetBinError(nbins, sqrt(recVB_pt->GetBinContent(nbins)));
		TGraphAsymmErrors* ptEff  = new TGraphAsymmErrors(recVB_pt,  genVB_pt,  "cp");
		//ptEff->GetYaxis()->SetRangeUser(0.,1.01);
		ptEff->Write();
	}
	if(genVB_E != nullptr && recVB_E != nullptr){
		int nbins = genVB_E->GetNbinsX();
		genVB_E->Fill(genVB_E->GetBinCenter(nbins), genVB_E->GetBinContent(nbins+1));
		genVB_E->SetBinError(nbins, sqrt(genVB_E->GetBinContent(nbins)));
		recVB_E->Fill(recVB_E->GetBinCenter(nbins), recVB_E->GetBinContent(nbins+1));
		recVB_E->SetBinError(nbins, sqrt(recVB_E->GetBinContent(nbins)));
		TGraphAsymmErrors* EEff   = new TGraphAsymmErrors(recVB_E,   genVB_E,   "cp");
		//EEff->GetYaxis()->SetRangeUser(0.,1.01);
		EEff->Write();
	}
	if(genVB_eta != nullptr && recVB_eta != nullptr){
		int nbins = genVB_eta->GetNbinsX();
		genVB_eta->Fill(genVB_eta->GetBinCenter(nbins), genVB_eta->GetBinContent(nbins+1));
		genVB_eta->SetBinError(nbins, sqrt(genVB_eta->GetBinContent(nbins)));
		recVB_eta->Fill(recVB_eta->GetBinCenter(nbins), recVB_eta->GetBinContent(nbins+1));
		recVB_eta->SetBinError(nbins, sqrt(recVB_eta->GetBinContent(nbins)));
		TGraphAsymmErrors* etaEff = new TGraphAsymmErrors(recVB_eta, genVB_eta, "cp");
		//etaEff->GetYaxis()->SetRangeUser(0.,1.01);
		etaEff->Write();
	}
	if(genVB_M != nullptr && recVB_M != nullptr){
		int nbins = genVB_M->GetNbinsX();
		genVB_M->Fill(genVB_M->GetBinCenter(nbins), genVB_M->GetBinContent(nbins+1));
		genVB_M->SetBinError(nbins, sqrt(genVB_M->GetBinContent(nbins)));
		recVB_M->Fill(recVB_M->GetBinCenter(nbins), recVB_M->GetBinContent(nbins+1));
		recVB_M->SetBinError(nbins, sqrt(recVB_M->GetBinContent(nbins)));
		TGraphAsymmErrors* MEff = new TGraphAsymmErrors(recVB_M, genVB_M, "cp");
		//etaEff->GetYaxis()->SetRangeUser(0.,1.01);
		MEff->Write();
	}
}


void VZZAnalyzer::simpleGraphs(){
  if(ZZ != nullptr){
		theHistograms.fill("ZZmass","ZZ mass;[GeV/c^2]", 200,0.,500., ZZ->mass(), theWeight);
		theHistograms.fill("ZZpt","ZZ pt;[GeV/c]", 200,0.,500., ZZ->pt(), theWeight);
  }
	/*else
		cout<<"evtN_: ZZ == nullptr\n";*/
	theHistograms.fill("genVBParticles","genVBParticles->size()", 5,-0.5,4.5, genVBParticles->size(), 1./*theWeight*/);
	foreach(const phys::Boson<phys::Particle>& gen, *genVBParticles){
    if(abs(gen.id()) == 24 && abs(gen.daughter(0).id()) < 10) {
      theHistograms.fill("genMass_W", "mass of gen W", 100,35.5,135.5, gen.mass(), 1./*theWeight*/);
    }
    if(abs(gen.id()) == 23 && abs(gen.daughter(0).id()) < 10) {
      theHistograms.fill("genMass_Z", "mass of gen Z", 100,35.5,135.5, gen.mass(), 1./*theWeight*/);
    }
  }
	
	theHistograms.fill("AK4_{rec} size","AK4_{rec} mass", 8,-0.5,7.5, jets->size(), 1./*theWeight*/);
	theHistograms.fill("AK8_{rec} size","AK8_{rec} mass", 8,-0.5,7.5, jetsAK8->size(), 1./*theWeight*/);
	theHistograms.fill("AK4_{gen} size","AK4_{gen} mass", 8,-0.5,7.5, genJets->size(), 1./*theWeight*/);
	theHistograms.fill("AK8_{gen} size","AK8_{gen} mass", 8,-0.5,7.5, genJetsAK8->size(), 1./*theWeight*/);
	
	/*foreach(const Jet& jet, *jets)
		theHistograms.fill("All AK4_{rec} mass", "All AK4_{rec} mass", 25,-5.,120., jet.corrPrunedMass(), theWeight);*/
	
	foreach(const Jet& jet, *jetsAK8){
		theHistograms.fill("All AK8_{rec} mass", "All AK8_{rec} mass", 25,-5.,120., jet.corrPrunedMass(), 1./*theWeight*/);
		if(ZZ != nullptr)
			if(ZZ->pt() > 1.)
				theHistograms.fill<TH2F>("All AK8_{rec} mass vs ZZ pt", "All AK8_{rec} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), jet.mass(), 1.);
	}
	
	/*foreach(const Particle& jet, *genJets)
		theHistograms.fill("All AK4_{gen} mass", "All AK4_{gen} mass", 25,-5.,120., jet.mass(), theWeight);*/
	
	foreach(const Particle& jet, *genJetsAK8){
		theHistograms.fill("All AK8_{gen} mass", "All AK8_{gen} mass", 25,-5.,120., jet.mass(), 1./*theWeight*/);
		if(ZZ != nullptr)
			if(ZZ->pt() > 1.)
				theHistograms.fill<TH2F>("All AK8_{gen} mass vs ZZ pt", "All AK8_{gen} mass vs ZZ pt", PT_2D_SIZE, MASS_2D_SIZE, ZZ->pt(), jet.mass(), 1./*theWeight*/);
	}
}


template<class P = phys::Jet>
bool VZZAnalyzer::findBestVFromSing(/*const*/ std::vector<P>* js, const P*& thisCandidate, VCandType& thisCandType){
	bool isAccurate = false;
	thisCandType = VCandType::None;
	if(js->size() < 1)
		return false;
	size_t indexZ = 0;
	size_t indexW = 0;
	float minDifZ = 1.;
	float minDifW = 1.;
	float tmpMass = 0.;
	for(size_t i = 0; i < js->size(); ++i){
		tmpMass = js->at(i).mass();
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
	
	if(minDifZ < minDifW){
		thisCandidate = &(js->at(indexZ));
		isAccurate = ZBosonDefinition(*thisCandidate);
		thisCandType = VCandType::Z;
	}
	else{
		thisCandidate = &(js->at(indexW));
		isAccurate = WBosonDefinition(*thisCandidate);
		thisCandType = VCandType::W;
	}
	return isAccurate;
}


template <class J = phys::Jet>
bool VZZAnalyzer::findBestVFromPair(const std::vector<J>* js, phys::Boson<J>*& thisCandidate, VCandType& thisCandType){
	bool isAccurate = false;
	thisCandType = VCandType::None;
	if(js->size() < 2)
		return false;
		
	pair<size_t, size_t> indicesZ(0,0);
	pair<size_t, size_t> indicesW(0,0);
	float minDifZ = 1.;
	float minDifW = 1.;
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
		
		if(minDifZ < minDifW){
			thisCandidate = new Boson<J>(js->at(indicesZ.first), js->at(indicesZ.second));
			isAccurate = ZBosonDefinition(*thisCandidate);
			thisCandType = VCandType::Z;
		}
		else{
			thisCandidate = new Boson<J>(js->at(indicesW.first), js->at(indicesW.second));
			isAccurate = WBosonDefinition(*thisCandidate);
			thisCandType = VCandType::W;
		}
	}
	return isAccurate;
}

template <class P = phys::Particle>
bool VZZAnalyzer::findBestVPoint(std::vector<const P*>* js, const P*& thisCandidate, VCandType& thisCandType){
	bool isAccurate = false;
	thisCandType = VCandType::None;
	if(js->size() == 0)
		return false;
	size_t indexZ = 0;
	size_t indexW = 0;
	float minDifZ = 1.;
	float minDifW = 1.;
	float tmpMass = 0.;
	for(size_t i = 0; i < js->size(); ++i){
		tmpMass = js->at(i)->mass();
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
	if(minDifZ < minDifW){
		thisCandidate = (js->at(indexZ));
		isAccurate = ZBosonDefinition(*thisCandidate); //calling the version for a Particle, not the one for a Boson<Particle>
		thisCandType = VCandType::Z;
	}
	else{
		thisCandidate = (js->at(indexW));
		isAccurate = WBosonDefinition(*thisCandidate);
		thisCandType = VCandType::W;
	}
	return isAccurate;
}

/*
void VZZAnalyzer::fillGenHadVBs(){
	if(genHadVBs_ == nullptr)
		genHadVBs_ = new vector<Boson<Particle>>;
	
	if(genHadVBs_->size() > 0)
		return;
	else{
		foreach(const Boson<Particle>& genVB, *genVBParticles){
			if( abs(genVB.daughter(0).id()) < 10 && abs(genVB.daughter(1).id()) < 10 )
				genHadVBs_->push_back(genVB); //The second condition is redundant
		}
	}
}*/
void VZZAnalyzer::fillGenHadVBs(){
	foreach(const Boson<Particle>& genVB, *genVBParticles)
		if( abs(genVB.daughter(0).id()) < 10 && abs(genVB.daughter(1).id()) < 10 )
			genHadVBs_->push_back(genVB); //The second condition is redundant
}

/*
void VZZAnalyzer::fillRecHadVBs(){
	if(AK4pairs_ == nullptr)
		AK4pairs_ = new vector<Boson<Jet>>;
	
	if(AK4pairs_->size() > 0)
		return;
	else
		for(size_t i = 0; i < jets->size(); ++i)
			for(size_t j = i+1; j < jets->size(); ++j)
				AK4pairs_->push_back( Boson<Jet>(jets->at(i), jets->at(j)) );
}*/
void VZZAnalyzer::fillRecHadVBs(){
	for(size_t i = 0; i < jets->size(); ++i)
		for(size_t j = i+1; j < jets->size(); ++j)
			AK4pairs_->push_back( Boson<Jet>(jets->at(i), jets->at(j)) );
}



void VZZAnalyzer::calcS(){
	TLorentzVector tot4g(ZZ->p4());
	foreach(const Particle& j, *genJets)
		tot4g += j.p4();
	sAK4g_ = tot4g.M();
	
	TLorentzVector tot8g(ZZ->p4());
	foreach(const Particle& j, *genJetsAK8)
		tot8g += j.p4();
	sAK8g_ = tot8g.M();
	
	TLorentzVector tot4r(ZZ->p4());
	foreach(const Jet& j, *jets)
		tot4r += j.p4();
	sAK4r_ = tot4r.M();
	
	TLorentzVector tot8r(ZZ->p4());
	foreach(const Jet& j, *jetsAK8)
		tot4g += j.p4();
	sAK8r_ = tot8r.M();
}




















