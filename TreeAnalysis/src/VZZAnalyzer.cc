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
#include "TGraphAsymmErrors.h"
#include "TString.h"

#include "boost/foreach.hpp"
#define foreach BOOST_FOREACH
#include "boost/assign/std/vector.hpp" 

#define PT_SIZE 40,0.,400.
#define E_SIZE 50,0.,1500.
#define ETA_SIZE 51,-5.05,5.05
#define M_SIZE 35,50.,120.
#define CAND_M_SIZE 40,50.,150.

using namespace boost::assign;

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;
using std::pair;

using namespace phys;

Int_t VZZAnalyzer::cut(){
	evtN_++;
	cout<<"\r\t\t"<<evtN_;
	return 1;
}

void VZZAnalyzer::begin(){
	cout<<"\n";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<" \tStart of VZZAnalyzer\t";
	for(char i=0; i<25; i++) cout<<"-";
	
	startTime_ = clock();
	
	cout<<"\nAnalyzed:\t      /"<<tree()->GetEntries()<<std::flush;
	return;
}

void VZZAnalyzer::analyze(){
	
	//bestCandidateAnalisys();
	//ptCutMVA();
	jetAnalisys();
	
	simpleGraphs();
	
	return;
}

void VZZAnalyzer::end(TFile & fout){
	cout<<'\n';
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
	
	
	float elapsedSec = (float)(clock()-startTime_)/CLOCKS_PER_SEC;
	int elapsedSecInt = (int)elapsedSec;
	cout<<"Candidate W: JetPair = "<<pairWFromJets_<<" \tJetSing = "<<singWFromJets_<<" \tAK8Pair = "<<pairWFromJetsAK8_<<" \tAK8Sing = "<<singWFromJetsAK8_;
	cout<<"\nCandidate Z: JetPair = "<<pairZFromJets_<<" \tJetSing = "<<singZFromJets_<<" \tAK8Pair = "<<pairZFromJetsAK8_<<" \tAK8Sing = "<<singZFromJetsAK8_;
	cout<<"\nEvents with a reconstructed VB: "<<recVBtot_<<" -> "<<recVBtot_*100./evtN_<<'%';
	cout<<"\nReconstructed near (deltaR < 0.4) a genVB: "<<goodRec_<<" \tevents with a genVB: "<<withGenVB_<<" \t--> "<<goodRec_*100./withGenVB_<<'%';
	cout<<"\nElapsed Time: "<<elapsedSec<<" s\t\t("<<elapsedSecInt/60<<"\' "<<elapsedSecInt%60<<"\")\n";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<" \tEnd of VZZAnalyzer\t";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<"\n\n";
}


void VZZAnalyzer::jetAnalisys(){
	//First: find the VB with hadronic deacay:
	vector<Boson<Particle>> genHadVBs;
	foreach(const Boson<Particle>& genVB, *genVBParticles){
		if( abs(genVB.daughter(0).id()) < 10 && abs(genVB.daughter(1).id()) < 10 )
			genHadVBs.push_back(genVB); //The second condition is redundant
	}
	if(genHadVBs.size() == 0) return; //No analisys can be done without a V-->JJ
	
	foreach(const Boson<Particle>& hadVB, genHadVBs){
		//Second.a: find the closest genAK8
		const Particle* closestGenAK8 = closestSing(genJetsAK8, hadVB);
		if(closestGenAK8){
			float dR = physmath::deltaR(*closestGenAK8, hadVB);
			if( dR < 0.4 ){
				float diffMass = hadVB.mass() - closestGenAK8->mass(); //genParticles --> mass()
				theHistograms.fill("#DeltaM (V_{gen}, AK8_{gen})", "#DeltaM (V_{gen}, AK8_{gen})", 40,-20,20, diffMass, 1.);
				theHistograms.fill("#DeltaR (V_{gen}, AK8_{gen})", "#DeltaR (V_{gen}, AK8_{gen})", 40,0.,0.4, dR, 1.);
			}
		}
	
		//Second.b: find the closest pair of genAK4
		Boson<Particle>* closestGenAK4s = closestPair(genJets, hadVB);
		if(closestGenAK4s){
			float dR = physmath::deltaR(*closestGenAK4s, hadVB);
			if( dR < 0.4 ){
				float diffMass = hadVB.mass() - closestGenAK4s->mass(); //genParticles --> mass()
				theHistograms.fill("#DeltaM (V_{gen}, AK4_{rec} pair)", "#DeltaM (V_{gen}, AK4_{rec} pair)", 40,-20,20, diffMass, 1.);
				theHistograms.fill("#DeltaR (V_{gen}, AK4_{rec} pair)", "#DeltaR (V_{gen}, AK4_{rec} pair)", 40,0.,0.4, dR, 1.);
			}
		}
		
		//Third.a: find the closest recAK8
		const Jet* closestAK8 = closestSing(jetsAK8, hadVB);
		if(closestAK8){
			float dR = physmath::deltaR(*closestAK8, hadVB);
			if( dR < 0.4 ){
				float diffMass = hadVB.mass() - closestAK8->corrPrunedMass(); //Jet -->corrPrunedmass()
				theHistograms.fill("#DeltaM_{corr-prun} (V_{gen}, AK8_{rec})", "#DeltaM_{corr-prun} (V_{gen}, AK8_{rec})", 40,-20,20, diffMass, 1.);
				theHistograms.fill("#DeltaR (V_{gen}, AK8_{rec})", "#DeltaR (V_{gen}, AK8_{rec})", 40,0.,0.4, dR, 1.);
			}
		}
		
		//Third.b: find the closest pair of recAK4
		Boson<Jet>* closestAK4s = closestPair(jets, hadVB);
		if(closestAK4s){
			float dR = physmath::deltaR(*closestAK4s, hadVB);
			if( dR < 0.4 ){
				float diffMass = hadVB.mass() - closestAK4s->mass(); //Boson is Particle --> mass()
				theHistograms.fill("#DeltaM (V_{gen}, AK4_{gen} pair)", "#DeltaM (V_{gen}, AK4_{gen} pair)", 40,-20,20, diffMass, 1.);
				theHistograms.fill("#DeltaR (V_{gen}, AK4_{gen} pair)", "#DeltaR (V_{gen}, AK4_{gen} pair)", 40,0.,0.4, dR, 1.);
			}
		}
		
		if(closestGenAK4s) delete closestGenAK4s; //Cleanup of allocated particles
		if(closestAK4s)    delete closestAK4s;
	}
	
	//End: How often AK8 are better? Are there some variables that discriminate?
}


template <class P, class R = Boson<Particle>> // P = Jet or Particle
P* VZZAnalyzer::closestSing(vector<P>* cands, const R& reference){
	if(cands->size() < 1) return nullptr;
	else{
		if(cands->size() > 1)
			std::sort( cands->begin(), cands->end(), phys::DeltaRComparator(reference) );
		return &(cands->at(0));
	}
}

template <class P, class R = Boson<Particle>> // P = Jet or Particle
Boson<P>* VZZAnalyzer::closestPair(vector<P>* cands, const R& reference){
	if(cands->size() < 2) return nullptr;
	else{
		if(cands->size() == 2)
			return new Boson<P>( cands->at(0), cands->at(1) );
		else{
			//Find the pair with the closest p4
			pair<size_t, size_t> indices(0,0);
			float minDR = 0.4; //starting value = the threshold we use
			for(size_t i = 0; i < cands->size(); i++)
				for(size_t j = i+1; j < cands->size(); j++){
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

void VZZAnalyzer::ptCutMVA(){
	if(jetsAK8->size() == 0)
		return; //Give up: can't analize mass functions for AK8s if there are none
	vector<float> vector_cuts;
	for(int i = 0; i < 10; i++){
		vector_cuts.push_back((float)(170.+10*i)); //there's no variation before 170 GeV
	}
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
			float massesVal[1] = {corrPrunedMass};//[6] = {simpleMass, secvtxMass, corrPrunedMass, prunedMass, softDropMass, puppiMass};
			const char* massesNames[1] = {"corrPrunedMass"};//[6] = {"simpleMass", "secvtxMass", "corrPrunedMass", "prunedMass", "softDropMass", "puppiMass"};
			
			foreach(const float& cut, vector_cuts){
				if(jetsAK8->at(0).pt() > cut){
					for(int i = 0; i<1/*6*/ /*massesVal's size*/; i++){
						char name[32];
						sprintf(name, "%s: pt > %d", massesNames[i], (int)(cut));
						theHistograms.fill(name, name, 30,0.,120., massesVal[i], theWeight);
					}
				}
			}
		}
	}
}


void VZZAnalyzer::bestCandidateAnalisys(){
	//Find the best candidate in Jets and JetsAK8 as single fat Jet or as a pair (Boson<Jet>)
	phys::Boson<phys::Jet>* candPairJets = nullptr; //initialization
	VCandType typePairJets = VCandType::None; //initialization
	bool pairFromJets = findBestVFromPair(jets, candPairJets, typePairJets); //modifies VCandidate and candType
	if(pairFromJets){
		( typePairJets == VCandType::W ?  pairWFromJets_++  : pairZFromJets_++  );
		theHistograms.fill("pairFromJets_M","pairFromJets_M",CAND_M_SIZE, candPairJets->mass(), theWeight);
	}
	
	const phys::Jet* candSingJets = nullptr;
	VCandType typeSingJets = VCandType::None;
	bool singFromJets = findBestVFromSing(jets, candSingJets, typeSingJets);
	if(singFromJets){
		( typeSingJets == VCandType::W ?  singWFromJets_++  : singZFromJets_++  );
		theHistograms.fill("singFromJets_M","singFromJets_M",CAND_M_SIZE, candSingJets->mass(), theWeight);
	}
	
	phys::Boson<phys::Jet>* candPairJetsAK8 = nullptr;
	VCandType typePairJetsAK8 = VCandType::None;
	bool pairFromJetsAK8 = findBestVFromPair(jetsAK8, candPairJetsAK8, typePairJetsAK8);
	if(pairFromJetsAK8){
		(typePairJetsAK8==VCandType::W ? pairWFromJetsAK8_++:pairZFromJetsAK8_++);
		theHistograms.fill("pairFromAK8_M","pairFromAK8_M",CAND_M_SIZE, candPairJetsAK8->mass(), theWeight);
	}
	
	const phys::Jet* candSingJetsAK8 = nullptr;
	VCandType typeSingJetsAK8 = VCandType::None;
	bool singFromJetsAK8 = findBestVFromSing(jetsAK8, candSingJetsAK8, typeSingJetsAK8);
	if(singFromJetsAK8){
		(typeSingJetsAK8==VCandType::W ? singWFromJetsAK8_++:singZFromJetsAK8_++);
		theHistograms.fill("singFromAK8_M","singFromAK8_M",CAND_M_SIZE, candSingJetsAK8->mass(), theWeight);
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
	
	/*if((pairFromJets || singFromJets || pairFromJetsAK8 || singFromJetsAK8) && !b){
		cout<<"pairJets "<<pairFromJets<<" \tsingFromJets "<<singFromJets<<" \tpairFromJetsAK8 "<<pairFromJetsAK8<<" \tsingFromJetsAK8 "<<singFromJetsAK8<<'\n';
		bool pJZ = typePairJets == VCandType::Z;
		float dMJ = fabs(candPairJets->mass() - (pJZ ? phys::ZMASS : phys::WMASS));
		cout<<"candPairJets "<<(pJZ ? 'Z' : 'W')<<' '<<dMJ<<'\n';
	}*/
	
	
	// ----- ----- PLOTS! ----- -----
	if(b && ZZ != nullptr){
		recVBtot_++; //Total VB reconstructed
		//if(ZZ->mass() > 1.){
		TLorentzVector p4Cand = bestV->p4();
		TLorentzVector p4Tot = p4Cand + ZZ->p4();
		float ptTot = p4Tot.Pt();
		float mTot = p4Tot.M(); //Doesn't make much sense, but let's try anyway
		float deltaEtaTot = fabs(p4Cand.Eta() - ZZ->eta());
		float deltaRTot = ZZ->p4().DeltaR(p4Cand); //TLorentzVector::DeltaR()
		theHistograms.fill("bestV_Mass", "bestV_Mass;[GeV]", CAND_M_SIZE,bestV->mass(), theWeight);
		theHistograms.fill("ptTot_r_", "ptTot;[GeV]", 200,0.,400., ptTot, theWeight);
		theHistograms.fill("mTot", "mTot;[GeV]", 200,0.,2000., mTot, theWeight);
		theHistograms.fill("deltaEtaTot_r_", "#DeltaEtaTot", 100,0.,5., deltaEtaTot, theWeight);
		theHistograms.fill("deltaRTot_r_", "#DeltaRTot;[GeV]", 140,0.,7., deltaRTot, theWeight);

		float mWJJ_norm = fabs(p4Cand.M() - phys::WMASS)/phys::WMASS;
		theHistograms.fill("mWJJ_norm_r_","M_W-M_JJ_norm;[m_W]", 200,0.,.8, mWJJ_norm,theWeight);
		theHistograms.fill("mWJJ_norm_sq_r_","(M_W-M_JJ_norm)^2;[m_W]", 300,0.,.6, mWJJ_norm*mWJJ_norm, theWeight);
		//}
	}
	
	//Let's find the best genV    TODO: define signal region(?)
	Boson<Particle>* genVB = nullptr;
	if(b && genVBParticles->size() >= 1){
		if(genVBParticles->size() >= 2){
			sort( genVBParticles->begin(), genVBParticles->end(), phys::DeltaRComparator(bestV->p4()) );
		}
		genVB = &(genVBParticles->at(0));
		withGenVB_++;
		theHistograms.fill("_N_genVB_pt",   "genVB_pt",   PT_SIZE,  genVB->pt(),  1.);
		theHistograms.fill("_N_genVB_E",    "genVB_E",    E_SIZE,   genVB->e(),   1.);
		theHistograms.fill("_N_genVB_#eta", "genVB_#eta", ETA_SIZE, genVB->eta(), 1.);
		theHistograms.fill("_N_genVB_M",    "genVB_M",    M_SIZE,   genVB->mass(),1.);
		
		//Match of reconstructed VB with the closest generated VB
		float deltaRmin = physmath::deltaR(bestV->p4(), genVB->p4());
		float deltaMclosest = fabs(bestV->mass() - genVB->mass());
		theHistograms.fill("#DeltaR_gen_rec_VB_n_", "#DeltaR_gen_rec_VB_n_", 120,0.,6., deltaRmin);
		theHistograms.fill("#DeltaM_gen_rec_VB_n_", "#DeltaM_gen_rec_VB_n_", 100,0.,50., deltaMclosest);
		if(deltaRmin < 0.4){//Filling with the corresponding generated VB properties, for efficiency analisys
			goodRec_++;
			theHistograms.fill("_N_recVB_pt",   "recVB_pt",   PT_SIZE,  genVB->pt(),  1.);
			theHistograms.fill("_N_recVB_E",    "recVB_E",    E_SIZE,   genVB->e(),   1.);
			theHistograms.fill("_N_recVB_#eta", "recVB_#eta", ETA_SIZE, genVB->eta(), 1.);
			theHistograms.fill("_N_recVB_M",    "recVB_M",    M_SIZE,   genVB->mass(),1.);
		}
	}
	
	
	if(candPairJets)    delete candPairJets;   //Cleaning the Bosons created in this event
	//if(candSingJets)    delete candSingJets; //It's an element of jet, so it is const
	if(candPairJetsAK8) delete candPairJetsAK8;
	//if(candSingJetsAK8) delete candSingJetsAK8;
	if(candidates) delete candidates;
}

void VZZAnalyzer::simpleGraphs(){
	if(ZZ != nullptr)
		theHistograms.fill("ZZmass","ZZ mass;[GeV/c^2]", 200,0.,500., ZZ->mass(), theWeight);
	/*else
		cout<<"evtN_: ZZ == nullptr\n";*/
	theHistograms.fill("genVBParticles","genVBParticles->size()", 5,-0.5,4.5, genVBParticles->size(), theWeight);
	foreach(const phys::Boson<phys::Particle>& gen, *genVBParticles){
    if(abs(gen.id()) == 24 && abs(gen.daughter(0).id()) < 10) {
      theHistograms.fill("genMass_W", "mass of gen W", 100,35.5,135.5, gen.mass(), theWeight);
    }
    if(abs(gen.id()) == 23 && abs(gen.daughter(0).id()) < 10) {
      theHistograms.fill("genMass_Z", "mass of gen Z", 100,35.5,135.5, gen.mass(), theWeight);
    }
  }
	/*
	if(jets->size()>=1){
		theHistograms.fill("jet1_pt","jet1_pt;[GeV/c]", 200,0.,200., jets->at(0).pt(), theWeight);
		theHistograms.fill("jet1_E","jet1_E;[GeV]", 200,0.,400., jets->at(0).e(), theWeight);
		if(jets->size()>=2){
			theHistograms.fill("jet2_pt","jet2_pt;[GeV/c]", 200,0.,200., jets->at(1).pt(), theWeight);
			theHistograms.fill("jet2_E","jet2_E;[GeV]", 200,0.,400., jets->at(1).e(), theWeight);
			theHistograms.fill("jets_deltaR","jets_deltaR", 200,0.,7., physmath::deltaR(jets->at(0),jets->at(1)), theWeight);
			theHistograms.fill("jets12_M","jets12_M;;[GeV/c^2]", 200,0.,400., (jets->at(0).p4()+jets->at(1).p4()).M(), theWeight);
		} else{
			theHistograms.fill("jet2_pt","jet2_pt;[GeV/c]", 200,0.,200., 0., theWeight);
			theHistograms.fill("jet2_E","jet2_E;[GeV]", 200,0.,400., 0., theWeight);
			theHistograms.fill("jets12_M","jets12_M;;[GeV/c^2]", 200,0.,400., 0., theWeight);
		}
	} else{
		theHistograms.fill("jet1_pt","jet1_pt;[GeV/c]", 200,0.,200., 0., theWeight);
		theHistograms.fill("jet1_E","jet1_E;[GeV]", 200,0.,400., 0., theWeight);
		theHistograms.fill("jets12_M","jets12_M;;[GeV/c^2]", 200,0.,400., 0., theWeight);
	}
	*/
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
	for(size_t i = 0; i < js->size(); i++){
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
	for(size_t i = 0; i < js->size(); i++){
		for(size_t j = i+1; j < js->size(); j++){
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
	for(size_t i = 0; i < js->size(); i++){
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
