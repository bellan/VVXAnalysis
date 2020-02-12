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

#include "boost/foreach.hpp"
#define foreach BOOST_FOREACH
#include "boost/assign/std/vector.hpp" 
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
	phys::Boson<phys::Jet> candPairJets; //initialization
	VCandType typePairJets = VCandType::None; //initialization
	bool pairFromJets = findBestVFromPair(jets, candPairJets, typePairJets); //modifies VCandidate and candType
	
	phys::Jet candSingJets;
	VCandType typeSingJets = VCandType::None;
	bool singFromJets = findBestVFromSing(jets, candSingJets, typeSingJets);
	
	phys::Boson<phys::Jet> candPairJetsAK8;
	VCandType typePairJetsAK8 = VCandType::None;
	bool pairFromJetsAK8 = findBestVFromPair(jetsAK8, candPairJetsAK8, typePairJetsAK8);
	
	phys::Jet candSingJetsAK8;
	VCandType typeSingJetsAK8 = VCandType::None;
	bool singFromJetsAK8 = findBestVFromSing(jetsAK8, candSingJetsAK8, typeSingJetsAK8);
	
	//TODO: chose best candidate
	if(pairFromJets && ZZ != nullptr){		
		TLorentzVector p4JJ = candPairJets.p4(); //TODO: candPairJets is temporary
		TLorentzVector p4Tot = p4JJ + ZZ->p4();
		float ptTot = p4Tot.Pt();
		float mTot = p4Tot.M(); //Doesn't make much sense, but let's try anyway
		float deltaEtaTot = fabs(p4JJ.Eta() - ZZ->eta());
		float deltaRTot = ZZ->p4().DeltaR(p4JJ); //TLorentzVector::DeltaR()
		theHistograms.fill("ptTot_r_", "ptTot;[GeV]", 200,0.,400., ptTot, theWeight);
		theHistograms.fill("mTot", "mTot;[GeV]", 200,0.,2000., mTot, theWeight);
		theHistograms.fill("deltaEtaTot_r_", "#DeltaEtaTot", 100,0.,5., deltaEtaTot, theWeight);
		theHistograms.fill("deltaRTot_r_", "#DeltaRTot;[GeV]", 140,0.,7., deltaRTot, theWeight);
		
		float mWJJ_norm = fabs(p4JJ.M() - phys::WMASS)/phys::WMASS;
		theHistograms.fill("mWJJ_norm_r_","M_W-M_JJ_norm;[m_W]", 200,0.,.8, mWJJ_norm, theWeight);
		theHistograms.fill("mWJJ_norm_sq_r_","(M_W-M_JJ_norm)^2;[m_W]", 300,0.,.6, mWJJ_norm*mWJJ_norm, theWeight);
	}
	
	simpleGraphs();
	
	return;
}

void VZZAnalyzer::end(TFile &){
	float elapsedSec = (float)(clock()-startTime_)/CLOCKS_PER_SEC;
	int elapsedSecInt = (int)elapsedSec;
	cout<<"\nElapsed Time: "<<elapsedSec<<" s\t\t("<<elapsedSecInt/60<<"\' "<<elapsedSecInt%60<<"\")\n";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<" \tEnd of VZZAnalyzer\t";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<"\n\n";
}


void VZZAnalyzer::simpleGraphs(){
	if(ZZ != nullptr)
		theHistograms.fill("ZZmass","ZZ mass;[GeV/c^2]", 200,0.,500., ZZ->mass(), theWeight);
	/*else
		cout<<"evtN_: ZZ == nullptr\n";*/
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
}

bool VZZAnalyzer::findBestVFromSing(const std::vector<phys::Jet>* js, phys::Jet& thisCandidate, VCandType& thisCandType){
	bool isAccurate = false;
	thisCandType = VCandType::None;
	if(js->size() < 1)
		return false;
	size_t index = 0;
	float minDifZ = 1.;
	float minDifW = 1.;
	float massZCand = 0.;
	float massWCand = 0.;
	float tmpMass = 0.;
	for(size_t i = 0; i < js->size(); i++){
		tmpMass = js->at(i).mass();
		float diffZ = (tmpMass - phys::ZMASS)/phys::ZMASS;
		float diffZa = fabs(diffZ);
		float diffW = (tmpMass - phys::WMASS)/phys::WMASS;
		float diffWa = fabs(diffW);
		if(diffZa < minDifZ){
			minDifZ = diffZa;
			massZCand = tmpMass;
			index = i;
		}
		if(diffWa < minDifW){
			minDifW = diffWa;
			massWCand = tmpMass;
			index = i;
		}
	}
	thisCandidate = js->at(index);
	if(minDifZ < minDifW){
			isAccurate = true; //TODO: define a good candidate in case of fat jet
			thisCandType = VCandType::Z;
	}
	else{
		isAccurate = true;
		thisCandType = VCandType::W;
	}
	return isAccurate;
}

bool VZZAnalyzer::findBestVFromPair(const std::vector<phys::Jet>* js, phys::Boson<phys::Jet>& thisCandidate, VCandType& thisCandType){
	bool isAccurate = false;
	thisCandType = VCandType::None;
	if(js->size() < 2)
		return false;
		
	pair<size_t, size_t> indices(0,0);
	float minDifZ = 1.;
	float minDifW = 1.;
	float massZCand = 0.;
	float massWCand = 0.;
	float tmpMass = 0.;
	for(size_t i = 0; i < js->size(); i++){
		for(size_t j = i+1; j < js->size(); j++){
			tmpMass = (js->at(i).p4() + js->at(j).p4()).M();
			float diffZ = (tmpMass - phys::ZMASS)/phys::ZMASS;
			float diffZa = fabs(diffZ);
			float diffW = (tmpMass - phys::WMASS)/phys::WMASS;
			float diffWa = fabs(diffW);
			if(diffZa < minDifZ){
				minDifZ = diffZa;
				massZCand = tmpMass;
				indices = std::make_pair(i,j);
			}
			if(diffWa < minDifW){
				minDifW = diffWa;
				massWCand = tmpMass;
				indices = std::make_pair(i,j);
			}
		}
		thisCandidate = Boson<Jet>(js->at(indices.first), js->at(indices.second));
		
		if(minDifZ < minDifW){
			isAccurate = ZBosonDefinition(thisCandidate);
			thisCandType = VCandType::Z;
		}
		else{
			isAccurate = WBosonDefinition(thisCandidate);
			thisCandType = VCandType::W;
		}
	}
	return isAccurate;
}


