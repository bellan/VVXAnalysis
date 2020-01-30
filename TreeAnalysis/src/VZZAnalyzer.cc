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
#include "TSystem.h"
#include "TTree.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;

using namespace phys;

Int_t VZZAnalyzer::cut(){
	evtN++;
	cout<<"\r\t\t"<<evtN;
	return 1;
}

void VZZAnalyzer::begin(){
	cout<<"\n";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<" \tStart of VZZAnalyzer\t";
	for(char i=0; i<25; i++) cout<<"-";
	
	startTime = clock();
	
	cout<<"\nAnalyzed:\t      /"<<tree()->GetEntries()<<std::flush;
	return;
}

void VZZAnalyzer::analyze(){
	if(ZZ != nullptr)
		theHistograms.fill("ZZmass","ZZ mass;[GeV/c^2]", 200,0.,500., ZZ->mass(), theWeight);
	/*else
		cout<<"evtN: ZZ == nullptr\n";*/
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
	
	return;
}

void VZZAnalyzer::end(TFile &){
	
	
	float elapsedSec = (float)(clock()-startTime)/CLOCKS_PER_SEC;
	int elapsedSecInt = (int)elapsedSec;
	cout<<"\nElapsed Time: "<<elapsedSec<<" s\t\t("<<elapsedSecInt/60<<"\' "<<elapsedSecInt%60<<"\")\n";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<" \tEnd of VZZAnalyzer\t";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<"\n\n";
}
