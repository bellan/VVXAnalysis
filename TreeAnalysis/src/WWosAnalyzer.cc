#include "VVXAnalysis/TreeAnalysis/interface/WWosAnalyzer.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <TF1.h>
#include <vector>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;


using namespace phys;

Int_t WWosAnalyzer::cut() {
	return 1;
}

void WWosAnalyzer::begin(){
	cout<<"\n";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<" \tBegin of WWos\t ";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<"\n";
}

void WWosAnalyzer::analyze(){
	initStatistics();
	std::vector<phys::Particle>* genElectrons = new std::vector<phys::Particle>();	//genElectrons in this event
	std::vector<phys::Particle>* genMuons 		= new std::vector<phys::Particle>();	//muons in this event
	std::vector<phys::Particle>* leptons			= new std::vector<phys::Particle>();	//every lepton in this event
	
	foreach(const phys::Particle &gen, *genParticles){
		tempStatisticParticles(gen);
		if(!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt))) {
			continue;
		}
		switch (abs(gen.id())){
			case 11:
				fillParticlePlots("genElectrons", gen);
				genElectrons->push_back(gen);
				leptons->push_back(gen);
			case 13:
				fillParticlePlots("genMuons", gen);
				genMuons->push_back(gen);
				leptons->push_back(gen);
			case 15:
				fillParticlePlots("taus", gen);
				leptons->push_back(gen);
			default:
				continue;
		}
		
	}
	tempStatisticEvents();
	//fillBasicPlots();		//Inherited from EventAnalyzer
	inTheLastEvent();	
	std::sort(genElectrons->begin(), genElectrons->end(), PtComparator());	//Descending order
	std::sort(genMuons->begin(), genMuons->end(), PtComparator());
	std::sort(leptons->begin(), leptons->end(), PtComparator());
  
  bool genElectronsignal = false;
  bool muonSignal = false;
  
  //eletrons
  if(genElectrons->size() >=2){
		if(
				genElectrons->at(0).pt() > 25 && 
				met->p4().E() > 20 &&
				genElectrons->at(1).pt() > 15 && 
				genElectrons->at(0).charge() * genElectrons->at(1).charge() == -1
			)
			genElectronsignal = true;
	}
	else genElectronsignal = false;
	//muons
	if(genMuons->size() >=2){
		if(
				genMuons->at(0).pt() > 20 && 
				met->p4().E() > 20 &&
				genMuons->at(1).pt() > 15 && 
				genMuons->at(0).charge() * genMuons->at(1).charge() == -1
			) 
			muonSignal = true;
	}
	else muonSignal = false;
	
	if(genElectronsignal && !muonSignal){ 
		electronEvents++;
	}
	if(muonSignal && !genElectronsignal){ 
		muonEvents++;
		theHistograms.fill("Muon_E_analyzed", "Muon_E_analyzed", 200, 0, 200, genMuons->at(0).p4().E());
		theHistograms.fill("Muon_pt_analyzed", "Muon_pt_analyzed", 200, 0, 200, genMuons->at(0).pt());
	}
	if((genElectronsignal && !muonSignal) || (muonSignal && !genElectronsignal)) passingSelection++;
	
	//Efficiency
	foreach(const phys::Particle & ele, *electrons){
		findElectronMatch(ele, genElectrons);
	}
	totalElectrons += genElectrons->size();
	
	/*
	delete genElectrons;
	delete genMuons;
	delete leptons;
	*/
	
	return;
}

void WWosAnalyzer::end(){
	//doSomeFits();
	cout<<"Events: "<<9055<<" \tPassing selection: "<<passingSelection<<" \tEfficiency: "<<(1.-(float)passingSelection/9055)*100<<" %\n";
	cout<<"Electron events: "<<electronEvents<<" \tMuon events: "<<muonEvents<<"\n";
	cout<<"Matched Electrons: "<<matchedElectrons<<" \tTotal Electrons: "<<totalElectrons<<" \tEfficiency: "<<(float)matchedElectrons/totalElectrons*100.<< "%" <<"\n";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<" \tEnd of WWos\t ";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<"\n\n";
}


//	Helper functions
void WWosAnalyzer::findElectronMatch(const phys::Particle & rec, std::vector<phys::Particle> * gen){
	foreach(const phys::Particle & genEle, *gen){
		if(electronMatch(rec, genEle)){
			matchedElectrons++;
			return;
		}
	}
}

bool WWosAnalyzer::electronMatch(const phys::Particle & reconstructed, const phys::Particle & generated){
	//if(reconstructed.charge() != generated.charge()) return false;
	cout<<physmath::deltaR(reconstructed, generated);
	//return physmath::deltaR(reconstructed, generated) < 10.;
	return true;
}


void WWosAnalyzer::inTheLastEvent(){	//For some reasons, end() is not executed at the end of the analysis
	static int event = 0;
	event++;
	if(event == 9055){
		end();
	}
	return;
}


void WWosAnalyzer::fillBasicPlots(){
  theHistograms.fill<TH1I>("nvtx"     , "Number of vertices" , 100, 0, 100, nvtx             , theWeight);
  theHistograms.fill      ("met"      , "Missing energy"     , 200, 0, 800, met->pt()        , theWeight);
  theHistograms.fill<TH1I>("nmuons"    ,"Number of muons"    ,  10, 0, 10 , muons->size()    , theWeight);
  theHistograms.fill<TH1I>("nelectrons","Number of electrons",  10, 0, 10 , electrons->size(), theWeight);
	
	theHistograms.get("met")->GetXaxis()->SetTitle("[GeV/c]");
  /*foreach(const phys::Lepton& mu , *muons)     fillLeptonPlots  ("mu",  mu  );
  foreach(const phys::Electron& e, *electrons) fillElectronPlots("e" ,  e   );
  foreach(const phys::Jet& jet   , *jets)      fillJetPlots     ("jet", jet );*/
}



void WWosAnalyzer::fillParticlePlots(const std::string &type, const phys::Particle & lepton){
  theHistograms.fill(type+"_pt" ,    "p_{T} spectrum", 100,   0   , 500   ,lepton.pt()    , theWeight);
  theHistograms.fill(type+"_eta",    "#eta spectrum" , 100,  -5 ,   5 ,lepton.eta()   , theWeight);
  theHistograms.fill(type+"_phi",    "#phi spectrum" ,  50,  -3.15,   3.15,lepton.phi()   , theWeight);
  //theHistograms.fill(type+"_charge", "charge"        ,  50,  -25  ,  25   ,lepton.charge(), theWeight);
  
  theHistograms.get(type+"_pt")->GetXaxis()->SetTitle("[GeV/c]");
}



void WWosAnalyzer::initStatistics(){
	/*	cout << "------------------------------------------------------------------"<<endl;
	cout << "Run: " << run << " event: " << event << endl;*/
	
	counter = 0;
	eCounter = 0;
	mCounter = 0;
	
	promptCounter = 0;
	peCounter = 0;
	pmCounter = 0;
}

void WWosAnalyzer::tempStatisticParticles(const phys::Particle &par){
	counter++;
		if(abs(par.id()) == 11) eCounter++;
		if(abs(par.id()) == 13) mCounter++;
		
		if(par.genStatusFlags().test(phys::GenStatusBit::isPrompt)){
			promptCounter++;
			if(abs(par.id()) == 11) peCounter++;
			if(abs(par.id()) == 13) pmCounter++;
			
			theHistograms.fill("ptPromptParticle", "p_t (prompt)", 25, 0, 50, par.pt());
			theHistograms.fill("massPromptParticle", "m (prompt)", 1000, 0, 100, par.p4().M());
		}
	/*	if(
				(abs(gen.id()) != 11 && abs(gen.id()) != 13) || 
				(!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || 
				!(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))
			) 				
			continue;
	*/
	//	cout << "id: " << gen.id() << " pt: " << gen.pt() << endl;
	theHistograms.get("ptPromptParticle")->GetXaxis()->SetTitle("[GeV/c]");
	theHistograms.get("massPromptParticle")->GetXaxis()->SetTitle("[GeV/c^2]");
	theHistograms.fill("ptAllAnalyzedParticle", "p_t", 25, 0, 50, par.pt());
	theHistograms.get("ptAllAnalyzedParticle")->GetXaxis()->SetTitle("[GeV]");
	theHistograms.fill("ParticlesIDs", "ParticlesIDs", 40, -20, 20, par.id());
}

void WWosAnalyzer::tempStatisticEvents(){
	theHistograms.fill("ParticlesPerEvent", "ParticlesPerEvent", 40, 0, 40, counter);
	theHistograms.fill("genElectronsPerEvent", "genElectronsPerEvent", 25, 0, 25, eCounter);
	theHistograms.fill("genMuonsPerEvent", "genMuonsPerEvent", 25, 0, 25, mCounter);
	
	theHistograms.fill("PromptParticlesPerEvent", "PromptParticlesPerEvent", 25, 0, 25, promptCounter);
	theHistograms.fill("PromptgenElectronsPerEvent", "PromptgenElectronsPerEvent", 25, 0, 25, peCounter);
	theHistograms.fill("PromptMuonsPerEvent", "PromptMuonsPerEvent", 25, 0, 25, pmCounter);
	
	theHistograms.fill("ElectronsPerEvent", "ElectronsPerEvent", 25, 0, 25, electrons->size());
}

void WWosAnalyzer::doSomeFits(){
	TF1* func1 = new TF1("func1","[0]*pow(x,[1])*exp(-x/2)",2,25);
	func1->SetParLimits(1,1.5,2.0);
	func1->SetLineColor(4);
	
	TF1* func2 = new TF1("func2","[0]*pow(x,[1])*exp(-x/2)",2,25);
	func2->SetParLimits(1,1.75,2.0);
	func2->SetLineColor(3);
	
	TF1* func3 = new TF1("func3","[2]+[0]*pow(x,[1])*exp(-x/2)",0,25);
	func3->SetParameter(1,1.5);
	func3->SetLineColor(1);
	
	TH1* ePlot = theHistograms.get("genElectronsPerEvent");
	ePlot->Fit(func1,"R");
	getFitInfo(func1);
	ePlot->Fit(func2,"R+");
	getFitInfo(func2);
	ePlot->Fit(func3,"R+");
	getFitInfo(func3);
		
	ePlot->Draw("same");
	func1->Draw("same");
	func2->Draw("same");
}

void WWosAnalyzer::getFitInfo(TF1* p){
	cout<<"Chi^2: "<<p->GetChisquare()<<"\tNumber of DoF: "<<p->GetNDF()<<"\t(Pobability: "<<p->GetProb();
	cout<<")\n\n";
	for(int i = 0; i<70; i++) cout<<"-";
	cout<<"\n";
}

  
