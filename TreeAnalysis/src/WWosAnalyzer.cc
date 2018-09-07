#include "VVXAnalysis/TreeAnalysis/interface/WWosAnalyzer.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <TF1.h>
#include <TH1.h>
#include <vector>
#include <algorithm>    // std::min
#include <time.h>
#include <string>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;
using std::vector;
using std::string;

using namespace phys;

Int_t WWosAnalyzer::cut() {
	static unsigned long evtN = 0;
	evtN++;
	cout<<evtN<<"\b\b\b\b\b\b";
	return 1;
}

void WWosAnalyzer::begin(){
	cout<<"\n";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<" \tBegin of WWos\t ";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<"\nTotal Events: "<<NUMBER_OF_EVENTS<<" \tAnalyzed:\n";
	startTime = clock();
}

void WWosAnalyzer::analyze(){
	#ifdef DO_STATISTICS_ON_PARTICLES
	initStatistics();
	#endif
	genElectrons 	= new std::vector<phys::Particle>();	//prompt genElectrons in this event
	genMuons 		= new std::vector<phys::Particle>();	//prompt genMuons in this event
	genLeptons 		= new std::vector<phys::Particle>();	//every prompt genLepton in this event
	
	foreach(const phys::Particle &gen, *genParticles){
		#ifdef DO_STATISTICS_ON_PARTICLES
		tempStatisticParticles(gen);
		#endif
		if(!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt) && abs(gen.eta()) < 2.4) ) {
			continue;
		}
		switch (abs(gen.id())){
			case 11:
				fillParticlePlots("genElectrons", gen);
				genElectrons->push_back(gen);
				genLeptons->push_back(gen);
				break;
			case 13:
				fillParticlePlots("genMuons", gen);
				genMuons->push_back(gen);
				genLeptons->push_back(gen);
				break;
			case 15:
				fillParticlePlots("genTaus", gen);
				genLeptons->push_back(gen);
				break;
			default:
				continue;
		}
		
	}
	#ifdef DO_STATISTICS_ON_EVENTS
	tempStatisticEvents();
	fillBasicPlots();		//Inherited from EventAnalyzer
	#endif
	std::sort(genElectrons->begin(), genElectrons->end(), PtComparator());	//Descending order
	std::sort(genMuons->begin(), genMuons->end(), PtComparator());
	std::sort(genLeptons->begin(), genLeptons->end(), PtComparator());
  	/*
  	bool genElectronsignal = false;
  	bool muonSignal = false;
  
  	//electrons
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
				met->p4().E() > 100 &&
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
	*/
	
	//Efficiency of particle reconstruction
	//analyzeEfficiency();
	analyzeEfficiency(genElectrons, electrons, string("Electrons"), matchedElectrons);
	analyzeEfficiency(genMuons, muons, string("Muons"), matchedMuons);
	
	totalElectrons += genElectrons->size();
	totalMuons += genMuons->size();
	
	delete genElectrons;
	delete genMuons;
	delete genLeptons;
	
	return;
}

void WWosAnalyzer::end(TFile &){
	//doSomeFits();
	//cout<<"Events: "<<NUMBER_OF_EVENTS<<" \tPassing selection: "<<passingSelection<<" \tEfficiency: "<<(1.-(float)passingSelection/NUMBER_OF_EVENTS)*100.<<" %\n";
	//cout<<"Electron events: "<<electronEvents<<" \tMuon events: "<<muonEvents<<"\n";
	cout<<"Total Electrons: "<<totalElectrons<<" \tMatched Electrons: "<<matchedElectrons<<" \tEfficiency: "<<(float)(100*matchedElectrons)/totalElectrons<<" %" <<"\n";
	cout<<"\t\t   Of wich with mismatched charge:   "<<wrongChargeE<<" \tRatio:      "<<(float)(100*wrongChargeE)/totalElectrons<<" %\n";
	cout<<"Total Muons:     "<<totalMuons<<" \tMatched Muons:     "<<matchedMuons<<" \tEfficiency: "<<(float)(100*matchedMuons)/totalMuons<<" %" <<"\n";
	cout<<"\t\t   Of wich with mismatched charge:    "<<wrongChargeM<<" \tRatio:      "<<(float)(100*wrongChargeM)/totalMuons<<" %\n";
	
	normalizeHistograms(string("Electrons"));
	normalizeHistograms(string("Muons"));
	
	cout<<"\nElapsed Time: "<< (float)(clock()-startTime)/CLOCKS_PER_SEC<<" s\n";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<" \tEnd of WWos\t";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<"\n\n";
}

//Efficiency analysis
template <class T, class P, typename C>
void WWosAnalyzer::analyzeEfficiency(vector<T>* genGroup, vector<P>* recGroup, std::string name, C& counter){
	foreach(const phys::Particle & gen, *genGroup){
		//fillParticlePlots("gen"+name, gen);
		/*theHistograms.fill("gen"+name+"_vs_eta","gen"+name+"_vs_eta", 26*8+1, -2.6,2.6,gen.eta(), 1);
		theHistograms.fill("gen"+name+"_vs_pt","gen"+name+"_vs_pt", 100, 0., 500., gen.pt(), 1);*/
		
		phys::Particle* cand = findMatchingParticle(gen, recGroup);	//candidate is a reconstructed particle
		if(cand != nullptr){
			counter++;
			theHistograms.fill(name+"Match_deltaR",name+"Match_deltaR", 200, 0, 2, physmath::deltaR(gen, *cand), 1.);
			if(checkMatch(gen, *cand, 0.2)){
				
				if(gen.charge() != (*cand).charge()){
					if(name == string("Electrons")) wrongChargeE++;
					if(name == string("Muons")) wrongChargeM++;
				}	
				//Efficiency vs Eta
				theHistograms.fill(name+"Matched_vs_eta", name+"Matched_vs_eta", 261, -2.6, 2.6, gen.eta(), 1);
				//Efficiency vs pt
				theHistograms.fill(name+"Matched_vs_pt",name+"Matched_vs_pt",100,0.,500.,gen.pt(),1);
			}
			//Efficiency vs Tolerance
			/*for(int cTolerance = 1; cTolerance <= 50; cTolerance++){ //I like to iterate on ints
				if(checkMatch(gen, *cand, (double)cTolerance/100.)){
					theHistograms.fill(name+"Efficiency_vs_tolerance",name+"Efficiency_vs_tolerance", 50, 0.01, 0.51, (float)cTolerance/100.+0.0001, 1);	//To be normalized in end()
				}
			}*/
		}
	}
}


//	Helper functions

phys::Particle* WWosAnalyzer::findMatchingParticle(const phys::Particle& rec, std::vector<phys::Lepton>* candidates){
	if(candidates->size() == 0) return nullptr;
	int minPos = 0;
	double deltaRMin = physmath::deltaR(rec, candidates->at(0));
	double temp = 999.;
	for(int i = 0; i < candidates->size(); i++){
		temp = physmath::deltaR(rec, candidates->at(i));
		if(temp < deltaRMin){
			deltaRMin = temp;
			minPos = i;
		}
	}
	// We've searched the candidate with the minmimum deltaR. Now we check if this makes sense
	/*if(checkMatch(rec, candidates->at(minPos), 0.1)){
		return & (candidates->at(minPos));
	} else return nullptr;	//Otherwise we can assume that rec is something faking an electron*/
	
	return & (candidates->at(minPos));
}

template <class P, class T>
bool WWosAnalyzer::checkMatch(const /*phys::Particle&*/P& reconstructed, const /*phys::Particle&*/ T& generated, const float& tolerance){
	//if(reconstructed.charge() != generated.charge()) return false;	//NO
	return physmath::deltaR(reconstructed, generated) < tolerance;
}

void WWosAnalyzer::normalizeHistograms(std::string name){
	theHistograms.clone(name+"Matched_vs_eta", name+"Efficiency_vs_eta");
	theHistograms.get(name+"Efficiency_vs_eta")->Divide(theHistograms.get("gen"+name+"_eta"));
	theHistograms.get(name+"Efficiency_vs_eta")->SetTitle((name+"Efficiency_vs_#eta").string::c_str());
	//There's no overload of SetTitle(const char*) with SetTitle(std::string)
	theHistograms.get(name+"Efficiency_vs_eta")->GetXaxis()->SetTitle("#eta");
	
	theHistograms.clone(name+"Matched_vs_pt", name+"Efficiency_vs_pt");
	theHistograms.get(name+"Efficiency_vs_pt")->Divide(theHistograms.get("gen"+name+"_pt"));
	theHistograms.get(name+"Efficiency_vs_pt")->SetTitle((name+"Efficiency_vs_pt").string::c_str());
	theHistograms.get(name+"Efficiency_vs_pt")->GetXaxis()->SetTitle("pt [GeV/c]");
	
	/*theHistograms.get(name+"Efficiency_vs_tolerance")->Scale(1./(float)totalElectrons);
	theHistograms.get(name+"Efficiency_vs_tolerance")->GetYaxis()->SetRangeUser(0., 1.);
	theHistograms.get(name+"Efficiency_vs_tolerance")->GetXaxis()->SetTitle("deltaR");*/
}


// Random statistics	------------------------------------------------------------------------

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
  theHistograms.fill(type+"_pt" ,    "p_{T} spectrum", 100,  0   , 500  ,lepton.pt() , 1/*theWeight*/);
  theHistograms.fill(type+"_eta",    "#eta spectrum" , 261, -2.6 ,  2.6 ,lepton.eta(), 1/*theWeight*/);
  theHistograms.fill(type+"_phi",    "#phi spectrum" ,  50, -3.15, 3.15 ,lepton.phi(), 1/*theWeight*/);
  //theHistograms.fill(type+"_charge", "charge"        ,  50,  -25  ,  25   ,lepton.charge(), theWeight);
  
  theHistograms.get(type+"_pt")->GetXaxis()->SetTitle("[GeV/c]");
}


#ifdef DO_STATISTICS_ON_PARTICLES
void WWosAnalyzer::initStatistics(){
	/*	cout << "------------------------------------------------------------------"<<endl;
	cout << "Run: " << run << " event: " << event << endl;*/
	
	particleCounter = 0;
	eCounter = 0;
	mCounter = 0;
	
	promptCounter = 0;
	peCounter = 0;
	pmCounter = 0;
}

void WWosAnalyzer::tempStatisticParticles(const phys::Particle &par){
	particleCounter++;
		if(abs(par.id()) == 11) eCounter++;
		if(abs(par.id()) == 13) mCounter++;
		
		if(par.genStatusFlags().test(phys::GenStatusBit::isPrompt)){
			promptCounter++;
			if(abs(par.id()) == 11) peCounter++;
			if(abs(par.id()) == 13) pmCounter++;
			
			theHistograms.fill("ptPromptParticle", "p_t (prompt)", 25, 0, 50, par.pt());
			//theHistograms.fill("massPromptParticle", "m (prompt)", 1000, 0, 100, par.p4().M());
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
	//theHistograms.get("massPromptParticle")->GetXaxis()->SetTitle("[GeV/c^2]");
	//theHistograms.fill("ptAllAnalyzedParticle", "p_t", 25, 0, 50, par.pt());
	//theHistograms.get("ptAllAnalyzedParticle")->GetXaxis()->SetTitle("[GeV]");
	//theHistograms.fill("ParticlesIDs", "ParticlesIDs", 40, -20, 20, par.id());
}
#endif
#ifdef DO_STATISTICS_ON_EVENTS
void WWosAnalyzer::tempStatisticEvents(){
	theHistograms.fill("genParticlesPerEvent", "genParticlesPerEvent", 40, 0, 40, particleCounter);
	theHistograms.fill("genElectronsPerEvent", "genElectronsPerEvent", 25, 0, 25, eCounter);
	theHistograms.fill("genMuonsPerEvent", "genMuonsPerEvent", 25, 0, 25, mCounter);
	
	theHistograms.fill("PromptGenParticlesPerEvent", "PromptGenParticlesPerEvent", 15, 0, 15, promptCounter);
	theHistograms.fill("PromptGenElectronsPerEvent", "PromptGenElectronsPerEvent", 10, 0, 10, peCounter);
	theHistograms.fill("PromptMuonsPerEvent", "PromptMuonsPerEvent", 10, 0, 10, pmCounter);
	
	theHistograms.fill("ElectronsPerEvent", "ElectronsPerEvent", 20, 0, 20, electrons->size());
}
#endif
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

void getFitInfo(TF1* p){
	cout<<"Chi^2: "<<p->GetChisquare()<<"\tNumber of DoF: "<<p->GetNDF()<<"\t(Pobability: "<<p->GetProb();
	cout<<")\n\n";
	for(int i = 0; i<70; i++) cout<<"-";
	cout<<"\n";
}

  
