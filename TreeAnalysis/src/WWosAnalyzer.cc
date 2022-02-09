#include "VVXAnalysis/TreeAnalysis/interface/WWosAnalyzer.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <TF1.h>
#include <TH1.h>
#include "TTree.h"
//#include <TEfficiency.h> // ClopperPearson(Double_t, Double_t, Double_t, Bool_t)
#include <vector>
#include <algorithm>    // std::min
#include <time.h>				//clock_t, clock()
#include <string>
#include <iomanip>      // std::setprecision

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;
using std::vector;
using std::string;

using namespace phys;

//Uniform configuration for histograms to avoid errors. Format: <nbins>,<ufirst>,<ulast>
#define HISTO_deltaR_Config 	200,0.,0.5
#define HISTO_Eta_Config 			261,-2.6,2.6
#define HISTO_JEta_Config 		150,-5.,5.
#define HISTO_Phi_Config 			100,-3.15,3.15
#define HISTO_Pt_Config 			100,0.,500.
#define HISTO_CSV_Config			52,-0.03,1.01//102,-0.015,1.005
#define HISTO_Cut_Config			4,0,4

Int_t WWosAnalyzer::cut() {
	evtN++;
	cout<<"\r\t\t"<<evtN;//"\b\b\b\b\b\b";
	
	theHistograms->fill("PreCut", "PreCut", 5,0,5, 0, theWeight);	//Total
	//#####		At least a pair of jets with mjj > 150 GeV/c^2		#####
	bool jetOK = false;

	if(jets->size() < 2) return -1;
	for(uint i = 0; i < jets->size()-1; ++i){
		if(jetOK) break;		//Saves computation time
    for(uint k = i+1; k < jets->size(); ++k){
    	double mjj = ( jets->at(i).p4() + jets->at(k).p4() ).M();
    	if (mjj > 150){
    		jetOK = true; 
    		break;
  		}
  	}
	}
	if(!jetOK || jets->size() < 2) return -1;	//Saves computation time, avoiding unnecessary lepton analysis
	theHistograms->fill("PreCut", "PreCut", 5,0,5, 1, theWeight);	//mjj
	
	//#####		At least a pair of leptons with pt > 20 GeV/c		#####
  bool leptOK = false;
  int countMu20 = 0;
  int countEle20 = 0;
  float maxPt = 20.;
  float secPt = 20.;
  thisEventType = WWosEventTypes::ee;
  foreach(const Lepton& e, *electrons){
  	countEle20 += (e.pt() > 20 ? 1 : 0);
  	if(e.pt() > maxPt){
	  	maxPt = e.pt();
  		leadLepton = e;		//leadLepton and trailLepton are used in analyze()
		}
		else if(e.pt() > secPt){
	  	secPt = e.pt();
  		trailLepton = e;
		}
	}
  	
  foreach(const Lepton& m, *muons){
  	countMu20 += (m.pt() > 20 ? 1 : 0);
  	if(m.pt() > maxPt){
	  	maxPt = m.pt();
  		leadLepton = m;
  	//	if(thisEventType == WWosEventTypes::ee) thisEventType = WWosEventTypes::em;
  	//	else if(thisEventType == WWosEventTypes::em) thisEventType = WWosEventTypes::mm;
		}
		else if(m.pt() > secPt){
	  	secPt = m.pt();
  		trailLepton = m;
  	//	if(thisEventType == WWosEventTypes::ee) thisEventType = WWosEventTypes::em;
  	//	else if(thisEventType == WWosEventTypes::em) thisEventType = WWosEventTypes::mm;
		}
	}
  if(countEle20 + countMu20 >= 2) leptOK = true;	
  if(countEle20 + countMu20 > 2) lostEvents++;	//Events we loose by cuttig events with a third lepton with pt > 20 GeV/c
  
  /*if(false){
  cout<<"\n\nType: "<<eventType<<" \t"<<;
	cout<<"\nElectrons pt: \t\t";
	foreach(const Lepton& el, *electrons) cout<<el.pt()<<" - ";
	cout<<"\nMuons pt: \t\t";
	foreach(const Lepton& mu, *muons) cout<<mu.pt()<<" - ";
	cout<<"\nmaxPt: "<<maxPt<<" \tsecPt: "<<secPt;
  }*/
	
	if(!leptOK)	return -1;
	
	theHistograms->fill("PreCut", "PreCut", 5,0,5, 2, theWeight);	//2 leptons pt>20
	
	if(leptonCut())	//opposite charge and, if they're same flavour, then mll>120
		return -1;

	return 1;
}

void WWosAnalyzer::begin(){
	cout<<"\n";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<" \tBegin of WWos\t ";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<"\nAnalyzed:\t      /"<<tree()->GetEntries()<<std::flush;
	startTime = clock();
	evtN = 0;
}

void WWosAnalyzer::analyze(){
	passingCut++;
	TH1::SetDefaultSumw2(true);
	#ifdef DO_GEN_PARTICLES_ANALYSIS
	genParticlesAnalysis();
	#endif
	
	#ifdef LEPTON_CUT
	//#####		pre-second cut		#####
	string eventType = "";/*
	//IMPORTANT: TO HAVE CUMULATIVE GRAPHS, I'M OVERRIDING THE EVENT TYPE (only the string)
	if(thisEventType == WWosEventTypes::ee){			eeEvents++;		eventType = "ee"; }
	else if(thisEventType == WWosEventTypes::em){	emEvents++;		eventType = "em"; }
	else if(thisEventType == WWosEventTypes::mm){	mmEvents++;		eventType = "mm"; }*/
		
	std::sort(jets->begin(), jets->end(), phys::EComparator());	//They were sorted by pt
	if(jetCsvtagCut(0.95, 0.95)) return;	//Some samples have a prouction cut at 0.95
	
	jetPlots(" all");
	miscPlots(&leadLepton, &trailLepton, " all", true);
	leptonPlots(&leadLepton, &trailLepton, " lep all"/*eventType*/, true);
		
	theHistograms->fill("Cuts", "Cuts", HISTO_Cut_Config, 0, theWeight);	//0th cut: Preselected 
	theHistograms->fill("Cuts unweighted","Cuts unweighted",HISTO_Cut_Config, 0,1);
	
	if(newCuts()) return;	//automatically fills the "Cuts" plot
	
	passingSelection++;
	//lepton2DGraph(&leadLepton, &trailLepton, eventType, true);
	jetPlots();
	miscPlots(&leadLepton, &trailLepton, eventType, true);
	leptonPlots(&leadLepton, &trailLepton, "lep"/*eventType*/, true);
	#endif
	
	return;
}

void WWosAnalyzer::end(TFile &){
	cout << std::setprecision(4);
	cout<<"\nEvents: "<<evtN<<" \tPassing cut: "<<passingCut<<" \t(of the total): "<<((float)passingCut/evtN)*100.<<" %\n";
	cout<<"Events\t\tee: "<<eeEvents<<" \tmm : "<<mmEvents<<" \tem: "<<emEvents<<'\n';
	cout<<"Passing second cut: "<<passingSelection<<" \t(of the analyzed): "<<((float)passingSelection/passingCut)*100.<<" % \t(of the total): "<<((float)passingSelection/evtN)*100.<<" %\n";
	cout<<"Events with >=3 leptons: "<<lostEvents<<" \t(of the selected): "<<((float)lostEvents/passingSelection)*100.<<" % \t(of the analyzed): "<<((float)lostEvents/passingCut)*100.<<" % \t(of the total): "<<((float)lostEvents/evtN)*100.<<" %\n";
	
	nameCutGraph();
	
	#ifdef DO_GEN_PARTICLES_ANALYSIS
	endGenParticleAnalysis();
	#endif
	
	float elapsedSec = (float)(clock()-startTime)/CLOCKS_PER_SEC;
	int elapsedSecInt = (int)elapsedSec;
	cout<<"\nElapsed Time: "<<elapsedSec<<" s\t\t("<<elapsedSecInt/60<<"\' "<<elapsedSecInt%60<<"\")\n";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<" \tEnd of WWos\t";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<"\n\n";
}


bool WWosAnalyzer::leptonCut(){
	if(leadLepton.charge() * trailLepton.charge() > 0) return true;	// true = cut this event
	theHistograms->fill("PreCut", "PreCut", 5,0,5, 3, theWeight);	//Opposite charge
	
	float mll = (leadLepton.p4() + trailLepton.p4()).M();
	if(mll < 120 && (leadLepton.id() == -1*trailLepton.id())) return true; //mLL
	
	theHistograms->fill("PreCut", "PreCut", 5,0,5, 4, theWeight);	//mll
	return false;
}

bool WWosAnalyzer::jetCsvtagCut(const float &leadMax, const float &trailMax){
	if(jets->at(0).csvtagger() > leadMax)  return true;
	if(jets->at(1).csvtagger() > trailMax) return true;
	return false;
}

bool WWosAnalyzer::newCuts(){
	UInt_t ncut = 0;
	
	if( fabs(jets->at(0).eta() - jets->at(1).eta()) < 5. ) return true;
	theHistograms->fill("Cuts", "Cuts", HISTO_Cut_Config, ++ncut, theWeight); //1st cut: deltaEta > 5
	theHistograms->fill("Cuts unweighted","Cuts unweighted",HISTO_Cut_Config, ncut,1);
	
	if((leadLepton.p4() + trailLepton.p4()).Pt() < 120) return true;
	theHistograms->fill("Cuts", "Cuts", HISTO_Cut_Config, ++ncut, theWeight); //2nd cut: ptVectLL > 120
	theHistograms->fill("Cuts unweighted","Cuts unweighted",HISTO_Cut_Config, ncut,1);
	
	if(leadLepton.pt() + trailLepton.pt() + met->pt() < 255) return true;
	theHistograms->fill("Cuts", "Cuts", HISTO_Cut_Config, ++ncut, theWeight); //3rd cut: ptScalarLLmet> 255
	theHistograms->fill("Cuts unweighted","Cuts unweighted",HISTO_Cut_Config, ncut,1);
	
	/*
	//if((leadLepton.p4() + trailLepton.p4() + met->p4()).Pt() < 100) return true;
	//theHistograms->fill("Cuts", "Cuts", HISTO_Cut_Config, ++ncut, theWeight);
	
	if(jetCsvtagCut(0.8484, 0.8484)) return true;	//0.9535 (Thight) or 0.8484 (Medium)
	theHistograms->fill("Cuts", "Cuts", HISTO_Cut_Config, ++ncut, theWeight); //4th cut: medium jet csvtag
	*/
/*	if( jetCsvtagCut(0.4, 0.4) ) return true;
	theHistograms->fill("Cuts", "Cuts", HISTO_Cut_Config, ++ncut, theWeight); //5th cut: csvtag<0.4
	*/
/*	if(jets->at(0).e() < 150) return true;
	theHistograms->fill("Cuts", "Cuts", HISTO_Cut_Config, ++ncut, theWeight); //6th cut: Lead jet E > 150
	*/
	return false;
}

void WWosAnalyzer::nameCutGraph(){
	theHistograms->get("Cuts")->SetMinimum(1.);
	UInt_t ncut = 0;
	TAxis* aXCuts = theHistograms->get("Cuts")->GetXaxis();
	TAxis* aXUnweighted = theHistograms->get("Cuts unweighted")->GetXaxis();
	/*
	aXCuts->SetBinLabel(0+1, "Total");
	aXCuts->SetBinLabel(1+1, "mJJ > 150");
	aXCuts->SetBinLabel(2+1, "lepton pt>20");
	aXCuts->SetBinLabel(3+1, "lepton charge");
	aXCuts->SetBinLabel(4+1, "mLL > 120");	//if not same flavour, in that case pass anyway
	aXCuts->SetBinLabel(5+1, "medium jet csvtag");*/

	aXCuts->SetBinLabel(++ncut, "Preselected");        //0
	aXUnweighted->SetBinLabel(ncut, "Preselected");
	
	aXCuts->SetBinLabel(++ncut, "deltaEta (jets) > 5");//1
	aXUnweighted->SetBinLabel(ncut, "deltaEta (jets) > 5");
	
	aXCuts->SetBinLabel(++ncut, "ptVectLL > 120");     //2
	aXUnweighted->SetBinLabel(ncut, "ptVectLL > 120");
	
	aXCuts->SetBinLabel(++ncut, "ptScalarLLmet > 255");//3
	aXUnweighted->SetBinLabel(ncut, "ptScalarLLmet > 255");
	
	TH1* precut = theHistograms->get("PreCut");
	precut->SetTitle("PreCut;;# Events");
	TAxis* aXPreCuts = precut->GetXaxis();
	aXPreCuts->SetBinLabel(1, "Total");
	aXPreCuts->SetBinLabel(2, ">2jet && mjj>150");
	aXPreCuts->SetBinLabel(3, "2lep (pt>20)");
	aXPreCuts->SetBinLabel(4, "charge");
	aXPreCuts->SetBinLabel(5, "mll>120");
}

//	####################### Plots #################################

void WWosAnalyzer::leptonPlots(const Particle* lead, const Particle* trail, const string& type, bool useWeight){
	double w = (useWeight ? theWeight : 1.);
	theHistograms->fill("Leading "+type+" pt", "Leading "+type+" pt;pt [GeV/c]", HISTO_Pt_Config, lead->pt(), w);
	theHistograms->fill("Second "+type+" pt", "Second "+type+" pt;pt [GeV/c]", HISTO_Pt_Config, trail->pt(), w);
	float deltaEta = lead->eta() - trail->eta();
	theHistograms->fill("#Delta #eta "+type, "#Delta #eta "+type+";#Delta#eta", HISTO_JEta_Config, deltaEta, w);
	float dPhi = physmath::deltaPhi(lead->p4(), trail->p4());//lead->phi() - trail->phi();
	theHistograms->fill("#Delta #phi "+type, "#Delta #phi "+type+";#phi [rad]", 70, 0., 3.15, dPhi, w);
	float deltaR = physmath::deltaR(*lead, *trail);
	theHistograms->fill("#Delta R "+type, "#Delta R "+type+"; #Delta R", 100, 0., 5., deltaR, w);
}

void WWosAnalyzer::lepton2DGraph(const Particle* lead, const Particle* trail, const string& type, bool useWeight){
	double w = (useWeight ? theWeight : 1.);
	theHistograms->fill<TH2F>(type+": pt", type+": pt", HISTO_Pt_Config, HISTO_Pt_Config, lead->pt(), trail->pt(), w);
}


void WWosAnalyzer::jetPlots(const string& type){
	double w = theWeight;
	theHistograms->fill("n Jets", "n Jets", 10, 0., 10., jets->size(), w);
	//if(jets->size() >= 1){	//we already asked jets->size() >= 2 in cut()
	theHistograms->fill("Lead jet E"+type, "Lead jet E"+type+";E [GeV]",   100,0.,800.,jets->at(0).e(), w);
	theHistograms->fill("Lead jet pt"+type,"Lead jet pt"+type+";E [GeV/c]",100,0.,800.,jets->at(0).pt(),w);
	//	if(jets->size() >= 2){
	theHistograms->fill("Trail jet E"+type, "Trail jet E"+type+";E [GeV]",   100,0.,700., jets->at(1).e(),  w);
	theHistograms->fill("Trail jet pt"+type,"Trail jet pt"+type+";E [GeV/c]",100,0.,700., jets->at(1).pt(), w);
	foreach(const Particle& jet, *jets){
		theHistograms->fill("Jet Eta"+type, "Jet Eta"+type+";#eta", 100, -5., 5.,jet.eta(), w);
	}
	float deltaR = physmath::deltaR(jets->at(0), jets->at(1));
	theHistograms->fill("delta R jets"+type, "delta R jets"+type+"#Delta R", 100, 0., 10., deltaR, w);
	float deltaEta = fabs(jets->at(0).eta() - jets->at(1).eta());
	theHistograms->fill("delta Eta jets"+type, "delta Eta jets"+type+"#Delta#eta", 100,0.,10., deltaEta,w);
	float jjMass = (jets->at(0).p4() + jets->at(1).p4()).M();
	theHistograms->fill("mjj"+type,"mjj"+type+";mass [GeV/c^2]", 100, 0., 2000., jjMass, w);
	float angSep = jets->at(0).p4().Vect().Angle( jets->at(1).p4().Vect() );
	theHistograms->fill("angSep (jet)"+type, "angSep (jet)"+type+";[rad]", 100,0.,3.15, angSep, w);
	
	float leadCSV = jets->at(0).csvtagger();
	if(leadCSV < 0.) leadCSV = -0.01;
	float trailCSV = jets->at(1).csvtagger();
	if(trailCSV <0.) trailCSV =-0.01;
	theHistograms->fill("lead jet csvtagger"+type, "lead jet csvtagger"+type, HISTO_CSV_Config, leadCSV, w);
	theHistograms->fill("trail jet csvtagger"+type,"trail jet csvtagger"+type,HISTO_CSV_Config,trailCSV, w);
	
	//cout<<'\n'<<type<<' '<<jets->at(0).csvtagger()<<' '<<jets->at(1).csvtagger()<<' '<<w;
	
	//theHistograms->fill<TH2F>("jets csvtagger"+type, "jets csvtagger"+type, 100,0.,1., 100,0.,1., jets->at(0).csvtagger(), jets->at(1).csvtagger(), w);
//}
//}
}

void WWosAnalyzer::miscPlots(const Particle* lead, const Particle* trail, const string& type, bool useWeight){
	double w = (useWeight ? theWeight : 1.);
	
	float ptScalarLL = lead->pt() + trail->pt();
	theHistograms->fill("ptScalarLL"+type, "ptScalarLL"+type+";pt [GeV/c]", 100, 0., 1500., ptScalarLL, w);
	
	TLorentzVector Lp4 = lead->p4();
	TLorentzVector Tp4 = trail->p4();
	float ptVectLL = (Lp4 + Tp4).Pt();
	theHistograms->fill("ptVectLL"+type, "ptVectLL"+type+";pt [GeV/c]", 100, 0., 1000., ptVectLL, w);
	
	float mll = (Lp4 + Tp4).M();
	theHistograms->fill("mll"+type, "mll"+type+";mass [GeV/c^2]", 100, 0., 1000., mll, w);
	
	float angularSeparation = Lp4.Vect().Angle(Tp4.Vect());
	theHistograms->fill("angSeparationLL"+type, "angSeparationLL"+type+";[rad]", 50,0.,3.14, angularSeparation, w);
	
	//met
	float ptScalarLLmet = ptScalarLL + met->pt();
	theHistograms->fill("ptScalarLLmet"+type,"ptScalarLLmet"+type+";pt [GeV/c]", 100,0.,1000.,
	ptScalarLLmet, w);
	
	float ptVectLLmet = (Lp4 + Tp4 + met->p4()).Pt();
	theHistograms->fill("ptVectLLmet"+type,"ptVectLLmet"+type+";pt [GeV/c]", 100,0.,1000.,ptVectLLmet, w);
}


#ifdef DO_GEN_PARTICLES_ANALYSIS
//	##################################	genParticles ##################################
void WWosAnalyzer::genParticlesAnalysis(){
	genElectrons 		= new std::vector<phys::Particle>();	//prompt genElectrons in this event
	genMuons 				= new std::vector<phys::Particle>();	//prompt genMuons in this event
	genLeptons 			= new std::vector<phys::Particle>();	//every genLepton in this event (no neutrinos)
	genNeutrinos 		= new std::vector<phys::Particle>();
	genCleanedJets	= new std::vector<phys::Particle>();	//every prompt Jet in this event
	fakeJets				= new std::vector<phys::Particle>();	//isolated particles that appear in genJets
	
	genParticlesCategorization();	//foreach(	, *genParticles)
	
	theHistograms->fill("n genNeutrinos", "n genNeutrinos", 10, 0., 10., genNeutrinos->size(), theWeight);
	
	totalElectrons += genElectrons->size();
	totalMuons += genMuons->size();
	totjets += jets->size();
	totgenJets += genJets->size();
	totgenCleanedJets += genCleanedJets->size();
	totgenParticles += genParticles->size();
	
	//Efficiency of particle reconstruction
	#ifdef DO_EFFICIENCY_ANALYSIS
	analyzeEfficiency(genElectrons, electrons, string("Electrons"), matchedElectrons);
	analyzeEfficiency(genMuons, muons, string("Muons"), matchedMuons);
	analyzeEfficiency(genCleanedJets, jets, string("Jets"), matchedJets, 0.4);
	#endif
	
	genElectrons->clear();
	genMuons->clear();
	genLeptons->clear();	
	genCleanedJets->clear();
	genNeutrinos->clear();

	delete genElectrons;
	delete genMuons;
	delete genLeptons;	
	delete genCleanedJets;
	delete genNeutrinos;
}


void WWosAnalyzer::genParticlesCategorization(){	//Divides the genParticles betwenn e, mu, jet, etc.
	#ifdef DO_STATISTICS_ON_PARTICLES
	initStatistics();
	#endif
	foreach(const phys::Particle &gen, *genParticles){
		#ifdef DO_STATISTICS_ON_PARTICLES
		tempStatisticParticles(gen);
		#endif
		if(
				abs(gen.eta()) > 2.4 && 
				(abs(gen.id()) == 11 || abs(gen.id()) == 13 || abs(gen.id()) == 15 || abs(gen.id()) == 22) &&
				gen.pt() >20.
			)
			fakeJets->push_back(gen);
		if(
				!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || 
				!(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)) || 
				abs(gen.eta()) > 2.4 /*|| gen.pt() < 25*/
			){
			
			continue;
		}
		switch (abs(gen.id())){
			case 11:
				fillParticlePlots("genElectrons", gen);
				genElectrons->push_back(gen);
				genLeptons->push_back(gen);
				break;
			case 12:
				genNeutrinos->push_back(gen);
				break;
			case 13:
				fillParticlePlots("genMuons", gen);
				genMuons->push_back(gen);
				genLeptons->push_back(gen);
				break;
			case 14:
				genNeutrinos->push_back(gen);
				break;
	/*	case 15:
				fillParticlePlots("genTaus", gen);
				genLeptons->push_back(gen);
				break;*/
	/*	case 16:
				genNeutrinos->push_back(gen);
				break;*/
			default:
				continue;
		}
	}
	
	foreach(const phys::Particle &jet, *genJets){
		#ifdef DO_EFFICIENCY_ANALYSIS
		if(abs(jet.eta()) < 4.7){
			phys::Particle* cand = findMatchingParticle(jet, fakeJets/*genParticles*/); //
			if(cand != nullptr){
				if(checkMatch(jet, *cand, 0.4) && (abs(cand->id()) >= 11 && abs(cand->id())<=16)) continue; //All particles are also jets, this is a duplicate
			}
		/*else{
				cout<<"\nFound a nullptr\n";
				cout<<"\n\tgenJet:\n";
				cout<<jet<<"\n";
				cout<<"\n\tgenParticles:" << genLeptons->size()<<endl;
				foreach(Particle& genP, *genParticles){
					cout<<genP<<" - "<<(genP.genStatusFlags().test(phys::GenStatusBit::isPrompt))<<" - "<<(genP.genStatusFlags().test(phys::GenStatusBit::fromHardProcess))<<"\n";
				}
				for(int i=0; i<50; i++) cout<<"-";
			}*/
		}
		#endif
		genCleanedJets->push_back(jet);	//genCleanedJets
		fillParticlePlots("genJets", jet);
	}/*
	foreach(const phys::Particle &jet, *centralGenJets){
		if(jet.genStatusFlags().test(phys::GenStatusBit::isPrompt) && abs(jet.eta()) < 5.){
			fillParticlePlots("genJets", jet);
			genJets->push_back(jet);
		}
	}*/
	#ifdef DO_STATISTICS_ON_EVENTS
	tempStatisticEvents();
	fillBasicPlots();		//Inherited from EventAnalyzer
	#endif
	std::sort(genElectrons->begin(), genElectrons->end(), PtComparator());	//Descending order
	std::sort(genMuons->begin(), genMuons->end(), PtComparator());
	std::sort(genLeptons->begin(), genLeptons->end(), PtComparator());
	std::sort(genCleanedJets->begin(), genCleanedJets->end(), PtComparator());
}

void WWosAnalyzer::endGenParticleAnalysis(){
	//doSomeFits();
	
	#ifdef DO_EFFICIENCY_ANALYSIS
	cout<<"\nTotal Electrons: "<<totalElectrons<<" \tMatched Electrons: "<<matchedElectrons<<" \tEfficiency: "<<(float)(100*matchedElectrons)/totalElectrons<<" %" <<"\n";
	cout<<"\t\t   Of wich with mismatched charge:   "<<wrongChargeE<<" \tRatio:      "<<(float)(100*wrongChargeE)/totalElectrons<<" %\n";
	cout<<"Total Muons:     "<<totalMuons<<" \tMatched Muons:     "<<matchedMuons<<" \tEfficiency: "<<(float)(100*matchedMuons)/totalMuons<<" %" <<"\n";
	cout<<"\t\t   Of wich with mismatched charge:    "<<wrongChargeM<<" \tRatio:      "<<(float)(100*wrongChargeM)/totalMuons<<" %\n";
	cout<<"jets: "<<totjets<<" \tgenJets: "<<totgenJets<<" \tgenCleanedJets: "<<totgenCleanedJets<<" \ttotgenParticles: "<<totgenParticles<<"\n";
	
	normalizeHistograms(string("Electrons"));
	normalizeHistograms(string("Muons"));
	normalizeHistograms(string("Jets"));
	
	cout<<"Resolution 1/Pt (Electrons)\n";
	fitResolutionPt(string("Electrons"));
	cout<<"Resolution 1/Pt (Muons)\n";
	fitResolutionPt(string("Muons"));
	cout<<"Resolution E (Electrons)\n";
	fitResolutionE(string("Electrons"));
	cout<<"Resolution E (Muons)\n";
	fitResolutionE(string("Muons"));
	
	#endif
}
#endif	//DO_GEN_PARTICLES_ANALYSIS

//Efficiency analysis
#ifdef DO_EFFICIENCY_ANALYSIS
template <class T, class P, typename C>
void WWosAnalyzer::analyzeEfficiency(vector<T>* genGroup, vector<P>* recGroup, std::string name, C& counter, Float_t maxDeltaR){
	foreach(const phys::Particle & gen, *genGroup){
		
		phys::Particle* cand = findMatchingParticle(gen, recGroup);	//candidate is a reconstructed particle
		if(cand != nullptr){
			counter++;
			theHistograms->fill(name+"Match_deltaR",name+"Match_deltaR;#Delta R", 100,0.,0.5 , physmath::deltaR(gen, *cand), 1.);
			if(checkMatch(gen, *cand, maxDeltaR)){
				
				if(gen.charge() != (*cand).charge()){
					if(name == "Electrons") wrongChargeE++;
					if(name == "Muons") wrongChargeM++;
				}
				//Efficiency vs Eta
				if(name == "Jets"){
					theHistograms->fill("JetsMatched_vs_eta","JetsMatched_vs_eta;#eta",HISTO_JEta_Config, gen.eta(), 1.);
					theHistograms->fill("JetsMatched_Charge","JetsCharge", 60, -10., 10., gen.charge(), 1.);
				}
				else 
					theHistograms->fill(name+"Matched_vs_eta", name+"Matched_vs_eta;#eta", HISTO_Eta_Config, gen.eta(),1);
				//Efficiency vs Phi
				theHistograms->fill(name+"Matched_vs_phi",name+"Matched_vs_phi;#phi [rad]", HISTO_Phi_Config, gen.phi(),1);
				//Efficiency vs pt
				theHistograms->fill(name+"Matched_vs_pt",name+"Matched_vs_pt;pt [GeV/c]", HISTO_Pt_Config, gen.pt(),1);
				
				//Let's analyze the resolution too!
				resolutionAnalysis(*cand, gen, name);
			}/*
			else if(name == "Jets"){
				cout<<"\n\tgenLeptons:\n";
				foreach(Particle& genP, *genLeptons){
					cout<<genP<<"\n";
				}
				cout<<"\n\tgenJets:\n";
				foreach(Particle& genJ, *genGroup){
					cout<<genJ<<"\n";
				}
				cout<<"\n\tjets:\n"<<*cand;
				foreach(Particle& recP, *recGroup){
					cout<<recP<<"\n";
				}
				for(int i=0; i<50; i++) cout<<"-";
			}*/
		}
	}
}


//	Helper functions
template <class P>
phys::Particle* WWosAnalyzer::findMatchingParticle(const phys::Particle& rec, std::vector<P>* candidates){
	if(candidates->size() == 0) return nullptr;
	int minPos = 0;
	double deltaRMin = physmath::deltaR(rec, candidates->at(0));
	double temp = 999.;
	for(uint i = 1; i < candidates->size(); i++){
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
	normalizeEta(name);	//problem with jets
	normalizePhi(name);
	normalizePt(name);
	
	//Resolution

	theHistograms->get("(#Delta 1/pt)/(1/pt) "+name)->GetXaxis()->SetTitle("");
}

void WWosAnalyzer::normalizeEta(std::string name){
	TH1* matchEta = theHistograms->get(name+"Matched_vs_eta");
	if(matchEta == nullptr){
		cout<<"\""<<name<<"Matched_vs_eta"<<"\" not found\n";
		return;
	}
	theHistograms->clone(name+"Matched_vs_eta", name+"Efficiency_vs_eta");
	TH1* effEta = theHistograms->get(name+"Efficiency_vs_eta");
	if(effEta == nullptr){
		cout<<"\""<<name<<"Efficiency_vs_eta"<<"\" not found\n";
		return;
	}
	/*TGraphAsymmErrors* a = myTGraphAsymmErrors(theHistograms->get(name+"Efficiency_vs_eta"), theHistograms->get("gen"+name+"_eta"));
	a->Draw("AP");*/
	effEta->Divide(theHistograms->get("gen"+name+"_eta"));
	effEta->SetTitle((name+"Efficiency_vs_#eta").string::c_str());
	
	//There's no overload of SetTitle(const char*) with SetTitle(std::string)
	effEta->GetXaxis()->SetTitle("#eta");
}

void WWosAnalyzer::normalizePhi(std::string name){
	theHistograms->clone(name+"Matched_vs_phi", name+"Efficiency_vs_phi");
	TH1* effPhi = theHistograms->get(name+"Efficiency_vs_phi");
	if(effPhi == nullptr){
		cout<<"\""<<name<<"Efficiency_vs_phi"<<"\" not found\n";
	}
	effPhi->Divide(theHistograms->get("gen"+name+"_phi"));
	effPhi->SetTitle((name+"Efficiency_vs_#phi").string::c_str());
	effPhi->GetXaxis()->SetTitle("#phi [rad]");
}

void WWosAnalyzer::normalizePt(std::string name){
	theHistograms->clone(name+"Matched_vs_pt", name+"Efficiency_vs_pt");
	TH1* effPt = theHistograms->get(name+"Efficiency_vs_pt");
	if(effPt == nullptr){
		cout<<"\""<<name<<"Efficiency_vs_pt"<<"\" not found\n";
	}
	effPt->Divide(theHistograms->get("gen"+name+"_pt"));
	effPt->SetTitle((name+"Efficiency_vs_pt").string::c_str());
	effPt->GetXaxis()->SetTitle("pt [GeV/c]");
}


template <class P, class T>
void WWosAnalyzer::resolutionAnalysis(const T& rec, const P& gen, std::string name){
	
	//#Delta(1/pt) is gaussian
	float delta1Pt = (float)((1/rec.pt()-1/gen.pt())*gen.pt());//(D(1/p))/(1/p) = D(1/p) * p
	float deltaEN = (float)((rec.e()-gen.e())/gen.e());
	
		theHistograms->fill("(#Delta 1/pt)/(1/pt) "+name,"(#Delta 1/pt)/(1/pt) "+name, 100,-0.1,0.1, delta1Pt, 1);
		theHistograms->fill("(#Delta E)/E "+name,"(#Delta E)/E "+name, 100,-0.1,0.1, deltaEN, 1);
}

void WWosAnalyzer::fitResolutionPt(std::string name){
	TH1* thisPlot = theHistograms->get("(#Delta 1/pt)/(1/pt) "+name);
	TAxis* thisXAxis = thisPlot->GetXaxis();
	Double_t xm = thisXAxis->GetBinLowEdge(1);
	Int_t xup = thisXAxis->GetNbins();
	Double_t xM = thisXAxis->GetBinLowEdge(xup);
	TF1* ptInv = new TF1("ptInv","[0]*exp(-(x*x)/(2*[1]))+[2]",xm/3.,xM/3.);
	if(name == "Electrons"){
		ptInv->SetParameters(5000., 0.007, 50.);
		ptInv->SetParLimits(1,0.0001,0.001);
		ptInv->SetParLimits(2,0.,1000.);
	}
	if(name == "Muons"){
		ptInv->SetParameters(5000., 0.1, 50.);
		ptInv->SetParLimits(1,0.0001,0.001);
		ptInv->SetParLimits(2,0.,1000.);
	}
	ptInv->SetLineColor(2);
	thisPlot->Fit(ptInv,"R");	//R->Fit in function's range
	getFitInfo(ptInv);
}

void WWosAnalyzer::fitResolutionE(std::string name){
	TH1* thisPlot = theHistograms->get("(#Delta E)/E "+name);
	TAxis* thisXAxis = thisPlot->GetXaxis();
	Double_t xm = thisXAxis->GetBinLowEdge(1);
	Int_t xup = thisXAxis->GetNbins();
	Double_t xM = thisXAxis->GetBinLowEdge(xup);
	TF1* funcE = new TF1("funcE","[0]*exp(-(x*x)/(2*[1]))+[2]",xm/3.,xM/3.);
	if(name == "Electrons"){
		funcE->SetParameters(5000., 0.007, 50.);
		funcE->SetParLimits(0,0.,100000.);
		funcE->SetParLimits(1,0.0001,0.01);
		funcE->SetParLimits(2,0.,2000.);
	}
	if(name == "Muons"){
		funcE->SetParameters(5000., 0.1, 50.);
		funcE->SetParLimits(0,0.,100000.);
		funcE->SetParLimits(1,0.0001,0.01);
		funcE->SetParLimits(2,0.,2000.);
	}
	funcE->SetLineColor(3);
	thisPlot->Fit(funcE,"R");
	getFitInfo(funcE);
}
#endif	//DO_EFFICIENCY_ANALYSIS


// Random statistics 	------------------------------------------------------------------------

void WWosAnalyzer::fillBasicPlots(){
  theHistograms->fill<TH1I>("nvtx"     , "Number of vertices" , 100, 0, 100, nvtx             , theWeight);
  theHistograms->fill      ("met"      , "Missing energy"     , 200, 0, 800, met->pt()        , theWeight);
  theHistograms->fill<TH1I>("nmuons"    ,"Number of muons"    ,  10, 0, 10 , muons->size()    , theWeight);
  theHistograms->fill<TH1I>("nelectrons","Number of electrons",  10, 0, 10 , electrons->size(), theWeight);
	
	theHistograms->get("met")->GetXaxis()->SetTitle("[GeV/c]");
  /*foreach(const phys::Lepton& mu , *muons)     fillLeptonPlots  ("mu",  mu  );
  foreach(const phys::Electron& e, *electrons) fillElectronPlots("e" ,  e   );
  foreach(const phys::Jet& jet   , *jets)      fillJetPlots     ("jet", jet );*/
}



void WWosAnalyzer::fillParticlePlots(const std::string &type, const phys::Particle & lepton){
  theHistograms->fill(type+"_pt" ,   "p_{T} spectrum",HISTO_Pt_Config	,lepton.pt(), 1/*theWeight*/);
  if(type == "genJets" || type == "allJets") 
  	theHistograms->fill(type+"_eta",	"#eta spectrum",HISTO_JEta_Config	,lepton.eta(), 1/*theWeight*/);
  else 
  	theHistograms->fill(type+"_eta", "#eta spectrum",HISTO_Eta_Config	,lepton.eta(), 1/*theWeight*/);
  theHistograms->fill(type+"_phi", 	"#phi spectrum",HISTO_Phi_Config	,lepton.phi(), 1/*theWeight*/);
  //theHistograms->fill(type+"_charge", "charge"    ,	50,		-25 ,  25 ,lepton.charge(), theWeight);
  
  //theHistograms->get(type+"_pt")->GetXaxis()->SetTitle("[GeV/c]");
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
			
			theHistograms->fill("ptPromptParticle", "p_t (prompt)", 25, 0, 50, par.pt());
			//theHistograms->fill("massPromptParticle", "m (prompt)", 1000, 0, 100, par.p4().M());
		}
	/*	if(
				(abs(gen.id()) != 11 && abs(gen.id()) != 13) || 
				(!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || 
				!(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))
			) 				
			continue;
	*/
	//	cout << "id: " << gen.id() << " pt: " << gen.pt() << endl;
	theHistograms->get("ptPromptParticle")->GetXaxis()->SetTitle("[GeV/c]");
	//theHistograms->get("massPromptParticle")->GetXaxis()->SetTitle("[GeV/c^2]");
	//theHistograms->fill("ptAllAnalyzedParticle", "p_t", 25, 0, 50, par.pt());
	//theHistograms->get("ptAllAnalyzedParticle")->GetXaxis()->SetTitle("[GeV]");
	//theHistograms->fill("ParticlesIDs", "ParticlesIDs", 40, -20, 20, par.id());
}
#endif
#ifdef DO_STATISTICS_ON_EVENTS
void WWosAnalyzer::tempStatisticEvents(){
	theHistograms->fill("genParticlesPerEvent", "genParticlesPerEvent", 40, 0, 40, particleCounter);
	theHistograms->fill("genElectronsPerEvent", "genElectronsPerEvent", 25, 0, 25, eCounter);
	theHistograms->fill("genMuonsPerEvent", "genMuonsPerEvent", 25, 0, 25, mCounter);
	
	theHistograms->fill("PromptGenParticlesPerEvent", "PromptGenParticlesPerEvent", 15, 0, 15, promptCounter);
	theHistograms->fill("PromptGenElectronsPerEvent", "PromptGenElectronsPerEvent", 10, 0, 10, peCounter);
	theHistograms->fill("PromptMuonsPerEvent", "PromptMuonsPerEvent", 10, 0, 10, pmCounter);
	
	theHistograms->fill("ElectronsPerEvent", "ElectronsPerEvent", 20, 0, 20, electrons->size());
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
	
	TH1* ePlot = theHistograms->get("genElectronsPerEvent");
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

/*void confidence_binomial_clopper_pearson(int n, int k, double &xlow, double &xhigh, double level=.68540158589942957) {
   double alpha = (1.0 - level)/2;
   xlow = (k == 0) ? 0.0 : ROOT::Math::beta_quantile(alpha, k, n-k + 1.0);
   xhigh = (k == n) ? 1.0 : ROOT::Math::beta_quantile(1.0 - alpha, k + 1.0, n-k);
}*/

/* temp stuff from analyze()
bool eeSignal = false;
  bool purEE = false;			//Only ee
  bool mmSignal = false;
  bool purMM = false;			//Only mm
  bool emSignal = false;
  bool purEM = false;			//Only em
  
  if(genCleanedJets->size() >= 2){
  	//electrons
  	if(genElectrons->size() >= 2){
			if(
					fabs(genCleanedJets->at(0).eta()) > 2.4 && 
					fabs(genCleanedJets->at(1).eta()) > 2.4 &&
					genElectrons->at(0).e() > 20. && 
					met->e() > 15. &&
					genElectrons->at(1).e() > 15. && 
					genElectrons->at(0).charge() * genElectrons->at(1).charge() == -1
				)
				eeSignal = true;
		}
		else eeSignal = false;
		//muons
		if(genMuons->size() >= 2){
			if(
					fabs(genCleanedJets->at(0).eta()) > 2.4 && 
					fabs(genCleanedJets->at(1).eta()) > 2.4 &&
					genMuons->at(0).e() > 20. && 
					met->e() > 15. &&
					genMuons->at(1).e() > 15. && 
					genMuons->at(0).charge() * genMuons->at(1).charge() == -1
				) 
				mmSignal = true;
		}
		else mmSignal = false;
		//mixed
		if((genElectrons->size() >= 1) && (genMuons->size() >= 1)){
			if(
					fabs(genCleanedJets->at(0).eta()) > 2.4 && 
					fabs(genCleanedJets->at(1).eta()) > 2.4 &&
					genElectrons->at(0).e() > 20. &&
					genMuons->at(0).e() > 15. &&
					met->e() > 15. &&
					genMuons->at(0).charge() * genElectrons->at(0).charge() == -1
				)
				emSignal = true;
		}
		else emSignal = false;
	}
	
	if(eeSignal && !mmSignal && !emSignal){ 
		eeEvents++;
		purEE = true;
	}
	if(mmSignal && !eeSignal && !emSignal){ 
		mmEvents++;
		purMM = true;
	}
	if(emSignal && !eeSignal && !mmSignal){ 
		emEvents++;
		purEM = true;
	}
	if(eeSignal || mmSignal || emSignal){
		if(purEE || purMM || purEM) passingSelection++;
		else multiSignEvents++;
	}
*/


 
