#include "VVXAnalysis/TreeAnalysis/interface/VVGammaAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include "VVXAnalysis/Commons/interface/GenVBHelper.h"

#include "TTree.h"

#include <map>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;
using std::vector;
using namespace colour;

using namespace phys;
using namespace physmath;

#define CUT_LAYOUT 7,0.5,7.5
#define DEBUG

void VVGammaAnalyzer::begin(){
	cout<<'\n';
	for(char i=0; i<25; ++i) cout<<'-';
	cout<<" Start of VVGammaAnalyzer ";
	for(char i=0; i<25; ++i) cout<<'-';
	cout<<'\n';
	
	initPhotonSF();
	
	std::string spaces( ceil(log10(tree()->GetEntries())), ' ' );
	cout<<"Analyzed:\t"<<spaces<<'/'<<tree()->GetEntries()<<std::flush;
	//cout<<'\n'; //TEMP
	
	return;
}


Int_t VVGammaAnalyzer::cut() {
  ++evtN_; totEvtW_ += theWeight;
	cout<<"\r\t\t"<<evtN_;  //TEMP
	
	// Cleanup
	kinPhotons_->clear();
	
	bool haveZVlep = false;
	bool have2l2j = false;
	bool haveKinPhoton = false;
	
	if(ZZ){
		regionWord = std::bitset<128>(ZZ->region());
		if(regionWord.test(26))
			return -1; //TEMP exclude signal region when running with -r MC
	}
	
	theHistograms.fill("AAA cuts w", "Cuts weighted", CUT_LAYOUT, 1, theWeight);
	theHistograms.fill("AAA cuts u", "Cuts unweighted", CUT_LAYOUT, 1);
	theHistograms.fill("C: n leptons", "n leptons", 7,-0.5,6.5, electrons->size()+muons->size(), theWeight);
	
	
	// ----- BASELINE SELECTION -----
	// ----- Cut2: require at least a ZZ or a WZ candidate
	if( (ZZ && ZZ->pt() > 1.) || (ZW && ZW->pt() > 1.) )
		haveZVlep = true;
	else
		haveZVlep = false;
	
	
	if(haveZVlep){
		theHistograms.fill("AAA cuts w", "Cuts weighted", CUT_LAYOUT, 2, theWeight);
		theHistograms.fill("AAA cuts u", "Cuts unweighted", CUT_LAYOUT, 2);
	}
	
	have2l2j = (muons->size()+electrons->size()==2) && (jets->size()==2 || jetsAK8->size()>=1);
	if(have2l2j){
		theHistograms.fill("AAA cuts w", "Cuts weighted", CUT_LAYOUT, 3, theWeight);
		theHistograms.fill("AAA cuts u", "Cuts unweighted", CUT_LAYOUT, 3);
	}
	
	// ----- Cut1: Require at least 1 loose photon with pt > 20 GeV
	vector<Photon> kinPh_pixEle;
	for(auto ph : *photons){
		//ID and electron veto
		//if(ph.hasPixelSeed() || !ph.passElectronVeto()) continue;
		//if(! ph.cutBasedIDLoose()) continue;
		
		//Kinematic selection
		if(ph.pt() < 20) continue;
		float ph_aeta = fabs(ph.eta());
		if(ph_aeta > 2.4) continue;
		if(ph_aeta > 1.4442 && ph_aeta < 1.566) continue;
		
		//Electrons and muons matching
		bool match = false;
		for(auto el : *electrons){
			if(deltaR(ph,el) < 0.3){
				match = true;
				break;
			}
		}
		if(match) continue;
		
		for(auto mu : *muons){
			if(deltaR(ph,mu) < 0.3){
				match = true;
				break;
			}
		}
		if(match) continue;
		
		kinPh_pixEle.push_back(ph);
		if(ph.hasPixelSeed() || !ph.passElectronVeto()) continue;
		kinPhotons_->push_back(ph);
	}
	bool eleVeto = false;
	bool pixSeed = false;
	for(auto ph : kinPh_pixEle){
		if(!ph.hasPixelSeed())    pixSeed = true;
		if(ph.passElectronVeto()) eleVeto = true;
	}
	if(pixSeed){
		theHistograms.fill("AAA cuts w", "Cuts weighted", CUT_LAYOUT, 4, theWeight);
		theHistograms.fill("AAA cuts u", "Cuts unweighted", CUT_LAYOUT, 4);
	}
	if(eleVeto){
		theHistograms.fill("AAA cuts w", "Cuts weighted", CUT_LAYOUT, 5, theWeight);
		theHistograms.fill("AAA cuts u", "Cuts unweighted", CUT_LAYOUT, 5);
	}
	
	haveKinPhoton = kinPhotons_->size() >= 1;
	
	theHistograms.fill("C: n good ph", "n good #gamma | ZZ/WZ exists", 7,-0.5,6.5, kinPhotons_->size(), theWeight);
	
	if(haveKinPhoton){
		theHistograms.fill("AAA cuts w", "Cuts weighted", CUT_LAYOUT, 6, theWeight);
		theHistograms.fill("AAA cuts u", "Cuts unweighted", CUT_LAYOUT, 6);
	}
	
	//if(!haveZVlep) //TEMP
	//	return -1;
	if(!haveKinPhoton) 
		return -1;
	else
		return 1;
}

void VVGammaAnalyzer::analyze(){
	++analyzedN_; analyzedW_ += theWeight;
	theHistograms.fill("AAA cuts w", "Cuts weighted", CUT_LAYOUT, 7, theWeight);
	theHistograms.fill("AAA cuts u", "Cuts unweighted", CUT_LAYOUT, 7);
	
	//int nPhLoose = 0;
	for(auto ph : *kinPhotons_){
		theHistograms.fill("A: ph kinem", "photons passing kinematic cut;p_{t} [GeV/c];|#eta|", pt_bins, aeta_bins, ph.pt(), fabs(ph.eta()), theWeight);
		theHistograms.fill("A:u ph kinem", "photons passing kinematic cut;p_{t} [GeV/c];|#eta|", pt_bins, aeta_bins, ph.pt(), fabs(ph.eta()));
		
		int nCuts_Loose = nCutsIDLoose(ph);
		theHistograms.fill("A:u nCuts loose ph", "Number of loose cuts passed by photons", 6,-0.5,5.5, nCuts_Loose);
		
		if(nCuts_Loose >= 3){
			theHistograms.fill("A: ph loose 3", "kin photons #geq3 loose;p_{t} [GeV/c];|#eta|", pt_bins, aeta_bins, ph.pt(), fabs(ph.eta()), theWeight);
			theHistograms.fill("A:u ph loose 3", "kin photons #geq3 loose (unweighted);p_{t} [GeV/c];|#eta|", pt_bins, aeta_bins, ph.pt(), fabs(ph.eta()));
		}
		if(nCuts_Loose >= 4){
			theHistograms.fill("A: ph loose 4", "kin photons #geq4 loose;p_{t} [GeV/c];|#eta|", pt_bins, aeta_bins, ph.pt(), fabs(ph.eta()), theWeight);
			theHistograms.fill("A:u ph loose 4", "kin photons #geq4 loose (unweighted);p_{t} [GeV/c];|#eta|", pt_bins, aeta_bins, ph.pt(), fabs(ph.eta()));
			if(ph.cutBasedIDLoose()){
				//nPhLoose++;
				theHistograms.fill("A: ph loose 5", "kin photons passing loose;p_{t} [GeV/c];|#eta|", pt_bins, aeta_bins, ph.pt(), fabs(ph.eta()), theWeight);
				theHistograms.fill("A:u ph loose 5", "kin photons passing loose (unweighted);p_{t} [GeV/c];|#eta|", pt_bins, aeta_bins, ph.pt(), fabs(ph.eta()));
			}
		}
		/*
		int nCuts_Medium = nCutsIDMedium(ph);
		theHistograms.fill("A:u nCuts medium ph", "Number of medium cuts passed by photons", 6,-0.5,5.5, nCuts_Medium);
		
		if(nCuts_Medium >= 4){
			theHistograms.fill("A: ph medium TOT", "kin photons #geq4 medium;p_{t} [GeV/c];|#eta|", pt_bins, aeta_bins, ph.pt(), fabs(ph.eta()), theWeight);
			theHistograms.fill("A:u ph medium TOT", "kin photons #geq4 medium (unweighted);p_{t} [GeV/c];|#eta|", pt_bins, aeta_bins, ph.pt(), fabs(ph.eta()));
			if(ph.cutBasedIDMedium()){
				theHistograms.fill("A: ph medium PASS", "kin photons passing medium;p_{t} [GeV/c];|#eta|", pt_bins, aeta_bins, ph.pt(), fabs(ph.eta()), theWeight);
				theHistograms.fill("A:u ph medium PASS", "kin photons passing medium (unweighted);p_{t} [GeV/c];|#eta|", pt_bins, aeta_bins, ph.pt(), fabs(ph.eta()));
			}
		}
		*/
	}
	//if(nPhLoose > 0){
	//	theHistograms.fill("AAA cuts w", "Cuts weighted", CUT_LAYOUT, 8, theWeight);
	//	theHistograms.fill("AAA cuts u", "Cuts unweighted", CUT_LAYOUT, 8);
	//}
	
	if(theMCInfo.isMC()){
		genEventSetup();
		
		for(size_t i = 0; i < kinPhotons_->size(); i++){
			Photon& rec = kinPhotons_->at(i);
			bool matched = false;
			for(auto gen: *genPhotons_){
				double dR = deltaR(gen, rec);
				if(dR < 0.1){
					matched = true;
					break;
				}
			}
			if(matched) continue;
			
			theHistograms.fill("FAKE G: loose", "Unmatched loose photons", pt_bins, aeta_bins, rec.pt(), fabs(rec.eta()), theWeight);
			if(rec.cutBasedIDMedium()){
				theHistograms.fill("FAKE G: medium", "Unmatched medium photons", pt_bins, aeta_bins, rec.pt(), fabs(rec.eta()), theWeight);
			}
		}
	}
	
	/*
	for(auto ph: *genPhotons_){
		auto flags = ph.genStatusFlags();
		cout<<"pt: " <<ph.pt()<<"   ";
		cout<<"eta: "<<ph.eta()<<"   ";
		cout<<"isPrompt: "       <<flags.test(phys::GenStatusBit::isPrompt)       <<"   ";
		cout<<"isHardProcess: "  <<flags.test(phys::GenStatusBit::isHardProcess)  <<"   ";
		cout<<"fromHardProcess: "<<flags.test(phys::GenStatusBit::fromHardProcess)<<"   ";
		cout<<"fromHardProcessBeforeFSR: "<<flags.test(phys::GenStatusBit::fromHardProcessBeforeFSR)<<"   ";
		cout<<'\n';
	}
	for(auto ph: *kinPhotons_){
		cout<<"pt: " <<ph.pt()<<"   ";
		cout<<"eta: "<<ph.eta()<<"   ";
		cout<<"loose: "<<ph.cutBasedIDLoose()<<"   ";
		cout<<"medium: "<<ph.cutBasedIDMedium()<<"   ";
		cout<<'\n';
	}
	*/
	// genEventHistos();
	
	// effPhotons();
	
	// SR: medium photon ID - CR: loose && !medium photon ID
	Photon& thePhoton = kinPhotons_->front();
	bool isSR_G = thePhoton.cutBasedIDLoose();  
	std::string region(isSR_G ? "G ID" : "G kin");
	
	theHistograms.fill(region+": G pt", Form("p_{t}^{#gamma} %s;GeV/c", region.c_str()), 50,0.,200., thePhoton.pt(), theWeight);
	if(ZZ && ZZ->pt() > 1.){
		theHistograms.fill(region+": ZZ mass", Form("m_{4l} %s;GeV/c^{2}", region.c_str()), 25,0.,500., ZZ->mass(), theWeight);
		theHistograms.fill(region+": Z0 mass", Form("m_{Z0} %s;GeV/c^{2}", region.c_str()), 35,55.,125., ZZ->first().mass(), theWeight);
		theHistograms.fill(region+": Z1 mass", Form("m_{Z1} %s;GeV/c^{2}", region.c_str()), 35,55.,125., ZZ->second().mass(), theWeight);
		theHistograms.fill(region+": ZZG mass", Form("m_{ZZ#gamma} %s;GeV/c", region.c_str()), 50,0.,500., (ZZ->p4()+thePhoton.p4()).M(), theWeight);
		theHistograms.fill(region+": Z0G mass", Form("m_{Z0#gamma} %s;GeV/c", region.c_str()), 50,0.,500., (ZZ->first().p4()+thePhoton.p4()).M(), theWeight);
		theHistograms.fill(region+": Z1G mass", Form("m_{Z1#gamma} %s;GeV/c", region.c_str()), 50,0.,500., (ZZ->second().p4()+thePhoton.p4()).M(), theWeight);
	}
	
	double leadElpt = electrons->size() > 0 ? electrons->at(0).pt() : 0.;
	double leadMupt = muons->size()     > 0 ? muons->at(0).pt()     : 0.;
	
	theHistograms.fill(region+": lead L pt", "Leading lepton p_{t};p_{t} [GeV/c]", 20,0.,400., std::max(leadMupt, leadElpt), theWeight);
	
	/*
  theHistograms.fill("nZtoChLep"    , "Number of Z->ll per event" , 7,0,7, genVBHelper_.ZtoChLep().size());
  theHistograms.fill("nZtoNeutrinos", "Number of Z->nn per event" , 7,0,7, genVBHelper_.ZtoNeutrinos().size());
  theHistograms.fill("nWtoLep"      , "Number of W->lnu per event", 7,0,7, genVBHelper_.WtoLep().size());
  theHistograms.fill("nZtoQ"        , "Number of Z->qq per event" , 7,0,7, genVBHelper_.ZtoQ().size());
  theHistograms.fill("nWtoQ"        , "Number of W->qq' per event", 7,0,7, genVBHelper_.WtoQ().size());
  */
  
  //theHistograms.fill("nPhotons"     , "Number of photons' per event", 7,0,7, photons->size());

  //int nVBs = genVBHelper_.ZtoChLep().size() + genVBHelper_.ZtoNeutrinos().size() + genVBHelper_.WtoLep().size() + genVBHelper_.ZtoQ().size() + genVBHelper_.WtoQ().size();
  //theHistograms.fill("nVBs", "Number of VB per event", 7,0,7, nVBs);
}


void VVGammaAnalyzer::end(TFile& fout){
	cout<<'\n';
	
	// Label names
	endNameHistos();
	
	cout<<Form("\nTotal events: %lu (weighted: %.2f)", evtN_, totEvtW_);
	cout<<Form("\nPassing cut:  %lu (weighted: %.2f)", analyzedN_, analyzedW_);
	cout<<Form("\nFraction:     %.1f %% (weighted: %.1f %%)", 100.*analyzedN_/evtN_, 100.*analyzedW_/totEvtW_)<<'\n';	
	
	
	float elapsedSec = (float)(clock()-startTime_)/CLOCKS_PER_SEC;
	int elapsedSecInt = (int)elapsedSec;
	cout<<"\nElapsed Time: "<<elapsedSec<<" s\t\t("<<elapsedSecInt/60<<"\' "<<elapsedSecInt%60<<"\")\n";
	for(char i=0; i<25; ++i) cout<<'-';
	cout<<" End of VZZAnalyzer ";
	for(char i=0; i<25; ++i) cout<<'-';
	cout<<"\n\n";
}


void VVGammaAnalyzer::endNameHistos(){
	std::map<int, std::string> cutsName = {
		{1, "All"},
		{2, "ZZ || ZW"},
		{3, "2l2j || 2l1J"},
		{4, "pixelSeed"},
		{5, "electronVeto"},
		{6, "#gamma kin"},
		{7, "analyzed"}
	};
	TH1* cuts = theHistograms.get("AAA cuts w");
	if(cuts){
		TAxis* x = cuts->GetXaxis();
		for(auto cn : cutsName)
			x->SetBinLabel(cn.first, cn.second.c_str());
	}
	TH1* cuts_u = theHistograms.get("AAA cuts u");
	if(cuts_u){
		TAxis* x = cuts_u->GetXaxis();
		for(auto cn : cutsName)
			x->SetBinLabel(cn.first, cn.second.c_str());
	}
}


void VVGammaAnalyzer::genEventSetup(){
	genQuarks_->clear();
	genChLeptons_->clear();
	genNeutrinos_->clear();
	genPhotons_->clear();
	
	genZlepCandidates_->clear();
	genWlepCandidates_->clear();
	genZhadCandidates_->clear();
	genWhadCandidates_->clear();
	
	genZZ_ = DiBoson<Particle, Particle>();
	genWZ_ = DiBoson<Particle, Particle>();
	
	// Sort gen particles
	for(auto p : *genParticles){
		if(abs(p.id()) < 9)
			genQuarks_->push_back(p);
		else if(abs(p.id()) == 11 || abs(p.id()) == 13){
			genChLeptons_->push_back(p);
		}
		else if(abs(p.id()) == 12 || abs(p.id()) == 14)
			genNeutrinos_->push_back(p);
		else if( p.id() == 22 && p.pt() > 18.)
			genPhotons_->push_back(p);
	}
	
	// Gen W --> l nu
	if(genNeutrinos_->size() > 0 && genChLeptons_->size() > 0){
		for(auto l : *genChLeptons_){
			for(auto v : *genNeutrinos_){
				if( abs(l.id() + v.id()) == 1 ){
					Boson<Particle> Wcand(l,v);
					if(GenWBosonDefinition(Wcand))
						genWlepCandidates_->push_back(Wcand);
				}
			}
		}
	}
	
	// Gen W --> q q'bar
	if(genQuarks_->size() >= 2){
		for(size_t i = 0  ; i < genQuarks_->size(); ++i){
		Particle& q1 = genQuarks_->at(i);
		if(q1.id() > 5) continue;
			for(size_t j = i+1; j < genQuarks_->size(); ++j){
				Particle& q2 = genQuarks_->at(j);
				if(q2.id() > 5) continue;
				
				if( (q1.id() * q2.id() < 0) && ( abs(q1.id()+q2.id()) % 2 ==1 ) ){
					Boson<Particle> Wcand(q1,q2);
					if(GenWBosonDefinition(Wcand))
						genWhadCandidates_->push_back(Wcand);
				}
			}
		}
	}
	
	// Gen Z --> q qbar
	if(genQuarks_->size() >= 2){
		for(size_t i = 0  ; i < genQuarks_->size(); ++i){
		Particle& q1 = genQuarks_->at(i);
		if(q1.id() > 5) continue;
			for(size_t j = i+1; j < genQuarks_->size(); ++j){
				Particle& q2 = genQuarks_->at(j);
				if(q2	.id() > 5) continue;
				
				if( q1.id() + q2.id() == 0 ){
					Boson<Particle> Zcand(q1,q2);
					if(ZBosonDefinition(Zcand))
						genZhadCandidates_->push_back(Zcand);
				}
			}
		}
	}
	
	// Gen Z --> l lbar
	if(genChLeptons_->size() >= 2){
		for(size_t i = 0 ; i < genChLeptons_->size(); ++i){
		Particle& l1 = genChLeptons_->at(i);
			for(size_t j = i+1; j < genChLeptons_->size(); ++j){
				Particle& l2 = genChLeptons_->at(j);
				
				if( l1.id() + l2.id() == 0 ){
					Boson<Particle> Zcand(l1,l2);
					if(ZBosonDefinition(Zcand))
						genZhadCandidates_->push_back(Zcand);
				}
			}
		}
	}
	
	// genZZ --> 4l
	if(genChLeptons_->size() >= 4 && genZlepCandidates_->size() >= 2){
		std::sort(genZlepCandidates_->begin(), genZlepCandidates_->end(), MassComparator(phys::ZMASS));
		Boson<Particle>& Z0 = genZlepCandidates_->front();
		
		// Vector containing the rest of the Zll candidates
		vector<Boson<Particle>> Zll(genZlepCandidates_->begin()+1, genZlepCandidates_->end());
		std::sort(Zll.begin(), Zll.end(), ScalarSumPtComparator());
		Boson<Particle>* pZ1 = nullptr;
		for(size_t i = 0; i < Zll.size(); ++i){
			if(! haveCommonDaughter(Z0, Zll.at(i))){
				pZ1 = &(Zll.at(i));
				break;
			}
		}
		if(pZ1)
			genZZ_ = DiBoson<Particle, Particle>(Z0, *pZ1);
	}
	
	// genZW --> 3l nu
	if(genChLeptons_->size() >= 3 && genZlepCandidates_->size() >= 1 && genWlepCandidates_->size() >= 1){	
		std::sort(genZlepCandidates_->begin(), genZlepCandidates_->end(), MassComparator(phys::ZMASS));
		Boson<Particle>& Z0 = genZlepCandidates_->front();
		
		std::sort(genWlepCandidates_->begin(), genWlepCandidates_->end(), MassComparator(phys::WMASS));
		Boson<Particle>& W0 = genWlepCandidates_->front();
		
		genWZ_ = DiBoson<Particle, Particle>(Z0, W0);
	}
}


void VVGammaAnalyzer::genEventHistos(){
	theHistograms.fill("GEN genZlepCandidates_", "# genZlepCandidates_", 4,-0.5,3.5, genZlepCandidates_->size());
	theHistograms.fill("GEN genWlepCandidates_", "# genWlepCandidates_", 4,-0.5,3.5, genWlepCandidates_->size());
	theHistograms.fill("GEN genZhadCandidates_", "# genZhadCandidates_", 4,-0.5,3.5, genZhadCandidates_->size());
	theHistograms.fill("GEN genWhadCandidates_", "# genWhadCandidates_", 4,-0.5,3.5, genWhadCandidates_->size());
	
	for(auto v : *genZlepCandidates_)
		theHistograms.fill("GEN genZlepCandidates_ mass", "mass genZlepCandidates_", 35.,50.,120., v.mass());
	for(auto v : *genWlepCandidates_)
		theHistograms.fill("GEN genWlepCandidates_ mass", "mass genWlepCandidates_", 35.,50.,120., v.mass());
	for(auto v : *genZhadCandidates_)
		theHistograms.fill("GEN genZhadCandidates_ mass", "mass genZhadCandidates_", 35.,50.,120., v.mass());
	for(auto v : *genWhadCandidates_)
		theHistograms.fill("GEN genWhadCandidates_ mass", "mass genWhadCandidates_", 35.,50.,120., v.mass());
}


void VVGammaAnalyzer::effPhotons(){
	// Photon gen-reco matching
	vector<Photon> goodPhotonsC(*kinPhotons_);
	
	for(auto gPh : *genPhotons_){
		if(gPh.pt() < 20.) continue;
		float ph_aeta = fabs(gPh.eta());
		if(ph_aeta > 2.4) continue;
		//if(ph_aeta > 1.4442 && ph_aeta < 1.566) continue;
		
		theHistograms.fill("eff: den G pt", "p_{t} gen #gamma", 25,0.,250., gPh.pt(), theWeight);
		//theHistograms.fill("eff: den G eta", "#eta gen #gamma",25,-2.5,2.5, gPh.eta(),theWeight);
		theHistograms.fill("eff: den G eta","#eta gen #gamma", eta_bins, gPh.eta(),theWeight);
		
		if(goodPhotonsC.size() == 0 ) continue;
		
		std::sort(goodPhotonsC.begin(), goodPhotonsC.end(), DeltaRComparator(gPh));
		Photon& rPh = goodPhotonsC.front();
		if(deltaR(rPh, gPh) < 0.1){
			theHistograms.fill("res: dR gen rec","#DeltaR(#gamma_{GEN}, #gamma_{REC})", 20,0.,0.1, deltaR(rPh,gPh), theWeight);
			theHistograms.fill("eff: num G pt","p_{t} gen #gamma", 25,0.,250., gPh.pt(), theWeight);
			theHistograms.fill("eff: num G eta","#eta gen #gamma",eta_bins, gPh.eta(),theWeight);
		}
	}
	
	
	unsigned int tightPh = 0;
	for(auto ph: *kinPhotons_){
		theHistograms.fill("A: goodPh pt","p_{t} loose #gamma", 25,0.,250., ph.pt(), theWeight);
		if(ph.cutBasedIDTight()){
			theHistograms.fill("A: goodPh tight pt","p_{t} tight #gamma", 25,0.,250., ph.pt(), theWeight);
			theHistograms.fill("A: goodPh tight eta","#eta tight #gamma",25,-2.5,2.5, ph.eta(),theWeight);
		}
	}
	
	theHistograms.fill("A: tightG num", "# tight #gamma", 5,-0.5,4.5, tightPh, theWeight);
}


void VVGammaAnalyzer::initPhotonSF(){
	if(! theMCInfo.isMC()){
		for(auto wp : {WPCutID::Loose, WPCutID::Medium} )
			photonSFhists_[wp] = nullptr;
		return;
	}
		
	int year = theMCInfo.setup();
	auto curDir = gDirectory;
	// Opening a TFile changes the TDirectory. When closing the file, all the histograms would be deleted
	for(auto wp : {WPCutID::Loose, WPCutID::Medium} ){
		const char* wp_str;
		switch(wp){
			case WPCutID::Loose : wp_str = "Loose" ; break;
			case WPCutID::Medium: wp_str = "Medium"; break;
		}
		auto fname = Form("../data/%d_Photons%s.root", year, wp_str);
		
		TFile* SFfile = TFile::Open(fname, "r");
		if(! (SFfile && SFfile->IsOpen()) ){
			//cout<<"Warning: photon SF file not found\n";
			photonSFhists_[wp] = nullptr;
		} else{
			auto hist = (TH2F*) SFfile->Get("EGamma_SF2D");//->Clone();
			hist->SetDirectory(curDir);
			photonSFhists_[wp] = hist;
			SFfile->Close();
		}
		delete SFfile;
	}
}


double VVGammaAnalyzer::getSF(const phys::Photon& ph, WPCutID wp) const{ 
	auto hist = photonSFhists_.at(wp);
	return ( hist ? hist->GetBinContent( hist->FindBin(ph.eta(), ph.pt()) ) : 1.);
}

double VVGammaAnalyzer::getSFerr(const phys::Photon& ph, WPCutID wp) const{ 
	auto hist = photonSFhists_.at(wp);
	return ( hist ? hist->GetBinError( hist->FindBin(ph.eta(), ph.pt()) ) : 0.);
}


const vector<double> VVGammaAnalyzer::pt_bins(
{20., 35., 50., 100., 200., 500.}
);

const vector<double> VVGammaAnalyzer::eta_bins(
{-2.4, -2., -1.566, -1.4442, -0.8, 0., 0.8, 1.4442, 1.566, 2., 2.4}
//{-2.4, -2.1, -1.8, -1.566, -1.4442, -1.2, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.2, 1.4442, 1.566, 1.8, 2.1, 2.4}
);

const vector<double> VVGammaAnalyzer::aeta_bins(
{0., 0.8, 1.4442, 1.566, 2., 2.4}
//{0., 0.3, 0.6, 0.9, 1.2, 1.4442, 1.566, 1.8, 2.1, 2.4}
);
