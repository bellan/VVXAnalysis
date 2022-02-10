#include "VVXAnalysis/TreeAnalysis/interface/VVGammaAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include "VVXAnalysis/Commons/interface/GenVBHelper.h"

#include "TTree.h"

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

#define CUT_LAYOUT 6,0.5,6.5
#define DEBUG

void VVGammaAnalyzer::begin(){
	cout<<'\n';
	for(char i=0; i<25; ++i) cout<<'-';
	cout<<" Start of VVGammaAnalyzer ";
	for(char i=0; i<25; ++i) cout<<'-';
	cout<<'\n';
	/*
	TFile* photonSFfile = TFile::Open("../data/2018_PhotonsMedium.root", "r");
	if(! (photonSFfile && photonSFfile->IsOpen()) ){
		cout<<"Warning: photon SF file not found\n";
	} else{
		photonSFhist = (TH2F*) photonSFfile->Get("EGamma_SF2D")->Clone();
		photonSFfile->Close();
	}
	delete photonSFfile;*/
	
	std::string spaces( ceil(log10(tree()->GetEntries())), ' ' );
	cout<<"Analyzed:\t"<<spaces<<'/'<<tree()->GetEntries()<<std::flush;
	//cout<<'\n'; //TEMP
	
	return;
}


Int_t VVGammaAnalyzer::cut() {
  ++evtN_; totEvtW_ += theWeight;
	cout<<"\r\t\t"<<evtN_;  //TEMP
	
	//cout<<"theWeight:"<<theWeight<<'\n';
	
	theHistograms->fill("AAA cuts w", "Cuts weighted", CUT_LAYOUT, 1, theWeight);
	theHistograms->fill("AAA cuts u", "Cuts unweighted", CUT_LAYOUT, 1);
	
	// Cleanup
	goodPhotons_->clear();
	
	bool haveZVlep = false;
	bool have2l2j = false;
	bool haveGoodPhoton = false;
	
	theHistograms->fill("C: n leptons", "n leptons", 7,-0.5,6.5, electrons->size()+muons->size(), theWeight);
	
	
	// ----- BASELINE SELECTION -----
	// ----- Cut2: require at least a ZZ or a WZ candidate
	if( (ZZ && ZZ->pt() > 1.) || (ZW && ZW->pt() > 1.) )
		haveZVlep = true;
	else
		haveZVlep = false;
	
	
	if(haveZVlep){
		theHistograms->fill("AAA cuts w", "Cuts weighted", CUT_LAYOUT, 2, theWeight);
		theHistograms->fill("AAA cuts u", "Cuts unweighted", CUT_LAYOUT, 2);
	}
	
	have2l2j = (muons->size()+electrons->size()==2) && (jets->size()==2 || jetsAK8->size()>=1);
	if(have2l2j){
		theHistograms->fill("AAA cuts w", "Cuts weighted", CUT_LAYOUT, 3, theWeight);
		theHistograms->fill("AAA cuts u", "Cuts unweighted", CUT_LAYOUT, 3);
	}
	
	// ----- Cut1: Require at least 1 loose photon with pt > 20 GeV
	for(auto ph : *photons){
		//ID and electron veto
		if(ph.hasPixelSeed() || !ph.passElectronVeto()) continue;
		if(! ph.cutBasedIDLoose()) continue;
		
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
		
		goodPhotons_->push_back(ph);
	}
	haveGoodPhoton = goodPhotons_->size() >= 1;
	
	theHistograms->fill("C: n good ph", "n good #gamma | ZZ/WZ exists", 7,-0.5,6.5, goodPhotons_->size(), theWeight);
	
	if(haveGoodPhoton){
		theHistograms->fill("AAA cuts w", "Cuts weighted", CUT_LAYOUT, 4, theWeight);
		theHistograms->fill("AAA cuts u", "Cuts unweighted", CUT_LAYOUT, 4);
	}
	
	if(!haveZVlep) 
		return -1;
	else if(!haveGoodPhoton) 
		return -1;
	else
		return 1;
}

void VVGammaAnalyzer::analyze(){
	++analyzedN_; analyzedW_ += theWeight;
	theHistograms->fill("AAA cuts w", "Cuts weighted", CUT_LAYOUT, 5, theWeight);
	theHistograms->fill("AAA cuts u", "Cuts unweighted", CUT_LAYOUT, 5);
	
	for(auto ph : *goodPhotons_){
		theHistograms->fill("A: gamma loose", "loose photons;p_{t} [GeV/c];|#eta|", pt_bins, aeta_bins, ph.pt(), fabs(ph.eta()), theWeight);
		if(ph.cutBasedIDMedium()){
			theHistograms->fill("A: gamma medium", "medium photons;p_{t} [GeV/c];|#eta|", pt_bins, aeta_bins, ph.pt(), fabs(ph.eta()), theWeight);
			theHistograms->fill("AAA cuts w", "Cuts weighted", CUT_LAYOUT, 6, theWeight);
			theHistograms->fill("AAA cuts u", "Cuts unweighted", CUT_LAYOUT, 6);
			break;
		}
	}
	
	if(theSampleInfo.isMC()){
		genEventSetup();
		
		for(size_t i = 0; i < goodPhotons_->size(); i++){
			Photon& rec = goodPhotons_->at(i);
			bool matched = false;
			for(auto gen: *genPhotons_){
				double dR = deltaR(gen, rec);
				if(dR < 0.1){
					matched = true;
					break;
				}
			}
			if(matched) continue;
			
			theHistograms->fill("FAKE G: loose", "Unmatched loose photons", pt_bins, aeta_bins, rec.pt(), fabs(rec.eta()), theWeight);
			if(rec.cutBasedIDMedium()){
				theHistograms->fill("FAKE G: medium", "Unmatched medium photons", pt_bins, aeta_bins, rec.pt(), fabs(rec.eta()), theWeight);
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
	for(auto ph: *goodPhotons_){
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
	Photon& thePhoton = goodPhotons_->front();
	
	theHistograms->fill("photon_pt", "p_{t}^{#gamma};GeV/c", 50,0.,200., thePhoton.pt(), theWeight);
	if((region_==CR2P2F || region_==CR3P1F || region_==SR4P) && ZZ && ZZ->pt() > 1.){
		theHistograms->fill("ZZ_mass", "m_{4l};GeV/c^{2}", 25,0.,500., ZZ->mass(), theWeight);
		theHistograms->fill("Z0_mass", "m_{Z0};GeV/c^{2}", 35,55.,125., ZZ->first().mass(), theWeight);
		theHistograms->fill("Z1_mass", "m_{Z1};GeV/c^{2}", 35,55.,125., ZZ->second().mass(), theWeight);
		theHistograms->fill("Z0_l0_pt", "m_{l00};GeV/c", 20,0.,400., ZZ->first().daughter(0).pt(), theWeight);
		theHistograms->fill("Z0_l1_pt", "m_{l00};GeV/c", 20,0.,400., ZZ->first().daughter(1).pt(), theWeight);
		theHistograms->fill("Z1_l0_pt", "m_{l00};GeV/c", 20,0.,400., ZZ->second().daughter(0).pt(), theWeight);
		theHistograms->fill("Z1_l1_pt", "m_{l00};GeV/c", 20,0.,400., ZZ->second().daughter(1).pt(), theWeight);
		// theHistograms->fill("ZZG mass", "m_{ZZ#gamma};GeV/c", 50,0.,500., (ZZ->p4()+thePhoton.p4()).M(), theWeight);
		// theHistograms->fill("Z0G mass", "m_{Z0#gamma};GeV/c", 50,0.,500., (ZZ->first().p4()+thePhoton.p4()).M(), theWeight);
		// theHistograms->fill("Z1G mass", "m_{Z1#gamma};GeV/c", 50,0.,500., (ZZ->second().p4()+thePhoton.p4()).M(), theWeight);
	}
	if( ((region_>=CR110 && region_<=CR000) || region_==SR3P ) && ZW && ZW->pt() > 1.){
		theHistograms->fill("ZW_massT", "m_{T,3l};GeV/c^{2}", 25,0.,500., ZW->p4().Mt(), theWeight);
		theHistograms->fill("Z_mass", "m_{Z};GeV/c^{2}", 35,55.,125., ZW->first().mass(), theWeight);
		theHistograms->fill("W_massT", "m_{T,W};GeV/c^{2}", 35,55.,125., ZW->second().p4().Mt(), theWeight);
		theHistograms->fill("Z_l0_pt", "m_{l00};GeV/c", 20,0.,400., ZW->first().daughter(0).pt(), theWeight);
		theHistograms->fill("Z_l1_pt", "m_{l00};GeV/c", 20,0.,400., ZW->first().daughter(1).pt(), theWeight);
		theHistograms->fill("W_l_pt", "m_{l00};GeV/c", 20,0.,400., ZW->second().daughter(0).pt(), theWeight);
		theHistograms->fill("W_MET_pt", "m_{l00};GeV/c", 20,0.,400., ZW->second().daughter(1).pt(), theWeight);
	}
	if(region_==CRLFR && ZL && ZL->first.pt() > 1.){
		theHistograms->fill("ZL_mass", "m_{3l};GeV/c^{2}", 25,0.,500., (ZL->first.p4()+ZL->second.p4()).M(), theWeight);
		theHistograms->fill("Z_mass", "m_{Z};GeV/c^{2}", 35,55.,125., ZL->first.mass(), theWeight);
		theHistograms->fill("Z_l0_pt", "m_{l00};GeV/c", 20,0.,400., ZL->first.daughter(0).pt(), theWeight);
		theHistograms->fill("Z_l1_pt", "m_{l00};GeV/c", 20,0.,400., ZL->first.daughter(1).pt(), theWeight);
		theHistograms->fill("L_pt", "m_{l00};GeV/c", 20,0.,400., ZL->second.pt(), theWeight);
	}
	return; //TEMP
	
	double leadElpt = electrons->size() > 0 ? electrons->at(0).pt() : 0.;
	double leadMupt = muons->size()     > 0 ? muons->at(0).pt()     : 0.;
	
	// electrons
	for(auto e : *electrons){
		theHistograms->fill("B: pt all e", "p_{t} all e; p_{t}", 50,0.,200., e.pt(), theWeight);
		theHistograms->fill("B: eta all e", "#eta all e; p_{t}", 25,-2.5,2.5, e.eta(), theWeight);
	}
	if(electrons->size() > 0){
		theHistograms->fill("B: pt e1", "p_{t} e1; p_{t}", 50,0.,200., electrons->at(0).pt(), theWeight);
		theHistograms->fill("B: eta e1", "#eta e1; p_{t}", 25,-2.5,2.5, electrons->at(0).eta(), theWeight);
		//theHistograms->fill("B: dR(e1, a1)", "#DeltaR(e1, #gamma1)", 50,0.,5., deltaR(electrons->at(0), goodPhotons_->at(0)), theWeight);
	}
	if(electrons->size() > 1){
		theHistograms->fill("B: pt e2", "p_{t} e2; p_{t}", 50,0.,200., electrons->at(1).pt(), theWeight);
		theHistograms->fill("B: eta e2", "#eta e2; p_{t}", 25,-2.5,2.5, electrons->at(1).eta(), theWeight);
	}
	
	// muons
	for(auto mu : *muons){
		theHistograms->fill("B: pt all mu", "p_{t} all #mu; p_{t}", 50,0.,200., mu.pt(), theWeight);
		theHistograms->fill("B: eta all mu", "#eta all #mu; eta", 25,-2.5,2.5, mu.eta(), theWeight);
	}
	if(muons->size() > 0){
		theHistograms->fill("B: pt mu1", "p_{t} #mu1; p_{t}", 50,0.,200., muons->at(0).pt(), theWeight);
		theHistograms->fill("B: eta mu1", "#eta #mu1; p_{t}", 25,-2.5,2.5, muons->at(0).eta(), theWeight);
		//theHistograms->fill("B: dR(mu1, a1)", "#DeltaR(#mu1, #gamma1)", 50,0.,5., deltaR(muons->at(0), goodPhotons_->at(0)), theWeight);
	}
	if(muons->size() > 1){
		theHistograms->fill("B: pt mu2", "p_{t} #mu2; p_{t}", 50,0.,200., muons->at(1).pt(), theWeight);
		theHistograms->fill("B: eta mu2", "#eta #mu2; p_{t}", 25,-2.5,2.5, muons->at(1).eta(), theWeight);
	}
	
	/*
  theHistograms->fill("nZtoChLep"    , "Number of Z->ll per event" , 7,0,7, genVBHelper_.ZtoChLep().size());
  theHistograms->fill("nZtoNeutrinos", "Number of Z->nn per event" , 7,0,7, genVBHelper_.ZtoNeutrinos().size());
  theHistograms->fill("nWtoLep"      , "Number of W->lnu per event", 7,0,7, genVBHelper_.WtoLep().size());
  theHistograms->fill("nZtoQ"        , "Number of Z->qq per event" , 7,0,7, genVBHelper_.ZtoQ().size());
  theHistograms->fill("nWtoQ"        , "Number of W->qq' per event", 7,0,7, genVBHelper_.WtoQ().size());
  */
  
  //theHistograms->fill("nPhotons"     , "Number of photons' per event", 7,0,7, photons->size());

  //int nVBs = genVBHelper_.ZtoChLep().size() + genVBHelper_.ZtoNeutrinos().size() + genVBHelper_.WtoLep().size() + genVBHelper_.ZtoQ().size() + genVBHelper_.WtoQ().size();
  //theHistograms->fill("nVBs", "Number of VB per event", 7,0,7, nVBs);
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
	TH1* cuts = theHistograms->get("AAA cuts w");
	if(cuts){
		TAxis* x = cuts->GetXaxis();
		x->SetBinLabel(1, "All");
		x->SetBinLabel(2, "ZZ || ZW");
		x->SetBinLabel(3, "2l2j || 2l1J");
		x->SetBinLabel(4, "#gamma loose");
		x->SetBinLabel(5, "analyzed");
		x->SetBinLabel(6, "#gamma medium");
	}
	TH1* cuts_u = theHistograms->get("AAA cuts u");
	if(cuts_u){
		TAxis* x = cuts_u->GetXaxis();
		x->SetBinLabel(1, "All");
		x->SetBinLabel(2, "ZZ || ZW");
		x->SetBinLabel(3, "2l2j || 2l1J");
		x->SetBinLabel(4, "#gamma loose");
		x->SetBinLabel(5, "analyzed");
		x->SetBinLabel(6, "#gamma medium");
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
	theHistograms->fill("GEN genZlepCandidates_", "# genZlepCandidates_", 4,-0.5,3.5, genZlepCandidates_->size());
	theHistograms->fill("GEN genWlepCandidates_", "# genWlepCandidates_", 4,-0.5,3.5, genWlepCandidates_->size());
	theHistograms->fill("GEN genZhadCandidates_", "# genZhadCandidates_", 4,-0.5,3.5, genZhadCandidates_->size());
	theHistograms->fill("GEN genWhadCandidates_", "# genWhadCandidates_", 4,-0.5,3.5, genWhadCandidates_->size());
	
	for(auto v : *genZlepCandidates_)
		theHistograms->fill("GEN genZlepCandidates_ mass", "mass genZlepCandidates_", 35.,50.,120., v.mass());
	for(auto v : *genWlepCandidates_)
		theHistograms->fill("GEN genWlepCandidates_ mass", "mass genWlepCandidates_", 35.,50.,120., v.mass());
	for(auto v : *genZhadCandidates_)
		theHistograms->fill("GEN genZhadCandidates_ mass", "mass genZhadCandidates_", 35.,50.,120., v.mass());
	for(auto v : *genWhadCandidates_)
		theHistograms->fill("GEN genWhadCandidates_ mass", "mass genWhadCandidates_", 35.,50.,120., v.mass());
}


void VVGammaAnalyzer::effPhotons(){
	// Photon gen-reco matching
	vector<Photon> goodPhotonsC(*goodPhotons_);
	
	for(auto gPh : *genPhotons_){
		if(gPh.pt() < 20.) continue;
		float ph_aeta = fabs(gPh.eta());
		if(ph_aeta > 2.4) continue;
		//if(ph_aeta > 1.4442 && ph_aeta < 1.566) continue;
		
		theHistograms->fill("eff: den G pt", "p_{t} gen #gamma", 25,0.,250., gPh.pt(), theWeight);
		//theHistograms->fill("eff: den G eta", "#eta gen #gamma",25,-2.5,2.5, gPh.eta(),theWeight);
		theHistograms->fill("eff: den G eta","#eta gen #gamma", eta_bins, gPh.eta(),theWeight);
		
		if(goodPhotonsC.size() == 0 ) continue;
		
		std::sort(goodPhotonsC.begin(), goodPhotonsC.end(), DeltaRComparator(gPh));
		Photon& rPh = goodPhotonsC.front();
		if(deltaR(rPh, gPh) < 0.1){
			theHistograms->fill("res: dR gen rec","#DeltaR(#gamma_{GEN}, #gamma_{REC})", 20,0.,0.1, deltaR(rPh,gPh), theWeight);
			theHistograms->fill("eff: num G pt","p_{t} gen #gamma", 25,0.,250., gPh.pt(), theWeight);
			theHistograms->fill("eff: num G eta","#eta gen #gamma",eta_bins, gPh.eta(),theWeight);
		}
	}
	
	
	unsigned int tightPh = 0;
	for(auto ph: *goodPhotons_){
		theHistograms->fill("A: goodPh pt","p_{t} loose #gamma", 25,0.,250., ph.pt(), theWeight);
		if(ph.cutBasedIDTight()){
			tightPh++;
			theHistograms->fill("A: goodPh tight pt","p_{t} tight #gamma", 25,0.,250., ph.pt(), theWeight);
			theHistograms->fill("A: goodPh tight eta","#eta tight #gamma",25,-2.5,2.5, ph.eta(),theWeight);
		}
	}
	
	theHistograms->fill("A: tightG num", "# tight #gamma", 5,-0.5,4.5, tightPh, theWeight);
	
	if(tightPh > 0){
		theHistograms->fill("AAA cuts w", "Cuts weighted", CUT_LAYOUT, 6, theWeight);
		theHistograms->fill("AAA cuts u", "Cuts unweighted", CUT_LAYOUT, 6);
	}
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
