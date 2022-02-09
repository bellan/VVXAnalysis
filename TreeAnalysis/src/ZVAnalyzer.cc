#include "VVXAnalysis/TreeAnalysis/interface/ZVAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/DataFormats/interface/Jet.h"
#include <TF1.h>
#include <TH1.h>
#include <TEfficiency.h>
#include <vector>
#include <boost/foreach.hpp>
#include <string>
#include <time.h>
#include "Math/Math.h"
#include "Math/QuantFuncMathCore.h"


#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp>

using namespace boost::assign;

using std::cout;
using std::endl;

using namespace phys;
using namespace physmath;

#define HISTO_deltaR_Config     200,0.,0.5
#define HISTO_Eta_Config             261,-2.6,2.6
#define HISTO_JEta_Config         250,-5.,5.
#define HISTO_Phi_Config             100,-3.15,3.15
#define HISTO_Pt_Config             100,0.,500.

Int_t ZVAnalyzer::cut() {
    totalEvents++;
    if (( (electrons-> size() == 2) xor (muons->size() == 2)) && (jets->size() > 3)) return 1;
    return -1;
}

void ZVAnalyzer::begin(){
    cout<<"\n";
    for(char i=0; i<25; i++) cout<<"-";
    cout<<" \tBegin of ZV\t ";
    for(char i=0; i<25; i++) cout<<"-";
    cout<<"\n";
    startTime = clock();
}

void ZVAnalyzer::GenAnalysis(){
    std::vector<phys::Particle> *genElectrons = new std::vector<phys::Particle>(); //genElectrons in this events
    std::vector<phys::Particle> *genMuons = new std::vector<phys::Particle>(); //genMuons in this events
    std::vector<phys::Particle> *genLeptons = new std::vector<phys::Particle>(); //genLeptons in this events
    std::vector<phys::Particle> *gencleanedjets = new std::vector<phys::Particle>();  //genJets in this events

    
    foreach(const phys::Particle &gen, *genParticles){
        tempStatisticParticles(gen);
        if(/*!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)) ||*/ !(abs(gen.eta())<2.4) || gen.pt()< 25) {
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
                fillParticlePlots("genLeptons", gen);
                genLeptons->push_back(gen);
                break;
            default:
                continue;
        }
    }
    
    foreach(const phys::Particle &jets, *genJets){
        phys::Particle* cand = findMatchingParticle(jets, genLeptons);
        if (cand != 0){
            if(checkMatch(jets, *cand, 0.4) && (abs(cand->id()) >= 11 && abs(cand->id()) <= 16)){
                continue;
            }
            else {
                if(fabs(jets.eta()) < 4.7){
                    fillParticlePlots("genJets",jets);
                    gencleanedjets->push_back(jets);
                }
            }
        }
    }
    
    
    tempStatisticEvents();
    
    std::sort(genElectrons->begin(), genElectrons->end(), PtComparator());    //Descending order
    std::sort(genMuons->begin(), genMuons->end(), PtComparator());
    std::sort(genLeptons->begin(), genLeptons->end(), PtComparator());
    std::sort(gencleanedjets->begin(), gencleanedjets->end(), PtComparator());
    
    
    totalElectrons += genElectrons->size();
    totalMuons += genMuons->size();
    totalJets += gencleanedjets->size();
    totalJets2 += genJets->size();
    
    analyzeEfficiency(genElectrons, electrons, std::string("Electrons"), matchedElectrons, 0.1);
    analyzeEfficiency(genMuons, muons, std::string("Muons"), matchedMuons, 0.1);
    analyzeEfficiency(gencleanedjets, jets, std::string("Jets"),matchedJets, 0.4);
    
    delete genElectrons;
    delete genMuons;
    delete genLeptons;
    delete gencleanedjets;
    
}

void ZVAnalyzer::endGenAnalysis(){
    normalizeHistograms(std::string("Electrons"));
    normalizeHistograms(std::string("Muons"));
    normalizeHistograms(std::string("Jets"));
    doSomeFits(std::string("Electrons"));
    doSomeFits(std::string("Muons"));
    
    cout << "The number of electrons is: " << totalElectrons << "and the number of matched electrons is: " << matchedElectrons << "\n";
    cout <<  "Efficiency: " << (matchedElectrons/totalElectrons)*100 << " %\n";
    
    cout << "The number of muons is: " << totalMuons << "and the number of matched muons is: " << matchedMuons << "\n";
    cout <<  "Efficiency: " << (matchedMuons/totalMuons)*100 << " %\n";
    
    cout << "The number of jets is: " << totalJets << "and the number of matched jets is: " << matchedJets << "\n";
    cout <<  "Efficiency: " << (matchedJets/totalJets)*100 << " %\n";

}

void ZVAnalyzer::analyze(){
    passingSelection++;
    
            switch (electrons->size()) {
                    float deltaEta;
                    float deltaPhi;
                    float deltaR;
                case 2:
                    if (!((electrons->at(0).charge())*(electrons->at(1).charge())<0)){
                        badevents++;
                        return;
                    }
                    else
                        theHistograms->fill("Leading e pt", "Leading e pt", HISTO_Pt_Config, electrons->at(0).pt(), theWeight);
                    theHistograms->fill("Second e pt", "Second e pt", HISTO_Pt_Config, electrons->at(1).pt(), theWeight);
                    deltaEta = electrons->at(0).eta() - electrons->at(1).eta();
                    theHistograms->fill("#Delta #eta ee", "#Delta #eta ee", HISTO_JEta_Config, deltaEta, theWeight);
                    deltaPhi = electrons->at(0).phi() - electrons->at(1).phi();
                    theHistograms->fill("#Delta #phi ee", "#Delta #phi ee", HISTO_Phi_Config, deltaPhi, theWeight);
                    deltaR = physmath::deltaR(electrons->at(0), electrons->at(1));
                    theHistograms->fill("#Delta R ee", "#Delta R ee", 100, 0., 2., deltaR, theWeight);
                    break;
                case 0:
                    if (!((muons->at(0).charge())*(muons->at(1).charge())<0)){
                        badevents++;
                        return;
                    }
                    else
                        theHistograms->fill("Leading m pt", "Leading m pt", HISTO_Pt_Config, muons->at(0).pt(), theWeight);
                    theHistograms->fill("Second m pt", "Second m pt", HISTO_Pt_Config, muons->at(1).pt(), theWeight);
                    deltaEta = muons->at(0).eta() - muons->at(1).eta();
                    theHistograms->fill("#Delta #eta mm", "#Delta #eta mm", HISTO_JEta_Config, deltaEta, theWeight);
                    deltaPhi = muons->at(0).phi() - muons->at(1).phi();
                    theHistograms->fill("#Delta #phi mm", "#Delta #phi mm", HISTO_Phi_Config, deltaPhi, theWeight);
                    deltaR = physmath::deltaR(muons->at(0), muons->at(1));
                    theHistograms->fill("#Delta R mm", "#Delta R mm", 100, 0., 2., deltaR, theWeight);
                    
                default:
                    break;
            }
    
    
    GenAnalysis(); //All the analysis about efficiency and resolution
    
    
    
    
    

    return;
    
}

void ZVAnalyzer::end(TFile &){
    
    endGenAnalysis();
    
    cout << "The number of good events: " << passingSelection << "\tFirst cut efficiency: " << (1.-passingSelection/totalEvents)*100. << "%\n";
    cout << "The number of bad events: " << badevents << "\tFraction of bad events that pass the first selection: " << (badevents/passingSelection)*100. << "%\n";
    cout << "The new number of good events: " << passingSelection-badevents << "\tSecond cut efficiency: " << (1.-(passingSelection-badevents)/totalEvents)*100. << "%\n";
    
    cout<<"\nElapsed Time: "<< (float)(clock()-startTime)/CLOCKS_PER_SEC<<" s\n";
    for(char i=0; i<25; i++) cout<<"-";
    cout<<" \tEnd of ZV\t";
    for(char i=0; i<25; i++) cout<<"-";
    cout<<"\n\n";
}

template<class T, class P>
void ZVAnalyzer::analyzeResolutionpt(const T& generated, const P& reconstructed, std::string name) {
    double ideltapt = 1./generated.pt() - 1./reconstructed.pt();
    theHistograms->fill(name+"pt_Resolution", name+"pt_Resolution", 200, -0.2, 0.2, ideltapt/(1./generated.pt()),1);
    theHistograms->get(name+"pt_Resolution")->GetXaxis()->SetTitle("#Delta(#frac{1}{pt})*pt");
};

template <class T, class P>
void ZVAnalyzer::analyzeResolutionEnergy(const T& generated, const P& reconstructed, std::string name) {
    double deltaE = generated.e() - reconstructed.e();
    theHistograms->fill(name+"Energy_Resolution", name+"Energy_Resolution", 200, -0.2, 0.2, deltaE/(generated.e()),1);
    theHistograms->get(name+"Energy_Resolution")->GetXaxis()->SetTitle("#frac{#DeltaE}{E}");
};

template <class T, class P, typename C>
void ZVAnalyzer::analyzeEfficiency(std::vector<T>* genGroup, std::vector<P>* recGroup, std::string name, C& counter, const float& tolerance){
    foreach(const phys::Particle & gen, *genGroup){
        phys::Particle* cand = findMatchingParticle(gen, recGroup);    //candidate is a reconstructed particle
        if(cand != nullptr){
            counter++;
            if (checkMatch(gen, *cand, tolerance)){
               if(name == "Jets"){
                    theHistograms->fill(name+"Matched_vs_pt",name+"Matched_vs_pt", HISTO_Pt_Config, gen.pt(), 1);
                    theHistograms->fill(name+"Matched_vs_eta",name+"Matched_vs_eta", HISTO_JEta_Config, gen.eta(), 1);
                }
                //Efficiency vs Eta
            theHistograms->fill(name+"Matched_vs_eta",name+"Matched_vs_eta", HISTO_Eta_Config, gen.eta(), 1);
                //Efficiency vs pt
            theHistograms->fill(name+"Matched_vs_pt",name+"Matched_vs_pt",HISTO_Pt_Config, gen.pt(), 1);
            }
        analyzeResolutionpt(gen, *cand, name);
        analyzeResolutionEnergy(gen, *cand, name);
        }
    }
}

template <class P>
phys::Particle* ZVAnalyzer::findMatchingParticle(const phys::Particle& gen, std::vector<P>* candidates){
    if(candidates->size() == 0) return nullptr;
    minPos=0;
    minDeltaR= physmath::deltaR(gen,candidates->at(0));
    double temp=999.;
        for(uint i=1;i<candidates->size();i++){
                temp=physmath::deltaR(gen,candidates->at(i));
                if (temp<minDeltaR){
                    minPos=i;
                    minDeltaR=temp;
                }
        }
    return & (candidates->at(minPos));
}

template <class P, class T>
bool ZVAnalyzer::checkMatch(const /*phys::Particle&*/P& generated, /*const phys::Particle&*/ T& reconstructed,  const float& tolerance){
    return physmath::deltaR(generated, reconstructed) < tolerance;
}


void ZVAnalyzer::normalizeHistograms(std::string name){
    theHistograms->clone(name+"Matched_vs_eta", name+"Efficiency_vs_eta");
    theHistograms->get(name+"Efficiency_vs_eta")->Divide(theHistograms->get("gen"+name+"_eta"));
    theHistograms->get(name+"Efficiency_vs_eta")->SetTitle((name+"Efficiency_vs_eta").std::string::c_str());
    //There's no overload of SetTitle(const char*) with SetTitle(std::string)
    theHistograms->get(name+"Efficiency_vs_eta")->GetXaxis()->SetTitle("#eta");
    
    theHistograms->clone(name+"Matched_vs_pt", name+"Efficiency_vs_pt");
    theHistograms->get(name+"Efficiency_vs_pt")->Divide(theHistograms->get("gen"+name+"_pt"));
    theHistograms->get(name+"Efficiency_vs_pt")->SetTitle((name+"Efficiency_vs_pt").std::string::c_str());
    theHistograms->get(name+"Efficiency_vs_pt")->GetXaxis()->SetTitle("pt [GeV/c]");
    
}

void ZVAnalyzer::fillBasicPlots(){
        theHistograms->fill<TH1I>("nvtx"     , "Number of vertices" , 100, 0, 100, nvtx             , theWeight);
        theHistograms->fill<TH1I>("nmuons"    ,"Number of muons"    ,  10, 0, 10 , muons->size()    , theWeight);
        theHistograms->fill<TH1I>("nelectrons","Number of electrons",  10, 0, 10 , electrons->size(), theWeight);
        
        theHistograms->get("met")->GetXaxis()->SetTitle("[GeV/c]");
    }
    
    
template <class P>
void ZVAnalyzer::fillParticlePlots(const std::string &type, const P& lepton){
    if (type == "genJets") {
        theHistograms->fill(type+"_eta",    "#eta spectrum" ,HISTO_JEta_Config ,lepton.eta()   , 1);
    }
        theHistograms->fill(type+"_pt" ,    "p_{T} spectrum", HISTO_Pt_Config,lepton.pt()    , 1);
        theHistograms->fill(type+"_eta",    "#eta spectrum" , HISTO_Eta_Config ,lepton.eta()   , 1);
        theHistograms->fill(type+"_phi",    "#phi spectrum" ,  HISTO_Phi_Config,lepton.phi()   , 1);
    
        
        theHistograms->get(type+"_pt")->GetXaxis()->SetTitle("[GeV/c]");
    }

void ZVAnalyzer::initStatistics(){
        counter = 0;
        eCounter = 0;
        mCounter = 0;
        
        promptCounter = 0;
        peCounter = 0;
        pmCounter = 0;
    }
    
void ZVAnalyzer::tempStatisticParticles(const phys::Particle &par){
        counter++;
        if(abs(par.id()) == 11) eCounter++;
        if(abs(par.id()) == 13) mCounter++;
        
        if(par.genStatusFlags().test(phys::GenStatusBit::isPrompt)){
            promptCounter++;
            if(abs(par.id()) == 11) peCounter++;
            if(abs(par.id()) == 13) pmCounter++;
            
            theHistograms->fill("ptPromptParticle", "p_t (prompt)", 25, 0, 50, par.pt());
            theHistograms->fill("massPromptParticle", "m (prompt)", 1000, 0, 100, par.p4().M());
        }
        theHistograms->get("ptPromptParticle")->GetXaxis()->SetTitle("[GeV/c]");
        theHistograms->get("massPromptParticle")->GetXaxis()->SetTitle("[GeV/c^2]");
        theHistograms->fill("ptAllAnalyzedParticle", "p_t", 25, 0, 50, par.pt());
        theHistograms->get("ptAllAnalyzedParticle")->GetXaxis()->SetTitle("[GeV]");
        theHistograms->fill("ParticlesIDs", "ParticlesIDs", 40, -20, 20, par.id());
    }
    
void ZVAnalyzer::tempStatisticEvents(){
        theHistograms->fill("ParticlesPerEvent", "ParticlesPerEvent", 40, 0, 40, counter);
        theHistograms->fill("genElectronsPerEvent", "genElectronsPerEvent", 25, 0, 25, eCounter);
        theHistograms->fill("genMuonsPerEvent", "genMuonsPerEvent", 25, 0, 25, mCounter);
        
        theHistograms->fill("PromptParticlesPerEvent", "PromptParticlesPerEvent", 25, 0, 25, promptCounter);
        theHistograms->fill("PromptgenElectronsPerEvent", "PromptgenElectronsPerEvent", 25, 0, 25, peCounter);
        theHistograms->fill("PromptMuonsPerEvent", "PromptMuonsPerEvent", 25, 0, 25, pmCounter);
        
        theHistograms->fill("ElectronsPerEvent", "ElectronsPerEvent", 25, 0, 25, electrons->size());
    }
    
void ZVAnalyzer::doSomeFits(std::string name){
    
    if (name == "Electrons") {
        TF1* func1 = new TF1("func1","[0]*exp(-x*x/(2*[1])) + [2]",-0.05,0.05);
        func1->SetParameters(2000,0.2,100);
        func1->SetParLimits(0,0,10000000);
        func1->SetParLimits(2,0,2000);
        func1->SetLineColor(4);
        theHistograms->get(name+"pt_Resolution")->Fit(func1,"R+");
        
        TF1* func2 = new TF1("func2","[0]*exp(-x*x/(2*[1])) + [2]",-0.05,0.05);
        func2->SetParameters(2000,0.2,100);
        func2->SetParLimits(0,0,10000000);
        func2->SetParLimits(2,0,2000);
        func2->SetLineColor(3);
        theHistograms->get(name+"Energy_Resolution")->Fit(func2, "R+");
    }
    if (name == "Muons"){
        TF1* func3 = new TF1("func3","[0]*exp(-x*x/(2*[1])) + [2]",-0.05,0.05);
        func3->SetParameters(2500,0.2,100);
        func3->SetParLimits(0,0,10000000);
        func3->SetParLimits(2,0,2500);
        func3->SetLineColor(4);
        theHistograms->get(name+"pt_Resolution")->Fit(func3,"R+");
        
        TF1* func4 = new TF1("func4","[0]*exp(-x*x/(2*[1])) + [2]",-0.05,0.05);
        func4->SetParameters(2500,0.2,100);
        func4->SetParLimits(0,0,10000000);
        func4->SetParLimits(2,0,2500);
        func4->SetLineColor(3);
        theHistograms->get(name+"Energy_Resolution")->Fit(func4, "R+");
    }
}
    
void ZVAnalyzer::getFitInfo(TF1* p){
        cout<<"Chi^2: "<<p->GetChisquare()<<"\tNumber of DoF: "<<p->GetNDF()<<"\t(Pobability: "<<p->GetProb();
        cout<<")\n\n";
        for(int i = 0; i<70; i++) cout<<"-";
        cout<<"\n";
    }

