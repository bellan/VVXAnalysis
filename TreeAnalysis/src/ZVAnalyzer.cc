#include "VVXAnalysis/TreeAnalysis/interface/ZVAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include <TF1.h>
#include <TH1.h>
#include <vector>
#include <boost/foreach.hpp>
#include <string>
#include <time.h>

#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp>

using namespace boost::assign;

using std::cout;
using std::endl;


using namespace phys;
using namespace physmath;

Int_t ZVAnalyzer::cut() {
  return 1;
}

void ZVAnalyzer::begin(){
    cout<<"\n";
    for(char i=0; i<25; i++) cout<<"-";
    cout<<" \tBegin of ZV\t ";
    for(char i=0; i<25; i++) cout<<"-";
    cout<<"\n";
    startTime = clock();
}



void ZVAnalyzer::analyze(){

    std::vector<phys::Particle> *genElectrons = new std::vector<phys::Particle>(); //genElectrons in this events
    std::vector<phys::Particle> *genMuons = new std::vector<phys::Particle>(); //genMuons in this events
    std::vector<phys::Particle> *leptons = new std::vector<phys::Particle>(); //genLeptons in this events
    
  /*cout << "------------------------------------------------------------------"<<endl;
  cout << "Run: " << run << " event: " << event << endl;*/
  
  foreach(const phys::Particle &gen, *genParticles){
      tempStatisticParticles(gen);
      if(!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt) && abs(gen.eta())<2.4)) {
          continue;
      }
      switch (abs(gen.id())){
          case 11:
              fillParticlePlots("genElectrons", gen);
              genElectrons->push_back(gen);
              leptons->push_back(gen);
              break;
          case 13:
              fillParticlePlots("genMuons", gen);
              genMuons->push_back(gen);
              leptons->push_back(gen);
              break;
          case 15:
              fillParticlePlots("taus", gen);
              leptons->push_back(gen);
              break;
          default:
              continue;
      }
  }
      
      tempStatisticEvents();
     
      std::sort(genElectrons->begin(), genElectrons->end(), PtComparator());    //Descending order
      std::sort(genMuons->begin(), genMuons->end(), PtComparator());
      std::sort(leptons->begin(), leptons->end(), PtComparator());
    
  
    totalElectrons += genElectrons->size();
    
    
    totalMuons += genMuons->size();
    
 
    analyzeEfficiency(genElectrons, electrons, std::string("Electrons"), matchedElectrons);
    analyzeEfficiency(genMuons, muons, std::string("Muons"), matchedMuons);
    

    delete genElectrons;
    delete genMuons;
    delete leptons;
    
    return;
    
}

template<class T, class P>
void ZVAnalyzer::analyzeResolutionpt(const T& generated, const P& reconstructed, std::string name) {
    double ideltapt = 1./generated.pt() - 1./reconstructed.pt();
    theHistograms.fill(name+"pt_Resolution", name+"pt_Resolution", 200, -0.2, 0.2, ideltapt/(1./generated.pt()),1);
};

template <class T, class P>
void ZVAnalyzer::analyzeResolutionEnergy(const T& generated, const P& reconstructed, std::string name) {
    double deltaE = generated.e() - reconstructed.e();
    theHistograms.fill(name+"Energy_Resolution", name+"Energy_Resolution", 200, -0.2, 0.2, deltaE/(generated.e()),1);
};




template <class T, class P, typename C>
void ZVAnalyzer::analyzeEfficiency(std::vector<T>* genGroup, std::vector<P>* recGroup, std::string name, C& counter){
    foreach(const phys::Particle & gen, *genGroup){
       // fillParticlePlots("gen"+name, gen);
       
        phys::Particle* cand = findMatchingParticle(gen, recGroup);    //candidate is a reconstructed particle
        if(cand != nullptr){
            counter++;
            if(checkMatch(gen, *cand, 0.1)){
                //Efficiency vs Eta
                theHistograms.fill(name+"Matched_vs_eta",name+"Matched_vs_eta", 100, -3, 3, gen.eta(), 1);
                //Efficiency vs pt
                theHistograms.fill(name+"Matched_vs_pt",name+"Matched_vs_pt",100,0.,500., gen.pt(), 1);
            }
        analyzeResolutionpt(gen, *cand, name);
        analyzeResolutionEnergy(gen, *cand, name);
        }
    }
}


phys::Particle* ZVAnalyzer::findMatchingParticle(const phys::Particle& gen, std::vector<phys::Lepton>* candidates){
    if(candidates->size() == 0) return nullptr;
    minPos=0;
    minDeltaR= physmath::deltaR(gen,candidates->at(0));
    double temp=999.;
        for(int i=1;i<candidates->size();i++){
           // if (rec.charge() == candidates->at(i).charge()){
                temp=physmath::deltaR(gen,candidates->at(i));
                if (temp<minDeltaR){
                    minPos=i;
                    minDeltaR=temp;
                }
            //else return nullptr;
           // }
        }
    return & (candidates->at(minPos));
    
}

template <class P, class T>
bool ZVAnalyzer::checkMatch(const /*phys::Particle&*/P& generated, /*const phys::Particle&*/ T& reconstructed,  const float& tolerance){
    //if(reconstructed.charge() != generated.charge()) return false;    //NO
    return physmath::deltaR(generated, reconstructed) < tolerance;
}

void ZVAnalyzer::end(TFile &){
    
    normalizeHistograms(std::string("Electrons"));
    normalizeHistograms(std::string("Muons"));
    doSomeFits(std::string("Electrons"));
    doSomeFits(std::string("Muons"));
    
    
    cout << "The number of electrons is: " << totalElectrons << "and the number of matched electrons is: " << matchedElectrons << "\n";
    cout <<  "Efficiency: " << (matchedElectrons/totalElectrons)*100 << " %\n";
    
    cout << "The number of muons is: " << totalMuons << "and the number of matched muons is: " << matchedMuons << "\n";
    cout <<  "Efficiency: " << (matchedMuons/totalMuons)*100 << " %\n";
    for(char i=0; i<25; i++) cout<<"-";
    cout<<" \tEnd of ZV\t ";
    for(char i=0; i<25; i++) cout<<"-";
    cout<<"\n\n";
    
    cout<<"\nElapsed Time: "<< (float)(clock()-startTime)/CLOCKS_PER_SEC<<" s\n";
    for(char i=0; i<25; i++) cout<<"-";
    cout<<" \tEnd of ZV\t";
    for(char i=0; i<25; i++) cout<<"-";
    cout<<"\n\n";
}

void ZVAnalyzer::normalizeHistograms(std::string name){
    theHistograms.clone(name+"Matched_vs_eta", name+"Efficiency_vs_eta");
    theHistograms.get(name+"Efficiency_vs_eta")->Divide(theHistograms.get("gen"+name+"_eta"));
    theHistograms.get(name+"Efficiency_vs_eta")->SetTitle((name+"Efficiency_vs_eta").std::string::c_str());
    //There's no overload of SetTitle(const char*) with SetTitle(std::string)
    theHistograms.get(name+"Efficiency_vs_eta")->GetXaxis()->SetTitle("#eta");
    
    
    theHistograms.clone(name+"Matched_vs_pt", name+"Efficiency_vs_pt");
    theHistograms.get(name+"Efficiency_vs_pt")->Divide(theHistograms.get("gen"+name+"_pt"));
    theHistograms.get(name+"Efficiency_vs_pt")->SetTitle((name+"Efficiency_vs_pt").std::string::c_str());
    theHistograms.get(name+"Efficiency_vs_pt")->GetXaxis()->SetTitle("pt [GeV/c]");
    
   /* theHistograms.get(name+"Efficiency_vs_tolerance")->Scale(1./(float)totalElectrons);
    theHistograms.get(name+"Efficiency_vs_tolerance")->GetYaxis()->SetRangeUser(0., 1.);
    theHistograms.get(name+"Efficiency_vs_tolerance")->GetXaxis()->SetTitle("deltaR");*/
}


    void ZVAnalyzer::fillBasicPlots(){
        theHistograms.fill<TH1I>("nvtx"     , "Number of vertices" , 100, 0, 100, nvtx             , theWeight);
        theHistograms.fill<TH1I>("nmuons"    ,"Number of muons"    ,  10, 0, 10 , muons->size()    , theWeight);
        theHistograms.fill<TH1I>("nelectrons","Number of electrons",  10, 0, 10 , electrons->size(), theWeight);
        
        theHistograms.get("met")->GetXaxis()->SetTitle("[GeV/c]");
        /*foreach(const phys::Lepton& mu , *muons)     fillLeptonPlots  ("mu",  mu  );
         foreach(const phys::Electron& e, *electrons) fillElectronPlots("e" ,  e   );
         foreach(const phys::Jet& jet   , *jets)      fillJetPlots     ("jet", jet );*/
    }
    
    
    
    void ZVAnalyzer::fillParticlePlots(const std::string &type, const phys::Particle & lepton){
        theHistograms.fill(type+"_pt" ,    "p_{T} spectrum", 100,   0   , 500   ,lepton.pt()    , 1);
        theHistograms.fill(type+"_eta",    "#eta spectrum" , 100,  -3 ,   3 ,lepton.eta()   , 1);
        theHistograms.fill(type+"_phi",    "#phi spectrum" ,  50,  -3.15,   3.15,lepton.phi()   , 1);
        //theHistograms.fill(type+"_charge", "charge"        ,  50,  -25  ,  25   ,lepton.charge(), theWeight);
        
        theHistograms.get(type+"_pt")->GetXaxis()->SetTitle("[GeV/c]");
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
            
            theHistograms.fill("ptPromptParticle", "p_t (prompt)", 25, 0, 50, par.pt());
            theHistograms.fill("massPromptParticle", "m (prompt)", 1000, 0, 100, par.p4().M());
        }
        theHistograms.get("ptPromptParticle")->GetXaxis()->SetTitle("[GeV/c]");
        theHistograms.get("massPromptParticle")->GetXaxis()->SetTitle("[GeV/c^2]");
        theHistograms.fill("ptAllAnalyzedParticle", "p_t", 25, 0, 50, par.pt());
        theHistograms.get("ptAllAnalyzedParticle")->GetXaxis()->SetTitle("[GeV]");
        theHistograms.fill("ParticlesIDs", "ParticlesIDs", 40, -20, 20, par.id());
    }
    
    void ZVAnalyzer::tempStatisticEvents(){
        theHistograms.fill("ParticlesPerEvent", "ParticlesPerEvent", 40, 0, 40, counter);
        theHistograms.fill("genElectronsPerEvent", "genElectronsPerEvent", 25, 0, 25, eCounter);
        theHistograms.fill("genMuonsPerEvent", "genMuonsPerEvent", 25, 0, 25, mCounter);
        
        theHistograms.fill("PromptParticlesPerEvent", "PromptParticlesPerEvent", 25, 0, 25, promptCounter);
        theHistograms.fill("PromptgenElectronsPerEvent", "PromptgenElectronsPerEvent", 25, 0, 25, peCounter);
        theHistograms.fill("PromptMuonsPerEvent", "PromptMuonsPerEvent", 25, 0, 25, pmCounter);
        
        theHistograms.fill("ElectronsPerEvent", "ElectronsPerEvent", 25, 0, 25, electrons->size());
    }
    
void ZVAnalyzer::doSomeFits(std::string name){
    
    if (name == "Electrons") {
        TF1* func1 = new TF1("func1","[0]*exp(-x*x/(2*[1])) + [2]",-0.05,0.05);
        func1->SetParameters(2000,0.2,100);
        func1->SetParLimits(0,0,10000000);
        func1->SetParLimits(2,0,2000);
        func1->SetLineColor(4);
        theHistograms.get(name+"pt_Resolution")->Fit(func1,"R+");
        
      
        
        TF1* func2 = new TF1("func2","[0]*exp(-x*x/(2*[1])) + [2]",-0.05,0.05);
        func2->SetParameters(2000,0.2,100);
        func2->SetParLimits(0,0,10000000);
        func2->SetParLimits(2,0,2000);
        func2->SetLineColor(3);
        theHistograms.get(name+"Energy_Resolution")->Fit(func2, "R+");
    }
    if (name == "Muons"){
        TF1* func3 = new TF1("func3","[0]*exp(-x*x/(2*[1])) + [2]",-0.05,0.05);
        func3->SetParameters(2500,0.2,100);
        func3->SetParLimits(0,0,10000000);
        func3->SetParLimits(2,0,2500);
        func3->SetLineColor(4);
        theHistograms.get(name+"pt_Resolution")->Fit(func3,"R+");
        
        TF1* func4 = new TF1("func4","[0]*exp(-x*x/(2*[1])) + [2]",-0.05,0.05);
        func4->SetParameters(2500,0.2,100);
        func4->SetParLimits(0,0,10000000);
        func4->SetParLimits(2,0,2500);
        func4->SetLineColor(3);
        theHistograms.get(name+"Energy_Resolution")->Fit(func4, "R+");
    }
    
        
    }
    
    void ZVAnalyzer::getFitInfo(TF1* p){
        cout<<"Chi^2: "<<p->GetChisquare()<<"\tNumber of DoF: "<<p->GetNDF()<<"\t(Pobability: "<<p->GetProb();
        cout<<")\n\n";
        for(int i = 0; i<70; i++) cout<<"-";
        cout<<"\n";
    }
