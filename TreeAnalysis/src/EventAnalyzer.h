#ifndef EventAnalyzer_h
#define EventAnalyzer_h

#include <vector>
#include <string>
#include <utility>
#include <map>

#include "Electron.h"
#include "Lepton.h"
#include "Jet.h"

class TFile;
class TTree;
class TBranch;
class TH1;

class EventAnalyzer {
public:
  enum METType {Std,NoMu,NoEl};

  EventAnalyzer(std::string filename, double lumi = 1., double externalXSection = -1.);
  virtual ~EventAnalyzer();

  struct PtComparator{
    template<typename LEP>
    bool operator()( const LEP & a , 
		     const LEP & b) const{ 
      return a.p4().Pt() > b.p4().Pt(); 
    }
  };
  
  // To steer the loop over all events. User is not supposed to change this.
  virtual void     loop(const std::string outputfile);

 protected:
  // Functions to be overloaded in the concrete instance of the EventAnalyzer class.
  virtual void  begin() {}
  virtual void  book()  {}
  virtual Int_t cut();
  virtual void  analyze() = 0;
  virtual void  end(TFile &) {};  
  
 private:
  // Structural functions to access the tree branches. User is not supposed to change these.
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  
  // Some basic plots. User may want to change these, thou they should be used only for very basic plots.
  void book(const std::string&);
  void bookExtra(const std::string&);
  virtual void     basicPlots();
  void fillLeptonPlots(const std::string&, const phys::Lepton&);
  void fillElectronExtraPlots(const std::string &type, const phys::Electron &electron);
  

 private:
  TTree *theTree;
  int fCurrent; 

 protected:
  // Histograms container
  std::map<std::string,TH1*> thePlots;
  
  double theWeight;
  double theSampleWeight;
  int    theCutCounter;

  // Access to the branches
  Int_t    nvtx  ; TBranch *b_nvtx   ;
  Double_t rho   ; TBranch *b_rho    ;
  Double_t weight; TBranch *b_weight ;
  Double_t puweight; TBranch *b_puweight ;
  Double_t xsec  ; TBranch *b_xsec   ;
  Int_t totalEvents; TBranch *b_totalEvents;

  //MET
  phys::Particle *met   ; TBranch *b_met    ;

  //Muons
  std::vector<phys::Lepton> *muons      ; TBranch *b_muons    ;    

  //Electrons
  std::vector<phys::Electron> *electrons; TBranch *b_electrons;

  // Jets  
  std::vector<phys::Jet> *jets; TBranch *b_jets;
  
};
#endif


