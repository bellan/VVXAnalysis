#ifndef FilterController_h
#define FilterController_h

/** \class FilterController
*
* No description available.
*
* $Date: 2013/01/30 21:54:21 $
* $Revision: 1.7 $
* \author N. Amapane - CERN
*/

#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/Common/interface/TriggerResults.h>
#include <FWCore/Common/interface/TriggerNames.h>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>

class FilterController {
 public:
  
  enum Trigger{DiMu, DiEle, MuEle, TriEle};
  typedef std::map<Trigger,bool> TriggerResults;

  /// Constructor
  FilterController(const edm::ParameterSet& pset);

  /// Destructor
  virtual ~FilterController(){}

  /// Pass skim
  bool passSkim(const edm::Event & event) { short bw=0; return passSkim(event,bw); }
  
  /// Pass skim (set bit in trigworld)
  bool passSkim(const edm::Event & event, short& trigworld, bool makeAnd = true);

  /// Get Simple trigger status
  TriggerResults getTriggerResults(const edm::Event & event, short& trigword);

  /// Pass trigger requests
  bool passTrigger(Channel channel, const edm::Event & event) { short bw=0; return passTrigger(channel,event,bw); }

  /// Pass trigger requests (and set bits in trigworld)
  bool passTrigger(Channel channel, const edm::Event & event, short& trigword);


  /// Pass MC filters specified in the card "MCFilterPath"
  bool passMCFilter(const edm::Event & event);
  
  /// Pass the specified filter
  bool passFilter(const edm::Event & event, const std::string& filterPath, bool fromHLT=false);

  bool isMC() const {return isMC_;};

  /// Running condition to be emulated (2011 or 2012)
  int setup() {return theSetup;};

  /// Type of MC sample (2011, 2012); can be different from setup()!
  int sampleType() {return theSampleType;}

  std::string PD;

 private:
  bool isMC_;
  int theSetup;
  int theSampleType;
  std::vector<std::string> skimPaths;
  std::string MCFilter;
  edm::EventID cachedEvtId;
  edm::Handle<edm::TriggerResults> triggerResults;
  edm::Handle<edm::TriggerResults> triggerResultsHLT;
  const edm::TriggerNames* triggerNames;
  const edm::TriggerNames* triggerNamesHLT;

  void eventInit(const edm::Event & event);

};
#endif
