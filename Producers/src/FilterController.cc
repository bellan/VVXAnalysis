#include "VVXAnalysis/Producers/interface/FilterController.h"
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;

FilterController::FilterController(const edm::ParameterSet& pset) :
  PD(pset.getParameter<std::string>("PD")),
  isMC_(pset.getUntrackedParameter<bool>("isMC")),
  theSetup(pset.getParameter<int>("setup")),
  theSampleType(pset.getParameter<int>("sampleType")),
  skimPaths(pset.getParameter<std::vector<std::string> >("skimPaths")),
  MCFilter(pset.getParameter<std::string>("MCFilterPath")){
  
  // Check for inconsistent configurations
  if ((theSampleType != 2011 && theSampleType != 2012) ||
      ((theSampleType != theSetup) && (!isMC_ || theSampleType!=2011))) {
    cout << "ERROR: FilterController: inconsistent setup" << theSampleType << " " << theSetup << " " <<isMC_ << endl;
    abort();
  }
  
  
  if ((isMC_&&PD!="") || (!isMC_ && (PD!="DoubleEle" && PD!="DoubleMu" && PD!="MuEG"))) {
    cout << "ERROR: FilterController: isMC: " << isMC_ << " PD: " << PD << endl;
    abort();
  }

  if (!isMC_&&MCFilter!="") {
    cout << "ERROR: FilterController: MCFilter= " << MCFilter << " when isMC=0"
	 << endl;
    abort();
  }
  
}

void
FilterController::eventInit(const edm::Event & event) {
  // Initialize trigger results table
  if (event.id()==cachedEvtId) return;
  if (event.getByLabel(edm::InputTag("TriggerResults"), triggerResults)) {
    triggerNames = &(event.triggerNames(*triggerResults));
  } else {
    cout << "ERROR: failed to get TriggerResults" << endl;
  }

// if (event.getByLabel(InputTag("TriggerResults","","HLT"), triggerResultsHLT)) {
// triggerNamesHLT = &(event.triggerNames(*triggerResultsHLT));
// } else {
// cout << "ERROR: failed to get TriggerResults:HLT" << endl;
// }

  cachedEvtId = event.id();
}

bool
FilterController::passMCFilter(const edm::Event & event){
  if (MCFilter=="") return true;
  return passFilter(event, MCFilter);
}

bool
FilterController::passSkim(const edm::Event & event, short& trigworld, bool makeAnd){

  bool evtPassSkim = makeAnd ? true : false;
  if (skimPaths.size()==0) evtPassSkim=true;
  else 
    for (vector<string>::const_iterator name = skimPaths.begin(); name!= skimPaths.end(); ++name) 
      evtPassSkim = makeAnd ? evtPassSkim && passFilter(event, *name) : evtPassSkim || passFilter(event, *name);
  
  if (evtPassSkim) set_bit_16(trigworld,8);
  return evtPassSkim;
}




short FilterController::getTriggerWord(const edm::Event & event){

  short trigword = 0;

  bool passDiMu = passFilter(event, "triggerDiMu");
  bool passDiEle = passFilter(event, "triggerDiEle");
  bool passMuEle = passFilter(event, "triggerMuEle");
  if ((!isMC_) && theSetup == 2011 ) { // follow changes in trigger menu in data 2011 (see wiki)
    int irun=event.id().run();
    if (irun>=175973) {
      passMuEle = passFilter(event, "triggerMuEle3");
    } else if (irun>=167914) {
      passMuEle = passFilter(event, "triggerMuEle2");
    }
  }
  bool passTriEle = false;
  if (theSetup == 2012) {
    passTriEle = passFilter(event, "triggerTriEle");
  }

  bool passAtLeastOneTrigger = passDiMu || passDiEle || passMuEle || passTriEle;
  
  if (passAtLeastOneTrigger)   set_bit_16(trigword,0);
  if (passDiMu)                set_bit_16(trigword,1);
  if (passDiEle)               set_bit_16(trigword,2);
  if (passMuEle)               set_bit_16(trigword,3);
  if (passTriEle)              set_bit_16(trigword,4);
 
  // To be matched with channel == EEEE final state
  if ((PD=="" || PD=="DoubleEle") && (passDiEle || passTriEle))            set_bit_16(trigword,5); 
  
  // To be matched with channel == EEMM final state
  if ((PD=="" && (passDiEle || passDiMu || passMuEle)) ||
      (PD=="DoubleEle" && passDiEle) ||
      (PD=="DoubleMu" && passDiMu && !passDiEle) ||
      (PD=="MuEG" && passMuEle && !passDiMu && !passDiEle ))               set_bit_16(trigword,6);
  
  // To be matched with channel == MMMM final state
  if ((PD=="" || PD=="DoubleMu") && passDiMu)                              set_bit_16(trigword,7); 
  

  // To be matched with channel == ZLL or ZL final states
  if ((PD=="" && (passDiEle || passDiMu || passMuEle || passTriEle)) ||
      (PD=="DoubleEle" && (passDiEle || passTriEle)) ||
      (PD=="DoubleMu" && passDiMu && !passDiEle && !passTriEle) ||
      (PD=="MuEG" && passMuEle && !passDiMu && !passDiEle && !passTriEle)) set_bit_16(trigword,8);
    
  // Note: for MMMM final state the requirement is passDiMu, so it needs to be matched with bit 1 (see few line above)

  return trigword;
}


bool
FilterController::passTrigger(Channel channel, const short& trigword) const{

  if(channel == NONE) return test_bit_16(trigword,0);
  if(channel == EEEE) return test_bit_16(trigword,5);
  if(channel == EEMM) return test_bit_16(trigword,6);
  if(channel == MMMM) return test_bit_16(trigword,7);
  if(channel == ZLL || channel == ZL) return test_bit_16(trigword,8);
  else{
    edm::LogWarning("FilterController") << "Unknown channel, do not know what to do.";
    return false;
  }
  
}


bool
FilterController::passFilter(const edm::Event & event, const string& filterPath, bool fromHLT) {
  eventInit(event);

  edm::Handle<edm::TriggerResults> myTriggerResults;
  const edm::TriggerNames* myTriggerNames = 0;

  if (fromHLT) {
// myTriggerNames = triggerNamesHLT;
// myTriggerResults = triggerResultsHLT;
  } else {
    myTriggerNames = triggerNames;
    myTriggerResults = triggerResults;
  }

  unsigned i = myTriggerNames->triggerIndex(filterPath);
  if (i== myTriggerNames->size()){
    cout << "ERROR: FilterController::isTriggerBit: path does not exist! " << filterPath << endl;
    abort();
  }
  return myTriggerResults->accept(i);

}

