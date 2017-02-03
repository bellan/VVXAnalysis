#include "VVXAnalysis/Producers/interface/FilterController.h"
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>

using namespace std;

FilterController::FilterController(const edm::ParameterSet& pset,  edm::ConsumesCollector && consumesCollector) :
  PD(pset.getParameter<std::string>("PD")),
  isMC_(pset.getUntrackedParameter<bool>("isMC")),
  theSetup(pset.getParameter<int>("setup")),
  theSampleType(pset.getParameter<int>("sampleType")),
  skimPaths(pset.getParameter<std::vector<std::string> >("skimPaths")),
  MCFilter(pset.getParameter<std::string>("MCFilterPath")),
  triggerToken_(consumesCollector.consumes<edm::TriggerResults>(edm::InputTag("TriggerResults"))){
  
  // Check for inconsistent configurations
  if ( ( theSampleType!=2011 && theSampleType!=2012 && theSampleType!=2015 && theSampleType!=2016 && theSampleType!=2017) ||
       ( theSetup!=2011 && theSetup!=2012 && theSetup!=2015 && theSetup!=2016 && theSetup!=2017) ||
       ( theSampleType!=theSetup ) // No sample rescaling supported as of now.
       // We may add exception for MC only when needed.
       ) {
    cout << "ERROR: FilterController: inconsistent setup " << theSampleType << " " << theSetup << " " <<isMC_ << endl;
    abort();
  }
  
  
  if ((isMC_&&PD!="") || (!isMC_ && (PD!="DoubleEle" && PD!="DoubleMu" && PD!="MuEG" && PD!= "SingleElectron" && PD!= "SingleMuon" ))) {
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
  if (event.getByToken(triggerToken_, triggerResults)) {
    triggerNames = &(event.triggerNames(*triggerResults));
  } else {
    cout << "ERROR: failed to get TriggerResults" << endl;
  }

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
  
  if (evtPassSkim) set_bit_16(trigworld,10); //not use
  return evtPassSkim;
}


short FilterController::getTriggerWord(const edm::Event & event){
  short trigword = 0;

  bool passDiMu  = passFilter(event, "triggerDiMu" );
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

  bool passTriEle    = false;
  bool passTriMu     = false;
  bool passSingleEle = false;
  bool passSingleMu  = false;

  if (theSetup >= 2012) {  passTriEle = passFilter(event,"triggerTriEle"); }
  if (theSetup >= 2015) {
    passTriMu     = passFilter(event,"triggerTriMu");
    passSingleEle = passFilter(event,"triggerSingleEle");
  }

  if (theSetup >= 2016) { passSingleMu  = passFilter(event,"triggerSingleMu"); }

  bool passAtLeastOneTrigger = passDiMu || passDiEle || passMuEle || passTriEle || passSingleEle || passTriMu || passSingleMu ;
  
  if (passAtLeastOneTrigger)   set_bit_16(trigword,0);
  if (passDiMu)                set_bit_16(trigword,1);
  if (passDiEle)               set_bit_16(trigword,2);
  if (passMuEle)               set_bit_16(trigword,3);
  if (passTriEle)              set_bit_16(trigword,4);
  if (passTriMu)               set_bit_16(trigword,5);
  if (passSingleEle)           set_bit_16(trigword,6);
  if (passSingleMu)            set_bit_16(trigword,7);

 
  // This is the trigger selection logic to select ZZ and ZLL events
  if( ( PD == ""                                 &&  passAtLeastOneTrigger)    || 
      ((PD == "DoubleEle" || PD == "DoubleEG"  ) && (passDiEle || passTriEle)) ||
      ((PD == "DoubleMu"  || PD == "DoubleMuon") && (passDiMu  || passTriMu)  && !passDiEle && !passTriEle) ||
      ((PD == "MuEG"      || PD == "MuonEG"    ) &&  passMuEle                && !passDiMu  && !passTriMu && !passDiEle && !passTriEle) ||
      ( PD == "SingleElectron"                   &&  passSingleEle            && !passMuEle && !passDiMu  && !passTriMu && !passDiEle && !passTriEle) ||
      ( PD == "SingleMuon"                       &&  passSingleMu             && !passMuEle && !passDiMu  && !passTriMu && !passDiEle && !passTriEle && !passSingleEle)) 
    set_bit_16(trigword,8); 

  // This is the trigger logic to select ZL events
  
  if( ( PD == ""                                 && (passDiEle || passDiMu )) ||
      //					     || (passSingleEle && !passMuEle && !passDiMu  && !passTriMu && !passDiEle && !passTriEle)  //Add single lepton only with a proper matching of triggered lepton and Z lepton 
      //					     || (passSingleMu  && !passMuEle && !passDiMu  && !passTriMu && !passDiEle && !passTriEle && !passSingleEle))) ||
      ((PD == "DoubleEle" || PD == "DoubleEG")   &&  passDiEle)              ||
      ((PD == "DoubleMu"  || PD == "DoubleMuon") &&  passDiMu && !passDiEle) ) 
      // ( PD == "SingleElectron"                   &&  passSingleEle            && !passMuEle && !passDiMu  && !passTriMu && !passDiEle && !passTriEle) ||
      // ( PD == "SingleMuon"                       &&  passSingleMu             && !passMuEle && !passDiMu  && !passTriMu && !passDiEle && !passTriEle && !passSingleEle)) 
      set_bit_16(trigword,9);
  return trigword;
}

bool
FilterController::passTrigger(Channel channel, const short& trigword) const{
  
  if(channel == NONE)                 return test_bit_16(trigword,0);
  if(channel == ZZ || channel == ZLL) return test_bit_16(trigword,8);
  if(channel == ZL)                   return test_bit_16(trigword,9);
  else{
    edm::LogWarning("FilterController") << "Unknown channel, do not know what to do.";
    return false;
  }
  
}


bool
FilterController::passFilter(const edm::Event & event, const string& filterPath) {
  eventInit(event);


  unsigned i = triggerNames->triggerIndex(filterPath);
  if (i== triggerNames->size()){
    cout << "ERROR: FilterController::isTriggerBit: path does not exist! " << filterPath << endl;
    abort();
  }
  return triggerResults->accept(i);
}

