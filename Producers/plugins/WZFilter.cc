#include "VVXAnalysis/Commons/interface/AriEle.h"
#include "VVXAnalysis/Commons/interface/Constants.h"




vector<DiBosonLepton> WZCand;

BosonLepton Wtemp;
DiBosonLepton WZtemp;
ZLCompositeCandidates Zltemp;

// filter on Z's mass
foreach(const ZLCompositeCandidate Zl, *ZLCand){
  if(Zl.first.mass() > 60. && Zl.first.mass() < 120.){
    Zltemp.push_back(Zl);
  }
}

if(Zltemp.size() == 0){
  return;
}

else if(Zltemp.size() > 0){
   
  // best Z is the one with closer mass to Z
  sort(Zltemp.begin(), Zltemp.end(), pairMassComparator(0, ZMASS));
  
  // best W is the one with the lepton with higher pt
  bool choosingW = kTRUE;
  
  for(int i = 0; choosingW; i++){
    Wtemp = BosonLepton(Zltemp[i].second, Lepton(met->p4()), copysign(24, Zltemp[i].second.charge()));
    
    // filter on W's trmass
    if(Wtemp.p4().Mt() > 30 && Wtemp.p4().Mt() < 500){
      WZCand.push_back(DiBosonLepton(Wtemp, Zltemp[i].first));
    }
    
    if(i == (int)Zltemp.size() || Zltemp[i].first.mass() != Zltemp[i + 1].first.mass()){
      choosingW = kFALSE;
    }
  }
  
  if(WZCand.size() == 0){
    return;
  }
  
  // best WZ couple is the first of WZCand
  sort(WZCand.begin(), WZCand.end(), WZPtComparator());
}
