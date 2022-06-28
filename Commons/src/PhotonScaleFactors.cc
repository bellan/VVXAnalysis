#include "VVXAnalysis/Commons/interface/PhotonScaleFactors.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include <TFile.h>
#include <iostream>

PhotonScaleFactors::PhotonScaleFactors(int year,
				       const std::string& Filename, bool preVFP) : phoSFHelper_(preVFP), year_(year){
  TFile *f = new TFile(Filename.c_str());

}

std::pair<double, double> PhotonScaleFactors::efficiencyScaleFactor(const phys::Photon &pho) const{
  /***
   * NEED TO ADD SUPERCLUSTER ETA FOR PHOTON
    float eta = abs(lep.id())==11 ? lep.scEta() : lep.eta(); !!! NEED TO ADD SC.ETA PHOTON
    ***/
  return pho.cutBasedIDLoose() ? efficiencyScaleFactor(pho.pt(), pho.eta() /*, pho.id()*/) : std::make_pair(1.,0.); 


}

double PhotonScaleFactors::weight(const phys::Photon &pho) const{
  return efficiencyScaleFactor(pho).first;
}


/*std::pair<double,double> PhotonScaleFactors::fakeRateScaleFactor(const double& lepPt, const double& lepEta, int lepId) const {
  double fakeRate    = 1.;
  double fakeRateUnc = 0.;

  //Distinguish channel Z->ee and Z->mumu
  if(abs(lepId) == 13){
    double pt  = lepPt < 200 ? lepPt : 199;
    fakeRate      =  hFRMu_->GetBinContent( hFRMu_->FindBin(fabs(lepEta),pt)); 
    fakeRateUnc   =  hFRMu_->GetBinError( hFRMu_->FindBin(fabs(lepEta),pt)); 
    //    std::cout<<"13  pt "<<pt<<" eta "<<fabs(lepEta)<<" "<<fakeRate/(1-fakeRate)<<std::endl;
  }
  
  else if (abs(lepId) == 11){
    double pt  = lepPt < 200 ? lepPt : 199;
    //    double pt  = lepPt < 80 ? lepPt : 79;
    fakeRate      =  hFREl_->GetBinContent( hFREl_->FindBin(fabs(lepEta),pt)); 
    fakeRateUnc   =  hFREl_->GetBinError( hFREl_->FindBin(fabs(lepEta),pt)); 
    //    std::cout<<"11  pt "<<pt<<" eta "<<fabs(lepEta)<<" "<<fakeRate/(1-fakeRate)<<std::endl;
  }

  else {
    abort();
  }
  
  if(fakeRate < 0.001 || fakeRate > 10.)
    std::cout << colour::Warning("Fake rate scale factor out of range") << " Lepton ID = " << lepId << ", pt =  " << lepPt << ", eta = " << lepEta << ", scale factor = " << fakeRate << std::endl;

  return std::make_pair(fakeRate/(1-fakeRate), fakeRateUnc/pow(1-fakeRate,2));
}
  std::pair<double,double> LeptonScaleFactors::fakeRateScaleFactor(const phys::Lepton& lep) const{
    return fakeRateScaleFactor(lep.pt(), lep.eta(), lep.id()); 
  }
*/
