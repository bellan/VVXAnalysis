#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include <TFile.h>
#include <iostream>

LeptonScaleFactors::LeptonScaleFactors(int year,
				       const std::string& muonFRFilename, const std::string& electronFRFilename):year_(year){


  TFile *fFRMu = new TFile(muonFRFilename.c_str());
  TFile *fFREl = new TFile(electronFRFilename.c_str());

  hFRMu_  = dynamic_cast<TH2D*>(fFRMu->Get("fakeRate_m"));                                                 
  hFREl_  = dynamic_cast<TH2D*>(fFREl->Get("fakeRate_e"));

}


std::pair<double, double> LeptonScaleFactors::efficiencyScaleFactor(const phys::Lepton& lep) const{
  float eta = abs(lep.id())==11 ? lep.scEta() : lep.eta();
  return lep.passFullSel() ? efficiencyScaleFactor(lep.pt(), eta, lep.id(), lep.isInCracks()) : std::make_pair(1.,0.); 
}


double LeptonScaleFactors::weight(const phys::Boson<phys::Lepton> &Z) const{
  return efficiencyScaleFactor(Z.daughter(0)).first * efficiencyScaleFactor(Z.daughter(1)).first;
}

double LeptonScaleFactors::weight(const phys::DiBoson<phys::Lepton,phys::Lepton> &ZZ) const{
  return weight(ZZ.first()) * weight(ZZ.second());
}

std::pair<double,double> LeptonScaleFactors::fakeRateScaleFactor(const double& lepPt, const double& lepEta, int lepId) const {
  
  double fakeRate    = 1.;
  double fakeRateUnc = 0.;

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
  
