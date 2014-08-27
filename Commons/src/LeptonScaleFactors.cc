#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
#include "VVXAnalysis/Commons/interface/Colours.h"

#include <TFile.h>
#include <iostream>

LeptonScaleFactors::LeptonScaleFactors(const std::string& muonEffFilename, const std::string& electronEffFilename,
				       const std::string& muonFRFilename, const std::string& electronFRFilename){

  TFile *fEffMu = new TFile(muonEffFilename.c_str());
  TFile *fEffEl = new TFile(electronEffFilename.c_str());
  
  hEffMu_ = dynamic_cast<TH2F*>(fEffMu->Get("TH2D_ALL_2012"));
  hEffEl_ = dynamic_cast<TH2F*>(fEffEl->Get("h_electronScaleFactor_RecoIdIsoSip"));

  TFile *fFRMu = new TFile(muonFRFilename.c_str());
  TFile *fFREl = new TFile(electronFRFilename.c_str());

  hFRMu_.first  = dynamic_cast<TH1D*>(fFRMu->Get("h1D_FRmu_EB"));
  hFRMu_.second = dynamic_cast<TH1D*>(fFRMu->Get("h1D_FRmu_EE"));
  hFREl_.first  = dynamic_cast<TH1D*>(fFREl->Get("h1D_FRel_EB"));
  hFREl_.second = dynamic_cast<TH1D*>(fFREl->Get("h1D_FRel_EE"));

}

double LeptonScaleFactors::efficiencyScaleFactor(const double& lepPt, const double& lepEta, int lepId) const {

  double sFactor = 1.;
  
  double pt  = lepPt;
  double eta = lepEta;
  int    id  = abs(lepId);
  
  //avoid to go out of the TH boundary
  if(id == 13 && pt > 99.) pt = 99.;
  else if(id == 11){
    eta = abs(lepEta);
    if(pt > 199.) pt = 199.;
  }
  
  if(id == 13){
    sFactor  = hEffMu_->GetBinContent(hEffMu_->GetXaxis()->FindBin(pt), hEffMu_->GetYaxis()->FindBin(eta));
    if(pt < 5.) sFactor = 0.;
  }
  else if(id == 11){
    sFactor  = hEffEl_->GetBinContent(hEffEl_->GetXaxis()->FindBin(pt), hEffEl_->GetYaxis()->FindBin(eta));
  }
  else {
    abort();
  }
  
  if(sFactor < 0.001 || sFactor > 10.){
    std::cout << colour::Warning("Efficiency scale factor out of range") << " Lepton ID = " << id << ", pt =  " << pt << ", eta = " << eta << ", scale factor = " << sFactor << std::endl;
  }
  return sFactor;
}

double LeptonScaleFactors::efficiencyScaleFactor(const phys::Lepton& lep) const{
  return efficiencyScaleFactor(lep.pt(), lep.eta(), lep.id()); 
}

double LeptonScaleFactors::weight(const phys::Boson<phys::Lepton> &Z) const{
  return efficiencyScaleFactor(Z.daughter(0)) * efficiencyScaleFactor(Z.daughter(1));
}

double LeptonScaleFactors::weight(const phys::DiBoson<phys::Lepton,phys::Lepton> &ZZ) const{
  return weight(ZZ.first()) * weight(ZZ.second());
}


std::pair<double,double> LeptonScaleFactors::fakeRateScaleFactor(const double& lepPt, const double& lepEta, int lepId) const {
  
  double fakeRate    = 1.;
  double fakeRateUnc = 0.;
  
  double pt  = lepPt < 80 ? pt : 79;

  
  if(abs(lepId) == 13){
    if(fabs(lepEta) <= 1.2){
      fakeRate    = hFRMu_.first->GetBinContent(hFRMu_.first->GetXaxis()->FindBin(pt));
      fakeRateUnc = hFRMu_.first->GetBinError(hFRMu_.first->GetXaxis()->FindBin(pt));
    }
    else{
      fakeRate    = hFRMu_.second->GetBinContent(hFRMu_.second->GetXaxis()->FindBin(pt));
      fakeRateUnc = hFRMu_.second->GetBinError(hFRMu_.second->GetXaxis()->FindBin(pt));
    }
  }
  else if(abs(lepId) == 11){
    if(fabs(lepEta) <= 1.45){
      fakeRate    = hFREl_.first->GetBinContent(hFREl_.first->GetXaxis()->FindBin(pt));
      fakeRateUnc = hFREl_.first->GetBinError(hFREl_.first->GetXaxis()->FindBin(pt));
    }
    else{
      fakeRate    = hFREl_.second->GetBinContent(hFREl_.second->GetXaxis()->FindBin(pt));
      fakeRateUnc = hFREl_.second->GetBinError(hFREl_.second->GetXaxis()->FindBin(pt));
    }
  }
  else {
    abort();
  }
  
  if(fakeRate < 0.001 || fakeRate > 10.){
    std::cout << colour::Warning("Efficiency scale factor out of range") << " Lepton ID = " << lepId << ", pt =  " << pt << ", eta = " << lepEta << ", scale factor = " << fakeRate << std::endl;
  }

  return std::make_pair(fakeRate/(1-fakeRate),fakeRateUnc/pow(1-fakeRate,2));
}


std::pair<double,double> LeptonScaleFactors::fakeRateScaleFactor(const phys::Lepton& lep) const{
  return fakeRateScaleFactor(lep.pt(), lep.eta(), lep.id()); 
}
