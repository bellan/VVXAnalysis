#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include <TFile.h>
#include <iostream>

LeptonScaleFactors::LeptonScaleFactors(const std::string& muonEffFilename, const std::string& electronEffFilename, const std::string& electronEffCracksfilename,
				       const std::string& electronEffRecoFilename,
				       const std::string& muonFRFilename, const std::string& electronFRFilename){

  TFile *fEffMu       = new TFile(muonEffFilename.c_str());
  TFile *fEffEl       = new TFile(electronEffFilename.c_str());
  TFile *fEffElCracks = new TFile(electronEffCracksfilename.c_str());  
  TFile *fEffElReco   = new TFile(electronEffRecoFilename.c_str());
  
  
  hEffMu_       = dynamic_cast<TH2F*>(fEffMu->Get("FINAL"));
  hEffEl_       = dynamic_cast<TH2F*>(fEffEl->Get("EGamma_SF2D"));
  hEffElCracks_ = dynamic_cast<TH2F*>(fEffElCracks->Get("EGamma_SF2D"));
  hEffElReco_   = dynamic_cast<TH2F*>(fEffElReco->Get("EGamma_SF2D"));

  TFile *fFRMu = new TFile(muonFRFilename.c_str());
  TFile *fFREl = new TFile(electronFRFilename.c_str());

  hFRMu_  = dynamic_cast<TH2D*>(fFRMu->Get("fakeRate_m"));                                                 
  hFREl_  = dynamic_cast<TH2D*>(fFREl->Get("fakeRate_e"));

}


std::pair<double, double> LeptonScaleFactors::efficiencyScaleFactor(const double& pt, const double& eta, int id, bool isInCracks) const {

  std::string  year = "2017"; //To be added in argument constructor

  const TH2F *hDataMCSF = 0;

  int xbin = -2;
  int ybin = -2;

  double sFactor    = 0; 
  double sFactorUnc = 0;

  double recoSfactor    = 0; 
  double recoSfactorUnc = 0; 


  if(abs(id) == 13){
    hDataMCSF = hEffMu_;
    xbin = hDataMCSF->GetXaxis()->FindBin(eta);
    ybin = hDataMCSF->GetYaxis()->FindBin(pt);
    if(pt >= hDataMCSF->GetYaxis()->GetXmax()) ybin = hDataMCSF->GetYaxis()->GetLast();
    else if (pt < hDataMCSF->GetYaxis()->GetXmin()) ybin = hDataMCSF->GetYaxis()->GetFirst();   // ...should never happen
    sFactor    = hDataMCSF->GetBinContent(xbin,ybin); 
    sFactorUnc = hDataMCSF->GetBinError(xbin,ybin);
  }

  else if (abs(id) == 11) {
    if(!isInCracks){
      hDataMCSF = hEffEl_;
    }
    else {
      hDataMCSF = hEffElCracks_;
    }
    if(year =="2017"){
      xbin = hDataMCSF->GetXaxis()->FindBin(eta);
      ybin = hDataMCSF->GetYaxis()->FindBin(pt);
      if(pt >= hDataMCSF->GetYaxis()->GetXmax())      ybin = hDataMCSF->GetYaxis()->GetLast();
      else if (pt < hDataMCSF->GetYaxis()->GetXmin()) ybin = hDataMCSF->GetYaxis()->GetFirst();   // ...should never happen
    }
      else if(year=="2016"){
      xbin = hDataMCSF->GetXaxis()->FindBin(abs(eta));
      ybin = hDataMCSF->GetYaxis()->FindBin(pt);
      if(pt >= hDataMCSF->GetYaxis()->GetXmax())      ybin = hDataMCSF->GetYaxis()->GetLast();
      else if (pt < hDataMCSF->GetYaxis()->GetXmin()) ybin = hDataMCSF->GetYaxis()->GetFirst();   // ...should never happen
    }

    sFactor        =  hDataMCSF->GetBinContent(xbin,ybin); 
    recoSfactor    =  hEffElReco_->GetBinContent(hEffElReco_->FindBin(eta,50)); // The histogram depend only on eta so 50 is just a value inside the range. 
    recoSfactorUnc =  hEffElReco_->GetBinError(hEffElReco_->FindBin(eta,50));

    if(pt < 20. || pt > 80.) recoSfactorUnc += 0.01;

    sFactorUnc   = TMath::Sqrt(TMath::Power(hDataMCSF->GetBinError(xbin,ybin)/sFactor,2)+TMath::Power(recoSfactorUnc/recoSfactor,2)); 
    sFactor      *= recoSfactor;

  }
  else{
    std::cout << colour::Warning("Efficiency scale factor asked for an unknown particle") << " ID = " << id << std::endl;
    abort();
  }

  if(sFactor < 0.001 || sFactor > 10.){
    std::cout << colour::Warning("Efficiency scale factor out of range") << " Lepton ID = " << id << ", pt =  " << pt << ", eta = " << eta << ", scale factor = " << sFactor << std::endl;
  }

  return std::make_pair(sFactor, sFactorUnc) ;
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
  
