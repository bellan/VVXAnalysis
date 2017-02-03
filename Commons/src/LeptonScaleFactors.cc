#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include <TFile.h>
#include <iostream>

//LeptonScaleFactors::LeptonScaleFactors(const std::string& muonEffFilename, const std::string& electronEffFilename, const std::string& electronEffCracksfilename,
LeptonScaleFactors::LeptonScaleFactors(const std::string& muonEffFilename, const std::string& electronEffFilename, const std::string& electronGSFEffFilename,
				       const std::string& muonFRFilename, const std::string& electronFRFilename){

  TFile *fEffMu   = new TFile(muonEffFilename.c_str());
  TFile *fEffEl   = new TFile(electronEffFilename.c_str());
  TFile *fGSEffEl = new TFile(electronGSFEffFilename.c_str());
  //  TFile *fEffElCracks = new TFile(electronEffCracksfilename.c_str());
  
  hEffMu_       = dynamic_cast<TH2F*>(fEffMu->Get("FINAL"));
  hEffEl_       = dynamic_cast<TH2F*>(fEffEl->Get("ele_scale_factors"));
  hEffElCracks_ = dynamic_cast<TH2F*>(fEffEl->Get("ele_scale_factors_gap"));

  hEffGSEl_         = dynamic_cast<TH2F*>(fGSEffEl->Get("EGamma_SF2D"));
  //2015
  // hEffEl_       = dynamic_cast<TH2F*>(fEffEl->Get("hScaleFactors_IdIsoSip"));
  // hEffElCracks_ = dynamic_cast<TH2F*>(fEffElCracks->Get("hScaleFactors_IdIsoSip_Cracks"));

  TFile *fFRMu = new TFile(muonFRFilename.c_str());
  TFile *fFREl = new TFile(electronFRFilename.c_str());

  hFRMu_.first  = dynamic_cast<TH1F*>(fFRMu->Get("h1D_FRmu_EB"));                                                 
  hFRMu_.second = dynamic_cast<TH1F*>(fFRMu->Get("h1D_FRmu_EE"));
  hFREl_.first  = dynamic_cast<TH1F*>(fFREl->Get("h1D_FRel_EB"));
  hFREl_.second = dynamic_cast<TH1F*>(fFREl->Get("h1D_FRel_EE"));

  //TGraphAsym for binomial error
  // grFRMu_.first  = dynamic_cast<TGraphAsymmErrors*>(fFRMu->Get("grFakeRate_NoWZ_h1D_FRmu_EB"));                                                 
  // grFRMu_.second = dynamic_cast<TGraphAsymmErrors*>(fFRMu->Get("grFakeRate_NoWZ_h1D_FRmu_EE"));
  // grFREl_.first  = dynamic_cast<TGraphAsymmErrors*>(fFREl->Get("grFakeRate_NoWZ_h1D_FRel_EB"));
  // grFREl_.second = dynamic_cast<TGraphAsymmErrors*>(fFREl->Get("grFakeRate_NoWZ_h1D_FRel_EE"));

}

std::pair<double, double> LeptonScaleFactors::efficiencyScaleFactor(const double& pt, const double& eta, int id, bool isInCracks) const {

  std::string  year = "2016"; //To be added in argument constructor

  const TH2F *hDataMCSF = 0;
  //  const TH2F *hDataMCSF_Unc = 0; del

  int xbin = -2;
  int ybin = -2;

  double sFactor    = 0; 
  double sFactorErr = 0;

  if(abs(id) == 13){
    hDataMCSF = hEffMu_;
    xbin = hDataMCSF->GetXaxis()->FindBin(eta);
    ybin = hDataMCSF->GetYaxis()->FindBin(pt);
    if(pt >= hDataMCSF->GetYaxis()->GetXmax()) ybin = hDataMCSF->GetYaxis()->GetLast();
    else if (pt < hDataMCSF->GetYaxis()->GetXmin()) ybin = hDataMCSF->GetYaxis()->GetFirst();   // ...should never happen
    sFactor    = hDataMCSF->GetBinContent(xbin,ybin); 
    sFactorErr = hDataMCSF->GetBinError(xbin,ybin);

  }

  else if (abs(id) == 11) {
    if(!isInCracks){
      hDataMCSF = hEffEl_;
      //      hDataMCSF_Unc = hEffEl_Unc_;
    }
    else {
      hDataMCSF = hEffElCracks_;
      //    hDataMCSF_Unc = hEffElCracks_Unc_;
    }
    if(year =="2015"){
      xbin = hDataMCSF->GetXaxis()->FindBin(pt);
      ybin = hDataMCSF->GetYaxis()->FindBin(eta);
      if(pt >= hDataMCSF->GetXaxis()->GetXmax()) xbin = hDataMCSF->GetXaxis()->GetLast();
      else if (pt < hDataMCSF->GetXaxis()->GetXmin()) xbin = hDataMCSF->GetXaxis()->GetFirst();   // ...should never happen
      
    }
    else if(year=="2016"){
      xbin = hDataMCSF->GetXaxis()->FindBin(abs(eta));
      ybin = hDataMCSF->GetYaxis()->FindBin(pt);
      if(pt >= hDataMCSF->GetYaxis()->GetXmax()) ybin = hDataMCSF->GetYaxis()->GetLast();
      else if (pt < hDataMCSF->GetYaxis()->GetXmin()) ybin = hDataMCSF->GetYaxis()->GetFirst();   // ...should never happen
    }

    sFactor    = hDataMCSF->GetBinContent(xbin,ybin); 
    //    sFactorErr = hDataMCSF_Unc->GetBinError(xbin,ybin); //in 2016 ele_scale_factors_uncertainties content and errors are equals 
    sFactorErr = TMath::Sqrt(TMath::Power(hDataMCSF->GetBinError(xbin,ybin),2)+TMath::Power(hEffGSEl_->GetBinError(eta,20),2)); 

    sFactor*=  hEffGSEl_->GetBinContent(hEffGSEl_->FindBin(eta,20)); // The histogram depend only on eta so 20 is just a value inside the range. 
  }
  else{
    std::cout << colour::Warning("Efficiency scale factor asked for an unknown particle") << " ID = " << id << std::endl;
    abort();
  }


  if(sFactor < 0.001 || sFactor > 10.){
    std::cout << colour::Warning("Efficiency scale factor out of range") << " Lepton ID = " << id << ", pt =  " << pt << ", eta = " << eta << ", scale factor = " << sFactor << std::endl;
  }

  return std::make_pair(sFactor, sFactorErr) ;
}

std::pair<double, double> LeptonScaleFactors::efficiencyScaleFactor(const phys::Lepton& lep) const{
  
  return lep.passFullSel() ? efficiencyScaleFactor(lep.pt(), lep.eta(), lep.id(), lep.isInCracks()) : std::make_pair(1.,0.); 
}


double LeptonScaleFactors::weight(const phys::Boson<phys::Lepton> &Z) const{
  return efficiencyScaleFactor(Z.daughter(0)).first * efficiencyScaleFactor(Z.daughter(1)).first;
}

double LeptonScaleFactors::weight(const phys::DiBoson<phys::Lepton,phys::Lepton> &ZZ) const{
  return weight(ZZ.first()) * weight(ZZ.second());
}


//std::pair<double,std::pair<double,double>> LeptonScaleFactors::fakeRateScaleFactor(const double& lepPt, const double& lepEta, int lepId) const {
std::pair<double,double> LeptonScaleFactors::fakeRateScaleFactor(const double& lepPt, const double& lepEta, int lepId) const {
  
  double fakeRate    = 1.;
  double fakeRateUnc = 0.;
  //  double  ptvalue = 0;  
  Int_t bin = -1;
  double pt  = lepPt < 200 ? lepPt : 199;

  
  if(abs(lepId) == 13){
    if(fabs(lepEta) <= 1.2){
      bin = hFRMu_.first->GetXaxis()->FindBin(pt);
      fakeRate    = hFRMu_.first->GetBinContent(bin);
      fakeRateUnc = hFRMu_.first->GetBinError(bin);
    }
    else{
      bin = hFRMu_.second->GetXaxis()->FindBin(pt);
      fakeRate    = hFRMu_.second->GetBinContent(bin);
      fakeRateUnc = hFRMu_.second->GetBinError(bin);
    }
  }
  else if(abs(lepId) == 11){
    if(fabs(lepEta) <= 1.45){
      bin = hFREl_.first->GetXaxis()->FindBin(pt);
      fakeRate    = hFREl_.first->GetBinContent(bin);
      fakeRateUnc = hFREl_.first->GetBinError(bin);
    }
    else{
      bin = hFREl_.second->GetXaxis()->FindBin(pt);

      fakeRate    = hFREl_.second->GetBinContent(bin);
      fakeRateUnc = hFREl_.second->GetBinError(bin);
    }
  }
  else {
    abort();
  }
  
  if(fakeRate < 0.001 || fakeRate > 10.)
    std::cout << colour::Warning("Fake rate scale factor out of range") << " Lepton ID = " << lepId << ", pt =  " << pt << ", eta = " << lepEta << ", scale factor = " << fakeRate << std::endl;

  return std::make_pair(fakeRate/(1-fakeRate), fakeRateUnc/pow(1-fakeRate,2));
}
  std::pair<double,double> LeptonScaleFactors::fakeRateScaleFactor(const phys::Lepton& lep) const{
    return fakeRateScaleFactor(lep.pt(), lep.eta(), lep.id()); 
  }
  
