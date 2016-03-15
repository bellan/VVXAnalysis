#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include <TFile.h>
#include <iostream>

LeptonScaleFactors::LeptonScaleFactors(const std::string& muonEffFilename, const std::string& electronEffFilename, const std::string& electronEffCracksfilename,
				       const std::string& muonFRFilename, const std::string& electronFRFilename){

  TFile *fEffMu = new TFile(muonEffFilename.c_str());
  TFile *fEffEl = new TFile(electronEffFilename.c_str());
  TFile *fEffElCracks = new TFile(electronEffCracksfilename.c_str());
  
  hEffMu_      = dynamic_cast<TH2F*>(fEffMu->Get("FINAL"));
  hEffEl_      = dynamic_cast<TH2F*>(fEffEl->Get("hScaleFactors_IdIsoSip"));
  hEffElCracks_ = dynamic_cast<TH2F*>(fEffElCracks->Get("hScaleFactors_IdIsoSip_Cracks"));

  TFile *fFRMu = new TFile(muonFRFilename.c_str());
  TFile *fFREl = new TFile(electronFRFilename.c_str());

  hFRMu_.first  = dynamic_cast<TH1F*>(fFRMu->Get("h1D_FRmu_EB"));                                                 
  hFRMu_.second = dynamic_cast<TH1F*>(fFRMu->Get("h1D_FRmu_EE"));
  hFREl_.first  = dynamic_cast<TH1F*>(fFREl->Get("h1D_FRel_EB"));
  hFREl_.second = dynamic_cast<TH1F*>(fFREl->Get("h1D_FRel_EE"));

  //TGraphAsym for binomial error
  grFRMu_.first  = dynamic_cast<TGraphAsymmErrors*>(fFRMu->Get("grFakeRate_NoWZ_h1D_FRmu_EB"));                                                 
  grFRMu_.second = dynamic_cast<TGraphAsymmErrors*>(fFRMu->Get("grFakeRate_NoWZ_h1D_FRmu_EE"));
  grFREl_.first  = dynamic_cast<TGraphAsymmErrors*>(fFREl->Get("grFakeRate_NoWZ_h1D_FRel_EB"));
  grFREl_.second = dynamic_cast<TGraphAsymmErrors*>(fFREl->Get("grFakeRate_NoWZ_h1D_FRel_EE"));

}

std::pair<double, double> LeptonScaleFactors::efficiencyScaleFactor(const double& pt, const double& eta, int id, bool isInCracks) const {

  const TH2F *hDataMCSF = 0;

  int xbin = -2;
  int ybin = -2;

  if(abs(id) == 13){
    hDataMCSF = hEffMu_;
    xbin = hDataMCSF->GetXaxis()->FindBin(eta);
    ybin = hDataMCSF->GetYaxis()->FindBin(pt);
    if(pt >= hDataMCSF->GetYaxis()->GetXmax()) ybin = hDataMCSF->GetYaxis()->GetLast();
    else if (pt < hDataMCSF->GetYaxis()->GetXmin()) ybin = hDataMCSF->GetYaxis()->GetFirst();   // ...should never happen
  }

  else if (abs(id) == 11) {
    hDataMCSF = hEffEl_;
    if(isInCracks) hDataMCSF = hEffElCracks_;
    xbin = hDataMCSF->GetXaxis()->FindBin(pt);
    ybin = hDataMCSF->GetYaxis()->FindBin(eta);
    if(pt >= hDataMCSF->GetXaxis()->GetXmax()) xbin = hDataMCSF->GetXaxis()->GetLast();
    else if (pt < hDataMCSF->GetXaxis()->GetXmin()) xbin = hDataMCSF->GetXaxis()->GetFirst();   // ...should never happen
  }
  else{
    std::cout << colour::Warning("Efficiency scale factor asked for an unknown particle") << " ID = " << id << std::endl;
    abort();
  }

  double sFactor    = hDataMCSF->GetBinContent(xbin,ybin); 
  double sFactorErr = hDataMCSF->GetBinError(xbin,ybin);

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


std::pair<double,std::pair<double,double>> LeptonScaleFactors::fakeRateScaleFactor(const double& lepPt, const double& lepEta, int lepId) const {
  
  double fakeRate    = 1.;
  double fakeRateUncUp = 0.;
  double fakeRateUncDown = 0.;
  double  ptvalue = 0;  
  Int_t bin = -1;
  double pt  = lepPt < 200 ? lepPt : 199;

  
  if(abs(lepId) == 13){
    if(fabs(lepEta) <= 1.2){
      bin = hFRMu_.first->GetXaxis()->FindBin(pt);
      grFRMu_.first->GetPoint(bin-1,ptvalue,fakeRate);
      fakeRateUncUp  = grFRMu_.first->GetErrorYhigh(bin-1);
      fakeRateUncDown  = grFRMu_.first->GetErrorYlow(bin-1);
      fakeRate    = hFRMu_.first->GetBinContent(bin);
    }
    else{

      bin = hFRMu_.second->GetXaxis()->FindBin(pt);
      grFRMu_.second->GetPoint(bin-1,ptvalue,fakeRate);
      fakeRateUncUp  = grFRMu_.second->GetErrorYhigh(bin-1);
      fakeRateUncDown  = grFRMu_.second->GetErrorYlow(bin-1);
      fakeRate    = hFRMu_.second->GetBinContent(bin);
    }
  }
  else if(abs(lepId) == 11){
    if(fabs(lepEta) <= 1.45){
      bin = hFREl_.first->GetXaxis()->FindBin(pt);
      grFREl_.first->GetPoint(bin-1,ptvalue,fakeRate);
      fakeRateUncUp  = grFREl_.first->GetErrorYhigh(bin-1);
      fakeRateUncDown  = grFREl_.first->GetErrorYlow(bin-1);
      fakeRate    = hFREl_.first->GetBinContent(bin);
    }
    else{
      bin = hFREl_.second->GetXaxis()->FindBin(pt);
      grFREl_.second->GetPoint(bin-1,ptvalue,fakeRate);
      fakeRateUncUp  = grFREl_.second->GetErrorYhigh(bin-1);
      fakeRateUncDown  = grFREl_.second->GetErrorYlow(bin-1);
      fakeRate    = hFREl_.second->GetBinContent(bin);
    }
  }
  else {
    abort();
  }
  
  if(fakeRate < 0.001 || fakeRate > 10.)
    std::cout << colour::Warning("Fake rate scale factor out of range") << " Lepton ID = " << lepId << ", pt =  " << pt << ", eta = " << lepEta << ", scale factor = " << fakeRate << std::endl;

  return std::make_pair(fakeRate/(1-fakeRate), std::make_pair(fakeRateUncUp/pow(1-fakeRate,2),(fakeRateUncDown/pow(1-fakeRate,2))));
}

    std::pair<double,std::pair<double,double>> LeptonScaleFactors::fakeRateScaleFactor(const phys::Lepton& lep) const{
  return fakeRateScaleFactor(lep.pt(), lep.eta(), lep.id()); 
}
