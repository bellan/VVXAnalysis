#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include <TFile.h>
#include <iostream>

LeptonScaleFactors::LeptonScaleFactors(int year, std::string data_tag) : lepSFHelper_(data_tag), year_(year){
  const char* filename = Form("$CMSSW_BASE/src/VVXAnalysis/Commons/data/leptonFakeRates_%d.root", year);

  TFile fFR(filename);
  std::cout << "[LeptonScaleFactors]: Getting FakeRate from \"" << fFR.GetName() << '\"' << std::endl;

  TH2D* hFRMu = dynamic_cast<TH2D*>(fFR.Get("fakeRate_m"));
  TH2D* hFREl = dynamic_cast<TH2D*>(fFR.Get("fakeRate_e"));
  hFRMu->SetDirectory(0);
  hFREl->SetDirectory(0);
  hFRMu_.reset(std::move(hFRMu));
  hFREl_.reset(std::move(hFREl));

  fFR.Close();
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

  TH2D *hFR = nullptr;
  switch(abs(lepId)){
  case 13:
    hFR = hFRMu_.get(); break;
  case 11:
    hFR = hFREl_.get(); break;
  default:
    abort();
  }
  
  double highest_pt = hFR->GetYaxis()->GetBinUpEdge(hFR->GetNbinsY());
  double pt = lepPt < highest_pt ? lepPt : highest_pt - 0.5;
  int bin = hFR->FindBin(fabs(lepEta), pt);
  fakeRate    = hFR->GetBinContent(bin);
  fakeRateUnc = hFR->GetBinError(bin);

  if(fakeRate < 0.001 || fakeRate > 10.)
    std::cout << colour::Warning("Fake rate scale factor out of range") << " Lepton ID = " << lepId << ", pt =  " << lepPt << ", eta = " << lepEta << ", scale factor = " << fakeRate << std::endl;

  return std::make_pair(fakeRate/(1-fakeRate), fakeRateUnc/pow(1-fakeRate,2));
}

  std::pair<double,double> LeptonScaleFactors::fakeRateScaleFactor(const phys::Lepton& lep) const{
    return fakeRateScaleFactor(lep.pt(), lep.eta(), lep.id()); 
  }
  
