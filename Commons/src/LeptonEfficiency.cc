#include "VVXAnalysis/Commons/interface/LeptonEfficiency.h"
#include "VVXAnalysis/Commons/interface/Colours.h"

#include <TFile.h>
#include <iostream>

LeptonEfficiency::LeptonEfficiency(){
  TFile *fMu = new TFile("../../VVXAnalysis/Commons/data/scale_factors_muons2012.root");
  TFile *fEl = new TFile("../../VVXAnalysis/Commons/data/scale_factors_ele2012.root");
  
  hMu_ = dynamic_cast<TH2F*>(fMu->Get("TH2D_ALL_2012"));
  hEl_ = dynamic_cast<TH2F*>(fEl->Get("h_electronScaleFactor_RecoIdIsoSip"));
}

double LeptonEfficiency::scaleFactor(const double& lepPt, const double& lepEta, int lepId) const {

  double sFactor = 1.;
  //double errCorr = 0.;
  //double errCorrSyst = 0.;
  
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
    sFactor  = hMu_->GetBinContent(hMu_->GetXaxis()->FindBin(pt), hMu_->GetYaxis()->FindBin(eta));
    //errCorr = hMu_->GetBinError  (hMu_->GetXaxis()->FindBin(pt), hMu_->GetYaxis()->FindBin(eta));
    //add the systematics on T&P corrections (for muons only, electrons have them already included)
    //if(pt >= 15.) errCorrSyst = 0.005;
    //else errCorrSyst = 0.015;
    if(pt < 5.) sFactor = 0.;
  }
  else if(id == 11){
    sFactor  = hEl_->GetBinContent(hEl_->GetXaxis()->FindBin(pt), hEl_->GetYaxis()->FindBin(eta));
    //errCorr = hEl_->GetBinError  (hEl_->GetXaxis()->FindBin(pt), hEl_->GetYaxis()->FindBin(eta));
  }
  else {
    abort();
  }
  
  if(sFactor < 0.001 || sFactor > 10.){
    std::cout << colour::Warning("Efficiency scale factor out of range") << " Lepton ID = " << id << ", pt =  " << pt << ", eta = " << eta << ", scale factor = " << sFactor << std::endl;
  }
  return sFactor;
}

double LeptonEfficiency::weight(const phys::Boson<phys::Lepton> &Z) const{
  return scaleFactor(Z.daughter(0).pt(),Z.daughter(0).eta(),Z.daughter(0).id()) * 
    scaleFactor(Z.daughter(1).pt(),Z.daughter(1).eta(),Z.daughter(1).id());
}

double LeptonEfficiency::weight(const phys::DiBoson<phys::Lepton,phys::Lepton> &ZZ) const{
  return weight(ZZ.first()) * weight(ZZ.second());
}


