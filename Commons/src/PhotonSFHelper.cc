#include <VVXAnalysis/Commons/interface/PhotonSFHelper.h>

#include <iostream>

#include "TFile.h"

PhotonSFHelper::PhotonSFHelper(int year, bool preVFP)
{
  const char dirpath[] = "$CMSSW_BASE/src/VVXAnalysis/Commons/data/PhotonEffScaleFactors/";
  const char *filename;
  
  switch(year){
    case 2016:
      filename = preVFP ? "PhotonSF_UL2016preVFP.root" : "PhotonSF_UL2016postVFP.root";
      break;
    case 2017:
      filename = "PhotonSF_UL2017.root";
      break;
    case 2018:
      filename = "PhotonSF_UL2018.root";
      break;
    default:
      edm::LogError("PhotonSFHelper::") << "Photon SFs for " << year << " are not supported!";
      throw std::logic_error(Form("Photon SFs for %d are not supported", year));  // abort();  // Let's not be so catastrophic
   }

  TFile root_file(Form("%s/%s", dirpath, filename), "READ");
  h_Pho.reset(std::move( (TH2F*) root_file.Get("EGamma_SF2D")->Clone() ));
  std::cout << "[PhotonSFHelper] SF map opened from root file \"" << root_file.GetName() << "\"." << std::endl;
  root_file.Close();
}


float PhotonSFHelper::getSF(float pt, float eta/*, float SCeta*/) const
{
  return h_Pho->GetBinContent(h_Pho->GetXaxis()->FindBin(eta), h_Pho->GetYaxis()->FindBin(pt));
}


float PhotonSFHelper::getSFError(float pt, float eta/*, float SCeta*/) const
{
  float SelSF = 1.0;
  float SelSF_Unc = 0.0;
  float SFError = 0.0;
  
  SelSF = h_Pho->GetBinContent(h_Pho->GetXaxis()->FindBin(eta), h_Pho->GetYaxis()->FindBin(pt));
  SelSF_Unc=h_Pho->GetBinError(h_Pho->GetXaxis()->FindBin(eta), h_Pho->GetYaxis()->FindBin(pt));
  
  SFError = SelSF != 0. ? SelSF_Unc/SelSF : 1.;
  /**
     SFError = sqrt( RecoSF_Unc*RecoSF_Unc/(RecoSF*RecoSF) + SelSF_Unc*SelSF_Unc/(SelSF*SelSF) ); // assume full correlation between different electrons (and uncorrelated reco and sel uncertainties)
  **/

   return SFError;
}
