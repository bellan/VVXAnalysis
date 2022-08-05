#include <VVXAnalysis/Commons/interface/PhotonSFHelper.h>

#include <iostream>

#include "TFile.h"
#include "TSystem.h"

PhotonSFHelper::PhotonSFHelper(int year, bool preVFP)
{
  if(year < 2016 || year > 2018){
    auto erroromsg = Form("Photon SFs for %d are not supported", year);
    edm::LogError("VVXAnalysis|Commons|PhotonSFHelper") << "[PhotonSFHelper] " << erroromsg;
    throw std::logic_error(erroromsg);
  }
  
  const char* filename = _getSFFilename(year, preVFP);
  if(filename == nullptr)
    edm::LogError("VVXAnalysis|Commons|PhotonSFHelper") << "[PhotonSFHelper] No SF file found for year " << year << (year==2016 ? Form(", preVFP: %d", preVFP) : "" );
  
  TFile root_file(filename, "READ");
  
  h_Pho.reset(std::move( (TH2F*) root_file.Get("EGamma_SF2D")->Clone() ));
  
  edm::LogInfo("VVXAnalysis|Commons|PhotonSFHelper") << "[PhotonSFHelper] SF map opened from root file \"" << root_file.GetName() << "\".";
  std::cout << "[PhotonSFHelper] SF map opened from root file \"" << root_file.GetName() << "\".\n";
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


const char* PhotonSFHelper::_getSFFilename(int year, bool preVFP) const{
  const char *filename;
  
  const char* thefile;
  switch(year){
    case 2016:
      thefile = preVFP ? "egammaEffi.txt_EGM2D_Pho_Loose_UL16.root" : "egammaEffi.txt_EGM2D_Pho_Loose_UL16_postVFP.root";
      break;
    case 2017:
      thefile = "egammaEffi.txt_EGM2D_Pho_Loose_UL17.root";
      break;
    case 2018:
      thefile = "egammaEffi.txt_EGM2D_Pho_Loose_UL18.root";
      break;
  }
  
  // First try the official egamma folder
  const char repository[] = "/eos/cms/store/group/phys_egamma/SF-Repository";
  
  if(!gSystem->AccessPathName(repository)){
    const char *folder = Form("UL%d", year%100);
    
    if(year == 2016)
      folder = Form("%s/%s", folder, preVFP ? "preVFP" : "postVFP");    
    
    const char *filename = Form("%s/%s/Photons/Loose/%s", repository, folder, thefile);
    
    if(!gSystem->AccessPathName(filename))
      return filename;
  }
  
  // If the official repository is unreachable, fall back to local copies
  filename = Form("$CMSSW_BASE/src/VVXAnalysis/Commons/data/%s", thefile);
  
  if(!gSystem->AccessPathName(filename))
    return filename;
  else
    return nullptr;
}
