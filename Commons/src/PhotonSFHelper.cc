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

  maxPt = h_Pho->GetYaxis()->GetBinUpEdge(h_Pho->GetNbinsY());
  
  edm::LogInfo("VVXAnalysis|Commons|PhotonSFHelper") << "[PhotonSFHelper] SF map opened from root file \"" << root_file.GetName() << "\".";
  std::cout << "[PhotonSFHelper] SF map opened from root file \"" << root_file.GetName() << "\".\n";
  root_file.Close();
}


int PhotonSFHelper::findPhotonBin(float pt, float eta) const
{
  // Check for pt greater than last bin
  if(pt > maxPt) pt = maxPt - 0.1;
  return h_Pho->FindFixBin(eta, pt);
}


float PhotonSFHelper::getSF(float pt, float eta/*, float SCeta*/) const
{
  return h_Pho->GetBinContent(findPhotonBin(pt, eta));
}


float PhotonSFHelper::getSFError(float pt, float eta/*, float SCeta*/) const
{
  return h_Pho->GetBinError(findPhotonBin(pt, eta));
}


const char* PhotonSFHelper::_getSFFilename(int year, bool preVFP) const{
  const char *path;          // Final path to the file
  const char POGrepo[] =     "/eos/cms/store/group/phys_egamma/SF-Repository";  // Path to the central EGamma repository
  const char *folderEGamma;  // Folder inside the EGamma repository
  const char *filename;      // Name of the SF file both in the EGamma folder and in our local "data/"

  switch(year){
    case 2016:
      filename = preVFP ?
	"egammaEffi.txt_EGM2D_Pho_Loose_UL16.root" :
	"egammaEffi.txt_EGM2D_Pho_Loose_UL16_postVFP.root";
      folderEGamma = preVFP ?
	"UL16/preVFP/Photons/Loose" :
	"UL16/postVFP/Photons/Loose";
      break;
    case 2017:
      filename = "egammaEffi.txt_EGM2D.root";
      folderEGamma = "UL17/Photons/Loose/passingLoose100XV2_lowpTChebychev_addGaus";
      break;
    case 2018:
      filename = "egammaEffi.txt_EGM2D.root";
      folderEGamma = "UL18/Photons/Loose/passingLoose100XV2";
      break;
  }
  
  // First try the official egamma folder
  if(!gSystem->AccessPathName(POGrepo)){
    path = Form("%s/%s/%s", POGrepo, folderEGamma, filename);

    if(!gSystem->AccessPathName(path))
      return path;
    else
      std::cout << Form("[PhotonSFHelper] SF file not found in EGamma repo: (%s)\n", path);
  }
  // If the official repository is unreachable, fall back to local copies
  else{
    path = Form("$CMSSW_BASE/src/VVXAnalysis/Commons/data/%s", filename);
    
    if(!gSystem->AccessPathName(path))
      return path;
    else
      std::cout << Form("[PhotonSFHelper] SF file not found local folder: (%s)\n", path);
  }
  
  return nullptr;
}
