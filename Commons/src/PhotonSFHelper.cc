#include <VVXAnalysis/Commons/interface/PhotonSFHelper.h>

using namespace std;

PhotonSFHelper::PhotonSFHelper(bool preVFP)
{
   // 2016 preVFP & postVFP
   TString    fipPhoton_2016 = Form("$CMSSW_BASE/src/VVXAnalysis/Commons/data/PhotonEffScaleFactors/PhotonSF_UL2016postVFP.root");
   if(preVFP) fipPhoton_2016 = Form("$CMSSW_BASE/src/VVXAnalysis/Commons/data/PhotonEffScaleFactors/PhotonSF_UL2016preVFP.root");
   root_file = TFile::Open(fipPhoton_2016,"READ");
   h_Pho_2016 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();

   // 2017  
    TString fipPhoton_2017 = Form("$CMSSW_BASE/src/VVXAnalysis/Commons/data/PhotonEffScaleFactors/PhotonSF_UL2017.root");
    root_file = TFile::Open(fipPhoton_2017,"READ");
    h_Pho_2017 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();
 
   // 2018
    TString fipPhoton_2018 = Form("$CMSSW_BASE/src/VVXAnalysis/Commons/data/PhotonEffScaleFactors/PhotonSF_UL2018.root");
    root_file = TFile::Open(fipPhoton_2018,"READ");
    h_Pho_2018 = (TH2F*) root_file->Get("EGamma_SF2D")->Clone();

   cout << "[PhotonSFHelper] SF maps opened from root files." << endl;
}

PhotonSFHelper::~PhotonSFHelper()
{
}

float PhotonSFHelper::getSF(int year, float pt, float eta/*, float SCeta*/) const
{
  //float RecoSF = 1.0;
   float SelSF = 1.0;
   float SF = 1.0;

   //cout << "year = " << year << " pt = " << pt << " eta = " << eta << " SCeta = " << SCeta  << endl;

   // Photon ID SFs
     if     (year == 2016)
      {
       SelSF = h_Pho_2016->GetBinContent(h_Pho_2016->GetXaxis()->FindBin(eta),h_Pho_2016->GetYaxis()->FindBin(pt));
      }
     else if (year == 2017)
      {
       SelSF = h_Pho_2017->GetBinContent(h_Pho_2017->GetXaxis()->FindBin(eta),h_Pho_2017->GetYaxis()->FindBin(pt));
      }
     else if (year == 2018)
      {
       SelSF = h_Pho_2018->GetBinContent(h_Pho_2016->GetXaxis()->FindBin(eta),h_Pho_2018->GetYaxis()->FindBin(pt));
      }
      else {
         edm::LogError("PhotonSFHelper::") << "Photon SFs for " << year << " is not supported!";
         abort();
      }

      SF = SelSF;

    return SF;
}


float PhotonSFHelper::getSFError(int year, float pt, float eta/*, float SCeta*/) const
{
  //float RecoSF = 1.0;
   float SelSF = 1.0;
   //   float RecoSF_Unc = 0.0;
   float SelSF_Unc = 0.0;
   float SFError = 0.0;

   // Photon ID SFs
     if     (year == 2016)
      {
       SelSF = h_Pho_2016->GetBinContent(h_Pho_2016->GetXaxis()->FindBin(eta),h_Pho_2016->GetYaxis()->FindBin(pt));
       SelSF_Unc=h_Pho_2016->GetBinError(h_Pho_2016->GetXaxis()->FindBin(eta),h_Pho_2016->GetYaxis()->FindBin(pt));
      }
     else if (year == 2017)
      {
       SelSF = h_Pho_2017->GetBinContent(h_Pho_2017->GetXaxis()->FindBin(eta),h_Pho_2017->GetYaxis()->FindBin(pt));
       SelSF_Unc=h_Pho_2017->GetBinError(h_Pho_2017->GetXaxis()->FindBin(eta),h_Pho_2017->GetYaxis()->FindBin(pt));
      }
     else if (year == 2018)
      {
       SelSF = h_Pho_2018->GetBinContent(h_Pho_2016->GetXaxis()->FindBin(eta),h_Pho_2018->GetYaxis()->FindBin(pt));
       SelSF_Unc=h_Pho_2018->GetBinError(h_Pho_2018->GetXaxis()->FindBin(eta),h_Pho_2018->GetYaxis()->FindBin(pt));
      }
      else {
         edm::LogError("PhotonSFHelper::") << "Photon SFs for " << year << " is not supported!";
         abort();
      }

      SFError = SelSF_Unc/SelSF;
      /**
SFError = sqrt( RecoSF_Unc*RecoSF_Unc/(RecoSF*RecoSF) + SelSF_Unc*SelSF_Unc/(SelSF*SelSF) ); // assume full correlation between different electrons (and uncorrelated reco and sel uncertainties)
      **/

   return SFError;
}
