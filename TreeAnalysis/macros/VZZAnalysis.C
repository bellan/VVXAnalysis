/**
 *	Macro to get analysis of Z and W hadronic decay 
 *	
 *	It may become necessary to change the path and/or to expand the samples list	
 *	
 *  $Date: 2018/09/18 11:20:03 $
 *  $Revision: 0.1 $
 *  \author C. Tarricone cristiano.tarrico@edu.unito.it
 */

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define TEST_MODE

using namespace std;

void VZZAnalysis(){
  
  TFile* resultWZZ = TFile::Open("~/VVXAnalysis/TreeAnalysis/results/WZZAnalyzer_MC/WZZ.root");
  TFile* resultZZZ = TFile::Open("~/VVXAnalysis/TreeAnalysis/results/WZZAnalyzer_MC/ZZZ.root");
  
  TH1F* hFakeZ = (TH1F*)resultWZZ->Get("genZmass");
  TH1F* hW = (TH1F*)resultWZZ->Get("genWmass");

  TH1F* hFakeW =  (TH1F*)resultZZZ->Get("genWmass");
  TH1F* hZ  = (TH1F*)resultZZZ->Get("genZmass");

 

  
  hW->Draw();
  hW->SetLineColor(kBlue);

  hFakeW->Draw("same");
  hFakeW->SetLineColor(kRed);

  new TCanvas();
  
  hZ->Draw();
  hZ->SetLineColor(kBlue);

  hFakeZ->Draw("same");
  hFakeZ->SetLineColor(kRed);


  TList *listW = new TList;
  listW->Add(hW);
  listW->Add(hFakeW);

  TH1F *hTotW = (TH1F*)hW->Clone("hTotW");
  hTotW->Reset();
  hTotW->Merge(listW);
  new TCanvas();
  hTotW->Draw();

  
  TList *listZ = new TList;
  listZ->Add(hZ);
  listZ->Add(hFakeZ);

  TH1F *hTotZ = (TH1F*)hZ->Clone("hTotZ");
  hTotZ->Reset();
  hTotZ->Merge(listZ);
  new TCanvas();
  hTotZ->Draw();


  








}
