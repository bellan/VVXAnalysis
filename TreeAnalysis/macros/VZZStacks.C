/**
 *  Macro to plot together VZZ signal and background, from results of WZZAnalyzer 
 *			
 *  $Date: 2019/10/17 11:56:47 
 *  $Revision: 0.2 $
 *
 *  \author C. Tarricone cristiano.tarrico@edu.unito.it
 */

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TMath.h>

#include "THStack.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define TEST_MODE

using namespace std;

void VZZStacks(){

  TString path =  "~/VVXAnalysis/TreeAnalysis/results/WZZAnalyzer_SR/";      



  THStack* hs = new THStack("hs","");


  TFile* resultWZ = TFile::Open(path + "WZ.root");
  TH1F* h2 = (TH1F*)resultWZ->Get("recoVMass_bckg0");

  h2->SetFillColor(kBlue);
  hs->Add(h2);


  TFile* resultggZZ2e2mu = TFile::Open(path + "ggZZ2e2mu.root");
  TH1F* h3 = (TH1F*)resultggZZ2e2mu->Get("recoVMass_bckg0");

  h3->SetFillColor(kGreen);
  hs->Add(h3);


  TFile* resultggZZ4mu = TFile::Open(path + "ggZZ4mu.root");
  TH1F* h4 = (TH1F*)resultggZZ4mu->Get("recoVMass_bckg0");

  h4->SetFillColor(kYellow);
  hs->Add(h4);


  TFile* resultggZZ4e = TFile::Open(path + "ggZZ4e.root");
  TH1F* h5 = (TH1F*)resultggZZ4e->Get("recoVMass_bckg0");

  h5->SetFillColor(kOrange);
  hs->Add(h5);


  TFile* resultWZZ = TFile::Open(path + "WZZ.root");
  TH1F* h6 = (TH1F*)resultWZZ->Get("recoVMass_bckg0");

  h6->SetFillColor(kPink);
  hs->Add(h6);


  TFile* resultZZZ = TFile::Open(path + "ZZZ.root");
  TH1F* h7 = (TH1F*)resultZZZ->Get("recoVMass_bckg0");

  h7->SetFillColor(kMagenta);
  hs->Add(h7);


  TFile* resultZZ = TFile::Open(path +"ZZTo4lamcatnlo.root" );
  TH1F* h1 = (TH1F*)resultZZ->Get("recoVMass_bckg0");

  h1->SetFillColor(kViolet);
  hs->Add(h1);

  
  TH1F* h8 = (TH1F*)resultZZZ->Get("recoVMass_sign0");

  h8->SetFillColor(kRed);
  hs->Add(h8);


  TH1F* h9 = (TH1F*)resultWZZ->Get("recoVMass_sign0");

  h9->SetFillColor(7);
  hs->Add(h9);

   
  TCanvas *cs = new TCanvas("cs","cs");
  TText T; T.SetTextFont(42); T.SetTextAlign(21);
  cs->Divide(2,2);
  cs->cd(4); hs->Draw("lego1"); T.DrawTextNDC(.5,.95,"Option \"lego1\"");
  return cs;

}
