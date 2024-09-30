#include <iostream>
#include <fstream>
#include <math.h>   
#include "TH1F.h"
#include "TF1.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TAxis.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include <TFile.h>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace std;
//using namespace boost::assign;
//using namespace colour;


void plotDrawer() {
  TFile* myFile = TFile::Open("results/2018/WlllnuAnalyzer_MC/WZTo3LNu.root");
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------- //
  // ------- W DECAY MODE 1: e+ e- mu nu -------- //
  
  // -- Gen leptons pt -- //
  TH1F* histoGenElPt = (TH1F *)myFile->Get("GEN_el_pt_mode1");
  TCanvas *c1 = new TCanvas("c1","GEN_el_pt_mode1",200,10,600,400);
  c1->SetFillColor(0);
  c1->cd();
  //histoGenElPt->SetLineColor(1);
  histoGenElPt->Draw();
  
  TH1F* histoGenMuPt = (TH1F *)myFile->Get("GEN_mu_pt_mode1");
  TCanvas *c2 = new TCanvas("c2","GEN_mu_pt_mode1",200,10,600,400);
  c2->SetFillColor(0);
  c2->cd();
  histoGenMuPt->Draw();
  
  // -- Invariant mass of the e+e- pair -- //
  TH1F* histoElPairInvariantMass = (TH1F *)myFile->Get("GEN_el_pair_invariant_mass_mode1");
  TCanvas *c3 = new TCanvas("c3","GEN_el_pair_invariant_mass_mode1",200,10,600,400);
  c3->SetFillColor(0);
  c3->cd();
  histoElPairInvariantMass->Draw();
  
  // -- Invariant mass of the 4 leptons e+ e- mu nu -- //
  TH1F* histoFourLepInvariantMass = (TH1F *)myFile->Get("GEN_four_leptons_invariant_mass_mode1");
  TCanvas *c4 = new TCanvas("c4","GEN_four_leptons_invariant_mass_mode1",200,10,600,400);
  c4->SetFillColor(0);
  c4->cd();
  histoFourLepInvariantMass->Draw();
  
  // -- Transverse mass of the 3 charged leptons e+ e- mu -- //
  TH1F* histoChLepTransverseMass = (TH1F *)myFile->Get("GEN_charged_leptons_transverse_mass_mode1");
  TCanvas *c5 = new TCanvas("c5","GEN_charged_leptons_transverse_mass_mode1",200,10,600,400);
  c5->SetFillColor(0);
  c5->cd();
  histoChLepTransverseMass->Draw();
  
  // -- RECONSTRUCTED LEPTONS EFFICIENCY e+ e- mu -- //
  
  // -- Gen Leptons & Rec Leptons match -- //
  TH1F* histoRecElCompatibility = (TH1F *)myFile->Get("REC_el_compatibility_mode1");
  TCanvas *c6 = new TCanvas("c6","REC_el_compatibility_mode1",200,10,600,400);
  c6->SetFillColor(0);
  c6->cd();
  histoRecElCompatibility->Draw();
  
  TH1F* histoRecMuCompatibility = (TH1F *)myFile->Get("REC_mu_compatibility_mode1");
  TCanvas *c7 = new TCanvas("c7","REC_mu_compatibility_mode1",200,10,600,400);
  c7->SetFillColor(0);
  c7->cd();
  histoRecMuCompatibility->Draw();
  
  // -- Gen Leptons & Rec Leptons pt match (deltaR < 0.01) -- //
  TH1F* histoRecElPtCompatibility = (TH1F *)myFile->Get("REC_el_pt_compatibility_mode1");
  TCanvas *c8 = new TCanvas("c8","REC_el_pt_compatibility_mode1",200,10,600,400);
  c8->SetFillColor(0);
  c8->cd();
  histoRecElPtCompatibility->Draw();
  
  TH1F* histoRecMuPtCompatibility = (TH1F *)myFile->Get("REC_mu_pt_compatibility_mode1");
  TCanvas *c9 = new TCanvas("c9","REC_mu_pt_compatibility_mode1",200,10,600,400);
  c9->SetFillColor(0);
  c9->cd();
  histoRecMuPtCompatibility->Draw();
  
  // -- Gen Muon & Rec Muon charge compatibility -- //
  TH1F* histoRecMuChargeCompatibility = (TH1F *)myFile->Get("REC_mu_charge_compatibility_mode1");
  TCanvas *c10 = new TCanvas("c10","REC_mu_charge_compatibility_mode1",200,10,600,400);
  c10->SetFillColor(0);
  c10->cd();
  histoRecMuChargeCompatibility->Draw();
  
  // -------------- EFFICIENCY --------------- //
  
  // -- Rec & Gen Electron Efficiency -- //
  TH1F* histoRecGenElEff = (TH1F *)histoRecElCompatibility->Clone("histoRecGenElEff");
  histoRecGenElEff->Divide(histoGenElPt);
  TCanvas *c11 = new TCanvas("c11","REC_el_efficiency_mode1",200,10,600,400);
  c11->SetFillColor(0);
  c11->cd();
  histoRecGenElEff->Draw();
  
  // -- Rec & Gen Muon Efficiency -- //
  TH1F* histoRecGenMuEff = (TH1F*)histoRecMuCompatibility->Clone("histoRecGenMuEff");
  histoRecGenMuEff->Divide(histoGenMuPt);
  TCanvas *c12 = new TCanvas("c12","REC_mu_efficiency_mode1",200,10,600,400);
  c12->SetFillColor(0);
  c12->cd();
  histoRecGenMuEff->Draw();
  
  // -- Rec & Gen Electron pt Efficiency -- //
  TH1F* histoRecGenElPtEff = (TH1F*)histoRecElPtCompatibility->Clone("histoRecGenElPtEff");
  histoRecGenElPtEff->Divide(histoGenElPt);
  TCanvas *c13 = new TCanvas("c13","REC_el_pt_efficiency_mode1",200,10,600,400);
  c13->SetFillColor(0);
  c13->cd();
  histoRecGenElPtEff->Draw();
  
  // -- Rec & Gen Muon pt Efficiency -- //
  TH1F* histoRecGenMuPtEff = (TH1F*)histoRecMuPtCompatibility->Clone("histoRecGenMuPtEff");
  histoRecGenMuPtEff->Divide(histoGenMuPt);
  TCanvas *c14 = new TCanvas("c14","REC_mu_pt_efficiency_mode1",200,10,600,400);
  c14->SetFillColor(0);
  c14->cd();
  histoRecGenMuPtEff->Draw();
  
  // -- Rec & Gen Muon charge Efficiency -- //
  TH1F* histoRecGenMuChargeEff = (TH1F*)histoRecMuChargeCompatibility->Clone("histoRecGenMuChargeEff");
  histoRecGenMuChargeEff->Divide(histoGenMuPt);
  TCanvas *c15 = new TCanvas("c15","REC_mu_charge_efficiency_mode1",200,10,600,400);
  c15->SetFillColor(0);
  c15->cd();
  histoRecGenMuChargeEff->Draw();
  
  
  
  
  
  /*
  // ----------------------------------------------------------------------------------------------------------------------------------------------------- //
  // ------- W DECAY MODE 2: mu+ mu- e nu -------- //
  
  // -- Gen leptons pt -- //
  TH1F* histoGenElPt_A = (TH1F *)myFile->Get("GEN_el_pt_mode2");
  TCanvas *cA1 = new TCanvas("cA1","GEN_el_pt_mode2",200,10,600,400);
  cA1->SetFillColor(0);
  cA1->cd();
  histoGenElPt_A->Draw();
  
  TH1F* histoGenMuPt_A = (TH1F *)myFile->Get("GEN_mu_pt_mode2");
  TCanvas *cA2 = new TCanvas("cA2","GEN_mu_pt_mode2",200,10,600,400);
  cA2->SetFillColor(0);
  cA2->cd();
  histoGenMuPt_A->Draw();
  
  // -- Invariant mass of the mu+mu- pair -- //
  TH1F* histoMuPairInvariantMass_A = (TH1F *)myFile->Get("GEN_mu_pair_invariant_mass_mode2");
  TCanvas *cA3 = new TCanvas("cA3","GEN_mu_pair_invariant_mass_mode2",200,10,600,400);
  cA3->SetFillColor(0);
  cA3->cd();
  histoMuPairInvariantMass_A->Draw();
  
  // -- Invariant mass of the 4 leptons mu+ mu- e nu -- //
  TH1F* histoFourLepInvariantMass_A = (TH1F *)myFile->Get("GEN_four_leptons_invariant_mass_mode2");
  TCanvas *cA4 = new TCanvas("cA4","GEN_four_leptons_invariant_mass_mode2",200,10,600,400);
  cA4->SetFillColor(0);
  cA4->cd();
  histoFourLepInvariantMass_A->Draw();
  
  // -- Transverse mass of the 3 charged leptons mu+ mu- e -- //
  TH1F* histoChLepTransverseMass_A = (TH1F *)myFile->Get("GEN_charged_leptons_transverse_mass_mode2");
  TCanvas *cA5 = new TCanvas("cA5","GEN_charged_leptons_transverse_mass_mode2",200,10,600,400);
  cA5->SetFillColor(0);
  cA5->cd();
  histoChLepTransverseMass_A->Draw(); 
  
  // -- RECONSTRUCTED LEPTONS EFFICIENCY e+ e- mu -- //
  
  // -- Gen Leptons & Rec Leptons match -- //
  TH1F* histoRecElCompatibility_A = (TH1F *)myFile->Get("REC_el_compatibility_mode2");
  TCanvas *cA6 = new TCanvas("cA6","REC_el_compatibility_mode2",200,10,600,400);
  cA6->SetFillColor(0);
  cA6->cd();
  histoRecElCompatibility_A->Draw();
  
  TH1F* histoRecMuCompatibility_A = (TH1F *)myFile->Get("REC_mu_compatibility_mode2");
  TCanvas *cA7 = new TCanvas("cA7","REC_mu_compatibility_mode2",200,10,600,400);
  cA7->SetFillColor(0);
  cA7->cd();
  histoRecMuCompatibility_A->Draw();
  
  // -- Gen Leptons & Rec Leptons pt match (deltaR < 0.01) -- //
  TH1F* histoRecElPtCompatibility_A = (TH1F *)myFile->Get("REC_el_pt_compatibility_mode2");
  TCanvas *cA8 = new TCanvas("cA8","REC_el_pt_compatibility_mode2",200,10,600,400);
  cA8->SetFillColor(0);
  cA8->cd();
  histoRecElPtCompatibility_A->Draw();
  
  TH1F* histoRecMuPtCompatibility_A = (TH1F *)myFile->Get("REC_mu_pt_compatibility_mode2");
  TCanvas *cA9 = new TCanvas("cA9","REC_mu_pt_compatibility_mode2",200,10,600,400);
  cA9->SetFillColor(0);
  cA9->cd();
  histoRecMuPtCompatibility_A->Draw();
  
  // -- Gen Electron & Rec Electron charge compatibility -- //
  TH1F* histoRecElChargeCompatibility_A = (TH1F *)myFile->Get("REC_el_charge_compatibility_mode2");
  TCanvas *cA10 = new TCanvas("cA10","REC_mel_charge_compatibility_mode2",200,10,600,400);
  cA10->SetFillColor(0);
  cA10->cd();
  histoRecElChargeCompatibility_A->Draw();
  */
  
  
  
  
  
  /*
  // -- Gen Electrons Efficiency -- //
  TH1* histoGenElPt = theHistograms->get("GEN_el_pt");
  TH1* histoGenElMatchPt = theHistograms->get("GEN_REC_el_pt");
  TH1* histoGenElEff = (TH1*)histoGenElPt->Clone("histoGenElEff");
  histoGenElEff->Divide(histoGenElPt);
  TCanvas *c1 = new TCanvas("c1","GEN_REC_el_pt",200,10,600,400);
  c1->SetFillColor(0);
  c1->cd();
  histoGenElEff->Draw();
  // -- Gen Muons Efficiency -- //
  TH1* histoGenMuPt = theHistograms->get("GEN_mu_pt");
  TH1* histoGenMuMatchPt = theHistograms->get("GEN_REC_mu_pt");
  TH1* histoGenMuEff = (TH1*)histoGenMuPt->Clone("histoGenMuEff");
  histoGenMuEff->Divide(histoGenMuPt);
  TCanvas *c2 = new TCanvas("c2","GEN_REC_mu_pt",200,10,600,400);
  c2->SetFillColor(0);
  c2->cd();
  histoGenMuEff->Draw();
  */
}
