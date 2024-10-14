#include <iostream>
#include <fstream>
#include <math.h>   
#include "TH1F.h"
#include "TF1.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TGraph.h"
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
  //TCanvas *c1 = new TCanvas("c1","GEN_el_pt_mode1",200,10,600,400);
  //c1->SetFillColor(0);
  //c1->cd();
  //histoGenElPt->Draw();
  
  TH1F* histoGenMuPt = (TH1F *)myFile->Get("GEN_mu_pt_mode1");
  //TCanvas *c2 = new TCanvas("c2","GEN_mu_pt_mode1",200,10,600,400);
  //c2->SetFillColor(0);
  //c2->cd();
  //histoGenMuPt->Draw();
  
  // -- Invariant mass of the e+e- pair -- //
  TH1F* histoElPairInvariantMass = (TH1F *)myFile->Get("GEN_el_pair_invariant_mass_mode1");
  //TCanvas *c3 = new TCanvas("c3","GEN_el_pair_invariant_mass_mode1",200,10,600,400);
  //c3->SetFillColor(0);
  //c3->cd();
  //histoElPairInvariantMass->Draw();
  
  // -- Invariant mass of the 4 leptons e+ e- mu nu -- //
  TH1F* histoFourLepInvariantMass = (TH1F *)myFile->Get("GEN_four_leptons_invariant_mass_mode1");
  TCanvas *c4 = new TCanvas("c4","GEN_four_leptons_invariant_mass_mode1",200,10,600,400);
  c4->SetFillColor(0);
  c4->cd();
  histoFourLepInvariantMass->Draw();
  
  // -- Transverse mass of the 3 charged leptons e+ e- mu -- //
  TH1F* histoChLepTransverseMass = (TH1F *)myFile->Get("GEN_charged_leptons_transverse_mass_mode1");
  //TCanvas *c5 = new TCanvas("c5","GEN_charged_leptons_transverse_mass_mode1",200,10,600,400);
  //c5->SetFillColor(0);
  //c5->cd();
  //histoChLepTransverseMass->Draw();
  
  
  // -- RECONSTRUCTED LEPTONS EFFICIENCY e+ e- mu -- //
  
  // -- Gen Leptons & Rec Leptons match -- //
  TH1F* histoRecElCompatibility = (TH1F *)myFile->Get("REC_el_compatibility_mode1");
  //TCanvas *c6 = new TCanvas("c6","REC_el_compatibility_mode1",200,10,600,400);
  //c6->SetFillColor(0);
  //c6->cd();
  //histoRecElCompatibility->Draw();
  
  TH1F* histoRecMuCompatibility = (TH1F *)myFile->Get("REC_mu_compatibility_mode1");
  //TCanvas *c7 = new TCanvas("c7","REC_mu_compatibility_mode1",200,10,600,400);
  //c7->SetFillColor(0);
  //c7->cd();
  //histoRecMuCompatibility->Draw();
  
  // -- Gen Leptons & Rec Leptons pt match (deltaR < 0.01) -- //
  TH1F* histoRecElPtCompatibility = (TH1F *)myFile->Get("REC_el_pt_compatibility_mode1");
  //TCanvas *c8 = new TCanvas("c8","REC_el_pt_compatibility_mode1",200,10,600,400);
  //c8->SetFillColor(0);
  //c8->cd();
  //histoRecElPtCompatibility->Draw();
  
  TH1F* histoRecMuPtCompatibility = (TH1F *)myFile->Get("REC_mu_pt_compatibility_mode1");
  //TCanvas *c9 = new TCanvas("c9","REC_mu_pt_compatibility_mode1",200,10,600,400);
  //c9->SetFillColor(0);
  //c9->cd();
  //histoRecMuPtCompatibility->Draw();
  
  // -- Gen Muon & Rec Muon charge compatibility -- //
  TH1F* histoRecMuChargeCompatibility = (TH1F *)myFile->Get("REC_mu_charge_compatibility_mode1");
  //TCanvas *c10 = new TCanvas("c10","REC_mu_charge_compatibility_mode1",200,10,600,400);
  //c10->SetFillColor(0);
  //c10->cd();
  //histoRecMuChargeCompatibility->Draw();
  
  // -- Correlation Factor GEN invariant mass GEN transverse mass -- //
  TH2F* histoGenFourLepTransverseInvariantMass_mode1 = (TH2F *)myFile->Get("GEN_four_lep_transv_mass_four_lep_inv_mass_mode1");
  cout << " " << endl;
  cout << "Correlation factor between GEN invariant & transverse mass (mode 1) = " << histoGenFourLepTransverseInvariantMass_mode1->GetCorrelationFactor() << endl;
  
  // -- Correlation Factor GEN invariant mass REC transverse mass -- //
  TH2F* histoRecGenFourLepTransverseMass_mode1 = (TH2F *)myFile->Get("REC_four_lep_transv_mass_GEN_four_lep_transv_mass_mode1");
  cout << "Correlation factor between REC & GEN transverse mass (mode 1) = " << histoRecGenFourLepTransverseMass_mode1->GetCorrelationFactor() << endl;
  cout << " " << endl;
  
  /*
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
  */
  
  
  // -------------- Signal effiency and background efficiency --------------- //
  TH1F* signalEffNum_mode1 = (TH1F *)myFile->Get("GEN_REC_signal_mode1");
  TH1F* signalEffDen_mode1 = (TH1F *)myFile->Get("GEN_signal_mode1");
  TH1F* signalEff_mode1 = (TH1F*)signalEffNum_mode1->Clone("signalEff_mode1");
  signalEff_mode1->Divide(signalEffDen_mode1);
  TCanvas *c16 = new TCanvas("c16","Signal_efficiency_mode1",200,10,600,400);
  c16->SetFillColor(0);
  c16->cd();
  signalEff_mode1->Draw();
  
  TH1F* backgroundEffNum_mode1 = (TH1F *)myFile->Get("not_GEN_REC_signal_mode1");
  TH1F* backgroundEffDen_mode1 = (TH1F *)myFile->Get("not_GEN_signal_mode1");
  TH1F* backgroundEff_mode1 = (TH1F*)backgroundEffNum_mode1->Clone("backgroundEff_mode1");
  TCanvas *c17 = new TCanvas("c17","Background_efficiency_mode1",200,10,600,400);
  backgroundEff_mode1->Divide(backgroundEffDen_mode1);
  c17->SetFillColor(0);
  c17->cd();
  backgroundEff_mode1->Draw();
  
  
  // -- Decay type events distribution -- //
  double Wtype_mode1 = histoFourLepInvariantMass->Integral(histoFourLepInvariantMass->FindBin(75),histoFourLepInvariantMass->FindBin(85));
  double Ztype_mode1 = histoFourLepInvariantMass->Integral(histoFourLepInvariantMass->FindBin(88),histoFourLepInvariantMass->FindBin(130));
  double WZtype_mode1 = histoFourLepInvariantMass->Integral(histoFourLepInvariantMass->FindBin(170),histoFourLepInvariantMass->FindBin(280));
  
  cout << "#events_Wtype (mode 1) = " << Wtype_mode1 << ";   #events_Ztype (mode 1) = " << Ztype_mode1 << ";   #events_WZtype (mode 1) = " << WZtype_mode1 << endl;
  cout << "(#events_Ztype)/(#events_Wtype) (mode 1) = " << Ztype_mode1/Wtype_mode1 << ";   (#events_WZtype)/(#events_Wtype) (mode 1) = " << WZtype_mode1/Wtype_mode1 << endl;
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------- //
  // ------- W DECAY MODE 2: mu+ mu- e nu -------- //
  
  // -- Gen leptons pt -- //
  TH1F* histoGenElPt_A = (TH1F *)myFile->Get("GEN_el_pt_mode2");
  //TCanvas *cA1 = new TCanvas("cA1","GEN_el_pt_mode2",200,10,600,400);
  //cA1->SetFillColor(0);
  //cA1->cd();
  //histoGenElPt_A->Draw();
  
  TH1F* histoGenMuPt_A = (TH1F *)myFile->Get("GEN_mu_pt_mode2");
  //TCanvas *cA2 = new TCanvas("cA2","GEN_mu_pt_mode2",200,10,600,400);
  //cA2->SetFillColor(0);
  //cA2->cd();
  //histoGenMuPt_A->Draw();
  
  // -- Invariant mass of the mu+mu- pair -- //
  TH1F* histoMuPairInvariantMass_A = (TH1F *)myFile->Get("GEN_mu_pair_invariant_mass_mode2");
  //TCanvas *cA3 = new TCanvas("cA3","GEN_mu_pair_invariant_mass_mode2",200,10,600,400);
  //cA3->SetFillColor(0);
  //cA3->cd();
  //histoMuPairInvariantMass_A->Draw();
  
  // -- Invariant mass of the 4 leptons mu+ mu- e nu -- //
  TH1F* histoFourLepInvariantMass_A = (TH1F *)myFile->Get("GEN_four_leptons_invariant_mass_mode2");
  TCanvas *cA4 = new TCanvas("cA4","GEN_four_leptons_invariant_mass_mode2",200,10,600,400);
  cA4->SetFillColor(0);
  cA4->cd();
  histoFourLepInvariantMass_A->Draw();
  
  // -- Transverse mass of the 3 charged leptons mu+ mu- e -- //
  TH1F* histoChLepTransverseMass_A = (TH1F *)myFile->Get("GEN_charged_leptons_transverse_mass_mode2");
  //TCanvas *cA5 = new TCanvas("cA5","GEN_charged_leptons_transverse_mass_mode2",200,10,600,400);
  //cA5->SetFillColor(0);
  //cA5->cd();
  //histoChLepTransverseMass_A->Draw(); 
  
  // -- RECONSTRUCTED LEPTONS EFFICIENCY e+ e- mu -- //
  
  // -- Gen Leptons & Rec Leptons match -- //
  TH1F* histoRecElCompatibility_A = (TH1F *)myFile->Get("REC_el_compatibility_mode2");
  //TCanvas *cA6 = new TCanvas("cA6","REC_el_compatibility_mode2",200,10,600,400);
  //cA6->SetFillColor(0);
  //cA6->cd();
  //histoRecElCompatibility_A->Draw();
  
  TH1F* histoRecMuCompatibility_A = (TH1F *)myFile->Get("REC_mu_compatibility_mode2");
  //TCanvas *cA7 = new TCanvas("cA7","REC_mu_compatibility_mode2",200,10,600,400);
  //cA7->SetFillColor(0);
  //cA7->cd();
  //histoRecMuCompatibility_A->Draw();
  
  // -- Gen Leptons & Rec Leptons pt match (deltaR < 0.01) -- //
  TH1F* histoRecElPtCompatibility_A = (TH1F *)myFile->Get("REC_el_pt_compatibility_mode2");
  //TCanvas *cA8 = new TCanvas("cA8","REC_el_pt_compatibility_mode2",200,10,600,400);
  //cA8->SetFillColor(0);
  //cA8->cd();
  //histoRecElPtCompatibility_A->Draw();
  
  TH1F* histoRecMuPtCompatibility_A = (TH1F *)myFile->Get("REC_mu_pt_compatibility_mode2");
  //TCanvas *cA9 = new TCanvas("cA9","REC_mu_pt_compatibility_mode2",200,10,600,400);
  //cA9->SetFillColor(0);
  //cA9->cd();
  //histoRecMuPtCompatibility_A->Draw();
  
  // -- Gen Electron & Rec Electron charge compatibility -- //
  TH1F* histoRecElChargeCompatibility_A = (TH1F *)myFile->Get("REC_el_charge_compatibility_mode2");
  //TCanvas *cA10 = new TCanvas("cA10","REC_mel_charge_compatibility_mode2",200,10,600,400);
  //cA10->SetFillColor(0);
  //cA10->cd();
  //histoRecElChargeCompatibility_A->Draw();
  
  // -- Correlation Factor GEN invariant mass GEN transverse mass -- //
  TH2F* histoGenFourLepTransverseInvariantMass_mode2 = (TH2F *)myFile->Get("GEN_four_lep_transv_mass_four_lep_inv_mass_mode2");
  cout << " " << endl;
  cout << "Correlation factor between GEN invariant & transverse mass (mode 2) = " << histoGenFourLepTransverseInvariantMass_mode2->GetCorrelationFactor() << endl;
  
  // -- Correlation Factor GEN invariant mass REC transverse mass -- //
  TH2F* histoRecGenFourLepTransverseMass_mode2 = (TH2F *)myFile->Get("REC_four_lep_transv_mass_GEN_four_lep_transv_mass_mode2");
  cout << "Correlation factor between REC & GEN transverse mass (mode 2) = " << histoRecGenFourLepTransverseMass_mode2->GetCorrelationFactor() << endl;
  cout << " " << endl;
  
  /*
  // -------------- EFFICIENCY --------------- //
  
  // -- Rec & Gen Electron Efficiency -- //
  TH1F* histoRecGenElEff_A = (TH1F *)histoRecElCompatibility_A->Clone("histoRecGenElEff_A");
  histoRecGenElEff_A->Divide(histoGenElPt_A);
  TCanvas *cA11 = new TCanvas("cA11","REC_el_efficiency_mode2",200,10,600,400);
  cA11->SetFillColor(0);
  cA11->cd();
  histoRecGenElEff_A->Draw();
  
  // -- Rec & Gen Muon Efficiency -- //
  TH1F* histoRecGenMuEff_A = (TH1F*)histoRecMuCompatibility_A->Clone("histoRecGenMuEff_A");
  histoRecGenMuEff_A->Divide(histoGenMuPt_A);
  TCanvas *cA12 = new TCanvas("cA12","REC_mu_efficiency_mode2",200,10,600,400);
  cA12->SetFillColor(0);
  cA12->cd();
  histoRecGenMuEff_A->Draw();
  
  // -- Rec & Gen Electron pt Efficiency -- //
  TH1F* histoRecGenElPtEff_A = (TH1F*)histoRecElPtCompatibility_A->Clone("histoRecGenElPtEff_A");
  histoRecGenElPtEff_A->Divide(histoGenElPt_A);
  TCanvas *cA13 = new TCanvas("cA13","REC_el_pt_efficiency_mode2",200,10,600,400);
  cA13->SetFillColor(0);
  cA13->cd();
  histoRecGenElPtEff_A->Draw();
  
  // -- Rec & Gen Muon pt Efficiency -- //
  TH1F* histoRecGenMuPtEff_A = (TH1F*)histoRecMuPtCompatibility_A->Clone("histoRecGenMuPtEff_A");
  histoRecGenMuPtEff_A->Divide(histoGenMuPt_A);
  TCanvas *cA14 = new TCanvas("cA14","REC_mu_pt_efficiency_mode2",200,10,600,400);
  cA14->SetFillColor(0);
  cA14->cd();
  histoRecGenMuPtEff_A->Draw();
  
  // -- Rec & Gen Muon charge Efficiency -- //
  TH1F* histoRecGenElChargeEff_A = (TH1F*)histoRecElChargeCompatibility_A->Clone("histoRecGenElChargeEff_A");
  histoRecGenElChargeEff_A->Divide(histoGenElPt_A);
  TCanvas *cA15 = new TCanvas("cA15","REC_el_charge_efficiency_mode2",200,10,600,400);
  cA15->SetFillColor(0);
  cA15->cd();
  histoRecGenElChargeEff_A->Draw();
  */
  
  
  // -------------- Signal efficiency and background efficiency --------------- //
  TH1F* signalEffNum_mode2 = (TH1F *)myFile->Get("GEN_REC_signal_mode2");
  TH1F* signalEffDen_mode2 = (TH1F *)myFile->Get("GEN_signal_mode2");
  TH1F* signalEff_mode2 = (TH1F*)signalEffNum_mode2->Clone("signalEff_mode2");
  signalEff_mode2->Divide(signalEffDen_mode2);
  TCanvas *cA16 = new TCanvas("cA16","Signal_efficiency_mode2",200,10,600,400);
  cA16->SetFillColor(0);
  cA16->cd();
  signalEff_mode2->Draw();
  
  TH1F* backgroundEffNum_mode2 = (TH1F *)myFile->Get("not_GEN_REC_signal_mode2");
  TH1F* backgroundEffDen_mode2 = (TH1F *)myFile->Get("not_GEN_signal_mode2");
  TH1F* backgroundEff_mode2 = (TH1F*)backgroundEffNum_mode2->Clone("backgroundEff_mode2");
  backgroundEff_mode2->Divide(backgroundEffDen_mode2);
  TCanvas *cA17 = new TCanvas("cA17","Background_efficiency_mode2",200,10,600,400);
  cA17->SetFillColor(0);
  cA17->cd();
  backgroundEff_mode2->Draw();
  
  
  // -- Decay type events distribution -- //
  double Wtype_mode2 = histoFourLepInvariantMass_A->Integral(histoFourLepInvariantMass_A->FindBin(75),histoFourLepInvariantMass_A->FindBin(85));                
  double Ztype_mode2 = histoFourLepInvariantMass_A->Integral(histoFourLepInvariantMass_A->FindBin(88),histoFourLepInvariantMass_A->FindBin(140));               
  double WZtype_mode2 = histoFourLepInvariantMass_A->Integral(histoFourLepInvariantMass_A->FindBin(170),histoFourLepInvariantMass_A->FindBin(280));             
  
  cout << "#events_Wtype (mode 2) = " << Wtype_mode2 << ";   #events_Ztype (mode 2) = " << Ztype_mode2 << ";   #events_WZtype (mode 2) = " << WZtype_mode2 << endl;
  cout << "(#events_Ztype)/(#events_Wtype) (mode 2) = " << Ztype_mode2/Wtype_mode2 << ";   (#events_WZtype)/(#events_Wtype) (mode 2) = " << WZtype_mode2/Wtype_mode2 << endl;
  
  
  
  
  
}
