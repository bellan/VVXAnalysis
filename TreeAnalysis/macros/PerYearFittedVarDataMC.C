#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TVector.h>
#include <TLatex.h>
#include <iostream>
#include <vector>
#include <cmath>

double minRangeSens = 0.;
double maxRangeSens = 2.0;
bool isLogScale = true;
int aggrBin = 1;

TString resDir="resultsPostZXRwgtFullRun2";
//TString resDir="resultsFullRun2";
TString varFitted="BDTScore";
//TString year="2018";
const std::vector<TString> years = {"Run2"};
const std::vector<TString> varNames = {"dPhiZG_","PhotonMVAID_", "recoZDaughter0Pt_","recoZDaughter1Pt_","DR_gammaClosestLept_","mllG_","recoVMass_", "recoVDaughter0Pt_", "recoVDaughter1Pt_", "recoZMass_", "ptGamma_","DR_gammaClosestJet_", "j0_p(g)","j1_p(g)", "j0_p(uds)","j1_p(uds)"};//"VZGMVAScore_"// {/*"recoZPt_", "ptGamma_", "recoVDaughter0Pt_", "recoVDaughter1Pt_", "recoZDaughter0Pt_", "recoZDaughter1Pt_",*/"recoZPt_","dPhiZG_","VZGMVAScore_","j0_p(bX)","j1_p(bX)","MET_","PhotonMVAID_", "ZepCorr_","recoZDaughter0Pt_","recoZDaughter1Pt_","DR_gammaClosestLept_","mllG_","mjjG_", "recoVMass_", "recoVDaughter0Pt_", "recoVDaughter1Pt_", "recoZMass_", "ptGamma_","DR_gammaClosestJet_", "j0_p(g)","j1_p(g)", "j0_p(uds)","j1_p(uds)","FWM_T0_fullSyst_", "FWM_T1_fullSyst_","jj_p(uds)Sum","jj_p(g)Sum"};//{"recoVMass","recoVDaughter0Pt"};// {"### HAD TOPO DJ cand mass"};//,"recoVDaughter0Pt","DR_gammaClosestLept", "DR_Lept","FWM_T0_fullSyst", "FWM_T1_fullSyst","FWM_T2_fullSyst", "FWM_T3_fullSyst", "FWM_T4_fullSyst","FWM_T5_fullSyst", "FWM_T6_jets", "FWM_T1_jets","FWM_T2_jets", "FWM_T3_jets","FWM_T4_jets", "FWM_T5_jets", "FWM_T6_jets", "mllG",  "recoVMass", "DR_gammaClosestJet", "recoVDaughter1Pt", "recoVDaughtersDeltaPhi", "recoVPt", "recoZDeltaPhi", "recoZEta", "recoZMass", "recoZPt"};//, "System_Pt"};//, "relativePT_G_vs_V","VBH0s","VBH0z","VBH0t","jjH0s","jjH0z","jjH0t"};
const std::vector<TString> regions = {"CRZOFF","CRFSRT","CRZON_FSRT", "CRZOFF_FSRT","CRZOFF_DIB","CR2P_1VL"};//


void PerYearFittedVarDataMC()
{
  gStyle->SetOptStat(0);
  for(int i=0; i<1; i++){
    TString year=years[i];
    TString path_to_rootOutput= "../"+resDir+"/"+year+"/VZGAnalyzer_SR2P/";
    //  TFile *WZGFile = new TFile(path_to_rootOutput+"WZGTo2L2jG.root", "READ");
    //TFile *ZZGFile = new TFile(path_to_rootOutput+"ZZGTo2L2jG.root", "READ");
    TFile *VZGFile = new TFile(path_to_rootOutput+"VZG.root", "READ");
    TFile *DYFile = new TFile(path_to_rootOutput+"DYJetsToLL_M50.root", "READ");
    TFile *ZGFile = new TFile(path_to_rootOutput+"ZGToLLG.root", "READ");
    //    TFile *FSRFile = new TFile(path_to_rootOutput+"FSR.root", "READ");
    /*  TFile *ttXFile = new TFile(path_to_rootOutput+"ttX.root", "READ");
	TFile *TZqFile = new TFile(path_to_rootOutput+"TZq.root", "READ");
	TFile *WXFile = new TFile(path_to_rootOutput+"WX.root", "READ");
	TFile *qqZZFile = new TFile(path_to_rootOutput+"qqZZ.root", "READ");
    */
    TFile *TTTo2L2NuFile = new TFile(path_to_rootOutput+"TTTo2L2Nu.root", "READ");
    TFile *DataFile = new TFile(path_to_rootOutput+"data_obs.root", "READ");


    TString varLabel="BDT Score";  //varFitted;
    TString varName="SYS_"+varFitted+"_central";

    // Leggi gli istogrammi di segnale e background dai file root
    TH1F *signalHist = (TH1F*)VZGFile->Get(varName);
    signalHist->GetXaxis()->SetRangeUser(-0.08,1.);

    //      TH1F *DYPromptHist = (TH1F*)DYFile->Get(varName+"prompt"+cutNb);
    TH1F *DYNonPromptHist = (TH1F*)DYFile->Get(varName);
    double kFactor=1.;
    if(year=="2016preVFP") kFactor=1./1.40;///1.28;
    if(year=="2016postVFP") kFactor=1./1.39;///1.28;
    if(year=="2017") kFactor=1./1.24;///1.14;
    if(year=="2018") kFactor=1./1.19;///1.08;
    kFactor=1.;
    DYNonPromptHist->Scale(kFactor);
    DYNonPromptHist->GetXaxis()->SetRangeUser(-0.08,1.);
    TH1F *ZGPromptHist = (TH1F*)ZGFile->Get(varName);
    ZGPromptHist->GetXaxis()->SetRangeUser(-0.08,1.);
    //      TH1F *ZGNonPromptHist = (TH1F*)ZGFile->Get(varName+"nonPrompt"+cutNb);

    //  TH1F *VZFSRHist = (TH1F*)VZGFile->Get(varName);
    TH1F *TTTo2L2NuHist = (TH1F*)TTTo2L2NuFile->Get(varName);
    TTTo2L2NuHist->GetXaxis()->SetRangeUser(-0.08,1.);
    //    TH1F *FSRHist = (TH1F*)FSRFile->Get(varName);
    TH1F *dataHist = (TH1F*)DataFile->Get(varName); //SR UNBLINDED
    dataHist->GetXaxis()->SetRangeUser(-0.08,1.);


    //      DYPromptHist->Rebin(aggrBin);
    ZGPromptHist->Rebin(aggrBin);

    DYNonPromptHist->Rebin(aggrBin);
    //      ZGNonPromptHist->Rebin(aggrBin);

    TTTo2L2NuHist->Rebin(aggrBin);
    //  VZFSRHist->Rebin(aggrBin);
    signalHist->Rebin(aggrBin);
    //    FSRHist->Rebin(aggrBin);

    dataHist->Rebin(aggrBin);
    /*
      signalHist->SetFillColor(kRed);
      signalHist->SetLineColor(kRed+1);
    */
    DYNonPromptHist->SetFillColor(kAzure+1);
    DYNonPromptHist->SetLineColor(kAzure+5);
      
    //      DYPromptHist->SetFillColor(kWhite);
    //      DYPromptHist->SetLineColor(kAzure+5);
      
    dataHist->SetMarkerStyle(20);
    dataHist->SetMarkerSize(.9);
    dataHist->SetMarkerColor(kBlack);
    dataHist->SetLineColor(kBlack);


    signalHist->SetFillColorAlpha(kRed+2, 3.0); // Colore blu, semi-trasparente
    signalHist->SetFillStyle(3004);// Linee diagonali
    signalHist->SetLineWidth(2);
    signalHist->SetLineColor(kRed);// Spessore del contorno aumentato
      
    ZGPromptHist->SetFillColor(kOrange);
    ZGPromptHist->SetLineColor(kOrange+3);

    //      ZGNonPromptHist->SetFillColor(kOrange+5);
    //      ZGNonPromptHist->SetLineColor(kOrange+7);

    /*
      DYHist->SetFillColor(kAzure+1);
      DYHist->SetLineColor(kAzure+5);
      ZGHist->SetFillColor(kOrange);
      ZGHist->SetLineColor(kOrange+3);
    */
    TTTo2L2NuHist->SetFillColor(kBlue);
    TTTo2L2NuHist->SetLineColor(kBlue+3);
    //    FSRHist->SetFillColor(kSpring+5);
    //    FSRHist->SetLineColor(kSpring+3);
    /*
      if (!signalHist) {
      std::cout << "Errore nella lettura dell'istogramma di segnale" << std::endl;
      return;
      }
      if (!DYHist) {
      std::cout << "Errore nella lettura dell'istogramma di background" << std::endl;
      return;
      }
    */
    // Crea il canvas per il plot
    TCanvas *canvas = new TCanvas("canvas", "BDTScoreForFit", 800, 1000);
    
    // Crea il pad principale per lo stack plot
    TPad *mainPad = new TPad("mainPad", "Main Pad", 0.0, 0.3, 0.95, 1.0);
    if (isLogScale) mainPad->SetLogy(); // Imposta l'asse Y in scala logaritmica

    mainPad->SetBottomMargin(0.1);
    mainPad->Draw();

    // Crea il pad inferiore per il rapporto segnale/fondo
    TPad *ratioPad = new TPad("ratioPad", "Ratio Pad", 0.0, 0.0, 0.95, 0.3);
    ratioPad->SetTopMargin(0.05);
    ratioPad->SetBottomMargin(0.3);
    ratioPad->Draw();

    // Disegna lo stack plot nel pad principale
    mainPad->cd();

    TH1F *MCHist = (TH1F*)ZGPromptHist->Clone("MCHist");
    //      if(regions[iRegion]=="CR2P_1VL"){
    //      	MCHist->Add(DYPromptHist);
    MCHist->Add(DYNonPromptHist);
    //      }
    MCHist->Add(signalHist);
    //MCHist->Add(ZGNonPromptHist);

    //    MCHist->Add(FSRHist);
    MCHist->Add(TTTo2L2NuHist);
    //	MCHist->Add(tZqHist);

    MCHist->SetFillStyle(3005);
    MCHist->SetMarkerStyle(1);
    MCHist->SetFillColor(kBlack);
    MCHist->SetLineColor(kBlack);// Spessore del contorno aumentato
    MCHist->SetLineWidth(1.0);

	
      
    THStack *stack = new THStack("stack", "");
    //      stack->Add(backgroundHist);
    stack->Add(signalHist);

    //	stack->Add(ttXHist);
    //	stack->Add(tZqHist);
      
    //    stack->Add(FSRHist);
    stack->Add(TTTo2L2NuHist);
    stack->Add(ZGPromptHist);
    stack->Add(DYNonPromptHist);
    //stack->Add(ZGNonPromptHist);
    //stack->Add(DYPromptHist);
    //	stack->Add(ZGHist_FSR);


    stack->Draw("HIST"); // opzione "nostack" per sovrapporre gli istogrammi

    //stack->GetXaxis()->SetTitle(varName);
    stack->GetYaxis()->SetTitle("Events");
    stack->GetXaxis()->SetTitle(varLabel);

    //      stack->SetMinimum(1.5); // Imposta il valore minimo dell'asse Y logaritmico
    //if(cutNb>10)	  stack->SetMaximum(4000);

    if(isLogScale){
      stack->SetMinimum(1.1);
    }
    stack->SetMaximum(aggrBin*50000);

    //      stack->SetMaximum(10000);
    //      signalHist->Draw("HIST SAME");
    //signalHist->SetMinimum(2.1);
    //      dataHist->Draw("E SAME"); 
    dataHist->Draw("E SAME"); 
    MCHist->Draw("sameE2");
	
    TLegend *legend = new TLegend(0.86, 0.68, 0.99, 0.89);
    legend->AddEntry(signalHist, "VZ#gamma", "f");
    //      legend->AddEntry(DYPromptHist, "Drell-Yan + #gamma prompt)", "f");
    legend->AddEntry(DYNonPromptHist, "Drell-Yan", "f");
    legend->AddEntry(ZGPromptHist, "Z#gamma", "f");	
    //      legend->AddEntry(ZGNonPromptHist, "Z#gamma non-prompt", "f");	
    //	legend->AddEntry(ZGHist_FSR, "Z#gamma(FSR)", "f");
    //	legend->AddEntry(ttXHist, "qqZZ", "f");
    //    legend->AddEntry(FSRHist, "VZ+FSR", "f");
    legend->AddEntry(TTTo2L2NuHist, "tt #rightarrow 2l2#nu", "f");
    legend->AddEntry(MCHist, "Pred. unc.", "f");
    legend->AddEntry(dataHist, "Data", "lep");
    legend->Draw();

    TLatex *cmsText = new TLatex();
    cmsText->SetTextSize(0.05);
    cmsText->SetTextFont(62); // Font bold per "CMS"
    cmsText->SetTextAlign(13); // Allineamento al centro orizzontale (sinistra, alto)
    cmsText->DrawLatexNDC(0.1, 0.95, "CMS #bf{Private work}");

    // Aggiunta della luminosità in alto a destra
    TLatex *luminosityText = new TLatex();
    luminosityText->SetTextSize(0.04);
    luminosityText->SetTextFont(42); // Font normale per luminosità
    luminosityText->SetTextAlign(32); // Allineamento al centro orizzontale (destra, alto)
    /*
    if(year=="2016preVFP")  luminosityText->DrawLatexNDC(0.9, 0.93, "19.5 fb^{-1} (13 TeV)");
    if(year=="2016postVFP")  luminosityText->DrawLatexNDC(0.9, 0.93, "16.8 fb^{-1} (13 TeV)");
    if(year=="2017")  luminosityText->DrawLatexNDC(0.9, 0.93, "41.5 fb^{-1} (13 TeV)");
    if(year=="2018")  luminosityText->DrawLatexNDC(0.9, 0.93, "59.8 fb^{-1} (13 TeV)");
    */
    luminosityText->DrawLatexNDC(0.9, 0.93, "137.6 fb^{-1} (13 TeV)");
      
    TH1F *ratioHist = (TH1F*)dataHist->Clone("ratioHist");
      
    ratioHist->Divide(MCHist);
    ratioPad->cd();

    ratioHist->GetYaxis()->SetTitleSize(0.08);	
    ratioHist->GetXaxis()->SetTitleSize(0.08);	
    ratioHist->GetYaxis()->SetTitleOffset(0.5);	
    ratioHist->GetYaxis()->SetTitle("Data/MC");
    ratioHist->GetYaxis()->SetLabelSize(0.06);	
    ratioHist->GetXaxis()->SetLabelSize(0.06);	
    ratioHist->GetXaxis()->SetTitle(varLabel);
    ratioHist->SetTitle("");
    ratioHist->Draw("E");
    ratioHist->GetYaxis()->SetRangeUser(minRangeSens,maxRangeSens);

    ratioHist->SetMarkerStyle(20);
    ratioHist->SetMarkerSize(.9);
    ratioHist->SetMarkerColor(kBlack);
    ratioHist->SetLineColor(kBlack);

    TLine* line = new TLine(ratioHist->GetXaxis()->GetXmin(), 1.0, ratioHist->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(2);
    line->Draw("same");
	
    /*
      TH1F *ratioHist = (TH1F*)DYNonPromptHist->Clone("ratioHist");
      ratioHist->Scale(1000000);
      ratioPad->cd();
      ratioHist->Draw("same");
      double xMax = ratioHist->GetXaxis()->GetXmax();
      double xMin = ratioHist->GetXaxis()->GetXmin();
      TText *blindText = new TText( (xMin + xMax)/2 - (xMax - xMin)/8  , 1.1, "BLINDED");
      blindText->SetTextSize(0.12);
      blindText->Draw("same");
	
      ratioHist->GetYaxis()->SetTitleSize(0.08);
      ratioHist->GetXaxis()->SetTitleSize(0.08);
      ratioHist->GetYaxis()->SetTitleOffset(0.5);
      ratioHist->GetYaxis()->SetTitle("Data/MC");
      ratioHist->GetYaxis()->SetLabelSize(0.06);	
      ratioHist->GetXaxis()->SetLabelSize(0.06);	
      ratioHist->GetXaxis()->SetTitle(varLabel);
      ratioHist->SetTitle("");
      ratioHist->GetYaxis()->SetRangeUser(minRangeSens,maxRangeSens);

      TLine* line = new TLine(ratioHist->GetXaxis()->GetXmin(), 1.0, ratioHist->GetXaxis()->GetXmax(), 1.0);
      line->SetLineStyle(2);
      line->Draw("same");
    */
    /*
    if (isLogScale){
      canvas->SaveAs("../DistributionPlots/LogScale/"+resDir+"/"+varFitted+"ForFit_"+year+".png");
      canvas->SaveAs("../DistributionPlots/LogScale/"+resDir+"/"+varFitted+"ForFit_"+year+".pdf");
    }
    else {
      canvas->SaveAs("../DistributionPlots/LinScale/"+resDir+"/"+varFitted+"ForFit_"+year+".png");
      canvas->SaveAs("../DistributionPlots/LinScale/"+resDir+"/"+varFitted+"ForFit_"+year+".pdf");
    }
    */

    canvas->SaveAs("../DistributionPlots/LogScale/"+resDir+"/"+varFitted+"ForFit_Run2.png");
    canvas->SaveAs("../DistributionPlots/LogScale/"+resDir+"/"+varFitted+"ForFit_Run2.pdf");

    // Chiudi i file root
    VZGFile->Close();
    //    FSRFile->Close();
    TTTo2L2NuFile->Close();
    DYFile->Close();
    ZGFile->Close();
    DataFile->Close();
    //    ttXFile->Close();
    //    tZqFile->Close();
    //ZGEnrichFile->Close();
  }
  
}
