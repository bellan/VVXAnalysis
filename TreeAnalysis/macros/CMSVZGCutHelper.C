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
double maxRangeSens = 1.;
bool isLogScale = true;
int aggrBin = 1;

const std::vector<TString> varNames = {"mllG_","recoZMass_"};//{"DR_gammaClosestLept_","recoVMass_","mllG_","j1_p(g)","PhotonMVAID_", "j1_Girth", "recoVDaughter1Pt_", "j0_p(uds)", "recoVPt_", "H_T\ ","j0_Girth","dRJ0Gamma_","dRJ1Gamma_", "ptGamma_","dPhiZG_", "recoZMass_"};

//, {"VZGMVAScore_"};// {, "recoZDaughter0Pt_","recoZDaughter1Pt_","mllG_","DR_gammaClosestJet_", };// {/*"recoZPt_", "ptGamma_", "recoVDaughter0Pt_", "recoVDaughter1Pt_", "recoZDaughter0Pt_", "recoZDaughter1Pt_",*/"recoZPt_","dPhiZG_","VZGMVAScore_","j0_p(bX)","j1_p(bX)","MET_","PhotonMVAID_", "ZepCorr_","recoZDaughter0Pt_","recoZDaughter1Pt_","DR_gammaClosestLept_","mllG_","mjjG_", "recoVMass_", "recoVDaughter0Pt_", "recoVDaughter1Pt_", "recoZMass_", "ptGamma_","DR_gammaClosestJet_", "j0_p(g)","j1_p(g)", "j0_p(uds)","j1_p(uds)","FWM_T0_fullSyst_", "FWM_T1_fullSyst_","jj_p(uds)Sum","jj_p(g)Sum"};//{"recoVMass","recoVDaughter0Pt"};// {"### HAD TOPO DJ cand mass"};//,"recoVDaughter0Pt","DR_gammaClosestLept", "DR_Lept","FWM_T0_fullSyst", "FWM_T1_fullSyst","FWM_T2_fullSyst", "FWM_T3_fullSyst", "FWM_T4_fullSyst","FWM_T5_fullSyst", "FWM_T6_jets", "FWM_T1_jets","FWM_T2_jets", "FWM_T3_jets","FWM_T4_jets", "FWM_T5_jets", "FWM_T6_jets", "mllG",  "recoVMass", "DR_gammaClosestJet", "recoVDaughter1Pt", "recoVDaughtersDeltaPhi", "recoVPt", "recoZDeltaPhi", "recoZEta", "recoZMass", "recoZPt"};//, "System_Pt"};//, "relativePT_G_vs_V","VBH0s","VBH0z","VBH0t","jjH0s","jjH0z","jjH0t"};
const std::vector<TString> regions = {"CRZON_FSRT", "CRZOFF", "CR2P_1VL","CRFSRT", "CRZOFF_FSRT"};//


void CMSVZGCutHelper()
{
  gStyle->SetOptStat(0);

  TString path_to_rootOutput= "../results_withLepUnc/Run2/VZGAnalyzer_SR2P/";
  //  TFile *WZGFile = new TFile(path_to_rootOutput+"WZGTo2L2jG.root", "READ");
  //TFile *ZZGFile = new TFile(path_to_rootOutput+"ZZGTo2L2jG.root", "READ");
  TFile *VZGFile = new TFile(path_to_rootOutput+"VZG.root", "READ");
  TFile *DYFile = new TFile(path_to_rootOutput+"DYJetsToLL_M50.root", "READ");
  TFile *ZGFile = new TFile(path_to_rootOutput+"ZGToLLG.root", "READ");
  /*  TFile *ttXFile = new TFile(path_to_rootOutput+"ttX.root", "READ");
  TFile *TZqFile = new TFile(path_to_rootOutput+"TZq.root", "READ");
  TFile *WXFile = new TFile(path_to_rootOutput+"WX.root", "READ");
  TFile *qqZZFile = new TFile(path_to_rootOutput+"qqZZ.root", "READ");
  */
  TFile *TTTo2L2NuFile = new TFile(path_to_rootOutput+"TTTo2L2Nu.root", "READ");
  TFile *DataFile = new TFile(path_to_rootOutput+"data_obs.root", "READ");
  /*
      if (isLogScale) {
	canvas->SaveAs("../DistributionPlots/Run2/LogScale/resultsForANv3plus_"+varName+"_"+regions[iRegion]+".png");
	canvas->SaveAs("../DistributionPlots/Run2/LogScale/resultsForANv3plus_"+varName+"_"+regions[iRegion]+".pdf");
      }/*
      else {
	canvas->SaveAs(path_to_rootOutput+"LinScale/"+varName+"_"+regions[iRegion]+".png");
	canvas->SaveAs(path_to_rootOutput+"LinScale/"+varName+"_"+regions[iRegion]+".pdf");
	}*/
  
  //Now for SR
  for (int cutNb = 4; cutNb < 5; cutNb++){
    for (int iVars=0; iVars<varNames.size(); ++iVars){
      TString varName=varNames[iVars];
      TString varLabel=varNames[iVars];
      if(varName=="mllG_") varLabel="m_{ll#gamma} [GeV]";
      if(varName=="VZGMVAScore_") varLabel="BDT Score";
      if(varName=="recoZMass_") varLabel="m_{ll} [GeV]";
      if(varName=="recoVMass_") varLabel="m_{jj} [GeV]";
      if(varName=="recoVDaughter0Pt_") varLabel="p_{T}^{j0} [GeV]";
      if(varName=="recoVDaughter1Pt_") varLabel="p_{T}^{j1} [GeV]";
      if(varName=="recoZDaughter0Pt_") varLabel="p_{T}^{l0} [GeV]";
      if(varName=="recoZDaughter1Pt_") varLabel="p_{T}^{l1} [GeV]";
      if(varName=="DR_gammaClosestLept_") varLabel="#Delta R_{l#gamma}";
      if(varName=="DR_gammaClosestJet_") varLabel="#Delta R_{j#gamma}";
      if(varName=="dPhiZG_") varLabel="#Delta #phi_{Z#gamma}";

      // Leggi gli istogrammi di segnale e background dai file root
      TH1F *signalHist = (TH1F*)VZGFile->Get(varName+"sign"+cutNb);
      TH1F *backgroundHist = (TH1F*)DYFile->Get(varName+"nonPrompt"+cutNb);

      //      TH1F *DYPromptHist = (TH1F*)DYFile->Get(varName+"prompt"+cutNb);
      TH1F *DYNonPromptHist = (TH1F*)DYFile->Get(varName+"nonPrompt"+cutNb);
      TH1F *ZGPromptHist = (TH1F*)ZGFile->Get(varName+"prompt"+cutNb);
      //      TH1F *ZGNonPromptHist = (TH1F*)ZGFile->Get(varName+"nonPrompt"+cutNb);

      TH1F *VZFSRHist = (TH1F*)VZGFile->Get(varName+"bckg"+cutNb);
      TH1F *TTTo2L2NuHist = (TH1F*)TTTo2L2NuFile->Get(varName+"all"+cutNb);
      //      TH1F *dataHist = (TH1F*)DataFile->Get(varName+"all"+cutNb); //SR UNBLINDED

      backgroundHist->Add(ZGPromptHist);
      backgroundHist->Add(TTTo2L2NuHist);
      backgroundHist->Add(VZFSRHist);


      //      DYPromptHist->Rebin(aggrBin);
      ZGPromptHist->Rebin(aggrBin);

      DYNonPromptHist->Rebin(aggrBin);
      //      ZGNonPromptHist->Rebin(aggrBin);

      TTTo2L2NuHist->Rebin(aggrBin);
      VZFSRHist->Rebin(aggrBin);
      signalHist->Rebin(aggrBin);
      backgroundHist->Rebin(aggrBin);
      /*
      signalHist->SetFillColor(kRed);
      signalHist->SetLineColor(kRed+1);
      */
      DYNonPromptHist->SetFillColor(kAzure+1);
      DYNonPromptHist->SetLineColor(kAzure+5);
      
      //      DYPromptHist->SetFillColor(kWhite);
      //      DYPromptHist->SetLineColor(kAzure+5);
      
      //      dataHist->Rebin(aggrBin);
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
      VZFSRHist->SetFillColor(kSpring+5);
      VZFSRHist->SetLineColor(kSpring+3);
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
      TCanvas *canvas = new TCanvas("canvas", varNames.at(iVars)+" after cut "+cutNb, 800, 1000);
    
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

      MCHist->Add(VZFSRHist);
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
      
      stack->Add(VZFSRHist);
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
      stack->SetMaximum(aggrBin*10000);

      //      stack->SetMaximum(10000);
      //      signalHist->Draw("HIST SAME");
      //signalHist->SetMinimum(2.1);
      //      dataHist->Draw("E SAME"); 
      MCHist->Draw("sameE2");
	
      TLegend *legend = new TLegend(0.86, 0.67, 0.99, 0.89);
      legend->AddEntry(signalHist, "VZ#gamma", "f");
      //      legend->AddEntry(DYPromptHist, "Drell-Yan + #gamma prompt)", "f");
      legend->AddEntry(DYNonPromptHist, "Drell-Yan", "f");
      legend->AddEntry(ZGPromptHist, "Z#gamma", "f");	
      //      legend->AddEntry(ZGNonPromptHist, "Z#gamma non-prompt", "f");	
      //	legend->AddEntry(ZGHist_FSR, "Z#gamma(FSR)", "f");
      //	legend->AddEntry(ttXHist, "qqZZ", "f");
      legend->AddEntry(VZFSRHist, "VZ+FSR", "f");
      legend->AddEntry(TTTo2L2NuHist, "tt #rightarrow 2l2#nu", "f");
      legend->AddEntry(MCHist, "Pred. unc.", "f");
      legend->Draw();

      TLatex *cmsText = new TLatex();
      cmsText->SetTextSize(0.05);
      cmsText->SetTextFont(62); // Font bold per "CMS"
      cmsText->SetTextAlign(13); // Allineamento al centro orizzontale (sinistra, alto)
      cmsText->DrawLatexNDC(0.1, 0.95, "CMS #bf{Work in progress}");

      // Aggiunta della luminosità in alto a destra
      TLatex *luminosityText = new TLatex();
      luminosityText->SetTextSize(0.04);
      luminosityText->SetTextFont(42); // Font normale per luminosità
      luminosityText->SetTextAlign(32); // Allineamento al centro orizzontale (destra, alto)
      luminosityText->DrawLatexNDC(0.9, 0.93, "137.6 fb^{-1} (13 TeV)");


      // Calcola l'istogramma cumulato di segnale
      TH1F *signalCumulativeFWD = (TH1F*)signalHist->GetCumulative();
      TH1F *backgroundCumulativeFWD = (TH1F*)backgroundHist->GetCumulative();
      TH1F *signalCumulativeBWD = (TH1F*)signalHist->GetCumulative(false);
      TH1F *backgroundCumulativeBWD = (TH1F*)backgroundHist->GetCumulative(false);


      std::cout << "__" << std::endl;
      std::cout << "CumulativeHistos created" << std::endl;
      std::cout << "__" << std::endl;



      // Calcola l'istogramma cumulato di segnale + fondo
      TH1F *total = (TH1F*)signalHist->Clone("total");
      total->Add(backgroundHist);
      TH1F *totalCumulativeFWD = (TH1F*)signalCumulativeFWD->Clone("totalCumulativeFWD");
      totalCumulativeFWD->Add(backgroundCumulativeFWD);
      TH1F *totalCumulativeBWD = (TH1F*)signalCumulativeBWD->Clone("totalCumulativeBWD");
      totalCumulativeBWD->Add(backgroundCumulativeBWD);


      std::cout << "__" << std::endl;
      std::cout << "Total Histos created" << std::endl;
      std::cout << "__" << std::endl;

      TH1F *sqrtTotal = (TH1F*)total->Clone("sqrtTotal");
      sqrtTotal->Reset();
      for (Int_t i=1; i<=total->GetNbinsX(); i++)
	{
	  Double_t binContent = total->GetBinContent(i);
	  Double_t sqrtContent;
	  if(binContent<=0) sqrtContent =0.000001;
	  else  sqrtContent = TMath::Sqrt(binContent);
	  std::cout << "Bin" << i  << "Den" << sqrtContent <<std::endl;
	  sqrtTotal->SetBinContent(i, sqrtContent);
	}

      std::cout << "__" << std::endl;
      std::cout << "sqrtTotal reset" << std::endl;
      std::cout << "__" << std::endl;

      // Calcola la radice quadrata della somma degli istogrammi cumulati di segnale e fondo
      TH1F *sqrtTotalCumulativeFWD = (TH1F*)totalCumulativeFWD->Clone("sqrtTotalCumulativeFWD");
      sqrtTotalCumulativeFWD->Reset();
      for (Int_t i=1; i<=totalCumulativeFWD->GetNbinsX(); i++)
	{
	  Double_t binContent = totalCumulativeFWD->GetBinContent(i);
	  Double_t sqrtContent;
	  if(binContent==0) sqrtContent =0.000001;
	  else  sqrtContent = TMath::Sqrt(binContent);
	  //	  std::cout << "Bin" << i  << "Den" << sqrtContent <<std::endl;
	  sqrtTotalCumulativeFWD->SetBinContent(i, sqrtContent);
	}
      TH1F *sqrtTotalCumulativeBWD = (TH1F*)totalCumulativeBWD->Clone("sqrtTotalCumulativeBWD");
      sqrtTotalCumulativeBWD->Reset();
      for (Int_t i=1; i<=totalCumulativeBWD->GetNbinsX(); i++)
	{
	  Double_t binContent = totalCumulativeBWD->GetBinContent(i);
	  Double_t sqrtContent;
	  if(binContent==0) sqrtContent =0.000001;
	  else  sqrtContent = TMath::Sqrt(binContent);
	  //	  std::cout << "Bin" << i  << "Den" << sqrtContent <<std::endl;
	  sqrtTotalCumulativeBWD->SetBinContent(i, sqrtContent);
	}


      std::cout << "__" << std::endl;
      std::cout << "sqrtTotal cumulative filled" << std::endl;
      std::cout << "__" << std::endl;


      // Calcola il rapporto segnale/fondo bin per bin
      TH1F *ratioHist = (TH1F*)signalHist->Clone("ratioHist");
      ratioHist->Divide(sqrtTotal);
      //ratioHist->Divide(total);//CT: modified to plot sig/bkg
      /*      
      TH1F *ratioHist = (TH1F*)DYNonPromptHist->Clone("ratioHist");
      ratioHist->Scale(1000000);
      ratioPad->cd();
      ratioHist->Draw("same");
      *//*      
      double xMax = ratioHist->GetXaxis()->GetXmax();
      double xMin = ratioHist->GetXaxis()->GetXmin();
      TText *blindText = new TText( (xMin + xMax)/2 - (xMax - xMin)/8  , 1.1, "BLINDED");
      blindText->SetTextSize(0.12);
      blindText->Draw("same");
      */
      TH1F *ratioHistFWD = (TH1F*)signalCumulativeFWD->Clone("ratioHistFWD");
      ratioHistFWD->Divide(sqrtTotalCumulativeFWD);
      TH1F *ratioHistBWD = (TH1F*)signalCumulativeBWD->Clone("ratioHistBWD");
      ratioHistBWD->Divide(sqrtTotalCumulativeBWD);
      
      // Disegna il rapporto segnale/fondo nel pad inferiore
      ratioPad->cd();

      ratioHist->GetYaxis()->SetTitleSize(0.08);
      ratioHist->GetXaxis()->SetTitleSize(0.08);
      ratioHist->GetYaxis()->SetTitleOffset(0.5);
      ratioHist->GetYaxis()->SetTitle("sig./#sqrt{tot.}");
      ratioHist->GetYaxis()->SetLabelSize(0.06);	
      ratioHist->GetXaxis()->SetLabelSize(0.06);	
      ratioHist->GetXaxis()->SetTitle(varLabel);
      ratioHist->SetTitle("");
      ratioHist->GetYaxis()->SetRangeUser(minRangeSens,maxRangeSens);
      
      ratioHist->SetStats(0);
      ratioHist->SetTitle("");
      ratioHist->GetXaxis()->SetTitle(varLabel);
      //      ratioHist->GetYaxis()->SetTitle("sensitivity");
      ratioHist->GetYaxis()->SetRangeUser(minRangeSens,maxRangeSens);
      ratioHist->SetMarkerStyle(8);
      ratioHist->SetLineColor(kBlack);
      ratioHist->Draw("SAMES");
      /*      
      ratioHistFWD->SetStats(0);
      ratioHistFWD->SetTitle("");
      ratioHistFWD->GetYaxis()->SetRangeUser(minRangeSens,maxRangeSens);
      ratioHistFWD->SetMarkerStyle(22);
      ratioHistFWD->SetLineColor(kBlack);
      ratioHistFWD->SetLineColor(kWhite);
      ratioHistFWD->Draw("SAMES");
      
      ratioHistBWD->SetStats(0);
      ratioHistBWD->SetTitle("");
      ratioHistBWD->GetYaxis()->SetRangeUser(minRangeSens,maxRangeSens);
      ratioHistBWD->SetMarkerStyle(18);
      ratioHistBWD->SetLineColor(kBlue);
      ratioHistFWD->SetLineColor(kWhite);
      ratioHistBWD->Draw("SAMES");
      */
      TLine* line = new TLine(ratioHist->GetXaxis()->GetXmin(), 1.0, ratioHist->GetXaxis()->GetXmax(), 1.0);
      line->SetLineStyle(2);
      line->Draw("same");

      if (isLogScale){
	canvas->SaveAs("../DistributionPlots/Run2/LogScale/SoverBForANv3plus_"+varName+"_SR"+cutNb+".png");
	canvas->SaveAs("../DistributionPlots/Run2/LogScale/SoverBForANv3plus_"+varName+"_SR"+cutNb+".pdf");
      }
      /*	
      if (isLogScale){
	canvas->SaveAs(path_to_rootOutput+"LogScale/"+varName+"_"+cutNb+".png");
	canvas->SaveAs(path_to_rootOutput+"LogScale/"+varName+"_"+cutNb+".pdf");
      }
      else {
	canvas->SaveAs(path_to_rootOutput+"LinScale/"+varName+"_"+cutNb+".png");
	canvas->SaveAs(path_to_rootOutput+"LinScale/"+varName+"_"+cutNb+".pdf");
      }
      */
    }
  }
  
  // Chiudi i file root
  VZGFile->Close();
  //  ZZGFile->Close();
  TTTo2L2NuFile->Close();
  DYFile->Close();
  ZGFile->Close();
  DataFile->Close();
  //    ttXFile->Close();
  //    tZqFile->Close();
  //ZGEnrichFile->Close();

  
}
