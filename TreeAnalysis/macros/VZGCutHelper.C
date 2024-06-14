#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TVector.h>

double minRangeSens = 1.8;
double maxRangeSens = 2.0;

const std::vector<TString> varNames = { "recoVMass", "FWM_T0_fullSyst", "DR_gammaClosestJet", "FWM_T2_fullSyst", "recoVDaughter0Pt", "recoVDaughter1Pt", "recoVDaughtersDeltaPhi", "recoVPt", "recoZDeltaPhi", "recoZEta", "recoZMass", "recoZPt"};//, "System_Pt"};//, "relativePT_G_vs_V","VBH0s","VBH0z","VBH0t","jjH0s","jjH0z","jjH0t"};

void VZGCutHelper()
{
  // Apri i file root contenenti gli istogrammi di segnale e background
  TString path_to_rootOutput= "../results/2016preVFP/VZGAnalyzer_SR2P/";
  TFile *WZGFile = new TFile(path_to_rootOutput+"WZGTo2L2jG.root", "READ");
  TFile *ZZGFile = new TFile(path_to_rootOutput+"ZZGTo2L2jG.root", "READ");
  TFile *backgroundFile = new TFile(path_to_rootOutput+"DYJetsToLL_M50.root", "READ");
  
  for (int cutNb = 7; cutNb < 8; cutNb++){
    for (int iVars=0; iVars<varNames.size(); ++iVars){
      TString varName=varNames[iVars];
      // Leggi gli istogrammi di segnale e background dai file root
      TH1F *signalHist = (TH1F*)WZGFile->Get(varName+"_sign"+cutNb);
      TH1F *ZZGHist = (TH1F*)ZZGFile->Get(varName+"_sign"+cutNb);
      signalHist->Add(ZZGHist);
      TH1F *backgroundHist = (TH1F*)backgroundFile->Get(varName+"_bckg"+cutNb);

      signalHist->SetFillColor(kRed);
      backgroundHist->SetFillColor(kBlue);
      // Verifica che gli istogrammi siano stati letti correttamente
      if (!signalHist) {
        std::cout << "Errore nella lettura dell'istogramma di segnale" << std::endl;
        return;
      }
      if (!backgroundHist) {
        std::cout << "Errore nella lettura dell'istogramma di background" << std::endl;
        return;
      }
      std::cout << "__" << std::endl;
      std::cout << "Histos exist" << std::endl;
      std::cout << "__" << std::endl;

      // Crea il canvas per il plot
      TCanvas *canvas = new TCanvas("canvas", varNames.at(iVars)+"  after cut "+cutNb, 800, 600);
    
      // Crea il pad principale per lo stack plot
      TPad *mainPad = new TPad("mainPad", "Main Pad", 0.0, 0.3, 1.0, 1.0);
      mainPad->SetLogy(); // Imposta l'asse Y in scala logaritmica
      mainPad->SetBottomMargin(0.1);
      mainPad->Draw();

      // Crea il pad inferiore per il rapporto segnale/fondo
      TPad *ratioPad = new TPad("ratioPad", "Ratio Pad", 0.0, 0.0, 1.0, 0.3);
      ratioPad->SetTopMargin(0.05);
      ratioPad->SetBottomMargin(0.3);
      ratioPad->Draw();

      // Disegna lo stack plot nel pad principale
      mainPad->cd();
      THStack *stack = new THStack("stack", varNames.at(iVars)+"  after cut "+cutNb);
      stack->Add(signalHist);
      stack->Add(backgroundHist);
      stack->Draw("HIST"); // opzione "nostack" per sovrapporre gli istogrammi
      stack->GetXaxis()->SetTitle(varName);
      stack->GetYaxis()->SetTitle("Events");
      stack->GetYaxis()->SetRange(1.0,1000);
      stack->SetMinimum(0.1); // Imposta il valore minimo dell'asse Y logaritmico
      TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
      legend->AddEntry(backgroundHist, "Background", "f");
      legend->AddEntry(signalHist, "Signal", "f");
      legend->Draw();

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
	  std::cout << "Bin" << i  << "Den" << sqrtContent <<std::endl;
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
	  std::cout << "Bin" << i  << "Den" << sqrtContent <<std::endl;
	  sqrtTotalCumulativeBWD->SetBinContent(i, sqrtContent);
	}


      std::cout << "__" << std::endl;
      std::cout << "sqrtTotal cumulative filled" << std::endl;
      std::cout << "__" << std::endl;


      // Calcola il rapporto segnale/fondo bin per bin
      TH1F *ratioHist = (TH1F*)signalHist->Clone("ratioHist");
      ratioHist->Divide(sqrtTotal);
      /*for (Int_t i=1; i<=ratioHist->GetNbinsX(); i++)
	{
        std::cout << "Bin" << i  << "Num" << signalCumulativeFWD->GetBinContent(i) <<std::endl;
        std::cout << "Bin" << i  << "Sens " << ratioHistFWD->GetBinContent(i) <<std::endl;

	}*/
      // Calcola il rapporto segnale/fondo bin per bin
      TH1F *ratioHistFWD = (TH1F*)signalCumulativeFWD->Clone("ratioHist");
      ratioHistFWD->Divide(sqrtTotalCumulativeFWD);


      std::cout << "__" << std::endl;
      std::cout << "ratio issued" << std::endl;
      std::cout << "__" << std::endl;

      for (Int_t i=1; i<=ratioHistFWD->GetNbinsX(); i++)
	{
	  std::cout << "Bin" << i  << "Num" << signalCumulativeFWD->GetBinContent(i) <<std::endl;
	  std::cout << "Bin" << i  << "Sens " << ratioHistFWD->GetBinContent(i) <<std::endl;

	}
      // Calcola il rapporto segnale/fondo bin per bin
      TH1F *ratioHistBWD = (TH1F*)signalCumulativeBWD->Clone("ratioHist");
      ratioHistBWD->Divide(sqrtTotalCumulativeBWD);
      for (Int_t i=1; i<=ratioHistBWD->GetNbinsX(); i++)
	{
	  std::cout << "Bin" << i  << "Num" << signalCumulativeBWD->GetBinContent(i) <<std::endl;
	  std::cout << "Bin" << i  << "Sens " << ratioHistBWD->GetBinContent(i) <<std::endl;

	}
    

      // Disegna il rapporto segnale/fondo nel pad inferiore
      ratioPad->cd();

      ratioHist->SetStats(0);
      ratioHist->SetTitle("");
      ratioHist->GetXaxis()->SetTitle(varName);
      ratioHist->GetYaxis()->SetTitle("sensitivity");
      ratioHist->GetYaxis()->SetRangeUser(minRangeSens,maxRangeSens);
      ratioHist->SetMarkerStyle(22);
      ratioHist->SetLineColor(kBlack);
      ratioHist->Draw("SAMES");

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


      // Salva il canvas in un file
      //canvas->Print("stack_plot.png");
      canvas->SaveAs(path_to_rootOutput+varName+cutNb+".png");
    }
  }
  // Chiudi i file root
  WZGFile->Close();
  ZZGFile->Close();
  backgroundFile->Close();
}
