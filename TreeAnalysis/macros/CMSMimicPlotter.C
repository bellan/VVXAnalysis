#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TVector.h>
#include <TLatex.h>

double minRangeSens = 1.5;
double maxRangeSens = 2.5;
int aggrBin=1;
bool isLogScale = true;
double sigZoom = 50.;


const std::vector<TString> varNames = {"recoZMass_"};//,"dPhiZG_", "j0_p(bX)","j1_p(bX)","MET_","PhotonMVAID_", "ZepCorr_","recoZDaughter0Pt_","recoZDaughter1Pt_","DR_gammaClosestLept_","mllG_","mjjG_", "recoVMass_", "recoVDaughter0Pt_", "recoVDaughter1Pt_", "recoZMass_", "ptGamma_","DR_gammaClosestJet_", "j0_p(g)","j1_p(g)", "j0_p(uds)","j1_p(uds)","FWM_T0_fullSyst_", "FWM_T1_fullSyst_","FWM_T2_fullSyst_","jj_p(uds)Sum","jj_p(g)Sum"};//{"recoVMass","recoVDaughter0Pt"};// {"### HAD TOPO DJ cand mass"};//,"recoVDaughter0Pt","DR_gammaClosestLept", "DR_Lept","FWM_T0_fullSyst", "FWM_T1_fullSyst","FWM_T2_fullSyst", "FWM_T3_fullSyst", "FWM_T4_fullSyst","FWM_T5_fullSyst", "FWM_T6_jets", "FWM_T1_jets","FWM_T2_jets", "FWM_T3_jets","FWM_T4_jets", "FWM_T5_jets", "FWM_T6_jets", "mllG",  "recoVMass", "DR_gammaClosestJet", "recoVDaughter1Pt", "recoVDaughtersDeltaPhi", "recoVPt", "recoZDeltaPhi", "recoZEta", "recoZMass", "recoZPt"};//, "System_Pt"};//, "relativePT_G_vs_V","VBH0s","VBH0z","VBH0t","jjH0s","jjH0z","jjH0t"};

void CMSMimicPlotter()
{
  double scaleLumi2017 = 1.;//7.035--> 2016post||3.32-->2017||2.29-->2018||8.16--> 2016post||7.035--> 2016pre  
   //double scaleLumi2016postVFP = 7.035;//7.035--> 2016post||3.32-->2017||2.29-->2018||8.16--> 2016post||7.035--> 2016pre  
  double scaleLumi2018 = 1.;//7.035--> 2016post||3.32-->2017||2.29-->2018||8.16--> 2016post||7.035--> 2016pre  

  //  double scaleStat = 4.646;
  // Apri i file root contenenti gli istogrammi di segnale e background
  //  TString path_to_rootOutput= "../results/2016preVFP/VZGAnalyzer_SR2P/";
  //  TString path_to_rootOutput= "../results/Run2/SR_cutBased/";
  TString path_to_rootOutput= "../resultsForThesisPlotsFullRun2/Run2/VZGAnalyzer_SR2P/";
  //  TFile *WZGFile = new TFile(path_to_rootOutput+"WZGTo2L2jG.root", "READ");
  //TFile *ZZGFile = new TFile(path_to_rootOutput+"ZZGTo2L2jG.root", "READ");
  TFile *VZGFile = new TFile(path_to_rootOutput+"VZG.root", "READ");
  TFile *DYFile = new TFile(path_to_rootOutput+"DYJetsToLL_M50.root", "READ");
  TFile *ZGFile = new TFile(path_to_rootOutput+"ZGToLLG.root", "READ");
  TFile *TTTo2L2NuFile = new TFile(path_to_rootOutput+"TTTo2L2Nu.root", "READ");
  /*
  int NJ=3;
  for(int VHCandRank==0; VHCandRank < 3; VHCandRank ++){
    TH1F *signalHist = (TH1F*)VZGFile->Get(NJ+"J_DJCand"+VHCandRank+"_jj_QGL");
    TH1F *backgroundHist = (TH1F*)DYFile->Get(NJ+"J_DJCand"+VHCandRank+"_jj_QGL");
    TH1F *ZGHist = (TH1F*)ZGFile->Get(NJ+"J_DJCand"+VHCandRank+"_jj_QGL");
  //  TH1F *VZFSRHist = (TH1F*)VZGFile->Get(varName+"bckg"+cutNb);
  }

  VZGFile->Close();
  //  ZZGFile->Close();
  DYFile->Close();
  ZGFile->Close();
  return;
  */
  for (int cutNb = 3; cutNb < 5; cutNb++){
    /*    if(cutNb>1)
      {
	minRangeSens = 0.8;
	maxRangeSens = 1.7;
      }
    */
    for (int iVars=0; iVars<varNames.size(); ++iVars){
      TString varName=varNames[iVars];
      // Leggi gli istogrammi di segnale e background dai file root
      TH1F *signalHist = (TH1F*)VZGFile->Get(varName+"sign"+cutNb);

      TH1F *DYPromptHist = (TH1F*)DYFile->Get(varName+"prompt"+cutNb);
      TH1F *DYNonPromptHist = (TH1F*)DYFile->Get(varName+"nonPrompt"+cutNb);
      TH1F *ZGPromptHist = (TH1F*)ZGFile->Get(varName+"prompt"+cutNb);
      TH1F *ZGNonPromptHist = (TH1F*)ZGFile->Get(varName+"nonPrompt"+cutNb);

      TH1F *VZFSRHist = (TH1F*)VZGFile->Get(varName+"bckg"+cutNb);
      TH1F *TTTo2L2NuHist = (TH1F*)TTTo2L2NuFile->Get(varName+"all"+cutNb);

      //      TH1F *ZZGHist = (TH1F*)ZZGFile->Get(varName+"sign"+cutNb);
      //      signalHist->Add(ZZGHist);
      //      std::cout << "printing var " << varName+"sign"+cutNb << std::endl;
      if(!isLogScale)signalHist->Scale(sigZoom*scaleLumi2018);

      //      std::cout << varNames[iVars] << std::endl;
      
      ZGPromptHist->Scale(scaleLumi2018);
      ZGNonPromptHist->Scale(scaleLumi2018);
      DYPromptHist->Scale(scaleLumi2018);
      DYNonPromptHist->Scale(scaleLumi2018);
      TTTo2L2NuHist->Scale(scaleLumi2018);
      VZFSRHist->Scale(scaleLumi2018);

      signalHist->Rebin(aggrBin);
      ZGPromptHist->Rebin(aggrBin);
      ZGNonPromptHist->Rebin(aggrBin);
      DYPromptHist->Rebin(aggrBin);
      DYNonPromptHist->Rebin(aggrBin);
      TTTo2L2NuHist->Rebin(aggrBin);
      VZFSRHist->Rebin(aggrBin);

      /*      
      if(varName=="mllG_" || varName=="recoVMass_" || varName=="VZGMVAScore_"){
	backgroundHist->Rebin(aggrBin);
	ZGHist->Rebin(aggrBin);
	VZFSRHist->Rebin(aggrBin);
	signalHist->Rebin(aggrBin); 
      }
      if(varName=="jj_p(g)Sum" || varName=="j0_p(g)" || varName=="j1_p(g)"){
	backgroundHist->Rebin(aggrBin);
	ZGHist->Rebin(aggrBin);
	VZFSRHist->Rebin(aggrBin);
	signalHist->Rebin(aggrBin); 
      }
      */


      //      stack->GetXaxis()->SetRangeUser(60,120);

      /*      if(varName=="recoVMass_"){
	backgroundHist->GetXaxis()->SetRangeUser(60,120);
	ZGHist->GetXaxis()->SetRangeUser(60,120);
	signalHist->GetXaxis()->SetRangeUser(60,120); 
      }
      */
      
      //      signalHist->SetFillColorAlpha(kRed+2, 3.0); // Colore blu, semi-trasparente
      //      signalHist->SetFillStyle(3004);// Linee diagonali
      signalHist->SetFillColor(0); // Colore blu, semi-trasparente
      signalHist->SetLineWidth(2);
      signalHist->SetLineColor(kRed);// Spessore del contorno aumentato
      //      signalHist->SetFillColor(kGreen);
      DYNonPromptHist->SetFillColor(kAzure+1);
      DYNonPromptHist->SetLineColor(kAzure+5);
      
      DYPromptHist->SetFillColor(kWhite);
      DYPromptHist->SetLineColor(kAzure+5);
      

      ZGPromptHist->SetFillColor(kOrange);
      ZGPromptHist->SetLineColor(kOrange+3);

      ZGNonPromptHist->SetFillColor(kOrange+5);
      ZGNonPromptHist->SetLineColor(kOrange+7);

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
      //      backgroundHist->SetLineColor(kAzure);
      // Verifica che gli istogrammi siano stati letti correttamente
      /*
      std::cout << "__" << std::endl;
      std::cout << "Histos exist" << std::endl;
      std::cout << "__" << std::endl;
      */
      // Crea il canvas per il plot
      TCanvas *canvas = new TCanvas("canvas", varNames.at(iVars)+"  after cut "+cutNb, 800, 1000);
    
      // Crea il pad principale per lo stack plot
      TPad *mainPad = new TPad("mainPad", "Main Pad", 0.0, 0.3, 1.0, 1.0);
      if(isLogScale) mainPad->SetLogy(); // Imposta l'asse Y in scala logaritmica
      mainPad->SetBottomMargin(0.1);
      mainPad->Draw();

      // Crea il pad inferiore per il rapporto segnale/fondo
      TPad *ratioPad = new TPad("ratioPad", "Ratio Pad", 0.0, 0.0, 1.0, 0.3);
      ratioPad->SetTopMargin(0.05);
      ratioPad->SetBottomMargin(0.3);
      ratioPad->Draw();

      // Disegna lo stack plot nel pad principale
      mainPad->cd();
      if(isLogScale)ZGPromptHist->SetMinimum(1.1);
	/*
      backgroundHist->SetMinimum(1.1); // Imposta il valore minimo dell'asse Y logaritmico
      backgroundHist->Draw("HIST");
      signalHist->Draw("HIST SAME");
      */

      TH1F *MCHist = (TH1F*)ZGPromptHist->Clone("MCHist");
      //      if(regions[iRegion]=="CR2P_1VL"){
      //      	MCHist->Add(DYPromptHist);
      	MCHist->Add(DYNonPromptHist);
	//      }
	//      MCHist->Add(signalHist);
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
      //stack->Add(signalHist);
      
      stack->Add(VZFSRHist);
      stack->Add(TTTo2L2NuHist);
      stack->Add(ZGPromptHist);
      stack->Add(DYNonPromptHist);
      //      stack->Add(signalHist);
      //      stack->Add(ZGNonPromptHist);
      //      stack->Add(DYPromptHist);

      stack->Draw("HIST"); // opzione "nostack" per sovrapporre gli istogrammi
      if(varName=="VZGMVAScore_")      stack->GetXaxis()->SetTitle("BDT Score");
      else      stack->GetXaxis()->SetTitle(varName);
      stack->GetYaxis()->SetTitle("Events");
      //      stack->SetMinimum(1.5); // Imposta il valore minimo dell'asse Y logaritmico
      //      if(cutNb>5)	  stack->SetMaximum(6000);
       
      
      signalHist->Draw("HIST SAME");
      signalHist->SetMinimum(2.1);

      MCHist->Draw("sameE2");

      double xLegend;
      if(varNames[iVars]=="VZGMVAScore_") xLegend = 0.58;
      else xLegend = 0.78;
      TLegend *legend = new TLegend(xLegend, 0.62, xLegend+0.19, 0.89);
      if(!isLogScale && sigZoom==50)legend->AddEntry(signalHist, "VZ#gamma x50", "f");
      else legend->AddEntry(signalHist, "VZ#gamma", "f");
      //      legend->AddEntry(DYPromptHist, "Drell-Yan + #gamma prompt)", "f");
      legend->AddEntry(DYNonPromptHist, "Drell-Yan + #gamma", "f");
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
      cmsText->SetTextFont(62); // 62 Font bold per "CMS"
      cmsText->SetTextAlign(13); // Allineamento al centro orizzontale (sinistra, alto)
      cmsText->DrawLatexNDC(0.1, 0.95, "CMS #bf{Private Work}");

      // Aggiunta della luminosità in alto a destra
      TLatex *luminosityText = new TLatex();
      luminosityText->SetTextSize(0.04);
      luminosityText->SetTextFont(42); // Font normale per luminosità
      luminosityText->SetTextAlign(32); // Allineamento al centro orizzontale (destra, alto)
      luminosityText->DrawLatexNDC(0.84, 0.93, "137.6 fb^{-1} (13 TeV)");


      // Calcola l'istogramma cumulato di segnale
      TH1F *signalCumulativeFWD = (TH1F*)signalHist->GetCumulative();
      TH1F *backgroundCumulativeFWD = (TH1F*)MCHist->GetCumulative();
      TH1F *signalCumulativeBWD = (TH1F*)signalHist->GetCumulative(false);
      TH1F *backgroundCumulativeBWD = (TH1F*)MCHist->GetCumulative(false);

      if(!isLogScale){
	signalCumulativeFWD->Scale(1./sigZoom);
	signalCumulativeBWD->Scale(1./sigZoom);
      }
      
      /*
      std::cout << "__" << std::endl;
      std::cout << "CumulativeHistos created" << std::endl;
      std::cout << "__" << std::endl;
      */
      for (Int_t i=0; i<=signalCumulativeBWD->GetNbinsX(); i++)
	{
	  Double_t binContent = signalCumulativeBWD->GetBinContent(i);
	  Double_t newBinContent;
	  
	  if(binContent<0.9)
	    signalCumulativeBWD->SetBinContent(i, 0.);
	}
      for (Int_t i=0; i<=signalCumulativeFWD->GetNbinsX(); i++)
	{
	  Double_t binContent = signalCumulativeFWD->GetBinContent(i);
	  Double_t newBinContent;
	  
	  if(binContent<0.9)
	    signalCumulativeFWD->SetBinContent(i, 0.);
	}




      // Calcola l'istogramma cumulato di segnale + fondo
      TH1F *total = (TH1F*)signalHist->Clone("total");
      if(!isLogScale)total->Scale(1./sigZoom);
      total->Add(MCHist);
      TH1F *totalCumulativeFWD = (TH1F*)signalCumulativeFWD->Clone("totalCumulativeFWD");
      totalCumulativeFWD->Add(backgroundCumulativeFWD);
      TH1F *totalCumulativeBWD = (TH1F*)signalCumulativeBWD->Clone("totalCumulativeBWD");
      totalCumulativeBWD->Add(backgroundCumulativeBWD);
      /*
      std::cout << "__" << std::endl;
      std::cout << "Total Histos created" << std::endl;
      std::cout << "__" << std::endl;
      */
      TH1F *sqrtTotal = (TH1F*)total->Clone("sqrtTotal");
      sqrtTotal->Reset();
      for (Int_t i=0; i<=total->GetNbinsX(); i++)
	{
	  Double_t binContent = total->GetBinContent(i);
	  Double_t sqrtContent;
	  
	  if(binContent<=0) sqrtContent =1.;
	  else  sqrtContent = TMath::Sqrt(binContent);
	  /*
	  std::cout << "Bin" << i  << "Den" << sqrtContent <<std::endl;
	  */
	  sqrtTotal->SetBinContent(i, sqrtContent);

	}
      /*
      std::cout << "__" << std::endl;
      std::cout << "sqrtTotal reset" << std::endl;
      std::cout << "__" << std::endl;
      */
      // Calcola la radice quadrata della somma degli istogrammi cumulati di segnale e fondo
      TH1F *sqrtTotalCumulativeFWD = (TH1F*)totalCumulativeFWD->Clone("sqrtTotalCumulativeFWD");
      sqrtTotalCumulativeFWD->Reset();
      for (Int_t i=0; i<=totalCumulativeFWD->GetNbinsX(); i++)
	{
	  Double_t binContent = totalCumulativeFWD->GetBinContent(i);
	  Double_t sqrtContent;
	  
	  if(binContent==0.) sqrtContent =0.001;
	  else  sqrtContent = TMath::Sqrt(binContent);
	  /*  
	  std::cout << "Bin" << i  << "Den" << sqrtContent <<std::endl;
	  */
	  sqrtTotalCumulativeFWD->SetBinContent(i, sqrtContent);
	  
	}
      TH1F *sqrtTotalCumulativeBWD = (TH1F*)totalCumulativeBWD->Clone("sqrtTotalCumulativeBWD");
      sqrtTotalCumulativeBWD->Reset();
      for (Int_t i=0; i<=totalCumulativeBWD->GetNbinsX(); i++)
	{
	  Double_t binContent = totalCumulativeBWD->GetBinContent(i);
	  Double_t sqrtContent;
	  
	  if(binContent<=0.)	    sqrtContent =1.;
	  else  sqrtContent = TMath::Sqrt(binContent);
	  
	  //	  std::cout << "Bin" << i  << "Den " << sqrtContent <<std::endl;
	  sqrtTotalCumulativeBWD->SetBinContent(i, sqrtContent);
	}
      /*

      std::cout << "__" << std::endl;
      std::cout << "sqrtTotal cumulative filled" << std::endl;
      std::cout << "__" << std::endl;
      */
      

      // Calcola il rapporto segnale/fondo bin per bin
      TH1F *ratioHist = (TH1F*)signalHist->Clone("ratioHist");
      ratioHist->Divide(sqrtTotal);
      /*for (Int_t i=0; i<=ratioHist->GetNbinsX(); i++)
	{
        std::cout << "Bin" << i  << "Num" << signalCumulativeFWD->GetBinContent(i) <<std::endl;
        std::cout << "Bin" << i  << "Sens " << ratioHistFWD->GetBinContent(i) <<std::endl;

	}*/
      // Calcola il rapporto segnale/fondo bin per bin
      TH1F *ratioHistFWD = (TH1F*)signalCumulativeFWD->Clone("ratioHist");
      ratioHistFWD->Divide(sqrtTotalCumulativeFWD);
      for (Int_t i=0; i<=ratioHistFWD->GetNbinsX(); i++)
	{
	  Double_t binContent = ratioHistFWD->GetBinContent(i);
	  
	  if(binContent>3.){
	    binContent =0.;
	    ratioHistFWD->SetBinContent(i, binContent);
	  }
	}

      /*      

      std::cout << "__" << std::endl;
      std::cout << "ratio issued" << std::endl;
      std::cout << "__" << std::endl;

      for (Int_t i=0; i<=ratioHistFWD->GetNbinsX(); i++)
	{
	  std::cout << "FWD Bin" << i  << "Num " << signalCumulativeFWD->GetBinContent(i) <<std::endl;
	  std::cout << "FWD Bin" << i  << "Sens " << ratioHistFWD->GetBinContent(i) <<std::endl;

	}
      */
      // Calcola il rapporto segnale/fondo bin per bin
      TH1F *ratioHistBWD = (TH1F*)signalCumulativeBWD->Clone("ratioHist");
      ratioHistBWD->Divide(sqrtTotalCumulativeBWD);
      for (Int_t i=0; i<=ratioHistBWD->GetNbinsX(); i++)
	{
	  Double_t binContent = ratioHistBWD->GetBinContent(i);
	  
	  if(binContent>3.){
	    binContent =0.;
	    ratioHistBWD->SetBinContent(i, binContent);
	  }
	}
      /*
      for (Int_t i=0; i<=ratioHistBWD->GetNbinsX(); i++)
	{
	  std::cout << "BWD Bin" << i  << "Num " << signalCumulativeBWD->GetBinContent(i) <<std::endl;
	  std::cout << "BWD Bin" << i  << "Sens " << ratioHistBWD->GetBinContent(i) <<std::endl;

	}
      */

      // Disegna il rapporto segnale/fondo nel pad inferiore
      ratioPad->cd();
      /*
      ratioHist->SetStats(0);
      ratioHist->SetTitle("");
      ratioHist->GetXaxis()->SetTitle(varName);
      ratioHist->GetYaxis()->SetTitle("sensitivity");
      ratioHist->GetYaxis()->SetRangeUser(minRangeSens,maxRangeSens);
      ratioHist->SetMarkerStyle(22);
      ratioHist->SetLineColor(kBlack);
      ratioHist->Draw("SAMES");
      */
      ratioHistFWD->SetStats(0);
      ratioHistFWD->SetTitle("");
      ratioHistFWD->GetYaxis()->SetRangeUser(minRangeSens,maxRangeSens);
      ratioHistFWD->SetMarkerStyle(22);
      ratioHistFWD->SetLineColor(kGreen);
      ratioHistFWD->SetFillColor(kWhite);
      ratioHistFWD->Draw("SAMES");

      ratioHistBWD->SetStats(0);
      ratioHistBWD->SetTitle("");
      ratioHistBWD->GetYaxis()->SetRangeUser(minRangeSens,maxRangeSens);
      ratioHistBWD->SetMarkerStyle(22);
      ratioHistBWD->SetLineColor(kViolet);
      ratioHistFWD->SetFillColor(kWhite);
      ratioHistBWD->Draw("SAMES");

      if (isLogScale){
	canvas->SaveAs("../DistributionPlots/Run2/LogScale/resultsForThesisPlots_"+varName+cutNb+"_.png");
	canvas->SaveAs("../DistributionPlots/Run2/LogScale/resultsForThesisPlots_"+varName+cutNb+"_.pdf");
      }

      // Salva il canvas in un file
      //canvas->Print("stack_plot.png");
      if(!isLogScale && sigZoom==50){
	canvas->SaveAs(path_to_rootOutput+"sensPlots/"+"50xSig/"+varName+cutNb+".png");
	canvas->SaveAs(path_to_rootOutput+"sensPlots/"+"50xSig/"+varName+cutNb+".pdf");
      }
      else if (!isLogScale) {
	canvas->SaveAs(path_to_rootOutput+"sensPlots/noZoomScale/"+varName+cutNb+".png");
	canvas->SaveAs(path_to_rootOutput+"sensPlots/noZoomScale/"+varName+cutNb+".pdf");
      }
    }
  }
  // Chiudi i file root
  VZGFile->Close();
  //  ZZGFile->Close();
  DYFile->Close();
  ZGFile->Close();
  TTTo2L2NuFile->Close();
  
}
