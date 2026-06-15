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

bool makeSRplots = true;
double minRangeSens = 0.;
double maxRangeSens = 2.0;
bool isLogScale = true;
int aggrBin = 2;

TString resDir= "resultsZGvectSum";//"resultsForNNLOshapeDisagreementFullRun2";//"resultsForNNLOshapeDisagreement";//
//TString year="2018";
const std::vector<TString> years = {"2017"};//"Run2"};//

const std::vector<TString> varNames = {"POST-RWGT_ptGamma_"};//"PRE-RWGT_ZGammaPtSum+_", "POST-RWGT_ZGammaPtSum_"};// ,"PRE-RWGT_J0bVetoVar_", "PRE-RWGT_J1bVetoVar_"};//, "PRE-RWGT_ptGamma_","PRE-RWGT_recoZPt_","BIN_PtG20-40_PRE-RWGT_recoZPt_", "BIN_PtG40-60_PRE-RWGT_recoZPt_", "BIN_PtG60-80_PRE-RWGT_recoZPt_", "BIN_PtG80-100_PRE-RWGT_recoZPt_", "BIN_PtG100-200_PRE-RWGT_recoZPt_",  };// {"DR_gammaClosestLept_", "recoVMass_", "j1_p(g)", "j0_p(uds)", "mllG_", "PhotonMVAID_", "recoVDaughter0Pt_", "recoVDaughter1Pt_","ptGamma_", "dPhiZG_"};// {"dPhiZG_","VZGMVAScore_","PhotonMVAID_", "recoZDaughter0Pt_","recoZDaughter1Pt_","DR_gammaClosestLept_","mllG_","recoVMass_", "recoVDaughter0Pt_", "recoVDaughter1Pt_", "recoZMass_", "ptGamma_","DR_gammaClosestJet_", "j0_p(g)","j1_p(g)", "j0_p(uds)","j1_p(uds)"};// {/*"recoZPt_", "ptGamma_", "recoVDaughter0Pt_", "recoVDaughter1Pt_", "recoZDaughter0Pt_", "recoZDaughter1Pt_",*/"recoZPt_","dPhiZG_","VZGMVAScore_","j0_p(bX)","j1_p(bX)","MET_","PhotonMVAID_", "ZepCorr_","recoZDaughter0Pt_","recoZDaughter1Pt_","DR_gammaClosestLept_","mllG_","mjjG_", "recoVMass_", "recoVDaughter0Pt_", "recoVDaughter1Pt_", "recoZMass_", "ptGamma_","DR_gammaClosestJet_", "j0_p(g)","j1_p(g)", "j0_p(uds)","j1_p(uds)","FWM_T0_fullSyst_", "FWM_T1_fullSyst_","jj_p(uds)Sum","jj_p(g)Sum"};//{"recoVMass","recoVDaughter0Pt"};// {"### HAD TOPO DJ cand mass"};//,"recoVDaughter0Pt","DR_gammaClosestLept", "DR_Lept","FWM_T0_fullSyst", "FWM_T1_fullSyst","FWM_T2_fullSyst", "FWM_T3_fullSyst", "FWM_T4_fullSyst","FWM_T5_fullSyst", "FWM_T6_jets", "FWM_T1_jets","FWM_T2_jets", "FWM_T3_jets","FWM_T4_jets", "FWM_T5_jets", "FWM_T6_jets", "mllG",  "recoVMass", "DR_gammaClosestJet", "recoVDaughter1Pt", "recoVDaughtersDeltaPhi", "recoVPt", "recoZDeltaPhi", "recoZEta", "recoZMass", "recoZPt"};//, "System_Pt"};//, "relativePT_G_vs_V","VBH0s","VBH0z","VBH0t","jjH0s","jjH0z","jjH0t"};
const std::vector<TString> regions = {"CR2P_1VL", "CRZON_FSRT"};//"CR2P_1VL_bContam", "CR2P_1VL_DJbContam","CRZON_FSRT_DJbContam" };// "CRZON_FSRT_bVetoed", "CRZON_FSRT_DJbVetoed"};//{"CR2P_1VL", "CRZON_FSRT_bVetoed", "CR2P_1VL_bVetoed", "CRZON_FSRT_DJbVetoed", "CR2P_1VL_DJbVetoed"};//, "CRZON_FSRT_bContam", "CR2P_1VL_bContam", "CRZON_FSRT_DJbContam", "CR2P_1VL_DJbContam"};// {"CR2P_1VL"};//{"CR2P_1VL_bVetoed"};// {"CRZOFF"};// {"CR2P_1VL","CRZON_FSRT","CRZOFF", "CRZON_FSRT", "CRZOFF_FSRT", "CRZOFF_DIB", "CRFSRT" };//


void CMSDataMC_testBkg()
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
    TFile *elseFile = new TFile(path_to_rootOutput+"else.root", "READ"); //WZTo3LNu.root
    /*
    TFile *WZFile = new TFile(path_to_rootOutput+"WZTo3LNu.root", "READ"); //WZTo3LNu.root
    TFile *ZZFile = new TFile(path_to_rootOutput+"ZZTo2Q2L.root", "READ"); //WZTo3LNu.root
    */
    /*  TFile *ttXFile = new TFile(path_to_rootOutput+"ttX.root", "READ");
	TFile *TZqFile = new TFile(path_to_rootOutput+"TZq.root", "READ");
	TFile *WXFile = new TFile(path_to_rootOutput+"WX.root", "READ");
	TFile *qqZZFile = new TFile(path_to_rootOutput+"qqZZ.root", "READ");
    */
    TFile *TTTo2L2NuFile = new TFile(path_to_rootOutput+"TTTo2L2Nu.root", "READ"); //TTTo2L2Nu
    TFile *DataFile = new TFile(path_to_rootOutput+"data_obs.root", "READ");

    for(int iRegion=0; iRegion<regions.size(); iRegion++){
      for (int iVars=0; iVars<varNames.size(); ++iVars){
	TString varName=varNames[iVars];
	TString varLabel=varNames[iVars];
	if(varName=="PRE-RWGT_ZGammaPtSum+_") varLabel="pre-rew., |#bf{p_{T}^{ll}}+#bf{p_{T}^{#gamma}}| [GeV]";
	if(varName=="POST-RWGT_ZGammaPtSum_") varLabel="post-rew., |#bf{p_{T}^{ll}}+#bf{p_{T}^{#gamma}}| [GeV]";
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
	if(varName=="ptGamma_") varLabel="p_{T}^{#gamma} [GeV]";
	if(varName=="PhotonMVAID_") varLabel="#gamma MVA ID";
	if(varName=="j1_p(g)") varLabel="sub-leading jet DeepFlavor(g) ";
	if(varName=="j0_p(uds)") varLabel="leading jet DeepFlavor(uds) ";
	if(varName=="PRE-RWGT_recoZPt_") varLabel="Z Cand. p_{T} [GeV]";
	if(varName=="POST-RWGT_recoZPt_") varLabel="Z Cand. p_{T} [GeV]";
	if(varName=="PRE-RWGT_ptGamma_") varLabel="#gamma Cand. p_{T} [GeV]";
	if(varName=="BIN_PtG20-40_PRE-RWGT_recoZPt_") varLabel="Z Cand. p_{T} [GeV]";
	if(varName=="BIN_PtG40-60_PRE-RWGT_recoZPt_") varLabel="Z Cand. p_{T} [GeV]";
	if(varName=="BIN_PtG60-80_PRE-RWGT_recoZPt_") varLabel="Z Cand. p_{T} [GeV]";
	if(varName=="BIN_PtG80-100_PRE-RWGT_recoZPt_") varLabel="Z Cand. p_{T} [GeV]";
	if(varName=="BIN_PtG100-200_PRE-RWGT_recoZPt_") varLabel="Z Cand. p_{T} [GeV]";
	// Leggi gli istogrammi di segnale e background dai file root

	cout<<varName<<"nonPrompt_"<<regions[iRegion]<<"2"<<endl;


	TH1F *signalHist = (TH1F*)VZGFile->Get(varName+"sign_"+regions[iRegion]+"2");
   
	//      TH1F *DYPromptHist = (TH1F*)DYFile->Get(varName+"prompt_"+regions[iRegion]+"2");
	cout<<varName<<"nonPrompt_"<<regions[iRegion]<<"2"<<endl;
	TH1F *DYNonPromptHist = (TH1F*)DYFile->Get(varName+"nonPrompt_"+regions[iRegion]+"2");
	TH1F *ZGPromptHist = (TH1F*)ZGFile->Get(varName+"prompt_"+regions[iRegion]+"2");
	//      TH1F *ZGNonPromptHist = (TH1F*)ZGFile->Get(varName+"nonPrompt_"+regions[iRegion]+"2");
	TH1F *VZFSRHist = (TH1F*)VZGFile->Get(varName+"bckg_"+regions[iRegion]+"2");
	TH1F *TTTo2L2NuHist = (TH1F*)TTTo2L2NuFile->Get(varName+"all_"+regions[iRegion]+"2");
	/*
	TH1F *ZZHist = (TH1F*)ZZFile->Get(varName+"all_"+regions[iRegion]+"2");
	TH1F *WZHist = (TH1F*)WZFile->Get(varName+"all_"+regions[iRegion]+"2");
	*/
	TH1F *elseHist = (TH1F*)elseFile->Get(varName+"all_"+regions[iRegion]+"2");
	TH1F *dataHist = (TH1F*)DataFile->Get(varName+"all_"+regions[iRegion]+"2");
	/*
	TH1F *signalHist = (TH1F*)VZGFile->Get(varName+"sign"+regions[iRegion]+"2");
   
	//      TH1F *DYPromptHist = (TH1F*)DYFile->Get(varName+"prompt"+regions[iRegion]+"2");
	cout<<varName<<"nonPrompt"<<regions[iRegion]<<"2"<<endl;
	TH1F *DYNonPromptHist = (TH1F*)DYFile->Get(varName+"nonPrompt"+regions[iRegion]+"2");
	TH1F *ZGPromptHist = (TH1F*)ZGFile->Get(varName+"prompt"+regions[iRegion]+"2");
	//      TH1F *ZGNonPromptHist = (TH1F*)ZGFile->Get(varName+"nonPrompt"+regions[iRegion]+"2");
	TH1F *VZFSRHist = (TH1F*)VZGFile->Get(varName+"bckg"+regions[iRegion]+"2");
	TH1F *TTTo2L2NuHist = (TH1F*)TTTo2L2NuFile->Get(varName+"all"+regions[iRegion]+"2");
	TH1F *dataHist = (TH1F*)DataFile->Get(varName+"all"+regions[iRegion]+"2");
	*/
	//      DYPromptHist->Rebin(aggrBin);
	DYNonPromptHist->Rebin(aggrBin);
	ZGPromptHist->Rebin(aggrBin);
	//      ZGNonPromptHist->Rebin(aggrBin);
	VZFSRHist->Rebin(aggrBin);
	signalHist->Rebin(aggrBin);
	dataHist->Rebin(aggrBin);
	TTTo2L2NuHist->Rebin(aggrBin);
	/*
	WZHist->Rebin(aggrBin);
        ZZHist->Rebin(aggrBin);
	*/
        elseHist->Rebin(aggrBin);
	/*            
		      if(varName=="mllG_" || varName=="recoVMass_" || varName=="VZGMVAScore_"){
		      backgroundHist->Rebin(2);
		      ZGHist->Rebin(2);
		      VZFSRHist->Rebin(2);
		      signalHist->Rebin(2);
		      dataHist->Rebin(2);
		      }
		      if(varName=="jj_p(g)Sum" || varName=="j0_p(g)" || varName=="j1_p(g)"){
		      backgroundHist->Rebin(2);
		      ZGHist->Rebin(2);
		      VZFSRHist->Rebin(2);
		      signalHist->Rebin(2);
		      dataHist->Rebin(2);
		      }
	*/

	dataHist->SetMarkerStyle(20);
	dataHist->SetMarkerSize(.9);
	dataHist->SetMarkerColor(kBlack);
	dataHist->SetLineColor(kBlack);
      
	signalHist->SetFillColorAlpha(kRed+2, 3.0); // Colore blu, semi-trasparente
	signalHist->SetFillStyle(3004);// Linee diagonali
	signalHist->SetLineWidth(2);
	signalHist->SetLineColor(kRed);// Spessore del contorno aumentato
	/*
	  signalHist->SetFillColor(kRed);
	  signalHist->SetLineColor(kRed+1);
	*/
	DYNonPromptHist->SetFillColor(kAzure+1);
	DYNonPromptHist->SetLineColor(kAzure+5);
      
	//      DYPromptHist->SetFillColor(kWhite);
	//      DYPromptHist->SetLineColor(kAzure+5);
      
	ZGPromptHist->SetFillColor(kOrange);
	ZGPromptHist->SetLineColor(kOrange+3);

	//      ZGNonPromptHist->SetFillColor(kOrange+5);
	//      ZGNonPromptHist->SetLineColor(kOrange+7);
	/*
	  ZGNonPromptHist->SetFillColorAlpha(kOrange, 3.0); // Colore blu, semi-trasparente
	  ZGNonPromptHist->SetFillStyle(3004);// Linee diagonali
	  ZGNonPromptHist->SetLineWidth(2);
	  ZGNonPromptHist->SetLineColor(kOrange+3);// Spessore del contorno aumentato
	*/
	VZFSRHist->SetFillColor(kSpring+5);
	VZFSRHist->SetLineColor(kSpring+3);
      
	TTTo2L2NuHist->SetFillColor(kBlue); // Colore blu, semi-trasparente
	//      TTTo2L2NuHist->SetFillStyle(3004);// Linee diagonali
	//      TTTo2L2NuHist->SetLineWidth(1);
	TTTo2L2NuHist->SetLineColor(kBlue);// Spessore del contorno aumentato

	elseHist->SetFillColor(kMagenta+2);// Spessore del contorno aumentato
	elseHist->SetLineColor(kMagenta+2);// Spessore del contorno aumentato
	/*
	WZHist->SetFillColor(kMagenta+2);// Spessore del contorno aumentato
	WZHist->SetLineColor(kMagenta+2);// Spessore del contorno aumentato

	ZZHist->SetFillColor(kMagenta);// Spessore del contorno aumentato
	ZZHist->SetLineColor(kMagenta);// Spessore del contorno aumentato
	*/

	/*
	  if (!signalHist) {
	  std::cout << "Errore nella lettura dell'istogramma di segnale" << std::endl;
	  return;
	  }
	  if (!DYHist) {
	  std::cout << "Errore nella lettura dell'istogramma di background" << std::endl;
	  return;
	  }
	*//*
	    std::cout << "__" << std::endl;
	    std::cout << "Histos exist" << std::endl;
	    std::cout << "__" << std::endl;
	  */
	// Crea il canvas per il plot
	TCanvas *canvas = new TCanvas("canvas", varNames.at(iVars)+" Region "+regions[iRegion], 800, 1000);
    
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
	/*
	  backgroundHist->SetMinimum(1.1); // Imposta il valore minimo dell'asse Y logaritmico
	  backgroundHist->Draw("HIST");
	  signalHist->Draw("HIST SAME");
	*/
	TH1F *MCHist = (TH1F*)ZGPromptHist->Clone("MCHist");
	//      if(regions[iRegion]=="CR2P_1VL"){
	//	MCHist->Add(DYPromptHist);
	MCHist->Add(DYNonPromptHist);
	//      }
	MCHist->Add(signalHist);
	//MCHist->Add(ZGPromptHist);
	//MCHist->Add(ZGNonPromptHist);
      
	MCHist->Add(VZFSRHist);
	MCHist->Add(TTTo2L2NuHist);
	//	MCHist->Add(ttXHist);
	//	MCHist->Add(tZqHist);
	//	MCHist->Add(WZHist);
	//	MCHist->Add(ZZHist);
	MCHist->Add(elseHist);

	
	MCHist->SetFillStyle(3005);
	MCHist->SetMarkerStyle(1);
	MCHist->SetFillColor(kBlack);
	MCHist->SetLineColor(kBlack);// Spessore del contorno aumentato
	MCHist->SetLineWidth(1.0);

	
      
	THStack *stack = new THStack("stack", "");
	//      stack->Add(backgroundHist);
	stack->Add(signalHist);

	stack->Add(VZFSRHist);
	/*
	stack->Add(WZHist);
	stack->Add(ZZHist);
	*/
	stack->Add(elseHist);
	stack->Add(TTTo2L2NuHist);

	
	//	stack->Add(ttXHist);
	//	stack->Add(tZqHist);
      
	stack->Add(ZGPromptHist);
	stack->Add(DYNonPromptHist);


	//stack->Add(ZGNonPromptHist);
	//      stack->Add(DYPromptHist);
	//stack->Add(ZGHist);
	//	stack->Add(ZGHist_FSR);


	stack->Draw("HIST"); // opzione "nostack" per sovrapporre gli istogrammi

	//stack->GetXaxis()->SetTitle(varName);
	stack->GetYaxis()->SetTitle("Events");
	stack->GetXaxis()->SetTitle(varLabel);

	//      stack->SetMinimum(1.5); // Imposta il valore minimo dell'asse Y logaritmico
	//if(cutNb>10)	  stack->SetMaximum(4000);
      
	if( regions[iRegion]=="CR2P_1F" || regions[iRegion]=="CR2P_1VL"){
	  if (varName == "recoZMass_" || varName == "recoVDaughter0Pt_" || varName == "recoVDaughter0Pt_") stack->SetMaximum(aggrBin*50000);
	  else stack->SetMaximum(aggrBin*50000);
	}
      
	if(isLogScale) stack->SetMinimum(1.1);

	/*if(varName=="recoZMass_" || varName=="recoVDaughter1Pt_" || varName=="VZGMVAScore_")	  stack->SetMaximum(8000);
	  else stack->SetMaximum(3000);
	*/
	//      signalHist->Draw("HIST SAME");
	//signalHist->SetMinimum(2.1);
	dataHist->Draw("E SAME"); 
	MCHist->Draw("sameE2");
	
	TLegend *legend = new TLegend(0.86, 0.65, 0.99, 0.89);
	legend->AddEntry(signalHist, "VZ#gamma", "f");
	//legend->AddEntry(DYPromptHist, "Drell-Yan + prompt #gamma", "f");
	if(varName=="PRE-RWGT_recoZPt_"){
	  legend->AddEntry(DYNonPromptHist, "Drell-Yan (pre-rew.)", "f");
	  legend->AddEntry(ZGPromptHist, "Z#gamma (pre-rew.)", "f");	
	}else if(varName=="POST-RWGT_recoZPt_"){
	  legend->AddEntry(DYNonPromptHist, "Drell-Yan (post-rew.)", "f");
	  legend->AddEntry(ZGPromptHist, "Z#gamma (post-rew.)", "f");	
	}else{
	  legend->AddEntry(DYNonPromptHist, "Drell-Yan", "f");
	  legend->AddEntry(ZGPromptHist, "Z#gamma", "f");	
	}
	//      legend->AddEntry(ZGNonPromptHist, "Z#gamma non-prompt", "f");	
	//	legend->AddEntry(ZGHist_FSR, "Z#gamma(FSR)", "f");
	//	legend->AddEntry(ttXHist, "qqZZ", "f");
	legend->AddEntry(VZFSRHist, "VZ+FSR", "f");
	legend->AddEntry(TTTo2L2NuHist, "tt #rightarrow 2l2#nu", "f");
	//legend->AddEntry(WZHist, "WZ", "f");
	legend->AddEntry(elseHist, "other", "f");
	legend->AddEntry(MCHist, "Pred. unc.", "f");
	legend->AddEntry(dataHist, "Data", "lep");
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
	if(year=="2016preVFP")  luminosityText->DrawLatexNDC(0.9, 0.93, "19.5 fb^{-1} (13 TeV)");
	if(year=="2016postVFP")  luminosityText->DrawLatexNDC(0.9, 0.93, "16.8 fb^{-1} (13 TeV)");
	if(year=="2017")  luminosityText->DrawLatexNDC(0.9, 0.93, "41.5 fb^{-1} (13 TeV)");
	if(year=="2018")  luminosityText->DrawLatexNDC(0.9, 0.93, "59.8 fb^{-1} (13 TeV)");
	if(year=="Run2")  luminosityText->DrawLatexNDC(0.9, 0.93, "137.6 fb^{-1} (13 TeV)");
	//luminosityText->DrawLatexNDC(0.9, 0.93, year+" (13 TeV)");



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

	if (isLogScale){
	  canvas->SaveAs("../Data_vs_MC/"+year+"/LogScale/"+resDir+"_"+varName+"_"+regions[iRegion]+".png");
	  canvas->SaveAs("../Data_vs_MC/"+year+"/LogScale/"+resDir+"_"+varName+"_"+regions[iRegion]+".pdf");
	}
	else {
	  canvas->SaveAs("../Data_vs_MC/"+year+"/LogScale/"+resDir+"_"+varName+"_"+regions[iRegion]+".png");
	  canvas->SaveAs("../Data_vs_MC/"+year+"/LogScale/"+resDir+"_"+varName+"_"+regions[iRegion]+".pdf");
	}
      }
    }

    //Now for SR
    if(makeSRplots){
      for (int cutNb = 2; cutNb < 3; cutNb++){
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

	  //      TH1F *DYPromptHist = (TH1F*)DYFile->Get(varName+"prompt"+cutNb);
	  TH1F *DYNonPromptHist = (TH1F*)DYFile->Get(varName+"nonPrompt"+cutNb);
	  TH1F *ZGPromptHist = (TH1F*)ZGFile->Get(varName+"prompt"+cutNb);
	  //      TH1F *ZGNonPromptHist = (TH1F*)ZGFile->Get(varName+"nonPrompt"+cutNb);

	  TH1F *VZFSRHist = (TH1F*)VZGFile->Get(varName+"bckg"+cutNb);
	  TH1F *TTTo2L2NuHist = (TH1F*)TTTo2L2NuFile->Get(varName+"all"+cutNb);
	  //      TH1F *dataHist = (TH1F*)DataFile->Get(varName+"all"+cutNb); //SR UNBLINDED
	  TH1F *elseHist = (TH1F*)elseFile->Get(varName+"all"+cutNb);


	  //      DYPromptHist->Rebin(aggrBin);
	  ZGPromptHist->Rebin(aggrBin);

	  DYNonPromptHist->Rebin(aggrBin);
	  //      ZGNonPromptHist->Rebin(aggrBin);

	  TTTo2L2NuHist->Rebin(aggrBin);
	  VZFSRHist->Rebin(aggrBin);
	  signalHist->Rebin(aggrBin);
	  elseHist->Rebin(aggrBin);
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
	  elseHist->SetFillColor(kMagenta+2);
	  elseHist->SetLineColor(kMagenta+2);
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
	  MCHist->Add(elseHist);

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
	  stack->Add(elseHist);
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
	  legend->AddEntry(elseHist, "other", "f");
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
	  if(year=="2016preVFP")  luminosityText->DrawLatexNDC(0.9, 0.93, "19.5 fb^{-1} (13 TeV)");
	  if(year=="2016postVFP")  luminosityText->DrawLatexNDC(0.9, 0.93, "16.8 fb^{-1} (13 TeV)");
	  if(year=="2017")  luminosityText->DrawLatexNDC(0.9, 0.93, "41.5 fb^{-1} (13 TeV)");
	  if(year=="2018")  luminosityText->DrawLatexNDC(0.9, 0.93, "59.8 fb^{-1} (13 TeV)");
	  if(year=="Run2")  luminosityText->DrawLatexNDC(0.9, 0.93, "137.6 fb^{-1} (13 TeV)");

	  //luminosityText->DrawLatexNDC(0.9, 0.93, year+" (13 TeV)");


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

	  if (isLogScale){
	    canvas->SaveAs("../Data_vs_MC/"+year+"/LogScale/"+resDir+"_"+varName+"_SR"+cutNb+".png");
	    canvas->SaveAs("../Data_vs_MC/"+year+"/LogScale/"+resDir+"_"+varName+"_SR"+cutNb+".pdf");
	  }
	  else {
	    canvas->SaveAs("../Data_vs_MC/"+year+"/LogScale/"+resDir+"_"+varName+"_SR"+cutNb+".png");
	    canvas->SaveAs("../Data_vs_MC/"+year+"/LogScale/"+resDir+"_"+varName+"_SR"+cutNb+".pdf");
	  }

	}
      }
    }  
    // Chiudi i file root
    VZGFile->Close();
    elseFile->Close();
    TTTo2L2NuFile->Close();
    DYFile->Close();
    ZGFile->Close();
    DataFile->Close();
    //    ttXFile->Close();
    //    tZqFile->Close();
    //ZGEnrichFile->Close();
  }
  
}
