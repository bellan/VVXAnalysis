#include "Efficiency.h"

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TEfficiency.h"

TFile* openTFile(const char* path){
	#ifdef TEST_MODE
	std::cout<<"\nOpening \""<<path<<"\" \t";
	#endif
	TFile* file = nullptr;
	file = TFile::Open(path, "READ");
	if(file == nullptr){
		std::cout<<"Error: file is nullptr\n";
		return nullptr;
	}
	if(file->IsZombie()){
		std::cout<<"\""<<path<<"\" does not exist\n";
		file->Close();
		delete file;
		return nullptr;
	}
	return file;
}

void doEfficiency(TFile* fin, const char* nTot, const char* nPas, const char* nName, const char* nTitle, bool extraFix){
	std::cout<<" ----- "<<nName<<" ----- 	\n";
	TH1* hTot = (TH1*)fin->Get(nTot);
	if(!hTot){
		std::cout<<Form("Could not find \"%s\".\n", nTot);
		return;
	}
	TH1* hPas = (TH1*)fin->Get(nPas);
	if(!hPas){
		std::cout<<Form("Could not find \"%s\".\n", nPas);
		return;
	}
	
	if(extraFix) doExtraFix(hTot, hPas);
	
	auto type = (hTot->InheritsFrom("TH1I") ? "cp" : "normal");
	TGraphAsymmErrors* eff = new TGraphAsymmErrors(hPas, hTot, type);
	eff->SetName(nName);
	eff->SetTitle(nTitle);
	eff->SetMinimum(0.);
	eff->SetMaximum(1.01);
	eff->Write();
}

void doExtraFix(TH1* total, TH1* pass){
	for(Int_t i = 0; i < total->GetNbinsX(); ++i) {
		if(pass->GetBinContent(i) > total->GetBinContent(i)) {
			float tmp = pass->GetBinContent(i);
			pass->SetBinContent(i, total->GetBinContent(i));
			std::cout<<"(Fixed bin "<<i<<":  old = "<<tmp<<"  new = "<<total->GetBinContent(i)<<")\n";
		}
	}
}


void doEfficiency(TFile* fin, const char* nTot, const char* nPas, const char* nName, bool extraFix){
	TH1* hTot = (TH1*)fin->Get(nTot);
	if(!hTot){
		std::cout<<Form("Could not find \"%s\".\n", nTot);
		return;
	}
	TString labels = getAxesLabels(hTot);
	const char* nTitle = Form("%s;%s", nName, labels.Data()); // No problem if they're empty
	doEfficiency(fin, nTot, nPas, nName, nTitle, extraFix);
}


TString getAxesLabels(const TH1* hist){
	if(!hist) return TString();
	const char* xTitle = hist->GetXaxis()->GetTitle();
	const char* yTitle = hist->GetYaxis()->GetTitle();
	return TString(Form("%s;%s", xTitle, yTitle));  // No problem if they're empty
}


void verboseEff(const TH1* hTot, const TH1* hPas, double cl){
	if(cl < 0. || 1.<cl){
		cerr << "Error: confidence level must be 0 < cl < 1. Given cl = "<<cl<<'\n';
		return;
	}
	for(int bin = 0; bin <= hTot->GetNbinsX()+1; bin++){
		double tot = hTot->GetBinContent(bin);
		double pas = hPas->GetBinContent(bin);
		cout<<'\t'<<" ["<<hTot->GetBinLowEdge(bin)<<", "<<hTot->GetBinLowEdge(bin+1) <<"]: ";
		if(tot > 0)
			cout << '[' << TEfficiency::ClopperPearson(tot,pas,cl,false) << ", " << pas/tot << ", " << TEfficiency::ClopperPearson(tot,pas,cl,true) << ']';
		else cout<<"[0, 0, 1]";
		cout<<'\n';
	}
}


void doEfficiency2D(TFile* fin, const char* nTot, const char* nPas, const char* nName, const char* nTitle){
	std::cout<<" ----- "<<nName<<" ----- 	\n";
	TH2* hTot = (TH2*)fin->Get(nTot);
	if(!hTot){
		std::cout<<Form("Could not find \"%s\".\n", nTot);
		return;
	}
	TH2* hPas = (TH2*)fin->Get(nPas);
	if(!hPas){
		std::cout<<Form("Could not find \"%s\".\n", nPas);
		return;
	}
	
	// 2D efficiency map
	TH2* eff = (TH2*) hPas->Clone(nName);
	eff->Divide(hTot);
	//TGraphAsymmErrors* eff = new TGraphAsymmErrors(hPas, hTot, type);
	
	eff->SetName(nName);
	eff->SetTitle(nTitle);
	eff->SetMinimum(0.);
	eff->SetMaximum(1.01);
	eff->Write();
	
	// Projections
	TH1D* hTot_X = hTot->ProjectionX();
	TH1D* hPas_X = hPas->ProjectionX();
	TH1D* hTot_Y = hTot->ProjectionY();
	TH1D* hPas_Y = hPas->ProjectionY();
	
	TGraphAsymmErrors* eff_X = new TGraphAsymmErrors(hPas_X, hTot_X, "cp");
	TGraphAsymmErrors* eff_Y = new TGraphAsymmErrors(hPas_Y, hTot_Y, "cp");
	eff_X->SetName(Form("%s_projX", eff->GetName()));
	eff_Y->SetName(Form("%s_projY", eff->GetName()));
	auto mainTitle = eff->GetTitle();
	auto label_Z = eff->GetZaxis()->GetTitle();
	eff_X->SetTitle(Form("%s;%s;%s", mainTitle, eff->GetXaxis()->GetTitle(), label_Z));
	eff_Y->SetTitle(Form("%s;%s;%s", mainTitle, eff->GetYaxis()->GetTitle(), label_Z));
	
	eff_X->Write();
	eff_Y->Write();
	
	//cout<<"Proj X"<<'\n';
	//verboseEff(hTot_X, hPas_X);
	//cout<<"Proj Y"<<'\n';
	//verboseEff(hTot_Y, hPas_Y);
}
