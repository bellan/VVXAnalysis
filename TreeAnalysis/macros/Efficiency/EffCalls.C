#include "Efficiency.h"

#include <iostream>

char defaultPath[] = "../../results";

void effJets(const char* path, const char* name);  //recJets with respect to genJets
void effQuarks(const char* path, const char* name);  // genJets with respect to quarks
//TH1* unwrap2D(const TH2* th2);

void effQuarks(const char* name){
	effQuarks(defaultPath, name);
}
void effPhotons(const char* path, const char* name);

// +++++ main +++++
void effCalls(){
	char sampleName[] = "2018/VVGammaAnalyzer_MC/data.root";
	//effQuarks(defaultPath, sampleName);
	//effJets(defaultPath, sampleName);
	effPhotons(defaultPath, sampleName);
}


void effPhotons(const char* path, const char* name){
	TFile* fileIn = openTFile(Form("%s/%s", path, name));
	if(!fileIn) return;
	
	TFile fileOut("Eff_Photons.root", "RECREATE");
	if(fileOut.IsZombie()){
		fileOut.Close();
		fileIn->Close();
		std::cout<<"Could not recreate file \"Eff_Photons.root\"";
		return;
	}
	fileOut.cd();
	
	//doEfficiency2D(fileIn, "FAKE G: loose","FAKE G: medium","Fake G","Fake rate #gamma;p_{t} [GeV/c];#eta;Efficiency");
	doEfficiency2D(fileIn, "A:u ph loose TOT","A:u ph loose PASS","Fake G loose","#gamma_{LOOSE} / #gamma_{#geq4};p_{t} [GeV/c];|#eta|;Fake Rate");
	doEfficiency2D(fileIn, "A:u ph medium TOT","A:u ph medium PASS","Fake G medium","#gamma_{MEDIUM} / #gamma_{#geq4};p_{t} [GeV/c];|#eta|;Fake Rate");
	doEfficiency2D(fileIn, "A:u ph kinem","A:u ph loose PASS","Fake G loose kin","#gamma_{LOOSE} / #gamma_{kin};p_{t} [GeV/c];|#eta|;Fake Rate");
	
	fileOut.Close();
	fileIn->Close();
	delete fileIn;
}


void effJets(const char* path, const char* name){
	TFile* fileIn = openTFile(Form("%s/%s", path, name));
	if(!fileIn) return;
	
	TFile fileOut("Eff_Jets.root", "RECREATE");
	if(fileOut.IsZombie()){
		fileOut.Close();
		fileIn->Close();
		std::cout<<"Could not recreate file \"Eff_Quarks.root\"";
		return;
	}
	fileOut.cd();
	
	doEfficiency(fileIn, "effJ genAK4: eta", "effJ recAK4: eta", "Eff recAK4 vs eta", "Efficiency recAK4s vs #eta(j);#eta;Efficiency");
	doEfficiency(fileIn, "effJ genAK8: eta", "effJ recAK8: eta", "Eff recAK8 vs eta", "Efficiency recAK8s vs #eta(j);#eta;Efficiency");
	doEfficiency(fileIn, "effJ genAK4: pt",  "effJ recAK4: pt",  "Eff recAK4 vs pt", "Efficiency recAK4s vs pt(j);pt [GeV/c];Efficiency");
	doEfficiency(fileIn, "effJ genAK8: pt",  "effJ recAK8: pt",  "Eff recAK8 vs pt", "Efficiency recAK8s vs pt(j);pt [GeV/c];Efficiency");
	
	fileOut.Close();
	fileIn->Close();
	delete fileIn;
}


void effQuarks(const char* path, const char* name){
	TFile* fileIn = openTFile(Form("%s/%s", path, name));
	if(!fileIn) return;
	
	TFile fileOut("Eff_Quarks.root", "RECREATE");
	if(fileOut.IsZombie()){
		fileOut.Close();
		fileIn->Close();
		std::cout<<"Could not recreate file \"Eff_Quarks.root\"";
		return;
	}
	fileOut.cd();
	
	doEfficiency(fileIn, "deltaR quark", "deltaR quark: AK4s", "Eff genAK4s vs dRqq", "Efficiency genAK4s vs #DeltaR(q,q);#DeltaR(q,q)");
	doEfficiency(fileIn, "deltaR quark", "deltaR quark: AK8", "Eff genAK8 vs dRqq", "Efficiency genAK8 vs #DeltaR(q,q);#DeltaR(q,q)");
	//doEfficiency(fileIn, "Eta quark", "Eta quark: AK4", "Eff genAK4 vs eta q", "Efficiency genAK4 vs #eta(q);#eta(q)");
	doEfficiency(fileIn, "Eta quarks", "Eta quarks: AK8", "Eff genAK8 vs eta qq", "Efficiency genAK8 vs #eta(qq);#eta(qq)");
	doEfficiency(fileIn, "pt quarks", "pt quarks: AK4s", "Eff genAK4s vs pt qq", "Efficiency genAK4 vs pt(qq);pt(qq) [GeV/c]");
	doEfficiency(fileIn, "pt quarks", "pt quarks: AK8", "Eff genAK8 vs pt qq", "Efficiency genAK8 vs pt(qq);pt(qq) [GeV/c]", true);
	
	// Exist = knowing the p4 of the original quark(s)
	// Found = chosen by the bestV algorithm
	doEfficiency(fileIn, "pt quarks: AK4s", "pt quarks: AK4s found good", "Eff (strict) choice genAK4s vs pt qq", "Efficiency (strict) choice of genAK4s vs pt(qq);pt(qq) [GeV/c]");
	doEfficiency(fileIn, "pt quarks: AK4s algo", "pt quarks: AK4s found good", "Purity (strict) choice genAK4s vs pt qq", "Purity choice of genAK4s vs pt(qq);pt(qq) [GeV/c]");
	doEfficiency(fileIn, "pt quarks: AK4s", "pt quarks: AK4s found", "Eff (loose) choice genAK4s vs pt qq", "Efficiency (loose) choice of genAK4s vs pt(qq);pt(qq) [GeV/c]");
	doEfficiency(fileIn, "pt quarks: AK4s algo", "pt quarks: AK4s found", "Purity (loose) choice genAK4s vs pt qq", "Purity choice of genAK4s vs pt(qq);pt(qq) [GeV/c]");
	
	doEfficiency(fileIn, "pt quarks: AK8", "pt quarks: AK8 found good", "Eff (strict) choice genAK8 vs pt qq", "Efficiency (strict) choice of genAK8 vs pt(qq);pt(qq) [GeV/c]");
	doEfficiency(fileIn, "pt quarks: AK8 algo", "pt quarks: AK8 found good", "Purity (strict) choice genAK8 vs pt qq", "Purity choice of genAK8 vs pt(qq);pt(qq) [GeV/c]", true);
	doEfficiency(fileIn, "pt quarks: AK8", "pt quarks: AK8 found", "Eff (loose) choice genAK8 vs pt qq", "Efficiency (loose) choice of genAK8 vs pt(qq);pt(qq) [GeV/c]");
	doEfficiency(fileIn, "pt quarks: AK8 algo", "pt quarks: AK8 found", "Purity (loose) choice genAK8 vs pt qq", "Purity choice of genAK8 vs pt(qq);pt(qq) [GeV/c]", true);
	
	fileOut.Close();
	fileIn->Close();
	delete fileIn;
}
/*
TH1* unwrap2DI(const TH2* th2){
	TAxis xaxis = hTot->GetXaxis();
	int nbinsx = xaxis->GetNbins();
	std::vector<double> xbins; xbins.reserve(nbinsx+1);
	for(Int_t bx = 1; bx <= nbinsx; bx++){
		xbins->push_back(xaxis->GetBinLowEdge(bx));
	}
	xbins->push_back(xaxis->GetBinUpEdge(nbinsx));
	
	TAxis yaxis = hTot->GetXaxis();
	int nbinsy = yaxis->GetNbins();
	std::vector<double> ybins; ybins.reserve(nbinsy+1);
	for(Int_t by = 1; by <= nbinsy; by++){
		ybins->push_back(yaxis->GetBinLowEdge(by));
	}
	ybins->push_back(yaxis->GetBinUpEdge(nbinsy));
	
	TH1* unwrapped = new TH1(Form("%s_unwrapped", th2->GetName()), nbinsx*nbinsy,1,nbinsx*nbinsy+1);
	
	for(Int_t by = 1; by <= nbinsy; by++){
		for(Int_t bx = 1; bx <= nbinsx; bx++){
			th2->GetBinContent(bx, by);
			th2->GetBinError(bx, by);
		}
	}
}*/


/*
TH1* unwrap2DF(const TH2* th2){
	TAxis xaxis = hTot->GetXaxis();
	int nbinsx = xaxis->GetNbins();
	std::vector<double> xbins; xbins.reserve(nbinsx+1);
	for(Int_t bx = 1; bx <= nbinsx; bx++){
		xbins->push_back(xaxis->GetBinLowEdge(bx));
	}
	xbins->push_back(xaxis->GetBinUpEdge(nbinsx));
	
	TAxis yaxis = hTot->GetXaxis();
	int nbinsy = yaxis->GetNbins();
	std::vector<double> ybins; ybins.reserve(nbinsy+1);
	for(Int_t by = 1; by <= nbinsy; by++){
		ybins->push_back(yaxis->GetBinLowEdge(by));
	}
	ybins->push_back(yaxis->GetBinUpEdge(nbinsy));
	
	for(Int_t by = 1; by <= nbinsy; by++){
		TH1F slice(Form("%s_slicey%d", th2->GetName(), by), xbins.data());
		for(Int_t bx = 1; bx <= nbinsx; bx++){
			slice.SetBinContent(bx, th2->GetBinContent(bx, by);
			slice.SetBinError(bx, th2->GetBinError(bx, by);
		}
	}
}*/
