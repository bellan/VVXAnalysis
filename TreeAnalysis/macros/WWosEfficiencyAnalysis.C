/**
 *	Macro to get efficiency from results of WWosAnalyzer and ZVAnalyzer
 *	
 *	Usage: 	root [-l] [-b] [-q] 'WWosEfficiencyAnalysis.C("<sample name>")'		-l do not display banner 	-b run in background		-q close after finishing
 *	e.g.: 	root -l 'WWosEfficiencyAnalysis.C("WWEW")'
 *	It may become necessary to change the path and/or to expand the samples list	
 *	
 *  $Date: 2018/09/12 11:37:23 $
 *  $Revision: 0.2 $
 *  \author A. Mecca alberto.mecca@edu.unito.it
 */

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define TEST_MODE

using namespace std;

void WWosEfficiencyAnalysis(TString requestedSample){
	
	TString path = "~/Workspace/VVXAnalysis/TreeAnalysis/results/WWosAnalyzer_MC/";
	vector<TString> samples = {"WWEW", "WWQCD", "WWEWQCD"};
	vector<TString> parNames = {"Electrons", "Muons", "Jets"};
	vector<TString> typeNames = {"eta", "phi", "pt"};


	TString sampleName;
	if(requestedSample != ""){
		for(int i = 0; i < samples.size(); i++){
			if(requestedSample == samples.at(i))
				sampleName = requestedSample;
		}
	}
	
	if(sampleName == ""){
		cout<<"Unknown sample \""<<sampleName<<"\"\n";
		return;
	}
	cout<<"Opening \""<<sampleName<<".root\"\n";
	TFile* result = TFile::Open(path + sampleName + ".root");
	
	foreach(TString& name, parNames){
		#ifdef TEST_MODE
		cout<<name<<"\n";
		#endif
		foreach(TString& type, typeNames){
			#ifdef TEST_MODE
			cout<<"\t"<<type<<"\n";
			#endif
			TH1F* hMatch = (TH1F*)result->Get(name+"Matched_vs_"+type);
			if(hMatch == nullptr){
				#ifdef TEST_MODE
				cout<<"Could not open \""<<name<<"Matched_vs_"<<type<<"\"\n";
				#endif
			}
			
			TH1F* hTot = (TH1F*)result->Get("gen"+name+"_"+type);
			if(hTot == nullptr){
				#
				cout<<"Could not open \""<<"gen"<<name<<"_"<<type<<"\"\n";
			}
			
			if(hMatch != nullptr && hTot != nullptr){
				TGraphAsymmErrors* a = new TGraphAsymmErrors(hMatch, hTot, "cp");
				a->SetTitle(name+"Efficiency_vs_"+type);
				a->GetYaxis()->SetRangeUser(0.,1.01);
				TCanvas *cDrawing = new TCanvas(name+"Efficiency_vs_"+type, name+"Efficiency_vs_"+type, 10,0,1280,1024);
				cDrawing->cd();
				a->Draw("AP");
			}
		}
	}
	result->Close("R");
}







