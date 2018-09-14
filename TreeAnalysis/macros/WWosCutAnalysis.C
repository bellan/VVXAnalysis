/**
 *	Macro to study cuts from results of WWosAnalyzer
 *	
 *	It may become necessary to change the path and/or signal and background	
 *	
 *  $Date: 2018/09/13 13:37:23 $
 *  $Revision: 0.1 $
 *  \author A. Mecca alberto.mecca@edu.unito.it
 */

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define TEST_MODE

using std::cout;
using std::vector;

void WWosCutAnalysis(){
	
	TString path = "~/Workspace/VVXAnalysis/TreeAnalysis/results/WWosAnalyzer_MC/";
	TString sigName = "WWEW";
	TString bkgName = "WWQCD";
	vector<TString> eventTypes = {"ee", "em", "mm"};
	vector<TString> cutConditions = {"lead pt >=", "second pt >=", "deltaR <="};
	
	
	cout<<"Opening \""<<sigName<<".root\"\n";
	TFile* signalFile = TFile::Open(path + sigName + ".root");
	cout<<"Opened \""<<sigName<<".root\"\n";
	
	cout<<"Opening \""<<bkgName<<".root\"\n";
	TFile* backgroundFile = TFile::Open(path + bkgName + ".root");
	cout<<"Opened \""<<bkgName<<".root\"\n";
	
	foreach(TString& type, eventTypes){
	cout<<type<<'\n';
		foreach(TString& condition, cutConditions){
			cout<<'\t'<<condition<<'\n';
			
			TString graphName = type+": events passing "+condition;
			TH1F* hSig = (TH1F*)signalFile->Get(graphName);
			/*cout<<"****"<<graphName<<endl;
			TCanvas *cDrawing = new TCanvas(type+" : "+condition, type+" : "+condition, 10,0,1280,1024);
			cDrawing->cd();
			hSig->Draw("hist");*/
			if(hSig == nullptr) 
				cout<<"Could not open \""<<graphName<<"\" from \""<<sigName<<".root"<<"\"\n";
			else
				cout<<"Retrieved \""<<graphName<<"\" from \""<<sigName<<".root"<<"\"\n";
			
			
			/*TCanvas *c = new TCanvas(type+" "+condition, type+" "+condition, 10,0,1280,1024);
			c->cd();
			hSig->Draw();*/
			
			TH1F* hBkg = (TH1F*)backgroundFile->Get(graphName);
			if(hBkg == nullptr) 
				cout<<"Could not open \""<<graphName<<"\" from \""<<bkgName<<".root"<<"\"\n";
			else
				cout<<"Retrieved \""<<graphName<<"\" from \""<<bkgName<<".root"<<"\"\n";
				
			int nbins = hSig->GetNbinsX();
			if(nbins != hBkg->GetNbinsX()) cout<<"Different number of bins\n";
			
			
			/*TAxis* thisXAxis = hSig->GetXaxis();
			Double_t xm = thisXAxis->GetBinLowEdge(1);
			Int_t xup = thisXAxis->GetNbins();
			Double_t xM = thisXAxis->GetBinLowEdge(xup);*/
			
			//hSig->Divide(hBkg);
			TH1F* hRatio = (TH1F*)hSig->Clone(type+" ratio "+condition);
			hRatio->Sumw2();
			//TH1F* hNewBkg = (TH1F*)hBkg->Clone(type+" bkg "+condition);
			//hNewBkg->Sumw2();
			
			TCanvas *c1 = new TCanvas(type+" :1 "+condition, type+" :1 "+condition, 10,0,1280,1024);
			c1->Divide(2,1);
			c1->cd(1);
			hSig->Draw();
			c1->cd(2);
			hBkg->Draw();
			
			
			hRatio->Divide(hBkg);
			
			TCanvas *cDrawing = new TCanvas(type+" : "+condition, type+" : "+condition, 10,0,1280,1024);
			cDrawing->cd();
			hRatio->Draw();
			
			/*TH1F* newGraph = new TH1F(type+" : "+condition, type+" : "+condition, nbins, xm, xM);
			for(float t = 0.; t <= xM; t++){
				newGraph->Fill(t, )
			}*/
		}
	}
}







