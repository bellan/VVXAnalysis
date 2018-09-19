/**
 *	Macro to study cuts from results of WWosAnalyzer
 *	
 *	It may become necessary to change the path and/or signal and background	
 *	
 *  $Date: 2018/09/13 13:37:23 $
 *  $Revision: 0.2 $
 *  \author A. Mecca alberto.mecca@edu.unito.it
 */

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <string>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

//#define TEST_MODE

using std::cout;
using std::vector;
using std::bitset;

TString nullTString = TString("");
TH1F* openGraph(TFile*, TString&, TString& = nullTString);	//points to the graph in the TFile
void printGraphs(TH1F*, TH1F*, const TString&);
TH1F* integralGraph(const TH1F* const origin);	//Creates a copy which is modified
TH1F* sqrtGraph(const TH1F* const origin);		//Creates a copy which is modified

void WWosCutAnalysis(){
	
	TString path = "~/Workspace/VVXAnalysis/TreeAnalysis/results/WWosAnalyzer_MC/";
	TString sigName = "WWEW";
	TString bkgName = "TTTo2L2Nu";
	vector<TString> eventTypes = {};//{"ee", "em", "mm"}; 
	vector<TString> cutConditions = {"lead pt >=", "second pt >="/*, "deltaR <="*/};
	
//graphs not categorized by event type
	vector<TString> extraGraphs = {"Lead jet E", "Lead jet pt", "delta R jets", "delta Eta jets", "mjj"};		
	vector<TString> leptonPlots = {"ptScalarLL", "ptVectLL", "mll", "angSeparationLL", "ptScalarLLmet", "ptVectLLmet"};
	
	TString graphNameToPrint;
	
	cout<<"Opening \""<<sigName<<".root\"\n";
	TFile* signalFile = TFile::Open(path + sigName + ".root");
	#ifdef TEST_MODE
	cout<<"Opened \""<<sigName<<".root\"\n";
	#endif
	
	cout<<"Opening \""<<bkgName<<".root\"\n";
	TFile* backgroundFile = TFile::Open(path + bkgName + ".root");
	#ifdef TEST_MODE
	cout<<"Opened \""<<bkgName<<".root\"\n";
	#endif
	
//################################	EVENT PLOTS ################################
	foreach(TString& type, eventTypes){
	cout<<type<<'\n';
	
		foreach(TString& condition, cutConditions){
			cout<<'\t'<<condition<<'\n';
			
			TString graphName = type+": events passing "+condition;
			
			TH1F* hSig = openGraph(signalFile, graphName, sigName);
			TH1F* hBkg = openGraph(backgroundFile, graphName, bkgName);
			
			graphNameToPrint = type+" "+condition;
			
			printGraphs(hSig, hBkg, graphNameToPrint);
		}
	}
	
	cout<<'\n';
	
//################################	JET PLOTS ################################
	foreach(TString& graphName, extraGraphs){
	
		TH1F* hSig = openGraph(signalFile, graphName, sigName);	
		TH1F* hBkg = openGraph(backgroundFile, graphName, bkgName);
		
		printGraphs(hSig, hBkg, graphName);
	}
	
//################################	LEPTON PLOTS ################################
	foreach(TString& graphName, leptonPlots){
	
		TH1F* hSig = openGraph(signalFile, graphName, sigName);
		TH1F* hBkg = openGraph(backgroundFile, graphName, bkgName);
		
		printGraphs(hSig, hBkg, graphName);
	}
}



// ################################################################################################
TH1F* openGraph(TFile* file, TString& graphName, TString& fileName){
	TH1F* hSig = (TH1F*)file->Get(graphName); 
	if(hSig == nullptr) 
		cout<<"Could not open \""<<graphName<<"\" from \""<<fileName<<".root"<<"\"\n";
	else{
		#ifdef TEST_MODE
		cout<<"Retrieved \""<<graphName<<"\" from \""<<fileName<<".root"<<"\"\n";
		#endif
	}
	return hSig;
}


void printGraphs(TH1F* hSig, TH1F* hBkg, const TString& graphName){
	int nbins = hSig->GetNbinsX();
	if(nbins != hBkg->GetNbinsX()) cout<<"Different number of bins\n";
	
	TString sigName = "Signal "+graphName;
	TString bkgName = "Background "+graphName;
	/*TAxis* thisXAxis = hSig->GetXaxis();
	Double_t xm = thisXAxis->GetBinLowEdge(1);
	Int_t xup = thisXAxis->GetNbins();
	Double_t xM = thisXAxis->GetBinLowEdge(xup);*/
	
	TH1F* hSigIntegral = integralGraph(hSig);
	TH1F* hBkbIntegral = integralGraph(hBkg);
	TH1F* hBkbIntegralSqrt = sqrtGraph(hBkbIntegral);
	
	//	Straigth comparison of signal and background
	TCanvas *c1 = new TCanvas(graphName+" comparison", graphName+" comparison", 10,0,1280,1024);
	c1->Divide(2,1);
	c1->cd(1);
	hSig->SetTitle(sigName);
	hSig->Draw();
	c1->cd(2);
	hBkg->SetTitle(bkgName);
	hBkg->Draw();
	
	//	Signal/background 
/*	TH1F* hRatio = (TH1F*)hSig->Clone(ratioGrName);
	hRatio->Divide(hBkg);
	TCanvas *c2 = new TCanvas(ratioGrName, ratioGrName, 10,0,1280,1024);
	c2->cd();
	hRatio->Draw();*/
	
	//	Integral plots comparison
	TCanvas *cIC = new TCanvas("Integral "+graphName+" comparison", "Integral "+graphName+" comparison", 10,0,1280,1024);
	cIC->Divide(2,1);
	cIC->cd(1);
	hSigIntegral->SetTitle(("Integral "+sigName).Data());
	hSigIntegral->Draw();
	cIC->cd(2);
	hBkbIntegral->SetTitle(("Integral "+bkgName).Data());
	hBkbIntegral->Draw();
	
	//	Integral plots SAME
	TCanvas *cIS = new TCanvas("Integral "+graphName, "Integral "+graphName, 10,0,1280,1024);
	cIS->cd();
	hSigIntegral->SetLineColor(2);
	hBkbIntegral->SetTitle(("Integral "+graphName).Data());
	hBkbIntegral->Draw();
	hSigIntegral->Draw("SAME");
	
	TLegend* lIS = new TLegend();
	lIS->AddEntry(hSigIntegral, "Signal Integral");
	lIS->AddEntry(hBkbIntegral, "Background Integral");
	lIS->Draw("SAME");
	
	//	Integral over sqrt(integral)
	TCanvas* cIoverSqrtOfI = new TCanvas("IoverSqrtOfI "+graphName, "IoverSqrtOfI "+graphName, 10,0,1280,1024);
	cIoverSqrtOfI->cd();
	TH1F* hIoverSqrtI = (TH1F*)hSigIntegral->Clone(graphName+"C");
	hIoverSqrtI->Divide(hBkbIntegralSqrt);
	hIoverSqrtI->SetTitle(("Integr(Sig)/sqrt(integr(bkg)) "+graphName).Data());
	hIoverSqrtI->SetLineColor(1);
	hIoverSqrtI->Draw();
	
	return;
}

TH1F* integralGraph(const TH1F* const origin){
	TString name = origin->GetName();
	TH1F* newGraph = (TH1F*)origin->Clone((name+"I").Data());
	int nbins = origin->GetNbinsX();
	for(int i = 1; i <= nbins; i++)
		newGraph->SetBinContent(i, origin->Integral(i,-1));
	
	return newGraph;
}

TH1F* sqrtGraph(const TH1F* const origin){
	TString name = origin->GetName();
	TH1F* newGraph = (TH1F*)origin->Clone((name+"Sqrt").Data());
	int nbins = origin->GetNbinsX();
	for(int i = 1; i <= nbins; i++){
		float s = sqrt(newGraph->GetBinContent(i));
		newGraph->SetBinContent(i, s);
	}
	return newGraph;
}







