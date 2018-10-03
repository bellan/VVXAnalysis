/**
 *	Macro to study cuts from results of WWosAnalyzer
 *	
 *	It may become necessary to change the path and/or signal and background	
 *	
 *	Usage: 	root [-l] 'WWosCutAnalysis.C("<signal sample name>", "<background sample name>", "<category of graphs>", "<requested graphs by name>", "<type of graphs>")'
 *	e.g.: 	root [-l] 'WWosCutAnalysis.C("WWEW", "TTTo2L2Nu", "jetlep", "all", "156")'
 *
 *	It's possible to use it with standard values, by simply calling:	root -l WWosCutAnalysis
 *	
 *	<category of graphs> = {"cut", "jet", "lepton"}
 *	<requested graphs by name> depends on the category
 *	<type of graphs> = {
 *		1 sig and bkg comparison
 *		2 sig to bkg ratio
 *		3 Integral graphs comparison
 *		4 Integral graps on same canvas
 *		5 sig/sqrt(bkg) (integrals)
 *		6 sig/bkg (integrals)
 *	}
 *	
 *  $Date: 2018/09/13 13:37:23 $
 *  $Revision: 0.3 $
 *  \author A. Mecca alberto.mecca@edu.unito.it
 */

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <THStack.h>
#include <string>
#include <algorithm>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

//#define TEST_MODE

using std::cout;
using std::vector;
using std::bitset;
using std::string;

TString nullTString = TString("");

struct GNames{
public:
	GNames(const TString& a,const TString& b,const vector<TString>& c) : graphName(a), sigGrName(b), bkgGrNames(c){};
	TString graphName;
	TString sigGrName;
	vector <TString> bkgGrNames;
};

template <class TH=TH1F>
struct MyGraphs{
	MyGraphs(TH* a, vector<TH*> b): hSig(a), vhBkgs(b) {};
	TH* hSig;
	vector<TH*> vhBkgs;
};

TFile* openRootFile(const TString& /*path*/, const TString& /*name*/);	//prints to cout the name of the file
template <class TH1Z = TH1F>
TH1Z* openGraph(TFile*, const TString&, const TString& = nullTString);//points to the graph in the TFile
template<size_t T>
void printGraphs(TH1F*, TH1F*, const TString&, const bitset<T>& grToDo, const TString& = nullTString, const TString& = nullTString);
TH1F* integralGraph(const TH1F* const origin); //Creates a copy which is modified
TH1F* reverseIntegrGraph(const TH1F* const origin);	//Integrates from the last bin
TH1F* sqrtGraph(const TH1F* const origin);		 //Creates a copy which is modified
void compareCutGraphs(const TH1I* const hSig, const TH1I* const hBkg, const TString& sigName,const TString& bkgName);

template <class TH>
void printGraphCompare(MyGraphs<TH>& graphs, const GNames& names);


void WWosCutAnalysis(string requestSig = string(""), string requestBkg = string(""), string requestType = string(""), string reqGrName = string(""), string sGrToDo = string("")){
	
	TString path = "~/Workspace/VVXAnalysis/TreeAnalysis/results/WWosAnalyzer_MC/";
	
	std::transform(reqGrName.begin(), reqGrName.end(), reqGrName.begin(), ::tolower);
	std::transform(requestType.begin(), requestType.end(), requestType.begin(), ::tolower);
	TString sigName = (requestSig == string("") ?   "WWEW"    : requestSig);
	TString bkgName = (requestBkg == string("") ? "TTTo2L2Nu" : requestBkg);
	
	vector<TString> eventTypes = {};//{"ee", "em", "mm"}; 
	vector<TString> cutConditions = {"lead pt >=", "second pt >="/*, "deltaR <="*/};
	
//graphs not categorized by event type
	vector<TString> jetPlots = {"Lead jet E", "Lead jet pt", "Trail jet E", "Trail jet pt", "delta R jets", "delta Eta jets", "mjj", "lead jet csvtagger", "trail jet csvtagger"};		
	vector<TString> leptonPlots = {"ptScalarLL", "ptVectLL", "mll", "angSeparationLL", "ptScalarLLmet", "ptVectLLmet"};
	
	
	TFile* signalFile = openRootFile(path, sigName);
	TFile* backgroundFile = openRootFile(path, bkgName);
	
	bitset<3> toDo(string("000"));
	
	if(requestType.find("all") != string::npos) toDo.set();  //All bits to 1
	if(requestType.find("cut") != string::npos) toDo.set(0); //set to 1
	if(requestType.find("jet") != string::npos) toDo.set(1);	
	if(requestType.find("lep") != string::npos) toDo.set(2);
	if(requestType == string(""))								toDo.set(2);
	
	bitset<6> grToDo;	
	if(sGrToDo == string("")){ 
		grToDo.set(0); 
		grToDo.set(4); 
	}
	else if(sGrToDo.find("all") != string::npos)
		grToDo.set();
	else{
		if(sGrToDo.find('1') != string::npos) grToDo.set(0); //set to 1
		if(sGrToDo.find('2') != string::npos) grToDo.set(1);	
		if(sGrToDo.find('3') != string::npos) grToDo.set(2);
		if(sGrToDo.find('4') != string::npos) grToDo.set(3);
		if(sGrToDo.find('5') != string::npos) grToDo.set(4);
		if(sGrToDo.find('6') != string::npos) grToDo.set(5);
	}
	
	cout<<"toDo: "<<toDo<<" \tgrToDo: "<<grToDo<<'\n';
	
//################################	EVENT PLOTS ################################
	if(toDo.test(0)){
		/*foreach(TString& type, eventTypes){
			#ifdef TEST_MODE
			cout<<type<<'\n';
			#endif
			foreach(TString& condition, cutConditions){
				#ifdef TEST_MODE
				cout<<'\t'<<condition<<'\n';
				#endif
				TString graphName = type+": events passing "+condition;
				
				TH1F* hSig = openGraph(signalFile, graphName, sigName);
				TH1F* hBkg = openGraph(backgroundFile, graphName, bkgName);
				
				printGraphs(hSig, hBkg, type+" "+condition, grToDo, sigName, bkgName);
			}
		}*/
		
		TH1I* hSig = openGraph<TH1I>(signalFile, "Cuts", sigName);
		TH1I* hBkg = openGraph<TH1I>(backgroundFile, "Cuts", bkgName);
		
		compareCutGraphs(hSig, hBkg, sigName, bkgName);
	}
	
//################################	JET PLOTS ################################
	if(toDo.test(1)){
		foreach(TString& graphName, jetPlots){
			if(reqGrName != string(""))
				if(reqGrName.find(graphName) == string::npos && reqGrName.find("all") == string::npos) continue;
			TH1F* hSig = openGraph(signalFile, graphName, sigName);	
			TH1F* hBkg = openGraph(backgroundFile, graphName, bkgName);
		
			printGraphs(hSig, hBkg, graphName, grToDo, sigName, bkgName);
		}
	}
	
//################################	LEPTON PLOTS ################################
	if(toDo.test(2)){
		foreach(TString& graphName, leptonPlots){
			TString copy = graphName;			copy.ToLower();
			if(reqGrName != string(""))
				if(reqGrName.find(graphName) == string::npos && reqGrName.find("all") == string::npos) continue;
			TH1F* hSig = openGraph(signalFile, graphName, sigName);
			TH1F* hBkg = openGraph(backgroundFile, graphName, bkgName);
				
			printGraphs(hSig, hBkg, graphName, grToDo, sigName, bkgName);
		}
	}
	cout<<'\n';
}



// ################################################################################################
TFile* openRootFile(const TString& path, const TString& name){
	cout<<"Opening \""<<name<<".root\"\n";
	TFile* file = TFile::Open(path + name + ".root");
	if(file == nullptr) cout<<"\""<<name<<".root\" does not exist\n";
	return file;
}

template <class TH1Z = TH1F>
TH1Z* openGraph(TFile* file, const TString& graphName, const TString& fileName){
	if(file == nullptr) return nullptr;
	TH1Z* hSig = (TH1Z*)file->Get(graphName); 
	if(hSig == nullptr) 
		cout<<"Could not open \""<<graphName<<"\" from \""<<fileName<<".root"<<"\"\n";
	else{
		#ifdef TEST_MODE
		cout<<"Retrieved \""<<graphName<<"\" from \""<<fileName<<".root"<<"\"\n";
		#endif
	}
	return hSig;
}

template<size_t T>
void printGraphs(TH1F* hSig, TH1F* hBkg, const TString& graphName, const bitset<T>& grToDo, const TString& sigName,const TString& bkgName){
	if(hSig == nullptr || hBkg == nullptr) return;
	int nbins = hSig->GetNbinsX();
	if(nbins != hBkg->GetNbinsX()) cout<<"Different number of bins\n";
	
	TString sigGrName = ((sigName == nullTString ?  "Signal"  : sigName) + " " + graphName);
	TString bkgGrName = ((bkgName == nullTString ?"Background": bkgName) + " " + graphName);
	/*TAxis* thisXAxis = hSig->GetXaxis();
	Double_t xm = thisXAxis->GetBinLowEdge(1);
	Int_t xup = thisXAxis->GetNbins();
	Double_t xM = thisXAxis->GetBinLowEdge(xup);*/
	
	//	#####	Graph creation #####
	TH1F* hSigIntegral;
	TH1F* hBkbIntegral;
	TH1F* hBkbIntegralSqrt;
	if(graphName.Contains("csvtagger", TString::kIgnoreCase)){
		hSigIntegral = reverseIntegrGraph(hSig);
		hBkbIntegral = reverseIntegrGraph(hBkg);
	} else{
		hSigIntegral = integralGraph(hSig);
		hBkbIntegral = integralGraph(hBkg);
	}
	hBkbIntegralSqrt = sqrtGraph(hBkbIntegral);
	
	
	//	-----	Straigth comparison of signal and background
	if(grToDo.test(0)){
		TCanvas *c1 = new TCanvas(graphName+" comparison", graphName+" comparison", 10,0,1280,1024);
		c1->Divide(2,1);
		c1->cd(1);
		hSig->SetTitle(sigGrName);
		hSig->Draw();
		c1->cd(2);
		hBkg->SetTitle(bkgGrName);
		hBkg->Draw();
	}
	
	//	-----	Signal/background ratio
	if(grToDo.test(1)){
		TH1F* hRatio = (TH1F*)hSig->Clone((graphName+" R").Data());
		hRatio->Divide(hBkg);
		TCanvas *c2 = new TCanvas(graphName+" ratio", graphName+" ratio", 10,0,1280,1024);
		c2->cd();
		hRatio->Draw();
	}
	
	//	-----	Integral plots comparison
	hSigIntegral->SetTitle((sigGrName+" Integral").Data());
	hSigIntegral->SetLineColor(2);
	hBkbIntegral->SetTitle((bkgGrName+" Integral").Data());
	if(grToDo.test(2)){
		TCanvas *cIC = new TCanvas(graphName+" Integral comparison", graphName+" Integral comparison", 10,0,1280,1024);
		cIC->Divide(2,1);
		cIC->cd(1);
		hSigIntegral->Draw();
		cIC->cd(2);
		hBkbIntegral->Draw();
	}
	
	//	-----	Integral plots SAME
	if(grToDo.test(3)){
		TCanvas *cIS = new TCanvas(graphName+" Integral", graphName+" Integral", 10,0,1280,1024);
		cIS->cd();
		hBkbIntegral->Draw();
		hSigIntegral->Draw("SAME");
	
		TLegend* lIS = new TLegend(0.2, 0.10);
		lIS->AddEntry(hSigIntegral, sigName+" Integral");
		lIS->AddEntry(hBkbIntegral, bkgName+" Integral");
		lIS->Draw("SAME");
	}
	
//	###########
	/*TStyle *MyStyle = new TStyle("MyStyle","My Style");
	MyStyle->SetTitleFontSize(0.04);
	gROOT->SetStyle("MyStyle");*/
//  ###########
	
	//	-----	Integral over sqrt(integral)
	if(grToDo.test(4)){
		//################ temp stuff
		/*
		TCanvas* cISandSqrtOfIB = new TCanvas(graphName+" I(S)andSqrt(I(B))", graphName+" I(S)andSqrt(I(B))", 10,0,1280,1024);
		cISandSqrtOfIB->cd();
		hSigIntegral->Draw();
		hBkbIntegralSqrt->Draw("SAME");
		
		TLegend* lIS = new TLegend(0.2, 0.10);
		lIS->AddEntry(hSigIntegral, "I("+sigName+")");
		lIS->AddEntry(hBkbIntegralSqrt, "#sqrt{I("+bkgName+")}");
		lIS->Draw("SAME");
		*/
		//################ end temp stuff
		TCanvas* cIoverSqrtOfI = new TCanvas(graphName+" I(S)/Sqrt(I(B))", graphName+" I(S)/Sqrt(I(B))", 10,0,1280,1024);
		cIoverSqrtOfI->cd();
		TH1F* hIoverSqrtI = (TH1F*)hSigIntegral->Clone(graphName+" S/#sqrt{B}");
		hIoverSqrtI->Divide(hBkbIntegralSqrt);
		TString laSIoSBI("I(S)/#sqrt{I(Bkg)} - "+graphName);
		hIoverSqrtI->SetTitle(laSIoSBI);//(("Integr(sig)/sqrt(integr(bkg)) "+graphName).Data());
		hIoverSqrtI->SetLineColor(2);
		hIoverSqrtI->SetMinimum(0.);
		hIoverSqrtI->Draw();
	}
	
	//	-----	signal kept over background kept
	if(grToDo.test(5)){
		TCanvas* cSIoBI = new TCanvas(graphName+" I(S)/I(B)", graphName+" I(S)/I(B)", 10,0,1280,1024);
		cSIoBI->cd(); 
		TH1F* hISoIB = (TH1F*)hSigIntegral->Clone(graphName+" I(S)/I(B)");
		hISoIB->Divide(hBkbIntegral);
		TString laSIoBI("I(S)/I(B) - "+graphName);
		hISoIB->SetTitle(laSIoBI);//(("Integr(sig)/Integr(bkg) "+graphName).Data());
		hISoIB->SetLineColor(2);
		hISoIB->SetMinimum(0.);
		hISoIB->Draw();
	}
	
	return;
}

void printCompareCanv(const TH1F* const hSig, const TH1F* const hBkg, const TString& graphName){
	
}

TH1F* integralGraph(const TH1F* const origin){
	TString name = origin->GetName();
	TH1F* newGraph = (TH1F*)origin->Clone((name+" I").Data());
	int nbins = origin->GetNbinsX();
	for(int i = 1; i <= nbins; i++)
		newGraph->SetBinContent(i, origin->Integral(i,-1));
	
	return newGraph;
}

TH1F* reverseIntegrGraph(const TH1F* const origin){
	TString name = origin->GetName();
	TH1F* newGraph = (TH1F*)origin->Clone((name+" rI").Data());
	int nbins = origin->GetNbinsX();
	for(int i = nbins; i >= 1; i--)
		newGraph->SetBinContent(i, origin->Integral(1, i));
	
	return newGraph;
}

TH1F* sqrtGraph(const TH1F* const origin){
	TString name = origin->GetName();
	TH1F* newGraph = (TH1F*)origin->Clone((name+" Sqrt").Data());
	int nbins = origin->GetNbinsX();
	for(int i = 1; i <= nbins; i++){
		float s = sqrt(newGraph->GetBinContent(i));
		newGraph->SetBinContent(i, s);
	}
	return newGraph;
}

void compareCutGraphs(const TH1I* const hSig, const TH1I* const hBkg, const TString& sigName,const TString& bkgName){
	unsigned int nbins = hSig->GetNbinsX();
	TH1F* hSigf = new TH1F("Cuts"+sigName,"Cuts"+sigName, nbins,0,nbins);
	TH1F* hBkgf = new TH1F("Cuts"+bkgName,"Cuts"+bkgName, nbins,0,nbins);
	const TAxis* const aSigAxis = hSig->GetXaxis();
	TAxis* aSigfAxis = hSigf->GetXaxis();
	TAxis* aBkgfAxis = hSigf->GetXaxis();
	for(int i = 1; i<= nbins; i++){
		hSigf->SetBinContent(i, (float)hSig->GetBinContent(i));
		aSigfAxis->SetBinLabel(i, aSigAxis->GetBinLabel(i)); 
		hBkgf->SetBinContent(i, (float)hBkg->GetBinContent(i));
		aBkgfAxis->SetBinLabel(i, aSigAxis->GetBinLabel(i)); 
	}
	TCanvas* cCut = new TCanvas("Cuts", "Cuts", 10,0,1280,1024);
	cCut->cd();
	hSigf->Scale(1./hSig->GetBinContent(1));
	hSigf->SetLineColor(2);
	hBkgf->Scale(1./hBkg->GetBinContent(1));
	hSigf->Draw();
	hBkgf->Draw("SAME");
	TLegend* lCut = new TLegend(0.2, 0.10);
	lCut->AddEntry(hSigf, sigName);
	lCut->AddEntry(hBkgf, bkgName);
	lCut->Draw("SAME");
}

template <class TH>
void printGraphCompare(MyGraphs<TH>& graphs, const GNames& names){
	TCanvas *c1 = new TCanvas(names.graphName+" comp", names.graphName+" comp", 10,0,1280,1024);
	THStack* stack1 = new THStack(names.graphName+" comp",names.graphName+" comp");
	
	TH* hSigC = (TH*)graphs.hSig->Clone(names.sigGrName);
	hSigC->SetTitle(names.sigGrName);
	stack1->Add(hSigC);
	
	vector<TH*> vhBkgCs;
	for(int i = 0; i < graphs.vhBkgs.size(); i++){
		TH* cg = (TH*)(graphs.vhBkgs.at(i)->Clone(names.bkgGrNames.at(i)));
		cg->SetTitle(names.bkgGrNames.at(i));
		vhBkgCs.push_back(cg);
		stack1->Add(cg);
	}
	
	stack1->Draw();
}









