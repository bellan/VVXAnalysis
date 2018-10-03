/**
 *	Macro to study cuts from results of WWosAnalyzer
 *	
 *	It may become necessary to change the path and/or signal and background	
 *	
 *	Usage: 	root [-l] 'WWosCutAnalysis.C("<category of graphs>", "<requested graphs by name>", "<type of graphs>")'
 *	e.g.: 	root [-l] 'WWosCutAnalysis.C("jetlep", "all", "13")'
 *
 *	It's possible to use it with standard values, by simply calling:	root -l WWosCutAnalysis
 *	
 *	<category of graphs> = {"cut", "jet", "lep", "jP", "lP"} //P = pre-second cut
 *	<requested graphs by name> depends on the category
 *	<type of graphs> = {
 *		0 sig and bkg same canvas (no stack)
 *		1 sig and bkg same canvas
 *		2 Integral graps on same canvas
 *		3 sig/sqrt(bkg) (integrals)
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
	TString toString() const{
		TString res = "graphName: "+graphName+"\nsigGrName: "+sigGrName+"\nbkgGrNames: ";
		foreach(const TString& b, bkgGrNames) res += b+" \t";
		return res;
	}
};

template <class TH=TH1F>
struct MyGraphs{
	MyGraphs(TH* a, vector<TH*> b): hSig(a), vhBkgs(b) {};
	MyGraphs(TH* a, vector<TH*> b, TH* c, vector<TH*> d): hSig(a), vhBkgs(b), hSigInt(c), vhBkgInt(d) {};
	TH* hSig;
	vector<TH*> vhBkgs;
	TH* hSigInt;
	vector<TH*> vhBkgInt;
	TString toString() const{
		return TString("hSig? ")+(hSig != nullptr)+" - vhBkgs.size(): "+vhBkgs.size();
	}
};

TFile* openRootFile(const TString&/*path*/,const TString&/*name*/);//prints to cout the name of the file
template <class TH1Z = TH1F>
TH1Z* openGraph(TFile*, const TString&, const TString& = nullTString);//points to the graph in the TFile
TH1F* integralGraph(const TH1F* const origin); //Creates a copy which is modified
TH1F* reverseIntegrGraph(const TH1F* const origin);	//Integrates from the last bin
TH1F* sqrtGraph(const TH1F* const origin);		 //Creates a copy which is modified

template <size_t N, class TH=TH1F>
void newPrintGraph(const MyGraphs<TH>& theGraphs, const GNames& theNames, const bitset<N>& grToDo);

// Specific types of graph
template <class TH>
void compareCutGraphs(const MyGraphs<TH>&, const GNames&);
template <class TH>
void printGraphSame(const MyGraphs<TH>&, const GNames&);
template <class TH>
void printGraphCompare(const MyGraphs<TH>&, const GNames&);
template <class TH>
void printGraphIntegr(const MyGraphs<TH>&, const GNames&);
template <class TH>
void printGraphSqrt(const MyGraphs<TH>&, const GNames&);

//Struct building
template <class TH = TH1F>
MyGraphs<TH>* buildMyGraphs(const GNames& names, TFile* const signalFile, vector<TFile*>& backgroundFiles);

//Colors
vector<Color_t> myColors ={kRed, kOrange-3, kGreen+1, kAzure+10, kBlue, kViolet+1, kMagenta, kOrange-2};

void WWosStackGraphs(string requestType = string(""), string reqGrName = string(""), string sGrToDo = string("")){
	
	TString path = "~/Workspace/VVXAnalysis/TreeAnalysis/results/WWosAnalyzer_MC/";
	//TString path = "~/Workspace/Private/Results/";
	cout<<"Source path: \""<<path<<"\"\n";
	vector<TString> backgrounds = {"WWQCD", "TTTo2L2Nu", "DYJetsToLL_M50", /*"TTJets"*/};
	
	std::transform(reqGrName.begin(), reqGrName.end(), reqGrName.begin(), ::tolower);
	std::transform(requestType.begin(), requestType.end(), requestType.begin(), ::tolower);
	TString sigName = "WWEW";
	
//graphs not categorized by event type
	vector<TString> jetPlots = {"Lead jet E", "Lead jet pt", "Trail jet E", "Trail jet pt", "delta R jets", "delta Eta jets", "mjj", "lead jet csvtagger", "trail jet csvtagger"};
	vector<TString> jetAll(jetPlots.size());
	for(size_t k = 0; k < jetPlots.size(); k++)
		jetAll.at(k) = (jetPlots.at(k)+" all");
		
	vector<TString> leptonPlots = {"ptScalarLL", "ptVectLL", "mll", "angSeparationLL", "ptScalarLLmet", "ptVectLLmet"};
	vector<TString> lepAll(leptonPlots.size());
	for(size_t k = 0; k < leptonPlots.size(); k++)
		lepAll.at(k) = (leptonPlots.at(k)+" all");
	
	
	TFile* signalFile = openRootFile(path, sigName);

	vector<TFile*> backgroundFiles;
	vector<size_t> notFound;
	for(size_t i = 0; i < backgrounds.size(); i++){
		//cout<<"backgrounds.at("<<i<<") = "<<backgrounds.at(i)<<'\n';
		backgroundFiles.push_back( openRootFile(path, backgrounds.at(i)) );
		cout<<" \tnullptr? "<<(backgroundFiles.at(i) == nullptr);
		if(backgroundFiles.at(i) == nullptr)
			notFound.push_back(i);
	}
	foreach(size_t i, notFound){
		backgrounds.erase(backgrounds.begin() + i);
		backgroundFiles.erase(backgroundFiles.begin() + i);
	}
	
	cout<<"\nbackgrounds:";
	foreach(TString& b, backgrounds) cout<<"\n\t"<<b;
	cout<<"\nbackgroundFiles:";
	foreach(TFile*& f, backgroundFiles) cout<<"\n\tfile found";
	cout<<'\n';
	
	bitset<5> toDo(string("00000"));
	if(requestType.find("all") != string::npos) toDo.set();  //All bits to 1
	if(requestType.find("cut") != string::npos) toDo.set(0); //set to 1
	if(requestType.find("jet") != string::npos) toDo.set(1);	
	if(requestType.find("lep") != string::npos) toDo.set(2);
	if(requestType.find("jp")  != string::npos) toDo.set(3);
	if(requestType.find("lp")  != string::npos) toDo.set(4);
	if(requestType == string(""))								toDo.set(2); //Default
	
	bitset<4> grToDo;	
	if(sGrToDo == string("")){ 	//Default
		grToDo.set(1); 
		grToDo.set(3); 
	}
	else if(sGrToDo.find("all") != string::npos)
		grToDo.set();
	else{
		if(sGrToDo.find('0') != string::npos) grToDo.set(0); //set to 1
		if(sGrToDo.find('1') != string::npos) grToDo.set(1);	
		if(sGrToDo.find('2') != string::npos) grToDo.set(2);
		if(sGrToDo.find('3') != string::npos) grToDo.set(3);
	/*if(sGrToDo.find('4') != string::npos) grToDo.set(4);
		if(sGrToDo.find('5') != string::npos) grToDo.set(5);*/
	}
	
	cout<<"toDo: "<<toDo<<" \tgrToDo: "<<grToDo<<'\n';
	
//################################	EVENT PLOTS ################################
	if(toDo.test(0)){
		
		TH1I* hSig = openGraph<TH1I>(signalFile, "Cuts", sigName);
		
		//	#####	Struct creation #####
		GNames theNames("Cuts", sigName, backgrounds);
		MyGraphs<TH1F> theGraphs = *buildMyGraphs(theNames, signalFile, backgroundFiles);
		//	#####	--------------- #####
			
		compareCutGraphs(theGraphs, theNames);	//Integration makes no sense here
	}
	
//################################	JET PLOTS ################################
	if(toDo.test(1)){
		foreach(TString& graphName, jetPlots){
			if(reqGrName != string(""))
				if(reqGrName.find(graphName) == string::npos && reqGrName.find("all") == string::npos) continue;
			
			//	#####	Struct creation #####
			GNames theNames(graphName, sigName, backgrounds);
			MyGraphs<TH1F> theGraphs = *buildMyGraphs(theNames, signalFile, backgroundFiles);
			//	#####	--------------- #####
			
			newPrintGraph(theGraphs, theNames, grToDo);
		}
	}
	
//################################	LEPTON PLOTS ################################
	if(toDo.test(2)){
		foreach(TString& graphName, leptonPlots){
			TString copy = graphName;			copy.ToLower();
			if(reqGrName != string(""))
				if(reqGrName.find(graphName) == string::npos && reqGrName.find("all") == string::npos) continue;
			
			//	#####	Struct creation #####
			GNames theNames(graphName, sigName, backgrounds);
			MyGraphs<TH1F> theGraphs = *buildMyGraphs(theNames, signalFile, backgroundFiles);
			//	#####	--------------- #####
			
			newPrintGraph(theGraphs, theNames, grToDo);
		}
	}
	
//################################	JET PRE-CUT ################################
	if(toDo.test(3)){
		foreach(TString& graphName, jetAll){
			TString copy = graphName;			copy.ToLower();
			if(reqGrName != string(""))
				if(reqGrName.find(graphName) == string::npos && reqGrName.find("all") == string::npos) continue;
			
			//	#####	Struct creation #####
			GNames theNames(graphName, sigName, backgrounds);
			MyGraphs<TH1F> theGraphs = *buildMyGraphs(theNames, signalFile, backgroundFiles);
			//	#####	--------------- #####
			
			newPrintGraph(theGraphs, theNames, grToDo);
		}
	}

//################################	LEPTON PRE-CUT ################################
	if(toDo.test(4)){
		foreach(TString& graphName, lepAll){
			TString copy = graphName;			copy.ToLower();
			if(reqGrName != string(""))
				if(reqGrName.find(graphName) == string::npos && reqGrName.find("all") == string::npos) continue;
			
			//	#####	Struct creation #####
			GNames theNames(graphName, sigName, backgrounds);
			MyGraphs<TH1F> theGraphs = *buildMyGraphs(theNames, signalFile, backgroundFiles);
			//	#####	--------------- #####
			
			newPrintGraph(theGraphs, theNames, grToDo);
		}
	}
	
	cout<<'\n';
}



// ################################################################################################
TFile* openRootFile(const TString& path, const TString& name){
	cout<<"\nOpening \""<<name<<".root\"";
	TFile* file;
	try{
		file = TFile::Open(path + name + ".root");
	} catch(exception e) {cout<<"Exception: "<<e.what()<<'\n'; }
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


// ########################## Utils ####################################

TH1F* integralGraph(const TH1F* const origin){
	TString name = origin->GetName();
	TH1F* newGraph = (TH1F*)origin->Clone((name+" I").Data());
	int nbins = origin->GetNbinsX();
	for(int i = 1; i <= nbins; i++)
		newGraph->SetBinContent(i, origin->Integral(i,-1));
	
	return newGraph;
}

TH1F* reverseIntegrGraph(const TH1F* const origin){
	if(origin == nullptr) return nullptr;
	TString name = origin->GetName();
	TH1F* newGraph = (TH1F*)origin->Clone((name+" rI").Data());
	int nbins = origin->GetNbinsX();
	for(int i = nbins; i >= 1; i--)
		newGraph->SetBinContent(i, origin->Integral(1, i));
	
	return newGraph;
}

TH1F* sqrtGraph(const TH1F* const origin){
	if(origin == nullptr) return nullptr;
	TString name = origin->GetName();
	TH1F* newGraph = (TH1F*)origin->Clone((name+" Sqrt").Data());
	int nbins = origin->GetNbinsX();
	for(int i = 1; i <= nbins; i++){
		float s = sqrt(newGraph->GetBinContent(i));
		newGraph->SetBinContent(i, s);
	}
	return newGraph;
}

// ########################## Graphs ####################################
template<class TH = TH1F>
void compareCutGraphs(const MyGraphs<TH>& theGraphs, const GNames& names){
	cout<<"\n----compareCutGraphs: "<<names.graphName;
	//Check range
	TAxis* sigXAxis = theGraphs.hSig->GetXaxis();
	Double_t xm = sigXAxis->GetBinLowEdge(1);
	Int_t nbins = sigXAxis->GetNbins();
	Double_t xM = sigXAxis->GetBinLowEdge(nbins+1);
	
	//1) cuts sig/sqrt(bkg)
	TCanvas* cCut = new TCanvas("Cuts sig/#sqrt{bkg}", "Cuts sig/#sqrt{bkg}", 10,0,1280,1024);
	TH* hSigC = (TH*)theGraphs.hSig->Clone("Cuts sig/#sqrt{bkg}");
	TH* hBkgSum = new TH("Cuts bkg","Cuts bkg", nbins, xm, xM);
	for(int i = 1; i<= nbins; i++)
		hBkgSum->GetXaxis()->SetBinLabel(i, sigXAxis->GetBinLabel(i));
	
	//2) cuts over same graph (not scaled)
	TCanvas* cCutSame = new TCanvas("Cuts same", "Cuts same", 10,0,1280,1024);
	//THStack* stackNotScaled = new THStack("Cuts (not scaled)","Cuts (not scaled)");
	TH* hSigC2 = (TH*)theGraphs.hSig->Clone("Cuts same");
	hSigC2->SetTitle("Cuts (not scaled)");
	hSigC2->SetLineColor(myColors.at(0));
	hSigC2->SetFillColor(kRed-4);
	hSigC2->Draw();//stackNotScaled->Add(hSigC2);
	
	TLegend* legendCut = new TLegend(0.78, 0.94, 0.98, 0.74);
	legendCut->AddEntry(hSigC2, names.sigGrName);
	
	//3) cuts over same graph (scaled)
	TCanvas* cCutScaled = new TCanvas("Cuts same norm", "Cuts same norm", 10,0,1280,1024);
	//THStack* stackScaled = new THStack("Cuts (scaled)","Cuts (scaled)");
	TH* hSigCNorm = (TH*)theGraphs.hSig->Clone("Cuts sig norm");
	hSigCNorm->Scale(1./hSigCNorm->GetBinContent(1));
	hSigCNorm->SetTitle("Cuts (scaled)");
	hSigCNorm->SetLineColor(myColors.at(0));
	hSigCNorm->SetFillColor(kRed-4);
	hSigCNorm->Draw();//stackScaled->Add(hSigCNorm);
	
	TLegend* legendCutNorm = new TLegend(0.78, 0.94, 0.98, 0.74);
	legendCutNorm->AddEntry(hSigCNorm, names.sigGrName);
	
	
	for(int i = 0; i < theGraphs.vhBkgs.size(); i++){
		//3: scaled
		TH* cgNorm = (TH*)(theGraphs.vhBkgs.at(i)->Clone(names.bkgGrNames.at(i)+" norm"));
		cgNorm->Scale(1./cgNorm->GetBinContent(1));
		cgNorm->SetLineColor(myColors.at(1+i));
		cgNorm->SetFillColor(1+i);//cgNorm->SetFillColorAlpha(1+i,35);
		cCutScaled->cd();
		cgNorm->Draw("SAME");//stackScaled->Add(cgNorm);
		legendCutNorm->AddEntry(cgNorm, names.bkgGrNames.at(i));
		
		//--------------------------------------------------------
		TH* cg     = (TH*)(theGraphs.vhBkgs.at(i)->Clone(names.bkgGrNames.at(i)+" C"));
		cg->SetTitle("Cuts same "+names.bkgGrNames.at(i));
		cg->SetLineColor(myColors.at(1+i));
		//2: not scaled
		cCutSame->cd();
		cg->Draw("SAME");//stackNotScaled->Add(cg);
		legendCut->AddEntry(cg, names.bkgGrNames.at(i));
		//1: sig/sqrt{bkg}
		hBkgSum->Add(cg);
		cout<<"\nAdded: \""<<names.bkgGrNames.at(i)<<'\"';
	}
	
	//1: sig/sqrt{bkg}
	cCut->cd();
	hSigC->Divide(sqrtGraph(hBkgSum));
	hSigC->SetTitle("Cuts sig/#sqrt{bkg}");
	hSigC->SetMinimum(0.);
	hSigC->Draw();
	
	//2: not scaled
	cCutSame->cd();
	//stackNotScaled->Draw();
	legendCut->Draw("SAME");
	
	//3: scaled
	cCutScaled->cd();
	//stackScaled->Draw();
	legendCutNorm->Draw("SAME");
}


template <class TH>
void printGraphSame(const MyGraphs<TH>& graphs, const GNames& names){
	cout<<"\n----printGraphSame: "<<names.graphName;
	TCanvas *c0 = new TCanvas(names.graphName+" same", names.graphName+" same", 10,0,1280,1024);
	TLegend* legend0 = new TLegend(0.78, 0.94, 0.98, 0.74);
	
	for(int i = 0; i < graphs.vhBkgs.size(); i++){
		TH* cg = (TH*)(graphs.vhBkgs.at(i)->Clone(names.bkgGrNames.at(i)));
		//vhBkgCs.push_back(cg);
		cg->SetTitle(names.graphName+" same");
		cg->SetLineColor(1+myColors.at(i));
		legend0->AddEntry(cg, names.bkgGrNames.at(i));
		cg->Draw("SAME");
		//cout<<"\nAdded: \""<<names.bkgGrNames.at(i)<<'\"';
	}
	
	TH* hSigC = (TH*)graphs.hSig->Clone(names.sigGrName);
	hSigC->SetTitle(names.sigGrName);
	hSigC->SetLineColor(myColors.at(0));
	hSigC->Draw("SAME");
	legend0->AddEntry(hSigC, names.sigGrName);
	
	legend0->Draw("SAME");
	//c0->BuildLegend(0.7,0.8,0.95,0.95);	//It's a kind of magic
}

template <class TH>
void printGraphCompare(const MyGraphs<TH>& graphs, const GNames& names){
	cout<<"\n----printGraphCompare: "<<names.graphName;
	TCanvas *c1 = new TCanvas(names.graphName+" comp", names.graphName+" comp", 10,0,1280,1024);
	THStack* stack1 = new THStack(names.graphName+" comp",names.graphName+" comp");
	
	/*cout<<'\n'<<graphs.toString();
	cout<<'\n'<<names.toString();*/
	
	//vector<TH*> vhBkgCs;
	for(int i = 0; i < graphs.vhBkgs.size(); i++){
		TH* cg = (TH*)(graphs.vhBkgs.at(i)->Clone(names.bkgGrNames.at(i)));
		//vhBkgCs.push_back(cg);
		cg->SetTitle(names.bkgGrNames.at(i));
		cg->SetLineColor(1+myColors.at(i));
		stack1->Add(cg);
		//cout<<"\nAdded: \""<<names.bkgGrNames.at(i)<<'\"';
	}
	
	TH* hSigC = (TH*)graphs.hSig->Clone(names.sigGrName);
	hSigC->SetTitle(names.sigGrName);
	hSigC->SetLineColor(myColors.at(0));
	stack1->Add(hSigC);
	
	stack1->Draw();
	c1->BuildLegend(0.7,0.8,0.95,0.95);	//It's a kind of magic
}

template <class TH>
void printGraphIntegr(const MyGraphs<TH>& graphs, const GNames& names){
	cout<<"\n----printGraphIntegr: "<<names.graphName;
	TCanvas *cI = new TCanvas(names.graphName+" Integral", names.graphName+" Integral", 10,0,1280,1024);
	THStack* stackIS = new THStack(names.graphName+" Integral",names.graphName+" Integral");
	
	for(int i = 0; i < graphs.vhBkgs.size(); i++){
		TH* cg = (TH*)( graphs.vhBkgInt.at(i)->Clone(names.bkgGrNames.at(i)+" Integral") ); //Copy integral of background
		cg->SetTitle(names.bkgGrNames.at(i)+" Integral");
		cg->SetLineColor(1+myColors.at(i));
		stackIS->Add(cg);
		//cout<<"\nAdded: \""<<names.bkgGrNames.at(i)<<'\"';
	}
	
	TH* hSigIntC = (TH*)graphs.hSig->Clone(names.sigGrName+" Integral"); //Copy integral of signal
	hSigIntC->SetTitle(names.sigGrName+" Integral");
	hSigIntC->SetLineColor(myColors.at(0));
	stackIS->Add(hSigIntC);
	
	cI->cd();
	stackIS->Draw();
	cI->BuildLegend(0.7,0.8,0.95,0.95);
}

template <class TH>
void printGraphSqrt(const MyGraphs<TH>& graphs, const GNames& names){
	cout<<"\n----printGraphSqrt: "<<names.graphName;
	TCanvas *cIS = new TCanvas(names.graphName+" I(S)/#sqrt{I(B)}", names.graphName+" I(S)/#sqrt{I(B)}", 10,0,1280,1024);
	
	TH* hSigIntC = (TH*)graphs.hSig->Clone(names.graphName+" I(S)/#sqrt{I(B)}"); //Copy integral of signal
	hSigIntC->SetTitle(names.graphName+" I(S)/#sqrt{I(B)}");
	hSigIntC->SetMinimum(0.);
	
	//Check range
	TAxis* sigXAxis = graphs.hSigInt->GetXaxis();
	Double_t xm = sigXAxis->GetBinLowEdge(1);
	Int_t nbins = sigXAxis->GetNbins();
	Double_t xM = sigXAxis->GetBinLowEdge(nbins+1);

	TH* hBkgSqrt = new TH(names.graphName+" #sqrt{I(allBkgs)}",names.graphName+" #sqrt{I(allBkgs)}", nbins, xm, xM);
	
	for(int i = 0; i < graphs.vhBkgs.size(); i++){
		TH* cg = (TH*)( graphs.vhBkgInt.at(i)->Clone(names.bkgGrNames.at(i)+" Integral") ); //Copy integral of background
		cg->SetTitle(names.bkgGrNames.at(i)+" Integral");
		cg->SetLineColor(1+i);
		hBkgSqrt->Add(cg);
		//cout<<"\nAdded: \""<<names.bkgGrNames.at(i)<<'\"';
	}
	
	hSigIntC->Divide(sqrtGraph(hBkgSqrt));
	
	cIS->cd();
	hSigIntC->SetFillColor(kBlue-7);
	hSigIntC->Draw();
	//cIS->BuildLegend(0.7,0.8,0.95,0.95);
}

// ########################## All the graphs in one place ####################################
template <size_t N, class TH=TH1F>
void newPrintGraph(const MyGraphs<TH>& theGraphs, const GNames& theNames, const bitset<N>& grToDo){
	
	if(grToDo.test(0))
		printGraphSame(theGraphs, theNames);
		
	if(grToDo.test(1))
		printGraphCompare(theGraphs, theNames);
		
	if(grToDo.test(2))
		printGraphIntegr(theGraphs,theNames);
		
	if(grToDo.test(3))
		printGraphSqrt(theGraphs,theNames);
	
	return;
}

template <class TH = TH1F>
MyGraphs<TH>* buildMyGraphs(const GNames& names, TFile* const signalFile, vector<TFile*>& backgroundFiles){
	if(signalFile == nullptr) return nullptr;
	TH* hSig = openGraph(signalFile, names.graphName, names.sigGrName);
	vector<TH*> vhBackgrounds;
	vector<TH*> vhBackgrIntegr;
	for(int i = 0; i < backgroundFiles.size(); i++){
		TH1F* pippo = openGraph(backgroundFiles.at(i), names.graphName, names.bkgGrNames.at(i));
		vhBackgrounds.push_back( pippo );
		vhBackgrIntegr.push_back( integralGraph(pippo) );
	}
	return new MyGraphs<TH>(hSig, vhBackgrounds, integralGraph(hSig), vhBackgrIntegr);
}

