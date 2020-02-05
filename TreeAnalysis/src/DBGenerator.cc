#include "VVXAnalysis/TreeAnalysis/interface/DBGenerator.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <iostream>
#include <fstream>			//open(), close(), <<
#include <string>				//find_last_of()
#include <time.h>				//clock_t, clock()
#include "TSystem.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;

using namespace phys;

Int_t DBGenerator::cut(){
	evtN++;
	cout<<"\r\t\t"<<evtN;
	return 1;
}

void DBGenerator::begin(){
	cout<<"\n";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<" \tStart of DBGenerator\t";
	for(char i=0; i<25; i++) cout<<"-";
	
	startTime = clock();
	
	cout<<"\n\tfilename = \""<<fileName<<"\"\n"; //fileName comes from EventAnalyzer.h
	std::string::size_type startOfName = fileName.find_last_of("/");
	std::string::size_type endOfName = fileName.find_last_of(".");
	std::string::size_type nameLength = fileName.length();
	cout<<"\tstart = "<<startOfName<<" \tend = "<<endOfName<<" \ttotal length = "<<nameLength<<"\n";
	string strippedName = "UnnamedTree";
	if(startOfName!=string::npos && endOfName!=string::npos){
		strippedName = fileName.substr(startOfName+1, endOfName-startOfName-1);
		cout<<"\tname = "<<strippedName<<' ';
	}
	
	if(gSystem->AccessPathName("databases")) //"Bizarre convention": returns false if path exists
		gSystem->MakeDirectory("databases");
	
	string outPath = "databases/"+strippedName+".txt";
	cout<<"\toutName = "<<outPath<<"\n";
	outputFile.open(outPath, std::ios::trunc);
	if(!outputFile.is_open())
		exit(0);
	
	cout<<"\nAnalyzed:\t      /"<<tree()->GetEntries()<<std::flush;
	return;
}

void DBGenerator::analyze(){
	//Here we write relevant variables to the outputFile
	outputFile<<met->pt()<<SEP_CHAR<<met->e();
	
	//Writing electrons variables
	outputFile<<SEP_CHAR<<electrons->size();
	if(electrons->size()>=1){
		outputFile<<SEP_CHAR<<electrons->at(0).pt()<<SEP_CHAR<<electrons->at(0).e();
		if(electrons->size()>=2)
			outputFile<<SEP_CHAR<<electrons->at(1).pt()<<SEP_CHAR<<electrons->at(1).e();
		else
			printZeros(2);
	} else
		printZeros(4);
	
	//Writing muons variables
	outputFile<<SEP_CHAR<<muons->size();
	if(muons->size()>=1){
		outputFile<<SEP_CHAR<<muons->at(0).pt()<<SEP_CHAR<<muons->at(0).e();
		if(muons->size()>=2)
			outputFile<<SEP_CHAR<<muons->at(1).pt()<<SEP_CHAR<<muons->at(1).e();
		else
			printZeros(2);
	} else
		printZeros(4);
	
	//Writing jets variables
	if(jets->size()>=1){
		outputFile<<SEP_CHAR<<jets->at(0).pt()<<SEP_CHAR<<jets->at(0).e();
		if(jets->size()>=2)
			outputFile<<SEP_CHAR<<jets->at(1).pt()<<SEP_CHAR<<jets->at(1).e();
		else
			printZeros(2);
	} else
		printZeros(4);
	
	//Other jets-related variables
	if(jets->size()>=2){
		TLorentzVector p4JJ = jets->at(0).p4() + jets->at(1).p4();
		TLorentzVector p4Tot = p4JJ + ZZ->p4();
		float ptTot = p4Tot.Pt();
		float mTot = p4Tot.M(); //Doesn't make much sense, but let's try anyway
		float deltaEtaTot = fabs(p4JJ.Eta() - ZZ->eta());
		float deltaRTot = ZZ->p4().DeltaR(p4JJ); //TLorentzVector::DeltaR()
		float mWJJ_norm = fabs(p4JJ.M() - phys::WMASS)/phys::WMASS;
		float mWJJ_norm_sq = mWJJ_norm*mWJJ_norm;
		outputFile<<SEP_CHAR<< ptTot<<SEP_CHAR<< mTot<<SEP_CHAR<< deltaEtaTot<<SEP_CHAR<< deltaRTot<<SEP_CHAR<< mWJJ_norm<<SEP_CHAR<<mWJJ_norm_sq;
	}
	else
		printZeros(6);
	
	outputFile<<std::endl;
}

void DBGenerator::end(TFile &){
	outputFile.close();
	
	float elapsedSec = (float)(clock()-startTime)/CLOCKS_PER_SEC;
	int elapsedSecInt = (int)elapsedSec;
	cout<<"\nElapsed Time: "<<elapsedSec<<" s\t\t("<<elapsedSecInt/60<<"\' "<<elapsedSecInt%60<<"\")\n";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<" \tEnd of DBGenerator\t";
	for(char i=0; i<25; i++) cout<<"-";
	cout<<"\n\n";
}

void DBGenerator::printZeros(size_t nzeros){
	for(size_t i = 0; i < nzeros; i++){
		outputFile<<SEP_CHAR<<0;
	}
}
