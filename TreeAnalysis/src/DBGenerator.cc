#include "VVXAnalysis/TreeAnalysis/interface/DBGenerator.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <iostream>
#include <fstream>			//open(), close(), <<
#include <string>				//find_last_of()
#include <time.h>				//clock_t, clock()
#include <utility>			//std::pair(), std::make_pair()
#include <stdarg.h>     //va_list, va_start, va_arg, va_end
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
using std::pair;

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
	
	string outPath = "databases/"+strippedName+".csv";
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
	/*
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
	*/
	//Writing jets variables
	if(jets->size()>=1){
		printVars(2, jets->at(0).pt(), jets->at(0).e());
		if(jets->size()>=2)
			printVars(2, jets->at(1).pt(), jets->at(1).e());
		else
			printZeros(2);
	} else
		printZeros(4);
	
	//Other jets-related variables
	phys::Boson<phys::Jet> VCandidate; //initialization
	VCandType candType = VCandType::None; //initialization
	bool weHaveACand = findBestVCandidate(jets, VCandidate, candType);
	if(weHaveACand){
		TLorentzVector p4JJ = VCandidate.p4();
		TLorentzVector p4Tot = p4JJ + ZZ->p4();
		float ptTot = p4Tot.Pt();
		float mTot = p4Tot.M(); //Doesn't make much sense, but let's try anyway
		float deltaEtaTot = fabs(p4JJ.Eta() - ZZ->eta());
		float deltaRTot = ZZ->p4().DeltaR(p4JJ); //TLorentzVector::DeltaR()
		float mWJJ_norm = fabs(p4JJ.M() - phys::WMASS)/phys::WMASS;
		float mWJJ_norm_sq = mWJJ_norm*mWJJ_norm;
		printVars(6, ptTot, mTot, deltaEtaTot, deltaRTot, mWJJ_norm, mWJJ_norm_sq);
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
		outputFile << SEP_CHAR << 0;
	}
}
void DBGenerator::printVars(int n, ...){
	double val;
  va_list vl;
  va_start(vl, n);
  for (int i = 0; i < n; i++){
		val = va_arg(vl, double); //For historic reasons floats are promoted to double anyway
		outputFile << SEP_CHAR << (float)val;
	}
  va_end(vl);
}


bool DBGenerator::findBestVCandidate(const std::vector<phys::Jet>* js, phys::Boson<phys::Jet>& thisCandidate, VCandType& thisCandType){
	bool isAccurate = false;
	if(js->size()>=2){
		pair<size_t, size_t> indices(0,0);
		float minDifZ = 1.;
		float minDifW = 1.;
		// float massZCand = 0.;
		// float massWCand = 0.;
		float tmpMass = 0.;
		for(size_t i = 0; i < js->size(); i++){
			for(size_t j = i+1; j < js->size(); j++){
				tmpMass = (js->at(i).p4() + js->at(j).p4()).M();
				float diffZ = (tmpMass - phys::ZMASS)/phys::ZMASS;
				float diffZa = fabs(diffZ);
				float diffW = (tmpMass - phys::WMASS)/phys::WMASS;
				float diffWa = fabs(diffW);
				if(diffZa < minDifZ){
					minDifZ = diffZa;
					//massZCand = tmpMass;
					indices = std::make_pair(i,j);
				}
				if(diffWa < minDifW){
					minDifW = diffWa;
					//massWCand = tmpMass;
					indices = std::make_pair(i,j);
				}
			}
		}
		thisCandidate = Boson<Jet>(js->at(indices.first), js->at(indices.second));
		
		if(minDifZ < minDifW){
			isAccurate = ZBosonDefinition(thisCandidate);
			thisCandType = VCandType::Z;
		}
		else{
			isAccurate = WBosonDefinition(thisCandidate);
			thisCandType = VCandType::W;
		}
	}
	else
		thisCandType = VCandType::None;
	return isAccurate;
}

