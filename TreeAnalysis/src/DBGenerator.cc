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


void DBGenerator::begin(){
	cout<<"\n";
	for(char i=0; i<25; ++i) cout<<"-";
	cout<<" \tStart of DBGenerator\t";
	for(char i=0; i<25; ++i) cout<<"-";
	
	startTime = clock();
	
	
	cout<<"\n\tfilename = \""<<fileName<<"\"\n"; //fileName comes from EventAnalyzer.h
	
	std::string::size_type startOfName = fileName.find_last_of("/");
	std::string::size_type endOfName = fileName.find_last_of(".");
	std::string::size_type nameLength = fileName.length();
	
	std::string::size_type startOfYear = fileName.find_first_of("/");
	std::string::size_type endOfYear = fileName.find_last_of("/");
	
	cout<<"\tstart = "<<startOfName<<" \tend = "<<endOfName<<" \ttotal length = "<<nameLength<<"\n";
	std::string strippedName = "UnnamedTree";
	if(startOfName!=string::npos && endOfName!=string::npos){
		strippedName = fileName.substr(startOfName+1, endOfName-startOfName-1);
		cout<<"\tname = "<<strippedName<<' ';
	}
	
	std::string year_str = "0";
	if(startOfYear!=string::npos && endOfYear!=string::npos){
		year_str = fileName.substr(startOfYear+1, endOfYear-startOfYear-1);
		cout<<"\tyear = "<<year_str<<' ';
	}
	
	
	if(gSystem->AccessPathName("databases")) //"Bizarre convention": returns false if path exists
		gSystem->MakeDirectory("databases");
	
	string outPath = "databases/"+strippedName+'_'+year_str+".csv";
	cout<<"\toutName = "<<outPath<<"\n";
	outputFile.open(outPath, std::ios::trunc);
	if(!outputFile.is_open())
		exit(0);
	
	cout<<"\nAnalyzed:\t      /"<<tree()->GetEntries()<<std::flush;
	
	//Memory allocation
	if(genHadVBs_ == nullptr)  genHadVBs_ = new vector<Boson<Particle>>;
	
	isSigFile_ = false;
	foreach(const std::string& name, signalNames_)
		if(TString(fileName).Contains(name.c_str())){
			isSigFile_ = true;
			break;
		}
	
	return;
}


Int_t DBGenerator::cut(){
	++evtN;
	cout<<"\r\t\t"<<evtN;
	
	if(genHadVBs_)  genHadVBs_->clear();
	fillGenHadVBs();
	
	return 1;
}


void DBGenerator::analyze(){
	int sigType = VZZAnalyzer::isSignal();
	if(/*isSigFile_ &&*/ topology.test(4) && topology.test(0) && sigType)
		outputFile<<1;  // 1  Is signal at gen level?
	else
		outputFile<<0;
		
	outputFile<<SEP_CHAR<<theWeight;  // 2
	
	//Here we write relevant variables to the outputFile
	outputFile<<SEP_CHAR<<met->pt();  // 3
	
	const Particle* genSigVB = nullptr;
	switch(sigType){
		case 1:
			genSigVB = &(genHadVBs_->front());
			break;
		case 2:
			genSigVB = &(genJetsAK8->front());
			break;
	}
	
	if(genSigVB){
		printVars(1, genSigVB->pt());  // 4
		
		if(ZZ && ZZ->pt() > 1.){
			float sHat = (ZZ->p4() + genSigVB->p4()).M();
			printVars(1, sHat);  // 5
			
			float dR1 = physmath::deltaR(*genSigVB,   ZZ->first());
			float dR2 = physmath::deltaR(*genSigVB,   ZZ->second());
			float dR3 = physmath::deltaR(ZZ->first(), ZZ->second());
			printVars(3, dR1, dR2, dR3);  // 6, 7, 8
		}
		else
			printZeroes(4);
	}
	else
		printZeroes(5);
		
	if(ZZ && ZZ->pt() > 1.){
		printVars(ZZ->mass(), ZZ->pt());
	}
	else
		printZeroes(2);
	
	
	/*	//Other jets-related variables
	VCandType candType = VCandType::None; //initialization
	phys::Boson<phys::Jet>* VCandidate = findBestVFromPair(jets, candType);
	if(VCandidate){
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
		printZeroes(6);*/
	
	outputFile<<std::endl;
}

void DBGenerator::end(TFile &){
	outputFile.close();
	
	if(genHadVBs_) delete genHadVBs_; //Deallocates memory
	
	float elapsedSec = (float)(clock()-startTime)/CLOCKS_PER_SEC;
	int elapsedSecInt = (int)elapsedSec;
	cout<<"\nElapsed Time: "<<elapsedSec<<" s\t\t("<<elapsedSecInt/60<<"\' "<<elapsedSecInt%60<<"\")\n";
	for(char i=0; i<25; ++i) cout<<"-";
	cout<<" \tEnd of DBGenerator\t";
	for(char i=0; i<25; ++i) cout<<"-";
	cout<<"\n\n";
}

/*
int DBGenerator::isSignal(){ //isHadSignal  //TODO: link to equivalent function in VZZAnalyzer
	bool twoAK4 = ( genHadVBs_->size() >= 1 );
	bool oneAK8 = ( genJetsAK8->size() >= 1 );
	if(!twoAK4 && !oneAK8)
		return 0;
	
	if(twoAK4){
		std::stable_sort(genHadVBs_->begin(), genHadVBs_->end(), Mass2Comparator(phys::ZMASS, phys::WMASS));
		float mass4 = genHadVBs_->front().mass();
		if(phys::WMASS-30. < mass4 && mass4 < phys::ZMASS+30.)
			return 1;
	}
	
	if(oneAK8){
		std::stable_sort(genJetsAK8->begin(), genJetsAK8->end(), Mass2Comparator(phys::ZMASS, phys::WMASS));
		float mass8 = genJetsAK8->front().mass();
		if(phys::WMASS-30. < mass8 && mass8 < phys::ZMASS+30.)
			return 2;
	}
		
	return 0;
}
*/

void DBGenerator::printZeroes(size_t nzeros){
	for(size_t i = 0; i < nzeros; ++i){
		outputFile << SEP_CHAR << 0;
	}
}
void DBGenerator::printVars(size_t n, ...){
	double val;
  va_list vl;
  va_start(vl, n);
  for (size_t i = 0; i < n; ++i){
		val = va_arg(vl, double); //For historic reasons floats are promoted to double anyway
		outputFile << SEP_CHAR << (float)val;
	}
  va_end(vl);
}

/*
void DBGenerator::fillGenHadVBs(){
	foreach(const Boson<Particle>& genVB, *genVBParticles)
		if( abs(genVB.daughter(0).id()) < 10 && abs(genVB.daughter(1).id()) < 10 )
			genHadVBs_->push_back(genVB); //The second condition is redundant
}


template <class J = phys::Jet>
phys::Boson<J>* DBGenerator::findBestVFromPair(const std::vector<J>* js, VCandType& thisCandType){
	thisCandType = VCandType::None;
	if(js->size() < 2)
		return nullptr;
		
	pair<size_t, size_t> indicesZ(0,0);
	pair<size_t, size_t> indicesW(0,0);
	float minDifZ = 50.;
	float minDifW = 50.;
	float tmpMass = 0.;
	for(size_t i = 0; i < js->size(); ++i){
		for(size_t j = i+1; j < js->size(); ++j){
			tmpMass = (js->at(i).p4() + js->at(j).p4()).M();
			float diffZa = fabs(tmpMass - phys::ZMASS);
			float diffWa = fabs(tmpMass - phys::WMASS);
			if(diffZa < minDifZ){
				minDifZ = diffZa;
				indicesZ = std::make_pair(i,j);
			}
			if(diffWa < minDifW){
				minDifW = diffWa;
				indicesW = std::make_pair(i,j);
			}
		}
	}
	
	phys::Boson<J>* thisCandidate = nullptr;
	if(minDifZ < minDifW){
		thisCandidate = new Boson<J>(js->at(indicesZ.first), js->at(indicesZ.second));
		if(!ZBosonDefinition(*thisCandidate))
			return nullptr;
		thisCandType = VCandType::Z;
	}
	else{
		thisCandidate = new Boson<J>(js->at(indicesW.first), js->at(indicesW.second));
		if(!WBosonDefinition(*thisCandidate))
			return nullptr;
		thisCandType = VCandType::W;
	}
	return thisCandidate;
}
*/
