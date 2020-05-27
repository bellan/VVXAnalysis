#include "VVXAnalysis/TreeAnalysis/interface/DBGenerator.h"
#include "VVXAnalysis/TreeAnalysis/interface/VZZAnalyzer.h"
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
	
	startTime_ = clock();
	
	//Memory allocation
	if(AK4pairs_ == nullptr)    AK4pairs_ = new vector<Boson<Jet>>;
	if(genHadVBs_ == nullptr)  genHadVBs_ = new vector<Boson<Particle>>;
	//if(AllGenVBjj_ == nullptr) AllGenVBjj_= new vector<Boson<Particle>>;
	//if(!AK4GenRec_) AK4GenRec_ = new vector<pair<Boson<Particle>, Boson<Jet>>>;
	
	cout<<"\n\tfilename = \""<<fileName<<"\"\n"; //fileName comes from EventAnalyzer.h
	
	// format: samples/<year>/<sample_name>.root
	
	std::string::size_type startOfName = fileName.find_last_of("/");
	std::string::size_type endOfName = fileName.find_last_of(".");
	std::string::size_type nameLength = fileName.length();
	
	std::string::size_type startOfYear = fileName.find_first_of("/");
	std::string::size_type endOfYear = fileName.find_last_of("/");
	
	//cout<<"\tstart = "<<startOfName<<" \tend = "<<endOfName<<" \ttotal length = "<<nameLength<<"\n";
	std::string strippedName = "UnnamedTree";
	if(startOfName != string::npos  &&  endOfName != string::npos){
		strippedName = fileName.substr(startOfName+1, endOfName-startOfName-1);
		cout<<"\tname = "<<strippedName<<' ';
	}
	
	std::string year_str = "0";
	if(startOfYear != string::npos  &&  endOfYear != string::npos){
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
	
	isSigFile_ = false;
	foreach(const std::string& name, signalNames_)
		if(TString(fileName).Contains(name.c_str())){
			isSigFile_ = true;
			break;
		}
	
	return;
}


Int_t DBGenerator::cut(){
	++evtN_;
	cout<<"\r\t\t"<<evtN_;
	
	//Cleanup of previous event
	sigVB_ = nullptr;
	if(genHadVBs_)  genHadVBs_->clear(); //Destroys objects but keeps memory allocated
	if(AK4pairs_)    AK4pairs_->clear();
	//if(AK4GenRec_)  AK4GenRec_->clear();
	//if(AllGenVBjj_)AllGenVBjj_->clear();
	
	//Preparation for this event
	fillGenHadVBs();  // -->genHadVBs_
	fillRecHadVBs();  // -->AK4pairs_
	//fillAK4GenRec();
	//fillGenVBtoAK4();
	
	sigType_ = VZZAnalyzer::isSignal(/*Uses genJets(AK4) and genJetsAK8*/);
	
	return 1;
}


void DBGenerator::analyze(){
	++analyzedN_; analyzedW_ += theWeight;	
	
	//Particle* candClosest =nullptr;//Contains information from the simulation. Probably useless
	const Particle* candBestV = nullptr;  // Obtainable looking only at data
	// As sigVB, it only points to the Particle, it doesn't own it. Maybe change this?
	
	switch(sigType_){
		case 1:
			sigVB_ = &(genHadVBs_->front());
			/*if(AK4pairs_->size() > 0){
				stable_sort(AK4pairs_->begin(), AK4pairs_->end(), DeltaRComparator(*sigVB_) );
				if(physmath::deltaR(AK4pairs_->front(), *sigVB_) < 0.4)
					candClosest = &(AK4pairs_->front());
			}*/ // else if(jetsAK8->size() > 0) goto ak8_label;
			break;
		case 2:
			sigVB_ = &(genJetsAK8->front());
			/*if(jetsAK8->size() != 0){
				stable_sort(jetsAK8->begin(), jetsAK8->end(), DeltaRComparator(*sigVB_) );
				if(physmath::deltaR(jetsAK8->front(), *sigVB_) < 0.4)
					candClosest = &(jetsAK8->front());
			}*/ // else if(AK4pairs_->size() > 0) goto ak4_label;
			break;
		case 0:
			break;
		default:
			cout<<"Error, sigType_ has an invalid value in analyze()\n";
			return;
	}
	
	//RECONSTUCTION: Move it to VZZ, make nonprivate member function
	VCandType temp = VCandType::None;
	int sigRecType = 0;
	Particle* cand4 = VZZAnalyzer::findBestVFromPair(jets, temp);  // = Particle -> Boson<Jet>*
	const Particle* cand8 = VZZAnalyzer::findBestVFromSing(jetsAK8, temp);//= Particle* -> Jet*
	if(!cand4 && !cand8) return;  // IMPORTANT: JUMP EVENT IF NO VB CAN BE RECONSTRUCTED!
	if(cand4 && cand8){
		if(VZZAnalyzer::minDM(cand4->mass()) < 30. || VZZAnalyzer::minDM(cand8->mass()) < 30.){
			if( VZZAnalyzer::minDM(cand4->mass()) < VZZAnalyzer::minDM(cand8->mass()) ){
				sigRecType = 1;
				candBestV = cand4;
			}
			else{
				sigRecType = 2;
				candBestV = cand8;
			}
		}
	}
	else if(cand4) { candBestV = cand4; 	sigRecType = 1; }
	else if(cand8) { candBestV = cand8; 	sigRecType = 2; }
	//END RECONSTRUCTION
	
	// IMPORTANT: JUMP EVENT IF NO VB CAN BE RECONSTRUCTED!
	if(!candBestV || !ZZ || ZZ->pt() < 1.) return;
	
	if(/*isSigFile_ &&*/ topology.test(4) && topology.test(0) && sigType_)
		outputFile<<1;  // 1  Is signal at gen level?
	else
		outputFile<<0;
		
	outputFile<<SEP_CHAR<<theWeight;  // 2
	
	//Here we write relevant variables to the outputFile
	outputFile<<SEP_CHAR<<met->pt();  // 3
	
	if(sigType_ && !sigVB_)
		cout<<evtN_<<"sigType_ && !sigVB_\n";
	
	mainEvtRec(sigRecType, candBestV);
	
	// GEN VARIABLES --> old
	/*
	const Particle* genSigVB = nullptr;
	switch(sigType){
		case 1:
			genSigVB = &(genHadVBs_->front());
			break;
		case 2:
			genSigVB = &(genJetsAK8->front());
			break;
	}
	
	if(ZZ && ZZ->pt() > 1.){  // 10out
		printVars(2, ZZ->mass(), ZZ->pt());  // 4, 5
		
		if(genSigVB){  // 8out
			printVars(3, genSigVB->pt(), genSigVB->e(), fabs(genSigVB->eta()));  // 6, 7, 8
			
			float sHat = (ZZ->p4() + genSigVB->p4()).M();
			printVars(1, sHat);  // 9
			
			float ptVZZ = (ZZ->p4() + genSigVB->p4()).Pt();
			printVars(1, ptVZZ);  // 10
			
			float dR1 = physmath::deltaR(*genSigVB,   ZZ->first());
			float dR2 = physmath::deltaR(*genSigVB,   ZZ->second());
			float dR3 = physmath::deltaR(ZZ->first(), ZZ->second());
			printVars(3, dR1, dR2, dR3);  // 11, 12, 13	
		}
		else{  // 8 out
			printZeroes(7);  // 6, ... 12
			printVars(1, physmath::deltaR(ZZ->first(), ZZ->second()));  // 13
		}
	}
	else
		printZeroes(10);  // 3, ... 13
	*/
	
	
	
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
	
	if(cand4) delete cand4;
	
	outputFile<<std::endl;
}

void DBGenerator::end(TFile &){
	outputFile.close();
	
	if(genHadVBs_) delete genHadVBs_; //Deallocates memory
	if(AK4pairs_)  delete AK4pairs_;
	
	cout<<"\nPassing cut: "<<Form("%lu (weighted: %.2f)", analyzedN_, analyzedW_)<<'\n';
	
	float elapsedSec = (float)(clock()-startTime_)/CLOCKS_PER_SEC;
	int elapsedSecInt = (int)elapsedSec;
	cout<<"\nElapsed Time: "<<elapsedSec<<" s\t\t("<<elapsedSecInt/60<<"\' "<<elapsedSecInt%60<<"\")\n";
	for(char i=0; i<25; ++i) cout<<"-";
	cout<<" \tEnd of DBGenerator\t";
	for(char i=0; i<25; ++i) cout<<"-";
	cout<<"\n\n";
}


void DBGenerator::mainEvtRec(int sigRecType, const Particle* recV){
	//recV is from findBestVfrom[...]()
	float dMHad = minDM(recV->mass());
	float ptHad = recV->pt();
	float aEtaHad = fabs(recV->eta());
	printVars(3, dMHad, ptHad, aEtaHad);
	
	//if(ZZ && ZZ->pt() > 1.){
	printVars(3, ZZ->mass(), ZZ->pt(), ZZ->eta());
	
		//if(recV){  // Either there's no rec VB or it's too far in dR
	TLorentzVector candp4(recV->p4());
	TLorentzVector Z1p4(ZZ->first().p4());
	TLorentzVector Z2p4(ZZ->second().p4());
	
	// Angles
	float ang0 = Z1p4.Angle(Z2p4.Vect());
	float ang1 = Z1p4.Angle(candp4.Vect());
	float ang2 = candp4.Angle(Z2p4.Vect());
	
	// phi
	float dPhi0 = physmath::deltaPhi(ZZ->first(), ZZ->second());
	float dPhi1 = physmath::deltaPhi(ZZ->first(), *recV);
	float dPhi2 = physmath::deltaPhi(*recV, ZZ->second());
	
	// eta
	float dEta0 = ZZ->first().eta() - ZZ->second().eta();
	float dEta1 = ZZ->first().eta() - recV->eta();
	float dEta2 = recV->eta() - ZZ->second().eta();
	
	// R
	float dR0 = physmath::deltaPhi(ZZ->first(), ZZ->second());
	float dR1 = physmath::deltaPhi(ZZ->first(), *recV);
	float dR2 = physmath::deltaPhi(*recV, ZZ->second());
	
	// projection - relative pt
	float pt0 = candp4.P() * sin( candp4.Angle(ZZ->p4().Vect()) );
	float pt1 = Z1p4.P() * sin( Z1p4.Angle((candp4 + Z2p4).Vect()) );
	float pt2 = Z2p4.P() * sin( Z2p4.Angle((candp4 + Z1p4).Vect()) );
	
	printVars(15, ang0,ang1,ang2, dPhi0,dPhi1,dPhi2, dEta0,dEta1,dEta2, dR0,dR1,dR2, pt0,pt1,pt2);
	
	printVars(2, (candp4+ZZ->p4()).M(), (candp4+ZZ->p4()).Pt());
	
		//}
		/*else{
			TLorentzVector Z1p4(ZZ->first().p4());
			TLorentzVector Z2p4(ZZ->second().p4());
			
			float ang0 = Z1p4.Angle(Z2p4.Vect());
			float dPhi0 = physmath::deltaPhi(ZZ->first(), ZZ->second());
			float dEta0 = ZZ->first().eta() - ZZ->second().eta();
			float dR0 = physmath::deltaPhi(ZZ->first(), ZZ->second());
			
			printVars(15, ang0,0.,0., dPhi0,0.,0., dEta0,0.,0., dR0,0.,0., 0.,0.,0.);
		}*/
	//}
	//else{
		//printZeroes(15);
	//}
}


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
