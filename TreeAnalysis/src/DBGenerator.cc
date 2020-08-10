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
	for(char i=0; i<25; ++i) cout<<'-';
	cout<<" \tStart of DBGenerator\t";
	for(char i=0; i<25; ++i) cout<<'-';
	
	startTime_ = clock();
	
	//Memory allocation
	if(AK4pairs_ == nullptr)    AK4pairs_ = new vector<Boson<Jet>>;
	if(genHadVBs_ == nullptr)  genHadVBs_ = new vector<Boson<Particle>>;
	//if(AllGenVBjj_ == nullptr) AllGenVBjj_= new vector<Boson<Particle>>;
	//if(!AK4GenRec_) AK4GenRec_ = new vector<pair<Boson<Particle>, Boson<Jet>>>;
	
	
	isSigFile_ = false;
	foreach(const std::string& name, signalNames_)
		if(TString(fileName).Contains(name.c_str())){
			isSigFile_ = true;
			break;
		}
	
	// Name parsing	
	cout<<"\n\tfilename = \""<<fileName<<"\"\n"; //fileName comes from EventAnalyzer.h
	// format: samples/<year>/<sample_name>.root
	
	std::string::size_type startOfName = fileName.find_last_of('/');
	std::string::size_type endOfName = fileName.find_last_of('.');
	std::string::size_type nameLength = fileName.length();
	
	std::string::size_type startOfYear = fileName.find_first_of('/');
	std::string::size_type endOfYear = fileName.find_last_of('/');
	
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
	cout<<"\toutName = "<<outPath<<"\"\n";
	outputFile_.open(outPath, std::ios::trunc);
	if(!outputFile_.is_open())
		exit(1);
	outputFile_.sync_with_stdio(false);
	
	if(isSigFile_){
		string outPathAK4 = "databases/"+strippedName+'_'+year_str+"_TagAK4.csv";
		cout<<"\toutName (Jet Tagger AK4) = "<<outPathAK4<<"\"\n";
		outputFileAK4_.open(outPathAK4, std::ios::trunc);
		if(!outputFileAK4_.is_open())
			exit(1);
		
		string outPathAK8 = "databases/"+strippedName+'_'+year_str+"_TagAK8.csv";
		cout<<"\toutName (Jet Tagger AK8) = "<<outPathAK8<<"\"\n";
		outputFileAK8_.open(outPathAK8, std::ios::trunc);
		if(!outputFileAK8_.is_open())
			exit(1);
	}
	
	cout<<"\nAnalyzed:\t      /"<<tree()->GetEntries()<<std::flush;
	
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
	if(!candBestV || !ZZ || ZZ->isValid() || ZZ->passFullSelection()) return;
	
	if(/*isSigFile_ &&*/ topology.test(4) && topology.test(0) && sigType_)
		outputFile_<<1;  // 1  Is signal at gen level?
	else
		outputFile_<<0;
		
	outputFile_<<SEP_CHAR<<theWeight;  // 2
	
	//Here we write relevant variables to the outputFile
	outputFile_<<SEP_CHAR<<met->pt();  // 3
	
	
	mainEvtRec(sigRecType, candBestV);
	
	if(cand4) delete cand4;
	
	outputFile_<<'\n';  // Let the buffer be flushed automatically
	
	writeTagger(outputFileAK4_, outputFileAK8_);  // contains the endline
	
}

void DBGenerator::end(TFile &){
	outputFile_.close();
	if(outputFileAK4_.is_open()) outputFileAK4_.close();
	if(outputFileAK8_.is_open()) outputFileAK8_.close();
	
	if(genHadVBs_) delete genHadVBs_; //Deallocates memory
	if(AK4pairs_)  delete AK4pairs_;
	
	cout<<"\nPassing cut: "<<Form("%lu (weighted: %.2f)", analyzedN_, analyzedW_)<<'\n';
	
	float elapsedSec = (float)(clock()-startTime_)/CLOCKS_PER_SEC;
	int elapsedSecInt = (int)elapsedSec;
	cout<<"\nElapsed Time: "<<elapsedSec<<" s\t\t("<<elapsedSecInt/60<<"\' "<<elapsedSecInt%60<<"\")\n";
	for(char i=0; i<25; ++i) cout<<'-';
	cout<<" \tEnd of DBGenerator\t";
	for(char i=0; i<25; ++i) cout<<'-';
	cout<<"\n\n";
}


void DBGenerator::writeTagger(std::ofstream& outFile4, std::ofstream& outFile8){
	vector<Particle> genQuarks;
	size_t nLep = 0;
	
	foreach(const Particle& p, *genParticles){
		unsigned int aID = abs(p.id());
		if(aID < 10)       // Is it a quark?
			genQuarks.push_back(Particle(p));
		else if(aID < 20)
			++nLep;
	}
	
	// 2q final state
	if(genQuarks.size() != 2) // || nLep != 4)
		return;
	
	Boson<Particle> qq(genQuarks.at(0), genQuarks.at(1));
	
	//GenJets AK4  IMPORTANT: eff/res ONLY for the 4l 2q channel (no 4q/6q)!
	Jet rec4_1, rec4_2;
	Boson<Jet> bestAK4;
	float dR0 = 99.;
	float dR1 = 99.;
	if(jets->size() >= 2){
		std::sort(jets->begin(), jets->end(), phys::DeltaRComparator(qq.daughter(0)));
		dR0 = physmath::deltaR(jets->front(), qq.daughter(0));
		if(dR0 < 0.4)
			rec4_1 = Jet(jets->front());
		
		std::sort(jets->begin(), jets->end(), phys::DeltaRComparator(qq.daughter(1)));
		dR1 = physmath::deltaR(jets->front(), qq.daughter(1));
		if(dR1 < 0.4)
			rec4_2 = Jet(jets->front());
			
		if(dR0 < 0.4 && dR1 < 0.4){
			bestAK4 = Boson<Jet>(rec4_1, rec4_2);
			writeInfoAK4(bestAK4, outFile4);
			outFile4<<'\n';  // '\n' instead of std::endl because it does not force buffer flush
		}
	}
	
	
	//GenJets AK8  IMPORTANT: eff/res ONLY for the 4l 2q channel (no 4q/6q)!
	if(jetsAK8->size() > 0){
		std::sort(jetsAK8->begin(), jetsAK8->end(), phys::DeltaRComparator(qq));
		
		if(physmath::deltaR(jetsAK8->front(), qq) > 0.4)
			return;
		else{
			writeInfoAK8(jetsAK8->front(), outFile8);
			outFile8<<'\n';  // '\n' instead of std::endl because it does not force buffer flush
		}
	}
}


void DBGenerator::writeInfoAK4(const Boson<Jet> bestAK4, std::ofstream& outFile){
	std::sort(AK4pairs_->begin(), AK4pairs_->end(), phys::DeltaRComparator(bestAK4));
	for(size_t i = 0; i < AK4pairs_->size(); ++i){
		Boson<Jet>& jj = AK4pairs_->at(i);
		if(i == 0)
			outFile<<1;
		else
			outFile<<0;
		printVars(1, theWeight);
		
		float dR = physmath::deltaR(jj.daughter(0), jj.daughter(1));
		float pt = jj.pt();
		float pt_scalar = jj.daughter(0).pt() + jj.daughter(1).pt();
		float m = jj.mass();
		float mDM = minDM(m);
		
		printVars(outFile, 5, dR, pt, pt_scalar, m, mDM);
	}
}


void DBGenerator::writeInfoAK8(const Jet bestAK8, std::ofstream& outFile){
	std::sort(jetsAK8->begin(), jetsAK8->end(), phys::DeltaRComparator(bestAK8)); // Just to be sure
	for(size_t i = 0; i < jetsAK8->size(); ++i){
		Jet& j = jetsAK8->at(i);
		if(i == 0)
			outFile<<1;
		else
			outFile<<0;
		printVars(1, theWeight);
		
		float tau21 = j.tau2()/j.tau1();  // This should be high
		float tau32 = j.tau3()/j.tau2();  // This should be low
		//float ptau21 = j.puppiTau2()/j.puppiTau1();
		//float ptau32 = j.puppiTau3()/j.puppiTau2();
		float pt = j.pt();
		float m = j.chosenAlgoMass();
		float mDM = minDM(m);
		
		//printVars(7, tau21, tau32, ptau21, ptau32, pt, m, mDM);
		printVars(outFile, 5, tau21, tau32, pt, m, mDM);
	}
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
		outputFile_ << SEP_CHAR << 0;
	}
}
void DBGenerator::printVars(size_t n, ...){
	double val;
  va_list vl;
  va_start(vl, n);
  for (size_t i = 0; i < n; ++i){
		val = va_arg(vl, double); //For historic reasons floats are promoted to double anyway
		outputFile_ << SEP_CHAR << (float)val;
	}
  va_end(vl);
}

void DBGenerator::printZeroes(std::ofstream& outFile, size_t nzeros){
	for(size_t i = 0; i < nzeros; ++i){
		outFile << SEP_CHAR << 0;
	}
}
void DBGenerator::printVars(std::ofstream& outFile, size_t n, ...){
	double val;
  va_list vl;
  va_start(vl, n);
  for (size_t i = 0; i < n; ++i){
		val = va_arg(vl, double); //For historic reasons floats are promoted to double anyway
		outFile << SEP_CHAR << (float)val;
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
