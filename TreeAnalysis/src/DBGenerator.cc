#include "VVXAnalysis/TreeAnalysis/interface/DBGenerator.h"
#include "VVXAnalysis/TreeAnalysis/interface/VZZAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/VZZAnalyzer_impl.h"
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

//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
//#include <numpy/arrayobject.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/assign/std/vector.hpp> 

//#define USE_PYTHON
//#define PRINT_ONLY_CAND

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
	//qq_    = (Boson<Particle>*) ::operator new (sizeof(Boson<Particle>));
	sigVB_ = (Boson<Particle>*) ::operator new (sizeof(Boson<Particle>));
	
	isSigFile_ = false;
	foreach(const std::string& name, signalNames_)
		if(TString(fileName).Contains(name.c_str())){
			isSigFile_ = true;
			break;
		}
	
	// Name parsing	
	cout<<"\n\tfilename = \""<<fileName<<"\"\n"; //fileName comes from EventAnalyzer.h
	// format: samples/<year>/<sample_name>.root
	
	//std::string::size_type nameLength = fileName.length();
	
	std::string::size_type startOfName = fileName.find_last_of('/');
	std::string::size_type endOfName = fileName.find_last_of('.');
	
	std::string::size_type startOfYear = fileName.find_first_of('/');
	std::string::size_type endOfYear = startOfName;  //fileName.find_last_of('/');
	
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
	
	// Output files for jet tagger
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

	#ifdef USE_PYTHON  // Scikit predictor test
	Py_Initialize();
	cout<<"Py interpreter initialized? "<<Py_IsInitialized()<<std::endl;
	
	PyRun_SimpleString(
   "import sys\n"
   "sys.path.append('./python')\n"
	);
	
	helper_module_ = PyImport_ImportModule("VZZhelper"); // import module
	if (!helper_module_){
		cout<<"Error: could not load \"VZZhelper\""<<std::endl;
		Py_Finalize();
		exit(2);
	}
	
	AK4_classifier_ = PyObject_CallMethod(helper_module_, (char*)"load_object", (char*)"s", (char*)"VZZ_AK4_tree.pkl");
	if(AK4_classifier_ == Py_None){
		cout<<"Error: could not load AK4_classifier_."<<std::endl;
		Py_DECREF(helper_module_);
		Py_Finalize();
		exit(3);
	}
	
	AK8_classifier_ = PyObject_CallMethod(helper_module_, (char*)"load_object", (char*)"s", (char*)"VZZ_AK8_tree.pkl");
	if(AK8_classifier_ == Py_None){
		cout<<"Error: could not load AK8_classifier_."<<std::endl;
		Py_DECREF(AK4_classifier_);
		Py_DECREF(helper_module_);
		Py_Finalize();
		exit(3);
	}
	#endif
	
	cout<<"\nAnalyzed:\t      /"<<tree()->GetEntries()<<std::flush;
	//cout<<'\n';
	return;
}


Int_t DBGenerator::cut(){
	++evtN_;
	cout<<"\r\t\t"<<evtN_;
	
	
	//Cleanup of previous event
	if(genHadVBs_)  genHadVBs_->clear(); //Destroys objects but keeps memory allocated
	if(AK4pairs_)    AK4pairs_->clear();
	//if(AK4GenRec_)  AK4GenRec_->clear();
	//if(AllGenVBjj_)AllGenVBjj_->clear();
	if(qq_.p() > 1.)      qq_    = Boson<Particle>();
	if(sigVB_->p() > 1.) *sigVB_ = Particle();
	
	//Preparation for this event
	fillGenHadVBs();  // -->genHadVBs_
	fillRecHadVBs();  // -->AK4pairs_
	//fillAK4GenRec();
	//fillGenVBtoAK4();
	makeQQ(true);
	
	//Hadronic signal
	sigType_ = VZZAnalyzer::isSignal(); //Uses genJets(AK4) and qq_. Modifies sigVB
	
	// REC BASELINE
	if(!ZZ || ZZ->p() < 1.) return -1;
	if(jets->size() < 2. && jetsAK8->size() == 0) return -2;
	
	return 1;
}


void DBGenerator::analyze(){
	++analyzedN_; analyzedW_ += theWeight;
	
	//vector<double> test(15, 0.);
	//getPyPrediction(test, AK8_classifier_);
	
	mainEvtRec(outputFile_); //TEMP
	
	writeTagger(outputFileAK4_, outputFileAK8_);  // contains the endline
}

void DBGenerator::end(TFile &){
	if(outputFile_.is_open())    outputFile_.close();
	if(outputFileAK4_.is_open()) outputFileAK4_.close();
	if(outputFileAK8_.is_open()) outputFileAK8_.close();
	
	if(genHadVBs_) delete genHadVBs_; //Deallocates memory
	if(AK4pairs_)  delete AK4pairs_;
	
	#ifdef USE_PYTHON
	//Python cleanup and shutdown of the interpreter
	Py_DECREF(AK4_classifier_);
	Py_DECREF(AK8_classifier_);
	Py_DECREF(helper_module_);
	Py_Finalize();
	#endif
	
	cout<<"\nPassing cut: "<<Form("%lu (weighted: %.2f)", analyzedN_, analyzedW_)<<'\n';
	
	float elapsedSec = (float)(clock()-startTime_)/CLOCKS_PER_SEC;
	int elapsedSecInt = (int)elapsedSec;
	cout<<"\nElapsed Time: "<<elapsedSec<<" s\t\t("<<elapsedSecInt/60<<"\' "<<elapsedSecInt%60<<"\")\n";
	for(char i=0; i<25; ++i) cout<<'-';
	cout<<" \tEnd of DBGenerator\t";
	for(char i=0; i<25; ++i) cout<<'-';
	cout<<"\n\n";
}


void DBGenerator::mainEvtRec(std::ofstream& outFile){
	// GEN SIGNAL
	bool sigDefinition = sigType_ > 0 && topology.test(0) && topology.test(4);
	
	//RECONSTUCTION: Move it to VZZ, make nonprivate member function
	const Particle* candBestV = nullptr;  // Obtainable looking only at data
	// As sigVB, it only points to the Particle, it doesn't own it
	
	int sigRecType = 0;
	if(AK4pairs_->size() > 0){
		std::sort(AK4pairs_->begin(), AK4pairs_->end(), Mass2Comparator(phys::ZMASS, phys::WMASS));
		if(physmath::minDM(AK4pairs_->front().mass()) < 30.){
			candBestV = &AK4pairs_->front();
			sigRecType = 1;
		}
	}
	if(sigRecType == 0 && jetsAK8->size() > 0){
		std::sort(jetsAK8->begin(), jetsAK8->end(), Mass2Comparator(phys::ZMASS, phys::WMASS));
		if(physmath::minDM(jetsAK8->front().chosenAlgoMass()) < 30.){
			candBestV = &jetsAK8->front();
			sigRecType = 2;
		}
	}
	
	if(!candBestV) return;  // IMPORTANT: JUMP EVENT IF NO VB CAN BE RECONSTRUCTED!
	
	//else return;
	//END RECONSTRUCTION
	
	// IMPORTANT: JUMP EVENT IF NO VB CAN BE RECONSTRUCTED!
	if(!candBestV || !ZZ->isValid() || !ZZ->passFullSelection()) return;
	
	if(isSigFile_ && sigDefinition)
		outFile<<1;  // 1  Is signal at gen level?
	else
		outFile<<0;
	
	outFile<<SEP_CHAR<<theWeight;  // 2
	
	//Here we write relevant variables to the outputFile
	writeMainEvt(outFile, sigRecType, candBestV);  // uses AK(4|8)classifier_
	
	//if(cand4) delete cand4;
	
	outFile<<'\n';  // Let the buffer be flushed automatically
}


void DBGenerator::writeTagger(std::ofstream& outFile4, std::ofstream& outFile8){
	
	if(qq_.p() < 1.){ // || nLep != 4)
		// this is not a signal event, but the jet algorithms may find a false candidate. This is the kind of error the classifier shall train against
		writeInfoAK4(nullptr, outFile4);
		writeInfoAK8(nullptr, outFile8);
		return;
	}
	
	//AK4:
	Boson<Jet>* bestAK4 = nullptr;
	const Jet* rec4_0;
	const Jet* rec4_1;
	size_t i4_0 = 99;
	size_t i4_1 = 99;
	vector<Jet> cpJets; //(*jets);
	
	if(jets->size() < 2)
		goto AK8;
	
	cpJets = vector<Jet>(*jets);
	
	rec4_0 = VZZAnalyzer::closestSing(jets, qq_.daughter(0), i4_0);
	rec4_1 = VZZAnalyzer::closestSing(jets, qq_.daughter(1), i4_1);
	
	if(rec4_0 && rec4_1){
		float dR0 = physmath::deltaR(*rec4_0, qq_.daughter(0));
		float dR1 = physmath::deltaR(*rec4_1, qq_.daughter(1));
		
		if(i4_0 == i4_1){ // If they're matched to the same jet
			if(dR0 < dR1){ // Search for a second, different jet
				cpJets.erase(cpJets.begin() + i4_0);
				rec4_1 = VZZAnalyzer::closestSing(&cpJets, qq_.daughter(1), i4_1);
			} else if(dR1 < dR0){
				cpJets.erase(cpJets.begin() + i4_1);
				rec4_0 = VZZAnalyzer::closestSing(&cpJets, qq_.daughter(1), i4_0);
			}
		}
	}

	if(rec4_0 && rec4_1){		
		bestAK4 = new Boson<Jet>(*rec4_0, *rec4_1);
	}
	writeInfoAK4(bestAK4, outFile4);
	
	AK8:
	if(jetsAK8->size() > 0){
		const Jet* rec8 = closestSing(jetsAK8, qq_);
		//TODO DOUBT a selection on the mass?
		writeInfoAK8(rec8, outFile8);
	}
	
	if(bestAK4) delete bestAK4;
}


void DBGenerator::writeInfoAK4(const Boson<Jet>* bestAK4, std::ofstream& outFile){
	// If it exists, the true good boson is always written.
	if(bestAK4){
		vector<double>* features = getAK4features(*bestAK4);
		outFile<<1<<SEP_CHAR<<theWeight;
		printVars(outFile, *features);
		outFile<<'\n';  // '\n' instead of std::endl because it does not force buffer flush
		delete features;
	}
	
	#ifdef PRINT_ONLY_CAND // If the algorithm found some other pair, that is the kind of error the classifier shall train against
	Boson<Jet>* cand = findBestVFromPair(jets);
	// Cand must exist (least we dereference a nullptr);
	if(cand && (!bestAK4 || *cand != *bestAK4)){ //see the overload of operator==() in "Particle.h"
		vector<double>* features = getAK4features(*cand);
		outFile<<0<<SEP_CHAR<<theWeight;
		printVars(outFile, *features);
		outFile<<'\n';  // '\n' instead of std::endl because it does not force buffer flush
		delete features;
	}
	
	#else  // Include all the other jet pairs
	size_t start = 0;
	if(bestAK4){
		std::sort(AK4pairs_->begin(), AK4pairs_->end(), phys::DeltaRComparator(*bestAK4));
		start = 1;  // The first pair was already written
	}
	float tmpM = 0.;
	for(size_t i = start; i < AK4pairs_->size(); ++i){
		tmpM = AK4pairs_->at(i).mass();
		if(tmpM < 50. || tmpM > 120.) 
			continue;  // We would exclude them anyway
		outFile<<0<<SEP_CHAR<<theWeight;
		vector<double>* features = getAK4features(AK4pairs_->at(i));
		printVars(outFile, *features);
		outFile<<'\n';
		delete features;
	}
	#endif
}


void DBGenerator::writeInfoAK8(const Jet* bestAK8, std::ofstream& outFile){
	// If it exists, the true good boson is always written.
	
	if(bestAK8){
		vector<double>* features = getAK8features(*bestAK8);
		outFile<<1<<SEP_CHAR<<theWeight;
		printVars(outFile, *features);
		outFile<<'\n';  // '\n' instead of std::endl because it does not force buffer flush
		delete features;
	}
	
	#ifdef PRINT_ONLY_CAND // If the algorithm found some other pair, that is the kind of error the classifier shall train against
	const Jet* cand = findBestVFromSing(jetsAK8);
	// Cand must exist (least we dereference a nullptr);
	if(cand && (!bestAK8 || *cand != *bestAK8)){ //see the overload of operator==() in "Particle.h"
		vector<double>* features = getAK8features(*cand);
		outFile<<0<<SEP_CHAR<<theWeight;
		printVars(outFile, *features);
		outFile<<'\n';  // '\n' instead of std::endl because it does not force buffer flush
		delete features;
	}
	
	#else  // Include all the other jet pairs
	
	size_t start = 0;
	if(bestAK8){
		std::sort(jetsAK8->begin(), jetsAK8->end(), phys::DeltaRComparator(*bestAK8));
		start = 1;  // The first jet was already written
	}
	float tmpM = 0.;
	for(size_t i = start; i < jetsAK8->size(); ++i){
		tmpM = jetsAK8->at(i).mass();
		if(tmpM < 50. || tmpM > 120.)
			continue;  // We would exclude them anyway	
		
		outFile<<0<<SEP_CHAR<<theWeight;
		vector<double>* features = getAK8features(jetsAK8->at(i));
		printVars(outFile, *features);
		outFile<<'\n';
		delete features;
	}
	#endif
}


void DBGenerator::writeMainEvt(std::ofstream& outFile, int sigRecType, const Particle* recV){
	//output of scikit
	#ifdef USE_PYTHON
	double pred = -10.;
	if(AK4_classifier_ && sigRecType==1){
		vector<double>* bufFeat = VZZAnalyzer::getAK4features( *((Boson<Jet>*)recV) );
		pred = VZZAnalyzer::getPyPrediction(*bufFeat, AK4_classifier_);
		delete bufFeat;
	}
	else if(AK8_classifier_ && sigRecType==2){
		vector<double>* bufFeat = VZZAnalyzer::getAK8features( *(Jet*)recV );
		pred = VZZAnalyzer::getPyPrediction(*bufFeat, AK8_classifier_);
		delete bufFeat;
	}
	
	printVars(outFile, 2, (double)sigRecType, pred);       // 2
	#endif
	// outFile<<SEP_CHAR<<met->pt();
	
	//recV is from findBestVfrom[...]()
	
	printVars(outFile, 3, ZZ->mass(), ZZ->pt(), ZZ->eta()); // 5
	
	TLorentzVector candp4(recV->p4());
	TLorentzVector Z1p4(ZZ->first().p4());
	TLorentzVector Z2p4(ZZ->second().p4());                 // 8 (+3)
	
	// Angles
	float ang0 = Z1p4.Angle(Z2p4.Vect());
	float ang1 = Z1p4.Angle(candp4.Vect());
	float ang2 = candp4.Angle(Z2p4.Vect());
	float ang3 = candp4.Angle(ZZ->p4().Vect());             // 12 (+7)
	
	// phi
	float dPhi0 = physmath::deltaPhi(ZZ->first(), ZZ->second());
	float dPhi1 = physmath::deltaPhi(ZZ->first(), *recV);
	float dPhi2 = physmath::deltaPhi(*recV, ZZ->second());
	float dPhi3 = physmath::deltaPhi(*recV, *ZZ);           // 16 (+11)
	
	// eta
	float dEta0 = ZZ->first().eta() - ZZ->second().eta();
	float dEta1 = ZZ->first().eta() - recV->eta();
	float dEta2 = recV->eta() - ZZ->second().eta();
	float dEta3 = recV->eta() - ZZ->eta();                  // 20 (+15)
	
	// R
	float dR0 = physmath::deltaR(ZZ->first(), ZZ->second());
	float dR1 = physmath::deltaR(ZZ->first(), *recV);
	float dR2 = physmath::deltaR(*recV, ZZ->second());
	float dR3 = physmath::deltaR(*recV, *ZZ);               // 24 (+19)
	
	// projection - relative pt
	float pt0 = candp4.P() * sin( candp4.Angle(ZZ->p4().Vect()) );
	float pt1 = Z1p4.P() * sin( Z1p4.Angle((candp4 + Z2p4).Vect()) );
	float pt2 = Z2p4.P() * sin( Z2p4.Angle((candp4 + Z1p4).Vect()) );
	
	printVars(outFile, 19, ang0,ang1,ang2,ang3, dPhi0,dPhi1,dPhi2,dPhi3, dEta0,dEta1,dEta2,dEta3, dR0,dR1,dR2,dR3, pt0,pt1,pt2);
	
	printVars(outFile, 2, (candp4+ZZ->p4()).M(), (candp4+ZZ->p4()).Pt());
	
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

void DBGenerator::printVars(std::ofstream& outFile, size_t n, const double* buf){
  for (size_t i = 0; i < n; ++i){
		outFile << SEP_CHAR << (float)(buf[i]);
	}
}

void DBGenerator::printVars(std::ofstream& outFile, const vector<double>& vect){
  for (size_t i = 0; i < vect.size(); ++i){
		outFile << SEP_CHAR << (float)(vect.at(i));
	}
}

