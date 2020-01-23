#include "VVXAnalysis/TreeAnalysis/interface/DBGenerator.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <iostream>
#include <fstream>
#include <string>
#include "TSystem.h"

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
	return 1;
}

void DBGenerator::begin(){
	cout<<"\n\tfilename = \""<<fileName<<"\"\n"; //fileName derives from EventAnalyzer.h
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
		gSystem->MakeDirectory("databases");//mkdir("Trees",S_IRWXU);
	
	string out = "databases/"+strippedName+".txt";
	cout<<"\toutName = "<<out<<"\n\n";
	outputFile.open(out, std::ios::trunc);
	if(!outputFile.is_open())
		exit(0);
	
	return;
}

void DBGenerator::analyze(){
	//Here we write relevant variables to the outputFile
}

void DBGenerator::end(TFile &){
	outputFile.close();
}
