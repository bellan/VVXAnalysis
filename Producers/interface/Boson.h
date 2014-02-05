#ifndef Boson_H
#define Boson_H

#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include "Math/GenVector/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

class Boson {

 public:

  Boson(LorentzVector daughter1 = LorentzVector(0., 0., 0., 0.), LorentzVector daughter2  = LorentzVector(0., 0., 0., 0.), 
	int bosonId = -1, int daughtersId = 1):fdaughter1(daughter1),fdaughter2(daughter2),fbosonId(bosonId),fdaughtersId(daughtersId){

  }

  void Setdaughter1(LorentzVector newdaughter1) {fdaughter1 = newdaughter1;}
  void Setdaughter2(LorentzVector newdaughter2) {fdaughter2 = newdaughter2;}
  void SetbosonId(int newbosonId) {fbosonId = newbosonId;}
  void SetdaughtersId(int newdaughtersId) {fdaughtersId = newdaughtersId;}
  LorentzVector p4() const {return fdaughter1+fdaughter2;}
  LorentzVector p4daughter1() const {return fdaughter1;}
  LorentzVector p4daughter2() const {return fdaughter2;}
  int bosonId() const {return fbosonId;}
  int daughtersId() const {return fdaughtersId;}

 private:

  LorentzVector fdaughter1;
  LorentzVector fdaughter2;
  int fbosonId;
  int fdaughtersId;
};

#endif
