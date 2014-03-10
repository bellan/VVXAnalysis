#ifndef Boson_H
#define Boson_H

#include "Math/GenVector/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

class Boson {

 public:

  Boson(LorentzVector daughter1 = LorentzVector(0., 0., 0., 0.), LorentzVector daughter2  = LorentzVector(0., 0., 0., 0.), 
	int bosonId = -1, int daughter1Id = -1, int daughter2Id = -1):fdaughter1(daughter1),fdaughter2(daughter2),fbosonId(bosonId),fdaughter1Id(daughter1Id),fdaughter2Id(daughter2Id){

  }

  void Setdaughter1(LorentzVector newdaughter1) {fdaughter1 = newdaughter1;}
  void Setdaughter2(LorentzVector newdaughter2) {fdaughter2 = newdaughter2;}
  void SetbosonId(int newbosonId) {fbosonId = newbosonId;}
  void Setdaughter1Id(int newdaughter1Id) {fdaughter1Id = newdaughter1Id;} 
  void Setdaughter2Id(int newdaughter2Id) {fdaughter2Id = newdaughter2Id;}
  LorentzVector p4() const {return fdaughter1+fdaughter2;}
  LorentzVector p4daughter1() const {return fdaughter1;}
  LorentzVector p4daughter2() const {return fdaughter2;}
  int bosonId() const {return fbosonId;}
  int daughter1Id() const {return fdaughter1Id;}
  int daughter2Id() const {return fdaughter2Id;}

 private:

  LorentzVector fdaughter1;
  LorentzVector fdaughter2;
  int fbosonId;
  int fdaughter1Id;
  int fdaughter2Id;
};

#endif
