#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"

TFile* openTFile(const char* path);

void doExtraFix(TH1* total, TH1* pass);  // Apparently floating point approximations can make the value in a bin of "pass" slightly greater than the corresponding one in "total". This fix consists in manually setting the value of the "pass" bin to be equal to the one in "total"

void doEfficiency(TFile* fin, const char* nTot, const char* nPas, const char* nName, const char* nTitle, bool extraFix=false);
void doEfficiency(TFile* fin, const char* nTot, const char* nPas, const char* nName, bool extraFix=false);

void doEfficiency2D(TFile* fin, const char* nTot, const char* nPas, const char* nName, const char* nTitle);

TString getAxesLabels(const TH1*);

void verboseEff(const TH1* hTot, const TH1* hPas, double cl=0.95);

#endif
