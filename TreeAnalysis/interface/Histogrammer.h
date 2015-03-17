#ifndef VVXAnalysis_TreeAnalysis_Histogrammer_H
#define VVXAnalysis_TreeAnalysis_Histogrammer_H

/** \class Histogrammer
 *  Class for easy histogramming
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include <map>

#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TString.h"

class Histogrammer{
  typedef std::map<std::string,TH1*> TH1map;

 public:
 Histogrammer()
   : profile_(false)
   , category_(0){
    TH1F::SetDefaultSumw2(kTRUE);
  }
  
  ~Histogrammer(){}

  // Methods for 1D histograms
  template<typename H>
    TH1* book(const std::string& name, const std::string& title, int bins, const double& min, const double& max){
    TH1map::iterator f = thePlots.find(name);
    if(f != thePlots.end()) return f->second; 
    else{
      thePlots[name] = new H(TString(name),TString(title),bins,min,max);
      return thePlots[name];
    }
  }

 template<typename H>
    TH1* book(const std::string& name, int bins, const double& min, const double& max){
   return book<H>(name, name, bins, min, max);
 }
 
 template<typename H>
   void fill(const std::string& name, const std::string& title, int bins, const double& min, const double& max, const double& value, const double& weight = 1){

   if(profile())
     // Profile by gen category, right now it is a word of 64 bit
     book<TH2F>(name, name, 64, 0, 64, bins, min, max)->Fill(*category_,value,weight);
   else
     book<H>(name, name, bins, min, max)->Fill(value,weight);

 }

 template<typename H>
   void fill(const std::string& name, int bins, const double& min, const double& max, const double& value, const double& weight = 1){
   fill<H>(name, name, bins, min, max, value, weight);
 }

 void fill(const std::string& name, const std::string& title, int bins, const double& min, const double& max, const double& value, const double& weight = 1){
   fill<TH1F>(name, name, bins, min, max, value, weight);
 }
 
 void fill(const std::string& name, int bins, const double& min, const double& max, const double& value, const double& weight = 1){
   fill(name, name, bins, min, max, value, weight);
 }

 // --------------------- Variable binning size ---------------------------

 template<typename H>
   TH1* book(const std::string& name, const std::string& title, const std::vector<double>& xbins){
   TH1map::iterator f = thePlots.find(name);
   if(f != thePlots.end()) return f->second;
    else{
      thePlots[name] = new H(TString(name),TString(title),xbins.size()-1,&xbins[0]);
      return thePlots[name];
    }
 }
 
 template<typename H>
   TH1* book(const std::string& name, const std::vector<double>& xbins){
   return book<H>(name, name, xbins);
 }
 
 template<typename H>
   void fill(const std::string& name, const std::string& title, const std::vector<double>& xbins, const double& value, const double& weight = 1){

   /* if(profile()) */
   /*   // Profile by gen category, right now it is a word of 64 bit */
   /*   book<TH2F>(name, name, 64, 0, 64, bins, xbins)->Fill(*category_,value,weight); */
   /* else */
     book<H>(name, name, xbins)->Fill(value,weight);

 }

 template<typename H>
   void fill(const std::string& name, const std::vector<double>& xbins, const double& value, const double& weight = 1){
   fill<H>(name, name, xbins, value, weight);
 }

 void fill(const std::string& name, const std::string& title, const std::vector<double>& xbins, const double& value, const double& weight = 1){
   fill<TH1F>(name, name, xbins, value, weight);
 }
 
 void fill(const std::string& name, const std::vector<double>& xbins, const double& value, const double& weight = 1){
   fill(name, name, xbins, value, weight);
 }

 // -----------------------------------------------------------------------


 void profile(const int& category) {
   category_ = &category;
   profile_  = true;
 }

 bool profile() const {
   // Do not profile data
   return profile_ && *category_ >=0;
 }

 
 // Methods for 2D histograms

 template<typename H>
   TH2* book(const std::string& name, const std::string& title, 
	     const double& xbins, const double& xmin, const double& xmax,
	     const double& ybins, const double& ymin, const double& ymax){
    TH1map::iterator f = thePlots.find(name);
    if(f != thePlots.end()) return dynamic_cast<TH2*>(f->second); 
    else{
      thePlots[name] = new H(TString(name),TString(title),xbins,xmin,xmax,ybins,ymin,ymax);
      return dynamic_cast<TH2*>(thePlots[name]);
    }
  }

 template<typename H>
   TH1* book(const std::string& name, 
	     const double& xbins, const double& xmin, const double& xmax,
	     const double& ybins, const double& ymin, const double& ymax){
   return book<H>(name, name, xbins, xmin, xmax, ybins, ymin, ymax);
 }
 
 template<typename H>
   void fill(const std::string& name, const std::string& title, 
	     const double& xbins, const double& xmin, const double& xmax, 
	     const double& ybins, const double& ymin, const double& ymax, 
	     const double& xvalue, const double& yvalue, const double& weight = 1){
   book<H>(name, name, xbins, xmin, xmax, ybins, ymin, ymax)->Fill(xvalue,yvalue,weight);
 }

 template<typename H>
   void fill(const std::string& name, 
	     const double& xbins, const double& xmin, const double& xmax, 
	     const double& ybins, const double& ymin, const double& ymax, 
	     const double& xvalue, const double& yvalue, const double& weight = 1){
   fill<H>(name, name, xbins, xmin, xmax, ybins, ymin, ymax, xvalue, yvalue, weight);
 }

 void fill(const std::string& name, const std::string& title, 
	   const double& xbins, const double& xmin, const double& xmax, 
	   const double& ybins, const double& ymin, const double& ymax, 
	   const double& xvalue, const double& yvalue, const double& weight = 1){
   book<TH2F>(name, name, xbins, xmin, xmax, ybins, ymin, ymax)->Fill(xvalue,yvalue,weight);
 }
 
 void fill(const std::string& name, 
	   const double& xbins, const double& xmin, const double& xmax, 
	   const double& ybins, const double& ymin, const double& ymax, 
	   const double& xvalue, const double& yvalue, const double& weight = 1){
   fill(name, name, xbins, xmin, xmax, ybins, ymin, ymax, xvalue, yvalue, weight);
 }




 // Methods for all histogram types
 void write(TFile& fout){
   fout.cd(); 
   for(TH1map::const_iterator h = thePlots.begin(); h != thePlots.end(); ++h)
     h->second->Write();
 }


  
 private:
  // Histograms container
  std::map<std::string,TH1*> thePlots;
  bool profile_;
  const int  *category_;

};



#endif
