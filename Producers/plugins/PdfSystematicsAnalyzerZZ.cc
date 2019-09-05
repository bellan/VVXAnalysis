////////// Header section /////////////////////////////////////////////
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH2F.h"
#include "ZZAnalysis/AnalysisStep/interface/bitops.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include <iostream>

#include "VVXAnalysis/Producers/interface/FilterController.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH


class PdfSystematicsAnalyzerZZ: public edm::EDAnalyzer {
public:
  
  PdfSystematicsAnalyzerZZ(const edm::ParameterSet& pset);
  virtual ~PdfSystematicsAnalyzerZZ();
  virtual void analyze(const edm::Event & ev, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void endJob() ;
  virtual void endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& setup);
  //virtual void endRun(const edm::Run& run, const edm::EventSetup& setup);

  void book(TFileService * fs);

private:
  
  FilterController filterController_;

  std::vector<std::string> filterNames_;
  Int_t genCategory_;
  std::vector<edm::InputTag> pdfWeightTags_;
  unsigned int originalEvents_;
  unsigned int selectedEvents_;
  Float_t theNumberOfEvents_;
  std::vector<int> pdfStart_;

  std::vector<double> weightedSelectedEvents_;
  std::vector<double> weighted2SelectedEvents_;
  std::vector<double> weightedEvents_;

  std::map<int,std::vector<double>> NewweightedSelectedEvents_;
  std::map<int,std::vector<double>> Newweighted2SelectedEvents_;
  std::map<int,std::vector<double>> NewweightedEvents_;

  // std::vector<double> weightedSelectedEvents2e2mu_;
  // std::vector<double> weighted2SelectedEvents2e2mu_;
  // std::vector<double> weightedEvents2e2mu_;

  edm::InputTag theGenCategoryLabel; 

 //int finstate_;
  TH2F * h_PdfSet1,  * h_PdfSet2, * h_PdfSet3;
  TH2F * h_Pdf2e2muSet1,  * h_Pdf2e2muSet2, * h_Pdf2e2muSet3;
  std::map<int,std::map<int,TH2 *> > _hPlots;
  // int FinStat;
};

////////// Source code ////////////////////////////////////////////////

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

/////////////////////////////////////////////////////////////////////////////////////

void PdfSystematicsAnalyzerZZ::book(TFileService * fs)
{
 

  //  std::vector<double> EmptyVec;

  for(int f=1; f<=4; f++){

    NewweightedSelectedEvents_[f];//=EmptyVec;
    Newweighted2SelectedEvents_[f];//=EmptyVec;
    NewweightedEvents_[f]; //=EmptyVec;
  }

  _hPlots[1][0] =  fs->make<TH2F>("hPdf4muSet1" , "Pdf4muSet1", 55,  0., 55.,8,1.,9.);
  _hPlots[1][1] =  fs->make<TH2F>("hPdf4muSet2" , "Pdf4muSet2", 55,  0., 55.,8,1.,9.);
  _hPlots[1][2] =  fs->make<TH2F>("hPdf4muSet3" , "Pdf4muSet3", 55,  0., 55.,8,1.,9.);

  _hPlots[2][0] =  fs->make<TH2F>("hPdf4eSet1" , "Pdf4eSet1", 55,  0., 55.,8,1.,9.);
  _hPlots[2][1] =  fs->make<TH2F>("hPdf4eSet2" , "Pdf4eSet2", 55,  0., 55.,8,1.,9.);
  _hPlots[2][2] =  fs->make<TH2F>("hPdf4eSet3" , "Pdf4eSet3", 55,  0., 55.,8,1.,9.);

  _hPlots[3][0] =  fs->make<TH2F>("hPdf2e2muSet1" , "Pdf2e2muSet1", 55,  0., 55.,8,1.,9.);
  _hPlots[3][1] =  fs->make<TH2F>("hPdf2e2muSet2" , "Pdf2e2muSet2", 55,  0., 55.,8,1.,9.);
  _hPlots[3][2] =  fs->make<TH2F>("hPdf2e2muSet3" , "Pdf2e2muSet3", 55,  0., 55.,8,1.,9.);

  _hPlots[4][0] =  fs->make<TH2F>("hPdfTotSet1" , "PdfTotSet1", 55,  0., 55.,8,1.,9.);
  _hPlots[4][1] =  fs->make<TH2F>("hPdfTotSet2" , "PdfTotSet2", 55,  0., 55.,8,1.,9.);
  _hPlots[4][2] =  fs->make<TH2F>("hPdfTotSet3" , "PdfTotSet3", 55,  0., 55.,8,1.,9.);

}

/////////////////////////////////////////////////////////////////////////////////////

PdfSystematicsAnalyzerZZ::PdfSystematicsAnalyzerZZ(const edm::ParameterSet& pset) :
  //finstate_ (pset.getParameter<int>("FinalState")),
  filterController_(pset, consumesCollector()),
  filterNames_(pset.getParameter<std::vector<std::string> > ("FilterNames")),
  pdfWeightTags_(pset.getUntrackedParameter<std::vector<edm::InputTag> > ("PdfWeightTags")) { 

  theGenCategoryLabel = pset.getUntrackedParameter<edm::InputTag>("GenCategory" , edm::InputTag("genCategory"));

  theNumberOfEvents_=0;  

  edm::Service<TFileService> fs;
  PdfSystematicsAnalyzerZZ::book(fs.operator->());

  // h_PdfSet1   = fs->make<TH2F>("hPdfSet1" , "PdfSet1", 55,  0., 55.,8,1.,9.);
  // h_PdfSet2   = fs->make<TH2F>("hPdfSet2" , "PdfSet2", 55,  0., 55.,8,1.,9.);
  // h_PdfSet3   = fs->make<TH2F>("hPdfSet3" , "PdfSet3", 55,  0., 55.,8,1.,9.);

}

/////////////////////////////////////////////////////////////////////////////////////
PdfSystematicsAnalyzerZZ::~PdfSystematicsAnalyzerZZ(){}

/////////////////////////////////////////////////////////////////////////////////////
void PdfSystematicsAnalyzerZZ::beginJob(){
      originalEvents_ = 0;
      selectedEvents_ = 0;

      edm::LogVerbatim("PDFAnalysis") << "PDF uncertainties will be determined for the following sets: ";
      for (unsigned int i=0; i<pdfWeightTags_.size(); ++i) {
            edm::LogVerbatim("PDFAnalysis") << "\t" << pdfWeightTags_[i].instance();
            pdfStart_.push_back(-1);
      }
}

/////////////////////////////////////////////////////////////////////////////////////
void PdfSystematicsAnalyzerZZ::endJob(){
  
  if (originalEvents_==0) {
    edm::LogVerbatim("PDFAnalysis") << "NO EVENTS => NO RESULTS";
    return;
  }
  if (selectedEvents_==0) {
    edm::LogVerbatim("PDFAnalysis") << "NO SELECTED EVENTS => NO RESULTS";
    return;
  }
  
  edm::LogVerbatim("PDFAnalysis") << "\n>>>> Begin of PDF weight systematics summary >>>>";
  edm::LogVerbatim("PDFAnalysis") << "Total number of analyzed data: " << originalEvents_ << " [events]";
  double originalAcceptance = double(selectedEvents_)/originalEvents_;
  edm::LogVerbatim("PDFAnalysis") << "Total number of selected data: " << selectedEvents_ << " [events], corresponding to acceptance: [" << originalAcceptance*100 << " +- " << 100*sqrt( originalAcceptance*(1.-originalAcceptance)/originalEvents_) << "] %";
  


  edm::LogVerbatim("PDFAnalysis") << "\n>>>>> PDF UNCERTAINTIES ON RATE >>>>>>";
  for (unsigned int i=0; i<pdfWeightTags_.size(); ++i) {
    bool nnpdfFlag = (pdfWeightTags_[i].instance().substr(0,5)=="NNPDF");
    unsigned int nmembers = weightedSelectedEvents_.size()-pdfStart_[i];
   
    if (i<pdfWeightTags_.size()-1) nmembers = pdfStart_[i+1] - pdfStart_[i];
    unsigned int npairs = (nmembers-1)/2;
    edm::LogVerbatim("PDFAnalysis") << "RATE Results for PDF set " << pdfWeightTags_[i].instance() << " ---->";
    
    
    double events_central = weightedSelectedEvents_[pdfStart_[i]]; 
    std::cout<<"Events before gen filter "<<theNumberOfEvents_<<std::endl<<"Original events after gen filter "<<originalEvents_<<std::endl;
    //       _hPlots[FinStat][i]->Fill(0.,8,originalEvents_);  
    for(int f=1; f<=4; f++){
      std::string hname= _hPlots[f][i]->GetName();
      _hPlots[f][i]->SetName((hname+"_"+pdfWeightTags_[i].instance()).c_str());
      //      _hPlots[f][i]->Fill(0.,8,originalEvents_);  
      _hPlots[f][i]->Fill(0.,8,theNumberOfEvents_);  
      std::cout<<" pdf "<<pdfWeightTags_[i].instance()<<" finstat "<<f<<" orig "<<originalEvents_<<" weighted events "<<NewweightedEvents_[f][pdfStart_[i]]<<std::endl;
      //   std::cout<<"f "<<f<<" "<<NewweightedEvents_[f][pdfStart_[i]]<<std::endl;
      _hPlots[f][i]->Fill(0.,7.,NewweightedEvents_[f][pdfStart_[i]]);
      _hPlots[f][i]->Fill(0.,6.,NewweightedSelectedEvents_[f][pdfStart_[i]]);
      _hPlots[f][i]->Fill(0.,5.,Newweighted2SelectedEvents_[f][pdfStart_[i]]);
      
    }    
        
    edm::LogVerbatim("PDFAnalysis") << "\tEstimate for central PDF member: " << int(events_central) << " [events]";
    double events2_central = weighted2SelectedEvents_[pdfStart_[i]];
    edm::LogVerbatim("PDFAnalysis") << "\ti.e. [" << std::setprecision(4) << 100*(events_central-selectedEvents_)/selectedEvents_ << " +- " <<
      100*sqrt(events2_central-events_central+selectedEvents_*(1-originalAcceptance))/selectedEvents_ 
				    << "] % relative variation with respect to original PDF";
    
    if (npairs>0) {
      edm::LogVerbatim("PDFAnalysis") << "\tNumber of eigenvectors for uncertainty estimation: " << npairs;
      double wplus = 0.;
      double wminus = 0.;
      unsigned int nplus = 0;
      unsigned int nminus = 0;

      for (unsigned int j=0; j<npairs; ++j) {
	double wa = weightedSelectedEvents_[pdfStart_[i]+2*j+1]/events_central-1.;
	double wb = weightedSelectedEvents_[pdfStart_[i]+2*j+2]/events_central-1.; 
	if (nnpdfFlag) {
	  if (wa>0.) {
	    wplus += wa*wa; 
	    nplus++;
	  } else {
	    wminus += wa*wa;
	    nminus++;
	  }
	  if (wb>0.) {
	    wplus += wb*wb; 
	    nplus++;
	  } else {
	    wminus += wb*wb;
	    nminus++;
	  }
	} else {
	  if (wa>wb) {
	    if (wa<0.) wa = 0.;
	    if (wb>0.) wb = 0.;
	    wplus += wa*wa;
	    wminus += wb*wb;
	  } else {
	    if (wb<0.) wb = 0.;
	    if (wa>0.) wa = 0.;
	    wplus += wb*wb;
	    wminus += wa*wa;
	  }
	}
      }
      if (wplus>0) wplus = sqrt(wplus);
      if (wminus>0) wminus = sqrt(wminus);
      if (nnpdfFlag) {
	if (nplus>0) wplus /= sqrt(nplus);
	if (nminus>0) wminus /= sqrt(nminus);
      }
	
	edm::LogVerbatim("PDFAnalysis") << "\tRelative uncertainty with respect to central member: +" << std::setprecision(4) << 100.*wplus << " / -" << std::setprecision(4) << 100.*wminus << " [%]";
      } 
      else {
	edm::LogVerbatim("PDFAnalysis") << "\tNO eigenvectors for uncertainty estimation";
      }
    }
    
  edm::LogVerbatim("PDFAnalysis") << "\n>>>>> PDF UNCERTAINTIES ON ACCEPTANCE >>>>>>";
  for (unsigned int i=0; i<pdfWeightTags_.size(); ++i) {
    bool nnpdfFlag = (pdfWeightTags_[i].instance().substr(0,5)=="NNPDF");
    unsigned int nmembers = weightedEvents_.size()-pdfStart_[i];
    if (i<pdfWeightTags_.size()-1) nmembers = pdfStart_[i+1] - pdfStart_[i];
    unsigned int npairs = (nmembers-1)/2;
    edm::LogVerbatim("PDFAnalysis") << "ACCEPTANCE Results for PDF set " << pdfWeightTags_[i].instance() << " ---->";
    
    double acc_central = 0.;
    double acc2_central = 0.;
    if (weightedEvents_[pdfStart_[i]]>0) {
      acc_central = weightedSelectedEvents_[pdfStart_[i]]/weightedEvents_[pdfStart_[i]]; 
      acc2_central = weighted2SelectedEvents_[pdfStart_[i]]/weightedEvents_[pdfStart_[i]]; 
    }
    double waverage = weightedEvents_[pdfStart_[i]]/originalEvents_;
    edm::LogVerbatim("PDFAnalysis") << "\tEstimate for central PDF member acceptance: [" << acc_central*100 << " +- " << 
      100*sqrt((acc2_central/waverage-acc_central*acc_central)/originalEvents_)
				    << "] %";
    double xi = acc_central-originalAcceptance;
    double deltaxi = (acc2_central-(originalAcceptance+2*xi+xi*xi))/originalEvents_;
    if (deltaxi>0) deltaxi = sqrt(deltaxi); //else deltaxi = 0.;
    edm::LogVerbatim("PDFAnalysis") << "\ti.e. [" << std::setprecision(4) << 100*xi/originalAcceptance << " +- " << std::setprecision(4) << 100*deltaxi/originalAcceptance << "] % relative variation with respect to the original PDF";
    
    if (npairs>0) {
      edm::LogVerbatim("PDFAnalysis") << "\tNumber of eigenvectors for uncertainty estimation: " << npairs;
      double wplus = 0.;
      double wminus = 0.;
      unsigned int nplus = 0;
      unsigned int nminus = 0;
      for (unsigned int j=0; j<npairs; ++j) {

	for(int f=1; f<=4; f++){
	_hPlots[f][i]->Fill(j,1,NewweightedEvents_[f][pdfStart_[i]+2*j+1]);
	_hPlots[f][i]->Fill(j,2,NewweightedSelectedEvents_[f][pdfStart_[i]+2*j+1]);
	_hPlots[f][i]->Fill(j,3,NewweightedEvents_[f][pdfStart_[i]+2*j+2]);	
	_hPlots[f][i]->Fill(j,4,NewweightedSelectedEvents_[f][pdfStart_[i]+2*j+2]);
	}
	
	double wa = 0.;
	if (weightedEvents_[pdfStart_[i]+2*j+1]>0) wa = (weightedSelectedEvents_[pdfStart_[i]+2*j+1]/weightedEvents_[pdfStart_[i]+2*j+1])/acc_central-1.;
	double wb = 0.;
	if (weightedEvents_[pdfStart_[i]+2*j+2]>0) wb = (weightedSelectedEvents_[pdfStart_[i]+2*j+2]/weightedEvents_[pdfStart_[i]+2*j+2])/acc_central-1.;

	if (nnpdfFlag) {
	  if (wa>0.) {
	    wplus += wa*wa; 
	    nplus++;
	  } else {
	    wminus += wa*wa;
	    nminus++;
	  }
	  if (wb>0.) {
	    wplus += wb*wb; 
	    nplus++;
	  } else {
	    wminus += wb*wb;
	    nminus++;
	  }
	} else {
	  if (wa>wb) {
	    if (wa<0.) wa = 0.;
	    if (wb>0.) wb = 0.;
	    wplus += wa*wa;
	    wminus += wb*wb;
	  } else {
	    if (wb<0.) wb = 0.;
	    if (wa>0.) wa = 0.;
	    wplus += wb*wb;
	    wminus += wa*wa;
	  }
	}
      }
      if (wplus>0) wplus = sqrt(wplus);
      if (wminus>0) wminus = sqrt(wminus);
      if (nnpdfFlag) {
	if (nplus>0) wplus /= sqrt(nplus);
	if (nminus>0) wminus /= sqrt(nminus);
      }
      edm::LogVerbatim("PDFAnalysis") << "\tRelative uncertainty with respect to central member: +" << std::setprecision(4) << 100.*wplus << " / -" << std::setprecision(4) << 100.*wminus << " [%]";
    } else {
      edm::LogVerbatim("PDFAnalysis") << "\tNO eigenvectors for uncertainty estimation";
    }
  }
  edm::LogVerbatim("PDFAnalysis") << ">>>> End of PDF weight systematics summary >>>>";
  
}

/////////////////////////////////////////////////////////////////////////////////////
void PdfSystematicsAnalyzerZZ::analyze(const edm::Event & ev, const edm::EventSetup&){
  
  if (!filterController_.passMCFilter(ev)) return;

  edm::Handle<std::vector<double> > weightHandle;
  for (unsigned int i=0; i<pdfWeightTags_.size(); ++i) {
    if (!ev.getByLabel(pdfWeightTags_[i], weightHandle)) {
      if (originalEvents_==0) {
	edm::LogError("PDFAnalysis") << ">>> WARNING: some weights not found!";
	edm::LogError("PDFAnalysis") << ">>> But maybe OK, if you are prefiltering!";
	edm::LogError("PDFAnalysis") << ">>> If things are OK, this warning should disappear after a while!";
      }

      return;     
    }
  }

  //Check final state  
  edm::Handle<int> genCategory;
  ev.getByLabel(theGenCategoryLabel, genCategory);
  genCategory_ = *genCategory;
  //if(genCategory_ == 0) return; // fixme!!
  // std::cout<<"genCategory "<<genCategory_<<std::endl;
  int FinStat=0;
  // RB: FIXME!!!!!!
  if     (test_bit(genCategory_,7) && !test_bit(genCategory_,8)) FinStat = 1; //4mu
  else if(test_bit(genCategory_,8) && !test_bit(genCategory_,7)) FinStat = 2; //4e
  else if(test_bit(genCategory_,7) && test_bit(genCategory_,8)) FinStat = 3; //2e2mu
  else {
    FinStat = 0;
  }


  originalEvents_++;
  
  bool selectedEvent = true;
  edm::Handle<edm::TriggerResults> triggerResults;
  if (!ev.getByLabel(edm::InputTag("TriggerResults"), triggerResults)) {
    edm::LogError("PDFAnalysis") << ">>> TRIGGER collection does not exist !!!";
    
    return;
  }
  
  const edm::TriggerNames & triggerNames = ev.triggerNames(*triggerResults);
  foreach(const std::string &filterName, filterNames_){
    unsigned i = triggerNames.triggerIndex(filterName);
    // Search also if a producer put the result inside the event, insteadof into the trigger results
    edm::Handle<bool> filterResult;
    bool foundIntoTheEvent = ev.getByLabel(filterName, filterResult);
    
    if (i == triggerNames.size() && !foundIntoTheEvent){
      std::cout << "ERROR: PdfSystematicsAnalyzerZZ::isTriggerBit: path does not exist anywhere! " << filterName << std::endl;
      abort();
    }
    
    if (i != triggerNames.size()) selectedEvent &=  triggerResults->accept(i);
    else                          selectedEvent &= *filterResult;
	
  }
  



      if (selectedEvent) selectedEvents_++;

      for (unsigned int i=0; i<pdfWeightTags_.size(); ++i) {
	if (!ev.getByLabel(pdfWeightTags_[i], weightHandle)) {

	  return;
	}
            std::vector<double> weights = (*weightHandle);
            unsigned int nmembers = weights.size();
	    // Set up arrays the first time wieghts are read
            if (pdfStart_[i]<0) {
                  pdfStart_[i] = weightedEvents_.size();
                  for (unsigned int j=0; j<nmembers; ++j) {
		    
		    for(int f=1; f<=4; f++){
		    NewweightedEvents_[f].push_back(0.);
		    NewweightedSelectedEvents_[f].push_back(0.);
		    Newweighted2SelectedEvents_[f].push_back(0.);  
		    }
		    weightedEvents_.push_back(0.);
		    weightedSelectedEvents_.push_back(0.);
		    weighted2SelectedEvents_.push_back(0.);
		  }
            }
            
            for (unsigned int j=0; j<nmembers; ++j) { 
	    
	      NewweightedEvents_[FinStat][pdfStart_[i]+j] += weights[j];
	      NewweightedEvents_[4][pdfStart_[i]+j] += weights[j];
	      weightedEvents_[pdfStart_[i]+j] += weights[j];

 
                  if (selectedEvent) {

                        NewweightedSelectedEvents_[FinStat][pdfStart_[i]+j] += weights[j];
                        Newweighted2SelectedEvents_[FinStat][pdfStart_[i]+j] += weights[j]*weights[j];

			NewweightedSelectedEvents_[4][pdfStart_[i]+j] += weights[j];
			Newweighted2SelectedEvents_[4][pdfStart_[i]+j] += weights[j]*weights[j];

                        weightedSelectedEvents_[pdfStart_[i]+j] += weights[j];
                        weighted2SelectedEvents_[pdfStart_[i]+j] += weights[j]*weights[j];
                  }
            }

            /*
            printf("\n>>>>>>>>> Run %8d Event %d, members %3d PDF set %s : Weights >>>> \n", ev.id().run(), ev.id().event(), nmembers, pdfWeightTags_[i].instance().data());
            for (unsigned int i=0; i<nmembers; i+=5) {
                  for (unsigned int j=0; ((j<5)&&(i+j<nmembers)); ++j) {
                        printf(" %2d: %7.4f", i+j, weights[i+j]);
                  }
                  safe_printf("\n");
            }
            */
      }
      //	    std::cout<<"Return true"<<std::endl;
}

// void PdfSystematicsAnalyzerZZ::endRun(const edm::Run& run, const edm::EventSetup& setup){
//   std::cout<<"end run"<<std::endl;
//}

//bool PdfSystematicsAnalyzerZZ::endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& setup)
void PdfSystematicsAnalyzerZZ::endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& setup)
{
  // std::cout<<" hei "<<std::endl;
  edm::Handle<edm::MergeableCounter> prePathCounter;
  lumi.getByLabel("genEventCounter", prePathCounter);       // Counter of input events in the input pattuple
  theNumberOfEvents_ += prePathCounter->value;
  std::cout<<" pre pat counter "<<prePathCounter->value;
  std::cout<<" numbers of event so far "<<theNumberOfEvents_<<std::endl;
}
DEFINE_FWK_MODULE(PdfSystematicsAnalyzerZZ);
