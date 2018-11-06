#ifndef VVXAnalysis_DataFormats_GenEventInfo_H
#define VVXAnalysis_DataFormats_GenEventInfo_H

/** \class GenEventInfo
 *  No description available.
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.3 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include <TObject.h>
#include <iostream>


class TreePlanter;

namespace phys {

  

  class GenEventInfo: public TObject {
    
    friend class ::TreePlanter;

  public:
    
    /// Constructor
    GenEventInfo(){
      
      mcprocweight_       = 1.;
      puweight_           = 1.; 
      puweightUp_        = 0.; 
      puweightDn_        = 0.; 
      
      genCategory_    = -1;
        
      kFactor_ggZZ_     = 1.; 
      kFactor_qqZZM_    = 1.; 
      kFactor_qqZZPt_   = 1.;
      kFactor_qqZZdPhi_ = 1.;
      kFactor_EWKqqZZ_  = 1.;
      
      
      LHEPDFScale_ = 0;
      LHEweight_QCDscale_muR1_muF1_ = 0;
      LHEweight_QCDscale_muR1_muF2_ = 0;
      LHEweight_QCDscale_muR1_muF0p5_ = 0;
      LHEweight_QCDscale_muR2_muF1_ = 0;
      LHEweight_QCDscale_muR2_muF2_ = 0;
      LHEweight_QCDscale_muR2_muF0p5_ = 0;
      LHEweight_QCDscale_muR0p5_muF1_ = 0;
      LHEweight_QCDscale_muR0p5_muF2_ = 0;
      LHEweight_QCDscale_muR0p5_muF0p5_ = 0;  
      LHEweight_PDFVariation_Up_ = 0;
      LHEweight_PDFVariation_Dn_ = 0;
      LHEweight_AsMZ_Up_ = 0;
      LHEweight_AsMZ_Dn_ = 0;
      
      p_JJVBF_BKG_MCFM_JECNominal_ = 0.;
      p_JJQCD_BKG_MCFM_JECNominal_ = 0.;
      p_JJVBF_BKG_MCFM_JECUp_ = 0.;     
      p_JJQCD_BKG_MCFM_JECUp_ = 0.;     
      p_JJVBF_BKG_MCFM_JECDn_ = 0.;     
      p_JJQCD_BKG_MCFM_JECDn_ = 0.;     
      p_JJEW_BKG_MCFM_JECNominal_ = 0.;  
      p_JJEW_BKG_MCFM_JECUp_ = 0.;     
      p_JJEW_BKG_MCFM_JECDn_ = 0.;     
    };
	
    /// Destructor
    virtual ~GenEventInfo(){};
    
    // Operations



    friend std::ostream&  operator<<(std::ostream& os, const GenEventInfo& ev){
      
      //os << endl;
      return os;
    }
    
  private:

    Double_t mcprocweight_;
    Double_t puweight_;
    Double_t puweightUp_;
    Double_t puweightDn_;
    
    Int_t genCategory_;

    float kFactor_ggZZ_;
    float kFactor_qqZZM_;
    float kFactor_qqZZPt_;
    float kFactor_qqZZdPhi_;
    float kFactor_EWKqqZZ_;
    
    
    
    Float_t LHEPDFScale_;
    Float_t LHEweight_QCDscale_muR1_muF1_ ;
    Float_t LHEweight_QCDscale_muR1_muF2_ ;
    Float_t LHEweight_QCDscale_muR1_muF0p5_ ;
    Float_t LHEweight_QCDscale_muR2_muF1_ ;
    Float_t LHEweight_QCDscale_muR2_muF2_ ;
    Float_t LHEweight_QCDscale_muR2_muF0p5_ ;
    Float_t LHEweight_QCDscale_muR0p5_muF1_ ;
    Float_t LHEweight_QCDscale_muR0p5_muF2_ ;
    Float_t LHEweight_QCDscale_muR0p5_muF0p5_ ;
    Float_t LHEweight_PDFVariation_Up_;
    Float_t LHEweight_PDFVariation_Dn_;
    Float_t LHEweight_AsMZ_Up_;
    Float_t LHEweight_AsMZ_Dn_;
    
    
    mutable Float_t p_JJVBF_BKG_MCFM_JECNominal_;
    mutable Float_t p_JJQCD_BKG_MCFM_JECNominal_;
    mutable Float_t p_JJVBF_BKG_MCFM_JECUp_;     
    mutable Float_t p_JJQCD_BKG_MCFM_JECUp_;     
    mutable Float_t p_JJVBF_BKG_MCFM_JECDn_;     
    mutable Float_t p_JJQCD_BKG_MCFM_JECDn_;     
    mutable Float_t p_JJEW_BKG_MCFM_JECNominal_;
    mutable Float_t p_JJEW_BKG_MCFM_JECUp_;
    mutable Float_t p_JJEW_BKG_MCFM_JECDn_;     
    
    

    ClassDef(GenEventInfo, 1) //     
  };

}

#endif

