#ifndef VVXAnalysis_DataFormats_GenEventWeights_H
#define VVXAnalysis_DataFormats_GenEventWeights_H

/** \class GenEventWeights
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

  

  class GenEventWeights: public TObject {
    
    friend class ::TreePlanter;

  public:
    
    /// Constructor
    GenEventWeights(){
      
      mcprocweight_       = 1.;
      puweight_           = 1.; 
      puweightUp_        = 0.; 
      puweightDn_        = 0.; 
      
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
      
      L1PrefiringWeight_   = 1.;
      L1PrefiringWeightUp_ = 1.;
      L1PrefiringWeightDn_ = 1.;
    };
	
    /// Destructor
    virtual ~GenEventWeights(){};
    
    // Operations
    double mcProcWeight()         const {return mcprocweight_;}
    double puWeight()             const {return puweight_;}
    double puWeightUncUp()        const {return puweightUp_/puweight_;}
    double puWeightUncDn()        const {return puweightDn_/puweight_;}

    double kF_ggZZ    () const {return kFactor_ggZZ_    ;} 
    double kF_qqZZM   () const {return kFactor_qqZZM_   ;}
    double kF_qqZZPt  () const {return kFactor_qqZZPt_  ;}
    double kF_qqZZdPhi() const {return kFactor_qqZZdPhi_;}
    double kF_EWKqqZZ () const {return kFactor_EWKqqZZ_ ;}

    float PDFScale           () const {return LHEPDFScale_;}
    float QCDscale_muR1F1    () const {return LHEweight_QCDscale_muR1_muF1_ ;}
    float QCDscale_muR1F2    () const {return LHEweight_QCDscale_muR1_muF2_ ;}
    float QCDscale_muR1F0p5  () const {return LHEweight_QCDscale_muR1_muF0p5_ ;}
    float QCDscale_muR2F1    () const {return LHEweight_QCDscale_muR2_muF1_ ;}
    float QCDscale_muR2F2    () const {return LHEweight_QCDscale_muR2_muF2_ ;}
    float QCDscale_muR2F0p5  () const {return LHEweight_QCDscale_muR2_muF0p5_ ;}
    float QCDscale_muR0p5F1  () const {return LHEweight_QCDscale_muR0p5_muF1_ ;}
    float QCDscale_muR0p5F2  () const {return LHEweight_QCDscale_muR0p5_muF2_ ;}
    float QCDscale_muR0p5F0p5() const {return LHEweight_QCDscale_muR0p5_muF0p5_ ;}
    float PDFVar_Up          () const {return LHEweight_PDFVariation_Up_;}
    float PDFVar_Down        () const {return LHEweight_PDFVariation_Dn_;}
    float alphas_MZ_Up       () const {return LHEweight_AsMZ_Up_;}
    float alphas_MZ_Down     () const {return LHEweight_AsMZ_Dn_;}

    float L1PrefiringWeight  () const {return L1PrefiringWeight_;  }
    float L1PrefiringWeightUp() const {return L1PrefiringWeightUp_;}
    float L1PrefiringWeightDn() const {return L1PrefiringWeightDn_;}



    friend std::ostream&  operator<<(std::ostream& os, const GenEventWeights& ev){
      
      //os << endl;
      return os;
    }
    
  private:

    Double_t mcprocweight_;
    Double_t puweight_;
    Double_t puweightUp_;
    Double_t puweightDn_;
    

    Float_t kFactor_ggZZ_;
    Float_t kFactor_qqZZM_;
    Float_t kFactor_qqZZPt_;
    Float_t kFactor_qqZZdPhi_;
    Float_t kFactor_EWKqqZZ_;
    
    
    // LHE weights
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

    Float_t L1PrefiringWeight_;  
    Float_t L1PrefiringWeightUp_;
    Float_t L1PrefiringWeightDn_;

    ClassDef(GenEventWeights, 1) //     
  };

}

#endif

