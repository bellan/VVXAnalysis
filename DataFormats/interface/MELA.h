#ifndef VVXAnalysis_DataFormats_MELA_H
#define VVXAnalysis_DataFormats_MELA_H

/** \class MELA
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

  

  class MELA: public TObject {
    
    friend class ::TreePlanter;

  public:
    
    /// Constructor
    MELA(){
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
    virtual ~MELA(){};
    
    // Operations



    friend std::ostream&  operator<<(std::ostream& os, const MELA& ev){
      
      //os << endl;
      return os;
    }
    
  private:

    mutable Float_t p_JJVBF_BKG_MCFM_JECNominal_;
    mutable Float_t p_JJQCD_BKG_MCFM_JECNominal_;
    mutable Float_t p_JJVBF_BKG_MCFM_JECUp_;     
    mutable Float_t p_JJQCD_BKG_MCFM_JECUp_;     
    mutable Float_t p_JJVBF_BKG_MCFM_JECDn_;     
    mutable Float_t p_JJQCD_BKG_MCFM_JECDn_;     
    mutable Float_t p_JJEW_BKG_MCFM_JECNominal_;
    mutable Float_t p_JJEW_BKG_MCFM_JECUp_;
    mutable Float_t p_JJEW_BKG_MCFM_JECDn_;     
    
    

    ClassDef(MELA, 1) //     
  };

}

#endif

