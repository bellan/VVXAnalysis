#ifndef VVXAnalysis_DataFormats_RegionsCounter_H
#define VVXAnalysis_DataFormats_RegionsCounter_H

/** \class RegionsCounter
 *  No description available.
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.3 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include <TObject.h>
#include <iostream>
#include <map>

#include "VVXAnalysis/DataFormats/interface/RegionTypes.h"
#include "VVXAnalysis/Commons/interface/Colours.h"

class TreePlanter;

namespace phys {
  class RegionsCounter: public TObject {
    
    friend class ::TreePlanter;

  public:
    
    /// Constructor
    RegionsCounter(){

      eventsInRegions_[phys::SR2P]      = 0;   
      eventsInRegions_[phys::SR2P_1L]   = 0; 
      
      eventsInRegions_[phys::SR3P]      = 0;   
      eventsInRegions_[phys::CR110]     = 0;  
      eventsInRegions_[phys::CR101]     = 0;  
      eventsInRegions_[phys::CR011]     = 0;  
      eventsInRegions_[phys::CR100]     = 0;  
      eventsInRegions_[phys::CR001]     = 0;  
      eventsInRegions_[phys::CR010]     = 0;  
      eventsInRegions_[phys::CR000]     = 0;  
      eventsInRegions_[phys::SR3P_1L]   = 0; 
      eventsInRegions_[phys::CRLFR]     = 0;  
      
      eventsInRegions_[phys::SR4P]      = 0;   
      eventsInRegions_[phys::CR2P2F]    = 0; 
      eventsInRegions_[phys::CR3P1F]    = 0; 
      eventsInRegions_[phys::SR4P_1L]   = 0; 
      eventsInRegions_[phys::CR2P2F_HZZ] = 0;
      eventsInRegions_[phys::CR3P1F_HZZ] = 0;
      eventsInRegions_[phys::SR_HZZ]     = 0;

      blinded_ = true;
    };        

	
    /// Destructor
    virtual ~RegionsCounter(){};
    
    // Getter functions
    bool blinded()   const {return blinded_;}
    bool unblinded() const {return !blinded_;}
    //
    void blind()   {blinded_ = true;}
    void unblind() {blinded_ = false;}

    
    Int_t& operator[](phys::RegionTypes rt)             { return eventsInRegions_[rt]; }
    const Int_t& operator[](phys::RegionTypes rt) const { return eventsInRegions_.at(rt); }

    RegionsCounter& operator+=(const RegionsCounter& rc){


      eventsInRegions_[phys::SR2P]      += rc[phys::SR2P]      ;   
      eventsInRegions_[phys::SR2P_1L]   += rc[phys::SR2P_1L]   ; 
      					  					
      eventsInRegions_[phys::SR3P]      += rc[phys::SR3P]      ;   
      eventsInRegions_[phys::CR110]     += rc[phys::CR110]     ;  
      eventsInRegions_[phys::CR101]     += rc[phys::CR101]     ;  
      eventsInRegions_[phys::CR011]     += rc[phys::CR011]     ;  
      eventsInRegions_[phys::CR100]     += rc[phys::CR100]     ;  
      eventsInRegions_[phys::CR001]     += rc[phys::CR001]     ;  
      eventsInRegions_[phys::CR010]     += rc[phys::CR010]     ;  
      eventsInRegions_[phys::CR000]     += rc[phys::CR000]     ;  
      eventsInRegions_[phys::SR3P_1L]   += rc[phys::SR3P_1L]   ; 
      eventsInRegions_[phys::CRLFR]     += rc[phys::CRLFR]     ;  
      					  					
      eventsInRegions_[phys::SR4P]      += rc[phys::SR4P]      ;   
      eventsInRegions_[phys::CR2P2F]    += rc[phys::CR2P2F]    ; 
      eventsInRegions_[phys::CR3P1F]    += rc[phys::CR3P1F]    ; 
      eventsInRegions_[phys::SR4P_1L]   += rc[phys::SR4P_1L]   ; 
      eventsInRegions_[phys::CR2P2F_HZZ]+= rc[phys::CR2P2F_HZZ];
      eventsInRegions_[phys::CR3P1F_HZZ]+= rc[phys::CR3P1F_HZZ];
      eventsInRegions_[phys::SR_HZZ]    += rc[phys::SR_HZZ]    ;

      blinded_ = blinded_ || rc.blinded();
      
      return *this; // return the result by reference
    }


  RegionsCounter& operator=(const RegionsCounter& rc){


      eventsInRegions_[phys::SR2P]      = rc[phys::SR2P]      ;   
      eventsInRegions_[phys::SR2P_1L]   = rc[phys::SR2P_1L]   ; 
      					  					
      eventsInRegions_[phys::SR3P]      = rc[phys::SR3P]      ;   
      eventsInRegions_[phys::CR110]     = rc[phys::CR110]     ;  
      eventsInRegions_[phys::CR101]     = rc[phys::CR101]     ;  
      eventsInRegions_[phys::CR011]     = rc[phys::CR011]     ;  
      eventsInRegions_[phys::CR100]     = rc[phys::CR100]     ;  
      eventsInRegions_[phys::CR001]     = rc[phys::CR001]     ;  
      eventsInRegions_[phys::CR010]     = rc[phys::CR010]     ;  
      eventsInRegions_[phys::CR000]     = rc[phys::CR000]     ;  
      eventsInRegions_[phys::SR3P_1L]   = rc[phys::SR3P_1L]   ; 
      eventsInRegions_[phys::CRLFR]     = rc[phys::CRLFR]     ;  
      					  					
      eventsInRegions_[phys::SR4P]      = rc[phys::SR4P]      ;   
      eventsInRegions_[phys::CR2P2F]    = rc[phys::CR2P2F]    ; 
      eventsInRegions_[phys::CR3P1F]    = rc[phys::CR3P1F]    ; 
      eventsInRegions_[phys::SR4P_1L]   = rc[phys::SR4P_1L]   ; 
      eventsInRegions_[phys::CR2P2F_HZZ]= rc[phys::CR2P2F_HZZ];
      eventsInRegions_[phys::CR3P1F_HZZ]= rc[phys::CR3P1F_HZZ];
      eventsInRegions_[phys::SR_HZZ]    = rc[phys::SR_HZZ]    ;

      return *this; // return the result by reference
    }





  friend std::ostream&  operator<<(std::ostream& os, const RegionsCounter& ev){      
    
    using namespace colour;
    os << std::endl
                           << "| Region     | Events   |"                         << std::endl
                           << "| ---------- | -------- |"                         << std::endl;
    if(ev.unblinded()){ os << "| SR2P       | " << Green(ev[phys::SR2P]      ) << "\t|" << std::endl
			   << "| SR2P_1L    | " << Green(ev[phys::SR2P_1L]   ) << "\t|" << std::endl
      
			   << "| SR3P       | " << Green(ev[phys::SR3P]      ) << "\t|" << std::endl;}
                        os << "| CR110      | " << Green(ev[phys::CR110]     ) << "\t|" << std::endl
			   << "| CR101      | " << Green(ev[phys::CR101]     ) << "\t|" << std::endl
			   << "| CR011      | " << Green(ev[phys::CR011]     ) << "\t|" << std::endl
			   << "| CR100      | " << Green(ev[phys::CR100]     ) << "\t|" << std::endl
			   << "| CR001      | " << Green(ev[phys::CR001]     ) << "\t|" << std::endl
			   << "| CR010      | " << Green(ev[phys::CR010]     ) << "\t|" << std::endl
			   << "| CR000      | " << Green(ev[phys::CR000]     ) << "\t|" << std::endl;
    if(ev.unblinded())  os << "| SR3P_1L    | " << Green(ev[phys::SR3P_1L]   ) << "\t|" << std::endl;
                        os << "| CRLFR      | " << Green(ev[phys::CRLFR]     ) << "\t|" << std::endl;
      
    if(ev.unblinded())  os << "| SR4P       | " << Green(ev[phys::SR4P]      ) << "\t|" << std::endl;
                        os << "| CR2P2F     | " << Green(ev[phys::CR2P2F]    ) << "\t|" << std::endl
			   << "| CR3P1F     | " << Green(ev[phys::CR3P1F]    ) << "\t|" << std::endl;
    if(ev.unblinded())  os << "| SR4P_1L    | " << Green(ev[phys::SR4P_1L]   ) << "\t|" << std::endl;
                        os << "| CR2P2F_HZZ | " << Green(ev[phys::CR2P2F_HZZ]) << "\t|" << std::endl
                           << "| CR3P1F_HZZ | " << Green(ev[phys::CR3P1F_HZZ]) << "\t|" << std::endl;
    if(ev.unblinded())  os << "| SR_HZZ     | " << Green(ev[phys::SR_HZZ]    ) << "\t|" << std::endl;							
    
    return os;
  }
    
  private:
      
    std::map<phys::RegionTypes,Int_t> eventsInRegions_;
    bool blinded_;

    ClassDef(RegionsCounter, 2) //     
  };

}

#endif

