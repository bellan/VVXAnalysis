#ifndef VVXAnalysis_Commons_RegionTypes_h
#define VVXAnalysis_Commons_RegionTypes_h

#include <string>

namespace phys{
  //enum RegionTypes {SR, CR, CR2P2F, CR3P1F, SR_HZZ, CR_HZZ, CR2P2F_HZZ, CR3P1F_HZZ, MC, MC_HZZ, SR3L};

  // * indicate a SR
  enum RegionTypes{SR4P,    // * 4 tight (pass ID) leptons (baseline SR for ZZ)
		   CR3P1F,  //   3 tight leptons + 1 loose but not tight lepton
		   CR2P2F,  //   2 tight leptons + 2 loose but not tight leptons
		   SR4P_1L, // * 4 tight leptons + 1 loose photon (baseline SR for ZZgamma)
		   SR4P_1P, // * 4 tight leptons + 1 tight photon (SR for ZZgamma)
		   CR4P_1F, //   4 tight leptons + 1 loose but not tight photon
		   
		   SR3P,    // * 3 tight (pass ID) leptons (baseline SR for WZ) 
		   CR110,   //   2 tight leptons + 1 loose but not tight lepton
		   CR101,   //   2 tight leptons + 1 loose but not tight lepton
		   CR011,   //   2 tight leptons + 1 loose but not tight lepton
		   CR100,   //   1 tight leptons + 2 loose but not tight lepton
		   CR001,   //   1 tight leptons + 2 loose but not tight lepton
		   CR010,   //   1 tight leptons + 2 loose but not tight lepton
		   CR000,   //   3 loose but not tight lepton
		   SR3P_1L, // * 3 tight leptons + 1 loose photon (baseline SR for WZgamma)
		   SR3P_1P, // * 3 tight leptons + 1 tight photon (SR for WZgamma)
		   CR3P_1F, //   3 tight leptons + 1 loose but not tight photon
		   CRLFR  , //   2 tight leptons + 1 loose lepton (CR for lepton fake rate measurement)
 
		   // Beware: 2P regions require the presence of either 2 AK4 jets or a AK8 jet
		   SR2P,    // * 2 tight (pass ID) leptons (baseline SR for ZV) 
		   SR2P_1L, // * 2 tight leptons + 1 loose photon (baseline SR for ZVgamma)
		   SR2P_1P, // * 2 tight leptons + 1 tight photon (SR for ZVgamma)
		   CR2P_1F, //   2 tight leptons + 1 loose but not tight photon

		   CR,
		   SR_HZZ, CR_HZZ, CR2P2F_HZZ, CR3P1F_HZZ, MC_HZZ,
		   MC};
  
  RegionTypes regionType(const std::string& input);
  std::string regionType(RegionTypes input);
  
  
  enum Channel {ZZ, ZW, ZL, ZV, UNDEF};
  Channel channelType(const std::string& input);
  std::string channelType(Channel input);
  

};

#endif
