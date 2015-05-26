#ifndef VVXAnalysis_DataFormats_GenStatusBit_H
#define VVXAnalysis_DataFormats_GenStatusBit_H

namespace phys{
  enum GenStatusBit {
    
    /////////////////////////////////////////////////////////////////////////////
    ///these are robust, generator-independent functions for categorizing
    //mainly final state particles, but also intermediate hadrons/taus

    //is particle prompt (not from hadron, muon, or tau decay) and final state
    promptFinalState,
     
    //is particle prompt (not from hadron, muon, or tau decay) and decayed
    //such as a prompt tau
    promptDecayed,
     
    //this particle is a direct decay product of a prompt tau and is final state
    //(eg an electron or muon from a leptonic decay of a prompt tau)
    directPromptTauDecayProductFinalState,
     

    /////////////////////////////////////////////////////////////////////////////
    //these are generator history-dependent functions for tagging particles
    //associated with the hard process
    //Currently implemented for Pythia 6 and Pythia 8 status codes and history   
    //and may not have 100% consistent meaning across all types of processes
    //Users are strongly encouraged to stick to the more robust flags above,
    //as well as the expanded set available in GenStatusFlags.h
   
    //this particle is part of the hard process
    hardProcess,
   
    //this particle is the final state direct descendant of a hard process particle  
    fromHardProcessFinalState,
     
    //this particle is the decayed direct descendant of a hard process particle
    //such as a tau from the hard process    
    fromHardProcessDecayed,
     
    //this particle is a direct decay product of a hardprocess tau and is final state
    //(eg an electron or muon from a leptonic decay of a tau from the hard process)
    directHardProcessTauDecayProductFinalState,
   
    //this particle is the direct descendant of a hard process particle of the same pdg id.
    //For outgoing particles the kinematics are those before QCD or QED FSR
    //This corresponds roughly to status code 3 in pythia 6
    //This is the most complex and error prone of all the flags and you are strongly encouraged
    //to consider using the others to fill your needs.
    fromHardProcessBeforeFSR,
     
    //provided for convenience.  Use this one if you were using status 3 before and didn't know or care what it exactly meant
    mostlyLikePythia6Status3,
     
    //this particle is the last copy of the particle in the chain  with the same pdg id
    //(and therefore is more likely, but not guaranteed, to carry the final physical momentum)    
    lastCopy,
   
    //this particle is the last copy of the particle in the chain with the same pdg id
    //before QED or QCD FSR
    //(and therefore is more likely, but not guaranteed, to carry the momentum after ISR)  
    lastCopyBeforeFSR
  };
}
#endif
