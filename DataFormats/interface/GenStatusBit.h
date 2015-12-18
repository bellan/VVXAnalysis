#ifndef VVXAnalysis_DataFormats_GenStatusBit_H
#define VVXAnalysis_DataFormats_GenStatusBit_H

namespace phys{
  enum GenStatusBit {
    

    /////////////////////////////////////////////////////////////////////////////
    //these are robust, generator-independent functions for categorizing
    //mainly final state particles, but also intermediate hadrons/taus    


    //is particle prompt (not from hadron, muon, or tau decay)
     isPrompt,

    
    //is particle a decayed hadron, muon, or tau (does not include resonance decays like W,Z,Higgs,top,etc)
    //This flag is equivalent to status 2 in the current HepMC standard
    //but older generators (pythia6, herwig6) predate this and use status 2 also for other intermediate
    //particles/states    
     isDecayedLeptonHadron,

    
    //this particle is a direct or indirect tau decay product
     isTauDecayProduct,

    
    //this particle is a direct or indirect decay product of a prompt tau
     isPromptTauDecayProduct,

    
    //this particle is a direct tau decay product
     isDirectTauDecayProduct,

    
    //this particle is a direct decay product from a prompt tau 
     isDirectPromptTauDecayProduct,

    
    //this particle is a direct decay product from a hadron
     isDirectHadronDecayProduct,

    
    /////////////////////////////////////////////////////////////////////////////
    //these are generator history-dependent functions for tagging particles
    //associated with the hard process
    //Currently implemented for Pythia 6 and Pythia 8 status codes and history   
    //and may not have 100% consistent meaning across all types of processes
    //Users are strongly encouraged to stick to the more robust flags above    
    
    //this particle is part of the hard process
     isHardProcess,

    
    //this particle is the direct descendant of a hard process particle of the same pdg id
     fromHardProcess,

    
    //this particle is a direct or indirect decay product of a tau
    //from the hard process
     isHardProcessTauDecayProduct,

    
    //this particle is a direct decay product of a tau
    //from the hard process
     isDirectHardProcessTauDecayProduct,

    
    //this particle is the direct descendant of a hard process particle of the same pdg id
    //For outgoing particles the kinematics are those before QCD or QED FSR
    //This corresponds roughly to status code 3 in pythia 6    
     fromHardProcessBeforeFSR,

    
    //this particle is the first copy of the particle in the chain with the same pdg id 
     isFirstCopy,

    
    //this particle is the last copy of the particle in the chain with the same pdg id
    //(and therefore is more likely, but not guaranteed, to carry the final physical momentum)    
     isLastCopy,

    
    //this particle is the last copy of the particle in the chain with the same pdg id
    //before QED or QCD FSR
    //(and therefore is more likely, but not guaranteed, to carry the momentum after ISR)  
     isLastCopyBeforeFSR


 
  };
}
#endif
