from ZZAnalysis.NanoAnalysis.tools import setConf, getConf, insertBefore, insertAfter

'''Specific configurations'''
setConf("PROCESS_CR"   ,True)
setConf("PROCESS_ZL"   , True)
setConf("FILTER_EVENTS","Z")
setConf("JES_SPLITTING", False)

#from ZZAnalysis.NanoAnalysis.nanoZZ4lAnalysis import *

print("Overriding ZZ4lAnalysis PostProcessor")




more_branchsel_in =['keep Photon*','keep *MET*']
more_branchsel_out=['keep Photon*','keep PFMET*', 'keep PuppiMET*', 'keep MET*', 'keep regionWord*',
                    'drop Jet_nSVs', 'drop Jet_sv*', 'drop Jet_hfa*', 'drop Jet_PNet*', 'drop Jet_UParTAK4*', 'drop Jet_rawFactor', 'drop HTXS*',
                    'drop Elec*HEEP',
                    'drop Elec*_chg',
                    'drop Elec*yErr',
                    'drop Electron_I*',
                    'drop Electron_PreshowerEnergy',
                    'drop Electron_cutBased',
                    'drop Electron_dr*',
                    'drop Electron_dzErr',
                    'drop Electron_eInvMinusPInv',
                    'drop Electron_ecalEnergyError',
                    'drop Electron_fbrem',
                    'drop Electron_gs*',
                    'drop Electron_isEcalDriven',
                    'drop Electron_jetDF',
                    'drop Electron_jetNDauCharged',
                    'drop Electron_jetPtRelv2',
                    'drop Electron_jetRelIso',
                    'drop Electron_lostHits',
                    'drop Electron_miniPFRelIso_all',
                    'drop Electron_r*',
                    'drop Electron_sc*',
                    'drop Electron_se*',
                    'drop Electron_sm*',
                    'drop Electron_superclusterEta',
                    'drop Electron_svIdx',
                    'drop Electron_tightCharge',
                    'drop Electron_uncorrected_pt',
                    'drop Electron_vidNestedWPBitmap',
                    'drop JetLeadingIdx',
                    'drop JetSubleadingIdx',
                    'drop Jet_*VJet',
                    'drop Jet_btagDeepFlavCvL',
                    'drop Jet_btagPNetC*',
                    'drop Jet_btagUParTAK4C*',
                    'drop Jet_btagUParTAK4Ele',
                    'drop Jet_btagUParTAK4Mu',
                    'drop Jet_btagUParTAK4S*',
                    'drop Jet_btagUParTAK4UDG',
                    'drop Jet_btagUParTAK4p*',
                    'drop Jet_hadronFlavour',
                    'drop Jet_hfcentralEtaStripSize',
                    'drop Jet_hfs*',
                    'drop Jet_muonS*',
                    'drop Jet_ptThreshold',
                    'drop Muon*soId',
                    'drop Muon_I*',
                    'drop Muon_V*',
                    'drop Muon_b*',
                    'drop Muon_dxyErr',
                    'drop Muon_dxybsErr',
                    'drop Muon_dzErr',
                    'drop Muon_h*',
                    'drop Muon_inTimeMuon',
                    'drop Muon_isGlobal',
                    'drop Muon_isStandalone',
                    'drop Muon_isTracker',
                    'drop Muon_jetDF',
                    'drop Muon_jetNDauCharged',
                    'drop Muon_jetRelIso',
                    'drop Muon_looseId',
                    'drop Muon_me*',
                    'drop Muon_mi*',
                    'drop Muon_mvaLowPt',
                    'drop Muon_n*',
                    'drop Muon_pfRelIso03_chg',
                    'drop Muon_pn*',
                    'drop Muon_promptMVA',
                    'drop Muon_ptErr',
                    'drop Muon_sc*',
                    'drop Muon_segmentComp',
                    'drop Muon_sm*',
                    'drop Muon_so*',
                    'drop Muon_svIdx',
                    'drop Muon_t*',
                    'drop Muon_uncorrected_pt',
                    'drop Phot*calo',
                    'drop Phot*idth',
                    'drop Photon_ecalPFClusterIso',
                    'drop Photon_en*',
                    'drop Photon_es*',
                    'drop Photon_haloTaggerMVAVal',
                    'drop Photon_r9',
                    'drop Photon_s4',
                    'drop Photon_se*',
                    'drop Photon_superclusterEta',
                    'drop Photon_t*',
                    'drop Photon_vidNestedWPBitmap',
                    'drop ZZCand_*_mass',
                    'drop nCle*']
                    




#setConf("branchsel_in_ext"  , 'keep Photon*' , append=True)
#setConf("branchsel_out_ext" , 'keep Photon*' , append=True)
#setConf("branchsel_out_ext" , 'keep regionWord*' , append=True) 

setConf("branchsel_in_ext" , more_branchsel_in ) 
setConf("branchsel_out_ext", more_branchsel_out) 

from VVXAnalysis.NanoAnalysis.FSEventTaggerAndFilter import FSEventTaggerAndFilter as FSTagger
#from VVXAnalysis.NanoAnalysis.VVXEventTaggerAndFilter import VVXEventTaggerAndFilter as VVXTagger
from VVXAnalysis.NanoAnalysis.Regions import flagDefinitions
from VVXAnalysis.NanoAnalysis.Regions import Flags

flagsSelection = [
     Flags.L2P_Z1TMass_J2,
     Flags.L3L,
     Flags.L4L
  
    # Flags.L4P, Flags.L3P1F, Flags.L2P2F,
    # Flags.L3P, Flags.L2P1F, Flags.L1P2F, Flags.L3F,
    # Flags.L2P1L,
    # Flags.L4P_P1P
    #Flags.L4P_P1mvaL
]

def customizeFSTagger_(p):
    insertAfter(p.modules, 'jetFiller', FSTagger(flagDefinitions, flagsSelection))

setConf("customizations",customizeFSTagger_,append=True)
    
# p = PostProcessor(".", fileNames,
#                   prefetch=True, longTermCache=False,
#                   cut=preselection, # pre-selection cuts (to speed up processing)
#                   branchsel=branchsel_in, # select branches to be read
#                   outputbranchsel=branchsel_out, # select branches to be written out
#                   jsonInput=jsonFile, # path of json file for data
#                   modules=ZZSequence,
#                   noOut=False, # True = do not write out skimmed nanoAOD file
#                   haddFileName="ZZ4lAnalysis.root", # name of output nanoAOD file
# #                  histFileName="histos.root", histDirName="plots", # file containing histograms
#                   maxEntries=0, # Number of events to be read
#                   firstEntry=0, # First event to be read
#                   provenance = False
#                   )

# for cf in customizations :
#     print(f"Applying process customization: {cf.__name__}")
#     cf(p)

# # Print sequence to be run:
# print("Sequence to be run:")
# for mod in p.modules:
#     print(" ", mod.__class__.__name__)
# print ("", flush=True)

