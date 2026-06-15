from ZZAnalysis.NanoAnalysis.tools import setConf, getConf, insertBefore, insertAfter

'''Specific configurations'''
setConf("PROCESS_CR"   ,True)
setConf("PROCESS_ZL"   , True)
setConf("FILTER_EVENTS","Z")

from ZZAnalysis.NanoAnalysis.nanoZZ4lAnalysis import *

print("Overriding ZZ4lAnalysis PostProcessor")


item = 'drop Photon*'
if item in branchsel_in: branchsel_in.remove(item)
branchsel_out.append('keep Photon*')
branchsel_out.append('keep regionWord*')


from VVXAnalysis.NanoAnalysis.FSEventTaggerAndFilter import FSEventTaggerAndFilter as FSTagger

#from VVXAnalysis.NanoAnalysis.VVXEventTaggerAndFilter import VVXEventTaggerAndFilter as VVXTagger

from VVXAnalysis.NanoAnalysis.Regions import flagDefinitions
from VVXAnalysis.NanoAnalysis.Regions import Flags

flagsSelection = [
    # Flags.L4P, Flags.L3P1F, Flags.L2P2F,
    # Flags.L3P, Flags.L2P1F, Flags.L1P2F, Flags.L3F,
    # Flags.L2P1L,
    # Flags.L2P_J2
    Flags.L4P_P1P
]
    
insertAfter(ZZSequence, 'jetFiller', FSTagger(flagDefinitions, flagsSelection))


p = PostProcessor(".", fileNames,
                  prefetch=True, longTermCache=False,
                  cut=preselection, # pre-selection cuts (to speed up processing)
                  branchsel=branchsel_in, # select branches to be read
                  outputbranchsel=branchsel_out, # select branches to be written out
                  jsonInput=jsonFile, # path of json file for data
                  modules=ZZSequence,
                  noOut=False, # True = do not write out skimmed nanoAOD file
                  haddFileName="ZZ4lAnalysis.root", # name of output nanoAOD file
#                  histFileName="histos.root", histDirName="plots", # file containing histograms
                  maxEntries=0, # Number of events to be read
                  firstEntry=0, # First event to be read
                  provenance = False
                  )

for cf in customizations :
    print(f"Applying process customization: {cf.__name__}")
    cf(p)

# Print sequence to be run:
print("Sequence to be run:")
for mod in p.modules:
    print(" ", mod.__class__.__name__)
print ("", flush=True)

