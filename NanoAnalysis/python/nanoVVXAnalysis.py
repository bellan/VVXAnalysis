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


from VVXAnalysis.NanoAnalysis.VVXEventTaggerAndFilter import VVXEventTaggerAndFilter as VVXTagger

VVX_Regions = [
    {'name'   : 'R4P',        
     'leptons'  : {
         'selection': [
             {
                 'name': 'pt10',
                 'cuts': {'pt': ('>', 10)},
                 'min_particles': 2,
                 'max_particles': float('inf')
             },
             {
                 'name': 'pt20',
                 'cuts': {'pt': ('>', 20)},
                 'min_particles': 1,
                 'max_particles': float('inf')
             },
          ]
     }
     }
]
    
insertAfter(ZZSequence, 'jetFiller', VVXTagger(VVX_Regions))


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

