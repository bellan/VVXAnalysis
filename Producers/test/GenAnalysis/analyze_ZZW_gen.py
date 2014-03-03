import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

####
sample = "QED6_0"

####

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring()
)

    

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.TFileService=cms.Service('TFileService',
                                fileName=cms.string('VVVAnalysis.root')
                                )

if (sample=="4ljj") :
    process.source.fileNames = cms.untracked.vstring('root://lxcms00//data/VVV/madgraph_pythia_4ljj.root')
    process.TFileService.fileName = "4ljj_3_GenAn.root"
elif (sample=="QED6_0") :
    process.source.fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/b/bellan/QGC/4ljj_QED6_run_0.root')
#    process.source.eventsToProcess = cms.untracked.VEventRange("1:9868")
    process.TFileService.fileName = "prova.root"
elif (sample=="QED6_1") :
    process.source.fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/b/bellan/QGC/4ljj_QED6_run_1.root')
    process.TFileService.fileName = "QED6_1_3_cat2.root"
elif (sample=="QED6_2") :
    process.source.fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/b/bellan/QGC/4ljj_QED6_run_2.root')
    process.TFileService.fileName = "QED6_2_3_cat2.root"
elif (sample=="QED6_3") :
    process.source.fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/b/bellan/QGC/4ljj_QED6_run_3.root')
    process.TFileService.fileName = "QED6_3_3_cat2.root"
elif (sample=="QED6_4") :
    process.source.fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/b/bellan/QGC/4ljj_QED6_run_4.root')
    process.TFileService.fileName = "QED6_4_3_cat2.root"
    ####### ZZW events #######
   # process.source.eventsToProcess = cms.untracked.VEventRange("1:1684" , "1:3700", "1:4981", "1:5260", "1:5295", "1:5308", "1:6480", "1:6595", "1:7564", "1:8676", )
##     "1:9660")
  #  process.TFileService.fileName = "QED_0signal.root"
    ####### ZZZ events #######
   ##  process.source.eventsToProcess = cms.untracked.VEventRange("1:2187", "1:2368", "1:2712", "1:4200", "1:7634", "1:8241")
##     process.TFileService.fileName = "4ljj_QED_ZZZ.root"
elif (sample=="W+ZZ") :
    process.source.fileNames = cms.untracked.vstring('root://lxcms00//data/VVV/W+ZZ_4ljj_PxD.root')
    process.TFileService.fileName = "W+ZZ_3_GenAn.root"
elif (sample=="W-ZZ") :
    process.source.fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/b/bellan/QGC/W-ZZ_4ljj_PxD_run_0.root')
    process.TFileService.fileName = "W-ZZ_3_GenAn.root"

 
# process.options.wantSummary = False

###
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
                                   maxEventsToPrint = cms.untracked.int32(100),
                                   printVertex = cms.untracked.bool(False),
                                   src = cms.InputTag("genParticles")
                                   )
process.genAnalyzer = cms.EDAnalyzer("ZZWCombinedGenAnalyzer")

    

## process.genCategory =  cms.EDFilter("GenFilterCategory",
##                                     Category = cms.int32(-1),
##                                     SignalDefinition = cms.int32(3))
## process.myAnalyzer = cms.EDAnalyzer("ZZWGenAnalyzer",
##                                    Category = cms.InputTag("genCategory")
##                                    )

#process.analysis = cms.Path(process.genCategory*process.myAnalyzer)
process.analysis = cms.Path(process.genAnalyzer) #*process.printTree)
#process.analysis = cms.Path(process.printTree*process.GenAnalyzer)

