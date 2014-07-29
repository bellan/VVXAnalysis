import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring()
)

    

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1



#process.TFileService=cms.Service('TFileService',
#                               fileName=cms.string('VVVAnalysis.root')
#                              )


    
process.source.fileNames = cms.untracked.vstring(
    # '/store/cmst3/group/cmgtools/CMG/WZZ_8TeV-aMCatNLO-herwig/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_13_1_YZF.root'    
    # 'root://lxcms00//data3/2013/HZZ_cmgTuple/BE539_H1258TeV.root' #533 V5_15_0 version
    #'/store/cmst3/group/cmgtools/CMG/WZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_10_1_nLP.root'
    #'/store/cmst3/group/cmgtools/CMG/WZZ_8TeV-aMCatNLO-herwig/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_100_1_GEb.root'
    #'/store/cmst3/group/cmgtools/CMG/ZZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_10_1_UV1.root'
    '/store/cmst3/user/cmgtools/CMG/ZZTo2e2mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_100_1_irQ.root'
    #'/store/cmst3/user/cmgtools/CMG//ZZTo4mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_100_1_UR6.root'
    #'/store/cmst3/user/cmgtools/CMG/VBF_phantom_8TeV/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/cmgTuple_42.root'
    )

# process.options.wantSummary = False

###
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
                                   maxEventsToPrint = cms.untracked.int32(-1),
                                   printVertex = cms.untracked.bool(False),
                                   src = cms.InputTag("genParticlesPruned")
                                    )    

option = int('0',2) 

process.zzGenCategory =  cms.EDFilter("ZZGenFilterCategory",
                                      Topology = cms.int32(option), # -1 means get everything. 1 means only ZZ,and so on...
                                      Option = cms.bool(True), # True for specific selection, False for every good signal (Two ZZ in correct mass range ecc...) 
                                      src = cms.InputTag("genParticlesPruned")
                                      )


process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string('output.root'),
                                  # put this if you have a filter
                                  SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('analysis')
    ),
                                  
)

#process.analysis = cms.Path(process.printTree*process.zzGenCategory)
process.analysis = cms.Path(process.zzGenCategory)

process.out = cms.EndPath(process.output)



