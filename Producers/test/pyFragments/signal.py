process.genCategory0   =  cms.EDFilter("ZZGenFilterCategory", Topology = cms.int32(SIGNALDEFINITION), ParticleStatus = cms.int32(1), GenJets = cms.InputTag("genJetSel"), src = cms.InputTag("genParticlesPruned"))
process.signalFilters += process.genCategory0
process.postSkimSignalCounter = cms.EDProducer("EventCountProducer")
process.signalFilters += process.postSkimSignalCounter
