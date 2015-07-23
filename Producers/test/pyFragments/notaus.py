process.genTaus             = cms.EDFilter("PdgIdAndStatusCandViewSelector", src = cms.InputTag("genParticlesPruned"), pdgId = cms.vint32( 15 ), status = cms.vint32( 3 ))
process.genTauCounterFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("genTaus"), minNumber = cms.uint32(1))
process.preSkimCounter      = cms.EDProducer("EventCountProducer")
process.signalFilters += process.genTaus
process.signalFilters += ~process.genTauCounterFilter
process.signalFilters += process.preSkimCounter
