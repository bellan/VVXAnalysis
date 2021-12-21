PT20_10     = "((daughter(0).pt > 20 && daughter(1).pt > 10) || (daughter(0).pt() > 10 && daughter(1).pt > 20))"
ZMASSWINDOW = "abs(mass -91.19) <= 10"
LEPTONSELECTION = "userFloat('isBestZ') && userFloat('Z1Presel') && userFloat('GoodIsoLeptons')"


process.selectedZCand = cms.EDFilter("PATCompositeCandidateSelector",
                                       src = cms.InputTag("ZCand"),
                                       cut = cms.string(LEPTONSELECTION + " && " + PT20_10 + " && " + ZMASSWINDOW)
                                   )



# For lepton-jet cleaning
process.muonsFromZV     = cms.EDProducer("PATMuonsFromCompositeCandidates"    , src =  cms.InputTag("selectedZCand"), SplitLevel = cms.int32(0))
process.electronsFromZV = cms.EDProducer("PATElectronsFromCompositeCandidates", src =  cms.InputTag("selectedZCand"), SplitLevel = cms.int32(0))
process.leptonsFromZV   = cms.Sequence(process.muonsFromZV + process.electronsFromZV)

process.pathFor2LeptonsAnalysis = cms.Path(process.selectedZCand + process.leptonsFromZV)



from VVXAnalysis.Producers.EventFilter_cfg import eventFilter
process.select2leptonsRegions = eventFilter.clone()
process.select2leptonsRegions.minTightLeptons = cms.int32(2)
process.select2leptonsRegions.minLooseLeptons = cms.int32(2)
process.select2leptonsRegions.maxTightLeptons = cms.int32(2)
process.select2leptonsRegions.maxLooseLeptons = cms.int32(2)
process.select2leptonsRegions.minAK4s = cms.int32(2)
process.select2leptonsRegions.minAK8s = cms.int32(1)

process.candSR2P       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("selectedZCand"), cut = cms.string(''))
process.candSR2PFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("candSR2P"), minNumber = cms.uint32(1))
process.SR2P           = cms.Path(process.select2leptonsRegions * process.candSR2P * process.candSR2PFilter)
process.SR2PCounter    = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("SR2P","pathFor2LeptonsAnalysis","triggerForZV"))

process.select2leptons1photonRegions = process.select2leptonsRegions.clone()
process.select2leptons1photonRegions.minPhotons = cms.int32(1)
process.select2leptons1photonRegions.maxPhotons = cms.int32(1)

process.SR2P_1L         = cms.Path(process.select2leptons1photonRegions * process.candSR2P * process.candSR2PFilter)
process.SR2P1LCounter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("SR2P_1L","pathFor2LeptonsAnalysis","triggerForZV"))





