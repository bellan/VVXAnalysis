### If it is MC, run also the signal definition path
if IsMC:
    # Empty sequence to attach the signal filter (if specified in the CSV file)
    process.mcWZSelectionCounter = cms.EDProducer("EventCountProducer") # not really needeed... it is mainly an hack to get the path executed
    process.WZsignalFilters = cms.Sequence(process.mcWZSelectionCounter) 
    process.mcWZSelection   = cms.Path(process.WZsignalFilters)
    MCFILTER = "mcWZSelection"

    genCategory =  cms.EDFilter("WZGenFilterCategory",
                                Topology       = cms.int32(SIGNALDEFINITION), 
                                src            = cms.InputTag("genParticlesFromHardProcess"),
                                GenJets        = cms.InputTag("selectedGenJets"),
                                GenJetsAK8     = cms.InputTag("selectedGenJetsAK8"),
                                )
    process.WZgenCategory = genCategory

    process.WZsignalCounter    = cms.EDProducer("EventCountProducer")
    process.WZsignalDefinition = cms.Path(process.WZgenCategory * process.WZsignalCounter)

