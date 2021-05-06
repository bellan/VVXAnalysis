### If it is MC, run also the signal definition path
# if IsMC:
#     # Empty sequence to attach the signal filter (if specified in the CSV file)
#     process.mcWZSelectionCounter = cms.EDProducer("EventCountProducer") # not really needeed... it is mainly an hack to get the path executed
#     process.WZsignalFilters = cms.Sequence(process.mcWZSelectionCounter) 
#     process.mcWZSelection   = cms.Path(process.WZsignalFilters)
#     MCFILTER = "mcWZSelection"

#     genCategory =  cms.EDFilter("WZGenFilterCategory",
#                                 Topology       = cms.int32(SIGNALDEFINITION), 
#                                 src            = cms.InputTag("genParticlesFromHardProcess"),
#                                 GenJets        = cms.InputTag("selectedGenJets"),
#                                 GenJetsAK8     = cms.InputTag("selectedGenJetsAK8"),
#                                 )
#     process.WZgenCategory = genCategory

#     process.WZsignalCounter    = cms.EDProducer("EventCountProducer")
#     process.WZsignalDefinition = cms.Path(process.WZgenCategory * process.WZsignalCounter)

process.ZlForWZ = cms.EDFilter("PATCompositeCandidateSelector",
                              src = cms.InputTag("ZlSelected"),
                              cut = cms.string("daughter(1).masterClone.userFloat('isGood')" + " && " +
                                               "daughter(1).masterClone.userFloat('passCombRelIsoPFFSRCorr')"
                                              ),
                              checkCharge = cms.bool(False)
)

process.ZlForWZCountFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("ZlForWZ"), 
                                          minNumber = cms.uint32(1) # maxNumber does not exist
                                          )


process.ZWCand = cms.EDProducer("PATCandViewShallowCloneCombiner",
                                decay = cms.string('ZlForWZ slimmedMETs'),
                                cut = cms.string("daughter(0).daughter(1).pt > 20" + " && " + "daughter(1).pt > 20" 
                                ),
                                checkCharge = cms.bool(False)
)


process.WZjjPath = cms.Path(process.ZlForWZ * process.ZlForWZCountFilter * process.ZWCand)
