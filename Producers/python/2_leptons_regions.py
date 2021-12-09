PT20_10     = "((daughter(0).pt > 20 && daughter(1).pt > 10) || (daughter(0).pt() > 10 && daughter(1).pt > 20))"
ZMASSWINDOW = "abs(mass -91.19) <= 10"
LEPTONSELECTION = "userFloat('isBestZ') && userFloat('Z1Presel') && userFloat('GoodIsoLeptons')"


process.selectedZCand = cms.EDFilter("PATCompositeCandidateSelector",
                                       src = cms.InputTag("ZCand"),
                                       cut = cms.string(LEPTONSELECTION + " && " + PT20_10 + " && " + ZMASSWINDOW)
                                   )


process.pathFor2LeptonsAnalysis = cms.Path(process.selectedZCand)



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
process.SR2PCounter    = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("SR2P","pathFor2LeptonsAnalysis","zzTrigger"))

process.select2leptons1photonRegions = process.select2leptonsRegions.clone()
process.select2leptons1photonRegions.minPhotons = cms.int32(1)
process.select2leptons1photonRegions.maxPhotons = cms.int32(1)

process.SR2P_1L         = cms.Path(process.select2leptons1photonRegions * process.candSR2P * process.candSR2PFilter)
process.SR2P_1LCounter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("SR2P_1L","pathFor2LeptonsAnalysis","zzTrigger"))



# First, we select the tight leptons
# TIGHTMUONPOG = "isGlobalMuon && isPFMuon" + " && " + "globalTrack.normalizedChi2 < 10 && globalTrack.hitPattern.numberOfValidMuonHits > 0 && numberOfMatchedStations > 1" + " && " + "dB < 0.2 && userFloat('dz') < 0.5" + " && " + "innerTrack.hitPattern.numberOfValidPixelHits > 0 && innerTrack.hitPattern.trackerLayersWithMeasurement > 5"


# process.tightMuonsPOG = cms.EDFilter("PATMuonRefSelector",
#                                     src = cms.InputTag("appendPhotons:muons"),
#                                     cut = cms.string(TIGHTMUONPOG + " && userFloat('passCombRelIsoPFFSRCorr')")
                                     #cut = cms.string(GOODMUON + " && userFloat('passCombRelIsoPFFSRCorr')")
#                                )

#TIGHTELEPOG = 





#process.bareZCandFromTighLeptons = cms.EDProducer("PATCandViewShallowCloneCombiner",
#                                                  decay = cms.string('tightMuonsPOG@+ tightMuonsPOG@-'),
#                                                  cut = cms.string('abs(daughter(0).pdgId()) == abs(daughter(1).pdgId())'),
#                                                  checkCharge = cms.bool(True)
#) 


# # then we rank them
# process.ZCandFromLooseL            = process.ZCand.clone()
# process.ZCandFromLooseL.src        = cms.InputTag("bareZCandFromLooseL")
# process.ZCandFromLooseL.bestZAmong = cms.string("mass > 40 && mass < 120")
# process.ZCandFromLooseL.flags      = cms.PSet(GoodLeptons = cms.string(""),
#                                               GoodIsoLeptons = cms.string(""),
#                                               Z1Presel = cms.string(ZMASSWINDOW+PT20_10),
#                                           )


# #  Just check that the trigger thresholds are passed
# process.selectedZCandFromLooseL = cms.EDFilter("PATCompositeCandidateSelector",
#                                                 src = cms.InputTag("ZlCandFromLooseL"),
#                                                 cut = cms.string("(daughter(0).daughter(0).pt > 20 && (daughter(0).daughter(1).pt > 10 || daughter(1).pt > 10))" + 
#                                                                  " || " +
#                                                                  "(daughter(0).daughter(1).pt > 20 && (daughter(0).daughter(0).pt > 10 || daughter(1).pt > 10))" + 
#                                                                  " || " +
#                                                                  "(daughter(1).pt > 20 && (daughter(0).daughter(0).pt > 10 || daughter(0).daughter(1).pt > 10))"
#                                                              ),
#                                                 checkCharge = cms.bool(False)
#                                             )

# # At this point we have all possible 3 leptons combinations. We now build the WZ bare candidate asking 20 GeV of MET
# process.bareZWCand = cms.EDProducer("PATCandViewShallowCloneCombiner",
#                                     decay = cms.string('selectedZlCandFromLooseL slimmedMETs'),
#                                     cut = cms.string("daughter(1).pt > 20"),
#                                     checkCharge = cms.bool(False)
# )



#process.pathFor2LeptonsAnalysis = cms.Path(process.tightMuonsPOG            +  # 
#                                           process.bareZCandFromTighLeptons 
                                           #process.ZCandFromLooseL          +  # 
                                           #process.ZlCandFromLooseL         +  # 
                                           #process.selectedZlCandFromLooseL +  # 
                                           #process.bareZWCand)                 # 
#                                           )

# from VVXAnalysis.Producers.EventFilter_cfg import eventFilter
# process.select3leptonsRegions = eventFilter.clone()
# process.select3leptonsRegions.minTightLeptons = cms.int32(0)
# process.select3leptonsRegions.minLooseLeptons = cms.int32(3)
# process.select3leptonsRegions.maxTightLeptons = cms.int32(3)
# process.select3leptonsRegions.maxLooseLeptons = cms.int32(3)


# PASSZW_1 = "daughter(0).daughter(0).masterClone.userFloat('isGood') && daughter(0).daughter(0).masterClone.userFloat('passCombRelIsoPFFSRCorr')"
# PASSZW_2 = "daughter(0).daughter(1).masterClone.userFloat('isGood') && daughter(0).daughter(1).masterClone.userFloat('passCombRelIsoPFFSRCorr')"
# PASSZW_3 = "daughter(1).daughter(0).masterClone.userFloat('isGood') && daughter(1).daughter(0).masterClone.userFloat('passCombRelIsoPFFSRCorr')"

# SEL111 =       PASSZW_1 + " && " +       PASSZW_2 + " && " +       PASSZW_3
# SEL110 =       PASSZW_1 + " && " +       PASSZW_2 + " && " + "!" + PASSZW_3
# SEL101 =       PASSZW_1 + " && " + "!" + PASSZW_2 + " && " +       PASSZW_3
# SEL011 = "!" + PASSZW_1 + " && " +       PASSZW_2 + " && " +       PASSZW_3
# SEL100 =       PASSZW_1 + " && " + "!" + PASSZW_2 + " && " + "!" + PASSZW_3
# SEL001 = "!" + PASSZW_1 + " && " + "!" + PASSZW_2 + " && " +       PASSZW_3
# SEL010 = "!" + PASSZW_1 + " && " +       PASSZW_2 + " && " + "!" + PASSZW_3
# SEL000 = "!" + PASSZW_1 + " && " + "!" + PASSZW_2 + " && " + "!" + PASSZW_3


# process.candSR3l       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("bareZWCand"), cut = cms.string(SEL111))
# process.candSR3lFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("candSR3l"), minNumber = cms.uint32(1))
# process.sr3l           = cms.Path(process.select3leptonsRegions * process.candSR3l * process.candSR3lFilter)
# process.sr3lCounter    = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("sr3l","pathFor3LeptonsAnalysis","zzTrigger"))

# process.cand110       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("bareZWCand"), cut = cms.string(SEL110))
# process.cand110Filter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand110"), minNumber = cms.uint32(1))
# process.cr110         = cms.Path(process.select3leptonsRegions * process.cand110 * process.cand110Filter)
# process.cr110Counter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("cr110","pathFor3LeptonsAnalysis","zzTrigger"))

# process.cand101       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("bareZWCand"), cut = cms.string(SEL101))
# process.cand101Filter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand101"), minNumber = cms.uint32(1))
# process.cr101         = cms.Path(process.select3leptonsRegions * process.cand101 * process.cand101Filter)
# process.cr101Counter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("cr101","pathFor3LeptonsAnalysis","zzTrigger"))

# process.cand011       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("bareZWCand"), cut = cms.string(SEL011))
# process.cand011Filter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand011"), minNumber = cms.uint32(1))
# process.cr011         = cms.Path(process.select3leptonsRegions * process.cand011 * process.cand011Filter)
# process.cr011Counter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("cr011","pathFor3LeptonsAnalysis","zzTrigger"))

# process.cand100       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("bareZWCand"), cut = cms.string(SEL100))
# process.cand100Filter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand100"), minNumber = cms.uint32(1))
# process.cr100         = cms.Path(process.select3leptonsRegions * process.cand100 * process.cand100Filter)
# process.cr100Counter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("cr100","pathFor3LeptonsAnalysis","zzTrigger"))

# process.cand001       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("bareZWCand"), cut = cms.string(SEL001))
# process.cand001Filter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand001"), minNumber = cms.uint32(1))
# process.cr001         = cms.Path(process.select3leptonsRegions * process.cand001 * process.cand001Filter)
# process.cr001Counter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("cr001","pathFor3LeptonsAnalysis","zzTrigger"))

# process.cand010       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("bareZWCand"), cut = cms.string(SEL010))
# process.cand010Filter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand010"), minNumber = cms.uint32(1))
# process.cr010         = cms.Path(process.select3leptonsRegions * process.cand010 * process.cand010Filter)
# process.cr010Counter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("cr010","pathFor3LeptonsAnalysis","zzTrigger"))

# process.cand000       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("bareZWCand"), cut = cms.string(SEL000))
# process.cand000Filter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand000"), minNumber = cms.uint32(1))
# process.cr000         = cms.Path(process.select3leptonsRegions * process.cand000 * process.cand000Filter)
# process.cr000Counter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("cr000","pathFor3LeptonsAnalysis","zzTrigger"))
