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





### We need to build up several CRs. For ref see CMS AN-2017/135 pag. 15-16


# First, we create al possible SFOS, both e and mu, pairs
process.bareZCandFromLooseL     = process.bareZCand.clone()
process.bareZCandFromLooseL.cut = cms.string('mass > 0 && abs(daughter(0).pdgId()) == abs(daughter(1).pdgId())')

# then we rank them
process.ZCandFromLooseL            = process.ZCand.clone()
process.ZCandFromLooseL.src        = cms.InputTag("bareZCandFromLooseL")
process.ZCandFromLooseL.bestZAmong = cms.string("mass > 40 && mass < 120")
process.ZCandFromLooseL.flags      = cms.PSet(GoodLeptons = cms.string(""),
                                              GoodIsoLeptons = cms.string(""),
                                              Z1Presel = cms.string("mass > 60 && mass < 120"),
                                          )

# combine the Z candidate with a free lepton.
process.ZlCandFromLooseL       = process.ZlCand.clone()
process.ZlCandFromLooseL.decay = cms.string('ZCandFromLooseL softLeptons')


#  Just check that the trigger thresholds are passed
process.selectedZlCandFromLooseL = cms.EDFilter("PATCompositeCandidateSelector",
                                                src = cms.InputTag("ZlCandFromLooseL"),
                                                cut = cms.string("((daughter(0).daughter(0).pt > 20 && (daughter(0).daughter(1).pt > 10 || daughter(1).pt > 10))" + 
                                                                 " || " +
                                                                 "(daughter(0).daughter(1).pt > 20 && (daughter(0).daughter(0).pt > 10 || daughter(1).pt > 10))" + 
                                                                 " || " +
                                                                 "(daughter(1).pt > 20 && (daughter(0).daughter(0).pt > 10 || daughter(0).daughter(1).pt > 10)))" +
                                                                 " && " + "daughter(0).masterClone.userFloat('isBestZ')" +
                                                                 " && " + "daughter(0).masterClone.userFloat('Z1Presel')"
                                                             ),
                                                checkCharge = cms.bool(False)
                                            )

# At this point we have all possible 3 leptons combinations. We now build the WZ bare candidate asking 20 GeV of MET
process.bareZWCand = cms.EDProducer("PATCandViewShallowCloneCombiner",
                                    decay = cms.string('selectedZlCandFromLooseL slimmedMETs'),
                                    cut = cms.string("daughter(1).pt > 20"),
                                    checkCharge = cms.bool(False)
)



process.pathFor3LeptonsAnalysis = cms.Path(process.bareZCandFromLooseL      +  # Z from loose leptons
                                           process.ZCandFromLooseL          +  # best Z from all loose leptons
                                           process.ZlCandFromLooseL         +  # best Z + a free loose lepton
                                           process.selectedZlCandFromLooseL +  # best Z + a free loose lepton w/ trigger requirements
                                           process.bareZWCand)                 # Zl + MET (WZ bare candidate)


from VVXAnalysis.Producers.EventFilter_cfg import eventFilter
process.select3leptonsRegions = eventFilter.clone()
process.select3leptonsRegions.minTightLeptons = cms.int32(0)
process.select3leptonsRegions.minLooseLeptons = cms.int32(3)
process.select3leptonsRegions.maxTightLeptons = cms.int32(3)
process.select3leptonsRegions.maxLooseLeptons = cms.int32(3)
process.select3leptonsRegions.minPhotons = cms.int32(0)
process.select3leptonsRegions.maxPhotons = cms.int32(1000)


PASSZW_1 = "daughter(0).daughter(0).daughter(0).masterClone.userFloat('isGood') && daughter(0).daughter(0).daughter(0).masterClone.userFloat('passCombRelIsoPFFSRCorr')"
PASSZW_2 = "daughter(0).daughter(0).daughter(1).masterClone.userFloat('isGood') && daughter(0).daughter(0).daughter(1).masterClone.userFloat('passCombRelIsoPFFSRCorr')"
PASSZW_3 = "daughter(0).daughter(1).masterClone.userFloat('isGood')             && daughter(0).daughter(1).masterClone.userFloat('passCombRelIsoPFFSRCorr')"

SEL111 =       PASSZW_1 + " && " +       PASSZW_2 + " && " +       PASSZW_3
SEL110 =       PASSZW_1 + " && " +       PASSZW_2 + " && " + "!" + PASSZW_3
SEL101 =       PASSZW_1 + " && " + "!" + PASSZW_2 + " && " +       PASSZW_3
SEL011 = "!" + PASSZW_1 + " && " +       PASSZW_2 + " && " +       PASSZW_3
SEL100 =       PASSZW_1 + " && " + "!" + PASSZW_2 + " && " + "!" + PASSZW_3
SEL001 = "!" + PASSZW_1 + " && " + "!" + PASSZW_2 + " && " +       PASSZW_3
SEL010 = "!" + PASSZW_1 + " && " +       PASSZW_2 + " && " + "!" + PASSZW_3
SEL000 = "!" + PASSZW_1 + " && " + "!" + PASSZW_2 + " && " + "!" + PASSZW_3


process.candSR3P       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("bareZWCand"), cut = cms.string(SEL111))
process.candSR3PFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("candSR3P"), minNumber = cms.uint32(1))
process.SR3P           = cms.Path(process.select3leptonsRegions * process.candSR3P * process.candSR3PFilter)
process.SR3PCounter    = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("SR3P","pathFor3LeptonsAnalysis","zzTrigger"))

process.cand110       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("bareZWCand"), cut = cms.string(SEL110))
process.cand110Filter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand110"), minNumber = cms.uint32(1))
process.CR110         = cms.Path(process.select3leptonsRegions * process.cand110 * process.cand110Filter)
process.CR110Counter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("CR110","pathFor3LeptonsAnalysis","zzTrigger"))

process.cand101       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("bareZWCand"), cut = cms.string(SEL101))
process.cand101Filter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand101"), minNumber = cms.uint32(1))
process.CR101         = cms.Path(process.select3leptonsRegions * process.cand101 * process.cand101Filter)
process.CR101Counter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("CR101","pathFor3LeptonsAnalysis","zzTrigger"))

process.cand011       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("bareZWCand"), cut = cms.string(SEL011))
process.cand011Filter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand011"), minNumber = cms.uint32(1))
process.CR011         = cms.Path(process.select3leptonsRegions * process.cand011 * process.cand011Filter)
process.CR011Counter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("CR011","pathFor3LeptonsAnalysis","zzTrigger"))

process.cand100       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("bareZWCand"), cut = cms.string(SEL100))
process.cand100Filter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand100"), minNumber = cms.uint32(1))
process.CR100         = cms.Path(process.select3leptonsRegions * process.cand100 * process.cand100Filter)
process.CR100Counter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("CR100","pathFor3LeptonsAnalysis","zzTrigger"))

process.cand001       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("bareZWCand"), cut = cms.string(SEL001))
process.cand001Filter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand001"), minNumber = cms.uint32(1))
process.CR001         = cms.Path(process.select3leptonsRegions * process.cand001 * process.cand001Filter)
process.CR001Counter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("CR001","pathFor3LeptonsAnalysis","zzTrigger"))

process.cand010       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("bareZWCand"), cut = cms.string(SEL010))
process.cand010Filter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand010"), minNumber = cms.uint32(1))
process.CR010         = cms.Path(process.select3leptonsRegions * process.cand010 * process.cand010Filter)
process.CR010Counter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("CR010","pathFor3LeptonsAnalysis","zzTrigger"))

process.cand000       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("bareZWCand"), cut = cms.string(SEL000))
process.cand000Filter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand000"), minNumber = cms.uint32(1))
process.CR000         = cms.Path(process.select3leptonsRegions * process.cand000 * process.cand000Filter)
process.CR000Counter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("CR000","pathFor3LeptonsAnalysis","zzTrigger"))


process.select3leptons1photonRegions = process.select3leptonsRegions.clone()
process.select3leptons1photonRegions.minPhotons = cms.int32(1)
process.select3leptons1photonRegions.maxPhotons = cms.int32(1)

process.SR3P_1L        = cms.Path(process.select3leptons1photonRegions * process.candSR3P * process.candSR3PFilter)
process.SR3P1LCounter = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("SR3P_1L","pathFor3LeptonsAnalysis","zzTrigger"))


# Path for lepton fake rate measuerement

PT20_10     = "((daughter(0).daughter(0).pt > 20 && daughter(0).daughter(1).pt > 10) || (daughter(0).daughter(0).pt() > 10 && daughter(0).daughter(1).pt > 20))"
ZMASSWINDOW = "abs(daughter(0).mass -91.19) <= 10"


### Basic object for the fake rate measurement CR
process.ZlSelected = cms.EDFilter("PATCompositeCandidateSelector",
                                  src = cms.InputTag("ZlCand"),
                                  cut = cms.string(PT20_10 + " && " + ZMASSWINDOW)
                                  )


from VVXAnalysis.Producers.EventFilter_cfg import eventFilter
process.selectLeptonFakeRateRegion = eventFilter.clone()
process.selectLeptonFakeRateRegion.minTightLeptons = cms.int32(2)
process.selectLeptonFakeRateRegion.minLooseLeptons = cms.int32(3)
process.selectLeptonFakeRateRegion.maxTightLeptons = cms.int32(3)
process.selectLeptonFakeRateRegion.maxLooseLeptons = cms.int32(3)
process.selectLeptonFakeRateRegion.minPhotons = cms.int32(0)
process.selectLeptonFakeRateRegion.maxPhotons = cms.int32(1000)


process.pathForLeptonFakeRateAnalysis = cms.Path(process.ZlSelected)
process.candZLFilter  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("ZlSelected"), minNumber = cms.uint32(1))
process.CRLFR         = cms.Path(process.selectLeptonFakeRateRegion * process.candZLFilter)
process.CRLFRCounter   = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("CRLFR","pathForLeptonFakeRateAnalysis","zzTrigger"))
