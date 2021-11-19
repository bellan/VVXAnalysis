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


PT20_10     = "((daughter(0).daughter(0).pt > 20 && daughter(0).daughter(1).pt > 10) || (daughter(0).daughter(0).pt() > 10 && daughter(0).daughter(1).pt > 20))"
ZMASSWINDOW = "abs(daughter(0).mass -91.19) <= 10"


### Basic object for the fake rate measurement CR
process.ZlSelected = cms.EDFilter("PATCompositeCandidateSelector",
                                  src = cms.InputTag("ZlCand"),
                                  cut = cms.string(PT20_10 + " && "+ZMASSWINDOW)
                                  )



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
                                              Z1Presel = cms.string("mass > 40 && mass < 120"),
                                          )

# combine the Z candidate with a free lepton.
process.ZlCandFromLooseL       = process.ZlCand.clone()
process.ZlCandFromLooseL.decay = cms.string('ZCandFromLooseL softLeptons')


#  Just check that the trigger thresholds are passed
process.selectedZlCandFromLooseL = cms.EDFilter("PATCompositeCandidateSelector",
                                                src = cms.InputTag("ZlCandFromLooseL"),
                                                cut = cms.string("(daughter(0).daughter(0).pt > 20 && (daughter(0).daughter(1).pt > 10 || daughter(1).pt > 10))" + 
                                                                 " || " +
                                                                 "(daughter(0).daughter(1).pt > 20 && (daughter(0).daughter(0).pt > 10 || daughter(1).pt > 10))" + 
                                                                 " || " +
                                                                 "(daughter(1).pt > 20 && (daughter(0).daughter(0).pt > 10 || daughter(0).daughter(1).pt > 10))"
                                                             ),
                                                checkCharge = cms.bool(False)
                                            )

# At this point we have all possible 3 leptons combinations. We now build the WZ bare candidate asking 20 GeV of MET
process.bareZWCand = cms.EDProducer("PATCandViewShallowCloneCombiner",
                                    decay = cms.string('selectedZlCandFromLooseL slimmedMETs'),
                                    cut = cms.string("daughter(1).pt > 20"),
                                    checkCharge = cms.bool(False)
)



# Filter to select the events
# it is a two stages filtering. One to reduce the computational time (preSelect), and another to do the proper selection (select)
from VVXAnalysis.Producers.EventFilter_cfg import eventFilter
process.preSelect3leptonsRegions = eventFilter.clone()
process.preSelect3leptonsRegions.minLooseLeptons = cms.int32(3)
process.preSelect3leptonsRegions.maxLooseLeptons = cms.int32(3)
process.preSelect3leptonsRegions.jetsAK4         = cms.InputTag("slimmedJets")
process.preSelect3leptonsRegions.muons           = cms.InputTag("slimmedMuons")
process.preSelect3leptonsRegions.electrons       = cms.InputTag("slimmedElectrons")


from VVXAnalysis.Producers.EventFilter_cfg import eventFilter
process.select3leptonsRegions = eventFilter.clone()
process.select3leptonsRegions.minTightLeptons = cms.int32(0)
process.select3leptonsRegions.minLooseLeptons = cms.int32(3)
process.select3leptonsRegions.maxTightLeptons = cms.int32(3)
process.select3leptonsRegions.maxLooseLeptons = cms.int32(3)


process.pathFor3LeptonsAnalysis = cms.Path(#process.preSelect3leptonsRegions *
                                           process.ZlSelected               +  # CR for fake rate mesurement 
                                           process.bareZCandFromLooseL      +  # Z from loose leptons
                                           process.ZCandFromLooseL          +  # best Z from all loose leptons
                                           process.ZlCandFromLooseL         +  # best Z + a free loose lepton
                                           process.selectedZlCandFromLooseL +  # best Z + a free loose lepton w/ trigger requirements
                                           process.bareZWCand)                 # Zl + MET (WZ bare candidate)



#process.ZlForWZCountFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("ZlForWZ"), 
#                                          minNumber = cms.uint32(1) # maxNumber does not exist
#                                          )


# process.ZlForWZ = cms.EDFilter("PATCompositeCandidateSelector",
#                               src = cms.InputTag("ZlSelected"),
#                               cut = cms.string("daughter(1).masterClone.userFloat('isGood')" + " && " +
#                                                "daughter(1).masterClone.userFloat('passCombRelIsoPFFSRCorr')"
#                                               ),
#                               checkCharge = cms.bool(False)
# )

# process.ZlForWZCountFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("ZlForWZ"), 
#                                           minNumber = cms.uint32(1) # maxNumber does not exist
#                                           )


# process.ZWCand = cms.EDProducer("PATCandViewShallowCloneCombiner",
#                                 decay = cms.string('ZlForWZ slimmedMETs'),
#                                 cut = cms.string("daughter(0).daughter(1).pt > 20" + " && " + "daughter(1).pt > 20" 
#                                 ),
#                                 checkCharge = cms.bool(False)
# )


#process.WZjjPath = cms.Path(process.ZlForWZ * process.ZlForWZCountFilter * process.ZWCand)

