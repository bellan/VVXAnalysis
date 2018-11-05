### ......................................................................... ###
### Build the W->jj candidate out of the previously disambiguated jet collection, restricted to the central regions
### ......................................................................... ###

process.centralJets = cms.EDFilter("EtaPtMinCandViewSelector", 
                                   src = cms.InputTag("disambiguatedJets"),
                                   ptMin   = cms.double(10),
                                   etaMin = cms.double(-2.4),
                                   etaMax = cms.double(2.4)
                                   )

process.bareVhadCand = cms.EDProducer("CandViewShallowCloneCombiner",
                                      decay = cms.string('centralJets centralJets'),
                                      cut = cms.string('mass > 20'), # protect against ghosts
                                      checkCharge = cms.bool(False))

process.VhadCand = cms.EDProducer("WCandidateFiller",
                                  src = cms.InputTag("bareVhadCand"))

process.VhadSequence = cms.Sequence(process.centralJets * process.bareVhadCand * process.VhadCand)

process.VZjjprePreselectionCounter       = cms.EDProducer("EventCountProducer")
process.VZjjpostPreselectionCounter      = cms.EDProducer("EventCountProducer")


process.VZjjPath = cms.Path(  process.VZjjprePreselectionCounter
                            * process.VhadSequence
                            * process.VZjjpostPreselectionCounter
                            )
