import FWCore.ParameterSet.Config as cms

eventFilter = cms.EDFilter("EventFilter",
                           muons = cms.InputTag("appendPhotons:muons"),
                           tightMuonSelection = cms.string("userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"),
                           electrons = cms.InputTag("appendPhotons:electrons"),
                           tightElectronSelection = cms.string("userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"),

                           minLooseLeptons = cms.int32(2),
                           maxLooseLeptons = cms.int32(4),
                           minTightLeptons = cms.int32(0),
                           maxTightLeptons = cms.int32(4),

                           photons = cms.InputTag("filteredPhotons"),
                           photonSelection = cms.string(""),
                           minPhotons = cms.int32(1),        
                           maxPhotons = cms.int32(1),        

                           jetsAK4 = cms.InputTag("dressedJets"),
                           jetAK4Selection = cms.string(""),
                           minAK4s = cms.int32(0),           

                           jetsAK8 = cms.InputTag("correctedJetsAK8:corrJets"),
                           jetAK8Selection = cms.string(""),
                           minAK8s = cms.int32(0),           

                           minAK4orMinAK8 = cms.bool(True),     

                           MET = cms.InputTag("slimmedMETs"),    
                           METSelection = cms.string("")
)

