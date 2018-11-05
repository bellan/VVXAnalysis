### ----------------------------------------------------------------------
### Based on ZZ->4l strategy.
###----------------------------------------------------------------------
#####################################################################################################################################################


### ------------------------------------------------------------------------- ###
### Activate some skimming and cleaning of the collections
### ------------------------------------------------------------------------- ###
SkimPaths.append("preselection")


### ......................................................................... ###
### Build collections of muons and electrons that pass a quality criteria (isGood + isolation) and that are NOT selected to form the ZZ best candidate that pass the full selection
### ......................................................................... ###


# Muons cleaning. First, create a muon collection from the best ZZ candidate grand daughters
process.muonsFromZZ = cms.EDProducer("PATMuonsFromCompositeCandidates", src =  cms.InputTag("ZZFiltered"), SplitLevel = cms.int32(1))

# Muons cleaning. Second, remove from the muon collection the muons that come from the best ZZ candidate that pass the full selection (previous collection).
# The newly produced collection is also filtered in muon quality and isolation
process.postCleaningMuons = cms.EDProducer("PATMuonCleaner",
                                           # pat electron input source
                                           src = cms.InputTag("appendPhotons:muons"),
                                           # preselection (any string-based cut for pat::Muons)
                                           preselection = cms.string("pt > 5 && userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"),
                                           # overlap checking configurables
                                           checkOverlaps = cms.PSet(
        muons = cms.PSet(
            src       = cms.InputTag("muonsFromZZ"), # Start from loose lepton def
            algorithm = cms.string("byDeltaR"),
            preselection        = cms.string(""), 
            deltaR              = cms.double(0.05),  
            checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
            pairCut             = cms.string(""),
            requireNoOverlaps   = cms.bool(True), # overlaps don't cause the electron to be discared
            )
        ),
                                           # finalCut (any string-based cut for pat::Muons)
                                           finalCut = cms.string(''),
                                           )


# Electrons cleaning. First, create a electron collection from the best ZZ candidate grand daughters
process.electronsFromZZ = cms.EDProducer("PATElectronsFromCompositeCandidates", src =  cms.InputTag("ZZFiltered"), SplitLevel = cms.int32(1))

# Electrons cleaning. Second, remove from the electron collection the electrons that come from the best ZZ candidate that pass the full selection (previous collection).
# The newly produced collection is also filtered in electron quality and isolation
process.postCleaningElectrons = cms.EDProducer("PATElectronCleaner",
                                               # pat electron input source
                                               src = cms.InputTag("appendPhotons:electrons"),
                                               # preselection (any string-based cut for pat::Electron)
                                              preselection = cms.string("pt > 7 && userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"),
                                               # overlap checking configurables
                                               checkOverlaps = cms.PSet(
        electrons = cms.PSet(
            src       = cms.InputTag("electronsFromZZ"), # Start from loose lepton def
            algorithm = cms.string("byDeltaR"),
            preselection        = cms.string(""), #
            deltaR              = cms.double(0.05),  
            checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
            pairCut             = cms.string(""),
            requireNoOverlaps   = cms.bool(True), # overlaps don't cause the electron to be discared
            )
        ),
                                               # finalCut (any string-based cut for pat::Electron)
                                               finalCut = cms.string(''),
                                               )


### ......................................................................... ###
# Remove from the event the jets that have leptons from the ZZ best candidate.
# Jets are also checked against other good isolated leptons not coming from the ZZ best candidate. 
# The jets, to be stored in the event, must pass the preselction (specified below by the user) AND the looseID + PU veto that is implemented in the code (same algo as for H->ZZ VBF selection)
### ......................................................................... ###

## FIXME: Logic need to be recheck as of new FSR strategy has been implemented
process.disambiguatedJets = cms.EDProducer("JetsWithLeptonsRemover",
                                           JetPreselection      = cms.string("pt > 20"),
                                           DiBosonPreselection  = cms.string(""),
                                           MuonPreselection     = cms.string(""),
                                           ElectronPreselection = cms.string(""),
                                           MatchingType         = cms.string("byDeltaR"), 
                                           Jets      = cms.InputTag("dressedJets"),
                                           Muons     = cms.InputTag("postCleaningMuons"),
                                           Electrons = cms.InputTag("postCleaningElectrons"),
                                           Diboson   = cms.InputTag("ZZFiltered"),
                                           cleanFSRFromLeptons = cms.bool(True),
                                           DebugPrintOuts      =  cms.untracked.bool(False),
                                           DebugPlots= cms.untracked.bool(False)
                                           )




### ......................................................................... ###
### Clean the 4 lepton candidates to select only the best possible candidate, requiring that passes the FULL selection
### ......................................................................... ###

Z1MASS            = "daughter('Z1').mass>60 && daughter('Z1').mass<120"
Z2MASS            = "daughter('Z2').mass>60 && daughter('Z2').mass<120"
ZZWITHONSHELLZS   = (BESTCAND_AMONG + "&&" + Z1MASS + "&&" + Z2MASS)



NOGHOST4lOF = (" ( ( (abs(daughter(0).daughter(0).pdgId) != abs(daughter(1).daughter(0).pdgId)) && "+
                    "( (deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).daughter(0).eta, daughter(1).daughter(0).phi)>0.05 ) && " +
                      "(deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).daughter(1).eta, daughter(1).daughter(1).phi)>0.05 ) && " +
                      "(deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).daughter(0).eta, daughter(1).daughter(0).phi)>0.05 ) && " +
                      "(deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).daughter(1).eta, daughter(1).daughter(1).phi)>0.05 ) ) ) || " + 
                  "(abs(daughter(0).daughter(0).pdgId)==abs(daughter(1).daughter(0).pdgId) )  )") 


process.ZZCand.flags.SR_ZZOnShell = cms.string(SR + "&&" + ZZWITHONSHELLZS + "&&" + NOGHOST4lOF) 
# Uncomment the lines below if you want a smaller finding region!
#process.ZZCand.bestCandAmong = cms.PSet(isBestCand = cms.string(ZZWITHONSHELLZS))


process.ZZSelectedCand = cms.EDFilter("PATCompositeCandidateSelector",
                                      src = cms.InputTag("ZZCand"),
                                      cut = cms.string("userFloat('isBestCand') && userFloat('SR')")
                                      #cut = cms.string("userFloat('isBestCand') && userFloat('SR_ZZOnShell')")
                                      )

### ......................................................................... ###
### Clean the Z+2 lepton candidates to select only the best possible candidate for that region, requiring that at least one lepton fails the FULL selection
### ......................................................................... ###


CR_Z1MASS_LARGE  = "daughter(0).mass>40 && daughter(0).mass<120"
CR_Z1MASS  = "daughter(0).mass>60 && daughter(0).mass<120"
CR_Z2MASS  = "daughter(1).mass>60 && daughter(1).mass<120"
# Value here below are the ones used for the H->ZZ analysis and here for cross-check for dedicated studies. 
#CR_Z1MASS  = "daughter(0).mass>40 && daughter(0).mass<120"
#CR_Z2MASS  = "daughter(1).mass>12 && daughter(1).mass<120"
CR_ZLLos_MASSWINDOW = (CR_Z1MASS    + "&&" + CR_Z2MASS)

# CR 3P1F
process.ZLLCand.flags.CRZLLos_3P1F_ZZOnShell = cms.string(CR_ZLLosSEL_3P1F + "&&" + CR_ZLLos_MASSWINDOW)
# Uncomment the lines below if you want a smaller finding region!
# process.ZLLCand.bestCandAmong.isBestCRZLLos_3P1F = cms.string(CR_BESTZLLos_3P1F + CR_ZLLos_MASSWINDOW)

# CR 2P2F
process.ZLLCand.flags.CRZLLos_2P2F_ZZOnShell = cms.string(CR_ZLLosSEL_2P2F + "&&" + CR_ZLLos_MASSWINDOW)
# Uncomment the lines below if you want a smaller finding region!
# process.ZLLCand.bestCandAmong.isBestCRZLLos_2P2F = cms.string(CR_BESTZLLos_2P2F + CR_ZLLos_MASSWINDOW)

process.ZLLFiltered3P1F = cms.EDFilter("PATCompositeCandidateSelector",
                                       src = cms.InputTag("ZLLCand"),
                                       cut = cms.string("userFloat('isBestCRZLLos_3P1F') && userFloat('CRZLLos_3P1F')")
                                       )

process.ZLLFiltered2P2F = cms.EDFilter("PATCompositeCandidateSelector",
                                       src = cms.InputTag("ZLLCand"),
                                       cut = cms.string("userFloat('isBestCRZLLos_2P2F') && userFloat('CRZLLos_2P2F')")
                                       )


# Merger of all ZZ final states.
process.ZZFiltered = cms.EDProducer("PATCompositeCandidateMergerWithPriority",
                                    src = cms.VInputTag(cms.InputTag("ZZSelectedCand"),
                                                        cms.InputTag("ZLLFiltered2P2F"), cms.InputTag("ZLLFiltered3P1F")),
                                    priority = cms.vint32(1,0,0)
                                    )

### ------------------------------------------------------------------------- ###
### Define the post reconstruction cleaning sequence
### ------------------------------------------------------------------------- ###

process.postRecoCleaning = cms.Sequence( process.ZZSelectedCand
                                         + process.ZLLFiltered2P2F
                                         + process.ZLLFiltered3P1F
                                         + process.ZZFiltered
                                         + process.muonsFromZZ*process.postCleaningMuons 
                                         + process.electronsFromZZ*process.postCleaningElectrons
                                         + process.disambiguatedJets
                                         )



### ------------------------------------------------------------------------- ###
### Run the preselection
### ------------------------------------------------------------------------- ###

# Skim counters
process.prePreselectionCounter       = cms.EDProducer("EventCountProducer")
#process.preSkimSignalCounter         = cms.EDProducer("EventCountProducer")
process.postPreselectionCounter      = cms.EDProducer("EventCountProducer")


### Some filters

# Number of disambiguated jets
process.jetCounterFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("disambiguatedJets"), minNumber = cms.uint32(0))

# Select only events with one such candidate
#process.zzCounterFilter  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("ZZFiltered"), minNumber = cms.uint32(0))

# Looser preselection: ask only for a at least a Z + 1 soft lepton
#process.zlCounterFilter  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("ZlCand"), minNumber = cms.uint32(1))







### If it is MC, run also the signal definition path
if IsMC:
    # Empty sequence to attach the signal filter (if specified in the CSV file)
    process.mcSelectionCounter = cms.EDProducer("EventCountProducer") # not really needeed... it is mainly an hack to get the path executed
    process.signalFilters = cms.Sequence(process.mcSelectionCounter) 
    process.mcSelection   = cms.Path(process.signalFilters)
    MCFILTER = "mcSelection"

    genCategory =  cms.EDFilter("ZZGenFilterCategory",
                                Topology       = cms.int32(SIGNALDEFINITION), 
                                ParticleStatus = cms.int32(1), 
                                src            = cms.InputTag("prunedGenParticles"),
                                GenJets        = cms.InputTag("slimmedGenJets"),
                                )
    process.genCategory = genCategory

    process.kFactor = cms.EDProducer('kfactorProducer',
                                     isMC  = cms.untracked.bool(IsMC),
                                     src            = cms.InputTag("prunedGenParticles"))
    
 
    process.signalCounter    = cms.EDProducer("EventCountProducer")
    process.signalDefinition = cms.Path(process.genCategory * process.kFactor * process.signalCounter)








process.zzAndzlFilterCombiner = cms.EDFilter("ZLFilter", ZLL = cms.InputTag("ZZFiltered"), ZL = cms.InputTag("ZlCand"),
                                             isMC         = cms.untracked.bool(IsMC),                                            
                                             ZLSelection = cms.string("((daughter(0).daughter(0).pt > 20 && daughter(0).daughter(1).pt > 10) || (daughter(0).daughter(0).pt() > 10 && daughter(0).daughter(1).pt > 20)) && abs(daughter(0).mass -91.19) <= 10") # Add SIP cut here and remove it in eventanalyzer.cc
                                             #is MC, topology,  SR, CRs,  
                                                               )


### Path that pre-select the higher level objects that will input the TreePlanter
process.preselection = cms.Path( process.prePreselectionCounter
                                 * process.CR
                                 * process.postRecoCleaning 
                                 * process.zzAndzlFilterCombiner * process.jetCounterFilter
                                 * process.postPreselectionCounter)


# Some counters and paths functional to the fake lepton background estimation and signal region check


process.zzTrigger = cms.EDFilter("ZZTriggerFilter", src = cms.InputTag("ZZFiltered"),
                                 isMC         = cms.untracked.bool(IsMC),
                                 setup        = cms.int32(LEPTON_SETUP),
                                 sampleType   = cms.int32(SAMPLE_TYPE),
                                 PD           = cms.string(PD),
                                 skimPaths    = cms.vstring(SkimPaths),
                                 MCFilterPath = cms.string(MCFILTER)
                                 )

process.cand2P2F       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("ZZFiltered"), cut = cms.string(BOTHFAIL))
process.cand2P2FFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand2P2F"), minNumber = cms.uint32(1))
process.cr2P2F         = cms.Path(process.cand2P2F * process.cand2P2FFilter)
process.cr2P2FCounter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("cr2P2F","preselection","zzTrigger"))

process.cand3P1F       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("ZZFiltered"), cut = cms.string(PASSD0_XOR_PASSD1))
process.cand3P1FFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("cand3P1F"), minNumber = cms.uint32(1))
process.cr3P1F         = cms.Path(process.cand3P1F * process.cand3P1FFilter)
process.cr3P1FCounter  = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("cr3P1F","preselection","zzTrigger"))

process.candSR       = cms.EDFilter("PATCompositeCandidateSelector", src = cms.InputTag("ZZFiltered"), cut = cms.string(BOTHPASS))
process.candSRFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("candSR"), minNumber = cms.uint32(1))
process.sr           = cms.Path(process.candSR * process.candSRFilter)
process.srCounter    = cms.EDProducer("SelectedEventCountProducer", names = cms.vstring("sr","preselection","zzTrigger"))



