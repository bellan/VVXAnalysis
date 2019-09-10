#process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
   
if IsMC:
    process.jec = cms.ESSource("PoolDBESSource",
        DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(1)
        ),
        timetype = cms.string('runnumber'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
                    tag    = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016V4_MC_AK4PFchs'), #for 80X/Moriond17
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
                    tag    = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016V4_MC_AK8PFchs'), #for 80X/Moriond17
                    label  = cms.untracked.string('AK8PFchs')
                    ),
                ),
             connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/JEC/Summer16_23Sep2016V4_MC.db'), #for 80X/Moriond17
            )
else:
   process.jec = cms.ESSource("PoolDBESSource",
         DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(1)
            ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
                    tag    = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016AllV4_DATA_AK4PFchs'), #for 80X/Moriond17
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
                    tag    = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016AllV4_DATA_AK8PFchs'), #for 80X/Moriond17
                    label  = cms.untracked.string('AK8PFchs')
                    ),
                ), 
            connect = cms.string('sqlite_fip:ZZAnalysis/AnalysisStep/data/JEC/Summer16_23Sep2016AllV4_DATA.db'), #for 80X/Moriond17
)

### reapply JEC
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector

process.patJetCorrFactorsReapplyJECFat = updatedPatJetCorrFactors.clone(
                                    src     = cms.InputTag("slimmedJetsAK8"),
                                    levels  = ['L1FastJet','L2Relative','L3Absolute'],
                                    payload = 'AK8PFchs')

process.patJetsReapplyJECFat = updatedPatJets.clone(
                                    jetSource = cms.InputTag("slimmedJetsAK8"),
                                    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJECFat") ))

### FIXME: should we be dressing the fat jets as well?
process.goodJetsFat = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                         filterParams = pfJetIDSelector.clone(),
                         src = cms.InputTag("patJetsReapplyJECFat"),
                         filter = cms.bool(True) )

# Clean jets wrt. good (preFSR-)isolated leptons
process.cleanJetsFat = cms.EDProducer("JetsWithLeptonsRemover",
                                   Jets      = cms.InputTag("goodJetsFat"),
                                   Muons     = cms.InputTag("appendPhotons:muons"),
                                   Electrons = cms.InputTag("appendPhotons:electrons"),
                                   Diboson   = cms.InputTag(""),
                                   JetPreselection      = cms.string(""),
                                   MuonPreselection     = cms.string("userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"),
                                   ElectronPreselection = cms.string("userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"),
                                   DiBosonPreselection  = cms.string(""),
                                   MatchingType = cms.string("byDeltaR"),
                                   DeltaRCut = cms.untracked.double(0.8),  
                                   cleanFSRFromLeptons = cms.bool(True),
                                   DebugPlots = cms.untracked.bool(False)
                                   )

process.corrJetsProducer = cms.EDProducer ( "CorrJetsProducer",
                                    jets    = cms.InputTag( "cleanJetsFat" ),
                                    vertex  = cms.InputTag( "goodPrimaryVertices" ), 
                                    rho     = cms.InputTag( "fixedGridRhoFastjetAll"   ),
                                    payload = cms.string  ( "AK8PFchs" ),
                                    isData  = cms.bool    (  False ))

process.fatJets = cms.Path( process.patJetCorrFactorsReapplyJECFat + process.patJetsReapplyJECFat + process.goodJetsFat )
