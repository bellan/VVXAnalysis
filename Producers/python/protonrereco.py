import FWCore.ParameterSet.Config as cms
from Configuration.EventContent.EventContent_cff import *

# load common code
import direct_simu_reco_cff as profile
process = cms.Process('CTPPSTestAcceptance', profile.era)
profile.LoadConfig(process)
profile.config.SetDefaults(process)

# minimal logger settings
process.MessageLogger = cms.Service("MessageLogger",
  statistics = cms.untracked.vstring(),
  destinations = cms.untracked.vstring('cerr'),
  cerr = cms.untracked.PSet(
    threshold = cms.untracked.string('WARNING')
  )
)

# number of events
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(100000)
)

# override LHCInfo source
process.load("CalibPPS.ESProducers.ctppsLHCInfoRandomXangleESSource_cfi")
process.ctppsLHCInfoRandomXangleESSource.generateEveryNEvents = 1
process.ctppsLHCInfoRandomXangleESSource.xangleHistogramFile = "CrossingAngles2018.root"
process.ctppsLHCInfoRandomXangleESSource.xangleHistogramObject = "hxang"
process.ctppsLHCInfoRandomXangleESSource.beamEnergy = 6500.
process.ctppsLHCInfoRandomXangleESSource.betaStar = 0.40
process.esPreferLHCInfo = cms.ESPrefer("CTPPSLHCInfoRandomXangleESSource", "ctppsLHCInfoRandomXangleESSource")


# override beam-parameter source
process.load("CalibPPS.ESProducers.ctppsBeamParametersFromLHCInfoESSource_cfi")

process.ctppsBeamParametersFromLHCInfoESSource.beamDivX45 = process.ctppsBeamParametersESSource.beamDivX45
process.ctppsBeamParametersFromLHCInfoESSource.beamDivX56 = process.ctppsBeamParametersESSource.beamDivX56
process.ctppsBeamParametersFromLHCInfoESSource.beamDivY45 = process.ctppsBeamParametersESSource.beamDivY45
process.ctppsBeamParametersFromLHCInfoESSource.beamDivY56 = process.ctppsBeamParametersESSource.beamDivY56

process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetX45 = process.ctppsBeamParametersESSource.vtxOffsetX45
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetX56 = process.ctppsBeamParametersESSource.vtxOffsetX56
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetY45 = process.ctppsBeamParametersESSource.vtxOffsetY45
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetY56 = process.ctppsBeamParametersESSource.vtxOffsetY56
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetZ45 = process.ctppsBeamParametersESSource.vtxOffsetZ45
process.ctppsBeamParametersFromLHCInfoESSource.vtxOffsetZ56 = process.ctppsBeamParametersESSource.vtxOffsetZ56

process.ctppsBeamParametersFromLHCInfoESSource.vtxStddevX = process.ctppsBeamParametersESSource.vtxStddevX
process.ctppsBeamParametersFromLHCInfoESSource.vtxStddevY = process.ctppsBeamParametersESSource.vtxStddevY
process.ctppsBeamParametersFromLHCInfoESSource.vtxStddevZ = process.ctppsBeamParametersESSource.vtxStddevZ

process.ctppsBeamParametersESSource.beamDivX45
process.ctppsBeamParametersESSource.beamDivX56
process.ctppsBeamParametersESSource.beamDivY45
process.ctppsBeamParametersESSource.beamDivY56

# do not apply vertex smearing again                                                                         
process.ctppsBeamParametersESSource.vtxStddevX = 0
process.ctppsBeamParametersESSource.vtxStddevY = 0
process.ctppsBeamParametersESSource.vtxStddevZ = 0

# undo CMS vertex shift                                                                                      
#process.ctppsBeamParametersESSource.vtxOffsetX45 = +0.2475 * 1E-1
#process.ctppsBeamParametersESSource.vtxOffsetY45 = -0.6924 * 1E-1
#process.ctppsBeamParametersESSource.vtxOffsetZ45 = -8.1100 * 1E-1

# event source                                                                                                                  
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
      #'file:miniAOD_test100ev.root'
      #'/store/group/phys_pps/MC/requests_2018/private/AAZZ_bSM/AODSIM/GGToZZ_bSM_A0Z_1E-5_ACZ_0E0_13TeV-fpmc_100kEVT/GGToZZ_bSM_A0Z_1E-5_ACZ_0E0_13TeV-fpmc_100kEVT_AODSIM_0.root'
      'file:stepAOD1.root'
  )
)

# update settings of beam-smearing module
process.beamDivergenceVtxGenerator.src = cms.InputTag("")
process.beamDivergenceVtxGenerator.srcGenParticle = cms.VInputTag(
    cms.InputTag("genPUProtons","genPUProtons"),
    cms.InputTag("genParticles")
    #cms.InputTag("prunedGenParticles")
)

# acceptance plotter
process.ctppsAcceptancePlotter = cms.EDAnalyzer("CTPPSAcceptancePlotter",
  tagHepMC = cms.InputTag("generator", "unsmeared"),
  tagTracks = cms.InputTag("ctppsLocalTrackLiteProducer"),

  rpId_45_F = process.rpIds.rp_45_F,
  rpId_45_N = process.rpIds.rp_45_N,
  rpId_56_N = process.rpIds.rp_56_N,
  rpId_56_F = process.rpIds.rp_56_F,

  outputFile = cms.string("test_acceptance.root")
)

# distribution plotter
process.ctppsTrackDistributionPlotter = cms.EDAnalyzer("CTPPSTrackDistributionPlotter",
  tagTracks = cms.InputTag("ctppsLocalTrackLiteProducer"),
  x_pitch_pixels = cms.untracked.double(80E-3), # to be synchronised with process.ctppsDirectProtonSimulation.pitchPixelsVer

  rpId_45_F = process.rpIds.rp_45_F,
  rpId_45_N = process.rpIds.rp_45_N,
  rpId_56_N = process.rpIds.rp_56_N,
  rpId_56_F = process.rpIds.rp_56_F,

  outputFile = cms.string("test_acceptance_xy.root")
)

# processing path
process.p = cms.Path(
  process.generator
  * process.beamDivergenceVtxGenerator
  * process.ctppsDirectProtonSimulation
  * process.reco_local
  * process.ctppsProtons
  * process.ctppsAcceptancePlotter
  * process.ctppsTrackDistributionPlotter
)


process.out = cms.OutputModule('PoolOutputModule',
    fileName = cms.untracked.string('testoutput.root'),
    outputCommands = AODSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.out.outputCommands.append('keep *_*_*_*')

process.outpath = cms.EndPath(process.out)

process.schedule = cms.Schedule(
    process.p,
    process.outpath
)
