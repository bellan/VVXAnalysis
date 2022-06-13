from datetime import datetime
import FWCore.ParameterSet.Config as cms
from Configuration.EventContent.EventContent_cff import *
from Configuration.Eras.Modifier_ctpps_2018_cff import ctpps_2018
from Configuration.ProcessModifiers.run2_miniAOD_UL_cff import run2_miniAOD_UL

# load common code
import direct_simu_reco_cff as profile
process = cms.Process('CTPPSTestAcceptance', profile.era)
profile.LoadConfig(process)
profile.config.SetDefaults(process)

process.load("Validation.CTPPS.simu_config.year_2018_cff")
process.ctppsRPAlignmentCorrectionsDataESSourceXML.MisalignedFiles = ["PPtoPPWWjets/PPtoPPWWjets/python/PPS_2018_Alignments/2018_postTS2.xml"]
process.ctppsRPAlignmentCorrectionsDataESSourceXML.RealFiles = ["PPtoPPWWjets/PPtoPPWWjets/python/PPS_2018_Alignments/2018_postTS2.xml"]

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
  optionalPSet = cms.untracked.bool(True),
  reportEvery = cms.untracked.int32(1),
    limit = cms.untracked.int32(-1),
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


from pileupSamples_cfi import *
print "Mixing PU samples: "
print "\n".join(pileupSamples(2018,"D")[:3])
print "...\n" 
process.load("protonPreMix.protonPreMix.ctppsProtonMixer_cfi")
process.ctppsProtonMixer.PUFilesList = cms.vstring(*pileupSamples(2018,"D"))
process.ctppsProtonMixer.Verbosity = 0
import FWCore.PythonUtilities.LumiList as LumiList
pu_jsonFile = "/eos/project-c/ctpps/Operations/DataExternalConditions/2018/CMSgolden_2RPGood_anyarms_Era"+"D"+".json"
process.ctppsProtonMixer.lumisToProcess = LumiList.LumiList(filename = pu_jsonFile).getVLuminosityBlockRange()
print "Using JSON file for PU: "+pu_jsonFile+"\n"
process.RandomNumberGeneratorService.ctppsProtonMixer = cms.PSet(initialSeed = cms.untracked.uint32(datetime.now().time().microsecond))

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
  fileNames = cms.untracked.vstring('file:../../../stepAOD1.root')
)

# update settings of beam-smearing module
process.beamDivergenceVtxGenerator.src = cms.InputTag("")
process.beamDivergenceVtxGenerator.srcGenParticle = cms.VInputTag(
    cms.InputTag("genParticles")
)

# rng service for efficiency
process.load("protonPreMix.protonPreMix.ppsEfficiencyProducer_cfi")

process.RandomNumberGeneratorService.ppsEfficiencyProducer = cms.PSet(initialSeed = cms.untracked.uint32(datetime.now().time().microsecond))
process.ppsEfficiencyProducer.year = cms.int32(2018)
process.ppsEfficiencyProducer.era = cms.string("D2")
process.ppsEfficiencyProducer.mixedProtonsSrc = cms.InputTag("ctppsProtonMixer")
process.ppsEfficiencyProducer.efficiencyFileName_Near = cms.string("file:/afs/cern.ch/user/a/abellora/public/forGiovanni/Efficiency_reMiniAOD/pixelEfficiencies_radiation_reMiniAOD.root")  
process.ppsEfficiencyProducer.efficiencyFileName_Far = cms.string("file:/afs/cern.ch/user/a/abellora/public/forGiovanni/Efficiency_reMiniAOD/pixelEfficiencies_multiRP_reMiniAOD.root")


# Processing path
process.p = cms.Path(
  process.generator
  * process.beamDivergenceVtxGenerator
  * process.ctppsDirectProtonSimulation
  * process.reco_local
  * process.ctppsProtons
  * process.ctppsProtonMixer
  * process.ppsEfficiencyProducer
)


process.out = cms.OutputModule('PoolOutputModule',
    fileName = cms.untracked.string('file:../../../rereco1mix2.root'),
    outputCommands = AODSIMEventContent.outputCommands,
    #overrideBranchesSplitLevel = cms.untracked.VPSet(
    #    cms.untracked.PSet(
    #        branch = cms.untracked.string('recoForwardProtons_ctppsProtonMixer_*'),
    #        splitLevel = cms.untracked.int32(99)
    #    )
                               #),
    splitLevel = cms.untracked.int32(0)
)

process.out.outputCommands.append('keep *_*_*_*')
process.out.outputCommands.append('keep recoForwardProtons_ctppsProtonMixer__*')

process.outpath = cms.EndPath(process.out)

process.schedule = cms.Schedule(
    process.p,
    process.outpath
)
