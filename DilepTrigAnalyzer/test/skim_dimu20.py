import FWCore.ParameterSet.Config as cms

process = cms.Process('SKIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
       '/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/730/00000/0C1BE802-B221-E611-87ED-02163E0144E5.root',
       '/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/730/00000/0E2F8FFA-B121-E611-9B01-02163E01339C.root',
       '/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/730/00000/1264A406-B221-E611-8934-02163E014602.root',
       '/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/730/00000/1AF27A2B-B221-E611-BB6E-02163E0137D0.root',
       '/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/730/00000/1CBC1100-B221-E611-BD93-02163E014594.root',
   ),
)


process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('RelVal nevts:100'),
    name = cms.untracked.string('Applications')
)

# JSON

import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_271036-273730_13TeV_PromptReco_Collisions16_JSON.txt').getVLuminosityBlockRange()

# skim definitions

process.highptMuons = cms.EDFilter("PtMinCandViewSelector",
                           src = cms.InputTag("slimmedMuons"),
                           ptMin = cms.double(20)
)

process.dimuonFilter = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag("highptMuons"),
  minNumber = cms.uint32(2),
)

process.dimu_filter_step = cms.Path(process.highptMuons + process.dimuonFilter)


# Output definition

process.Out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string ("skim_5files.root"),
                               outputCommands = cms.untracked.vstring('keep *'),
                               SelectEvents = cms.untracked.PSet(
                                   SelectEvents = cms.vstring('dimu_filter_step')
                               ),
)
  
process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v8"

# Path and EndPath definitions
process.output_step = cms.EndPath(process.Out)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
