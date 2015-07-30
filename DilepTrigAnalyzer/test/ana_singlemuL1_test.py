import FWCore.ParameterSet.Config as cms

process = cms.Process('HLTANALYZER')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load('TrigAnalyzer.DilepTrigAnalyzer.singleMuTrigAnalyzerL1RECO_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/o/olivito/hltsub/CMSSW_7_4_7/src/pickevents_251883_failtkmu.root'),
)


process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('RelVal nevts:100'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('histos_test.root')
                                   )

### analyzer configuration

process.singleMuTrigAnalyzerL1RECO.triggerName = cms.string("HLT_IsoTkMu20_v2")
process.singleMuTrigAnalyzerL1RECO.printL1MuonScales = cms.bool(False)
process.singleMuTrigAnalyzerL1RECO.verbose = cms.bool(False)

process.GlobalTag.globaltag = "74X_dataRun2_Prompt_v0"

# Path and EndPath definitions
process.HLTanalyzers = cms.Path(process.singleMuTrigAnalyzerL1RECO)
