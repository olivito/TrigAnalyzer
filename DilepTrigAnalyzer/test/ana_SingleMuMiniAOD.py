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

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
#   fileNames = cms.untracked.vstring('/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/730/00000/00378BFA-B121-E611-93ED-02163E0146D2.root'),
   fileNames = cms.untracked.vstring('file:skim.root'),
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
                                       fileName = cms.string('histos_100k_skim.root')
                                   )

### analyzer configuration

process.singleMuTrigAnalyzerMiniAOD = cms.EDAnalyzer("SingleMuTrigAnalyzerMiniAOD")
process.singleMuTrigAnalyzerMiniAOD.triggerName = cms.untracked.string("HLT_IsoTkMu22_v2")
process.singleMuTrigAnalyzerMiniAOD.verbose = cms.untracked.bool(False)

process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v8"

# Path and EndPath definitions
process.HLTanalyzers = cms.Path(process.singleMuTrigAnalyzerMiniAOD)
