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

process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
#   fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/o/olivito/public/hlt/CMSSW_8_0_14/src/output_doubletkmu_noiter2_v20_reco_Run2016B_2e33_PU17_1file.root'),
#   fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/o/olivito/public/hlt/CMSSW_8_0_14/src/output_doubletkmu_2losthits_v4_reco_Run2016B_2e33_PU17_1file.root'),
#   fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/o/olivito/public/hlt/CMSSW_8_0_14/src/output_doubletkmu_noiter2_v16_newgt_reco_Run2016D_8e33_PU28_1file.root'),
#   fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/o/olivito/public/hlt/CMSSW_8_0_14/src/output_doubletkmu_noiter2_v16_newgt_reco_Run2016F_11e33_PU38_1file.root'),
#   fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/o/olivito/public/hlt/CMSSW_8_0_14/src/output_doubletkmu_noiter2_v16_newgt_2losthits_iter0_reco_Run2016F_11e33_PU38_1file.root'),
#   fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/o/olivito/public/hlt/CMSSW_8_0_14/src/output_doubletkmu_noiter2_v20_reco_Run2016F_11e33_PU38_1file.root'),
#   fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/o/olivito/public/hlt/CMSSW_8_0_14/src/output_doubletkmu_2losthits_v4_reco_Run2016F_11e33_PU38_1file.root'),
   fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/o/olivito/public/hlt/CMSSW_8_0_14/src/output_doubletkmu_oldtkmu_v2_reco_Run2016F_HIPfix_10p5e33_PU35_1file.root'),
#   fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/o/olivito/public/hlt/CMSSW_8_0_14/src/output_doubletkmu_2losthits_v8_reco_Run2016F_HIPfix_10p5e33_PU35_1file.root'),
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
#                                       fileName = cms.string('histos_2losthits_v4_L1denom_TkMu24_2e33_PU17_1file.root')
#                                       fileName = cms.string('histos_v16_newgt_reHLT_match_MuMuDZ_8e33_PU28_1file.root')
#                                       fileName = cms.string('histos_2losthits_v4_L1denom_TkMu24_11e33_PU38_1file.root')
                                       fileName = cms.string('histos_oldtkmu_v2_L1denom_TkMu24_HIPfix_10p5e33_PU35_1file.root')
#                                       fileName = cms.string('histos_2losthits_v8_L1denom_TkMu24_HIPfix_10p5e33_PU35_1file.root')
                                   )

### analyzer configuration

process.singleMuTrigTnPAnalyzerRECO = cms.EDAnalyzer("SingleMuTrigTnPAnalyzerRECO")
process.singleMuTrigTnPAnalyzerRECO.processName = cms.untracked.string("reHLT")
process.singleMuTrigTnPAnalyzerRECO.triggerResultsTag = cms.untracked.InputTag("TriggerResults", "", "reHLT")

#process.singleMuTrigTnPAnalyzerRECO.processName = cms.untracked.string("HLT")
#process.singleMuTrigTnPAnalyzerRECO.triggerResultsTag = cms.untracked.InputTag("TriggerResults", "", "HLT")

process.singleMuTrigTnPAnalyzerRECO.refTriggerName = cms.untracked.string("HLT_IsoMu24_v2")
process.singleMuTrigTnPAnalyzerRECO.sigTriggerName = cms.untracked.string("HLT_TkMu24_v1")
process.singleMuTrigTnPAnalyzerRECO.denomFilterName = cms.untracked.string("hltL1fL1sMu22L1Filtered0")
process.singleMuTrigTnPAnalyzerRECO.verbose = cms.untracked.bool(False)

process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v8"

# Path and EndPath definitions
process.HLTanalyzers = cms.Path(process.singleMuTrigTnPAnalyzerRECO)
