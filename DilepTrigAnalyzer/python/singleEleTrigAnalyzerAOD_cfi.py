import FWCore.ParameterSet.Config as cms

singleEleTrigAnalyzerAOD = cms.EDAnalyzer("SingleEleTrigAnalyzerAOD",
    processName = cms.string("reHLT"),
    triggerResults = cms.InputTag("TriggerResults","","reHLT"),
    triggerEvent = cms.InputTag("hltTriggerSummaryAOD","","reHLT"),
    genParticles = cms.InputTag("genParticles"),
    genPt = cms.double(25.),
    genEta = cms.double(2.1),
    verbose = cms.bool(False)
)
