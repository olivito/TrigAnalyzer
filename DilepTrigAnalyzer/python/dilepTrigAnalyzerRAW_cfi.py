import FWCore.ParameterSet.Config as cms

dilepTrigAnalyzerRAW = cms.EDAnalyzer("DilepTrigAnalyzerRAW",
    processName = cms.string("HLT"),
    triggerName = cms.string("@"),
    triggerResults = cms.InputTag("TriggerResults","","HLT"),
    triggerEventWithRefs = cms.InputTag("hltTriggerSummaryRAW","","HLT")
)
