import FWCore.ParameterSet.Config as cms

dilepTrigAnalyzerRECO = cms.EDAnalyzer("DilepTrigAnalyzerRECO",
    processName = cms.string("HLT1"),
    triggerName = cms.string("@"),
    triggerResults = cms.InputTag("TriggerResults","","HLT1"),
    triggerEventWithRefs = cms.InputTag("hltTriggerSummary","","HLT1")
)
