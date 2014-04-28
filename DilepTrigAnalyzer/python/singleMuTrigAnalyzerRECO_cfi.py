import FWCore.ParameterSet.Config as cms

singleMuTrigAnalyzerRECO = cms.EDAnalyzer("SingleMuTrigAnalyzerRECO",
    processName = cms.string("reHLT"),
    triggerEventWithRefs = cms.InputTag("hltTriggerSummary","","reHLT"),
    muonsInputTag = cms.InputTag("muons"),
    vtxInputTag = cms.InputTag("offlinePrimaryVertices"),
    verbose = cms.bool(False)
)
