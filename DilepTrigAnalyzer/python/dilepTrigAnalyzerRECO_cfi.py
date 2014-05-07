import FWCore.ParameterSet.Config as cms

dilepTrigAnalyzerRECO = cms.EDAnalyzer("DilepTrigAnalyzerRECO",
    processName = cms.string("HLT1"),
    triggerName = cms.string("@"),
    triggerResults = cms.InputTag("TriggerResults","","HLT1"),
    triggerEventWithRefs = cms.InputTag("hltTriggerSummary","","HLT1"),
    electronsInputTag = cms.InputTag("gsfElectrons"),
    muonsInputTag = cms.InputTag("muons"),
    vtxInputTag = cms.InputTag("offlinePrimaryVertices"),
    doOffGenMatch = cms.bool(False),
    genParticles = cms.InputTag("genParticles","","SIM"),
    verbose = cms.bool(False)
)
