import FWCore.ParameterSet.Config as cms

trigAnalyzerRECORef = cms.EDAnalyzer("TrigAnalyzerRECORef",
    processName = cms.string("reHLT"),
    triggerNames = cms.vstring(""),
    triggerNamesShort = cms.vstring(""),
    triggerResults = cms.InputTag("TriggerResults","","reHLT"),
    triggerEventWithRefs = cms.InputTag("hltTriggerSummary","","reHLT"),
    electronsInputTag = cms.InputTag("gsfElectrons"),
    muonsInputTag = cms.InputTag("muons"),
    vtxInputTag = cms.InputTag("offlinePrimaryVertices"),
    doOffGenMatch = cms.bool(False),
    genParticles = cms.InputTag("genParticles","","SIM"),
    verbose = cms.bool(False)
)
