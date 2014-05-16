#ifndef HLTcore_TrigAnalyzerRECORef_h
#define HLTcore_TrigAnalyzerRECORef_h

/** \class TrigAnalyzerRECORef
 *
 *  
 *  This class is an EDAnalyzer analyzing the combined HLT information for RAW
 *    for dilepton triggers (based on HLTEventAnalyzerRAW)
 *
 *  $Date: 2012/01/30 09:40:35 $
 *  $Revision: 1.1 $
 *
 *  \author Dominick Olivito
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"

//
// class declaration
//
class TrigAnalyzerRECORef : public edm::EDAnalyzer {
  
  typedef math::XYZTLorentzVectorF LorentzVector;

 public:
  explicit TrigAnalyzerRECORef(const edm::ParameterSet&);
  ~TrigAnalyzerRECORef();

  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  trigger::VRmuon getTrigMuons(const edm::Event&, const edm::EventSetup&, const std::string& triggerName);

  int chargedHadronVertex(const reco::VertexCollection& vertices, 
			const reco::PFCandidate& pfcand ) const;
  int trackVertex(const reco::VertexCollection& vertices, 
			const reco::TrackRef& trackRef ) const;

 private:

  void bookHists(edm::Service<TFileService>& fs, const std::string& suffix = "");
  void fillHists(const reco::Muon& muon, const std::string& suffix = "");

  float muonPFiso(const reco::Muon& mu);

  /// module config parameters
  std::string   processName_;
  std::vector<std::string>   triggerNames_;
  std::vector<std::string>   triggerNamesShort_;
  edm::InputTag triggerResultsTag_;
  edm::InputTag triggerEventTag_;
  edm::InputTag triggerEventWithRefsTag_;
  edm::InputTag electronsInputTag_;
  edm::InputTag muonsInputTag_;
  edm::InputTag vtxInputTag_;
  bool reqTrigMatch_;
  std::vector<double> offPt_;
  bool offTight_;
  double offDxy_;
  bool doOffGenMatch_;
  edm::InputTag genParticlesTag_;
  bool verbose_;

  /// additional class data memebers
  edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
  edm::Handle<trigger::TriggerEventWithRefs> triggerEventWithRefsHandle_;
  HLTConfigProvider hltConfig_;
  TH1F* h_results_;
  edm::Handle<reco::VertexCollection> vertexHandle_;
  edm::Handle<reco::GsfElectronCollection> elsHandle_;
  edm::Handle<reco::MuonCollection> musHandle_;
  edm::Handle<reco::GenParticleCollection> genParticlesHandle_;

  std::map<std::string,TH1F*> hists_1d_;
  std::map<std::string,TH2F*> hists_2d_;

  unsigned int trigpass_results_;

};
#endif
