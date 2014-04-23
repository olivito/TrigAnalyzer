#ifndef HLTcore_SingleMuTrigAnalyzerRECO_h
#define HLTcore_SingleMuTrigAnalyzerRECO_h

/** \class SingleMuTrigAnalyzerRECO
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
class SingleMuTrigAnalyzerRECO : public edm::EDAnalyzer {
  
  typedef math::XYZTLorentzVectorF LorentzVector;

  struct StudyLepton {
    LorentzVector lv;
    float trkiso;
    float pfiso;
    int type;
    bool isHLT;
    bool isGen;
    float vz;
    int charge;
    int npixhits;
    int algo;
    int nhits;
    float dxy;
    float dz_hlt;
  };

 public:
  explicit SingleMuTrigAnalyzerRECO(const edm::ParameterSet&);
  ~SingleMuTrigAnalyzerRECO();

  enum hltTrigs {tkmu, mu, notrig};

  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  bool analyzeTrigger(const edm::Event&, const edm::EventSetup&, const hltTrigs triggerEnum);

  int chargedHadronVertex(const reco::VertexCollection& vertices, 
			const reco::PFCandidate& pfcand ) const;
  int trackVertex(const reco::VertexCollection& vertices, 
			const reco::TrackRef& trackRef ) const;

 private:

  void bookHists(edm::Service<TFileService>& fs, const std::string& suffix, bool hlt = false);
  void fillHists(const StudyLepton& lead, const StudyLepton& subl, const std::string& suffix, bool isHLT = false);
  void fillHistsRecoHLT(const StudyLepton& off_lead, const StudyLepton& off_subl, const StudyLepton& hlt_lead, const StudyLepton& hlt_subl, const std::string& suffix);

  void bookHistsGen(edm::Service<TFileService>& fs, const std::string& suffix);
  void fillHistsGen(const StudyLepton& mu, const std::string& suffix);

  float muonPFiso(const reco::Muon& mu);

  /// module config parameters
  std::string   processName_;
  std::string   baseTriggerName_;
  std::string   tkTriggerName_;
  edm::InputTag triggerResultsTag_;
  edm::InputTag triggerEventWithRefsTag_;
  edm::InputTag muonsInputTag_;
  edm::InputTag vtxInputTag_;
  edm::InputTag hltVtxInputTag_;
  bool dumpHLTTracks_;
  edm::InputTag hltTracksTag_;
  bool reqTrigMatch_;
  float offPt_;
  bool doOffGenMatch_;
  edm::InputTag genParticlesTag_;
  bool compareToGen_;
  bool verbose_;

  /// additional class data memebers
  edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
  edm::Handle<trigger::TriggerEventWithRefs> triggerEventWithRefsHandle_;
  HLTConfigProvider hltConfig_;
  std::vector<std::string> hltTriggerNames_;
  std::vector<std::string> hltShortNames_;
  TH1F* h_results_;
  edm::Handle<reco::VertexCollection> vertexHandle_;
  edm::Handle<reco::VertexCollection> hltVertexHandle_;
  edm::Handle<reco::MuonCollection> musHandle_;
  edm::Handle<reco::TrackCollection> hltTracksHandle_;
  edm::Handle<reco::GenParticleCollection> genParticlesHandle_;

  std::map<std::string,TH1F*> hists_1d_;
  std::map<std::string,TH2F*> hists_2d_;

  unsigned int trigpass_results_;

  /// payload extracted from TriggerEventWithRefs

  trigger::Vids        muonIds_;
  trigger::VRmuon      muonRefs_;

};
#endif
