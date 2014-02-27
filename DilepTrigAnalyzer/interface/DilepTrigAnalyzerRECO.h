#ifndef HLTcore_DilepTrigAnalyzerRECO_h
#define HLTcore_DilepTrigAnalyzerRECO_h

/** \class DilepTrigAnalyzerRECO
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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"

//
// class declaration
//
class DilepTrigAnalyzerRECO : public edm::EDAnalyzer {
  
  typedef math::XYZTLorentzVectorF LorentzVector;

  struct StudyLepton {
    LorentzVector lv;
    float trkiso;
    float pfiso;
    int type;
    bool isHLT;
    float vz;
    int charge;
  };

 public:
  explicit DilepTrigAnalyzerRECO(const edm::ParameterSet&);
  ~DilepTrigAnalyzerRECO();

  enum hltTrigs {mmi, mm, mmitk, mmtk, emi, em, mei, me, notrig};

  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  bool analyzeTrigger(const edm::Event&, const edm::EventSetup&, const hltTrigs triggerEnum);

  int chargedHadronVertex(const reco::VertexCollection& vertices, 
			const reco::PFCandidate& pfcand ) const;
  int trackVertex(const reco::VertexCollection& vertices, 
			const reco::TrackRef& trackRef ) const;

 private:

  void bookHists(edm::Service<TFileService>& fs, const std::string& suffix, bool hlt = false);
  // void fillHists(const LorentzVector& lead, const LorentzVector& subl, const std::string& suffix, bool hlt);
  // void fillMuonIsoHists(const reco::MuonCollection& col, const int& lead_idx, const int& subl_idx, const std::string& suffix);
  // void fillHistsRecoHLT(const LorentzVector& off_lead, const LorentzVector& off_subl, const LorentzVector& hlt_lead, const LorentzVector& hlt_subl, const std::string& suffix, bool isem = false);
  void fillHists(const StudyLepton& lead, const StudyLepton& subl, const std::string& suffix, bool isHLT = false);
  void fillHistsRecoHLT(const StudyLepton& off_lead, const StudyLepton& off_subl, const StudyLepton& hlt_lead, const StudyLepton& hlt_subl, const std::string& suffix);

  float muonPFiso(const reco::Muon& mu);

  /// module config parameters
  std::string   processName_;
  std::string   triggerName_;
  std::string   mmBaseTriggerName_;
  std::string   mmIsoTriggerName_;
  std::string   mmtkBaseTriggerName_;
  std::string   mmtkIsoTriggerName_;
  edm::InputTag triggerResultsTag_;
  edm::InputTag triggerEventTag_;
  edm::InputTag triggerEventWithRefsTag_;
  edm::InputTag electronsInputTag_;
  edm::InputTag muonsInputTag_;
  edm::InputTag vtxInputTag_;
  edm::InputTag hltVtxInputTag_;
  bool getHLTIsoVals_;
  edm::InputTag isoValMapGlbTag_;
  edm::InputTag isoValMapTrkTag_;
  bool dumpHLTPFCands_;
  edm::InputTag hltTracksGlbTag_;
  edm::InputTag hltPFCandsGlbTag_;
  edm::InputTag hltPFCandsTrkTag_;
  edm::InputTag offPFCandsTag_;
  float offLeadPt_;
  float offSublPt_;
  bool verbose_;

  /// additional class data memebers
  edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
  edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
  edm::Handle<trigger::TriggerEventWithRefs> triggerEventWithRefsHandle_;
  HLTConfigProvider hltConfig_;
  std::vector<std::string> hltTriggerNames_;
  std::vector<std::string> hltShortNames_;
  TH1F* h_results_mm_;
  TH1F* h_results_offdilep_mm_;
  TH1F* h_results_em_;
  edm::Handle<reco::VertexCollection> vertexHandle_;
  edm::Handle<reco::VertexCollection> hltVertexHandle_;
  edm::Handle<reco::GsfElectronCollection> elsHandle_;
  edm::Handle<reco::MuonCollection> musHandle_;
  edm::Handle<edm::ValueMap<float> > isoValMapGlbHandle_;
  edm::Handle<edm::ValueMap<float> > isoValMapTrkHandle_;
  edm::Handle<reco::TrackCollection> hltTracksGlbHandle_;
  edm::Handle<reco::PFCandidateCollection> hltPFCandsGlbHandle_;
  edm::Handle<reco::PFCandidateCollection> hltPFCandsTrkHandle_;
  edm::Handle<reco::PFCandidateCollection> offPFCandsHandle_;

  std::map<std::string,TH1F*> hists_1d_;
  std::map<std::string,TH2F*> hists_2d_;

  unsigned int trigpass_results_;
  unsigned int trigpass_results_offdilep_;

  /// payload extracted from TriggerEventWithRefs

  trigger::Vids        electronIds_;
  trigger::VRelectron  electronRefs_;
  trigger::Vids        muonIds_;
  trigger::VRmuon      muonRefs_;

};
#endif
