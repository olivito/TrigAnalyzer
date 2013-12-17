#ifndef HLTcore_TrackCompAnalyzer_h
#define HLTcore_TrackCompAnalyzer_h

/** \class TrackCompAnalyzer
 *
 *  
 *  This class is an EDAnalyzer comparing online and offline iterative tracks
 *    regionally around HLT muons
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

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"

//
// class declaration
//
class TrackCompAnalyzer : public edm::EDAnalyzer {
  
  typedef math::XYZTLorentzVectorF LorentzVector;

 public:
  explicit TrackCompAnalyzer(const edm::ParameterSet&);
  ~TrackCompAnalyzer();

  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  //  bool analyzeTrigger(const edm::Event&, const edm::EventSetup&, const hltTrigs triggerEnum);

 private:

  void bookHists(edm::Service<TFileService>& fs, const std::string& suffix);
  void bookHistsRecoHLT(edm::Service<TFileService>& fs, const std::string& suffix);
  void fillHists(const reco::Track& trk, const std::string& suffix);
  void fillHistsRecoHLT(const reco::Track& off_trk, const reco::Track& hlt_trk, const std::string& suffix);

  float delta_phi(float phi1, float phi2);

  /// module config parameters
  std::string   processName_;
  std::string   triggerName_;
  edm::InputTag triggerResultsTag_;
  edm::InputTag triggerEventTag_;
  edm::InputTag triggerEventWithRefsTag_;
  edm::InputTag hltTracksInputTag_;
  edm::InputTag offlineTracksInputTag_;
  edm::InputTag beamSpotInputTag_;
  bool verbose_;

  /// additional class data memebers
  edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
  edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
  edm::Handle<trigger::TriggerEventWithRefs> triggerEventWithRefsHandle_;
  HLTConfigProvider hltConfig_;
  // std::vector<std::string> hltTriggerNames_;
  // std::vector<std::string> hltShortNames_;
  edm::Handle<reco::TrackCollection> hltTracksHandle_;
  edm::Handle<reco::TrackCollection> offlineTracksHandle_;
  edm::Handle<reco::BeamSpot> beamSpotHandle_;

  std::map<std::string,TH1F*> hists_1d_;
  std::map<std::string,TH2F*> hists_2d_;

  /// payload extracted from TriggerEventWithRefs

  // trigger::Vids        electronIds_;
  // trigger::VRelectron  electronRefs_;
  trigger::Vids        muonIds_;
  trigger::VRmuon      muonRefs_;

};
#endif
