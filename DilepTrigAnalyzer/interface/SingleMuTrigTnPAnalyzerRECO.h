#ifndef HLTcore_SingleMuTrigTnPAnalyzerRECO_h
#define HLTcore_SingleMuTrigTnPAnalyzerRECO_h

/** \class SingleMuTrigTnPAnalyzerRECO
 *
 *  
 *  This class is an EDAnalyzer analyzing the HLT information for AOD
 *    for single electron triggers (based on HLTEventAnalyzerAOD)
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

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"

typedef math::XYZTLorentzVectorF LorentzVector;

//
// class declaration
//
class SingleMuTrigTnPAnalyzerRECO : public edm::EDAnalyzer {

 public:
  explicit SingleMuTrigTnPAnalyzerRECO(const edm::ParameterSet&);
  ~SingleMuTrigTnPAnalyzerRECO();

  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  bool analyzeTrigger(const edm::Event&, const edm::EventSetup&, const std::string& triggerName);

 private:

  void bookHists(edm::Service<TFileService>& fs, const std::string& suffix);
  void fillHists(const LorentzVector& lv, const std::string& suffix);
  std::vector<LorentzVector> getTriggerObjects(const unsigned int& triggerIndex, const std::string& targetModuleLabel = "");
  float muonPFiso(const reco::Muon& mu);

  /// module config parameters
  std::string   processName_;
  std::string   refTriggerName_;
  std::string   sigTriggerName_;
  std::string   denomFilterName_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerEventToken_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muonsToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  
  edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
  edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
  HLTConfigProvider hltConfig_;
  edm::Handle<reco::VertexCollection> vertexHandle_;
  edm::Handle<edm::View<reco::Muon> > musHandle_;
  
  double tagPt_;
  double tagEta_;
  double probePt_;
  double probeEta_;
  double probePtForEta_;
  bool verbose_;

  std::map<std::string,TH1F*> hists_1d_;

};
#endif
