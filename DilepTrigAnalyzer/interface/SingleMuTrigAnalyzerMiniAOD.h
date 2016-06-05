#ifndef HLTcore_SingleMuTrigAnalyzerMiniAOD_h
#define HLTcore_SingleMuTrigAnalyzerMiniAOD_h

/** \class SingleMuTrigAnalyzerMiniAOD
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
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"

//
// class declaration
//
class SingleMuTrigAnalyzerMiniAOD : public edm::EDAnalyzer {
  
  typedef math::XYZTLorentzVectorF LorentzVector;

 public:
  explicit SingleMuTrigAnalyzerMiniAOD(const edm::ParameterSet&);
  ~SingleMuTrigAnalyzerMiniAOD();

  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  bool analyzeTrigger(const edm::Event&, const edm::EventSetup&, const std::string& triggerName);

 private:

  void bookHists(edm::Service<TFileService>& fs, const std::string& suffix);
  void fillHists(const LorentzVector& lv, const std::string& suffix);
  float muonPFiso(const pat::Muon& mu);

  /// module config parameters
  std::string   processName_;
  std::string   triggerName_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectStandAloneToken_;
  edm::EDGetTokenT<edm::View<pat::Muon> > muonsToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  double tagPt_;
  double tagEta_;
  double probePt_;
  double probeEta_;
  bool verbose_;

  /// additional class data memebers
  edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
  HLTConfigProvider hltConfig_;

  std::map<std::string,TH1F*> hists_1d_;

};
#endif
