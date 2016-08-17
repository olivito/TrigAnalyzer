#ifndef HLTcore_DiMuTrigAnalyzerRECO_h
#define HLTcore_DiMuTrigAnalyzerRECO_h

/** \class DiMuTrigAnalyzerRECO
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

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"

typedef math::XYZTLorentzVectorF LorentzVector;

//
// class declaration
//
class DiMuTrigAnalyzerRECO : public edm::EDAnalyzer {

 public:
  explicit DiMuTrigAnalyzerRECO(const edm::ParameterSet&);
  ~DiMuTrigAnalyzerRECO();

  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

 private:

  std::vector<LorentzVector> getTriggerObjects(const unsigned int& triggerIndex);
  float muonPFiso(const reco::Muon& mu);

  /// module config parameters
  std::string   processName_;
  std::string   refTriggerName_;
  std::string   sigTriggerName_;
  bool doTrigObjectMatching_;
  bool verbose_;

  /// additional class data memebers
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerEventToken_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muonsToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  
  edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
  edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
  HLTConfigProvider hltConfig_;
  
  edm::Handle<reco::VertexCollection> vertexHandle_;
  edm::Handle<edm::View<reco::Muon> > musHandle_;

  std::map<std::string,TH1F*> hists_1d_;
  std::map<std::string,TH2F*> hists_2d_;

};
#endif
