#ifndef HLTcore_DilepTrigAnalyzerRAW_h
#define HLTcore_DilepTrigAnalyzerRAW_h

/** \class DilepTrigAnalyzerRAW
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
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"

#include "TH1.h"

//
// class declaration
//
class DilepTrigAnalyzerRAW : public edm::EDAnalyzer {
  
 public:
  explicit DilepTrigAnalyzerRAW(const edm::ParameterSet&);
  ~DilepTrigAnalyzerRAW();

  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  bool analyzeTrigger(const edm::Event&, const edm::EventSetup&, const std::string& triggerName);

 private:

  /// module config parameters
  std::string   processName_;
  std::string   triggerName_;
  edm::InputTag triggerResultsTag_;
  edm::InputTag triggerEventWithRefsTag_;
  edm::InputTag isoDepMapTag_;

  /// additional class data memebers
  edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
  edm::Handle<trigger::TriggerEventWithRefs> triggerEventWithRefsHandle_;
  HLTConfigProvider hltConfig_;
  std::vector<std::string> dimuTriggerNames_;
  std::vector<std::string> emuTriggerNames_;
  TH1F* h_results_mm_;
  TH1F* h_results_em_;

  std::vector<TH1F*> h_lead_pt_;
  std::vector<TH1F*> h_subl_pt_;
  std::vector<TH1F*> h_lead_eta_;
  std::vector<TH1F*> h_subl_eta_;

  /// payload extracted from TriggerEventWithRefs

  trigger::Vids        electronIds_;
  trigger::VRelectron  electronRefs_;
  trigger::Vids        muonIds_;
  trigger::VRmuon      muonRefs_;
  trigger::Vids        jetIds_;
  trigger::VRjet       jetRefs_;
  trigger::Vids        pixtrackIds_;
  trigger::VRpixtrack  pixtrackRefs_;

  trigger::Vids        l1emIds_;
  trigger::VRl1em      l1emRefs_;
  trigger::Vids        l1muonIds_;
  trigger::VRl1muon    l1muonRefs_;
  trigger::Vids        l1jetIds_;
  trigger::VRl1jet     l1jetRefs_;

  trigger::Vids        pfjetIds_;
  trigger::VRpfjet     pfjetRefs_;

};
#endif
