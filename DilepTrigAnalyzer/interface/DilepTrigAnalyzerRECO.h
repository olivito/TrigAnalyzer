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

#include "TH1.h"

//
// class declaration
//
class DilepTrigAnalyzerRECO : public edm::EDAnalyzer {
  
 public:
  explicit DilepTrigAnalyzerRECO(const edm::ParameterSet&);
  ~DilepTrigAnalyzerRECO();

  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  bool analyzeTrigger(const edm::Event&, const edm::EventSetup&, const std::string& triggerName);

  enum mmTrigs {mm, mmi, mmtk, mmitk};
  enum emTrigs {em, emi, me, mei};

 private:

  bool isTriggerPassed(const std::string& triggerName);
  std::string lastModuleLabel(const std::string& triggerName);

  /// module config parameters
  std::string   processName_;
  std::string   triggerName_;
  edm::InputTag triggerResultsTag_;
  edm::InputTag triggerEventTag_;

  /// additional class data memebers
  edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
  edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
  HLTConfigProvider hltConfig_;
  std::vector<std::string> dimuTriggerNames_;
  std::vector<std::string> emuTriggerNames_;
  TH1F* h_results_mm_;
  TH1F* h_results_em_;

  std::vector<TH1F*> h_lead_pt_;
  std::vector<TH1F*> h_subl_pt_;
  std::vector<TH1F*> h_lead_eta_;
  std::vector<TH1F*> h_subl_eta_;

};
#endif
