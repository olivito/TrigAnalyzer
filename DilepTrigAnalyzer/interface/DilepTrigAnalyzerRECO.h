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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"

//
// class declaration
//
class DilepTrigAnalyzerRECO : public edm::EDAnalyzer {
  
  typedef math::XYZTLorentzVectorF LorentzVector;

 public:
  explicit DilepTrigAnalyzerRECO(const edm::ParameterSet&);
  ~DilepTrigAnalyzerRECO();

  enum hltTrigs {mm, mmi, mmtk, mmitk, em, emi, me, mei};

  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  bool analyzeTrigger(const edm::Event&, const edm::EventSetup&, const hltTrigs triggerEnum);

 private:

  // bool isTriggerPassed(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& triggerName);
  // std::string lastModuleLabel(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& triggerName);

  void bookHists(edm::Service<TFileService>& fs, const std::string& suffix);
  void fillHists(const LorentzVector& lead, const LorentzVector& subl, const std::string& suffix, bool hlt);

  /// module config parameters
  std::string   processName_;
  std::string   triggerName_;
  edm::InputTag triggerResultsTag_;
  edm::InputTag triggerEventTag_;
  edm::InputTag electronsInputTag_;
  edm::InputTag muonsInputTag_;
  bool verbose_;

  /// additional class data memebers
  edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
  edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
  HLTConfigProvider hltConfig_;
  std::vector<std::string> hltTriggerNames_;
  std::vector<std::string> hltShortNames_;
  TH1F* h_results_mm_;
  TH1F* h_results_em_;

  std::map<std::string,TH1F*> hists_1d_;

};
#endif
