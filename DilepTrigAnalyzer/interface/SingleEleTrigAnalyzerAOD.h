#ifndef HLTcore_SingleEleTrigAnalyzerAOD_h
#define HLTcore_SingleEleTrigAnalyzerAOD_h

/** \class SingleEleTrigAnalyzerAOD
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

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"

//
// class declaration
//
class SingleEleTrigAnalyzerAOD : public edm::EDAnalyzer {
  
  typedef math::XYZTLorentzVectorF LorentzVector;

 public:
  explicit SingleEleTrigAnalyzerAOD(const edm::ParameterSet&);
  ~SingleEleTrigAnalyzerAOD();

  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  bool analyzeTrigger(const edm::Event&, const edm::EventSetup&, const std::string& triggerName);

 private:

  void bookHistsGen(edm::Service<TFileService>& fs, const std::string& suffix);
  void fillHistsGen(const LorentzVector& lv, const std::string& suffix);

  /// module config parameters
  std::string   processName_;
  std::string   triggerName_;
  edm::InputTag triggerResultsTag_;
  edm::InputTag triggerEventTag_;
  edm::InputTag genParticlesTag_;
  double genPt_;
  double genEta_;
  bool verbose_;

  /// additional class data memebers
  edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
  edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
  HLTConfigProvider hltConfig_;
  edm::Handle<reco::GenParticleCollection> genParticlesHandle_;

  std::map<std::string,TH1F*> hists_1d_;

};
#endif
