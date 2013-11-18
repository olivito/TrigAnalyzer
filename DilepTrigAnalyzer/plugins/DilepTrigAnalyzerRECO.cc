/** \class DilepTrigAnalyzerRECO
 *
 * See header file for documentation
 *
 *  $Date: 2012/01/30 09:40:35 $
 *  $Revision: 1.1 $
 *
 *  \author Dominick Olivito
 *
 */

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DilepTrigAnalyzer/DilepTrigAnalyzer/interface/DilepTrigAnalyzerRECO.h"

// need access to class objects being referenced to get their content!
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/HcalIsolatedTrack/interface/IsolatedPixelTrackCandidate.h"
#include "DataFormats/L1Trigger/interface/L1HFRings.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/TauReco/interface/PFTau.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <cassert>

//
// constructors and destructor
//
//____________________________________________________________________________
DilepTrigAnalyzerRECO::DilepTrigAnalyzerRECO(const edm::ParameterSet& ps) : 
  processName_(ps.getParameter<std::string>("processName")),
  triggerName_(ps.getParameter<std::string>("triggerName")),
  triggerResultsTag_(ps.getParameter<edm::InputTag>("triggerResults")),
  triggerEventTag_(ps.getParameter<edm::InputTag>("triggerEvent"))
{
  using namespace std;
  using namespace edm;

  cout << "DilepTrigAnalyzerRECO configuration: " << endl
       << "   ProcessName = " << processName_ << endl
       << "   TriggerName = " << triggerName_ << endl
       << "   TriggerResultsTag = " << triggerResultsTag_.encode() << endl
       << "   TriggerEventTag = " << triggerEventTag_.encode() << endl;

  dimuTriggerNames_.reserve(4);
  dimuTriggerNames_.at(mm) = "HLT_Mu17_Mu8_v23";
  dimuTriggerNames_.at(mmi) = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1";
  dimuTriggerNames_.at(mmtk) = "HLT_Mu17_TkMu8_v15";
  dimuTriggerNames_.at(mmitk) = "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1";

  emuTriggerNames_.reserve(4);
  emuTriggerNames_.at(em) = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10";
  emuTriggerNames_.at(emi) = "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1";
  emuTriggerNames_.at(me) = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10";
  emuTriggerNames_.at(mei) = "HLT_Mu17_TrkIsoVVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1";

  // histogram setup
  edm::Service<TFileService> fs;
  h_results_mm_ = fs->make<TH1F>("h_results_mm" , ";Trigger Results" , 16 , -0.5 , 15.5 );
  h_results_em_ = fs->make<TH1F>("h_results_em" , ";Trigger Results" , 16 , -0.5 , 15.5 );

}

//____________________________________________________________________________
DilepTrigAnalyzerRECO::~DilepTrigAnalyzerRECO()
{
}

//
// member functions
//
//____________________________________________________________________________
void
DilepTrigAnalyzerRECO::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  using namespace std;
  using namespace edm;

  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
      // check if trigger names in (new) config
      if (triggerName_=="dimu") {
	for (unsigned int itrig = 0; itrig < dimuTriggerNames_.size(); ++itrig) {
	  unsigned int triggerIndex(hltConfig_.triggerIndex(dimuTriggerNames_.at(itrig)));
	  if (triggerIndex>=n) {
	    cout << "DilepTrigAnalyzerRECO::analyze:"
		 << " TriggerName " << dimuTriggerNames_.at(itrig) 
		 << " not available in (new) config!" << endl;
	  }
	} // loop over dimu triggers
      } // if dimu triggers 
      else if (triggerName_=="emu") {
	for (unsigned int itrig = 0; itrig < emuTriggerNames_.size(); ++itrig) {
	  unsigned int triggerIndex(hltConfig_.triggerIndex(emuTriggerNames_.at(itrig)));
	  if (triggerIndex>=n) {
	    cout << "DilepTrigAnalyzerRECO::analyze:"
		 << " TriggerName " << emuTriggerNames_.at(itrig) 
		 << " not available in (new) config!" << endl;
	  }
	} // loop over emu triggers
      } // if emu triggers 
      else {
	cout << "DilepTrigAnalyzerRECO::analyze: TriggerName " << triggerName_ 
	     << " not recognized!!" << endl;
      }
    } // if changed
  } else {
    cout << "DilepTrigAnalyzerRECO::analyze:"
	 << " config extraction failure with process name "
	 << processName_ << endl;
  }

}

//____________________________________________________________________________
// ------------ method called to produce the data  ------------
void
DilepTrigAnalyzerRECO::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  
  cout << endl;

  // get event products
  iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
    cout << "DilepTrigAnalyzerRECO::analyze: Error in getting TriggerResults product from Event!" << endl;
    return;
  }
  iEvent.getByLabel(triggerEventTag_,triggerEventHandle_);
  if (!triggerEventHandle_.isValid()) {
    cout << "DilepTrigAnalyzerRECO::analyze: Error in getting TriggerEvent product from Event!" << endl;
    return;
  }
  // sanity check
  assert(triggerResultsHandle_->size()==hltConfig_.size());
  
  // analyze this event for the triggers requested
  if (triggerName_=="@") {
    const unsigned int n(hltConfig_.size());
    for (unsigned int i=0; i!=n; ++i) {
      analyzeTrigger(iEvent,iSetup,hltConfig_.triggerName(i));
    }
  } else if (triggerName_=="dimu") {
    unsigned int results(0);
    for (unsigned int itrig=0; itrig < dimuTriggerNames_.size(); ++itrig) {
      bool pass = analyzeTrigger(iEvent,iSetup,dimuTriggerNames_.at(itrig));
      if (pass) results |= 1 << itrig;
    }
    h_results_mm_->Fill(results);
  } else if (triggerName_=="emu") {
    unsigned int results(0);
    for (unsigned int itrig=0; itrig < emuTriggerNames_.size(); ++itrig) {
      bool pass = analyzeTrigger(iEvent,iSetup,emuTriggerNames_.at(itrig));
      if (pass) results |= 1 << itrig;
    }
    h_results_em_->Fill(results);
  } else {
    analyzeTrigger(iEvent,iSetup,triggerName_);
  }

  return;

}

//____________________________________________________________________________

bool DilepTrigAnalyzerRECO::analyzeTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& triggerName) {
  
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  cout << endl;

  const unsigned int n(hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));

  // abort on invalid trigger name
  if (triggerIndex>=n) {
    cout << "DilepTrigAnalyzerRECO::analyzeTrigger: path "
	 << triggerName << " - not found!" << endl;
    return false;
  }
  
  cout << "DilepTrigAnalyzerRECO::analyzeTrigger: path "
       << triggerName << " [" << triggerIndex << "]" << endl;
  // modules on this trigger path
  const unsigned int m(hltConfig_.size(triggerIndex));
  const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));

  bool wasRun = triggerResultsHandle_->wasrun(triggerIndex);
  bool accept = triggerResultsHandle_->accept(triggerIndex);
  bool error = triggerResultsHandle_->error(triggerIndex);
  // Results from TriggerResults product
  cout << " Trigger path status:"
       << " WasRun=" << wasRun
       << " Accept=" << accept
       << " Error =" << error
       << endl;
  const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
  cout << " Last active module - label/type: "
       << moduleLabels[moduleIndex] << "/" << hltConfig_.moduleType(moduleLabels[moduleIndex])
       << " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]"
       << endl;
  assert (moduleIndex<m);

  // loop over trigger and reco objects, match, make plots?

  // Results from TriggerEvent product - Attention: must look only for
  // modules actually run in this path for this event!
  for (unsigned int j=0; j<=moduleIndex; ++j) {
    const string& moduleLabel(moduleLabels[j]);
    const string  moduleType(hltConfig_.moduleType(moduleLabel));
    // check whether the module is packed up in TriggerEvent product
    const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(moduleLabel,"",processName_)));
    if (filterIndex<triggerEventHandle_->sizeFilters()) {
      cout << " 'L3' filter in slot " << j << " - label/type " << moduleLabel << "/" << moduleType << endl;
      const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
      const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
      const size_type nI(VIDS.size());
      const size_type nK(KEYS.size());
      assert(nI==nK);
      const size_type n(max(nI,nK));
      cout << "   " << n  << " accepted 'L3' objects found: " << endl;
      const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
      for (size_type i=0; i!=n; ++i) {
        const TriggerObject& TO(TOC[KEYS[i]]);
	// just need to apply cuts to trigger objects to find correct ones
        cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
             << TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()
             << endl;
      }
    }
  }


  return (wasRun && accept && !error);
}

//____________________________________________________________________________
bool DilepTrigAnalyzerRECO::isTriggerPassed(const std::string& triggerName) {

  const unsigned int n(hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));

  // abort on invalid trigger name
  if (triggerIndex>=n) {
    cout << "DilepTrigAnalyzerRECO::isTriggerPassed: path "
	 << triggerName << " - not found!" << endl;
    return false;
  }
  
  // modules on this trigger path
  const unsigned int m(hltConfig_.size(triggerIndex));
  const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));

  bool wasRun = triggerResultsHandle_->wasrun(triggerIndex);
  bool accept = triggerResultsHandle_->accept(triggerIndex);
  bool error = triggerResultsHandle_->error(triggerIndex);

  return (wasRun && accept && !error);
}

//____________________________________________________________________________
std::string DilepTrigAnalyzerRECO::lastModuleLabel(const std::string& triggerName) {

  const unsigned int n(hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));

  // abort on invalid trigger name
  if (triggerIndex>=n) {
    cout << "DilepTrigAnalyzerRECO::lastModuleLabel: path "
	 << triggerName << " - not found!" << endl;
    return false;
  }
  
  // modules on this trigger path
  const unsigned int m(hltConfig_.size(triggerIndex));
  const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));

  const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
  assert (moduleIndex<m);

  return moduleLabels[moduleIndex];
}
