/** \class DilepTrigAnalyzerRAW
 *
 * See header file for documentation
 *
 *  $Date: 2012/01/30 09:40:35 $
 *  $Revision: 1.14 $
 *
 *  \author Martin Grunewald
 *
 */

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "TrigAnalyzer/DilepTrigAnalyzer/interface/DilepTrigAnalyzerRAW.h"

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
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <cassert>

//
// constructors and destructor
//
//____________________________________________________________________________
DilepTrigAnalyzerRAW::DilepTrigAnalyzerRAW(const edm::ParameterSet& ps) : 
  processName_(ps.getParameter<std::string>("processName")),
  triggerName_(ps.getParameter<std::string>("triggerName")),
  triggerResultsTag_(ps.getParameter<edm::InputTag>("triggerResults")),
  triggerEventWithRefsTag_(ps.getParameter<edm::InputTag>("triggerEventWithRefs")),
  isoDepMapTag_(ps.getParameter<edm::InputTag>("isoDepMap"))
{
  using namespace std;
  using namespace edm;

  cout << "DilepTrigAnalyzerRAW configuration: " << endl
       << "   ProcessName = " << processName_ << endl
       << "   TriggerName = " << triggerName_ << endl
       << "   TriggerResultsTag = " << triggerResultsTag_.encode() << endl
       << "   TriggerEventWithRefsTag = " << triggerEventWithRefsTag_.encode() << endl
       << "   IsoDepMapTag = " << isoDepMapTag_.encode() << endl;

  dimuTriggerNames_.push_back("HLT_Mu17_Mu8_v23");
  dimuTriggerNames_.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1");
  dimuTriggerNames_.push_back("HLT_Mu17_TkMu8_v15");
  dimuTriggerNames_.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1");

  emuTriggerNames_.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10");
  emuTriggerNames_.push_back("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1");
  emuTriggerNames_.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10");
  emuTriggerNames_.push_back("HLT_Mu17_TrkIsoVVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1");

  // histogram setup
  edm::Service<TFileService> fs;
  h_results_mm_ = fs->make<TH1F>("h_results_mm" , ";Trigger Results" , 16 , -0.5 , 15.5 );
  h_results_em_ = fs->make<TH1F>("h_results_em" , ";Trigger Results" , 16 , -0.5 , 15.5 );

}

//____________________________________________________________________________
DilepTrigAnalyzerRAW::~DilepTrigAnalyzerRAW()
{
}

//
// member functions
//
//____________________________________________________________________________
void
DilepTrigAnalyzerRAW::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  using namespace std;
  using namespace edm;

  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
      // check if trigger name in (new) config
      if (triggerName_!="@") { // "@" means: analyze all triggers in config
	const unsigned int n(hltConfig_.size());
	if (triggerName_=="dimu") {
	  for (unsigned int itrig = 0; itrig < dimuTriggerNames_.size(); ++itrig) {
	    unsigned int triggerIndex(hltConfig_.triggerIndex(dimuTriggerNames_.at(itrig)));
	    if (triggerIndex>=n) {
	      cout << "DilepTrigAnalyzerRAW::analyze:"
		   << " TriggerName " << dimuTriggerNames_.at(itrig) 
		   << " not available in (new) config!" << endl;
	    }
	  } // loop over dimu triggers
	} // if dimu triggers 
	else if (triggerName_=="emu") {
	  for (unsigned int itrig = 0; itrig < emuTriggerNames_.size(); ++itrig) {
	    unsigned int triggerIndex(hltConfig_.triggerIndex(emuTriggerNames_.at(itrig)));
	    if (triggerIndex>=n) {
	      cout << "DilepTrigAnalyzerRAW::analyze:"
		   << " TriggerName " << emuTriggerNames_.at(itrig) 
		   << " not available in (new) config!" << endl;
	    }
	  } // loop over emu triggers
	} // if emu triggers 
	else { // one specific trigger
	  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
	  if (triggerIndex>=n) {
	    cout << "DilepTrigAnalyzerRAW::analyze:"
		 << " TriggerName " << triggerName_ 
		 << " not available in (new) config!" << endl;
	    cout << "Available TriggerNames are: " << endl;
	    hltConfig_.dump("Triggers");
	  }
	} // if one specific trigger
      } // not using all triggers
    } // if changed
  } else {
    cout << "DilepTrigAnalyzerRAW::analyze:"
	 << " config extraction failure with process name "
	 << processName_ << endl;
  }

}

//____________________________________________________________________________
// ------------ method called to produce the data  ------------
void
DilepTrigAnalyzerRAW::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  
  cout << endl;

  // get event products
  iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
    cout << "DilepTrigAnalyzerRAW::analyze: Error in getting TriggerResults product from Event!" << endl;
    return;
  }
  iEvent.getByLabel(triggerEventWithRefsTag_,triggerEventWithRefsHandle_);
  if (!triggerEventWithRefsHandle_.isValid()) {
    cout << "DilepTrigAnalyzerRAW::analyze: Error in getting TriggerEventWithRefs product from Event!" << endl;
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

bool DilepTrigAnalyzerRAW::analyzeTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& triggerName) {
  
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
    cout << "DilepTrigAnalyzerRAW::analyzeTrigger: path "
	 << triggerName << " - not found!" << endl;
    return false;
  }
  
  cout << "DilepTrigAnalyzerRAW::analyzeTrigger: path "
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

  // Results from TriggerEventWithRefs product
  electronIds_.clear();
  electronRefs_.clear();
  muonIds_.clear();
  muonRefs_.clear();
  // jetIds_.clear();
  // jetRefs_.clear();
  // pixtrackIds_.clear();
  // pixtrackRefs_.clear();
  // l1emIds_.clear();
  // l1emRefs_.clear();
  // l1muonIds_.clear();
  // l1muonRefs_.clear();
  // l1jetIds_.clear();
  // l1jetRefs_.clear();
  // l1hfringsIds_.clear();
  // l1hfringsRefs_.clear();
  // pfjetIds_.clear();
  // pfjetRefs_.clear();

  // Attention: must look only for modules actually run in this path
  // for this event!
  for (unsigned int j=0; j<=moduleIndex; ++j) {
    const string& moduleLabel(moduleLabels[j]);
    const string  moduleType(hltConfig_.moduleType(moduleLabel));
    // check whether the module is packed up in TriggerEventWithRef product
    const unsigned int filterIndex(triggerEventWithRefsHandle_->filterIndex(InputTag(moduleLabel,"",processName_)));
    if (filterIndex<triggerEventWithRefsHandle_->size()) {
      cout << " Filter in slot " << j << " - label/type " << moduleLabel << "/" << moduleType << endl;
      cout << " Filter packed up at: " << filterIndex << endl;
      cout << "  Accepted objects:" << endl;

      triggerEventWithRefsHandle_->getObjects(filterIndex,electronIds_,electronRefs_);
      const unsigned int nElectrons(electronIds_.size());
      if (nElectrons>0) {
  	cout << "   Electrons: " << nElectrons << "  - the objects: # id pt eta phi vz id key" << endl;
  	for (unsigned int i=0; i!=nElectrons; ++i) {
  	  cout << "   " << i << " " << electronIds_[i]
  	       << " " << electronRefs_[i]->pt()
  	       << " " << electronRefs_[i]->eta()
  	       << " " << electronRefs_[i]->phi()
  	       << " " << electronRefs_[i]->vz()
  	       << " " << electronRefs_[i].id()
  	       << " " << electronRefs_[i].key()
  	       << endl;
  	}
      }

      triggerEventWithRefsHandle_->getObjects(filterIndex,muonIds_,muonRefs_);
      const unsigned int nMuons(muonIds_.size());
      if (nMuons>0) {

	bool isovals = false;
	//	Handle<edm::ValueMap<reco::IsoDeposit> > depMap;
	Handle<edm::ValueMap<float> > depMap;

	//	if (moduleLabel == "hltDiMuonGlb17Trk8DzFiltered0p2TrkIsoFiltered0p4") {
	if (moduleLabel == "hltDiMuonGlb17Trk8DzFiltered0p2") {
	  cout << "found the module before iso! hurray!" << endl;
	  iEvent.getByLabel(isoDepMapTag_,depMap);
	  isovals = true;
	}

  	cout << "   Muons: " << nMuons << "  - the objects: # id pt eta phi vz id key" << endl;
  	for (unsigned int i=0; i!=nMuons; ++i) {
  	  cout << "   " << i << " " << muonIds_[i]
  	       << " " << muonRefs_[i]->pt()
  	       << " " << muonRefs_[i]->eta()
  	       << " " << muonRefs_[i]->phi()
  	       << " " << muonRefs_[i]->vz()
  	       << " " << muonRefs_[i].id()
  	       << " " << muonRefs_[i].key();
	  if (isovals) {
	    const edm::ValueMap<float> ::value_type & muonDeposit = (*(depMap))[muonRefs_[i]];
	    cout << ", isoval: " << muonDeposit << endl;
	    // const edm::ValueMap<reco::IsoDeposit> ::value_type & muonDeposit = (*(depMap))[muonRefs_[i]];
	    // cout << ", isoval: " << muonDeposit.sumWithin(0.3) 
	    // 	 << ", muonDeposit: " << muonDeposit.print() << endl;
	  }
  	  cout << endl;
  	}
      }

  //     triggerEventWithRefsHandle_->getObjects(filterIndex,jetIds_,jetRefs_);
  //     const unsigned int nJets(jetIds_.size());
  //     if (nJets>0) {
  // 	cout << "   Jets: " << nJets << "  - the objects: # id pt" << endl;
  // 	for (unsigned int i=0; i!=nJets; ++i) {
  // 	  cout << "   " << i << " " << jetIds_[i]
  // 	       << " " << jetRefs_[i]->pt()
  // 	       << endl;
  // 	}
  //     }

  //     triggerEventWithRefsHandle_->getObjects(filterIndex,pixtrackIds_,pixtrackRefs_);
  //     const unsigned int nPixTracks(pixtrackIds_.size());
  //     if (nPixTracks>0) {
  // 	cout << "   PixTracks: " << nPixTracks << "  - the objects: # id pt" << endl;
  // 	for (unsigned int i=0; i!=nPixTracks; ++i) {
  // 	  cout << "   " << i << " " << pixtrackIds_[i]
  // 	       << " " << pixtrackRefs_[i]->pt()
  // 	       << endl;
  // 	}
  //     }

  //     triggerEventWithRefsHandle_->getObjects(filterIndex,l1emIds_,l1emRefs_);
  //     const unsigned int nL1EM(l1emIds_.size());
  //     if (nL1EM>0) {
  // 	cout << "   L1EM: " << nL1EM << "  - the objects: # id pt" << endl;
  // 	for (unsigned int i=0; i!=nL1EM; ++i) {
  // 	  cout << "   " << i << " " << l1emIds_[i]
  // 	       << " " << l1emRefs_[i]->pt()
  // 	       << endl;
  // 	}
  //     }

  //     triggerEventWithRefsHandle_->getObjects(filterIndex,l1muonIds_,l1muonRefs_);
  //     const unsigned int nL1Muon(l1muonIds_.size());
  //     if (nL1Muon>0) {
  // 	cout << "   L1Muon: " << nL1Muon << "  - the objects: # id pt" << endl;
  // 	for (unsigned int i=0; i!=nL1Muon; ++i) {
  // 	  cout << "   " << i << " " << l1muonIds_[i]
  // 	       << " " << l1muonRefs_[i]->pt()
  // 	       << endl;
  // 	}
  //     }

  //     triggerEventWithRefsHandle_->getObjects(filterIndex,l1jetIds_,l1jetRefs_);
  //     const unsigned int nL1Jet(l1jetIds_.size());
  //     if (nL1Jet>0) {
  // 	cout << "   L1Jet: " << nL1Jet << "  - the objects: # id pt" << endl;
  // 	for (unsigned int i=0; i!=nL1Jet; ++i) {
  // 	  cout << "   " << i << " " << l1jetIds_[i]
  // 	       << " " << l1jetRefs_[i]->pt()
  // 	       << endl;
  // 	}
  //     }

  //     triggerEventWithRefsHandle_->getObjects(filterIndex,pfjetIds_,pfjetRefs_);
  //     const unsigned int nPFJets(pfjetIds_.size());
  //     if (nPFJets>0) {
  // 	cout << "   PFJets: " << nPFJets << "  - the objects: # id pt" << endl;
  // 	for (unsigned int i=0; i!=nPFJets; ++i) {
  // 	  cout << "   " << i << " " << pfjetIds_[i]
  // 	       << " " << pfjetRefs_[i]->pt()
  // 	       << endl;
  // 	}
  //     }

    }
  }

  return (wasRun && accept && !error);
}
