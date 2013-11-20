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
#include "TrigAnalyzer/DilepTrigAnalyzer/interface/DilepTrigAnalyzerRECO.h"

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

// ROOT includes
#include "Math/VectorUtil.h"

#include <cassert>

//
// constructors and destructor
//
//____________________________________________________________________________
DilepTrigAnalyzerRECO::DilepTrigAnalyzerRECO(const edm::ParameterSet& ps) : 
  processName_(ps.getParameter<std::string>("processName")),
  triggerName_(ps.getParameter<std::string>("triggerName")),
  triggerResultsTag_(ps.getParameter<edm::InputTag>("triggerResults")),
  triggerEventTag_(ps.getParameter<edm::InputTag>("triggerEvent")),
  electronsInputTag_(ps.getParameter<edm::InputTag>("electronsInputTag")),
  muonsInputTag_(ps.getParameter<edm::InputTag>("muonsInputTag")),
  verbose_(ps.getParameter<bool>("verbose"))
{
  using namespace std;
  using namespace edm;

  cout << "DilepTrigAnalyzerRECO configuration: " << endl
       << "   ProcessName = " << processName_ << endl
       << "   TriggerName = " << triggerName_ << endl
       << "   TriggerResultsTag = " << triggerResultsTag_.encode() << endl
       << "   TriggerEventTag = " << triggerEventTag_.encode() << endl
       << "   ElectronsInputTag = " << electronsInputTag_.encode() << endl
       << "   MuonsInputTag = " << muonsInputTag_.encode() << endl
       << "   Verbose = " << verbose_ << endl;

  // hltTriggerNames_.reserve(8);
  // hltTriggerNames_.at(mm) = "HLT_Mu17_Mu8_v23";
  // hltTriggerNames_.at(mmi) = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1";
  // hltTriggerNames_.at(mmtk) = "HLT_Mu17_TkMu8_v15";
  // hltTriggerNames_.at(mmitk) = "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1";
  // hltTriggerNames_.at(em) = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10";
  // hltTriggerNames_.at(emi) = "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1";
  // hltTriggerNames_.at(me) = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10";
  // hltTriggerNames_.at(mei) = "HLT_Mu17_TrkIsoVVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1";

  // hltShortNames_.reserve(8);
  // hltShortNames_.at(mm) = "mm";
  // hltShortNames_.at(mmi) = "mmi";
  // hltShortNames_.at(mmtk) = "mmtk";
  // hltShortNames_.at(mmitk) = "mmitk";
  // hltShortNames_.at(em) = "em";
  // hltShortNames_.at(emi) = "emi";
  // hltShortNames_.at(me) = "me";
  // hltShortNames_.at(mei) = "mei";

  hltTriggerNames_.push_back("HLT_Mu17_Mu8_v23");
  hltTriggerNames_.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1");
  hltTriggerNames_.push_back("HLT_Mu17_TkMu8_v15");
  hltTriggerNames_.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1");
  hltTriggerNames_.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10");
  hltTriggerNames_.push_back("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1");
  hltTriggerNames_.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10");
  hltTriggerNames_.push_back("HLT_Mu17_TrkIsoVVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1");

  hltShortNames_.push_back("mm");
  hltShortNames_.push_back("mmi");
  hltShortNames_.push_back("mmtk");
  hltShortNames_.push_back("mmitk");
  hltShortNames_.push_back("em");
  hltShortNames_.push_back("emi");
  hltShortNames_.push_back("me");
  hltShortNames_.push_back("mei");

  // histogram setup
  edm::Service<TFileService> fs;
  if (triggerName_ == "dimu") {
    h_results_mm_ = fs->make<TH1F>("h_results_mm" , ";Trigger Results" , 16 , -0.5 , 15.5 );
    for (unsigned int itrig=0; itrig<=(unsigned int)mmitk; ++itrig) {
      bookHists(fs,hltShortNames_.at(itrig));
    }
  } else if (triggerName_ == "emu") {
    h_results_em_ = fs->make<TH1F>("h_results_em" , ";Trigger Results" , 16 , -0.5 , 15.5 );
    for (unsigned int itrig=(unsigned int)em; itrig<=(unsigned int)mei; ++itrig) {
      bookHists(fs,hltShortNames_.at(itrig));
    }
  }

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
      const unsigned int n(hltConfig_.size());
      // check if trigger names in (new) config
      if (triggerName_=="dimu") {
	for (unsigned int itrig = 0; itrig <= (unsigned int)mmitk; ++itrig) {
	  unsigned int triggerIndex(hltConfig_.triggerIndex(hltTriggerNames_.at(itrig)));
	  if (triggerIndex>=n) {
	    cout << "DilepTrigAnalyzerRECO::analyze:"
		 << " TriggerName " << hltTriggerNames_.at(itrig) 
		 << " not available in (new) config!" << endl;
	  }
	} // loop over dimu triggers
      } // if dimu triggers 
      else if (triggerName_=="emu") {
	for (unsigned int itrig = (unsigned int)em; itrig <= (unsigned int)mei; ++itrig) {
	  unsigned int triggerIndex(hltConfig_.triggerIndex(hltTriggerNames_.at(itrig)));
	  if (triggerIndex>=n) {
	    cout << "DilepTrigAnalyzerRECO::analyze:"
		 << " TriggerName " << hltTriggerNames_.at(itrig) 
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

  if (verbose_) cout << endl;

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
  if (triggerName_=="dimu") {
    unsigned int results(0);
    for (unsigned int itrig=0; itrig <= (unsigned int)mmitk; ++itrig) {
      bool pass = analyzeTrigger(iEvent,iSetup,(hltTrigs)itrig);
      if (pass) results |= 1 << itrig;
    }
    h_results_mm_->Fill(results);
  } else if (triggerName_=="emu") {
    unsigned int results(0);
    for (unsigned int itrig=(unsigned int)em; itrig <= (unsigned int)mei ; ++itrig) {
      bool pass = analyzeTrigger(iEvent,iSetup,(hltTrigs)itrig);
      if (pass) results |= 1 << (itrig - (unsigned int)em);
    }
    h_results_em_->Fill(results);
  } else {
	cout << "DilepTrigAnalyzerRECO::analyze: TriggerName " << triggerName_ 
	     << " not recognized!!" << endl;
      }

  if (verbose_) cout << endl;

  return;
}

//____________________________________________________________________________

bool DilepTrigAnalyzerRECO::analyzeTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup, const hltTrigs triggerEnum) {
  
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  if (verbose_) cout << endl;

  const std::string triggerName = hltTriggerNames_.at(triggerEnum);
  const std::string triggerShort = hltShortNames_.at(triggerEnum);

  const unsigned int ntrigs(hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));

  // abort on invalid trigger name
  if (triggerIndex>=ntrigs) {
    cout << "DilepTrigAnalyzerRECO::analyzeTrigger: path "
	 << triggerName << " - not found!" << endl;
    return false;
  }

  // check what kind of trigger we have
  bool ismm = false;
  bool isem = false;
  bool isme = false;
  if ( triggerEnum <= mmitk ) {
    ismm = true;
  } else if ( triggerEnum <= emi ) {
    isem = true;
  } else if (triggerEnum <= mei) {
    isme = true;
  }

  if (!ismm && !isem && !isme) {
    cout << "DilepTrigAnalyzerRECO::analyzeTrigger: triggerEnum: " << triggerEnum
	 << " not recognized, aborting.." << endl;
    return false;
  }
  
  if (verbose_) {
    cout << "DilepTrigAnalyzerRECO::analyzeTrigger: path "
	 << triggerName << " [" << triggerIndex << "]" << endl;
  }
  // modules on this trigger path
  const unsigned int m(hltConfig_.size(triggerIndex));
  const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));

  bool wasRun = triggerResultsHandle_->wasrun(triggerIndex);
  bool accept = triggerResultsHandle_->accept(triggerIndex);
  bool error = triggerResultsHandle_->error(triggerIndex);
  const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
  // Results from TriggerResults product
  if (verbose_) {
    cout << " Trigger path status:"
	 << " WasRun=" << wasRun
	 << " Accept=" << accept
	 << " Error =" << error
	 << endl;
    cout << " Last active module - label/type: "
	 << moduleLabels[moduleIndex] << "/" << hltConfig_.moduleType(moduleLabels[moduleIndex])
	 << " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]"
	 << endl;
  }
  assert (moduleIndex<m);

  if (!wasRun || !accept || error) return false;

  // loop over trigger and reco objects, match, make plots

  // first, get trigger objects from last filter

  // Results from TriggerEvent product - Attention: must look only for
  // modules actually run in this path for this event!
  const string& moduleLabel(moduleLabels[moduleIndex-1]);
  const string  moduleType(hltConfig_.moduleType(moduleLabel));
  // check whether the module is packed up in TriggerEvent product
  const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(moduleLabel,"",processName_)));
  if (!(filterIndex<triggerEventHandle_->sizeFilters())) {
    cout << "DilepTrigAnalyzerRECO::analyzeTrigger: filterIndex out of range! triggerName: " << triggerName
	 << ", filterIndex: " << filterIndex << ", sizeFilters(): " << triggerEventHandle_->sizeFilters() << std::endl;
    return true;
  }

  if (verbose_) {
    cout << " 'L3' filter in slot " << moduleIndex-1 << " - label/type " << moduleLabel << "/" << moduleType << endl;
  }
  const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
  const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
  const size_type nI(VIDS.size());
  const size_type nK(KEYS.size());
  assert(nI==nK);
  const size_type n(max(nI,nK));
  if (verbose_) {
    cout << "   " << n  << " accepted 'L3' objects found: " << endl;
  }

  LorentzVector hlt_el_lead;
  LorentzVector hlt_el_subl;
  int hlt_el_lead_idx = -1;
  int hlt_el_subl_idx = -1;

  LorentzVector hlt_mu_lead;
  LorentzVector hlt_mu_subl;
  int hlt_mu_lead_idx = -1;
  int hlt_mu_subl_idx = -1;

  const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
  for (size_type i=0; i!=n; ++i) {
    const TriggerObject& TO(TOC[KEYS[i]]);
    LorentzVector lv( TO.particle().p4() );
    // electrons
    if (VIDS[i] == 82) {
      if ( (hlt_el_lead_idx == -1) || (lv.pt() > hlt_el_lead.pt()) ) {
	hlt_el_subl_idx = hlt_el_lead_idx;
	hlt_el_subl = hlt_el_lead;
	hlt_el_lead_idx = (int)i;
	hlt_el_lead = lv;
      } else if ( (hlt_el_subl_idx == -1) || (lv.pt() > hlt_el_subl.pt()) ) {
	// check dR to remove exact duplicates
	if (ROOT::Math::VectorUtil::DeltaR(hlt_el_lead,lv) > 0.001) {
	  hlt_el_subl_idx = (int)i;
	  hlt_el_subl = lv;
	}
      }
    }
    // muons
    else if (VIDS[i] == 83) {
      if ( (hlt_mu_lead_idx == -1) || (lv.pt() > hlt_mu_lead.pt()) ) {
	hlt_mu_subl_idx = hlt_mu_lead_idx;
	hlt_mu_subl = hlt_mu_lead;
	hlt_mu_lead_idx = (int)i;
	hlt_mu_lead = lv;
      } else if ( (hlt_mu_subl_idx == -1) || (lv.pt() > hlt_mu_subl.pt()) ) {
	if (ROOT::Math::VectorUtil::DeltaR(hlt_mu_lead,lv) > 0.001) {
	  hlt_mu_subl_idx = (int)i;
	  hlt_mu_subl = lv;
	}
      }
    }

    // just need to apply cuts to trigger objects to find correct ones
    if (verbose_) {
      cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
	   << TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()
	   << endl;
    }
  }

  if ((hlt_el_lead_idx == -1) && (hlt_mu_lead_idx == -1)) {
    cout << "DilepTrigAnalyzerRECO::analyzeTrigger: no valid trigger leptons! n(obj) = " << n << endl;
    return true;
  }

  // make plots based on trigger path
  if (ismm) {
    fillHists(hlt_mu_lead,hlt_mu_subl,triggerShort,true);
  } else if (isem) {
    fillHists(hlt_el_lead,hlt_mu_lead,triggerShort,true);
  } else if (isme) {
    fillHists(hlt_mu_lead,hlt_el_lead,triggerShort,true);
  }

  // now loop on reco electrons and muons

  //-------------------------------------
  //   reco electrons 
  //-------------------------------------

  Handle<View<GsfElectron> > els_h;
  iEvent.getByLabel(electronsInputTag_, els_h);

  // needed???
  //  View<GsfElectron> gsfElColl = *(els_h.product());
  // Handle<GsfElectronCollection> els_coll_h;
  // iEvent.getByLabel(electronsInputTag_, els_coll_h);  

  unsigned int elsIndex = 0;
  View<GsfElectron>::const_iterator els_end = els_h->end();  // Iterator
  for( View<GsfElectron>::const_iterator el = els_h->begin(); el != els_end; ++el, ++elsIndex ) {

  }

  //-------------------------------------
  //   reco muons 
  //-------------------------------------

  Handle<View<Muon> > muon_h;
  iEvent.getByLabel( muonsInputTag_ , muon_h );

  unsigned int muonIndex = 0;
  View<Muon>::const_iterator muons_end = muon_h->end();  // Iterator
  for ( View<Muon>::const_iterator muon = muon_h->begin(); muon != muons_end; ++muon, ++muonIndex ) {

  }


  return true;
}

//____________________________________________________________________________
void DilepTrigAnalyzerRECO::bookHists(edm::Service<TFileService>& fs, const std::string& suffix) {

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  std::string hlt_suf("_hlt");

  hists_1d_["h_lead_pt"+suf+hlt_suf] = fs->make<TH1F>(Form("h_lead_pt%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Leading p_{T} [GeV]" , 100 , 0. , 100. );
  hists_1d_["h_subl_pt"+suf+hlt_suf] = fs->make<TH1F>(Form("h_subl_pt%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Subleading p_{T} [GeV]" , 100 , 0. , 100. );
  hists_1d_["h_lead_eta"+suf+hlt_suf] = fs->make<TH1F>(Form("h_lead_eta%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Leading #eta" , 100 , -3. , 3. );
  hists_1d_["h_subl_eta"+suf+hlt_suf] = fs->make<TH1F>(Form("h_subl_eta%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Subleading #eta" , 100 , -3. , 3. );
  hists_1d_["h_mll"+suf+hlt_suf] = fs->make<TH1F>(Form("h_mll%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT M_{ll} [GeV]" , 150 , 0. , 150. );
  hists_1d_["h_dr"+suf+hlt_suf] = fs->make<TH1F>(Form("h_dr%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT #DeltaR" , 600 , 0. , 6. );

  hists_1d_["h_lead_pt"+suf] = fs->make<TH1F>(Form("h_lead_pt%s",suf.c_str()) , "; Leading p_{T} [GeV]" , 100 , 0. , 100. );
  hists_1d_["h_subl_pt"+suf] = fs->make<TH1F>(Form("h_subl_pt%s",suf.c_str()) , "; Subleading p_{T} [GeV]" , 100 , 0. , 100. );
  hists_1d_["h_lead_eta"+suf] = fs->make<TH1F>(Form("h_lead_eta%s",suf.c_str()) , "; Leading #eta" , 100 , -3. , 3. );
  hists_1d_["h_subl_eta"+suf] = fs->make<TH1F>(Form("h_subl_eta%s",suf.c_str()) , "; Subleading #eta" , 100 , -3. , 3. );
  hists_1d_["h_mll"+suf] = fs->make<TH1F>(Form("h_mll%s",suf.c_str()) , "; M_{ll} [GeV]" , 150 , 0. , 150. );
  hists_1d_["h_dr"+suf] = fs->make<TH1F>(Form("h_dr%s",suf.c_str()) , "; #DeltaR" , 600 , 0. , 6. );

  return;
}

//____________________________________________________________________________
void DilepTrigAnalyzerRECO::fillHists(const LorentzVector& lead, const LorentzVector& subl, const std::string& suffix, bool hlt) {

  if (lead.pt() <= 0.) {
    std::cout << "DilepTrigAnalyzerRECO::fillHists: invalid lead pt: " << lead.pt() << std::endl;
    return;
  }

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  std::string hlt_suf("");
  if (hlt) hlt_suf = "_hlt";

  LorentzVector dilep = lead+subl;

  hists_1d_["h_lead_pt"+suf+hlt_suf]->Fill(lead.pt());
  hists_1d_["h_lead_eta"+suf+hlt_suf]->Fill(lead.eta());
  if (subl.pt() > 0) {
    hists_1d_["h_subl_pt"+suf+hlt_suf]->Fill(subl.pt());
    hists_1d_["h_subl_eta"+suf+hlt_suf]->Fill(subl.eta());
    hists_1d_["h_mll"+suf+hlt_suf]->Fill(dilep.M());
    hists_1d_["h_dr"+suf+hlt_suf]->Fill(ROOT::Math::VectorUtil::DeltaR(lead,subl));
  }

  return;
}

// //____________________________________________________________________________
// bool DilepTrigAnalyzerRECO::isTriggerPassed(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& triggerName) {

//   const unsigned int n(hltConfig_.size());
//   const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
//   assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));

//   // abort on invalid trigger name
//   if (triggerIndex>=n) {
//     std::cout << "DilepTrigAnalyzerRECO::isTriggerPassed: path "
// 	      << triggerName << " - not found!" << std::endl;
//     return false;
//   }
  
//   bool wasRun = triggerResultsHandle_->wasrun(triggerIndex);
//   bool accept = triggerResultsHandle_->accept(triggerIndex);
//   bool error = triggerResultsHandle_->error(triggerIndex);

//   return (wasRun && accept && !error);
// }

// //____________________________________________________________________________
// std::string DilepTrigAnalyzerRECO::lastModuleLabel(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& triggerName) {

//   const unsigned int n(hltConfig_.size());
//   const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
//   assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));

//   // abort on invalid trigger name
//   if (triggerIndex>=n) {
//     std::cout << "DilepTrigAnalyzerRECO::lastModuleLabel: path "
// 	      << triggerName << " - not found!" << std::endl;
//     return false;
//   }
  
//   // modules on this trigger path
//   const unsigned int m(hltConfig_.size(triggerIndex));
//   const std::vector<std::string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));

//   const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
//   assert (moduleIndex<m);

//   return moduleLabels[moduleIndex];
// }
