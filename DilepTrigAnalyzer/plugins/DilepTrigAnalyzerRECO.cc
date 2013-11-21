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
// #include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
// #include "DataFormats/JetReco/interface/CaloJet.h"
// #include "DataFormats/Candidate/interface/CompositeCandidate.h"
// #include "DataFormats/METReco/interface/MET.h"
// #include "DataFormats/METReco/interface/CaloMET.h"
// #include "DataFormats/HcalIsolatedTrack/interface/IsolatedPixelTrackCandidate.h"
// #include "DataFormats/L1Trigger/interface/L1HFRings.h"
// #include "DataFormats/L1Trigger/interface/L1EmParticle.h"
// #include "DataFormats/L1Trigger/interface/L1JetParticle.h"
// #include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
// #include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
// #include "DataFormats/JetReco/interface/PFJet.h"
// #include "DataFormats/TauReco/interface/PFTau.h"

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
  vtxInputTag_(ps.getParameter<edm::InputTag>("vtxInputTag")),
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
       << "   VtxInputTag = " << vtxInputTag_.encode() << endl
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
      bookHists(fs,hltShortNames_.at(itrig)+"_tight");
      bookHists(fs,hltShortNames_.at(itrig)+"_trigiso");
      bookHists(fs,hltShortNames_.at(itrig)+"_tiso");
    }
  } else if (triggerName_ == "emu") {
    h_results_em_ = fs->make<TH1F>("h_results_em" , ";Trigger Results" , 16 , -0.5 , 15.5 );
    for (unsigned int itrig=(unsigned int)em; itrig<=(unsigned int)mei; ++itrig) {
      bookHists(fs,hltShortNames_.at(itrig));
      bookHists(fs,hltShortNames_.at(itrig)+"_tight");
      bookHists(fs,hltShortNames_.at(itrig)+"_trigiso");
      bookHists(fs,hltShortNames_.at(itrig)+"_tiso");
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

  // retrieve necessary containers
  iEvent.getByLabel(vtxInputTag_, vertexHandle_);
  iEvent.getByLabel(electronsInputTag_, elsHandle_);
  iEvent.getByLabel( muonsInputTag_ , musHandle_ );

  
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

  //------------------------------------
  //  hlt objects
  //------------------------------------

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
  //   reco vertices
  //-------------------------------------

  // find vertex 0 in vertex container
  const VertexCollection* vertexCollection = vertexHandle_.product();
  VertexCollection::const_iterator firstGoodVertex = vertexCollection->end();
  for ( VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); ++vtx ) {
    if (  !vtx->isFake() && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0 ) {
      firstGoodVertex = vtx;
      break;
    }
  } // loop on vertices

  if (firstGoodVertex == vertexCollection->end()) {
    cout << "DilepTrigAnalyzerRECO::analyzeTrigger: didn't find any good offline vertices!! size: " 
	 << vertexCollection->size() << std::endl;
    return true;
  }

  //-------------------------------------
  //   reco electrons 
  //-------------------------------------

  const float dr_trigmatch = 0.2;

  LorentzVector off_el_lead;
  LorentzVector off_el_subl;
  int off_el_lead_idx = -1;
  int off_el_subl_idx = -1;

  // needed???
  //  View<GsfElectron> gsfElColl = *(elsHandle_.product());
  // Handle<GsfElectronCollection> els_coll_h;
  // iEvent.getByLabel(electronsInputTag_, els_coll_h);  

  unsigned int elsIndex = 0;
  GsfElectronCollection::const_iterator els_end = elsHandle_->end();  // Iterator
  for( GsfElectronCollection::const_iterator el = elsHandle_->begin(); el != els_end; ++el, ++elsIndex ) {
    LorentzVector lv(el->p4());

    // check for match to trigger objects
    bool match = false;
    if ( (hlt_el_lead_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(lv,hlt_el_lead) < dr_trigmatch) ) match = true;
    else if ( (hlt_el_subl_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(lv,hlt_el_subl) < dr_trigmatch) ) match = true;
    if (!match) continue;

    // pt ordering
    if ( (off_el_lead_idx == -1) || (lv.pt() > off_el_lead.pt()) ) {
      off_el_subl_idx = off_el_lead_idx;
      off_el_subl = off_el_lead;
      off_el_lead_idx = (int)elsIndex;
      off_el_lead = lv;
    } else if ( (off_el_subl_idx == -1) || (lv.pt() > off_el_subl.pt()) ) {
      // check dR to remove exact duplicates
      if (ROOT::Math::VectorUtil::DeltaR(off_el_lead,lv) > 0.001) {
	off_el_subl_idx = (int)elsIndex;
	off_el_subl = lv;
      }
    }

  } // loop on reco gsf electrons

  //-------------------------------------
  //   reco muons 
  //-------------------------------------

  LorentzVector off_mu_lead;
  LorentzVector off_mu_subl;
  int off_mu_lead_idx = -1;
  int off_mu_subl_idx = -1;

  LorentzVector off_mu_tight_lead;
  LorentzVector off_mu_tight_subl;
  int off_mu_tight_lead_idx = -1;
  int off_mu_tight_subl_idx = -1;

  LorentzVector off_mu_trigiso_lead;
  LorentzVector off_mu_trigiso_subl;
  int off_mu_trigiso_lead_idx = -1;
  int off_mu_trigiso_subl_idx = -1;

  LorentzVector off_mu_tiso_lead;
  LorentzVector off_mu_tiso_subl;
  int off_mu_tiso_lead_idx = -1;
  int off_mu_tiso_subl_idx = -1;

  MuonCollection muons_good;
  MuonCollection muons_dup;

  // loop first to remove duplicate muons..
  unsigned int muonIndex = 0;
  MuonCollection::const_iterator muons_end = musHandle_->end();  // Iterator
  for ( MuonCollection::const_iterator muon = musHandle_->begin(); muon != muons_end; ++muon, ++muonIndex ) {
    LorentzVector lv(muon->p4());

    bool duplicate = false;
    // check if this muon is already flagged as a duplicate
    for ( MuonCollection::const_iterator muon2 = muons_dup.begin(); muon2 != muons_dup.end(); ++muon2 ) {
      // !!!!!!!! this may not work
      if (muon == muon2) {
	duplicate = true;
	break;
      }
    } // loop on found duplicates
    if (duplicate) continue;

    // require match to trigger objects
    bool match = false;
    if ( (hlt_mu_lead_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(lv,hlt_mu_lead) < dr_trigmatch) ) match = true;
    else if ( (hlt_mu_subl_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(lv,hlt_mu_subl) < dr_trigmatch) ) match = true;
    if (!match) continue;

    // check tight muon ID, to break ambiguity
    bool pass_tight_1 = muon::isTightMuon(*muon,*firstGoodVertex);

    // // check for duplicate muons using dR 0.1
    duplicate = false;
    for ( MuonCollection::const_iterator muon2 = muon+1; muon2 != muons_end; ++muon2 ) {
      LorentzVector lv2(muon2->p4());
      if (ROOT::Math::VectorUtil::DeltaR(lv,lv2) < 0.1) {
	// check whether each passes tight ID
	bool pass_tight_2 = muon::isTightMuon(*muon2,*firstGoodVertex);
	// mu1 passes, mu2 doesn't: keep mu1
	if (pass_tight_1 && !pass_tight_2) {
	  muons_dup.push_back(*muon2);
	  continue;
	} 
	// mu2 passes, mu1 doesn't: keep mu2
	else if (!pass_tight_1 && pass_tight_2) {
	  muons_dup.push_back(*muon);
	  duplicate = true;
	  break;
	}
	// if both or neither pass tight, take the highest pt muon
	else {
	  if (lv.pt() > lv2.pt()) {
	    muons_dup.push_back(*muon2);
	    continue;
	  } else {
	    muons_dup.push_back(*muon);
	    duplicate = true;
	    break;
	  }
	}
      } // if dR < 0.1 
    } // loop over mu2
    if (duplicate) continue;

    muons_good.push_back(*muon);

  } // loop on reco muons


  // loop again on duplicate-cleaned muon collection


  muonIndex = 0;
  muons_end = muons_good.end();  // Iterator
  for ( MuonCollection::const_iterator muon = muons_good.begin(); muon != muons_end; ++muon, ++muonIndex ) {
    LorentzVector lv(muon->p4());

    // pt ordering
    if ( (off_mu_lead_idx == -1) || (lv.pt() > off_mu_lead.pt()) ) {
      off_mu_subl_idx = off_mu_lead_idx;
      off_mu_subl = off_mu_lead;
      off_mu_lead_idx = (int)muonIndex;
      off_mu_lead = lv;
    } else if ( (off_mu_subl_idx == -1) || (lv.pt() > off_mu_subl.pt()) ) {
      off_mu_subl_idx = (int)muonIndex;
      off_mu_subl = lv;
    }

    bool pass_tight = muon::isTightMuon(*muon,*firstGoodVertex);
    float trkiso = muon->pfIsolationR03().sumChargedHadronPt;
    bool pass_trkiso_trig = bool(trkiso/lv.pt() < 0.4);
    float reliso = muonPFiso(*muon);
    bool pass_iso_tight = bool(reliso < 0.15);

    if (pass_tight) {
      // pt ordering
      if ( (off_mu_tight_lead_idx == -1) || (lv.pt() > off_mu_tight_lead.pt()) ) {
	off_mu_tight_subl_idx = off_mu_tight_lead_idx;
	off_mu_tight_subl = off_mu_tight_lead;
	off_mu_tight_lead_idx = (int)muonIndex;
	off_mu_tight_lead = lv;
      } else if ( (off_mu_tight_subl_idx == -1) || (lv.pt() > off_mu_tight_subl.pt()) ) {
	off_mu_tight_subl_idx = (int)muonIndex;
	off_mu_tight_subl = lv;
      }
    }

    if (pass_tight && pass_trkiso_trig) {
      // pt ordering
      if ( (off_mu_trigiso_lead_idx == -1) || (lv.pt() > off_mu_trigiso_lead.pt()) ) {
	off_mu_trigiso_subl_idx = off_mu_trigiso_lead_idx;
	off_mu_trigiso_subl = off_mu_trigiso_lead;
	off_mu_trigiso_lead_idx = (int)muonIndex;
	off_mu_trigiso_lead = lv;
      } else if ( (off_mu_trigiso_subl_idx == -1) || (lv.pt() > off_mu_trigiso_subl.pt()) ) {
	off_mu_trigiso_subl_idx = (int)muonIndex;
	off_mu_trigiso_subl = lv;
      }
    }

    if (pass_tight && pass_iso_tight) {
      // pt ordering
      if ( (off_mu_tiso_lead_idx == -1) || (lv.pt() > off_mu_tiso_lead.pt()) ) {
	off_mu_tiso_subl_idx = off_mu_tiso_lead_idx;
	off_mu_tiso_subl = off_mu_tiso_lead;
	off_mu_tiso_lead_idx = (int)muonIndex;
	off_mu_tiso_lead = lv;
      } else if ( (off_mu_tiso_subl_idx == -1) || (lv.pt() > off_mu_tiso_subl.pt()) ) {
	off_mu_tiso_subl_idx = (int)muonIndex;
	off_mu_tiso_subl = lv;
      }
    }

  } // loop on duplicate-cleaned muons


  //-------------------------------------
  //   reco lepton plots
  //-------------------------------------

  // make plots based on trigger path
  if (ismm) {
    fillHists(off_mu_lead,off_mu_subl,triggerShort,false);
    fillMuonIsoHists(muons_good,off_mu_lead_idx,off_mu_subl_idx,triggerShort);
    if (off_mu_tight_lead_idx >= 0) {
      fillHists(off_mu_tight_lead,off_mu_tight_subl,triggerShort+"_tight",false);
      fillMuonIsoHists(muons_good,off_mu_tight_lead_idx,off_mu_tight_subl_idx,triggerShort);
    }
    if (off_mu_trigiso_lead_idx >= 0) {
      fillHists(off_mu_trigiso_lead,off_mu_trigiso_subl,triggerShort+"_trigiso",false);
      fillMuonIsoHists(muons_good,off_mu_trigiso_lead_idx,off_mu_trigiso_subl_idx,triggerShort);
    }
    if (off_mu_tiso_lead_idx >= 0) {
      fillHists(off_mu_tiso_lead,off_mu_tiso_subl,triggerShort+"_tiso",false);
      fillMuonIsoHists(muons_good,off_mu_tiso_lead_idx,off_mu_tiso_subl_idx,triggerShort);
    }
  } else if (isem) {
    fillHists(off_el_lead,off_mu_lead,triggerShort,false);
    fillMuonIsoHists(muons_good,off_mu_lead_idx,-1,triggerShort);
    if (off_mu_tight_lead_idx >= 0) {
      fillHists(off_el_lead,off_mu_tight_lead,triggerShort+"_tight",false);
      fillMuonIsoHists(muons_good,-1,off_mu_tight_lead_idx,triggerShort);
    }
    if (off_mu_trigiso_lead_idx >= 0) {
      fillHists(off_el_lead,off_mu_trigiso_lead,triggerShort+"_trigiso",false);
      fillMuonIsoHists(muons_good,-1,off_mu_trigiso_lead_idx,triggerShort);
    }
    if (off_mu_tiso_lead_idx >= 0) {
      fillHists(off_el_lead,off_mu_tiso_lead,triggerShort+"_tiso",false);
      fillMuonIsoHists(muons_good,-1,off_mu_tiso_lead_idx,triggerShort);
    }
  } else if (isme) {
    fillHists(off_mu_lead,off_el_lead,triggerShort,false);
    fillMuonIsoHists(muons_good,off_mu_lead_idx,-1,triggerShort);
    if (off_mu_tight_lead_idx >= 0) {
      fillHists(off_mu_tight_lead,off_el_lead,triggerShort+"_tight",false);
      fillMuonIsoHists(muons_good,off_mu_tight_lead_idx,-1,triggerShort);
    }
    if (off_mu_trigiso_lead_idx >= 0) {
      fillHists(off_mu_trigiso_lead,off_el_lead,triggerShort+"_trigiso",false);
      fillMuonIsoHists(muons_good,off_mu_trigiso_lead_idx,-1,triggerShort);
    }
    if (off_mu_tiso_lead_idx >= 0) {
      fillHists(off_mu_tiso_lead,off_el_lead,triggerShort+"_tiso",false);
      fillMuonIsoHists(muons_good,off_mu_tiso_lead_idx,-1,triggerShort);
    }
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

  hists_1d_["h_lead_reltrkiso"+suf] = fs->make<TH1F>(Form("h_lead_reltrkiso%s",suf.c_str()) , "; Leading trkiso / p_{T}" , 200 , 0. , 2. );
  hists_1d_["h_subl_reltrkiso"+suf] = fs->make<TH1F>(Form("h_subl_reltrkiso%s",suf.c_str()) , "; Subleading trkiso / p_{T}" , 200 , 0. , 2. );
  hists_1d_["h_lead_relpfiso"+suf] = fs->make<TH1F>(Form("h_lead_relpfiso%s",suf.c_str()) , "; Leading pfiso / p_{T}" , 200 , 0. , 2. );
  hists_1d_["h_subl_relpfiso"+suf] = fs->make<TH1F>(Form("h_subl_relpfiso%s",suf.c_str()) , "; Subleading pfiso / p_{T}" , 200 , 0. , 2. );

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

//____________________________________________________________________________
void DilepTrigAnalyzerRECO::fillMuonIsoHists(const reco::MuonCollection& col, const int& lead_idx, const int& subl_idx, const std::string& suffix) {

  if (col.size() == 0) {
    std::cout << "DilepTrigAnalyzerRECO::fillMuonIsoHists: MuonCollection has size 0" << std::endl;
    return;
  }

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  if (lead_idx >= 0) {
    float reltrkiso = col.at(lead_idx).pfIsolationR03().sumChargedHadronPt/col.at(lead_idx).pt();
    hists_1d_["h_lead_reltrkiso"+suf]->Fill(reltrkiso);
    hists_1d_["h_lead_reltrkiso"+suf]->Fill(muonPFiso(col.at(lead_idx)));
  }
  if (subl_idx >= 0) {
    float reltrkiso = col.at(subl_idx).pfIsolationR03().sumChargedHadronPt/col.at(subl_idx).pt();
    hists_1d_["h_subl_reltrkiso"+suf]->Fill(reltrkiso);
    hists_1d_["h_subl_reltrkiso"+suf]->Fill(muonPFiso(col.at(subl_idx)));
  }

  return;
}

//____________________________________________________________________________
float DilepTrigAnalyzerRECO::muonPFiso(const reco::Muon& muon) {

  reco::MuonPFIsolation pfStructR03 = muon.pfIsolationR03();
  float chiso = pfStructR03.sumChargedHadronPt;
  float nhiso = pfStructR03.sumNeutralHadronEt;
  float emiso = pfStructR03.sumPhotonEt;
  float deltaBeta = pfStructR03.sumPUPt;
  float pt = muon.pt();

  //  float absiso = chiso + nhiso + emiso;
  float absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return (absiso / pt);

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
