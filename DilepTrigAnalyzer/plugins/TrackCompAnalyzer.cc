/** \class TrackCompAnalyzer
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
#include "TrigAnalyzer/DilepTrigAnalyzer/interface/TrackCompAnalyzer.h"

#include "DataFormats/Common/interface/ValueMap.h"
// #include "DataFormats/VertexReco/interface/Vertex.h"
// #include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
// #include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
// #include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

// ROOT includes
#include "Math/VectorUtil.h"

#include <cassert>

//
// constructors and destructor
//
//____________________________________________________________________________
TrackCompAnalyzer::TrackCompAnalyzer(const edm::ParameterSet& ps) : 
  processName_(ps.getParameter<std::string>("processName")),
  triggerName_(ps.getParameter<std::string>("triggerName")),
  triggerResultsTag_(ps.getParameter<edm::InputTag>("triggerResults")),
  triggerEventTag_(ps.getParameter<edm::InputTag>("triggerEvent")),
  triggerEventWithRefsTag_(ps.getParameter<edm::InputTag>("triggerEventWithRefs")),
  hltTracksInputTag_(ps.getParameter<edm::InputTag>("hltTracksInputTag")),
  offlineTracksInputTag_(ps.getParameter<edm::InputTag>("offlineTracksInputTag")),
  offlinePFCandsInputTag_(ps.getParameter<edm::InputTag>("offlinePFCandsInputTag")),
  offlineVerticesInputTag_(ps.getParameter<edm::InputTag>("offlineVerticesInputTag")),
  muonsInputTag_(ps.getParameter<edm::InputTag>("muonsInputTag")),
  beamSpotInputTag_(ps.getParameter<edm::InputTag>("beamSpotInputTag")),
  verbose_(ps.getParameter<bool>("verbose"))
{
  using namespace std;
  using namespace edm;

  cout << "TrackCompAnalyzer configuration: " << endl
       << "   processName = " << processName_ << endl
       << "   triggerName = " << triggerName_ << endl
       << "   triggerResultsTag = " << triggerResultsTag_.encode() << endl
       << "   triggerEventTag = " << triggerEventTag_.encode() << endl
       << "   triggerEventWithRefsTag = " << triggerEventWithRefsTag_.encode() << endl
       << "   hltTracksInputTag = " << hltTracksInputTag_.encode() << endl
       << "   offlineTracksInputTag = " << offlineTracksInputTag_.encode() << endl
       << "   offlinePFCandsInputTag = " << offlinePFCandsInputTag_.encode() << endl
       << "   offlineVerticesInputTag = " << offlineVerticesInputTag_.encode() << endl
       << "   muonsInputTag = " << muonsInputTag_.encode() << endl
       << "   beamSpotInputTag = " << beamSpotInputTag_.encode() << endl
       << "   verbose = " << verbose_ << endl;

  // hltTriggerNames_.push_back("HLT_Mu17_Mu8_v23");
  // hltTriggerNames_.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1");
  // hltTriggerNames_.push_back("HLT_Mu17_TkMu8_v15");
  // hltTriggerNames_.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1");
  // hltTriggerNames_.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10");
  // hltTriggerNames_.push_back("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1");
  // hltTriggerNames_.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10");
  // hltTriggerNames_.push_back("HLT_Mu17_TrkIsoVVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1");

  // histogram setup
  edm::Service<TFileService> fs;
  bookHists(fs,"offall_iter0");
  bookHists(fs,"offreg_iter0");
  bookHists(fs,"offonly_iter0");
  bookHists(fs,"offall_others");
  bookHists(fs,"offreg_others");
  bookHists(fs,"offonly_others");

  bookHists(fs,"hltreg");
  bookHistsMuHLT(fs,"hltreg");

  bookHists(fs,"hltmatch");
  bookHistsMuHLT(fs,"hltmatch");
  bookHistsPFHLT(fs,"hltmatch");

  bookHists(fs,"hltpfnonhad");
  bookHistsMuHLT(fs,"hltpfnonhad");
  bookHistsPFHLT(fs,"hltpfnonhad");

  bookHists(fs,"hltpfnonvtx");
  bookHistsMuHLT(fs,"hltpfnonvtx");
  bookHistsPFHLT(fs,"hltpfnonvtx");

  bookHists(fs,"hltnopf");
  bookHistsMuHLT(fs,"hltnopf");
  bookHistsRecoHLT(fs,"hltnopf");

  bookHists(fs,"hltonly");
  bookHistsMuHLT(fs,"hltonly");

  hists_1d_["h_drmin_hltonly"] = fs->make<TH1F>(Form("h_drmin_hltonly") , "; min #DeltaR(off,HLT)" , 600 , 0. , 6. );
  hists_1d_["h_offvtx_hltmatch"] = fs->make<TH1F>(Form("h_offvtx_hltmatch") , "; off vtx" , 32 , -2.5 , 29.5 );
  hists_1d_["h_offvtx_hltpfnonhad"] = fs->make<TH1F>(Form("h_offvtx_hltpfnonhad") , "; off vtx" , 32 , -2.5 , 29.5 );
  hists_1d_["h_offvtx_hltpfnonvtx"] = fs->make<TH1F>(Form("h_offvtx_hltpfnonvtx") , "; off vtx" , 32 , -2.5 , 29.5 );
  //  hists_1d_["h_dphi_offmu"] = fs->make<TH1F>(Form("h_dphi_offmu") , "; #Delta#phi(off,#mu)" , 600 , -2*ROOT::Math::Pi() , 2*ROOT::Math::Pi() );
  bookHistsRecoHLT(fs,"offhlt");

}

//____________________________________________________________________________
TrackCompAnalyzer::~TrackCompAnalyzer()
{
}

//
// member functions
//
//____________________________________________________________________________
void
TrackCompAnalyzer::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  using namespace std;
  using namespace edm;

  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
      const unsigned int n(hltConfig_.size());
      // check if trigger names in (new) config
      unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
      if (triggerIndex>=n) {
	cout << "TrackCompAnalyzer::analyze:"
	     << " TriggerName " << triggerName_ 
	     << " not available in (new) config!" << endl;
      }
    } // if changed
  } else {
    cout << "TrackCompAnalyzer::analyze:"
	 << " config extraction failure with process name "
	 << processName_ << endl;
  }

}

//____________________________________________________________________________
// ------------ method called to produce the data  ------------
void
TrackCompAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  // bool skip = true;
  // // toggle verbosity based on event number
  // if (iEvent.id().event() == 123672589 || iEvent.id().event() == 124630217 || iEvent.id().event() == 123927053 
  //     || iEvent.id().event() == 125711062 || iEvent.id().event() == 126045757 || iEvent.id().event() == 127045570
  //     || iEvent.id().event() == 125252732 || iEvent.id().event() == 125368560 || iEvent.id().event() == 147040431) {
  //   verbose_ = true;
  //   skip = false;
  //   cout << endl << "- Event with HLT inefficiency: run: " << iEvent.id().run() << ", event: " 
  // 	 << iEvent.id().event() << endl;
  // } else {
  //   verbose_ = false;
  // }

  // if (skip) return;

  if (verbose_) cout << endl;

  // get event products
  iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
    cout << "TrackCompAnalyzer::analyze: Error in getting TriggerResults product from Event!" << endl;
    return;
  }
  iEvent.getByLabel(triggerEventTag_,triggerEventHandle_);
  if (!triggerEventHandle_.isValid()) {
    cout << "TrackCompAnalyzer::analyze: Error in getting TriggerEvent product from Event!" << endl;
    return;
  }
  iEvent.getByLabel(triggerEventWithRefsTag_,triggerEventWithRefsHandle_);
  if (!triggerEventWithRefsHandle_.isValid()) {
    cout << "TrackCompAnalyzer::analyze: Error in getting TriggerEventWithRefs product from Event!" << endl;
    return;
  }
  // sanity check
  assert(triggerResultsHandle_->size()==hltConfig_.size());

  // retrieve necessary containers
  iEvent.getByLabel(hltTracksInputTag_, hltTracksHandle_);
  iEvent.getByLabel(offlineTracksInputTag_, offlineTracksHandle_);
  iEvent.getByLabel(offlinePFCandsInputTag_, offlinePFCandsHandle_);
  iEvent.getByLabel(offlineVerticesInputTag_, offlineVerticesHandle_);
  iEvent.getByLabel(muonsInputTag_, musHandle_ );
  iEvent.getByLabel(beamSpotInputTag_,beamSpotHandle_);

  if (verbose_) cout << endl;

  const unsigned int ntrigs(hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName_));

  // abort on invalid trigger name
  if (triggerIndex>=ntrigs) {
    cout << "TrackCompAnalyzer::analyze: path "
	 << triggerName_ << " - not found!" << endl;
    return;
  }

  // check what kind of trigger we have
  bool ismm = false;
  bool ismmtk = false;
  bool isem = false;
  bool isme = false;
  string moduleLabelPreIso;
  if ( triggerName_.find("HLT_Mu17_Mu8") != string::npos ) {
    ismm = true;
    moduleLabelPreIso = "hltDiMuonGlb17Glb8DzFiltered0p2";
  } else if ( triggerName_.find("HLT_Mu17_TkMu8") != string::npos ) {
    ismmtk = true;
    moduleLabelPreIso = "hltDiMuonGlb17Trk8DzFiltered0p2";
  } else if ( triggerName_.find("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL") != string::npos ) {
    isem = true;
    moduleLabelPreIso = "hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8";
  } else if ( triggerName_.find("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL") != string::npos ) {
    isme = true;
    moduleLabelPreIso = "hltL1Mu12EG7L3MuFiltered17";
  }

  if (!ismm && !ismmtk && !isem && !isme) {
    cout << "TrackCompAnalyzer::analyze: triggerName: " << triggerName_
	 << " not recognized, aborting.." << endl;
    return;
  }
  
  if (verbose_) {
    cout << "TrackCompAnalyzer::analyze: path "
	 << triggerName_ << " [" << triggerIndex << "]" << endl;
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

  if (!wasRun || !accept || error) return;

  // loop over trigger and reco objects, match, make plots

  // first, get trigger objects from last filter

  //------------------------------------
  //  hlt objects
  //------------------------------------

  bool foundMuons = false;

  vector<RecoChargedCandidateRef> hlt_mus;

  // Results from TriggerEvent product - Attention: must look only for
  // modules actually run in this path for this event!
  // -- loop backwards through modules to find last filter for each object type
  //  for (unsigned int j=0; j<=moduleIndex; ++j) {
  for (unsigned int j=moduleIndex; j!=0; --j) {
    const string& moduleLabel(moduleLabels[j]);
    // check for the module before the isolation filter
    if (moduleLabel != moduleLabelPreIso) continue;
    const string  moduleType(hltConfig_.moduleType(moduleLabel));
    // check whether the module is packed up in TriggerEvent product
    const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(moduleLabel,"",processName_)));
    //    if (filterIndex>=triggerEventHandle_->sizeFilters()) continue;
    if (filterIndex>=triggerEventWithRefsHandle_->size()) continue;
    if (verbose_) {
      cout << " 'L3' filter in slot " << j << " - label/type " << moduleLabel << "/" << moduleType << endl
	   << " Filter packed up at: " << filterIndex << endl;
    }
    if (moduleLabel == "hltBoolEnd") continue;

    muonIds_.clear();
    muonRefs_.clear();
    triggerEventWithRefsHandle_->getObjects(filterIndex,muonIds_,muonRefs_);
    const unsigned int nMuons(muonIds_.size());
    if (nMuons>0) {

      if (verbose_) {
	cout << "   Muons: " << nMuons << ", MuonRefs: " << muonRefs_.size()
	     << "  - the objects: # id pt eta phi vz id key" << endl;
      }
      for (unsigned int i=0; i!=nMuons; ++i) {

	if (!foundMuons) {
	  LorentzVector lv(muonRefs_.at(i)->p4());

	  bool overlap = false;
	  for (unsigned int k=0; k < hlt_mus.size(); ++k) {
	    if (ROOT::Math::VectorUtil::DeltaR(hlt_mus.at(k)->p4(),lv) < 0.001) {
	      overlap = true;
	      break;
	    }
	  }

	  if (overlap) continue;
	  hlt_mus.push_back(muonRefs_.at(i));

	} // if !foundMuons

	if (verbose_) {
	  cout << "   " << i
	       << " " << muonIds_.at(i)
	       << " " << muonRefs_.at(i)->pt()
	       << " " << muonRefs_.at(i)->eta()
	       << " " << muonRefs_.at(i)->phi()
	       << " " << muonRefs_.at(i)->vz()
	       << " " << muonRefs_.at(i).id()
	       << " " << muonRefs_.at(i).key();
	  cout << endl;
	}

      } // loop over hlt muons
    } // if muons in TriggerEventWithRefs

    // --------- code for triggerEvent, instead of triggerEventWithRefs -------------------

    // const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
    // const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
    // const size_type nI(VIDS.size());
    // const size_type nK(KEYS.size());
    // assert(nI==nK);
    // const size_type n(max(nI,nK));
    // if (verbose_) cout << "   " << n  << " accepted 'L3' objects found: " << endl;
    // const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
    // for (size_type i=0; i!=n; ++i) {
    //   const TriggerObject& TO(TOC[KEYS[i]]);
    //   LorentzVector lv( TO.particle().p4() );

    //   // muons
    //   if (!foundMuons && (VIDS[i] == 83)) {
    // 	if ( (hlt_mu_lead_idx == -1) || (lv.pt() > hlt_mu_lead.pt()) ) {
    // 	  hlt_mu_subl_idx = hlt_mu_lead_idx;
    // 	  hlt_mu_subl = hlt_mu_lead;
    // 	  hlt_mu_lead_idx = (int)i;
    // 	  hlt_mu_lead = lv;
    // 	} else if ( (hlt_mu_subl_idx == -1) || (lv.pt() > hlt_mu_subl.pt()) ) {
    // 	  if (ROOT::Math::VectorUtil::DeltaR(hlt_mu_lead,lv) > 0.001) {
    // 	    hlt_mu_subl_idx = (int)i;
    // 	    hlt_mu_subl = lv;
    // 	  }
    // 	}
    //   } // hlt muons

    //   if (verbose_) {
    // 	cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
    // 	     << TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()
    // 	     << endl;
    //   }
    // } // loop on trig objects

    if (hlt_mus.size() > 0) {
      foundMuons = true;
      //      hlt_mu_filter_idx = filterIndex;
    }

    if (foundMuons) break;
  } // backwards loop on modules

  if (hlt_mus.size() == 0) {
    cout << "TrackCompAnalyzer::analyze: no valid trigger leptons!" << endl;
    return;
  }

  const TrackCollection* hltTracksCollection = hltTracksHandle_.product();
  const TrackCollection* offlineTracksCollection = offlineTracksHandle_.product();
  const PFCandidateCollection* offlinePFCandsCollection = offlinePFCandsHandle_.product();
  const VertexCollection* offlineVerticesCollection = offlineVerticesHandle_.product();
  const TrackCollection::const_iterator hltEnd = hltTracksCollection->end();
  const TrackCollection::const_iterator offlineEnd = offlineTracksCollection->end();
  const PFCandidateCollection::const_iterator pfcandEnd = offlinePFCandsCollection->end();

  if (verbose_) cout << endl;

  //-------------------------------------
  //   reco vertices
  //-------------------------------------

  // find vertex 0 in vertex container
  VertexCollection::const_iterator firstGoodVertex = offlineVerticesCollection->end();
  int vtx_idx = 0;
  for ( VertexCollection::const_iterator vtx = offlineVerticesCollection->begin(); vtx != offlineVerticesCollection->end(); ++vtx ) {
    if (  !vtx->isFake() && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0 ) {
      firstGoodVertex = vtx;
      break;
    }
    ++vtx_idx;
  } // loop on vertices

  if (firstGoodVertex == offlineVerticesCollection->end()) {
    cout << "TrackCompAnalyzer::analyze: didn't find any good offline vertices!! size: " 
  	 << offlineVerticesCollection->size() << endl;
    return;
  } else if (verbose_) {
    cout << " - first good vertex: " << vtx_idx << ", z: " << firstGoodVertex->position().Z() << ", ntracks: " 
	 << firstGoodVertex->tracksSize() << endl;
  }

  if (verbose_) cout << endl;

  // -------------------------------------------------------------------
  // --------- loop on offline muons -----------------------------------
  // -------------------------------------------------------------------

  const float dr_trigmatch = 0.2;

  // loop first to remove duplicate muons..
  MuonCollection::const_iterator muons_end = musHandle_->end();  // Iterator
  for ( MuonCollection::const_iterator muon = musHandle_->begin(); muon != muons_end; ++muon ) {

    // require match to trigger objects
    bool trigmatch = false;
    for (unsigned int ihlt=0; ihlt < hlt_mus.size(); ++ihlt) {
      float dr = ROOT::Math::VectorUtil::DeltaR(muon->p4(),hlt_mus.at(ihlt)->p4());
      if (dr < dr_trigmatch) {
	trigmatch = true;
	break;
      }
    }
    if (!trigmatch) continue;

    // some issue with muon trackRefs sometimes..
    int vtx = -2;
    if (!muon->track().isNull()) vtx = trackVertex(*offlineVerticesCollection,muon->track());
    if (verbose_) {
      cout << " - reco muon, pt: " << muon->pt() << ", eta: " << muon->eta()
	   << ", phi: " << muon->phi() << ", vz: " << muon->vz() 
           << ", vtx: " << vtx << endl;
    }
  }

  if (verbose_) cout << endl;

  // -------------------------------------------------------------------
  // --------- loop on HLT tracks --------------------------------------
  // -------------------------------------------------------------------

  // first loop on HLT tracks and remove muons, make plots
  TrackCollection hltTracksCollectionNoMuons;
  for ( TrackCollection::const_iterator hlt_trk = hltTracksCollection->begin(); hlt_trk != hltEnd; ++hlt_trk ) {

    // check dR from hlt muons: those tracks aren't used in track iso
    bool muon_overlap = false;
    bool muon_iso_region = false;
    float mindr_mu = 99.;
    int mindr_mu_idx = -1;
    for (unsigned int ihlt=0; ihlt < hlt_mus.size(); ++ihlt) {
      float dr = ROOT::Math::VectorUtil::DeltaR(hlt_trk->momentum(),hlt_mus.at(ihlt)->p4());
      float dz = hlt_mus.at(ihlt)->vz() - hlt_trk->vz();
      if (dr < 0.01) {
	muon_overlap = true;
	break;
      }
      else if (dr < 0.3 && fabs(dz) < 0.2) {
	muon_iso_region = true;
	if (dr < mindr_mu) {
	  mindr_mu = dr;
	  mindr_mu_idx = ihlt;
	}
      }
    }
    float mindr_mu_dz = 99.;
    if (muon_iso_region && mindr_mu_idx >= 0) {
      mindr_mu_dz = hlt_mus.at(mindr_mu_idx)->vz() - hlt_trk->vz();
    }
    if (muon_overlap) {
      if (verbose_) {
	cout << " * hlt muon track, track pt: " << hlt_trk->pt() << ", eta: " << hlt_trk->eta()
	     << ", phi: " << hlt_trk->phi() << ", vz: " << hlt_trk->vz()
	     << ", nhits: " << hlt_trk->numberOfValidHits() 
	     << ", algo: " <<hlt_trk->algo() << endl;
      }
      continue;
    }
    else if (!muon_iso_region) continue;
    else if (verbose_) {
      cout << " - hlt track. Closest muon: " << mindr_mu_idx 
	   << ", dr: " << mindr_mu << ", dz: " << mindr_mu_dz << endl
	   << "   pt: " << hlt_trk->pt() << ", eta: " << hlt_trk->eta()
	   << ", phi: " << hlt_trk->phi() << ", vz: " << hlt_trk->vz()
	   << ", nhits: " << hlt_trk->numberOfValidHits() 
	   << ", algo: " <<hlt_trk->algo() << endl;
    }

    hltTracksCollectionNoMuons.push_back(*hlt_trk);

    fillHists(*hlt_trk,"hltreg");
    TrackRef mu_trk = hlt_mus.at(mindr_mu_idx)->track();
    fillHistsMuHLT(*mu_trk,*hlt_trk,"hltreg");

    float mindr_pf = 99.;
    PFCandidateCollection::const_iterator mindr_pfcand = pfcandEnd;
    for ( PFCandidateCollection::const_iterator off_pf = offlinePFCandsCollection->begin(); off_pf != pfcandEnd; ++off_pf ) {
      // require charged cand
      if (abs(off_pf->charge()) != 1) continue;
      //      if (abs(off_pf->pdgId()) != 211) continue;
      //      if (vtx != 0 && vtx != -1) continue;
      float dr = ROOT::Math::VectorUtil::DeltaR(off_pf->p4(),hlt_trk->momentum());
      if (dr < mindr_pf) {
	mindr_pf = dr;
	mindr_pfcand = off_pf;
      }
    } // off pf cand loop

    // check for a match within dr < 0.02
    if (mindr_pf < 0.02) {
      int vtx = -2;
      if (!mindr_pfcand->trackRef().isNull()) vtx = chargedHadronVertex(*offlineVerticesCollection,*mindr_pfcand);

      if ((mindr_pfcand->particleId() == 1) && (vtx == 0 || vtx == -1)) {
	fillHists(*hlt_trk,"hltmatch");
	fillHistsMuHLT(*mu_trk,*hlt_trk,"hltmatch");
	fillHistsPFHLT(*mindr_pfcand,*hlt_trk,"hltmatch");
	hists_1d_["h_offvtx_hltmatch"]->Fill(vtx);
	if (verbose_) {
	  cout << "   - match to pf cand: dr: " << mindr_pf << ", vtx: " << vtx 
	       << ", vz: " << mindr_pfcand->vertex().z() << endl;
	}
      } else if (mindr_pfcand->particleId() != 1) {
	fillHists(*hlt_trk,"hltpfnonhad");
	fillHistsMuHLT(*mu_trk,*hlt_trk,"hltpfnonhad");
	fillHistsPFHLT(*mindr_pfcand,*hlt_trk,"hltpfnonhad");
	hists_1d_["h_offvtx_hltpfnonhad"]->Fill(vtx);
	if (verbose_) {
	  cout << "   @ match to non-had pf cand: dr: " << mindr_pf << ", vtx: " << vtx 
	       << ", vz: " << mindr_pfcand->vertex().z() << ", id: " << mindr_pfcand->particleId() << endl;
	}
      } else if (vtx != 0 && vtx != -1) {
	fillHists(*hlt_trk,"hltpfnonvtx");
	fillHistsMuHLT(*mu_trk,*hlt_trk,"hltpfnonvtx");
	fillHistsPFHLT(*mindr_pfcand,*hlt_trk,"hltpfnonvtx");
	hists_1d_["h_offvtx_hltpfnonvtx"]->Fill(vtx);
	if (verbose_) {
	  cout << "   & match to unused pf cand: dr: " << mindr_pf << ", vtx: " << vtx 
	       << ", vz: " << mindr_pfcand->vertex().z() << ", id: " << mindr_pfcand->particleId() << endl;
	}
      }
    } // pf cand match

    // only loop over offline tracks if we didn't find a charged pf cand
    else if (mindr_pf >= 0.02) {
      if (verbose_) {
	cout << "   ! no match to pf cands!! mindr: " << mindr_pf << endl;
      }

      float mindr_off = 99.;
      TrackCollection::const_iterator mindr_offtrk = offlineEnd;
      for ( TrackCollection::const_iterator off_trk = offlineTracksCollection->begin(); off_trk != offlineEnd; ++off_trk ) {
	float dr = ROOT::Math::VectorUtil::DeltaR(off_trk->momentum(),hlt_trk->momentum());
	if (dr < mindr_off) {
	  mindr_off = dr;
	  mindr_offtrk = off_trk;
	}
      } // off track loop

      // check for a match within dr < 0.02
      if (mindr_off > 0.02) {
	fillHists(*hlt_trk,"hltonly");
	fillHistsMuHLT(*mu_trk,*hlt_trk,"hltonly");
	hists_1d_["h_drmin_hltonly"]->Fill(mindr_off);
	if (verbose_) {
	  cout << "   ! no match to offline!! mindr: " << mindr_off << endl;
	}
      } else {
	fillHists(*hlt_trk,"hltnopf");
	fillHistsMuHLT(*mu_trk,*hlt_trk,"hltnopf");
	fillHistsRecoHLT(*mindr_offtrk,*hlt_trk,"hltnopf");
	if (verbose_) {
	  cout << "   - match to offline track: dr: " << mindr_off << endl;
	}
      }
    } // no match to pf


  } // loop over hlt tracks

  if (verbose_) cout << endl;

  const TrackCollection::const_iterator hltNoMuonsEnd = hltTracksCollectionNoMuons.end();

  // -------------------------------------------------------------------
  // --------- loop on offline tracks ----------------------------------
  // -------------------------------------------------------------------

  // loop on offline tracks, select those within the online tracking regions
  //  make plots for offline tracks, check for online matches
  //  TrackCollection offlineTracksCollectionReg;
  for ( TrackCollection::const_iterator off_trk = offlineTracksCollection->begin(); off_trk != offlineEnd; ++off_trk ) {
    // check for iter0 (algo 4) and high purity
    //    if (off_trk->algo() != 4) continue;
    if (!off_trk->quality(TrackBase::highPurity)) continue;

    if (off_trk->algo() == 4) fillHists(*off_trk,"offall_iter0");
    else fillHists(*off_trk,"offall_others");

    // check for region around hlt muons
    bool region = false;
    for (unsigned int ihlt=0; ihlt < hlt_mus.size(); ++ihlt) {
      float dr = ROOT::Math::VectorUtil::DeltaR(off_trk->momentum(),hlt_mus.at(ihlt)->p4());
      float dz = off_trk->vz() - hlt_mus.at(ihlt)->vz();
      if (dr > 0.01 && dr < 0.3 && fabs(dz) < 0.2) {
	region = true;
	break;
      }
      // float dphi = delta_phi(hlt_mus.at(ihlt)->phi(),off_trk->phi());
      // if ( (fabs(hlt_mus.at(ihlt)->eta() - off_trk->eta()) < 0.5) && (fabs(dphi) < 0.5) ) {
      // 	region = true;
      // 	break;
      // }
    } // loop over hlt muons
    //    if (!region) continue;

    if (region) {
      if (off_trk->algo() == 4) fillHists(*off_trk,"offreg_iter0");
      else fillHists(*off_trk,"offreg_others");
    }
    // offlineTracksCollectionReg.push_back(*off_trk);

    if (verbose_ && region) {
      cout << " - offline track: pt: " << off_trk->pt() << ", eta: " << off_trk->eta()
	   << ", phi: " << off_trk->phi() << ", nhits: " << off_trk->numberOfValidHits() 
	   << ", algo: " << off_trk->algo() << endl;
    }

    // next match to hlt tracks
    float mindr = 99.;
    TrackCollection::const_iterator hlt_match = hltNoMuonsEnd;
    for ( TrackCollection::const_iterator hlt_trk = hltTracksCollectionNoMuons.begin(); hlt_trk != hltNoMuonsEnd; ++hlt_trk ) {
      float dr = ROOT::Math::VectorUtil::DeltaR(off_trk->momentum(),hlt_trk->momentum()); 
      if (dr < mindr) {
	mindr = dr;
	hlt_match = hlt_trk;
      }
    } // hlt track loop

    // check for HLT track within dR < 0.02
    if (mindr < 0.02) {
      //      fillHistsRecoHLT(*off_trk,*hlt_match,"offhlt");
      if (verbose_ && region) {
	cout << "   - hlt match, dr = " << mindr << ", track pt: " << hlt_match->pt() << ", eta: " << hlt_match->eta()
	     << ", phi: " << hlt_match->phi() << ", nhits: " << hlt_match->numberOfValidHits() 
	     << ", algo: " << hlt_match->algo() << endl;
      }
    } 
    else if (region) {
      if (off_trk->algo() == 4) fillHists(*off_trk,"offonly_iter0");
      else fillHists(*off_trk,"offonly_others");
    }

  } // off track loop

  if (verbose_) cout << endl;

  return;
}

//____________________________________________________________________________
void TrackCompAnalyzer::bookHists(edm::Service<TFileService>& fs, const std::string& suffix) {

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  hists_1d_["h_pt"+suf] = fs->make<TH1F>(Form("h_pt%s",suf.c_str()) , "; p_{T} [GeV]" , 500 , 0. , 50. );
  hists_1d_["h_eta"+suf] = fs->make<TH1F>(Form("h_eta%s",suf.c_str()) , "; #eta" , 100 , -3. , 3. );
  hists_1d_["h_nhits"+suf] = fs->make<TH1F>(Form("h_nhits%s",suf.c_str()) , "; N(hits)" , 31 , -0.5 , 30.5 );
  hists_1d_["h_chi2"+suf] = fs->make<TH1F>(Form("h_chi2%s",suf.c_str()) , "; #chi^{2}" , 500 , 0. , 500. );
  hists_1d_["h_normchi2"+suf] = fs->make<TH1F>(Form("h_normchi2%s",suf.c_str()) , "; Normalized #chi^{2}" , 300 , 0. , 30. );
  hists_1d_["h_dz"+suf] = fs->make<TH1F>(Form("h_dz%s",suf.c_str()) , "; d_{z} [cm]" , 500 , -10. , 10. );
  hists_1d_["h_dxy_bs"+suf] = fs->make<TH1F>(Form("h_dxy_bs%s",suf.c_str()) , "; d_{xy} wrt BeamSpot [cm]" , 100 , -0.5 , 0.5 );
  hists_1d_["h_algo"+suf] = fs->make<TH1F>(Form("h_algo%s",suf.c_str()) , "; Algo" , 15 , -0.5 , 14.5 );


  return;
}

//____________________________________________________________________________
void TrackCompAnalyzer::bookHistsRecoHLT(edm::Service<TFileService>& fs, const std::string& suffix) {

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  hists_1d_["h_dr"+suf] = fs->make<TH1F>(Form("h_dr%s",suf.c_str()) , "; #DeltaR(off,HLT)" , 600 , 0. , 6. );
  hists_1d_["h_dpt"+suf] = fs->make<TH1F>(Form("h_dpt%s",suf.c_str()) , "; (p_{T}^{off} - p_{T}^{HLT}) / p_{T}^{off}" , 500 , -5. , 5. );
  hists_1d_["h_dzoffhlt"+suf] = fs->make<TH1F>(Form("h_dzoffhlt%s",suf.c_str()) , "; #Deltaz(off,HLT) [cm]" , 500 , -10. , 10. );
  hists_1d_["h_dnhits"+suf] = fs->make<TH1F>(Form("h_dnhits%s",suf.c_str()) , "; N(hits,off) - N(hits,HLT)" , 41 , -10.5 , 30.5 );
  hists_1d_["h_algooff"+suf] = fs->make<TH1F>(Form("h_algooff%s",suf.c_str()) , "; Offline Algo" , 15 , -0.5 , 14.5 );
  hists_1d_["h_qualoff"+suf] = fs->make<TH1F>(Form("h_qualoff%s",suf.c_str()) , "; Offline Quality" , 3 , -0.5 , 2.5 );

  return;
}

//____________________________________________________________________________
void TrackCompAnalyzer::bookHistsPFHLT(edm::Service<TFileService>& fs, const std::string& suffix) {

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  hists_1d_["h_pfid"+suf] = fs->make<TH1F>(Form("h_pfid%s",suf.c_str()) , "; PF Particle Id" , 3 , 0.5 , 3.5 );

  bookHistsRecoHLT(fs,suffix);

  return;
}

//____________________________________________________________________________
void TrackCompAnalyzer::bookHistsMuHLT(edm::Service<TFileService>& fs, const std::string& suffix) {

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  hists_1d_["h_drmuhlt"+suf] = fs->make<TH1F>(Form("h_drmuhlt%s",suf.c_str()) , "; #DeltaR(#mu,HLT)" , 600 , 0. , 6. );
  hists_1d_["h_dzmuhlt"+suf] = fs->make<TH1F>(Form("h_dzmuhlt%s",suf.c_str()) , "; #Deltaz(#mu,HLT) [cm]" , 500 , -10. , 10. );

  return;
}

//____________________________________________________________________________
void TrackCompAnalyzer::fillHists(const reco::Track& trk, const std::string& suffix) {

  if (trk.pt() <= 0.) {
    std::cout << "TrackCompAnalyzer::fillHists: invalid trk pt: " << trk.pt() << std::endl;
    return;
  }

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  hists_1d_["h_pt"+suf]->Fill(trk.pt());
  hists_1d_["h_eta"+suf]->Fill(trk.eta());
  hists_1d_["h_nhits"+suf]->Fill(trk.numberOfValidHits());
  hists_1d_["h_chi2"+suf]->Fill(trk.chi2());
  hists_1d_["h_normchi2"+suf]->Fill(trk.normalizedChi2());
  hists_1d_["h_dz"+suf]->Fill(trk.dz());
  hists_1d_["h_dxy_bs"+suf]->Fill(trk.dxy(*(beamSpotHandle_.product())));
  hists_1d_["h_algo"+suf]->Fill(trk.algo());

  return;
}

//____________________________________________________________________________
void TrackCompAnalyzer::fillHistsRecoHLT(const reco::Track& off_trk, const reco::Track& hlt_trk, const std::string& suffix) {

  using namespace reco;

  if (off_trk.pt() <= 0.) {
    std::cout << "TrackCompAnalyzer::fillHistsRecoHLT: invalid off_trk pt: " << off_trk.pt() << std::endl;
    return;
  }

  if (hlt_trk.pt() <= 0.) {
    std::cout << "TrackCompAnalyzer::fillHistsRecoHLT: invalid hlt_trk pt: " << hlt_trk.pt() << std::endl;
    return;
  }

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  float dr = ROOT::Math::VectorUtil::DeltaR(off_trk.momentum(),hlt_trk.momentum());
  hists_1d_["h_dpt"+suf]->Fill( (off_trk.pt() - hlt_trk.pt()) / off_trk.pt() );
  hists_1d_["h_dr"+suf]->Fill(dr);
  hists_1d_["h_dzoffhlt"+suf]->Fill(off_trk.vz() - hlt_trk.vz());
  hists_1d_["h_dnhits"+suf]->Fill(off_trk.numberOfValidHits() - hlt_trk.numberOfValidHits());
  hists_1d_["h_algooff"+suf]->Fill(off_trk.algo());

  int quality = int(TrackBase::loose);
  if (off_trk.quality(TrackBase::highPurity)) quality = int(TrackBase::highPurity);
  else if (off_trk.quality(TrackBase::tight)) quality = int(TrackBase::tight);
  hists_1d_["h_qualoff"+suf]->Fill(quality);

  return;
}

//____________________________________________________________________________
void TrackCompAnalyzer::fillHistsPFHLT(const reco::PFCandidate& off_pf, const reco::Track& hlt_trk, const std::string& suffix) {

  using namespace reco;

  if (off_pf.pt() <= 0.) {
    std::cout << "TrackCompAnalyzer::fillHistsPFHLT: invalid off_pf pt: " << off_pf.pt() << std::endl;
    return;
  }

  if (hlt_trk.pt() <= 0.) {
    std::cout << "TrackCompAnalyzer::fillHistsPFHLT: invalid hlt_trk pt: " << hlt_trk.pt() << std::endl;
    return;
  }

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  hists_1d_["h_pfid"+suf]->Fill(off_pf.particleId());

  if (!off_pf.trackRef().isNull())  fillHistsRecoHLT(*(off_pf.trackRef()),hlt_trk,suffix);

  return;
}

//____________________________________________________________________________
void TrackCompAnalyzer::fillHistsMuHLT(const reco::Track& hlt_mu, const reco::Track& hlt_trk, const std::string& suffix) {

  if (hlt_mu.pt() <= 0.) {
    std::cout << "TrackCompAnalyzer::fillHistsMuHLT: invalid hlt_mu pt: " << hlt_mu.pt() << std::endl;
    return;
  }

  if (hlt_trk.pt() <= 0.) {
    std::cout << "TrackCompAnalyzer::fillHistsMuHLT: invalid hlt_trk pt: " << hlt_trk.pt() << std::endl;
    return;
  }

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  float dr = ROOT::Math::VectorUtil::DeltaR(hlt_mu.momentum(),hlt_trk.momentum());
  hists_1d_["h_drmuhlt"+suf]->Fill(dr);
  hists_1d_["h_dzmuhlt"+suf]->Fill(hlt_mu.vz() - hlt_trk.vz());

  return;
}

//____________________________________________________________________________
// returns dphi in range [-pi,pi[
float TrackCompAnalyzer::delta_phi(float phi1, float phi2) {
  float dphi = phi1 - phi2;
  if (dphi < -ROOT::Math::Pi()) dphi += 2*ROOT::Math::Pi();
  else if (dphi >= ROOT::Math::Pi()) dphi -= 2*ROOT::Math::Pi();

  return dphi;
}

//____________________________________________________________________________
// returns best match vertex.  Taken from PFPileUpAlgo, chargedHadronVertex
int TrackCompAnalyzer::chargedHadronVertex( const reco::VertexCollection& vertices, const reco::PFCandidate& pfcand ) const {
  if (pfcand.trackRef().isNull()) return -2;
  return trackVertex(vertices,pfcand.trackRef());
}

//____________________________________________________________________________
// returns best match vertex.  Taken from PFPileUpAlgo, chargedHadronVertex
int TrackCompAnalyzer::trackVertex( const reco::VertexCollection& vertices, const reco::TrackRef& trackRef ) const {

  reco::TrackBaseRef trackBaseRef( trackRef );
  
  size_t  iVertex = 0;
  unsigned index=0;
  unsigned nFoundVertex = 0;
  typedef reco::VertexCollection::const_iterator IV;
  typedef reco::Vertex::trackRef_iterator IT;
  float bestweight=0;
  for(IV iv=vertices.begin(); iv!=vertices.end(); ++iv, ++index) {

    const reco::Vertex& vtx = *iv;
    
    // loop on tracks in vertices
    for(IT iTrack=vtx.tracks_begin(); 
	iTrack!=vtx.tracks_end(); ++iTrack) {
	 
      const reco::TrackBaseRef& baseRef = *iTrack;

      // one of the tracks in the vertex is the same as 
      // the track considered in the function
      if(baseRef == trackBaseRef ) {
	float w = vtx.trackWeight(baseRef);
	//select the vertex for which the track has the highest weight
	if (w > bestweight){
	  bestweight=w;
	  iVertex=index;
	  nFoundVertex++;
	}	 	
      }
    }
  }

  if (nFoundVertex>0){
    if (nFoundVertex!=1)
      edm::LogWarning("TrackOnTwoVertex")<<"a track is shared by at least two verteces. Used to be an assert";
    return iVertex;
  }
  // no vertex found with this track. 

  const bool checkClosestZVertex_ = true;
  // optional: as a secondary solution, associate the closest vertex in z
  if ( checkClosestZVertex_ ) {

    double dzmin = 10000;
    double ztrack = trackRef->vertex().z();
    bool foundVertex = false;
    index = 0;
    for(IV iv=vertices.begin(); iv!=vertices.end(); ++iv, ++index) {

      double dz = fabs(ztrack - iv->z());
      if(dz<dzmin) {
	dzmin = dz; 
	iVertex = index;
	foundVertex = true;
      }
    }

    if( foundVertex ) 
      return iVertex;  

  }


  return -1 ;
}

