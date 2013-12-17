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
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
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
  bookHists(fs,"hlt");
  bookHists(fs,"hltonly");
  hists_1d_["h_drmin_hltonly"] = fs->make<TH1F>(Form("h_drmin_hltonly") , "; min #DeltaR(off,HLT)" , 600 , 0. , 6. );
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

  vector<LorentzVector> hlt_mus;

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
	    if (ROOT::Math::VectorUtil::DeltaR(hlt_mus.at(k),lv) < 0.001) {
	      overlap = true;
	      break;
	    }
	  }

	  if (overlap) continue;
	  hlt_mus.push_back(lv);

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
  const TrackCollection::const_iterator hltEnd = hltTracksCollection->end();
  const TrackCollection::const_iterator offlineEnd = offlineTracksCollection->end();

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
      float dphi = delta_phi(hlt_mus.at(ihlt).phi(),off_trk->phi());
      if ( (fabs(hlt_mus.at(ihlt).eta() - off_trk->eta()) < 0.5) && (fabs(dphi) < 0.5) ) {
      	region = true;
      	break;
      }
    } // loop over hlt muons
    //    if (!region) continue;

    if (region) {
      if (off_trk->algo() == 4) fillHists(*off_trk,"offreg_iter0");
      else fillHists(*off_trk,"offreg_others");
    }
    // offlineTracksCollectionReg.push_back(*off_trk);

    if (verbose_) {
      cout << " - offline track: pt: " << off_trk->pt() << ", eta: " << off_trk->eta()
	   << ", phi: " << off_trk->phi() << ", nhits: " << off_trk->numberOfValidHits() 
	   << ", algo: " << off_trk->algo() << endl;
    }

    // next match to hlt tracks
    float mindr = 99.;
    TrackCollection::const_iterator hlt_match = hltEnd;
    for ( TrackCollection::const_iterator hlt_trk = hltTracksCollection->begin(); hlt_trk != hltEnd; ++hlt_trk ) {
      float dr = ROOT::Math::VectorUtil::DeltaR(off_trk->momentum(),hlt_trk->momentum()); 
      if (dr < mindr) {
	mindr = dr;
	hlt_match = hlt_trk;
      }
    } // hlt track loop

    // check for HLT track within dR < 0.02
    if (mindr < 0.02) {
      fillHistsRecoHLT(*off_trk,*hlt_match,"offhlt");
      if (verbose_) {
	cout << "   - hlt match, dr = " << mindr << ", track pt: " << hlt_match->pt() << ", eta: " << hlt_match->eta()
	     << ", phi: " << hlt_match->phi() << ", nhits: " << hlt_match->numberOfValidHits() << endl;
      }
    } 
    else if (region) {
      if (off_trk->algo() == 4) fillHists(*off_trk,"offonly_iter0");
      else fillHists(*off_trk,"offonly_others");
    }

  } // off track loop

  if (verbose_) cout << endl;

  // loop on hlt tracks, make plots
  //  also check for offline matches, plot online only tracks
  //  const TrackCollection::const_iterator offlineRegEnd = offlineTracksCollectionReg.end();
  for ( TrackCollection::const_iterator hlt_trk = hltTracksCollection->begin(); hlt_trk != hltEnd; ++hlt_trk ) {

    // check dR from hlt muons: those tracks aren't used in track iso
    bool muon_overlap = false;
    for (unsigned int ihlt=0; ihlt < hlt_mus.size(); ++ihlt) {
      if (ROOT::Math::VectorUtil::DeltaR(hlt_trk->momentum(),hlt_mus.at(ihlt)) < 0.01) {
	muon_overlap = true;
	break;
      }
    }
    if (muon_overlap) continue;

    fillHists(*hlt_trk,"hlt");

    float mindr = 99.;
    for ( TrackCollection::const_iterator off_trk = offlineTracksCollection->begin(); off_trk != offlineEnd; ++off_trk ) {
      float dr = ROOT::Math::VectorUtil::DeltaR(off_trk->momentum(),hlt_trk->momentum());
      if (dr < mindr) mindr = dr;
    } // off reg track loop

    // check for a match within dr < 0.02
    if (mindr > 0.02) {
      fillHists(*hlt_trk,"hltonly");
      hists_1d_["h_drmin_hltonly"]->Fill(mindr);
      if (verbose_) {
	cout << " - unmatched hlt, mindr: " << mindr << ", track pt: " << hlt_trk->pt() << ", eta: " << hlt_trk->eta()
	     << ", phi: " << hlt_trk->phi() << ", nhits: " << hlt_trk->numberOfValidHits() << endl;
      }
    }

  } // hlt track loop

  //-------------------------------------
  //   reco vertices
  //-------------------------------------

  // // find vertex 0 in vertex container
  // const VertexCollection* vertexCollection = vertexHandle_.product();
  // VertexCollection::const_iterator firstGoodVertex = vertexCollection->end();
  // for ( VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); ++vtx ) {
  //   if (  !vtx->isFake() && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0 ) {
  //     firstGoodVertex = vtx;
  //     break;
  //   }
  // } // loop on vertices

  // if (firstGoodVertex == vertexCollection->end()) {
  //   cout << "TrackCompAnalyzer::analyze: didn't find any good offline vertices!! size: " 
  // 	 << vertexCollection->size() << std::endl;
  //   return true;
  // }

  if (verbose_) cout << endl;

  return;
}

//____________________________________________________________________________
void TrackCompAnalyzer::bookHists(edm::Service<TFileService>& fs, const std::string& suffix) {

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  hists_1d_["h_pt"+suf] = fs->make<TH1F>(Form("h_pt%s",suf.c_str()) , "; p_{T} [GeV]" , 500 , 0. , 50. );
  hists_1d_["h_eta"+suf] = fs->make<TH1F>(Form("h_eta%s",suf.c_str()) , "; #eta" , 100 , -3. , 3. );
  hists_1d_["h_nhits"+suf] = fs->make<TH1F>(Form("h_nhits%s",suf.c_str()) , "; N(hits)" , 30 , -0.5 , 30.5 );
  hists_1d_["h_chi2"+suf] = fs->make<TH1F>(Form("h_chi2%s",suf.c_str()) , "; #chi^{2}" , 500 , 0. , 500. );
  hists_1d_["h_normchi2"+suf] = fs->make<TH1F>(Form("h_normchi2%s",suf.c_str()) , "; Normalized #chi^{2}" , 300 , 0. , 30. );
  hists_1d_["h_dxy"+suf] = fs->make<TH1F>(Form("h_dxy%s",suf.c_str()) , "; d_{xy} [cm]" , 100 , -0.5 , 0.5 );
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
  hists_1d_["h_dnhits"+suf] = fs->make<TH1F>(Form("h_dnhits%s",suf.c_str()) , "; N(hits,off) - N(hits,HLT)" , 40 , -10.5 , 30.5 );

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
  hists_1d_["h_dxy"+suf]->Fill(trk.dxy());
  hists_1d_["h_dz"+suf]->Fill(trk.dz());
  hists_1d_["h_dxy_bs"+suf]->Fill(trk.dxy(*(beamSpotHandle_.product())));
  hists_1d_["h_algo"+suf]->Fill(trk.algo());

  return;
}

//____________________________________________________________________________
void TrackCompAnalyzer::fillHistsRecoHLT(const reco::Track& off_trk, const reco::Track& hlt_trk, const std::string& suffix) {

  if (off_trk.pt() <= 0.) {
    std::cout << "TrackCompAnalyzer::fillHists: invalid off_trk pt: " << off_trk.pt() << std::endl;
    return;
  }

  if (hlt_trk.pt() <= 0.) {
    std::cout << "TrackCompAnalyzer::fillHists: invalid hlt_trk pt: " << hlt_trk.pt() << std::endl;
    return;
  }

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  float dr = ROOT::Math::VectorUtil::DeltaR(off_trk.momentum(),hlt_trk.momentum());
  hists_1d_["h_dpt"+suf]->Fill( (off_trk.pt() - hlt_trk.pt()) / off_trk.pt() );
  hists_1d_["h_dr"+suf]->Fill(dr);
  hists_1d_["h_dnhits"+suf]->Fill(off_trk.numberOfValidHits() - hlt_trk.numberOfValidHits());

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
