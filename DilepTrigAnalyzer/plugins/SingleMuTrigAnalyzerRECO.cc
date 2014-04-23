/** \class SingleMuTrigAnalyzerRECO
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
#include "TrigAnalyzer/DilepTrigAnalyzer/interface/SingleMuTrigAnalyzerRECO.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

// ROOT includes
#include "Math/VectorUtil.h"

#include <cassert>

using namespace reco;
using namespace edm;

//
// constructors and destructor
//
//____________________________________________________________________________
SingleMuTrigAnalyzerRECO::SingleMuTrigAnalyzerRECO(const edm::ParameterSet& ps) : 
  processName_(ps.getParameter<std::string>("processName")),
  baseTriggerName_(ps.getParameter<std::string>("baseTriggerName")),
  tkTriggerName_(ps.getParameter<std::string>("tkTriggerName")),
  triggerResultsTag_(ps.getParameter<edm::InputTag>("triggerResults")),
  triggerEventWithRefsTag_(ps.getParameter<edm::InputTag>("triggerEventWithRefs")),
  muonsInputTag_(ps.getParameter<edm::InputTag>("muonsInputTag")),
  vtxInputTag_(ps.getParameter<edm::InputTag>("vtxInputTag")),
  hltVtxInputTag_(ps.getParameter<edm::InputTag>("hltVtxInputTag")),
  dumpHLTTracks_(ps.getParameter<bool>("dumpHLTTracks")),
  hltTracksTag_(ps.getParameter<edm::InputTag>("hltTracks")),
  reqTrigMatch_(ps.getParameter<bool>("reqTrigMatch")),
  offPt_(ps.getParameter<double>("offPt")),
  doOffGenMatch_(ps.getParameter<bool>("doOffGenMatch")),
  genParticlesTag_(ps.getParameter<edm::InputTag>("genParticles")),
  compareToGen_(ps.getParameter<bool>("compareToGen")),
  verbose_(ps.getParameter<bool>("verbose"))
{
  using namespace std;
  using namespace edm;

  cout << "SingleMuTrigAnalyzerRECO configuration: " << endl
       << "   ProcessName = " << processName_ << endl
       << "   baseTriggerName = " << baseTriggerName_ << endl
       << "   tkTriggerName = " << tkTriggerName_ << endl
       << "   TriggerResultsTag = " << triggerResultsTag_.encode() << endl
       << "   TriggerEventWithRefsTag = " << triggerEventWithRefsTag_.encode() << endl
       << "   MuonsInputTag = " << muonsInputTag_.encode() << endl
       << "   HLTVtxInputTag = " << hltVtxInputTag_.encode() << endl
       << "   VtxInputTag = " << vtxInputTag_.encode() << endl
       << "   dumpHLTTracks = " << dumpHLTTracks_ << endl
       << "   HLTTracksTag = " << hltTracksTag_.encode() << endl
       << "   ReqTrigMatch = " << reqTrigMatch_ << endl
       << "   OffPt = " << offPt_ << endl
       << "   DoOffGenMatch = " << doOffGenMatch_ << endl
       << "   GenParticlesTag = " << genParticlesTag_.encode() << endl
       << "   CompareToGen = " << compareToGen_ << endl
       << "   Verbose = " << verbose_ << endl;

  hltTriggerNames_.push_back(tkTriggerName_);
  hltTriggerNames_.push_back(baseTriggerName_);

  hltShortNames_.push_back("tkmu");
  hltShortNames_.push_back("mu");

  // histogram setup
  edm::Service<TFileService> fs;
  h_results_ = fs->make<TH1F>("h_results" , ";Trigger Results" , 4 , -0.5 , 3.5 );
  for (unsigned int itrig=0; itrig<=(unsigned int)mu; ++itrig) {
    bookHists(fs,hltShortNames_.at(itrig),false);
    bookHistsGen(fs,hltShortNames_.at(itrig)+"_gen");
    bookHistsGen(fs,hltShortNames_.at(itrig)+"_gen_match");
    bookHistsGen(fs,hltShortNames_.at(itrig)+"_gen_nomatch");
    bookHists(fs,hltShortNames_.at(itrig)+"_tight");
    bookHists(fs,hltShortNames_.at(itrig)+"_tight_match");
    bookHists(fs,hltShortNames_.at(itrig)+"_tight_nomatch");
  }

}

//____________________________________________________________________________
SingleMuTrigAnalyzerRECO::~SingleMuTrigAnalyzerRECO()
{
}

//
// member functions
//
//____________________________________________________________________________
void
SingleMuTrigAnalyzerRECO::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  using namespace std;
  using namespace edm;

  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
      const unsigned int n(hltConfig_.size());
      // check if trigger names in (new) config
      for (unsigned int itrig = 0; itrig <= (unsigned int)mu; ++itrig) {
	unsigned int triggerIndex(hltConfig_.triggerIndex(hltTriggerNames_.at(itrig)));
	if (triggerIndex>=n) {
	  cout << "SingleMuTrigAnalyzerRECO::analyze:"
	       << " TriggerName " << hltTriggerNames_.at(itrig) 
	       << " not available in (new) config!" << endl;
	}
      } // loop over triggers
    } // if changed
  } else {
    cout << "SingleMuTrigAnalyzerRECO::analyze:"
	 << " config extraction failure with process name "
	 << processName_ << endl;
  }

}

//____________________________________________________________________________
// ------------ method called to produce the data  ------------
void
SingleMuTrigAnalyzerRECO::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;

  if (verbose_) cout << endl;

  // get event products
  iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
    cout << "SingleMuTrigAnalyzerRECO::analyze: Error in getting TriggerResults product from Event!" << endl;
    return;
  }
  iEvent.getByLabel(triggerEventWithRefsTag_,triggerEventWithRefsHandle_);
  if (!triggerEventWithRefsHandle_.isValid()) {
    cout << "SingleMuTrigAnalyzerRECO::analyze: Error in getting TriggerEventWithRefs product from Event!" << endl;
    return;
  }
  // sanity check
  assert(triggerResultsHandle_->size()==hltConfig_.size());

  // retrieve necessary containers
  iEvent.getByLabel(vtxInputTag_, vertexHandle_);
  iEvent.getByLabel( muonsInputTag_ , musHandle_ );

  // for printing out HLT PF Cand info
  if (dumpHLTTracks_) {
    iEvent.getByLabel(hltTracksTag_, hltTracksHandle_);
    iEvent.getByLabel(hltVtxInputTag_, hltVertexHandle_);
  }

  if (doOffGenMatch_ || compareToGen_) {
    iEvent.getByLabel(genParticlesTag_, genParticlesHandle_);
  }
  
  trigpass_results_ = 0;
  for (unsigned int itrig=0; itrig <= (unsigned int)mu; ++itrig) {
    bool pass = analyzeTrigger(iEvent,iSetup,(hltTrigs)itrig);
    if (pass) trigpass_results_ |= 1 << itrig;
  }
  h_results_->Fill(trigpass_results_);

  if (verbose_) cout << endl;

  return;
}

//____________________________________________________________________________

bool SingleMuTrigAnalyzerRECO::analyzeTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup, const hltTrigs triggerEnum) {
  
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
    cout << "SingleMuTrigAnalyzerRECO::analyzeTrigger: path "
	 << triggerName << " - not found!" << endl;
    return false;
  }

  if (verbose_) {
    cout << "SingleMuTrigAnalyzerRECO::analyzeTrigger: path "
	 << triggerName << " [" << triggerIndex << "]" << endl;
  }
  // modules on this trigger path
  const unsigned int m(hltConfig_.size(triggerIndex));
  const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));

  bool wasRun = triggerResultsHandle_->wasrun(triggerIndex);
  bool accept = triggerResultsHandle_->accept(triggerIndex);
  bool error = triggerResultsHandle_->error(triggerIndex);
  const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
  bool passedLevel1 = (moduleLabels[moduleIndex].find("hltL1sMu16") == string::npos);
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

  //  if (!wasRun || !accept || error) return false;
  if (!wasRun || error || !passedLevel1) return false;

  // loop over trigger and reco objects, match, make plots

  // first, get trigger objects from last filter

  //------------------------------------
  //  hlt objects
  //------------------------------------

  // // Results from TriggerEventWithRefs product
  // muonIds_.clear();
  // muonRefs_.clear();
  bool foundMuons = false;

  StudyLepton hlt_mu_lead;
  StudyLepton hlt_mu_subl;
  StudyLepton hlt_mu_third;
  hlt_mu_lead.type = 13;
  hlt_mu_subl.type = 13;
  hlt_mu_third.type = 13;
  hlt_mu_lead.isHLT = true;
  hlt_mu_subl.isHLT = true;
  hlt_mu_third.isHLT = true;
  hlt_mu_lead.isGen = false;
  hlt_mu_subl.isGen = false;
  hlt_mu_third.isGen = false;
  int hlt_mu_lead_idx = -1;
  int hlt_mu_subl_idx = -1;
  int hlt_mu_third_idx = -1;
  //  int hlt_mu_filter_idx = -1;

  // Results from TriggerEvent product - Attention: must look only for
  // modules actually run in this path for this event!
  // -- loop backwards through modules to find last filter for each object type
  //  for (unsigned int j=0; j<=moduleIndex; ++j) {
  for (unsigned int j=moduleIndex; j!=0; --j) {
    const string& moduleLabel(moduleLabels[j]);
    const string  moduleType(hltConfig_.moduleType(moduleLabel));
    // check whether the module is packed up in TriggerEvent product
    const unsigned int filterIndex(triggerEventWithRefsHandle_->filterIndex(InputTag(moduleLabel,"",processName_)));
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
	     << "  - the objects: # id pt eta phi vz id key eta_trkref phi_trkref" << endl;
      }
      for (unsigned int i=0; i!=nMuons; ++i) {

	float iso = -99.;
	if (!foundMuons) {
	  LorentzVector lv(muonRefs_.at(i)->p4());

	  if ( (hlt_mu_lead_idx == -1) || (lv.pt() > hlt_mu_lead.lv.pt()) ) {
	    hlt_mu_third_idx = hlt_mu_subl_idx;
	    hlt_mu_third.lv = hlt_mu_subl.lv;
	    hlt_mu_third.trkiso = hlt_mu_subl.trkiso;
	    hlt_mu_third.vz = hlt_mu_subl.vz;
	    hlt_mu_third.charge = hlt_mu_subl.charge;
	    hlt_mu_subl_idx = hlt_mu_lead_idx;
	    hlt_mu_subl.lv = hlt_mu_lead.lv;
	    hlt_mu_subl.trkiso = hlt_mu_lead.trkiso;
	    hlt_mu_subl.vz = hlt_mu_lead.vz;
	    hlt_mu_subl.charge = hlt_mu_lead.charge;
	    hlt_mu_lead_idx = (int)i;
	    hlt_mu_lead.lv = lv;
	    hlt_mu_lead.trkiso = iso;
	    hlt_mu_lead.vz = muonRefs_.at(i)->vz();
	    hlt_mu_lead.charge = muonRefs_.at(i)->charge();
	  } else if ( (hlt_mu_subl_idx == -1) || (lv.pt() > hlt_mu_subl.lv.pt()) ) {
	    if (ROOT::Math::VectorUtil::DeltaR(hlt_mu_lead.lv,lv) > 0.001) {
	      hlt_mu_third_idx = hlt_mu_subl_idx;
	      hlt_mu_third.lv = hlt_mu_subl.lv;
	      hlt_mu_third.trkiso = hlt_mu_subl.trkiso;
	      hlt_mu_third.vz = hlt_mu_subl.vz;
	      hlt_mu_third.charge = hlt_mu_subl.charge;
	      hlt_mu_subl_idx = (int)i;
	      hlt_mu_subl.lv = lv;
	      hlt_mu_subl.trkiso = iso;
	      hlt_mu_subl.vz = muonRefs_.at(i)->vz();
	      hlt_mu_subl.charge = muonRefs_.at(i)->charge();
	    }
	  } else if ( (hlt_mu_third_idx == -1) || (lv.pt() > hlt_mu_third.lv.pt()) ) {
	    if ( (ROOT::Math::VectorUtil::DeltaR(hlt_mu_lead.lv,lv) > 0.001) &&
		 (ROOT::Math::VectorUtil::DeltaR(hlt_mu_subl.lv,lv) > 0.001) ) {
	      hlt_mu_third_idx = (int)i;
	      hlt_mu_third.lv = lv;
	      hlt_mu_third.trkiso = iso;
	      hlt_mu_third.vz = muonRefs_.at(i)->vz();
  	      hlt_mu_third.charge = muonRefs_.at(i)->charge();
	    }
	  }
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

    if (hlt_mu_lead_idx >= 0) {
      foundMuons = true;
      //      hlt_mu_filter_idx = filterIndex;
    }

    if (foundMuons) break;
  } // backwards loop on modules

  if (hlt_mu_lead_idx == -1) {
    cout << "SingleMuTrigAnalyzerRECO::analyzeTrigger: no valid trigger leptons!" << endl;
    //    return true;
  }

  // make plots based on trigger path, if the trigger was passed
  if (accept) {
    //    fillHists(hlt_mu_lead,hlt_mu_subl,triggerShort,true);
  }

  //-------------------------------------
  //   compare to gen level
  //-------------------------------------

  const float dr_trigmatch = 0.2;

  StudyLepton gen_mu;
  gen_mu.type = 13;
  gen_mu.isHLT = false;
  gen_mu.isGen = true;
  bool gen_mu_trigmatch = false;
  bool found_gen_mu = false;

  if (compareToGen_) {
      GenParticleCollection::const_iterator genps_end = genParticlesHandle_->end();  // Iterator
      for ( GenParticleCollection::const_iterator genp = genParticlesHandle_->begin(); genp != genps_end; ++genp ) {
	//      if (genp->status() != 3) continue;
	// -- pythia 8..
	if ((genp->status() != 1) && (genp->status() != 23)) continue;
	if (abs(genp->pdgId()) != 13) continue;
	// W mother
	if (abs(genp->mother()->pdgId()) != 24) continue;

	found_gen_mu = true;

	// pt, eta acceptance cuts
	if ((genp->pt() > offPt_) && (fabs(genp->eta()) < 2.4)) {
	  gen_mu.lv = genp->p4();

	  gen_mu_trigmatch = false;
	  if (accept) {
	    if ( (hlt_mu_lead_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(genp->p4(),hlt_mu_lead.lv) < dr_trigmatch) ) gen_mu_trigmatch = true;
	    else if ( (hlt_mu_subl_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(genp->p4(),hlt_mu_subl.lv) < dr_trigmatch) ) gen_mu_trigmatch = true;
	    else if ( (hlt_mu_third_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(genp->p4(),hlt_mu_third.lv) < dr_trigmatch) ) gen_mu_trigmatch = true;
	  }
	}

	break;
      } // loop over genps

      if (!found_gen_mu) cout << "!! WARNING! didn't find gen muon.." << endl;

      if (gen_mu.lv.pt() > 0.) {
	fillHistsGen(gen_mu,triggerShort+"_gen");
	if (gen_mu_trigmatch) {
	  fillHistsGen(gen_mu,triggerShort+"_gen_match");
	} else {
	  fillHistsGen(gen_mu,triggerShort+"_gen_nomatch");
	}
      }

  } // if compareToGen

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
    cout << "SingleMuTrigAnalyzerRECO::analyzeTrigger: didn't find any good offline vertices!! size: " 
	 << vertexCollection->size() << std::endl;
    return accept;
  }

  //-------------------------------------
  //   reco muons 
  //-------------------------------------

  StudyLepton off_mu_lead;
  StudyLepton off_mu_subl;
  off_mu_lead.type = 13;
  off_mu_subl.type = 13;
  off_mu_lead.isHLT = false;
  off_mu_subl.isHLT = false;
  off_mu_lead.isGen = false;
  off_mu_subl.isGen = false;
  int off_mu_lead_idx = -1;
  int off_mu_subl_idx = -1;
  bool off_mu_lead_trigmatch = false;
  bool off_mu_subl_trigmatch = false;

  StudyLepton off_mu_tight_lead;
  StudyLepton off_mu_tight_subl;
  off_mu_tight_lead.type = 13;
  off_mu_tight_subl.type = 13;
  off_mu_tight_lead.isHLT = false;
  off_mu_tight_subl.isHLT = false;
  off_mu_tight_lead.isGen = false;
  off_mu_tight_subl.isGen = false;
  int off_mu_tight_lead_idx = -1;
  int off_mu_tight_subl_idx = -1;
  bool off_mu_tight_lead_trigmatch = false;
  bool off_mu_tight_subl_trigmatch = false;

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
    if (accept) {
      if ( (hlt_mu_lead_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(lv,hlt_mu_lead.lv) < dr_trigmatch) ) match = true;
      else if ( (hlt_mu_subl_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(lv,hlt_mu_subl.lv) < dr_trigmatch) ) match = true;
      else if ( (hlt_mu_third_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(lv,hlt_mu_third.lv) < dr_trigmatch) ) match = true;
    }
    if (!match && reqTrigMatch_) {
      if (verbose_) cout << ", FAILS trig match" << std::endl;
      continue;
    }

    // match to gen object: status 3 muons
    if (doOffGenMatch_) {
      match = false;
      GenParticleCollection::const_iterator genps_end = genParticlesHandle_->end();  // Iterator
      for ( GenParticleCollection::const_iterator genp = genParticlesHandle_->begin(); genp != genps_end; ++genp ) {
	//      if (genp->status() != 3) continue;
	// -- pythia 8..
	if ((genp->status() != 1) && (genp->status() != 23)) continue;
	if (abs(genp->pdgId()) != 13) continue;
	// W mother
	if (abs(genp->mother()->pdgId()) != 24) continue;

	float dr = ROOT::Math::VectorUtil::DeltaR(lv,genp->p4());
	if (dr < 0.1) {
	  match = true;
	  break;
	}
      } // loop over genps
      if (!match) continue;
    } // do gen match

    // basic dz cut to remove muons from large z
    bool pass_dz = bool(fabs(muon->muonBestTrack()->dz(firstGoodVertex->position())) < 0.5);
    if (!pass_dz) continue;

    // min pt cut
    if (muon->pt() < offPt_) continue;

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

  if (!muons_good.size()) return accept;

  // loop again on duplicate-cleaned muon collection
  muonIndex = 0;
  muons_end = muons_good.end();  // Iterator
  for ( MuonCollection::const_iterator muon = muons_good.begin(); muon != muons_end; ++muon, ++muonIndex ) {
    LorentzVector lv(muon->p4());

    //    bool pass_loose = muon::isLooseMuon(*muon);
    bool pass_tight = muon::isTightMuon(*muon,*firstGoodVertex);
    float trkiso = muon->pfIsolationR03().sumChargedHadronPt;
    // bool pass_trkiso_trig = bool(trkiso/lv.pt() < 0.4);
    float pfiso = muonPFiso(*muon);
    // bool pass_iso_loose = bool(pfiso/lv.pt() < 0.4);
    // bool pass_iso_tight = bool(pfiso/lv.pt() < 0.15);
    TrackRef innerTrack = muon->innerTrack();
    int algo = 0;
    int nhits = 0;
    int npixhits = 0;
    float dxy = -999.;
    if (innerTrack.isNonnull()) {
      npixhits = innerTrack->hitPattern().numberOfValidPixelHits();
      algo = innerTrack->algo();
      nhits = innerTrack->numberOfValidHits();
      dxy = innerTrack->dxy(firstGoodVertex->position());
    }
    float dz_hlt = -999.;
    if (dumpHLTTracks_) {
      const VertexCollection* hltVertexCollection = hltVertexHandle_.product();
      if (hltVertexCollection->size()) {
	dz_hlt = muon->vz() - hltVertexCollection->at(0).position().Z();
      }
    }

    // check match to trigger objects
    bool match = false;
    if (accept) {
      if ( (hlt_mu_lead_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(lv,hlt_mu_lead.lv) < dr_trigmatch) ) match = true;
      else if ( (hlt_mu_subl_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(lv,hlt_mu_subl.lv) < dr_trigmatch) ) match = true;
      else if ( (hlt_mu_third_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(lv,hlt_mu_third.lv) < dr_trigmatch) ) match = true;
    }

    // pt ordering
    if ( ((off_mu_lead_idx == -1) || (lv.pt() > off_mu_lead.lv.pt())) && (lv.pt() > offPt_) ) {
      if (off_mu_lead_idx != -1) {
	off_mu_subl_idx = off_mu_lead_idx;
        off_mu_subl_trigmatch = match;
	off_mu_subl = off_mu_lead;
      }

      off_mu_lead_idx = (int)muonIndex;
      off_mu_lead_trigmatch = match;
      off_mu_lead.lv = lv;
      off_mu_lead.trkiso = trkiso;
      off_mu_lead.pfiso = pfiso;
      off_mu_lead.vz = muon->vz();
      off_mu_lead.charge = muon->charge();
      off_mu_lead.npixhits = npixhits;
      off_mu_lead.algo = algo;
      off_mu_lead.nhits = nhits;
      off_mu_lead.dxy = dxy;
      off_mu_lead.dz_hlt = dz_hlt;
    } else if ( ((off_mu_subl_idx == -1) || (lv.pt() > off_mu_subl.lv.pt())) && (lv.pt() > offPt_) ) {
      off_mu_subl_idx = (int)muonIndex;
      off_mu_subl_trigmatch = match;
      off_mu_subl.lv = lv;
      off_mu_subl.trkiso = trkiso;
      off_mu_subl.pfiso = pfiso;
      off_mu_subl.vz = muon->vz();
      off_mu_subl.charge = muon->charge();
      off_mu_subl.npixhits = npixhits;
      off_mu_subl.algo = algo;
      off_mu_subl.nhits = nhits;
      off_mu_subl.dxy = dxy;
      off_mu_subl.dz_hlt = dz_hlt;
    }

    if (pass_tight) {
      // pt ordering
      if ( (((off_mu_tight_lead_idx == -1) || (lv.pt() > off_mu_tight_lead.lv.pt())))  && (lv.pt() > offPt_) ) {
	if (off_mu_tight_lead_idx != -1) {
	  off_mu_tight_subl_idx = off_mu_tight_lead_idx;
          off_mu_tight_subl_trigmatch = match;
	  off_mu_tight_subl = off_mu_tight_lead;
	}

	off_mu_tight_lead_idx = (int)muonIndex;
        off_mu_tight_lead_trigmatch = match;
	off_mu_tight_lead.lv = lv;
	off_mu_tight_lead.trkiso = trkiso;
	off_mu_tight_lead.pfiso = pfiso;
	off_mu_tight_lead.vz = muon->vz();
	off_mu_tight_lead.charge = muon->charge();
        off_mu_tight_lead.npixhits = npixhits;
        off_mu_tight_lead.algo = algo;
        off_mu_tight_lead.nhits = nhits;
        off_mu_tight_lead.dxy = dxy;
        off_mu_tight_lead.dz_hlt = dz_hlt;
      } else if ( ((off_mu_tight_subl_idx == -1) || (lv.pt() > off_mu_tight_subl.lv.pt()))  && (lv.pt() > offPt_) ) {
	off_mu_tight_subl_idx = (int)muonIndex;
        off_mu_tight_subl_trigmatch = match;
	off_mu_tight_subl.lv = lv;
	off_mu_tight_subl.trkiso = trkiso;
	off_mu_tight_subl.pfiso = pfiso;
	off_mu_tight_subl.vz = muon->vz();
        off_mu_tight_subl.charge = muon->charge();
        off_mu_tight_subl.npixhits = npixhits;
        off_mu_tight_subl.algo = algo;
        off_mu_tight_subl.nhits = nhits;
        off_mu_tight_subl.dxy = dxy;
        off_mu_tight_subl.dz_hlt = dz_hlt;
      }
    }


  } // loop on duplicate-cleaned muons

  // dummy line to prevent errors..
  if (off_mu_tight_subl_trigmatch && off_mu_lead_trigmatch && off_mu_subl_trigmatch) {}

  //-------------------------------------
  //   reco lepton plots
  //-------------------------------------

  // no offline sel
  fillHists(off_mu_lead,off_mu_subl,triggerShort,false);
  // offline tight
  if (off_mu_tight_lead_idx >= 0) {
    fillHists(off_mu_tight_lead,off_mu_tight_subl,triggerShort+"_tight",false);
    if (off_mu_tight_lead_trigmatch) {
      fillHists(off_mu_tight_lead,off_mu_tight_subl,triggerShort+"_tight_match",false);
    } else {
      fillHists(off_mu_tight_lead,off_mu_tight_subl,triggerShort+"_tight_nomatch",false);
    }
  } // tight
  // fillHistsRecoHLT(off_mu_tight_lead,off_mu_tight_subl,hlt_mu_lead,hlt_mu_subl,triggerShort+"_tight");

  // printout for cases where trigger is inefficient
  if ( verbose_ && (triggerEnum == tkmu) && (off_mu_tight_lead_idx >= 0) && !off_mu_tight_lead_trigmatch ) {
    cout << "++ HLT inefficiency, run: " << iEvent.id().run()
	 << ", event: " << iEvent.id().event() << endl
	 << "     off mu pt eta phi: " << off_mu_tight_lead.lv.pt() << " " << off_mu_tight_lead.lv.eta() << " " << off_mu_tight_lead.lv.phi() 
	 << ", iso: " << off_mu_tight_lead.trkiso << ", nhits: " << off_mu_tight_lead.nhits << ", algo: " << off_mu_tight_lead.algo 
	 << ", npixhits: " << off_mu_tight_lead.npixhits << ", vz: " << off_mu_tight_lead.vz  << ", dxy: " << off_mu_tight_lead.dxy << endl;
    // loop over HLT tracks near offline muon and print info
    if (dumpHLTTracks_) {
      edm::Handle<reco::TrackCollection> hltTracksHandle = hltTracksHandle_;
      const VertexCollection* hltVertexCollection = hltVertexHandle_.product();

      cout << "   ++++ online vertex z: ";
      if (hltVertexCollection->size()) {
	cout << hltVertexCollection->at(0).position().Z() << endl;
      } else {
	cout << "none found!!!" << endl;
      }
      cout << "   ++++ hlt tracks near leading offline muon:" << endl;
      // find hlt tracks near offline muons
      TrackCollection::const_iterator hlt_trks_end = hltTracksHandle->end();  // Iterator
      //	for ( TrackCollection::const_iterator hlt_trk = hltTracksHandle->begin(); hlt_trk != hlt_trks_end; ++hlt_trk ) {
      for ( unsigned int hlt_trk_idx = 0; hlt_trk_idx < hltTracksHandle->size(); ++hlt_trk_idx ) {
	reco::TrackRef hlt_trk_ref(hltTracksHandle,hlt_trk_idx);
	float dR = ROOT::Math::VectorUtil::DeltaR(off_mu_tight_lead.lv,hlt_trk_ref->momentum());
	if (dR > 0.3) continue;
	//	int vtx = trackVertex(*hltVertexCollection,hlt_trk_ref);
	float dz = hlt_mu_lead.vz - hlt_trk_ref->vz();
	//	  if (fabs(dz) > 1.0) continue;
	if ((fabs(dz) < 0.2) && (dR < 0.01)) cout << "    * ";
	else cout << "      ";
	cout << "pt: " << (*hlt_trk_ref).pt() << ", eta: " << (*hlt_trk_ref).eta() << ", phi: " << (*hlt_trk_ref).phi()
	     << ", vz: " << (*hlt_trk_ref).vz() // << ", vtx: " << vtx
	     << ", nhits: " << (*hlt_trk_ref).numberOfValidHits() << ", algo: " << (*hlt_trk_ref).algo() 
	     << ", dR: " << dR << endl;
      } // loop on hlt tracks

    } // if dumpHLTTracks

  } // if verbose etc

  return accept;
}

//____________________________________________________________________________
void SingleMuTrigAnalyzerRECO::bookHists(edm::Service<TFileService>& fs, const std::string& suffix, bool hlt) {

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  std::string hlt_suf("_hlt");

  //  std::cout << "SingleMuTrigAnalyzerRECO::bookHists: called with suffix: " << suffix << ", isHLT: " << hlt << std::endl;

  // if (hlt) {
  //   hists_1d_["h_lead_pt"+suf+hlt_suf] = fs->make<TH1F>(Form("h_lead_pt%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Leading p_{T} [GeV]" , 100 , 0. , 100. );
  //   hists_1d_["h_subl_pt"+suf+hlt_suf] = fs->make<TH1F>(Form("h_subl_pt%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Subleading p_{T} [GeV]" , 100 , 0. , 100. );
  //   hists_1d_["h_lead_eta"+suf+hlt_suf] = fs->make<TH1F>(Form("h_lead_eta%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Leading #eta" , 100 , -3. , 3. );
  //   hists_1d_["h_subl_eta"+suf+hlt_suf] = fs->make<TH1F>(Form("h_subl_eta%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Subleading #eta" , 100 , -3. , 3. );
  //   hists_1d_["h_mll"+suf+hlt_suf] = fs->make<TH1F>(Form("h_mll%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT M_{ll} [GeV]" , 150 , 0. , 150. );
  //   hists_1d_["h_dr"+suf+hlt_suf] = fs->make<TH1F>(Form("h_dr%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT #DeltaR" , 600 , 0. , 6. );

  // hists_1d_["h_lead_abstrkiso"+suf+hlt_suf] = fs->make<TH1F>(Form("h_lead_abstrkiso%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Leading trkiso [GeV]" , 200 , 0. , 10. );
  // hists_1d_["h_subl_abstrkiso"+suf+hlt_suf] = fs->make<TH1F>(Form("h_subl_abstrkiso%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Subleading trkiso [GeV]" , 200 , 0. , 10. );
  // hists_1d_["h_lead_reltrkiso"+suf+hlt_suf] = fs->make<TH1F>(Form("h_lead_reltrkiso%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Leading trkiso / p_{T}" , 200 , 0. , 2. );
  // hists_1d_["h_subl_reltrkiso"+suf+hlt_suf] = fs->make<TH1F>(Form("h_subl_reltrkiso%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Subleading trkiso / p_{T}" , 200 , 0. , 2. );

  // }

  hists_1d_["h_nlep_off"+suf] = fs->make<TH1F>(Form("h_nlep_off%s",suf.c_str()) , "; Offline N(lep)" , 3 , -0.5 , 2.5 );

  hists_1d_["h_lead_pt"+suf] = fs->make<TH1F>(Form("h_lead_pt%s",suf.c_str()) , "; Leading p_{T} [GeV]" , 100 , 0. , 100. );
  hists_1d_["h_subl_pt"+suf] = fs->make<TH1F>(Form("h_subl_pt%s",suf.c_str()) , "; Subleading p_{T} [GeV]" , 100 , 0. , 100. );
  hists_1d_["h_lead_eta"+suf] = fs->make<TH1F>(Form("h_lead_eta%s",suf.c_str()) , "; Leading #eta" , 100 , -3. , 3. );
  hists_1d_["h_subl_eta"+suf] = fs->make<TH1F>(Form("h_subl_eta%s",suf.c_str()) , "; Subleading #eta" , 100 , -3. , 3. );
  hists_1d_["h_lead_npixhits"+suf] = fs->make<TH1F>(Form("h_lead_npixhits%s",suf.c_str()) , "; Leading N(pix hits)" , 5 , -0.5 , 4.5 );
  hists_1d_["h_subl_npixhits"+suf] = fs->make<TH1F>(Form("h_subl_npixhits%s",suf.c_str()) , "; Subleading N(pix hits)" , 5 , -0.5 , 4.5 );
  hists_1d_["h_lead_algo"+suf] = fs->make<TH1F>(Form("h_lead_algo%s",suf.c_str()) , "; Leading track algo" , 15 , -0.5 , 14.5 );
  hists_1d_["h_subl_algo"+suf] = fs->make<TH1F>(Form("h_subl_algo%s",suf.c_str()) , "; Subleading track algo" , 15 , -0.5 , 14.5 );
  hists_1d_["h_lead_nhits"+suf] = fs->make<TH1F>(Form("h_lead_nhits%s",suf.c_str()) , "; Leading N(hits)" , 25 , -0.5 , 24.5 );
  hists_1d_["h_subl_nhits"+suf] = fs->make<TH1F>(Form("h_subl_nhits%s",suf.c_str()) , "; Subleading N(hits)" , 25 , -0.5 , 24.5 );
  hists_1d_["h_lead_dxy"+suf] = fs->make<TH1F>(Form("h_lead_dxy%s",suf.c_str()) , "; Leading dxy wrt vertex" , 100 , -0.2 , 0.2 );
  hists_1d_["h_subl_dxy"+suf] = fs->make<TH1F>(Form("h_subl_dxy%s",suf.c_str()) , "; Subleading dxy wrt vertex" , 100 , -0.2 , 0.2 );
  hists_1d_["h_lead_dz_hlt"+suf] = fs->make<TH1F>(Form("h_lead_dz_hlt%s",suf.c_str()) , "; Leading dz wrt hlt vertex" , 100 , -20 , 20 );
  hists_1d_["h_subl_dz_hlt"+suf] = fs->make<TH1F>(Form("h_subl_dz_hlt%s",suf.c_str()) , "; Subleading dz wrt hlt vertex" , 100 , -20 , 20 );
  hists_1d_["h_mll"+suf] = fs->make<TH1F>(Form("h_mll%s",suf.c_str()) , "; M_{ll} [GeV]" , 150 , 0. , 150. );
  hists_1d_["h_dr"+suf] = fs->make<TH1F>(Form("h_dr%s",suf.c_str()) , "; #DeltaR" , 600 , 0. , 6. );

  // hists_1d_["h_lead_pt_os"+suf] = fs->make<TH1F>(Form("h_lead_pt_os%s",suf.c_str()) , "; Leading p_{T} [GeV]" , 100 , 0. , 100. );
  // hists_1d_["h_subl_pt_os"+suf] = fs->make<TH1F>(Form("h_subl_pt_os%s",suf.c_str()) , "; Subleading p_{T} [GeV]" , 100 , 0. , 100. );
  // hists_1d_["h_lead_eta_os"+suf] = fs->make<TH1F>(Form("h_lead_eta_os%s",suf.c_str()) , "; Leading #eta" , 100 , -3. , 3. );
  // hists_1d_["h_subl_eta_os"+suf] = fs->make<TH1F>(Form("h_subl_eta_os%s",suf.c_str()) , "; Subleading #eta" , 100 , -3. , 3. );
  // hists_1d_["h_mll_os"+suf] = fs->make<TH1F>(Form("h_mll_os%s",suf.c_str()) , "; M_{ll} [GeV]" , 150 , 0. , 150. );
  // hists_1d_["h_dr_os"+suf] = fs->make<TH1F>(Form("h_dr_os%s",suf.c_str()) , "; #DeltaR" , 600 , 0. , 6. );

  // hists_1d_["h_lead_pt_ss"+suf] = fs->make<TH1F>(Form("h_lead_pt_ss%s",suf.c_str()) , "; Leading p_{T} [GeV]" , 100 , 0. , 100. );
  // hists_1d_["h_subl_pt_ss"+suf] = fs->make<TH1F>(Form("h_subl_pt_ss%s",suf.c_str()) , "; Subleading p_{T} [GeV]" , 100 , 0. , 100. );
  // hists_1d_["h_lead_eta_ss"+suf] = fs->make<TH1F>(Form("h_lead_eta_ss%s",suf.c_str()) , "; Leading #eta" , 100 , -3. , 3. );
  // hists_1d_["h_subl_eta_ss"+suf] = fs->make<TH1F>(Form("h_subl_eta_ss%s",suf.c_str()) , "; Subleading #eta" , 100 , -3. , 3. );
  // hists_1d_["h_mll_ss"+suf] = fs->make<TH1F>(Form("h_mll_ss%s",suf.c_str()) , "; M_{ll} [GeV]" , 150 , 0. , 150. );
  // hists_1d_["h_dr_ss"+suf] = fs->make<TH1F>(Form("h_dr_ss%s",suf.c_str()) , "; #DeltaR" , 600 , 0. , 6. );

  // hists_1d_["h_lead_abstrkiso"+suf] = fs->make<TH1F>(Form("h_lead_abstrkiso%s",suf.c_str()) , "; Leading trkiso [GeV]" , 200 , 0. , 10. );
  // hists_1d_["h_subl_abstrkiso"+suf] = fs->make<TH1F>(Form("h_subl_abstrkiso%s",suf.c_str()) , "; Subleading trkiso [GeV]" , 200 , 0. , 10. );
  // hists_1d_["h_lead_reltrkiso"+suf] = fs->make<TH1F>(Form("h_lead_reltrkiso%s",suf.c_str()) , "; Leading trkiso / p_{T}" , 200 , 0. , 2. );
  // hists_1d_["h_subl_reltrkiso"+suf] = fs->make<TH1F>(Form("h_subl_reltrkiso%s",suf.c_str()) , "; Subleading trkiso / p_{T}" , 200 , 0. , 2. );

  // hists_1d_["h_lead_abspfiso"+suf] = fs->make<TH1F>(Form("h_lead_abspfiso%s",suf.c_str()) , "; Leading pfiso [GeV]" , 200 , 0. , 10. );
  // hists_1d_["h_subl_abspfiso"+suf] = fs->make<TH1F>(Form("h_subl_abspfiso%s",suf.c_str()) , "; Subleading pfiso [GeV]" , 200 , 0. , 10. );
  // hists_1d_["h_lead_relpfiso"+suf] = fs->make<TH1F>(Form("h_lead_relpfiso%s",suf.c_str()) , "; Leading pfiso / p_{T}" , 200 , 0. , 2. );
  // hists_1d_["h_subl_relpfiso"+suf] = fs->make<TH1F>(Form("h_subl_relpfiso%s",suf.c_str()) , "; Subleading pfiso / p_{T}" , 200 , 0. , 2. );

  // hists_1d_["h_lead_offhlt_dpt"+suf] = fs->make<TH1F>(Form("h_lead_offhlt_dpt%s",suf.c_str()) , "; Leading (p_{T}^{off} - p_{T}^{HLT}) / p_{T}^{off}" , 250 , -5. , 5. );
  // hists_1d_["h_subl_offhlt_dpt"+suf] = fs->make<TH1F>(Form("h_subl_offhlt_dpt%s",suf.c_str()) , "; Subleading (p_{T}^{off} - p_{T}^{HLT}) / p_{T}^{off}" , 250 , -5. , 5. );
  // hists_1d_["h_lead_offhlt_dr"+suf] = fs->make<TH1F>(Form("h_lead_offhlt_dr%s",suf.c_str()) , "; Leading #DeltaR(off,HLT)" , 600 , 0. , 6. );
  // hists_1d_["h_subl_offhlt_dr"+suf] = fs->make<TH1F>(Form("h_subl_offhlt_dr%s",suf.c_str()) , "; Subleading #DeltaR(off,HLT)" , 600 , 0. , 6. );

  // hists_2d_["h_lead_abstrkiso_hlt_vs_off"+suf] = fs->make<TH2F>(Form("h_lead_abstrkiso_hlt_vs_off%s",suf.c_str()) , "; Leading offline trkiso [GeV]; Leading HLT trkiso [GeV]" , 50 , 0. , 10. , 50 , 0. , 10. );
  // hists_2d_["h_subl_abstrkiso_hlt_vs_off"+suf] = fs->make<TH2F>(Form("h_subl_abstrkiso_hlt_vs_off%s",suf.c_str()) , "; Subleading offline trkiso [GeV]; Subleading HLT trkiso [GeV]" , 50 , 0. , 10. , 50 , 0. , 10. );
  // hists_2d_["h_lead_reltrkiso_hlt_vs_off"+suf] = fs->make<TH2F>(Form("h_lead_reltrkiso_hlt_vs_off%s",suf.c_str()) , "; Leading offline trkiso / p_{T}; Leading HLT trkiso / p_{T}" , 50 , 0. , 2. , 50 , 0. , 2. );
  // hists_2d_["h_subl_reltrkiso_hlt_vs_off"+suf] = fs->make<TH2F>(Form("h_subl_reltrkiso_hlt_vs_off%s",suf.c_str()) , "; Subleading offline trkiso / p_{T}; Subleading HLT trkiso / p_{T}" , 50 , 0. , 2. , 50 , 0. , 2. );


  return;
}

//____________________________________________________________________________
void SingleMuTrigAnalyzerRECO::fillHists(const StudyLepton& lead, const StudyLepton& subl, const std::string& suffix, bool isHLT) {

  // if (lead.lv.pt() <= 0.) {
  //   std::cout << "SingleMuTrigAnalyzerRECO::fillHists: invalid lead pt: " << lead.lv.pt() << std::endl;
  //   return;
  // }

  //  std::cout << "SingleMuTrigAnalyzerRECO::fillHists: called with suffix: " << suffix << ", isHLT: " << isHLT << std::endl;

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  std::string hlt_suf("");
  if (isHLT) hlt_suf = "_hlt";

  int nlep = 0;
  if (lead.lv.pt() > 0) ++nlep;
  if (subl.lv.pt() > 0) ++nlep;
  if (!isHLT) hists_1d_["h_nlep_off"+suf+hlt_suf]->Fill(nlep);
  if (nlep == 0) return;

  if (lead.lv.pt() > 0) {
    hists_1d_["h_lead_pt"+suf+hlt_suf]->Fill(lead.lv.pt());
    hists_1d_["h_lead_eta"+suf+hlt_suf]->Fill(lead.lv.eta());
    hists_1d_["h_lead_npixhits"+suf+hlt_suf]->Fill(lead.npixhits);
    hists_1d_["h_lead_algo"+suf+hlt_suf]->Fill(lead.algo);
    hists_1d_["h_lead_nhits"+suf+hlt_suf]->Fill(lead.nhits);
    hists_1d_["h_lead_dxy"+suf+hlt_suf]->Fill(lead.dxy);
    hists_1d_["h_lead_dz_hlt"+suf+hlt_suf]->Fill(lead.dz_hlt);
    // if (lead.type == 13) {
    //   hists_1d_["h_lead_abstrkiso"+suf+hlt_suf]->Fill(lead.trkiso);
    //   hists_1d_["h_lead_reltrkiso"+suf+hlt_suf]->Fill(lead.trkiso/lead.lv.pt());
    //   if (!isHLT) {
    // 	hists_1d_["h_lead_abspfiso"+suf+hlt_suf]->Fill(lead.pfiso);
    // 	hists_1d_["h_lead_relpfiso"+suf+hlt_suf]->Fill(lead.pfiso/lead.lv.pt());
    //   }
    // }
  } // valid lead pt

  if (subl.lv.pt() > 0) {
    hists_1d_["h_subl_pt"+suf+hlt_suf]->Fill(subl.lv.pt());
    hists_1d_["h_subl_eta"+suf+hlt_suf]->Fill(subl.lv.eta());
    hists_1d_["h_subl_npixhits"+suf+hlt_suf]->Fill(subl.npixhits);
    hists_1d_["h_subl_algo"+suf+hlt_suf]->Fill(subl.algo);
    hists_1d_["h_subl_nhits"+suf+hlt_suf]->Fill(subl.nhits);
    hists_1d_["h_subl_dxy"+suf+hlt_suf]->Fill(subl.dxy);
    hists_1d_["h_subl_dz_hlt"+suf+hlt_suf]->Fill(subl.dz_hlt);

    if (lead.lv.pt() > 0) {
      LorentzVector dilep = lead.lv+subl.lv;
      hists_1d_["h_mll"+suf+hlt_suf]->Fill(dilep.M());
      float dr = ROOT::Math::VectorUtil::DeltaR(lead.lv,subl.lv);
      hists_1d_["h_dr"+suf+hlt_suf]->Fill(dr);
      // if (!isHLT && (lead.charge == subl.charge)) {
      // 	hists_1d_["h_lead_pt_ss"+suf+hlt_suf]->Fill(lead.lv.pt());
      // 	hists_1d_["h_lead_eta_ss"+suf+hlt_suf]->Fill(lead.lv.eta());
      // 	hists_1d_["h_subl_pt_ss"+suf+hlt_suf]->Fill(subl.lv.pt());
      // 	hists_1d_["h_subl_eta_ss"+suf+hlt_suf]->Fill(subl.lv.eta());
      // 	hists_1d_["h_mll_ss"+suf+hlt_suf]->Fill(dilep.M());
      // 	hists_1d_["h_dr_ss"+suf+hlt_suf]->Fill(dr);
      // } else if (!isHLT) {
      // 	hists_1d_["h_lead_pt_os"+suf+hlt_suf]->Fill(lead.lv.pt());
      // 	hists_1d_["h_lead_eta_os"+suf+hlt_suf]->Fill(lead.lv.eta());
      // 	hists_1d_["h_subl_pt_os"+suf+hlt_suf]->Fill(subl.lv.pt());
      // 	hists_1d_["h_subl_eta_os"+suf+hlt_suf]->Fill(subl.lv.eta());
      // 	hists_1d_["h_mll_os"+suf+hlt_suf]->Fill(dilep.M());
      // 	hists_1d_["h_dr_os"+suf+hlt_suf]->Fill(dr);
      // }
    }
    // if (subl.type == 13) {
    //   hists_1d_["h_subl_abstrkiso"+suf+hlt_suf]->Fill(subl.trkiso);
    //   hists_1d_["h_subl_reltrkiso"+suf+hlt_suf]->Fill(subl.trkiso/subl.lv.pt());
    //   if (!isHLT) {
    // 	hists_1d_["h_subl_abspfiso"+suf+hlt_suf]->Fill(subl.pfiso);
    // 	hists_1d_["h_subl_relpfiso"+suf+hlt_suf]->Fill(subl.pfiso/subl.lv.pt());
    //   }
    // }
  } // valid subl pt

  return;
}

//____________________________________________________________________________
void SingleMuTrigAnalyzerRECO::fillHistsRecoHLT(const StudyLepton& off_lead, const StudyLepton& off_subl, const StudyLepton& hlt_lead, const StudyLepton& hlt_subl, const std::string& suffix) {

  // if (off_lead.pt() <= 0.) {
  //   std::cout << "SingleMuTrigAnalyzerRECO::fillHists: invalid lead pt: " << off_lead.pt() << std::endl;
  //   return;
  // }

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  int nlep = 0;
  if (off_lead.lv.pt() > 0) ++nlep;
  if (off_subl.lv.pt() > 0) ++nlep;
  if (nlep == 0) return;

  // if (off_lead.lv.pt() > 0.) {
  //   bool match_lead = false;
  //   if (hlt_lead.lv.pt() > 0. && (off_lead.type == hlt_lead.type)) {
  //     float dr_hlt_lead = ROOT::Math::VectorUtil::DeltaR(off_lead.lv,hlt_lead.lv);
  //     if (dr_hlt_lead < 0.2) {
  // 	hists_1d_["h_lead_offhlt_dpt"+suf]->Fill( (off_lead.lv.pt() - hlt_lead.lv.pt()) / off_lead.lv.pt() );
  // 	hists_1d_["h_lead_offhlt_dr"+suf]->Fill(dr_hlt_lead);
  // 	if (off_lead.type == 13) {
  // 	  hists_2d_["h_lead_abstrkiso_hlt_vs_off"+suf]->Fill(off_lead.trkiso,hlt_lead.trkiso);
  // 	  hists_2d_["h_lead_reltrkiso_hlt_vs_off"+suf]->Fill(off_lead.trkiso/off_lead.lv.pt(),hlt_lead.trkiso/hlt_lead.lv.pt());
  // 	}
  // 	match_lead = true;
  //     }
  //   } // valid hlt_lead
  //   // if we didn't match the lead, and the leptons are the same type, then check subleading
  //   else if (!match_lead && hlt_subl.lv.pt() > 0. && (off_lead.type == hlt_subl.type)) {
  //     float dr_hlt_subl = ROOT::Math::VectorUtil::DeltaR(off_lead.lv,hlt_subl.lv);
  //     if (dr_hlt_subl < 0.2) {
  // 	hists_1d_["h_lead_offhlt_dpt"+suf]->Fill( (off_lead.lv.pt() - hlt_subl.lv.pt()) / off_lead.lv.pt() );
  // 	hists_1d_["h_lead_offhlt_dr"+suf]->Fill(dr_hlt_subl);
  // 	if (off_lead.type == 13) {
  // 	  hists_2d_["h_lead_abstrkiso_hlt_vs_off"+suf]->Fill(off_lead.trkiso,hlt_subl.trkiso);
  // 	  hists_2d_["h_lead_reltrkiso_hlt_vs_off"+suf]->Fill(off_lead.trkiso/off_lead.lv.pt(),hlt_subl.trkiso/hlt_subl.lv.pt());
  // 	}
  //     }
  //   } // hlt_subl
  // } // valid off_lead

  // if (off_subl.lv.pt() > 0.) {
  //   bool match_lead = false;
  //   // if this isn't an emu trigger, check the leading lepton at the HLT
  //   if (hlt_lead.lv.pt() > 0. && (off_subl.type == hlt_lead.type)) {
  //     float dr_hlt_lead = ROOT::Math::VectorUtil::DeltaR(off_subl.lv,hlt_lead.lv);
  //     if (dr_hlt_lead < 0.2) {
  // 	hists_1d_["h_subl_offhlt_dpt"+suf]->Fill( (off_subl.lv.pt() - hlt_lead.lv.pt()) / off_subl.lv.pt() );
  // 	hists_1d_["h_subl_offhlt_dr"+suf]->Fill(dr_hlt_lead);
  // 	if (off_subl.type == 13) {
  // 	  hists_2d_["h_subl_abstrkiso_hlt_vs_off"+suf]->Fill(off_subl.trkiso,hlt_lead.trkiso);
  // 	  hists_2d_["h_subl_reltrkiso_hlt_vs_off"+suf]->Fill(off_subl.trkiso/off_subl.lv.pt(),hlt_lead.trkiso/hlt_lead.lv.pt());
  // 	}
  // 	match_lead = true;
  //     }
  //   } // valid hlt_lead
  //   else if (!match_lead && hlt_subl.lv.pt() > 0. && (off_subl.type == hlt_subl.type)) {
  //     float dr_hlt_subl = ROOT::Math::VectorUtil::DeltaR(off_subl.lv,hlt_subl.lv);
  //     if (dr_hlt_subl < 0.2) {
  // 	hists_1d_["h_subl_offhlt_dpt"+suf]->Fill( (off_subl.lv.pt() - hlt_subl.lv.pt()) / off_subl.lv.pt() );
  // 	hists_1d_["h_subl_offhlt_dr"+suf]->Fill(dr_hlt_subl);
  // 	if (off_subl.type == 13) {
  // 	  hists_2d_["h_subl_abstrkiso_hlt_vs_off"+suf]->Fill(off_subl.trkiso,hlt_subl.trkiso);
  // 	  hists_2d_["h_subl_reltrkiso_hlt_vs_off"+suf]->Fill(off_subl.trkiso/off_subl.lv.pt(),hlt_subl.trkiso/hlt_subl.lv.pt());
  // 	}
  //     }
  //   } // hlt_subl
  // } // valid off_subl

  return;
}

//____________________________________________________________________________
void SingleMuTrigAnalyzerRECO::bookHistsGen(edm::Service<TFileService>& fs, const std::string& suffix) {

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  hists_1d_["h_pt"+suf] = fs->make<TH1F>(Form("h_pt%s",suf.c_str()) , "; p_{T} [GeV]" , 100 , 0. , 100. );
  hists_1d_["h_eta"+suf] = fs->make<TH1F>(Form("h_eta%s",suf.c_str()) , "; #eta" , 100 , -3. , 3. );

  return;
}

//____________________________________________________________________________
void SingleMuTrigAnalyzerRECO::fillHistsGen(const StudyLepton& mu, const std::string& suffix) {

  //  std::cout << "SingleMuTrigAnalyzerRECO::fillHists: called with suffix: " << suffix << ", isHLT: " << isHLT << std::endl;

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  hists_1d_["h_pt"+suf]->Fill(mu.lv.pt());
  hists_1d_["h_eta"+suf]->Fill(mu.lv.eta());

  return;
}

//____________________________________________________________________________
float SingleMuTrigAnalyzerRECO::muonPFiso(const reco::Muon& muon) {

  reco::MuonPFIsolation pfStructR03 = muon.pfIsolationR03();
  float chiso = pfStructR03.sumChargedHadronPt;
  float nhiso = pfStructR03.sumNeutralHadronEt;
  float emiso = pfStructR03.sumPhotonEt;
  float deltaBeta = pfStructR03.sumPUPt;

  //  float absiso = chiso + nhiso + emiso;
  float absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return absiso;

}

//____________________________________________________________________________
// returns best match vertex.  Taken from PFPileUpAlgo, chargedHadronVertex
int SingleMuTrigAnalyzerRECO::chargedHadronVertex( const reco::VertexCollection& vertices, const reco::PFCandidate& pfcand ) const {
  if (pfcand.trackRef().isNull()) return -2;
  return trackVertex(vertices,pfcand.trackRef());
}

//____________________________________________________________________________
// returns best match vertex.  Taken from PFPileUpAlgo, chargedHadronVertex
int SingleMuTrigAnalyzerRECO::trackVertex( const reco::VertexCollection& vertices, const reco::TrackRef& trackRef ) const {

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


