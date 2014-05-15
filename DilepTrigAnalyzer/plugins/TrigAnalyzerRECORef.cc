/** \class TrigAnalyzerRECORef
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
#include "TrigAnalyzer/DilepTrigAnalyzer/interface/TrigAnalyzerRECORef.h"

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

#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"

// ROOT includes
#include "Math/VectorUtil.h"

#include <cassert>

using namespace reco;
using namespace edm;

//
// constructors and destructor
//
//____________________________________________________________________________
TrigAnalyzerRECORef::TrigAnalyzerRECORef(const edm::ParameterSet& ps) : 
  processName_(ps.getParameter<std::string>("processName")),
  triggerNames_(ps.getParameter<std::vector<std::string> >("triggerNames")),
  triggerNamesShort_(ps.getParameter<std::vector<std::string> >("triggerNamesShort")),
  triggerResultsTag_(ps.getParameter<edm::InputTag>("triggerResults")),
  triggerEventWithRefsTag_(ps.getParameter<edm::InputTag>("triggerEventWithRefs")),
  electronsInputTag_(ps.getParameter<edm::InputTag>("electronsInputTag")),
  muonsInputTag_(ps.getParameter<edm::InputTag>("muonsInputTag")),
  vtxInputTag_(ps.getParameter<edm::InputTag>("vtxInputTag")),
  offLeadPt_(ps.getParameter<double>("offLeadPt")),
  offSublPt_(ps.getParameter<double>("offSublPt")),
  doOffGenMatch_(ps.getParameter<bool>("doOffGenMatch")),
  genParticlesTag_(ps.getParameter<edm::InputTag>("genParticles")),
  verbose_(ps.getParameter<bool>("verbose"))
{
  using namespace std;
  using namespace edm;

  cout << "TrigAnalyzerRECORef configuration: " << endl
       << "   ProcessName = " << processName_ << endl
       << "   TriggerNames = "; 
  for (unsigned int i=0; i<triggerNames_.size(); ++i) {
    cout << triggerNames_.at(i) << ", ";
  }
  cout << endl 
       << "   TriggerNamesShort = "; 
  for (unsigned int i=0; i<triggerNamesShort_.size(); ++i) {
    cout << triggerNamesShort_.at(i) << ", ";
  }
  cout << endl
       << "   TriggerResultsTag = " << triggerResultsTag_.encode() << endl
       << "   TriggerEventWithRefsTag = " << triggerEventWithRefsTag_.encode() << endl
       << "   ElectronsInputTag = " << electronsInputTag_.encode() << endl
       << "   MuonsInputTag = " << muonsInputTag_.encode() << endl
       << "   VtxInputTag = " << vtxInputTag_.encode() << endl
       << "   OffLeadPt = " << offLeadPt_ << endl
       << "   OffSublPt = " << offSublPt_ << endl
       << "   DoOffGenMatch = " << doOffGenMatch_ << endl
       << "   GenParticlesTag = " << genParticlesTag_.encode() << endl
       << "   Verbose = " << verbose_ << endl;

  if (triggerNames_.size() == 0) {
    cout << "ERROR: must specify at least one trigger!!" << endl;
  } else if (triggerNames_.size() != triggerNamesShort_.size()) {
    cout << "ERROR: triggerNames and triggerNamesShort must have same size!!" << endl;
  }

  // histogram setup
  edm::Service<TFileService> fs;
  // unsigned int bitsize = 2**(triggerNames_.size());
  // h_results_ = fs->make<TH1F>("h_results_" , ";Trigger Results" , bitsize , -0.5 , float(bitsize)-0.5 );
  bookHists(fs);
  for (unsigned int itrig=0; itrig < triggerNames_.size(); ++itrig) {
    bookHists(fs,triggerNamesShort_.at(itrig)+"_match");
    bookHists(fs,triggerNamesShort_.at(itrig)+"_match_onZ");
    bookHists(fs,triggerNamesShort_.at(itrig)+"_match_offZ");
    bookHists(fs,triggerNamesShort_.at(itrig)+"_nomatch");
    bookHists(fs,triggerNamesShort_.at(itrig)+"_nomatch_onZ");
    bookHists(fs,triggerNamesShort_.at(itrig)+"_nomatch_offZ");
  }

}

//____________________________________________________________________________
TrigAnalyzerRECORef::~TrigAnalyzerRECORef()
{
}

//
// member functions
//
//____________________________________________________________________________
void
TrigAnalyzerRECORef::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  using namespace std;
  using namespace edm;

  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
      const unsigned int n(hltConfig_.size());
      // check if trigger names in (new) config
      for (unsigned int itrig = 0; itrig < triggerNames_.size(); ++itrig) {
	unsigned int triggerIndex(hltConfig_.triggerIndex(triggerNames_.at(itrig)));
	if (triggerIndex>=n) {
	  cout << "TrigAnalyzerRECORef::analyze:"
	       << " TriggerName " << triggerNames_.at(itrig) 
	       << " not available in (new) config!" << endl;
	}
      } // loop over triggers
    } // if changed
  } else {
    cout << "TrigAnalyzerRECORef::analyze:"
	 << " config extraction failure with process name "
	 << processName_ << endl;
  }

}

//____________________________________________________________________________
// ------------ method called to produce the data  ------------
void
TrigAnalyzerRECORef::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;

  if (verbose_) cout << endl;

  // get event products
  iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
    cout << "TrigAnalyzerRECORef::analyze: Error in getting TriggerResults product from Event!" << endl;
    return;
  }
  iEvent.getByLabel(triggerEventWithRefsTag_,triggerEventWithRefsHandle_);
  if (!triggerEventWithRefsHandle_.isValid()) {
    cout << "TrigAnalyzerRECORef::analyze: Error in getting TriggerEventWithRefs product from Event!" << endl;
    return;
  }
  // sanity check
  assert(triggerResultsHandle_->size()==hltConfig_.size());

  // retrieve necessary containers
  iEvent.getByLabel(vtxInputTag_, vertexHandle_);
  iEvent.getByLabel(electronsInputTag_, elsHandle_);
  iEvent.getByLabel( muonsInputTag_ , musHandle_ );

  if (doOffGenMatch_) {
    iEvent.getByLabel(genParticlesTag_, genParticlesHandle_);
  }

  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------
  //   analyze event
  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------
  
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
    cout << "TrigAnalyzerRECORef::analyzeTrigger: didn't find any good offline vertices!! size: " 
	 << vertexCollection->size() << std::endl;
    return;
  }

  //-------------------------------------
  //   reco muons
  //-------------------------------------

  // 1. select good muons for study
  // 2. determine if event is on/offZ etc
  // 3. loop on triggers to collect trig muons
  // 4. loop back over good muons to check trig decisions and make plots

  MuonCollection muons_good;
  MuonCollection muons_dup;

  // loop first to remove duplicate muons and select good muons for study
  unsigned int muonIndex = 0;
  MuonCollection::const_iterator muons_end = musHandle_->end();  // Iterator
  for ( MuonCollection::const_iterator muon = musHandle_->begin(); muon != muons_end; ++muon, ++muonIndex ) {
    LorentzVector lv(muon->p4());

    //    if (verbose_) cout << " - reco muon: pt: " << lv.pt() << ", eta: " << lv.eta() << ", phi: " << lv.phi(); 

    bool duplicate = false;
    // check if this muon is already flagged as a duplicate
    for ( MuonCollection::const_iterator muon2 = muons_dup.begin(); muon2 != muons_dup.end(); ++muon2 ) {
      // !!!!!!!! this may not work
      if (muon == muon2) {
	duplicate = true;
	break;
      }
    } // loop on found duplicates
    if (duplicate) {
      // if (verbose_) cout << ", DUPLICATE" << std::endl;
      continue;
    }

    // min pt cut
    if (muon->pt() < offSublPt_) {
      // if (verbose_) cout << ", FAILS pt" << std::endl;
      continue;
    }

    // basic dz cut to remove muons from large z
    bool pass_dz = bool(fabs(muon->muonBestTrack()->dz(firstGoodVertex->position())) < 0.5);
    if (!pass_dz) {
      // if (verbose_) cout << ", FAILS dz" << std::endl;
      continue;
    }

    // check tight muon ID, require here
    bool pass_tight_1 = muon::isTightMuon(*muon,*firstGoodVertex);
    if (!pass_tight_1) continue;

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
    if (duplicate) {
      // if (verbose_) cout << ", DUPLICATE" << std::endl;
      continue;
    }

    // if (verbose_) cout << ", passes presel cuts" << std::endl;

    muons_good.push_back(*muon);

  } // loop on reco muons

  // require at least one muon to continue
  if (!muons_good.size()) return;

  // check duplicate-cleaned muon collection to see if any pairs are on the z mass
  muonIndex = 0;
  muons_end = muons_good.end();  // Iterator
  bool onZ = false;
  for ( MuonCollection::const_iterator muon = muons_good.begin(); muon != muons_end; ++muon, ++muonIndex ) {
    LorentzVector lv(muon->p4());

    for ( MuonCollection::const_iterator muon2 = muons_good.begin(); muon2 != muons_end; ++muon2 ) {
      if (muon == muon2) continue;
      LorentzVector lv2(muon2->p4());
      float mass = (lv+lv2).M();
      if (mass > 81. && mass < 101.) {
	onZ = true;
	break;
      }
    } // 2nd loop on muons
  } // loop on good muons

  std::vector<trigger::VRmuon> triggersMuons;
  // loop on triggers to get collections of trigger muons
  for (unsigned int itrig = 0; itrig < triggerNames_.size(); ++itrig) {
    trigger::VRmuon trigMuons = getTrigMuons(iEvent,iSetup,triggerNames_.at(itrig));
    triggersMuons.push_back(trigMuons);
  }

  // loop to compare with triggers and make plots
  muonIndex = 0;
  const float dr_trigmatch = 0.2;
  for ( MuonCollection::const_iterator muon = muons_good.begin(); muon != muons_end; ++muon, ++muonIndex ) {
    LorentzVector lv(muon->p4());

    fillHists(*muon);

    for (unsigned int itrig = 0; itrig < triggerNames_.size(); ++itrig) {
      bool match = false;
      for (unsigned int itmu = 0; itmu < triggersMuons.at(itrig).size(); ++itmu) {
	if (ROOT::Math::VectorUtil::DeltaR(lv,triggersMuons.at(itrig).at(itmu)->p4()) < dr_trigmatch) {
	  match = true;
	  break;
	}
      } // loop over trig muons
      const std::string nameShort(triggerNamesShort_.at(itrig));
      if (match) {
	fillHists(*muon,nameShort+"_match");
	if (onZ) fillHists(*muon,nameShort+"_match_onZ");
	else fillHists(*muon,nameShort+"_match_offZ");
      } else {
	fillHists(*muon,nameShort+"_nomatch");
	if (onZ) fillHists(*muon,nameShort+"_nomatch_onZ");
	else fillHists(*muon,nameShort+"_nomatch_offZ");
      }
    } // loop over trigs

    // bool pass_loose = muon::isLooseMuon(*muon);
    // bool pass_tight = muon::isTightMuon(*muon,*firstGoodVertex);
    // float trkiso = muon->pfIsolationR03().sumChargedHadronPt;
    // bool pass_trkiso_trig = bool(trkiso/lv.pt() < 0.4);
    // float pfiso = muonPFiso(*muon);
    // bool pass_iso_loose = bool(pfiso/lv.pt() < 0.4);
    // bool pass_iso_tight = bool(pfiso/lv.pt() < 0.15);
    // const TrackRef siTrack  = muon->innerTrack();
    // float dxy = -999.;
    // if ( siTrack.isNonnull() && firstGoodVertex != vertexCollection->end() ) {
    //   dxy = siTrack->dxy(firstGoodVertex->position());
    // } 
    // bool pass_dxy = bool(fabs(dxy) < 0.02);

  } // end muon loop for plots

  if (verbose_) cout << endl;

  return;
}

//____________________________________________________________________________

trigger::VRmuon TrigAnalyzerRECORef::getTrigMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& triggerName) {
  
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  if (verbose_) cout << endl;

  trigger::Vids muonIds;
  trigger::VRmuon muonRefs;

  const unsigned int ntrigs(hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));

  // abort on invalid trigger name
  if (triggerIndex>=ntrigs) {
    cout << "TrigAnalyzerRECORef::analyzeTrigger: path "
	 << triggerName << " - not found!" << endl;
    return muonRefs;
  }

  if (verbose_) {
    cout << "TrigAnalyzerRECORef::analyzeTrigger: path "
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

  if (!wasRun || !accept || error) return muonRefs;

  // get trigger objects from last filter

  //------------------------------------
  //  hlt objects
  //------------------------------------

  bool foundMuons = false;

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

    // electronIds_.clear();
    // electronRefs_.clear();
    // triggerEventWithRefsHandle_->getObjects(filterIndex,electronIds_,electronRefs_);
    // const unsigned int nElectrons(electronIds_.size());
    // if (nElectrons>0) {
    //   if (verbose_) { 
    //   cout << "   Electrons: " << nElectrons << "  - the objects: # id pt eta phi vz id key" << endl;
    //   for (unsigned int i=0; i!=nElectrons; ++i) {

    // 	  cout << "   " << i << " " << electronIds_.at(i)
    // 	       << " " << electronRefs_.at(i)->pt()
    // 	       << " " << electronRefs_.at(i)->eta()
    // 	       << " " << electronRefs_.at(i)->phi()
    // 	       << " " << electronRefs_.at(i)->vz()
    // 	       << " " << electronRefs_.at(i).id()
    // 	       << " " << electronRefs_.at(i).key()
    // 	       << endl;
    //   } // electrons loop
    //  }
    // } // if nElectrons>0

    muonIds.clear();
    muonRefs.clear();
    triggerEventWithRefsHandle_->getObjects(filterIndex,muonIds,muonRefs);
    const unsigned int nMuons(muonIds.size());
    if (nMuons>0) {
      foundMuons = true;

      if (verbose_) {
	cout << "   Muons: " << nMuons << ", MuonRefs: " << muonRefs.size()
	     << "  - the objects: # id pt eta phi vz id key eta_trkref phi_trkref" << endl;
	for (unsigned int i=0; i!=nMuons; ++i) {
	  cout << "   " << i
	       << " " << muonIds.at(i)
	       << " " << muonRefs.at(i)->pt()
	       << " " << muonRefs.at(i)->eta()
	       << " " << muonRefs.at(i)->phi()
	       << " " << muonRefs.at(i)->vz()
	       << " " << muonRefs.at(i).id()
	       << " " << muonRefs.at(i).key();
	  cout << endl;
	}
      } // verbose

    } // if muons in TriggerEventWithRefs

    if (foundMuons) break;
  } // backwards loop on modules

  if (accept && !foundMuons) {
    cout << "TrigAnalyzerRECORef::getTrigMuons: no valid trigger leptons!  trigger: " << triggerName << endl;
  }

  return muonRefs;
}

//   if ((hlt_el_lead_idx == -1) && (hlt_mu_lead_idx == -1)) {
//     cout << "TrigAnalyzerRECORef::analyzeTrigger: no valid trigger leptons!" << endl;
//     return true;
//   }

//   // make plots based on trigger path
//   if (ismm) {
//     fillHists(hlt_mu_lead,hlt_mu_subl,triggerShort,true);
//   } else if (isem) {
//     fillHists(hlt_el_lead,hlt_mu_lead,triggerShort,true);
//   } else if (isme) {
//     fillHists(hlt_mu_lead,hlt_el_lead,triggerShort,true);
//   }

//   //-------------------------------------
//   //   pt thresholds - depend on trigger
//   //-------------------------------------

//   float leadPtThresh_e = offLeadPt_;
//   float sublPtThresh_e = offSublPt_;
//   float leadPtThresh_m = offLeadPt_;
//   float sublPtThresh_m = offSublPt_;
//   if (isem) leadPtThresh_m = offSublPt_;
//   else if (isme) leadPtThresh_e = offSublPt_;

//   //-------------------------------------
//   //   reco electrons 
//   //-------------------------------------

//   const float dr_trigmatch = 0.2;

//   StudyLepton off_el_lead;
//   StudyLepton off_el_subl;
//   off_el_lead.type = 11;
//   off_el_subl.type = 11;
//   off_el_lead.isHLT = false;
//   off_el_subl.isHLT = false;
//   int off_el_lead_idx = -1;
//   int off_el_subl_idx = -1;

//   // needed???
//   //  View<GsfElectron> gsfElColl = *(elsHandle_.product());
//   // Handle<GsfElectronCollection> els_coll_h;
//   // iEvent.getByLabel(electronsInputTag_, els_coll_h);  

//   unsigned int elsIndex = 0;
//   GsfElectronCollection::const_iterator els_end = elsHandle_->end();  // Iterator
//   for( GsfElectronCollection::const_iterator el = elsHandle_->begin(); el != els_end; ++el, ++elsIndex ) {
//     LorentzVector lv(el->p4());

//     // if (verbose_) cout << " - reco ele: pt: " << lv.pt() << ", eta: " << lv.eta() << ", phi: " << lv.phi(); 

//     if (lv.pt() < sublPtThresh_e) {
//       // if (verbose_) cout << ", FAILS pt" << std::endl;
//       continue;
//     }

//     // check for match to trigger objects
//     bool match = false;
//     if ( (hlt_el_lead_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(lv,hlt_el_lead.lv) < dr_trigmatch) ) match = true;
//     else if ( (hlt_el_subl_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(lv,hlt_el_subl.lv) < dr_trigmatch) ) match = true;
//     if (!match && reqTrigMatch_) {
//       // if (verbose_) cout << ", FAILS trig match" << std::endl;
//       continue;
//     }

//     // check loose ID (probably 2012 version)
// //     float chpfiso = el->pfIsolationVariables().sumChargedHadronPt;
// //     float nhpfiso = el->pfIsolationVariables().sumNeutralHadronPt;
// //     float phpfiso = el->pfIsolationVariables().sumPhotonEt;
// // bool pass_loose = EgammaCutBasedEleId::TestWP(EgammaCutBasedEleId::LOOSE , el, convs_h, beamSpotreco, vertexHandle, pfiso_ch, pfiso_em, pfiso_nh, rhoIso);

//     // if (verbose_) cout << ", passes presel cuts" << std::endl;

//     // pt ordering
//     if ( ((off_el_lead_idx == -1) || (lv.pt() > off_el_lead.lv.pt())) && (lv.pt() > leadPtThresh_e) ) {
//       off_el_subl_idx = off_el_lead_idx;
//       off_el_subl.lv = off_el_lead.lv;
//       off_el_lead_idx = (int)elsIndex;
//       off_el_lead.lv = lv;
//     } else if ( ((off_el_subl_idx == -1) || (lv.pt() > off_el_subl.lv.pt())) && (lv.pt() > sublPtThresh_e)) {
//       // check dR to remove exact duplicates
//       if (ROOT::Math::VectorUtil::DeltaR(off_el_lead.lv,lv) > 0.001) {
// 	off_el_subl_idx = (int)elsIndex;
// 	off_el_subl.lv = lv;
//       }
//     }

//   } // loop on reco gsf electrons

//   if (verbose_ && (isem || isme) && (off_el_lead_idx < 0)) {
//     cout << "TrigAnalyzerRECORef::analyzeTrigger: no good reco electrons: N(gsfele) =  " << elsHandle_->size() << std::endl;
//   }


//   // printout for cases where trigger is inefficient
//   //  if ( verbose_ && ismm && (off_mu_trigiso_lead_idx >= 0) && (off_mu_trigiso_subl_idx >= 0) ) {
//   if ( verbose_ && (triggerEnum == mm) && (off_mu_trigiso_lead_idx >= 0) && (off_mu_trigiso_subl_idx >= 0) ) {
//     // if ( verbose_ && (triggerEnum == em) && (off_mu_tliso_lead_idx >= 0) && (off_el_lead_idx >= 0) ) {
//     // check that this is the noniso trigger, and that the iso trigger failed
//     LorentzVector dilep = off_mu_trigiso_lead.lv + off_mu_trigiso_subl.lv;
//     float deltaR = ROOT::Math::VectorUtil::DeltaR(off_mu_trigiso_lead.lv,off_mu_trigiso_subl.lv);
//     cout << "-- event passes noniso trigger: " << triggerName;
//     cout << ", run: " << iEvent.id().run()
// 	 << ", event: " << iEvent.id().event() << ", mass: " << dilep.M() << ", deltaR: " << deltaR << endl
// 	 << "     hlt mu1 pt: " << hlt_mu_lead.lv.pt() << ", iso: " << hlt_mu_lead.trkiso
// 	 << ", hlt mu2 pt: " << hlt_mu_subl.lv.pt() << ", iso: " << hlt_mu_subl.trkiso << endl
// 	 << "     off mu1 pt: " << off_mu_trigiso_lead.lv.pt() << ", iso: " << off_mu_trigiso_lead.trkiso
// 	 << ", off mu2 pt: " << off_mu_trigiso_subl.lv.pt() << ", iso: " << off_mu_trigiso_subl.trkiso << endl;
//     cout << "     trig enum: " << triggerEnum << ", iso trigger enum: " << isoTriggerEnum
// 	 << ", iso trig decision: " << (1 << isoTriggerEnum) << endl;
//     // LorentzVector dilep = off_mu_tliso_lead.lv + off_el_lead.lv;
//     // float deltaR = ROOT::Math::VectorUtil::DeltaR(off_mu_tliso_lead.lv,off_el_lead.lv);
//     bool ineff = false;
//     if ( (isoTriggerEnum != notrig) && ((trigpass_results_ & (1 << isoTriggerEnum)) == 0) ) {
//       // if ( (isoTriggerEnum != notrig) && ((trigpass_results_ & (1 << (isoTriggerEnum - (unsigned int)emi) )) == 0) ) {
//       cout << "-- HLT inefficiency, event level! Non-iso trigger: " << triggerName;
//       ineff = true;
//     }
//     // else if ( (isoTriggerEnum != notrig) && ((trigpass_results_offdilep_ & (1 << isoTriggerEnum)) == 0) ) {
//     //   cout << "++ HLT inefficiency, offline object level! Non-iso trigger: " << triggerName;
//     //   ineff = true;
//     // }
//     if (ineff) {
//       cout << ", run: " << iEvent.id().run()
//       	   << ", event: " << iEvent.id().event() << ", mass: " << dilep.M() << ", deltaR: " << deltaR << endl
//       	   << "     hlt mu1 pt: " << hlt_mu_lead.lv.pt() << ", iso: " << hlt_mu_lead.trkiso
//       	   << ", hlt mu2 pt: " << hlt_mu_subl.lv.pt() << ", iso: " << hlt_mu_subl.trkiso << endl
//       	   << "     off mu1 pt: " << off_mu_trigiso_lead.lv.pt() << ", iso: " << off_mu_trigiso_lead.trkiso
//       	   << ", off mu2 pt: " << off_mu_trigiso_subl.lv.pt() << ", iso: " << off_mu_trigiso_subl.trkiso << endl;
//       // cout << ", run: " << iEvent.id().run()
//       // 	   << ", event: " << iEvent.id().event() << ", mass: " << dilep.M() << ", deltaR: " << deltaR << endl
//       // 	   << "     hlt mu pt: " << hlt_mu_lead.lv.pt() << ", iso: " << hlt_mu_lead.trkiso
//       // 	   << ", hlt el pt: " << hlt_el_lead.lv.pt()  << endl
//       // 	   << "     off mu pt: " << off_mu_tliso_lead.lv.pt() << ", iso: " << off_mu_tliso_lead.trkiso
//       // 	   << ", off el pt: " << off_el_lead.lv.pt() << endl;
//       // loop over pf cands near the hlt muons and print them out..
//       if (dumpHLTPFCands_) {
// 	// edm::Handle<reco::PFCandidateCollection> hltPFCandsHandle = hltPFCandsGlbHandle_;
// 	edm::Handle<reco::TrackCollection> hltTracksHandle = hltTracksGlbHandle_;
// 	// if (triggerEnum == mm) hltPFCandsHandle = hltPFCandsGlbHandle_;
// 	// else hltPFCandsHandle = hltPFCandsTrkHandle_;
// 	const VertexCollection* hltVertexCollection = hltVertexHandle_.product();

// 	// find hlt tracks near hlt muons
//         TrackCollection::const_iterator hlt_trks_end = hltTracksHandle->end();  // Iterator
// 	std::vector<reco::TrackRef> hlt_trks_lead;
// 	std::vector<reco::TrackRef> hlt_trks_subl;
// 	//	for ( TrackCollection::const_iterator hlt_trk = hltTracksHandle->begin(); hlt_trk != hlt_trks_end; ++hlt_trk ) {
// 	for ( unsigned int hlt_trk_idx = 0; hlt_trk_idx < hltTracksHandle->size(); ++hlt_trk_idx ) {
// 	  reco::TrackRef hlt_trk_ref(hltTracksHandle,hlt_trk_idx);
// 	  if (ROOT::Math::VectorUtil::DeltaR(hlt_mu_lead.lv,hlt_trk_ref->momentum()) < 0.3) {
// 	    hlt_trks_lead.push_back(hlt_trk_ref);
// 	  }
// 	  if (ROOT::Math::VectorUtil::DeltaR(hlt_mu_subl.lv,hlt_trk_ref->momentum()) < 0.3) {
// 	    hlt_trks_subl.push_back(hlt_trk_ref);
// 	  }
// 	} // loop on hlt tracks

// 	// find hlt pf cands near hlt muons
//         // PFCandidateCollection::const_iterator hlt_pfcands_end = hltPFCandsHandle->end();  // Iterator
// 	// PFCandidateCollection hlt_pfcands_lead;
// 	// PFCandidateCollection hlt_pfcands_subl;
// 	// for ( PFCandidateCollection::const_iterator hlt_pfcand = hltPFCandsHandle->begin(); hlt_pfcand != hlt_pfcands_end; ++hlt_pfcand ) {
// 	//   if (ROOT::Math::VectorUtil::DeltaR(hlt_mu_lead.lv,hlt_pfcand->p4()) < 0.3) {
// 	//     hlt_pfcands_lead.push_back(*hlt_pfcand);
// 	//   }
// 	//   // if (ROOT::Math::VectorUtil::DeltaR(hlt_mu_subl.lv,hlt_pfcand->p4()) < 0.3) {
// 	//   //   hlt_pfcands_subl.push_back(*hlt_pfcand);
// 	//   // }
// 	// } // loop on hlt pf cands

// 	// find off pf cands near hlt muons
//         PFCandidateCollection::const_iterator off_pfcands_end = offPFCandsHandle_->end();  // Iterator
// 	PFCandidateCollection off_pfcands_lead;
// 	PFCandidateCollection off_pfcands_subl;
// 	for ( PFCandidateCollection::const_iterator off_pfcand = offPFCandsHandle_->begin(); off_pfcand != off_pfcands_end; ++off_pfcand ) {
// 	  if (abs(off_pfcand->charge()) != 1) continue;

// 	  if (ROOT::Math::VectorUtil::DeltaR(hlt_mu_lead.lv,off_pfcand->p4()) < 0.3) {
// 	    off_pfcands_lead.push_back(*off_pfcand);
// 	  }
// 	  if (ROOT::Math::VectorUtil::DeltaR(hlt_mu_subl.lv,off_pfcand->p4()) < 0.3) {
// 	    off_pfcands_subl.push_back(*off_pfcand);
// 	  }
// 	} // loop on off pf cands

// 	// dump hlt pf cands near lead hlt muon
// 	cout << "   ---- hlt tracks near leading hlt muon:" << endl;
//         std::vector<reco::TrackRef>::const_iterator hlt_trks_lead_end = hlt_trks_lead.end();  // Iterator
// 	for ( std::vector<reco::TrackRef>::const_iterator hlt_trk = hlt_trks_lead.begin(); hlt_trk != hlt_trks_lead_end; ++hlt_trk ) {
// 	  float dR = ROOT::Math::VectorUtil::DeltaR(hlt_mu_lead.lv,(**hlt_trk).momentum());
// 	  int vtx = trackVertex(*hltVertexCollection,*hlt_trk);
// 	  float dz = hlt_mu_lead.vz - (**hlt_trk).vz();
// 	  //	  if (fabs(dz) > 1.0) continue;
// 	  if ((fabs(dz) < 0.2) && (dR > 0.01)) cout << "    * ";
// 	  else cout << "      ";
// 	  cout << "pt: " << (**hlt_trk).pt() << ", eta: " << (**hlt_trk).eta() << ", phi: " << (**hlt_trk).phi()
// 	       << ", vz: " << (**hlt_trk).vz() << ", vtx: " << vtx
// 	       << ", nhits: " << (**hlt_trk).numberOfValidHits() << ", algo: " << (**hlt_trk).algo() 
// 	       << ", dR: " << dR << endl;
// 	}

// 	// dump hlt pf cands near lead hlt muon
// 	// cout << "   ++++ hlt pf cands near leading hlt muon:" << endl;
//         // PFCandidateCollection::const_iterator hlt_pfcands_lead_end = hlt_pfcands_lead.end();  // Iterator
// 	// for ( PFCandidateCollection::const_iterator hlt_pfcand = hlt_pfcands_lead.begin(); hlt_pfcand != hlt_pfcands_lead_end; ++hlt_pfcand ) {
// 	//   float dR = ROOT::Math::VectorUtil::DeltaR(hlt_mu_lead.lv,hlt_pfcand->p4());
// 	//   int vtx = chargedHadronVertex(*hltVertexCollection,*hlt_pfcand);
// 	//   float dz = hlt_mu_lead.vz - hlt_pfcand->vz();
// 	//   //	  if (fabs(dz) > 1.0) continue;
// 	//   if (fabs(dz) < 0.2) cout << "    * ";
// 	//   else cout << "      ";
// 	//   cout << "pt: " << hlt_pfcand->pt() << ", eta: " << hlt_pfcand->eta() << ", phi: " << hlt_pfcand->phi()
// 	//        << ", id: " << hlt_pfcand->particleId() << ", vz: " << hlt_pfcand->vz() << ", vtx: " << vtx
// 	//        << ", nhits: " << hlt_pfcand->trackRef()->numberOfValidHits() << ", algo: " << hlt_pfcand->trackRef()->algo() 
// 	//        << ", dR: " << dR << endl;
// 	// }

// 	// dump off pf cands near lead hlt muon
// 	cout << "   ++++ off pf cands near leading hlt muon:" << endl;
//         PFCandidateCollection::const_iterator off_pfcands_lead_end = off_pfcands_lead.end();  // Iterator
// 	for ( PFCandidateCollection::const_iterator off_pfcand = off_pfcands_lead.begin(); off_pfcand != off_pfcands_lead_end; ++off_pfcand ) {
// 	  float dR = ROOT::Math::VectorUtil::DeltaR(hlt_mu_lead.lv,off_pfcand->p4());
// 	  int vtx = chargedHadronVertex(*vertexCollection,*off_pfcand);
// 	  float dz = hlt_mu_lead.vz - off_pfcand->vz();
// 	  if (fabs(dz) > 1.0) continue;
// 	  if (((vtx == 0) || (vtx == -1)) && (off_pfcand->particleId() == 1)) cout << "    * ";
// 	  else cout << "      ";
// 	  cout << "pt: " << off_pfcand->pt() << ", eta: " << off_pfcand->eta() << ", phi: " << off_pfcand->phi()
// 	       << ", id: " << off_pfcand->particleId() << ", vz: " << off_pfcand->vz() << ", vtx: " << vtx
// 	       << ", dR: " << dR << endl;
// 	  //	       << ", nhits: " << off_pfcand->trackRef()->numberOfValidHits() << ", dR: " << dR << endl;
// 	}

// 	// dump hlt tracks near subl hlt muon
// 	cout << "   ---- hlt tracks near subleading hlt muon:" << endl;
// 	std::vector<reco::TrackRef>::const_iterator hlt_trks_subl_end = hlt_trks_subl.end();  // Iterator
// 	for ( std::vector<reco::TrackRef>::const_iterator hlt_trk = hlt_trks_subl.begin(); hlt_trk != hlt_trks_subl_end; ++hlt_trk ) {
// 	  float dR = ROOT::Math::VectorUtil::DeltaR(hlt_mu_subl.lv,(**hlt_trk).momentum());
// 	  int vtx = trackVertex(*hltVertexCollection,*hlt_trk);
// 	  float dz = hlt_mu_subl.vz - (**hlt_trk).vz();
// 	  //	  if (fabs(dz) > 1.0) continue;
// 	  if ((fabs(dz) < 0.2) && (dR > 0.01)) cout << "    * ";
// 	  else cout << "      ";
// 	  cout << "pt: " << (**hlt_trk).pt() << ", eta: " << (**hlt_trk).eta() << ", phi: " << (**hlt_trk).phi()
// 	       << ", vz: " << (**hlt_trk).vz() << ", vtx: " << vtx
// 	       << ", nhits: " << (**hlt_trk).numberOfValidHits() << ", algo: " << (**hlt_trk).algo() 
// 	       << ", dR: " << dR << endl;
// 	}

// 	//   // dump hlt pf cands near subl hlt muon
// 	//   cout << "   ++++ hlt pf cands near subleading hlt muon:" << endl;
// 	//   PFCandidateCollection::const_iterator hlt_pfcands_subl_end = hlt_pfcands_subl.end();  // Iterator
// 	//   for ( PFCandidateCollection::const_iterator hlt_pfcand = hlt_pfcands_subl.begin(); hlt_pfcand != hlt_pfcands_subl_end; ++hlt_pfcand ) {
// 	//     float dR = ROOT::Math::VectorUtil::DeltaR(hlt_mu_subl.lv,hlt_pfcand->p4());
// 	//     int vtx = chargedHadronVertex(*hltVertexCollection,*hlt_pfcand);
// 	//     float dz = hlt_mu_lead.vz - hlt_pfcand->vz();
// 	//     //if (fabs(dz) > 1.0) continue;
// 	//     if (fabs(dz) < 0.2) cout << "    * ";
// 	//     else cout << "      ";
// 	//     cout << "pt: " << hlt_pfcand->pt() << ", eta: " << hlt_pfcand->eta() << ", phi: " << hlt_pfcand->phi()
// 	// 	 << ", id: " << hlt_pfcand->particleId() << ", vz: " << hlt_pfcand->vz() << ", vtx: " << vtx
// 	// 	 << ", nhits: " << hlt_pfcand->trackRef()->numberOfValidHits() << ", algo: " << hlt_pfcand->trackRef()->algo() 
// 	// 	 << ", dR: " << dR << endl;

// 	// dump off pf cands near subl hlt muon
// 	cout << "   ++++ off pf cands near subleading hlt muon:" << endl;
// 	PFCandidateCollection::const_iterator off_pfcands_subl_end = off_pfcands_subl.end();  // Iterator
// 	for ( PFCandidateCollection::const_iterator off_pfcand = off_pfcands_subl.begin(); off_pfcand != off_pfcands_subl_end; ++off_pfcand ) {
// 	  float dR = ROOT::Math::VectorUtil::DeltaR(hlt_mu_subl.lv,off_pfcand->p4());
// 	  int vtx = chargedHadronVertex(*vertexCollection,*off_pfcand);
// 	  float dz = hlt_mu_lead.vz - off_pfcand->vz();
// 	  if (fabs(dz) > 1.0) continue;
// 	  if (((vtx == 0) || (vtx == -1)) && (off_pfcand->particleId() == 1)) cout << "    * ";
// 	  else cout << "      ";
// 	  cout << "pt: " << off_pfcand->pt() << ", eta: " << off_pfcand->eta() << ", phi: " << off_pfcand->phi()
// 	       << ", id: " << off_pfcand->particleId() << ", vz: " << off_pfcand->vz() << ", vtx: " << vtx
// 	       << ", nhits: " << off_pfcand->trackRef()->numberOfValidHits() << ", dR: " << dR << endl;
// 	}


//       } // if dumpHLTPFCands
//     } // if ineff

//   }

//   return true;
// }

//____________________________________________________________________________
void TrigAnalyzerRECORef::bookHists(edm::Service<TFileService>& fs, const std::string& suffix) {

  //  std::cout << "TrigAnalyzerRECORef::bookHists: called with suffix: " << suffix << std::endl;

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  hists_1d_["h_pt"+suf] = fs->make<TH1F>(Form("h_pt%s",suf.c_str()) , ";  p_{T} [GeV]" , 100 , 0. , 100. );
  hists_1d_["h_eta"+suf] = fs->make<TH1F>(Form("h_eta%s",suf.c_str()) , ";  #eta" , 100 , -3. , 3. );

  hists_1d_["h_abstrkiso"+suf] = fs->make<TH1F>(Form("h_abstrkiso%s",suf.c_str()) , ";  charged pfiso [GeV]" , 200 , 0. , 10. );
  hists_1d_["h_reltrkiso"+suf] = fs->make<TH1F>(Form("h_reltrkiso%s",suf.c_str()) , ";  charged pfiso / p_{T}" , 200 , 0. , 2. );
  hists_1d_["h_abspfiso"+suf] = fs->make<TH1F>(Form("h_abspfiso%s",suf.c_str()) , ";  pfiso [GeV]" , 200 , 0. , 10. );
  hists_1d_["h_relpfiso"+suf] = fs->make<TH1F>(Form("h_relpfiso%s",suf.c_str()) , ";  pfiso / p_{T}" , 200 , 0. , 2. );

  return;
}

//____________________________________________________________________________
void TrigAnalyzerRECORef::fillHists(const reco::Muon& muon, const std::string& suffix) {

  //  std::cout << "TrigAnalyzerRECORef::fillHists: called with suffix: " << suffix << std::endl;

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  hists_1d_["h_pt"+suf]->Fill(muon.pt());
  hists_1d_["h_eta"+suf]->Fill(muon.eta());
  float trkiso = muon.pfIsolationR03().sumChargedHadronPt;
  float pfiso = muonPFiso(muon);
  hists_1d_["h_abstrkiso"+suf]->Fill(trkiso);
  hists_1d_["h_reltrkiso"+suf]->Fill(trkiso/muon.pt());
  hists_1d_["h_abspfiso"+suf]->Fill(pfiso);
  hists_1d_["h_relpfiso"+suf]->Fill(pfiso/muon.pt());

  return;
}

//____________________________________________________________________________
float TrigAnalyzerRECORef::muonPFiso(const reco::Muon& muon) {

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
int TrigAnalyzerRECORef::chargedHadronVertex( const reco::VertexCollection& vertices, const reco::PFCandidate& pfcand ) const {
  if (pfcand.trackRef().isNull()) return -2;
  return trackVertex(vertices,pfcand.trackRef());
}

//____________________________________________________________________________
// returns best match vertex.  Taken from PFPileUpAlgo, chargedHadronVertex
int TrigAnalyzerRECORef::trackVertex( const reco::VertexCollection& vertices, const reco::TrackRef& trackRef ) const {

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


