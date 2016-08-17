/** \class SingleMuTrigTnPAnalyzerRECO
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
#include "TrigAnalyzer/DilepTrigAnalyzer/interface/SingleMuTrigTnPAnalyzerRECO.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

// ROOT includes
#include "Math/VectorUtil.h"

#include <cassert>

using namespace reco;
using namespace edm;

//
// constructors and destructor
//
//____________________________________________________________________________
SingleMuTrigTnPAnalyzerRECO::SingleMuTrigTnPAnalyzerRECO(const edm::ParameterSet& ps) 
{
  using namespace std;
  using namespace edm;

  processName_ = ps.getUntrackedParameter<std::string>("processName","HLT");
  refTriggerName_ = ps.getUntrackedParameter<std::string>("refTriggerName","HLT_IsoMu24_v3");
  sigTriggerName_ = ps.getUntrackedParameter<std::string>("sigTriggerName","HLT_TkMu24_v1");
  denomFilterName_ = ps.getUntrackedParameter<std::string>("denomFilterName","hltL1fL1sMu22L1Filtered0");
  triggerResultsToken_ = consumes<edm::TriggerResults> (ps.getUntrackedParameter<edm::InputTag>("triggerResultsTag", edm::InputTag("TriggerResults", "", "reHLT")));
  triggerEventToken_ = consumes<trigger::TriggerEvent> (ps.getUntrackedParameter<edm::InputTag>("triggerEventTag", edm::InputTag("hltTriggerSummaryAOD")));
  muonsToken_ = consumes<View<reco::Muon> > (ps.getUntrackedParameter<edm::InputTag>("muonsInputTag",edm::InputTag("muons")));
  vtxToken_ = consumes<reco::VertexCollection> (ps.getUntrackedParameter<edm::InputTag>("vtxInputTag",edm::InputTag("offlinePrimaryVertices")));
  tagPt_ = ps.getUntrackedParameter<double>("tagPt",25.);
  tagEta_ = ps.getUntrackedParameter<double>("tagEta",2.4);
  probePt_ = ps.getUntrackedParameter<double>("probePt",20.);
  probeEta_ = ps.getUntrackedParameter<double>("probeEta",2.4);
  probePtForEta_ = ps.getUntrackedParameter<double>("probePt",30.);
  verbose_ = ps.getUntrackedParameter<bool>("verbose",false);
    
  // histogram setup
  edm::Service<TFileService> fs;
  hists_1d_["h_passtrig"] = fs->make<TH1F>("h_passtrig" , "; passed trigger" , 2 , 0. , 2. );
  hists_1d_["h_mll_allpairs"] = fs->make<TH1F>("h_mll_allpairs" , "; m_{ll} [GeV]" , 75 , 0. , 150. );
  hists_1d_["h_mll_cut"] = fs->make<TH1F>("h_mll_cut" , "; m_{ll} [GeV]" , 75 , 0. , 150. );
  bookHists(fs,"probe_all");
  bookHists(fs,"probe_pass");
  bookHists(fs,"probe_fail");

}

//____________________________________________________________________________
SingleMuTrigTnPAnalyzerRECO::~SingleMuTrigTnPAnalyzerRECO()
{
}

//
// member functions
//
//____________________________________________________________________________
void
SingleMuTrigTnPAnalyzerRECO::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  using namespace std;
  using namespace edm;

  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
      const unsigned int n(hltConfig_.size());
      // check if trigger names in (new) config
      unsigned int refTriggerIndex(hltConfig_.triggerIndex(refTriggerName_));
      if (refTriggerIndex>=n) {
	cout << "SingleMuTrigAnalyzerMiniAOD::analyze:"
	     << " RefTriggerName " << refTriggerName_ 
	     << " not available in config!" << endl;
      }
      unsigned int sigTriggerIndex(hltConfig_.triggerIndex(sigTriggerName_));
      if (sigTriggerIndex>=n) {
	cout << "SingleMuTrigAnalyzerMiniAOD::analyze:"
	     << " SigTriggerName " << sigTriggerName_ 
	     << " not available in config!" << endl;
      }
    } // if changed
  } else {
    cout << "DiMuTrigAnalyzerRECO::analyze:"
	 << " config extraction failure with process name "
	 << processName_ << endl;
  }

}

//____________________________________________________________________________
// ------------ method called to produce the data  ------------
void
SingleMuTrigTnPAnalyzerRECO::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  // get event products
  iEvent.getByToken(triggerResultsToken_,triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
    cout << "SingleMuTrigAnalyzerMiniAOD::analyze: Error in getting TriggerResults product from Event!" << endl;
    return;
  }
  iEvent.getByToken(triggerEventToken_,triggerEventHandle_);
  if (!triggerEventHandle_.isValid()) {
    cout << "Error in getting TriggerEvent product from Event!" << endl;
    return;
  }

  // sanity check
  assert(triggerResultsHandle_->size()==hltConfig_.size());

  // retrieve necessary containers
  Handle<reco::VertexCollection> vertexHandle_;
  iEvent.getByToken(vtxToken_, vertexHandle_);
  Handle<View<reco::Muon> > musHandle_;
  iEvent.getByToken( muonsToken_ , musHandle_ );

  if (verbose_) cout << endl;

  //-------------------------------------
  //   check trigger results
  //-------------------------------------

  const unsigned int ntrigs(hltConfig_.size());
  const unsigned int refTriggerIndex(hltConfig_.triggerIndex(refTriggerName_));
  assert(refTriggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(refTriggerName_));
  const unsigned int sigTriggerIndex(hltConfig_.triggerIndex(sigTriggerName_));
  assert(sigTriggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(sigTriggerName_));

  // abort on invalid trigger name
  if (refTriggerIndex>=ntrigs) {
    cout << "HTTrigAnalyzerMiniAOD::analyzeTrigger: path "
	 << refTriggerName_ << " - not found!" << endl;
    return;
  }
  if (sigTriggerIndex>=ntrigs) {
    cout << "HTTrigAnalyzerMiniAOD::analyzeTrigger: path "
	 << sigTriggerName_ << " - not found!" << endl;
    return;
  }

  if (verbose_) {
    cout << "HTTrigAnalyzerMiniAOD::analyzeTrigger: reference path "
	 << refTriggerName_ << " [" << refTriggerIndex << "]" << endl;
    cout << "HTTrigAnalyzerMiniAOD::analyzeTrigger: signal path "
	 << sigTriggerName_ << " [" << sigTriggerIndex << "]" << endl;
  }
  
  // modules on this trigger path
  bool refAccept = triggerResultsHandle_->accept(refTriggerIndex);
  bool sigAccept = triggerResultsHandle_->accept(sigTriggerIndex);

  if (refAccept) hists_1d_["h_passtrig"]->Fill(1);
  else {  
    // don't consider event if ref trigger didn't fire
    hists_1d_["h_passtrig"]->Fill(0);
    return;
  }

  //------------------------------------
  //  hlt objects
  //------------------------------------

  std::vector<LorentzVector> refTriggerObjects;
  if (refAccept) {
    refTriggerObjects = getTriggerObjects(refTriggerIndex);
    if (refTriggerObjects.size() == 0) {
      cout << "DiMuTrigAnalyzerRECO::analyze: WARNING!! no valid trigger objects for ref path: " << refTriggerName_ << endl;
    }
  }

  std::vector<LorentzVector> denomTriggerObjects;
  if (denomFilterName_.size()) denomTriggerObjects = getTriggerObjects(sigTriggerIndex,denomFilterName_);
  
  std::vector<LorentzVector> sigTriggerObjects;
  if (sigAccept) {
    sigTriggerObjects = getTriggerObjects(sigTriggerIndex);
    if (sigTriggerObjects.size() == 0) {
      cout << "DiMuTrigAnalyzerRECO::analyze: WARNING!! no valid trigger objects for sig path: " << sigTriggerName_ << endl;
    }
  }
  
  //-------------------------------------
  //   reco vertices -- need for muon ID
  //-------------------------------------

  unsigned int nvtx = 0;
  
  // find vertex 0 in vertex container
  const VertexCollection* vertexCollection = vertexHandle_.product();
  VertexCollection::const_iterator firstGoodVertex = vertexCollection->end();
  for ( VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); ++vtx ) {
    if (  !vtx->isFake() && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0 ) {
      if (firstGoodVertex == vertexCollection->end()) firstGoodVertex = vtx;
      ++nvtx;
    }
  } // loop on vertices

  if (firstGoodVertex == vertexCollection->end()) {
    cout << "SingleMuTrigAnalyzerMiniAOD::analyze: didn't find any good offline vertices!! size: " 
	 << vertexCollection->size() << std::endl;
    return;
  }


  //-------------------------------------
  //   reco muons 
  //-------------------------------------

  const float dr_trigmatch = 0.1;
  const float dr_trigdenommatch = 0.3;
  
  std::vector<LorentzVector> offTagMuons;
  std::vector<LorentzVector> offProbeMuons;
  std::vector<LorentzVector> offProbeMatchedMuons;

  if (verbose_) cout << "found offline muons:" << endl;
  View<reco::Muon>::const_iterator muons_end = musHandle_->end();  // Iterator
  for ( View<reco::Muon>::const_iterator muon_tag = musHandle_->begin(); muon_tag != muons_end; ++muon_tag ) {

    if (verbose_) cout << " - tag cand: pt: " << muon_tag->pt() << ", eta: " << muon_tag->eta()
		       << ", phi: " << muon_tag->phi() << ", iso: " << muonPFiso(*muon_tag) << endl;
    
    // check if muon passes tag cuts
    if (muon_tag->pt() < tagPt_) continue;
    if (fabs(muon_tag->eta()) > tagEta_) continue;
    if (!muon::isTightMuon(*muon_tag,*firstGoodVertex)) continue;
    if (verbose_) cout << "   - passes tight" << endl;
    if (muonPFiso(*muon_tag) > 0.15) continue;
    if (verbose_) cout << "   - passes iso" << endl;
    
    // check if muon matches trigger
    bool trigmatch_tag = false;
    for (unsigned int itrig=0; itrig < refTriggerObjects.size(); ++itrig) {
      if (ROOT::Math::VectorUtil::DeltaR(muon_tag->p4(),refTriggerObjects.at(itrig)) < dr_trigmatch) trigmatch_tag = true;
    }
    if (!trigmatch_tag) continue;
    if (verbose_) cout << "   - matched to trigger" << endl;

    // good tag muon: look for any probe muons in the event
    for ( View<reco::Muon>::const_iterator muon_probe = musHandle_->begin(); muon_probe != muons_end; ++muon_probe ) {
      
      // make sure probe isn't the same as the tag
      if (ROOT::Math::VectorUtil::DeltaR(muon_tag->p4(),muon_probe->p4()) < 0.01) continue;

      // check probe cuts
      if (muon_probe->pt() < probePt_) continue;
      if (fabs(muon_probe->eta()) > probeEta_) continue;
      if (!muon::isTightMuon(*muon_probe,*firstGoodVertex)) continue;
      if (muonPFiso(*muon_probe) > 0.15) continue;

      if (denomFilterName_.size()) {
	// check if probe muon matches trigger denom
	bool trigdenommatch_probe = false;
	for (unsigned int itrig=0; itrig < denomTriggerObjects.size(); ++itrig) {
	  if (ROOT::Math::VectorUtil::DeltaR(muon_probe->p4(),denomTriggerObjects.at(itrig)) < dr_trigdenommatch) trigdenommatch_probe = true;
	}
	if (!trigdenommatch_probe) continue;
      }
      
      if (verbose_) cout << " - probe cand: pt: " << muon_probe->pt() << ", eta: " << muon_probe->eta()
			 << ", phi: " << muon_probe->phi() << ", iso: " << muonPFiso(*muon_probe) << endl;
    
      // check dimuon mass: must be in range 81-101
      LorentzVector dimuon = LorentzVector(muon_tag->p4() + muon_probe->p4());
      hists_1d_["h_mll_allpairs"]->Fill(dimuon.M());
      if (verbose_) cout << " - probe dimuon mass: " << dimuon.M()  << endl;
      if (dimuon.M() < 81. || dimuon.M() > 101.) continue;
      hists_1d_["h_mll_cut"]->Fill(dimuon.M());

      fillHists(LorentzVector(muon_probe->p4()),"probe_all");

      // check if probe muon matches trigger
      bool trigmatch_probe = false;
      for (unsigned int itrig=0; itrig < sigTriggerObjects.size(); ++itrig) {
	if (ROOT::Math::VectorUtil::DeltaR(muon_probe->p4(),sigTriggerObjects.at(itrig)) < dr_trigmatch) trigmatch_probe = true;
      }

      if (trigmatch_probe) fillHists(LorentzVector(muon_probe->p4()),"probe_pass");
      else fillHists(LorentzVector(muon_probe->p4()),"probe_fail");

    } // loop over probes
  } // loop over tags (offline muons)

  if (verbose_) cout << endl;
  return;
}

//____________________________________________________________________________
void SingleMuTrigTnPAnalyzerRECO::bookHists(edm::Service<TFileService>& fs, const std::string& suffix) {

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  hists_1d_["h_pt"+suf] = fs->make<TH1F>(Form("h_pt%s",suf.c_str()) , "; p_{T} [GeV]" , 100 , 0. , 100. );
  hists_1d_["h_eta"+suf] = fs->make<TH1F>(Form("h_eta%s",suf.c_str()) , "; #eta" , 100 , -3. , 3. );
  hists_1d_["h_phi"+suf] = fs->make<TH1F>(Form("h_phi%s",suf.c_str()) , "; #phi" , 100 , -3.14 , 3.14 );

  return;
}

//____________________________________________________________________________
void SingleMuTrigTnPAnalyzerRECO::fillHists(const LorentzVector& lv, const std::string& suffix) {

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  hists_1d_["h_pt"+suf]->Fill(lv.pt());
if (lv.pt() > probePtForEta_) {
    hists_1d_["h_eta"+suf]->Fill(lv.eta());
    hists_1d_["h_phi"+suf]->Fill(lv.phi());
  }

  return;
}

//____________________________________________________________________________
float SingleMuTrigTnPAnalyzerRECO::muonPFiso(const reco::Muon& muon) {

  reco::MuonPFIsolation pfStructR03 = muon.pfIsolationR03();
  float chiso = pfStructR03.sumChargedHadronPt;
  float nhiso = pfStructR03.sumNeutralHadronEt;
  float emiso = pfStructR03.sumPhotonEt;
  float deltaBeta = pfStructR03.sumPUPt;

  //  float absiso = chiso + nhiso + emiso;
  float absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return absiso/muon.pt();

}

//____________________________________________________________________________
std::vector<LorentzVector> SingleMuTrigTnPAnalyzerRECO::getTriggerObjects(const unsigned int& triggerIndex, const std::string& targetModuleLabel) {

  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;
  
  std::vector<LorentzVector> triggerObjects;
  
  // modules on this trigger path
  const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));
  const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
  
  for (unsigned int j=moduleIndex; j!=0; --j) {
    const string& moduleLabel(moduleLabels[j]);
    const string  moduleType(hltConfig_.moduleType(moduleLabel));
    // check whether the module is packed up in TriggerEvent product
    const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(moduleLabel,"",processName_)));
    if (filterIndex>=triggerEventHandle_->sizeFilters()) continue;
    if (targetModuleLabel.size() > 0 && (targetModuleLabel != moduleLabel)) continue;
    //if (filterIndex>=triggerEventHandle_->size()) continue;
    if (verbose_) {
      cout << " 'L3' filter in slot " << j << " - label/type " << moduleLabel << "/" << moduleType << endl
	   << " Filter packed up at: " << filterIndex << endl;
    }
    if (moduleLabel == "hltBoolEnd") continue;
    const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
    const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
    const size_type nI(VIDS.size());
    const size_type nK(KEYS.size());
    assert(nI==nK);
    const size_type n(max(nI,nK));
    if (verbose_) cout << "   " << n  << " accepted 'L3' objects found: " << endl;
    const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
    for (size_type i=0; i!=n; ++i) {
      const TriggerObject& TO(TOC[KEYS[i]]);
      LorentzVector lv( TO.particle().p4() );

      triggerObjects.push_back(lv);

      if (verbose_) {
    	cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
    	     << TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()
    	     << endl;
      }
    } // loop on trig objects

    if (triggerObjects.size() > 0) break;
  } // backwards loop on modules

  return triggerObjects;
}

