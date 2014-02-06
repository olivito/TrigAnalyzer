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
DilepTrigAnalyzerRECO::DilepTrigAnalyzerRECO(const edm::ParameterSet& ps) : 
  processName_(ps.getParameter<std::string>("processName")),
  triggerName_(ps.getParameter<std::string>("triggerName")),
  mmIsoTriggerName_(ps.getParameter<std::string>("mmIsoTriggerName")),
  mmtkIsoTriggerName_(ps.getParameter<std::string>("mmtkIsoTriggerName")),
  triggerResultsTag_(ps.getParameter<edm::InputTag>("triggerResults")),
  triggerEventTag_(ps.getParameter<edm::InputTag>("triggerEvent")),
  triggerEventWithRefsTag_(ps.getParameter<edm::InputTag>("triggerEventWithRefs")),
  electronsInputTag_(ps.getParameter<edm::InputTag>("electronsInputTag")),
  muonsInputTag_(ps.getParameter<edm::InputTag>("muonsInputTag")),
  vtxInputTag_(ps.getParameter<edm::InputTag>("vtxInputTag")),
  hltVtxInputTag_(ps.getParameter<edm::InputTag>("hltVtxInputTag")),
  getHLTIsoVals_(ps.getParameter<bool>("getHLTIsoVals")),
  isoValMapGlbTag_(ps.getParameter<edm::InputTag>("isoValMapGlb")),
  isoValMapTrkTag_(ps.getParameter<edm::InputTag>("isoValMapTrk")),
  dumpHLTPFCands_(ps.getParameter<bool>("dumpHLTPFCands")),
  hltPFCandsGlbTag_(ps.getParameter<edm::InputTag>("hltPFCandsGlb")),
  hltPFCandsTrkTag_(ps.getParameter<edm::InputTag>("hltPFCandsTrk")),
  offPFCandsTag_(ps.getParameter<edm::InputTag>("offPFCands")),
  offLeadPt_(ps.getParameter<double>("offLeadPt")),
  offSublPt_(ps.getParameter<double>("offSublPt")),
  verbose_(ps.getParameter<bool>("verbose"))
{
  using namespace std;
  using namespace edm;

  cout << "DilepTrigAnalyzerRECO configuration: " << endl
       << "   ProcessName = " << processName_ << endl
       << "   TriggerName = " << triggerName_ << endl
       << "   mmIsoTriggerName = " << mmIsoTriggerName_ << endl
       << "   mmtkIsoTriggerName = " << mmtkIsoTriggerName_ << endl
       << "   TriggerResultsTag = " << triggerResultsTag_.encode() << endl
       << "   TriggerEventTag = " << triggerEventTag_.encode() << endl
       << "   TriggerEventWithRefsTag = " << triggerEventWithRefsTag_.encode() << endl
       << "   ElectronsInputTag = " << electronsInputTag_.encode() << endl
       << "   MuonsInputTag = " << muonsInputTag_.encode() << endl
       << "   HLTVtxInputTag = " << hltVtxInputTag_.encode() << endl
       << "   VtxInputTag = " << vtxInputTag_.encode() << endl
       << "   GetHLTIsoVals = " << getHLTIsoVals_ << endl
       << "   IsoValMapGlbTag = " << isoValMapGlbTag_.encode() << endl
       << "   IsoValMapTrkTag = " << isoValMapTrkTag_.encode() << endl
       << "   DumpHLTPFCands = " << dumpHLTPFCands_ << endl
       << "   HLTPFCandsGlbTag = " << hltPFCandsGlbTag_.encode() << endl
       << "   HLTPFCandsTrkTag = " << hltPFCandsTrkTag_.encode() << endl
       << "   OffPFCandsTrkTag = " << offPFCandsTag_.encode() << endl
       << "   OffLeadPt = " << offLeadPt_ << endl
       << "   OffSublPt = " << offSublPt_ << endl
       << "   Verbose = " << verbose_ << endl;

  //  hltTriggerNames_.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1");
  //  hltTriggerNames_.push_back("HLT_Mu17_IterTrkIsoVVL_Mu8_IterTrkIsoVVL_v1");
  //  hltTriggerNames_.push_back("HLT_Mu17_ChPFIsoVVL_Mu8_ChPFIsoVVL_v1");
  //  hltTriggerNames_.push_back("HLT_Mu17_NoMuChPFIsoVVL_Mu8_NoMuChPFIsoVVL_v1");
  hltTriggerNames_.push_back(mmIsoTriggerName_);
  hltTriggerNames_.push_back("HLT_Mu17_Mu8_v23");
  //  hltTriggerNames_.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1");
  //  hltTriggerNames_.push_back("HLT_Mu17_IterTrkIsoVVL_TkMu8_IterTrkIsoVVL_v1");
  //  hltTriggerNames_.push_back("HLT_Mu17_ChPFIsoVVL_TkMu8_ChPFIsoVVL_v1");
  //  hltTriggerNames_.push_back("HLT_Mu17_NoMuChPFIsoVVL_TkMu8_NoMuChPFIsoVVL_v1");
  hltTriggerNames_.push_back(mmtkIsoTriggerName_);
  hltTriggerNames_.push_back("HLT_Mu17_TkMu8_v16");
  hltTriggerNames_.push_back("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1");
  hltTriggerNames_.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10");
  hltTriggerNames_.push_back("HLT_Mu17_TrkIsoVVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1");
  hltTriggerNames_.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10");

  hltShortNames_.push_back("mmi");
  hltShortNames_.push_back("mm");
  hltShortNames_.push_back("mmitk");
  hltShortNames_.push_back("mmtk");
  hltShortNames_.push_back("emi");
  hltShortNames_.push_back("em");
  hltShortNames_.push_back("mei");
  hltShortNames_.push_back("me");

  // histogram setup
  edm::Service<TFileService> fs;
  if (triggerName_ == "dimu") {
    h_results_mm_ = fs->make<TH1F>("h_results_mm" , ";Trigger Results" , 16 , -0.5 , 15.5 );
    h_results_offdilep_mm_ = fs->make<TH1F>("h_results_offdilep_mm" , ";Trigger Results" , 16 , -0.5 , 15.5 );
    for (unsigned int itrig=0; itrig<=(unsigned int)mmtk; ++itrig) {
      bookHists(fs,hltShortNames_.at(itrig),true);
      bookHists(fs,hltShortNames_.at(itrig)+"_tight");
      bookHists(fs,hltShortNames_.at(itrig)+"_trigiso");
      bookHists(fs,hltShortNames_.at(itrig)+"_liso");
      bookHists(fs,hltShortNames_.at(itrig)+"_tliso");
      bookHists(fs,hltShortNames_.at(itrig)+"_tiso");
    }
  } else if (triggerName_ == "emu") {
    h_results_em_ = fs->make<TH1F>("h_results_em" , ";Trigger Results" , 16 , -0.5 , 15.5 );
    for (unsigned int itrig=(unsigned int)emi; itrig<=(unsigned int)me; ++itrig) {
      bookHists(fs,hltShortNames_.at(itrig),true);
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
	for (unsigned int itrig = 0; itrig <= (unsigned int)mmtk; ++itrig) {
	  unsigned int triggerIndex(hltConfig_.triggerIndex(hltTriggerNames_.at(itrig)));
	  if (triggerIndex>=n) {
	    cout << "DilepTrigAnalyzerRECO::analyze:"
		 << " TriggerName " << hltTriggerNames_.at(itrig) 
		 << " not available in (new) config!" << endl;
	  }
	} // loop over dimu triggers
      } // if dimu triggers 
      else if (triggerName_=="emu") {
	for (unsigned int itrig = (unsigned int)emi; itrig <= (unsigned int)me; ++itrig) {
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
  // iEvent.getByLabel(triggerEventTag_,triggerEventHandle_);
  // if (!triggerEventHandle_.isValid()) {
  //   cout << "DilepTrigAnalyzerRECO::analyze: Error in getting TriggerEvent product from Event!" << endl;
  //   return;
  // }
  iEvent.getByLabel(triggerEventWithRefsTag_,triggerEventWithRefsHandle_);
  if (!triggerEventWithRefsHandle_.isValid()) {
    cout << "DilepTrigAnalyzerRECO::analyze: Error in getting TriggerEventWithRefs product from Event!" << endl;
    return;
  }
  // sanity check
  assert(triggerResultsHandle_->size()==hltConfig_.size());

  // get trigger iso val collections, if requested
  if ( getHLTIsoVals_ ) {
    iEvent.getByLabel(isoValMapGlbTag_,isoValMapGlbHandle_);
    if (!isoValMapGlbHandle_.isValid()) {
      cout << "DilepTrigAnalyzerRECO::analyzeTrigger: Error in getting isoValMapGlb product from Event!" << endl;
    }

    iEvent.getByLabel(isoValMapTrkTag_,isoValMapTrkHandle_);
    if (!isoValMapTrkHandle_.isValid()) {
      cout << "DilepTrigAnalyzerRECO::analyzeTrigger: Error in getting isoValMapTrk product from Event!" << endl;
    }
  } 

  // retrieve necessary containers
  iEvent.getByLabel(vtxInputTag_, vertexHandle_);
  iEvent.getByLabel(electronsInputTag_, elsHandle_);
  iEvent.getByLabel( muonsInputTag_ , musHandle_ );

  // for printing out HLT PF Cand info
  if (dumpHLTPFCands_) {
    iEvent.getByLabel(hltPFCandsGlbTag_, hltPFCandsGlbHandle_);
    iEvent.getByLabel(hltPFCandsTrkTag_, hltPFCandsTrkHandle_);
    iEvent.getByLabel(offPFCandsTag_, offPFCandsHandle_);
    iEvent.getByLabel(hltVtxInputTag_, hltVertexHandle_);
  }
  
  // analyze this event for the triggers requested
  if (triggerName_=="dimu") {
    trigpass_results_ = 0;
    trigpass_results_offdilep_ = 0;
    for (unsigned int itrig=0; itrig <= (unsigned int)mmtk; ++itrig) {
      bool pass = analyzeTrigger(iEvent,iSetup,(hltTrigs)itrig);
      if (pass) trigpass_results_ |= 1 << itrig;
    }
    h_results_mm_->Fill(trigpass_results_);
    h_results_offdilep_mm_->Fill(trigpass_results_offdilep_);
  } else if (triggerName_=="emu") {
    trigpass_results_ = 0;
    trigpass_results_offdilep_ = 0;
    for (unsigned int itrig=(unsigned int)emi; itrig <= (unsigned int)me ; ++itrig) {
      bool pass = analyzeTrigger(iEvent,iSetup,(hltTrigs)itrig);
      if (pass) trigpass_results_ |= 1 << (itrig - (unsigned int)emi);
    }
    h_results_em_->Fill(trigpass_results_);
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
  if ( triggerEnum <= mmtk ) {
    ismm = true;
  } else if ( triggerEnum <= em ) {
    isem = true;
  } else if (triggerEnum <= me) {
    isme = true;
  }

  if (!ismm && !isem && !isme) {
    cout << "DilepTrigAnalyzerRECO::analyzeTrigger: triggerEnum: " << triggerEnum
	 << " not recognized, aborting.." << endl;
    return false;
  }

  // if non-iso trigger, get iso trigger result to compare to
  hltTrigs isoTriggerEnum = notrig;
  if ( triggerEnum == mm || triggerEnum == mmtk || triggerEnum == em || triggerEnum == me ) {
    isoTriggerEnum = hltTrigs(int(triggerEnum)-1);
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

  // // Results from TriggerEventWithRefs product
  // electronIds_.clear();
  // electronRefs_.clear();
  // muonIds_.clear();
  // muonRefs_.clear();

  bool foundMuons = false;
  bool foundElectrons = ismm; // don't look for electrons on dimu trigger

  StudyLepton hlt_el_lead;
  StudyLepton hlt_el_subl;
  hlt_el_lead.type = 11;
  hlt_el_subl.type = 11;
  hlt_el_lead.isHLT = true;
  hlt_el_subl.isHLT = true;
  int hlt_el_lead_idx = -1;
  int hlt_el_subl_idx = -1;
  //  int hlt_el_filter_idx = -1;

  StudyLepton hlt_mu_lead;
  StudyLepton hlt_mu_subl;
  StudyLepton hlt_mu_third;
  hlt_mu_lead.type = 13;
  hlt_mu_subl.type = 13;
  hlt_mu_third.type = 13;
  hlt_mu_lead.isHLT = true;
  hlt_mu_subl.isHLT = true;
  hlt_mu_third.isHLT = true;
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

    electronIds_.clear();
    electronRefs_.clear();
    triggerEventWithRefsHandle_->getObjects(filterIndex,electronIds_,electronRefs_);
    const unsigned int nElectrons(electronIds_.size());
    if (nElectrons>0) {
      if (verbose_) cout << "   Electrons: " << nElectrons << "  - the objects: # id pt eta phi vz id key" << endl;
      for (unsigned int i=0; i!=nElectrons; ++i) {

	if (!foundElectrons) {
	  LorentzVector lv(electronRefs_.at(i)->p4());

	  if ( (hlt_el_lead_idx == -1) || (lv.pt() > hlt_el_lead.lv.pt()) ) {
	    hlt_el_subl_idx = hlt_el_lead_idx;
	    hlt_el_subl.lv = hlt_el_lead.lv;
	    hlt_el_lead_idx = (int)i;
	    hlt_el_lead.lv = lv;
	  } else if ( (hlt_el_subl_idx == -1) || (lv.pt() > hlt_el_subl.lv.pt()) ) {
	    // check dR to remove exact duplicates
	    if (ROOT::Math::VectorUtil::DeltaR(hlt_el_lead.lv,lv) > 0.001) {
	      hlt_el_subl_idx = (int)i;
	      hlt_el_subl.lv = lv;
	    }
	  }
	} // if !foundElectrons

	if (verbose_) {
  	  cout << "   " << i << " " << electronIds_.at(i)
  	       << " " << electronRefs_.at(i)->pt()
  	       << " " << electronRefs_.at(i)->eta()
  	       << " " << electronRefs_.at(i)->phi()
  	       << " " << electronRefs_.at(i)->vz()
  	       << " " << electronRefs_.at(i).id()
  	       << " " << electronRefs_.at(i).key()
  	       << endl;
	}
      } // electrons loop
    } // if nElectrons>0


    muonIds_.clear();
    muonRefs_.clear();
    triggerEventWithRefsHandle_->getObjects(filterIndex,muonIds_,muonRefs_);
    const unsigned int nMuons(muonIds_.size());
    if (nMuons>0) {

      bool isovals = false;

      // mm / em / me / mmtk cases
      if ( getHLTIsoVals_ && 
	   ( (moduleLabel.find("hltDiMuonGlb17Glb8DzFiltered0p2") != string::npos) ||
	     (moduleLabel.find("hltL1Mu12EG7L3MuFiltered17") != string::npos) ||
	     (moduleLabel.find("hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8") != string::npos) || 
	     (moduleLabel.find("hltDiMuonGlb17Trk8DzFiltered0p2") != string::npos) ) ) {
	//	if (verbose_) cout << "found the module before/cutting on iso! hurray!" << endl;
	//	iEvent.getByLabel(isoValMapGlbTag_,isoValMapGlbHandle_);
	bool success = true;
	if (!isoValMapGlbHandle_.isValid() && (triggerEnum == mm || triggerEnum == mmi)) {
	  cout << "DilepTrigAnalyzerRECO::analyzeTrigger: Error in getting isoValMapGlb product from Event!" << endl;
	  success = false;
	}
	// else if (isoValMapGlbHandle_->idSize() != nMuons) {
	//   cout << "DilepTrigAnalyzerRECO::analyzeTrigger: WARNING: isoValMapGlb size == " << isoValMapGlbHandle_->idSize()
	//        << ", nMuons == " << nMuons << ". Contents: " << endl;
	//   for (unsigned int j=0; j < isoValMapGlbHandle_->ids().size(); ++j) {
	//     cout << "   " << isoValMapGlbHandle_->ids().at(j).first << endl;
	//   }
	// }

	//	iEvent.getByLabel(isoValMapTrkTag_,isoValMapTrkHandle_);
	if (!isoValMapTrkHandle_.isValid() && (triggerEnum == mmtk || triggerEnum == mmitk)) {
	  cout << "DilepTrigAnalyzerRECO::analyzeTrigger: Error in getting isoValMapTrk product from Event!" << endl;
	  success = false;
	}
	// else if (isoValMapTrkHandle_->idSize() != nMuons) {
	//   cout << "DilepTrigAnalyzerRECO::analyzeTrigger: WARNING: isoValMapTrk size == " << isoValMapTrkHandle_->idSize()
	//        << ", nMuons == " << nMuons << ". Contents: " << endl;
	//   for (unsigned int j=0; j < isoValMapTrkHandle_->ids().size(); ++j) {
	//     cout << "   " << isoValMapTrkHandle_->ids().at(j).first << endl;
	//   }
	// }

	if (success) isovals = true;
      } // filter before/after iso

      if (verbose_) {
	cout << "   Muons: " << nMuons << ", MuonRefs: " << muonRefs_.size()
	     << "  - the objects: # id pt eta phi vz id key" << endl;
      }
      for (unsigned int i=0; i!=nMuons; ++i) {

	float iso = -99.;
	if (!foundMuons) {
	  LorentzVector lv(muonRefs_.at(i)->p4());
	  if (isovals) {
	    // check first that id is in map
	    //   if not in map, try other map
	    //   if not in either map, raise hell
	    if ( (triggerEnum == mm || triggerEnum == mmi) && isoValMapGlbHandle_->contains(muonRefs_.at(i).id()) ) {
	      const edm::ValueMap<float> ::value_type & muonDeposit = (*(isoValMapGlbHandle_))[muonRefs_[i]];
	      iso = muonDeposit;
	    }
	    else if ( (triggerEnum == mmtk || triggerEnum == mmitk) && isoValMapTrkHandle_->contains(muonRefs_.at(i).id()) ) {
	      const edm::ValueMap<float> ::value_type & muonDeposit = (*(isoValMapTrkHandle_))[muonRefs_[i]];
	      iso = muonDeposit;
	    } else {
	      cout << "DilepTrigAnalyzerRECO::analyzeTrigger: ERROR: couldn't find id in iso val maps."
		   << " id = " << muonRefs_.at(i).id() << endl
		   << "  Glb Map contains:" << endl;
	      for (unsigned int j=0; j < isoValMapGlbHandle_->ids().size(); ++j) {
		cout << "   " << isoValMapGlbHandle_->ids().at(j).first << endl;
	      }
	      cout << "  GlbTrk Map contains:" << endl;
	      for (unsigned int j=0; j < isoValMapTrkHandle_->ids().size(); ++j) {
		cout << "   " << isoValMapTrkHandle_->ids().at(j).first << endl;
	      }
	      cout << "  GlbTrk Map contains:" << endl;
	    }
	  }

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
	  if (isovals) cout << ", isoval: " << iso;
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

    //   // electrons
    //   if (!foundElectrons && (VIDS[i] == 82)) {
    // 	if ( (hlt_el_lead_idx == -1) || (lv.pt() > hlt_el_lead.pt()) ) {
    // 	  hlt_el_subl_idx = hlt_el_lead_idx;
    // 	  hlt_el_subl = hlt_el_lead;
    // 	  hlt_el_lead_idx = (int)i;
    // 	  hlt_el_lead = lv;
    // 	} else if ( (hlt_el_subl_idx == -1) || (lv.pt() > hlt_el_subl.pt()) ) {
    // 	  // check dR to remove exact duplicates
    // 	  if (ROOT::Math::VectorUtil::DeltaR(hlt_el_lead,lv) > 0.001) {
    // 	    hlt_el_subl_idx = (int)i;
    // 	    hlt_el_subl = lv;
    // 	  }
    // 	}
    //   } // hlt electrons

    //   // muons
    //   else if (!foundMuons && (VIDS[i] == 83)) {
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

    if (hlt_el_lead_idx >= 0) {
      foundElectrons = true;
      //      hlt_el_filter_idx = filterIndex;
    }

    if (hlt_mu_lead_idx >= 0) {
      foundMuons = true;
      //      hlt_mu_filter_idx = filterIndex;
    }

    if (foundElectrons && foundMuons) break;
  } // backwards loop on modules

  if ((hlt_el_lead_idx == -1) && (hlt_mu_lead_idx == -1)) {
    cout << "DilepTrigAnalyzerRECO::analyzeTrigger: no valid trigger leptons!" << endl;
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

  StudyLepton off_el_lead;
  StudyLepton off_el_subl;
  off_el_lead.type = 11;
  off_el_subl.type = 11;
  off_el_lead.isHLT = false;
  off_el_subl.isHLT = false;
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
    if ( (hlt_el_lead_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(lv,hlt_el_lead.lv) < dr_trigmatch) ) match = true;
    else if ( (hlt_el_subl_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(lv,hlt_el_subl.lv) < dr_trigmatch) ) match = true;
    if (!match) continue;

    // pt ordering
    if ( (off_el_lead_idx == -1) || (lv.pt() > off_el_lead.lv.pt()) ) {
      off_el_subl_idx = off_el_lead_idx;
      off_el_subl.lv = off_el_lead.lv;
      off_el_lead_idx = (int)elsIndex;
      off_el_lead.lv = lv;
    } else if ( (off_el_subl_idx == -1) || (lv.pt() > off_el_subl.lv.pt()) ) {
      // check dR to remove exact duplicates
      if (ROOT::Math::VectorUtil::DeltaR(off_el_lead.lv,lv) > 0.001) {
	off_el_subl_idx = (int)elsIndex;
	off_el_subl.lv = lv;
      }
    }

  } // loop on reco gsf electrons

  //-------------------------------------
  //   reco muons 
  //-------------------------------------

  StudyLepton off_mu_lead;
  StudyLepton off_mu_subl;
  off_mu_lead.type = 13;
  off_mu_subl.type = 13;
  off_mu_lead.isHLT = false;
  off_mu_subl.isHLT = false;
  int off_mu_lead_idx = -1;
  int off_mu_subl_idx = -1;

  StudyLepton off_mu_tight_lead;
  StudyLepton off_mu_tight_subl;
  off_mu_tight_lead.type = 13;
  off_mu_tight_subl.type = 13;
  off_mu_tight_lead.isHLT = false;
  off_mu_tight_subl.isHLT = false;
  int off_mu_tight_lead_idx = -1;
  int off_mu_tight_subl_idx = -1;

  StudyLepton off_mu_trigiso_lead;
  StudyLepton off_mu_trigiso_subl;
  off_mu_trigiso_lead.type = 13;
  off_mu_trigiso_subl.type = 13;
  off_mu_trigiso_lead.isHLT = false;
  off_mu_trigiso_subl.isHLT = false;
  int off_mu_trigiso_lead_idx = -1;
  int off_mu_trigiso_subl_idx = -1;

  StudyLepton off_mu_liso_lead;
  StudyLepton off_mu_liso_subl;
  off_mu_liso_lead.type = 13;
  off_mu_liso_subl.type = 13;
  off_mu_liso_lead.isHLT = false;
  off_mu_liso_subl.isHLT = false;
  int off_mu_liso_lead_idx = -1;
  int off_mu_liso_subl_idx = -1;

  StudyLepton off_mu_tliso_lead;
  StudyLepton off_mu_tliso_subl;
  off_mu_tliso_lead.type = 13;
  off_mu_tliso_subl.type = 13;
  off_mu_tliso_lead.isHLT = false;
  off_mu_tliso_subl.isHLT = false;
  int off_mu_tliso_lead_idx = -1;
  int off_mu_tliso_subl_idx = -1;

  StudyLepton off_mu_tiso_lead;
  StudyLepton off_mu_tiso_subl;
  off_mu_tiso_lead.type = 13;
  off_mu_tiso_subl.type = 13;
  off_mu_tiso_lead.isHLT = false;
  off_mu_tiso_subl.isHLT = false;
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
    if ( (hlt_mu_lead_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(lv,hlt_mu_lead.lv) < dr_trigmatch) ) match = true;
    else if ( (hlt_mu_subl_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(lv,hlt_mu_subl.lv) < dr_trigmatch) ) match = true;
    else if ( (hlt_mu_third_idx >= 0) && (ROOT::Math::VectorUtil::DeltaR(lv,hlt_mu_third.lv) < dr_trigmatch) ) match = true;
    if (!match) continue;

    // basic dz cut to remove muons from large z
    bool pass_dz = bool(fabs(muon->muonBestTrack()->dz(firstGoodVertex->position())) < 0.5);
    if (!pass_dz) continue;

    // min pt cut
    if (muon->pt() < offSublPt_) continue;

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

  ///////////////////////////
  // TransientTrackBuilder //
  ///////////////////////////
  ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

  // loop again on duplicate-cleaned muon collection
  muonIndex = 0;
  muons_end = muons_good.end();  // Iterator
  for ( MuonCollection::const_iterator muon = muons_good.begin(); muon != muons_end; ++muon, ++muonIndex ) {
    LorentzVector lv(muon->p4());

    bool pass_loose = muon::isLooseMuon(*muon);
    bool pass_tight = muon::isTightMuon(*muon,*firstGoodVertex);
    float trkiso = muon->pfIsolationR03().sumChargedHadronPt;
    bool pass_trkiso_trig = bool(trkiso/lv.pt() < 0.4);
    float pfiso = muonPFiso(*muon);
    bool pass_iso_loose = bool(pfiso/lv.pt() < 0.4);
    bool pass_iso_tight = bool(pfiso/lv.pt() < 0.15);
    // get 3d impact parameter significance, for ZZ sel
    const TrackRef siTrack  = muon->innerTrack();
    float ip3d_val = -999.;
    float ip3d_err = -999.;
    float ip3d_sig = -999.;
    float dxy = -999.;
    if ( siTrack.isNonnull() && firstGoodVertex != vertexCollection->end() ) {
      TransientTrack tt       = theTTBuilder->build( siTrack );
      Measurement1D ip3D      = IPTools::absoluteImpactParameter3D( tt, *firstGoodVertex ).second;
      ip3d_val = ip3D.value();
      ip3d_err = ip3D.error();
      ip3d_sig = ip3d_val/ip3d_err;
      dxy = siTrack->dxy(firstGoodVertex->position());
    } 
    bool pass_ip3d = bool(fabs(ip3d_sig) < 4.);
    bool pass_dxy = bool(fabs(dxy) < 0.02);

    // pt ordering
    if ( ((off_mu_lead_idx == -1) || (lv.pt() > off_mu_lead.lv.pt())) && (lv.pt() > offLeadPt_) ) {
      if (off_mu_lead_idx != -1) {
	off_mu_subl_idx = off_mu_lead_idx;
	off_mu_subl = off_mu_lead;
      }

      off_mu_lead_idx = (int)muonIndex;
      off_mu_lead.lv = lv;
      off_mu_lead.trkiso = trkiso;
      off_mu_lead.pfiso = pfiso;
      off_mu_lead.vz = muon->vz();
      off_mu_lead.charge = muon->charge();
    } else if ( ((off_mu_subl_idx == -1) || (lv.pt() > off_mu_subl.lv.pt())) && (lv.pt() > offSublPt_) ) {
      off_mu_subl_idx = (int)muonIndex;
      off_mu_subl.lv = lv;
      off_mu_subl.trkiso = trkiso;
      off_mu_subl.pfiso = pfiso;
      off_mu_subl.vz = muon->vz();
      off_mu_subl.charge = muon->charge();
    }

    if (pass_tight) {
      // pt ordering
      if ( (((off_mu_tight_lead_idx == -1) || (lv.pt() > off_mu_tight_lead.lv.pt())))  && (lv.pt() > offLeadPt_) ) {
	if (off_mu_tight_lead_idx != -1) {
	  off_mu_tight_subl_idx = off_mu_tight_lead_idx;
	  off_mu_tight_subl = off_mu_tight_lead;
	}

	off_mu_tight_lead_idx = (int)muonIndex;
	off_mu_tight_lead.lv = lv;
	off_mu_tight_lead.trkiso = trkiso;
	off_mu_tight_lead.pfiso = pfiso;
	off_mu_tight_lead.vz = muon->vz();
	off_mu_tight_lead.charge = muon->charge();
      } else if ( ((off_mu_tight_subl_idx == -1) || (lv.pt() > off_mu_tight_subl.lv.pt()))  && (lv.pt() > offSublPt_) ) {
	off_mu_tight_subl_idx = (int)muonIndex;
	off_mu_tight_subl.lv = lv;
	off_mu_tight_subl.trkiso = trkiso;
	off_mu_tight_subl.pfiso = pfiso;
	off_mu_tight_subl.vz = muon->vz();
        off_mu_tight_subl.charge = muon->charge();
      }
    }

    //    if (pass_tight && pass_trkiso_trig) {
    if (pass_loose && pass_trkiso_trig) {
      // pt ordering
      if ( ((off_mu_trigiso_lead_idx == -1) || (lv.pt() > off_mu_trigiso_lead.lv.pt()))  && (lv.pt() > offLeadPt_)) {
	if (off_mu_trigiso_lead_idx != -1) {
	  off_mu_trigiso_subl_idx = off_mu_trigiso_lead_idx;
	  off_mu_trigiso_subl = off_mu_trigiso_lead;
	}

	off_mu_trigiso_lead_idx = (int)muonIndex;
	off_mu_trigiso_lead.lv = lv;
	off_mu_trigiso_lead.trkiso = trkiso;
	off_mu_trigiso_lead.pfiso = pfiso;
	off_mu_trigiso_lead.vz = muon->vz();
        off_mu_trigiso_lead.charge = muon->charge();
      } else if ( ((off_mu_trigiso_subl_idx == -1) || (lv.pt() > off_mu_trigiso_subl.lv.pt()))  && (lv.pt() > offSublPt_)) {
	off_mu_trigiso_subl_idx = (int)muonIndex;
	off_mu_trigiso_subl.lv = lv;
	off_mu_trigiso_subl.trkiso = trkiso;
	off_mu_trigiso_subl.pfiso = pfiso;
	off_mu_trigiso_subl.vz = muon->vz();
        off_mu_trigiso_subl.charge = muon->charge();
      }
    }

    if (pass_loose && pass_ip3d && pass_iso_loose) {
      // pt ordering
      if ( ((off_mu_liso_lead_idx == -1) || (lv.pt() > off_mu_liso_lead.lv.pt()))  && (lv.pt() > offLeadPt_)) {
	if (off_mu_liso_lead_idx != -1) {
	  off_mu_liso_subl_idx = off_mu_liso_lead_idx;
	  off_mu_liso_subl = off_mu_liso_lead;
	}

	off_mu_liso_lead_idx = (int)muonIndex;
	off_mu_liso_lead.lv = lv;
	off_mu_liso_lead.trkiso = trkiso;
	off_mu_liso_lead.pfiso = pfiso;
	off_mu_liso_lead.vz = muon->vz();
        off_mu_liso_lead.charge = muon->charge();
      } else if ( ((off_mu_liso_subl_idx == -1) || (lv.pt() > off_mu_liso_subl.lv.pt()))  && (lv.pt() > offSublPt_)) {
	off_mu_liso_subl_idx = (int)muonIndex;
	off_mu_liso_subl.lv = lv;
	off_mu_liso_subl.trkiso = trkiso;
	off_mu_liso_subl.pfiso = pfiso;
	off_mu_liso_subl.vz = muon->vz();
        off_mu_liso_subl.charge = muon->charge();
      }
    }

    if (pass_tight && pass_iso_loose) {
      // pt ordering
      if ( ((off_mu_tliso_lead_idx == -1) || (lv.pt() > off_mu_tliso_lead.lv.pt()))  && (lv.pt() > offLeadPt_)) {
	if (off_mu_tliso_lead_idx != -1) {
	  off_mu_tliso_subl_idx = off_mu_tliso_lead_idx;
	  off_mu_tliso_subl = off_mu_tliso_lead;
	}

	off_mu_tliso_lead_idx = (int)muonIndex;
	off_mu_tliso_lead.lv = lv;
	off_mu_tliso_lead.trkiso = trkiso;
	off_mu_tliso_lead.pfiso = pfiso;
	off_mu_tliso_lead.vz = muon->vz();
        off_mu_tliso_lead.charge = muon->charge();
      } else if ( ((off_mu_tliso_subl_idx == -1) || (lv.pt() > off_mu_tliso_subl.lv.pt()))  && (lv.pt() > offSublPt_)) {
	off_mu_tliso_subl_idx = (int)muonIndex;
	off_mu_tliso_subl.lv = lv;
	off_mu_tliso_subl.trkiso = trkiso;
	off_mu_tliso_subl.pfiso = pfiso;
	off_mu_tliso_subl.vz = muon->vz();
        off_mu_tliso_subl.charge = muon->charge();
      }
    }

    if (pass_tight && pass_iso_tight && pass_dxy) {
      // pt ordering
      if ( ((off_mu_tiso_lead_idx == -1) || (lv.pt() > off_mu_tiso_lead.lv.pt()))  && (lv.pt() > offLeadPt_)) {
	if (off_mu_tiso_lead_idx != -1) {
	  off_mu_tiso_subl_idx = off_mu_tiso_lead_idx;
	  off_mu_tiso_subl = off_mu_tiso_lead;
	}

	off_mu_tiso_lead_idx = (int)muonIndex;
	off_mu_tiso_lead.lv = lv;
	off_mu_tiso_lead.trkiso = trkiso;
	off_mu_tiso_lead.pfiso = pfiso;
	off_mu_tiso_lead.vz = muon->vz();
        off_mu_tiso_lead.charge = muon->charge();
      } else if ( ((off_mu_tiso_subl_idx == -1) || (lv.pt() > off_mu_tiso_subl.lv.pt()))  && (lv.pt() > offSublPt_)) {
	off_mu_tiso_subl_idx = (int)muonIndex;
	off_mu_tiso_subl.lv = lv;
	off_mu_tiso_subl.trkiso = trkiso;
	off_mu_tiso_subl.pfiso = pfiso;
	off_mu_tiso_subl.vz = muon->vz();
        off_mu_tiso_subl.charge = muon->charge();
      }
    }

  } // loop on duplicate-cleaned muons


  //-------------------------------------
  //   reco lepton plots
  //-------------------------------------

  // make plots based on trigger path
  if (ismm) {
    fillHists(off_mu_lead,off_mu_subl,triggerShort,false);
    fillHistsRecoHLT(off_mu_lead,off_mu_subl,hlt_mu_lead,hlt_mu_subl,triggerShort);
    if (off_mu_tight_lead_idx >= 0) {
      fillHists(off_mu_tight_lead,off_mu_tight_subl,triggerShort+"_tight",false);
      fillHistsRecoHLT(off_mu_tight_lead,off_mu_tight_subl,hlt_mu_lead,hlt_mu_subl,triggerShort+"_tight");
    }
    if (off_mu_trigiso_lead_idx >= 0) {
      fillHists(off_mu_trigiso_lead,off_mu_trigiso_subl,triggerShort+"_trigiso",false);
      fillHistsRecoHLT(off_mu_trigiso_lead,off_mu_trigiso_subl,hlt_mu_lead,hlt_mu_subl,triggerShort+"_trigiso");
      if (off_mu_trigiso_subl_idx >= 0) trigpass_results_offdilep_ |= 1 << int(triggerEnum);
    }
    if (off_mu_liso_lead_idx >= 0) {
      fillHists(off_mu_liso_lead,off_mu_liso_subl,triggerShort+"_liso",false);
      fillHistsRecoHLT(off_mu_liso_lead,off_mu_liso_subl,hlt_mu_lead,hlt_mu_subl,triggerShort+"_liso");
    }
    if (off_mu_tliso_lead_idx >= 0) {
      fillHists(off_mu_tliso_lead,off_mu_tliso_subl,triggerShort+"_tliso",false);
      fillHistsRecoHLT(off_mu_tliso_lead,off_mu_tliso_subl,hlt_mu_lead,hlt_mu_subl,triggerShort+"_tliso");
    }
    if (off_mu_tiso_lead_idx >= 0) {
      fillHists(off_mu_tiso_lead,off_mu_tiso_subl,triggerShort+"_tiso",false);
      fillHistsRecoHLT(off_mu_tiso_lead,off_mu_tiso_subl,hlt_mu_lead,hlt_mu_subl,triggerShort+"_tiso");
    }
  } // ismm

  else if (isem) {
    fillHists(off_el_lead,off_mu_lead,triggerShort,false);
    fillHistsRecoHLT(off_el_lead,off_mu_lead,hlt_el_lead,hlt_mu_lead,triggerShort);
    if (off_mu_tight_lead_idx >= 0) {
      fillHists(off_el_lead,off_mu_tight_lead,triggerShort+"_tight",false);
      fillHistsRecoHLT(off_el_lead,off_mu_tight_lead,hlt_el_lead,hlt_mu_lead,triggerShort+"_tight");
    }
    if (off_mu_trigiso_lead_idx >= 0) {
      fillHists(off_el_lead,off_mu_trigiso_lead,triggerShort+"_trigiso",false);
      fillHistsRecoHLT(off_el_lead,off_mu_trigiso_lead,hlt_el_lead,hlt_mu_lead,triggerShort+"_trigiso");
    }
    if (off_mu_tiso_lead_idx >= 0) {
      fillHists(off_el_lead,off_mu_tiso_lead,triggerShort+"_tiso",false);
      fillHistsRecoHLT(off_el_lead,off_mu_tiso_lead,hlt_el_lead,hlt_mu_lead,triggerShort+"_tiso");
    }
  } // isem

  else if (isme) {
    fillHists(off_mu_lead,off_el_lead,triggerShort,false);
    fillHistsRecoHLT(off_mu_lead,off_el_lead,hlt_mu_lead,hlt_el_lead,triggerShort);
    if (off_mu_tight_lead_idx >= 0) {
      fillHists(off_mu_tight_lead,off_el_lead,triggerShort+"_tight",false);
      fillHistsRecoHLT(off_mu_tight_lead,off_el_lead,hlt_mu_lead,hlt_el_lead,triggerShort+"_tight");
    }
    if (off_mu_trigiso_lead_idx >= 0) {
      fillHists(off_mu_trigiso_lead,off_el_lead,triggerShort+"_trigiso",false);
      fillHistsRecoHLT(off_mu_trigiso_lead,off_el_lead,hlt_mu_lead,hlt_el_lead,triggerShort+"_trigiso");
    }
    if (off_mu_tiso_lead_idx >= 0) {
      fillHists(off_mu_tiso_lead,off_el_lead,triggerShort+"_tiso",false);
      fillHistsRecoHLT(off_mu_tiso_lead,off_el_lead,hlt_mu_lead,hlt_el_lead,triggerShort+"_tiso");
    }
  } // isme


  // printout for cases where trigger is inefficient
  if ( verbose_ && ismm && (off_mu_trigiso_lead_idx >= 0) && (off_mu_trigiso_subl_idx >= 0) ) {
    // check that this is the noniso trigger, and that the iso trigger failed
    LorentzVector dilep = off_mu_trigiso_lead.lv + off_mu_trigiso_subl.lv;
    float deltaR = ROOT::Math::VectorUtil::DeltaR(off_mu_trigiso_lead.lv,off_mu_trigiso_subl.lv);
    bool ineff = false;
    if ( (isoTriggerEnum != notrig) && ((trigpass_results_ & (1 << isoTriggerEnum)) == 0) ) {
      cout << "-- HLT inefficiency, event level! Non-iso trigger: " << triggerName;
      ineff = true;
    }
    else if ( (isoTriggerEnum != notrig) && ((trigpass_results_offdilep_ & (1 << isoTriggerEnum)) == 0) ) {
      cout << "++ HLT inefficiency, offline object level! Non-iso trigger: " << triggerName;
      ineff = true;
    }
    if (ineff) {
      cout << ", run: " << iEvent.id().run()
	   << ", event: " << iEvent.id().event() << ", mass: " << dilep.M() << ", deltaR: " << deltaR << endl
	   << "     hlt mu1 pt: " << hlt_mu_lead.lv.pt() << ", iso: " << hlt_mu_lead.trkiso
	   << ", hlt mu2 pt: " << hlt_mu_subl.lv.pt() << ", iso: " << hlt_mu_subl.trkiso << endl
	   << "     off mu1 pt: " << off_mu_trigiso_lead.lv.pt() << ", iso: " << off_mu_trigiso_lead.trkiso
	   << ", off mu2 pt: " << off_mu_trigiso_subl.lv.pt() << ", iso: " << off_mu_trigiso_subl.trkiso << endl;
      // loop over pf cands near the hlt muons and print them out..
      if (dumpHLTPFCands_) {
	edm::Handle<reco::PFCandidateCollection> hltPFCandsHandle;
	if (triggerEnum == mm) hltPFCandsHandle = hltPFCandsGlbHandle_;
	else hltPFCandsHandle = hltPFCandsTrkHandle_;
	const VertexCollection* hltVertexCollection = hltVertexHandle_.product();

	// find hlt pf cands near hlt muons
        PFCandidateCollection::const_iterator hlt_pfcands_end = hltPFCandsHandle->end();  // Iterator
	PFCandidateCollection hlt_pfcands_lead;
	PFCandidateCollection hlt_pfcands_subl;
	for ( PFCandidateCollection::const_iterator hlt_pfcand = hltPFCandsHandle->begin(); hlt_pfcand != hlt_pfcands_end; ++hlt_pfcand ) {
	  if (ROOT::Math::VectorUtil::DeltaR(hlt_mu_lead.lv,hlt_pfcand->p4()) < 0.3) {
	    hlt_pfcands_lead.push_back(*hlt_pfcand);
	  }
	  if (ROOT::Math::VectorUtil::DeltaR(hlt_mu_subl.lv,hlt_pfcand->p4()) < 0.3) {
	    hlt_pfcands_subl.push_back(*hlt_pfcand);
	  }
	} // loop on hlt pf cands

	// find off pf cands near hlt muons
        PFCandidateCollection::const_iterator off_pfcands_end = offPFCandsHandle_->end();  // Iterator
	PFCandidateCollection off_pfcands_lead;
	PFCandidateCollection off_pfcands_subl;
	for ( PFCandidateCollection::const_iterator off_pfcand = offPFCandsHandle_->begin(); off_pfcand != off_pfcands_end; ++off_pfcand ) {
	  if (abs(off_pfcand->charge()) != 1) continue;

	  if (ROOT::Math::VectorUtil::DeltaR(hlt_mu_lead.lv,off_pfcand->p4()) < 0.3) {
	    off_pfcands_lead.push_back(*off_pfcand);
	  }
	  if (ROOT::Math::VectorUtil::DeltaR(hlt_mu_subl.lv,off_pfcand->p4()) < 0.3) {
	    off_pfcands_subl.push_back(*off_pfcand);
	  }
	} // loop on off pf cands

	// dump hlt pf cands near lead hlt muon
	cout << "   ---- hlt pf cands near leading hlt muon:" << endl;
        PFCandidateCollection::const_iterator hlt_pfcands_lead_end = hlt_pfcands_lead.end();  // Iterator
	for ( PFCandidateCollection::const_iterator hlt_pfcand = hlt_pfcands_lead.begin(); hlt_pfcand != hlt_pfcands_lead_end; ++hlt_pfcand ) {
	  float dR = ROOT::Math::VectorUtil::DeltaR(hlt_mu_lead.lv,hlt_pfcand->p4());
	  int vtx = chargedHadronVertex(*hltVertexCollection,*hlt_pfcand);
	  float dz = hlt_mu_lead.vz - hlt_pfcand->vz();
	  if (fabs(dz) > 1.0) continue;
	  else if (fabs(dz) < 0.2) cout << "    * ";
	  else cout << "      ";
	  cout << "pt: " << hlt_pfcand->pt() << ", eta: " << hlt_pfcand->eta() << ", phi: " << hlt_pfcand->phi()
	       << ", id: " << hlt_pfcand->particleId() << ", vz: " << hlt_pfcand->vz() << ", vtx: " << vtx
	       << ", nhits: " << hlt_pfcand->trackRef()->numberOfValidHits() << ", dR: " << dR << endl;
	}

	// dump off pf cands near lead hlt muon
	cout << "   ++++ off pf cands near leading hlt muon:" << endl;
        PFCandidateCollection::const_iterator off_pfcands_lead_end = off_pfcands_lead.end();  // Iterator
	for ( PFCandidateCollection::const_iterator off_pfcand = off_pfcands_lead.begin(); off_pfcand != off_pfcands_lead_end; ++off_pfcand ) {
	  float dR = ROOT::Math::VectorUtil::DeltaR(hlt_mu_lead.lv,off_pfcand->p4());
	  int vtx = chargedHadronVertex(*vertexCollection,*off_pfcand);
	  float dz = hlt_mu_lead.vz - off_pfcand->vz();
	  if (fabs(dz) > 1.0) continue;
	  if (((vtx == 0) || (vtx == -1)) && (off_pfcand->particleId() == 1)) cout << "    * ";
	  else cout << "      ";
	  cout << "pt: " << off_pfcand->pt() << ", eta: " << off_pfcand->eta() << ", phi: " << off_pfcand->phi()
	       << ", id: " << off_pfcand->particleId() << ", vz: " << off_pfcand->vz() << ", vtx: " << vtx
	       << ", nhits: " << off_pfcand->trackRef()->numberOfValidHits() << ", dR: " << dR << endl;
	}

	// dump hlt pf cands near subl hlt muon
	cout << "   ---- hlt pf cands near subleading hlt muon:" << endl;
        PFCandidateCollection::const_iterator hlt_pfcands_subl_end = hlt_pfcands_subl.end();  // Iterator
	for ( PFCandidateCollection::const_iterator hlt_pfcand = hlt_pfcands_subl.begin(); hlt_pfcand != hlt_pfcands_subl_end; ++hlt_pfcand ) {
	  float dR = ROOT::Math::VectorUtil::DeltaR(hlt_mu_subl.lv,hlt_pfcand->p4());
	  int vtx = chargedHadronVertex(*hltVertexCollection,*hlt_pfcand);
	  float dz = hlt_mu_lead.vz - hlt_pfcand->vz();
	  if (fabs(dz) > 1.0) continue;
	  else if (fabs(dz) < 0.2) cout << "    * ";
	  else cout << "      ";
	  cout << "pt: " << hlt_pfcand->pt() << ", eta: " << hlt_pfcand->eta() << ", phi: " << hlt_pfcand->phi()
	       << ", id: " << hlt_pfcand->particleId() << ", vz: " << hlt_pfcand->vz() << ", vtx: " << vtx
	       << ", nhits: " << hlt_pfcand->trackRef()->numberOfValidHits() << ", dR: " << dR << endl;
	}

	// dump off pf cands near subl hlt muon
	cout << "   ++++ off pf cands near subleading hlt muon:" << endl;
        PFCandidateCollection::const_iterator off_pfcands_subl_end = off_pfcands_subl.end();  // Iterator
	for ( PFCandidateCollection::const_iterator off_pfcand = off_pfcands_subl.begin(); off_pfcand != off_pfcands_subl_end; ++off_pfcand ) {
	  float dR = ROOT::Math::VectorUtil::DeltaR(hlt_mu_subl.lv,off_pfcand->p4());
	  int vtx = chargedHadronVertex(*vertexCollection,*off_pfcand);
	  float dz = hlt_mu_lead.vz - off_pfcand->vz();
	  if (fabs(dz) > 1.0) continue;
	  if (((vtx == 0) || (vtx == -1)) && (off_pfcand->particleId() == 1)) cout << "    * ";
	  else cout << "      ";
	  cout << "pt: " << off_pfcand->pt() << ", eta: " << off_pfcand->eta() << ", phi: " << off_pfcand->phi()
	       << ", id: " << off_pfcand->particleId() << ", vz: " << off_pfcand->vz() << ", vtx: " << vtx
	       << ", nhits: " << off_pfcand->trackRef()->numberOfValidHits() << ", dR: " << dR << endl;
	}

      } // if dumpHLTPFCands
    } // if ineff

  }

  return true;
}

//____________________________________________________________________________
void DilepTrigAnalyzerRECO::bookHists(edm::Service<TFileService>& fs, const std::string& suffix, bool hlt) {

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  std::string hlt_suf("_hlt");

  if (hlt) {
    hists_1d_["h_lead_pt"+suf+hlt_suf] = fs->make<TH1F>(Form("h_lead_pt%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Leading p_{T} [GeV]" , 100 , 0. , 100. );
    hists_1d_["h_subl_pt"+suf+hlt_suf] = fs->make<TH1F>(Form("h_subl_pt%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Subleading p_{T} [GeV]" , 100 , 0. , 100. );
    hists_1d_["h_lead_eta"+suf+hlt_suf] = fs->make<TH1F>(Form("h_lead_eta%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Leading #eta" , 100 , -3. , 3. );
    hists_1d_["h_subl_eta"+suf+hlt_suf] = fs->make<TH1F>(Form("h_subl_eta%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Subleading #eta" , 100 , -3. , 3. );
    hists_1d_["h_mll"+suf+hlt_suf] = fs->make<TH1F>(Form("h_mll%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT M_{ll} [GeV]" , 150 , 0. , 150. );
    hists_1d_["h_dr"+suf+hlt_suf] = fs->make<TH1F>(Form("h_dr%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT #DeltaR" , 600 , 0. , 6. );

  hists_1d_["h_lead_abstrkiso"+suf+hlt_suf] = fs->make<TH1F>(Form("h_lead_abstrkiso%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Leading trkiso [GeV]" , 200 , 0. , 10. );
  hists_1d_["h_subl_abstrkiso"+suf+hlt_suf] = fs->make<TH1F>(Form("h_subl_abstrkiso%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Subleading trkiso [GeV]" , 200 , 0. , 10. );
  hists_1d_["h_lead_reltrkiso"+suf+hlt_suf] = fs->make<TH1F>(Form("h_lead_reltrkiso%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Leading trkiso / p_{T}" , 200 , 0. , 2. );
  hists_1d_["h_subl_reltrkiso"+suf+hlt_suf] = fs->make<TH1F>(Form("h_subl_reltrkiso%s%s",suf.c_str(),hlt_suf.c_str()) , "; HLT Subleading trkiso / p_{T}" , 200 , 0. , 2. );

  }

  hists_1d_["h_lead_pt"+suf] = fs->make<TH1F>(Form("h_lead_pt%s",suf.c_str()) , "; Leading p_{T} [GeV]" , 100 , 0. , 100. );
  hists_1d_["h_subl_pt"+suf] = fs->make<TH1F>(Form("h_subl_pt%s",suf.c_str()) , "; Subleading p_{T} [GeV]" , 100 , 0. , 100. );
  hists_1d_["h_lead_eta"+suf] = fs->make<TH1F>(Form("h_lead_eta%s",suf.c_str()) , "; Leading #eta" , 100 , -3. , 3. );
  hists_1d_["h_subl_eta"+suf] = fs->make<TH1F>(Form("h_subl_eta%s",suf.c_str()) , "; Subleading #eta" , 100 , -3. , 3. );
  hists_1d_["h_mll"+suf] = fs->make<TH1F>(Form("h_mll%s",suf.c_str()) , "; M_{ll} [GeV]" , 150 , 0. , 150. );
  hists_1d_["h_dr"+suf] = fs->make<TH1F>(Form("h_dr%s",suf.c_str()) , "; #DeltaR" , 600 , 0. , 6. );

  hists_1d_["h_lead_pt_os"+suf] = fs->make<TH1F>(Form("h_lead_pt_os%s",suf.c_str()) , "; Leading p_{T} [GeV]" , 100 , 0. , 100. );
  hists_1d_["h_subl_pt_os"+suf] = fs->make<TH1F>(Form("h_subl_pt_os%s",suf.c_str()) , "; Subleading p_{T} [GeV]" , 100 , 0. , 100. );
  hists_1d_["h_lead_eta_os"+suf] = fs->make<TH1F>(Form("h_lead_eta_os%s",suf.c_str()) , "; Leading #eta" , 100 , -3. , 3. );
  hists_1d_["h_subl_eta_os"+suf] = fs->make<TH1F>(Form("h_subl_eta_os%s",suf.c_str()) , "; Subleading #eta" , 100 , -3. , 3. );
  hists_1d_["h_mll_os"+suf] = fs->make<TH1F>(Form("h_mll_os%s",suf.c_str()) , "; M_{ll} [GeV]" , 150 , 0. , 150. );
  hists_1d_["h_dr_os"+suf] = fs->make<TH1F>(Form("h_dr_os%s",suf.c_str()) , "; #DeltaR" , 600 , 0. , 6. );

  hists_1d_["h_lead_pt_ss"+suf] = fs->make<TH1F>(Form("h_lead_pt_ss%s",suf.c_str()) , "; Leading p_{T} [GeV]" , 100 , 0. , 100. );
  hists_1d_["h_subl_pt_ss"+suf] = fs->make<TH1F>(Form("h_subl_pt_ss%s",suf.c_str()) , "; Subleading p_{T} [GeV]" , 100 , 0. , 100. );
  hists_1d_["h_lead_eta_ss"+suf] = fs->make<TH1F>(Form("h_lead_eta_ss%s",suf.c_str()) , "; Leading #eta" , 100 , -3. , 3. );
  hists_1d_["h_subl_eta_ss"+suf] = fs->make<TH1F>(Form("h_subl_eta_ss%s",suf.c_str()) , "; Subleading #eta" , 100 , -3. , 3. );
  hists_1d_["h_mll_ss"+suf] = fs->make<TH1F>(Form("h_mll_ss%s",suf.c_str()) , "; M_{ll} [GeV]" , 150 , 0. , 150. );
  hists_1d_["h_dr_ss"+suf] = fs->make<TH1F>(Form("h_dr_ss%s",suf.c_str()) , "; #DeltaR" , 600 , 0. , 6. );

  hists_1d_["h_lead_abstrkiso"+suf] = fs->make<TH1F>(Form("h_lead_abstrkiso%s",suf.c_str()) , "; Leading trkiso [GeV]" , 200 , 0. , 10. );
  hists_1d_["h_subl_abstrkiso"+suf] = fs->make<TH1F>(Form("h_subl_abstrkiso%s",suf.c_str()) , "; Subleading trkiso [GeV]" , 200 , 0. , 10. );
  hists_1d_["h_lead_reltrkiso"+suf] = fs->make<TH1F>(Form("h_lead_reltrkiso%s",suf.c_str()) , "; Leading trkiso / p_{T}" , 200 , 0. , 2. );
  hists_1d_["h_subl_reltrkiso"+suf] = fs->make<TH1F>(Form("h_subl_reltrkiso%s",suf.c_str()) , "; Subleading trkiso / p_{T}" , 200 , 0. , 2. );

  hists_1d_["h_lead_abspfiso"+suf] = fs->make<TH1F>(Form("h_lead_abspfiso%s",suf.c_str()) , "; Leading pfiso [GeV]" , 200 , 0. , 10. );
  hists_1d_["h_subl_abspfiso"+suf] = fs->make<TH1F>(Form("h_subl_abspfiso%s",suf.c_str()) , "; Subleading pfiso [GeV]" , 200 , 0. , 10. );
  hists_1d_["h_lead_relpfiso"+suf] = fs->make<TH1F>(Form("h_lead_relpfiso%s",suf.c_str()) , "; Leading pfiso / p_{T}" , 200 , 0. , 2. );
  hists_1d_["h_subl_relpfiso"+suf] = fs->make<TH1F>(Form("h_subl_relpfiso%s",suf.c_str()) , "; Subleading pfiso / p_{T}" , 200 , 0. , 2. );

  hists_1d_["h_lead_offhlt_dpt"+suf] = fs->make<TH1F>(Form("h_lead_offhlt_dpt%s",suf.c_str()) , "; Leading (p_{T}^{off} - p_{T}^{HLT}) / p_{T}^{off}" , 250 , -5. , 5. );
  hists_1d_["h_subl_offhlt_dpt"+suf] = fs->make<TH1F>(Form("h_subl_offhlt_dpt%s",suf.c_str()) , "; Subleading (p_{T}^{off} - p_{T}^{HLT}) / p_{T}^{off}" , 250 , -5. , 5. );
  hists_1d_["h_lead_offhlt_dr"+suf] = fs->make<TH1F>(Form("h_lead_offhlt_dr%s",suf.c_str()) , "; Leading #DeltaR(off,HLT)" , 600 , 0. , 6. );
  hists_1d_["h_subl_offhlt_dr"+suf] = fs->make<TH1F>(Form("h_subl_offhlt_dr%s",suf.c_str()) , "; Subleading #DeltaR(off,HLT)" , 600 , 0. , 6. );

  hists_2d_["h_lead_abstrkiso_hlt_vs_off"+suf] = fs->make<TH2F>(Form("h_lead_abstrkiso_hlt_vs_off%s",suf.c_str()) , "; Leading offline trkiso [GeV]; Leading HLT trkiso [GeV]" , 50 , 0. , 10. , 50 , 0. , 10. );
  hists_2d_["h_subl_abstrkiso_hlt_vs_off"+suf] = fs->make<TH2F>(Form("h_subl_abstrkiso_hlt_vs_off%s",suf.c_str()) , "; Subleading offline trkiso [GeV]; Subleading HLT trkiso [GeV]" , 50 , 0. , 10. , 50 , 0. , 10. );
  hists_2d_["h_lead_reltrkiso_hlt_vs_off"+suf] = fs->make<TH2F>(Form("h_lead_reltrkiso_hlt_vs_off%s",suf.c_str()) , "; Leading offline trkiso / p_{T}; Leading HLT trkiso / p_{T}" , 50 , 0. , 2. , 50 , 0. , 2. );
  hists_2d_["h_subl_reltrkiso_hlt_vs_off"+suf] = fs->make<TH2F>(Form("h_subl_reltrkiso_hlt_vs_off%s",suf.c_str()) , "; Subleading offline trkiso / p_{T}; Subleading HLT trkiso / p_{T}" , 50 , 0. , 2. , 50 , 0. , 2. );


  return;
}

//____________________________________________________________________________
void DilepTrigAnalyzerRECO::fillHists(const StudyLepton& lead, const StudyLepton& subl, const std::string& suffix, bool isHLT) {

  // if (lead.lv.pt() <= 0.) {
  //   std::cout << "DilepTrigAnalyzerRECO::fillHists: invalid lead pt: " << lead.lv.pt() << std::endl;
  //   return;
  // }

  //  std::cout << "DilepTrigAnalyzerRECO::fillHists: called with suffix: " << suffix << ", isHLT: " << isHLT << std::endl;

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  std::string hlt_suf("");
  if (isHLT) hlt_suf = "_hlt";

  if (lead.lv.pt() > 0) {
    hists_1d_["h_lead_pt"+suf+hlt_suf]->Fill(lead.lv.pt());
    hists_1d_["h_lead_eta"+suf+hlt_suf]->Fill(lead.lv.eta());
    if (lead.type == 13) {
      hists_1d_["h_lead_abstrkiso"+suf+hlt_suf]->Fill(lead.trkiso);
      hists_1d_["h_lead_reltrkiso"+suf+hlt_suf]->Fill(lead.trkiso/lead.lv.pt());
      if (!isHLT) {
	hists_1d_["h_lead_abspfiso"+suf+hlt_suf]->Fill(lead.pfiso);
	hists_1d_["h_lead_relpfiso"+suf+hlt_suf]->Fill(lead.pfiso/lead.lv.pt());
      }
    }
  } // valid lead pt

  if (subl.lv.pt() > 0) {
    hists_1d_["h_subl_pt"+suf+hlt_suf]->Fill(subl.lv.pt());
    hists_1d_["h_subl_eta"+suf+hlt_suf]->Fill(subl.lv.eta());

    if (lead.lv.pt() > 0) {
      LorentzVector dilep = lead.lv+subl.lv;
      hists_1d_["h_mll"+suf+hlt_suf]->Fill(dilep.M());
      float dr = ROOT::Math::VectorUtil::DeltaR(lead.lv,subl.lv);
      hists_1d_["h_dr"+suf+hlt_suf]->Fill(dr);
      if (!isHLT && (lead.charge == subl.charge)) {
	hists_1d_["h_lead_pt_ss"+suf+hlt_suf]->Fill(lead.lv.pt());
	hists_1d_["h_lead_eta_ss"+suf+hlt_suf]->Fill(lead.lv.eta());
	hists_1d_["h_subl_pt_ss"+suf+hlt_suf]->Fill(subl.lv.pt());
	hists_1d_["h_subl_eta_ss"+suf+hlt_suf]->Fill(subl.lv.eta());
	hists_1d_["h_mll_ss"+suf+hlt_suf]->Fill(dilep.M());
	hists_1d_["h_dr_ss"+suf+hlt_suf]->Fill(dr);
      } else if (!isHLT) {
	hists_1d_["h_lead_pt_os"+suf+hlt_suf]->Fill(lead.lv.pt());
	hists_1d_["h_lead_eta_os"+suf+hlt_suf]->Fill(lead.lv.eta());
	hists_1d_["h_subl_pt_os"+suf+hlt_suf]->Fill(subl.lv.pt());
	hists_1d_["h_subl_eta_os"+suf+hlt_suf]->Fill(subl.lv.eta());
	hists_1d_["h_mll_os"+suf+hlt_suf]->Fill(dilep.M());
	hists_1d_["h_dr_os"+suf+hlt_suf]->Fill(dr);
      }
    }
    if (subl.type == 13) {
      hists_1d_["h_subl_abstrkiso"+suf+hlt_suf]->Fill(subl.trkiso);
      hists_1d_["h_subl_reltrkiso"+suf+hlt_suf]->Fill(subl.trkiso/subl.lv.pt());
      if (!isHLT) {
	hists_1d_["h_subl_abspfiso"+suf+hlt_suf]->Fill(subl.pfiso);
	hists_1d_["h_subl_relpfiso"+suf+hlt_suf]->Fill(subl.pfiso/subl.lv.pt());
      }
    }
  } // valid subl pt

  return;
}

//____________________________________________________________________________
void DilepTrigAnalyzerRECO::fillHistsRecoHLT(const StudyLepton& off_lead, const StudyLepton& off_subl, const StudyLepton& hlt_lead, const StudyLepton& hlt_subl, const std::string& suffix) {

  // if (off_lead.pt() <= 0.) {
  //   std::cout << "DilepTrigAnalyzerRECO::fillHists: invalid lead pt: " << off_lead.pt() << std::endl;
  //   return;
  // }

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  if (off_lead.lv.pt() > 0.) {
    bool match_lead = false;
    if (hlt_lead.lv.pt() > 0. && (off_lead.type == hlt_lead.type)) {
      float dr_hlt_lead = ROOT::Math::VectorUtil::DeltaR(off_lead.lv,hlt_lead.lv);
      if (dr_hlt_lead < 0.2) {
	hists_1d_["h_lead_offhlt_dpt"+suf]->Fill( (off_lead.lv.pt() - hlt_lead.lv.pt()) / off_lead.lv.pt() );
	hists_1d_["h_lead_offhlt_dr"+suf]->Fill(dr_hlt_lead);
	if (off_lead.type == 13) {
	  hists_2d_["h_lead_abstrkiso_hlt_vs_off"+suf]->Fill(off_lead.trkiso,hlt_lead.trkiso);
	  hists_2d_["h_lead_reltrkiso_hlt_vs_off"+suf]->Fill(off_lead.trkiso/off_lead.lv.pt(),hlt_lead.trkiso/hlt_lead.lv.pt());
	}
	match_lead = true;
      }
    } // valid hlt_lead
    // if we didn't match the lead, and the leptons are the same type, then check subleading
    else if (!match_lead && hlt_subl.lv.pt() > 0. && (off_lead.type == hlt_subl.type)) {
      float dr_hlt_subl = ROOT::Math::VectorUtil::DeltaR(off_lead.lv,hlt_subl.lv);
      if (dr_hlt_subl < 0.2) {
	hists_1d_["h_lead_offhlt_dpt"+suf]->Fill( (off_lead.lv.pt() - hlt_subl.lv.pt()) / off_lead.lv.pt() );
	hists_1d_["h_lead_offhlt_dr"+suf]->Fill(dr_hlt_subl);
	if (off_lead.type == 13) {
	  hists_2d_["h_lead_abstrkiso_hlt_vs_off"+suf]->Fill(off_lead.trkiso,hlt_subl.trkiso);
	  hists_2d_["h_lead_reltrkiso_hlt_vs_off"+suf]->Fill(off_lead.trkiso/off_lead.lv.pt(),hlt_subl.trkiso/hlt_subl.lv.pt());
	}
      }
    } // hlt_subl
  } // valid off_lead

  if (off_subl.lv.pt() > 0.) {
    bool match_lead = false;
    // if this isn't an emu trigger, check the leading lepton at the HLT
    if (hlt_lead.lv.pt() > 0. && (off_subl.type == hlt_lead.type)) {
      float dr_hlt_lead = ROOT::Math::VectorUtil::DeltaR(off_subl.lv,hlt_lead.lv);
      if (dr_hlt_lead < 0.2) {
	hists_1d_["h_subl_offhlt_dpt"+suf]->Fill( (off_subl.lv.pt() - hlt_lead.lv.pt()) / off_subl.lv.pt() );
	hists_1d_["h_subl_offhlt_dr"+suf]->Fill(dr_hlt_lead);
	if (off_subl.type == 13) {
	  hists_2d_["h_subl_abstrkiso_hlt_vs_off"+suf]->Fill(off_subl.trkiso,hlt_lead.trkiso);
	  hists_2d_["h_subl_reltrkiso_hlt_vs_off"+suf]->Fill(off_subl.trkiso/off_subl.lv.pt(),hlt_lead.trkiso/hlt_lead.lv.pt());
	}
	match_lead = true;
      }
    } // valid hlt_lead
    else if (!match_lead && hlt_subl.lv.pt() > 0. && (off_subl.type == hlt_subl.type)) {
      float dr_hlt_subl = ROOT::Math::VectorUtil::DeltaR(off_subl.lv,hlt_subl.lv);
      if (dr_hlt_subl < 0.2) {
	hists_1d_["h_subl_offhlt_dpt"+suf]->Fill( (off_subl.lv.pt() - hlt_subl.lv.pt()) / off_subl.lv.pt() );
	hists_1d_["h_subl_offhlt_dr"+suf]->Fill(dr_hlt_subl);
	if (off_subl.type == 13) {
	  hists_2d_["h_subl_abstrkiso_hlt_vs_off"+suf]->Fill(off_subl.trkiso,hlt_subl.trkiso);
	  hists_2d_["h_subl_reltrkiso_hlt_vs_off"+suf]->Fill(off_subl.trkiso/off_subl.lv.pt(),hlt_subl.trkiso/hlt_subl.lv.pt());
	}
      }
    } // hlt_subl
  } // valid off_subl

  return;
}

//____________________________________________________________________________
float DilepTrigAnalyzerRECO::muonPFiso(const reco::Muon& muon) {

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
int DilepTrigAnalyzerRECO::chargedHadronVertex( const reco::VertexCollection& vertices, const reco::PFCandidate& pfcand ) const {
  if (pfcand.trackRef().isNull()) return -2;
  return trackVertex(vertices,pfcand.trackRef());
}

//____________________________________________________________________________
// returns best match vertex.  Taken from PFPileUpAlgo, chargedHadronVertex
int DilepTrigAnalyzerRECO::trackVertex( const reco::VertexCollection& vertices, const reco::TrackRef& trackRef ) const {

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


