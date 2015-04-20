/** \class SingleEleTrigAnalyzerAOD
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
#include "TrigAnalyzer/DilepTrigAnalyzer/interface/SingleEleTrigAnalyzerAOD.h"

#include "DataFormats/HLTReco/interface/TriggerObject.h"

// ROOT includes
#include "Math/VectorUtil.h"

#include <cassert>

using namespace reco;
using namespace edm;

//
// constructors and destructor
//
//____________________________________________________________________________
SingleEleTrigAnalyzerAOD::SingleEleTrigAnalyzerAOD(const edm::ParameterSet& ps) : 
  processName_(ps.getParameter<std::string>("processName")),
  triggerName_(ps.getParameter<std::string>("triggerName")),
  triggerResultsTag_(ps.getParameter<edm::InputTag>("triggerResults")),
  triggerEventTag_(ps.getParameter<edm::InputTag>("triggerEvent")),
  genParticlesTag_(ps.getParameter<edm::InputTag>("genParticles")),
  genPt_(ps.getParameter<double>("genPt")),
  genEta_(ps.getParameter<double>("genEta")),
  verbose_(ps.getParameter<bool>("verbose"))
{
  using namespace std;
  using namespace edm;

  cout << "SingleEleTrigAnalyzerAOD configuration: " << endl
       << "   ProcessName = " << processName_ << endl
       << "   triggerName = " << triggerName_ << endl
       << "   TriggerResultsTag = " << triggerResultsTag_.encode() << endl
       << "   TriggerEventTag = " << triggerEventTag_.encode() << endl
       << "   GenParticlesTag = " << genParticlesTag_.encode() << endl
       << "   GenPt = " << genPt_ << endl
       << "   GenEta = " << genEta_ << endl
       << "   Verbose = " << verbose_ << endl;

  // histogram setup
  edm::Service<TFileService> fs;
  bookHistsGen(fs,"gen");
  bookHistsGen(fs,"gen_match");
  bookHistsGen(fs,"gen_nomatch");

}

//____________________________________________________________________________
SingleEleTrigAnalyzerAOD::~SingleEleTrigAnalyzerAOD()
{
}

//
// member functions
//
//____________________________________________________________________________
void
SingleEleTrigAnalyzerAOD::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
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
	cout << "SingleEleTrigAnalyzerAOD::analyze:"
	     << " TriggerName " << triggerName_ 
	     << " not available in (new) config!" << endl;
      }
    } // if changed
  } else {
    cout << "SingleEleTrigAnalyzerAOD::analyze:"
	 << " config extraction failure with process name "
	 << processName_ << endl;
  }

}

//____________________________________________________________________________
// ------------ method called to produce the data  ------------
void
SingleEleTrigAnalyzerAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;

  if (verbose_) cout << endl;

  // get event products
  iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
    cout << "SingleEleTrigAnalyzerAOD::analyze: Error in getting TriggerResults product from Event!" << endl;
    return;
  }
  iEvent.getByLabel(triggerEventTag_,triggerEventHandle_);
  if (!triggerEventHandle_.isValid()) {
    cout << "SingleEleTrigAnalyzerAOD::analyze: Error in getting TriggerEvent product from Event!" << endl;
    return;
  }

  // sanity check
  assert(triggerResultsHandle_->size()==hltConfig_.size());

  iEvent.getByLabel(genParticlesTag_, genParticlesHandle_);

  analyzeTrigger(iEvent, iSetup, triggerName_);

  if (verbose_) cout << endl;

  return;
}

//____________________________________________________________________________

bool SingleEleTrigAnalyzerAOD::analyzeTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& triggerName) {
  
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  if (verbose_) cout << endl;

  const unsigned int ntrigs(hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));

  // abort on invalid trigger name
  if (triggerIndex>=ntrigs) {
    cout << "SingleEleTrigAnalyzerAOD::analyzeTrigger: path "
	 << triggerName << " - not found!" << endl;
    return false;
  }

  if (verbose_) {
    cout << "SingleEleTrigAnalyzerAOD::analyzeTrigger: path "
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

  //  if (!wasRun || !accept || error) return false;
  if (!wasRun || error) return false;

  // loop over trigger and reco objects, match, make plots

  // first, get trigger objects from last filter

  //------------------------------------
  //  hlt objects
  //------------------------------------

  bool foundElectrons = false;
  std::vector<LorentzVector> trigElectrons;

    // --------- code for triggerEvent, instead of triggerEventWithRefs -------------------

  for (unsigned int j=moduleIndex; j!=0; --j) {
    const string& moduleLabel(moduleLabels[j]);
    const string  moduleType(hltConfig_.moduleType(moduleLabel));
    // check whether the module is packed up in TriggerEvent product
    const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(moduleLabel,"",processName_)));
    if (filterIndex>=triggerEventHandle_->sizeFilters()) continue;
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

      // electrons
      if (!foundElectrons && (VIDS[i] == 81)) {
	trigElectrons.push_back(lv);
      } // hlt electrons

      if (verbose_) {
    	cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
    	     << TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()
    	     << endl;
      }
    } // loop on trig objects

    if (trigElectrons.size() > 0) {
      foundElectrons = true;
    }

    if (foundElectrons) break;
  } // backwards loop on modules

  if (accept && trigElectrons.size() == 0) {
    cout << "SingleEleTrigAnalyzerAOD::analyzeTrigger: WARNING!! no valid trigger leptons!" << endl;
  }

  //-------------------------------------
  //   compare to gen level
  //-------------------------------------

  const float dr_trigmatch = 0.2;

  GenParticleCollection::const_iterator genps_end = genParticlesHandle_->end();  // Iterator
  for ( GenParticleCollection::const_iterator genp = genParticlesHandle_->begin(); genp != genps_end; ++genp ) {
    // can use status 3 for pythia 6
    //      if (genp->status() != 3) continue;
    // allow status 23 for pythia 8
    if ((genp->status() != 1) && (genp->status() != 23)) continue;
    if (abs(genp->pdgId()) != 11) continue;
    // W mother
    if (abs(genp->mother()->pdgId()) != 24) continue;

    // pt, eta acceptance cuts
    if ((genp->pt() > genPt_) && (fabs(genp->eta()) < genEta_)) {

      LorentzVector lv(genp->p4());

      bool trigmatch = false;
      if (accept) {
	for (unsigned int itrig=0; itrig < trigElectrons.size(); ++itrig) {
	  if (ROOT::Math::VectorUtil::DeltaR(lv,trigElectrons.at(itrig)) < dr_trigmatch) trigmatch = true;
	}
      }

      fillHistsGen(lv,"gen");
      if (trigmatch) {
	fillHistsGen(lv,"gen_match");
      } else {
	fillHistsGen(lv,"gen_nomatch");
      }

    } // in acceptance

  } // loop over genps

  return accept;
}

//____________________________________________________________________________
void SingleEleTrigAnalyzerAOD::bookHistsGen(edm::Service<TFileService>& fs, const std::string& suffix) {

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  hists_1d_["h_pt"+suf] = fs->make<TH1F>(Form("h_pt%s",suf.c_str()) , "; p_{T} [GeV]" , 100 , 0. , 100. );
  hists_1d_["h_eta"+suf] = fs->make<TH1F>(Form("h_eta%s",suf.c_str()) , "; #eta" , 100 , -3. , 3. );

  return;
}

//____________________________________________________________________________
void SingleEleTrigAnalyzerAOD::fillHistsGen(const LorentzVector& lv, const std::string& suffix) {

  //  std::cout << "SingleEleTrigAnalyzerAOD::fillHists: called with suffix: " << suffix << std::endl;

  std::string suf(suffix);
  if (suffix.size()) suf = "_"+suffix;

  hists_1d_["h_pt"+suf]->Fill(lv.pt());
  hists_1d_["h_eta"+suf]->Fill(lv.eta());

  return;
}
