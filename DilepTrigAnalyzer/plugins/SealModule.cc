#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "TrigAnalyzer/DilepTrigAnalyzer/interface/DilepTrigAnalyzerRAW.h"
#include "TrigAnalyzer/DilepTrigAnalyzer/interface/DilepTrigAnalyzerRECO.h"
#include "TrigAnalyzer/DilepTrigAnalyzer/interface/SingleMuTrigAnalyzerRECO.h"
#include "TrigAnalyzer/DilepTrigAnalyzer/interface/TrigAnalyzerRECORef.h"
#include "TrigAnalyzer/DilepTrigAnalyzer/interface/TrackCompAnalyzer.h"

DEFINE_FWK_MODULE(DilepTrigAnalyzerRAW);
DEFINE_FWK_MODULE(DilepTrigAnalyzerRECO);
DEFINE_FWK_MODULE(SingleMuTrigAnalyzerRECO);
DEFINE_FWK_MODULE(TrigAnalyzerRECORef);
DEFINE_FWK_MODULE(TrackCompAnalyzer);
