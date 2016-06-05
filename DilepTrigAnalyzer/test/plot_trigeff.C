#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TString.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TStyle.h"

void plot_trigeff (const TString& infile = "histos_100k_skim.root") {

  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);

  TH1::SetDefaultSumw2();
  
  TFile* f_in = new TFile(infile);

  // ---- efficiency vs pt
  
  TH1F* h_pt_denom = (TH1F*) f_in->Get("singleMuTrigAnalyzerMiniAOD/h_pt_probe_all");
  TH1F* h_pt_num = (TH1F*) f_in->Get("singleMuTrigAnalyzerMiniAOD/h_pt_probe_pass");
  
  TCanvas* c_pt = new TCanvas("c_pt","c_pt");
  c_pt->SetGrid(1,1);
  c_pt->cd();

  TH2F* h_pt_axis = new TH2F("h_pt_axis",";p_{T} [GeV];Efficiency of HLT_IsoMu22",100,0,100,20,0,1);
  h_pt_axis->GetYaxis()->SetTitleOffset(0.98);
  h_pt_axis->Draw();
  
  TEfficiency* h_pt_eff = new TEfficiency(*h_pt_num, *h_pt_denom);
  h_pt_eff->SetLineColor(kRed);
  h_pt_eff->SetMarkerColor(kRed);
  
  h_pt_eff->Draw("pe same");
  
  // ---- efficiency vs eta
  
  TH1F* h_eta_denom = (TH1F*) f_in->Get("singleMuTrigAnalyzerMiniAOD/h_eta_probe_all");
  TH1F* h_eta_num = (TH1F*) f_in->Get("singleMuTrigAnalyzerMiniAOD/h_eta_probe_pass");
  
  TCanvas* c_eta = new TCanvas("c_eta","c_eta");
  c_eta->SetGrid(1,1);
  c_eta->cd();

  TH2F* h_eta_axis = new TH2F("h_eta_axis",";#eta;Efficiency of HLT_IsoMu22",100,-3.,3.,20,0,1);
  h_eta_axis->GetYaxis()->SetTitleOffset(0.98);
  h_eta_axis->Draw();
  
  TEfficiency* h_eta_eff = new TEfficiency(*h_eta_num, *h_eta_denom);
  h_eta_eff->SetLineColor(kRed);
  h_eta_eff->SetMarkerColor(kRed);
  
  h_eta_eff->Draw("pe same");

  // ---- efficiency vs phi
  
  TH1F* h_phi_denom = (TH1F*) f_in->Get("singleMuTrigAnalyzerMiniAOD/h_phi_probe_all");
  TH1F* h_phi_num = (TH1F*) f_in->Get("singleMuTrigAnalyzerMiniAOD/h_phi_probe_pass");
  
  TCanvas* c_phi = new TCanvas("c_phi","c_phi");
  c_phi->SetGrid(1,1);
  c_phi->cd();

  TH2F* h_phi_axis = new TH2F("h_phi_axis",";#phi;Efficiency of HLT_IsoMu22",100,-3.14,3.14,20,0,1);
  h_phi_axis->GetYaxis()->SetTitleOffset(0.98);
  h_phi_axis->Draw();
  
  TEfficiency* h_phi_eff = new TEfficiency(*h_phi_num, *h_phi_denom);
  h_phi_eff->SetLineColor(kRed);
  h_phi_eff->SetMarkerColor(kRed);
  
  h_phi_eff->Draw("pe same");
 
  return;
}
