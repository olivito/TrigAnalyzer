#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TString.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TStyle.h"

//const TString& infile = "histos_noiter2_v20_L1denom_TkMu24_2e33_PU17_1file.root";
//const TString& infile = "histos_2losthits_v4_L1denom_TkMu24_2e33_PU17_1file.root";
//const TString& infile = "histos_noiter2_v20_L1denom_TkMu24_11e33_PU38_1file.root";
//const TString& infile = "histos_2losthits_v4_L1denom_TkMu24_11e33_PU38_1file.root";
//const TString& infile = "histos_oldtkmu_v2_L1denom_TkMu24_HIPfix_10p5e33_PU35_1file.root";
const TString& infile = "histos_2losthits_v8_L1denom_TkMu24_HIPfix_10p5e33_PU35_1file.root";

void plot_trigeff_singlemu () {

  //TString label = "L1denom_TkMu24_noiter2_v20_PU17";
  //TString label = "L1denom_TkMu24_2losthits_v4_PU17";
  //TString label = "L1denom_TkMu24_noiter2_v20_PU38";
  //    TString label = "L1denom_TkMu24_2losthits_v4_PU38";
  //  TString label = "L1denom_TkMu24_oldtkmu_v2_HIPfix_PU35";
    TString label = "L1denom_TkMu24_2losthits_v8_HIPfix_PU35";
  
  TString analyzer_name = "singleMuTrigTnPAnalyzerRECO";
  
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  
  TFile* f_in = new TFile(infile);

  // ---- efficiency vs pt
  
  TH1F* h_pt_denom = (TH1F*) f_in->Get(analyzer_name+"/h_pt_probe_all");
  TH1F* h_pt_num = (TH1F*) f_in->Get(analyzer_name+"/h_pt_probe_pass");

  // value for eff, pt > 30
  float n_denom_ptcut = h_pt_denom->Integral(31,-1);
  float n_num_ptcut = h_pt_num->Integral(31,-1);
  float eff_ptcut = n_num_ptcut / n_denom_ptcut;
  float eff_ptcut_UP = TEfficiency::ClopperPearson(n_denom_ptcut, n_num_ptcut, 0.68, 1);
  float eff_ptcut_DN = TEfficiency::ClopperPearson(n_denom_ptcut, n_num_ptcut, 0.68, 0);
  std::cout << "effieciency for pt > 30: " << eff_ptcut
	    << ", + " << eff_ptcut_UP - eff_ptcut << " - " << eff_ptcut - eff_ptcut_DN << std::endl;
  
  h_pt_denom->Rebin(5);
  h_pt_num->Rebin(5);
  
  TCanvas* c_pt = new TCanvas("c_pt","c_pt");
  c_pt->SetGrid(1,1);
  c_pt->cd();

  TH2F* h_pt_axis = new TH2F("h_pt_axis",";p_{T} [GeV];Efficiency of HLT_TkMu24",100,0,100,20,0,1);
  h_pt_axis->GetYaxis()->SetTitleOffset(0.98);
  h_pt_axis->Draw();
  
  TEfficiency* h_pt_eff = new TEfficiency(*h_pt_num, *h_pt_denom);
  h_pt_eff->SetLineColor(kRed);
  h_pt_eff->SetMarkerColor(kRed);
  
  h_pt_eff->Draw("pe same");

  c_pt->SaveAs(Form("eff_pt_%s.pdf",label.Data()));
  
  // ---- efficiency vs eta
  
  TH1F* h_eta_denom = (TH1F*) f_in->Get(analyzer_name+"/h_eta_probe_all");
  TH1F* h_eta_num = (TH1F*) f_in->Get(analyzer_name+"/h_eta_probe_pass");
  h_eta_denom->Rebin(5);
  h_eta_num->Rebin(5);
  
  TCanvas* c_eta = new TCanvas("c_eta","c_eta");
  c_eta->SetGrid(1,1);
  c_eta->cd();

  TH2F* h_eta_axis = new TH2F("h_eta_axis",";#eta;Efficiency of HLT_TkMu24",100,-3.,3.,20,0,1);
  h_eta_axis->GetYaxis()->SetTitleOffset(0.98);
  h_eta_axis->Draw();
  
  TEfficiency* h_eta_eff = new TEfficiency(*h_eta_num, *h_eta_denom);
  h_eta_eff->SetLineColor(kRed);
  h_eta_eff->SetMarkerColor(kRed);
  
  h_eta_eff->Draw("pe same");

  c_eta->SaveAs(Form("eff_eta_%s.pdf",label.Data()));

  // ---- efficiency vs phi
  
  TH1F* h_phi_denom = (TH1F*) f_in->Get(analyzer_name+"/h_phi_probe_all");
  TH1F* h_phi_num = (TH1F*) f_in->Get(analyzer_name+"/h_phi_probe_pass");
  h_phi_denom->Rebin(5);
  h_phi_num->Rebin(5);
  
  TCanvas* c_phi = new TCanvas("c_phi","c_phi");
  c_phi->SetGrid(1,1);
  c_phi->cd();

  TH2F* h_phi_axis = new TH2F("h_phi_axis",";#phi;Efficiency of HLT_TkMu24",100,-3.14,3.14,20,0,1);
  h_phi_axis->GetYaxis()->SetTitleOffset(0.98);
  h_phi_axis->Draw();
  
  TEfficiency* h_phi_eff = new TEfficiency(*h_phi_num, *h_phi_denom);
  h_phi_eff->SetLineColor(kRed);
  h_phi_eff->SetMarkerColor(kRed);
  
  h_phi_eff->Draw("pe same");
 
  c_phi->SaveAs(Form("eff_phi_%s.pdf",label.Data()));

  return;
}
