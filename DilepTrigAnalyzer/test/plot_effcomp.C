#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TString.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TEfficiency.h"
#include "TStyle.h"

void plot_effcomp () {
  
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);

  TH1::SetDefaultSumw2();

  // TString suffix = "8e33_PU28_trigcomp";

  // TFile* f1 = new TFile("histos_MuMuDZ_8e33_PU28_1file.root");
  // TFile* f2 = new TFile("histos_MuTkDZ_8e33_PU28_1file.root");
  // TFile* f3 = new TFile("histos_TkTkDZ_8e33_PU28_1file.root");

  TString suffix = "HIPfix_10p5e33_PU35_trigcomp";

  TFile* f1 = new TFile("histos_oldtkmu_v2_notrigmatch_MuMuDZ_10p5e33_PU35_1file.root");
  TFile* f2 = new TFile("histos_oldtkmu_v2_notrigmatch_MuTkDZ_10p5e33_PU35_1file.root");
  TFile* f3 = new TFile("histos_oldtkmu_v2_notrigmatch_iso_TkTkDZ10_10p5e33_PU35_1file.root");

  // TString suffix = "TkTkDZ10_match_lumicomp";

  // TFile* f1 = new TFile("histos_v16_newgt_reHLT_match_iso_TkTkDZ10_11e33_PU38_1file.root");
  // TFile* f2 = new TFile("histos_v16_newgt_reHLT_match_iso_TkTkDZ10_8e33_PU28_1file.root");
  // TFile* f3 = new TFile("histos_v16_newgt_reHLT_match_iso_TkTkDZ10_2e33_PU17_1file.root");

  TH1F* h_denom1 = (TH1F*) ((TH1F*) f1->Get("diMuTrigAnalyzerRECO/h_nvtx_denom"))->Clone("h_denom1");
  TH1F* h_num1 = (TH1F*) ((TH1F*) f1->Get("diMuTrigAnalyzerRECO/h_nvtx_num"))->Clone("h_num1");
  TH1F* h_denom2 = (TH1F*) ((TH1F*) f2->Get("diMuTrigAnalyzerRECO/h_nvtx_denom"))->Clone("h_denom2");
  TH1F* h_num2 = (TH1F*) ((TH1F*) f2->Get("diMuTrigAnalyzerRECO/h_nvtx_num"))->Clone("h_num2");
  TH1F* h_denom3 = (TH1F*) ((TH1F*) f3->Get("diMuTrigAnalyzerRECO/h_nvtx_denom"))->Clone("h_denom3");
  TH1F* h_num3 = (TH1F*) ((TH1F*) f3->Get("diMuTrigAnalyzerRECO/h_nvtx_num"))->Clone("h_num3");

  TEfficiency* h_eff1 = new TEfficiency(*h_num1, *h_denom1);
  h_eff1->SetName("h_nvtx_eff1");
  h_eff1->SetLineColor(kRed);
  h_eff1->SetMarkerColor(kRed);

  TEfficiency* h_eff2 = new TEfficiency(*h_num2, *h_denom2);
  h_eff2->SetName("h_nvtx_eff2");
  h_eff2->SetLineColor(kGreen+2);
  h_eff2->SetMarkerColor(kGreen+2);

  TEfficiency* h_eff3 = new TEfficiency(*h_num3, *h_denom3);
  h_eff3->SetName("h_nvtx_eff3");
  h_eff3->SetLineColor(kBlue);
  h_eff3->SetMarkerColor(kBlue);

  TCanvas* c = new TCanvas("c","c");
  c->SetGrid(1,1);
  c->cd();

  //  TH2F* h_axis = new TH2F("h_axis",";N(vtx);DZ Efficiency",10,0,50,22,0,1.05);
  TH2F* h_axis = new TH2F("h_axis",";N(vtx);DZ Efficiency",10,0,50,22,0.6,1.05);
  h_axis->GetYaxis()->SetTitleOffset(0.98);
  h_axis->Draw();

  h_eff3->Draw("pe same");
  h_eff2->Draw("pe same");
  h_eff1->Draw("pe same");
  
  //  TLegend* leg = new TLegend(0.2,0.16,0.63,0.47);
  TLegend* leg = new TLegend(0.2,0.16,0.84,0.36);
  
  // //  leg->SetHeader("Run 276315, 8e33, PU28");
  // leg->SetHeader("Run 278017, 1.1e34, PU38");
  leg->SetHeader("Run 278803, 1.05e34, PU35, HIP fix");
  leg->AddEntry(h_eff3,"TkMu17_TkMu8_DZ1p0","pe");
  leg->AddEntry(h_eff2,"Mu17_TkMu8_DZ","pe");
  leg->AddEntry(h_eff1,"Mu17_Mu8_DZ","pe");

  //leg->SetHeader("Mu17_Mu8_DZ");
  //leg->SetHeader("Mu17_TkMu8_DZ");
  //  leg->SetHeader("TkMu17_TkMu8_DZ");
  // leg->SetHeader("TkMu17_TkMu8_DZ10");
  // leg->AddEntry(h_eff3,"Run 273503, 2e33, PU17","pe");
  // leg->AddEntry(h_eff2,"Run 276315, 8e33, PU28","pe");
  // leg->AddEntry(h_eff1,"Run 278017, 1.1e34, PU38","pe");

  leg->Draw("same");

  c->SaveAs(Form("effcomp_%s.pdf",suffix.Data()));
  c->SaveAs(Form("effcomp_%s.eps",suffix.Data()));
  c->SaveAs(Form("effcomp_%s.png",suffix.Data()));

  // f1->Close();
  // f2->Close();
  // f3->Close();
  
  return;
}
