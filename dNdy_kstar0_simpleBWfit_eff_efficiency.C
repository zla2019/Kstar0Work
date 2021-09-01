#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TProfile.h"
#include "TStyle.h"

#include "../../func.h"

double IMcut[2] = {
    0.75,
    1.05}; // K*0 mass from PDG = 895.94 +- 0.22 [MeV/c2], width 48.7 +- 0.8 MeV
int colorIndex[] = {1, 2, 3, 4, 6, 7, 8};
double sigma = 0.049;
float br = 1.00;

TGraphErrors *convertToGraph(TH1F *h);

void dNdy_kstar0_simpleBWfit_eff_efficiency(TString cutName = "Mix",
                                            int isAntiKstar = 0) {

  bool debug = true;

  // 0-10%, 10-40%, 40-60%, 0-80%

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadTopMargin(0.05);

  TH3F::SetDefaultSumw2();
  TH2F::SetDefaultSumw2();
  TH1F::SetDefaultSumw2();

  TH1F *ax;
  TF1 *zero = new TF1("zero", "[0]", 0, 110);
  zero->SetLineStyle(2);
  zero->SetLineWidth(1);
  zero->SetLineColor(1);

  TLatex lt;
  lt.SetNDC();
  lt.SetTextFont(43);
  lt.SetTextSize(14);

  //==================================================================

  const int nCentFlow = 4;
  TString centName[] = {"0-10%", "10-40%", "40-60%", "0-80%"};
  TString etaName[] = {"-0.9 < y < -0.6", "-0.6 < y < -0.3", "-0.3 < y < 0",
                       "0 < y < 0.3"};
  TString etaNameFow[] = {"-1.0 < y < -0.8"}; // used at begining
  TString etaNameMid[] = {"-0.2 < y < 0.2"};

  // combined rap bin
  TString etaNameTot[] = {"-1.0 < y < -0.8", "-0.9 < y < -0.6",
                          "-0.6 < y < -0.3", "-0.3 < y < 0",
                          "0 < y < 0.3",     "-0.2 < y < 0.2"};
  TString etaNameTotScale[] = {"-1.0<y<-0.8, x2^{-1}", "-0.9<y<-0.6",
                               "-0.6<y<-0.3, x2",      "-0.3<y<0, x2^{2}",
                               "0<y<0.3,    x2^{3}",   "-0.2<y<0.2,  x2^{4}"};

  double fctetaXaxis[] = {-0.9, -0.6, -0.3, 0, 0.3};
  double fctetaXaxisFow[] = {-1.0, -0.8};
  double fctetaXaxisMid[] = {-0.2, 0.2};

  float rapiditywidth = 0.3;
  float rapiditywidthFow = 0.2;
  float rapiditywidthMid = 0.4;

  //-1.0 < y < -0.8
  const int nPtFow = 3;
  // double fctpTFowXaxis[nPtFow+1] = {0.4, 0.7, 1.0, 1.5, 2.0};
  double fctpTFowXaxis[nPtFow + 1] = {0.4, 0.9, 1.4, 2.0};

  //-0.8 < y < 0.2
  const int nEta = 4;
  const int nPt = 4;
  double fctpTXaxis[nEta][nPt + 1] = {
      //{0.4, 0.6, 0.9, 1.2, 1.5, 2.0},
      {0.4, 0.8, 1.2, 1.6, 2.0},
      {0.4, 0.8, 1.2, 1.6, 2.0},
      {0.4, 0.8, 1.2, 1.6, 2.0},
      {0.4, 0.8, 1.2, 1.6, 2.0},
      //        {0.4, 0.8, 1.2, 1.6, 2.0} // 0, 0.2
  };

  const int nPtMid = 4;
  double fctpTMidXaxis[nPtMid + 1] = {0.4, 0.8, 1.2, 1.6, 2.0};
  // double fctpTMidXaxis[nPtMid+1] = {0.4, 0.7, 1.0, 1.3, 1.6, 2.0};

  // nTotEta begin to include all !!!
  const int nTotEta = 6;
  // float scaler4Show[nTotEta] = {0.5, 1., 2., 4., 8., 16., 32}; //scale for
  // clarity
  float scaler4Show[nTotEta] = {1., 1., 1., 1., 1., 1}; // scale for clarity
  const Double_t RapiditypTLow[nTotEta][nCentFlow] = {
      {0.4, 0.4, 0.4}, //-1.0, -0.8
      {0.4, 0.4, 0.4}, //-0.9, -0.6
      {0.4, 0.4, 0.4}, //-0.6, -0.3
      {0.4, 0.4, 0.4}, //-0.3, -0.0
      {0.4, 0.4, 0.4}, // 0, 0.3
      {0.4, 0.4, 0.4}  // -0.2, 0.2
  };
  const Double_t RapiditypTHigh[nTotEta][nCentFlow] = {
      {2.0, 2.0, 2.0}, {2.0, 2.0, 2.0}, {2.0, 2.0, 2.0},
      {2.0, 2.0, 2.0}, {2.0, 2.0, 2.0}, {2.0, 2.0, 2.0},
  };
  //==================================================================
  // centrality
  //    TFile *pol_file =
  //    TFile::Open("kstar0_3gev_TPCandTOF_spinYmp5to0_flowYmp8to0_ana_hist.root");
  //    TFile *pol_file =
  //    TFile::Open("19154032_5000010_kstar0Tree_TPCorTOF_hist_Aug12.root");
  //    TFile *pol_file =
  //    TFile::Open("1.6MMuDataTree_New_TPCorTOF_hist_Aug12.root"); TFile
  //    *pol_file = TFile::Open("Kstar0_20210823_embedding.root");
  TFile *pol_file =
      TFile::Open("../../Kstar0_20210824_RealData_TPCandTOF.root");
  TFile *efficiency = TFile::Open("../../KstarEfficiency.root");
  TFile *efficiency_bigRange =
      TFile::Open("../../KstarEfficiency_binning2.root");
  TFile *efficiency_fow_mid = TFile::Open("../../KstarEfficiency_fow_mid.root");
  TH2F *Efficiency[4];
  TH2F *EfficiencyBigRange[4];
  TH2F *EfficiencyFow[4];
  TH2F *EfficiencyMid[4];
  for (int icent = 0; icent < 4; ++icent) {
    Efficiency[icent] = new TH2F("", "", 20, -1, 1, 100, 0, 5);
    TH2F *htmp = (TH2F *)efficiency->Get(Form("Efficiency_cent%i", icent));
    EfficiencyBigRange[icent] =
        (TH2F *)efficiency_bigRange->Get(Form("Efficiency_cent%i", icent));
    EfficiencyFow[icent] =
        (TH2F *)efficiency_fow_mid->Get(Form("Efficiency_fow_cent%i", icent));
    EfficiencyMid[icent] =
        (TH2F *)efficiency_fow_mid->Get(Form("Efficiency_mid_cent%i", icent));
    for (int ix = 6; ix < 26; ++ix) {
      for (int iy = 0; iy < Efficiency[icent]->GetYaxis()->GetNbins(); ++iy) {
        Efficiency[icent]->SetBinContent(ix - 5, iy + 1,
                                         htmp->GetBinContent(ix, iy + 1));
      }
    }
  }
  TH1F *hCentrality_event = (TH1F *)pol_file->Get("hCentrality");
  int nEvents = hCentrality_event->GetEntries();
  int nEvents_cent[nCentFlow];
  nEvents_cent[0] = hCentrality_event->GetBinContent(9) +
                    hCentrality_event->GetBinContent(8); // 0-10%
  nEvents_cent[1] = hCentrality_event->GetBinContent(7) +
                    hCentrality_event->GetBinContent(6) +
                    hCentrality_event->GetBinContent(5); // 10-40%
  nEvents_cent[2] = hCentrality_event->GetBinContent(4) +
                    hCentrality_event->GetBinContent(3); // 40-60%
  nEvents_cent[3] = hCentrality_event->GetBinContent(2) +
                    hCentrality_event->GetBinContent(1) + nEvents_cent[2] +
                    nEvents_cent[1] + nEvents_cent[0]; // 0-80%
  cout << "Number of events: " << nEvents_cent[0] << " (0-10%), "
       << nEvents_cent[1] << " (10-40%), " << nEvents_cent[2] << " (40-60%), "
       << nEvents_cent[3] << " (0-80%)." << endl;

  // TFile *pol_file   = TFile::Open(Form("macro/output/sys_%s_sig.root",
  // cutName.Data())); TFile *polBg_file =
  // TFile::Open(Form("macro/output/sys_%s_bkg.root", cutName.Data()));
  /*
  // efficiency
  TH1F *hEff_accReg_we[nEta][nCentFlow];
  TH1F *hEff_accReg_we_fow[nCentFlow];
  TH1F *hEff_accReg_we_mid[nCentFlow];
  auto feff = new TFile(Form("macro/output/sys_eff_%s.root", cutName.Data()),
  "r"); for(int icent=0; icent<nCentFlow; icent++){ for(int ieta=0; ieta<nEta;
  ieta++) hEff_accReg_we[ieta][icent] =
  (TH1F*)feff->Get(Form("hEff_accReg_we_%d_%d", ieta, icent));
      hEff_accReg_we_fow[icent] = (TH1F*)feff->Get(Form("hEff_accReg_we_fow_%d",
  icent)); hEff_accReg_we_mid[icent] =
  (TH1F*)feff->Get(Form("hEff_accReg_we_mid_%d", icent));
  }
   */

  TString histName = "NumKstar0InvMassvsPtY";
  TString particleName = "Kstar0";

  if (isAntiKstar == 1) {
    histName = "NumAntiKstar0InvMassvsPtY";
    particleName = "antiKstar0";
  }
  // if(ksigma == 4) cutName = "nsigma_2";

  system(Form("mkdir -p ../../files/plots_%s", cutName.Data()));
  //==========================================================================================================================================
  // get historgram
  TH3F *h3InvMassXiYvsPt_tot[nCentFlow]; // 0-10%, 10-40%, 40-60%, 0-80%
  TH3F *h3InvMassXiYvsPt_bg[nCentFlow];

  // eta bins: -0.8， 0.2
  TH1F *MinvLam_EPflow_tot[nCentFlow][nEta][nPt];
  TH1F *MinvLam_EPflow_bg[nCentFlow][nEta][nPt];

  // eta bin: -1.0， -0.8
  TH1F *MinvLamFow_EPflow_tot[nCentFlow][nPtFow];
  TH1F *MinvLamFow_EPflow_bg[nCentFlow][nPtFow];

  // eta bin: -0.5， 0
  TH1F *MinvLamMid_EPflow_tot[nCentFlow][nPtMid];
  TH1F *MinvLamMid_EPflow_bg[nCentFlow][nPtMid];

  // 0-10%
  h3InvMassXiYvsPt_tot[0] =
      (TH3F *)pol_file->Get(Form("%s_Cent0", histName.Data())); // 0-5%
  h3InvMassXiYvsPt_tot[0]->Add(
      (TH3F *)pol_file->Get(Form("%s_Cent1", histName.Data()))); // 5-10%

  h3InvMassXiYvsPt_bg[0] = (TH3F *)pol_file->Get(
      Form("%s_%s_Cent0", histName.Data(), cutName.Data())); // 0-5%
  h3InvMassXiYvsPt_bg[0]->Add((TH3F *)pol_file->Get(
      Form("%s_%s_Cent1", histName.Data(), cutName.Data()))); // 5-10%

  // 10-40%
  h3InvMassXiYvsPt_tot[1] =
      (TH3F *)pol_file->Get(Form("%s_Cent2", histName.Data())); // 10-20%
  h3InvMassXiYvsPt_tot[1]->Add(
      (TH3F *)pol_file->Get(Form("%s_Cent3", histName.Data()))); // 20-30%
  h3InvMassXiYvsPt_tot[1]->Add(
      (TH3F *)pol_file->Get(Form("%s_Cent4", histName.Data()))); // 30-40%

  h3InvMassXiYvsPt_bg[1] = (TH3F *)pol_file->Get(
      Form("%s_%s_Cent2", histName.Data(), cutName.Data())); // 10-20%
  h3InvMassXiYvsPt_bg[1]->Add((TH3F *)pol_file->Get(
      Form("%s_%s_Cent3", histName.Data(), cutName.Data()))); // 20-30%
  h3InvMassXiYvsPt_bg[1]->Add((TH3F *)pol_file->Get(
      Form("%s_%s_Cent4", histName.Data(), cutName.Data()))); // 30-40%

  // 40-60%
  h3InvMassXiYvsPt_tot[2] =
      (TH3F *)pol_file->Get(Form("%s_Cent5", histName.Data())); // 40-50%
  h3InvMassXiYvsPt_tot[2]->Add(
      (TH3F *)pol_file->Get(Form("%s_Cent6", histName.Data()))); // 50-60%

  h3InvMassXiYvsPt_bg[2] = (TH3F *)pol_file->Get(
      Form("%s_%s_Cent5", histName.Data(), cutName.Data())); // 40-50%
  h3InvMassXiYvsPt_bg[2]->Add((TH3F *)pol_file->Get(
      Form("%s_%s_Cent6", histName.Data(), cutName.Data()))); // 50-60%

  // 0-80%
  h3InvMassXiYvsPt_tot[3] =
      (TH3F *)pol_file->Get(Form("%s_Cent7", histName.Data())); // 60-70%
  h3InvMassXiYvsPt_tot[3]->Add(
      (TH3F *)pol_file->Get(Form("%s_Cent8", histName.Data()))); // 70-80%
  h3InvMassXiYvsPt_tot[3]->Add(h3InvMassXiYvsPt_tot[0]);         // 0-10%
  h3InvMassXiYvsPt_tot[3]->Add(h3InvMassXiYvsPt_tot[1]);         // 10-40%
  h3InvMassXiYvsPt_tot[3]->Add(h3InvMassXiYvsPt_tot[2]);         // 40-60%

  h3InvMassXiYvsPt_bg[3] = (TH3F *)pol_file->Get(
      Form("%s_%s_Cent7", histName.Data(), cutName.Data())); // 60-70%
  h3InvMassXiYvsPt_bg[3]->Add((TH3F *)pol_file->Get(
      Form("%s_%s_Cent8", histName.Data(), cutName.Data()))); // 70-80%
  h3InvMassXiYvsPt_bg[3]->Add(h3InvMassXiYvsPt_bg[0]);        // 0-10%
  h3InvMassXiYvsPt_bg[3]->Add(h3InvMassXiYvsPt_bg[1]);        // 10-40%
  h3InvMassXiYvsPt_bg[3]->Add(h3InvMassXiYvsPt_bg[2]);        // 40-60%

  for (int icent = 0; icent < nCentFlow; icent++) {
    // rapidity bins: (-0.8， 0.2)
    for (int ieta = 0; ieta < nEta; ieta++) {
      int eta_low =
          h3InvMassXiYvsPt_tot[icent]->GetYaxis()->FindBin(fctetaXaxis[ieta]);
      int eta_high = h3InvMassXiYvsPt_tot[icent]->GetYaxis()->FindBin(
                         fctetaXaxis[ieta + 1]) -
                     1;
      if (ieta == 0)
        eta_low = 3;
      cout << "Rapidity bins (-0.8， 0.2): " << ieta << " -- " << eta_low
           << ", " << eta_high << endl; // 3~4, 5~6, 7~8, 9~10, 11~12

      for (int ipt = 0; ipt < nPt; ipt++) {
        int ptBinIdx_Low = h3InvMassXiYvsPt_tot[icent]->GetXaxis()->FindBin(
            fctpTXaxis[ieta][ipt]);
        int ptBinIdx_High = h3InvMassXiYvsPt_tot[icent]->GetXaxis()->FindBin(
                                fctpTXaxis[ieta][ipt + 1]) -
                            1;
        cout << "  pT bin " << ipt << " -- " << ptBinIdx_Low << ", "
             << ptBinIdx_High << endl;

        MinvLam_EPflow_tot[icent][ieta][ipt] =
            (TH1F *)h3InvMassXiYvsPt_tot[icent]->ProjectionZ(
                Form("sig_cent%d_eta%d_pt%d", icent, ieta, ipt), ptBinIdx_Low,
                ptBinIdx_High, eta_low, eta_high);
        MinvLam_EPflow_bg[icent][ieta][ipt] =
            (TH1F *)h3InvMassXiYvsPt_bg[icent]->ProjectionZ(
                Form("bkg_cent%d_eta%d_pt%d", icent, ieta, ipt), ptBinIdx_Low,
                ptBinIdx_High, eta_low, eta_high);
        MinvLam_EPflow_tot[icent][ieta][ipt]->Sumw2();
        MinvLam_EPflow_bg[icent][ieta][ipt]->Sumw2();
      }
    }

    // rapidity bin: (-1.0， -0.8)
    int eta_low =
        h3InvMassXiYvsPt_tot[icent]->GetYaxis()->FindBin(fctetaXaxisFow[0]);
    int eta_high =
        h3InvMassXiYvsPt_tot[icent]->GetYaxis()->FindBin(fctetaXaxisFow[1]);
    cout << "Rapidity bin (-1.0， -0.8): FOW -- " << eta_low << ", " << eta_high
         << endl; // 1~2
    for (int ipt = 0; ipt < nPtFow; ipt++) {
      int ptBinFowIdx_Low =
          h3InvMassXiYvsPt_tot[icent]->GetXaxis()->FindBin(fctpTFowXaxis[ipt]);
      int ptBinFowIdx_High = h3InvMassXiYvsPt_tot[icent]->GetXaxis()->FindBin(
                                 fctpTFowXaxis[ipt + 1]) -
                             1;
      cout << "  pT bin " << ipt << " -- " << ptBinFowIdx_Low << ", "
           << ptBinFowIdx_High << endl;

      MinvLamFow_EPflow_tot[icent][ipt] =
          (TH1F *)h3InvMassXiYvsPt_tot[icent]->ProjectionZ(
              Form("sig_cent%d_fow_pt%d", icent, ipt), ptBinFowIdx_Low,
              ptBinFowIdx_High, eta_low, eta_high);
      MinvLamFow_EPflow_bg[icent][ipt] =
          (TH1F *)h3InvMassXiYvsPt_bg[icent]->ProjectionZ(
              Form("bkg_cent%d_fow_pt%d", icent, ipt), ptBinFowIdx_Low,
              ptBinFowIdx_High, eta_low, eta_high);
      MinvLamFow_EPflow_tot[icent][ipt]->Sumw2();
      MinvLamFow_EPflow_bg[icent][ipt]->Sumw2();
    }

    // mid rapidity bin: (-0.2, 0.2)
    eta_low =
        h3InvMassXiYvsPt_tot[icent]->GetYaxis()->FindBin(fctetaXaxisMid[0]);
    eta_high =
        h3InvMassXiYvsPt_tot[icent]->GetYaxis()->FindBin(fctetaXaxisMid[1]) - 1;
    cout << "Rapidity bin (-0.2， 0.2): MID -- " << eta_low << ", " << eta_high
         << endl; // 9~12
    for (int ipt = 0; ipt < nPtMid; ipt++) {
      int ptBinMidIdx_Low =
          h3InvMassXiYvsPt_tot[icent]->GetXaxis()->FindBin(fctpTMidXaxis[ipt]);
      int ptBinMidIdx_High = h3InvMassXiYvsPt_tot[icent]->GetXaxis()->FindBin(
                                 fctpTMidXaxis[ipt + 1]) -
                             1;
      cout << "  pT bin " << ipt << " -- " << ptBinMidIdx_Low << ", "
           << ptBinMidIdx_High << endl;

      MinvLamMid_EPflow_tot[icent][ipt] =
          (TH1F *)h3InvMassXiYvsPt_tot[icent]->ProjectionZ(
              Form("sig_cent%d_mid_pt%d", icent, ipt), ptBinMidIdx_Low,
              ptBinMidIdx_High, eta_low, eta_high);
      MinvLamMid_EPflow_bg[icent][ipt] =
          (TH1F *)h3InvMassXiYvsPt_bg[icent]->ProjectionZ(
              Form("bkg_cent%d_mid_pt%d", icent, ipt), ptBinMidIdx_Low,
              ptBinMidIdx_High, eta_low, eta_high);
      MinvLamMid_EPflow_tot[icent][ipt]->Sumw2();
      MinvLamMid_EPflow_bg[icent][ipt]->Sumw2();
    }
  }

  // get acceptance histograms
  TH2F *h2XiYvsPt_tot[nCentFlow]; // 0-10%, 10-40%, 40-60%
  TH2F *h2XiYvsPt_bg[nCentFlow];
  for (int icent = 0; icent < nCentFlow; icent++) {
    h3InvMassXiYvsPt_tot[icent]->GetZaxis()->SetRange(16,
                                                      45); // invMass 0.75~1.05
    h3InvMassXiYvsPt_bg[icent]->GetZaxis()->SetRange(16, 45);

    h2XiYvsPt_tot[icent] = (TH2F *)h3InvMassXiYvsPt_tot[icent]->Project3D("xy");
    h2XiYvsPt_bg[icent] = (TH2F *)h3InvMassXiYvsPt_bg[icent]->Project3D("xy");
  }

  TCanvas *ca_acceptance =
      new TCanvas("ca_acceptance", "Acceptance", 450 * nCentFlow, 400);
  ca_acceptance->Divide(nCentFlow, 1);
  TCanvas *efficiencyCheck =
      new TCanvas("efficiency_Check", "Efficiency", 450 * nCentFlow, 400);
  efficiencyCheck->Divide(nCentFlow, 1);

  for (int icent = 0; icent < nCentFlow; icent++) {
    ca_acceptance->cd(icent + 1)->SetLogz();
    //	Efficiency[icent]->GetXaxis()->SetRangeUser(-1.0, 0.4);
    //	Efficiency[icent]->GetYaxis()->SetRangeUser(0, 3.0);
    //	h2XiYvsPt_tot[icent]->Divide(Efficiency[icent]);
    h2XiYvsPt_tot[icent]->Draw("COLZ");
    h2XiYvsPt_tot[icent]->GetXaxis()->SetRangeUser(-1.0, 0.4);
    h2XiYvsPt_tot[icent]->GetYaxis()->SetRangeUser(0, 3.0);
    h2XiYvsPt_tot[icent]->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    h2XiYvsPt_tot[icent]->GetXaxis()->SetTitle("Rapidity");
    efficiencyCheck->cd(icent + 1)->SetLogz();
    Efficiency[icent]->GetXaxis()->SetRangeUser(-1.0, 0.4);
    Efficiency[icent]->GetYaxis()->SetRangeUser(0, 3.0);
    Efficiency[icent]->Draw("COLZ");
  }
  ca_acceptance->Print(
      Form("../../files/plots_%s/acceptance_%s_%s_EPmethod_fxt3GeV.png",
           cutName.Data(), particleName.Data(), cutName.Data()));
  efficiencyCheck->Print(
      Form("../../files/plots_%s/efficency_%s_%s_EPmethod_fxt3GeV.png",
           cutName.Data(), particleName.Data(), cutName.Data()));
  //    return;
  //==========================================================================================================================================
  // extract Xi signal

  TF1 *f1_cent;
  TF1 *f2_bg;
  TF1 *f2_sig;

  // rapidity bins: (-0.8， 0.2)
  double IMcut_3Sigma[nCentFlow][nEta][nPt][2];
  double FitQuality[nCentFlow][nEta][nPt]
                   [6]; // Chi^2, NDF, m0, m0Err, Gamma, GammaErr

  TH1F *MinvLam_EPflow_sig[nCentFlow][nEta][nPt];
  TH1F *MinvLam_EPflowReal_sig[nCentFlow][nEta][nPt];
  TF1 *fMinvLam_EPflow_tot[nCentFlow][nEta][nPt];
  TF1 *fMinvLam_EPflow_bg[nCentFlow][nEta][nPt];
  TF1 *fMinvLam_EPflow_sig[nCentFlow][nEta][nPt];

  // rapidity bin: (-1.0, -0.8)
  double IMcutFow_3Sigma[nCentFlow][nPtFow][2];
  double FitQualityFow[nCentFlow][nPtFow]
                      [6]; // Chi^2, NDF, m0, m0Err, Gamma, GammaErr

  TH1F *MinvLamFow_EPflow_sig[nCentFlow][nPtFow];
  TH1F *MinvLamFow_EPflowReal_sig[nCentFlow][nPtFow];
  TF1 *fMinvLamFow_EPflow_tot[nCentFlow][nPtFow];
  TF1 *fMinvLamFow_EPflow_bg[nCentFlow][nPtFow];
  TF1 *fMinvLamFow_EPflow_sig[nCentFlow][nPtFow];

  // seperate Mid and other rapidity in case of different pT binnning
  // rapidity bin: (-0.2, 0.2)
  double IMcutMid_3Sigma[nCentFlow][nPtMid][2];
  double FitQualityMid[nCentFlow][nPtMid]
                      [6]; // Chi^2, NDF, m0, m0Err, Gamma, GammaErr

  TH1F *MinvLamMid_EPflow_sig[nCentFlow][nPtMid];
  TH1F *MinvLamMid_EPflowReal_sig[nCentFlow][nPtMid];
  TF1 *fMinvLamMid_EPflow_tot[nCentFlow][nPtMid];
  TF1 *fMinvLamMid_EPflow_bg[nCentFlow][nPtMid];
  TF1 *fMinvLamMid_EPflow_sig[nCentFlow][nPtMid];

  double lowFitRange = 0.75;
  double highFitRange = 1.20;

  // from inclusive K*0 invariant mass spectrum with 0.4 < pT < 2.0 and -0.5 < y
  // < 0 for 4 centrality bins
  double initMassMean[] = {0.883, 0.884, 0.887, 0.885};
  double initMassMeanErr[] = {0.006, 0.003, 0.002, 0.002};
  double initMassWidth[] = {0.044, 0.037, 0.037, 0.041};
  double initMassWidthErr[] = {0.014, 0.008, 0.008, 0.007};

  // rapidity bins: (-0.8， 0.2)
  for (int icent = 0; icent < nCentFlow; icent++) {
    for (int ieta = 0; ieta < nEta; ieta++) {
      cout << "######################## FITTING INVARIANT MASS BIN: "
           << centName[icent] << " & " << etaName[ieta]
           << " ########################" << endl;
      for (int ipt = 0; ipt < nPt; ipt++) {
        double nInclLambda_r = MinvLam_EPflow_tot[icent][ieta][ipt]->Integral(
            MinvLam_EPflow_tot[icent][ieta][ipt]->FindBin(1.2),
            MinvLam_EPflow_tot[icent][ieta][ipt]->FindBin(1.5));
        double nBgLambda_r = MinvLam_EPflow_bg[icent][ieta][ipt]->Integral(
            MinvLam_EPflow_bg[icent][ieta][ipt]->FindBin(1.2),
            MinvLam_EPflow_bg[icent][ieta][ipt]->FindBin(1.5));
        double nInclLambda_l = MinvLam_EPflow_tot[icent][ieta][ipt]->Integral(
            MinvLam_EPflow_tot[icent][ieta][ipt]->FindBin(0.6),
            MinvLam_EPflow_tot[icent][ieta][ipt]->FindBin(0.75));
        double nBgLambda_l = MinvLam_EPflow_bg[icent][ieta][ipt]->Integral(
            MinvLam_EPflow_bg[icent][ieta][ipt]->FindBin(0.6),
            MinvLam_EPflow_bg[icent][ieta][ipt]->FindBin(0.75));
        double scaler = 0;
        // if((nBgLambda_l + nBgLambda_r) !=0) scaler= (nInclLambda_l +
        // nInclLambda_r) / (nBgLambda_l + nBgLambda_r);
        if (nBgLambda_r != 0)
          scaler = nInclLambda_r / nBgLambda_r;
        MinvLam_EPflow_bg[icent][ieta][ipt]->Scale(scaler);

        // combinatorial method: mixing or rotation
        MinvLam_EPflow_sig[icent][ieta][ipt] =
            SubtractBG(MinvLam_EPflow_tot[icent][ieta][ipt],
                       MinvLam_EPflow_bg[icent][ieta][ipt]);

        f1_cent = new TF1("f1_cent", fitFuncbw, lowFitRange, highFitRange, 6);
        f1_cent->SetParNames("Scale1", "Gamma", "Mean", "p0", "p1", "p2");
        // f1_cent->SetParLimits(0, 0., 50000);
        f1_cent->SetParameter(1, initMassWidth[icent]);
        f1_cent->SetParLimits(
            1, initMassWidth[icent] - 1.0 * initMassWidthErr[icent],
            initMassWidth[icent] + 1.0 * initMassWidthErr[icent]);
        // f1_cent->SetParLimits(1, 0.03, 0.07);
        f1_cent->SetParameter(2, initMassMean[icent]);
        f1_cent->SetParLimits(
            2, initMassMean[icent] - 5.0 * initMassMeanErr[icent],
            initMassMean[icent] + 5.0 * initMassMeanErr[icent]);
        // f1_cent->SetParLimits(2, 0.85, 0.93);
        f1_cent->SetLineColor(kBlack);
        // f1_cent->SetNpx(50);

        f2_bg = new TF1("f2_bg", background, lowFitRange, highFitRange, 3);
        // f2_bg = new TF1("f2_bg", background, 1.3, 1.35, 2);
        f2_bg->SetLineColor(kBlue);
        f2_bg->SetLineStyle(2);
        // f2_bg->SetNpx(50);

        MinvLam_EPflow_sig[icent][ieta][ipt]->Fit("f1_cent", "MO0", "",
                                                  lowFitRange, highFitRange);
        f1_cent->SetParameters(
            f1_cent->GetParameter("Scale1"), f1_cent->GetParameter("Gamma"),
            f1_cent->GetParameter("Mean"), f1_cent->GetParameter("p0"),
            f1_cent->GetParameter("p1"), f1_cent->GetParameter("p2"));
        MinvLam_EPflow_sig[icent][ieta][ipt]->Fit("f1_cent", "MO0", "",
                                                  lowFitRange, highFitRange);

        double *par;
        const double *parErr;
        par = f1_cent->GetParameters();
        parErr = f1_cent->GetParErrors();
        f2_bg->SetParameters(&par[3]);
        f2_bg->SetParErrors(&parErr[3]);

        f2_sig = new TF1("f2_sig", breitwigner, lowFitRange, highFitRange,
                         3); // order
        f2_sig->SetLineColor(kRed);
        f2_sig->SetLineStyle(1);
        // f2_sig->SetNpx(50);
        f2_sig->SetParameters(&par[0]);
        f2_sig->SetParErrors(&parErr[0]);

        fMinvLam_EPflow_bg[icent][ieta][ipt] = (TF1 *)f2_bg;
        fMinvLam_EPflow_sig[icent][ieta][ipt] = (TF1 *)f2_sig;
        fMinvLam_EPflow_tot[icent][ieta][ipt] = (TF1 *)f1_cent;

        FitQuality[icent][ieta][ipt][0] = f1_cent->GetChisquare();
        FitQuality[icent][ieta][ipt][1] = f1_cent->GetNDF();
        FitQuality[icent][ieta][ipt][2] = par[2]; // mean
        FitQuality[icent][ieta][ipt][3] = parErr[2];
        FitQuality[icent][ieta][ipt][4] = par[1]; // gamma
        FitQuality[icent][ieta][ipt][5] = parErr[1];

        IMcut_3Sigma[icent][ieta][ipt][0] = par[2] - 2.0 * par[1];
        IMcut_3Sigma[icent][ieta][ipt][1] = par[2] + 2.0 * par[1];
        // if(f1_cent->GetParameter("Width") < 0.02) {
        IMcut_3Sigma[icent][ieta][ipt][0] =
            initMassMean[icent] - 2.0 * initMassWidth[icent];
        IMcut_3Sigma[icent][ieta][ipt][1] =
            initMassMean[icent] + 2.0 * initMassWidth[icent];
        //}
      }
    }
    //========================================================================

    // rapdity bin: (-1.0, -0.8)
    cout << "######################## FITTING INVARIANT MASS BIN: "
         << centName[icent] << " & " << etaNameFow[0]
         << " ########################" << endl;
    for (int ipt = 0; ipt < nPtFow; ipt++) {
      double nInclLambda_r = MinvLamFow_EPflow_tot[icent][ipt]->Integral(
          MinvLamFow_EPflow_tot[icent][ipt]->FindBin(1.2),
          MinvLamFow_EPflow_tot[icent][ipt]->FindBin(1.5));
      double nBgLambda_r = MinvLamFow_EPflow_bg[icent][ipt]->Integral(
          MinvLamFow_EPflow_bg[icent][ipt]->FindBin(1.2),
          MinvLamFow_EPflow_bg[icent][ipt]->FindBin(1.5));
      double nInclLambda_l = MinvLamFow_EPflow_tot[icent][ipt]->Integral(
          MinvLamFow_EPflow_tot[icent][ipt]->FindBin(0.6),
          MinvLamFow_EPflow_tot[icent][ipt]->FindBin(0.75));
      double nBgLambda_l = MinvLamFow_EPflow_bg[icent][ipt]->Integral(
          MinvLamFow_EPflow_bg[icent][ipt]->FindBin(0.6),
          MinvLamFow_EPflow_bg[icent][ipt]->FindBin(0.75));
      double scaler = 0;
      // if((nBgLambda_l + nBgLambda_r) !=0) scaler= (nInclLambda_l +
      // nInclLambda_r) / (nBgLambda_l + nBgLambda_r);
      if (nBgLambda_r != 0)
        scaler = nInclLambda_r / nBgLambda_r;
      MinvLamFow_EPflow_bg[icent][ipt]->Scale(scaler);

      // Rotation method
      MinvLamFow_EPflow_sig[icent][ipt] = SubtractBG(
          MinvLamFow_EPflow_tot[icent][ipt], MinvLamFow_EPflow_bg[icent][ipt]);

      f1_cent = new TF1("f1_cent", fitFuncbw, lowFitRange, highFitRange, 6);
      f1_cent->SetParNames("Scale1", "Gamma", "Mean", "p0", "p1", "p2");
      // f1_cent->SetParLimits(0, 0., 50000);
      f1_cent->SetParameter(1, initMassWidth[icent]);
      f1_cent->SetParLimits(
          1, initMassWidth[icent] - 1.0 * initMassWidthErr[icent],
          initMassWidth[icent] + 1.0 * initMassWidthErr[icent]);
      // f1_cent->SetParLimits(1, 0.03, 0.07);
      f1_cent->SetParameter(2, initMassMean[icent]);
      f1_cent->SetParLimits(2,
                            initMassMean[icent] - 5.0 * initMassMeanErr[icent],
                            initMassMean[icent] + 5.0 * initMassMeanErr[icent]);
      // f1_cent->SetParLimits(2, 0.85, 0.93);
      f1_cent->SetLineColor(kBlack);
      // f1_cent->SetNpx(50);

      f2_bg = new TF1("f2_bg", background, lowFitRange, highFitRange, 3);
      f2_bg->SetLineColor(kBlue);
      f2_bg->SetLineStyle(2);
      // f2_bg->SetNpx(50);

      MinvLamFow_EPflow_sig[icent][ipt]->Fit("f1_cent", "MO0", "", lowFitRange,
                                             highFitRange);
      f1_cent->SetParameters(
          f1_cent->GetParameter("Scale1"), f1_cent->GetParameter("Gamma"),
          f1_cent->GetParameter("Mean"), f1_cent->GetParameter("p0"),
          f1_cent->GetParameter("p1"), f1_cent->GetParameter("p2"));
      MinvLamFow_EPflow_sig[icent][ipt]->Fit("f1_cent", "MO0", "", lowFitRange,
                                             highFitRange);

      double *par;
      const double *parErr;
      par = f1_cent->GetParameters();
      parErr = f1_cent->GetParErrors();
      f2_bg->SetParameters(&par[3]);
      f2_bg->SetParErrors(&parErr[3]);

      f2_sig =
          new TF1("f2_sig", breitwigner, lowFitRange, highFitRange, 3); // order
      f2_sig->SetLineColor(kRed);
      f2_sig->SetLineStyle(1);
      // f2_sig->SetNpx(50);
      f2_sig->SetParameters(&par[0]);
      f2_sig->SetParErrors(&parErr[0]);

      fMinvLamFow_EPflow_bg[icent][ipt] = (TF1 *)f2_bg;
      fMinvLamFow_EPflow_sig[icent][ipt] = (TF1 *)f2_sig;
      fMinvLamFow_EPflow_tot[icent][ipt] = (TF1 *)f1_cent;

      FitQualityFow[icent][ipt][0] = f1_cent->GetChisquare();
      FitQualityFow[icent][ipt][1] = f1_cent->GetNDF();
      FitQualityFow[icent][ipt][2] = par[2];
      FitQualityFow[icent][ipt][3] = parErr[2];
      FitQualityFow[icent][ipt][4] = par[1];
      FitQualityFow[icent][ipt][5] = parErr[1];

      IMcutFow_3Sigma[icent][ipt][0] = par[2] - 2.0 * par[1];
      IMcutFow_3Sigma[icent][ipt][1] = par[2] + 2.0 * par[1];
      // if(f1_cent->GetParameter("Width") < 0.02) {
      IMcutFow_3Sigma[icent][ipt][0] =
          initMassMean[icent] - 2.0 * initMassWidth[icent];
      IMcutFow_3Sigma[icent][ipt][1] =
          initMassMean[icent] + 2.0 * initMassWidth[icent];
      //}
    }
    //========================================================================

    // mid rapdity bin
    cout << "######################## FITTING INVARIANT MASS BIN: "
         << centName[icent] << " & " << etaNameMid[0]
         << " ########################" << endl;
    for (int ipt = 0; ipt < nPtMid; ipt++) {
      double nInclLambda_r = MinvLamMid_EPflow_tot[icent][ipt]->Integral(
          MinvLamMid_EPflow_tot[icent][ipt]->FindBin(1.2),
          MinvLamMid_EPflow_tot[icent][ipt]->FindBin(1.5));
      double nBgLambda_r = MinvLamMid_EPflow_bg[icent][ipt]->Integral(
          MinvLamMid_EPflow_bg[icent][ipt]->FindBin(1.2),
          MinvLamMid_EPflow_bg[icent][ipt]->FindBin(1.5));
      double nInclLambda_l = MinvLamMid_EPflow_tot[icent][ipt]->Integral(
          MinvLamMid_EPflow_tot[icent][ipt]->FindBin(0.6),
          MinvLamMid_EPflow_tot[icent][ipt]->FindBin(0.75));
      double nBgLambda_l = MinvLamMid_EPflow_bg[icent][ipt]->Integral(
          MinvLamMid_EPflow_bg[icent][ipt]->FindBin(0.6),
          MinvLamMid_EPflow_bg[icent][ipt]->FindBin(0.75));
      double scaler = 0;
      // if((nBgLambda_l + nBgLambda_r) !=0) scaler= (nInclLambda_l +
      // nInclLambda_r) / (nBgLambda_l + nBgLambda_r);
      if (nBgLambda_r != 0)
        scaler = nInclLambda_r / nBgLambda_r;
      MinvLamMid_EPflow_bg[icent][ipt]->Scale(scaler);

      // Rotation method
      MinvLamMid_EPflow_sig[icent][ipt] = SubtractBG(
          MinvLamMid_EPflow_tot[icent][ipt], MinvLamMid_EPflow_bg[icent][ipt]);

      f1_cent = new TF1("f1_cent", fitFuncbw, lowFitRange, highFitRange, 6);
      f1_cent->SetParNames("Scale1", "Gamma", "Mean", "p0", "p1", "p2");
      // f1_cent->SetParLimits(0, 0., 50000);
      f1_cent->SetParameter(1, initMassWidth[icent]);
      f1_cent->SetParLimits(
          1, initMassWidth[icent] - 1.0 * initMassWidthErr[icent],
          initMassWidth[icent] + 1.0 * initMassWidthErr[icent]);
      // f1_cent->SetParLimits(1, 0.03, 0.09);
      f1_cent->SetParameter(2, initMassMean[icent]);
      f1_cent->SetParLimits(2,
                            initMassMean[icent] - 5.0 * initMassMeanErr[icent],
                            initMassMean[icent] + 5.0 * initMassMeanErr[icent]);
      // f1_cent->SetParLimits(2, 0.85, 0.93);
      f1_cent->SetLineColor(kBlack);
      // f1_cent->SetNpx(50);

      f2_bg = new TF1("f2_bg", background, lowFitRange, highFitRange, 3);
      f2_bg->SetLineColor(kBlue);
      f2_bg->SetLineStyle(2);
      // f2_bg->SetNpx(50);

      MinvLamMid_EPflow_sig[icent][ipt]->Fit("f1_cent", "MO0", "", lowFitRange,
                                             highFitRange);
      f1_cent->SetParameters(
          f1_cent->GetParameter("Scale1"), f1_cent->GetParameter("Gamma"),
          f1_cent->GetParameter("Mean"), f1_cent->GetParameter("p0"),
          f1_cent->GetParameter("p1"), f1_cent->GetParameter("p2"));
      MinvLamMid_EPflow_sig[icent][ipt]->Fit("f1_cent", "MO0", "", lowFitRange,
                                             highFitRange);

      double *par;
      const double *parErr;
      par = f1_cent->GetParameters();
      parErr = f1_cent->GetParErrors();
      f2_bg->SetParameters(&par[3]);
      f2_bg->SetParErrors(&parErr[3]);

      f2_sig =
          new TF1("f2_sig", breitwigner, lowFitRange, highFitRange, 3); // order
      f2_sig->SetLineColor(kRed);
      f2_sig->SetLineStyle(1);
      // f2_sig->SetNpx(50);
      f2_sig->SetParameters(&par[0]);
      f2_sig->SetParErrors(&parErr[0]);

      fMinvLamMid_EPflow_bg[icent][ipt] = (TF1 *)f2_bg;
      fMinvLamMid_EPflow_sig[icent][ipt] = (TF1 *)f2_sig;
      fMinvLamMid_EPflow_tot[icent][ipt] = (TF1 *)f1_cent;

      FitQualityMid[icent][ipt][0] = f1_cent->GetChisquare();
      FitQualityMid[icent][ipt][1] = f1_cent->GetNDF();
      FitQualityMid[icent][ipt][2] = par[2];
      FitQualityMid[icent][ipt][3] = parErr[2];
      FitQualityMid[icent][ipt][4] = par[1];
      FitQualityMid[icent][ipt][5] = parErr[1];

      IMcutMid_3Sigma[icent][ipt][0] = par[2] - 2.0 * par[1];
      IMcutMid_3Sigma[icent][ipt][1] = par[2] + 2.0 * par[1];
      // if(f1_cent->GetParameter("Width") < 0.02) {
      IMcutMid_3Sigma[icent][ipt][0] =
          initMassMean[icent] - 2.0 * initMassWidth[icent];
      IMcutMid_3Sigma[icent][ipt][1] =
          initMassMean[icent] + 2.0 * initMassWidth[icent];
      //}
    }
  }
  //==========================================================================================================================================
  TLine *lineL[nPtMid];
  TLine *lineR[nPtMid];
  float eff_weight;

  double nXiYield[nCentFlow][nEta][nPt] = {0.};
  double nXiYieldErr[nCentFlow][nEta][nPt] = {0.};

  TCanvas *ca_invMass[nCentFlow][nEta];
  for (int icent = 0; icent < nCentFlow; icent++) {
    for (int ieta = 0; ieta < nEta; ieta++) {
      ca_invMass[icent][ieta] = new TCanvas(
          Form("ca_invMass_cent%d_eta%d", icent, ieta),
          Form("ca_invMass_cent%d_eta%d", icent, ieta), 450 * nPt, 800);
      ca_invMass[icent][ieta]->Divide(nPt, 2);
    }
  }
  for (int icent = 0; icent < nCentFlow; icent++) {
    for (int ieta = 0; ieta < nEta; ieta++) {
      for (int ipt = 0; ipt < nPt; ipt++) {
        ca_invMass[icent][ieta]->cd(ipt + 1);

        float maxCounts = MinvLam_EPflow_tot[icent][ieta][ipt]->GetBinContent(
            MinvLam_EPflow_tot[icent][ieta][ipt]->GetMaximumBin());
        float maxCountsErr = MinvLam_EPflow_tot[icent][ieta][ipt]->GetBinError(
            MinvLam_EPflow_tot[icent][ieta][ipt]->GetMaximumBin());

        float lowerEdgeY = -0.15 * (maxCounts + maxCountsErr);
        float upperEdgeY = 1.15 * (maxCounts + maxCountsErr);

        ax = gPad->DrawFrame(0.6, lowerEdgeY, 1.5, upperEdgeY); // no EPcorr
        SetAxis(ax, 1.5, 1.8);
        ax->SetXTitle("M_{inv} [GeV/c^{2}] ");
        ax->SetYTitle("Counts");
        ax->Draw();

        TLine *unitaty = new TLine(0.6, 0., 1.5, 0.);
        unitaty->Draw();

        MinvLam_EPflow_tot[icent][ieta][ipt]->Draw("ESAME");
        MinvLam_EPflow_tot[icent][ieta][ipt]->SetMarkerStyle(20);
        MinvLam_EPflow_tot[icent][ieta][ipt]->SetMarkerSize(1.0);
        MinvLam_EPflow_tot[icent][ieta][ipt]->SetMarkerColor(kBlack);

        MinvLam_EPflow_bg[icent][ieta][ipt]->Draw("ESAME");
        MinvLam_EPflow_bg[icent][ieta][ipt]->SetMarkerStyle(24);
        MinvLam_EPflow_bg[icent][ieta][ipt]->SetMarkerSize(1.0);
        MinvLam_EPflow_bg[icent][ieta][ipt]->SetMarkerColor(kBlue);

        lineL[ipt] = new TLine(IMcut[0], lowerEdgeY, IMcut[0], upperEdgeY);
        lineL[ipt]->SetLineColor(kBlack);
        lineL[ipt]->SetLineStyle(2);
        lineL[ipt]->SetLineWidth(2);
        lineL[ipt]->Draw("SAME");

        lineR[ipt] = new TLine(IMcut[1], lowerEdgeY, IMcut[1], upperEdgeY);
        lineR[ipt]->SetLineColor(kBlack);
        lineR[ipt]->SetLineStyle(2);
        lineR[ipt]->SetLineWidth(2);
        lineR[ipt]->Draw("SAME");

        TLatex *tex_sys =
            new TLatex(1.1, 0.85 * upperEdgeY,
                       Form("%s 3GeV Au+Au ", centName[icent].Data()));
        tex_sys->SetTextFont(42);
        tex_sys->SetTextSize(0.04);
        tex_sys->Draw("same");

        TLatex *tex_eta = new TLatex(1.1, 0.75 * upperEdgeY,
                                     Form("%s", etaName[ieta].Data()));
        tex_eta->SetTextFont(42);
        tex_eta->SetTextSize(0.04);
        tex_eta->Draw("same");

        TLatex *tex_phi =
            new TLatex(1.1, 0.65 * upperEdgeY,
                       Form("%4.2f < p_{T} < %4.2f GeV/c",
                            fctpTXaxis[ieta][ipt], fctpTXaxis[ieta][ipt + 1]));
        tex_phi->SetTextFont(42);
        tex_phi->SetTextSize(0.04);
        tex_phi->Draw("same");

        ca_invMass[icent][ieta]->cd(nPt + ipt + 1);
        float minCounts, minCountsErr;
        float minTemp = 0., maxTemp = 0.;
        int invMBinIdx_Low =
            MinvLam_EPflow_sig[icent][ieta][ipt]->GetXaxis()->FindBin(
                lowFitRange);
        int invMBinIdx_High =
            MinvLam_EPflow_sig[icent][ieta][ipt]->GetXaxis()->FindBin(
                highFitRange);
        for (int iinvM = invMBinIdx_Low; iinvM < invMBinIdx_High; iinvM++) {
          maxCounts =
              MinvLam_EPflow_sig[icent][ieta][ipt]->GetBinContent(iinvM);
          maxCountsErr =
              MinvLam_EPflow_sig[icent][ieta][ipt]->GetBinError(iinvM);
          if ((maxCounts + maxCountsErr) > maxTemp)
            maxTemp = maxCounts + maxCountsErr;

          minCounts =
              MinvLam_EPflow_sig[icent][ieta][ipt]->GetBinContent(iinvM);
          minCountsErr =
              MinvLam_EPflow_sig[icent][ieta][ipt]->GetBinError(iinvM);
          if ((minCounts + minCountsErr) < minTemp)
            minTemp = minCounts + minCountsErr;
        }

        if (minTemp >= 0)
          lowerEdgeY = -0.15 * maxTemp;
        else
          lowerEdgeY = 1.15 * minTemp;
        if (maxTemp > 0)
          upperEdgeY = 1.15 * maxTemp;
        else
          upperEdgeY = maxTemp - minTemp;

        ax = gPad->DrawFrame(lowFitRange, lowerEdgeY, highFitRange,
                             upperEdgeY); // no EPcorr
        SetAxis(ax, 1.5, 1.8);
        ax->SetXTitle("M_{inv} [GeV/c^{2}] ");
        ax->SetYTitle("Counts");
        ax->Draw();

        MinvLam_EPflow_sig[icent][ieta][ipt]->Draw("E");
        MinvLam_EPflow_sig[icent][ieta][ipt]->GetXaxis()->SetRangeUser(
            lowFitRange, highFitRange);
        MinvLam_EPflow_sig[icent][ieta][ipt]->GetYaxis()->SetRangeUser(
            lowerEdgeY, upperEdgeY);
        MinvLam_EPflow_sig[icent][ieta][ipt]->GetXaxis()->SetTitle(
            "M_{inv} [GeV/c^{2}]");
        MinvLam_EPflow_sig[icent][ieta][ipt]->GetYaxis()->SetTitle("Counts");
        MinvLam_EPflow_sig[icent][ieta][ipt]->SetMarkerStyle(20);
        MinvLam_EPflow_sig[icent][ieta][ipt]->SetMarkerSize(0.5);
        MinvLam_EPflow_sig[icent][ieta][ipt]->SetMarkerColor(kRed);

        TLine *unitaty1 = new TLine(lowFitRange, 0., highFitRange, 0.);
        unitaty1->Draw();

        TLine *Kstar0mass = new TLine(initMassMean[icent], lowerEdgeY,
                                      initMassMean[icent], upperEdgeY);
        Kstar0mass->SetLineColor(kGreen);
        Kstar0mass->SetLineStyle(2);
        Kstar0mass->SetLineWidth(2);
        Kstar0mass->Draw("SAME");

        fMinvLam_EPflow_tot[icent][ieta][ipt]->Draw("SAME");
        fMinvLam_EPflow_sig[icent][ieta][ipt]->Draw("SAME");
        fMinvLam_EPflow_bg[icent][ieta][ipt]->Draw("SAME");

        MinvLam_EPflowReal_sig[icent][ieta][ipt] = SubtractBGFn(
            MinvLam_EPflow_sig[icent][ieta][ipt],
            fMinvLam_EPflow_bg[icent][ieta][ipt], lowFitRange, highFitRange);
        MinvLam_EPflowReal_sig[icent][ieta][ipt]->SetLineWidth(2.0);
        MinvLam_EPflowReal_sig[icent][ieta][ipt]->SetLineColor(kRed);
        MinvLam_EPflowReal_sig[icent][ieta][ipt]->Draw("HIST SAME");

        TLine *lineLFit =
            new TLine(IMcut_3Sigma[icent][ieta][ipt][0], lowerEdgeY,
                      IMcut_3Sigma[icent][ieta][ipt][0], upperEdgeY);
        lineLFit->SetLineColor(kBlack);
        lineLFit->SetLineStyle(2);
        lineLFit->SetLineWidth(2);
        lineLFit->Draw("SAME");

        TLine *lineRFit =
            new TLine(IMcut_3Sigma[icent][ieta][ipt][1], lowerEdgeY,
                      IMcut_3Sigma[icent][ieta][ipt][1], upperEdgeY);
        lineRFit->SetLineColor(kBlack);
        lineRFit->SetLineStyle(2);
        lineRFit->SetLineWidth(2);
        lineRFit->Draw("SAME");

        int lowerMassBin_Kstar0 =
            MinvLam_EPflowReal_sig[icent][ieta][ipt]->FindBin(
                IMcut_3Sigma[icent][ieta][ipt][0]);
        int upperMassBin_Kstar0 =
            MinvLam_EPflowReal_sig[icent][ieta][ipt]->FindBin(
                IMcut_3Sigma[icent][ieta][ipt][1]) -
            1;
        double nSignalErr = 0., nSignal = 0.;
        nSignal = MinvLam_EPflowReal_sig[icent][ieta][ipt]->IntegralAndError(
            lowerMassBin_Kstar0, upperMassBin_Kstar0, nSignalErr);

        if (isnan(nSignal)) {

          // fMinvLam_EPflow_bg[icent][ieta][ipt]->GetXaxis()->SetRangeUser(1.25,1.40);
          TF1 *fBkgCounts = fMinvLam_EPflow_bg[icent][ieta][ipt];
          fBkgCounts->SetNpx(45);
          TH1F *h1_BkgCounts = (TH1F *)fBkgCounts->GetHistogram()->Clone();

          for (int ix = 0;
               ix < MinvLam_EPflow_sig[icent][ieta][ipt]->GetNbinsX(); ix++) {
            float xL =
                MinvLam_EPflow_sig[icent][ieta][ipt]->GetBinLowEdge(ix + 1);
            float xR =
                MinvLam_EPflow_sig[icent][ieta][ipt]->GetBinLowEdge(ix + 2);
            float bgFn =
                fMinvLam_EPflow_bg[icent][ieta][ipt]->Integral(xL, xR) /
                MinvLam_EPflow_sig[icent][ieta][ipt]->GetBinWidth(ix + 1);
            //            cout<<bgFn<<" :   2nd:
            //            "<<h1_BkgCounts->GetBinContent(ix+1)<<"   :
            //            "<<MinvLam_EPflow_sig[icent][ieta][ipt]->GetNbinsX()<<endl;
          }

          TH1F *hBkgCounts =
              SubtractBG(MinvLam_EPflow_sig[icent][ieta][ipt], h1_BkgCounts);
          int lowerMassBinT =
              hBkgCounts->FindBin(IMcut_3Sigma[icent][ieta][ipt][0]);
          int upperMassBinT =
              hBkgCounts->FindBin(IMcut_3Sigma[icent][ieta][ipt][1]) - 1;
          nSignal = hBkgCounts->IntegralAndError(lowerMassBinT, upperMassBinT,
                                                 nSignalErr);
          /*
           double nSignalT = 0., nSignalErrT = 0.;
           for(int i=lowerMassBin_Kstar0; i<upperMassBin_Kstar0+1; i++) {
           nSignalT += MinvLam_EPflow_sig[icent][ieta][ipt]->GetBinContent(i);
           nSignalErrT += MinvLam_EPflow_sig[icent][ieta][ipt]->GetBinError(i) *
           MinvLam_EPflow_sig[icent][ieta][ipt]->GetBinError(i);
           }

           nSignal = nSignalT;
           nSignalErr = sqrt(nSignalErrT);
           */
          //          MinvLam_EPflowReal_sig[icent][ieta][ipt] =
          //          (TH1F*)hBkgCounts->Clone(Form("Sig_cent%d_eta%d_pt%d",
          //          icent,ieta,ipt));
        }

        // if(isnan(nSignal)) nSignal = 0;
        lowerMassBin_Kstar0 = MinvLam_EPflow_tot[icent][ieta][ipt]->FindBin(
            IMcut_3Sigma[icent][ieta][ipt][0]);
        upperMassBin_Kstar0 = MinvLam_EPflow_tot[icent][ieta][ipt]->FindBin(
                                  IMcut_3Sigma[icent][ieta][ipt][1]) -
                              1;
        float nTotal = MinvLam_EPflow_tot[icent][ieta][ipt]->Integral(
            lowerMassBin_Kstar0, upperMassBin_Kstar0);
        float nBackground = nTotal - nSignal;

        // eff_weight = hEff_accReg_we[ieta][icent]->GetBinContent(ipt+1);
        // if(eff_weight!=0)  eff_weight = 1./eff_weight;
        //                eff_weight = 1.0;
        eff_weight =
            1. / EfficiencyBigRange[icent]->GetBinContent(ieta + 1, ipt + 1);

        nXiYield[icent][ieta][ipt] = nSignal * eff_weight;
        nXiYieldErr[icent][ieta][ipt] =
            nSignalErr * eff_weight; // sqrt(nSignal);
        if (icent == 1) {
          if (debug) {
            std::cout << "zla ipt == " << ipt << " ieta == " << ieta
                      << " eff_weight == "
                      << EfficiencyBigRange[icent]->GetBinContent(ieta + 1,
                                                                  ipt + 1)
                      << " nsignal == " << nSignal
                      << " yield == " << nXiYield[icent][ieta][ipt]
                      << std::endl;
          }
        }

        TLatex *tex_chi2 = new TLatex(
            0.91, 0.85 * upperEdgeY,
            Form("#chi^{2}/NDF = %4.2f / %d", FitQuality[icent][ieta][ipt][0],
                 (int)FitQuality[icent][ieta][ipt][1]));
        tex_chi2->SetTextFont(42);
        tex_chi2->SetTextSize(0.04);
        tex_chi2->Draw("same");

        TLatex *tex_signif =
            new TLatex(0.91, 0.75 * upperEdgeY,
                       Form("S / #sqrt{S+B} = %6.2f", nSignal / sqrt(nTotal)));
        tex_signif->SetTextFont(42);
        tex_signif->SetTextSize(0.04);
        tex_signif->Draw("same");

        TLatex *tex_SBR = new TLatex(
            0.91, 0.65 * upperEdgeY,
            Form("S=%5.1f#pm%5.1f, B=%5.1f", nSignal, nSignalErr, nBackground));
        tex_SBR->SetTextFont(42);
        tex_SBR->SetTextSize(0.04);
        tex_SBR->Draw("same");

        TLatex *tex_mean = new TLatex(0.77, 0.45 * upperEdgeY,
                                      Form("#mu=%5.4f#pm%5.4f",
                                           FitQuality[icent][ieta][ipt][2],
                                           FitQuality[icent][ieta][ipt][3]));
        tex_mean->SetTextFont(42);
        tex_mean->SetTextSize(0.04);
        tex_mean->Draw("same");

        TLatex *tex_gamma = new TLatex(0.77, 0.35 * upperEdgeY,
                                       Form("#Gamma=%5.4f#pm%5.4f",
                                            FitQuality[icent][ieta][ipt][4],
                                            FitQuality[icent][ieta][ipt][5]));
        tex_gamma->SetTextFont(42);
        tex_gamma->SetTextSize(0.04);
        tex_gamma->Draw("same");
      }

      ca_invMass[icent][ieta]->cd();
      ca_invMass[icent][ieta]->Print(Form(
          "../../files/plots_%s/"
          "invMassVsYCent%d_dNdY%d_%s_%s_EPmethod_fxt3GeV.png",
          cutName.Data(), icent, ieta, particleName.Data(), cutName.Data()));
    }
  }
  //==========================================================================================================================================
  // rapidity bin: (-1.0, -0.8)
  double nXiYieldFow[nCentFlow][nPtFow] = {0.};
  double nXiYieldErrFow[nCentFlow][nPtFow] = {0.};
  TCanvas *ca_invMassFow[nCentFlow];
  for (int icent = 0; icent < nCentFlow; icent++) {
    ca_invMassFow[icent] =
        new TCanvas(Form("ca_invMass_cent%d_Fow", icent),
                    Form("ca_invMass_cent%d_Fow", icent), 450 * nPtFow, 800);
    ca_invMassFow[icent]->Divide(nPtFow, 2);
  }

  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  gPad->SetBottomMargin(0.1);
  gPad->SetLeftMargin(0.1);

  for (int icent = 0; icent < nCentFlow; icent++) {
    for (int ipt = 0; ipt < nPtFow; ipt++) {
      ca_invMassFow[icent]->cd(ipt + 1);

      float maxCounts = MinvLamFow_EPflow_tot[icent][ipt]->GetBinContent(
          MinvLamFow_EPflow_tot[icent][ipt]->GetMaximumBin());
      float maxCountsErr = MinvLamFow_EPflow_tot[icent][ipt]->GetBinError(
          MinvLamFow_EPflow_tot[icent][ipt]->GetMaximumBin());

      float lowerEdgeY = -0.35 * (maxCounts + maxCountsErr);
      float upperEdgeY = 1.05 * (maxCounts + maxCountsErr);

      ax = gPad->DrawFrame(0.6, lowerEdgeY, 1.5, upperEdgeY); // no EPcorr
      SetAxis(ax, 3.0, 3.0);
      ax->SetXTitle("M_{inv} [GeV/c^{2}] ");
      ax->SetYTitle("Counts");
      ax->Draw();

      TLine *unitaty = new TLine(0.6, 0., 1.5, 0.);
      unitaty->Draw();

      MinvLamFow_EPflow_tot[icent][ipt]->Draw("ESAME");
      MinvLamFow_EPflow_tot[icent][ipt]->SetMarkerStyle(20);
      MinvLamFow_EPflow_tot[icent][ipt]->SetMarkerSize(1.0);
      MinvLamFow_EPflow_tot[icent][ipt]->SetMarkerColor(kBlack);

      MinvLamFow_EPflow_bg[icent][ipt]->Draw("ESAME");
      MinvLamFow_EPflow_bg[icent][ipt]->SetMarkerStyle(24);
      MinvLamFow_EPflow_bg[icent][ipt]->SetMarkerSize(1.0);
      MinvLamFow_EPflow_bg[icent][ipt]->SetMarkerColor(kBlue);

      lineL[ipt] = new TLine(IMcut[0], lowerEdgeY, IMcut[0], upperEdgeY);
      lineL[ipt]->SetLineColor(kBlack);
      lineL[ipt]->SetLineStyle(2);
      lineL[ipt]->SetLineWidth(2);
      lineL[ipt]->Draw("SAME");

      lineR[ipt] = new TLine(IMcut[1], lowerEdgeY, IMcut[1], upperEdgeY);
      lineR[ipt]->SetLineColor(kBlack);
      lineR[ipt]->SetLineStyle(2);
      lineR[ipt]->SetLineWidth(2);
      lineR[ipt]->Draw("SAME");

      TLatex *tex_sys =
          new TLatex(1.1, 0.85 * upperEdgeY,
                     Form("%s 3GeV Au+Au ", centName[icent].Data()));
      tex_sys->SetTextFont(42);
      tex_sys->SetTextSize(0.04);
      tex_sys->Draw("same");

      TLatex *tex_eta =
          new TLatex(1.1, 0.75 * upperEdgeY, Form("%s", etaNameFow[0].Data()));
      tex_eta->SetTextFont(42);
      tex_eta->SetTextSize(0.04);
      tex_eta->Draw("same");

      TLatex *tex_phi =
          new TLatex(1.1, 0.65 * upperEdgeY,
                     Form("%4.2f < p_{T} < %4.2f GeV/c", fctpTFowXaxis[ipt],
                          fctpTFowXaxis[ipt + 1]));
      tex_phi->SetTextFont(42);
      tex_phi->SetTextSize(0.04);
      tex_phi->Draw("same");

      ca_invMassFow[icent]->cd(nPtFow + ipt + 1);
      float minCounts, minCountsErr;
      float minTemp = 0., maxTemp = 0.;
      int invMBinIdx_Low =
          MinvLamFow_EPflow_sig[icent][ipt]->GetXaxis()->FindBin(lowFitRange);
      int invMBinIdx_High =
          MinvLamFow_EPflow_sig[icent][ipt]->GetXaxis()->FindBin(highFitRange);
      for (int iinvM = invMBinIdx_Low; iinvM < invMBinIdx_High; iinvM++) {
        maxCounts = MinvLamFow_EPflow_sig[icent][ipt]->GetBinContent(iinvM);
        maxCountsErr = MinvLamFow_EPflow_sig[icent][ipt]->GetBinError(iinvM);
        if ((maxCounts + maxCountsErr) > maxTemp)
          maxTemp = maxCounts + maxCountsErr;

        minCounts = MinvLamFow_EPflow_sig[icent][ipt]->GetBinContent(iinvM);
        minCountsErr = MinvLamFow_EPflow_sig[icent][ipt]->GetBinError(iinvM);
        if ((minCounts + minCountsErr) < minTemp)
          minTemp = minCounts + minCountsErr;
      }

      if (minTemp >= 0)
        lowerEdgeY = -0.15 * maxTemp;
      else
        lowerEdgeY = 1.15 * minTemp;
      if (maxTemp > 0)
        upperEdgeY = 1.15 * maxTemp;
      else
        upperEdgeY = maxTemp - minTemp;

      ax = gPad->DrawFrame(lowFitRange, lowerEdgeY, highFitRange,
                           upperEdgeY); // no EPcorr
      SetAxis(ax, 1.5, 1.8);
      ax->SetXTitle("M_{inv} [GeV/c^{2}] ");
      ax->SetYTitle("Counts");
      ax->Draw();

      MinvLamFow_EPflow_sig[icent][ipt]->Draw("E");
      MinvLamFow_EPflow_sig[icent][ipt]->GetXaxis()->SetRangeUser(lowFitRange,
                                                                  highFitRange);
      MinvLamFow_EPflow_sig[icent][ipt]->GetYaxis()->SetRangeUser(lowerEdgeY,
                                                                  upperEdgeY);
      MinvLamFow_EPflow_sig[icent][ipt]->GetXaxis()->SetTitle(
          "M_{inv} [GeV/c^{2}]");
      MinvLamFow_EPflow_sig[icent][ipt]->GetYaxis()->SetTitle("Counts");
      MinvLamFow_EPflow_sig[icent][ipt]->SetMarkerStyle(20);
      MinvLamFow_EPflow_sig[icent][ipt]->SetMarkerSize(0.5);
      MinvLamFow_EPflow_sig[icent][ipt]->SetMarkerColor(kRed);

      TLine *unitaty1 = new TLine(lowFitRange, 0., highFitRange, 0.);
      unitaty1->Draw();

      TLine *Kstar0mass = new TLine(initMassMean[icent], lowerEdgeY,
                                    initMassMean[icent], upperEdgeY);
      Kstar0mass->SetLineColor(kGreen);
      Kstar0mass->SetLineStyle(2);
      Kstar0mass->SetLineWidth(2);
      Kstar0mass->Draw("SAME");

      fMinvLamFow_EPflow_tot[icent][ipt]->Draw("SAME");
      fMinvLamFow_EPflow_sig[icent][ipt]->Draw("SAME");
      fMinvLamFow_EPflow_bg[icent][ipt]->Draw("SAME");

      MinvLamFow_EPflowReal_sig[icent][ipt] = SubtractBGFn(
          MinvLamFow_EPflow_sig[icent][ipt], fMinvLamFow_EPflow_bg[icent][ipt],
          lowFitRange, highFitRange);
      MinvLamFow_EPflowReal_sig[icent][ipt]->SetLineWidth(2.0);
      MinvLamFow_EPflowReal_sig[icent][ipt]->SetLineColor(kRed);
      MinvLamFow_EPflowReal_sig[icent][ipt]->Draw("HIST SAME");

      TLine *lineLFit = new TLine(IMcutFow_3Sigma[icent][ipt][0], lowerEdgeY,
                                  IMcutFow_3Sigma[icent][ipt][0], upperEdgeY);
      lineLFit->SetLineColor(kBlack);
      lineLFit->SetLineStyle(2);
      lineLFit->SetLineWidth(2);
      lineLFit->Draw("SAME");

      TLine *lineRFit = new TLine(IMcutFow_3Sigma[icent][ipt][1], lowerEdgeY,
                                  IMcutFow_3Sigma[icent][ipt][1], upperEdgeY);
      lineRFit->SetLineColor(kBlack);
      lineRFit->SetLineStyle(2);
      lineRFit->SetLineWidth(2);
      lineRFit->Draw("SAME");

      int lowerMassBin_Kstar0 = MinvLamFow_EPflowReal_sig[icent][ipt]->FindBin(
          IMcutFow_3Sigma[icent][ipt][0]);
      int upperMassBin_Kstar0 = MinvLamFow_EPflowReal_sig[icent][ipt]->FindBin(
                                    IMcutFow_3Sigma[icent][ipt][1]) -
                                1;
      double nSignalErr = 0., nSignal = 0.;
      nSignal = MinvLamFow_EPflowReal_sig[icent][ipt]->IntegralAndError(
          lowerMassBin_Kstar0, upperMassBin_Kstar0, nSignalErr);
      // nSignal =
      // MinvLamFow_EPflow_sig[icent][ipt]->IntegralAndError(lowerMassBin_Kstar0,
      // upperMassBin_Kstar0, nSignalErr); for(int i=lowerMassBin_Kstar0;
      // i<upperMassBin_Kstar0+1; i++) {
      //    nSignal += MinvLamFow_EPflowReal_sig[icent][ipt]->GetBinContent(i);
      //}

      // eff_weight = hEff_accReg_we_fow[icent]->GetBinContent(ipt+1);
      // if(eff_weight !=0) eff_weight = 1./ eff_weight;
      //            eff_weight = 1.;
      eff_weight = 1. / EfficiencyFow[icent]->GetBinContent(1, ipt + 1);

      nXiYieldFow[icent][ipt] = nSignal * eff_weight;
      nXiYieldErrFow[icent][ipt] = nSignalErr * eff_weight; // sqrt(nSignal);

      lowerMassBin_Kstar0 = MinvLamFow_EPflow_tot[icent][ipt]->FindBin(
          IMcutFow_3Sigma[icent][ipt][0]);
      upperMassBin_Kstar0 = MinvLamFow_EPflow_tot[icent][ipt]->FindBin(
                                IMcutFow_3Sigma[icent][ipt][1]) -
                            1;
      float nTotal = MinvLamFow_EPflow_tot[icent][ipt]->Integral(
          lowerMassBin_Kstar0, upperMassBin_Kstar0);
      float nBackground = nTotal - nSignal;

      TLatex *tex_chi2 = new TLatex(0.91, 0.85 * upperEdgeY,
                                    Form("#chi^{2}/NDF = %4.2f / %d",
                                         FitQualityFow[icent][ipt][0],
                                         (int)FitQualityFow[icent][ipt][1]));
      tex_chi2->SetTextFont(42);
      tex_chi2->SetTextSize(0.04);
      tex_chi2->Draw("same");

      TLatex *tex_signif =
          new TLatex(0.91, 0.75 * upperEdgeY,
                     Form("S / #sqrt{S+B} = %6.2f", nSignal / sqrt(nTotal)));
      tex_signif->SetTextFont(42);
      tex_signif->SetTextSize(0.04);
      tex_signif->Draw("same");

      TLatex *tex_SBR = new TLatex(
          0.91, 0.65 * upperEdgeY,
          Form("S=%5.1f#pm%5.1f, B=%5.1f", nSignal, nSignalErr, nBackground));
      tex_SBR->SetTextFont(42);
      tex_SBR->SetTextSize(0.04);
      tex_SBR->Draw("same");

      TLatex *tex_mean =
          new TLatex(0.77, 0.45 * upperEdgeY,
                     Form("#mu=%5.4f#pm%5.4f", FitQualityFow[icent][ipt][2],
                          FitQualityFow[icent][ipt][3]));
      tex_mean->SetTextFont(42);
      tex_mean->SetTextSize(0.04);
      tex_mean->Draw("same");

      TLatex *tex_gamma =
          new TLatex(0.77, 0.35 * upperEdgeY,
                     Form("#Gamma=%5.4f#pm%5.4f", FitQualityFow[icent][ipt][4],
                          FitQualityFow[icent][ipt][5]));
      tex_gamma->SetTextFont(42);
      tex_gamma->SetTextSize(0.04);
      tex_gamma->Draw("same");
    }

    ca_invMassFow[icent]->cd();
    ca_invMassFow[icent]->Print(
        Form("../../files/plots_%s/"
             "invMassVsYCent%d_dNdYFow_%s_%s_EPmethod_fxt3GeV.png",
             cutName.Data(), icent, particleName.Data(), cutName.Data()));
  }
  //==========================================================================================================================================

  // mid rapidity bin
  double nXiYieldMid[nCentFlow][nPtMid] = {0.};
  double nXiYieldErrMid[nCentFlow][nPtMid] = {0.};
  TCanvas *ca_invMassMid[nCentFlow];
  for (int icent = 0; icent < nCentFlow; icent++) {
    ca_invMassMid[icent] =
        new TCanvas(Form("ca_invMass_cent%d_mid", icent),
                    Form("ca_invMass_cent%d_mid", icent), 450 * nPtMid, 800);
    ca_invMassMid[icent]->Divide(nPtMid, 2);
  }

  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  gPad->SetBottomMargin(0.1);
  gPad->SetLeftMargin(0.1);

  for (int icent = 0; icent < nCentFlow; icent++) {
    for (int ipt = 0; ipt < nPtMid; ipt++) {
      ca_invMassMid[icent]->cd(ipt + 1);

      float maxCounts = MinvLamMid_EPflow_tot[icent][ipt]->GetBinContent(
          MinvLamMid_EPflow_tot[icent][ipt]->GetMaximumBin());
      float maxCountsErr = MinvLamMid_EPflow_tot[icent][ipt]->GetBinError(
          MinvLamMid_EPflow_tot[icent][ipt]->GetMaximumBin());

      float lowerEdgeY = -0.35 * (maxCounts + maxCountsErr);
      float upperEdgeY = 1.05 * (maxCounts + maxCountsErr);

      ax = gPad->DrawFrame(0.6, lowerEdgeY, 1.5, upperEdgeY); // no EPcorr
      SetAxis(ax, 3.0, 3.0);
      ax->SetXTitle("M_{inv} [GeV/c^{2}] ");
      ax->SetYTitle("Counts");
      ax->Draw();

      TLine *unitaty = new TLine(0.6, 0., 1.5, 0.);
      unitaty->Draw();

      MinvLamMid_EPflow_tot[icent][ipt]->Draw("ESAME");
      MinvLamMid_EPflow_tot[icent][ipt]->SetMarkerStyle(20);
      MinvLamMid_EPflow_tot[icent][ipt]->SetMarkerSize(1.0);
      MinvLamMid_EPflow_tot[icent][ipt]->SetMarkerColor(kBlack);

      MinvLamMid_EPflow_bg[icent][ipt]->Draw("ESAME");
      MinvLamMid_EPflow_bg[icent][ipt]->SetMarkerStyle(24);
      MinvLamMid_EPflow_bg[icent][ipt]->SetMarkerSize(1.0);
      MinvLamMid_EPflow_bg[icent][ipt]->SetMarkerColor(kBlue);

      lineL[ipt] = new TLine(IMcut[0], lowerEdgeY, IMcut[0], upperEdgeY);
      lineL[ipt]->SetLineColor(kBlack);
      lineL[ipt]->SetLineStyle(2);
      lineL[ipt]->SetLineWidth(2);
      lineL[ipt]->Draw("SAME");

      lineR[ipt] = new TLine(IMcut[1], lowerEdgeY, IMcut[1], upperEdgeY);
      lineR[ipt]->SetLineColor(kBlack);
      lineR[ipt]->SetLineStyle(2);
      lineR[ipt]->SetLineWidth(2);
      lineR[ipt]->Draw("SAME");

      TLatex *tex_sys =
          new TLatex(1.1, 0.85 * upperEdgeY,
                     Form("%s 3GeV Au+Au ", centName[icent].Data()));
      tex_sys->SetTextFont(42);
      tex_sys->SetTextSize(0.04);
      tex_sys->Draw("same");

      TLatex *tex_eta =
          new TLatex(1.1, 0.75 * upperEdgeY, Form("%s", etaNameMid[0].Data()));
      tex_eta->SetTextFont(42);
      tex_eta->SetTextSize(0.04);
      tex_eta->Draw("same");

      TLatex *tex_phi =
          new TLatex(1.1, 0.65 * upperEdgeY,
                     Form("%4.2f < p_{T} < %4.2f GeV/c", fctpTMidXaxis[ipt],
                          fctpTMidXaxis[ipt + 1]));
      tex_phi->SetTextFont(42);
      tex_phi->SetTextSize(0.04);
      tex_phi->Draw("same");

      ca_invMassMid[icent]->cd(nPtMid + ipt + 1);
      float minCounts, minCountsErr;
      float minTemp = 0., maxTemp = 0.;
      int invMBinIdx_Low =
          MinvLamMid_EPflow_sig[icent][ipt]->GetXaxis()->FindBin(lowFitRange);
      int invMBinIdx_High =
          MinvLamMid_EPflow_sig[icent][ipt]->GetXaxis()->FindBin(highFitRange);
      for (int iinvM = invMBinIdx_Low; iinvM < invMBinIdx_High; iinvM++) {
        maxCounts = MinvLamMid_EPflow_sig[icent][ipt]->GetBinContent(iinvM);
        maxCountsErr = MinvLamMid_EPflow_sig[icent][ipt]->GetBinError(iinvM);
        if ((maxCounts + maxCountsErr) > maxTemp)
          maxTemp = maxCounts + maxCountsErr;

        minCounts = MinvLamMid_EPflow_sig[icent][ipt]->GetBinContent(iinvM);
        minCountsErr = MinvLamMid_EPflow_sig[icent][ipt]->GetBinError(iinvM);
        if ((minCounts + minCountsErr) < minTemp)
          minTemp = minCounts + minCountsErr;
      }

      if (minTemp >= 0)
        lowerEdgeY = -0.15 * maxTemp;
      else
        lowerEdgeY = 1.15 * minTemp;
      if (maxTemp > 0)
        upperEdgeY = 1.15 * maxTemp;
      else
        upperEdgeY = maxTemp - minTemp;

      ax = gPad->DrawFrame(lowFitRange, lowerEdgeY, highFitRange,
                           upperEdgeY); // no EPcorr
      SetAxis(ax, 1.5, 1.8);
      ax->SetXTitle("M_{inv} [GeV/c^{2}] ");
      ax->SetYTitle("Counts");
      ax->Draw();

      MinvLamMid_EPflow_sig[icent][ipt]->Draw("E");
      MinvLamMid_EPflow_sig[icent][ipt]->GetXaxis()->SetRangeUser(lowFitRange,
                                                                  highFitRange);
      MinvLamMid_EPflow_sig[icent][ipt]->GetYaxis()->SetRangeUser(lowerEdgeY,
                                                                  upperEdgeY);
      MinvLamMid_EPflow_sig[icent][ipt]->GetXaxis()->SetTitle(
          "M_{inv} [GeV/c^{2}]");
      MinvLamMid_EPflow_sig[icent][ipt]->GetYaxis()->SetTitle("Counts");
      MinvLamMid_EPflow_sig[icent][ipt]->SetMarkerStyle(20);
      MinvLamMid_EPflow_sig[icent][ipt]->SetMarkerSize(0.5);
      MinvLamMid_EPflow_sig[icent][ipt]->SetMarkerColor(kRed);

      TLine *unitaty1 = new TLine(lowFitRange, 0., highFitRange, 0.);
      unitaty1->Draw();

      TLine *Kstar0mass = new TLine(initMassMean[icent], lowerEdgeY,
                                    initMassMean[icent], upperEdgeY);
      Kstar0mass->SetLineColor(kGreen);
      Kstar0mass->SetLineStyle(2);
      Kstar0mass->SetLineWidth(2);
      Kstar0mass->Draw("SAME");

      fMinvLamMid_EPflow_tot[icent][ipt]->Draw("SAME");
      fMinvLamMid_EPflow_sig[icent][ipt]->Draw("SAME");
      fMinvLamMid_EPflow_bg[icent][ipt]->Draw("SAME");

      MinvLamMid_EPflowReal_sig[icent][ipt] = SubtractBGFn(
          MinvLamMid_EPflow_sig[icent][ipt], fMinvLamMid_EPflow_bg[icent][ipt],
          lowFitRange, highFitRange);
      MinvLamMid_EPflowReal_sig[icent][ipt]->SetLineWidth(2.0);
      MinvLamMid_EPflowReal_sig[icent][ipt]->SetLineColor(kRed);
      MinvLamMid_EPflowReal_sig[icent][ipt]->Draw("HIST SAME");

      TLine *lineLFit = new TLine(IMcutMid_3Sigma[icent][ipt][0], lowerEdgeY,
                                  IMcutMid_3Sigma[icent][ipt][0], upperEdgeY);
      lineLFit->SetLineColor(kBlack);
      lineLFit->SetLineStyle(2);
      lineLFit->SetLineWidth(2);
      lineLFit->Draw("SAME");

      TLine *lineRFit = new TLine(IMcutMid_3Sigma[icent][ipt][1], lowerEdgeY,
                                  IMcutMid_3Sigma[icent][ipt][1], upperEdgeY);
      lineRFit->SetLineColor(kBlack);
      lineRFit->SetLineStyle(2);
      lineRFit->SetLineWidth(2);
      lineRFit->Draw("SAME");

      int lowerMassBin_Kstar0 = MinvLamMid_EPflowReal_sig[icent][ipt]->FindBin(
          IMcutMid_3Sigma[icent][ipt][0]);
      int upperMassBin_Kstar0 = MinvLamMid_EPflowReal_sig[icent][ipt]->FindBin(
                                    IMcutMid_3Sigma[icent][ipt][1]) -
                                1;
      double nSignalErr = 0., nSignal = 0.;
      nSignal = MinvLamMid_EPflowReal_sig[icent][ipt]->IntegralAndError(
          lowerMassBin_Kstar0, upperMassBin_Kstar0, nSignalErr);
      // nSignal =
      // MinvLamMid_EPflow_sig[icent][ipt]->IntegralAndError(lowerMassBin_Kstar0,
      // upperMassBin_Kstar0, nSignalErr); for(int i=lowerMassBin_Kstar0;
      // i<upperMassBin_Kstar0+1; i++) {
      //    nSignal += MinvLamMid_EPflowReal_sig[icent][ipt]->GetBinContent(i);
      //}

      // eff_weight = hEff_accReg_we_mid[icent]->GetBinContent(ipt+1);
      // if(eff_weight !=0) eff_weight = 1./ eff_weight;
      //            eff_weight = 1.0;
      eff_weight = 1. / EfficiencyMid[icent]->GetBinContent(1, ipt + 1);

      nXiYieldMid[icent][ipt] = nSignal * eff_weight;
      nXiYieldErrMid[icent][ipt] = nSignalErr * eff_weight; // sqrt(nSignal);

      lowerMassBin_Kstar0 = MinvLamMid_EPflow_tot[icent][ipt]->FindBin(
          IMcutMid_3Sigma[icent][ipt][0]);
      upperMassBin_Kstar0 = MinvLamMid_EPflow_tot[icent][ipt]->FindBin(
                                IMcutMid_3Sigma[icent][ipt][1]) -
                            1;
      float nTotal = MinvLamMid_EPflow_tot[icent][ipt]->Integral(
          lowerMassBin_Kstar0, upperMassBin_Kstar0);
      float nBackground = nTotal - nSignal;

      TLatex *tex_chi2 = new TLatex(0.91, 0.85 * upperEdgeY,
                                    Form("#chi^{2}/NDF = %4.2f / %d",
                                         FitQualityMid[icent][ipt][0],
                                         (int)FitQualityMid[icent][ipt][1]));
      tex_chi2->SetTextFont(42);
      tex_chi2->SetTextSize(0.04);
      tex_chi2->Draw("same");

      TLatex *tex_signif =
          new TLatex(0.91, 0.75 * upperEdgeY,
                     Form("S / #sqrt{S+B} = %6.2f", nSignal / sqrt(nTotal)));
      tex_signif->SetTextFont(42);
      tex_signif->SetTextSize(0.04);
      tex_signif->Draw("same");

      TLatex *tex_SBR = new TLatex(
          0.91, 0.65 * upperEdgeY,
          Form("S=%5.1f#pm%5.1f, B=%5.1f", nSignal, nSignalErr, nBackground));
      tex_SBR->SetTextFont(42);
      tex_SBR->SetTextSize(0.04);
      tex_SBR->Draw("same");

      TLatex *tex_mean =
          new TLatex(0.77, 0.45 * upperEdgeY,
                     Form("#mu=%5.4f#pm%5.4f", FitQualityMid[icent][ipt][2],
                          FitQualityMid[icent][ipt][3]));
      tex_mean->SetTextFont(42);
      tex_mean->SetTextSize(0.04);
      tex_mean->Draw("same");

      TLatex *tex_gamma =
          new TLatex(0.77, 0.35 * upperEdgeY,
                     Form("#Gamma=%5.4f#pm%5.4f", FitQualityMid[icent][ipt][4],
                          FitQualityMid[icent][ipt][5]));
      tex_gamma->SetTextFont(42);
      tex_gamma->SetTextSize(0.04);
      tex_gamma->Draw("same");
    }

    ca_invMassMid[icent]->cd();
    ca_invMassMid[icent]->Print(
        Form("../../files/plots_%s/"
             "invMassVsYCent%d_dNdYmid_%s_%s_EPmethod_fxt3GeV.png",
             cutName.Data(), icent, particleName.Data(), cutName.Data()));
  }

  //==========================================================================================================================================
  // get spectrum, merge all rapidity bins !!!
  // mid rapidity may have different pT bins

  float ptBinArr[nEta][nPt + 1];
  float mtBinArr[nEta][nPt + 1];
  for (int ieta = 0; ieta < nEta; ieta++) {
    for (int ipt = 0; ipt < nPt + 1; ipt++) {
      ptBinArr[ieta][ipt] = fctpTXaxis[ieta][ipt];
      mtBinArr[ieta][ipt] =
          sqrt(ptBinArr[ieta][ipt] * ptBinArr[ieta][ipt] + ximass * ximass) -
          ximass;
      cout << "mT bin edge is :::::" << mtBinArr[ieta][ipt] << "    ";
    }
    cout << endl;
  }
  // Fow rapidity
  float ptBinArrFow[nPtFow + 1];
  float mtBinArrFow[nPtFow + 1];
  for (int ipt = 0; ipt < nPtFow + 1; ipt++) {
    ptBinArrFow[ipt] = fctpTFowXaxis[ipt];
    mtBinArrFow[ipt] =
        sqrt(ptBinArrFow[ipt] * ptBinArrFow[ipt] + ximass * ximass) - ximass;
    cout << mtBinArrFow[ipt] << ", mtBinArrFow[ipt]    ";
  }
  cout << endl;
  // Mid rapidity
  float ptBinArrMid[nPtMid + 1];
  float mtBinArrMid[nPtMid + 1];
  for (int ipt = 0; ipt < nPtMid + 1; ipt++) {
    ptBinArrMid[ipt] = fctpTMidXaxis[ipt];
    mtBinArrMid[ipt] =
        sqrt(ptBinArrMid[ipt] * ptBinArrMid[ipt] + ximass * ximass) - ximass;
    cout << mtBinArrMid[ipt] << ", mtBinArrMid[ipt]    ";
  }
  cout << endl;

  TH1F *hXiYieldVsPt[nCentFlow][nTotEta]; // dN/pT/dpT/dy
  TH1F *hXiYieldVsMt[nCentFlow][nTotEta];
  TH1F *hXiYieldVsPt_raw[nCentFlow][nTotEta]; // dN/dpT
  TH1F *hXiYieldVsMt_raw[nCentFlow][nTotEta];
  TH1F *hXidNdpT[nTotEta][nCentFlow]; // dN/dpT/dy

  for (int icent = 0; icent < nCentFlow; icent++) {
    cout << "$$$$$$$$$$$$$$$$$$$$$$ " << centName[icent]
         << " $$$$$$$$$$$$$$$$$$$$$$" << endl;
    for (int ieta = 0; ieta < nTotEta; ieta++) {
      cout << "$$$$$$$$$$$$$$$$$$$$$$ " << etaNameTot[ieta]
           << " $$$$$$$$$$$$$$$$$$$$$$" << endl;

      if (ieta == 0) { //-1.0 < y < -0.8
        hXiYieldVsPt[icent][ieta] =
            new TH1F(Form("hXiYieldVsPt_cent%d_Y%d", icent, ieta), "", nPtFow,
                     ptBinArrFow);
        hXiYieldVsPt[icent][ieta]->Sumw2();
        hXiYieldVsMt[icent][ieta] =
            new TH1F(Form("hXiYieldVsMt_cent%d_Y%d", icent, ieta), "", nPtFow,
                     mtBinArrFow);
        hXiYieldVsMt[icent][ieta]->Sumw2();

        hXiYieldVsPt_raw[icent][ieta] =
            new TH1F(Form("hXiYieldVsPt_raw_cent%d_Y%d", icent, ieta), "",
                     nPtFow, ptBinArrFow);
        hXiYieldVsPt_raw[icent][ieta]->Sumw2();
        hXiYieldVsMt_raw[icent][ieta] =
            new TH1F(Form("hXiYieldVsMt_raw_cent%d_Y%d", icent, ieta), "",
                     nPtFow, mtBinArrFow);
        hXiYieldVsMt_raw[icent][ieta]->Sumw2();

        hXidNdpT[ieta][icent] = new TH1F(
            Form("hXidNdpT_Y%d_cent%d", ieta, icent), "", nPtFow, ptBinArrFow);
        hXidNdpT[ieta][icent]->Sumw2();

        for (int ipt = 0; ipt < nPtFow; ipt++) {
          hXiYieldVsPt_raw[icent][ieta]->SetBinContent(
              ipt + 1, nXiYieldFow[icent][ipt] /
                           (ptBinArrFow[ipt + 1] - ptBinArrFow[ipt])); // dN/dpT
          hXiYieldVsPt_raw[icent][ieta]->SetBinError(
              ipt + 1, nXiYieldErrFow[icent][ipt] /
                           (ptBinArrFow[ipt + 1] - ptBinArrFow[ipt]));

          float binCenterFow = 0.5 * (ptBinArrFow[ipt] + ptBinArrFow[ipt + 1]);
          float yieldPtFow =
              nXiYieldFow[icent][ipt] / br / 2 / rapiditywidthFow /
              TMath::Pi() / nEvents_cent[icent] / binCenterFow /
              (ptBinArrFow[ipt + 1] - ptBinArrFow[ipt]) * scaler4Show[ieta];
          float yieldPtErrFow =
              yieldPtFow * nXiYieldErrFow[icent][ipt] / nXiYieldFow[icent][ipt];
          hXiYieldVsPt[icent][ieta]->SetBinContent(
              ipt + 1, yieldPtFow); // 1/(B.R. N_{ev}2#piP_{T})dN^{2}/dP_{T}dy
                                    // (GeV/c)^{-2}
          hXiYieldVsPt[icent][ieta]->SetBinError(ipt + 1, yieldPtErrFow);

          hXiYieldVsMt_raw[icent][ieta]->SetBinContent(
              ipt + 1, nXiYieldFow[icent][ipt] /
                           (mtBinArrFow[ipt + 1] - mtBinArrFow[ipt]));
          hXiYieldVsMt_raw[icent][ieta]->SetBinError(
              ipt + 1, nXiYieldErrFow[icent][ipt] /
                           (mtBinArrFow[ipt + 1] - mtBinArrFow[ipt]));

          float binCenterMtFow =
              sqrt(binCenterFow * binCenterFow + ximass * ximass);
          float yieldMtFow =
              nXiYieldFow[icent][ipt] / br / 2 / rapiditywidthFow /
              TMath::Pi() / nEvents_cent[icent] / binCenterMtFow /
              (mtBinArrFow[ipt + 1] - mtBinArrFow[ipt]) * scaler4Show[ieta];
          float yieldMtErrFow =
              yieldMtFow * nXiYieldErrFow[icent][ipt] / nXiYieldFow[icent][ipt];
          hXiYieldVsMt[icent][ieta]->SetBinContent(
              ipt + 1, yieldMtFow); // 1/(B.R. N_{ev}2#pim_{T})dN^{2}/dm_{T}dy
                                    // (GeV/c)^{-2}
          hXiYieldVsMt[icent][ieta]->SetBinError(ipt + 1, yieldMtErrFow);

          hXidNdpT[ieta][icent]->SetBinContent(
              ipt + 1, nXiYieldFow[icent][ipt] /
                           (ptBinArrFow[ipt + 1] - ptBinArrFow[ipt]) /
                           rapiditywidthFow / nEvents_cent[icent] / 2. /
                           TMath::Pi() / br); // dN/dpT/dy
          hXidNdpT[ieta][icent]->SetBinError(
              ipt + 1, nXiYieldErrFow[icent][ipt] /
                           (ptBinArrFow[ipt + 1] - ptBinArrFow[ipt]) /
                           rapiditywidthFow / nEvents_cent[icent] / 2. /
                           TMath::Pi() / br);

          cout << "Xi counts at pT bin " << ipt << ": "
               << nXiYieldFow[icent][ipt] << "+-" << nXiYieldErrFow[icent][ipt]
               << endl;
          cout << "Xi counts at pT bin " << ipt
               << " after weight: " << yieldMtFow << "+-" << yieldMtErrFow
               << endl;
        }
      } else if (ieta == nTotEta - 1) { //-0.2 < y < 0.2
        hXiYieldVsPt[icent][ieta] =
            new TH1F(Form("hXiYieldVsPt_cent%d_Y%d", icent, ieta), "", nPtMid,
                     ptBinArrMid);
        hXiYieldVsPt[icent][ieta]->Sumw2();
        hXiYieldVsMt[icent][ieta] =
            new TH1F(Form("hXiYieldVsMt_cent%d_Y%d", icent, ieta), "", nPtMid,
                     mtBinArrMid);
        hXiYieldVsMt[icent][ieta]->Sumw2();

        hXiYieldVsPt_raw[icent][ieta] =
            new TH1F(Form("hXiYieldVsPt_raw_cent%d_Y%d", icent, ieta), "",
                     nPtMid, ptBinArrMid);
        hXiYieldVsPt_raw[icent][ieta]->Sumw2();
        hXiYieldVsMt_raw[icent][ieta] =
            new TH1F(Form("hXiYieldVsMt_raw_cent%d_Y%d", icent, ieta), "",
                     nPtMid, mtBinArrMid);
        hXiYieldVsMt_raw[icent][ieta]->Sumw2();

        hXidNdpT[ieta][icent] = new TH1F(
            Form("hXidNdpT_Y%d_cent%d", ieta, icent), "", nPtMid, ptBinArrMid);
        hXidNdpT[ieta][icent]->Sumw2();

        for (int ipt = 0; ipt < nPtMid; ipt++) {
          hXiYieldVsPt_raw[icent][ieta]->SetBinContent(
              ipt + 1, nXiYieldMid[icent][ipt] /
                           (ptBinArrMid[ipt + 1] - ptBinArrMid[ipt])); // dN/dpT
          hXiYieldVsPt_raw[icent][ieta]->SetBinError(
              ipt + 1, nXiYieldErrMid[icent][ipt] /
                           (ptBinArrMid[ipt + 1] - ptBinArrMid[ipt]));

          float binCenterMid = 0.5 * (ptBinArrMid[ipt] + ptBinArrMid[ipt + 1]);
          float yieldPtMid =
              nXiYieldMid[icent][ipt] / br / 2 / rapiditywidthMid /
              TMath::Pi() / nEvents_cent[icent] / binCenterMid /
              (ptBinArrMid[ipt + 1] - ptBinArrMid[ipt]) * scaler4Show[ieta];
          float yieldPtErrMid =
              yieldPtMid * nXiYieldErrMid[icent][ipt] / nXiYieldMid[icent][ipt];
          hXiYieldVsPt[icent][ieta]->SetBinContent(
              ipt + 1, yieldPtMid); // 1/(B.R. N_{ev}2#piP_{T})dN^{2}/dP_{T}dy
                                    // (GeV/c)^{-2}
          hXiYieldVsPt[icent][ieta]->SetBinError(ipt + 1, yieldPtErrMid);

          hXiYieldVsMt_raw[icent][ieta]->SetBinContent(
              ipt + 1, nXiYieldMid[icent][ipt] /
                           (mtBinArrMid[ipt + 1] - mtBinArrMid[ipt]));
          hXiYieldVsMt_raw[icent][ieta]->SetBinError(
              ipt + 1, nXiYieldErrMid[icent][ipt] /
                           (mtBinArrMid[ipt + 1] - mtBinArrMid[ipt]));

          float binCenterMtMid =
              sqrt(binCenterMid * binCenterMid + ximass * ximass);
          float yieldMtMid =
              nXiYieldMid[icent][ipt] / br / 2 / rapiditywidthMid /
              TMath::Pi() / nEvents_cent[icent] / binCenterMtMid /
              (mtBinArrMid[ipt + 1] - mtBinArrMid[ipt]) * scaler4Show[ieta];
          float yieldMtErrMid =
              yieldMtMid * nXiYieldErrMid[icent][ipt] / nXiYieldMid[icent][ipt];
          hXiYieldVsMt[icent][ieta]->SetBinContent(
              ipt + 1, yieldMtMid); // 1/(B.R. N_{ev}2#pim_{T})dN^{2}/dm_{T}dy
                                    // (GeV/c)^{-2}
          hXiYieldVsMt[icent][ieta]->SetBinError(ipt + 1, yieldMtErrMid);

          hXidNdpT[ieta][icent]->SetBinContent(
              ipt + 1, nXiYieldMid[icent][ipt] /
                           (ptBinArrMid[ipt + 1] - ptBinArrMid[ipt]) /
                           rapiditywidthMid / nEvents_cent[icent] / 2. /
                           TMath::Pi() / br); // dN/dpT/dy
          hXidNdpT[ieta][icent]->SetBinError(
              ipt + 1, nXiYieldErrMid[icent][ipt] /
                           (ptBinArrMid[ipt + 1] - ptBinArrMid[ipt]) /
                           rapiditywidthMid / nEvents_cent[icent] / 2. /
                           TMath::Pi() / br);

          cout << "Xi counts at pT bin " << ipt << ": "
               << nXiYieldMid[icent][ipt] << "+-" << nXiYieldErrMid[icent][ipt]
               << endl;
          cout << "Xi counts at pT bin " << ipt
               << " after weight: " << yieldMtMid << "+-" << yieldMtErrMid
               << endl;
        }
      } else {
        int nieta = ieta - 1;
        hXiYieldVsPt[icent][ieta] =
            new TH1F(Form("hXiYieldVsPt_cent%d_Y%d", icent, ieta), "", nPt,
                     ptBinArr[nieta]);
        hXiYieldVsPt[icent][ieta]->Sumw2();
        hXiYieldVsMt[icent][ieta] =
            new TH1F(Form("hXiYieldVsMt_cent%d_Y%d", icent, ieta), "", nPt,
                     mtBinArr[nieta]);
        hXiYieldVsMt[icent][ieta]->Sumw2();

        hXiYieldVsPt_raw[icent][ieta] =
            new TH1F(Form("hXiYieldVsPt_raw_cent%d_Y%d", icent, ieta), "", nPt,
                     ptBinArr[nieta]);
        hXiYieldVsPt_raw[icent][ieta]->Sumw2();
        hXiYieldVsMt_raw[icent][ieta] =
            new TH1F(Form("hXiYieldVsMt_raw_cent%d_Y%d", icent, ieta), "", nPt,
                     mtBinArr[nieta]);
        hXiYieldVsMt_raw[icent][ieta]->Sumw2();

        hXidNdpT[ieta][icent] = new TH1F(
            Form("hXidNdpT_Y%d_cent%d", ieta, icent), "", nPt, ptBinArr[nieta]);
        hXidNdpT[ieta][icent]->Sumw2();

        for (int ipt = 0; ipt < nPt; ipt++) {
          hXiYieldVsPt_raw[icent][ieta]->SetBinContent(
              ipt + 1,
              nXiYield[icent][nieta][ipt] /
                  (ptBinArr[nieta][ipt + 1] - ptBinArr[nieta][ipt])); // dN/dpT
          hXiYieldVsPt_raw[icent][ieta]->SetBinError(
              ipt + 1, nXiYieldErr[icent][nieta][ipt] /
                           (ptBinArr[nieta][ipt + 1] - ptBinArr[nieta][ipt]));

          float binCenter =
              0.5 * (ptBinArr[nieta][ipt] + ptBinArr[nieta][ipt + 1]);
          float yieldPt = nXiYield[icent][nieta][ipt] / br / 2 / rapiditywidth /
                          TMath::Pi() / nEvents_cent[icent] / binCenter /
                          (ptBinArr[nieta][ipt + 1] - ptBinArr[nieta][ipt]) *
                          scaler4Show[ieta];
          float yieldPtErr = yieldPt * nXiYieldErr[icent][nieta][ipt] /
                             nXiYield[icent][nieta][ipt];
          hXiYieldVsPt[icent][ieta]->SetBinContent(
              ipt + 1,
              yieldPt); // 1/(B.R. N_{ev}2#piP_{T})dN^{2}/dP_{T}dy (GeV/c)^{-2}
          hXiYieldVsPt[icent][ieta]->SetBinError(ipt + 1, yieldPtErr);

          hXiYieldVsMt_raw[icent][ieta]->SetBinContent(
              ipt + 1,
              nXiYield[icent][nieta][ipt] /
                  (mtBinArr[nieta][ipt + 1] - mtBinArr[nieta][ipt])); // dN/dmT
          hXiYieldVsMt_raw[icent][ieta]->SetBinError(
              ipt + 1, nXiYieldErr[icent][nieta][ipt] /
                           (mtBinArr[nieta][ipt + 1] - mtBinArr[nieta][ipt]));

          float binCenterMt = sqrt(binCenter * binCenter + ximass * ximass);
          float yieldMt = nXiYield[icent][nieta][ipt] / br / 2 / rapiditywidth /
                          TMath::Pi() / nEvents_cent[icent] / binCenterMt /
                          (mtBinArr[nieta][ipt + 1] - mtBinArr[nieta][ipt]) *
                          scaler4Show[ieta];
          cout << "mT and pT, check:     " << yieldMt << " : " << yieldPt
               << endl;
          cout << "check:  " << nXiYield[icent][nieta][ipt] << " : "
               << nEvents_cent[icent] << endl;
          float yieldMtErr = yieldMt * nXiYieldErr[icent][nieta][ipt] /
                             nXiYield[icent][nieta][ipt];
          hXiYieldVsMt[icent][ieta]->SetBinContent(
              ipt + 1,
              yieldPt); // 1/(B.R. N_{ev}2#pim_{T})dN^{2}/dm_{T}dy (GeV/c)^{-2}
          hXiYieldVsMt[icent][ieta]->SetBinError(ipt + 1, yieldMtErr);

          hXidNdpT[ieta][icent]->SetBinContent(
              ipt + 1, nXiYield[icent][nieta][ipt] /
                           (ptBinArr[nieta][ipt + 1] - ptBinArr[nieta][ipt]) /
                           rapiditywidth / nEvents_cent[icent] / 2. /
                           TMath::Pi() / br); // dN/dpTdy
          if (debug) {
            std::cout << "zla icent == " << icent << " ieta == " << nieta
                      << " ipt == " << ipt << " bincontent == "
                      << nXiYield[icent][nieta][ipt] / rapiditywidth /
                             nEvents_cent[icent] / br
                      << "and nXiYield == " << nXiYield[icent][nieta][ipt]
                      << std::endl;
          }
          hXidNdpT[ieta][icent]->SetBinError(
              ipt + 1, nXiYieldErr[icent][nieta][ipt] /
                           (ptBinArr[nieta][ipt + 1] - ptBinArr[nieta][ipt]) /
                           rapiditywidth / nEvents_cent[icent] / 2. /
                           TMath::Pi() / br);
          cout << "Xi counts at pT bin " << ipt << " : "
               << nXiYield[icent][nieta][ipt] << "+-"
               << nXiYieldErr[icent][nieta][ipt] << endl;
          cout << "Xi counts at pT bin " << ipt << " after weight: " << yieldPt
               << "+-" << yieldPtErr << endl;
        }
      }
    }
  }
  //	return;
  //==========================================================================================================================================
  // fitting using expmT, extract pt integrated dNdy

  TF1 *f1_exp;
  double xiYield_int[nCentFlow][nTotEta];
  double xiYieldErr_int[nCentFlow][nTotEta];

  double *par;
  const double *parErr;
  double temperature[nCentFlow][nTotEta];
  double temperatureErr[nCentFlow][nTotEta];

  TF1 *ptmptexp[nCentFlow][nTotEta];
  TF1 *ptdptexp[nCentFlow][nTotEta];
  TF1 *ptexp[nCentFlow][nTotEta]; // just for plots

  TCanvas *ca_dNdpTdy = new TCanvas("ca_dNdpTdy", "", 450 * nCentFlow, 400);
  ca_dNdpTdy->Divide(nCentFlow, 1);

  for (int icent = 0; icent < nCentFlow; icent++) {
    ca_dNdpTdy->cd(icent + 1)->SetLogy();
    ax = gPad->DrawFrame(
        0,
        1e-2 * hXidNdpT[nTotEta - 1][icent]->GetBinContent(
                   hXidNdpT[nTotEta - 1][icent]->GetMinimumBin()),
        2.5,
        2 * hXidNdpT[nTotEta - 1][icent]->GetBinContent(
                hXidNdpT[nTotEta - 1][icent]->GetMaximumBin()));
    SetAxis(ax, 1, 1);
    ax->SetYTitle("K^{*0} dN/dp_{T}/dy");
    ax->SetXTitle("p_{T} [GeV/c]");
    ax->Draw();

    TLegend *legT = new TLegend(0.25, 0.2, 0.45, 0.6);
    legT->SetFillColor(10);
    legT->SetLineStyle(3004);
    legT->SetLineColor(10);
    legT->SetLineWidth(0.);
    legT->SetTextFont(42);
    legT->SetTextSize(0.03);
    legT->SetHeader(centName[icent].Data());

    for (int ieta = 0; ieta < nTotEta; ieta++) {
      f1_exp = new TF1("f1_exp", fitFuncExp, 0, 2.5, 2); //[0]*exp(-x/[1])*x
      f1_exp->SetParNames("Scale1", "Temperature");
      f1_exp->SetParameters(1.5e-03, 0.12);
      f1_exp->SetLineColor(colorIndex[ieta]);
      f1_exp->SetLineWidth(2);
      cout << "**** fit parameter **** "
           << " : " << etaNameTot[ieta].Data() << endl;
      hXidNdpT[ieta][icent]->Fit("f1_exp", "IRMN", "",
                                 RapiditypTLow[ieta][icent],
                                 RapiditypTHigh[ieta][icent]);
      f1_exp->Draw("SAME");

      hXidNdpT[ieta][icent]->SetYTitle("K^{*0} dN/dp_{T}/dy");
      hXidNdpT[ieta][icent]->SetXTitle("p_{T} [GeV/c]");

      par = f1_exp->GetParameters();
      parErr = f1_exp->GetParErrors();
      temperature[icent][ieta] = par[1];
      temperatureErr[icent][ieta] = parErr[1];
      cout << "chcking the fit parameter T in eta: " << etaNameTot[ieta].Data()
           << " : icent: " << icent << " :" << par[1]
           << "  ; error: " << parErr[1] << endl;

      ptmptexp[icent][ieta] =
          new TF1(Form("ptmptexp_cent%d_y%d", icent, ieta), "[0]*exp(-x/[1])",
                  0, sqrt(2.5 * 2.5 + ximass * ximass) - ximass);
      ptmptexp[icent][ieta]->SetParameters(
          par[0] * scaler4Show[ieta], par[1],
          ximass); // here 0.2 represents rapidity width

      ptdptexp[icent][ieta] = new TF1(
          Form("ptdptexp_cent%d_y%d", icent, ieta),
          "[0]*exp(-(sqrt(x*x + 0.89594*0.89594)-0.89594)/[1])", 0, 2.5);
      ptdptexp[icent][ieta]->SetParameters(
          par[0] * scaler4Show[ieta], par[1],
          ximass); // here 0.2 represents rapidity width

      ptexp[icent][ieta] =
          new TF1(Form("ptexp_cent%d_y%d", icent, ieta), fitFuncExp, 0, 2.5, 2);
      //      ptexp[icent][ieta]->SetParameters(par[0]*scaler4Show[ieta],
      //      par[1], ximass);
      ptexp[icent][ieta]->SetParameters(par[0], par[1], ximass);
      ptexp[icent][ieta]->SetChisquare(f1_exp->GetChisquare());
      ptexp[icent][ieta]->SetNDF(f1_exp->GetNDF());

      float dNdy = 0.;
      float dNdyErr = 0.;
      if (ieta == 0) {
        for (int ipt = 0; ipt < nPtFow; ipt++) {
          float ptbinwidth = ptBinArrFow[ipt + 1] - ptBinArrFow[ipt];
          dNdy += hXidNdpT[ieta][icent]->GetBinContent(ipt + 1) * ptbinwidth *
                  2 * TMath::Pi();
          dNdyErr += pow(hXidNdpT[ieta][icent]->GetBinError(ipt + 1), 2) *
                     pow(ptbinwidth, 2) * pow(2 * TMath::Pi(), 2);
        }
      } else if (ieta == nTotEta - 1) {
        for (int ipt = 0; ipt < nPtMid; ipt++) {
          float ptbinwidth = ptBinArrMid[ipt + 1] - ptBinArrMid[ipt];
          dNdy += hXidNdpT[ieta][icent]->GetBinContent(ipt + 1) * ptbinwidth *
                  2 * TMath::Pi();
          dNdyErr += pow(hXidNdpT[ieta][icent]->GetBinError(ipt + 1), 2) *
                     pow(ptbinwidth, 2) * pow(2 * TMath::Pi(), 2);
        }
      } else {
        int nieta = ieta - 1;
        for (int ipt = 0; ipt < nPt; ipt++) {
          float ptbinwidth = ptBinArr[nieta][ipt + 1] - ptBinArr[nieta][ipt];
          dNdy += hXidNdpT[ieta][icent]->GetBinContent(ipt + 1) * ptbinwidth *
                  2 * TMath::Pi();
          dNdyErr += pow(hXidNdpT[ieta][icent]->GetBinError(ipt + 1), 2) *
                     pow(ptbinwidth, 2) * pow(2 * TMath::Pi(), 2);
        }
      }
      cout << "dNdpTdy $$$$$$$$$$$$$$$$$$$ dNdy : " << dNdy << "+-"
           << sqrt(dNdyErr) << endl;

      xiYield_int[icent][ieta] =
          dNdy + (f1_exp->Integral(0, RapiditypTLow[ieta][icent]) +
                  f1_exp->Integral(RapiditypTHigh[ieta][icent], 2.5)) *
                     2 * TMath::Pi();
      if (debug) {
        std::cout << "zla 2 icent == " << icent << " ieta == " << ieta - 1
                  << " dNdy == " << dNdy
                  << " xiYield_int == " << xiYield_int[icent][ieta];
      }
      xiYieldErr_int[icent][ieta] =
          xiYield_int[icent][ieta] / dNdy * sqrt(dNdyErr);
      xiYield_int[icent][ieta] = xiYield_int[icent][ieta];
      cout << "total xiYield_int: " << xiYield_int[icent][ieta] << "+-"
           << xiYieldErr_int[icent][ieta] << endl;

      hXidNdpT[ieta][icent]->Draw("ESAME");
      hXidNdpT[ieta][icent]->SetMarkerStyle(24);
      hXidNdpT[ieta][icent]->SetMarkerSize(1.0);
      hXidNdpT[ieta][icent]->SetMarkerColor(colorIndex[ieta]);
      hXidNdpT[ieta][icent]->SetLineColor(colorIndex[ieta]);
      f1_exp->Draw("SAME");
      f1_exp->SetLineColor(colorIndex[ieta]);
      legT->AddEntry(hXidNdpT[ieta][icent], Form("%s", etaNameTot[ieta].Data()),
                     "lp");
    }
    legT->Draw("same");
  }
  ca_dNdpTdy->cd();
  ca_dNdpTdy->Print(
      Form("../../files/plots_%s/dNdpTdy_%s_%s_EPmethod_fxt3GeV.png",
           cutName.Data(), particleName.Data(), cutName.Data()));
  //================================================================================================
  // fit to dNdpT, to check if it's consistent with fitting to dNdpTdy !!!

  TCanvas *ca_dNdpT = new TCanvas("ca_dNdpT", "", 450 * nCentFlow, 400);
  ca_dNdpT->Divide(nCentFlow, 1);

  for (int icent = 0; icent < nCentFlow; icent++) {
    ca_dNdpT->cd(icent + 1)->SetLogy();

    ax = gPad->DrawFrame(0, 1, 2.5, 2.0e5);
    SetAxis(ax, 1, 1);
    ax->SetYTitle("K^{*0} dN/dp_{T}");
    ax->SetXTitle("p_{T} [GeV/c]");
    ax->Draw();

    TLegend *legT = new TLegend(0.2, 0.2, 0.4, 0.6);
    legT->SetFillColor(10);
    legT->SetLineStyle(3004);
    legT->SetLineColor(10);
    legT->SetLineWidth(0.);
    legT->SetTextFont(42);
    legT->SetTextSize(0.03);
    legT->SetHeader(centName[icent].Data());

    for (int ieta = 0; ieta < nTotEta; ieta++) {

      f1_exp = new TF1("f1_exp", fitFuncExp, 0, 2.5, 2); //[0]*exp(-x/[1])*x
      f1_exp->SetParNames("Scale1", "Temperature");
      f1_exp->SetParameter(0, 1.2 * 1.0e5);
      f1_exp->SetParameter(1, 0.5);
      if (ieta == 0) {
        f1_exp->SetParameter(0, 8.94775e+04);
        f1_exp->SetParameter(1, 1.23059e-01);
      }
      f1_exp->SetLineColor(colorIndex[ieta]);
      f1_exp->SetLineWidth(2);
      cout << "**** fit parameter **** " << endl;
      hXiYieldVsPt_raw[icent][ieta]->Fit("f1_exp", "IRMN", "",
                                         RapiditypTLow[ieta][icent],
                                         RapiditypTHigh[ieta][icent]);
      f1_exp->Draw("SAME");

      par = f1_exp->GetParameters();
      parErr = f1_exp->GetParErrors();
      cout << "chcking the fit parameter T in eta: " << etaNameTot[ieta].Data()
           << " : icent: " << icent << " :" << par[1]
           << "  ; error: " << parErr[1] << endl;

      float dNdy = 0.;
      float dNdyErr = 0.;
      float dNdyTot = 0;
      float dNdyErrTot = 0;

      double rapiditywidthTot;
      if (ieta == 0) {
        rapiditywidthTot = rapiditywidthFow;
        for (int ipt = 0; ipt < nPtFow; ipt++) {
          float ptbinwidth = ptBinArrFow[ipt + 1] - ptBinArrFow[ipt];
          dNdy += hXiYieldVsPt_raw[icent][ieta]->GetBinContent(ipt + 1) *
                  ptbinwidth / rapiditywidthTot;
          dNdyErr +=
              pow(hXiYieldVsPt_raw[icent][ieta]->GetBinError(ipt + 1), 2) *
              pow(ptbinwidth, 2) / pow(rapiditywidthTot, 2);
        }
      } else if (ieta == nTotEta - 1) {
        rapiditywidthTot = rapiditywidthMid;
        for (int ipt = 0; ipt < nPtMid; ipt++) {
          float ptbinwidth = ptBinArrMid[ipt + 1] - ptBinArrMid[ipt];
          dNdy += hXiYieldVsPt_raw[icent][ieta]->GetBinContent(ipt + 1) *
                  ptbinwidth / rapiditywidthTot;
          dNdyErr +=
              pow(hXiYieldVsPt_raw[icent][ieta]->GetBinError(ipt + 1), 2) *
              pow(ptbinwidth, 2) / pow(rapiditywidthTot, 2);
        }
      } else {
        rapiditywidthTot = rapiditywidth;
        int nieta = ieta - 1;
        for (int ipt = 0; ipt < nPt; ipt++) {
          float ptbinwidth = ptBinArr[nieta][ipt + 1] - ptBinArr[nieta][ipt];
          dNdy += hXiYieldVsPt_raw[icent][ieta]->GetBinContent(ipt + 1) *
                  ptbinwidth / rapiditywidthTot;
          dNdyErr +=
              pow(hXiYieldVsPt_raw[icent][ieta]->GetBinError(ipt + 1), 2) *
              pow(ptbinwidth, 2) / pow(rapiditywidthTot, 2);
        }
      }
      cout << "dNdpT $$$$$$$$$$$$$$$$$$$ dNdy : " << dNdy << "+-"
           << sqrt(dNdyErr) << endl;

      dNdyTot = dNdy + (f1_exp->Integral(0, RapiditypTLow[ieta][icent]) +
                        f1_exp->Integral(RapiditypTHigh[ieta][icent], 2.5)) /
                           rapiditywidthTot;
      dNdyErrTot =
          dNdyTot / dNdy * sqrt(dNdyErr) / nEvents_cent[icent] / br; // zyj
      dNdyTot = dNdyTot / nEvents_cent[icent] / br;
      cout << "xiYield_int total: " << dNdyTot << "+-" << dNdyErrTot << endl;
      cout << "================================================================"
              "========"
           << endl;

      hXiYieldVsPt_raw[icent][ieta]->Draw("ESAME");
      hXiYieldVsPt_raw[icent][ieta]->SetMarkerStyle(24);
      hXiYieldVsPt_raw[icent][ieta]->SetMarkerSize(1.0);
      hXiYieldVsPt_raw[icent][ieta]->SetMarkerColor(colorIndex[ieta]);
      hXiYieldVsPt_raw[icent][ieta]->SetLineColor(colorIndex[ieta]);
      f1_exp->Draw("SAME");
      f1_exp->SetLineColor(colorIndex[ieta]);
      legT->AddEntry(hXiYieldVsPt_raw[icent][ieta],
                     Form("%s", etaNameTot[ieta].Data()), "lp");
    }
    legT->Draw("same");
  }

  ca_dNdpT->cd();
  ca_dNdpT->Print(Form("../../files/plots_%s/dNdpT_%s_%s_EPmethod_fxt3GeV.png",
                       cutName.Data(), particleName.Data(), cutName.Data()));
  //=================================================================================================================

  TCanvas *ca_dNdy =
      new TCanvas("ca_dNdy", "K^{*0} dNdy", 450 * nCentFlow, 800);
  ca_dNdy->Divide(nCentFlow, 2);

  for (int icent = 0; icent < nCentFlow; icent++) {
    ca_dNdy->cd(icent + 1)->SetLogy();

    //    for(int ieta=0; ieta<nTotEta; ieta++)
    //    hXidNdpT[ieta][icent]->Scale(scaler4Show[ieta]);

    ax = gPad->DrawFrame(0, 1.0e-8, 2.5,
                         hXidNdpT[nTotEta - 1][icent]->GetBinContent(1) * 5);
    SetAxis(ax, 2.5, 2.5);
    ax->SetYTitle("d^{2}N/(B.R. 2#piN_{evt}p_{T}dy) [(GeV/c)^{-2}]");
    ax->SetXTitle("p_{T} [GeV/c]");
    ax->Draw();

    //    TLegend *legT = new TLegend(0.75, 0.2, 0.9, 0.5);
    TLegend *legT = new TLegend(0.2, 0.2, 0.35, 0.5);
    legT->SetFillColor(10);
    legT->SetLineStyle(3004);
    legT->SetLineColor(10);
    legT->SetLineWidth(0.);
    legT->SetTextFont(42);
    legT->SetTextSize(0.03);
    legT->SetHeader(centName[icent].Data());
    for (int ieta = 0; ieta < nTotEta; ieta++) {
      //      if(ieta==nTotEta-2) continue;

      hXidNdpT[ieta][icent]->Draw("ESAME");
      hXidNdpT[ieta][icent]->SetMarkerStyle(24);
      hXidNdpT[ieta][icent]->SetMarkerSize(1.0);
      hXidNdpT[ieta][icent]->SetMarkerColor(colorIndex[ieta]);
      hXidNdpT[ieta][icent]->SetLineColor(colorIndex[ieta]);
      //      legT->AddEntry(hXidNdpT[ieta][icent], Form("%s",
      //      etaNameTotScale[ieta].Data()), "lp");
      legT->AddEntry(hXidNdpT[ieta][icent],
                     Form("%s,#chi^{2}/ndf=%2.3f", etaNameTotScale[ieta].Data(),
                          ptexp[icent][ieta]->GetChisquare() /
                              ptexp[icent][ieta]->GetNDF()),
                     "lp");

      ptexp[icent][ieta]->SetLineColor(colorIndex[ieta]);
      ptexp[icent][ieta]->Draw("same");
    }
    legT->Draw();
  }
  //=================================================================================================================
  for (int icent = 0; icent < nCentFlow; icent++) {
    ca_dNdy->cd(icent + 1 + nCentFlow)->SetLogy();
    ax =
        gPad->DrawFrame(0, 1.0e-8, sqrt(2.5 * 2.5 + ximass * ximass) - ximass,
                        hXiYieldVsMt[icent][nTotEta - 1]->GetBinContent(1) * 5);
    //    if(icent==0) {
    //      ca_dNdy->cd(icent+2)->SetLogy();
    //      ax = gPad->DrawFrame( 0, 1.0e-9,
    //      sqrt(2.3*2.3+ximass*ximass)-ximass, 2.0e-1 );
    //    }
    //    else {
    //      ca_dNdy->cd(icent+3)->SetLogy();
    //      ax = gPad->DrawFrame( 0, 1.0e-9,
    //      sqrt(2.3*2.3+ximass*ximass)-ximass, 1.0e-1 );
    //    }
    //
    SetAxis(ax, 2.5, 2.5);
    // ax->SetLabelSize(0.1,"XY");
    ax->SetYTitle("d^{2}N/(B.R. 2#piN_{evt}m_{T}dydm_{T}) [(GeV/c)^{-2}]");
    ax->SetXTitle("m_{T} - m_{0} [GeV/c^{2}]");
    ax->Draw();
    TLegend *legT = new TLegend(0.2, 0.2, 0.35, 0.5);
    legT->SetFillColor(10);
    legT->SetLineStyle(3004);
    legT->SetLineColor(10);
    legT->SetLineWidth(0.);
    legT->SetTextFont(42);
    legT->SetTextSize(0.03);
    legT->SetHeader(centName[icent].Data());

    for (int ieta = 0; ieta < nTotEta; ieta++) {
      hXiYieldVsMt[icent][ieta]->Draw("ESAME");
      hXiYieldVsMt[icent][ieta]->SetMarkerStyle(24);
      hXiYieldVsMt[icent][ieta]->SetMarkerSize(1.0);
      hXiYieldVsMt[icent][ieta]->SetMarkerColor(colorIndex[ieta]);
      hXiYieldVsMt[icent][ieta]->SetLineColor(colorIndex[ieta]);
      ptmptexp[icent][ieta]->Draw("SAME");
      ptmptexp[icent][ieta]->SetLineColor(colorIndex[ieta]);

      legT->AddEntry(
          hXiYieldVsMt[icent][ieta],
          Form("%s: T=%4.3f#pm%4.3f GeV", etaNameTotScale[ieta].Data(),
               temperature[icent][ieta], temperatureErr[icent][ieta]),
          "lp");
    }
    legT->Draw();

    TLegend *leg = new TLegend(0.6, 0.8, 0.85, 0.92);
    leg->SetFillColor(10);
    leg->SetLineStyle(3004);
    leg->SetLineColor(10);
    leg->SetLineWidth(0.);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->SetHeader("K^{*0} @ 3 GeV Au+Au collision");
    leg->AddEntry(f1_exp, "m_{T} exponential fit", "lp");
    leg->Draw();
  }

  ca_dNdy->cd();
  ca_dNdy->Print(Form("../../files/plots_%s/dNdmT_%s_%s_EPmethod_fxt3GeV.png",
                      cutName.Data(), particleName.Data(), cutName.Data()));
  //================================================================================================
  // Temperature vs rapidity

  TH1F *hTeff = new TH1F("hTeff", "hTeff", nCentFlow, 0.5,
                         nCentFlow + 0.5); // save the Teff result
  TH1F *hdTdy[nCentFlow];
  TH1F *hdTdy_rf[nCentFlow];
  TH1F *hdTdy_Mid[nCentFlow];
  for (int icent = 0; icent < nCentFlow; icent++) {
    //    hdTdy[icent] = new TH1F(Form("hdTdy_icent%d", icent), "", 5, -0.8,
    //    0.2);
    hdTdy[icent] = new TH1F(Form("hdTdy_icent%d", icent), "", 8, -0.8, 0.8);
    hdTdy_rf[icent] = new TH1F(Form("hdTdy_rf_icent%d", icent), "", 4, 0, 0.8);
    hdTdy_Mid[icent] =
        new TH1F(Form("hdTdy_Mid_icent%d", icent), "", 1, fctetaXaxisMid);
  }
  for (int icent = 0; icent < nCentFlow; icent++) {
    for (int ieta = 0; ieta < 4; ieta++) {
      hdTdy[icent]->SetBinContent(ieta + 1, temperature[icent][ieta]);
      hdTdy[icent]->SetBinError(ieta + 1, temperatureErr[icent][ieta]);
    }
    for (int ieta = 0; ieta < 4; ieta++) {
      hdTdy_rf[icent]->SetBinContent(ieta + 1, temperature[icent][3 - ieta]);
      hdTdy_rf[icent]->SetBinError(ieta + 1, temperatureErr[icent][3 - ieta]);

      hdTdy[icent]->SetBinContent(ieta + 5, temperature[icent][3 - ieta]);
      hdTdy[icent]->SetBinError(ieta + 5, temperatureErr[icent][3 - ieta]);
    }
    hdTdy_Mid[icent]->SetBinContent(1, temperature[icent][nTotEta - 1]);
    hdTdy_Mid[icent]->SetBinError(1, temperatureErr[icent][nTotEta - 1]);
  }

  TF1 *f1_coshy_cent[nCentFlow];

  TCanvas *ca_T = new TCanvas("ca_T", "Xi dTdy", 450 * nCentFlow, 400);
  ca_T->Divide(nCentFlow, 1);
  for (int icent = 0; icent < nCentFlow; icent++) {
    ca_T->cd(icent + 1);
    ax = gPad->DrawFrame(-1.0, 0, 1.0, hdTdy[icent]->GetBinContent(4) * 1.7);
    SetAxis(ax, 1, 1);
    ax->SetYTitle("T [GeV]");
    ax->SetXTitle("y");
    ax->Draw();

    hdTdy[icent]->Draw("ESAME");
    hdTdy[icent]->SetMarkerStyle(20);
    hdTdy[icent]->SetMarkerSize(1.0);
    hdTdy[icent]->SetMarkerColor(kBlue);
    hdTdy[icent]->SetLineColor(kBlue);

    hdTdy_rf[icent]->Draw("ESAME");
    hdTdy_rf[icent]->SetMarkerStyle(24);
    hdTdy_rf[icent]->SetMarkerSize(1.0);
    hdTdy_rf[icent]->SetMarkerColor(kBlue);
    hdTdy_rf[icent]->SetLineColor(kBlue);

    hdTdy_Mid[icent]->Draw("ESAME");
    hdTdy_Mid[icent]->SetMarkerStyle(24);
    hdTdy_Mid[icent]->SetMarkerSize(1.0);
    hdTdy_Mid[icent]->SetMarkerColor(kRed);
    hdTdy_Mid[icent]->SetLineColor(kRed);

    f1_coshy_cent[icent] =
        new TF1(Form("f1_coshy_icent_%d", icent), coshy, -1.0, 1.0, 1);
    f1_coshy_cent[icent]->SetParNames("Teff");
    f1_coshy_cent[icent]->SetParameter(0, 0.150);
    f1_coshy_cent[icent]->SetLineColor(kBlue);
    f1_coshy_cent[icent]->SetLineWidth(2);
    f1_coshy_cent[icent]->SetLineStyle(2);

    TGraphErrors *gtmp = new TGraphErrors(8);
    for (int i = 0; i < 8; i++)
      gtmp->SetPoint(i, hdTdy[icent]->GetBinCenter(i + 1),
                     hdTdy[icent]->GetBinContent(i + 1));
    for (int i = 0; i < 8; i++)
      gtmp->SetPointError(i, 0.1, hdTdy[icent]->GetBinError(i + 1));
    gtmp->Fit(Form("f1_coshy_icent_%d", icent), "EX0+", "", -1.0, 1.0);
    //    hdTdy[icent]->Fit(Form("f1_coshy_icent_%d", icent),"0+","",-1.0,1.0);
    f1_coshy_cent[icent]->Draw("SAME");

    hTeff->SetBinContent(icent + 1, f1_coshy_cent[icent]->GetParameter(0));
    hTeff->SetBinError(icent + 1, f1_coshy_cent[icent]->GetParError(0));

    TLatex *tex_Teff_cent =
        new TLatex(-0.75, hdTdy[icent]->GetBinContent(4) * 1.5,
                   Form("T_{eff} = %4.3f #pm %4.3f GeV",
                        f1_coshy_cent[icent]->GetParameter(0),
                        f1_coshy_cent[icent]->GetParError(0)));
    tex_Teff_cent->SetTextFont(42);
    tex_Teff_cent->SetTextSize(0.05);
    tex_Teff_cent->Draw("same");

    TLegend *legTcent = new TLegend(0.3, 0.2, 0.75, 0.5);
    legTcent->SetFillColor(10);
    legTcent->SetLineStyle(3004);
    legTcent->SetLineColor(10);
    legTcent->SetLineWidth(0.);
    legTcent->SetTextFont(42);
    legTcent->SetTextSize(0.04);
    legTcent->SetHeader(
        Form("K^{*0}, %s Au+Au #sqrt{s_{NN}} = 3 GeV", centName[icent].Data()));
    legTcent->AddEntry(hdTdy[icent], "-0.8 < y < 0.2 (data)", "lp");
    legTcent->AddEntry(hdTdy_rf[icent], "0 < y < 0.8 (reflection)", "lp");
    legTcent->AddEntry(f1_coshy_cent[icent], "CosH(y) Fit: T_{eff}/cosh(y)",
                       "lp");
    legTcent->AddEntry(hdTdy_Mid[icent],
                       Form("%s (data)", etaNameMid[0].Data()), "lp");
    legTcent->Draw("same");
  }
  ca_T->cd();
  ca_T->Print(Form("../../files/plots_%s/dTdY_%s_%s_EPmethod_fxt3GeV.png",
                   cutName.Data(), particleName.Data(), cutName.Data()));

  //================================================================================================
  // dN/dy plots
  TH1F *hXidNdy[nCentFlow];
  TH1F *hXidNdy_rf[nCentFlow];
  TH1F *hXidNdy_Mid[nCentFlow];
  TF1 *f1_gaus_cent[nCentFlow];

  for (int icent = 0; icent < nCentFlow; icent++) {
    //    hXidNdy[icent] = new TH1F(Form("hXidNdy_icent%d", icent), "", 5, -0.8,
    //    0.2);
    hXidNdy[icent] = new TH1F(Form("hXidNdy_icent%d", icent), "", 6, -0.9, 0.9);
    hXidNdy_rf[icent] =
        new TH1F(Form("hXidNdy_rf_icent%d", icent), "", 3, 0, 0.9);
    hXidNdy_Mid[icent] =
        new TH1F(Form("hXidNdy_Mid_icent%d", icent), "", 1, fctetaXaxisMid);
  }

  for (int icent = 0; icent < nCentFlow; icent++) {
    //    for(int ieta=0; ieta<5; ieta++) {
    //      hXidNdy[icent]->SetBinContent(ieta+1, xiYield_int[icent][ieta]);
    //      hXidNdy[icent]->SetBinError(ieta+1, xiYieldErr_int[icent][ieta]);
    //    }
    for (int ieta = 0; ieta < 3; ieta++) {
      hXidNdy[icent]->SetBinContent(ieta + 1, xiYield_int[icent][ieta + 1]);
      if (debug) {
        std::cout << "zla check icent == " << icent << " ieta == " << ieta
                  << " yield == " << xiYield_int[icent][ieta] << std::endl;
      }
      hXidNdy[icent]->SetBinError(ieta + 1, xiYieldErr_int[icent][ieta + 1]);
    }
    for (int ieta = 0; ieta < 3; ieta++) {
      hXidNdy_rf[icent]->SetBinContent(ieta + 1,
                                       xiYield_int[icent][2 - ieta + 1]);
      hXidNdy_rf[icent]->SetBinError(ieta + 1,
                                     xiYieldErr_int[icent][2 - ieta + 1]);

      // hXidNdy[icent]->SetBinContent(ieta+5, xiYield_int[icent][3-ieta]);
      // hXidNdy[icent]->SetBinError(ieta+5, xiYieldErr_int[icent][3-ieta]);
    }

    hXidNdy_Mid[icent]->SetBinContent(1, xiYield_int[icent][nTotEta - 1]);
    hXidNdy_Mid[icent]->SetBinError(1, xiYieldErr_int[icent][nTotEta - 1]);
  }

  TCanvas *ca = new TCanvas("ca", "K*0 dNdy", 450 * nCentFlow, 400);
  ca->Divide(nCentFlow, 1);
  for (int icent = 0; icent < nCentFlow; icent++) {
    ca->cd(icent + 1)->SetGridy();

    ax = gPad->DrawFrame(-1.0, 0, 1.0, hXidNdy[icent]->GetBinContent(3) * 1.7);
    SetAxis(ax, 1, 1);
    ax->SetYTitle("K^{*0} dN/dy");
    ax->SetXTitle("y");
    ax->Draw();

    hXidNdy[icent]->Draw("ESAME");
    hXidNdy[icent]->SetMarkerStyle(20);
    hXidNdy[icent]->SetMarkerSize(1.0);
    hXidNdy[icent]->SetMarkerColor(kBlue);
    hXidNdy[icent]->SetLineColor(kBlue);

    hXidNdy_rf[icent]->Draw("ESAME");
    hXidNdy_rf[icent]->SetMarkerStyle(24);
    hXidNdy_rf[icent]->SetMarkerSize(1.0);
    hXidNdy_rf[icent]->SetMarkerColor(kBlue);
    hXidNdy_rf[icent]->SetLineColor(kBlue);

    hXidNdy_Mid[icent]->Draw("ESAME");
    hXidNdy_Mid[icent]->SetMarkerStyle(24);
    hXidNdy_Mid[icent]->SetMarkerSize(1.0);
    hXidNdy_Mid[icent]->SetMarkerColor(kRed);
    hXidNdy_Mid[icent]->SetLineColor(kRed);

    f1_gaus_cent[icent] =
        new TF1(Form("f1_gaus_cent%d", icent), gaus, -1.0, 1.0, 3);
    f1_gaus_cent[icent]->SetParNames("Scale1", "Mean", "Width");
    f1_gaus_cent[icent]->FixParameter(1, 0);
    f1_gaus_cent[icent]->SetParameter(0, 1e-3);
    f1_gaus_cent[icent]->SetParameter(2, 0.5);
    f1_gaus_cent[icent]->SetLineColor(kBlue);
    f1_gaus_cent[icent]->SetLineWidth(2);
    // TGraphErrors* gtmp = new TGraphErrors(8);
    // for(int i=0; i<8; i++)gtmp->SetPoint(i,
    // hXidNdy[icent]->GetBinCenter(i+1), hXidNdy[icent]->GetBinContent(i+1));
    // for(int i=0; i<8; i++)gtmp->SetPointError(i, 0.1,
    // hXidNdy[icent]->GetBinError(i+1)); gtmp->Fit(Form("f1_gaus_cent%d",
    // icent),"EX0","",-1.0,1.0);
    hXidNdy[icent]->Fit(Form("f1_gaus_cent%d", icent), "I+", "", -1.0, 1.);
    f1_gaus_cent[icent]->Draw("SAME");

    TLatex *tex_sigma_cent = new TLatex(
        -0.3, 0.2 * hXidNdy[icent]->GetBinContent(4),
        Form("#sigma = %3.3f#pm%3.3f", f1_gaus_cent[icent]->GetParameter(2),
             f1_gaus_cent[icent]->GetParError(2)));
    tex_sigma_cent->SetTextFont(42);
    tex_sigma_cent->SetTextSize(0.05);
    tex_sigma_cent->Draw("same");

    TLegend *leg = new TLegend(0.2, 0.75, 0.45, 0.9);
    leg->SetFillColor(10);
    leg->SetLineStyle(3004);
    leg->SetLineColor(10);
    leg->SetLineWidth(0.);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->SetHeader(
        Form("K^{*0}, %s Au+Au #sqrt{s_{NN}} = 3 GeV", centName[icent].Data()));
    leg->AddEntry(hXidNdy[icent], "-0.9 < y < 0 (data)", "lp");
    leg->AddEntry(hXidNdy_rf[icent], "0 < y < 0.9 (reflection)", "lp");
    leg->AddEntry(hXidNdy_Mid[icent], Form("%s (data)", etaNameMid[0].Data()),
                  "lp");
    leg->Draw("same");
  }
  ca->cd();
  ca->Print(Form("../../files/plots_%s/dNdY_%s_%s_EPmethod_fxt3GeV.png",
                 cutName.Data(), particleName.Data(), cutName.Data()));

  //================================================================================================
  auto fs = new TFile(Form("../../files/%s.root", cutName.Data()), "RECREATE");
  fs->cd();

  // dN/dpT/dy
  for (int ieta = 0; ieta < nTotEta; ieta++) {
    for (int icent = 0; icent < nCentFlow; icent++) {
      hXidNdpT[ieta][icent]->Write();
      ptexp[icent][ieta]->Write();
    }
  }

  for (int icent = 0; icent < nCentFlow; icent++) {
    for (int ieta = 0; ieta < nTotEta; ieta++) {
      hXiYieldVsPt_raw[icent][ieta]->Write();
    }
  }

  for (int icent = 0; icent < nCentFlow; icent++) {
    for (int ieta = 0; ieta < nTotEta; ieta++) {
      hXiYieldVsPt[icent][ieta]->Write(
          Form("hXiYieldVsPt_cent%d_Y%d", icent, ieta));
    }
  }
  for (int icent = 0; icent < nCentFlow; icent++) {
    for (int ieta = 0; ieta < nTotEta; ieta++) {
      hXiYieldVsMt[icent][ieta]->Write(
          Form("hXiYieldVsMt_cent%d_Y%d", icent, ieta));
    }
  }

  TGraphErrors *gXiYieldVsPt[nCentFlow][nTotEta];
  TGraphErrors *gXiYieldVsMt[nCentFlow][nTotEta];

  for (int icent = 0; icent < nCentFlow; icent++) {
    for (int ieta = 0; ieta < nTotEta; ieta++) {
      gXiYieldVsPt[icent][ieta] = convertToGraph(hXiYieldVsPt[icent][ieta]);
      gXiYieldVsMt[icent][ieta] = convertToGraph(hXiYieldVsMt[icent][ieta]);

      gXiYieldVsPt[icent][ieta]->Write(
          Form("gXiYieldVsPt_cent%d_Y%d", icent, ieta));
      gXiYieldVsMt[icent][ieta]->Write(
          Form("gXiYieldVsMt_cent%d_Y%d", icent, ieta));

      ptdptexp[icent][ieta]->Write();
      ptmptexp[icent][ieta]->Write();
    }
  }

  for (int icent = 0; icent < nCentFlow; icent++) {
    hdTdy[icent]->Write();
    hdTdy_rf[icent]->Write();
    hdTdy_Mid[icent]->Write();
    hXidNdy[icent]->Write();
    hXidNdy_rf[icent]->Write();
    hXidNdy_Mid[icent]->Write();
    f1_coshy_cent[icent]->Write();
  }
  hTeff->Write("hTeff");

  fs->Close();
  return;
}

TGraphErrors *convertToGraph(TH1F *h) {
  int nbins = h->GetNbinsX();

  TGraphErrors *g = new TGraphErrors(nbins);
  for (int i = 0; i < nbins; i++) {
    g->SetPoint(i, h->GetBinCenter(i + 1), h->GetBinContent(i + 1));
    g->SetPointError(i, 0.1, h->GetBinError(i + 1));
  }

  return g;
}
