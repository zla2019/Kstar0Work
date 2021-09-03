#include <TLegend.h>
#include <iostream>

void getEfficiencyRap() {
  static const int nPtBin = 16;
  static const int nEtaBin = 8;
  static const int nTotEtaBin = 7;
  const float ptBinBoundary[nPtBin + 1] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                           1.0, 1.1, 1.2, 1.3, 1.4, 1.5,
                                           1.6, 1.7, 1.8, 1.9, 2.0};
  const float etaBinBoundary[nEtaBin + 1] = {-0.8, -0.7, -0.6, -0.5, -0.4,
                                             -0.3, -0.2, -0.1, 0};

  TFile *fRec = TFile::Open("../../ReconstructionKstarInfo_v3.root");
  TFile *fMC = TFile::Open("../../InputKstarInfo_cent4.root");
  TFile *fEfficiency = new TFile("../../KstarEfficiency_Rap.root", "RECREATE");
	const std::string centRange[4] = {"0-10%", "10-40%", "40-80%", "0-80%"};

  TH2F *hRecAcc[4];
  TH2F *hMCAcc[4];
  TH2F *hEfficiency[4];
  TH1F *hEfficiencyRap[4][nEtaBin];

  for (int icent = 0; icent < 4; ++icent) {
      hRecAcc[icent] =
          (TH2F *)fRec->Get(Form("hKstar0RapidityvsPt_cent%i", icent));
      if (icent == 2) {
        hRecAcc[icent]->Reset();
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 6)));
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 5)));
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 7)));
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 8)));
      } else if (icent == 0) {
        hRecAcc[icent]->Reset();
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 0)));
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 1)));
      } else if (icent == 1) {
        hRecAcc[icent]->Reset();
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 2)));
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 3)));
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 4)));
      } else if (icent == 3) {
        hRecAcc[icent]->Reset();
        for (int icent2 = 0; icent2 < 9; ++icent2) {
          hRecAcc[icent]->Add(
              (TH2F *)fRec->Get(Form("hKstarAcc_cent%i", icent2)));
        }
      }
      hMCAcc[icent] = (TH2F *)fMC->Get(Form("hKstarAcc_cent%i", icent));
    hMCAcc[icent] = (TH2F *)fMC->Get(Form("hKstarAcc_cent%i", icent));
    hEfficiency[icent] =
        new TH2F(Form("Efficiency_cent%i", icent), "Efficiency", nEtaBin,
                 etaBinBoundary, nPtBin, ptBinBoundary);

    for (int ieta = 0; ieta < nEtaBin; ++ieta) {
      for (int ipt = 0; ipt < nPtBin; ++ipt) {
        int integralXLimit[2] = {hMCAcc[icent]->ProjectionX()->FindBin(
                                     etaBinBoundary[ieta] + 0.000001),
                                 hMCAcc[icent]->ProjectionX()->FindBin(
                                     etaBinBoundary[ieta + 1] - 0.000001)};
        int integralYLimit[2] = {hMCAcc[icent]->ProjectionY()->FindBin(
                                     ptBinBoundary[ipt] + 0.000001),
                                 hMCAcc[icent]->ProjectionY()->FindBin(
                                     ptBinBoundary[ipt + 1] - 0.000001)};
        float inputCount =
            hMCAcc[icent]->Integral(integralXLimit[0], integralXLimit[1],
                                    integralYLimit[0], integralYLimit[1]);
        float reconstructCount =
            hRecAcc[icent]->Integral(integralXLimit[0], integralXLimit[1],
                                     integralYLimit[0], integralYLimit[1]);
        float binContent = (float)reconstructCount / (float)inputCount;
        if (binContent > 1 || binContent == 0 || inputCount == 0) {
          continue;
        } else {
          hEfficiency[icent]->SetBinContent(ieta + 1, ipt + 1, binContent);
        }
      }
      //	hEfficiencyRap[icent][ieta] = new
      //TH1F(Form("Efficiency_cent%i_rap%i", icent, ieta), Form("Efficiency cent
      //%i, y %i", icent, ieta), nEtaBin, etaBinBoundary[0],
      //etaBinBoundary[nEtaBin]);
      hEfficiencyRap[icent][ieta] = (TH1F *)hEfficiency[icent]->ProjectionY(
          Form("Efficiency_cent%i_rap%i", icent, ieta), ieta + 1, ieta + 1);
    }
  }
  fEfficiency->cd();
  TCanvas *cTmp = new TCanvas("", "", 700 * 4, 600);
  cTmp->Divide(4, 1);
  for (int icent = 0; icent < 4; ++icent) {
    cTmp->cd(icent + 1);
    TLegend *leg = new TLegend(0.11, 0.5, 0.5, 0.87);
    leg->SetLineWidth(0);
    for (int ieta = 0; ieta < nEtaBin; ++ieta) {
      hEfficiencyRap[icent][ieta]->Write();
      if (ieta == 0) {
        hEfficiencyRap[icent][ieta]->SetLineColor(1);
        hEfficiencyRap[icent][ieta]->GetYaxis()->SetRangeUser(0, 0.75);
        hEfficiencyRap[icent][ieta]->SetStats(0);
        leg->AddEntry(hEfficiencyRap[icent][ieta],
                      Form("%0.2f < y < %0.2f", etaBinBoundary[ieta],
                           etaBinBoundary[ieta + 1]),
                      "lep");
	hEfficiencyRap[icent][ieta]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hEfficiencyRap[icent][ieta]->GetYaxis()->SetTitle("Efficiency");
	hEfficiencyRap[icent][ieta]->SetTitle(Form("Efficiency, %s", centRange[icent].c_str()));
        hEfficiencyRap[icent][ieta]->Draw();
        leg->Draw("same");
      } else {
        leg->AddEntry(hEfficiencyRap[icent][ieta],
                      Form("%0.2f < y < %0.2f", etaBinBoundary[ieta],
                           etaBinBoundary[ieta + 1]),
                      "lep");
        hEfficiencyRap[icent][ieta]->SetLineColor(1 + ieta);
        hEfficiencyRap[icent][ieta]->Draw("same");
      }
    }
    //    hEfficiency[icent]->Write();
  }
  //  hEfficiency[3]->Draw("colztext");
}
