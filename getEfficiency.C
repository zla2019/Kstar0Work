#include <iostream>

void getEfficiency() {
  static const int nPtBin = 16;
  static const int nEtaBin = 8;
  static const int nTotEtaBin = 7;
	static const int nCent = 9;
  const float ptBinBoundary[nPtBin + 1] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
  const float etaBinBoundary[nEtaBin + 1] = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0};

  TFile *fRec = TFile::Open("../../ReconstructionKstarInfo_v3.root");
  TFile *fMC = TFile::Open("../../InputKstarInfo_v3.root");
  TFile *fEfficiency = new TFile("../../KstarEfficiency_test.root", "RECREATE");

  TH2F *hRecAcc[nCent];
  TH2F *hMCAcc[nCent];
  TH2F *hEfficiency[nCent];

  for (int icent = 0; icent < 4; ++icent) {
    hRecAcc[icent] =
        (TH2F *)fRec->Get(Form("hKstar0RapidityvsPt_cent%i", icent));
    hMCAcc[icent] = (TH2F *)fMC->Get(Form("hKstarAcc_cent%i", icent));
/*
    hEfficiency[icent] =
        new TH2F(Form("Efficiency_cent%i", icent), "Efficiency", nEtaBin,
                 etaBinBoundary[0], etaBinBoundary[nEtaBin], nPtBin,
                 ptBinBoundary[0], ptBinBoundary[nPtBin]);
*/
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
    }
  }
  fEfficiency->cd();
  for (int icent = 0; icent < 4; ++icent) {
    hEfficiency[icent]->Write();
  }
  hEfficiency[3]->Draw("colztext");
}
