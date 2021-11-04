#include <iostream>

void getEfficiency() {
  static const int nPtBin = 5;
  static const int nEtaBin = 5;
  static const int nTotEtaBin = 7;
  static const int nCent = 4;
  const float ptBinBoundary[nPtBin + 1] = {0.2, 0.4, 0.8, 
                                           1.2, 1.6, 
                                           2.0};
  const float etaBinBoundary[nEtaBin + 1] = {-0.8, -0.6, -0.4,
                                             -0.2, 0, 0.2};

  TFile *fRec;
  TFile *fMC;
  if (nCent == 4) {
//    fRec = TFile::Open("../../ReconstructionKstarInfo_v3.root");
//    fMC = TFile::Open("../../InputKstarInfo_v4.root");
//    fRec = TFile::Open("../../ReconstructionKstarInfo_v6.root");
    fRec = TFile::Open("../../embedding_data/InputRecKstar_cent4_20211101.root");
//    fMC = TFile::Open("../../InputKstarInfo_cent9.root");
    fMC = TFile::Open("../../embedding_data/InputKstarInfo_cent4_20211101.root");
  } else if (nCent == 9) {
    fRec = TFile::Open("../../embedding_data/embedding_Rec_20211101.root");
    fMC = TFile::Open("../../embedding_data/InputKstarInfo_cent4_20211101.root");
  }
  TFile *fEfficiency = new TFile("../../efficiency/KstarEfficiency_20211101.root", "RECREATE");

  TH2F *hRecAcc[nCent];
  TH2F *hMCAcc[nCent];
  TH2F *hEfficiency[nCent];
  TH2F *hRecAccRebin[nCent];
  TH2F *hMCAccRebin[nCent];

  for (int icent = 0; icent < nCent; ++icent) {
    if (nCent == 4) {
	hRecAcc[icent] = (TH2F*)fRec->Get(Form("hKstarAcc_cent%i", icent));
	hMCAcc[icent] = (TH2F*)fMC->Get(Form("hKstarAcc_cent%i", icent));
    } else if (nCent == 9) {
      hRecAcc[icent] = (TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 8 - icent));
      hMCAcc[icent] = (TH2F *)fMC->Get(Form("hKstarAcc_cent%i", icent));
    }
    /*
        hEfficiency[icent] =
            new TH2F(Form("Efficiency_cent%i", icent), "Efficiency", nEtaBin,
                     etaBinBoundary[0], etaBinBoundary[nEtaBin], nPtBin,
                     ptBinBoundary[0], ptBinBoundary[nPtBin]);
    */
    hEfficiency[icent] =
        new TH2F(Form("Efficiency_cent%i", icent), "Efficiency", nEtaBin,
                 etaBinBoundary, nPtBin, ptBinBoundary);
    hRecAccRebin[icent] =
        new TH2F(Form("hRecAcc_cent%i", icent), "hRecAcc", nEtaBin,
                 etaBinBoundary, nPtBin, ptBinBoundary);
    hMCAccRebin[icent] =
        new TH2F(Form("hMCAcc_cent%i", icent), "hMCAcc", nEtaBin,
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
        int integralRecXLimit[2] = {hRecAcc[icent]->ProjectionX()->FindBin(
                                     etaBinBoundary[ieta] + 0.000001),
                                 hRecAcc[icent]->ProjectionX()->FindBin(
                                     etaBinBoundary[ieta + 1] - 0.000001)};
        int integralRecYLimit[2] = {hRecAcc[icent]->ProjectionY()->FindBin(
                                     ptBinBoundary[ipt] + 0.000001),
                                 hRecAcc[icent]->ProjectionY()->FindBin(
                                     ptBinBoundary[ipt + 1] - 0.000001)};
        float inputCount =
            hMCAcc[icent]->Integral(integralXLimit[0], integralXLimit[1],
                                    integralYLimit[0], integralYLimit[1]);
        float reconstructCount =
            hRecAcc[icent]->Integral(integralRecXLimit[0], integralRecXLimit[1],
                                     integralRecYLimit[0], integralRecYLimit[1]);

        hRecAccRebin[icent]->SetBinContent(ieta + 1, ipt + 1, reconstructCount);
        hMCAccRebin[icent]->SetBinContent(ieta + 1, ipt + 1, inputCount);
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
  for (int icent = 0; icent < nCent; ++icent) {
	hRecAccRebin[icent]->SetEntries(hRecAccRebin[icent]->Integral());
	hMCAccRebin[icent]->SetEntries(hMCAccRebin[icent]->Integral());
    hRecAccRebin[icent]->Write();
    hMCAccRebin[icent]->Write();
    hEfficiency[icent]->Write();
  }
  hEfficiency[2]->Draw("colztext");
}


float getEfficiency(int icent, float etaLower, float etaHigher, float ptLower, float ptHigher) {

  TFile *fRec;
  TFile *fMC;
  fRec = TFile::Open("../../embedding_data/InputRecKstar_cent4.root");
  fMC = TFile::Open("../../embedding_data/InputKstarInfo_cent4.root");

  TH2F *hRecAcc;
  TH2F *hMCAcc;

    hRecAcc = (TH2F*)fRec->Get(Form("hKstarAcc_cent%i", icent));
    hMCAcc = (TH2F*)fMC->Get(Form("hKstarAcc_cent%i", icent));
    int integralXLimit[2] = {hMCAcc->ProjectionX()->FindBin(
                                 etaLower + 0.000001),
                             hMCAcc->ProjectionX()->FindBin(
                                 etaHigher - 0.000001)};
    int integralYLimit[2] = {hMCAcc->ProjectionY()->FindBin(
                                 ptLower + 0.000001),
                             hMCAcc->ProjectionY()->FindBin(
                                 ptHigher - 0.000001)};
    int integralRecXLimit[2] = {hRecAcc->ProjectionX()->FindBin(
                                 etaLower + 0.000001),
                             hRecAcc->ProjectionX()->FindBin(
                                 etaHigher - 0.000001)};
    int integralRecYLimit[2] = {hRecAcc->ProjectionY()->FindBin(
                                 ptLower + 0.000001),
                             hRecAcc->ProjectionY()->FindBin(
                                 ptHigher - 0.000001)};
    float inputCount =
        hMCAcc->Integral(integralXLimit[0], integralXLimit[1],
                                integralYLimit[0], integralYLimit[1]);
    float reconstructCount =
        hRecAcc->Integral(integralRecXLimit[0], integralRecXLimit[1],
                                integralRecYLimit[0], integralRecYLimit[1]);
    float binContent = (float)reconstructCount / (float)inputCount;
	return binContent;
}
