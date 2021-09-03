#include <TLegend.h>
#include <iostream>

void getEfficiencyCent() {
  static const int nPtBin = 16;
  static const int nEtaBin = 1;
  static const int nTotEtaBin = 7;
  static const int nCent = 4;
  const float ptBinBoundary[nPtBin + 1] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                           1.0, 1.1, 1.2, 1.3, 1.4, 1.5,
                                           1.6, 1.7, 1.8, 1.9, 2.0};
  const float etaBinBoundary[nEtaBin + 1] = {-0.8, 0};
  const std::string centBinBoundary[10] = {"80%", "70%", "60%", "50%", "40%",
                                           "30%", "20%", "10%", "5%",  "0%"};
  const std::string centBinBoundary2[5] = {"0%", "10%", "40%", "80%", "80%"};

  TFile *fRec;
  TFile *fMC;
  TFile *fEfficiency;
  if (nCent == 9) {
    fRec = TFile::Open("../../ReconstructionKstarInfo_v4.root");
    fMC = TFile::Open("../../InputKstarInfo_v4.root");
    fEfficiency = new TFile("../../KstarEfficiency_Cent.root", "RECREATE");
  } else if (nCent == 4) {
    fRec = TFile::Open("../../ReconstructionKstarInfo_v3.root");
    fMC = TFile::Open("../../InputKstarInfo_v4.root");
    fEfficiency = new TFile("../../KstarEfficiency_Cent_test.root", "RECREATE");
  }

  TH2F *hRecAcc[nCent];
  TH2F *hMCAcc[nCent];
  TH2F *hEfficiency[nCent];
  TH1F *hEfficiencyRap[nCent];

  for (int icent = 0; icent < nCent; ++icent) {
    if (nCent == 9) {
      hRecAcc[icent] = (TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 8 - icent));
      hMCAcc[icent] = (TH2F *)fMC->Get(Form("hKstarAcc_cent%i", icent));
    } else if (nCent == 4) {
      hRecAcc[icent] =
          (TH2F *)fRec->Get(Form("hKstar0RapidityvsPt_cent%i", icent));
	hMCAcc[icent] = (TH2F*)(((TH2F *)fMC->Get(Form("hKstarAcc_cent%i", icent)))->Clone());
      if (icent == 2) {
        hRecAcc[icent]->Reset();
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 6)));
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 5)));
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 7)));
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 8)));

        hMCAcc[icent]->Reset();
        hMCAcc[icent]->Add((TH2F *)fMC->Get(Form("hKstarAcc_cent%i", 0)));
        hMCAcc[icent]->Add((TH2F *)fMC->Get(Form("hKstarAcc_cent%i", 1)));
        hMCAcc[icent]->Add((TH2F *)fMC->Get(Form("hKstarAcc_cent%i", 2)));
        hMCAcc[icent]->Add((TH2F *)fMC->Get(Form("hKstarAcc_cent%i", 3)));

      } else if (icent == 0) {
        hRecAcc[icent]->Reset();
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 0)));
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 1)));
        hMCAcc[icent]->Reset();
        hMCAcc[icent]->Add((TH2F *)fMC->Get(Form("hKstarAcc_cent%i", 7)));
        hMCAcc[icent]->Add((TH2F *)fMC->Get(Form("hKstarAcc_cent%i", 8)));
      } else if (icent == 1) {
        hRecAcc[icent]->Reset();
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 2)));
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 3)));
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 4)));
        hMCAcc[icent]->Reset();
        hMCAcc[icent]->Add((TH2F *)fMC->Get(Form("hKstarAcc_cent%i", 4)));
        hMCAcc[icent]->Add((TH2F *)fMC->Get(Form("hKstarAcc_cent%i", 5)));
        hMCAcc[icent]->Add((TH2F *)fMC->Get(Form("hKstarAcc_cent%i", 6)));
      } else if (icent == 3) {
        hRecAcc[icent]->Reset();
        hMCAcc[icent]->Reset();
        for (int icent2 = 0; icent2 < 9; ++icent2) {
          hRecAcc[icent]->Add(
              (TH2F *)fRec->Get(Form("hKstarAcc_cent%i", icent2)));
          hMCAcc[icent]->Add(
              (TH2F *)fMC->Get(Form("hKstarAcc_cent%i", icent2)));
        }
      }
//      hMCAcc[icent] = (TH2F *)fMC->Get(Form("hKstarAcc_cent%i", icent));
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
      // TH1F(Form("Efficiency_cent%i_rap%i", icent, ieta), Form("Efficiency
      // cent %i, y %i", icent, ieta), nEtaBin, etaBinBoundary[0],
      // etaBinBoundary[nEtaBin]);
      hEfficiencyRap[icent] = (TH1F *)hEfficiency[icent]->ProjectionY(
          Form("Efficiency_y_cent%i", icent), ieta + 1, ieta + 1);
    }
  }
  fEfficiency->cd();
  TCanvas *cTmp = new TCanvas("", "", 600, 400);
  cTmp->cd();
  TLegend *leg = new TLegend(0.11, 0.5, 0.5, 0.87);
  for (int icent = 0; icent < nCent; ++icent) {
    leg->SetLineWidth(0);
    for (int ieta = 0; ieta < nEtaBin; ++ieta) {
      hEfficiencyRap[icent]->Write();
      hEfficiency[icent]->Write();
      if (icent == 0) {
        hEfficiencyRap[icent]->SetLineColor(1);
        hEfficiencyRap[icent]->GetYaxis()->SetRangeUser(0, 0.75);
        hEfficiencyRap[icent]->SetStats(0);
        if (nCent == 9) {
          leg->AddEntry(hEfficiencyRap[icent],
                        Form("%s < cent < %s", centBinBoundary[icent].c_str(),
                             centBinBoundary[icent + 1].c_str()),
                        "lep");
        } else if (nCent == 4) {
          leg->AddEntry(hEfficiencyRap[icent],
                        Form("%s < cent < %s", centBinBoundary2[icent - (icent / 3) * 3].c_str(),
                             centBinBoundary2[icent + 1].c_str()),
                        "lep");
        }
	hEfficiencyRap[icent]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hEfficiencyRap[icent]->GetYaxis()->SetTitle("Efficiency");
        hEfficiencyRap[icent]->Draw();
        leg->Draw("same");
      } else {
        if (nCent == 9) {
          leg->AddEntry(hEfficiencyRap[icent],
                        Form("%s < cent < %s", centBinBoundary[icent].c_str(),
                             centBinBoundary[icent + 1].c_str()),
                        "lep");
        } else if (nCent == 4) {
          leg->AddEntry(hEfficiencyRap[icent],
                        Form("%s < cent < %s", centBinBoundary2[icent - (icent / 3) * 3].c_str(),
                             centBinBoundary2[icent + 1].c_str()),
                        "lep");
        }
        hEfficiencyRap[icent]->SetLineColor(1 + icent);
	hEfficiencyRap[icent]->SetMarkerColor(1 + icent);
        hEfficiencyRap[icent]->Draw("same");
      }
    }
    //    hEfficiency[icent]->Write();
  }
  //  hEfficiency[3]->Draw("colztext");
}
