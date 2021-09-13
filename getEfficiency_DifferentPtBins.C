#include <iostream>
#include <TH2Poly.h>

void getEfficiency_DifferentPtBins() {
  static const int nPtBin = 6;
  static const int nEtaBin = 6;
  static const int nTotEtaBin = 7;
  static const int nCent = 4;
  const float ptBinBoundary[nPtBin + 1] = {0.2, 0.4, 0.8, 
                                           1.2, 1.6, 
                                           2.0, 3.0};
  const float etaBinBoundary[nEtaBin + 1] = {-0.8, -0.6, -0.4,
                                             -0.2, 0, 0.2, 0.4};

        float pTBoundaryBuffer[nCent][nEtaBin][nPtBin + 1] = {
                { { 0.2, 0.8, 1.6, 2.0, 0.0, 0.0, 0.0 }, { 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 0.0 }, { 0.2, 0.8, 1.2, 1.6, 3.0, 0.0, 0.0 }, { 0.2, 1.0, 1.6, 3.0, 0.0, 0.0, 0.0 }, { 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0 }, { 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0 } }, //cent 0
                { { 0.2, 0.8, 1.6, 2.0, 0.0, 0.0, 0.0 }, { 0.2, 0.4, 0.8, 1.2, 2.0, 3.0, 0.0 }, { 0.2, 0.4, 1.2, 1.6, 3.0, 0.0, 0.0 }, { 0.2, 0.8, 1.6, 3.0, 0.0, 0.0, 0.0 }, { 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0 }, { 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0 } }, //cent 1
                { { 0.2, 0.4, 1.2, 2.0, 0.0, 0.0, 0.0 }, { 0.2, 0.4, 0.8, 1.2, 1.6, 3.0, 0.0 }, { 0.2, 0.4, 0.8, 1.2, 1.6, 3.0, 0.0 }, { 0.2, 0.4, 0.8, 2.0, 3.0, 0.0, 0.0 }, { 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0 }, { 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0 } }, //cent 2
                { { 0.2, 0.4, 0.8, 1.6, 3.0, 0.0, 0.0 }, { 0.2, 0.4, 0.8, 1.2, 1.6, 3.0, 0.0 }, { 0.2, 0.4, 0.8, 1.2, 1.6, 3.0, 0.0 }, { 0.2, 0.4, 0.8, 1.6, 3.0, 0.0, 0.0 }, { 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0 }, { 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0 } }  //cent 3
        };
        const int nPTBinsArr[nCent][nEtaBin] = { { 3, 5, 4, 3, 6, 6 }, { 3, 5, 4, 3, 6, 6 }, { 3, 5, 5, 4, 6, 6 }, { 4, 5, 5, 4, 6, 6 } };
        const float* pTBoundaryArr[nCent][nEtaBin];
        for(int icent = 0; icent < nCent; ++icent) {
                for(int ieta = 0; ieta < nEtaBin; ++ieta) {
                        pTBoundaryArr[icent][ieta] = &pTBoundaryBuffer[icent][ieta][0];
                }
        }



//	const float ptBinBoundary[nPtBin + 1] = { 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8 };
//	const float etaBinBoundary[nEtaBin + 1] = { -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6 };
	const std::string centName[4] = { "0-10%", "10-40%", "40-60%", "0-60%" };

  TFile *fRec;
  TFile *fMC;
  if (nCent == 4) {
    fRec = TFile::Open("../../ReconstructionKstarInfo_v6.root");
    fMC = TFile::Open("../../InputKstarInfo_cent9.root");
  } else if (nCent == 9) {
    fRec = TFile::Open("../../ReconstructionKstarInfo_v6.root");
    fMC = TFile::Open("../../InputKstarInfo_cent9.root");
  }
  TFile *fEfficiency = new TFile("../../KstarEfficiency_binning2_0.root", "RECREATE");

  TH2F *hRecAcc[nCent];
  TH2F *hMCAcc[nCent];
  TH2F *hEfficiency[nCent];
  TH2F *hRecAccRebin[nCent];
  TH2F *hMCAccRebin[nCent];

	TH2Poly *hEfficiency_poly[nCent];

  for (int icent = 0; icent < nCent; ++icent) {
    if (nCent == 4) {
      hRecAcc[icent] =
          (TH2F *)fRec->Get(Form("hKstar0RapidityvsPt_cent%i", icent));
        hMCAcc[icent] = (TH2F*)(((TH2F *)fMC->Get(Form("hKstarAcc_cent%i", icent)))->Clone());
      if (icent == 2) {
        hRecAcc[icent]->Reset();
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 6)));
        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 5)));
//        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 7)));
//        hRecAcc[icent]->Add((TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 8)));

        hMCAcc[icent]->Reset();
//        hMCAcc[icent]->Add((TH2F *)fMC->Get(Form("hKstarAcc_cent%i", 0)));
//        hMCAcc[icent]->Add((TH2F *)fMC->Get(Form("hKstarAcc_cent%i", 1)));
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
        for (int icent2 = 0; icent2 < 7; ++icent2) {
          hRecAcc[icent]->Add(
              (TH2F *)fRec->Get(Form("hKstarAcc_cent%i", 9 - icent2)));
          hMCAcc[icent]->Add(
              (TH2F *)fMC->Get(Form("hKstarAcc_cent%i", icent2)));
        }
      }
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
        new TH2F(Form("Efficiency_old_cent%i", icent), "Efficiency", nEtaBin,
                 etaBinBoundary, nPtBin, ptBinBoundary);


	hEfficiency_poly[icent] = new TH2Poly();
	hEfficiency_poly[icent]->SetName(Form("Efficiency_cent%i", icent));
	for(int ieta = 0; ieta < 6; ++ieta) {
		for(int ipt = 0; ipt < nPTBinsArr[icent][ieta]; ++ipt) {
			double xArr2[2] = { etaBinBoundary[ieta], etaBinBoundary[ieta + 1] };
			double yArr2[2] = { pTBoundaryArr[icent][ieta][ipt + 1], pTBoundaryArr[icent][ieta][ipt] };
			hEfficiency_poly[icent]->AddBin(etaBinBoundary[ieta], pTBoundaryArr[icent][ieta][ipt], etaBinBoundary[ieta + 1], pTBoundaryArr[icent][ieta][ipt + 1]);
		}
	}


    hRecAccRebin[icent] =
        new TH2F(Form("hRecAcc_cent%i", icent), "hRecAcc", nEtaBin,
                 etaBinBoundary, nPtBin, ptBinBoundary);
    hMCAccRebin[icent] =
        new TH2F(Form("hMCAcc_cent%i", icent), "hMCAcc", nEtaBin,
                 etaBinBoundary, nPtBin, ptBinBoundary);

    for (int ieta = 0; ieta < nEtaBin; ++ieta) {
      for (int ipt = 0; ipt < nPTBinsArr[icent][ieta]; ++ipt) {
        int integralXLimit[2] = {hMCAcc[icent]->ProjectionX()->FindBin(
                                     etaBinBoundary[ieta] + 0.000001),
                                 hMCAcc[icent]->ProjectionX()->FindBin(
                                     etaBinBoundary[ieta + 1] - 0.000001)};
        int integralYLimit[2] = {hMCAcc[icent]->ProjectionY()->FindBin(
                                     pTBoundaryArr[icent][ieta][ipt] + 0.000001),
                                 hMCAcc[icent]->ProjectionY()->FindBin(
                                     pTBoundaryArr[icent][ieta][ipt + 1] - 0.000001)};
        float inputCount =
            hMCAcc[icent]->Integral(integralXLimit[0], integralXLimit[1],
                                    integralYLimit[0], integralYLimit[1]);
        float reconstructCount =
            hRecAcc[icent]->Integral(integralXLimit[0], integralXLimit[1],
                                     integralYLimit[0], integralYLimit[1]);

        hRecAccRebin[icent]->SetBinContent(ieta + 1, ipt + 1, reconstructCount);
        hMCAccRebin[icent]->SetBinContent(ieta + 1, ipt + 1, inputCount);
        float binContent = (float)reconstructCount / (float)inputCount;
        if (binContent >= 1 || binContent == 0 || inputCount == 0) {
          continue;
        } else {
//		double xArrTmp[2] = { etaBinBoundary[ieta] + 0.000001, etaBinBoundary[ieta + 1] - 0.000001 };
//		double yArrTmp[2] = { pTBoundaryArr[icent][ieta][ipt + 1] - 0.000001, pTBoundaryArr[icent][ieta][ipt] + 0.000001 };
//		hEfficiency_poly[icent]->AddBin(2, xArrTmp, yArrTmp);
	
/*                double xArr2[2] = { etaBinBoundary[ieta], etaBinBoundary[ieta + 1] };
                double yArr2[2] = { pTBoundaryArr[icent][ieta][ipt + 1], pTBoundaryArr[icent][ieta][ipt] };
                int binNumber = hEfficiency_poly[icent]->AddBin(xArr2[0], yArr2[1], xArr2[1], yArr2[0]);*/
		int binNumber = hEfficiency_poly[icent]->FindBin((etaBinBoundary[ieta] + etaBinBoundary[ieta + 1]) / 2, (pTBoundaryArr[icent][ieta][ipt] + pTBoundaryArr[icent][ieta][ipt + 1]) / 2);
		hEfficiency_poly[icent]->SetBinContent(binNumber, binContent);
//          hEfficiency[icent]->SetBinContent(ieta + 1, ipt + 1, binContent);
        }
      }
    }
  }
  fEfficiency->cd();
        TCanvas *cTmp = new TCanvas("", "", 1280 * 2, 1080);
        cTmp->Divide(2, 2);
  for (int icent = 0; icent < nCent; ++icent) {
	cTmp->cd(icent + 1);
    hRecAccRebin[icent]->Write();
    hMCAccRebin[icent]->Write();
	hEfficiency_poly[icent]->GetXaxis()->SetTitle("y");
	hEfficiency_poly[icent]->GetYaxis()->SetTitle("p_{T} [GeV/c]");
	hEfficiency_poly[icent]->SetTitle(Form("Efficiency, %s", centName[icent].c_str()));
	hEfficiency_poly[icent]->SetStats(0);
	hEfficiency_poly[icent]->SetMarkerSize(2);
	hEfficiency_poly[icent]->Draw("colztext");
    hEfficiency_poly[icent]->Write();
  }
/*	TCanvas *cEfficiencySingle = new TCanvas("", "", 2560, 1440);
	hEfficiency[3]->GetXaxis()->SetTitle("y");
	hEfficiency[3]->GetYaxis()->SetTitle("p_{T} [GeV/c]");
	hEfficiency[3]->SetTitle(Form("Efficiency, %s", centName[3].c_str()));
	hEfficiency[3]->SetStats(0);*/
//	hEfficiency_poly[3]->Draw("colztext");
//	hEfficiency[3]->SetMarkerSize(2);
//	hEfficiency[3]->Draw("colztext");
//  hEfficiency[2]->Draw("colztext");
}
