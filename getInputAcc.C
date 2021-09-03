#include <iostream>

void getInputAcc()
{
	static const int nCent = 4;
	TFile *fInput = TFile::Open("../../InputMCTrackTree_v1.root");
	TFile *fOutput = new TFile(Form("../../InputKstarInfo_cent%i.root", nCent), "RECREATE");

	TH2F *hKstarAcc[nCent];
	for(int icent = 0; icent < nCent; ++icent) {
		hKstarAcc[icent] = new TH2F(Form("hKstarAcc_cent%i", icent), "Input Kstar Acc", 300, -1.5, 1.5, 100, 0, 5);
	}

	int cent;
	float px[5], py[5], pz[5], pt[5], eta[5], y[5];
	TTree *tInput = (TTree*)fInput->Get("mMcKstar");

	tInput->SetBranchAddress("mPx", px);
	tInput->SetBranchAddress("mPy", py);
	tInput->SetBranchAddress("mPz", pz);
	tInput->SetBranchAddress("mY", y);
	tInput->SetBranchAddress("mCent", &cent);

	Long64_t nEvent = tInput->GetEntries();
	for(Long64_t ievent = 0; ievent < nEvent; ++ievent) {
		tInput->GetEntry(ievent);
		for(int ikstar = 0; ikstar < 5; ++ikstar) {
			if(nCent == 9) {
				hKstarAcc[cent]->Fill(y[ikstar], sqrt(px[ikstar]*px[ikstar] + py[ikstar]*py[ikstar]));
			} else if(nCent == 4) {
				if(cent >= 7) {
					hKstarAcc[0]->Fill(y[ikstar], sqrt(px[ikstar]*px[ikstar] + py[ikstar]*py[ikstar]));
				} else if(cent >= 4) {
					hKstarAcc[1]->Fill(y[ikstar], sqrt(px[ikstar]*px[ikstar] + py[ikstar]*py[ikstar]));
				} else if(cent >= 0) {
					hKstarAcc[2]->Fill(y[ikstar], sqrt(px[ikstar]*px[ikstar] + py[ikstar]*py[ikstar]));
				}
				hKstarAcc[3]->Fill(y[ikstar], sqrt(px[ikstar]*px[ikstar] + py[ikstar]*py[ikstar]));
			}
		}
	}
	fOutput->cd();
	for(int icent = 0; icent < nCent; ++icent) {
		hKstarAcc[icent]->Write();
	}
}
