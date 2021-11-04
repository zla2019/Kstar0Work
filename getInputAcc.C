#include <iostream>

void getInputAcc()
{
	TFile *fInput = TFile::Open("../../embedding_data/embedding_Mc_20211101.root");
	TFile *fOutput = new TFile("../../embedding_data/InputKstarInfo_cent4_20211101.root", "RECREATE");
	TFile *fOutput2 = new TFile("../../embedding_data/InputKstarInfo_cent9_20211101.root", "RECREATE");

	TH2F *hKstarAcc[4];
	for(int icent = 0; icent < 4; ++icent) {
		hKstarAcc[icent] = new TH2F(Form("hKstarAcc_cent%i", icent), "Input Kstar Acc", 20, -1, 1, 100, 0, 5);
	}

	TH2F *hKstarAccCent9[9];
	for(int icent = 0; icent < 9; ++icent) {
		hKstarAccCent9[icent] = new TH2F(Form("hKstarAcc_cent%i", icent), "Input Kstar Acc", 20, -1, 1, 100, 0, 5);
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
			if(cent >= 7) {
				hKstarAcc[0]->Fill(y[ikstar], sqrt(px[ikstar]*px[ikstar] + py[ikstar]*py[ikstar]));
			} else if(cent >= 4) {
				hKstarAcc[1]->Fill(y[ikstar], sqrt(px[ikstar]*px[ikstar] + py[ikstar]*py[ikstar]));
			} else if(cent >= 2) {
				hKstarAcc[2]->Fill(y[ikstar], sqrt(px[ikstar]*px[ikstar] + py[ikstar]*py[ikstar]));
			}
			if(cent >= 2) {
				hKstarAcc[3]->Fill(y[ikstar], sqrt(px[ikstar]*px[ikstar] + py[ikstar]*py[ikstar]));
			}
			hKstarAccCent9[8 - cent]->Fill(y[ikstar], sqrt(px[ikstar]*px[ikstar] + py[ikstar]*py[ikstar]));
		}
	}

	fOutput->cd();
	for(int icent = 0; icent < 4; ++icent) {
		hKstarAcc[icent]->Write();
	}
	fOutput2->cd();
	for(int icent = 0; icent < 9; ++icent) {
		hKstarAccCent9[icent]->Write();
	}
}
