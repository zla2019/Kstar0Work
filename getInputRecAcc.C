#include <iostream>

void getInputRecAcc()
{
	TFile *fInput = TFile::Open("/home/zla/Kstar0Workdir/embedding_data/embedding_Rec_20211101.root");
	TH3F *h3Dplot[9];
	TH2F *hAcc[9];
	for(int icent = 0; icent < 9; ++icent) {
		h3Dplot[icent] = (TH3F*)fInput->Get(Form("NumKstar0InvMassvsPtY_Cent%i", icent));
		std::cout << "entries: " << h3Dplot[icent]->GetEntries() << std::endl;
		hAcc[icent] = (TH2F*)h3Dplot[icent]->Project3D("xy");
	}
	TFile *fOutput = new TFile("/home/zla/Kstar0Workdir/embedding_data/InputRecKstar_cent9_20211101.root", "RECREATE");
	fOutput->cd();
	for(int icent = 0; icent < 9; ++icent) {
		hAcc[icent]->SetName(Form("hKstarAcc_cent%i", icent));
		hAcc[icent]->Write();
	}

	TFile *fOutput2 = new TFile("/home/zla/Kstar0Workdir/embedding_data/InputRecKstar_cent4_20211101.root", "RECREATE");
	fOutput2->cd();
	TH2F *hAccCent4[4];
	hAccCent4[0] = (TH2F*)hAcc[1]->Clone();
	hAccCent4[0]->Add(hAcc[0]);
	hAccCent4[1] = (TH2F*)hAcc[2]->Clone();
	hAccCent4[1]->Add(hAcc[3]);
	hAccCent4[1]->Add(hAcc[4]);
	hAccCent4[2] = (TH2F*)hAcc[5]->Clone();
	hAccCent4[2]->Add(hAcc[6]);
	hAccCent4[3] = (TH2F*)hAcc[0]->Clone();
	for(int icent = 0; icent < 7; ++icent) {
		hAccCent4[3]->Add(hAcc[icent]);
	}
	for(int icent = 0; icent < 4; ++icent) {
		hAccCent4[icent]->SetName(Form("hKstarAcc_cent%i", icent));
		hAccCent4[icent]->Write();
	}
}
