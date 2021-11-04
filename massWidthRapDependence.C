#include <iostream>
#include "plot_Kstar0InvMass_CheckRapDepend.C"

void massWidthRapDependence() {
	gSystem->Load("plot_Kstar0InvMass_CheckRapDepend.C");
	const int rapBins = 5;
	float massWidth[rapBins] = { 0 };
	float massWidthErr[rapBins] = { 0 };
	const float raplimit[rapBins + 1] = { -1.0, -0.8, -0.6, -0.4, -0.2, 0};
	float rapMid[rapBins] = { 0 };
	float rapWidth[rapBins] = { 0 };
	for(int i = 0; i < rapBins; ++i) {
		rapMid[i] = (raplimit[i] + raplimit[i + 1]) / 2;
		rapWidth[i] = (raplimit[i + 1] - raplimit[i]) / 2;
	}
	for(int i = 0; i < rapBins; ++i) {
		plot_Kstar0InvMass_CheckRapDepend(1, &massWidth[i], &massWidthErr[i], &raplimit[i]);
	}
	for(int i = 0; i < rapBins; ++i) {
		std::cout << massWidth[i] << " ";
	}
	std::cout << std::endl;
	for(int i = 0; i < rapBins; ++i) {
		std::cout << massWidthErr[i] << " ";
	}
	std::cout << std::endl;

	TGraphErrors *gr = new TGraphErrors(rapBins, rapMid, massWidth, rapWidth, massWidthErr);
	gr->SetTitle("Mass Width Rap Dependence, 0-60% centrality");
	gr->GetXaxis()->SetTitle("y");
	gr->GetYaxis()->SetTitle("Mass Width");
	float massWidth2[1] = { 0 };
	float massWidthErr2[1] = { 0 };
	const float raplimit2[2] = { -1., 0 };
	plot_Kstar0InvMass_CheckRapDepend(1, &massWidth2[0], &massWidthErr2[0], &raplimit2[0]);
	gStyle->SetPadTopMargin(0.115);
	gStyle->SetPadBottomMargin(0.085);
	gStyle->SetPadRightMargin(0.07);
	gStyle->SetPadLeftMargin(0.1);
	gStyle->SetOptTitle(1);
	gStyle->SetOptStat(1);
        TCanvas *cTmp = new TCanvas("", "", 1280, 720);
        cTmp->cd();
	gr->Draw();
	TLine *massLine = new TLine(-0.8, massWidth2[0], 0, massWidth2[0]);
	massLine->SetLineStyle(2);
	massLine->SetLineColor(kRed);
	massLine->SetLineWidth(3);
	massLine->Draw("same");
	TText *test1 = new TText(-0.4, 0.4, "Mass width at y:-0.8~0");
//	test1->AddText("Mass width at y:-0.8~0");
	test1->Draw("same");
}
