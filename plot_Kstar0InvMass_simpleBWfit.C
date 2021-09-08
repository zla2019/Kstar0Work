#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TProfile.h"

double breitwigner(Double_t *x,Double_t *par) // Breit-Wigner function
{
    Double_t fitsig = (par[0]*par[1]/(2.0*3.14159265358979323846)) / (par[1]*par[1]/4.0 + (x[0] - par[2])*(x[0] - par[2])) ;
    return fitsig;
}

double background(Double_t *x,Double_t *par) //pol2
{
    Double_t fitbg = par[0]+par[1]*x[0]+par[2]*x[0]*x[0];//+par[3]*x[0]*x[0]*x[0];
    return fitbg;
}

double fitFuncbw(Double_t *x,Double_t *par) //bw + pol2
{
    Double_t fitval = background(x,&par[3])  + breitwigner(x,par);
    return fitval;
}

TH1F* SubtractBGFn( TH1F* hin, TF1* fbg, float xminFn, float xmaxFn ){
  
  TH1F *hout = new TH1F( *hin );
  for( int ix=0; ix<hout->GetNbinsX(); ix++ ){
    
    float fy = hin->GetBinContent(ix+1);
    float fyerr = hin->GetBinError(ix+1);
    float fyw= hin->Integral(ix+1,ix+1,"width");
    
    float xL = hin->GetBinLowEdge(ix+1);
    float xR = hin->GetBinLowEdge(ix+2);
    
    float sg = 0.0;
    float sgerr = 0.0;
    if( xL>=xminFn && xR<xmaxFn ){
      float bgFn = fbg->Integral(xL,xR);
      //      float bg = bgFn*fy/fyw;
      float bg = bgFn/(float)hin->GetBinWidth(ix+1); // zyj
      float bgerr = sqrt(abs(bg)); //simple estimation
      
      sg = fy - bg;
      sgerr = sqrt(fyerr*fyerr + bgerr*bgerr);
      if( sg<0 ) sg = 0;
    }
    hout->SetBinContent( ix+1, sg );
    hout->SetBinError( ix+1, sgerr );
  }
  
  return hout;
}

void plot_Kstar0InvMass_simpleBWfit(const bool isKstar0 = 1) {
    //EPD + TPC alignment along Z direction
    //    Inner      EPD       Outer                         TPC
    //-5.3 A -4.05 B -3.3 C -2.9 D -2.6         -2.0 A -1.1 eta-gap -1.0 B 0
    
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetPadLeftMargin(0.16);
    gStyle->SetPadRightMargin(0.03);
    gStyle->SetPadBottomMargin(0.18);
    gStyle->SetPadTopMargin(0.05);
    
    TString cutName = "Mix";
    TString particle = "Kstar0";
    TString histName = "NumKstar0InvMassvsPtY";
    
    if(!isKstar0) { particle = "AntiKstar0"; histName = "NumAntiKstar0InvMassvsPtY"; }
    
    //TFile *pol_file   = TFile::Open("kstar0_3gev_TPCorTOF_spinYmp5to0_flowYmp8to0_ana_hist.root");
//    TFile *pol_file   = TFile::Open("kstar0_3gev_TPCandTOF_spinYmp5to0_flowYmp8to0_ana_hist.root");
//    TFile *pol_file   = TFile::Open("19154032_5000010_kstar0Tree_TPCorTOF_hist_Aug12.root");
//    TFile *pol_file   = TFile::Open("Kstar0_20210823_embedding.root");
    TFile *pol_file   = TFile::Open("../../Kstar0_realData_binning1_20210901.root");
//    TFile *pol_file   = TFile::Open("../../Kstar0_realData_binning1_20210831.root");
//    TFile *pol_file   = TFile::Open("../../Kstar0_20210824_RealData_TPCandTOF.root");
    
    const int nCent = 4;
    TString centMergedName[] = {"0-10%", "10-40%", "40-80%", "0-80%"};
    TString ptName[]   = {"0.2<p_{T}<0.4 GeV/c", "0.4<p_{T}<0.8 GeV/c", "0.8<p_{T}<1.2 GeV/c", "1.2<p_{T}<1.6 GeV/c", "1.6<p_{T}<2.0 GeV/c", "2.0<p_{T}<3.0 GeV/c"};
    TString etaName[]  = {"-1.0<#eta<-0.8", "-0.8<#eta<-0.6", "-0.6<#eta<-0.4", "-0.4<#eta<-0.2", "-0.2<#eta<0", "0<#eta<0.2", "0.2<#eta<0.4", "0.4<#eta<1.0"};
    TString phiName[]  = { "-1", "-2/3", "-1/3", "0", "1/3", "2/3", "1"}; //unit: TMath::Pi()
    
    TH1F *ax;
    TF1 *zero = new TF1("zero","[0]",-10,100);
    zero->SetLineStyle(2);
    zero->SetLineWidth(1);
    zero->SetLineColor(1);
    
    TLatex lt;
    lt.SetNDC();
    lt.SetTextFont(43);
    lt.SetTextSize(14);
    
    //invariant mass spectrum of lambda
    TH1F     *invMassVsCent_tot[nCent];
    TH1F     *invMassVsCent_bg[nCent];
        
    TH3F *h3InvMassXiYvsPt_tot[nCent];
    TH3F *h3InvMassXiYvsPt_bg[nCent];
    
    //0-10%
    h3InvMassXiYvsPt_tot[0] = (TH3F*)  pol_file->Get(Form("%s_Cent0", histName.Data())); //0-5%
    h3InvMassXiYvsPt_tot[0]->Add( (TH3F*)pol_file->Get(Form("%s_Cent1", histName.Data())) ); //5-10%
    
    h3InvMassXiYvsPt_bg[0]  = (TH3F*)  pol_file->Get(Form("%s_%s_Cent0", histName.Data(), cutName.Data())); //0-5%
    h3InvMassXiYvsPt_bg[0]->Add( (TH3F*)pol_file->Get(Form("%s_%s_Cent1", histName.Data(), cutName.Data())) ); //5-10%
    
    //10-40%
    h3InvMassXiYvsPt_tot[1] = (TH3F*)  pol_file->Get(Form("%s_Cent2", histName.Data())); //10-20%
    h3InvMassXiYvsPt_tot[1]->Add( (TH3F*)pol_file->Get(Form("%s_Cent3", histName.Data())) ); //20-30%
    h3InvMassXiYvsPt_tot[1]->Add( (TH3F*)pol_file->Get(Form("%s_Cent4", histName.Data())) ); //30-40%
    
    h3InvMassXiYvsPt_bg[1]  = (TH3F*)  pol_file->Get(Form("%s_%s_Cent2", histName.Data(), cutName.Data())); //10-20%
    h3InvMassXiYvsPt_bg[1]->Add( (TH3F*)pol_file->Get(Form("%s_%s_Cent3", histName.Data(), cutName.Data())) ); //20-30%
    h3InvMassXiYvsPt_bg[1]->Add( (TH3F*)pol_file->Get(Form("%s_%s_Cent4", histName.Data(), cutName.Data())) ); //30-40%
    
    //40-80%
    h3InvMassXiYvsPt_tot[2] = (TH3F*)    pol_file->Get(Form("%s_Cent5", histName.Data())); //40-50%
    h3InvMassXiYvsPt_tot[2]->Add( (TH3F*)pol_file->Get(Form("%s_Cent6", histName.Data())) ); //50-60%
    h3InvMassXiYvsPt_tot[2]->Add( (TH3F*)pol_file->Get(Form("%s_Cent7", histName.Data())) ); //60-70%
    h3InvMassXiYvsPt_tot[2]->Add( (TH3F*)pol_file->Get(Form("%s_Cent8", histName.Data())) ); //70-80%
    
    h3InvMassXiYvsPt_bg[2]  = (TH3F*)   pol_file->Get(Form("%s_%s_Cent5", histName.Data(), cutName.Data())); //40-50%
    h3InvMassXiYvsPt_bg[2]->Add( (TH3F*)pol_file->Get(Form("%s_%s_Cent6", histName.Data(), cutName.Data())) ); //50-60%
    h3InvMassXiYvsPt_bg[2]->Add( (TH3F*)pol_file->Get(Form("%s_%s_Cent7", histName.Data(), cutName.Data())) ); //60-70%
    h3InvMassXiYvsPt_bg[2]->Add( (TH3F*)pol_file->Get(Form("%s_%s_Cent8", histName.Data(), cutName.Data())) ); //70-80%
    
    //0-80%
    h3InvMassXiYvsPt_tot[3] = (TH3F*)    pol_file->Get(Form("%s_Cent7", histName.Data())); //60-70%
    h3InvMassXiYvsPt_tot[3]->Add( (TH3F*)pol_file->Get(Form("%s_Cent8", histName.Data())) ); //70-80%
    h3InvMassXiYvsPt_tot[3]->Add( h3InvMassXiYvsPt_tot[0] ); //0-10%
    h3InvMassXiYvsPt_tot[3]->Add( h3InvMassXiYvsPt_tot[1] ); //10-40%
    h3InvMassXiYvsPt_tot[3]->Add( h3InvMassXiYvsPt_tot[2] ); //40-60%
    
    h3InvMassXiYvsPt_bg[3]  = (TH3F*)   pol_file->Get(Form("%s_%s_Cent7", histName.Data(), cutName.Data())); //60-70%
    h3InvMassXiYvsPt_bg[3]->Add( (TH3F*)pol_file->Get(Form("%s_%s_Cent8", histName.Data(), cutName.Data())) ); //70-80%
    h3InvMassXiYvsPt_bg[3]->Add( h3InvMassXiYvsPt_bg[0] ); //0-10%
    h3InvMassXiYvsPt_bg[3]->Add( h3InvMassXiYvsPt_bg[1] ); //10-40%
    h3InvMassXiYvsPt_bg[3]->Add( h3InvMassXiYvsPt_bg[2] ); //40-60%
    
    int lowEtaBin = ((TH1F*)h3InvMassXiYvsPt_tot[0]->ProjectionY())->FindBin(-0.5), highEtaBin = ((TH1F*)h3InvMassXiYvsPt_tot[0]->ProjectionY())->FindBin(0); //-0.5 < y < 0
    int ptBinIdx_Low  = ((TH1F*)h3InvMassXiYvsPt_tot[0]->ProjectionX())->FindBin(0.4), ptBinIdx_High = ((TH1F*)h3InvMassXiYvsPt_tot[0]->ProjectionX())->FindBin(2.0); // 0.4 < pT < 2.0
    for(int icent=0; icent<nCent; icent++) {
        invMassVsCent_tot[icent]   = (TH1F*)h3InvMassXiYvsPt_tot[icent]->ProjectionZ(Form("sig_Cent%d", icent), ptBinIdx_Low, ptBinIdx_High, lowEtaBin, highEtaBin);
        
        invMassVsCent_bg[icent]    = (TH1F*)h3InvMassXiYvsPt_bg[icent]->ProjectionZ(Form("bg_Cent%d", icent), ptBinIdx_Low, ptBinIdx_High, lowEtaBin, highEtaBin);
    }
    /*
    TCanvas *ca_Accep = new TCanvas("ca_Accep","tot. sig vs comb. bkg", 1200, 600);
    ca_Accep->Divide(2,2);
    
    ca_Accep->cd(1)->SetGrid();
    TH2F *hKstar0RapidityvsPt_sig = (TH2F*)hNumKstar0InvMassvsPtY_sig->Project3D("xy");
    hKstar0RapidityvsPt_sig->Draw("COLZ");
    hKstar0RapidityvsPt_sig->GetXaxis()->SetTitle("Rapidity (y_{C.M.})");
    hKstar0RapidityvsPt_sig->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    hKstar0RapidityvsPt_sig->GetXaxis()->SetRangeUser(-1., 1.);
    hKstar0RapidityvsPt_sig->GetYaxis()->SetRangeUser(0, 3.0);
    
    ca_Accep->cd(3);
    TH1F *hInvMassAcceptance_sig = (TH1F*)hNumKstar0InvMassvsPtY_sig->ProjectionZ("invM_sig",9,60,1,10); //0.4<pT<3.0, -1.0<y<0
    hInvMassAcceptance_sig->Draw("ESAME");
    hInvMassAcceptance_sig->SetMarkerStyle(20);
    hInvMassAcceptance_sig->SetMarkerSize(0.5);
    hInvMassAcceptance_sig->SetMarkerColor(kBlack);
    hInvMassAcceptance_sig->GetXaxis()->SetTitle("M_{inv} [GeV/c^{2}]");
    hInvMassAcceptance_sig->GetYaxis()->SetTitle("Counts");
    
    ca_Accep->cd(2)->SetGrid();
    TH3F *hNumKstar0InvMassvsPtY_bg = (TH3F*)pol_file->Get(Form("NumKstar0InvMassvsPtY_%s_Cent9",cutName.Data()));
    TH2F *hKstar0RapidityvsPt_bg = (TH2F*)hNumKstar0InvMassvsPtY_bg->Project3D("xy");
    hKstar0RapidityvsPt_bg->Draw("COLZ");
    hKstar0RapidityvsPt_bg->GetXaxis()->SetTitle("Rapidity (y_{C.M.})");
    hKstar0RapidityvsPt_bg->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    hKstar0RapidityvsPt_bg->GetXaxis()->SetRangeUser(-1., 1.);
    hKstar0RapidityvsPt_bg->GetYaxis()->SetRangeUser(0, 3.0);
    
    ca_Accep->cd(4);
    TH1F *hInvMassAcceptance_bg = (TH1F*)hNumKstar0InvMassvsPtY_bg->ProjectionZ("invM_bg",9,60,1,10); //0.4<pT<3.0, -1.0<y<0
    hInvMassAcceptance_bg->Draw("ESAME");
    hInvMassAcceptance_bg->SetMarkerStyle(20);
    hInvMassAcceptance_bg->SetMarkerSize(0.5);
    hInvMassAcceptance_bg->SetMarkerColor(kBlue);
    hInvMassAcceptance_bg->GetXaxis()->SetTitle("M_{inv} [GeV/c^{2}]");
    hInvMassAcceptance_bg->GetYaxis()->SetTitle("Counts");
    
    ca_Accep->cd();
    ca_Accep->Print(Form("%s/Acceptance_%s_fxt3GeV.png", cutName.Data(), particle.Data()));
    ca_Accep->Print(Form("%s/Acceptance_%s_fxt3GeV.pdf", cutName.Data(), particle.Data()));
    */
    
    //plotting
    TF1 *f1_cent;
    TF1 *f2_bg;
    TF1 *f2_sig;
    
    TH1F *h1_cent_tot;
    TH1F *h1_cent_bkg;
    TH1F *h1_cent_sig;
    
    TH1F *hKstar0Mass_sig[nCent];
    TH1F *hKstar0MassReal_sig[nCent];
    float scaler[nCent];
    TLine *lineMass[nCent];
    TLine *lineL[nCent];
    TLine *lineR[nCent];
    float maxCounts=0, maxCountsErr=0;
    float upperCounts[] = {32000, 15000};
    TCanvas *ca_invMass = new TCanvas("ca_invMass","invMass of Kstar0", 450*nCent, 800);
    ca_invMass->Divide(nCent,2);
    
    gPad->SetTopMargin(0.01);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gPad->SetLeftMargin(0.1);
    
    double lowFitRange  = 0.75;
    double highFitRange = 1.20;
    
    for(int icent=0; icent<nCent; icent++) {
        ca_invMass->cd(icent+1);
        
        double nInclKstar0_r = invMassVsCent_tot[icent]->Integral(invMassVsCent_tot[icent]->FindBin(1.2), invMassVsCent_tot[icent]->FindBin(1.5));
        double nBgKstar0_r   = invMassVsCent_bg[icent]->Integral(invMassVsCent_bg[icent]->FindBin(1.2), invMassVsCent_bg[icent]->FindBin(1.5));
        double nInclKstar0_l = invMassVsCent_tot[icent]->Integral(invMassVsCent_tot[icent]->FindBin(0.6), invMassVsCent_tot[icent]->FindBin(0.8));
        double nBgKstar0_l   = invMassVsCent_bg[icent]->Integral(invMassVsCent_bg[icent]->FindBin(0.6), invMassVsCent_bg[icent]->FindBin(0.8));
        //scaler[icent] = nInclKstar0_r / nBgKstar0_r;
        scaler[icent] = (nInclKstar0_r + nInclKstar0_l) / (nBgKstar0_r + nBgKstar0_l);
        invMassVsCent_bg[icent]->Scale(scaler[icent]);
        
        maxCounts       = invMassVsCent_tot[icent]->GetBinContent(invMassVsCent_tot[icent]->GetMaximumBin());
        maxCountsErr    = invMassVsCent_tot[icent]->GetBinError(invMassVsCent_tot[icent]->GetMaximumBin());
        
        ax = gPad->DrawFrame( 0.6, -0.05*(maxCounts + maxCountsErr), 1.5, 1.05*(maxCounts + maxCountsErr) ); // no EPcorr
        //SetAxis(ax, 3.0, 3.0);
        ax->SetXTitle("M_{inv} [GeV/c^{2}]");
        ax->SetYTitle("Counts");
        ax->Draw();
        
        invMassVsCent_tot[icent]->Draw("E1SAME");
        invMassVsCent_tot[icent]->SetMarkerStyle(20);
        invMassVsCent_tot[icent]->SetMarkerSize(0.5);
        invMassVsCent_tot[icent]->SetMarkerColor(kBlack);
        
        invMassVsCent_bg[icent]->Draw("ESAME");
        invMassVsCent_bg[icent]->SetMarkerStyle(24);
        invMassVsCent_bg[icent]->SetMarkerSize(0.5);
        invMassVsCent_bg[icent]->SetMarkerColor(kBlue);
        
        TLegend *legT;
        legT = new TLegend(0.4, 0.82, 0.85, 0.9);
        legT->SetFillColor(10);
        legT->SetLineStyle(3004);
        legT->SetLineColor(10);
        legT->SetLineWidth(0.);
        legT->SetTextFont(42);
        legT->SetTextSize(0.04);
        if(isKstar0)
            legT->SetHeader(Form("K^{*0}, %s, 0.4<p_{T}<2.0, -0.5<y<0", centMergedName[icent].Data()), "p");
        else
            legT->SetHeader(Form("#bar{K}^{*0}, %s, 0.4<p_{T}<3.0, -0.5<y<0", centMergedName[icent].Data()), "p");
        
        //legT->AddEntry(invMassVsCent_tot[icent], "total", "lep");
        //legT->AddEntry(invMassVsCent_bg[icent], "Comb. BG (mixing #pi)", "lep");
        legT->Draw();
        
        TLegend *legT_t;
        legT_t = new TLegend(0.6, 0.66, 0.85, 0.82);
        legT_t->SetFillColor(10);
        legT_t->SetLineStyle(3004);
        legT_t->SetLineColor(10);
        legT_t->SetLineWidth(0.);
        legT_t->SetTextFont(42);
        legT_t->SetTextSize(0.04);
        legT_t->AddEntry(invMassVsCent_tot[icent], "Same Event (SE)", "lep");
        legT_t->AddEntry(invMassVsCent_bg[icent], "Mixed Event (ME)", "lep");
        legT_t->Draw();
        
        ca_invMass->cd(icent+nCent+1);
        hKstar0Mass_sig[icent] = (TH1F*)invMassVsCent_tot[icent]->Clone();
        hKstar0Mass_sig[icent]->Add(invMassVsCent_bg[icent], -1.0);
        
        float lowerEdgeY, upperEdgeY;
        float minCounts, minCountsErr;
        float minTemp = 0., maxTemp = 0.;
        int invMBinIdx_Low  = hKstar0Mass_sig[icent]->GetXaxis()->FindBin(lowFitRange);
        int invMBinIdx_High = hKstar0Mass_sig[icent]->GetXaxis()->FindBin(highFitRange);
        for(int iinvM=invMBinIdx_Low; iinvM<invMBinIdx_High; iinvM++) {
            maxCounts    = hKstar0Mass_sig[icent]->GetBinContent(iinvM);
            maxCountsErr = hKstar0Mass_sig[icent]->GetBinError(iinvM);
            if( (maxCounts+maxCountsErr)>maxTemp )  maxTemp = maxCounts+maxCountsErr;
            
            minCounts    = hKstar0Mass_sig[icent]->GetBinContent(iinvM);
            minCountsErr = hKstar0Mass_sig[icent]->GetBinError(iinvM);
            if( (minCounts+minCountsErr)<minTemp )  minTemp = minCounts+minCountsErr;
        }
        
        if(minTemp>=0) lowerEdgeY = -0.15*maxTemp;
        else lowerEdgeY = 1.15*minTemp;
        if(maxTemp>0) upperEdgeY =  1.15*maxTemp;
        else upperEdgeY = maxTemp - minTemp;
        
        ax = gPad->DrawFrame( lowFitRange, lowerEdgeY, highFitRange, upperEdgeY ); // no EPcorr
        ax->SetXTitle("M_{inv} [GeV/c^{2}] ");
        ax->SetYTitle("Counts");
        ax->Draw();
        
        hKstar0Mass_sig[icent]->Draw("ESAME");
        hKstar0Mass_sig[icent]->GetXaxis()->SetRangeUser(lowFitRange, highFitRange);
        hKstar0Mass_sig[icent]->GetYaxis()->SetRangeUser(lowerEdgeY, upperEdgeY);
        hKstar0Mass_sig[icent]->GetXaxis()->SetTitle("M_{inv} [GeV/c^{2}]");
        hKstar0Mass_sig[icent]->GetYaxis()->SetTitle("Counts");
        hKstar0Mass_sig[icent]->SetMarkerStyle(20);
        hKstar0Mass_sig[icent]->SetMarkerSize(0.5);
        hKstar0Mass_sig[icent]->SetMarkerColor(kRed);
        
        TLine *unitaty = new TLine(lowFitRange, 0., highFitRange, 0.);
        unitaty->Draw();
        
        f1_cent = new TF1("f1_cent", fitFuncbw, lowFitRange, highFitRange, 6);//7); //order
        f1_cent->SetParNames("Scale1", "Width", "Mean","p0", "p1", "p2");//, "p3");
        //f1_cent->SetParLimits(0, 0, 50000);
        f1_cent->SetParLimits(1, 0.03, 0.07);
        f1_cent->SetParLimits(2, 0.85, 0.93);
        f1_cent->SetLineColor(kBlack);
        //f1_cent->SetNpx(hKstar0Mass_sig[i]->GetNbinsX());
        //f1_cent->SetNpx(45);   // 0.45/0.01 = 45
        
        f2_bg = new TF1("f2_bg", background, lowFitRange, highFitRange, 3);//4); //order
        f2_bg->SetLineColor(kBlue);
        f2_bg->SetLineStyle(2);
        //f2_bg->SetNpx(hKstar0Mass_sig[i]->GetNbinsX());
        //f2_bg->SetNpx(45);
        
        hKstar0Mass_sig[icent]->Fit("f1_cent", "MO0", "", lowFitRange, highFitRange);
        f1_cent->SetParameters(f1_cent->GetParameter("Scale1"),
                f1_cent->GetParameter("Width"),
                f1_cent->GetParameter("Mean"),
                f1_cent->GetParameter("p0"),
                f1_cent->GetParameter("p1"),
                f1_cent->GetParameter("p2")//,
                //f1_cent->GetParameter("p3")
                );
        hKstar0Mass_sig[icent]->Fit("f1_cent", "MO0", "", lowFitRange, highFitRange);
        
        f1_cent->Draw("SAME");

        double *par;
        const double *parErr;
        par = f1_cent->GetParameters();
        parErr = f1_cent->GetParErrors();
        f2_bg->SetParameters(&par[3]);
        f2_bg->Draw("SAME");
        
        f2_sig = new TF1("f2_sig", breitwigner, lowFitRange, highFitRange, 3); //order
        f2_sig->SetLineColor(kRed);
        f2_sig->SetLineStyle(2);
        //f2_sig->SetNpx(hKstar0Mass_sig[i]->GetNbinsX());
        //f2_sig->SetNpx(45);
        f2_sig->SetParameters(&par[0]);
        f2_sig->Draw("SAME");
        
        h1_cent_tot = (TH1F *)f1_cent->GetHistogram()->Clone();
        h1_cent_sig = (TH1F *)f1_cent->GetHistogram()->Clone();
        h1_cent_bkg = (TH1F *)f2_bg->GetHistogram()->Clone();
        h1_cent_sig->Add(h1_cent_bkg, -1.0);
        //h1_cent_sig->Draw("SAME");
        //h1_cent_sig->SetLineColor(kRed);
        
        hKstar0MassReal_sig[icent] = SubtractBGFn(hKstar0Mass_sig[icent], f2_bg, lowFitRange, highFitRange);
        hKstar0MassReal_sig[icent]->SetLineWidth(2.0);
        hKstar0MassReal_sig[icent]->SetLineColor(kRed);
        hKstar0MassReal_sig[icent]->Draw("HIST SAME");
        
        int lowerMassBin_Kstar0 = h1_cent_tot->FindBin(par[2]-2.0*par[1]);
        int upperMassBin_Kstar0 = h1_cent_tot->FindBin(par[2]+2.0*par[1])-1;
        //float nSignal = h1_cent_sig->Integral(lowerMassBin_Kstar0, upperMassBin_Kstar0);
        //float nTotal  = h1_cent_tot->Integral(lowerMassBin_Kstar0, upperMassBin_Kstar0);
        //float nBackground  = h1_cent_bkg->Integral(lowerMassBin_Kstar0, upperMassBin_Kstar0);
        
        lineMass[icent] = new TLine(par[2], lowerEdgeY, par[2], upperEdgeY);
        lineMass[icent]->SetLineColor(kGreen);
        lineMass[icent]->SetLineStyle(2);
        lineMass[icent]->SetLineWidth(2);
        lineMass[icent]->Draw("SAME");
        
        lineL[icent] = new TLine(par[2]-2.0*par[1], lowerEdgeY, par[2]-2.0*par[1], upperEdgeY);
        lineL[icent]->SetLineColor(kBlack);
        lineL[icent]->SetLineStyle(2);
        lineL[icent]->SetLineWidth(2);
        lineL[icent]->Draw("SAME");
        
        lineR[icent] = new TLine(par[2]+2.0*par[1],lowerEdgeY, par[2]+2.0*par[1], upperEdgeY);
        lineR[icent]->SetLineColor(kBlack);
        lineR[icent]->SetLineStyle(2);
        lineR[icent]->SetLineWidth(2);
        lineR[icent]->Draw("SAME");
        
        float nSignal = hKstar0MassReal_sig[icent]->Integral(lowerMassBin_Kstar0, upperMassBin_Kstar0);
        float nTotal  = invMassVsCent_tot[icent]->Integral(lowerMassBin_Kstar0, upperMassBin_Kstar0);
        float nBackground  = nTotal - nSignal; //invMassVsCent_bg[icent]->Integral(lowerMassBin_Kstar0, upperMassBin_Kstar0);
        cout << nTotal << "; " << nSignal << "; " << nBackground <<  endl;
        
        TLatex *tex_chi2 = new TLatex(0.76, 0.45*upperEdgeY, Form("#chi^{2}/NDF = %4.2f / %d", f1_cent->GetChisquare(), f1_cent->GetNDF()));
        tex_chi2->SetTextFont(42);
        tex_chi2->SetTextSize(0.04);
        tex_chi2->Draw("same");
        
        TLatex *tex_mean = new TLatex(0.76, 0.35*upperEdgeY, Form("m_{0} = %4.3f #pm %4.3f GeV/c^{2}", par[2], parErr[2]));
        tex_mean->SetTextFont(42);
        tex_mean->SetTextSize(0.04);
        tex_mean->Draw("same");
        
        TLatex *tex_gamma = new TLatex(0.76, 0.25*upperEdgeY, Form("#Gamma = %4.3f #pm %4.3f GeV/c^{2}", par[1], parErr[1]));
        tex_gamma->SetTextFont(42);
        tex_gamma->SetTextSize(0.04);
        tex_gamma->Draw("same");
        
        TLatex *tex_signif = new TLatex(0.76, 0.15*upperEdgeY, Form("S / #sqrt{S+B} = %6.0f / %4.0f = %6.2f", nSignal, sqrt(nTotal), nSignal/sqrt(nTotal)));
        tex_signif->SetTextFont(42);
        tex_signif->SetTextSize(0.04);
        tex_signif->Draw("same");
        
        TLatex *tex_SBR = new TLatex(0.76, 0.05*upperEdgeY, Form("S / B = %6.5f", nSignal/(nBackground)));
        tex_SBR->SetTextFont(42);
        tex_SBR->SetTextSize(0.04);
        tex_SBR->Draw("same");
        
        TLegend *legT1;
        legT1 = new TLegend(0.57, 0.6, 0.92, 0.9);
        legT1->SetFillColor(10);
        legT1->SetLineStyle(3004);
        legT1->SetLineColor(10);
        legT1->SetLineWidth(0.);
        legT1->SetTextFont(42);
        legT1->SetTextSize(0.04);
        if(isKstar0)
            legT->SetHeader(Form("K^{*0}, %s, 0.4<p_{T}<2.0, -0.5<y<0", centMergedName[icent].Data()), "p");
        else
            legT->SetHeader(Form("#bar{K}^{*0}, %s, 0.4<p_{T}<2.0, -0.5<y<0", centMergedName[icent].Data()), "p");
        legT1->AddEntry(hKstar0Mass_sig[icent], "SE - ME", "lep");
        legT1->AddEntry(f1_cent, "BW + pol2", "lp");
        legT1->AddEntry(f2_bg, "residual background", "lp");
        legT1->AddEntry(hKstar0MassReal_sig[icent], "K^{*0} signal", "l");
        legT1->Draw();
    }
    
    ca_invMass->cd();
    ca_invMass->Print(Form("%s/invMass%s_SimpleBWfit_fxt3GeV.png", cutName.Data(), particle.Data()));
    ca_invMass->Print(Form("%s/invMass%s_SimpleBWfit_fxt3GeV.pdf", cutName.Data(), particle.Data()));
}
