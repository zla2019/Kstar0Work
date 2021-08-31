#include <math.h>
#include "TMath.h"
#include <cmath>
#include "TMinuit.h"
using namespace TMath;

float ximass = 0.497611;


TLatex* drawLatex(Double_t x, Double_t y, char* text, Int_t textFont, Double_t textSize, Int_t colorIndex)
{
  TLatex *latex = new TLatex(x, y, text);
  latex->SetNDC();
  latex->SetTextFont(textFont);
  latex->SetTextSize(textSize);
  latex->SetTextColor(colorIndex);
  latex->Draw("same");
  return latex;
}


TF1 *fsin;
TF1 *fLastFunc;

double gaus(Double_t *x,Double_t *par) // Gaus function
{
  Double_t fitsig = (par[0]/(par[2]*sqrt(2*3.14159265358979323846)))*TMath::Exp(-0.5* ((x[0]-par[1])/par[2])* ((x[0]-par[1])/par[2]) ) ;
  return fitsig;
}

double doublegaus(Double_t *x,Double_t *par)
{
        Double_t fitsig = par[0]*exp(-0.5*(((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]))) + par[3]*exp(-0.5*(((x[0]-par[1])/par[4])*((x[0]-par[1])/par[4])));
        return fitsig;
}

double coshy(Double_t *x,Double_t *par) //cosh(y)
{
  Double_t teff = par[0]/TMath::CosH(x[0]);
  return teff;
}
/*
double background(Double_t *x,Double_t *par) //pol1
{
  Double_t fitbg = par[0]+par[1]*x[0];
  return fitbg;
}
*/
double background(Double_t *x,Double_t *par) //pol2
{
  Double_t fitbg = par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
  return fitbg;
}

double fline(Double_t *x, Double_t *par)
{
  if (x[0] > 0.75 && x[0] < 1.05) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0];
}

double fitFuncGaus(Double_t *x,Double_t *par) //gaus + pol2
{
  Double_t fitval = background(x,&par[3])  + gaus(x,par);
  return fitval;
}

double fitFuncDoubleGaus(Double_t *x,Double_t *par) //bw + pol3
{
    Double_t fitval = background(x,&par[5])  + doublegaus(x,par);
    return fitval;
}

double breitwigner(Double_t *x,Double_t *par) // Breit-Wigner function
{
    Double_t fitsig = (par[0]*par[1]/(2.0*3.14159265358979323846)) / (par[1]*par[1]/4.0 + (x[0] - par[2])*(x[0] - par[2])) ;
    return fitsig;
}

double fitFuncbw(Double_t *x,Double_t *par) //bw + pol2
{
    Double_t fitval = background(x,&par[3])  + breitwigner(x,par);
    return fitval;
}

double fitFuncExp(Double_t *x,Double_t *par) //mT exponential
{
  Double_t fitdNdy = par[0]*exp(-(sqrt(x[0]*x[0] + 0.497611*0.497611)-0.497611)/par[1])*x[0];
//  Double_t fitdNdy = par[0]*exp(-x[0]/par[1]);
  return fitdNdy;
}

double fitFuncLevy(Double_t *x,Double_t *par) //Levy function
{
  Double_t fitdNdy = par[0]*pow(1.0 + x[0]/(par[1]*par[2]), -par[2]);
  return fitdNdy;
}

double fitFuncBoltzmann(Double_t *x,Double_t *par) //Boltzmann function
{
  Double_t fitdNdy = par[0]*(x[0] + 0.89594)*exp(-(x[0] + 0.89594)/par[1]);
  return fitdNdy;
}

//Taken from https://gitlab.cern.ch/mkrzewic/AliRoot/-/blob/fef6f2b32ddd823bd1c11b130f321a585288b032/PWG/Tools/AliPWGFunc.cxx
Double_t IntegrandBG(const double * x, const double* p){
  // integrand for boltzman-gibbs blast wave
  // x[0] -> r (radius)
  // p[0] -> mass
  // p[1] -> pT (transverse momentum)
  // p[2] -> beta_max (surface velocity)
  // p[3] -> T (freezout temperature)
  // p[4] -> n (velocity profile)
  
  double x0 = x[0];
  
  double mass     = p[0];
  double pT       = p[1];
  double beta_max = p[2];
  double temp     = p[3];
  Double_t n      = p[4];
  
  // Keep beta within reasonable limits
  Double_t beta = beta_max * TMath::Power(x0, n);
  if (beta > 0.9999999999999999) beta = 0.9999999999999999;
  
  double mT      = TMath::Sqrt(mass*mass+pT*pT);
  
  double rho0   = TMath::ATanH(beta);
  double arg00 = pT*TMath::SinH(rho0)/temp;
  if (arg00 > 700.) arg00 = 700.; // avoid FPE
  double arg01 = mT*TMath::CosH(rho0)/temp;
  double f0 = x0*mT*TMath::BesselI0(arg00)*TMath::BesselK1(arg01);
  
  //  printf("r=%f, pt=%f, beta_max=%f, temp=%f, n=%f, mt=%f, beta=%f, rho=%f, argI0=%f, argK1=%f\n", x0, pT, beta_max, temp, n, mT, beta, rho0, arg00, arg01);
  
  return f0;
}

Double_t StaticBGdNdPt(const double * x, const double* p) {
  
  // implementation of BGBW (1/pt dNdpt)
  
  double pT = x[0];;
  
  
  double mass    = p[0];
  double beta    = p[1];
  double temp    = p[2];
  double n       = p[3];
  double norm    = p[4];
  
  static TF1 * fIntBG = 0;
  if(!fIntBG)
    fIntBG = new TF1 ("fIntBG", IntegrandBG, 0, 1, 5);
  
  fIntBG->SetParameters(mass, pT, beta, temp,n);
  double result = fIntBG->Integral(0,1);
  //  printf ("[%4.4f], Int :%f\n", pT, result);
  return result*norm;//*1e30;;
  
}

Double_t StaticBGdNdPtTimesPt(const double * x, const double* p) {
  // BGBW dNdpt implementation
  return x[0]*StaticBGdNdPt(x,p);
}

TF1 *GetBGBWdNdpt(Double_t mass, Double_t beta, Double_t temp,
                  Double_t n, Double_t norm, const char * name){
  
  // BGBW 1/pt dNdpt
  
  fLastFunc = new TF1 (name, StaticBGdNdPt, 0.0, 10, 5);
  fLastFunc->SetParameters(mass,beta,temp,n,norm);
  fLastFunc->FixParameter(0,mass);
  
  fLastFunc->SetParNames("mass", "#beta", "T", "n", "norm");
  fLastFunc->SetLineWidth(1);
  return fLastFunc;
  
}

TF1 *GetBGBWdNdptTimesPt(Double_t mass, Double_t beta, Double_t temp,
                         Double_t n, Double_t norm, const char * name){
  
  
  fLastFunc = new TF1 (name, StaticBGdNdPtTimesPt, 0.0, 10, 5);
  fLastFunc->SetParameters(mass,beta,temp,n,norm);
  fLastFunc->FixParameter(0,mass);
  
  fLastFunc->SetParNames("mass", "#beta", "T", "n", "norm");
  fLastFunc->SetLineWidth(1);
  return fLastFunc;
  
}

TH1F* SubtractBG( TH1F* hin, TH1F* hbg){
  
  TH1F *hout = new TH1F( *hin );
  for( int ix=0; ix<hout->GetNbinsX(); ix++ ){
    
    float fy = hin->GetBinContent(ix+1);
    float bg = hbg->GetBinContent(ix+1);
    
    float fyerr = hin->GetBinError(ix+1);
    float bgerr = hbg->GetBinError(ix+1);
    
    float sg = fy - bg;
    float sgerr = sqrt(fyerr*fyerr + bgerr*bgerr);
    
    hout->SetBinContent( ix+1, sg );
    hout->SetBinError( ix+1, sgerr );
    if(fy == 0) hout->SetBinError( ix+1, 1);
  }
  
  return hout;
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

void GetRawSignal( TGraphErrors *gphi, TGraphErrors *gout, const int &fix, double &fx, double fstat[2], int &fsgn ){
  
  Double_t min, edm, errdef;
  Int_t npari, nparx, stat;
  
  gphi->Fit( fsin, "MOR", "", 0, TMath::Pi());
  gMinuit->mnstat( min, edm, errdef, npari, nparx, stat );
  
  fstat[0] = fsin->GetChisquare();
  fstat[1] = fsin->GetNDF();
  
  if( stat==3 ) {
    gout->SetPoint     ( fix, fx, fsgn*fsin->GetParameter(1) );
    gout->SetPointError( fix, 0., fsin->GetParError(1) );
  }
  else {
    cout <<"GetRawSignal():: fix="<<fix<<" fx="<<fx<<" : Fitting status is NOT good. Check the fitting result !!" << endl;
    gout->SetPoint     ( fix, fx, -999. );
  }
  
}

void SetAxis( TH1F *ax, float xff, float yoff ){
  ax->SetTitleFont(43,"X");
  ax->SetTitleSize(16,"X");
  ax->SetTitleOffset(xff,"X");
  ax->SetTitleFont(43,"Y");
  ax->SetTitleSize(16,"Y");
  ax->SetTitleOffset(yoff,"Y");
  ax->SetNdivisions(505,"X");
  ax->SetNdivisions(505,"Y");
  ax->SetLabelFont(43,"X");
  ax->SetLabelFont(43,"Y");
  ax->SetLabelSize(12,"X");
  ax->SetLabelSize(12,"Y");
}

void SetAxis( TGraphErrors *g, float xff, float yoff ){
  g->GetXaxis()->SetTitleFont(43);
  g->GetXaxis()->SetTitleSize(16);
  g->GetXaxis()->SetTitleOffset(xff);
  g->GetYaxis()->SetTitleFont(43);
  g->GetYaxis()->SetTitleSize(16);
  g->GetYaxis()->SetTitleOffset(yoff);
  g->GetXaxis()->SetNdivisions(505);
  g->GetYaxis()->SetNdivisions(505);
  g->GetXaxis()->SetLabelFont(43);
  g->GetYaxis()->SetLabelFont(43);
  g->GetXaxis()->SetLabelSize(12);
  g->GetYaxis()->SetLabelSize(12);
}

void SetLegend( TLegend *l ){
  
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->SetTextFont(43);
  l->SetTextSize(13);
}

/////////////////////////////////////////////////////////////////////////////
// Levy function for dN/ptdpt
/////////////////////////////////////////////////////////////////////////////
Double_t LevyFcn(Double_t *x, Double_t *par)// x is pT
{
  Double_t A  = par[0];
  Double_t n  = par[1];
  Double_t T  = par[2];
  Double_t m0 = par[3];
  Double_t mT = sqrt(x[0] * x[0] + m0 * m0);
  
  Double_t a1 = A * (n - 1) * (n - 2);
  Double_t a2 = n * T * (n * T + m0 * (n - 2));
  Double_t a3 = pow(1 + (mT - m0) / n / T, -n);
  
  return  a1 / a2 * a3;
}

/////////////////////////////////////////////////////////////////////////////
// Levy function for dN/dpt
/////////////////////////////////////////////////////////////////////////////
Double_t LevyFcnPt(Double_t *x, Double_t *par)// x is pT
{
  Double_t A  = par[0];
  Double_t n  = par[1];
  Double_t T  = par[2];
  Double_t m0 = par[3];
  Double_t mT = sqrt(x[0] * x[0] + m0 * m0);
  
  Double_t a1 = A * (n - 1) * (n - 2);
  Double_t a2 = n * T * (n * T + m0 * (n - 2));
  Double_t a3 = pow(1 + (mT - m0) / n / T, -n);
  
  return  a1 / a2 * a3 * x[0];
}
/////////////////////////////////////////////////////////////////////////////
// Levy function for dN/dpt * pT
/////////////////////////////////////////////////////////////////////////////
Double_t LevyFcnPt2(Double_t *x, Double_t *par)//x is pT
{
  Double_t A  = par[0];
  Double_t n  = par[1];
  Double_t T  = par[2];
  Double_t m0 = par[3];
  Double_t mT = sqrt(x[0] * x[0] + m0 * m0);
  
  Double_t a1 = A * (n - 1) * (n - 2);
  Double_t a2 = n * T * (n * T + m0 * (n - 2));
  Double_t a3 = pow(1 + (mT - m0) / n / T, -n);
  
  return  a1 / a2 * a3 * x[0] * x[0];
}
/////////////////////////////////////////////////////////////////////////////
// powerlaw function for dN/ptdpt // dN/mTdmT
/////////////////////////////////////////////////////////////////////////////
Double_t powerlaw(Double_t *x, Double_t *par)// x is mT-m0
{
  const Double_t m0 = ximass;
  double mTm0 = x[0];
  double pT = sqrt((mTm0 + m0) * (mTm0 + m0) - m0 * m0);
  double dNdy = par[0];
  double meanPt = par[1];
  double n = par[2];
  double p0 = meanPt * (n - 3.) / 2.;
  double A = dNdy * 4.*(n - 1) * (n - 2) / (n - 3) / (n - 3) / meanPt / meanPt;
  return A * TMath::Power(1 + pT / p0, -n);
}

/////////////////////////////////////////////////////////////////////////////
// powerlaw function for dN/dpt // dN/dmT
/////////////////////////////////////////////////////////////////////////////
Double_t powerlawmT(Double_t *x, Double_t *par)// x is mT// y is * mT
{
  const Double_t m0 = ximass;
  double mTm0 = x[0] - m0;
  double mT   = x[0];
  double pT = sqrt(mT * mT - m0 * m0);
  double dNdy = par[0];
  double meanPt = par[1];
  double n = par[2];
  double p0 = meanPt * (n - 3.) / 2.;
  double A = dNdy * 4.*(n - 1) * (n - 2) / (n - 3) / (n - 3) / meanPt / meanPt;
  return A * TMath::Power(1 + pT / p0, -n) * mT;
}
/////////////////////////////////////////////////////////////////////////////
// exponential function for dN/ptdpt // dN/mTdmT // mTexponential
/////////////////////////////////////////////////////////////////////////////
Double_t exponential(Double_t *x, Double_t *par) // x is mT-m0
{
  return par[0] * exp(-x[0] / par[1]);
}

/////////////////////////////////////////////////////////////////////////////
// exponential*mT function for dN/dmT // mTexponential
/////////////////////////////////////////////////////////////////////////////
Double_t exponentialmT(Double_t *x, Double_t *par)// x is mT //y is *mT
{
  const Double_t m0 = ximass;
  return par[0] * exp(-(x[0] - m0) / par[1]) * x[0];
}

// // Next is dN/pTdpT and so on, changed to pT
// // x[0]  -> sqrt(m0*m0+x[0]*x[0]) - m0
/////////////////////////////////////////////////////////////////////////////
// powerlaw function for dN/ptdpt
/////////////////////////////////////////////////////////////////////////////
Double_t powerlaw2(Double_t *x, Double_t *par)// x is pT
{
  const Double_t m0 = ximass;
  double mTm0 = sqrt(x[0] * x[0] + m0 * m0) - m0;
  double pT = sqrt((mTm0 + m0) * (mTm0 + m0) - m0 * m0);
  double dNdy = par[0];
  double meanPt = par[1];
  double n = par[2];
  double p0 = meanPt * (n - 3.) / 2.;
  double A = dNdy * 4.*(n - 1) * (n - 2) / (n - 3) / (n - 3) / meanPt / meanPt;
  return A * TMath::Power(1 + pT / p0, -n);
}

/////////////////////////////////////////////////////////////////////////////
// powerlaw function for dN/dpt
/////////////////////////////////////////////////////////////////////////////
Double_t powerlawmT2(Double_t *x, Double_t *par) // x is pT // y is * pT
{
  const Double_t m0 = ximass;
  double mTm0 = sqrt(x[0] * x[0] + m0 * m0) - m0;
  double mT   = mTm0 + m0;
  double pT = sqrt((mTm0 + m0) * (mTm0 + m0) - m0 * m0);
  double dNdy = par[0];
  double meanPt = par[1];
  double n = par[2];
  double p0 = meanPt * (n - 3.) / 2.;
  double A = dNdy * 4.*(n - 1) * (n - 2) / (n - 3) / (n - 3) / meanPt / meanPt;
  return A * TMath::Power(1 + pT / p0, -n) * pT;
}

// // Next is dN/pTdpT and so on, changed to pT
// // x[0]  -> sqrt(m0*m0+x[0]*x[0]) - m0
/////////////////////////////////////////////////////////////////////////////
// exponential function for dN/ptdpt  // mTexponential
/////////////////////////////////////////////////////////////////////////////
Double_t exponential2(Double_t *x, Double_t *par) // x is pT
{
  const Double_t m0 = ximass;
  return par[0] * exp(-(sqrt(x[0] * x[0] + m0 * m0) - m0) / par[1]);
}

/////////////////////////////////////////////////////////////////////////////
// exponential*mT function for dN/dmT  // mTexponential
/////////////////////////////////////////////////////////////////////////////
Double_t exponentialmT2(Double_t *x, Double_t *par)// x is pT // y is *pT
{
  const Double_t m0 = ximass;
  return par[0] * exp(-(sqrt(x[0] * x[0] + m0 * m0) - m0) / par[1]) * x[0];
}

// // Next is dN/pTdpT and so on, changed to pT
// // x[0]  -> sqrt(m0*m0+x[0]*x[0]) - m0
/////////////////////////////////////////////////////////////////////////////
// mT Boltzmann function for dN/ptdpt  // mT Boltzmann
/////////////////////////////////////////////////////////////////////////////
Double_t BoltzmannmT(Double_t *x, Double_t *par)// x is pT // y is * pT
{
  const Double_t m0 = ximass;
  return par[0] * exp(-(sqrt(x[0] * x[0] + m0 * m0) - m0) / par[1]) * (sqrt(x[0] * x[0] + m0 * m0) - m0) * x[0];
}

/////////////////////////////////////////////////////////////////////////////
// exponential*pT function for dN/dpT  // pTexponential
/////////////////////////////////////////////////////////////////////////////
Double_t exponentialpT(Double_t *x, Double_t *par)// x is pT // y is *pT
{
  const Double_t m0 = ximass;
  // return par[0] * exp(-x[0] / par[1]) * x[0];
  return par[0] * exp(-pow(x[0], 1.5) / pow(par[1], 1.5)) * x[0];
}

/////////////////////////////////////////////////////////////////////////////
// exponential*pTgaus function for dN/dpT  // pT gauss exponential
/////////////////////////////////////////////////////////////////////////////
Double_t exponentialpTGaus(Double_t *x, Double_t *par)// x is pT // y is *pT
{
  const Double_t m0 = ximass;
  return par[0] * exp(-pow(x[0], 2) / pow(par[1], 2)) * x[0];
}

/////////////////////////////////////////////////////////////////////////////
// exponential*pTtrip function for dN/dpT  // pT trip exponential
/////////////////////////////////////////////////////////////////////////////
Double_t exponentialpTTrip(Double_t *x, Double_t *par)// x is pT // y is *pT
{
  const Double_t m0 = ximass;
  return par[0] * exp(-pow(x[0], 2.5) / pow(par[1], 2.5)) * x[0];
}

/////////////////////////////////////////////////////////////////////////////
// TBW function for dN/ptdpt
/////////////////////////////////////////////////////////////////////////////
double IntegrandBW(const double *x, const double *p)
{
  
  double radius = 1.;
  double Y = 1;
  double pi = TMath::Pi();
  
  double r   = x[0];
  double phi = x[1];
  double y = x[2];
  double mass = p[5];
  double pT   = p[0];
  double beta = p[1];
  double T    = p[2];
  double mT   = sqrt(mass * mass + pT * pT);
  double n    = p[3];
  double q    = p[4];
  double rho = ATanH(beta * pow(r / radius, n));
  double a1 = mT * CosH(y) * CosH(rho);
  double a2 = pT * SinH(rho) * Cos(phi);
  double b1 = 1 + (q - 1) / T * (a1 - a2);
  double b2 = 1 / (1 - q);
  
  return mT * CosH(y) * r * pow(b1, b2);
}

double BWFunc(const double *x, const double *p)
{
  
  double radius = 1.;
  double Y = 1;
  double pi = TMath::Pi();
  
  double pT   = x[0];
  double mass = p[5];
  double beta = p[0];
  double T    = p[1];
  double n    = p[2];
  double q    = p[3];
  double beta_s = (n + 2) / (n + 1) * beta;
  
  TF3 *fIntBG = 0;
  if (!fIntBG) fIntBG = new TF3("fIntBG", IntegrandBW, 0., radius, -TMath::Pi(), TMath::Pi(), -Y, Y, 6);
  
  fIntBG->SetParameters(pT, beta_s, T, n, q, mass);
  fIntBG->FixParameter(0, pT);
  fIntBG->FixParameter(5, mass);
  
  return fIntBG->Integral(0., radius, -TMath::Pi(), TMath::Pi(), -Y, Y, 1E-3) * p[4];
}

/////////////////////////////////////////////////////////////////////////////
// BW function for dN/ptdpt
/////////////////////////////////////////////////////////////////////////////
Double_t IntegrandBGBW(const double *x, const double *p)
{
  double radius = 1.;
  double Y = 1;
  double pi = TMath::Pi();
  double x0   = x[0];
  double mass = p[4];
  double pT   = p[0];
  double beta = p[1];
  double T    = p[2];
  double mT   = Sqrt(mass * mass + pT * pT);
  double n    = p[3];
  double beta_s = (n + 2) / 2 * beta;
  double rho0 = ATanH(beta_s * pow(x0 / radius, n));
  double a0   = pT * SinH(rho0) / T;
  double a1   = mT * CosH(rho0) / T;
  
  return x0 * mT * BesselI0(a0) * BesselK1(a1);
}

Double_t BGBW(const double *x, const double *p)
{
  double radius = 1.;
  double Y = 1;
  double pi = TMath::Pi();
  double pT   = x[0];
  double mass = p[4];
  double beta = p[0];
  double T    = p[1];
  double n    = p[2];
  
  TF1 *fIntBG = 0;
  if (!fIntBG) fIntBG = new TF1("fIntBG", IntegrandBGBW, 0., radius, 5);
  
  fIntBG->SetParameters(pT, beta, T, n, mass);
  fIntBG->FixParameter(0, pT);
  fIntBG->FixParameter(4, mass);
  
  double params[5] = {pT, beta, T, n, mass};
//  return fIntBG->Integral(0., radius, params) * p[3];
  return fIntBG->Integral(0., radius) * p[3];
}

/////////////////////////////////////////////////////////////////////////////
// BW function for dN/dpt
/////////////////////////////////////////////////////////////////////////////
//Taken from https://gitlab.cern.ch/mkrzewic/AliRoot/-/blob/fef6f2b32ddd823bd1c11b130f321a585288b032/PWG/Tools/AliPWGFunc.cxx
Double_t IntegrandBGBWPt(const double *x, const double *p)
{
  double radius = 1.;
  double Y = 1;
  double pi = TMath::Pi();
  double x0   = x[0];
  double mass = p[4];
  double pT   = p[0];
  double beta = p[1];
  double T    = p[2];
  double mT   = Sqrt(mass * mass + pT * pT);
  double n    = p[3];
  double beta_s = (n + 2) / 2 * beta;
  
  // Keep beta within reasonable limits
  Double_t beta_t = beta_s * pow(x0 / radius, n);
  if (beta_t > 0.9999999999999999) beta_t = 0.9999999999999999;
  
  double rho0 = ATanH(beta_s * pow(x0 / radius, n));
  double a0   = pT * SinH(rho0) / T;
  if (a0 > 700.) a0 = 700.; // avoid FPE
  
  double a1   = mT * CosH(rho0) / T;
  
  return x0 * mT * BesselI0(a0) * BesselK1(a1) * p[0];
}

Double_t BGBWPt(const double *x, const double *p)
{
  double radius = 1.;
  double Y = 1;
  double pi = TMath::Pi();
  double pT   = x[0];
  double mass = p[4];
  double beta = p[0];
  double T    = p[1];
  double n    = p[2];
  
  TF1 *fIntBG = 0;
  // if (!fIntBG) fIntBG = new TF1("fIntBG", IntegrandBGBWPt, 0., radius, 5);
  if (!fIntBG) fIntBG = new TF1("fIntBG", IntegrandBGBW, 0., radius, 5);
  
  fIntBG->SetParameters(pT, beta, T, n, mass);
  fIntBG->FixParameter(0, pT);
  fIntBG->FixParameter(4, mass);
  
  double params[5] = {pT, beta, T, n, mass};
//  return fIntBG->Integral(0., radius, params) * p[3] * x[0];
  return fIntBG->Integral(0., radius) * p[3] * x[0];
}

/////////////////////////////////////////////////////////////////////////////
// BW function for dN/dpt*pT
/////////////////////////////////////////////////////////////////////////////
Double_t IntegrandBGBWPt2(const double *x, const double *p)
{
  double radius = 1.;
  double Y = 1;
  double pi = TMath::Pi();
  double x0   = x[0];
  double mass = p[4];
  double pT   = p[0];
  double beta = p[1];
  double T    = p[2];
  double mT   = Sqrt(mass * mass + pT * pT);
  double n    = p[3];
  double beta_s = (n + 2) / 2 * beta;
  double rho0 = ATanH(beta_s * pow(x0 / radius, n));
  double a0   = pT * SinH(rho0) / T;
  double a1   = mT * CosH(rho0) / T;
  
  return x0 * mT * BesselI0(a0) * BesselK1(a1) * p[0] * p[0];
}

Double_t BGBWPt2(const double *x, const double *p)
{
  double radius = 1.;
  double Y = 1;
  double pi = TMath::Pi();
  double pT   = x[0];
  double mass = p[4];
  double beta = p[0];
  double T    = p[1];
  double n    = p[2];
  
  TF1 *fIntBG = 0;
  // if (!fIntBG) fIntBG = new TF1("fIntBG", IntegrandBGBWPt2, 0., radius, 5);
  if (!fIntBG) fIntBG = new TF1("fIntBG", IntegrandBGBW, 0., radius, 5);
  
  fIntBG->SetParameters(pT, beta, T, n, mass);
  fIntBG->FixParameter(0, pT);
  fIntBG->FixParameter(4, mass);
  
  double params[5] = {pT, beta, T, n, mass};
//  return fIntBG->Integral(0., radius, params) * p[3] * x[0] * x[0];
  return fIntBG->Integral(0., radius) * p[3] * x[0] * x[0];
}
