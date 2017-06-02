#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include "TRandom3.h"


void plotTimingResolution(){

  int iColors[3] = {kBlack, kRed, kBlue};
  int iColors05[3] = {kGray+1, kRed+2, kAzure-2};


  //from TB
  float paramA[3];
  float parErrA[3];
  float paramC[3];
  float parErrC[3];
  float noise[3];


  float paramA_MIP[3];
  float parErrA_MIP[3];
  float paramC_MIP[3];
  float parErrC_MIP[3];

  //100um
  paramA[0] = 1.00;
  parErrA[0] = 0.01;
  paramC[0] = 0.009;
  parErrC[0] = 0.001;

  paramA_MIP[0] = 0.69;
  parErrA_MIP[0] = 0.01;
  paramC_MIP[0] = 0.010;
  parErrC_MIP[0] = 0.001;

  //200um  
  paramA[1] = 1.06;
  parErrA[1] = 0.02;
  paramC[1] = 0.008;
  parErrC[1] = 0.001;

  paramA_MIP[1] = 0.38;
  parErrA_MIP[1] = 0.01;
  paramC_MIP[1] = 0.009;
  parErrC_MIP[1] = 0.001;

  //300um 
  paramA[2] = 1.11;
  parErrA[2] = 0.02;
  paramC[2] = 0.010;
  parErrC[2] = 0.001;

  paramA_MIP[2] = 0.34;
  parErrA_MIP[2] = 0.01;
  paramC_MIP[2] = 0.010;
  parErrC_MIP[2] = 0.001;

  //noise
  noise[0] = 0.269;
  noise[1] = 0.131;
  noise[2] = 0.066;
  //SoverN
  // noise[0] = 1.42;
  // noise[1] = 2.65;
  // noise[2] = 3.2;



  TF1* timeResolution_100 = new TF1("timeSi100", "sqrt(pow([0]/x/sqrt(2.), 2) + pow([1], 2) )", 1., 1000.);
  TF1* timeResolution_200 = new TF1("timeSi200", "sqrt(pow([0]/x/sqrt(2.), 2) + pow([1], 2) )", 1., 1000.);
  TF1* timeResolution_300 = new TF1("timeSi300", "sqrt(pow([0]/x/sqrt(2.), 2) + pow([1], 2) )", 1., 1000.);

  TF1* timeMIP_100 = new TF1("timeMIP_100", "sqrt(pow([0]/x*1.42/sqrt(2.), 2) + pow([1], 2) )", 1., 1000.);
  TF1* timeMIP_200 = new TF1("timeMIP_200", "sqrt(pow([0]/x*2.65/sqrt(2.), 2) + pow([1], 2) )", 1., 1000.);
  TF1* timeMIP_300 = new TF1("timeMIP_300", "sqrt(pow([0]/x*3.2/sqrt(2.), 2) + pow([1], 2) )", 1., 1000.);


  //new option
  TF1* timeResolution05_100 = new TF1("timeSi100_05", "sqrt(pow([0]/x/sqrt(2.), 2) + pow([1], 2) )", 1., 1000.);
  TF1* timeResolution05_200 = new TF1("timeSi200_05", "sqrt(pow([0]/x/sqrt(2.), 2) + pow([1], 2) )", 1., 1000.);
  TF1* timeResolution05_300 = new TF1("timeSi300_05", "sqrt(pow([0]/x/sqrt(2.), 2) + pow([1], 2) )", 1., 1000.);

  TF1* timeMIP05_100 = new TF1("timeMIP05_100", "sqrt(pow([0]/x*1.42/sqrt(2.), 2) + pow([1], 2) )", 1., 1000.);
  TF1* timeMIP05_200 = new TF1("timeMIP05_200", "sqrt(pow([0]/x*2.65/sqrt(2.), 2) + pow([1], 2) )", 1., 1000.);
  TF1* timeMIP05_300 = new TF1("timeMIP05_300", "sqrt(pow([0]/x*3.2/sqrt(2.), 2) + pow([1], 2) )", 1., 1000.);


  TCanvas* cTB = new TCanvas(); 
  gPad->SetLogx();
  gPad->SetLogy();

  timeResolution_100->SetParameters(paramA[0], paramC[0]);
  timeResolution_200->SetParameters(paramA[1], paramC[1]);
  timeResolution_300->SetParameters(paramA[2], paramC[2]);

  timeResolution_100->SetLineColor(iColors[0]);
  timeResolution_200->SetLineColor(iColors[1]);
  timeResolution_300->SetLineColor(iColors[2]);

  timeResolution_100->SetRange(1, 200);

  timeMIP_100->SetParameters(paramA_MIP[0], paramC_MIP[0]);
  timeMIP_200->SetParameters(paramA_MIP[1], paramC_MIP[1]);
  timeMIP_300->SetParameters(paramA_MIP[2], paramC_MIP[2]);

  timeMIP_100->SetLineColor(iColors[0]);
  timeMIP_200->SetLineColor(iColors[1]);
  timeMIP_300->SetLineColor(iColors[2]);
  timeMIP_100->SetLineStyle(2);
  timeMIP_200->SetLineStyle(2);
  timeMIP_300->SetLineStyle(2);

  //////////////////////////////////////////////////
  timeResolution05_100->SetParameters(paramA[0]/2., paramC[0]);
  timeResolution05_200->SetParameters(paramA[1]/2., paramC[1]);
  timeResolution05_300->SetParameters(paramA[2]/2., paramC[2]);

  timeResolution05_100->SetLineColor(iColors05[0]);
  timeResolution05_200->SetLineColor(iColors05[1]);
  timeResolution05_300->SetLineColor(iColors05[2]);

  timeResolution05_100->SetRange(1, 200);

  timeMIP05_100->SetParameters(paramA_MIP[0]/2., paramC_MIP[0]);
  timeMIP05_200->SetParameters(paramA_MIP[1]/2., paramC_MIP[1]);
  timeMIP05_300->SetParameters(paramA_MIP[2]/2., paramC_MIP[2]);

  timeMIP05_100->SetLineColor(iColors05[0]);
  timeMIP05_200->SetLineColor(iColors05[1]);
  timeMIP05_300->SetLineColor(iColors05[2]);
  timeMIP05_100->SetLineStyle(2);
  timeMIP05_200->SetLineStyle(2);
  timeMIP05_300->SetLineStyle(2);

  cTB->cd();
  timeResolution_100->Draw();
  timeResolution_200->Draw("same");
  timeResolution_300->Draw("same");
  timeMIP_100->Draw("same");
  timeMIP_200->Draw("same");
  timeMIP_300->Draw("same");

  //  cTB_05->cd();
  timeResolution05_100->Draw("same");
  timeResolution05_200->Draw("same");
  timeResolution05_300->Draw("same");
  timeMIP05_100->Draw("same");
  timeMIP05_200->Draw("same");
  timeMIP05_300->Draw("same");

  cTB->Print("timeResolution_TB.png", "png");

  // TCanvas* cTB_05 = new TCanvas(); 
  // gPad->SetLogx();
  // gPad->SetLogy();

  //  cTB_05->Print("timeResolution_TB_half.png", "png");



}
