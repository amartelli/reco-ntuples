#include "TLegend.h"
#include "TLatex.h"
#include <TH2.h>
#include <TH1.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TTree.h"
#include "TChain.h"
#include <vector>
#include <fstream>
#include <string>
#include "TROOT.h"
#include "TSystem.h"



void compareAxis(){
  gROOT->Macro("/afs/cern.ch/user/a/amartell/public/setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);



  TFile* inFM[2][4];
  inFM[0][0] = TFile::Open("../test/clustering_AxisAnalysis_Ele11_Pt5_def.root");
  inFM[0][1] = TFile::Open("../test/clustering_AxisAnalysis_Ele11_Pt5_def.root");
  inFM[0][2] = TFile::Open("../test/clustering_AxisAnalysis_Ele11_Pt5_def.root");
  inFM[0][3] = TFile::Open("../test/clustering_AxisAnalysis_Ele11_Pt5_def.root");
  inFM[1][0] = TFile::Open("../test/clustering_AxisAnalysis_Ele11_Pt5_ax.root");
  inFM[1][1] = TFile::Open("../test/clustering_AxisAnalysis_Ele11_Pt5_ax.root");
  inFM[1][2] = TFile::Open("../test/clustering_AxisAnalysis_Ele11_Pt5_ax.root");
  inFM[1][3] = TFile::Open("../test/clustering_AxisAnalysis_Ele11_Pt5_ax.root");
  

  /*
  //  inFM[0][0] = TFile::Open("../test/clustering_AxisAnalysis_Pho22_E5_def_VtxSM.root");
  inFM[0][0] = TFile::Open("../test/clustering_AxisAnalysis_Ele11_Pt5_def.root");
  inFM[0][1] = TFile::Open("../test/clustering_AxisAnalysis_Pho22_E60_def_VtxSM.root");
  inFM[0][2] = TFile::Open("../test/clustering_AxisAnalysis_Pho211_E5_def_VtxSM.root");
  inFM[0][3] = TFile::Open("../test/clustering_AxisAnalysis_Pho211_E30_def_VtxSM.root");
  //  inFM[1][0] = TFile::Open("../test/clustering_AxisAnalysis_Pho22_E5_ax_VtxSM.root");
  inFM[1][0] = TFile::Open("../test/clustering_AxisAnalysis_Ele11_Pt5_ax.root");
  inFM[1][1] = TFile::Open("../test/clustering_AxisAnalysis_Pho22_E60_ax_VtxSM.root");
  inFM[1][2] = TFile::Open("../test/clustering_AxisAnalysis_Pho211_E5_ax_VtxSM.root");
  inFM[1][3] = TFile::Open("../test/clustering_AxisAnalysis_Pho211_E30_ax_VtxSM.root");
  */
  /*
  inFM[0][0] = TFile::Open("../test/clustering_AxisAnalysis_Pho22_E5_def_Vtx000.root");
  inFM[0][1] = TFile::Open("../test/clustering_AxisAnalysis_Pho22_E60_def_Vtx000.root");
  inFM[0][2] = TFile::Open("../test/clustering_AxisAnalysis_Pho211_E5_def_Vtx000.root");
  inFM[0][3] = TFile::Open("../test/clustering_AxisAnalysis_Pho211_E30_def_Vtx000.root");
  inFM[1][0] = TFile::Open("../test/clustering_AxisAnalysis_Pho22_E5_ax_Vtx000.root");
  inFM[1][1] = TFile::Open("../test/clustering_AxisAnalysis_Pho22_E60_ax_Vtx000.root");
  inFM[1][2] = TFile::Open("../test/clustering_AxisAnalysis_Pho211_E5_ax_Vtx000.root");
  inFM[1][3] = TFile::Open("../test/clustering_AxisAnalysis_Pho211_E30_ax_Vtx000.root");
  */  
  TH1F* h_2DClSumE[2][4];
  TProfile* tp_2DClSumE_vsRadius_wrtGen[2][4];
  TH1F* h_dR_2DClSel_Gen[2][4];
  TH1F* h_mtcRHSumE[4];
  TProfile* tp_mtcRHSumE_vsRadius_wrtGen[4];
  TH1F* h_dR_MatchHit_Gen[4];

  std::vector<std::string> type;
  type.push_back("ID22_E5");
  type.push_back("ID22_E60");
  type.push_back("ID211_E5");
  type.push_back("ID211_E30");
  type.push_back("ID22_E5_v1");
  type.push_back("ID22_E60_v1");
  type.push_back("ID211_E5_v1");
  type.push_back("ID211_E30_v1");

  int iColor[9] = {kRed, kOrange-3, kOrange-2, kBlue, kBlue-9, kCyan, kGreen+1, kCyan-2, kYellow+2}; //kGray+1};
  int iColors[2] = {kRed, kBlue};

  for(int iA=0; iA<2; ++iA){
    for(int iP=0; iP<4; ++iP){
      std::cout << "iA = " << iA << " iP = " << iP << std::endl;
      //h_2DClSumE[iA][iP] = (TH1F*)inFM[iA][iP]->Get("ana/h_2DClSelSumE");
      h_2DClSumE[iA][iP] = (TH1F*)inFM[iA][iP]->Get("ana/h_2DClSelConsSumE");
      //      h_2DClSumE[iA][iP] = (TH1F*)inFM[iA][iP]->Get("ana/h_2DClSumE");
      h_2DClSumE[iA][iP]->SetName(("h_2DClSumE_"+type.at(iA*4 + iP)).c_str());

      //tp_2DClSumE_vsRadius_wrtGen[iA][iP] = (TProfile*)inFM[iA][iP]->Get("ana/tp_2DClSumE_vsRadius_wrtGen");
      tp_2DClSumE_vsRadius_wrtGen[iA][iP] = (TProfile*)inFM[iA][iP]->Get("ana/tp_2DClSelConsSumE_vsRadius_wrtGen");
      //tp_2DClSumE_vsRadius_wrtGen[iA][iP] = (TProfile*)inFM[iA][iP]->Get("ana/tp_2DClSelSumE_vsRadius_wrtGen");
      tp_2DClSumE_vsRadius_wrtGen[iA][iP]->SetName(("tp_2DClSumE_vsRadius_wrtGen_"+type.at(iA*4 + iP)).c_str());

      h_dR_2DClSel_Gen[iA][iP] = (TH1F*)inFM[iA][iP]->Get("ana/h_dR_2DClSel_Gen");
      h_dR_2DClSel_Gen[iA][iP]->SetName(("h_dR_2DClSel_Gen_"+type.at(iA*4 + iP)).c_str());

      h_2DClSumE[iA][iP]->SetLineColor(iColors[iA]);
      h_2DClSumE[iA][iP]->SetLineWidth(2);
      h_2DClSumE[iA][iP]->SetMarkerColor(iColors[iA]);
      tp_2DClSumE_vsRadius_wrtGen[iA][iP]->SetLineColor(iColors[iA]);
      tp_2DClSumE_vsRadius_wrtGen[iA][iP]->SetMarkerColor(iColors[iA]);
      tp_2DClSumE_vsRadius_wrtGen[iA][iP]->SetMarkerStyle(iA==0 ? 20 : 4);


      h_dR_2DClSel_Gen[iA][iP]->SetLineColor(iColors[iA]);
      h_dR_2DClSel_Gen[iA][iP]->SetLineWidth(2);
      h_dR_2DClSel_Gen[iA][iP]->SetMarkerColor(iColors[iA]);


      h_2DClSumE[iA][iP]->Rebin(4);
      h_dR_2DClSel_Gen[iA][iP]->Rebin(4);
      h_dR_2DClSel_Gen[iA][iP]->Scale(1./h_dR_2DClSel_Gen[iA][iP]->Integral());

      if(iA == 0){
	h_mtcRHSumE[iP] = (TH1F*)inFM[iA][iP]->Get("ana/h_mtcRHSumE");
	h_mtcRHSumE[iP]->SetName(("h_mtcRHSumE_"+type.at(iP)).c_str());
	tp_mtcRHSumE_vsRadius_wrtGen[iP] = (TProfile*)inFM[iA][iP]->Get("ana/tp_mtcRHSumE_vsRadius_wrtGen");
	tp_mtcRHSumE_vsRadius_wrtGen[iP]->SetName(("tp_mtcRHSumE_vsRadius_wrtGen_"+type.at(iP)).c_str());
	h_dR_MatchHit_Gen[iP] = (TH1F*)inFM[iA][iP]->Get("ana/h_dR_MatchHit_Gen");
	h_dR_MatchHit_Gen[iP]->SetName(("h_dR_MatchHit_Gen_"+type.at(iP)).c_str());

	h_mtcRHSumE[iP]->SetLineWidth(2);
	h_mtcRHSumE[iP]->SetLineColor(kBlack);
	h_mtcRHSumE[iP]->SetMarkerColor(kBlack);

	tp_mtcRHSumE_vsRadius_wrtGen[iP]->SetLineColor(kBlack);
	tp_mtcRHSumE_vsRadius_wrtGen[iP]->SetMarkerColor(kBlack);
	tp_mtcRHSumE_vsRadius_wrtGen[iP]->SetMarkerStyle(20);    

	h_dR_MatchHit_Gen[iP]->SetLineWidth(2);
	h_dR_MatchHit_Gen[iP]->SetLineColor(kBlack);
	h_dR_MatchHit_Gen[iP]->SetMarkerColor(kBlack);

	h_mtcRHSumE[iP]->Rebin(4);
	h_dR_MatchHit_Gen[iP]->Rebin(4);
	h_dR_MatchHit_Gen[iP]->Scale(1./h_dR_MatchHit_Gen[iP]->Integral());
      }
    }
  }

  std::cout << " ci sono " << std::endl;

  TLegend *legTGM = new TLegend(0.70,0.70,0.90,0.90,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetTextSize(0.02);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  legTGM->AddEntry(h_2DClSumE[0][0], "def.", "pl");
  legTGM->AddEntry(h_2DClSumE[1][0], "baryc_v1", "pl");
  legTGM->AddEntry(h_mtcRHSumE[0], "gen", "pl");

  std::cout << "leg ok  " << std::endl;

  //  std::string folder = "plots"+std::string(Form("_E%d",ene));
  std::string folder = "plots_Axis";

  TCanvas* ch_1[4];
  for(int iP=0; iP<4; ++iP){
    ch_1[iP] = new TCanvas();
    ch_1[iP]->cd();
    if(iP == 0) h_2DClSumE[0][iP]->GetYaxis()->SetRangeUser(0., 1000.);
    if(iP == 1) h_2DClSumE[0][iP]->GetYaxis()->SetRangeUser(0., 2000.);
    if(iP == 2) h_2DClSumE[0][iP]->GetYaxis()->SetRangeUser(0., 2000.);
    if(iP == 3) h_2DClSumE[0][iP]->GetYaxis()->SetRangeUser(0., 1000.);
    h_2DClSumE[0][iP]->GetXaxis()->SetTitle("reco/gen energy");
    h_2DClSumE[0][iP]->Draw();
    h_2DClSumE[1][iP]->Draw("same");
    h_mtcRHSumE[iP]->Draw("same");
    legTGM->SetHeader(type.at(iP).c_str());
    legTGM->Draw("same");
    ch_1[iP]->Print((folder+"/h_2DClSumE_"+type.at(iP)+".png").c_str(), "png");
  }


  TCanvas* ch_2[4];
  for(int iP=0; iP<4; ++iP){
    ch_2[iP] = new TCanvas();
    ch_2[iP]->cd();
    if(iP == 0 || iP == 1) tp_2DClSumE_vsRadius_wrtGen[0][iP]->GetYaxis()->SetRangeUser(0., 1.);
    if(iP == 2 || iP == 3) tp_2DClSumE_vsRadius_wrtGen[0][iP]->GetYaxis()->SetRangeUser(0., 0.2);
    tp_2DClSumE_vsRadius_wrtGen[0][iP]->GetXaxis()->SetTitle("radius wrt gen");
    tp_2DClSumE_vsRadius_wrtGen[0][iP]->GetYaxis()->SetTitle("reco/gen energy");
    tp_2DClSumE_vsRadius_wrtGen[0][iP]->Draw();
    tp_2DClSumE_vsRadius_wrtGen[1][iP]->Draw("same");
    tp_mtcRHSumE_vsRadius_wrtGen[iP]->Draw("same");
    legTGM->SetHeader(type.at(iP).c_str());
    legTGM->Draw("same");
    ch_2[iP]->Print((folder+"/tp_2DClSumE_vsRadius_wrtGen_"+type.at(iP)+".png").c_str(), "png");
  }


  TCanvas* ch_3[4];
  for(int iP=0; iP<4; ++iP){
    ch_3[iP] = new TCanvas();
    ch_3[iP]->cd();
    h_dR_2DClSel_Gen[0][iP]->GetXaxis()->SetTitle("radius wrt gen");
    h_dR_2DClSel_Gen[0][iP]->Draw();
    h_dR_2DClSel_Gen[1][iP]->Draw("same");
    h_dR_MatchHit_Gen[iP]->Draw("same");
    legTGM->SetHeader(type.at(iP).c_str());
    legTGM->Draw("same");
    ch_3[iP]->Print((folder+"/h_dR_2DClSel_Gen_"+type.at(iP)+".png").c_str(), "png");
  }

}
