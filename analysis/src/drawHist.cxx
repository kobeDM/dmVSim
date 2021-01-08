// STL
#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>

// ROOT
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TPad.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>
#include <TLegend.h>
#include <TF1.h>
#include <TPaletteAxis.h>
#include <TLine.h>
#include <TMath.h>

using namespace std;

string IntToString(int number)
{
  stringstream ss;
  ss << number;
  return ss.str();
}

// drawHist
int main(int argc, char *argv[]){

  gStyle->SetPalette(1);
  gStyle->SetStatY(0.96);

  if(argc != 6) {
    cout << "<ERROR>    Specify input file name as first argument." << endl;
    return 1;
  }
  string filename = argv[1];
  string::size_type index = filename.find(".dat");
  if( index == string::npos ) {
    cout << "<ERROR>    Inputfile is not found." << std::endl;
    return 1;
  }
  string infilename = "RESULT_" + filename.substr(0, index) + "_Nucleus_" + argv[2] + "_EneThr_"  + argv[3] + "keV_HeadTail_" + argv[4]
                                + "_" + argv[5] + "_Samples" + "_anal.root";
  TFile* file = new TFile(infilename.c_str());

  TH1D* hr             = (TH1D*)file->Get("r");
  TH1D* hNucleus       = (TH1D*)file->Get("Nucleus");
  TH1D* hCosTheta      = (TH1D*)file->Get("CosTheta");
  TH1D* hSinPhi        = (TH1D*)file->Get("SinPhi");
  TH1D* hEr            = (TH1D*)file->Get("Er");
  TH1D* hIsotropicFlag = (TH1D*)file->Get("IsotropicFlag");
  TH1D* hNDF_CosTheta  = (TH1D*)file->Get("hNDF_CosTheta");
  TH1D* hNDF           = (TH1D*)file->Get("hNDF");

  TH2D** hEr_vs_CosTheta_r;
  TH2D** hEr_vs_CosTheta_scaled_r;
  TH2D** Clone_hEr_vs_CosTheta_r;
  TH2D** hEr_vs_CosTheta_pseudo_r;
  TH1D** hCosTheta_r;
  TH1D** hCosTheta_scaled_r;
  TH1D** hCosTheta_pseudo_r;
  TGraph2DErrors** hg2DError_Er_vs_CosTheta_pseudo_r;
  TGraph** g_Chi2_CosTheta_r;
  TGraph** g_Delta_Chi2_CosTheta_r;
  TGraph** g_KS_Test_pValue_CosTheta_r;
  TGraph** g_Chi2_r;
  TGraph** g_Delta_Chi2_r;
  TGraph** g_KS_Test_pValue_r;

  hEr_vs_CosTheta_r = new TH2D*[11];
  hEr_vs_CosTheta_scaled_r = new TH2D*[11];
  Clone_hEr_vs_CosTheta_r = new TH2D*[11];
  hEr_vs_CosTheta_pseudo_r = new TH2D*[11];
  hCosTheta_r = new TH1D*[11];
  hCosTheta_scaled_r = new TH1D*[11];
  hCosTheta_pseudo_r = new TH1D*[11];
  hg2DError_Er_vs_CosTheta_pseudo_r = new TGraph2DErrors*[11];
  g_Chi2_CosTheta_r = new TGraph*[11];
  g_Delta_Chi2_CosTheta_r = new TGraph*[11];
  g_KS_Test_pValue_CosTheta_r = new TGraph*[11];
  g_Chi2_r = new TGraph*[11];
  g_Delta_Chi2_r = new TGraph*[11];
  g_KS_Test_pValue_r = new TGraph*[11];

  for(int iRatio=0;iRatio<11;iRatio++){
    string histname = "Er_vs_CosTheta_r" + IntToString(iRatio);
    hEr_vs_CosTheta_r[iRatio] = (TH2D*)file->Get(histname.c_str());
    histname = "Er_vs_CosTheta_scaled_r" + IntToString(iRatio);
    hEr_vs_CosTheta_scaled_r[iRatio] = (TH2D*)file->Get(histname.c_str());
    histname = "Clone_Er_vs_CosTheta_r" + IntToString(iRatio);
    Clone_hEr_vs_CosTheta_r[iRatio] = (TH2D*)hEr_vs_CosTheta_r[iRatio]->Clone(histname.c_str());
    histname = "Er_vs_CosTheta_pseudo_r" + IntToString(iRatio);
    hEr_vs_CosTheta_pseudo_r[iRatio] = (TH2D*)file->Get(histname.c_str());
    histname = "CosTheta_r" + IntToString(iRatio);
    hCosTheta_r[iRatio] = (TH1D*)file->Get(histname.c_str());
    histname = "CosTheta_scaled_r" + IntToString(iRatio);
    hCosTheta_scaled_r[iRatio] = (TH1D*)file->Get(histname.c_str());
    histname = "CosTheta_pseudo_r" + IntToString(iRatio);
    hCosTheta_pseudo_r[iRatio] = (TH1D*)file->Get(histname.c_str());
    string graphname = "g2DError_Er_vs_CosTheta_pseudo_r" + IntToString(iRatio);
    hg2DError_Er_vs_CosTheta_pseudo_r[iRatio] = (TGraph2DErrors*)file->Get(graphname.c_str());
    graphname = "g_Chi2_r" + IntToString(iRatio);
    g_Chi2_r[iRatio] = (TGraph*)file->Get(graphname.c_str());
    graphname = "g_Delta_Chi2_r" + IntToString(iRatio);
    g_Delta_Chi2_r[iRatio] = (TGraph*)file->Get(graphname.c_str());
    graphname = "g_KS_Test_pValue_r" + IntToString(iRatio);
    g_KS_Test_pValue_r[iRatio] = (TGraph*)file->Get(graphname.c_str());
  }

  cout << "Drawing Histogram and Graph in Canvas..." << endl;
  gStyle->SetOptStat("e");
  TCanvas** canvas;
  canvas = new TCanvas*[21];

  canvas[0] = new TCanvas("canvas_0", "canvas_0", 0, 0, 1200, 800);
  canvas[0]->SetBorderMode(0);
  canvas[0]->Divide(3,3);
  canvas[0]->cd(1);
  gPad->SetBorderMode(0);
  gPad->SetFillColor(0);
  gPad->GetFrame()->SetBorderMode(0);
  gPad->GetFrame()->SetLineColor(1);
  gPad->SetFrameBorderMode(0);
  hr->SetXTitle("r");
  hr->SetYTitle("[Entries / bin]");
  hr->GetYaxis()->SetTitleOffset(1.4);
  hr->SetFillColor(38);
  hr->Draw();
  hr->Draw("sametext0");
  canvas[0]->cd(2);
  gPad->SetBorderMode(0);
  gPad->SetFillColor(0);
  gPad->GetFrame()->SetBorderMode(0);
  gPad->SetFrameBorderMode(0);
  hNucleus->SetXTitle("Nucleus");
  hNucleus->SetYTitle("[Entries / bin]");
  hNucleus->GetYaxis()->SetTitleOffset(1.4);
  hNucleus->SetFillColor(1);
  hNucleus->Draw();
  canvas[0]->cd(3);
  canvas[0]->cd(3)->SetLogy();
  gPad->SetBorderMode(0);
  gPad->SetFillColor(0);
  gPad->GetFrame()->SetBorderMode(0);
  gPad->SetFrameBorderMode(0);
  hEr->SetXTitle("Er [keV]");
  hEr->SetYTitle("[Entries / bin]");
  hEr->GetYaxis()->SetTitleOffset(1.4);
  hEr->SetFillColor(5);
  hEr->Draw();
  gPad->SetFillColor(0);
  canvas[0]->cd(4);
  gPad->SetBorderMode(0);
  gPad->SetFillColor(0);
  gPad->GetFrame()->SetBorderMode(0);
  gPad->SetFrameBorderMode(0);
  hCosTheta->SetXTitle("cos#theta");
  hCosTheta->SetYTitle("[Entries / bin]");
  hCosTheta->GetYaxis()->SetTitleOffset(1.4);
  hCosTheta->SetFillColor(2);
  hCosTheta->Draw();
  //TF1* f1 = new TF1("f1", "pol3",-0.9,0.9);
  //f1->SetLineColor(5);
  //hCosTheta->Fit("f1");
  //f1->Draw("same");
  //cout << f1->GetChisquare() << f1->GetNDF() << endl;
  canvas[0]->cd(5);
  gPad->SetBorderMode(0);
  gPad->SetFillColor(0);
  gPad->GetFrame()->SetBorderMode(0);
  gPad->SetFrameBorderMode(0);
  hSinPhi->SetXTitle("sin#phi");
  hSinPhi->SetYTitle("[Entries / bin]");
  hSinPhi->GetYaxis()->SetTitleOffset(1.4);
  hSinPhi->SetFillColor(4);
  hSinPhi->Draw();
  canvas[0]->cd(6);
  gPad->SetBorderMode(0);
  gPad->SetFillColor(0);
  gPad->GetFrame()->SetBorderMode(0);
  gPad->SetFrameBorderMode(0);
  hEr->SetXTitle("Er [keV]");
  hEr->SetYTitle("[Entries / bin]");
  hEr->GetYaxis()->SetTitleOffset(1.4);
  hEr->SetFillColor(5);
  hEr->Draw();
  canvas[0]->cd(7);
  gPad->SetBorderMode(0);
  gPad->SetFillColor(0);
  gPad->GetFrame()->SetBorderMode(0);
  gPad->SetFrameBorderMode(0);
  hIsotropicFlag->SetXTitle("IsotropicFlag");
  hIsotropicFlag->SetYTitle("[Entries / bin]");
  hIsotropicFlag->GetYaxis()->SetTitleOffset(1.4);
  hIsotropicFlag->GetXaxis()->SetNdivisions(10);
  hIsotropicFlag->SetFillColor(6);
  hIsotropicFlag->Draw();
  hIsotropicFlag->Draw("same text0");

  canvas[1] = new TCanvas("canvas_1", "canvas_1", 0, 0, 1200, 800);
  canvas[1]->SetBorderMode(0);
  canvas[1]->Divide(4,3);
  for(int iRatio=0;iRatio<11;iRatio++){
    canvas[1]->cd(iRatio+1);
    gPad->SetBorderMode(0);
    gPad->SetFillColor(0);
    gPad->GetFrame()->SetBorderMode(0);
    gPad->SetFrameBorderMode(0);
    hCosTheta_r[iRatio]->SetXTitle("cos#theta");
    hCosTheta_r[iRatio]->SetYTitle("[Entries / bin]");
    hCosTheta_r[iRatio]->GetYaxis()->SetTitleOffset(1.4);
    hCosTheta_r[iRatio]->SetFillColor(3);
    hCosTheta_r[iRatio]->Draw();
  }

  canvas[2] = new TCanvas("canvas_2", "canvas_2", 0, 0, 1200, 800);
  canvas[2]->SetBorderMode(0);
  canvas[2]->Divide(4,3);
  for(int iRatio=0;iRatio<11;iRatio++){
    canvas[2]->cd(iRatio+1);
    gPad->SetBorderMode(0);
    gPad->SetFillColor(0);
    gPad->GetFrame()->SetBorderMode(0);
    gPad->SetFrameBorderMode(0);
    hCosTheta_scaled_r[iRatio]->SetXTitle("cos#theta");
    hCosTheta_scaled_r[iRatio]->SetYTitle("[Entries / bin]");
    hCosTheta_scaled_r[iRatio]->GetYaxis()->SetTitleOffset(1.4);
    hCosTheta_scaled_r[iRatio]->SetFillColor(8);
    hCosTheta_scaled_r[iRatio]->Draw();
  }

  canvas[3] = new TCanvas("canvas_3", "canvas_3", 0, 0, 1200, 800);
  canvas[3]->SetBorderMode(0);
  canvas[3]->Divide(4,3);
  for(int iRatio=0;iRatio<11;iRatio++){
    canvas[3]->cd(iRatio+1);
    gPad->SetBorderMode(0);
    gPad->SetFillColor(0);
    gPad->GetFrame()->SetBorderMode(0);
    gPad->SetFrameBorderMode(0);
    hCosTheta_pseudo_r[iRatio]->SetXTitle("cos#theta");
    hCosTheta_pseudo_r[iRatio]->SetYTitle("[Entries / bin]");
    hCosTheta_pseudo_r[iRatio]->GetYaxis()->SetTitleOffset(1.4);
    hCosTheta_pseudo_r[iRatio]->SetFillColor(29);
    hCosTheta_pseudo_r[iRatio]->Draw();
  }

  canvas[4] = new TCanvas("canvas_4", "canvas_4", 0, 0, 1200, 800);
  canvas[4]->SetBorderMode(0);
  canvas[4]->Divide(4,3);
  for(int iRatio=0;iRatio<11;iRatio++){
    canvas[4]->cd(iRatio+1);
    gPad->SetBorderMode(0);
    gPad->SetFillColor(0);
    gPad->GetFrame()->SetBorderMode(0);
    gPad->SetFrameBorderMode(0);
    hEr_vs_CosTheta_r[iRatio]->SetXTitle("cos#theta");
    hEr_vs_CosTheta_r[iRatio]->SetYTitle("Er [keV]");
    hEr_vs_CosTheta_r[iRatio]->GetYaxis()->SetTitleOffset(1.4);
    hEr_vs_CosTheta_r[iRatio]->Draw("colz");
    //canvas[4]->Update();
    //TPaletteAxis *palette = (TPaletteAxis*)hEr_vs_CosTheta_r[iRatio]->GetListOfFunctions()->FindObject("palette");
    //palette->SetX1NDC(0.9);
    //palette->SetX2NDC(0.92);
    //palette->SetY1NDC(0.1);
    //palette->SetY2NDC(0.9);
    //palette->GetAxis()->SetLabelSize(0.02);
    //canvas[4]->Modified();
    //canvas[4]->Update();
  }

  for(int iRatio=0;iRatio<11;iRatio++){
    string canvasname = "canvas_" + IntToString(iRatio+5);
    canvas[iRatio+5] = new TCanvas(canvasname.c_str(), canvasname.c_str(), 0, 0, 1200, 800);
    canvas[iRatio+5]->SetBorderMode(0);
    canvas[iRatio+5]->Divide(1,1);
    canvas[iRatio+5]->cd(1);
    //canvas[iRatio+5]->cd(1)->SetPhi(180);
    //canvas[iRatio+5]->cd(1)->SetTheta(0);
    canvas[iRatio+5]->cd(1)->SetPhi(200);
    canvas[iRatio+5]->cd(1)->SetTheta(30);
    gPad->SetBorderMode(0);
    gPad->SetFillColor(0);
    gPad->GetFrame()->SetBorderMode(0);
    gPad->SetFrameBorderMode(0);
    hEr_vs_CosTheta_scaled_r[iRatio]->GetXaxis()->SetTitle("cos#theta");
    hEr_vs_CosTheta_scaled_r[iRatio]->GetYaxis()->SetTitle("Er [keV]");
    hEr_vs_CosTheta_scaled_r[iRatio]->GetZaxis()->SetTitle("[Entries / bin]");
    hEr_vs_CosTheta_scaled_r[iRatio]->GetXaxis()->SetLabelOffset(0.01);
    hEr_vs_CosTheta_scaled_r[iRatio]->GetXaxis()->SetTitleOffset(1.6);
    hEr_vs_CosTheta_scaled_r[iRatio]->GetYaxis()->SetTitleOffset(1.8);
    hEr_vs_CosTheta_scaled_r[iRatio]->GetXaxis()->SetLabelSize(0.045);
    hEr_vs_CosTheta_scaled_r[iRatio]->Draw("LEGO2");
    hg2DError_Er_vs_CosTheta_pseudo_r[iRatio]->SetFillColor(29);
    hg2DError_Er_vs_CosTheta_pseudo_r[iRatio]->SetMarkerSize(0.9);
    hg2DError_Er_vs_CosTheta_pseudo_r[iRatio]->SetMarkerStyle(20);
    hg2DError_Er_vs_CosTheta_pseudo_r[iRatio]->SetMarkerColor(kRed);
    hg2DError_Er_vs_CosTheta_pseudo_r[iRatio]->SetLineColor(1);
    hg2DError_Er_vs_CosTheta_pseudo_r[iRatio]->SetLineColor(1);
    hg2DError_Er_vs_CosTheta_pseudo_r[iRatio]->SetLineWidth(2);
    hg2DError_Er_vs_CosTheta_pseudo_r[iRatio]->Draw("same err p0");
    TLegend *leg = new TLegend(0.78,0.79,0.98,0.89);
    leg->AddEntry(hEr_vs_CosTheta_scaled_r[iRatio],"sim_template_data","f");
    leg->AddEntry(hg2DError_Er_vs_CosTheta_pseudo_r[iRatio],"sim_pseudo_data","ep");
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->Draw("same");
  }

  canvas[16] = new TCanvas("canvas_16", "canvas_16", 0, 0, 1200, 800);
  canvas[16]->SetBorderMode(0);
  canvas[16]->Divide(4,3);
  int ndf = hNDF->GetBinLowEdge(2);
  double reduced_chi2_68 = TMath::ChisquareQuantile(0.68, ndf)/ndf;
  double reduced_chi2_90 = TMath::ChisquareQuantile(0.90, ndf)/ndf;
  double reduced_chi2_95 = TMath::ChisquareQuantile(0.95, ndf)/ndf;
  TLine *line_chi2_68 = new TLine(0., reduced_chi2_68, 110., reduced_chi2_68);
  line_chi2_68->SetLineColor(80);
  line_chi2_68->SetLineWidth(1);
  line_chi2_68->SetLineStyle(2);
  TLine *line_chi2_90 = new TLine(0., reduced_chi2_90, 110., reduced_chi2_90);
  line_chi2_90->SetLineColor(7);
  line_chi2_90->SetLineWidth(1);
  line_chi2_90->SetLineStyle(2);
  TLine *line_chi2_95 = new TLine(0., reduced_chi2_95, 110., reduced_chi2_95);
  line_chi2_95->SetLineColor(6);
  line_chi2_95->SetLineWidth(1);
  line_chi2_95->SetLineStyle(2);
  for(int iRatio=0;iRatio<11;iRatio++){
    canvas[16]->cd(iRatio+1)->SetGridy();
    canvas[16]->cd(iRatio+1)->SetTickx();
    canvas[16]->cd(iRatio+1)->SetTicky();
    gPad->SetBorderMode(0);
    gPad->SetFillColor(0);
    gPad->GetFrame()->SetBorderMode(0);
    gPad->SetFrameBorderMode(0);
    g_Chi2_r[iRatio]->GetXaxis()->SetTitle("r");
    g_Chi2_r[iRatio]->GetYaxis()->SetTitle("chi2");
    g_Chi2_r[iRatio]->GetXaxis()->SetTitleSize(0.05);
    g_Chi2_r[iRatio]->GetYaxis()->SetTitleSize(0.05);
    g_Chi2_r[iRatio]->GetXaxis()->SetLabelSize(0.045);
    g_Chi2_r[iRatio]->SetMarkerColor(4);
    g_Chi2_r[iRatio]->Draw("AC*");
    line_chi2_68->Draw("same");
    line_chi2_90->Draw("same");
    line_chi2_95->Draw("same");
  }
  TLegend *legend_chi2 = new TLegend(0.1,0.7,0.3,0.9);
  legend_chi2->AddEntry(line_chi2_68,"68%","L");
  legend_chi2->AddEntry(line_chi2_90,"90%","L");
  legend_chi2->AddEntry(line_chi2_95,"95%","L");
  canvas[16]->cd(12);
  legend_chi2->Draw();

  canvas[17] = new TCanvas("canvas_17", "canvas_17", 0, 0, 1200, 800);
  canvas[17]->SetBorderMode(0);
  canvas[17]->Divide(4,3);
  for(int iRatio=0;iRatio<11;iRatio++){
    canvas[17]->cd(iRatio+1)->SetGridy();
    canvas[17]->cd(iRatio+1)->SetTickx();
    canvas[17]->cd(iRatio+1)->SetTicky();
    canvas[17]->cd(iRatio+1)->SetLogy();
    gPad->SetBorderMode(0);
    gPad->SetFillColor(0);
    gPad->GetFrame()->SetBorderMode(0);
    gPad->SetFrameBorderMode(0);
    g_Chi2_r[iRatio]->GetXaxis()->SetTitle("r");
    g_Chi2_r[iRatio]->GetYaxis()->SetTitle("chi2");
    g_Chi2_r[iRatio]->GetXaxis()->SetTitleSize(0.05);
    g_Chi2_r[iRatio]->GetYaxis()->SetTitleSize(0.05);
    g_Chi2_r[iRatio]->GetXaxis()->SetLabelSize(0.045);
    g_Chi2_r[iRatio]->SetMarkerColor(4);
    g_Chi2_r[iRatio]->Draw("AC*");
    line_chi2_68->Draw("same");
    line_chi2_90->Draw("same");
    line_chi2_95->Draw("same");
  }
  canvas[17]->cd(12);
  legend_chi2->Draw();

  canvas[18] = new TCanvas("canvas_18", "canvas_18", 0, 0, 1200, 800);
  canvas[18]->SetBorderMode(0);
  canvas[18]->Divide(4,3);
  TLine *line_delta_chi2_68 = new TLine(0., 1., 110., 1.);
  line_delta_chi2_68->SetLineColor(80);
  line_delta_chi2_68->SetLineWidth(1);
  line_delta_chi2_68->SetLineStyle(2);
  for(int iRatio=0;iRatio<11;iRatio++){
    canvas[18]->cd(iRatio+1)->SetGridy();
    canvas[18]->cd(iRatio+1)->SetTickx();
    canvas[18]->cd(iRatio+1)->SetTicky();
    gPad->SetBorderMode(0);
    gPad->SetFillColor(0);
    gPad->GetFrame()->SetBorderMode(0);
    gPad->SetFrameBorderMode(0);
    g_Delta_Chi2_r[iRatio]->GetXaxis()->SetTitle("r");
    g_Delta_Chi2_r[iRatio]->GetYaxis()->SetTitle("delta_chi2");
    g_Delta_Chi2_r[iRatio]->GetXaxis()->SetTitleSize(0.05);
    g_Delta_Chi2_r[iRatio]->GetYaxis()->SetTitleSize(0.05);
    g_Delta_Chi2_r[iRatio]->GetXaxis()->SetLabelSize(0.045);
    g_Delta_Chi2_r[iRatio]->SetMarkerColor(2);
    g_Delta_Chi2_r[iRatio]->Draw("A*");
    line_delta_chi2_68->Draw("same");
  }
  TLegend *legend_delta_chi2 = new TLegend(0.1,0.7,0.3,0.9);
  legend_delta_chi2->AddEntry(line_chi2_68,"68%","L");
  canvas[18]->cd(12);
  legend_delta_chi2->Draw();

  canvas[19] = new TCanvas("canvas_19", "canvas_19", 0, 0, 1200, 800);
  canvas[19]->SetBorderMode(0);
  canvas[19]->Divide(4,3);
  for(int iRatio=0;iRatio<11;iRatio++){
    canvas[19]->cd(iRatio+1)->SetGridy();
    canvas[19]->cd(iRatio+1)->SetTickx();
    canvas[19]->cd(iRatio+1)->SetTicky();
    canvas[19]->cd(iRatio+1)->SetLogy();
    gPad->SetBorderMode(0);
    gPad->SetFillColor(0);
    gPad->GetFrame()->SetBorderMode(0);
    gPad->SetFrameBorderMode(0);
    g_Delta_Chi2_r[iRatio]->GetXaxis()->SetTitle("r");
    g_Delta_Chi2_r[iRatio]->GetYaxis()->SetTitle("delta_chi2");
    g_Delta_Chi2_r[iRatio]->GetXaxis()->SetTitleSize(0.05);
    g_Delta_Chi2_r[iRatio]->GetYaxis()->SetTitleSize(0.05);
    g_Delta_Chi2_r[iRatio]->GetXaxis()->SetLabelSize(0.045);
    g_Delta_Chi2_r[iRatio]->SetMarkerColor(2);
    g_Delta_Chi2_r[iRatio]->Draw("A*");
    line_delta_chi2_68->Draw("same");
  }
  canvas[19]->cd(12);
  legend_delta_chi2->Draw();

  canvas[20] = new TCanvas("canvas_20", "canvas_20", 0, 0, 1200, 800);
  canvas[20]->SetBorderMode(0);
  canvas[20]->Divide(4,3);
  double kstest68CL = 1 - TMath::KolmogorovProb(0.96);
  TLine *line_68 = new TLine(0., kstest68CL, 110., kstest68CL);
  line_68->SetLineColor(80);
  line_68->SetLineWidth(1);
  line_68->SetLineStyle(2);
  double kstest90CL = 1 - TMath::KolmogorovProb(1.23);
  TLine *line_90 = new TLine(0., kstest90CL, 110., kstest90CL);
  line_90->SetLineColor(7);
  line_90->SetLineWidth(1);
  line_90->SetLineStyle(2);
  double kstest95CL = 1 - TMath::KolmogorovProb(1.36);
  TLine *line_95 = new TLine(0., kstest95CL, 110., kstest95CL);
  line_95->SetLineColor(6);
  line_95->SetLineWidth(1);
  line_95->SetLineStyle(2);
  for(int iRatio=0;iRatio<11;iRatio++){
    canvas[20]->cd(iRatio+1)->SetGridy();
    canvas[20]->cd(iRatio+1)->SetTickx();
    canvas[20]->cd(iRatio+1)->SetTicky();
    gStyle->SetNdivisions(211);
    gPad->SetBorderMode(0);
    gPad->SetFillColor(0);
    gPad->GetFrame()->SetBorderMode(0);
    gPad->SetFrameBorderMode(0);
    g_KS_Test_pValue_r[iRatio]->SetMaximum(1.1);
    g_KS_Test_pValue_r[iRatio]->GetXaxis()->SetTitle("r");
    g_KS_Test_pValue_r[iRatio]->GetYaxis()->SetTitle("ks test p-Value");
    g_KS_Test_pValue_r[iRatio]->GetXaxis()->SetTitleSize(0.05);
    g_KS_Test_pValue_r[iRatio]->GetYaxis()->SetTitleSize(0.05);
    g_KS_Test_pValue_r[iRatio]->GetXaxis()->SetLabelSize(0.045);
    g_KS_Test_pValue_r[iRatio]->SetMarkerColor(95);
    g_KS_Test_pValue_r[iRatio]->Draw("A*");
    line_68->Draw("same");
    line_90->Draw("same");
    line_95->Draw("same");
  }
  TLegend *legend_ks_test = new TLegend(0.1,0.7,0.3,0.9);
  legend_ks_test->AddEntry(line_68,"68%","L");
  legend_ks_test->AddEntry(line_90,"90%","L");
  legend_ks_test->AddEntry(line_95,"95%","L");
  canvas[18]->cd(12);
  legend_ks_test->Draw();

  string savefile1 = "RESULT_" + filename.substr(0, index) + "_Nucleus_" + argv[2] + "_EneThr_" + argv[3] + "keV_HeadTail_" + argv[4] + "_" + argv[5] + "_Samples"+ ".pdf(";
  string savefile2 = "RESULT_" + filename.substr(0, index) + "_Nucleus_" + argv[2] + "_EneThr_" + argv[3] + "keV_HeadTail_" + argv[4] + "_" + argv[5] + "_Samples"+ ".pdf";
  string savefile3 = "RESULT_" + filename.substr(0, index) + "_Nucleus_" + argv[2] + "_EneThr_" + argv[3] + "keV_HeadTail_" + argv[4] + "_" + argv[5] + "_Samples"+ ".pdf)";
  canvas[0]->Print(savefile1.c_str(),"pdf Portrait");
  canvas[1]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[2]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[3]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[4]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[5]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[6]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[7]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[8]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[9]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[10]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[11]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[12]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[13]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[14]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[15]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[16]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[17]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[18]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[19]->Print(savefile2.c_str(),"pdf Portrait");
  canvas[20]->Print(savefile3.c_str(),"pdf Portrait");

  return 0;
}
