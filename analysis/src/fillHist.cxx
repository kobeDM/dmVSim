// STL
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// ROOT
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TFile.h>

#include "param.h"

using namespace std;


int main(int argc, char *argv[]){
  if(argc != 2) {
    cout << "<ERROR>    Specify input file name as first argument." << endl;
    return 1;
  }
  string filename = argv[1];
  string::size_type index = filename.find(".dat");
  if( index == string::npos ) { 
    cout << "<ERROR>    Inputfile is not found." << endl;
    return 1;
  }
  string outfilename = "RESULT_" + filename.substr(0, index) + "_hist.root";
  TFile *f = new TFile(outfilename.c_str(),"recreate");

  Int_t    r = 0;
  Int_t    Nucleus = 0;
  Double_t CosTheta = 0;
  Double_t SinPhi = 0;
  Double_t Er = 0;
  Int_t    IsotropicFlag = 0;

  //-------- HISTGRAMS -------//
  TH2D* hEr_vs_CosTheta_r= new TH2D("hEr_vs_CosTheta","hEr_vs_CosTheta", nbin_CosTheta, bin_min_CosTheta, bin_max_CosTheta, nbin_Er, bin_min_Er, bin_max_Er);
  TH1D* hCosTheta      = new TH1D("CosTheta", "CosTheta", nbin_CosTheta, bin_min_CosTheta, bin_max_CosTheta);
  TH1D* hSinPhi        = new TH1D("SinPhi", "SinPhi", 240, -1.2, 1.2);
  TH1D* hEr            = new TH1D("Er", "Er", nbin_Er, bin_min_Er, bin_max_Er);
  TH1D* hIsotropicFlag = new TH1D("IsotropicFlag", "IsotropicFlag", 5, -0.75, 1.75);
  TH1D* hr       = new TH1D("r", "r", 60, 0, 120);
  TH1D* hNucleus = new TH1D("Nucleus", "Nucleus", 20, 0, 20);


  //--- open file -----
  ifstream inFile(filename.c_str());
  Long64_t n_evt = 0;
  while(inFile >> r >> Nucleus >> CosTheta >> SinPhi >> Er >> IsotropicFlag)
  {
    n_evt += 1;
    if(n_evt%10000==0) std::cout << "\rFill Data into Tree ... : Evevt # " << n_evt << " " << flush;

    hr->Fill(r);
    hNucleus->Fill(Nucleus);
    hCosTheta->Fill(CosTheta);
    hSinPhi->Fill(SinPhi);
    hEr->Fill(Er);
    hIsotropicFlag->Fill(IsotropicFlag);
    hEr_vs_CosTheta_r->Fill(CosTheta, Er);

  }
  cout << "\rFill Data into Tree ... : Evevt # " << n_evt << " " << flush;

  hr->Write();
  hNucleus->Write();
  hCosTheta->Write();
  hSinPhi->Write();
  hEr->Write();
  hIsotropicFlag->Write();
  hEr_vs_CosTheta_r->Write();


  f->Close();
  delete f;
  cout << endl;
  return 0;
}
