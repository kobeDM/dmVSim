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

using namespace std;

// fillTree
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
  string outfilename = "RESULT_" + filename.substr(0, index) + "_tree.root";
  TFile *f = new TFile(outfilename.c_str(),"recreate");

  Int_t    r = 0;
  Int_t    Nucleus = 0;
  Double_t CosTheta = 0;
  Double_t SinPhi = 0;
  Double_t Er = 0;
  Int_t    IsotropicFlag = 0;

  TTree* tree= new TTree("DMVD", "DMVD");
  tree->Branch("r", &r, "r/I");
  tree->Branch("Nucleus", &Nucleus, "Nucleus/I");
  tree->Branch("CosTheta", &CosTheta, "CosTheta/D");
  tree->Branch("SinPhi", &SinPhi, "SinPhi/D");
  tree->Branch("Er", &Er, "Er/D");
  tree->Branch("IsotropicFlag", &IsotropicFlag, "IsotropicFlag/I");

  //--- open file -----
  ifstream inFile(filename.c_str());
  Long64_t n_evt = 0;
  while(inFile >> r >> Nucleus >> CosTheta >> SinPhi >> Er >> IsotropicFlag)
  {
    n_evt += 1;
    if(n_evt%10000==0) std::cout << "\rFill Data into Tree ... : Evevt # " << n_evt << " " << flush;
    tree->Fill();
  }
  cout << "\rFill Data into Tree ... : Evevt # " << n_evt << " " << flush;
  tree->Write("", TObject::kOverwrite);
  f->Close();
  delete f;
  cout << endl;
  return 0;
}
