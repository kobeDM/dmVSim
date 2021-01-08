// STL
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cstdlib> // rand()
#include <ctime>   // time()

// ROOT
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>
#include <TMath.h>
#include <TApplication.h>

// My Header
#include "MyKolmogorovTest1D.h"
#include "MyKolmogorovTest2D.h"
#include "param.h"

using namespace std;
//bool DEBUG = true;
  bool DEBUG = false;

string IntToString(int integer)
{
  stringstream ss;
  ss << integer;
  return ss.str();
}

int StringToInt(string str)
{
   istringstream iss(str);
   int integer;
   iss >> integer;
   return integer;
}

double myRound(double d_in, int n_len)
{
  double d_out;
  d_out = d_in * pow(10.0, n_len);
  d_out = (double)(int)(d_out + 0.5);
  return d_out * pow(10.0, -n_len);
}

//==============================================
// runAnalysis
//==============================================
int main(int argc, char *argv[]){


  
  // Input Parameter
  if(argc != 7) {
    cerr << "<ERROR>    Specify input file name as first argument." << endl;
    cerr << "./main [r0.hist] [r25.root] [nucleus number] [energy threshold(keV)] [Headtail (ON or OFF] [sample number] " << endl;
    return 1;
  }
  string filename_r0 = argv[1];
  string filename_r25 = argv[2];
  string nucleus  = argv[3];
  string ene_thr  = argv[4];
  string headtail = argv[5];
  string n_sample = argv[6];
  double f_ene_thr  = atof(argv[4]);
  int    f_n_sample = (unsigned long)atof(argv[6]);

  TApplication app("app",&argc,argv);

  // Search .dat Filename and Set ROOT filename
  string::size_type index = filename_r0.find(".root");
  if( index == string::npos ) {
    cout << "<ERROR>    Inputfile is not found." << endl;
    return 1;
  }
  index = filename_r25.find(".root");
  if( index == string::npos ) {
    cout << "<ERROR>    Inputfile is not found." << endl;
    return 1;
  }


  //=======================================//
  // Open Tree ROOT file (r=0)
  //=======================================//
  TFile* pfile_r0 = new TFile(filename_r0.c_str());
  // HEADTAIL
  /*
    if(headtail == "OFF"){
      CosTheta = abs(CosTheta);
    }
  */
  //TH1D* hCosTheta_r0      = (TH1D*)pfile_r0->Get("CosTheta");
  TH1D* hEr_r0            = (TH1D*)pfile_r0->Get("Er");
  TH2D* hEr_vs_CosTheta_r0 = (TH2D*)pfile_r0->Get("hEr_vs_CosTheta");




  //-------Energy CUT------//
  for(int i=0;i<nbin_Er;i++){
    double value=hEr_vs_CosTheta_r0->GetYaxis()->GetBinCenter(i+1);
    if(value < f_ene_thr){
      for(int j=0;j<nbin_CosTheta;j++){
	int ent=hEr_vs_CosTheta_r0->GetBinContent(j+1,i+1);
	hEr_vs_CosTheta_r0->SetEntries(hEr_vs_CosTheta_r0->GetEntries()-ent);
	hEr_vs_CosTheta_r0->SetBinContent(j+1,i+1,0.0);
	hEr_vs_CosTheta_r0->SetBinError(j+1,i+1,0);

      }
    }
  }

  for(int i=0;i<nbin_Er;i++){
    double value=hEr_r0->GetXaxis()->GetBinCenter(i+1);
    if(value < f_ene_thr){
      int ent=hEr_r0->GetBinContent(i+1);
      hEr_r0->SetEntries(hEr_r0->GetEntries()-ent);
      hEr_r0->SetBinContent(i+1,0.0);
    }
  }

  //  TH1D* hCosTheta_r0=(TH1D*)hEr_vs_CosTheta_r0->ProjectionX("CosTheta_r0",bin_min_CosTheta,bin_max_CosTheta);

  // please use this FILL methode
  TH1D* hCosTheta_r0=new TH1D("CosTheta_r0","",nbin_CosTheta,bin_min_CosTheta,bin_max_CosTheta);
  for(int j=0;j<nbin_CosTheta;j++){
    double val=0;
    for(int i=0;i<nbin_Er;i++){
      val+=hEr_vs_CosTheta_r0->GetBinContent(j+1,i+1);
    }
    hCosTheta_r0->SetBinContent(j+1,val);
    hCosTheta_r0->SetEntries(hCosTheta_r0->GetEntries()+val);
  }



  //=======================================//
  // Open Tree ROOT file (r=25)
  //=======================================//
  TFile* pfile_r25 = new TFile(filename_r25.c_str());
  // HEADTAIL
  /*
    if(headtail == "OFF"){
      CosTheta = abs(CosTheta);
    }
  */
  TH1D* hEr_r25            = (TH1D*)pfile_r25->Get("Er");
  TH2D* hEr_vs_CosTheta_r25 = (TH2D*)pfile_r25->Get("hEr_vs_CosTheta");

  //--------- energy CUT ---------//
  for(int i=0;i<nbin_Er;i++){
    double value=hEr_vs_CosTheta_r25->GetYaxis()->GetBinCenter(i+1);
    if(value < f_ene_thr){
      for(int j=0;j<nbin_CosTheta;j++){
	int ent=hEr_vs_CosTheta_r25->GetBinContent(j+1,i+1);
	hEr_vs_CosTheta_r25->SetEntries(hEr_vs_CosTheta_r25->GetEntries()-ent);
	hEr_vs_CosTheta_r25->SetBinContent(j+1,i+1,0.0);
	hEr_vs_CosTheta_r25->SetBinError(j+1,i+1,0);
      }
    }
  }

  for(int i=0;i<nbin_Er;i++){
    double value=hEr_r25->GetXaxis()->GetBinCenter(i+1);
    if(value < f_ene_thr) {
      int ent=hEr_r25->GetBinContent(i+1);
      hEr_r25->SetEntries(hEr_r25->GetEntries()-ent);
      hEr_r25->SetBinContent(i+1,0);
    }
  }
  //  TH1D* hCosTheta_r25=(TH1D*)hEr_vs_CosTheta_r25->ProjectionX("CosTheta_r25",bin_min_CosTheta,bin_max_CosTheta);

  TH1D* hCosTheta_r25=new TH1D("CosTheta_r25","",nbin_CosTheta,bin_min_CosTheta,bin_max_CosTheta);
  for(int j=0;j<nbin_CosTheta;j++){
    double val=0;
    for(int i=0;i<nbin_Er;i++){
      val+=hEr_vs_CosTheta_r25->GetBinContent(j+1,i+1);
    }
    hCosTheta_r25->SetBinContent(j+1,val);
    hCosTheta_r25->SetEntries(hCosTheta_r25->GetEntries()+val);
  }




  //=================================//
  // Extract events for pseudo data
  //=================================//
  //Define Histogram for Pseudo-Experiment Data
  TH1D* hCosTheta_pseudo_r = new TH1D("hCosTheta_pseudo","hCosTheta_pseudo", nbin_CosTheta, bin_min_CosTheta, bin_max_CosTheta);
  TH2D* hEr_vs_CosTheta_pseudo_r = new TH2D("hEr_vs_CosTheta_pseudo","hEr_vs_CosTheta_pseudo", nbin_CosTheta, bin_min_CosTheta, bin_max_CosTheta, nbin_Er, bin_min_Er, bin_max_Er);


  int sim_pseudo_data_n_r = 0;
  while(sim_pseudo_data_n_r < f_n_sample){
    double cos=hCosTheta_r25->GetRandom();
    double ene=hEr_r25->GetRandom();
    hEr_vs_CosTheta_pseudo_r->Fill(cos,ene);
    hCosTheta_pseudo_r->Fill(cos);
    sim_pseudo_data_n_r +=1;
  }



  //=================//
  //--- Main Loop ---//
  //=================//
  // Search common empty bin of template simulation data for r : 0 ~ 10
  //-- 1-D CosTheta
  vector<int> empty_bin_cos;
  for(int xBin=1;xBin<=nbin_CosTheta;xBin++){
    bool empty_flag = false;
    if(hCosTheta_r0->GetBinContent(xBin)!=0) { empty_flag = false;}
    else                                    { empty_flag = true;        }

    if(empty_flag==true){
      empty_bin_cos.push_back(xBin);
      if(DEBUG) cout << "CosTheta empty bin : x = (" << xBin << ")" << endl;
    }
  }
  int n_empty_bin_cos = empty_bin_cos.size();
  if(DEBUG) cout << "number of CosTheta empty bin : " << n_empty_bin_cos << endl;

  //-- 2-D Er vs. CosTheta
  vector<int> empty_bin_x;
  vector<int> empty_bin_y;
  for(int xBin=1;xBin<=nbin_CosTheta;xBin++){
    for(int yBin=1;yBin<=nbin_Er;yBin++){
      bool empty_flag = false;
      if(hEr_vs_CosTheta_r0->GetBinContent(xBin, yBin)!=0) { empty_flag = false;}
      else                                                        { empty_flag = true;        }

      if(empty_flag==true){
        empty_bin_x.push_back(xBin);
        empty_bin_y.push_back(yBin);
        if(DEBUG) cout << "empty bin : (x, y) = (" << xBin <<  ", " << yBin << ")" << endl;
      }
    }
  }
  int n_empty_bin = empty_bin_x.size();
  if(DEBUG) cout << "number of empty bin : " << n_empty_bin << endl;



  // Define Histogram for Scaled Template Data
  TH1D* hCosTheta_scaled_r= (TH1D*)hCosTheta_r0->Clone("CosTheta_scaled");
  hCosTheta_scaled_r->SetNameTitle("CosTheta_scaled","CosTheta_scaled");
  double scale_factor_cos = (double)hCosTheta_pseudo_r->Integral()/hCosTheta_r0->Integral();
  hCosTheta_scaled_r->Scale(scale_factor_cos);

  TH2D* hEr_vs_CosTheta_scaled_r = new TH2D("Er_vs_CosTheta_scaled","Er_vs_CosTheta_scaled", nbin_CosTheta, bin_min_CosTheta, bin_max_CosTheta, nbin_Er, bin_min_Er, bin_max_Er);
  // Set Scale Factor for Template Data
  double scale_factor_e_cos_r = (double)hEr_vs_CosTheta_pseudo_r->Integral()/hEr_vs_CosTheta_r0->Integral();




  //=================================================//
  // Chi Squre Test
  //=================================================//
  vector<double> sim_pseudo_data_Er_r;        
  vector<double> sim_pseudo_data_CosTheta_r;  
  vector<double> sim_pseudo_data_BinContent_r;
  vector<double> sim_pseudo_data_BinError_r;

  double chi2_cos_r=0;
  double delta_chi2_cos_r=0;
  double chi2_r=0;
  double delta_chi2_r=0;
                                           //n :      0,     1,     2,     3,     4,     5,     6,      7,      8,      9,     10
  double poisson_twosided_68CI_error_lower[11] = {0.000, 0.174, 0.712, 1.373, 2.093, 2.849, 3.630,  4.429,  5.243,  6.069,  6.905};
  double poisson_twosided_68CI_error_upper[11] = {1.833, 3.289, 4.625, 5.904, 7.147, 8.365, 9.566, 10.751, 11.925, 13.089, 14.245};

  //------------ 1-D CosTheta ------------//
  for(int xBin=1;xBin<=nbin_CosTheta;xBin++){
    int tmp_bin_content_pseudo_r=0;
    double tmp_bin_content_scaled_r=0;
    tmp_bin_content_pseudo_r = hCosTheta_pseudo_r->GetBinContent(xBin);
    tmp_bin_content_scaled_r = hCosTheta_scaled_r->GetBinContent(xBin);

    // Set Error for Scaled Template Data
    hCosTheta_scaled_r->SetBinError(xBin, 0);// Notice that error of scaled template data equal zero.

    bool skip_chi2_cos_flag = false;
    for(int iEmptyBin=0;iEmptyBin<n_empty_bin_cos;iEmptyBin++){
      if(xBin==empty_bin_cos[iEmptyBin]){
      	skip_chi2_cos_flag = true;
      	break;
      }
    }
    if(skip_chi2_cos_flag==true) continue;
    // Calulate chi square // Notice that error of pseudo-experiment data histgram equal sqrt(n) for n>10, poisson error for 0<=n<=10.
    int    mes = tmp_bin_content_pseudo_r;
    double exp = tmp_bin_content_scaled_r;
    if(mes<0){
      cout << "<WARNING>    Negative entry bin is found." << endl;
    } else if(mes>=0 && mes<=10){
      if(exp>=mes){
	chi2_cos_r += (pow((double)mes-exp, 2.0)/pow(poisson_twosided_68CI_error_upper[mes]-mes, 2.0));
      } else {
	chi2_cos_r += (pow((double)mes-exp, 2.0)/pow(mes-poisson_twosided_68CI_error_lower[mes], 2.0));
      }
    } else if(mes>10){
      chi2_cos_r += (pow((double)mes-exp, 2.0)/mes);
    }
  }// End of loop for xBin



  //-------------- 2-D Er vs. CosTheta ---------//
  for(int xBin=1;xBin<=nbin_CosTheta;xBin++){
    for(int yBin=1;yBin<=nbin_Er;yBin++){
      int tmp_bin_content_pseudo_r=0;
      double tmp_bin_content_scaled_r=0;
      tmp_bin_content_pseudo_r = hEr_vs_CosTheta_pseudo_r->GetBinContent(xBin, yBin);
      tmp_bin_content_scaled_r  = scale_factor_e_cos_r*hEr_vs_CosTheta_r0->GetBinContent(xBin, yBin);

      
      // Set Entry and Error for Scaled Template Data
      hEr_vs_CosTheta_scaled_r->SetBinContent(xBin, yBin, tmp_bin_content_scaled_r);
      hEr_vs_CosTheta_scaled_r->SetBinError(xBin, yBin, 0);// Notice that error of scaled template data equal zero.

      bool skip_chi2_flag = false;
      for(int iEmptyBin=0;iEmptyBin<n_empty_bin;iEmptyBin++){
      	if(xBin==empty_bin_x[iEmptyBin] && yBin==empty_bin_y[iEmptyBin]){
      	  skip_chi2_flag = true;
      	  break;
        }
      }
      if(skip_chi2_flag==true) continue;

      /*
      sim_pseudo_data_Er_r.push_back((hEr_vs_CosTheta_pseudo_r->GetXaxis())->GetBinCenter(xBin));
      sim_pseudo_data_CosTheta_r.push_back((hEr_vs_CosTheta_pseudo_r->GetYaxis())->GetBinCenter(yBin));
      sim_pseudo_data_BinContent_r.push_back(tmp_bin_content_pseudo_r);
      sim_pseudo_data_BinError_r.push_back(sqrt(tmp_bin_content_pseudo_r));// Notice that error of pseudo-experiment data graph(not but histgram) equal sqrt(n).
      */

      double tmp_scaled_bin_content_r=0;
      tmp_scaled_bin_content_r = hEr_vs_CosTheta_scaled_r->GetBinContent(xBin, yBin);


      // Calulate chi square // Notice that error of pseudo-experiment data histgram equal sqrt(n) for n>10, poisson error for 0<=n<=10.
      int    mes = tmp_bin_content_pseudo_r;
      double exp = tmp_scaled_bin_content_r;
      if(mes<0){
	cout << "<WARNING>    Negative entry bin is found." << endl;
      } else if(mes>=0 && mes<=10){
	if(exp>=mes){
	  chi2_r += (pow((double)mes-exp, 2.0)/pow(poisson_twosided_68CI_error_upper[mes]-mes, 2.0));
	} else {
	  chi2_r += (pow((double)mes-exp, 2.0)/pow(mes-poisson_twosided_68CI_error_lower[mes], 2.0));
	}
      } else if(mes>10){
	chi2_r += (pow((double)mes-exp, 2.0)/mes);
      }
    }// End of loop for yBin
  }// End of loop for xBin


  // Calculate reduced chi square
  //1-D CosTheta
  double redu_chi2_cos_r=chi2_cos_r/(nbin_CosTheta-n_empty_bin_cos);
  //2-D Er vs. CosTheta
  double redu_chi2_r=chi2_r/((nbin_CosTheta*nbin_Er)-n_empty_bin);



  // Calculate delta chi square
  //1-D CosTheta
  double tmp_min_chi2_cos = chi2_cos_r;
  if(chi2_cos_r < tmp_min_chi2_cos) tmp_min_chi2_cos = chi2_cos_r;
  delta_chi2_cos_r = chi2_cos_r - tmp_min_chi2_cos;
  //2-D Er vs. CosTheta
  double tmp_min_chi2 = chi2_r;


  if(chi2_r < tmp_min_chi2) tmp_min_chi2 = chi2_r;
  delta_chi2_r = chi2_r - tmp_min_chi2;


  //========================================//
  // OUTPUT
  //========================================//
  cout << f_n_sample << "\t" << f_ene_thr << "\t" << (nbin_CosTheta-n_empty_bin_cos) << "\t" 
       << chi2_cos_r << "\t" << ((nbin_CosTheta*nbin_Er)-n_empty_bin) << "\t" << chi2_r << endl;

  if(DEBUG){
  TCanvas *c=new TCanvas("c","",800,800);
  c->Divide(2,2);
  c->cd(1);
  hEr_vs_CosTheta_r0->Draw("colz");
  c->cd(2);
  hEr_vs_CosTheta_r25->Draw("colz");
  c->cd(3);
  hCosTheta_r0->Draw("E");
  //  hCosTheta_pseudo_r->Draw();
  //  hCosTheta_scaled_r->Draw("same");
  //  hEr_r0->Draw();
  app.Run();
  }
  return 0;
}
