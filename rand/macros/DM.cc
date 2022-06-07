//20120524 nakamura debug
//20110820 nakamura Copy from CoGeNT
//20191011 ikeda debug for ROOT
#include "func.h"

const double M_D=100; // [GeV]

void DM( )
{


  //#include <fstream.h>
  //#include "/home/miuchi/root_macro/cpalette.cc"
  //cpalette();

  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  int canvas_flag = 1;


  Double_t v_E[13] = { 244.0 , 233.4 /*Jan*/, 240.0 , 247.4 , 253.7 , 257.2 , 257.4 /*Jun*/,
		       254.3 , 248.5 , 241.4 , 234.6 , 230.0 , 229.5 }; // Earth velocity [km/sec]
  Double_t k_0     = pow( PI*v_0*v_0 , 1.5 );
  Double_t k_1     = k_0 * ( TMath::Erf(v_esc/v_0) - 2/sqrt(PI) * v_esc/v_0 * exp(-v_esc*v_esc/v_0/v_0));
  Double_t const0  = k_0/k_1;
  Double_t const1[13];
  Double_t const2[13];
  Double_t const3  = exp(-v_esc*v_esc/v_0/v_0);
  for(int i=0;i<13;i++){
    const1[i] = sqrt(PI)/4.0*v_0/v_E[i];
    const2[i] = v_E[i]/v_0;
  }
  std::cout << "test" << std::endl;

  //**********************************************************************
  //** Normalized Energy Spectrum ****************************************
  //**********************************************************************
  TCanvas *c1;
  c1 = new TCanvas("c1", "c1", 25, 25, 700, 600);
  c1->SetFillColor(10);
  c1->SetGridx(1);
  c1->SetGridy(1);
  c1->SetLogx(0);
  c1->SetLogy(0);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);

  TF1 *spec_nor_ave = new TF1("spec_nor_ave","func_spec_nor(x,[0],[1],[2],[3])",0,10);
  spec_nor_ave->SetLineColor(2);
  spec_nor_ave->SetParameters(const0, const1[0], const2[0], const3);
  TF1 *spec_nor_Jun = new TF1("spec_nor_Lun","func_spec_nor(x,[0],[1],[2],[3])",0,10);
  spec_nor_Jun->SetLineColor(3);
  spec_nor_Jun->SetParameters(const0, const1[6], const2[6], const3);
  TF1 *spec_nor_Dec = new TF1("spec_nor_Dec","func_spec_nor(x,[0],[1],[2],[3])",0,10);
  spec_nor_Dec->SetLineColor(4);
  spec_nor_Dec->SetParameters(const0, const1[12], const2[12], const3);

  TH1F *waku1 = new TH1F("waku1","normalized spectrum",1,0,10);
  waku1->SetMinimum(0);
  waku1->SetMaximum(1);
  waku1->GetXaxis()->SetTitle("#frac{E}{E_{0}r} (normalized energy)");
  waku1->GetYaxis()->SetTitle("#frac{E_{0}r}{R_{0}} #frac{dR}{dE} (normalized event rate)");
  waku1->GetXaxis()->SetTitleOffset(1.3);
  waku1->GetYaxis()->SetTitleOffset(1.4);

  TLegend *legend1 = new TLegend(0.6,0.6,0.88,0.88);
  legend1->AddEntry(spec_nor_ave,"average","l");
  legend1->AddEntry(spec_nor_Jun,"June","l");
  legend1->AddEntry(spec_nor_Dec,"December","l");

  waku1->Draw("");
  spec_nor_ave->Draw("same");
  spec_nor_Jun->Draw("same");
  spec_nor_Dec->Draw("same");
  legend1->Draw("same");
  // c1->Print("energy_spectrum_kikakuka.eps");


  //**********************************************************************
  //** Form Factor *******************************************************
  //**********************************************************************



  /**********************************************************************/
  /** Nuclear rms charge radii ******************************************/
  /**********************************************************************/
  /*
	TCanvas *c3;
  c3 = new TCanvas("c3", "c3", 75, 75, 1000, 650);
  c3->SetFillColor(10);
  c3->SetGridx(1);
  c3->SetGridy(1);
  c3->SetLogx(0);
  c3->SetLogy(0);

  TF1 *radii1 = new TF1("radii1","sqrt(3.0/5.0*pow((1.23*x**(1.0/3.0)-0.6),2)+7.0/5.0*[1]*[1]*[0]*[0])",0,250);
  radii1->SetLineColor(2);
  double a_3 = 0.52;
  radii1->SetParameters(a_3,PI);

  TF1 *radii2 = new TF1("radii2","sqrt(3.0/5.0*[0]*[0]*x**(2.0/3.0)+3.0*[1]*[1])",0,250);
  radii2->SetLineColor(4);
  double r_n_3 = 1.0;//1.19; <- to match the shape of graph
  double s_3 = 1.0;//0.3; <- to match the shape of graph
  radii2->SetParameters(r_n_3,s_3);

  TF1 *radii3 = new TF1("radii3","0.89*x**(1.0/3.0)+0.3",0,250);
  radii3->SetLineColor(3);

  TH1F *waku3 = new TH1F("waku3","Nuclear rms charge radii",1,0,250);
  waku3->SetMinimum(0);
  waku3->SetMaximum(6);
  waku3->GetXaxis()->SetTitle("A [atomic number]");
  waku3->GetYaxis()->SetTitle("r_{rms}");

  TLegend *legend3 = new TLegend(0.6,0.3,0.88,0.39);
  legend3->AddEntry(radii1,"least-squares fit to c","l");
  legend3->AddEntry(radii2,"Engal [15] fit","l");
  legend3->AddEntry(radii3,"Eder [23] fit","l");

  waku3->Draw("");
  radii1->Draw("same");
  radii2->Draw("same");
  radii3->Draw("same");
  legend3->Draw("same");
	*/





  //**********************************************************************
  //** Cross Section  ****************************************************
  //**********************************************************************
  /*
	TCanvas *c4;
  c4 = new TCanvas("c4", "c4", 100, 100, 850, 450);
  c4->SetFillColor(10);

  TPad *pad41 = new TPad("pad41","pad41",0.0,0.0,0.5,1.0);
  pad41->SetFillColor(10);
  pad41->SetGridx(1);
  pad41->SetGridy(1);
  pad41->SetLogx(1);
  pad41->SetLogy(1);
  pad41->SetLeftMargin(0.15);
  pad41->SetRightMargin(0.05);
  pad41->SetBottomMargin(0.11);
  pad41->Draw();

  TPad *pad42 = new TPad("pad42","pad42",0.5,0.0,1.0,1.0);
  pad42->SetFillColor(10);
  pad42->SetGridx(1);
  pad42->SetGridy(1);
  pad42->SetLogx(1);
  pad42->SetLogy(1);
  pad42->SetLeftMargin(0.15);
  pad42->SetRightMargin(0.05);
  pad42->SetBottomMargin(0.11);
  pad42->Draw();

  // SI cross section
  TF1 *cross_SI = new TF1("cross_SI","func_cross_SI(x,[0])",10,1000);

  // SD cross section
  TF1 *cross_SD = new TF1("cross_SD","func_cross_SD(x,[0],[1])",10,1000);

  // material
  TF1* cross_SI_Ge = (TF1*)cross_SI->Clone("cross_SI_Ge");
  TF1* cross_SI_Xe = (TF1*)cross_SI->Clone("cross_SI_Xe");
  cross_SI->SetParameter(0,A);
  cross_SI_Ge->SetParameter(0,73);
  cross_SI_Xe->SetParameter(0,128);

  TF1* cross_SD_Ge = (TF1*)cross_SD->Clone("cross_SD_Ge");
  TF1* cross_SD_Xe = (TF1*)cross_SD->Clone("cross_SD_Xe");
  TF1* cross_SD_Xe2= (TF1*)cross_SD->Clone("cross_SD_Xe2");
  TF1* cross_SD_Li = (TF1*)cross_SD->Clone("cross_SD_Li");
  TF1* cross_SD_Na = (TF1*)cross_SD->Clone("cross_SD_Na");
  TF1* cross_SD_I  = (TF1*)cross_SD->Clone("cross_SD_I");
  TF1* cross_SD_Si = (TF1*)cross_SD->Clone("cross_SD_Si");
  TF1* cross_SD_W  = (TF1*)cross_SD->Clone("cross_SD_W");
  cross_SD->SetParameters(A, Ja);
  cross_SD_Ge->SetParameters(73, Ja_Ge);
  cross_SD_Xe->SetParameters(129,Ja_Xe);
  cross_SD_Xe2->SetParameters(131,Ja_Xe2);
  cross_SD_Li->SetParameters(7,  Ja_Li);
  cross_SD_Na->SetParameters(23, Ja_Na);
  cross_SD_I->SetParameters(127, Ja_I);
  cross_SD_Si->SetParameters(29, Ja_Si);
  cross_SD_W->SetParameters(183, Ja_W);

  // waku & draw
  TH1F *waku41 = new TH1F("waku41","",1,10,1000);
  waku41->SetMinimum(10000);
  waku41->SetMaximum(1000000000);
  waku41->GetXaxis()->SetTitle("M_{D} [GeV/c^{2}]");
  waku41->GetYaxis()->SetTitle("#sigma^{SI}_{DM-N} / #sigma^{SI}_{DM-p}");
  waku41->GetXaxis()->SetTitleOffset(1.2);
  waku41->GetYaxis()->SetTitleOffset(1.4);
  TLegend *legend41 = new TLegend(0.2,0.66,0.32,0.84);
  legend41->AddEntry(cross_SI,   "^{19}F","l");
  legend41->AddEntry(cross_SI_Ge,"^{73}Ge","l");
  legend41->AddEntry(cross_SI_Xe,"^{128}Xe","l");
  pad41->cd();
  waku41->Draw("");
  legend41->Draw("same");
  cross_SI->SetLineColor(2);
  cross_SI->Draw("same");
  cross_SI_Ge->SetLineColor(4);
  cross_SI_Ge->Draw("same");
  cross_SI_Xe->SetLineColor(3);
  cross_SI_Xe->Draw("same");

  TH1F *waku42 = new TH1F("waku42","",1,10,1000);
  waku42->SetMinimum(0.1);
  waku42->SetMaximum(1000);
  waku42->GetXaxis()->SetTitle("M_{D} [GeV/c^{2}]");
  waku42->GetYaxis()->SetTitle("#sigma^{SD}_{DM-N} / #sigma^{SD}_{DM-p}");
  waku42->GetXaxis()->SetTitleOffset(1.2);
  waku42->GetYaxis()->SetTitleOffset(1.4);
  TLegend *legend42 = new TLegend(0.8,0.12,0.93,0.5);
  legend42->AddEntry(cross_SD_Li,"^{7}Li","l");
  legend42->AddEntry(cross_SD   ,"^{19}F","l");
  legend42->AddEntry(cross_SD_Na,"^{23}Na","l");
  legend42->AddEntry(cross_SD_Si,"^{29}Si","l");
  legend42->AddEntry(cross_SD_Ge,"^{73}Ge","l");
  legend42->AddEntry(cross_SD_I, "^{127}I","l");
  legend42->AddEntry(cross_SD_Xe,"^{129}Xe","l");
  legend42->AddEntry(cross_SD_Xe2,"^{131}Xe","l");
  pad42->cd();
  waku42->Draw("");
  legend42->Draw("same");
  cross_SD->SetLineColor(2);
  cross_SD->Draw("same");
  cross_SD_Li->SetLineColor(1);
  cross_SD_Li->Draw("same");
  cross_SD_Na->SetLineColor(4);
  cross_SD_Na->Draw("same");
  cross_SD_I->SetLineColor(3);
  cross_SD_I->Draw("same");
  cross_SD_Ge->SetLineColor(2);
  cross_SD_Ge->SetLineStyle(2);
  cross_SD_Ge->Draw("same");
  cross_SD_Xe->SetLineColor(4);
  cross_SD_Xe->SetLineStyle(2);
  cross_SD_Xe->Draw("same");
  cross_SD_Xe2->SetLineColor(5);
  cross_SD_Xe2->SetLineStyle(2);
  cross_SD_Xe2->Draw("same");
  cross_SD_Si->SetLineColor(1);
  cross_SD_Si->SetLineStyle(2);
  cross_SD_Si->Draw("same");
  cross_SD_W->SetLineColor(3);
  cross_SD_W->SetLineStyle(2);
  cross_SD_W->Draw("same");
	*/




  //**********************************************************************
  //** Energy Spectrum *************************************
  //**********************************************************************
  TCanvas *c5;
  c5 = new TCanvas("c5", "c5", 125, 125, 1000, 600);
  c5->SetFillColor(10);

  // pad
  TPad *pad51 = new TPad("pad51","pad51",0.0,0.5,0.33,1.0);
  pad51->SetFillColor(10);
  pad51->SetGridx(1);
  pad51->SetGridy(1);
  pad51->SetLogx(0);
  pad51->SetLogy(1);
  pad51->Draw();

  TPad *pad52 = new TPad("pad52","pad52",0.0,0.0,0.33,0.5);
  pad52->SetFillColor(10);
  pad52->SetGridx(1);
  pad52->SetGridy(1);
  pad52->SetLogx(0);
  pad52->SetLogy(1);
  pad52->Draw();

  TPad *pad53 = new TPad("pad53","pad53",0.33,0.5,0.66,1.0);
  pad53->SetFillColor(10);
  pad53->SetGridx(1);
  pad53->SetGridy(1);
  pad53->SetLogx(0);
  pad53->SetLogy(1);
  pad53->Draw();

  TPad *pad54 = new TPad("pad54","pad54",0.33,0.0,0.66,0.5);
  pad54->SetFillColor(10);
  pad54->SetGridx(1);
  pad54->SetGridy(1);
  pad54->SetLogx(0);
  pad54->SetLogy(1);
  pad54->Draw();

  TPad *pad55 = new TPad("pad55","pad55",0.66,0.5,1.0,1.0);
  pad55->SetFillColor(10);
  pad55->SetGridx(1);
  pad55->SetGridy(1);
  pad55->SetLogx(0);
  pad55->SetLogy(1);
  pad55->Draw();

  TPad *pad56 = new TPad("pad56","pad56",0.66,0.0,1.0,0.5);
  pad56->SetFillColor(10);
  pad56->SetGridx(1);
  pad56->SetGridy(1);
  pad56->SetLogx(0);
  pad56->SetLogy(1);
  pad56->Draw();

  TF1 *spec_SI = new TF1("spec_SI","const5([4],[5],[6],[7])*func_spec_nor(x/const4([4],[5],[6]),[0],[1],[2],[3])*func_FSI_keV(x,[4])*[8]*func_cross_SI([5],[4])",0,200);

  TF1 *spec_SD = new TF1("spec_SD","const5([4],[5],[6],[7])*func_spec_nor(x/const4([4],[5],[6]),[0],[1],[2],[3])*func_FSD_keV(x,[4])*[8]*func_cross_SD([5],[4],[9])",0,1000);


  // *******SI 100GeV*******
  TF1* spec_SI_100_F  = (TF1*)spec_SI->Clone("spec_SI_100_F");
  TF1* spec_SI_100_Ge = (TF1*)spec_SI->Clone("spec_SI_100_Ge");
  TF1* spec_SI_100_Xe = (TF1*)spec_SI->Clone("spec_SI_100_Xe");
  TF1* spec_SI_100_Na = (TF1*)spec_SI->Clone("spec_SI_100_Na");
  TF1* spec_SI_100_Ar = (TF1*)spec_SI->Clone("spec_SI_100_Ar");
  TF1* spec_SI_100_I  = (TF1*)spec_SI->Clone("spec_SI_100_I");
  spec_SI_100_F->SetParameters(  const0, const1[0] ,const2[0] ,const3 , A,  M_D, v_0, rho, sigma_SI_0);
  spec_SI_100_Ge->SetParameters( const0, const1[0] ,const2[0] ,const3 , 74, M_D, v_0, rho, sigma_SI_0);
  spec_SI_100_Xe->SetParameters( const0, const1[0] ,const2[0] ,const3 , 132,M_D, v_0, rho, sigma_SI_0);
  spec_SI_100_Na->SetParameters( const0, const1[0] ,const2[0] ,const3 , 22, M_D, v_0, rho, sigma_SI_0);
  spec_SI_100_Ar->SetParameters( const0, const1[0] ,const2[0] ,const3 , 39, M_D, v_0, rho, sigma_SI_0);
  spec_SI_100_I->SetParameters(  const0, const1[0] ,const2[0] ,const3 , 127,M_D, v_0, rho, sigma_SI_0);

  // *******SD 100GeV*******
  TF1* spec_SD_100_F =  (TF1*)spec_SD->Clone("spec_SD_100_F");
  TF1* spec_SD_100_Ge = (TF1*)spec_SD->Clone("spec_SD_100_Ge");
  TF1* spec_SD_100_Xe = (TF1*)spec_SD->Clone("spec_SD_100_Xe");
  TF1* spec_SD_100_Li = (TF1*)spec_SD->Clone("spec_SD_100_Li");
  TF1* spec_SD_100_Na = (TF1*)spec_SD->Clone("spec_SD_100_Na");
  TF1* spec_SD_100_I =  (TF1*)spec_SD->Clone("spec_SD_100_I");
  TF1* spec_SD_100_F_Jun = (TF1*)spec_SD->Clone("spec_SD_100_F_Jun");
  TF1* spec_SD_100_F_Dec = (TF1*)spec_SD->Clone("spec_SD_100_F_Dec");
  spec_SD_100_F->SetParameters(  const0, const1[0] ,const2[0] ,const3 , A,  M_D, v_0, rho, sigma_SD_0, Ja);
  spec_SD_100_Ge->SetParameters( const0, const1[0] ,const2[0] ,const3 , 73, M_D, v_0, rho, sigma_SD_0, Ja_Ge);
  spec_SD_100_Xe->SetParameters( const0, const1[0] ,const2[0] ,const3 , 131,M_D, v_0, rho, sigma_SD_0, Ja_Xe);
  spec_SD_100_Li->SetParameters( const0, const1[0] ,const2[0] ,const3 , 7,  M_D, v_0, rho, sigma_SD_0, Ja_Li);
  spec_SD_100_Na->SetParameters( const0, const1[0] ,const2[0] ,const3 , 23, M_D, v_0, rho, sigma_SD_0, Ja_Na);
  spec_SD_100_I->SetParameters(  const0, const1[0] ,const2[0] ,const3 , 127,M_D, v_0, rho, sigma_SD_0, Ja_I);
  spec_SD_100_F_Jun->SetParameters( const0, const1[6] ,const2[6] ,const3   , A, M_D, v_0, rho, sigma_SD_0, Ja);
  spec_SD_100_F_Dec->SetParameters( const0, const1[12] ,const2[12] ,const3 , A, M_D, v_0, rho, sigma_SD_0, Ja);

  TF1* spec_SD_100_F_annual[12];
  double month[12];
  double rate_annual[12];

  for(int i=0;i<12;i++){
    spec_SD_100_F_annual[i]=(TF1*)spec_SD->Clone(Form("spec_SD_100_F_%d",i));
    spec_SD_100_F_annual[i]->SetParameters( const0, const1[i+1] ,const2[i+1] ,const3 , A, M_D, v_0, rho, sigma_SD_0, Ja);
    spec_SD_100_F_annual[i]->SetNpx(1000);
    rate_annual[i]=spec_SD_100_F_annual[i]->Eval(55);
    month[i]=i+1;
    cerr << rate_annual[i] << endl;
  }



  // thesis
  TGraph *g_annual=new TGraph(12,month,rate_annual);
  TCanvas *c_annual=new TCanvas("c_annual","c_annual",800,800);
  TH1F *frame_annual = new TH1F("frame_annual","",12,0,12);
  frame_annual->GetXaxis()->SetBinLabel(0,"Jan");
  frame_annual->GetXaxis()->SetBinLabel(5,"Jun");
  frame_annual->GetXaxis()->SetBinLabel(11,"Dec");
  frame_annual->GetXaxis()->SetLabelOffset(0.01);
  frame_annual->SetMinimum(-0.1);
  frame_annual->SetMaximum(0.1);
  frame_annual->GetXaxis()->SetNdivisions(505);
  //  frame_annual->GetXaxis()->SetTitle("Recoil energy (keV)");
  frame_annual->GetYaxis()->SetTitle("Event rate (counts/kg/days)");
  frame_annual->GetXaxis()->CenterTitle();
  frame_annual->GetYaxis()->CenterTitle();
  frame_annual->GetXaxis()->SetTitleFont(43);
  frame_annual->GetXaxis()->SetTitleSize(40);
  frame_annual->GetXaxis()->SetLabelFont(43);
  frame_annual->GetXaxis()->SetLabelSize(40);
  frame_annual->GetYaxis()->SetTitleFont(43);
  frame_annual->GetYaxis()->SetTitleSize(40);
  frame_annual->GetYaxis()->SetTitleOffset(1.5);
  frame_annual->GetYaxis()->SetLabelFont(43);
  frame_annual->GetYaxis()->SetLabelSize(40);
  frame_annual->Draw();

  //  g_annual->Draw("pl");
  TF1 *fit_annual=new TF1("fit_annual","[0]+[1]*cos(2*3.14*(x-[2])/[3])",0,13);
  fit_annual->SetParameter(0,0.21);
  fit_annual->SetParameter(1,0.01);
  fit_annual->SetParameter(2,0);
  fit_annual->SetParameter(3,12);
  g_annual->Fit(fit_annual);

  double mean=fit_annual->GetParameter(0);
  fit_annual->SetParameter(0,0);
  fit_annual->Draw("same l");

  //######################################
  // for D thesis
  //######################################
  TCanvas *c_spec_SD_annual=new TCanvas("c_spec_SD_annual","c_spec_SD_annual",800,800);
  //  c_spec_SD->SetLogy(1);
  TH1F *frame_spec_SD_annual = new TH1F("frame_spec_SD_annual","",1,0,200);
  frame_spec_SD_annual->SetMinimum(1e-3);
  frame_spec_SD_annual->SetMaximum(10);
  frame_spec_SD_annual->GetXaxis()->SetNdivisions(505);
  frame_spec_SD_annual->GetXaxis()->SetTitle("Recoil energy (keV)");
  frame_spec_SD_annual->GetYaxis()->SetTitle("Event rate (counts/keV/kg/days)");
  frame_spec_SD_annual->GetXaxis()->CenterTitle();
  frame_spec_SD_annual->GetYaxis()->CenterTitle();
  frame_spec_SD_annual->GetXaxis()->SetTitleFont(43);
  frame_spec_SD_annual->GetXaxis()->SetTitleSize(40);
  frame_spec_SD_annual->GetXaxis()->SetLabelFont(43);
  frame_spec_SD_annual->GetXaxis()->SetLabelSize(40);
  frame_spec_SD_annual->GetYaxis()->SetTitleFont(43);
  frame_spec_SD_annual->GetYaxis()->SetTitleSize(40);
  frame_spec_SD_annual->GetYaxis()->SetTitleOffset(1.5);
  frame_spec_SD_annual->GetYaxis()->SetLabelFont(43);
  frame_spec_SD_annual->GetYaxis()->SetLabelSize(40);
  frame_spec_SD_annual->Draw();

  spec_SD_100_F_Dec->SetNpx(1000);
  spec_SD_100_F_Dec->SetLineColor(kRed);
  spec_SD_100_F_Dec->SetLineWidth(3);
  spec_SD_100_F_Dec->Draw("same");
  spec_SD_100_F_Jun->SetNpx(1000);
  spec_SD_100_F_Jun->SetLineWidth(3);
  spec_SD_100_F_Jun->SetLineColor(kBlue);
  spec_SD_100_F_Jun->Draw("same");

  TLegend *leg_spec_SD_annual=new TLegend(0.5,0.6,0.8,0.8);
  leg_spec_SD_annual->SetBorderSize(0);
  leg_spec_SD_annual->SetTextFont(43);
  leg_spec_SD_annual->SetTextSize(30);
  leg_spec_SD_annual->AddEntry(spec_SD_100_F_Jun,"m_{#chi}=100 GeV","l");
  leg_spec_SD_annual->AddEntry(spec_SD_100_F_Dec,"m_{#chi}=100 GeV","l");
  leg_spec_SD_annual->Draw();
  //c_spec_SD_annual->Print("c_spec_SD_annual.eps");


  // *******SI 50GeV*******
  /*
  rho        = 0.3;
  M_D        = 50; // [GeV/c^2]
  sigma_SI_0 = 0.000001; //[pb]
  TF1* spec_SI_50_F = (TF1*)spec_SI->Clone("spec_SI_50_F");
  spec_SI_50_F->SetParameters( const0, const1[0] ,const2[0] ,const3 , A, M_D, v_0, rho, sigma_SI_0);
  */

  // *******SD 50GeV*******
  /*
  rho        = 0.3;
  M_D        = 50; // [GeV/c^2]
  sigma_SD_0 = 1; //[pb]
  */
  TF1* spec_SD_50_F = (TF1*)spec_SD->Clone("spec_SD_50_F");
  spec_SD_50_F->SetParameters( const0, const1[0] ,const2[0] ,const3 , A, 50, v_0, rho, sigma_SD_0, Ja);


  // *******SI 200GeV*******
  /*
  rho        = 0.3;
  M_D        = 200; // [GeV/c^2]
  sigma_SI_0 = 0.000001; //[pb]
  TF1* spec_SI_200_F = (TF1*)spec_SI->Clone("spec_SI_200_F");
  spec_SI_200_F->SetParameters( const0, const1[0] ,const2[0] ,const3 , A, M_D, v_0, rho, sigma_SI_0);
  */

  // *******SD 200GeV*******
  /*
  rho        = 0.3;
  M_D        = 200; // [GeV/c^2]
  sigma_SD_0 = 1; //[pb]
  */
  TF1* spec_SD_200_F = (TF1*)spec_SD->Clone("spec_SD_200_F");
  spec_SD_200_F->SetParameters( const0, const1[0] ,const2[0] ,const3 , A, 200, v_0, rho, sigma_SD_0, Ja);


  // *******SI 7GeV*******
  /*
  rho        = 0.3;
  M_D        = 7; // [GeV/c^2]
  sigma_SI_0 = 0.000001; //[pb]
  TF1* spec_SI_7_F = (TF1*)spec_SI->Clone("spec_SI_7_F");
  TF1* spec_SI_7_Ge = (TF1*)spec_SI->Clone("spec_SI_7_Ge");
  TF1* spec_SI_7_Xe = (TF1*)spec_SI->Clone("spec_SI_7_Xe");
  TF1* spec_SI_7_Ge_Jun = (TF1*)spec_SI->Clone("spec_SI_7_Ge_Jun");
  TF1* spec_SI_7_Ge_Dec = (TF1*)spec_SI->Clone("spec_SI_7_Ge_Dec");
  spec_SI_7_F->SetParameters( const0, const1[0] ,const2[0] ,const3 , A, M_D, v_0, rho, sigma_SI_0);
  spec_SI_7_Ge->SetParameters( const0, const1[0] ,const2[0] ,const3 , 73, M_D, v_0, rho, sigma_SI_0);
  spec_SI_7_Xe->SetParameters( const0, const1[0] ,const2[0] ,const3 , 131, M_D, v_0, rho, sigma_SI_0);
  spec_SI_7_Ge_Jun->SetParameters( const0, const1[6] ,const2[6] ,const3 , 73, M_D, v_0, rho, sigma_SI_0);
  spec_SI_7_Ge_Dec->SetParameters( const0, const1[12] ,const2[12] ,const3 , 73, M_D, v_0, rho, sigma_SI_0);
  */

  // *******SD 7GeV*******
  /*
  rho        = 0.3;
  M_D        = 7; // [GeV/c^2]
  sigma_SD_0 = 80; //[pb]
  TF1* spec_SD_7_F = (TF1*)spec_SD->Clone("spec_SD_7_F");
  TF1* spec_SD_7_Ge = (TF1*)spec_SD->Clone("spec_SD_7_Ge");
  TF1* spec_SD_7_Xe = (TF1*)spec_SD->Clone("spec_SD_7_Xe");
  TF1* spec_SD_7_Ge_Jun = (TF1*)spec_SD->Clone("spec_SD_7_Ge_Jun");
  TF1* spec_SD_7_Ge_Dec = (TF1*)spec_SD->Clone("spec_SD_7_Ge_Dec");
  spec_SD_7_F->SetParameters( const0, const1[0] ,const2[0] ,const3 , A, M_D, v_0, rho, sigma_SD_0, Ja);
  spec_SD_7_Ge->SetParameters( const0, const1[0] ,const2[0] ,const3 , 73, M_D, v_0, rho, sigma_SD_0, Ja_Ge);
  spec_SD_7_Xe->SetParameters( const0, const1[0] ,const2[0] ,const3 , 131, M_D, v_0, rho, sigma_SD_0, Ja_Xe);
  spec_SD_7_Ge_Jun->SetParameters( const0, const1[6] ,const2[6] ,const3 , 73, M_D, v_0, rho, sigma_SD_0, Ja_Ge);
  spec_SD_7_Ge_Dec->SetParameters( const0, const1[12] ,const2[12] ,const3 , 73, M_D, v_0, rho, sigma_SD_0, Ja_Ge);
  */

  /*
  cout<<"hoge"<<endl;
  // quenching
  Double_t fac_que[40], keV_que[40], spec_que[40];
  Double_t integ50,integ100;
  ifstream data("quenching.dat");
  Int_t index=0;
  while(data >> fac_que[index] >> keV_que[index]){ index++; }
  data.close();
  for(Int_t i=0; i<40; i++){
    spec_que[i]=(spec_SD_100_F->Eval(fac_que[i]));
    //cout<<i<<" "<<fac_que[i]<<" "<<spec_que[i]<<endl;
    if(spec_que[i]<0)spec_que[i]=0;
    if(keV_que[i]>50) integ50+=spec_que[i];
    if(keV_que[i]>100)integ100+=spec_que[i];
  }
  TGraph *graph = new TGraph(40,keV_que,spec_que);
  graph->SetLineColor(4);
  graph->SetLineWidth(3);
  cout<<"50keV:"<<graph->Eval(50)<<" 100keV:"<<graph->Eval(100)<<endl;;
  cout<<"50-200keV:"<<integ50<<" 100-200keV:"<<integ100<<endl;
  */

  // waku & draw
  TH1F *waku51 = new TH1F("waku51","SI, ^{19}F,  #sigma=1e-6 [pb],  average",1,0,200);
  waku51->SetMinimum(0.000001);
  waku51->SetMaximum(0.003);
  waku51->GetXaxis()->SetTitle("recoil energy [keV]");
  waku51->GetYaxis()->SetTitle("event rate [counts/keV/kg/days]");
  waku51->GetYaxis()->SetTitleOffset(1.3);
  TLegend *legend51 = new TLegend(0.6,0.6,0.88,0.88);
  //  legend51->AddEntry(spec_SI_50_F, "^{19}F  M_{D}=50GeV","l");
  legend51->AddEntry(spec_SI_100_F,"^{19}F  M_{D}=100GeV","l");
  //  legend51->AddEntry(spec_SI_200_F,"^{19}F  M_{D}=200GeV","l");
  pad51->cd();
  waku51->Draw("");
  legend51->Draw("same");
  //  spec_SI_50_F->SetLineColor(2);
  //  spec_SI_50_F->Draw("same");
  spec_SI_100_F->SetLineColor(3);
  spec_SI_100_F->Draw("same");
  //  spec_SI_200_F->SetLineColor(4);
  //  spec_SI_200_F->Draw("same");

  TH1F *waku52 = new TH1F("waku52","SD, ^{19}F,  #sigma=1 [pb],  average",1,0,200);
  waku52->SetMinimum(0.001);
  waku52->SetMaximum(10);
  waku52->GetXaxis()->SetTitle("recoil energy [keV]");
  waku52->GetYaxis()->SetTitle("event rate [counts/keV/kg/days]");
  waku52->GetYaxis()->SetTitleOffset(1.3);
  TLegend *legend52 = new TLegend(0.6,0.6,0.88,0.88);
  //  legend52->AddEntry(spec_SD_50_F, "^{19}F  M_{D}=50GeV","l");
  legend52->AddEntry(spec_SD_100_F,"^{19}F  M_{D}=100GeV","l");
  //  legend52->AddEntry(spec_SD_200_F,"^{19}F  M_{D}=200GeV","l");
  pad52->cd();
  waku52->Draw("");
  legend52->Draw("same");
  //graph->Draw("PC");//expected energy spectrum with quanching
  spec_SD_50_F->SetLineColor(2);
  spec_SD_50_F->Draw("same");
  spec_SD_100_F->SetLineColor(3);
  spec_SD_100_F->Draw("same");
  spec_SD_200_F->SetLineColor(4);
  spec_SD_200_F->Draw("same");

  //######################################
  // for D thesis
  //######################################
  TCanvas *c_spec_SD=new TCanvas("c_spec_SD","c_spec_SD",800,800);
  c_spec_SD->SetLogy(1);
  TH1F *frame_spec_SD = new TH1F("frame_spec_SD","",1,0,200);
  frame_spec_SD->SetMinimum(1e-3);
  frame_spec_SD->SetMaximum(10);
  frame_spec_SD->GetXaxis()->SetNdivisions(505);
  frame_spec_SD->GetXaxis()->SetTitle("Recoil energy (keV)");
  frame_spec_SD->GetYaxis()->SetTitle("Event rate (counts/keV/kg/days)");
  frame_spec_SD->GetXaxis()->CenterTitle();
  frame_spec_SD->GetYaxis()->CenterTitle();
  frame_spec_SD->GetXaxis()->SetTitleFont(43);
  frame_spec_SD->GetXaxis()->SetTitleSize(40);
  frame_spec_SD->GetXaxis()->SetLabelFont(43);
  frame_spec_SD->GetXaxis()->SetLabelSize(40);
  frame_spec_SD->GetYaxis()->SetTitleFont(43);
  frame_spec_SD->GetYaxis()->SetTitleSize(40);
  frame_spec_SD->GetYaxis()->SetTitleOffset(1.5);
  frame_spec_SD->GetYaxis()->SetLabelFont(43);
  frame_spec_SD->GetYaxis()->SetLabelSize(40);
  frame_spec_SD->Draw();
  spec_SD_50_F->SetNpx(1000);
  spec_SD_50_F->SetLineColor(kViolet+1);
  spec_SD_50_F->SetLineWidth(3);
  spec_SD_50_F->Draw("same");
  spec_SD_100_F->SetNpx(1000);
  spec_SD_100_F->SetLineColor(kRed);
  spec_SD_100_F->SetLineWidth(3);
  spec_SD_100_F->Draw("same");
  spec_SD_200_F->SetNpx(1000);
  spec_SD_200_F->SetLineWidth(3);
  spec_SD_200_F->SetLineColor(kBlue);
  spec_SD_200_F->Draw("same");

  TLegend *leg_spec_SD=new TLegend(0.5,0.6,0.8,0.8);
  leg_spec_SD->SetBorderSize(0);
  leg_spec_SD->SetTextFont(43);
  leg_spec_SD->SetTextSize(30);
  leg_spec_SD->AddEntry(spec_SD_50_F,"m_{#chi}=50 GeV","l");
  leg_spec_SD->AddEntry(spec_SD_100_F,"m_{#chi}=100 GeV","l");
  leg_spec_SD->AddEntry(spec_SD_200_F,"m_{#chi}=200 GeV","l");
  leg_spec_SD->Draw();
  //c_spec_SD->Print("c_spec_SD.eps");



  /*
  TH1F *waku53 = new TH1F("waku53","SI, M_{D}=7 [GeV],  #sigma=1e-6 [pb],  average",1,0,20);
  waku53->SetMinimum(0.00001);
  waku53->SetMaximum(1);
  waku53->GetXaxis()->SetTitle("recoil energy [keV]");
  waku53->GetYaxis()->SetTitle("event rate [counts/keV/kg/days]");
  waku53->GetYaxis()->SetTitleOffset(1.3);
  TLegend *legend53 = new TLegend(0.6,0.6,0.88,0.88);
  legend53->AddEntry(spec_SI_7_F, "^{19}F  M_{D}=7GeV","l");
  legend53->AddEntry(spec_SI_7_Ge,"^{73}Ge  M_{D}=7GeV","l");
  legend53->AddEntry(spec_SI_7_Xe,"^{131}Xe  M_{D}=7GeV","l");
  pad53->cd();
  waku53->Draw("");
  legend53->Draw("same");
  spec_SI_7_F->SetLineColor(2);
  spec_SI_7_F->Draw("same");
  spec_SI_7_Ge->SetLineColor(3);
  spec_SI_7_Ge->Draw("same");
  spec_SI_7_Xe->SetLineColor(4);
  spec_SI_7_Xe->Draw("same");

  TH1F *waku54 = new TH1F("waku54","SD, M_{D}=7 [GeV],  #sigma=1 [pb],  average",1,0,20);
  waku54->SetMinimum(1);
  waku54->SetMaximum(10000);
  waku54->GetXaxis()->SetTitle("recoil energy [keV]");
  waku54->GetYaxis()->SetTitle("event rate [counts/keV/kg/days]");
  waku54->GetYaxis()->SetTitleOffset(1.3);
  TLegend *legend54 = new TLegend(0.6,0.6,0.88,0.88);
  legend54->AddEntry(spec_SD_7_F, "^{19}F  M_{D}=7GeV","l");
  legend54->AddEntry(spec_SD_7_Ge,"^{73}Ge  M_{D}=7GeV","l");
  legend54->AddEntry(spec_SD_7_Xe,"^{131}Xe  M_{D}=7GeV","l");
  pad54->cd();
  waku54->Draw("");
  legend54->Draw("same");
  spec_SD_7_F->SetLineColor(2);
  spec_SD_7_F->Draw("same");
  spec_SD_7_Ge->SetLineColor(3);
  spec_SD_7_Ge->Draw("same");
  spec_SD_7_Xe->SetLineColor(4);
  spec_SD_7_Xe->Draw("same");
  */
  TH1F *waku55 = new TH1F("waku55","SI, M_{D}=100 [GeV],  #sigma=1e-6 [pb]",1,0,125);
  waku55->SetMinimum(0.00005);
  waku55->SetMaximum(0.04);
  waku55->GetXaxis()->SetTitle("recoil energy [keV]");
  waku55->GetYaxis()->SetTitle("event rate [counts/keV/kg/days]");
  waku55->GetXaxis()->SetTitleSize(0.04);
  waku55->GetYaxis()->SetTitleSize(0.04);
  waku55->GetXaxis()->SetLabelSize(0.04);
  waku55->GetYaxis()->SetLabelSize(0.04);
  waku55->GetYaxis()->SetTitleOffset(1.3);
  TLegend *legend55 = new TLegend(0.7,0.5,0.88,0.88);
  legend55->AddEntry(spec_SI_100_F, "^{19}F","l");
  legend55->AddEntry(spec_SI_100_Na,"^{22}Na","l");
  legend55->AddEntry(spec_SI_100_Ar,"^{39}Ar","l");
  legend55->AddEntry(spec_SI_100_Ge,"^{74}Ge","l");
  legend55->AddEntry(spec_SI_100_I, "^{127}I","l");
  legend55->AddEntry(spec_SI_100_Xe,"^{132}Xe","l");
  pad55->cd();
  waku55->Draw("");
  legend55->Draw("same");
  spec_SI_100_F->SetLineColor(2);
  spec_SI_100_F->Draw("same");
  spec_SI_100_Na->SetLineColor(4);
  spec_SI_100_Na->Draw("same");
  spec_SI_100_Ar->SetLineColor(1);
  spec_SI_100_Ar->Draw("same");
  spec_SI_100_Ge->SetLineColor(6);
  spec_SI_100_Ge->Draw("same");
  spec_SI_100_I->SetLineColor(3);
  spec_SI_100_I->Draw("same");
  spec_SI_100_Xe->SetLineColor(5);
  spec_SI_100_Xe->Draw("same");

  TH1F *waku56 = new TH1F("waku56","SD, M_{D}=100 [GeV],  #sigma=1 [pb]",1,0,125);
  waku56->SetMinimum(0.001);
  waku56->SetMaximum(4);
  waku56->GetXaxis()->SetTitle("recoil energy [keV]");
  waku56->GetYaxis()->SetTitle("event rate [counts/keV/kg/days]");
  waku56->GetXaxis()->SetTitleSize(0.04);
  waku56->GetYaxis()->SetTitleSize(0.04);
  waku56->GetXaxis()->SetLabelSize(0.04);
  waku56->GetYaxis()->SetLabelSize(0.04);
  waku56->GetYaxis()->SetTitleOffset(1.3);
  TLegend *legend56 = new TLegend(0.7,0.5,0.88,0.88);
  legend56->AddEntry(spec_SD_100_F, "^{19}F","l");
  legend56->AddEntry(spec_SD_100_Li,"^{7}Li","l");
  legend56->AddEntry(spec_SD_100_Na,"^{23}Na","l");
  legend56->AddEntry(spec_SD_100_Ge,"^{73}Ge","l");
  legend56->AddEntry(spec_SD_100_I, "^{127}I","l");
  legend56->AddEntry(spec_SD_100_Xe,"^{131}Xe","l");
  pad56->cd();
  waku56->Draw("");
  legend56->Draw("same");
  spec_SD_100_F->SetLineColor(2);
  spec_SD_100_F->Draw("same");
  spec_SD_100_Li->SetLineColor(1);
  spec_SD_100_Li->Draw("same");
  spec_SD_100_Na->SetLineColor(4);
  spec_SD_100_Na->Draw("same");
  spec_SD_100_Ge->SetLineColor(6);
  spec_SD_100_Ge->Draw("same");
  spec_SD_100_I->SetLineColor(3);
  spec_SD_100_I->Draw("same");
  spec_SD_100_Xe->SetLineColor(5);
  spec_SD_100_Xe->Draw("same");


  TCanvas *c51;
  c51 = new TCanvas("c51", "c51", 125, 125, 500, 500);
  c51->SetFillColor(10); c51->SetGrid(1); c51->SetLogy(1); c51->SetLeftMargin(0.15); c51->SetRightMargin(0.05);
  waku55->Draw("");
  legend55->Draw("same");
  spec_SI_100_F->SetLineColor(2);
  spec_SI_100_F->Draw("same");
  spec_SI_100_Na->SetLineColor(4);
  spec_SI_100_Na->Draw("same");
  spec_SI_100_Ar->SetLineColor(1);
  spec_SI_100_Ar->Draw("same");
  spec_SI_100_Ge->SetLineColor(6);
  spec_SI_100_Ge->Draw("same");
  spec_SI_100_I->SetLineColor(3);
  spec_SI_100_I->Draw("same");
  spec_SI_100_Xe->SetLineColor(5);
  spec_SI_100_Xe->Draw("same");

  TCanvas *c52;
  c52 = new TCanvas("c52", "c52", 125, 125, 500, 500);
  c52->SetFillColor(10); c52->SetGrid(1); c52->SetLogy(1); c52->SetLeftMargin(0.15); c52->SetRightMargin(0.05);
  waku56->Draw("");
  legend56->Draw("same");
  spec_SD_100_F->SetLineColor(2);
  spec_SD_100_F->Draw("same");
  spec_SD_100_Li->SetLineColor(1);
  spec_SD_100_Li->Draw("same");
  spec_SD_100_Na->SetLineColor(4);
  spec_SD_100_Na->Draw("same");
  spec_SD_100_Ge->SetLineColor(6);
  spec_SD_100_Ge->Draw("same");
  spec_SD_100_I->SetLineColor(3);
  spec_SD_100_I->Draw("same");
  spec_SD_100_Xe->SetLineColor(5);
  spec_SD_100_Xe->Draw("same");

  // c51->Print("energy_spectrum_SI.eps");
  // c52->Print("energy_spectrum_SD.eps");



















  //**********************************************************************
  //** Annual Modulation *************************************************
  //**********************************************************************
  TCanvas *c6;
  c6 = new TCanvas("c6", "c6", 150, 150, 1000, 600);
  c6->SetFillColor(10);

  // pad
  TPad *pad61 = new TPad("pad61","pad61",0.0,0.5,0.33,1.0);
  pad61->SetFillColor(10);
  pad61->SetGridx(1);
  pad61->SetGridy(1);
  pad61->SetLogx(0);
  pad61->SetLogy(0);
  pad61->Draw();

  TPad *pad62 = new TPad("pad62","pad62",0.0,0.0,0.33,0.5);
  pad62->SetFillColor(10);
  pad62->SetGridx(1);
  pad62->SetGridy(1);
  pad62->SetLogx(0);
  pad62->SetLogy(0);
  pad62->Draw();

  TPad *pad63 = new TPad("pad63","pad63",0.33,0.5,0.66,1.0);
  pad63->SetFillColor(10);
  pad63->SetGridx(1);
  pad63->SetGridy(1);
  pad63->SetLogx(0);
  pad63->SetLogy(0);
  pad63->Draw();

  TPad *pad64 = new TPad("pad64","pad64",0.33,0.0,0.66,0.5);
  pad64->SetFillColor(10);
  pad64->SetGridx(1);
  pad64->SetGridy(1);
  pad64->SetLogx(0);
  pad64->SetLogy(0);
  pad64->Draw();

  TPad *pad65 = new TPad("pad65","pad65",0.66,0.5,1.0,1.0);
  pad65->SetFillColor(10);
  pad65->SetGridx(1);
  pad65->SetGridy(1);
  pad65->SetLogx(0);
  pad65->SetLogy(0);
  pad65->Draw();

  TPad *pad66 = new TPad("pad66","pad66",0.66,0.0,1.0,0.5);
  pad66->SetFillColor(10);
  pad66->SetGridx(1);
  pad66->SetGridy(1);
  pad66->SetLogx(0);
  pad66->SetLogy(0);
  pad66->Draw();

  // Jun Dec difference
  /*
  TF1* spec_SI_7_Ge_dif = new TF1("spec_SI_7_Ge_dif","(spec_SI_7_Ge_Jun)-(spec_SI_7_Ge_Dec)",0,5);
  TF1* spec_SD_7_Ge_dif = new TF1("spec_SD_7_Ge_dif","(spec_SD_7_Ge_Jun)-(spec_SD_7_Ge_Dec)",0,5);
  */
  TF1* spec_SD_100_F_dif = new TF1("spec_SD_100_F_dif","(spec_SD_100_F_Jun)-(spec_SD_100_F_Dec)",0,5);


  // waku & draw
  TH1F *waku61 = new TH1F("waku61","SI, M_{D}=7 [GeV],  #sigma=1e-6 [pb],  annual",1,0,5);
  waku61->SetMinimum(0.001);
  waku61->SetMaximum(0.3);
  waku61->GetXaxis()->SetTitle("recoil energy [keV]");
  waku61->GetYaxis()->SetTitle("rate [counts/keV/kg/days]");
  waku61->GetYaxis()->SetTitleOffset(1.3);
  pad61->cd();
  waku61->Draw("");
  //  spec_SI_7_Ge->Draw("same");
  //  spec_SI_7_Ge_Jun->SetLineColor(2);
  //  spec_SI_7_Ge_Jun->Draw("same");
  //  spec_SI_7_Ge_Dec->SetLineColor(4);
  //  spec_SI_7_Ge_Dec->Draw("same");

  /*
  TH1F *waku62 = new TH1F("waku62","SD, M_{D}=7 [GeV],  #sigma=1 [pb],  annual",1,0,5);
  waku62->SetMinimum(0.001);
  waku62->SetMaximum(30);
  waku62->GetXaxis()->SetTitle("recoil energy [keV]");
  waku62->GetYaxis()->SetTitle("rate [counts/keV/kg/days]");
  waku62->GetYaxis()->SetTitleOffset(1.3);
  pad62->cd();
  waku62->Draw("");
  spec_SD_7_Ge->Draw("same");
  spec_SD_7_Ge_Jun->SetLineColor(2);
  spec_SD_7_Ge_Jun->Draw("same");
  spec_SD_7_Ge_Dec->SetLineColor(4);
  spec_SD_7_Ge_Dec->Draw("same");
  */

  /*
  TH1F *waku63 = new TH1F("waku63","SI, ^{73}Ge, M_{D}=7 [GeV],  #sigma=1e-6 [pb],  Jun-Dec",1,0,5);
  waku63->SetMinimum(0);
  waku63->SetMaximum(0.01);
  waku63->GetXaxis()->SetTitle("recoil energy [keV]");
  waku63->GetYaxis()->SetTitle("rate [counts/keV/kg/days]");
  waku63->GetYaxis()->SetTitleOffset(1.3);
  pad63->cd();
  waku63->Draw("");
  spec_SI_7_Ge_dif->SetLineColor(3);
  spec_SI_7_Ge_dif->Draw("same");
  */
  /*
  TH1F *waku64 = new TH1F("waku64","SD, ^{73}Ge, M_{D}=7 [GeV],  #sigma=1 [pb],  Jun-Dec",1,0,5);
  waku64->SetMinimum(0);
  waku64->SetMaximum(1);
  waku64->GetXaxis()->SetTitle("recoil energy [keV]");
  waku64->GetYaxis()->SetTitle("rate [counts/keV/kg/days]");
  waku64->GetYaxis()->SetTitleOffset(1.3);
  pad64->cd();
  waku64->Draw("");
  spec_SD_7_Ge_dif->SetLineColor(3);
  spec_SD_7_Ge_dif->Draw("same");
  */

  TH1F *waku65 = new TH1F("waku65","SD, ^{19}F, M_{D}=100 [GeV],  #sigma=1 [pb]",1,0,200);
  waku65->SetMinimum(0);
  waku65->SetMaximum(2);
  waku65->GetXaxis()->SetTitle("recoil energy [keV]");
  waku65->GetYaxis()->SetTitle("rate [counts/keV/kg/days]");
  waku65->GetYaxis()->SetTitleOffset(1.3);
  pad65->cd();
  waku65->Draw("");
  //  spec_SD_100_F->Draw("same");
  spec_SD_100_F_Jun->SetLineColor(2);
  spec_SD_100_F_Jun->Draw("same");
  spec_SD_100_F_Dec->SetLineColor(4);
  spec_SD_100_F_Dec->Draw("same");

  //######################################
  // for D thesis
  //######################################
  TCanvas *c_spec_SD_an=new TCanvas("c_spec_SD_an","c_spec_SD_an",800,800);
  c_spec_SD_an->SetLogy(1);
  TH1F *frame_spec_SD_an = new TH1F("frame_spec_SD_an","",1,0,200);
  frame_spec_SD_an->SetMinimum(1e-3);
  frame_spec_SD_an->SetMaximum(10);
  frame_spec_SD_an->GetXaxis()->SetNdivisions(505);
  frame_spec_SD_an->GetXaxis()->SetTitle("Recoil energy (keV)");
  frame_spec_SD_an->GetYaxis()->SetTitle("Event rate (counts/keV/kg/days)");
  frame_spec_SD_an->GetXaxis()->CenterTitle();
  frame_spec_SD_an->GetYaxis()->CenterTitle();
  frame_spec_SD_an->GetXaxis()->SetTitleFont(43);
  frame_spec_SD_an->GetXaxis()->SetTitleSize(40);
  frame_spec_SD_an->GetXaxis()->SetLabelFont(43);
  frame_spec_SD_an->GetXaxis()->SetLabelSize(40);
  frame_spec_SD_an->GetYaxis()->SetTitleFont(43);
  frame_spec_SD_an->GetYaxis()->SetTitleSize(40);
  frame_spec_SD_an->GetYaxis()->SetTitleOffset(1.5);
  frame_spec_SD_an->GetYaxis()->SetLabelFont(43);
  frame_spec_SD_an->GetYaxis()->SetLabelSize(40);
  frame_spec_SD_an->Draw();
  spec_SD_100_F_Jun->SetNpx(1000);
  spec_SD_100_F_Jun->SetLineColor(kRed);
  spec_SD_100_F_Jun->SetLineWidth(3);
  spec_SD_100_F_Jun->Draw("same");
  spec_SD_100_F_Dec->SetNpx(1000);
  spec_SD_100_F_Dec->SetLineColor(kBlue);
  spec_SD_100_F_Dec->SetLineWidth(3);
  spec_SD_100_F_Dec->Draw("same");


  TLegend *leg_spec_SD_an=new TLegend(0.55,0.6,0.8,0.8);
  leg_spec_SD_an->SetBorderSize(0);
  leg_spec_SD_an->SetTextFont(43);
  leg_spec_SD_an->SetTextSize(30);
  leg_spec_SD_an->AddEntry(spec_SD_100_F_Jun,"June","l");
  leg_spec_SD_an->AddEntry(spec_SD_100_F_Dec,"December","l");
  leg_spec_SD_an->Draw();
  //c_spec_SD_an->Print("c_spec_SD_annual.eps");




  /**********************************************************************/
  /** Cygnus direction **************************************************/
  /**********************************************************************/
  TCanvas *c7;
  c7 = new TCanvas("c7", "c7", 175, 175, 700, 600);
  c7->SetFillColor(10);
  c7->SetGridx(0);
  c7->SetGridy(1);
  c7->SetLogx(0);
  c7->SetLogy(0);

  double cygbin = 90;

  TH2F *cygnus = new TH2F("cygnus","",cygbin,-90,90,cygbin,-90,90);
  cygnus->GetXaxis()->SetTitle("Azimuth [degree]");
  cygnus->GetYaxis()->SetTitle("Altitude [degree]");

  TH2F *axis= new TH2F("axis","axis",1800,-90,90,1800,-90,90);
  for(int i=-900;i<900;i++){
    for(int j=-90;j<=90;j+=30){
      axis->Fill(j*cos(PI*i/1800),i/10);
    }
  }
  // cyg_phi=20.62h0m?
  double cyg_phi_ji  = 20.62;
  double cyg_phi_fun = 0;
  // cyg_theta=42.03deg
  double cyg_theta_deg   = 42.03;
  // kam_phi=137deg18min
  double kam_phi_deg = 137;
  double kam_phi_min = 18;
  // kam_theta=36deg25min
  double kam_theta_deg   = 36;
  double kam_theta_min   = 25;

  double cyg_phi_seki = 2*PI*(cyg_phi_ji/24+cyg_phi_fun/24/60);
  double cyg_theta_seki   = 2*PI*(cyg_theta_deg/360);
  double kam_phi_seki = 2*PI*(kam_phi_deg/360+kam_phi_min/360/60);
  double kam_theta_seki   = 2*PI*(kam_theta_deg/360+kam_theta_min/360/60);
  cout<<"*  cyg_theta_seki = "<<cyg_theta_seki<<"  cyg_phi_seki = "<<cyg_phi_seki<<endl;
  cout<<"*  kam_theta_seki = "<<kam_theta_seki<<"  kam_phi_seki = "<<kam_phi_seki<<endl;

  double day,hour,date_phi_seki,hour_phi_seki,phi_seki,x_seki,y_seki,z_seki,x,y,z,theta,phi,phi2,bincon,phibin,thetabin;
  int i=0;
  while(1){ // line
    day           = 170;//7(Mar)+30(Apr)+31(May)+30(Jun)+31(Jul)+31(Aug)+10(Sep)
    hour          = 21-(double)i/100;
    date_phi_seki = 2*PI*(day/365);
    hour_phi_seki = 2*PI*(hour/24);
    phi_seki      = cyg_phi_seki - (/*kam_phi_seki +*/-PI + date_phi_seki + hour_phi_seki);
    //cout<<"*  date_phi_seki = "<<date_phi_seki<<"  hour_phi_seki = "<<hour_phi_seki<<"  phi_seki = "<<phi_seki<<endl;
    x_seki = cos(cyg_theta_seki)*cos(phi_seki);
    y_seki = cos(cyg_theta_seki)*sin(phi_seki);
    z_seki = sin(cyg_theta_seki);
    x =  x_seki*cos(PI/2-kam_theta_seki) - z_seki*sin(PI/2-kam_theta_seki);
    y =  y_seki;
    z =  z_seki*cos(PI/2-kam_theta_seki) + x_seki*sin(PI/2-kam_theta_seki);
    theta = 180/PI*asin(z/sqrt(x*x+y*y+z*z));
    phi   = 180/PI*asin(y/sqrt(x*x+y*y));
    phi2 = phi*fabs(cos(PI/180*theta));
    //cout<<"*  theta = "<<theta<<"  phi = "<<phi<<"  phi2 = "<<phi2<<endl;

    phibin = (int)((phi2+90)/(180/cygbin));
    thetabin = (int)((theta+90)/(180/cygbin));
    //cout<<"*  thetabin = "<<thetabin<<"  phibin = "<<phibin<<endl;
    bincon = cygnus->GetBinContent(phibin+1,thetabin+1);
    if( bincon == 0 ){
      cygnus->Fill(phi2,theta);
      //cout<<"bincon"<<endl;
    }else{
      //cout<<"else"<<endl;
    }
    if(i>2400)break;
    i += 1;
  }
  i=0;
  while(1){ // point
    day           = 170;//7(Mar)+30(Apr)+31(May)+30(Jun)+31(Jul)+31(Aug)+10(Sep)
    hour          = 21-(double)i*3;
    date_phi_seki = 2*PI*(day/365);
    hour_phi_seki = 2*PI*(hour/24);
    phi_seki      = cyg_phi_seki - (-PI + date_phi_seki + hour_phi_seki);
    x_seki = cos(cyg_theta_seki)*cos(phi_seki);
    y_seki = cos(cyg_theta_seki)*sin(phi_seki);
    z_seki = sin(cyg_theta_seki);
    x =  x_seki*cos(PI/2-kam_theta_seki) - z_seki*sin(PI/2-kam_theta_seki);
    y =  y_seki;
    z =  z_seki*cos(PI/2-kam_theta_seki) + x_seki*sin(PI/2-kam_theta_seki);
    theta = 180/PI*asin(z/sqrt(x*x+y*y+z*z));
    phi   = 180/PI*asin(y/sqrt(x*x+y*y));
    phi2 = phi*fabs(cos(PI/180*theta));

    for(int tmp=0;tmp<20;tmp++){
      cygnus->Fill(phi2,theta);
    }
    if(i>6)break;
    i += 1;
  }
  //cygnus->SetDraw("COLZ");
  cygnus->Draw("COLZ");
  axis->Draw("same");

  TText *t1 = new TText(10,85,"9/10 21:00");
  TText *t2 = new TText(47,50,"9/10 18:00");
  TText *t3 = new TText(53,15,"9/10 15:00");
  TText *t4 = new TText(30,-12,"9/10 12:00");
  TText *t5 = new TText(-15,-20,"9/10 9:00");
  TText *t6 = new TText(-62,-5,"9/10 6:00");
  TText *t7 = new TText(-77,30,"9/10 3:00");
  TText *t8 = new TText(-60,60,"9/10 0:00");
  t1->SetTextSize(0.03);
  t2->SetTextSize(0.03);
  t3->SetTextSize(0.03);
  t4->SetTextSize(0.03);
  t5->SetTextSize(0.03);
  t6->SetTextSize(0.03);
  t7->SetTextSize(0.03);
  t8->SetTextSize(0.03);
  t1->Draw("same");
  t2->Draw("same");
  t3->Draw("same");
  t4->Draw("same");
  t5->Draw("same");
  t6->Draw("same");
  t7->Draw("same");
  t8->Draw("same");

  TBox *box = new TBox(90.2,-100,110,100);
  box->SetFillColor(10);
  box->Draw("same");




  /**********************************************************************/
  /** Angle *************************************************************/
  /**********************************************************************/
  TCanvas *c8;
  c8 = new TCanvas("c8", "c8", 200, 200, 1100, 700);
  c8->SetFillColor(10);

  // pad
  TPad *pad81 = new TPad("pad81","pad81",0.0,0.5,0.33,1.0);
  pad81->SetFillColor(10);
  pad81->SetGridx(1);
  pad81->SetGridy(1);
  pad81->SetLogx(0);
  pad81->SetLogy(0);
  pad81->SetLeftMargin(0.12);
  pad81->SetRightMargin(0.15);
  pad81->Draw();

  TPad *pad82 = new TPad("pad82","pad82",0.0,0.0,0.33,0.5);
  pad82->SetFillColor(10);
  pad82->SetGridx(1);
  pad82->SetGridy(1);
  pad82->SetLogx(0);
  pad82->SetLogy(0);
  pad82->SetLeftMargin(0.12);
  pad82->SetRightMargin(0.15);
  pad82->Draw();

  TPad *pad83 = new TPad("pad83","pad83",0.33,0.5,0.66,1.0);
  pad83->SetFillColor(10);
  pad83->SetGridx(1);
  pad83->SetGridy(1);
  pad83->SetLogx(0);
  pad83->SetLogy(0);
  pad83->SetLeftMargin(0.15);
  pad83->SetRightMargin(0.06);
  pad83->Draw();

  TPad *pad84 = new TPad("pad84","pad84",0.33,0.0,0.66,0.5);
  pad84->SetFillColor(10);
  pad84->SetGridx(1);
  pad84->SetGridy(1);
  pad84->SetLogx(0);
  pad84->SetLogy(0);
  pad84->SetLeftMargin(0.15);
  pad84->SetRightMargin(0.06);
  pad84->Draw();

  // Recoil angle-energy distribution
  double Emax_draw_SI = 400;
  double Emax_draw_SD = 400;

  TF2 *angle_SI = new TF2("angle_SI","0.5*const5([1],[2],[3],[4])*exp(-([0]*x-sqrt(y/const4([1],[2],[3])))**2)*func_FSI_keV(y,[1])*[5]*func_cross_SI([2],[1])",-1,1,0,Emax_draw_SI);

  TF2 *angle_SD = new TF2("angle_SD","0.5*const5([1],[2],[3],[4])*exp(-([0]*x-sqrt(y/const4([1],[2],[3])))**2)*func_FSD_keV(y,[1])*[5]*func_cross_SD([2],[1],[6])",-1,1,0,Emax_draw_SD);

  // parameter
  /*
  rho        = 0.3;
  M_D        = 100; // [GeV/c^2]
  sigma_SI_0 = 0.000001; //[pb]
  sigma_SD_0 = 1; //[pb]
  */

  angle_SI->SetParameters(const2[0] , A, M_D, v_0, rho, sigma_SI_0);
  angle_SD->SetParameters(const2[0] , A, M_D, v_0, rho, sigma_SD_0, Ja);

  // bin
  int enebin=200;
  int angbin=200;
  angle_SI->SetNpx(angbin);
  angle_SI->SetNpy(enebin);
  angle_SD->SetNpx(angbin);
  angle_SD->SetNpy(enebin);

  // waku & draw
  pad81->cd();
  TH2F *waku81 = new TH2F("waku81","SI, F, M_{D}=100 [GeV],  #sigma=1e-6 [pb]",100,-1,1,100,0,Emax_draw_SI);
  waku81->GetXaxis()->SetTitle("cos#theta");
  waku81->GetYaxis()->SetTitle("recoil energy [keV]");
  waku81->GetYaxis()->SetTitleOffset(1.2);
  //  waku81->Draw("");
  angle_SI->Draw("COLZ");

  pad82->cd();
  TH2F *waku82 = new TH2F("waku82","SD, F, M_{D}=100 [GeV],  #sigma=1 [pb]",100,-1,1,100,0,Emax_draw_SD);
  waku82->GetXaxis()->SetTitle("cos#theta");
  waku82->GetYaxis()->SetTitle("recoil energy [keV]");
  waku82->GetYaxis()->SetTitleOffset(1.2);
  //  waku82->Draw("");
  angle_SD->Draw("COLZ");



  TH2D *cosin_SI_2D = (TH2D*)angle_SI->GetHistogram();
  TH1D *cosin_SI_01 = (TH1D*)cosin_SI_2D->ProjectionX("cosin_SI_01",25, 50); // 50 -100 keV
  TH1D *cosin_SI_02 = (TH1D*)cosin_SI_2D->ProjectionX("cosin_SI_02",50,100); // 100-200 keV
  TH1D *cosin_SI_03 = (TH1D*)cosin_SI_2D->ProjectionX("cosin_SI_03",50, 60); // 100-120 keV
  cosin_SI_01->Scale(0.02);
  cosin_SI_02->Scale(0.01);
  cosin_SI_03->Scale(0.05);

  TH2D *cosin_SD_2D = (TH2D*)angle_SD->GetHistogram();


  TH1D *cosin_SD_01 = (TH1D*)cosin_SD_2D->ProjectionX("cosin_SD_01",cosin_SD_2D->GetYaxis()->FindBin(50),cosin_SD_2D->GetYaxis()->FindBin(60)); // 50 -100 keV
  //  double bin_num=(double) ( cosin_SD_2D->GetXaxis()->FindBin(60)- cosin_SD_2D->GetXaxis()->FindBin(50)+1);
  cosin_SD_01->Scale(1/50.);


  // for D thesis
  TH1D *cosin_SD_03 = (TH1D*)cosin_SD_2D->ProjectionX("cosin_SD_03",cosin_SD_2D->GetYaxis()->FindBin(50.), cosin_SD_2D->GetYaxis()->FindBin(60.)); // 50-60
  //bin_num=(double) ( cosin_SD_2D->GetXaxis()->FindBin(120)- cosin_SD_2D->GetXaxis()->FindBin(100)+1);
  cosin_SD_03->Scale(1/10.);
  cosin_SD_03->Scale(10.);  // to counts/cos/kg/days

  TH1D *cosin_SD_03_abs=(TH1D*)cosin_SD_03->Clone("cosin_SD_03_abs");
  int bin_num=cosin_SD_03_abs->GetNbinsX();
  for(int i=0;i<bin_num;i++){
    double this_value=cosin_SD_03_abs->GetBinContent(i+1);
    double cos=cosin_SD_03_abs->GetBinCenter(i+1);
    if(cos<0) continue;
    double ops_cos=-cos;
    double ops_value=cosin_SD_03_abs->GetBinContent(cosin_SD_03->GetXaxis()->FindBin(ops_cos));
    cosin_SD_03_abs->SetBinContent(i+1,this_value+ops_value);
  }



  TH1D *cosin_SD_02 = (TH1D*)cosin_SD_2D->ProjectionX("cosin_SD_02",50,100); // 100-200 keV
  //  cosin_SD_02->Scale(0.01);




  pad83->cd();
  cosin_SI_01->SetTitle("SI, F, M_{D}=100 [GeV],  #sigma=1e-6 [pb]");
  cosin_SI_01->GetXaxis()->SetTitle("cos#theta");
  cosin_SI_01->GetYaxis()->SetTitle("rate [count/keV/kg/days/cos#theta]");
  cosin_SI_01->GetYaxis()->SetTitleOffset(1.4);
  cosin_SI_01->SetLineColor(2);
  cosin_SI_02->SetLineColor(3);
  cosin_SI_03->SetLineColor(4);
  cosin_SI_01->Draw("same");
  cosin_SI_03->Draw("same");
  cosin_SI_02->Draw("same");
  TLegend *legend83 = new TLegend(0.17,0.6,0.4,0.88);
  legend83->AddEntry(cosin_SI_01,"50-100keV","l");
  legend83->AddEntry(cosin_SI_02,"100-200keV","l");
  legend83->AddEntry(cosin_SI_03,"100-120keV","l");
  legend83->Draw("same");

  pad84->cd();
  cosin_SD_01->SetTitle("SD, F, M_{D}=100 [GeV],  #sigma=1 [pb]");
  cosin_SD_01->GetXaxis()->SetTitle("cos#theta");
  cosin_SD_01->GetYaxis()->SetTitle("rate [count/keV/kg/days/cos#theta]");
  cosin_SD_01->GetYaxis()->SetTitleOffset(1.9);
  cosin_SD_01->SetLineColor(2);
  cosin_SD_02->SetLineColor(3);
  cosin_SD_03->SetLineColor(4);
  //  cosin_SD_01->Draw("same");
  cosin_SD_03->Draw("same HIST");
  //  cosin_SD_02->Draw("same");
  TLegend *legend84 = new TLegend(0.17,0.6,0.4,0.88);
  legend84->AddEntry(cosin_SD_01,"50-100keV","l");
  legend84->AddEntry(cosin_SD_02,"100-200keV","l");
  legend84->AddEntry(cosin_SD_03,"100-120keV","l");
  legend84->Draw("same");



  //######################################
  // for D thesis
  //######################################
  TCanvas *c_cos_SD=new TCanvas("c_cos_SD","c_cos_SD",800,800);
  TH1F *frame_cos_SD = new TH1F("frame_cos_SD","",1,-1,1);
  frame_cos_SD->SetMinimum(0);
  frame_cos_SD->SetMaximum(5.0);
  frame_cos_SD->GetXaxis()->SetNdivisions(505);
  frame_cos_SD->GetYaxis()->SetNdivisions(505);
  frame_cos_SD->GetXaxis()->SetTitle("cos#theta");
  frame_cos_SD->GetYaxis()->SetTitle("Event rate (counts/cos#theta/kg/days)");
  frame_cos_SD->GetXaxis()->CenterTitle();
  frame_cos_SD->GetYaxis()->CenterTitle();
  frame_cos_SD->GetXaxis()->SetTitleFont(43);
  frame_cos_SD->GetXaxis()->SetTitleSize(40);
  frame_cos_SD->GetXaxis()->SetLabelFont(43);
  frame_cos_SD->GetXaxis()->SetLabelSize(40);
  frame_cos_SD->GetYaxis()->SetTitleFont(43);
  frame_cos_SD->GetYaxis()->SetTitleSize(40);
  frame_cos_SD->GetYaxis()->SetTitleOffset(1.5);
  frame_cos_SD->GetYaxis()->SetLabelFont(43);
  frame_cos_SD->GetYaxis()->SetLabelSize(40);
  frame_cos_SD->Draw();
  //  cosin_SD_01->Draw("same HIST");
  cosin_SD_03->SetFillColor(kWhite);
  cosin_SD_03->SetLineWidth(3);
  cosin_SD_03->SetLineColor(kBlue);

  cosin_SD_03->Rebin(8);
  cosin_SD_03->Scale(1/8.);
  cosin_SD_03->Draw("same HIST");

  double iso_value=cosin_SD_03->Integral()/cosin_SD_03->GetNbinsX();
  TLine *line=new TLine(-1,iso_value,1,iso_value);
  line->SetLineWidth(3);
  line->SetLineStyle(9);
  line->Draw();

  TLegend *leg_cos_SD=new TLegend(0.2,0.6,0.5,0.8);
  leg_cos_SD->SetBorderSize(0);
  leg_cos_SD->SetTextFont(43);
  leg_cos_SD->SetTextSize(30);
  leg_cos_SD->AddEntry(cosin_SD_03,"50 keV #leq E #leq 60 keV","l");
  leg_cos_SD->AddEntry(line,"isotropy","l");
  leg_cos_SD->Draw();
  //c_cos_SD->Print("F_cos_SD.eps");

  TCanvas *c_cos_SD_abs=new TCanvas("c_cos_SD_abs","c_cos_SD_abs",800,800);
  TH1F *frame_cos_SD_abs = new TH1F("frame_cos_SD_abs","",1,0,1);
  frame_cos_SD_abs->SetMinimum(0);
  frame_cos_SD_abs->SetMaximum(5);
  frame_cos_SD_abs->GetXaxis()->SetNdivisions(505);
  frame_cos_SD_abs->GetYaxis()->SetNdivisions(505);
  frame_cos_SD_abs->GetXaxis()->SetTitle("|cos#theta|");
  frame_cos_SD_abs->GetYaxis()->SetTitle("Event rate (counts/cos#theta/kg/days)");
  frame_cos_SD_abs->GetXaxis()->CenterTitle();
  frame_cos_SD_abs->GetYaxis()->CenterTitle();
  frame_cos_SD_abs->GetXaxis()->SetTitleFont(43);
  frame_cos_SD_abs->GetXaxis()->SetTitleSize(40);
  frame_cos_SD_abs->GetXaxis()->SetLabelFont(43);
  frame_cos_SD_abs->GetXaxis()->SetLabelSize(40);
  frame_cos_SD_abs->GetYaxis()->SetTitleFont(43);
  frame_cos_SD_abs->GetYaxis()->SetTitleSize(40);
  frame_cos_SD_abs->GetYaxis()->SetTitleOffset(1.5);
  frame_cos_SD_abs->GetYaxis()->SetLabelFont(43);
  frame_cos_SD_abs->GetYaxis()->SetLabelSize(40);
  frame_cos_SD_abs->Draw();

  cosin_SD_03_abs->SetFillColor(kWhite);
  cosin_SD_03_abs->SetLineWidth(3);
  cosin_SD_03_abs->SetLineColor(kBlue);

  // calc forward-backward ratio
  //  int tmp_bin_num=cosin_SD_03_abs->GetNbinsX();
  double backward=cosin_SD_03_abs->Integral(cosin_SD_03_abs->GetXaxis()->FindBin(0.0),
					   cosin_SD_03_abs->GetXaxis()->FindBin(0.5));
  double forward=cosin_SD_03_abs->Integral(cosin_SD_03_abs->GetXaxis()->FindBin(0.5),
					    cosin_SD_03_abs->GetXaxis()->FindBin(1.0));
  cerr << "forward-backward ratio of abs cos = " << forward/backward << endl;



  cosin_SD_03_abs->Rebin(4);
  cosin_SD_03_abs->Scale(1/4.);
  cosin_SD_03_abs->Draw("same HIST");

  double iso_value_abs=cosin_SD_03_abs->Integral(cosin_SD_03_abs->GetXaxis()->FindBin(0.),cosin_SD_03_abs->GetXaxis()->FindBin(1.))/(cosin_SD_03_abs->GetNbinsX()/2);
  TLine *line_abs=new TLine(0,iso_value_abs,1,iso_value_abs);
  line_abs->SetLineWidth(3);
  line_abs->SetLineStyle(9);
  line_abs->Draw();

  line_abs->Draw();

  TLegend *leg_cos_SD_abs=new TLegend(0.2,0.6,0.5,0.8);
  leg_cos_SD_abs->SetBorderSize(0);
  leg_cos_SD_abs->SetTextFont(43);
  leg_cos_SD_abs->SetTextSize(30);
  leg_cos_SD_abs->AddEntry(cosin_SD_03_abs,"50 keV #leq E #leq 60 keV","l");
  leg_cos_SD_abs->AddEntry(line,"isotropy","l");
  leg_cos_SD_abs->Draw();
  //c_cos_SD_abs->Print("F_cos_SD_abs.eps");

  TCanvas *c_2d=new TCanvas("c_2d","",800,800);
  cosin_SD_2D->SetContour(1000);
  gStyle->SetPalette(55);
  cosin_SD_2D->GetXaxis()->CenterTitle();
  cosin_SD_2D->GetYaxis()->CenterTitle();
  cosin_SD_2D->GetXaxis()->SetTitleFont(43);
  cosin_SD_2D->GetXaxis()->SetTitleSize(40);
  cosin_SD_2D->GetXaxis()->SetLabelFont(43);
  cosin_SD_2D->GetXaxis()->SetLabelSize(40);
  cosin_SD_2D->GetYaxis()->SetTitleFont(43);
  cosin_SD_2D->GetYaxis()->SetTitleSize(40);
  cosin_SD_2D->GetYaxis()->SetTitleOffset(1.5);
  cosin_SD_2D->GetYaxis()->SetLabelFont(43);
  cosin_SD_2D->GetYaxis()->SetLabelSize(40);
  cosin_SD_2D->SetTitle("");
  cosin_SD_2D->GetYaxis()->SetTitle("Recoil energy (keV)");
  cosin_SD_2D->GetXaxis()->SetTitle("cos#theta");
  cosin_SD_2D->GetYaxis()->SetNdivisions(505);
  cosin_SD_2D->GetXaxis()->SetNdivisions(505);
  cosin_SD_2D->GetZaxis()->CenterTitle();
  cosin_SD_2D->GetZaxis()->SetTitleOffset(1.3);
  cosin_SD_2D->GetZaxis()->SetTitle(" (counts/keV/kg/days/cos#theta)");
  cosin_SD_2D->Draw("colz");
  //c_2d->Print("event_rate_2D.eps");

  // save
  // TFile *file = new TFile("expected_signal.root","RECREATE");
  // cosin_SD_03->Write();
  // cosin_SD_03_abs->Write();
  // file->Close();



  TF1* spec_SD_300_F = (TF1*)spec_SD->Clone("spec_SD_300_F");
  spec_SD_300_F->SetParameters( const0, const1[0] ,const2[0] ,const3 , A, 300, v_0, rho, sigma_SD_0, Ja);
  spec_SD_300_F->SetNpx( 1000 );
  TF1* spec_SD_400_F = (TF1*)spec_SD->Clone("spec_SD_400_F");
  spec_SD_400_F->SetParameters( const0, const1[0] ,const2[0] ,const3 , A, 400, v_0, rho, sigma_SD_0, Ja);
  spec_SD_400_F->SetNpx( 1000 );
  TF1* spec_SD_25_F = (TF1*)spec_SD->Clone("spec_SD_25_F");
  spec_SD_25_F->SetParameters( const0, const1[0] ,const2[0] ,const3 , A, 25, v_0, rho, sigma_SD_0, Ja);
  spec_SD_25_F->SetNpx( 1000 );
  TF1* spec_SD_10_F = (TF1*)spec_SD->Clone("spec_SD_10_F");
  spec_SD_10_F->SetParameters( const0, const1[0] ,const2[0] ,const3 , A, 10, v_0, rho, sigma_SD_0, Ja);
  spec_SD_10_F->SetNpx( 1000 );
  TF1* spec_SD_5_F = (TF1*)spec_SD->Clone("spec_SD_5_F");
  spec_SD_5_F->SetParameters( const0, const1[0] ,const2[0] ,const3 , A, 5, v_0, rho, sigma_SD_0, Ja);
  spec_SD_5_F->SetNpx( 1000 );


  TFile enFile( "energy.root", "RECREATE" );

  TH1F hist_SD_5_F  ( "hist_SD_5_F",   "hist_SD_5_F",   200, 0, 200 );
  TH1F hist_SD_10_F ( "hist_SD_10_F",  "hist_SD_10_F",  200, 0, 200 );
  TH1F hist_SD_25_F ( "hist_SD_25_F",  "hist_SD_25_F",  200, 0, 200 );
  TH1F hist_SD_50_F ( "hist_SD_50_F",  "hist_SD_50_F",  200, 0, 200 );
  TH1F hist_SD_100_F( "hist_SD_100_F", "hist_SD_100_F", 200, 0, 200 );
  TH1F hist_SD_200_F( "hist_SD_200_F", "hist_SD_200_F", 200, 0, 200 );
  TH1F hist_SD_300_F( "hist_SD_300_F", "hist_SD_300_F", 200, 0, 200 );
  TH1F hist_SD_400_F( "hist_SD_400_F", "hist_SD_400_F", 200, 0, 200 );
  hist_SD_5_F.SetDirectory( &enFile );
  hist_SD_10_F.SetDirectory( &enFile );
  hist_SD_25_F.SetDirectory( &enFile );
  hist_SD_50_F.SetDirectory( &enFile );
  hist_SD_100_F.SetDirectory( &enFile );
  hist_SD_200_F.SetDirectory( &enFile );
  hist_SD_300_F.SetDirectory( &enFile );
  hist_SD_400_F.SetDirectory( &enFile );
  for( int i = 0; i < 10000000; ++i ) {
      hist_SD_5_F.Fill( spec_SD_5_F->GetRandom( 0, 200 ) );
      hist_SD_10_F.Fill( spec_SD_10_F->GetRandom( 0, 200 ) );
      hist_SD_25_F.Fill( spec_SD_25_F->GetRandom( 0, 200 ) );
      hist_SD_50_F.Fill( spec_SD_50_F->GetRandom( 0, 200 ) );
      hist_SD_100_F.Fill( spec_SD_100_F->GetRandom( 0, 200 ) );
      hist_SD_200_F.Fill( spec_SD_200_F->GetRandom( 0, 200 ) );
      hist_SD_300_F.Fill( spec_SD_200_F->GetRandom( 0, 200 ) );
      hist_SD_400_F.Fill( spec_SD_200_F->GetRandom( 0, 200 ) );
  }
  hist_SD_5_F.Write( );
  hist_SD_10_F.Write( );
  hist_SD_25_F.Write( );
  hist_SD_50_F.Write( );
  hist_SD_100_F.Write( );
  hist_SD_200_F.Write( );
  hist_SD_300_F.Write( );
  hist_SD_400_F.Write( );
  enFile.Close( );
}
