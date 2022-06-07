#define v_c 299792.458 // [km/sec] light speed
#define M_p 0.932// [Gev/c^2] proton mass

//-------------------------------//
//------- some parameters -------//
//-------------------------------//
const double PI = TMath::Pi();
const double v_0     = 220; // typical DM velocity [km/sec]   230 for paper, 220 for NishimuraPhD
const double v_esc   = 650; // escape velocity from galaxy[km/sec]   600 for paper, 650 for NishimuraPhD
const double v_E[13] = { 244.0 , 233.4 , 240.0 , 247.4 , 253.7 , 257.2 , 257.4 , 
                     254.3 , 248.5 , 241.4 , 234.6 , 230.0 , 229.5 }; // Earth velocity [km/sec]

//----------------------------------//
//------ Dark Matter Constants -----//
//----------------------------------//
const double rho        = 0.3;		// [GeV/c^2/cm^3]
const double sigma_SI_0 = 0.000001;	// [pb]
const double sigma_SD_0 = 1;		// 100; //[pb]

//---------------------------------//
//------- Target Constants --------//
//---------------------------------//
// F
const double A  = 19;
const double Ja = 0.647*1.0;
// others
Double_t Ja_Ge = 0.065*0.078;
Double_t Ja_Xe = 0.124*0.264;
Double_t Ja_Xe2= 0.055*0.212;
Double_t Ja_Li = 0.244*0.925;
Double_t Ja_Na = 0.041*1.0;
Double_t Ja_I  = 0.007*1.0;
Double_t Ja_Si = 0.063*0.047;
Double_t Ja_W  = 0.003*0.143;

Double_t k_0     = pow( PI*v_0*v_0 , 1.5 );
Double_t k_1     = k_0 * ( TMath::Erf(v_esc/v_0) - 2./sqrt(PI) * v_esc/v_0 * exp(-v_esc*v_esc/v_0/v_0));


//---------------------------------//
//------ some function ------------//
//---------------------------------//
Double_t rate(Double_t M_D, Double_t A){
  return 4*M_D*A*M_p/(M_D+A*M_p)/(M_D+A*M_p);
}
Double_t E_0(Double_t M_D, Double_t v_0){
  return 0.5*M_D*v_0*v_0/v_c/v_c*1000000;
}
Double_t R_0(Double_t M_D, Double_t A, Double_t rho, Double_t v_0){
  return 361/A/M_p/M_D*(rho/0.3)*(v_0/220); // R_0[real] = R_0*sigma
}
Double_t const0(){
  return k_0/k_1;;
}
Double_t const4(Double_t A, Double_t M_D, Double_t v_0){
  return E_0(M_D,v_0)*rate(M_D,A);
}
Double_t const5(Double_t A, Double_t M_D, Double_t v_0, Double_t rho){
  return R_0(M_D,A,rho,v_0)/E_0(M_D,v_0)/rate(M_D,A);
}

//------- Form Factor Functions
Double_t s(void){ return 0.9; }
Double_t r_n_SI(Double_t A){ return 1.14*pow(A,(1.0/3.0)); }
Double_t r_n_SD(Double_t A){
  return 1.0*pow(A,1.0/3.0);
}
Double_t qE(Double_t A){
  return 6.92*0.001*sqrt(A);
}
Double_t func_FSI_qr(Double_t x){ return pow((3*(sin(x)-x*cos(x))/pow(x,3.)),2.); }
Double_t func_FSD_qr(Double_t x){
  return pow(sin(x)/x,2.0);
}
Double_t func_FSI_q (Double_t x, Double_t A){ return func_FSI_qr(x*r_n_SI(A)) *exp(-x*x*s()*s()); }
Double_t func_FSD_q (Double_t x, Double_t A){
  return func_FSD_qr(x*r_n_SD(A));
}
Double_t func_FSI_keV (Double_t x, Double_t A){ return func_FSI_q(sqrt(x)*qE(A), A); }
Double_t func_FSD_keV (Double_t x, Double_t A){
  return func_FSD_q(sqrt(x)*qE(A), A);
}

//---------- Cross Section Function
Double_t func_cross_SI(Double_t x, Double_t A){ return pow(A * A*(x+M_p)/(x+M_p*A), 2); }
Double_t func_cross_SD(Double_t x, Double_t A, Double_t Ja){
  return Ja/0.75*pow(A*(x+M_p)/(x+M_p*A),2);
}

Double_t func_cos_ene_SD(Double_t x, Double_t y, Double_t C2, Double_t A, Double_t M_D,
			 Double_t v_0, Double_t rho, Double_t sigma_SD_0, Double_t Ja){
  double val =  0.5*const5(A,M_D,v_0,rho)*exp(-pow(C2*x-sqrt(y/const4(A,M_D,v_0)),2.0))*func_FSD_keV(y,A)
    *sigma_SD_0*func_cross_SD(M_D,A,Ja);
  val*=76./88.; //for F of CF4
  return val;
}


Double_t func_diff_que(Double_t x){
  //return 0.01*96.2536*exp(-72.9137/(x+65.9335));//TRIM_F_in_CF4_1.0atm.txt
  return 0.01*97.0806*exp(-44.7858/(x+33.7443));//TRIM_F_in_CF4_0.1atm.txt
}

Double_t func_Eres(Double_t ene, Double_t Eres_keV, Double_t Eres_percent){
  return sqrt( pow(Eres_keV,2) +pow(ene*Eres_percent/100.0,2) );
}

/*
  Double_t func_cos_ene_SD_abs(Double_t x, Double_t y, Double_t C2, Double_t A, Double_t M_D,
  Double_t v_0, Double_t rho, Double_t sigma_SD_0, Double_t Ja){
  double val = func_cos_ene_SD(x,y,C2,A,M_D,v_0,rho,sigma_SD_0,Ja)
  +func_cos_ene_SD(-x,y,C2,A,M_D,v_0,rho,sigma_SD_0,Ja);
  if(x<0) val = 0;
  return val;
  }
*/

Double_t func_spec_nor(Double_t x, Double_t p0, Double_t p1, Double_t p2, Double_t p3){
  return p0*(p1*(TMath::Erf(sqrt(x)+p2)-TMath::Erf(sqrt(x)-p2))-p3);
}



// func
void processing(double i, double maxi){
  char chartime[32];
  system("date +%Y/%m/%d/%H:%M:%S> nowtime");
  ifstream timedata("nowtime");
  timedata>>chartime;
  cout<<" i = "<<i<<" doing. Upto "<<maxi<<" "<<chartime<<endl;
  system("rm nowtime");
  return;
}
