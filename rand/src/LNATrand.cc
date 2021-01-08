/*this code generates the energy density for LNAT profile.*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <string>
#include <sstream>

#include "func.h"

// My utility
#include "inc/shinclude.h"

int flg2; /*global variable*/


using namespace std;
int main(int argc,char **argv)
{

  if(argc != 4){
    std::cerr << "INPUT ERROR" << std::endl;
    std::cerr << "./main [DM mass (GeV)] [atom number] [ratio of anisotropic gauss]" << endl;
    return 0;
  }

  // dark matter mass
  double mdm=atof(argv[1]); // GeV

  // target atom
  int atom=atoi(argv[2]);
  if( atom !=10 && atom !=11 ){
    std::cerr << "INPUT ATOM ERROR" << std::endl;
    cerr << " 10=F, 11=Ag " << endl;
    return 0;
  }

  int r=atoi(argv[3]); //ratio of anisotropic gaussian


  //---------------------------------------//
  int snum=1e07; //number of signal 2720 511 //~
  //---------------------------------------//
  time_t t = time(NULL);
  std::cout <<  ctime(&t) << std::endl;
  std::cout << "atom is " << atom << " (0=C, 1=S, 2=Br, 3=I, 10=F, 11=Ag)" << std::endl;
  std::cout << "mdm=" << mdm <<  "signal=" << log10f(snum) << std::endl;
  std::cout <<  "r=" << r << endl;

  srand((unsigned) time(NULL));

  int isgnl,i1,i2,i3;
  double imax=4000.; // max x axis
  double *f3,*ff3,*g3, *h3,  *udmlist, *Ctheta2list;
  double *x; // x bin
  double vdm;
  double mnl[15]; // nuclear mass matrics
  double clight=299792.458; // km/s
  double delta1,Ctheta1,Ctheta2,Sphi1,Sphi2,un,qn,track,yzl,udm,tmp1, FF;

  FILE *fp1; // output file
  FILE *fp2; // output info file
  double CT,V,Udm,Un;
  double alpha,beta,gamma,Ct2E,Sp2E,vdmE0, boost;
  vec vdm0,vdmE,un0,unE;
  vec4 vtmp;
  int flgFF=0,vflg=0;//~修正

  double vdm1,vvdm1,vdm2,vdm3;

  x=(double*)malloc(100000*sizeof(double));
  f3=(double*)malloc(100000*sizeof(double));
  ff3=(double*)malloc(100000*sizeof(double));
  g3=(double*)malloc(100000*sizeof(double));
  h3=(double*)malloc(100000*sizeof(double));
  udmlist=(double*)malloc(200000*sizeof(double));
  Ctheta2list=(double*)malloc(200000*sizeof(double));


  //--------- Nucelar Mass Matrix -----------//
  //SI
  mnl[0]=0.932*(0.989*12.+0.011*13);		//C
  mnl[1]=0.932*(0.950*32.+0.008*33.+0.042*34.);	//S
  mnl[2]=0.932*(0.5069*79.+0.4931*81.);		//Br
  mnl[3]=0.932*127.;				//I
  mnl[11]=0.932*(0.5184*107.+0.48161*109.);	//Ag
  //SD
  mnl[10]=0.932*19.;				//F

  //---------------------------------------//
  // calculation of DM distribution
  //---------------------------------------//
  xlist(x,imax);
  f2list(f3,imax);//calculating invese funtion of LNATrand
  ff2list(ff3,imax);
  g2list(g3,imax);
  h2list(h3,imax);
	

  //---------------------------------------//
  // Create Output File
  //---------------------------------------//
  std::string filename_head="LNAT1e7_";
  std::string filename_ext=".dat";
  std::stringstream s_mdm;
  std::stringstream s_atom;
  std::stringstream s_r;
  s_mdm << mdm; s_atom << atom; s_r << r;
  std::string filename=filename_head + s_mdm.str() + "_" + s_atom.str() +"_" + s_r.str()+ filename_ext;
  std::string filename_ext_info=".info";
  std::string filename_info=filename_head + s_mdm.str() + "_" + s_atom.str() +"_" + s_r.str()+ filename_ext_info;
  fp1 = fopen(filename.c_str(), "wt");//~
  fp2 = fopen(filename_info.c_str(), "wt");//~



  //=================================//
  // ROOP START
  //=================================//
  double mn=mnl[atom];
  int j=0;
  int total_count=0;
  while( j < snum ){
    total_count++;

    if(j%1000000==0.) printf("event number = %d \r",j);

    for(vdmE0=100000.;vdmE0>650.;){
      vdm1=linear(f3,x,rand()/(RAND_MAX+1.0));
      vvdm1=linear(ff3,x,rand()/(RAND_MAX+1.0));
      vdm2=linear(g3,x,rand()/(RAND_MAX+1.0));
      vdm3=linear(h3,x,rand()/(RAND_MAX+1.0));
      vtmp=vdmE_rand(vdm1,vvdm1,vdm2,vdm3, r);
      vdmE.x=vtmp.x;
      vdmE.y=vtmp.y;
      vdmE.z=vtmp.z;
      if(vtmp.w > 3.) vflg=1;
      vdmE0=norm(vdmE);
    }
		
    V=mdm/(mdm+mn)*vdmE0;
    CT= randompm()*rand()/(RAND_MAX+1.0);
    Ctheta1=fCtheta1(mdm,mn,CT);
    Ctheta2=fCtheta2(CT);
    Udm=mn/(mdm+mn)*vdmE0;
    Un=mdm/mn*Udm;
    udm=fudm(V,Udm,CT,Ctheta1);
    un=fun(vdmE0,udm,mdm,mn);
    //Sphi1=0.;
    //Sphi2=-Sphi1;
    Sphi1=fSphi1();
    Sphi2=-Sphi1;
      
    tmp=randompm();
    un0.x=un*tmp*sqrt(1-sq(Ctheta2))*randompm()*sqrt(1-sq(Sphi2));
    un0.y=un*tmp*sqrt(1-sq(Ctheta2))*Sphi2;
    un0.z=un*Ctheta2;
    alpha=rotanglea(vdmE);
    beta=rotangleb(vdmE);
    unE=rotation(un0,alpha,beta);

    Ct2E=unE.z/norm(unE);
    Sp2E=randompm()*unE.y/norm(unE)/sqrt(1-sq(Ct2E));

    if((long double) (rand()/(RAND_MAX+1.0)) <= FF2(ER(mn,mdm,norm(unE)),atom)){

      fprintf(fp1,"%d %d %f %f %f %d\n",
	      r,atom,Ct2E,Sp2E,ER(mn,mdm,norm(unE)),vflg);
      j++;
      if(vflg!=0&&vflg!=1)
	printf("%d %d %f %f %f %d\n",
	       r,atom,Ct2E,Sp2E,ER(mn,mdm,norm(unE)),vflg);
    }
    //1   2   3    4      5
    vflg=0;
    flgFF=0;

  } //end of signal number loop
 


  // end 
  t=time(NULL);
  printf("%s", ctime(&t));
  fprintf(fp2,"%d\n",total_count);
  fclose(fp1);
  fclose(fp2);

  return 0;

}
