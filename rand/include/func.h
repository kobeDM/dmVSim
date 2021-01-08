typedef struct TAGvec {
    double x;
    double y;
    double z;
} vec;

typedef struct TAGvec4 {
  double x;
  double y;
  double z;
  double w;
} vec4;

double tmp;
double *vvdmz,*ffz;
double EE=2.7182818284;
double gamma_f=0.886227;

/*subroutines......................*/
/*subroutines......................*/

/*make square*/
double sq(double x){
  return pow(x,2.);
    }

/*......interpolate list f3[] by linear approx*/
double linear(double *f3,double *x, double f)
{
  double a,b;
  int j=0;
  while(f3[j]<f){
    j=j+1;
  };
  a=(x[j]-x[j-1])/(f3[j]-f3[j-1]);
  b=x[j]-a*f3[j];
   if(j==0){
    a=0;
    b=0;
  }
  return a*f+b;
}

/*.................................*/
double f1(double v)//vr
{
	double v0=240.4;
	return exp(-sq(v/v0))/(2*v0*gamma_f);
}

/*.................................*/
double f2(double x1,double x2)/*integrate from x1>x2*/
{
    int j;
    double fsum;
    double step=.5;
    fsum=0.;
    for(j=1;j*step<(x2-x1)+step;j++){
        fsum=fsum+(f1(x1+(j-1)*step)+f1(x1+j*step))/2.*step;
    };
    return fsum;
}

/*.................................*/
double ff1(double v)//vphi_1
{
    double v0=250.;
    return exp(-sq(v/v0))/(2*v0*gamma_f);
}

/*.................................*/
long double ff2(double x1,double x2)/*integrate from x1>x2*/
{
    int j;
    double fsum;
    double step=0.5;
    fsum=0.;
    for(j=1;j*step<(x2-x1)+step;j++){
        fsum=fsum+(ff1(x1+(j-1)*step)+ff1(x1+j*step))/2.*step;
    };
    return fsum;
}

/*.................................*/
double g1(double v)//vphi_2
{
    double v0=120.;
	double mu=150.;
	return exp(-sq(v-mu)/sq(v0))/(2*v0*gamma_f);
}

/*.................................*/
double g2(double x1,double x2)/*integrate from x1>x2*/
{
    int j;
    double gsum;
    double step=.5;
    gsum=0.;
    for(j=1;j*step<(x2-x1)+step;j++){
        gsum=gsum+(g1(x1+(j-1)*step)+g1(x1+j*step))/2.*step;
    };
    return gsum;
}

/*.................................*/
double h1(double v)//vz
{
    double v0=214.6;
	return exp(-sq(v)/sq(v0))/(2*v0*gamma_f);
}

/*.................................*/
double h2(double x1,double x2)/*integrate from x1>x2*/
{
    int j;
    double hsum;
    double step=.5;
    hsum=0.;
    for(j=1;j*step<(x2-x1)+step;j++){
        hsum=hsum+(h1(x1+(j-1)*step)+h1(x1+j*step))/2.*step;
    };
    return hsum;
}

/*.................................*/
int xlist(double x[],double xmax)
{
  int j;
  double step=.5;
  j=0;
  while(-1000.+j*step<xmax+step){
    x[j]=-1000.+step*j;
	j++;
  }
  return 0;
}

int f2list(double f3[],double xmax)
{
	int j;
	double step=.5;
	double s;
	j=0;
	s=f2(-1000.,1000.);
	while(-1000.+j*step<xmax+step){
	f3[j]=f2(-1000.,-1000.+step*j)/s;
	j++; 
	}
	return 0;
}

int ff2list(double ff3[],double xmax)
{
    int j;
    double step=.5;
    double s;
    j=0;
    s=ff2(-1000.,1000.);
    while(-1000.+j*step<xmax+step){
        ff3[j]=ff2(-1000.,-1000.+step*j)/s;
        j++; 
    }
    return 0;
}

int g2list(double g3[],double xmax)
{
	int j;
	double step=.5;
	double s;
	j=0;
	s=g2(-1000.,1000.);
	while(-1000.+j*step<xmax+step){
	g3[j]=g2(-1000.,-1000.+step*j)/s;
	j++; 
	}
	return 0;
}

int h2list(double h3[],double xmax)
{
	int j;
	double step=.5;
	double s;
	j=0;
	s=h2(-1000.,1000.);
	while(-1000.+j*step<xmax+step){
        h3[j]=h2(-1000.,-1000.+step*j)/s;
        //printf("%lf\n",h3[j]);
        j++; 
	}
	return 0;
}

//EOM**********************************************************************
//EOM**********************************************************************
//EOM**********************************************************************
/*make plus or minus 1 randamly*/
double randompm(){
  if (rand()%2==0){return 1.;}
  else{return -1.;}
}

double fSphi1(){
  //return sin((rand()%360+1.)*3.1415/180.);
	double pi=4.*atan(1);
	return sin((double)(rand()*(360 - 0 + 1.0)/(1.0 + RAND_MAX)*pi/180.));
}

double fCtheta1(double mdm,double mn,double CT){
  double ct1,x;

  x=mdm/mn;
  ct1=(1-sq(CT))/sq(CT+x)+1;
  ct1=1/sqrt(ct1);
  return ct1;//always positive since DM is heavy
}

double fudm(double V,double Udm,double CT,double Ct1){
  return (Udm*CT+V)/Ct1;
}

double fct2(double ct1,double mdm,double mn){
  double x=mdm/mn;
  double b,c,c2t2;
  b=-x*(1-sq(ct1));
  c=sq(x)*(1-sq(ct1))-sq(ct1);

  if(sq(b)-c<0){
    printf("root is negative...\n");
  }

  if(sq(-b+sqrt(sq(b)-c))>1){
  c2t2=-b-sqrt(sq(b)-c);
  //  printf("take - \n");
  return sqrt((c2t2+1)/2);}
  else{
      c2t2=-b+sqrt(sq(b)-c);
      //printf("take +\n");
      return sqrt((c2t2+1)/2);}
}

double fCtheta2(double CT){
  double clight=299792.458;//km/s
  double theta2,pi=4.*atan(1);

  theta2=(pi-acos(CT))/2;
  return cos(theta2);//pm is automatically fixed.
  }

double fun(double vdm,double udm,double mdm,double mn){
  return sqrt(mdm/mn*(sq(vdm)-sq(udm)));
}
/*double fun(double V,double Un,double CT,double Ct2){
  double clight=299792.458;//km/s

  return (-Un*CT+V)/Ct2;
}
*/

double ER(double mn, double mdm, double un){
    double clight=299792.458;//km/s
    //  printf("mdm=%e,un=%e,ER=%e\n",mdm,un,mdm*sq(un)/2);
    return pow(10,6)*mn*sq(un/clight)/2;}//keV

vec plus(vec a, vec b)
{
    vec res;
    res.x = a.x + b.x;
    res.y = a.y + b.y;
    res.z = a.z + b.z;
    return res;
}

double norm(vec a){
  return sqrt(sq(a.x)+sq(a.y)+sq(a.z));
}

double rotanglea(vec vdmE){
    //returns rotation angle alpha, which rotates vdmE to (0,0,Norm(vdmE))
    double alpha;
    double pi=4.*atan(1);
    
    alpha=-asin(vdmE.y/norm(vdmE));
    return alpha;
}

double rotangleb(vec vdmE){
    //returns rotation angle beta, which rotates vdmE to (0,0,Norm(vdmE))
    double beta;
    double pi=4.*atan(1);
    
    beta=atan2(vdmE.x,vdmE.z);
    return beta;
}

vec rotation(vec u,double alpha,double beta){//rotation from COM to Earth sys
    vec ans;
    double al=alpha,be=beta;

    ans.x = u.x*cos(be) + u.y*sin(al)*sin(be)+u.z*cos(al)*sin(be);
    ans.y =               u.y*cos(al)           - u.z*sin(al);
    ans.z = -u.x*sin(be) +u.y*sin(al)*cos(be) + u.z*cos(al)*cos(be);
    return ans;
}

vec4 vdmE_rand(double vdm1,double vvdm2, double vdm2,double vdm3, int r){//vdm1;vr vvdm1;vphi-1 vdm2;vphi-2 vdm3;vz /r corresponds to f (ratio of anisotropy)
    vec v1,ve;
    vec4 vans;
    int i,vflg=0;
    double normv,v1x,v1y,phi,ct;
    
    ve.x=0;
    ve.y=0;
    ve.z=232;

    v1.x=vdm1;  
    v1.y=vdm3;
	if((rand()%100+1)<(100-r)+1){//r=0のとき完全等方的
        v1.z=-vvdm2;
    } //1~100の乱数がr以上ならFig.5(b)の赤を
    else{v1.z=-vdm2;
        vflg=5;
        //printf("anisotropy\n");
    } //それ以外ならFig.5(b)の緑を使う

    if(norm(v1)>650.){
        v1.x=100000.;
        v1.y=100000.;
        v1.z=100000.;}
    else{
        //at the earth
		v1=plus(v1,ve);
      }
  vans.x = v1.x;
  vans.y = v1.y;
  vans.z = v1.z;
  vans.w = vflg*1.;
    return vans;
}


/*Form factor******************************/
/*Form factor******************************/
/*Form factor******************************/

long double FF2(double ER, int atom){//ER in keV unit
    double A, c, a, rn, s, qr, qs; //qr=q*rn
    double pi=4.*atan(1);
    double f,ftest=0.;

    if (atom == 0) A= 12.011;//C
    if (atom == 1) A= 32.06;//S
    if (atom == 2) A=79.904;//Br
    if (atom == 3) A=126.90 ;//I //form factorはSD, SIどちらを使うの？？？？？？？？？？？？？？？？
    if (atom == 10) A=18.998;//F
    if (atom == 11) A=107.8682;//Ag
  

    a=0.52;
    s=0.9;
    c=1.23*pow(A, 1./3.) - 0.60;
    rn=sqrt(c*c+7./3.*pi*pi*a*a-5.*s*s);//4.11
    qr = 6.92*pow(10,-3)*sqrt(A*ER)*rn;
    qs = sqrt(2*1.932*A*ER)*0.9/197.3;
    if(atom == 0||atom == 2||atom == 11){
    f=3.*(sin(qr)-qr*cos(qr))/(pow(qr,3.))*exp(-qs*qs/2.); //(4.3) of Lewin, for SI
    }else{
        if(qr<2.55||4.5<qr){f=sin(qr)/(qr);}
        else{f=sqrt(0.047);}
    }//(4.2-4.5) of Lewin, for SD
    //printf("FF=%.30lf\n",f*f);
    return f*f;}
/*end of subroutine*/
