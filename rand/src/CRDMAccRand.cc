#include "inc/shinclude.h"
#include "TF3.h"
#include "TRandom3.h"
#include "TTree.h"
#include "Math/IntegratorOptions.h"

const double V_LIGHT = 299792.458; //km/s
// const double PI = TMath::Pi( );

double randomPM()
{
    return gRandom->Integer( 2 ) ? 1.0 : -1.0;
}

/*make square*/
double sq( const double& x )
{
    return pow( x, 2.0 );
}

double gamFact( const double& velo )
{
    return 1.0 / sqrt( 1 + sq( velo/V_LIGHT ) );
}

double getSinDmScatPhi( ){
    return sin( gRandom->Rndm( ) * 2.0 * PI );
}

double getDmCosScatTheta( const double& dmM, const double& nuM, const double& cosScatAngleCom ){
    double x = dmM / nuM;
    double dmCosScatTheta = ( 1.0 - sq( cosScatAngleCom ) ) / sq( cosScatAngleCom + x ) + 1.0;
    return 1.0 / sqrt( dmCosScatTheta );//always positive if DM is heavy
}

double getNuCosScatTheta( const double& cosScatThetaCom ){
    double theta2 = ( PI - acos(cosScatThetaCom) ) * 0.5;
    return cos( theta2 );//pm is automatically fixed.
}

TVector3 getDmVVec( const double& velocity,
                    const double& theta,
                    const double& phi )
{
    TVector3 vec;
    vec.SetX( -1.0 * velocity * cos(theta) * cos(phi) );
    vec.SetY( -1.0 * velocity * cos(theta) * sin(phi) );
    vec.SetZ( -1.0 * velocity * sin(theta)            );
    return vec;
}

double getDmVExp( const double& dmInitVCom,
                  const double& dmFinVCom,
                  const double& cosScatThetaCom,
                  const double& cosScatThetaTmpExp )
{
    return ( cosScatThetaCom * dmFinVCom * gamFact( dmFinVCom ) + dmInitVCom * gamFact( dmInitVCom ) ) / cosScatThetaTmpExp;
}


double getNuFinVExp( const double& dmInitVExp,
                     const double& dmFinVExp,
                     const double& dmM,
                     const double& nuM )
{
    return sqrt( dmM / nuM * ( sq( dmInitVExp ) - sq( dmFinVExp ) ) );
}

// recoil energy [keV]
double getRecoilEnergy( const double& nuM,
                        const double& nuVExp )
{
    return pow( 10.0, 6.0 ) * nuM * sq( nuVExp/V_LIGHT ) * 0.5;
}

double getRotAngleAlpha( TVector3 dmInitVExpVec ){
    //returns rotation angle alpha, which rotates vdmE to (0,0,Norm(vdmE))
    return -1.0 * asin( dmInitVExpVec.Y( ) / dmInitVExpVec.Mag( ) );
}

double getRotAngleBeta( TVector3 dmInitVExpVec ){
    //returns rotation angle beta, which rotates vdmE to (0,0,Norm(vdmE))
    return atan2( dmInitVExpVec.X( ), dmInitVExpVec.Z( ) );
}

//rotation from COM to Earth sys
TVector3 rotation( TVector3      u,
                   const double& al,
                   const double& be )
{
    TVector3 ans;
    ans.SetX(  u.X( ) * cos( be ) + u.Y( ) * sin( al ) * sin( be ) + u.Z( ) * cos( al ) * sin( be ) );
    ans.SetY(                       u.Y( ) * cos( al )             - u.Z( ) * sin( al ) );
    ans.SetZ( -u.X( ) * sin( be ) + u.Y( ) * sin( al ) * cos( be ) + u.Z( ) * cos( al ) * cos( be ) );
    return ans;
}

/*Form factor******************************/

//ER in keV unit
double getFormFactorSq( const double& ER, const int& atom )
{
    double formFactor = 0.0;

    double A = 0.0;
    if( atom == 0  ) A = 12.011;   //C
    if( atom == 1  ) A = 32.06;    //S
    if( atom == 2  ) A = 79.904;   //Br
    if( atom == 3  ) A = 126.90;   //I //both for SD, SI
    if( atom == 10 ) A = 18.998;   //F
    if( atom == 11 ) A = 107.8682; //Ag

    double a  = 0.52;
    double s  = 0.9;
    double c  = 1.23 * pow( A, 0.333333 ) - 0.60;
    double rn = sqrt(c*c+7./3.*PI*PI*a*a-5.*s*s);//4.11
    double qr = 6.92*pow(10,-3)*sqrt(A*ER)*rn;
    double qs = sqrt(2*1.932*A*ER)*0.9/197.3;
    if( atom == 0 || atom == 2 || atom == 11 ) { //(4.3) of Lewin Smith, for SI
        formFactor = 3.0 * ( sin( qr ) - qr * cos( qr ) ) / ( pow( qr, 3.0 ) ) * exp( -0.5 *sq( qs ) );
    }
    else { //(4.2-4.5) of Lewin, for SD
        if( qr < 2.55 || 4.5 < qr ) formFactor = sin( qr ) / qr;
        else                        formFactor = sqrt( 0.047 );
    }

    return sq( formFactor );
}

int main( int argc, char** argv )
{
    time_t t = time( nullptr );
    printf("%s", ctime(&t));//time to start

    // set TRandom3
    if( gRandom != nullptr && gRandom->GetName( ) != "Random3" ) {
        delete gRandom;
        gRandom == nullptr;
    }

    if( gRandom == nullptr )
        gRandom = new TRandom3( 0 );

    std::cout << "Select randomizer: " << gRandom->GetName( ) << std::endl;

    if( argc != 3 ) {
        std::cerr << "INPUT ERROR" << std::endl;
        std::cerr << "./CRDMAccRand [input filename] [output filename]" << std::endl;
        abort( );
    }
    String input  = argv[1];
    String output = argv[2];
    
    double sysRelativeV = 0.0;
    double dmFinVCom  = 0.0;
    double nuFinVCom  = 0.0;
    double dmInitVExp = 0.0;
    double dmFinVExp  = 0.0;
    double nuFinVExp     = 0.0;

    double cosScatThetaCom = 0.0;
    double dmSinScatTheta  = 0.0, nuSinScatTheta = 0.0;
    double dmCosScatTheta  = 0.0, nuCosScatTheta = 0.0;
    double dmSinScatPhi    = 0.0, nuSinScatPhi   = 0.0;
    double dmCosScatPhi    = 0.0, nuCosScatPhi   = 0.0;

    double nuRecoilE = 0.0;

    double alpha = 0.0, beta = 0.0;

    TVector3 dmInitVExpVec;
    TVector3 nuFinVTmpExpVec;
    TVector3 nuFinVExpVec;
    
    int flgFF=0,vflg=0,group;//~修正

    double nuMList[15] = { };

    //SI
    nuMList[0]=0.932*(0.989*12.+0.011*13);//C
    nuMList[1]=0.932*(0.950*32.+0.008*33.+0.042*34.);//S
    nuMList[2]=0.932*(0.5069*79.+0.4931*81.);//Br
    nuMList[3]=0.932*127.;//I
    nuMList[11]=0.932*(0.5184*107.+0.48161*109.);//Ag
    //SD
    nuMList[10]=0.932*19.;//F

    FILE* fp1 = fopen("data/dist_blue10GeV.dat", "wt");//output file
    FILE* fp2 = fopen("data/output2_dist_blue10GeV.dat", "wt");//output file
    int atom = 10;//here
    printf("atom is %d (0=C, 1=S, 2=Br, 3=I, 10=F, 11=Ag)\n",atom);

    double nuM = nuMList[ atom ];
    double dmM = 10.;//DM mass
    // printf("dmM=%5.1f #event=10^%2.0lf\n", dmM, log10f(snum));

    // open input file
    TFile inputFile( input.c_str( ) );
    TTree* pInTree = dynamic_cast< TTree* >( inputFile.Get( "tree" ) );
    if( pInTree == nullptr ) {
        std::cerr << "Failed to read input file..." << std::endl;
        abort( );
    }

    double velocity = 0.0, theta = 0.0, phi = 0.0;
    pInTree->SetBranchAddress( "velocity", &velocity );
    pInTree->SetBranchAddress( "theta",    &theta    );
    pInTree->SetBranchAddress( "phi",      &phi      );

    // open output file
    TFile outputFile( output.c_str( ), "RECREATE" );
    TTree* pOutTree = new TTree( "tree", "tree" );
    pOutTree->SetDirectory( &outputFile );
    pOutTree->Branch( "dmM",           &dmM            );
    pOutTree->Branch( "nuM",           &nuM            );
    pOutTree->Branch( "atom",          &atom           );
    pOutTree->Branch( "dmInjV",        &velocity       );
    pOutTree->Branch( "dmInjTheta",    &theta          );
    pOutTree->Branch( "dmInjPhi",      &phi            );
    pOutTree->Branch( "nuRecCosTheta", &nuCosScatTheta );
    pOutTree->Branch( "nuRecSinPhi",   &nuSinScatPhi   );
    pOutTree->Branch( "nuRecV",        &nuFinVExp      );
    pOutTree->Branch( "nuRecE",        &nuRecoilE      );

    int totEvt = pInTree->GetEntries( );
    for( int evt = 0; evt < totEvt; ++evt ) {
        pInTree->GetEntry( evt );
        
        dmInitVExpVec = getDmVVec( velocity, theta, phi );
        dmInitVExp = dmInitVExpVec.Mag( );
        
        sysRelativeV = dmM / (dmM + nuM) * dmInitVExp;
        cosScatThetaCom = gRandom->Uniform( -1.0, 1.0 );

        dmCosScatTheta = getDmCosScatTheta( dmM, nuM, cosScatThetaCom );
        nuCosScatTheta = getNuCosScatTheta( cosScatThetaCom );

        dmFinVCom = nuM / ( dmM + nuM ) * dmInitVExp * gamFact( dmInitVExp );

        dmFinVExp    = getDmVExp( sysRelativeV, dmFinVCom, cosScatThetaCom, dmCosScatTheta );
        nuFinVExp    = getNuFinVExp( dmInitVExp, dmFinVExp, dmM, nuM );
        dmSinScatPhi = getSinDmScatPhi( );
        nuSinScatPhi = -dmSinScatPhi;
      
        nuSinScatTheta = sqrt( 1.0 - sq( nuCosScatTheta ) );
        nuCosScatPhi   = sqrt( 1.0 - sq( nuSinScatPhi   ) );

        nuFinVTmpExpVec.SetX( randomPM( ) * nuFinVExp * nuSinScatTheta * nuCosScatPhi );
        nuFinVTmpExpVec.SetY( randomPM( ) * nuFinVExp * nuSinScatTheta * nuSinScatPhi );
        nuFinVTmpExpVec.SetZ(               nuFinVExp * nuCosScatTheta                );

        alpha        = getRotAngleAlpha( dmInitVExpVec );
        beta         = getRotAngleBeta ( dmInitVExpVec );
        nuFinVExpVec = rotation( nuFinVTmpExpVec, alpha, beta );

        nuCosScatTheta = nuFinVExpVec.Z( ) / nuFinVExp;
        nuSinScatTheta = sqrt( 1.0 - sq( nuCosScatTheta ) );
        nuSinScatPhi   = randomPM( ) * nuFinVExpVec.Y( ) / nuFinVExp / nuSinScatTheta;

        nuRecoilE = getRecoilEnergy( nuM, nuFinVExp );
        if( gRandom->Rndm( ) <= getFormFactorSq( nuRecoilE, atom ) ){
            group = 1;
            fprintf(fp1,"%d %d %lf %lf %lf %lf\n",
                    group, atom, dmM, nuCosScatTheta, nuSinScatPhi, nuRecoilE );
            //1     2     3    4    5                    6
            fprintf(fp2,"%lf %lf %lf %lf %lf\n",
                    dmInitVExpVec.X( ), dmInitVExpVec.Y( ), dmInitVExpVec.Z( ), nuM, nuFinVExp );

            pOutTree->Fill( );
        }
    }//end of event number loop

    pOutTree->Write( );

    printf("\n");
    printf("%d events rejected.\n", flgFF);
    printf("%s", ctime(&t));//time to finish calc.
    fclose(fp1);
    fclose(fp2);

    return 0;
}


