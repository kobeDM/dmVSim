#include "CRDMFunc.h"
#include "TF3.h"
#include "TRandom3.h"
#include "TTree.h"
#include "Math/IntegratorOptions.h"

double   randomPM         ( );
double   sq               ( const double& x );
double   gamFact          ( const double& velo );

double   getSinDmScatPhi  ( );

double   getDmCosScatTheta( const double& dmM,
                            const double& nuM,
                            const double& cosScatAngleCom );

double   getNuCosScatTheta( const double& cosScatThetaCom );

TVector3 getDmVVec        ( const double& velocity,
                            const double& theta,
                            const double& phi );

double   getDmVExp        ( const double& dmInitVCom,
                            const double& dmFinVCom,
                            const double& cosScatThetaCom,
                            const double& cosScatThetaTmpExp );

double   getNuFinVExp     ( const double& dmInitVExp,
                            const double& dmFinVExp,
                            const double& dmM,
                            const double& nuM );

double   getRecoilEnergy  ( const double& nuM,
                            const double& nuVExp );

// double   getRotAngleAlpha ( TVector3      dmInitVExpVec );

// double   getRotAngleBeta  ( TVector3      dmInitVExpVec );

double   getRotAngleBeta  ( TVector3      dmInitVExpVec );

double   getRotAngleGamma ( TVector3      dmInitVExpVec );

// TVector3 rotation         ( TVector3      u,
//                             const double& al,
//                             const double& be );

TVector3 rotation         ( TVector3      u,
                            const double& be,
                            const double& ga );

double   getFormFactorSq  ( const double& ER,
                            const int&    atom );

//////////////////////////////////////////////////////////////////
//
// main
//
//////////////////////////////////////////////////////////////////
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
        std::cerr << "./SimNuclRecoil [input filename] [output filename]" << std::endl;
        abort( );
    }
    String input  = argv[1];
    String output = argv[2];
    
    double sysRelativeV = 0.0;
    double dmFinVCom  = 0.0;
    double nuFinVCom  = 0.0;
    double dmInitVExp = 0.0, dmInitVExpX = 0.0, dmInitVExpY = 0.0 , dmInitVExpZ = 0.0;
    double dmFinVExp  = 0.0, dmFinVExpX  = 0.0, dmFinVExpY  = 0.0 , dmFinVExpZ  = 0.0;
    double nuFinVExp  = 0.0, nuFinVExpX  = 0.0, nuFinVExpY  = 0.0 , nuFinVExpZ  = 0.0;

    double nuFinVTmpExp  = 0.0, nuFinVTmpExpX  = 0.0, nuFinVTmpExpY  = 0.0 , nuFinVTmpExpZ  = 0.0;

    double scatThetaCom    = 0.0;
    double cosScatThetaCom = 0.0;

    double dmScatThetaExp     = 0.0, nuScatThetaExp    = 0.0;
    double dmScatPhiExp       = 0.0, nuScatPhiExp      = 0.0;
    double dmSinScatThetaExp  = 0.0, nuSinScatThetaExp = 0.0;
    double dmCosScatThetaExp  = 0.0, nuCosScatThetaExp = 0.0;
    double dmSinScatPhiExp    = 0.0, nuSinScatPhiExp   = 0.0;
    double dmCosScatPhiExp    = 0.0, nuCosScatPhiExp   = 0.0;

    double dmScatThetaTmpExp     = 0.0, nuScatThetaTmpExp    = 0.0;
    double dmScatPhiTmpExp       = 0.0, nuScatPhiTmpExp      = 0.0;
    double dmSinScatThetaTmpExp  = 0.0, nuSinScatThetaTmpExp = 0.0;
    double dmCosScatThetaTmpExp  = 0.0, nuCosScatThetaTmpExp = 0.0;
    double dmSinScatPhiTmpExp    = 0.0, nuSinScatPhiTmpExp   = 0.0;
    double dmCosScatPhiTmpExp    = 0.0, nuCosScatPhiTmpExp   = 0.0;

    double nuRecoilE    = 0.0;
    double formFactorSq = 0.0;

    // double alpha = 0.0, beta = 0.0;
    double beta = 0.0, gamma = 0.0;

    TVector3 dmInitVExpVec;
    TVector3 nuFinVTmpExpVec;
    TVector3 dmFinVExpVec;
    TVector3 nuFinVExpVec;

    double rndm = 0.0;
    
    int group = 0;
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
    pOutTree->Branch( "dmM",              &dmM               );
    pOutTree->Branch( "nuM",              &nuM               );
    pOutTree->Branch( "atom",             &atom              );
    pOutTree->Branch( "dmInjV_debug",     &velocity          );
    pOutTree->Branch( "dmInjV",           &dmInitVExp        );
    pOutTree->Branch( "dmInjVX",          &dmInitVExpX       );
    pOutTree->Branch( "dmInjVY",          &dmInitVExpY       );
    pOutTree->Branch( "dmInjVZ",          &dmInitVExpZ       );
    pOutTree->Branch( "dmInjTheta",       &theta             );
    pOutTree->Branch( "dmInjPhi",         &phi               );

    pOutTree->Branch( "dmOutPhiTmpExp",   &dmScatPhiTmpExp   );
    pOutTree->Branch( "dmOutThetaTmpExp", &dmScatThetaTmpExp );
    pOutTree->Branch( "dmOutCosPhiTmpExp",   &dmCosScatPhiTmpExp   );
    pOutTree->Branch( "dmOutCosThetaTmpExp", &dmCosScatThetaTmpExp );

    pOutTree->Branch( "nuRecPhiTmpExp",   &nuScatPhiTmpExp   );
    pOutTree->Branch( "nuRecThetaTmpExp", &nuScatThetaTmpExp );
    pOutTree->Branch( "nuRecCosPhiTmpExp",   &nuCosScatPhiTmpExp   );
    pOutTree->Branch( "nuRecCosThetaTmpExp", &nuCosScatThetaTmpExp );
    pOutTree->Branch( "nuRecVTmpExpX",    &nuFinVTmpExpX  );
    pOutTree->Branch( "nuRecVTmpExpY",    &nuFinVTmpExpY  );
    pOutTree->Branch( "nuRecVTmpExpZ",    &nuFinVTmpExpZ  );

    pOutTree->Branch( "dmOutV",           &dmFinVExp         );
    pOutTree->Branch( "dmOutVX",          &dmFinVExpX        );
    pOutTree->Branch( "dmOutVY",          &dmFinVExpY        );
    pOutTree->Branch( "dmOutVZ",          &dmFinVExpZ        );
    pOutTree->Branch( "dmOutTheta",       &dmScatThetaExp    );
    pOutTree->Branch( "dmOutPhi",         &dmScatPhiExp      );
    pOutTree->Branch( "dmOutCosTheta",    &dmCosScatThetaExp );
    pOutTree->Branch( "dmOutSinPhi",      &dmSinScatPhiExp   );
                                        
    pOutTree->Branch( "nuRecTheta",       &nuScatThetaExp    );
    pOutTree->Branch( "nuRecPhi",         &nuScatPhiExp      );
    pOutTree->Branch( "nuRecCosTheta",    &nuCosScatThetaExp );
    pOutTree->Branch( "nuRecSinPhi",      &nuSinScatPhiExp   );
    pOutTree->Branch( "nuRecV",           &nuFinVExp         );
    pOutTree->Branch( "nuRecVX",          &nuFinVExpX        );
    pOutTree->Branch( "nuRecVY",          &nuFinVExpY        );
    pOutTree->Branch( "nuRecVZ",          &nuFinVExpZ        );
    pOutTree->Branch( "nuRecE",           &nuRecoilE         );
                                        
    pOutTree->Branch( "sysRelV",          &sysRelativeV      );

    pOutTree->Branch( "scatThetaCom",     &scatThetaCom      );
    pOutTree->Branch( "cosScatThetaCom",  &cosScatThetaCom   );
                                                              
    // pOutTree->Branch( "alpha",            &alpha             );
    // pOutTree->Branch( "beta",             &beta              );
    pOutTree->Branch( "beta",             &beta              );
    pOutTree->Branch( "gamma",            &gamma             );
    pOutTree->Branch( "formFactorSq",     &formFactorSq      );

    pOutTree->Branch( "rndm",             &rndm              );

    int totEvt = pInTree->GetEntries( );
    for( int evt = 0; evt < totEvt; ++evt ) {
        printProgressBar( evt, totEvt );
        pInTree->GetEntry( evt );
        
        dmInitVExpVec = getDmVVec( velocity, theta, phi );
        dmInitVExpX   = dmInitVExpVec.X( );
        dmInitVExpY   = dmInitVExpVec.Y( );
        dmInitVExpZ   = dmInitVExpVec.Z( );
        dmInitVExp    = dmInitVExpVec.Mag( );
        
        sysRelativeV = dmM / (dmM + nuM) * dmInitVExp;
        scatThetaCom = gRandom->Rndm( ) * 1.0 * PI;
        cosScatThetaCom = cos( scatThetaCom );

        dmCosScatThetaTmpExp = getDmCosScatTheta( dmM, nuM, cosScatThetaCom );
        dmSinScatThetaTmpExp = sqrt( 1.0 - sq( dmCosScatThetaTmpExp ) );
        nuCosScatThetaTmpExp = getNuCosScatTheta( cosScatThetaCom );
        nuSinScatThetaTmpExp = sqrt( 1.0 - sq( nuCosScatThetaTmpExp ) );
        dmScatThetaTmpExp    = asin( dmSinScatThetaTmpExp );
        nuScatThetaTmpExp    = asin( nuSinScatThetaTmpExp );

        dmSinScatPhiTmpExp = getSinDmScatPhi( );
        dmCosScatPhiTmpExp = sqrt( 1.0 - sq( dmSinScatPhiTmpExp ) );
        dmScatPhiTmpExp    = asin( dmSinScatPhiTmpExp );

        nuSinScatPhiTmpExp = sin( PI - dmScatPhiTmpExp );
        nuCosScatPhiTmpExp = sqrt( 1.0 - sq( nuSinScatPhiTmpExp ) );
        nuScatPhiTmpExp    = asin( nuSinScatPhiTmpExp );

        dmFinVCom = nuM / ( dmM + nuM ) * dmInitVExp * gamFact( dmInitVExp );
        dmFinVExp = getDmVExp( sysRelativeV, dmFinVCom, cosScatThetaCom, dmCosScatThetaTmpExp );
        nuFinVExp = getNuFinVExp( dmInitVExp, dmFinVExp, dmM, nuM );

        // nuFinVTmpExpVec.SetX( randomPM( ) * nuFinVExp * nuSinScatThetaTmpExp * nuCosScatPhiTmpExp );
        // nuFinVTmpExpVec.SetY( randomPM( ) * nuFinVExp * nuSinScatThetaTmpExp * nuSinScatPhiTmpExp );
        // nuFinVTmpExpVec.SetZ( randomPM( ) * nuFinVExp * nuCosScatThetaTmpExp                      );

        // nuFinVTmpExpVec.SetX( nuFinVExp * nuSinScatThetaTmpExp * nuCosScatPhiTmpExp );
        // nuFinVTmpExpVec.SetY( nuFinVExp * nuSinScatThetaTmpExp * nuSinScatPhiTmpExp );
        // nuFinVTmpExpVec.SetZ( nuFinVExp * nuCosScatThetaTmpExp                      );

        nuFinVTmpExpVec.SetY( randomPM( ) * nuFinVExp * nuSinScatThetaTmpExp * nuCosScatPhiTmpExp );
        nuFinVTmpExpVec.SetZ( randomPM( ) * nuFinVExp * nuSinScatThetaTmpExp * nuSinScatPhiTmpExp );
        nuFinVTmpExpVec.SetX( nuFinVExp * nuCosScatThetaTmpExp                      );


        nuFinVTmpExpX = nuFinVTmpExpVec.X( );
        nuFinVTmpExpY = nuFinVTmpExpVec.Y( );
        nuFinVTmpExpZ = nuFinVTmpExpVec.Z( );

        // alpha        = getRotAngleAlpha( dmInitVExpVec );
        // beta         = getRotAngleBeta ( dmInitVExpVec );
        // nuFinVExpVec = rotation( nuFinVTmpExpVec, alpha, beta );

        beta         = getRotAngleBeta ( dmInitVExpVec );
        gamma        = getRotAngleGamma( dmInitVExpVec );
        nuFinVExpVec = rotation( nuFinVTmpExpVec, beta, gamma );


        nuFinVExpX = nuFinVExpVec.X( );
        nuFinVExpY = nuFinVExpVec.Y( );
        nuFinVExpZ = nuFinVExpVec.Z( );
        nuSinScatThetaExp = nuFinVExpVec.Z( ) / nuFinVExp;
        nuCosScatThetaExp = sqrt( 1.0 - sq( nuSinScatThetaExp ) );
        nuSinScatPhiExp   = nuFinVExpVec.Y( ) / nuFinVExp / nuCosScatThetaExp;

        nuScatThetaExp = asin( nuSinScatThetaExp );
        nuScatPhiExp   = asin( nuSinScatPhiExp   );
        if( nuFinVExpVec.X( ) < 0.0 ) nuScatPhiExp = PI - nuScatPhiExp;
        if( nuScatPhiExp < 0.0 ) nuScatPhiExp = nuScatPhiExp + 2.0 * PI;

        nuRecoilE = getRecoilEnergy( nuM, nuFinVExp );
        
        rndm = gRandom->Rndm( );
        formFactorSq = getFormFactorSq( nuRecoilE, atom );
        if( rndm <= formFactorSq ) {
            group = 1;
            fprintf(fp1,"%d %d %lf %lf %lf %lf\n",
                    group, atom, dmM, nuCosScatThetaExp, nuSinScatPhiExp, nuRecoilE );
            //1     2     3    4    5                    6
            fprintf(fp2,"%lf %lf %lf %lf %lf\n",
                    dmInitVExpVec.X( ), dmInitVExpVec.Y( ), dmInitVExpVec.Z( ), nuM, nuFinVExp );

            // pOutTree->Fill( );
        }
        pOutTree->Fill( );

    }//end of event number loop

    pOutTree->Write( );

    printf("\n");
    printf("%s", ctime(&t));//time to finish calc.
    fclose(fp1);
    fclose(fp2);

    return 0;
}


//////////////////////////////////////////////////////////////////
//
// sub functions
//
//////////////////////////////////////////////////////////////////
double randomPM()
{
    return gRandom->Integer( 2 ) ? 1.0 : -1.0;
}

double sq( const double& x )
{
    return pow( x, 2.0 );
}

double gamFact( const double& velo )
{
    // return 1.0 / sqrt( 1 + sq( velo/V_LIGHT ) );
    return 1.0;
}

double getSinDmScatPhi( )
{
    return sin( gRandom->Rndm( ) * 2.0 * PI );
}
double getDmCosScatTheta( const double& dmM, const double& nuM, const double& cosScatAngleCom )
{
    double x = dmM / nuM;
    double dmCosScatTheta = ( 1.0 - sq( cosScatAngleCom ) ) / sq( cosScatAngleCom + x ) + 1.0;
    return 1.0 / sqrt( dmCosScatTheta );//always positive if DM is heavy
}

double getNuCosScatTheta( const double& cosScatThetaCom )
{
    double theta2 = ( PI - acos(cosScatThetaCom) ) * 0.5;
    return cos( theta2 );//pm is automatically fixed.
}

TVector3 getDmVVec( const double& velocity,
                    const double& theta,
                    const double& phi )
{
    TVector3 vec;
    // vec.SetX( -1.0 * velocity * cos(theta) * cos(phi) );
    // vec.SetY( -1.0 * velocity * cos(theta) * sin(phi) );
    // vec.SetZ( -1.0 * velocity * sin(theta)            );
    vec.SetX( velocity * cos(theta) * cos(phi) );
    vec.SetY( velocity * cos(theta) * sin(phi) );
    vec.SetZ( velocity * sin(theta)            );
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

// double getRotAngleAlpha( TVector3 dmInitVExpVec ){
//     //returns rotation angle alpha, which rotates vdmE to (0,0,Norm(vdmE))
//     return -1.0 * asin( dmInitVExpVec.Y( ) / dmInitVExpVec.Mag( ) );
// }

// double getRotAngleBeta( TVector3 dmInitVExpVec ){
//     //returns rotation angle beta, which rotates vdmE to (0,0,Norm(vdmE))
//     return atan2( dmInitVExpVec.X( ), dmInitVExpVec.Z( ) );
// }


double getRotAngleBeta( TVector3 dmInitVExpVec ){
    //returns rotation angle beta, which rotates vdmE to (0,0,Norm(vdmE))
    return -1.0 * asin( dmInitVExpVec.Z( ) / dmInitVExpVec.Mag( ) );
}

double getRotAngleGamma( TVector3 dmInitVExpVec ){
    //returns rotation angle alpha, which rotates vdmE to (0,0,Norm(vdmE))
    return atan2( dmInitVExpVec.Y( ), dmInitVExpVec.X( ) );
}

// //rotation from COM to Earth sys
// TVector3 rotation( TVector3      u,
//                    const double& al,
//                    const double& be )
// {
//     TVector3 ans;
//     ans.SetX(  u.X( ) * cos( be ) + u.Y( ) * sin( al ) * sin( be ) + u.Z( ) * cos( al ) * sin( be ) );
//     ans.SetY(                       u.Y( ) * cos( al )             - u.Z( ) * sin( al )             );
//     ans.SetZ( -u.X( ) * sin( be ) + u.Y( ) * sin( al ) * cos( be ) + u.Z( ) * cos( al ) * cos( be ) );
//     return ans;
// }

//rotation from COM to Earth sys
TVector3 rotation( TVector3      u,
                   const double& be,
                   const double& ga )
{
    TVector3 ans;
    ans.SetX(  u.X( ) * cos( be ) * cos( ga ) - u.Y( ) * sin( ga ) + u.Z( ) * sin( be ) * cos( ga ) );
    ans.SetY(  u.X( ) * cos( be ) * sin( ga ) + u.Y( ) * cos( ga ) + u.Z( ) * sin( be ) * sin( ga ) );
    ans.SetZ( -u.X( ) * sin( be )                                  + u.Z( ) * cos( be )             );
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

