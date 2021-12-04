#include "CRDMFunc.h"
#include "TF3.h"
#include "TRandom3.h"
#include "Math/IntegratorOptions.h"
#include "Math/GaussIntegrator.h"

const String CR_P_INPUT   = "./data/GALPROP_p.txt";

//////////////////////////////////////////////////////////////////
//
// main
//
//////////////////////////////////////////////////////////////////
int main( int argc, char** argv )
{
    // clock_t start = clock( );

    // time_t t = time( nullptr );
    // printf("%s", ctime(&t));//time to start

    if( argc != 6 ) {
        std::cerr << "INPUT ERROR" << std::endl;
        std::cerr << "./SimDMVelocity [l.o.s. (kpc)] [DM mass (GeV)] [profile (NFW or IT or EIN)] [The number of events] [output filename]" << std::endl;
        abort( );
    }

    double los      = std::stod(argv[1]);
    double dmM      = std::stod(argv[2]);
    String profile  = argv[3];
    int    evtNum   = std::stoi(argv[4]);
    String fileName = argv[5];

    if( profile != "NFW" && profile != "IT" && profile != "EIN" ) {
        std::cerr << "The dark matter profile: " << profile << " is not supported by this program..." << std::endl;
        std::cerr << "Choose NFW or IT (isothermal) or EIN (Einasto)." << std::endl;
        abort( );
    }

    // get proton flux distribution
    if( readGalprop( CR_P_INPUT ) == false ) {
        std::cerr << "Failed to read proton flux from file: " << CR_P_INPUT << std::endl;
        abort( );
    }
    
    double xsection = PC2CM*PC2CM*1e-30;
    double rhoScaleKPC = DM_RHO_SCALE / ( PC2CM*PC2CM*PC2CM ); // converted to [GeV / kpc^3]

    TF2* pFuncDir = nullptr;
    if( profile == "NFW" ) {
        pFuncDir = new TF2( "fluxDir", getDMNFWFluxDir, -0.5*TMath::Pi( ), 0.5*TMath::Pi( ), 0.0, 2.0 * TMath::Pi( ), 4, 2 );
    }
    else if( profile == "IT" ) {
        pFuncDir = new TF2( "fluxDir", getDMIsoThermalFluxDir, -0.5*TMath::Pi( ), 0.5*TMath::Pi( ), 0.0, 2.0 * TMath::Pi( ), 4, 2 );
    }
    else if( profile == "EIN" ) {
        pFuncDir = new TF2( "fluxDir", getDMEinastoFluxDir, -0.5*TMath::Pi( ), 0.5*TMath::Pi( ), 0.0, 2.0 * TMath::Pi( ), 4, 2 );
    }

    TF1* pFuncEn = new TF1( "fluxEn", getDMFluxEn, 0.0, V_LIGHT, 4, 1 );

    pFuncEn->SetParameter( 0, PROTON_MASS       );
    pFuncEn->SetParameter( 1, dmM               );
    pFuncEn->SetParameter( 2, xsection          ); // assume sigma_DM = 10^-30 [1/cm^2]
    pFuncEn->SetParameter( 3, LAMBDA_P          );

    pFuncDir->SetParameter( 0, los               );
    pFuncDir->SetParameter( 1, rhoScaleKPC       );
    pFuncDir->SetParameter( 2, DM_R_SCALE        );
    pFuncDir->SetParameter( 3, SUN_DISTANCE      );

    pFuncDir->SetNpx( 36 );
    pFuncDir->SetNpy( 72 );

    // create CDF
    double x[100000], y[100000];
    double xCDF[100000], yCDF[100000];
    double xCDFinv[100000], yCDFinv[100000];

    std::cout << "Create CDF" << std::endl;
    for( int i = 0; i < 100000; ++i ) {
        printProgressBar( i, 100000 );
        double val = 0.0;
        // x[i] = V_LIGHT * (double)i * 0.0001 * 0.1; // 0 < v < 0.1c
        x[i] = V_LIGHT * (double)i * 0.00001; // 0 < v < c
        xCDF[i] = x[i];
        y[i] =  pFuncEn->Eval( x[i] );
        if( i != 0 )
            yCDF[i] =  yCDF[i-1] + y[i];
        else
            yCDF[i] =  y[i];

        // if( i % 100 == 0 ) DEBUG(i);
    }

    for( int i = 100000; i > 0; --i ) {
        xCDFinv[i] = xCDF[i];
        if( i != 100000 )
            yCDFinv[i-1] = yCDFinv[i] + y[i-1];
        else
            yCDFinv[i-1] = y[i-1];
    }

    // normalize CDF
    double yCDFMax    = yCDF[99999];
    double yCDFinvMax = yCDFinv[0];
    for( int i = 0; i < 100000; ++i ) {
        yCDF[i]    = yCDF[i]    / yCDFMax;
        yCDFinv[i] = yCDFinv[i] / yCDFinvMax;

        yCDFinv[i] = 1.0 - yCDFinv[i];
    }

    double totValue = 1.0;

    gRandom->SetSeed( 0 );

    TFile file( fileName.c_str( ), "RECREATE" );
    TTree* pTree = new TTree( "tree", "tree" );
    double theta = 0.0, phi = 0.0, velocity = 0.0, velocityInv = 0.0;
    double vLight      = V_LIGHT;
    double pM          = PROTON_MASS;
    double rScaleKPC   = DM_R_SCALE;
    double sunDistance = SUN_DISTANCE;
    double lambdaP     = LAMBDA_P;

    double rnd         = 0.0;

    pTree->Branch( "theta",       &theta       );
    pTree->Branch( "phi",         &phi         );
    pTree->Branch( "velocity",    &velocity    );
    pTree->Branch( "velocityInv", &velocityInv );
    pTree->Branch( "invWeight",   &totValue    );
    pTree->Branch( "vLight",      &vLight      );

    pTree->Branch( "dmM",         &dmM         );
    pTree->Branch( "los",         &los         );

    pTree->Branch( "pM",          &pM          );
    pTree->Branch( "xsection",    &xsection    );
    pTree->Branch( "rhoScaleKPC", &rhoScaleKPC );
    pTree->Branch( "rScaleKPC",   &rScaleKPC   );
    pTree->Branch( "sunDistance", &sunDistance );
    pTree->Branch( "lambdaP",     &lambdaP     );

    pTree->Branch( "rnd",         &rnd         );
    
    for( int i = 0; i < evtNum; ++i ) {
        printProgressBar( i, evtNum );
        rnd = gRandom->Rndm( );
        pFuncDir->GetRandom2( theta, phi );
        velocity    = getVelocity( xCDF, yCDF, rnd );
        velocityInv = getVelocity( xCDFinv, yCDFinv, rnd );
        
        pTree->Fill( );
    }

    pTree->Write( );

    // TCanvas cvs( "cvs", "cvs", 800, 600 );
    // cvs.SetLogx(1);
    // cvs.SetLogy(1);

    TGraph graph( 30000, x, y );
    graph.SetName("vDist");
    graph.SetMarkerStyle( 20 );
    graph.Write();

    TGraph graphCDF( 30000, xCDF, yCDF );
    graphCDF.SetName("vCDF");
    graphCDF.SetMarkerStyle( 20 );
    graphCDF.Write();

    TGraph graphCDFinv( 30000, xCDFinv, yCDFinv );
    graphCDFinv.SetName("vCDFnv");
    graphCDFinv.SetMarkerStyle( 20 );
    graphCDFinv.Write();

    file.Close( );

    return 0;
}

