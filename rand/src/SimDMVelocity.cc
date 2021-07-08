#include "CRDMFunc.h"
#include "TF3.h"
#include "TRandom3.h"
#include "Math/IntegratorOptions.h"

const String CR_P_INPUT   = "./data/GALPROP_p.txt";

//////////////////////////////////////////////////////////////////
//
// main
//
//////////////////////////////////////////////////////////////////
int main( int argc, char** argv )
{
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

    TF3* pFunc = nullptr;
    if( profile == "NFW" ) {
        pFunc = new TF3( "flux", getDMNFWFluxV, -0.5*TMath::Pi( ), 0.5*TMath::Pi( ), 0.0, 2.0 * TMath::Pi( ), 0.0, V_LIGHT, 8, 3 );
    }
    else if( profile == "IT" ) {
        pFunc = new TF3( "flux", getDMIsoThermalFluxV, -0.5*TMath::Pi( ), 0.5*TMath::Pi( ), 0.0, 2.0 * TMath::Pi( ), 0.0, V_LIGHT, 8, 3 );
    }
    else if( profile == "EIN" ) {
        pFunc = new TF3( "flux", getDMEinastoFluxV, -0.5*TMath::Pi( ), 0.5*TMath::Pi( ), 0.0, 2.0 * TMath::Pi( ), 0.0, V_LIGHT, 8, 3 );
    }

    pFunc->SetParameter( 0, PROTON_MASS       );
    pFunc->SetParameter( 1, los               );
    pFunc->SetParameter( 2, dmM               );
    pFunc->SetParameter( 3, xsection          ); // assume sigma_DM = 10^-30 [1/cm^2]
    pFunc->SetParameter( 4, rhoScaleKPC       );
    pFunc->SetParameter( 5, DM_R_SCALE        );
    pFunc->SetParameter( 6, SUN_DISTANCE      );
    pFunc->SetParameter( 7, LAMBDA_P          );
    double totValue = PC2CM*PC2CM*pFunc->Integral( -0.5*TMath::Pi( ), 0.5*TMath::Pi( ), 0.0, 2.0 * TMath::Pi( ), 0.0, V_LIGHT );

    pFunc->SetNpx(200);
    pFunc->SetNpy(200);
    pFunc->SetNpz(5000);

    gRandom->SetSeed( 0 );

    TFile file( fileName.c_str( ), "RECREATE" );
    TTree* pTree = new TTree( "tree", "tree" );
    double theta = 0.0, phi = 0.0, velocity = 0.0;
    double vLight      = V_LIGHT;
    double pM          = PROTON_MASS;
    double rScaleKPC   = DM_R_SCALE;
    double sunDistance = SUN_DISTANCE;
    double lambdaP     = LAMBDA_P;
    pTree->Branch( "theta",       &theta       );
    pTree->Branch( "phi",         &phi         );
    pTree->Branch( "velocity",    &velocity    );
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
    
    for( int i = 0; i < evtNum; ++i ) {
        printProgressBar( i, evtNum );
        pFunc->GetRandom3( theta, phi, velocity );
        pTree->Fill( );
    }

    pTree->Write( );

    return 0;
}

