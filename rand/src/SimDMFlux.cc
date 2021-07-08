#include "CRDMFunc.h"
#include "TF3.h"
#include "TRandom3.h"

const String CR_P_INPUT   = "./data/GALPROP_p.txt";

int main( int argc, char** argv )
{
    if( argc != 5 ) {
        std::cerr << "INPUT ERROR" << std::endl;
        std::cerr << "./SimDMFlux [DM energy:Tx (GeV)] [DM mass (GeV)] [profile (NFW or IT or EIN)] [events]" << std::endl;
        abort( );
    }

    double dmE     = std::stod(argv[1]);
    double dmM     = std::stod(argv[2]);
    String profile = argv[3];
    int    evtNum  = std::stoi(argv[4]);
    
    if( profile != "NFW" && profile != "IT" ) {
        std::cerr << "The dark matter profile: " << profile << " is not supported by this program..." << std::endl;
        std::cerr << "Choose NFW or IT (isothermal) EIN (Einasto)." << std::endl;
        abort( );
    }
    
    double rhoScaleKPC = DM_RHO_SCALE / ( PC2CM*PC2CM*PC2CM ); // converted to [GeV / kpc^3]

    // get proton flux distribution
    if( readGalprop( CR_P_INPUT ) == false ) {
        std::cerr << "Failed to read proton flux from file: " << CR_P_INPUT << std::endl;
        abort( );
    }

    // DEBUG
    DEBUG(getTIntegral(PROTON_MASS, dmE, dmM));
    // return 0;

    TF3* pFunc = nullptr;
    if( profile == "NFW" ) {
        pFunc = new TF3( "flux", getDMNFWFlux, -0.5*TMath::Pi( ), 0.5*TMath::Pi( ), 0.0, 2.0 * TMath::Pi( ), 0.0, LOS_LIMIT, 8, 3 );
    }
    else if( profile == "IT" ) {
        pFunc = new TF3( "flux", getDMIsoThermalFlux, -0.5*TMath::Pi( ), 0.5*TMath::Pi( ), 0.0, 2.0 * TMath::Pi( ), 0.0, LOS_LIMIT, 8, 3 );
    }
    else if( profile == "EIN" ) {
        pFunc = new TF3( "flux", getDMEinastoFlux, -0.5*TMath::Pi( ), 0.5*TMath::Pi( ), 0.0, 2.0 * TMath::Pi( ), 0.0, LOS_LIMIT, 8, 3 );         
    }

    pFunc->SetParameter( 0, PROTON_MASS       );
    pFunc->SetParameter( 1, dmE               );
    pFunc->SetParameter( 2, dmM               );
    pFunc->SetParameter( 3, PC2CM*PC2CM*1e-30 ); // assume sigma_DM = 10^-30 [1/cm^2]
    pFunc->SetParameter( 4, rhoScaleKPC       );
    pFunc->SetParameter( 5, DM_R_SCALE        );
    pFunc->SetParameter( 6, SUN_DISTANCE      );
    pFunc->SetParameter( 7, LAMBDA_P          );
    
    double totValue = PC2CM*PC2CM*pFunc->Integral( -0.5*TMath::Pi( ), 0.5*TMath::Pi( ), 0.0, 2.0 * TMath::Pi( ), 0.0, LOS_LIMIT );
    // DEBUG(totValue);

    // pFunc->SetNpx(500);
    // pFunc->SetNpy(500);
    pFunc->SetNpz(100);

    gRandom->SetSeed( 0 );

    TFile file( "output.root", "RECREATE" );
    TTree* pTree = new TTree( "tree", "tree" );
    double theta = 0.0, phi = 0.0, los = 0.0;
    pTree->Branch( "theta",  &theta );
    pTree->Branch( "phi",    &phi );
    pTree->Branch( "los",    &los );
    pTree->Branch( "invWeight", &totValue );
    
    for( int i = 0; i < evtNum; ++i ) {
        printProgressBar( i, evtNum );
        pFunc->GetRandom3( theta, phi, los );
        pTree->Fill( );
    }

    pTree->Write( );
    
    return 0;
}

