#include "CRDMFunc.h"
#include "TF3.h"
#include "TRandom3.h"

const String CR_P_INPUT   = "./data/GALPROP_p.txt";

int main( int argc, char** argv )
{
    if( argc != 4 ) {
        std::cerr << "INPUT ERROR" << std::endl;
        std::cerr << "./SimDMFlux [DM energy:Tx (GeV)] [DM mass (GeV)]" << std::endl;
        abort( );
    }

    double dmE    = std::stod(argv[1]);
    double dmM    = std::stod(argv[2]);
    int    evtNum = std::stoi(argv[3]);
    
    double rhoScaleKPC = DM_RHO_SCALE / ( PC2CM*PC2CM*PC2CM ); // converted to [GeV / kpc^3]

    // get proton flux distribution
    if( readGalprop( CR_P_INPUT ) == false ) {
        std::cerr << "Failed to read proton flux from file: " << CR_P_INPUT << std::endl;
        abort( );
    }

    // DEBUG
    DEBUG(getTIntegral(PROTON_MASS, dmE, dmM));
    // return 0;

    // TF3 func( "flux", getDMFlux, -0.5*TMath::Pi( ), 0.5*TMath::Pi( ), 0.0, 2.0 * TMath::Pi( ), 0.0, 10.0, 8, 3 );
    TF3 func( "flux", getDMFlux, -0.5*TMath::Pi( ), 0.5*TMath::Pi( ), 0.0, 2.0 * TMath::Pi( ), 0.0, LOS_LIMIT, 8, 3 );
    func.SetParameter( 0, PROTON_MASS       );
    func.SetParameter( 1, dmE               );
    func.SetParameter( 2, dmM               );
    func.SetParameter( 3, PC2CM*PC2CM*1e-30 ); // assume sigma_DM = 10^-30 [1/cm^2]
    func.SetParameter( 4, rhoScaleKPC       );
    func.SetParameter( 5, DM_R_SCALE        );
    func.SetParameter( 6, SUN_DISTANCE      );
    func.SetParameter( 7, LAMBDA_P          );
    
    double totValue = PC2CM*PC2CM*func.Integral( -0.5*TMath::Pi( ), 0.5*TMath::Pi( ), 0.0, 2.0 * TMath::Pi( ), 0.0, LOS_LIMIT );
    // DEBUG(totValue);

    // func.SetNpx(500);
    // func.SetNpy(500);
    func.SetNpz(100);

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
        func.GetRandom3( theta, phi, los );
        pTree->Fill( );
    }

    pTree->Write( );
    
    return 0;
}

