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
    if( argc != 5 ) {
        std::cerr << "INPUT ERROR" << std::endl;
        std::cerr << "./SimDMVelocity [l.o.s. (kpc)] [DM mass (GeV)] [The number of events] [output filename]" << std::endl;
        abort( );
    }

    double los      = std::stod(argv[1]);
    double dmM      = std::stod(argv[2]);
    int    evtNum   = std::stoi(argv[3]);
    String fileName = argv[4];

    // get proton flux distribution
    if( readGalprop( CR_P_INPUT ) == false ) {
        std::cerr << "Failed to read proton flux from file: " << CR_P_INPUT << std::endl;
        abort( );
    }
    
    double xsection = PC2CM*PC2CM*1e-30;
    double rhoScaleKPC = DM_RHO_SCALE / ( PC2CM*PC2CM*PC2CM ); // converted to [GeV / kpc^3]

    TF3 func( "flux", getDMFluxV, -0.5*TMath::Pi( ), 0.5*TMath::Pi( ), 0.0, 2.0 * TMath::Pi( ), 0.0, V_LIGHT, 8, 3 );
    func.SetParameter( 0, PROTON_MASS       );
    func.SetParameter( 1, los               );
    func.SetParameter( 2, dmM               );
    func.SetParameter( 3, xsection          ); // assume sigma_DM = 10^-30 [1/cm^2]
    func.SetParameter( 4, rhoScaleKPC       );
    func.SetParameter( 5, DM_R_SCALE        );
    func.SetParameter( 6, SUN_DISTANCE      );
    func.SetParameter( 7, LAMBDA_P          );
    double totValue = PC2CM*PC2CM*func.Integral( -0.5*TMath::Pi( ), 0.5*TMath::Pi( ), 0.0, 2.0 * TMath::Pi( ), 0.0, V_LIGHT );

    // func.SetNpx(100);
    // func.SetNpy(100);
    // func.SetNpz(200);

    // Test
    // TCanvas cvs( "cvs", "cvs", 800, 600 );
    // func.Draw("lego");
    // cvs.SaveAs("test.png");

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
        // printProgressBar( i, evtNum );
        func.GetRandom3( theta, phi, velocity );
        pTree->Fill( );
    }

    pTree->Write( );

    return 0;
}

