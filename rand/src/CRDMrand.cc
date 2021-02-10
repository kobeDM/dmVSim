#include "inc/shinclude.h"
#include "TF3.h"
#include "TRandom3.h"
#include "TTree.h"
#include "Math/IntegratorOptions.h"

const double LOS_LIMIT    = 10.0;   // line-of-sight kpc
const double PROTON_MASS  = 0.938; // GeV
const double DM_XS        = 1e-35; // 10^-35 cm^2

const double DM_RHO_SCALE = 0.403; // GeV / cm^3
const double DM_R_SCALE   = 12.53; // kpc
const double SUN_DISTANCE = 8.0;   // kpc
const double LAMBDA_P     = 0.77;  // GeV

const double PC2CM        = 3.086*1e21; // kpc->cm

const double V_LIGHT      = 299792.458; // km/s

const String CR_P_INPUT   = "./data/GALPROP_p.txt";

std::vector< double > gEneAray;
std::vector< double > gFlxAray;

double getDMFlux       ( );

double getTMax         ( double pE,
                         double pM,
                         double dmM );

double getTMaxInv      ( double pE,
                         double pM,
                         double dmM );

double getTMin         ( double pM,
                         double dmE,
                         double dmM );

double getTIntegral    ( double pM,
                         double dmE,
                         double dmM );

double getDMFlux       ( double theta,
                         double phi,
                         double los,
                         double pM,
                         double dmE,
                         double dmM,
                         double dmXS,
                         double dmDScale,
                         double dmRScale,
                         double sunDist,
                         double lambdaP );

double getDMFlux       ( double* x,
                         double* par );

double getDMFluxV      ( double* x,
                         double* par );

void   printProgressBar( const int& index, const int& total );

bool   readGalprop     ( const String& input );

double getDiffFlux     ( const double& energy );


//////////////////////////////////////////////////////////////////
//
// main
//
//////////////////////////////////////////////////////////////////
int main( int argc, char** argv )
{
    if( argc != 5 ) {
        std::cerr << "INPUT ERROR" << std::endl;
        std::cerr << "./CRDMrand [l.o.s. (kpc)] [DM mass (GeV)] [The number of events] [output filename]" << std::endl;
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

    func.SetNpx(100);
    func.SetNpy(100);
    func.SetNpz(100);

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
        func.GetRandom3( theta, phi, velocity );
        pTree->Fill( );
    }

    pTree->Write( );

    return 0;
}


//////////////////////////////////////////////////////////////////
//
// sub functions
//
//////////////////////////////////////////////////////////////////
double getTMax( double pE,
                double pM,
                double dmM )
{
    if( pE < 0.0 || pM < 0.0 || dmM < 0.0 ) return -100.0;

    double num = pE*pE + 2*pM*pE;
    double den = pE + ((pM+dmM)*(pM+dmM)/(2.0*dmM));
    return num / den;
}

double getTMaxInv( double pE,
                   double pM,
                   double dmM )
{
    if( pE < 0.0 || pM < 0.0 || dmM < 0.0 ) return -100.0;
    return getDiffFlux( pE ) / getTMax( pE, pM, dmM );
}

double getTMaxInv( double* x,
                   double* par )
{
    double pE  = x[0];
    double pM  = par[0];
    double dmM = par[1];
    
    return getTMaxInv( pE, pM, dmM );
}


double getTMin( double pM,
                double dmE,
                double dmM )
{
    if( pM < 0.0 || dmE < 0.0 || dmM < 0.0 ) return -100.0;

    double first = dmE / 2.0 - pM;
    double tmp = 1.0 + (2.0*dmE/dmM)*pow((pM+dmM)/(2.0*pM-dmE) , 2 );
    
    double second = 0.0;
    if( first > 0.0 ) second = 1.0 + sqrt(tmp);
    else              second = 1.0 - sqrt(tmp);

    return first * second;
}


double getTIntegral( double pM,
                     double dmE,
                     double dmM )
{
    if( pM < 0.0 || dmE < 0.0 || dmM < 0.0 ) return -100.0;
    
    TF1 func( "TInt", getTMaxInv, getTMin( pM, dmE, dmM ), 10000000.0, 2, 1 );
    func.SetParameter( 0, pM );
    func.SetParameter( 1, dmM );

    ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6);
    ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1.E-6);
    ROOT::Math::IntegratorMultiDimOptions::SetDefaultRelTolerance(1.E-6);
    ROOT::Math::IntegratorMultiDimOptions::SetDefaultAbsTolerance(1.E-6);

    return func.Integral( getTMin( pM, dmE, dmM ), 10000000.0 );
}


double getDMFlux( double theta,
                  double phi,
                  double los,
                  double pM,
                  double dmE,
                  double dmM,
                  double dmXS,
                  double dmDScale,
                  double dmRScale,
                  double sunDist,
                  double lambdaP )
{
    if( pM < 0.0 || dmE < 0.0 || dmM < 0.0 || dmXS < 0.0 || dmDScale < 0.0 || dmRScale < 0.0 || sunDist < 0.0 )
        return -100.0;

    double TIntegral = getTIntegral( pM, dmE, dmM );

    double formFactor = 1.0 / pow( 1.0 + (2.0 * dmM * dmE) / (lambdaP*lambdaP), 2 );
    
    double sinTheta = sin( theta );
    double cosTheta = cos( theta );
    double cosPhi   = cos( phi );
    double rsq      = los*los + sunDist*sunDist - 2*los*sunDist*cosTheta*cosPhi;
    double r        = sqrt( fabs( rsq ) );

    double relR     = r / dmRScale;
    double rhoNFW   = dmDScale / (relR*pow( 1 + relR, 2 ));
    double jacobian = los*los*cosTheta;
    double fluxCorr = 1.0 / (4.0*TMath::Pi()*los*los);

    return TIntegral * fluxCorr * jacobian * dmXS / dmM * rhoNFW * formFactor * formFactor;
}


double getDMFlux( double* x,
                  double* par )
{
    double theta    = x[0];
    double phi      = x[1];
    double los      = x[2];
    double pM       = par[0];
    double dmE      = par[1];
    double dmM      = par[2];
    double dmXS     = par[3];
    double dmDScale = par[4];
    double dmRScale = par[5];
    double sunDist  = par[6];
    double lambdaP  = par[7];

    return getDMFlux( theta, phi, los, pM, dmE, dmM, dmXS, dmDScale, dmRScale, sunDist, lambdaP );
}


double getDMFluxV( double* x,
                   double* par )
{
    double theta    = x[0];
    double phi      = x[1];
    double velo     = x[2];
    double pM       = par[0];
    double los      = par[1];
    double dmM      = par[2];
    double dmXS     = par[3];
    double dmDScale = par[4];
    double dmRScale = par[5];
    double sunDist  = par[6];
    double lambdaP  = par[7];

    double gamma    = 1.0 / sqrt( 1 - pow( velo / V_LIGHT , 2) );
    double dmE      = dmM * (gamma - 1.0);

    double flux = getDMFlux( theta, phi, los, pM, dmE, dmM, dmXS, dmDScale, dmRScale, sunDist, lambdaP );

    // return velo * velo * flux * dmM * dmM * gamma * gamma * gamma; // Note!!: the local dark matter density is ignored at this stage... 
    return flux * dmM * dmM * gamma * gamma * gamma; // Note!!: the local dark matter density is ignored at this stage... 
}


void printProgressBar( const int& index, const int& total )
{
    if( index % 100 == 0 ) {
        String printBar = " [";
        double progress = static_cast< double >( index ) / static_cast< double >( total );
        for( int bar = 0; bar < 20; ++bar ) {
            double currentFraction = static_cast< double >( bar ) * 0.05;
            if( progress > currentFraction ) printBar += "/";
            else printBar += ".";
        }
        printBar += "] ";
        double percent = 100.0 * progress;
        StringStream percentSS;
        percentSS << std::setprecision( 2 ) << percent;
        String text = printBar + " ";
        text += percentSS.str( );
        std::cout << std::flush; 
        std::cout << text << "%\r" << std::flush; 
    }
    return;
}
    
bool readGalprop( const String& input )
{
    std::ifstream ifs;
    ifs.open( input );
    if( ifs.is_open( ) == false ) return false;

    double energy = 0.0, flux = 0.0;
    while( !ifs.eof( ) ) {
        String line = "";
        std::getline( ifs, line );
        if( line.length( ) <= 0 || strncmp( line.c_str( ), "#", 1 ) == 0 ) continue;
        
        StringStream ss( line );
        ss >> energy >> flux;
        energy *= 0.001; // MeV -> GeV
        gEneAray.push_back( energy );
        gFlxAray.push_back( flux / energy / energy );
    }

    return true;
}

double getDiffFlux( const double& energy )
{
    double flux = 0.0;
    double eneMinBin = 0.0, eneMaxBin = 0.0;
    double flxMinBin = 0.0, flxMaxBin = 0.0;

    int idx = 0;
    for( auto eneBin : gEneAray ) {
        if( idx == 0 ) {
            if( eneBin > energy ) {
                DEBUG("error: energy is too low.");
                break;
            }
            eneMinBin = eneBin;
            idx++;
            continue;
        }

        if( eneBin > energy ) {
            eneMaxBin = energy;
            flxMinBin = gFlxAray.at( idx - 1 );
            flxMaxBin = gFlxAray.at( idx );
            break;
        }
        
        idx++;
    }

    flux = (flxMaxBin - flxMinBin) / (eneMaxBin - eneMinBin) * (energy - eneMinBin) + flxMinBin;
    return flux;
}
