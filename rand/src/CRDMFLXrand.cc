#include "inc/shinclude.h"
#include "TF3.h"
#include "TRandom3.h"
#include "TTree.h"

const double LOS_LIMIT    = 10.0;   // line-of-sight kpc
const double PROTON_MASS  = 0.938; // GeV
const double DM_XS        = 1e-35; // 10^-35 cm^2

const double DM_RHO_SCALE = 0.403; // GeV / cm^3
const double DM_R_SCALE   = 12.53; // kpc
const double SUN_DISTANCE = 8.0;   // kpc
const double LAMBDA_P     = 0.77;  // GeV

const double PC2CM        = 3.086*1e21; // kpc->cm

const String CR_P_INPUT   = "/home/higashino/work/repo/dmvsim/rand/GALPROP_p.txt";

double getDMFlux( );

double getTMax( double pE,
                double pM,
                double dmM );

double getTMaxInv( double pE,
                   double pM,
                   double dmM );

double getTMin( double pM,
                double dmE,
                double dmM );

double getTIntegral( double pM,
                     double dmE,
                     double dmM );

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
                  double lambdaP );

double getDMFlux( double* x,
                  double* par );

int main( int argc, char** argv )
{
    if( argc != 4 ) {
        std::cerr << "INPUT ERROR" << std::endl;
        std::cerr << "./CRDMFLXrand [DM energy:Tx (GeV)] [DM mass (GeV)]" << std::endl;
        abort( );
    }

    double dmE    = std::stod(argv[1]);
    double dmM    = std::stod(argv[2]);
    int    evtNum = std::stoi(argv[3]);
    
    double rhoScaleKPC = DM_RHO_SCALE / ( PC2CM*PC2CM*PC2CM ); // converted to [GeV / kpc^3]

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
    DEBUG(totValue);

    func.SetNpx(500);
    func.SetNpy(500);
    func.SetNpz(500);

    gRandom->SetSeed( 0 );

    TFile file( "output.root", "RECREATE" );
    TTree* pTree = new TTree( "tree", "tree" );
    double theta = 0.0, phi = 0.0, los = 0.0;
    pTree->Branch( "theta",  &theta );
    pTree->Branch( "phi",    &phi );
    pTree->Branch( "los",    &los );
    pTree->Branch( "invWeight", &totValue );
    
    for( int i = 0; i < evtNum; ++i ) {
        func.GetRandom3( theta, phi, los );
        pTree->Fill( );
    }

    pTree->Write( );
    
    return 0;
}


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
    return 1.0 / getTMax( pE, pM, dmM );
    // return getTMax( pE, pM, dmM );
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
    
    TF1 func( "TInt", getTMaxInv, getTMin( pM, dmE, dmM ), 1000.0, 2, 1 );
    func.SetParameter( 0, pM );
    func.SetParameter( 1, dmM );

    // TCanvas cvs( "cvsTInt", "cvsTInt", 800, 600 );
    // func.Draw( );
    // cvs.SaveAs("cvsTInt.png");

    return func.Integral( getTMin( pM, dmE, dmM ), 1000.0 );
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

