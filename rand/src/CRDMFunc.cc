#include "CRDMFunc.h"

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
        gEneAray.push_back( energy * 0.001);  // GeV
        gFlxAray.push_back( flux / energy / energy * 1000.0 ); // /cm^2/s/sr/GeV
    }

    return true;
}

double getDiffFlux( const double& energy )
{
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

    return (flxMaxBin - flxMinBin) / (eneMaxBin - eneMinBin) * (energy - eneMinBin) + flxMinBin; // /cm^2/s/sr/GeV
}


double getTMax( double pE,
                double pM,
                double dmM )
{
    if( pE < 0.0 || pM < 0.0 || dmM < 0.0 ) return 0.0;

    double num = pE*pE + 2*pM*pE; // GeV^2
    double den = pE + ((pM+dmM)*(pM+dmM)/(2.0*dmM)); // GeV
    return num / den; // GeV
}

double getTMaxInv( double pE,
                   double pM,
                   double dmM )
{
    if( pE < 0.0 || pM < 0.0 || dmM < 0.0 ) return 0.0;

    return getDiffFlux( pE ) / getTMax( pE, pM, dmM ); // /cm^2/s/sr/GeV^2
}

double getTMaxInv( double* x,
                   double* par )
{
    double pE  = x[0];
    double pM  = par[0];
    double dmM = par[1];
    return getTMaxInv( pE, pM, dmM ); // /cm^2/s/sr/GeV^2
}

double getTMin( double pM,
                double dmE,
                double dmM )
{
    if( pM < 0.0 || dmE < 0.0 || dmM < 0.0 ) return 0.0;

    double first = dmE * 0.5 - pM; // GeV
    double tmp = 1.0 + (2.0*dmE/dmM)*pow((pM+dmM)/(2.0*pM-dmE) , 2 ); 
    
    double second = 0.0;
    if( first > 0.0 ) second = 1.0 + sqrt(tmp);
    else              second = 1.0 - sqrt(tmp);

    return first * second;  // GeV
}

double getTIntegral( double pM,
                     double dmE,
                     double dmM )
{
    if( pM < 0.0 || dmE < 0.0 || dmM < 0.0 ) return 0.0;
    
    TF1 func( "TInt", getTMaxInv, 1e-8, 1e+12, 2, 1 );
    func.SetParameter( 0, pM );
    func.SetParameter( 1, dmM );

    return func.Integral( getTMin( pM, dmE, dmM ), 1e+12, 1e-30 ); // /cm^2/s/sr/GeV
}


// energy component
double getDMFluxEn( double pM,
                    double dmE,
                    double dmM,
                    double dmXS,
                    double lambdaP )
{
    if( pM < 0.0 || dmE < 0.0 || dmM < 0.0 || dmXS < 0.0 ) return 0.0;

    double TIntegral = getTIntegral( pM, dmE, dmM ); // /cm^2/s/sr/GeV
    double formFactor = 1.0 / pow( 1.0 + (2.0 * dmM * dmE) / (lambdaP*lambdaP), 2 );

    return TIntegral * dmXS / dmM * formFactor * formFactor; // /s/sr/GeV2
}

double getDMFluxEn( double* x,
                    double* par )
{
    double velo     = x[0];
    double pM       = par[0];
    double dmM      = par[1];
    double dmXS     = par[2];
    double lambdaP  = par[3];

    double gamma    = 1.0 / sqrt( 1 - pow( velo / V_LIGHT , 2) );
    double dmE      = dmM * (gamma - 1.0);

    double flux = getDMFluxEn( pM, dmE, dmM, dmXS, lambdaP ); // /s/sr/GeV2
    return flux * dmM * dmM * gamma * gamma * gamma; // Note!!: the local dark matter density is ignored at this stage... 
}



/////////////////////////////////////////////////////////////////////////
// NFW profile
/////////////////////////////////////////////////////////////////////////
double getDMNFWFluxDir( double theta,
                        double phi,
                        double los,      // kpc
                        double dmDScale, // GeV/cm^3
                        double dmRScale, // kpc
                        double sunDist ) // kpc
{
    if( dmDScale < 0.0 || dmRScale < 0.0 || sunDist < 0.0 ) return -100.0;

    double cosTheta = cos( theta );
    double cosPhi   = cos( phi );
    double rsq      = los*los + sunDist*sunDist - 2*los*sunDist*cosTheta*cosPhi; // kpc^2
    double r        = sqrt( fabs( rsq ) ); // kpc
    if( r < 0.1 ) r = 0.1; // r constraint

    double relR     = r / dmRScale;
    double rhoNFW   = dmDScale / (relR*pow( 1 + relR, 2 )); // GeV/cm^3

    double jacobian = los*los; // Note: This jacobian should only affect the los. Do not add cosTheta!!!
    double fluxCorr = 1.0 / (4.0*TMath::Pi()*los*los);

    return fluxCorr * jacobian * rhoNFW * PC2CM * PC2CM * PC2CM; // GeV/kpc^3
}

double getDMNFWFluxDir( double* x,
                        double* par )
{
    double theta    = x[0];
    double phi      = x[1];
    double los      = par[0];
    double dmDScale = par[1];
    double dmRScale = par[2];
    double sunDist  = par[3];

    return getDMNFWFluxDir( theta, phi, los, dmDScale, dmRScale, sunDist );
}

double getDMNFWFluxDirLOS( double* x,
                           double* par )
{
    double los      = x[0];
    double theta    = par[0];
    double phi      = par[1];
    double dmDScale = par[2];
    double dmRScale = par[3];
    double sunDist  = par[4];

    return getDMNFWFluxDir( theta, phi, los, dmDScale, dmRScale, sunDist );
}

double getDMNFWFluxDirInt(  double theta,
                            double phi,
                            double los,
                            double dmDScale,
                            double dmRScale,
                            double sunDist )
{
    TF1 func( "NFWDirInt", getDMNFWFluxDirLOS, 0.0, los, 5, 1 );
    func.SetParameter( 0, theta    );
    func.SetParameter( 1, phi      );
    func.SetParameter( 2, dmDScale );
    func.SetParameter( 3, dmRScale );
    func.SetParameter( 4, sunDist  );

    return func.Integral( 0.0, los ) / PC2CM / PC2CM; // GeV/kpc^2 -> GeV/cm^2
}


double getDMNFWFluxDirInt(  double* x,
                            double* par )
{
    double theta    = x[0];
    double phi      = x[1];
    double los      = par[0];
    double dmDScale = par[1];
    double dmRScale = par[2];
    double sunDist  = par[3];

    return getDMNFWFluxDirInt( theta, phi, los, dmDScale, dmRScale, sunDist );
}

double getDMNFWFlux( double theta,
                     double phi,
                     double los,
                     double pM,
                     double dmE,
                     double dmM,
                     double dmXS,
                     double dmDScale, // GeV/cm^3
                     double dmRScale, // kpc
                     double sunDist,
                     double lambdaP )
{
    if( pM < 0.0 || dmE < 0.0 || dmM < 0.0 || dmXS < 0.0 || dmDScale < 0.0 || dmRScale < 0.0 || sunDist < 0.0 )
        return -100.0;
    
    double fluxDir = getDMNFWFluxDirInt( theta, phi, los, dmDScale, dmRScale, sunDist ); // GeV/cm^2
    double fluxEn  = getDMFluxEn( pM, dmE, dmM, dmXS, lambdaP ); // /s/sr/GeV^2

    if( fluxDir < 0 ) DEBUG("errorDir");
    if( fluxEn < 0 ) DEBUG("errorEn");

    return fluxDir * fluxEn; // /cm^2/GeV/s/sr
}

double getDMNFWFlux( double* x,
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

    return getDMNFWFlux( theta, phi, los, pM, dmE, dmM, dmXS, dmDScale, dmRScale, sunDist, lambdaP );
}

double getDMNFWFluxInt( double* x,
                        double* par )
{
    double theta    = x[0];
    double phi      = x[1];
    double dmE      = x[2];
    double los      = par[0];
    double pM       = par[1];
    double dmM      = par[2];
    double dmXS     = par[3];
    double dmDScale = par[4];
    double dmRScale = par[5];
    double sunDist  = par[6];
    double lambdaP  = par[7];

    return getDMNFWFlux( theta, phi, los, pM, dmE, dmM, dmXS, dmDScale, dmRScale, sunDist, lambdaP );
}

double getDMNFWFluxV( double* x,
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

    double flux = getDMNFWFlux( theta, phi, los, pM, dmE, dmM, dmXS, dmDScale, dmRScale, sunDist, lambdaP );
    return flux * dmM * dmM * gamma * gamma * gamma; // Note!!: the local dark matter density is ignored at this stage... 
}

/////////////////////////////////////////////////////////////////////////
// Iso-thermal profile
/////////////////////////////////////////////////////////////////////////
double getDMIsoThermalFluxDir( double theta,
                               double phi,
                               double los,
                               double dmDScale,
                               double dmRScale,
                               double sunDist )
{
    if( dmDScale < 0.0 || dmRScale < 0.0 || sunDist < 0.0 ) return -100.0;

    double cosTheta = cos( theta );
    double cosPhi   = cos( phi );
    double rsq      = los*los + sunDist*sunDist - 2*los*sunDist*cosTheta*cosPhi; // kpc^2
    double r        = sqrt( fabs( rsq ) ); // kpc

    double relR     = r / dmRScale;
    double rhoIsoThermal   = dmDScale / ( 1.0 + pow( relR, 2 ) ); // GeV/cm^3
    double jacobian = los*los; // Note: This jacobian should only affect the los. Do not add cosTheta!!!
    double fluxCorr = 1.0 / (4.0*TMath::Pi()*los*los);

    return fluxCorr * jacobian * rhoIsoThermal * PC2CM * PC2CM * PC2CM; // GeV/kpc^3
}

double getDMIsoThermalFluxDir( double* x,
                               double* par )
{
    double theta    = x[0];
    double phi      = x[1];
    double los      = par[0];
    double dmDScale = par[1];
    double dmRScale = par[2];
    double sunDist  = par[3];

    return getDMIsoThermalFluxDir( theta, phi, los, dmDScale, dmRScale, sunDist );
}


double getDMIsoThermalFluxDirLOS( double* x,
                                  double* par )
{
    double los      = x[0];
    double theta    = par[0];
    double phi      = par[1];
    double dmDScale = par[2];
    double dmRScale = par[3];
    double sunDist  = par[4];

    return getDMIsoThermalFluxDir( theta, phi, los, dmDScale, dmRScale, sunDist );
}


double getDMIsoThermalFluxDirInt( double theta,
                                  double phi,
                                  double los,
                                  double dmDScale,
                                  double dmRScale,
                                  double sunDist )
{
    TF1 func( "NFWDirInt", getDMIsoThermalFluxDirLOS, 0.0, los, 5, 1 );
    func.SetParameter( 0, theta    );
    func.SetParameter( 1, phi      );
    func.SetParameter( 2, dmDScale );
    func.SetParameter( 3, dmRScale );
    func.SetParameter( 4, sunDist  );

    return func.Integral( 0.0, los ) / PC2CM / PC2CM; // GeV/kpc^2 -> GeV/cm^2
}

double getDMIsoThermalFluxDirInt(  double* x,
                                   double* par )
{
    double theta    = x[0];
    double phi      = x[1];
    double los      = par[0];
    double dmDScale = par[1];
    double dmRScale = par[2];
    double sunDist  = par[3];

    return getDMIsoThermalFluxDirInt( theta, phi, los, dmDScale, dmRScale, sunDist );
}


double getDMIsoThermalFlux( double theta,
                            double phi,
                            double los,
                            double pM,
                            double dmE,
                            double dmM,
                            double dmXS,
                            double dmDScale, // GeV/cm^3
                            double dmRScale, // kpc
                            double sunDist,
                            double lambdaP )
{
    if( pM < 0.0 || dmE < 0.0 || dmM < 0.0 || dmXS < 0.0 || dmDScale < 0.0 || dmRScale < 0.0 || sunDist < 0.0 )
        return -100.0;

    double fluxDir = getDMIsoThermalFluxDirInt( theta, phi, los, dmDScale, dmRScale, sunDist );
    double fluxEn  = getDMFluxEn( pM, dmE, dmM, dmXS, lambdaP );

    return fluxDir * fluxEn; // /cm^2/GeV/s/sr
}

double getDMIsoThermalFlux( double* x,
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

    return getDMIsoThermalFlux( theta, phi, los, pM, dmE, dmM, dmXS, dmDScale, dmRScale, sunDist, lambdaP );
}


double getDMIsoThermalFluxInt( double* x,
                               double* par )
{
    double theta    = x[0];
    double phi      = x[1];
    double dmE      = x[2];
    double los      = par[0];
    double pM       = par[1];
    double dmM      = par[2];
    double dmXS     = par[3];
    double dmDScale = par[4];
    double dmRScale = par[5];
    double sunDist  = par[6];
    double lambdaP  = par[7];

    return getDMIsoThermalFlux( theta, phi, los, pM, dmE, dmM, dmXS, dmDScale, dmRScale, sunDist, lambdaP );
}

double getDMIsoThermalFluxV( double* x,
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

    double flux = getDMIsoThermalFlux( theta, phi, los, pM, dmE, dmM, dmXS, dmDScale, dmRScale, sunDist, lambdaP );
    return flux * dmM * dmM * gamma * gamma * gamma; // Note!!: the local dark matter density is ignored at this stage... 
}



/////////////////////////////////////////////////////////////////////////
// Einasto profile
/////////////////////////////////////////////////////////////////////////
double getDMEinastoFluxDir( double theta,
                            double phi,
                            double los,
                            double dmDScale,
                            double dmRScale,
                            double sunDist )
{
    if( dmDScale < 0.0 || dmRScale < 0.0 || sunDist < 0.0 ) return -100.0;

    double cosTheta = cos( theta );
    double cosPhi   = cos( phi );
    double rsq      = los*los + sunDist*sunDist - 2*los*sunDist*cosTheta*cosPhi;
    double r        = sqrt( fabs( rsq ) );

    double relR     = r / dmRScale;
    double rhoEinasto   = dmDScale * exp( -2.0 * ( pow( relR, EIN_ALPHA ) - 1.0 ) / EIN_ALPHA ); // GeV/cm^3
    double jacobian = los*los; // Note: This jacobian should only affect the los. Do not add cosTheta!!!
    double fluxCorr = 1.0 / (4.0*TMath::Pi()*los*los);

    return fluxCorr * jacobian * rhoEinasto * PC2CM * PC2CM * PC2CM; // GeV/kpc^3
}

double getDMEinastoFluxDir( double* x,
                            double* par )
{
    double theta    = x[0];
    double phi      = x[1];
    double los      = par[0];
    double dmDScale = par[1];
    double dmRScale = par[2];
    double sunDist  = par[3];

    return getDMEinastoFluxDir( theta, phi, los, dmDScale, dmRScale, sunDist );
}

double getDMEinastoFluxDirLOS( double* x,
                               double* par )
{
    double los      = x[0];
    double theta    = par[0];
    double phi      = par[1];
    double dmDScale = par[2];
    double dmRScale = par[3];
    double sunDist  = par[4];

    return getDMEinastoFluxDir( theta, phi, los, dmDScale, dmRScale, sunDist );
}

double getDMEinastoFluxDirInt(  double theta,
                                double phi,
                                double los,
                                double dmDScale,
                                double dmRScale,
                                double sunDist )
{
    TF1 func( "NFWDirInt", getDMEinastoFluxDirLOS, 0.0, los, 5, 1 );
    func.SetParameter( 0, theta    );
    func.SetParameter( 1, phi      );
    func.SetParameter( 2, dmDScale );
    func.SetParameter( 3, dmRScale );
    func.SetParameter( 4, sunDist  );

    return func.Integral( 0.0, los ) / PC2CM / PC2CM; // GeV/kpc^2 -> GeV/cm^2
}


double getDMEinastoFluxDirInt(  double* x,
                                double* par )
{
    double theta    = x[0];
    double phi      = x[1];
    double los      = par[0];
    double dmDScale = par[1];
    double dmRScale = par[2];
    double sunDist  = par[3];

    return getDMEinastoFluxDirInt( theta, phi, los, dmDScale, dmRScale, sunDist );
}

double getDMEinastoFlux( double theta,
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

    double fluxDir = getDMEinastoFluxDirInt( theta, phi, los, dmDScale, dmRScale, sunDist );
    double fluxEn  = getDMFluxEn( pM, dmE, dmM, dmXS, lambdaP );

    return fluxDir * fluxEn; // /cm^2/GeV/s/sr
}

double getDMEinastoFlux( double* x,
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

    return getDMEinastoFlux( theta, phi, los, pM, dmE, dmM, dmXS, dmDScale, dmRScale, sunDist, lambdaP );
}

double getDMEinastoFluxInt( double* x,
                            double* par )
{
    double theta    = x[0];
    double phi      = x[1];
    double dmE      = x[2];
    double los      = par[0];
    double pM       = par[1];
    double dmM      = par[2];
    double dmXS     = par[3];
    double dmDScale = par[4];
    double dmRScale = par[5];
    double sunDist  = par[6];
    double lambdaP  = par[7];

    return getDMEinastoFlux( theta, phi, los, pM, dmE, dmM, dmXS, dmDScale, dmRScale, sunDist, lambdaP );
}


double getDMEinastoFluxV( double* x,
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

    double flux = getDMEinastoFlux( theta, phi, los, pM, dmE, dmM, dmXS, dmDScale, dmRScale, sunDist, lambdaP );
    return flux * dmM * dmM * gamma * gamma * gamma; // Note!!: the local dark matter density is ignored at this stage... 
}




void   printProgressBar( const int& index, const int& total )
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

double getVelocity( double* xCDF, double* yCDF, const double& rndUni )
{
    if( xCDF == nullptr || yCDF == nullptr ) return -1.0;
    
    double retVal = 0.0;
    double x1 = 0.0, x2 = 0.0;
    double y1 = 0.0, y2 = 0.0;
    if( rndUni < yCDF[0]    ) return retVal;
    if( rndUni > yCDF[99999] ) return xCDF[99999];

    for( int i = 0; i < 100000; ++i ) {
        if( rndUni < yCDF[i] ) {
            // linear compensation
            x1 = xCDF[i-1];
            x2 = xCDF[i];
            y1 = yCDF[i-1];
            y2 = yCDF[i];

            retVal = x1 + (rndUni - y1) * (x2 - x1) / (y2 - y1);
            break;
        }
    }
    
    return retVal;
}
