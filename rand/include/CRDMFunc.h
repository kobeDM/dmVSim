#ifndef CRDM_FUNC_H
#define CRDM_FUNC_H

#include "inc/shinclude.h"

const double LOS_LIMIT    = 10.0;  // line-of-sight kpc
const double PROTON_MASS  = 0.938; // GeV
const double DM_XS        = 1e-35; // 10^-35 cm^2

const double DM_RHO_SCALE = 0.403; // GeV / cm^3, Publications of the Astronomical Society of Japan 64, 75 (2012).
const double DM_R_SCALE   = 12.53; // kpc,        Publications of the Astronomical Society of Japan 64, 75 (2012).
const double SUN_DISTANCE = 8.0;   // kpc
const double LAMBDA_P     = 0.77;  // GeV

const double EIN_ALPHA    = 0.17;  // arxiv:1103.2377

const double PC2CM        = 3.086*1e21; // kpc->cm

const double V_LIGHT      = 299792.458; // km/s

static std::vector< double > gEneAray;
static std::vector< double > gFlxAray;

bool   readGalprop     ( const String& input );
double getDiffFlux     ( const double& energy );
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
double getDMNFWFlux    ( double theta,
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
double getDMNFWFlux    ( double* x,
                         double* par );
double getDMNFWFluxV   ( double* x,
                         double* par );

double getDMIsoThermalFlux    ( double theta,
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
double getDMIsoThermalFlux    ( double* x,
                                double* par );
double getDMIsoThermalFluxV   ( double* x,
                                double* par );

double getDMEinastoFlux    ( double theta,
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
double getDMEinastoFlux    ( double* x,
                             double* par );
double getDMEinastoFluxV   ( double* x,
                             double* par );

void   printProgressBar( const int& index, const int& total );

#endif // CRDM_FUNC_H
