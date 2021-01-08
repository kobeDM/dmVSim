//////////////////////////////////////////////////////////////////
//
// Physics Utility (based on ROOT application)
//
// 2017/8/4
// Satoshi Higashino
// satoshi.higashino@cern.ch
//
//////////////////////////////////////////////////////////////////
#include "inc/physutil.h"

// ---------------------------------------------------------------
// Calc pseudo continuous b-tag score
//
float ShPUtil::GetPCScore( ShJet* pJet, const bool& isMV2c10 ) {
    if( pJet != nullptr ) return 0.0;
    return CalcPCScore( pJet->score, isMV2c10 );
}

void ShPUtil::SetPCScore( ShJet* pJet, const bool& isMV2c10 ) {
    if( pJet != nullptr ) return;
    pJet->PCscore = CalcPCScore( pJet->score, isMV2c10 );
    return;
}

float ShPUtil::CalcPCScore( const float& score, const bool& isMV2c10 ) {
    // note: working point was set for Rel21 ftag recommendation
       
    // for MV2c10
    if( isMV2c10 == true ) {
        if     ( score > 0.94 ) return 5;
        else if( score > 0.83 ) return 4;
        else if( score > 0.64 ) return 3;
        else if( score > 0.11 ) return 2;
        else                    return 1;
    }
    // for DL1
    else {
        if     ( score > 2.74 ) return 5;
        else if( score > 2.02 ) return 4;
        else if( score > 1.45 ) return 3;
        else if( score > 0.46 ) return 2;
        else                    return 1;
    }

    // it shouldn't be passed!
    return -1.0;
}
    
//////////////////////////////////////////////////////////////////
//
// Kinematics
//
//////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------
// Calc dR
double ShPUtil::GetDR( const TLorentzVector& p1, const TLorentzVector& p2 )
{
    return p1.DeltaR( p2 );
}
double ShPUtil::GetDR( const double& p1Phi, const double& p1Eta,
                       const double& p2Phi, const double& p2Eta )
{
    return ShUtil::GetSumOfSqrt( p1Phi - p2Phi, p1Eta - p2Eta );
}
double ShPUtil::GetDRrap( const TLorentzVector& p1, const TLorentzVector& p2 )
{
    return sqrt( pow( p1.DeltaPhi( p2 ), 2 ) + pow( p1.Rapidity( ) - p2.Rapidity( ), 2 ) );
}

//////////////////////////////////////////////////////////////////
//
// Create Fitting Function
//
//////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------
// Breit-Wigner
// 
TF1* ShPUtil::GetBreitWignerFunction( const String& name )
{
    return new TF1( name.c_str( ), "[0] / ( pow( x*x - [1]*[1], 2 ) + [1]*[1]*[2]*[2] )" );
}

// ---------------------------------------------------------------
// Power Law
// 
TF1* ShPUtil::GetPowLawFunction( const String& name )
{
    return new TF1( name.c_str( ), "[0] * pow(x, [1])" );
}

// ---------------------------------------------------------------
// Exponential
// 
TF1* ShPUtil::GetExpFunction( const String& name )
{
    return new TF1( name.c_str( ), "[0] * exp( [1] * x )" );
}

// ---------------------------------------------------------------
// Exponential of 2nd degree polynomial (ExpPoly2)
// 
TF1* ShPUtil::GetExpPoly2Function( const String& name )
{
    return new TF1( name.c_str( ), "exp( [0] + [1] * x + [2] * x * x )" );
}

// ---------------------------------------------------------------
// Exponential of 2nd degree polynomial (ExpPoly2)
// note: the function is only suited for the ttHyy analysis
// 
TF1* ShPUtil::GetExpPoly2FuncTTHYY( const String& name )
{
    return new TF1( name.c_str( ), "[0] * exp( (x - 100.0)/100.0 * ([1] + [2]*(x - 100.0)/100.0 ) )" );
}

// ---------------------------------------------------------------
// Exponential of 3nd degree polynomial (ExpPoly3)
// 
TF1* ShPUtil::GetExpPoly3Function( const String& name )
{
    return new TF1( name.c_str( ), "exp( [0] + [1] * x + [2] * x * x + [3] * x * x * x )" );
}

// ---------------------------------------------------------------
// Poisson distribution
//
// [0]: normalization factor
// [1]/[2]: mean 
// 
TF1* ShPUtil::GetPoissonFunction( const String& name )
{
    return new TF1( name.c_str( ), "[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1)" );
}

// ---------------------------------------------------------------
// Double sided crystal ball (DoubleCB)
//
// [0] : alphaH
// [1] : nH
// [2] : alphaL
// [3] : nL
// [4] : mu
// [5] : sigma
// [6] : norm
// 
TF1* ShPUtil::GetDoubleCBFunction( const String& name )
{
    return new TF1( name.c_str( ), "(exp( -0.5*[2]*[2] )/pow( ([2]/[3])*(([3]/[2])-[2]+([4]/[5])-(x/[5])), [1] ) * ( (0.5*ROOT::Math::erf(-100*(x - [4]+[2]*[5])) ) + 0.5 ) + exp( -0.5*[0]*[0] )/pow( ([0]/[1])*(([1]/[0])-[0]-([4]/[5])+(x/[5])), [1] ) * ( (0.5*ROOT::Math::erf(100*(x - [4]-[0]*[5])) ) + 0.5 ) + exp( -0.5*pow((x-[4])/[5],2) )*( (0.5*ROOT::Math::erf(100*(x - [4]+[2]*[5])) ) + 0.5 )*( (0.5*ROOT::Math::erf(-100*(x - [4]-[0]*[5])) ) + 0.5 ))*[6]" );
}

double ShPUtil::GetDoubleCBFunction( double x,
                                     double alphaH,
                                     double nH,
                                     double alphaL,
                                     double nL,
                                     double mu,
                                     double sigma,
                                     double norm )
{
    double t = ( x - mu ) / sigma;
    if( t < -alphaL ) {
        double a = exp( -0.5 * alphaL * alphaL );
        double b = nL / alphaL - alphaL;
        return norm * a / TMath::Power( alphaL / nL * (b - t), nL );
    }
    else if( t > alphaH ) {
        double a = exp( -0.5 * alphaH * alphaH );
        double b = nH / alphaH - alphaH;
        return norm * a / TMath::Power( alphaH / nH * ( b + t ), nH );
    }
    return norm * exp( -0.5 * t * t );
}

double ShPUtil::GetDoubleCBFunction( double* x,
                                     double* par )
{
    double var    = x[0];
    double alphaH = par[0];
    double nH     = par[1];
    double alphaL = par[2];
    double nL     = par[3];
    double mu     = par[4];
    double sigma  = par[5];
    double norm   = par[6];
    
    return GetDoubleCBFunction( var, alphaH, nH, alphaL, nL, mu, sigma, norm );
}

double ShPUtil::GetBernstein2( double x,
                               double a0,
                               double a1,
                               double a2 )
{
    double t = 1.0 - x;
    double s = x;
    double result0 = a0 * s * s * TMath::Binomial( 2, 0 );
    double result1 = a1 * t * s * TMath::Binomial( 2, 1 );
    double result2 = a2 * t * t * TMath::Binomial( 2, 2 );

    return result0 + result1 + result2;
}

double ShPUtil::GetBernstein3( double x,
                               double a0,
                               double a1,
                               double a2,
                               double a3 )
{
    double t = 1.0 - x;
    double s = x;
    double result0 = a0 * s * s * s * TMath::Binomial( 3, 0 );
    double result1 = a1 * t * s * s * TMath::Binomial( 3, 1 );
    double result2 = a2 * t * t * s * TMath::Binomial( 3, 2 );
    double result3 = a3 * t * t * t * TMath::Binomial( 3, 3 );

    return result0 + result1 + result2 + result3; 
}

double ShPUtil::GetBernstein4( double x,
                               double a0,
                               double a1,
                               double a2,
                               double a3,
                               double a4,
                               double norm )
{
    double t = 1.0 - x;
    double s = x;
    double result0 = a0 * s * s * s * s * TMath::Binomial( 4, 0 );
    double result1 = a1 * t * s * s * s * TMath::Binomial( 4, 1 );
    double result2 = a2 * t * t * s * s * TMath::Binomial( 4, 2 );
    double result3 = a3 * t * t * t * s * TMath::Binomial( 4, 3 );
    double result4 = a4 * t * t * t * t * TMath::Binomial( 4, 4 );

    return result0 + result1 + result2 + result3 + result4;
}

double ShPUtil::GetBernstein5( double x,
                               double a0,
                               double a1,
                               double a2,
                               double a3,
                               double a4,
                               double a5,
                               double norm )
{
    double t = 1.0 - x;
    double s = x;
    double result0 = a0 * s * s * s * s * s * TMath::Binomial( 5, 0 );
    double result1 = a1 * t * s * s * s * s * TMath::Binomial( 5, 1 );
    double result2 = a2 * t * t * s * s * s * TMath::Binomial( 5, 2 );
    double result3 = a3 * t * t * t * s * s * TMath::Binomial( 5, 3 );
    double result4 = a4 * t * t * t * t * s * TMath::Binomial( 5, 4 );
    double result5 = a5 * t * t * t * t * t * TMath::Binomial( 5, 5 );

    return result0 + result1 + result2 + result3 + result4 + result5;
}

// ---------------------------------------------------------------
// Gaussian 
//
double ShPUtil::GetGauss( double t,
                          double norm,
                          double mean,
                          double sigma )
{
    return norm * TMath::Gaus( t, mean, sigma, false );
}

// ---------------------------------------------------------------
// Double Gaussian (common mean value)
//
double ShPUtil::GetDGauss( double t,
                           double norm1,
                           double mean,
                           double sigma1,
                           double norm2,
                           double sigma2 )
{
    return norm1 * TMath::Gaus( t, mean, sigma1, false ) + norm2 * TMath::Gaus( t, mean, sigma2, false );
}

// ---------------------------------------------------------------
// Gaussian (normalized)
//
double ShPUtil::GetGaussN( double t,
                           double mean,
                           double sigma )
{
    return TMath::Gaus( t, mean, sigma, true );
}

// ---------------------------------------------------------------
// Double Gaussian (normalized, common mean value)
//
double ShPUtil::GetDGaussN( double t,
                            double fraction,
                            double mean,
                            double sigma1,
                            double sigma2 )
{
    return fraction * TMath::Gaus( t, mean, sigma1, true ) + ( 1.0 - fraction ) * TMath::Gaus( t, mean, sigma2, true );
}

    
// ---------------------------------------------------------------
// Flavor mixing probabirity 
//
// note: purity should be 1/sqrt(2) when decaying to a CP eigenstate
// 
double ShPUtil::GetMixingCPV( double t,
                              double norm,
                              double tau,
                              double dm,
                              double sin2phi1,
                              double cp,
                              double flv,
                              double purity )
{
    double decay = exp(-1.0 * fabs(t) / tau );

    double cosTerm = 2.0 * purity * purity - 1.0;
    double sinTerm = 2.0 * purity * sqrt( 1.0 - purity * purity ) * cp * flv * sin2phi1;

    return norm * decay * ( 1.0 + cosTerm * cos( dm * t ) + sinTerm * sin( dm * t ) );
}

double ShPUtil::GetMixingCPV( double* x,
                              double* par )
{
    double t        = x[0];
    double norm     = par[0];
    double tau      = par[1];
    double dm       = par[2];
    double sin2phi1 = par[3];
    double cp       = par[4];
    double flv      = par[5];
    double purity   = par[6];

    return ShPUtil::GetMixingCPV( t, norm, tau, dm, sin2phi1, cp, flv, purity );
}

double ShPUtil::GetMixingCPVPhase( double t,
                                   double norm,
                                   double tau,
                                   double dm,
                                   double cpPhase,
                                   double selPhase,
                                   double cp,
                                   double flv,
                                   double purity )
{
    double decay = exp(-1.0 * fabs(t) / tau );

    // double sinPhase = sin( selPhase - flv*cpPhase );
    // double sinPhase = 2.0 * 0.096 / 4.78 * sin( selPhase - flv*cpPhase );
    double sinPhase = 0.1 * sin( selPhase - flv*cpPhase );
    double cosTerm = 2.0 * purity * purity - 1.0;
    double sinTerm = 2.0 * purity * sqrt( 1.0 - purity * purity ) * cp * sinPhase;

    double normUnit = tau * (1.0 + cosTerm * 1.0 / ( 1.0 + dm*dm*tau*tau ) );
        
    return norm * decay * ( 1.0 + cosTerm * cos( dm * t ) + sinTerm * sin( dm * t ) ) / normUnit;
}

double ShPUtil::GetMixingCPVPhase( double* x,
                                   double* par )
{
    double t        = x[0];
    double norm     = par[0];
    double tau      = par[1];
    double dm       = par[2];
    double cpPhase  = par[3];
    double selPhase = par[4];
    double cp       = par[5];
    double flv      = par[6];
    double purity   = par[7];

    return ShPUtil::GetMixingCPVPhase( t, norm, tau, dm, cpPhase, selPhase, cp, flv, purity );
}



//////////////////////////////////////////////////////////////////
//
// Statistics
//
//////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------
// Get Binomial Error from efficiency and sample size
// 
// return: Binomial error
double ShPUtil::GetBinomialError( const double& eff, const int& num )
{
    return sqrt( eff * ( 1 - eff ) / (double)num );
}

// ---------------------------------------------------------------
// Get Binomial Error from numerator and denominator
// 
// return: Binomial error
double ShPUtil::GetBinomialError( const double& numerator, const double& denominator, const int& num )
{
    if( ShUtil::IsZero( denominator ) == true ) return -1.0;
    double eff = numerator / denominator;

    return GetBinomialError( eff, num );
}

// ---------------------------------------------------------------
// Get Simple significance
// Z = S / sqrt( S + B )
// 
double ShPUtil::SimpleSignificance( const double& signal, const double& bg )
{
    if( ShUtil::IsZero( signal + bg ) == true ) return 0.0;
    return signal / sqrt( signal + bg );    
}

// ---------------------------------------------------------------
// Get Simple significance
// Z = sqrt( 2 * ( (S + B) * log( 1 + S / B ) - S ) )
// 
double ShPUtil::Significance( const double& signal, const double& bg )
{
    if( ShUtil::IsZero( bg ) == true ) return 0.0;
    return sqrt( 2 * ( (signal + bg) * log( 1 + signal/bg ) - signal ) );
}

// ---------------------------------------------------------------
// Get Simple significance error
// Z = S / sqrt( S + B )
// 
double ShPUtil::SimpleSignificanceError( const double& signal,    const double& bg,
                                         const double& signalErr, const double& bgErr )
{
    if( ShUtil::IsZero( signal + bg ) == true ) return 0.0;
    double z0 = SimpleSignificance( signal, bg );

    double paramA = z0 / signal;
    double paramB = 0.5 * z0 * paramA;

    return paramA * sqrt( ( 1.0 - paramB ) * ( 1.0 - paramB ) * signalErr * signalErr + paramB * paramB * bgErr * bgErr );
}

// ---------------------------------------------------------------
// Get significance error
// Z = sqrt( 2 * ( (S + B) * log( 1 + S / B ) - S ) )
// 
double ShPUtil::SignificanceError( const double& signal,    const double& bg,
                                   const double& signalErr, const double& bgErr )
{
    if( ShUtil::IsZero( bg ) == true ) return 0.0;
    double z0 = Significance( signal, bg );
    double paramA = signal / bg;
    double paramB = log( 1 + paramA );
    double paramC = paramB * paramB * signalErr * signalErr;
    double paramD = ( paramB - paramA ) * ( paramB - paramA ) * bgErr * bgErr;

    return sqrt( paramC + paramD ) / z0;
}
