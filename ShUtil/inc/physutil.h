//////////////////////////////////////////////////////////////////
//
// Physics Utility (based on ROOT application)
//
// 2017/8/4
// Satoshi Higashino
// satoshi.higashino@cern.ch
//
//////////////////////////////////////////////////////////////////
#ifndef SH_PHYS_UTIL_H
#define SH_PHYS_UTIL_H

// COMMON
#include "rootutil.h"

// ROOT
#include "TLorentzVector.h"

class ShPUtil
{
public:

    //////////////////////////////////////////////////////////////////
    //
    // Jet
    //
    //////////////////////////////////////////////////////////////////
    // ---------------------------------------------------------------
    // Define ShJet instance
    struct ShJet
    {
        float pt;
        float phi;
        float eta;
        float e;
        float score;
        float PCscore;
    };

    // ---------------------------------------------------------------
    // Calc pseudo continuous b-tag score
    //
    static float GetPCScore( ShJet* pJet, const bool& isMV2c10 = true );
    static void  SetPCScore( ShJet* pJet, const bool& isMV2c10 = true );

    static float CalcPCScore( const float& score, const bool& isMV2c10 = true );
    
    //////////////////////////////////////////////////////////////////
    //
    // Kinematics
    //
    //////////////////////////////////////////////////////////////////
    // ---------------------------------------------------------------
    // Calc dR
    static double GetDR( const TLorentzVector& p1, const TLorentzVector& p2 );
    static double GetDR( const double& p1Phi,      const double& p1Eta,
                         const double& p2Phi,      const double& p2Eta );
    static double GetDRrap( const TLorentzVector& p1, const TLorentzVector& p2 );

    //////////////////////////////////////////////////////////////////
    //
    // Create Fitting Function
    //
    //////////////////////////////////////////////////////////////////
    // ---------------------------------------------------------------
    // Breit-Wigner
    // 
    static TF1* GetBreitWignerFunction( const String& name );

    // ---------------------------------------------------------------
    // Power Law
    // 
    static TF1* GetPowLawFunction( const String& name );

    // ---------------------------------------------------------------
    // Exponential
    // 
    static TF1* GetExpFunction( const String& name );

    // ---------------------------------------------------------------
    // Exponential of 2nd degree polynomial (ExpPoly2)
    // 
    static TF1* GetExpPoly2Function( const String& name );

    // ---------------------------------------------------------------
    // Exponential of 2nd degree polynomial (ExpPoly2)
    // note: the function is only suited for the ttHyy analysis
    // 
    static TF1* GetExpPoly2FuncTTHYY( const String& name );

    // ---------------------------------------------------------------
    // Exponential of 3nd degree polynomial (ExpPoly3)
    // 
    static TF1* GetExpPoly3Function( const String& name );

    // ---------------------------------------------------------------
    // Poisson distribution
    //
    // [0]: normalization factor
    // [1]/[2]: mean 
    // 
    static TF1* GetPoissonFunction( const String& name );

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
    static TF1* GetDoubleCBFunction( const String& name );

    static double GetDoubleCBFunction( double x,
                                       double alphaH,
                                       double nH,
                                       double alphaL,
                                       double nL,
                                       double mu,
                                       double sigma,
                                       double norm );

    static double GetDoubleCBFunction( double* x, double* par );
    
    static double GetBernstein2( double x,
                                 double a0,
                                 double a1,
                                 double a2 );
    static double GetBernstein3( double x,
                                 double a0,
                                 double a1,
                                 double a2,
                                 double a3 );
    static double GetBernstein4( double x,
                                 double a0,
                                 double a1,
                                 double a2,
                                 double a3,
                                 double a4,
                                 double norm );
    static double GetBernstein5( double x,
                                 double a0,
                                 double a1,
                                 double a2,
                                 double a3,
                                 double a4,
                                 double a5,
                                 double norm );

    // ---------------------------------------------------------------
    // Gaussian 
    //
    static double GetGauss( double t,
                            double norm,
                            double mean,
                            double sigma );

    // ---------------------------------------------------------------
    // Double Gaussian (common mean value)
    //
    static double GetDGauss( double t,
                             double norm1,
                             double mean,
                             double sigma1,
                             double norm2,
                             double sigma2 );

    // ---------------------------------------------------------------
    // Gaussian (normalized)
    //
    static double GetGaussN( double t,
                             double mean,
                             double sigma );

    // ---------------------------------------------------------------
    // Double Gaussian (normalized, common mean value)
    //
    static double GetDGaussN( double t,
                              double fraction,
                              double mean,
                              double sigma1,
                              double sigma2 );

    
    // ---------------------------------------------------------------
    // Flavor mixing probabirity 
    //
    // note: purity should be 1/sqrt(2) when decaying into a CP eigenstate
    // 
    static double GetMixingCPV( double t,
                                double norm,
                                double tau,
                                double dm,
                                double sin2phi1,
                                double cp,
                                double flv,
                                double purity );

    static double GetMixingCPV( double* x,
                                double* par );

    static double GetMixingCPVPhase( double t,
                                     double norm,
                                     double tau,
                                     double dm,
                                     double cpPhase,
                                     double selPhase,
                                     double cp,
                                     double flv,
                                     double purity );

    static double GetMixingCPVPhase( double* x,
                                     double* par );

    //////////////////////////////////////////////////////////////////
    //
    // Statistics
    //
    //////////////////////////////////////////////////////////////////
    // ---------------------------------------------------------------
    // Get Binomial Error from efficiency and sample size
    // 
    // return: Binomial error
    static double GetBinomialError( const double& eff, const int& num );

    // ---------------------------------------------------------------
    // Get Binomial Error from numerator and denominator
    // 
    // return: Binomial error
    static double GetBinomialError( const double& numerator, const double& denominator, const int& num );

    // ---------------------------------------------------------------
    // Get Simple significance
    // Z = S / sqrt( S + B )
    // 
    static double SimpleSignificance( const double& signal, const double& bg );

    // ---------------------------------------------------------------
    // Get Simple significance
    // Z = sqrt( 2 * ( (S + B) * log( 1 + S / B ) - S ) )
    // 
    static double Significance( const double& signal, const double& bg );

    // ---------------------------------------------------------------
    // Get Simple significance error
    // Z = S / sqrt( S + B )
    // 
    static double SimpleSignificanceError( const double& signal, const double& bg, const double& signalErr, const double& bgErr );

    // ---------------------------------------------------------------
    // Get significance error
    // Z = sqrt( 2 * ( (S + B) * log( 1 + S / B ) - S ) )
    // 
    static double SignificanceError( const double& signal, const double& bg, const double& signalErr, const double& bgErr );

};

#endif // SH_PHYS_UTIL_H
