//////////////////////////////////////////////////////////////////
//
// Efficiency plot maker
//
// 2018/12/24
// Satoshi Higashino
// satoshi.higashino@cern.ch
//
//////////////////////////////////////////////////////////////////
#ifndef SH_EFF_MAKER_H
#define SH_EFF_MAKER_H

// COMMON
#include "inc/shinclude.h"

// ---------------------------------------------------------------
// EffMake (1D)
// ---------------------------------------------------------------
class EffMaker
{
private:
    String m_name;
    double m_rangeMin;
    double m_rangeMax;
    int    m_nBin;

    std::vector< double > m_denArr;
    std::vector< double > m_numArr;
    
public:
    EffMaker( const String& name, const int& nBin, const double& rangeMin, const double& rangeMax );
    ~EffMaker( );

    String                 getName           ( ) { return m_name;     }
    double                 getRangeMin       ( ) { return m_rangeMin; }
    double                 getRangeMax       ( ) { return m_rangeMax; }
    int                    getNBin           ( ) { return m_nBin;     }
    
    std::vector< double >* getDenominator    ( );
    void                   getDenominator    ( std::vector< double >* pDenArr );
    std::vector< double >* getNumerator      ( );
    void                   getNumerator      ( std::vector< double >* pNumArr );

    bool                   setDenominator    ( const int&             bin,
                                               const double&          val );
    bool                   setNumerator      ( const int&             bin,
                                               const double&          val );

    bool                   countDenEvt       ( const double&          val );
    bool                   countNumEvt       ( const double&          val );

    int                    findBin           ( const double&          val );
    double                 findValue         ( const int&             bin );

    double                 getEff            ( const int&             bin );
    double                 getEff            ( const double&          val );
    double                 getEffErr         ( const int&             bin );
    double                 getEffErr         ( const double&          val );
    
    TH1F*                  createHist        ( const String&          name = "" );
    TH1F*                  createVarBinHist  ( const int&             nBin,
                                               const double*          pBins,
                                               const String&          name = "" );
                                               
};


// ---------------------------------------------------------------
// EffMake2D
//
// note: 
//             ...          ...         ...  
//             Y5           Y5          Y5   
//             Y4           Y4          Y4   
//        X1 { Y3 },   X2 { Y3 },  X3 { Y3 },
//             Y2           Y2          Y2   
//             Y1           Y1          Y1   
//
// ---------------------------------------------------------------
class EffMaker2D
{
private:
    String m_name;
    double m_rangeMinX;
    double m_rangeMaxX;
    int    m_nBinX;
    double m_rangeMinY;
    double m_rangeMaxY;
    int    m_nBinY;

    std::vector< std::vector< double > > m_denArr;
    std::vector< std::vector< double > > m_numArr;
    
public:
    EffMaker2D( const String& name,
                const int&    nBinX, const double& rangeMinX, const double& rangeMaxX,
                const int&    nBinY, const double& rangeMinY, const double& rangeMaxY );

    ~EffMaker2D( );

    String                 getName           ( ) { return m_name;      }
    double                 getRangeMinX      ( ) { return m_rangeMinX; }
    double                 getRangeMaxX      ( ) { return m_rangeMaxX; }
    int                    getNBinX          ( ) { return m_nBinX;     }
    double                 getRangeMinY      ( ) { return m_rangeMinY; }
    double                 getRangeMaxY      ( ) { return m_rangeMaxY; }
    int                    getNBinY          ( ) { return m_nBinY;     }
    
    std::vector< std::vector< double > >* getDenominator2D  ( );
    void                                  getDenominator2D  ( std::vector< std::vector< double > >* pDenArr );


    std::vector< double >*                getDenominatorX   ( const double&          val );
    void                                  getDenominatorX   ( std::vector< double >* pDenArr,
                                                              const double&          val);
    std::vector< double >*                getDenominatorY   ( const double&          val );
    void                                  getDenominatorY   ( std::vector< double >* pDenArr,
                                                              const double&          val);

    std::vector< std::vector< double > >* getNumerator2D    ( );
    void                                  getNumerator2D    ( std::vector< std::vector< double > >* pDenArr );

    std::vector< double >*                getNumeratorX     ( const double&          val );
    void                                  getNumeratorX     ( std::vector< double >* pNumArr,
                                                              const double&          val);
    std::vector< double >*                getNumeratorY     ( const double&          val );
    void                                  getNumeratorY     ( std::vector< double >* pNumArr,
                                                              const double&          val);

    bool                   setDenominator    ( const int&             binX,
                                               const int&             binY,
                                               const double&          val );

    bool                   setNumerator      ( const int&             binX,
                                               const int&             binY,
                                               const double&          val );

    bool                   countDenEvt       ( const double&          valX,
                                               const double&          valY );
    bool                   countNumEvt       ( const double&          valX,
                                               const double&          valY );

    int                    findBinX          ( const double&          val );
    int                    findBinY          ( const double&          val );
    double                 findValueX        ( const int&             bin );
    double                 findValueY        ( const int&             bin );

    double                 getEff            ( const int&             binX,
                                               const int&             binY );
    double                 getEff            ( const double&          valX,
                                               const double&          valY );
    double                 getEffErr         ( const int&             binX,
                                               const int&             binY );
    double                 getEffErr         ( const double&          valX,
                                               const double&          valY );
    
    TH1F*                  createHistX       ( const String&          name = "" );
    TH1F*                  createHistY       ( const String&          name = "" );
    TH2F*                  createHist2D      ( const String&          name = "" );

    TH2F*                  createVarBinHist2DX ( const int&             nBinX, 
                                                 const double*          pBinsX,
                                                 const String&          name = "" );
    TH2F*                  createVarBinHist2DY ( const int&             nBinY, 
                                                 const double*          pBinsY,
                                                 const String&          name = "" );
    TH2F*                  createVarBinHist2D  ( const int&             nBinX, 
                                                 const double*          pBinsX,
                                                 const int&             nBinY, 
                                                 const double*          pBinsY,
                                                 const String&          name = "" );
};
    

#endif // SH_EFF_MAKER_H
