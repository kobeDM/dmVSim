#include "inc/effMaker.h"

// ---------------------------------------------------------------
// EffMake (1D)
// ---------------------------------------------------------------
EffMaker::EffMaker( const String& name, const int& nBin, const double& rangeMin, const double& rangeMax )
    : m_name      ( name     )
    , m_nBin      ( nBin     )
    , m_rangeMin  ( rangeMin )
    , m_rangeMax  ( rangeMax )
{
    m_denArr.resize( nBin + 1 );
    m_numArr.resize( nBin + 1 );
}

EffMaker::~EffMaker( )
{
}

std::vector< double >*
EffMaker::getDenominator( )
{
    return &m_denArr;
}

void
EffMaker::getDenominator( std::vector< double >* pDenArr )
{
    pDenArr = &m_denArr;
    return;
}

std::vector< double >*
EffMaker::getNumerator( )
{
    return &m_numArr;
}

void
EffMaker::getNumerator( std::vector< double >* pNumArr )
{
    pNumArr = &m_numArr;
    return;
}

bool
EffMaker::setDenominator( const int&    bin,
                          const double& val )
{
    if( m_denArr.size( ) >= bin ) {
        ShUtil::Cerr( "EffMaker::setDenominator( ): Invalid access to the denominator." );
        return false;
    }

    m_denArr.at( bin ) = val;
    return true;
}

bool
EffMaker::setNumerator( const int&    bin,
                        const double& val )
{
    if( m_numArr.size( ) >= bin ) {
        ShUtil::Cerr( "EffMaker::setNumerator( ): Invalid access to the numerator." );
        return false;
    }

    m_numArr.at( bin ) = val;
    return true;
}

bool
EffMaker::countDenEvt( const double& val )
{
    if( findBin( val ) < 0 || findBin( val ) >= m_denArr.size( ) ) return false;
    m_denArr.at( findBin( val ) ) += 1.0;
    return true;
}
    
bool
EffMaker::countNumEvt( const double& val )
{
    if( findBin( val ) < 0 || findBin( val ) >= m_numArr.size( ) ) return false;
    m_numArr.at( findBin( val ) ) += 1.0;
    return true;
}

int
EffMaker::findBin( const double& val )
{
    if( val < m_rangeMin || val > m_rangeMax ) {
        ShUtil::Cerr( "EffMaker::findBin( ): value out of range" );
        return -1;
    }
    
    double tmpA   = val        - m_rangeMin;
    double tmpB   = m_rangeMax - m_rangeMin;
    double retVal = static_cast< double >( m_nBin ) * tmpA / tmpB;
    
    return static_cast< int >( retVal );
}

double
EffMaker::findValue( const int& bin )
{
    // Note: the "bin" should be count from 0.
    if( bin < 0 || bin > m_nBin ) {
        ShUtil::Cerr( "EffMaker::findValue( ): value out of range" );
        return -100.0;
    }

    double tmpC   = static_cast< double >( bin ) / static_cast< double >( m_nBin );
    double tmpB   = m_rangeMax - m_rangeMin;
    double retVal = tmpC * tmpB + m_rangeMin;

    return retVal;
}

double
EffMaker::getEff( const int& bin )
{
    if( bin < 0 || bin > m_nBin ) {
        ShUtil::Cerr( "EffMaker::getEff( ): value out of range" );
        return -100.0;
    }

    if( m_denArr.at( bin ) < 1.0 || m_numArr.at( bin ) < 1.0 ) return 0.0;
    return m_numArr.at( bin ) / m_denArr.at( bin );
}

double
EffMaker::getEff( const double& val )
{
    int bin = findBin( val );
    return getEff( bin );
}

double
EffMaker::getEffErr( const int& bin )
{
    if( bin < 0 || bin > m_nBin ) {
        ShUtil::Cerr( "EffMaker::getEffErr( ): value out of range" );
        return -100.0;
    }

    if( m_denArr.at( bin ) < 1.0 ) return 0.0;
    return ShPUtil::GetBinomialError( getEff( bin ), static_cast< int >( m_denArr.at( bin ) ) );
}

double
EffMaker::getEffErr( const double& val )
{
    int bin = findBin( val );
    return getEffErr( bin );
}

TH1F*
EffMaker::createHist( const String& name )
{
    String histName = ( name == "" ) ? m_name : name;
    TH1F* pHist = new TH1F( histName.c_str( ), histName.c_str( ), m_nBin, m_rangeMin, m_rangeMax );
    for( int bin = 0; bin < m_nBin; ++bin ) {
        double eff = getEff( bin );
        double effErr = getEffErr( bin );
        double val = findValue( bin ) + 0.00001;
        
        pHist->Fill( val, eff );
        pHist->SetBinError( bin + 1, effErr );
    }

    return pHist;
}

TH1F*
EffMaker::createVarBinHist( const int& nBin, const double* pBins, const String& name )
{
    if( pBins == nullptr ) return nullptr;

    String histName = ( name == "" ) ? m_name : name;

    // create new numerator & denominator counter
    std::vector< double > newDenArr;
    std::vector< double > newNumArr;
    newDenArr.resize( nBin );
    newNumArr.resize( nBin );
    TH1F* pHist = new TH1F( histName.c_str( ), histName.c_str( ), nBin, pBins );
    for( int primBin = 0; primBin < m_nBin; ++primBin ) {
        double val = findValue( primBin ) + 0.00001;
        int newBin = -1;
        for( int bin = 0; bin < nBin; ++bin )
            if( val > pBins[ bin ] && val < pBins[ bin + 1 ] ) newBin = bin;

        if( newBin < 0 || newBin > nBin ) continue;
        newDenArr.at( newBin ) += m_denArr.at( primBin );
        newNumArr.at( newBin ) += m_numArr.at( primBin );
    }

    for( int newBin = 0; newBin < nBin; ++newBin ) {
        double eff = 0.0;
        double effErr = 0.0001;
        double val = pBins[ newBin ] + 0.0001;
        if( newDenArr.at( newBin ) > 0.1 && newNumArr.at( newBin ) > 0.1 )
            eff = newNumArr.at( newBin ) / newDenArr.at( newBin );
        if( newDenArr.at( newBin ) > 0.1 && newNumArr.at( newBin ) > 0.1 )
            effErr = ShPUtil::GetBinomialError( eff, static_cast< int >( newDenArr.at( newBin ) ) );
        
        pHist->Fill( val, eff );
        pHist->SetBinError( newBin + 1, effErr );
    }

    return pHist;
}


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
EffMaker2D::EffMaker2D( const String& name,
                        const int&    nBinX, const double& rangeMinX, const double& rangeMaxX,
                        const int&    nBinY, const double& rangeMinY, const double& rangeMaxY )
    : m_name      ( name      )
    , m_nBinX     ( nBinX     )
    , m_rangeMinX ( rangeMinX )
    , m_rangeMaxX ( rangeMaxX )
    , m_nBinY     ( nBinY     )
    , m_rangeMinY ( rangeMinY )
    , m_rangeMaxY ( rangeMaxY )
{
    m_denArr.reserve( nBinX );
    m_numArr.reserve( nBinX );
    for( int i = 0; i < nBinX; ++i ) {
        std::vector< double > tmpVector;
        tmpVector.resize( nBinY );
        m_denArr.push_back( tmpVector );
        m_numArr.push_back( tmpVector );
    }
}

EffMaker2D::~EffMaker2D( )
{
}

std::vector< std::vector< double > >* 
EffMaker2D::getDenominator2D( )
{
    return &m_denArr;
}

void
EffMaker2D::getDenominator2D( std::vector< std::vector< double > >* pDenArr )
{
    pDenArr = &m_denArr;
    return;
}

std::vector< double >*
EffMaker2D::getDenominatorX( const double&          val )
{
    // not yet impremented
    return nullptr;
}

void
EffMaker2D::getDenominatorX( std::vector< double >* pDenArr,
                             const double&          val)
{
    // not yet impremented
    return;
}

std::vector< double >*
EffMaker2D::getDenominatorY( const double&          val )
{
    // not yet impremented
    return nullptr;
}

void
EffMaker2D::getDenominatorY( std::vector< double >* pDenArr,
                             const double&          val )
{
    // not yet impremented
    return;
}

std::vector< std::vector< double > >*
EffMaker2D::getNumerator2D( )
{
    // not yet impremented
    return nullptr;
}

void
EffMaker2D::getNumerator2D( std::vector< std::vector< double > >* pDenArr )
{
    // not yet impremented
    return;
}

std::vector< double >*
EffMaker2D::getNumeratorX( const double&          val )
{
    // not yet impremented
    return nullptr;
}

void
EffMaker2D::getNumeratorX( std::vector< double >* pNumArr,
                           const double&          val)
{
    // not yet impremented
    return;
}

std::vector< double >*
EffMaker2D::getNumeratorY( const double&          val )
{
    // not yet impremented
    return nullptr;
}

void
EffMaker2D::getNumeratorY( std::vector< double >* pNumArr,
                           const double&          val)
{
    // not yet impremented
    return;
}

bool
EffMaker2D::setDenominator( const int&             binX,
                            const int&             binY,
                            const double&          val )
{
    if( m_denArr.size( ) >= binX ) {
        ShUtil::Cerr( "EffMaker::setDenominator( ): Invalid access to the denominator X." );
        return false;
        if( m_denArr.at( binX ).size( ) >= binY ) {
            ShUtil::Cerr( "EffMaker::setDenominator( ): Invalid access to the denominator Y." );
            return false;
        }
    }

    m_denArr.at( binX ).at( binY ) = val;
    return true;
}

bool
EffMaker2D::setNumerator( const int&             binX,
                          const int&             binY,
                          const double&          val )
{
    if( m_numArr.size( ) >= binX ) {
        ShUtil::Cerr( "EffMaker::setNumerator( ): Invalid access to the numerator X." );
        return false;
        if( m_numArr.at( binX ).size( ) >= binY ) {
            ShUtil::Cerr( "EffMaker::setNumerator( ): Invalid access to the numerator Y." );
            return false;
        }
    }

    m_numArr.at( binX ).at( binY ) = val;
    return true;
}

bool
EffMaker2D::countDenEvt( const double&          valX,
                         const double&          valY )
{
    if( findBinX( valX ) < 0                 || findBinY( valY ) < 0                 ||
        findBinX( valX ) >= m_denArr.size( ) || findBinY( valY ) >= m_denArr.size( ) ) return false;
    m_denArr.at( findBinX( valX ) ).at( findBinY( valY ) ) += 1.0;
    return true;
}



bool
EffMaker2D::countNumEvt( const double&          valX,
                         const double&          valY )
{
    if( findBinX( valX ) < 0                 || findBinY( valY ) < 0                 ||
        findBinX( valX ) >= m_numArr.size( ) || findBinY( valY ) >= m_numArr.size( ) ) return false;
    m_numArr.at( findBinX( valX ) ).at( findBinY( valY ) ) += 1.0;
    return true;
}

int
EffMaker2D::findBinX( const double& val )
{
    if( val < m_rangeMinX || val > m_rangeMaxX ) {
        ShUtil::Cerr( "EffMaker::findBinX( ): value out of range" );
        return -1;
    }
    double tmpA   = val         - m_rangeMinX;
    double tmpB   = m_rangeMaxX - m_rangeMinX;
    double retVal = static_cast< double >( m_nBinX ) * tmpA / tmpB;
    
    return static_cast< int >( retVal );
}

int
EffMaker2D::findBinY( const double& val )
{
    if( val < m_rangeMinY || val > m_rangeMaxY ) {
        ShUtil::Cerr( "EffMaker::findBinY( ): value out of range" );
        return -1;
    }
    double tmpA   = val         - m_rangeMinY;
    double tmpB   = m_rangeMaxY - m_rangeMinY;
    double retVal = static_cast< double >( m_nBinY ) * tmpA / tmpB;
    
    return static_cast< int >( retVal );
}

double
EffMaker2D::findValueX( const int&             bin )
{
    // Note: the "bin" should be count from 0.
    if( bin < 0 || bin > m_nBinX ) {
        ShUtil::Cerr( "EffMaker::findValueX( ): value out of range" );
        return -100.0;
    }

    double tmpC   = static_cast< double >( bin ) / static_cast< double >( m_nBinX );
    double tmpB   = m_rangeMaxX - m_rangeMinX;
    double retVal = tmpC * tmpB + m_rangeMinX;

    return retVal;
}

double
EffMaker2D::findValueY( const int&             bin )
{
    // Note: the "bin" should be count from 0.
    if( bin < 0 || bin > m_nBinY ) {
        ShUtil::Cerr( "EffMaker::findValueY( ): value out of range" );
        return -100.0;
    }

    double tmpC   = static_cast< double >( bin ) / static_cast< double >( m_nBinY );
    double tmpB   = m_rangeMaxY - m_rangeMinY;
    double retVal = tmpC * tmpB + m_rangeMinY;

    return retVal;
}

double
EffMaker2D::getEff( const int&             binX,
                    const int&             binY )
{
    if( binX < 0 || binX > m_nBinX ||
        binX < 0 || binY > m_nBinY ) {
        ShUtil::Cerr( "EffMaker::getEff( ): value out of range" );
        return -100.0;
    }

    if( m_denArr.at( binX ).at( binY ) < 1.0 || m_numArr.at( binX ).at( binY ) < 1.0 ) return 0.0;
    return m_numArr.at( binX ).at( binY ) / m_denArr.at( binX ).at( binY );
}

double
EffMaker2D::getEff( const double&          valX,
                    const double&          valY )
{
    return getEff( findBinX( valX ), findBinY( valY ) );
}


double
EffMaker2D::getEffErr( const int&             binX,
                       const int&             binY )
{
    if( binX < 0 || binX > m_nBinX ||
        binY < 0 || binY > m_nBinY ) {
        ShUtil::Cerr( "EffMaker::getEffErr( ): value out of range" );
        return -100.0;
    }

    if( m_denArr.at( binX ).at( binY ) < 1.0 ) return 0.0;
    return ShPUtil::GetBinomialError( getEff( binX, binY ), static_cast< int >( m_denArr.at( binX ).at( binY ) ) );

}

double
EffMaker2D::getEffErr( const double&          valX,
                       const double&          valY )
{
    return getEffErr( findBinX( valX ), findBinY( valY ) );
}

TH1F*
EffMaker2D::createHistX( const String&          name )
{
    String histName = ( name == "" ) ? m_name : name;
    TH1F* pHist = new TH1F( histName.c_str( ), histName.c_str( ), m_nBinX, m_rangeMinX, m_rangeMaxX );
    for( int binX = 0; binX < m_nBinX; ++binX ) {
        double den = 0.0;
        double num = 0.0;
        for( int binY = 0; binY < m_nBinY; ++binY ) {
            den += m_denArr.at( binX ).at( binY );
            num += m_numArr.at( binX ).at( binY );
        }
        double eff = num / den;
        double effErr = ShPUtil::GetBinomialError( eff, static_cast< int >( den ) );
        double val = findValueX( binX ) + 0.00001;
        
        pHist->Fill( val, eff );
        pHist->SetBinError( binX + 1, effErr );
    }

    return pHist;    
}

TH1F*
EffMaker2D::createHistY( const String&          name )
{
    String histName = ( name == "" ) ? m_name : name;
    TH1F* pHist = new TH1F( histName.c_str( ), histName.c_str( ), m_nBinY, m_rangeMinY, m_rangeMaxY );
    for( int binY = 0; binY < m_nBinY; ++binY ) {
        double den = 0.0;
        double num = 0.0;
        for( int binX = 0; binX < m_nBinX; ++binX ) {
            den += m_denArr.at( binX ).at( binY );
            num += m_numArr.at( binX ).at( binY );
        }
        double eff = num / den;
        double effErr = ShPUtil::GetBinomialError( eff, static_cast< int >( den ) );
        double val = findValueY( binY ) + 0.00001;
        
        pHist->Fill( val, eff );
        pHist->SetBinError( binY + 1, effErr );
    }

    return pHist;

}

TH2F*
EffMaker2D::createHist2D( const String&          name )
{
    String histName = ( name == "" ) ? m_name : name;
    TH2F* pHist = new TH2F( histName.c_str( ), histName.c_str( ),
                            m_nBinX, m_rangeMinX, m_rangeMaxX,
                            m_nBinY, m_rangeMinY, m_rangeMaxY );
    for( int binX = 0; binX < m_nBinX; ++binX ) {
        for( int binY = 0; binY < m_nBinY; ++binY ) {
            double eff = getEff( binX, binY );
            double valX = findValueX( binX ) + 0.00001;
            double valY = findValueY( binY ) + 0.00001;
            pHist->Fill( valX, valY, eff );
        }
    }

    return pHist;
}

TH2F*
EffMaker2D::createVarBinHist2DX( const int&             nBinX, 
                                 const double*          pBinsX,
                                 const String&          name )
{
    if( pBinsX == nullptr ) return nullptr;
    String histName = ( name == "" ) ? m_name : name;
    int nBinY = m_nBinY;

    // create new numerator & denominator counter
    std::vector< std::vector< double > > newDenArr;
    std::vector< std::vector< double > > newNumArr;
    newDenArr.reserve( nBinX );
    newNumArr.reserve( nBinX );
    for( int i = 0; i < nBinX; ++i ) {
        std::vector< double > tmpVector;
        tmpVector.resize( nBinY );
        newDenArr.push_back( tmpVector );
        newNumArr.push_back( tmpVector );
    }

    TH2F* pHist = new TH2F( histName.c_str( ), histName.c_str( ), nBinX, pBinsX, nBinY, m_rangeMinY, m_rangeMaxY );
    for( int primBinX = 0; primBinX < m_nBinX; ++primBinX ) {
        for( int primBinY = 0; primBinY < m_nBinY; ++primBinY ) {
            double valX = findValueX( primBinX ) + 0.00001;
            double valY = findValueY( primBinY ) + 0.00001;
            int newBinX = -1, newBinY = primBinY;
            for( int binX = 0; binX < nBinX; ++binX )
                if( valX > pBinsX[ binX ] && valX < pBinsX[ binX + 1 ] ) newBinX = binX;
            
            if( newBinX < 0 || newBinY < 0 || newBinX > nBinX || newBinY > nBinY ) continue;
            newDenArr.at( newBinX ).at( newBinY ) += m_denArr.at( primBinX ).at( primBinY );
            newNumArr.at( newBinX ).at( newBinY ) += m_numArr.at( primBinX ).at( primBinY );
        }
    }

    for( int newBinX = 0; newBinX < nBinX; ++newBinX ) {
        for( int newBinY = 0; newBinY < nBinY; ++newBinY ) {
            double eff = -1.0;
            double valX = pBinsX[ newBinX ] + 0.00001;
            double valY = findValueY( newBinY ) + 0.00001;
            if( newDenArr.at( newBinX ).at( newBinY ) > 0.1 && newNumArr.at( newBinX ).at( newBinY ) > 0.1 )
                eff = newNumArr.at( newBinX ).at( newBinY ) / newDenArr.at( newBinX ).at( newBinY );
            
            pHist->Fill( valX, valY, eff );
        }
    }

    pHist->GetZaxis( )->SetRangeUser( 0.0, 1.0 );
    return pHist;
}


TH2F*
EffMaker2D::createVarBinHist2DY( const int&             nBinY, 
                                 const double*          pBinsY,
                                 const String&          name )
{
    if( pBinsY == nullptr ) return nullptr;
    String histName = ( name == "" ) ? m_name : name;
    int nBinX = m_nBinX;

    // create new numerator & denominator counter
    std::vector< std::vector< double > > newDenArr;
    std::vector< std::vector< double > > newNumArr;
    newDenArr.reserve( nBinX );
    newNumArr.reserve( nBinX );
    for( int i = 0; i < nBinX; ++i ) {
        std::vector< double > tmpVector;
        tmpVector.resize( nBinY );
        newDenArr.push_back( tmpVector );
        newNumArr.push_back( tmpVector );
    }

    TH2F* pHist = new TH2F( histName.c_str( ), histName.c_str( ), nBinX, m_rangeMinX, m_rangeMaxX, nBinY,pBinsY );
    for( int primBinX = 0; primBinX < m_nBinX; ++primBinX ) {
        for( int primBinY = 0; primBinY < m_nBinY; ++primBinY ) {
            double valX = findValueX( primBinX ) + 0.00001;
            double valY = findValueY( primBinY ) + 0.00001;
            int newBinX = primBinX, newBinY = -1;
            for( int binY = 0; binY < nBinY; ++binY )
                if( valY > pBinsY[ binY ] && valY < pBinsY[ binY + 1 ] ) newBinY = binY;
            
            if( newBinX < 0 || newBinY < 0 || newBinX > nBinX || newBinY > nBinY ) continue;
            newDenArr.at( newBinX ).at( newBinY ) += m_denArr.at( primBinX ).at( primBinY );
            newNumArr.at( newBinX ).at( newBinY ) += m_numArr.at( primBinX ).at( primBinY );
        }
    }

    for( int newBinX = 0; newBinX < nBinX; ++newBinX ) {
        for( int newBinY = 0; newBinY < nBinY; ++newBinY ) {
            double eff = -1.0;
            double valX = findValueX( newBinX ) + 0.00001;
            double valY = pBinsY[ newBinY ] + 0.00001;
            if( newDenArr.at( newBinX ).at( newBinY ) > 0.1 && newNumArr.at( newBinX ).at( newBinY ) > 0.1 )
                eff = newNumArr.at( newBinX ).at( newBinY ) / newDenArr.at( newBinX ).at( newBinY );
            
            pHist->Fill( valX, valY, eff );
        }
    }

    pHist->GetZaxis( )->SetRangeUser( 0.0, 1.0 );
    return pHist;
}

TH2F*
EffMaker2D::createVarBinHist2D( const int&             nBinX,
                                const double*          pBinsX,
                                const int&             nBinY,
                                const double*          pBinsY,
                                const String&          name )
{
    if( pBinsX == nullptr || pBinsY == nullptr ) return nullptr;
    String histName = ( name == "" ) ? m_name : name;

    // create new numerator & denominator counter
    std::vector< std::vector< double > > newDenArr;
    std::vector< std::vector< double > > newNumArr;
    newDenArr.reserve( nBinX );
    newNumArr.reserve( nBinX );
    for( int i = 0; i < nBinX; ++i ) {
        std::vector< double > tmpVector;
        tmpVector.resize( nBinY );
        newDenArr.push_back( tmpVector );
        newNumArr.push_back( tmpVector );
    }

    TH2F* pHist = new TH2F( histName.c_str( ), histName.c_str( ), nBinX, pBinsX, nBinY, pBinsY );
    for( int primBinX = 0; primBinX < m_nBinX; ++primBinX ) {
        for( int primBinY = 0; primBinY < m_nBinY; ++primBinY ) {
            double valX = findValueX( primBinX ) + 0.00001;
            double valY = findValueY( primBinY ) + 0.00001;
            int newBinX = -1, newBinY = -1;
            for( int binX = 0; binX < nBinX; ++binX )
                if( valX > pBinsX[ binX ] && valX < pBinsX[ binX + 1 ] ) newBinX = binX;
            for( int binY = 0; binY < nBinY; ++binY )
                if( valY > pBinsY[ binY ] && valY < pBinsY[ binY + 1 ] ) newBinY = binY;
            
            if( newBinX < 0 || newBinY < 0 || newBinX > nBinX || newBinY > nBinY ) continue;
            newDenArr.at( newBinX ).at( newBinY ) += m_denArr.at( primBinX ).at( primBinY );
            newNumArr.at( newBinX ).at( newBinY ) += m_numArr.at( primBinX ).at( primBinY );
        }
    }

    for( int newBinX = 0; newBinX < nBinX; ++newBinX ) {
        for( int newBinY = 0; newBinY < nBinY; ++newBinY ) {
            double eff = -1.0;
            double valX = pBinsX[ newBinX ] + 0.0001;
            double valY = pBinsY[ newBinY ] + 0.0001;
            if( newDenArr.at( newBinX ).at( newBinY ) > 0.1 && newNumArr.at( newBinX ).at( newBinY ) > 0.1 )
                eff = newNumArr.at( newBinX ).at( newBinY ) / newDenArr.at( newBinX ).at( newBinY );
            
            pHist->Fill( valX, valY, eff );
        }
    }

    pHist->GetZaxis( )->SetRangeUser( 0.0, 1.0 );
    return pHist;    
}
