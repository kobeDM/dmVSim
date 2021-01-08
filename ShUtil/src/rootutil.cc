//////////////////////////////////////////////////////////////////
//
// ROOT Utility
//
// 2016/4/19
// Satoshi Higashino
// satoshi.higashino@cern.ch
//
//////////////////////////////////////////////////////////////////
#include "inc/rootutil.h"

//////////////////////////////////////////////////////////////////
//
// Histogram
//
//////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------
// Set title range
void ShTUtil::SetXTitle( TH1* pHist, const String& title )
{
    if( pHist != nullptr || pHist->GetXaxis( ) != nullptr ) {
        pHist->GetXaxis( )->SetTitle( title.c_str( ) );
    }
}
void ShTUtil::SetYTitle( TH1* pHist, const String& title )
{
    if( pHist != nullptr || pHist->GetYaxis( ) != nullptr ) {
        pHist->GetYaxis( )->SetTitle( title.c_str( ) );
    }
}
void ShTUtil::SetZTitle( TH1* pHist, const String& title )
{
    if( pHist != nullptr || pHist->GetZaxis( ) != nullptr ) {
        pHist->GetZaxis( )->SetTitle( title.c_str( ) );
    }
}

// ---------------------------------------------------------------
// Set range
void ShTUtil::SetXRange( TH1* pHist, const Double_t& begin, const Double_t& end )
{
    if( pHist != nullptr || pHist->GetXaxis( ) != nullptr ) {
        if( begin < end )
            pHist->GetXaxis( )->SetRangeUser( begin, end );
    }
}
void ShTUtil::SetYRange( TH1* pHist, const Double_t& begin, const Double_t& end )
{
    if( pHist != nullptr || pHist->GetYaxis( ) != nullptr ) {
        if( begin < end )
            pHist->GetYaxis( )->SetRangeUser( begin, end );
    }
}
void ShTUtil::SetZRange( TH1* pHist, const Double_t& begin, const Double_t& end )
{
    if( pHist != nullptr || pHist->GetZaxis( ) != nullptr ) {
        if( begin < end )
            pHist->GetZaxis( )->SetRangeUser( begin, end );
    }
}

// ---------------------------------------------------------------
// Set draw style
void ShTUtil::SetStyles( TH1*    pHist,
                         Color_t lineColor,
                         Width_t lineWidth,
                         Color_t markColor,
                         Size_t  markSize, 
                         Style_t markStyle )
{
    if( pHist != nullptr ) {
        pHist->SetLineColor  ( lineColor );
        pHist->SetLineWidth  ( lineWidth );
        pHist->SetMarkerColor( markColor );
        pHist->SetMarkerSize ( markSize  );
        pHist->SetMarkerStyle( markStyle );
    }
}

// ---------------------------------------------------------------
// Set fill options
void ShTUtil::SetFillOption( TH1*    pHist,
                             Color_t fillColor,
                             Style_t fillStyle,
                             Bool_t  editLine,
                             Float_t transparency,
                             Color_t lineColor,
                             Width_t lineWidth,
                             Style_t lineStyle )
{
    if( pHist != nullptr ) {
        pHist->SetFillColorAlpha( fillColor, transparency );
        pHist->SetFillStyle     ( fillStyle );

        Color_t color = fillColor;
        Style_t style = 20;
        Width_t width = 2;
        if( editLine == true ) {
            color = lineColor;
            style = lineStyle;
            width = lineWidth;
        }

        pHist->SetLineColor     ( color );
        pHist->SetLineWidth     ( width );
        pHist->SetLineStyle     ( style );
    }
}

// ---------------------------------------------------------------
// Get histogram from TDirectory
TObject* ShTUtil::ReadObjFromDir( TDirectory* pDir, const String& keyName )
{
    TObject* pObj = nullptr;
    if( pDir == nullptr ) return pObj;

    TKey* pKey = pDir->GetKey( keyName.c_str( ) );
    if( pKey == nullptr ) return pObj;

    pObj = pKey->ReadObj( );
    return pObj;
}

bool ShTUtil::ReadObjListFromDir( TDirectory* pDir, std::list< TObject* >* pList )
{
    if( pDir == nullptr || pList == nullptr ) return false;

    TList* pTList = pDir->GetListOfKeys( );
    if( pTList == nullptr ) return false;
    TIter next( pTList );
    TKey* pKey;
    while( ( pKey = (TKey*)next( ) ) ) {
        if( pKey == nullptr ) continue;
        if( pKey->ReadObj( ) != nullptr )
            pList->push_back( pKey->ReadObj( ) );
    }
    return true;
}

// ---------------------------------------------------------------
// Get histogram from TDirectory
bool ShTUtil::SetPointPoissonError( TGraphAsymmErrors* pGraph,
                                    const int&         index,
                                    const double&      valx,
                                    const double&      valy,
                                    const double&      weight )
{
    if( pGraph == nullptr ) return false;

    double errUp = 0.0, errDown = 0.0;
    double y1 = 0.0, y2 = 0.0, d1 = 0.0, d2 = 0.0;
    if( valy > 0.0 ) {
        y1 = valy + 1.0;
        d1 = 1.0 - 1.0 / ( 9.0 * y1 ) + 1.0 / (3.0 * sqrt( y1 ) );
        errUp = y1 * d1 * d1 * d1 - valy;
        y2 = valy;
        d2 = 1.0 - 1.0 / ( 9.0 * y2 ) - 1.0 / (3.0 * sqrt( y2 ) );
        errDown = valy - y2 * d2 * d2 * d2;

        pGraph->SetPoint( index, valx, valy / weight );
        pGraph->SetPointError( index, 0.0, 0.0, errDown / weight, errUp / weight );
    }
    else {
        pGraph->SetPoint( index, valx, -999.0 );
        pGraph->SetPointError( index, 0.0, 0.0, 0.0, 0.0 );
    }
        
    return true;
}

//////////////////////////////////////////////////////////////////
//
// THStack
//
//////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------
// Set title range
void ShTUtil::SetXTitle( THStack* pStack, const String& title )
{
    if( pStack != nullptr || pStack->GetXaxis( ) != nullptr ) {
        pStack->GetXaxis( )->SetTitle( title.c_str( ) );
    }
}
// ---------------------------------------------------------------
// Set title range
void ShTUtil::SetYTitle( THStack* pStack, const String& title )
{
    if( pStack != nullptr || pStack->GetYaxis( ) != nullptr ) {
        pStack->GetYaxis( )->SetTitle( title.c_str( ) );
    }
}

//////////////////////////////////////////////////////////////////
//
// Function
//
//////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------
// Set draw style
void ShTUtil::SetStyles( TF1*    pFunc,
                         Color_t lineColor,
                         Width_t lineWidth )
{
    if( pFunc != nullptr ) {
        pFunc->SetLineColor  ( lineColor );
        pFunc->SetLineWidth  ( lineWidth );
    }
}

//////////////////////////////////////////////////////////////////
//
// ROOT Visual Objetcs
//
//////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------
// Create text in your plots
//   note: Detail is shown in : https://root.cern.ch/doc/master/classTLatex.html
//
TLatex* ShTUtil::CreateDrawText( const double&  x,
                                 const double&  y,
                                 const String&  text,
                                 const double&  size,
                                 const Color_t& color )
{
    if( text.size( ) <= 0 ) return nullptr;
    TLatex l;
    l.SetNDC( );
    l.SetTextColor( color );
    l.SetTextSize( size );
    return l.DrawLatex( x, y, text.c_str( ) );
}

// ---------------------------------------------------------------
// Create legend in your plots
//   note: Detail is shown in : https://root.cern.ch/doc/master/classTLegend.html
//
TLegend* ShTUtil::CreateLegend( const double& x1,
                                const double& y1,
                                const double& x2,
                                const double& y2 )
{
    TLegend* pLegend = new TLegend( x1, y1, x2, y2 );
    pLegend->SetFillStyle( 0 );
    pLegend->SetBorderSize( 0 );
    pLegend->SetTextFont( 42 );
    return pLegend;
}

