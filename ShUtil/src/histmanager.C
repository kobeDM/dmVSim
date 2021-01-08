#include "inc/histmanager.h"

// Create and store TH1F histogram
bool
ShTHistManager::createTH1F( const String& name,
                            const int&    Nbins,
                            const double& xmin,
                            const double& xmax,
                            const String& title )
{
    if( name.length( ) <= 0 ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH1F: invalid histogram name : " + name );
        return false;
    }
    if( hasTH1F( name ) == true ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH1F: Attempt to create second histogram named " + name );
        return false;
    }
    
    m_histTH1F[ name ] = new TH1F( name.c_str( ), title.c_str( ), Nbins, xmin, xmax );
    m_histTH1F[ name ]->Sumw2( );
    return true;
}

// Create and Store TH1F histogram
bool
ShTHistManager::createTH1F( const String&                name,
                            const std::vector< double >& bins,
                            const String&                title )
{
    if( name.length( ) <= 0 ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH1F: invalid histogram name : " + name );
        return false;
    }
    if( hasTH1F( name ) == true ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH1F: Attempt to create second histogram named " + name );
        return false;
    }
    m_histTH1F[ name ] = new TH1F( name.c_str( ), title.c_str( ), -1 + bins.size( ), &bins[ 0 ] );
    m_histTH1F[ name ]->Sumw2( );
    return true;
}

// Create and Store TH1F histogram
bool
ShTHistManager::createTH1F( const String&                name,
                            const int&                   nBin,
                            const double*                pBins,
                            const String&                title )
{
    if( pBins == nullptr ) return false;
    if( name.length( ) <= 0 ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH1F: invalid histogram name : " + name );
        return false;
    }
    if( hasTH1F( name ) == true ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH1F: Attempt to create second histogram named " + name );
        return false;
    }

    m_histTH1F[ name ] = new TH1F( name.c_str( ), title.c_str( ), nBin, pBins );
    m_histTH1F[ name ]->Sumw2( );

    return true;
}

// Create and store TH2F histogram
bool
ShTHistManager::createTH2F( const String& name,
                            const int&    NbinsX,
                            const double& xmin,
                            const double& xmax,
                            const int&    NBinsY,
                            const double& ymin,
                            const double& ymax,
                            const String& title)
{
    if( name.length( ) <= 0 ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH2F: invalid histogram name : " + name );
        return false;
    }
    if( hasTH2F( name ) == true ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH2F: Attempt to create second histogram named " + name );
        return false;
    }
    m_histTH2F[ name ] = new TH2F( name.c_str( ), title.c_str( ), NbinsX, xmin, xmax, NBinsY, ymin, ymax );
    m_histTH2F[ name ]->Sumw2( );
    return true;
}

// Create and store TH2F histogram
bool
ShTHistManager::createTH2F( const String&                name,
                            const std::vector< double >& xbins,
                            const std::vector< double >& ybins,
                            const String&                title)
{
    if( name.length( ) <= 0 ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH2F: invalid histogram name : " + name );
        return false;
    }
    if( hasTH2F( name ) == true ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH2F: Attempt to create second histogram named " + name );
        return false;
    }
    m_histTH2F[ name ] = new TH2F( name.c_str( ), title.c_str( ), xbins.size( ) - 1, &xbins[ 0 ],  ybins.size( ) - 1, &ybins[ 0 ] );
    m_histTH2F[ name ]->Sumw2( );
    return true;
}

// Create and Store TH2F histogram
bool
ShTHistManager::createTH2F( const String&                name,
                            const int&                   nBinX,
                            const double*                pBinsX,
                            const int&                   nBinY,
                            const double*                pBinsY,
                            const String&                title )
{
    if( pBinsX == nullptr || pBinsY == nullptr ) return false;
    if( name.length( ) <= 0 ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH2F: invalid histogram name : " + name );
        return false;
    }
    if( hasTH2F( name ) == true ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH2F: Attempt to create second histogram named " + name );
        return false;
    }

    m_histTH2F[ name ] = new TH2F( name.c_str( ), title.c_str( ), nBinX, pBinsX, nBinY, pBinsY );
    m_histTH2F[ name ]->Sumw2( );

    return true;
}


bool
ShTHistManager::createTH2F( const String&              name,
                            const int&                 nBinX,
                            const double*              pBinsX,
                            const int&                 NBinsY,
                            const double&              ymin,
                            const double&              ymax,
                            const String&              title )
{
    if( pBinsX == nullptr ) return false;
    if( name.length( ) <= 0 ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH2F: invalid histogram name : " + name );
        return false;
    }
    if( hasTH2F( name ) == true ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH2F: Attempt to create second histogram named " + name );
        return false;
    }

    m_histTH2F[ name ] = new TH2F( name.c_str( ), title.c_str( ), nBinX, pBinsX, NBinsY, ymin, ymax );
    m_histTH2F[ name ]->Sumw2( );

    return true;
}


bool
ShTHistManager::createTH2F( const String&              name,
                            const int&                 NBinsX,
                            const double&              xmin,
                            const double&              xmax,
                            const int&                 nBinY,
                            const double*              pBinsY,
                            const String&              title )
{
    if( pBinsY == nullptr ) return false;
    if( name.length( ) <= 0 ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH2F: invalid histogram name : " + name );
        return false;
    }
    if( hasTH2F( name ) == true ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH2F: Attempt to create second histogram named " + name );
        return false;
    }

    m_histTH2F[ name ] = new TH2F( name.c_str( ), title.c_str( ), NBinsX, xmin, xmax, nBinY, pBinsY );
    m_histTH2F[ name ]->Sumw2( );

    return true;
}


// Create and store TH3F histogram
bool
ShTHistManager::createTH3F( const String& name,
                            const int&    NbinsX,
                            const double& xmin,
                            const double& xmax,
                            const int&    NBinsY,
                            const double& ymin,
                            const double& ymax,
                            const int&    NBinsZ,
                            const double& zmin,
                            const double& zmax,
                            const String& title )
{
    if( name.length( ) <= 0 ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH3F: invalid histogram name : " + name );
        return false;
    }
    if( hasTH3F( name ) == true ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH3F: Attempt to create second histogram named " + name );
        return false;
    }
    m_histTH3F[ name ] = new TH3F( name.c_str( ), title.c_str( ), NbinsX, xmin, xmax, NBinsY, ymin, ymax, NBinsZ, zmin, zmax );
    m_histTH3F[ name ]->Sumw2( );
    return true;
}

// Create and store TH3F histogram
bool
ShTHistManager::createTH3F( const String&                name,
                            const std::vector< double >& xbins,
                            const std::vector< double >& ybins,
                            const std::vector< double >& zbins,
                            const String&                title)
{
    if( name.length( ) <= 0 ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH3F: invalid histogram name : " + name );
        return false;
    }
    if( hasTH3F( name ) == true ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH3F: Attempt to create second histogram named " + name );
        return false;
    }
    m_histTH3F[ name ] = new TH3F( name.c_str( ), title.c_str( ),
                                   xbins.size( ) - 1, &xbins[ 0 ],
                                   ybins.size( ) - 1, &ybins[ 0 ],
                                   zbins.size( ) - 1, &zbins[ 0 ] );
    m_histTH3F[ name ]->Sumw2( );
    return true;
}

// Create and Store TH3F histogram
bool
ShTHistManager::createTH3F( const String&                name,
                            const int&                   nBinX,
                            const double*                pBinsX,
                            const int&                   nBinY,
                            const double*                pBinsY,
                            const int&                   nBinZ,
                            const double*                pBinsZ,
                            const String&                title )
{
    if( pBinsX == nullptr || pBinsY == nullptr || pBinsZ == nullptr ) return false;
    if( name.length( ) <= 0 ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH1F: invalid histogram name : " + name );
        return false;
    }
    if( hasTH1F( name ) == true ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTH1F: Attempt to create second histogram named " + name );
        return false;
    }

    m_histTH3F[ name ] = new TH3F( name.c_str( ), title.c_str( ), nBinX, pBinsX, nBinY, pBinsY, nBinZ, pBinsZ );
    m_histTH3F[ name ]->Sumw2( );

    return true;
}

// Create and store TProfile histogram
bool
ShTHistManager::createTProfile( const String& name,
                                const int&    NbinsX,
                                const double& xmin,
                                const double& xmax,
                                const String& title )
{
    if( name.length( ) <= 0 ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTProfile: invalid histogram name : " + name );
        return false;
    }
    if( hasTProfile( name ) == true ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTProfile: Attempt to create second histogram named " + name );
        return false;
    }
    m_histTProfile[ name ] = new TProfile( name.c_str( ), title.c_str( ), NbinsX, xmin, xmax );
    m_histTProfile[ name ]->Sumw2( );
    return true;
}

// Create and store TProfile histogram
bool
ShTHistManager::createTProfile( const String&                name,
                                const std::vector< double >& xbins,
                                const String&                title)
{
    if( name.length( ) <= 0 ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTProfile: invalid histogram name : " + name );
        return false;
    }
    if( hasTProfile( name ) == true ) {
        ShUtil::Cerr( "ShTHistManager::createHistoTProfile: Attempt to create second histogram named " + name );
        return false;
    }
    m_histTProfile[ name ] = new TProfile( name.c_str( ), title.c_str( ), xbins.size( ) - 1, &xbins[ 0 ] );
    m_histTProfile[ name ]->Sumw2( );
    return true;
}

// Retrieve a TH1F histogram from the internal store
TH1F*
ShTHistManager::getTH1F( const String& name )
{
    if( hasTH1F( name ) == false ) {
        ShUtil::Cerr( "HistoStore::getTH1F requested histogram " + name + " cannot be accessed (did you forget to declare it?)" );
        return nullptr;
    }

    return m_histTH1F[ name ];
}

// Retrieve a TH2F histogram from the internal store
TH2F*
ShTHistManager::getTH2F( const String& name )
{
    if( hasTH2F( name ) == false ) {
        ShUtil::Cerr( "HistoStore::getTH2F requested histogram " + name + " cannot be accessed (did you forget to declare it?)" );
        return nullptr;
    }
    return m_histTH2F[ name ];
}

// Retrieve a TH3F histogram from the internal store
TH3F*
ShTHistManager::getTH3F( const String& name )
{
    if( hasTH3F( name ) == false ) {
        ShUtil::Cerr( "HistoStore::getTH3F requested histogram " + name + " cannot be accessed (did you forget to declare it?)" );
        return nullptr;
    }
    return m_histTH3F[ name ];
}

// Retrieve a TH2F histogram from the internal store
TProfile*
ShTHistManager::getTProfile( const String& name )
{
    if( hasTProfile( name ) == false ) {
        ShUtil::Cerr( "HistoStore::getTProfile " + name + " cannot be accessed (did you forget to declare it?)" );
        return nullptr;
    }
    return m_histTProfile[ name ];
}

// Retrieve a TF1 function from the internal store
TF1*
ShTHistManager::getTF1( const String& name )
{
    if( hasTF1( name ) == false ) {
        ShUtil::Cerr( "HistoStore::getTF1 requested function " + name + " cannot be accessed (did you forget to declare it?)" );
        return nullptr;
    }

    return m_funcTF1[ name ];
}

// Return vector of all histograms in internal store
bool
ShTHistManager::getListOfHistograms( std::vector< TH1* >* pArray )
{
    if( pArray == nullptr ) return false;

    pArray->reserve( m_histTH1F.size( ) + m_histTH2F.size( ) + m_histTH3F.size( ) + m_histTProfile.size( ) );
    // an iterator of a map is pair<keyType,valueType>
    for( auto h : m_histTH1F )     pArray->push_back( h.second );
    for( auto h : m_histTH2F )     pArray->push_back( h.second );
    for( auto h : m_histTH3F )     pArray->push_back( h.second );
    for( auto h : m_histTProfile ) pArray->push_back( h.second );
    return true;
}

// Save histograms
// String
// _getSaveOptionStr( const ShTSaveOption& saveAs )
// {
//     String retVal = "pdf";
//     switch( saveAs ) {
//     case pdf:
//         retVal = "pdf";
//         break;
//     case png:
//         retVal = "png";
//         break;
//     case ps:
//         retVal = "ps";
//         break;
//     case eps:
//         retVal = "eps";
//         break;
//     case jpg:
//         retVal = "jpg";
//         break;
//     default:
//         ShUtil::Cwarn( "Invalid save option. Save as default format (pdf)" );
//         break;
//     }

//     return retVal;
// }

bool
ShTHistManager::saveTH1F( const String& outputDir, const String& option )
{
    ShUtil::ExistCreateDir( outputDir );
    TCanvas cvs( "cvs", "cvs", 800, 600 );
    for( auto histPair : m_histTH1F ) {
        if( histPair.second != nullptr )
            histPair.second->Draw( option.c_str( ) );
        double max = histPair.second->GetMaximum( ) * 2.0;
        double min = 0.0;
        histPair.second->GetYaxis()->SetRangeUser( min, max );
        // ATLASLabel( 0.18, 0.88, "Work in Progress", 1, 0.05 );
        // ShTUtil::CreateDrawText( 0.18, 0.8, "#sqrt{s} = 13 TeV, #intLdt = 79.8 fb^{-1}", 0.05 );
        ShTUtil::CreateDrawText( 0.18, 0.88, "#sqrt{s} = 13 TeV", 0.05 );
        cvs.SaveAs( Form( "%s/%s.pdf", outputDir.c_str( ), histPair.first.c_str( ) ) );
        cvs.SaveAs( Form( "%s/%s.eps", outputDir.c_str( ), histPair.first.c_str( ) ) );
        cvs.SaveAs( Form( "%s/%s.png", outputDir.c_str( ), histPair.first.c_str( ) ) );
    }

    return true;
}

bool
ShTHistManager::saveTH2F( const String& outputDir, const String& option )
{
    ShUtil::ExistCreateDir( outputDir );
    TCanvas cvs( "cvs", "cvs", 800, 600 );
    gPad->SetRightMargin( 0.2 );
    cvs.SetGridx( 1 ); cvs.SetGridy( 1 );

    const Int_t NRGBs = 5; const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 1.00 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.80, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable( NRGBs, stops, red, green, blue, NCont );
    gStyle->SetNumberContours( NCont );

    for( auto histPair : m_histTH2F ) {
        if( histPair.second == nullptr ) continue;
        histPair.second->Draw( option.c_str( ) );
        TLine lineCorr( histPair.second->GetXaxis( )->GetXmin( ), histPair.second->GetYaxis( )->GetXmin( ),
                        histPair.second->GetXaxis( )->GetXmax( ), histPair.second->GetYaxis( )->GetXmax( ) );
        lineCorr.SetLineStyle( 2 );
        // lineCorr.Draw( );
        // ATLASLabel( 0.18, 0.88, "Work in Progress", 1, 0.05 );
        // ShTUtil::CreateDrawText( 0.18, 0.8, "#sqrt{s} = 13 TeV, #intLdt = 79.8 fb^{-1}", 0.05 );
        ShTUtil::CreateDrawText( 0.18, 0.88, "#sqrt{s} = 13 TeV", 0.05 );
        cvs.SaveAs( Form( "%s/%s.pdf", outputDir.c_str( ), histPair.first.c_str( ) ) );
        cvs.SaveAs( Form( "%s/%s.png", outputDir.c_str( ), histPair.first.c_str( ) ) );
        cvs.SaveAs( Form( "%s/%s.eps", outputDir.c_str( ), histPair.first.c_str( ) ) );
    }

    return true;
}

bool
ShTHistManager::saveTH3F( const String& outputDir, const String& option )
{
    ShUtil::ExistCreateDir( outputDir );
    TCanvas cvs( "cvs", "cvs", 800, 600 );
    for( auto histPair : m_histTH3F ) {
        if( histPair.second != nullptr )
            histPair.second->Draw( option.c_str( ) );
        cvs.SaveAs( Form( "%s/%s.pdf", outputDir.c_str( ), histPair.first.c_str( ) ) );
        cvs.SaveAs( Form( "%s/%s.png", outputDir.c_str( ), histPair.first.c_str( ) ) );
        cvs.SaveAs( Form( "%s/%s.eps", outputDir.c_str( ), histPair.first.c_str( ) ) );
    }

    return true;
}

bool
ShTHistManager::saveTProfile( const String& outputDir, const String& option )
{
    ShUtil::ExistCreateDir( outputDir );
    TCanvas cvs( "cvs", "cvs", 800, 600 );
    for( auto histPair : m_histTProfile ) {
        if( histPair.second != nullptr )
            histPair.second->Draw( option.c_str( ) );
        cvs.SaveAs( Form( "%s/%s.pdf", outputDir.c_str( ), histPair.first.c_str( ) ) );
        cvs.SaveAs( Form( "%s/%s.png", outputDir.c_str( ), histPair.first.c_str( ) ) );
        cvs.SaveAs( Form( "%s/%s.eps", outputDir.c_str( ), histPair.first.c_str( ) ) );
    }

    return true;
}

bool
ShTHistManager::saveTF1( const String& outputDir, const String& option )
{
    ShUtil::ExistCreateDir( outputDir );
    TCanvas cvs( "cvs", "cvs", 800, 600 );
    for( auto funcPair : m_funcTF1 ) {
        if( funcPair.second != nullptr )
            funcPair.second->Draw( option.c_str( ) );
        cvs.SaveAs( Form( "%s/%s.pdf", outputDir.c_str( ), funcPair.first.c_str( ) ) );
        cvs.SaveAs( Form( "%s/%s.png", outputDir.c_str( ), funcPair.first.c_str( ) ) );
        cvs.SaveAs( Form( "%s/%s.eps", outputDir.c_str( ), funcPair.first.c_str( ) ) );
    }

    return true;
}

bool
ShTHistManager::saveAll( const String& outputDir, const String& option )
{
    return ( saveTH1F( outputDir, option ) && 
             saveTH2F( outputDir, option ) && 
             saveTH3F( outputDir, option ) && 
             saveTProfile( outputDir, option ) );
}

bool
ShTHistManager::saveSameTH1F( const String&        outputDir,
                              const String&        keyword,
                              const String&        option,
                              const bool&          inheritColor,
                              const bool&          isLogY,
                              const bool&          isNorm,
                              const bool&          poissonErr,
                              const bool&          drawFitFunc,
                              const String&        signalName,
                              const String&        dataName )
{
    ShUtil::ExistCreateDir( outputDir );
    if( keyword.length( ) <= 0 ) return false;
    String keywordUnderBar = keyword + "_";
    double max = 0.0;
    double min = 100000.0;
    std::list< TH1F* > savePlots;
    for( auto histPair : m_histTH1F ) {
        String plotName = histPair.first;
        if( plotName.find( keywordUnderBar ) == String::npos ) continue;

        if( histPair.second == nullptr ) continue;
        TH1F* pHist = static_cast< TH1F* >( histPair.second->Clone( histPair.second->GetName( ) ) );
        if( plotName.find("Data") == String::npos )
            savePlots.push_front( histPair.second );
        else
            savePlots.push_back( histPair.second );
        // savePlots.push_back( pHist );
        double weight = 1.0;
        if( isNorm == true ) weight = histPair.second->GetSumOfWeights( );
        if( max < histPair.second->GetMaximum( ) / weight )
            max = histPair.second->GetMaximum( ) / weight;
        if( min > histPair.second->GetMinimum( ) / weight )
            min = histPair.second->GetMinimum( ) / weight;
    }
    // if( min < 0.0001 ) min = 0.0001;
    min = 0.008;


    if( savePlots.size( ) <= 0 ) {
        ShUtil::Cwarn( "no target plots in histmanager" );
        return true;
    }

    // tthyy specific (need to fix)
    int color = 3;
    TCanvas cvs( "cvs", "cvs", 800, 800 );
    if( isLogY == true ) cvs.SetLogy( 1 );
    TLegend* pLeg = ShTUtil::CreateLegend( );
    // pLeg->SetX1(0.6); pLeg->SetY1(0.65); pLeg->SetX2(0.93); pLeg->SetY2(0.92); 
    // pLeg->SetX1(0.25); pLeg->SetY1(0.4); pLeg->SetX2(0.75); pLeg->SetY2(0.7); 
    pLeg->SetX1(0.5); pLeg->SetY1(0.6); pLeg->SetX2(0.95); pLeg->SetY2(0.8); 
    double scale = 2.3;
    if( keyword.find("had1") != String::npos || keyword.find("had2") != String::npos ) scale = 2.5;
    if( keyword.find("had3") != String::npos || keyword.find("had4") != String::npos ) scale = 1.8;
    
    for( auto pHist : savePlots ) {
        TGraphAsymmErrors* pGraph = new TGraphAsymmErrors( );
        pGraph->SetLineWidth( 2 );
        pGraph->SetMarkerSize( 1 );
        if( pHist != nullptr ) {
            String histName = pHist->GetName( );
            if( isLogY == true )
                pHist->GetYaxis( )->SetRangeUser( min, max * 10.0 );
            else
                pHist->GetYaxis( )->SetRangeUser( 0.0, max * scale );
                
            // color
            if( inheritColor == true ) {
                pHist->SetMarkerColor( pHist->GetLineColor( ) );
                if( poissonErr == true ) pHist->SetMarkerSize( 0 );
                else                     pHist->SetMarkerSize( 1 );
                pHist->SetFillColor( 0 );
            }
            else {
                pHist->SetLineColor( color );
                pHist->SetMarkerColor( color );
                pHist->SetFillColor( 0 );

                if( histName.find( signalName ) != String::npos ) {
                    pHist->SetLineColor( kRed );
                    pHist->SetMarkerColor( kRed );
                }
                else if( histName.find( dataName ) != String::npos ) {
                    pHist->SetLineColor( kBlack );
                    pHist->SetMarkerColor( kBlack );
                }
                else {
                    ++color;
                }
            }

            if( isNorm == true ) pHist->GetYaxis( )->SetTitle( "A.U." );

            size_t pos = histName.find( keywordUnderBar );
            String sample = histName.substr( pos+keywordUnderBar.size(), histName.size( ) );
            double weight = isNorm ? pHist->GetSumOfWeights( ) : 1.0;
            if( histName.find( dataName.c_str( ) ) != String::npos ) {
                pHist->SetLineWidth( 1.0 );

                if( poissonErr == true ) {
                    pGraph->SetLineColor( pHist->GetLineColor( ) );
                    for( int i = 0; i < pHist->GetNbinsX( ); ++i ) {
                        ShTUtil::SetPointPoissonError( pGraph,
                                                       i,
                                                       pHist->GetBinCenter( i + 1 ),
                                                       pHist->GetBinContent( i + 1 ),
                                                       weight );
                        pHist->SetBinError( i + 1, 0.0001 );
                    }
                }
                if( isNorm == true ) pHist->Scale( 1.0 / weight );
                if( isLogY == true )
                    pHist->GetYaxis( )->SetRangeUser( min, max * 10.0 );
                else
                    pHist->GetYaxis( )->SetRangeUser( 0.0, max * scale );

                pHist->Draw( "same" );
                if( poissonErr == true ) pLeg->AddEntry( pGraph, sample.c_str( ), "lep" );
                else                     pLeg->AddEntry( pHist,  sample.c_str( ), "lep" );
            }
            else {
                if( isNorm == true ) pHist->Scale( 1.0 / weight );
                if( isLogY == true )
                    pHist->GetYaxis( )->SetRangeUser( min, max * 10.0 );
                else
                    pHist->GetYaxis( )->SetRangeUser( 0.0, max * scale );
                pHist->Draw( Form("same%s", option.c_str( ) ) );
                pLeg->AddEntry( pHist, sample.c_str( ), "f" );
            }
        }

        if( poissonErr == true ) pGraph->Draw("Psame");
    }

    pLeg->Draw();

    // ATLASLabel( 0.18, 0.88, "Work in Progress", 1, 0.05 );
    ShTUtil::CreateDrawText( 0.19, 0.86, "#sqrt{s} = 13 TeV, #intLdt = 79.8 fb^{-1}", 0.05 );
    // ShTUtil::CreateDrawText( 0.19, 0.78, "Leptonic region", 0.05 );
    // ShTUtil::CreateDrawText( 0.19, 0.78, "Hadronic region", 0.05 );
    // ShTUtil::CreateDrawText( 0.19, 0.78, "Leptonic category", 0.045 );
    // ShTUtil::CreateDrawText( 0.18, 0.83, "#sqrt{s} = 13 TeV, #intLdt = 79.8 fb^{-1}", 0.05 );

    String cat = "";
    if     ( keyword.find( "had1" ) != String::npos ) cat = "Hadronic Cat.1";
    else if( keyword.find( "had2" ) != String::npos ) cat = "Hadronic Cat.2";
    else if( keyword.find( "had3" ) != String::npos ) cat = "Hadronic Cat.3";
    else if( keyword.find( "had4" ) != String::npos ) cat = "Hadronic Cat.4";
    else if( keyword.find( "lep1" ) != String::npos ) cat = "Leptonic Cat.1";
    else if( keyword.find( "lep2" ) != String::npos ) cat = "Leptonic Cat.2";
    else if( keyword.find( "lep3" ) != String::npos ) cat = "Leptonic Cat.3";
    else if( keyword.find( "all"  ) != String::npos ) cat = "All Categories";
    ShTUtil::CreateDrawText( 0.19, 0.78, cat.c_str( ), 0.045 );

    if( keyword.find( "bdtscore_allTI" ) != String::npos ) {
        // TLine cut( 0, 0, 0, 0 );
        // double cutVal = -1.0;
        // if     ( keyword.find( "had" ) != String::npos ) cutVal = 0.911;
        // else if( keyword.find( "lep" ) != String::npos ) cutVal = 0.705;
        
        // cut.SetLineWidth( 2 );
        // cut.SetLineStyle( 2 );
        // cut.SetLineColor( kBlack );
        // cut.DrawLine( cutVal, 0.0, cutVal, max );

        // TArrow direction( cutVal, max * 4.0 / 5.0, cutVal + 0.5, max * 4.0 / 5.0, 0.05, "|>" );
        // direction.SetLineWidth( 2 );
        // direction.SetLineStyle( 2 );
        // direction.SetLineColor( kBlack );
        // direction.DrawArrow( cutVal, max * 4.0 / 5.0, cutVal + 0.04, max * 4.0 / 5.0, 0.015, "|>" );
    }
    
    if( drawFitFunc == true ) {
        TF1* SigFunc = getTF1( "Sig_" + keyword );
        TF1* BkgFunc = getTF1( "Bkg_" + keyword );
        TF1* TotFunc = getTF1( "Tot_" + keyword );

        if( SigFunc != nullptr && BkgFunc != nullptr && TotFunc != nullptr ) {
            SigFunc->Draw( "same" );
            BkgFunc->Draw( "same" );
            TotFunc->Draw( "same" );
        }
    }
    String filePathPdf = Form( "%s/%s.pdf", outputDir.c_str( ), keyword.c_str( ) );
    cvs.SaveAs( filePathPdf.c_str( ) );
    String filePathPng = Form( "%s/%s.png", outputDir.c_str( ), keyword.c_str( ) );
    cvs.SaveAs( filePathPng.c_str( ) );
    String filePathEps = Form( "%s/%s.eps", outputDir.c_str( ), keyword.c_str( ) );
    cvs.SaveAs( filePathEps.c_str( ) );
    if( isLogY == true ) cvs.SetLogy( 0 );

    return true;
}

bool
ShTHistManager::saveStackTH1F( const String&             outputDir,
                               const String&             keyword,
                               std::vector< String >*    pExcludeList,
                               const bool&               drawExcludeAsSame,
                               const String&             option,
                               const bool&               inheritColor )
{
    ShUtil::ExistCreateDir( outputDir );
    if( keyword.length( ) <= 0 ) return false;
    
    // ### note: pExcludeList could be nullptr!!! ###
    // if( pExcludeList == nullptr ) return false;

    std::list< TH1F* > savePlots;
    std::list< TH1F* > excludePlots;
    for( auto histPair : m_histTH1F ) {
        if( histPair.first.find( keyword ) == String::npos ) continue;
        if( histPair.second == nullptr ) continue;
        String plotName = histPair.second->GetName( );

        bool exclude = false;
        if( pExcludeList != nullptr && pExcludeList->size( ) > 0 ) {
            for( auto name : *pExcludeList )
                if( plotName.find( name ) != String::npos ) exclude = true;
        }

        if( exclude == true ) {
            if( drawExcludeAsSame == true ) excludePlots.push_back( histPair.second );
        }
        else {
            if( plotName.find("Data") == String::npos )
                savePlots.push_front( histPair.second );
            else
                savePlots.push_back( histPair.second );
        }
    }

    int color = 3;
    TCanvas cvs( "cvs", "cvs", 800, 600 );
    TLegend* pLeg = ShTUtil::CreateLegend( );
    pLeg->SetX1(0.2); pLeg->SetY1(0.7); pLeg->SetX2(0.8); pLeg->SetY2(0.9); 
    pLeg->SetNColumns( 2 );
    THStack hs( Form( "stack_%s", keyword.c_str( ) ), Form( "stack_%s", keyword.c_str( ) ) );
    String xtitle = "";
    String ytitle = "";
    bool setTitle = false;
    for( auto pHist: savePlots ){
        if( pHist == nullptr ) continue;

        if( inheritColor == true ) {
            pHist->SetMarkerColor( pHist->GetLineColor( ) );
            pHist->SetFillColor( pHist->GetLineColor( ) );
        }
        else {
            pHist->SetLineColor( color );
            pHist->SetMarkerColor( color );
            pHist->SetFillColor( color );
            ++color;
        }

        if( setTitle == false ) {
            xtitle = pHist->GetXaxis( )->GetTitle( );
            ytitle = pHist->GetYaxis( )->GetTitle( );
            setTitle = true;
        }

        hs.Add( pHist );

        String histName = pHist->GetName( );            
        size_t pos = histName.find( keyword );
        String sample = histName.substr( 0, pos-1 );
        pLeg->AddEntry( pHist, sample.c_str( ), "f" );
    }

    // tthyy specific (need to fix)
    double max = 0.0;
    if( drawExcludeAsSame == true && excludePlots.size( ) >0 ) {
        for( auto pHist : excludePlots ) {
            if( pHist == nullptr ) continue;

            String histName = pHist->GetName( );
            size_t pos = histName.find( keyword );
            String sample = histName.substr( 0, pos-1 );
            int colorSigBG = 0;
            if( histName.find( "tthyy" ) != String::npos )     colorSigBG = kRed;
            else if( histName.find( "data" ) != String::npos ) colorSigBG = kBlack;

            pHist->SetLineColor( colorSigBG );
            pHist->SetMarkerColor( colorSigBG );

            if( max < pHist->GetMaximum( ) )
                max = pHist->GetMaximum( );

            if( histName.find( "tthyy" ) != String::npos )
                pLeg->AddEntry( pHist, sample.c_str( ), "f" );
            else if( histName.find( "data" ) != String::npos )
                pLeg->AddEntry( pHist, sample.c_str( ), "lp" );

        }
    }

    if( max < hs.GetMaximum( ) )
        max = hs.GetMaximum( );
    
    hs.Draw( option.c_str( ) );
    ShTUtil::SetXTitle( &hs, xtitle.c_str( ) );
    ShTUtil::SetYTitle( &hs, ytitle.c_str( ) );
    // hs.SetMinimum( 0.0 );
    hs.SetMaximum( max * 2.0 );
    // if( hs.GetYaxis( ) != nullptr ) hs.GetYaxis( )->SetRangeUser( 0.0, max * 2.0 );
    hs.Draw( option.c_str( ) );
    if( drawExcludeAsSame == true && excludePlots.size( ) >0 ) {
        for( auto pHist : excludePlots ) {
            
            String histName = pHist->GetName( );
            if( histName.find( "tthyy" ) != String::npos ) {
                pHist->Draw("samehist");
            }
            else if( histName.find( "data" ) != String::npos ) {
                pHist->Draw("sameP");
            }
        }
    }

    pLeg->Draw();
    cvs.SaveAs( Form( "%s/stack_%s.pdf", outputDir.c_str( ), keyword.c_str( ) ) );
    cvs.SaveAs( Form( "%s/stack_%s.png", outputDir.c_str( ), keyword.c_str( ) ) );
    cvs.SaveAs( Form( "%s/stack_%s.eps", outputDir.c_str( ), keyword.c_str( ) ) );

    return true;
}



bool
ShTHistManager::saveToRootFile( TFile* pFile )
{
    if( pFile == nullptr ) return false;
    setDirectoryAll( pFile );

    for( auto histPair : m_histTH1F )
        if( histPair.second != nullptr ) histPair.second->Write( );
    for( auto histPair : m_histTH2F )
        if( histPair.second != nullptr ) histPair.second->Write( );
    for( auto histPair : m_histTH3F )
        if( histPair.second != nullptr ) histPair.second->Write( );
    for( auto histPair : m_histTProfile )
        if( histPair.second != nullptr ) histPair.second->Write( );

    setDirectoryAll( nullptr );
    return true;
}

void
ShTHistManager::setDirectoryTH1F( TDirectory* pDir )
{
    // ### note: pDir could be nullptr!!! ###
    // if( pDir == nullptr ) return false;
    for( auto histPair : m_histTH1F ) {
        if( histPair.second != nullptr )
            histPair.second->SetDirectory( pDir );
    }
    return;
}

void
ShTHistManager::setDirectoryTH2F( TDirectory* pDir )
{
    // ### note: pDir could be nullptr!!! ###
    // if( pDir == nullptr ) return false;
    for( auto histPair : m_histTH2F ) {
        if( histPair.second != nullptr )
            histPair.second->SetDirectory( pDir );
    }
    return;
}

void
ShTHistManager::setDirectoryTH3F( TDirectory* pDir )
{
    // ### note: pDir could be nullptr!!! ###
    // if( pDir == nullptr ) return false;
    for( auto histPair : m_histTH3F ) {
        if( histPair.second != nullptr )
            histPair.second->SetDirectory( pDir );
    }
    return;
}

void
ShTHistManager::setDirectoryTProfile( TDirectory* pDir )
{
    // ### note: pDir could be nullptr!!! ###
    // if( pDir == nullptr ) return false;
    for( auto histPair : m_histTProfile ) {
        if( histPair.second != nullptr )
            histPair.second->SetDirectory( pDir );
    }
    return;
}

void
ShTHistManager::setDirectoryAll( TDirectory* pDir )
{
    // ### note: pDir could be nullptr!!! ###
    // if( pDir == nullptr ) return false;
    setDirectoryTH1F( pDir );
    setDirectoryTH2F( pDir );
    setDirectoryTH3F( pDir );
    setDirectoryTProfile( pDir );
    return;
}

bool
ShTHistManager::normalizeTH1F( const bool& fitMax )
{
    bool retVal = true;
    for( auto histPair : m_histTH1F ) {
        if( normalizeTH1F( histPair.first, fitMax ) == false )
            retVal = false;
    }

    return true;
}

bool
ShTHistManager::normalizeTH1F( const String& name, const bool& fitMax )
{
    if( hasTH1F( name ) == false ) return false;
    TH1F* pHist = getTH1F( name );
    
    if( pHist == nullptr || pHist->GetEntries( ) <= 0 ) return false;
    double scaleVal = 1.0;
    if( fitMax == true )
        scaleVal = 1.0 / static_cast< double >( pHist->GetMaximum( ) );
    else
        scaleVal = 1.0 / static_cast< double >( pHist->Integral( ) );
    pHist->Scale( scaleVal );
    ShTUtil::SetYTitle( pHist, "A.U." );
    
    return true;
}

bool
ShTHistManager::normalizeTH1F( const double& normFactor )
{
    bool retVal = true;
    for( auto histPair : m_histTH1F ) {
        if( normalizeTH1F( histPair.first, normFactor ) == false )
            retVal = false;
    }

    return true;
}

bool
ShTHistManager::normalizeTH1F( const String& name, const double& normFactor )
{
    if( hasTH1F( name ) == false ) return false;
    TH1F* pHist = getTH1F( name );
    
    if( pHist == nullptr || pHist->GetEntries( ) <= 0 ) return false;
    if( normFactor < 0.0 ) return false;

    pHist->Scale( normFactor );
    ShTUtil::SetYTitle( pHist, "A.U." );
    return true;
}

bool
ShTHistManager::normalizeTH2F( const double& normFactor )
{
    bool retVal = true;
    for( auto histPair : m_histTH2F ) {
        if( normalizeTH2F( histPair.first, normFactor ) == false )
            retVal = false;
    }

    return true;
}

bool
ShTHistManager::normalizeTH2F( const String& name, const double& normFactor )
{
    if( hasTH2F( name ) == false ) return false;
    TH2F* pHist = getTH2F( name );
    
    if( pHist == nullptr || pHist->GetEntries( ) <= 0 ) return false;
    if( normFactor < 0.0 ) return false;

    pHist->Scale( normFactor );
    ShTUtil::SetZTitle( pHist, "A.U." );
    return true;
}

double
ShTHistManager::getNormFactorTH1F( const String& name, const bool& fitMax )
{
    if( hasTH1F( name ) == false ) return -1.0;
    TH1F* pHist = getTH1F( name );
    if( pHist == nullptr || pHist->GetEntries( ) <= 0 ) return -1.0;

    double scaleVal = 1.0;
    if( fitMax == true )
        scaleVal = 1.0 / static_cast< double >( pHist->GetMaximum( ) );
    else
        scaleVal = 1.0 / static_cast< double >( pHist->Integral( ) );
    pHist->Scale( scaleVal );
    ShTUtil::SetYTitle( pHist, "A.U." );
    
    return scaleVal;
}

TH1F*
ShTHistManager::createRatioPlot( TH1F* pNum, TH1F* pDen )
{
    if( pNum == nullptr || pDen == nullptr ) return nullptr;

    TH1F* pRatio = static_cast< TH1F* >( pNum->Clone( Form("Ratio_%s_%s", pNum->GetName( ), pDen->GetName( ) ) ) );
    if( pRatio == nullptr ) return nullptr;

    pRatio->Divide( pDen );

    return pRatio;
}

bool
ShTHistManager::saveRatioPlot( const String&        outputDir,
                               const String&        keyword,
                               const String&        numName,
                               const String&        denName,
                               const String&        option,
                               const bool&          isLogY )
{
    ShUtil::ExistCreateDir( outputDir );
    if( keyword.length( ) <= 0 ) return false;

    // check histograms
    if( hasTH1F( Form( "%s_%s", keyword.c_str( ), numName.c_str( ) ) ) == false ) {
        ShUtil::Cerr( "Numerator is not found" );
        return false;
    }
    if( hasTH1F( Form( "%s_%s", keyword.c_str( ), denName.c_str( ) ) ) == false ) {
        ShUtil::Cerr( "Denominator is not found" );
        return false;
    }

    TH1F* pNum = getTH1F( Form( "%s_%s", keyword.c_str( ), numName.c_str( ) ) );
    TH1F* pDen = getTH1F( Form( "%s_%s", keyword.c_str( ), denName.c_str( ) ) );
    double max = pNum->GetMaximum( ) > pDen->GetMaximum( ) ? pNum->GetMaximum( ) : pDen->GetMaximum( );
    double min = pNum->GetMinimum( ) > pDen->GetMinimum( ) ? pNum->GetMinimum( ) : pDen->GetMinimum( );
    if( min < 0.1 ) min = 0.1;

    TCanvas cvs( "cvs", "cvs", 800, 600 );
    cvs.Divide( 1, 2 );
    cvs.cd( 1 );
    gPad->SetPad( 0.0, 1.0, 1.0, 0.3 );
    gPad->SetBottomMargin( 0.02 );
    cvs.cd( 2 );
    gPad->SetPad( 0.0, 0.3, 1.0, 0.0 );
    gPad->SetTopMargin( 0.02 );
    gPad->SetBottomMargin( 0.45 );
    gPad->SetGridy( );

    cvs.cd( 1 );
    pNum->GetXaxis( )->SetLabelSize( 0.001 );
    pNum->SetMaximum( max * 2.0 );
    pNum->SetMinimum( 0.0 );
    if( isLogY == true ) {
        cvs.SetLogy( 1 );
        pNum->SetMinimum( min );
    }
    pNum->Draw( option.c_str( ) );
    pDen->Draw( ( option + "same" ).c_str( ) );
    ATLASLabel( 0.18, 0.8, "Work in Progress", 1, 0.05 );
    ShTUtil::CreateDrawText( 0.18, 0.7, "#sqrt{s} = 13 TeV, #intLdt = 79.8 fb^{-1}", 0.07 );

    cvs.cd( 2 );
    TH1F* pRatio = createRatioPlot( pNum, pDen );
    pRatio->GetXaxis( )->SetTitleSize( 0.12 );
    pRatio->GetXaxis( )->SetLabelSize( 0.12 );

    pRatio->GetYaxis( )->SetRangeUser( 0.5, 1.5 );
    pRatio->GetYaxis( )->SetTitleOffset( 0.55 );
    pRatio->GetYaxis( )->SetTitleSize( 0.12 );
    pRatio->GetYaxis( )->SetLabelSize( 0.12 );
    pRatio->GetYaxis( )->SetNdivisions( 505 );
    pRatio->Draw( option.c_str( ) );

    // String saveFormat = _getSaveOptionStr( saveAs );
    String filePathPdf = Form( "%s/ratio_%s.pdf", outputDir.c_str( ), keyword.c_str( ) );
    cvs.SaveAs( filePathPdf.c_str( ) );
    String filePathPng = Form( "%s/ratio_%s.png", outputDir.c_str( ), keyword.c_str( ) );
    cvs.SaveAs( filePathPng.c_str( ) );
    String filePathEps = Form( "%s/ratio_%s.eps", outputDir.c_str( ), keyword.c_str( ) );
    cvs.SaveAs( filePathEps.c_str( ) );

    
    return true;
}

bool
ShTHistManager::saveRatioPlot( const String&        outputDir,
                               const String&        keyword,
                               const String&        denName,
                               const String&        option,
                               const bool&          isLogY )
{
    ShUtil::ExistCreateDir( outputDir );
    if( keyword.length( ) <= 0 ) return false;

    // check histograms
    if( hasTH1F( Form( "%s_%s", keyword.c_str( ), denName.c_str( ) ) ) == false ) {
        ShUtil::Cerr( "Denominator is not found" );
        return false;
    }

    double max = 0.0;
    double min = 10000.0;
    String keywordUnderBar = keyword + "_";
    for( auto histPair : m_histTH1F ) {
        if( histPair.first.find( keywordUnderBar ) == String::npos ) continue;
        if( max < histPair.second->GetMaximum( ) )
            max = histPair.second->GetMaximum( );
        if( min > histPair.second->GetMinimum( ) )
            min = histPair.second->GetMinimum( );
    }
    if( min < 0.1 ) min = 0.1;
    
    TH1F* pDen = getTH1F( Form( "%s_%s", keyword.c_str( ), denName.c_str( ) ) );
    
    TCanvas cvs( "cvs", "cvs", 800, 600 );
    cvs.Divide( 1, 2 );

    TLegend* pLeg = ShTUtil::CreateLegend( );
    pLeg->SetX1( 0.7 ); pLeg->SetY1( 0.5 ); pLeg->SetX2( 0.92 ); pLeg->SetY2( 0.92 ); 

    cvs.cd( 1 );
    gPad->SetPad( 0.0, 1.0, 1.0, 0.3 );
    gPad->SetBottomMargin( 0.02 );
    cvs.cd( 2 );
    gPad->SetPad( 0.0, 0.3, 1.0, 0.0 );
    gPad->SetTopMargin( 0.02 );
    gPad->SetBottomMargin( 0.45 );
    gPad->SetGridy( );

    cvs.cd( 1 );
    pDen->GetXaxis( )->SetLabelSize( 0.001 );
    pDen->GetYaxis( )->SetRangeUser( 0.0, max * 2.0 );
    if( isLogY == true ) {
        cvs.SetLogy( 1 );
        pDen->GetYaxis( )->SetRangeUser( min, max * 10.0 );
    }
    pDen->Draw( option.c_str( ) );

    String histName = pDen->GetName( );
    size_t pos = histName.find( keywordUnderBar );
    String sample = histName.substr( pos + keywordUnderBar.size( ), histName.size( ) );
    pLeg->AddEntry( pDen, sample.c_str( ), "lp" );

    for( auto histPair : m_histTH1F ) {
        TH1F* pHist = histPair.second;
        if( pHist != nullptr ) {
            String name = histPair.first;
            if( name.find( keywordUnderBar ) != String::npos ) {
                if( name.find( denName ) != String::npos ) continue;
                pHist->Draw( ( option + "same" ).c_str( ) );
                histName = pHist->GetName( );
                pos = histName.find( keywordUnderBar );
                sample = histName.substr( pos + keywordUnderBar.size( ), histName.size( ) );
                pLeg->AddEntry( pHist, sample.c_str( ), "lp" );
            }
        }
    }
    if( isLogY == true ) cvs.SetLogy( 0 );
    pLeg->Draw( );
    ATLASLabel( 0.18, 0.85, "Internal", 1, 0.07 );
    ShTUtil::CreateDrawText( 0.18, 0.73, "#sqrt{s} = 13 TeV, #intLdt = 79.8 fb^{-1}", 0.07 );

    cvs.cd( 2 );
    TH1F* pRatio = createRatioPlot( pDen, pDen ); // ratio == 1
    pRatio->GetXaxis( )->SetTitleSize( 0.12 );
    pRatio->GetXaxis( )->SetLabelSize( 0.12 );

    // double ratioRangeMax = 1.0;
    // if( pRatio->GetMaximum( ) - 1.0 > 1.0 - pRatio->GetMinimum( ) )
    //     ratioRangeMax = pRatio->GetMaximum( ) - 1.0;
    // else
    //     ratioRangeMax = 1.0 - pRatio->GetMinimum( );
    // pRatio->GetYaxis( )->SetRangeUser( 1.0 - ratioRangeMax, 1.0 + ratioRangeMax );

    pRatio->GetYaxis( )->SetRangeUser( 0.75, 1.25 );
    pRatio->GetYaxis( )->SetTitleOffset( 0.55 );
    pRatio->GetYaxis( )->SetTitleSize( 0.12 );
    pRatio->GetYaxis( )->SetLabelSize( 0.12 );
    pRatio->GetYaxis( )->SetNdivisions( 505 );

    pRatio->SetLineWidth( 0.0 );
    pRatio->SetMarkerSize( 0.0 );

    pRatio->Draw( option.c_str( ) );
    for( auto histPair : m_histTH1F ) {
        TH1F* pHist = histPair.second;
        if( pHist != nullptr ) {
            String histName = histPair.first;
            if( histName.find( keywordUnderBar ) != String::npos ) {
                if( histName.find( denName ) != String::npos ) continue;
                if( createRatioPlot( pHist, pDen ) != nullptr )
                    createRatioPlot( pHist, pDen )->Draw( ( option + "same" ).c_str( ) );
            }
        }
    }
    
    String filePathPdf = Form( "%s/ratio_%s.pdf", outputDir.c_str( ), keyword.c_str( ) );
    cvs.SaveAs( filePathPdf.c_str( ) );
    String filePathPng = Form( "%s/ratio_%s.png", outputDir.c_str( ), keyword.c_str( ) );
    cvs.SaveAs( filePathPng.c_str( ) );
    String filePathEps = Form( "%s/ratio_%s.eps", outputDir.c_str( ), keyword.c_str( ) );
    cvs.SaveAs( filePathEps.c_str( ) );
    
    return true;
}
