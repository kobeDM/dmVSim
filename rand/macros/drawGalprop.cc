#include "inc/shinclude.h"

static std::vector< double > gEneAray;
static std::vector< double > gFlxAray;
static std::vector< double > gFlxE2Aray;

bool readGalprop( const String& input );

void drawGalprop( const String& inputFile, const String& outputDir )
{
    ShUtil::ExistCreateDir( outputDir );
    SetAtlasStyle( );
    
    if( readGalprop( inputFile ) == false ) return;
    
    double x[1000] = {};
    double y[1000] = {};
    double yE2[1000] = {};
    int arrSize = gEneAray.size( );
    for( int i = 0; i < arrSize; ++i ) {
        x[i] = gEneAray[i];
        y[i] = gFlxAray[i];
        yE2[i] = gFlxE2Aray[i];
    }

    TGraph gr( arrSize, x, y );
    TGraph grE2( arrSize, x, yE2 );
    TCanvas cvs( "cvs", "cvs", 800, 600 );
    cvs.SetLogx(1);
    cvs.SetLogy(1);
    gr.Draw( "AP" );
    cvs.SaveAs( Form( "%s/out.png", outputDir.c_str( ) ) );

    grE2.Draw( "AP" );
    cvs.SaveAs( Form( "%s/outE2.png", outputDir.c_str( ) ) );

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
        gEneAray.push_back( energy * 0.001);  // GeV
        gFlxAray.push_back( flux / energy / energy * 1000.0 ); // /cm^2/s/sr/GeV
        gFlxE2Aray.push_back( flux * 1000.0 ); // MeV^2/cm^2/s/sr/GeV
    }

    return true;
}
