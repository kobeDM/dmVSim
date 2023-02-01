#include "inc/shinclude.h"

void compareFlux( const String& inputFileList, const String& outputDir )
{
    ShUtil::ExistCreateDir( outputDir );

    SetAtlasStyle( );

    std::list< String > fileList;
    if( ShUtil::GetLines( inputFileList, &fileList ) == false ) return;

    std::vector< TH1D* > histThetaArr;
    std::vector< TH1D* > histPhiArr;
    std::vector< TH1D* > histCosGammaArr;
    std::vector< int   > losIntArr;
    for( auto filename : fileList ) {
        TFile file( filename.c_str( ) );
        TTree* pTree = dynamic_cast< TTree* >( file.Get( "tree" ) );
        if( pTree == nullptr ) continue;

        double phi, theta = -100.0;
        double flux = 0.0;
        double los = 0.0;
        pTree->SetBranchAddress( "phi", &phi );
        pTree->SetBranchAddress( "theta", &theta );
        pTree->SetBranchAddress( "invWeight", &flux );
        pTree->SetBranchAddress( "los", &los );
        pTree->GetEntry( 0 );

        int losInt = static_cast< int >( los );

        TH1D* pHistTheta    = new TH1D( Form("histThetaLos%d",    losInt), Form("histThetaLos%d",    losInt), 100, -90,  90  );
        TH1D* pHistPhi      = new TH1D( Form("histPhiLos%d",      losInt), Form("histPhiLos%d",      losInt), 100, -180, 180 );
        TH1D* pHistCosGamma = new TH1D( Form("histCosGammaLos%d", losInt), Form("histCosGammaLos%d", losInt), 100, -1,   1   );
        pHistTheta->SetDirectory( nullptr );
        pHistPhi->SetDirectory( nullptr );
        pHistCosGamma->SetDirectory( nullptr );
        
        int totEvt = pTree->GetEntries( );
        for( int i = 0; i < totEvt; ++i ) {
            ShUtil::PrintProgressBar( i, totEvt );
            pTree->GetEntry( i );
            double phiCorr = ( phi > TMath::Pi( ) ) ? (phi - 2.0 * TMath::Pi( )) : phi;
            pHistTheta->Fill( theta * 180.0 / TMath::Pi( ), 1.0 / (double)totEvt );
            pHistPhi->Fill( phiCorr * 180.0 / TMath::Pi( ), 1.0 / (double)totEvt );
            pHistCosGamma->Fill( cos( theta ) *cos( phiCorr ), 1.0 / (double)totEvt );
        }
        
        histThetaArr.push_back( pHistTheta );
        histPhiArr.push_back( pHistPhi );
        histCosGammaArr.push_back( pHistCosGamma);
        losIntArr.push_back( losInt );
    }

    int arrSize = histThetaArr.size( );
    if( histPhiArr.size( )      != arrSize ||
        histCosGammaArr.size( ) != arrSize ||
        losIntArr.size( )       != arrSize ) {
        return;
    }

    // set histograms
    for( int idx = 0; idx < arrSize; ++idx ) {
        TH1D* pHistTheta    = histThetaArr.at( idx );
        TH1D* pHistPhi      = histPhiArr.at( idx );
        TH1D* pHistCosGamma = histCosGammaArr.at( idx );
        if( pHistTheta == nullptr || pHistPhi == nullptr || pHistCosGamma == nullptr ) continue;
        DEBUG("test1");
        
        pHistTheta->SetLineColor( idx + 1 );
        pHistTheta->SetXTitle( "#it{#theta}_{DM}" );
        pHistTheta->SetYTitle( "d#Phi/d#it{#theta}_{DM} [/cm^{2}/s]" );
        pHistTheta->GetYaxis( )->SetRangeUser( 0.0, pHistTheta->GetMaximum( ) * 1.5 );

        pHistPhi->SetLineColor( idx + 1 );
        pHistPhi->SetXTitle( "#it{#phi}_{DM}" );
        pHistPhi->SetYTitle( "d#Phi/d#it{#phi}_{DM} [/cm^{2}/s]" );
        pHistPhi->GetYaxis( )->SetRangeUser( 0.0, pHistPhi->GetMaximum( ) * 1.5 );

        pHistCosGamma->SetLineColor( idx + 1 );
        pHistCosGamma->SetXTitle( "cos#it{#gamma}_{GC}" );
        pHistCosGamma->SetYTitle( "d#Phi/dcos#it{#gamma}_{GC} [/cm^{2}/s]" );
        pHistCosGamma->GetYaxis( )->SetRangeUser( 0.0, pHistCosGamma->GetMaximum( ) * 1.5 );

    }
    DEBUG("test2");

    TCanvas cvs( "cvs","cvs", 800, 600 );
    TLegend* pLegTheta = ShTUtil::CreateLegend( 0.7, 0.7, 0.9, 0.9 );
    for( int idx = 0; idx < arrSize; ++idx ) {
        TH1D* pHistTheta  = histThetaArr.at( idx );
        int   losInt      = losIntArr.at( idx );
        if( pHistTheta == nullptr ) continue;
        
        if( idx == 0 ) pHistTheta->Draw( "hist" );
        else           pHistTheta->Draw( "samehist" );

        pLegTheta->AddEntry( pHistTheta, Form("l.o.s. = %d", losInt), "l" );
    }
    pLegTheta->Draw( );
    cvs.SaveAs( Form( "%s/fluxCompTheta.pdf", outputDir.c_str( ) ) );
    cvs.SaveAs( Form( "%s/fluxCompTheta.png", outputDir.c_str( ) ) );
    cvs.SaveAs( Form( "%s/fluxCompTheta.eps", outputDir.c_str( ) ) );

    TLegend* pLegPhi = ShTUtil::CreateLegend( 0.7, 0.7, 0.9, 0.9 );
    for( int idx = 0; idx < arrSize; ++idx ) {
        TH1D* pHistPhi  = histPhiArr.at( idx );
        int   losInt      = losIntArr.at( idx );
        if( pHistPhi == nullptr ) continue;
        
        if( idx == 0 ) pHistPhi->Draw( "hist" );
        else           pHistPhi->Draw( "samehist" );

        pLegPhi->AddEntry( pHistPhi, Form("l.o.s. = %d", losInt), "l" );
    }
    pLegPhi->Draw( );
    cvs.SaveAs( Form( "%s/fluxCompPhi.pdf", outputDir.c_str( ) ) );
    cvs.SaveAs( Form( "%s/fluxCompPhi.png", outputDir.c_str( ) ) );
    cvs.SaveAs( Form( "%s/fluxCompPhi.eps", outputDir.c_str( ) ) );

    TLegend* pLegCosGamma = ShTUtil::CreateLegend( 0.7, 0.7, 0.9, 0.9 );
    for( int idx = 0; idx < arrSize; ++idx ) {
        TH1D* pHistCosGamma  = histCosGammaArr.at( idx );
        int   losInt      = losIntArr.at( idx );
        if( pHistCosGamma == nullptr ) continue;
        
        if( idx == 0 ) pHistCosGamma->Draw( "hist" );
        else           pHistCosGamma->Draw( "samehist" );

        pLegCosGamma->AddEntry( pHistCosGamma, Form("l.o.s. = %d", losInt), "l" );
    }
    pLegCosGamma->Draw( );
    cvs.SaveAs( Form( "%s/fluxCompCosGamma.pdf", outputDir.c_str( ) ) );
    cvs.SaveAs( Form( "%s/fluxCompCosGamma.png", outputDir.c_str( ) ) );
    cvs.SaveAs( Form( "%s/fluxCompCosGamma.eps", outputDir.c_str( ) ) );

    return;
}
