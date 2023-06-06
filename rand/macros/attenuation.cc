#include "inc/shinclude.h"
#include "CRDMFunc.h"
#include "CRDMFunc.cc"

void attenuation( )
{
    SetAtlasStyle( );
    // double dmM = 0.1; // 100 MeV
    // double dmM = 0.01; // 10 MeV
    double dmM = 0.0001; // 100 keV
    double xsection = 1e-32;
    
    double dmTArr[12]     = { 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0 };
    double dmTCorrArr[12] = { };

    for( int idx = 0; idx < 12; ++idx ) {
        double dmT = dmTArr[idx];

        // mantle1 (2900 km)
        for( int i = 0; i < 2900; i++ ) {
            dmT = dmT + attenuate( dmM, dmT, 100000.0, xsection, false );
        }
        DEBUG( dmT );
    
        // core (6940 km)
        for( int i = 0; i < 6940; i++ ) {
            dmT = dmT + attenuate( dmM, dmT, 100000.0, xsection, true );
        }
        DEBUG( dmT );

        // mantle2 (2900 km)
        for( int i = 0; i < 2900; i++ ) {
            dmT = dmT + attenuate( dmM, dmT, 100000.0, xsection, false );
        }
        dmTCorrArr[idx] = dmT;
        DEBUG( dmT );
    }

    TGraph g( 12, dmTArr, dmTCorrArr );
    TCanvas cvs( "cvs", "cvs", 800, 600 );
    cvs.SetGridx( 1 );
    cvs.SetGridy( 1 );
    cvs.SetLogx( 1 );
    cvs.SetLogy( 1 );
    g.Draw( "AP" );
    g.GetXaxis()->SetTitle( "Initial T_{#chi} [GeV]" );
    g.GetYaxis()->SetTitle( "Attenuated T_{#chi} [GeV]" );

    if     ( dmM < 0.001 ) ShTUtil::CreateDrawText( 0.2, 0.7, Form("#it{m}_{#chi} = %1.1lf keV", dmM*1000000.0 ) );
    else if( dmM < 1.0   ) ShTUtil::CreateDrawText( 0.2, 0.7, Form("#it{m}_{#chi} = %1.1lf MeV", dmM*1000.0 ) );
    else                   ShTUtil::CreateDrawText( 0.2, 0.7, Form("#it{m}_{#chi} = %1.1lf GeV", dmM ) );

    cvs.SaveAs( "attenuation.pdf" );

    TH1F histT( "hT", "hT", 1000, 0.0, 0.1 );
    TH1F histV( "hV", "hV", 1000, 0.0, 1.0 );

    TFile file( "fv_20230214/fv_NFW_los1dm1.root" );
    // TFile file( "fv_20230214/fv_NFW_los8dm1.root" );
    // TFile file( "fv_20230214/fv_NFW_los8dm00001.root" );
    TTree* pTree = dynamic_cast< TTree* >( file.Get( "tree" ) );
    if( pTree == nullptr ) return;

    double tree_dmM = 0.0, tree_dmV = 0.0;
    pTree->SetBranchAddress( "dmM", &tree_dmM );
    pTree->SetBranchAddress( "velocity", &tree_dmV );
    pTree->GetEntry( 0 );
    DEBUG(tree_dmV);
    DEBUG(tree_dmM);
    int tot = pTree->GetEntries( );
    for( int i = 0; i < tot; ++i ) {
        ShUtil::PrintProgressBar( i, tot );
        pTree->GetEntry( i );
        double dmGamma = 1.0 / sqrt(1.0 - tree_dmV*tree_dmV / V_LIGHT/V_LIGHT );
        double dmMom   = dmGamma * tree_dmM * tree_dmV / V_LIGHT;
        double dmE     = sqrt( dmMom * dmMom + tree_dmM * tree_dmM );
        double dmT     = dmE - tree_dmM;

        dmT = tree_dmM * (dmGamma - 1.0);
        // dmT = 0.5 * tree_dmM * tree_dmV*tree_dmV / V_LIGHT/V_LIGHT;

        histT.Fill( dmT, dmT );
        histV.Fill( tree_dmV / V_LIGHT, tree_dmV * tree_dmV );
    }

    histT.GetXaxis()->SetRangeUser( 1e-6, 0.1 );
    histT.GetXaxis()->SetLimits( 1e-5, 0.1 );
    histT.GetXaxis()->SetTitle( "Initial T_{#chi} [GeV]" );
    histT.GetYaxis()->SetTitle( "MC Events" );
    histT.Draw( );

    if     ( tree_dmM < 0.001 ) ShTUtil::CreateDrawText( 0.7, 0.8, Form("#it{m}_{#chi} = %1.1lf keV", tree_dmM*1000000.0 ) );
    else if( tree_dmM < 1.0   ) ShTUtil::CreateDrawText( 0.7, 0.8, Form("#it{m}_{#chi} = %1.1lf MeV", tree_dmM*1000.0 ) );
    else                        ShTUtil::CreateDrawText( 0.7, 0.8, Form("#it{m}_{#chi} = %1.1lf GeV", tree_dmM ) );

    cvs.SaveAs( "histT.pdf" );

    // histV.GetXaxis()->SetRangeUser( 1e-6, 0.1 );
    histV.GetXaxis()->SetLimits( 1e-1, 1.0 );
    histV.GetXaxis()->SetTitle( "#beta" );
    histV.GetYaxis()->SetTitle( "MC Events" );
    histV.Draw( );

    if     ( tree_dmM < 0.001 ) ShTUtil::CreateDrawText( 0.7, 0.8, Form("#it{m}_{#chi} = %1.1lf keV", tree_dmM*1000000.0 ) );
    else if( tree_dmM < 1.0   ) ShTUtil::CreateDrawText( 0.7, 0.8, Form("#it{m}_{#chi} = %1.1lf MeV", tree_dmM*1000.0 ) );
    else                        ShTUtil::CreateDrawText( 0.7, 0.8, Form("#it{m}_{#chi} = %1.1lf GeV", tree_dmM ) );

    cvs.SaveAs( "histV.pdf" );

    return;
}
