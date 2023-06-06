#include "inc/shinclude.h"
#include "CRDMFunc.h"
#include "CRDMFunc.cc"

void recoilEnergy( )
{
    SetAtlasStyle( );

    TFile file( "recoil_20230417/recoil_NFW_los8dm001_Ag.root" );
    TTree* pTree = dynamic_cast< TTree* >( file.Get( "tree" ) );
    if( pTree == nullptr ) return;

    double dmInjV = 0.0, nuRecE = 0.0, dmM = 0.0;
    pTree->SetBranchAddress( "dmM",   &dmM    );
    pTree->SetBranchAddress( "dmInjV",&dmInjV );
    pTree->SetBranchAddress( "nuRecE",&nuRecE );

    TH2F hist2D( "hist", "hist", 100,0,200.0,100,0,10);
    int tot = pTree->GetEntries( );
    for( int i = 0; i < tot; ++i ) {
        ShUtil::PrintProgressBar( i, tot );
        pTree->GetEntry( i );
        double dmGamma = 1.0 / sqrt(1.0 - dmInjV*dmInjV / V_LIGHT/V_LIGHT );
        double dmMom   = dmGamma * dmM * dmInjV / V_LIGHT;
        double dmE     = sqrt( dmMom * dmMom + dmM * dmM );
        double dmT     = dmE - dmM;

        dmT = dmM * (dmGamma - 1.0);
        hist2D.Fill( dmT*1000.0, nuRecE*1000 ); // GeV -> MeV
    }

    TCanvas cvs( "cvs", "cvs", 800, 600 ); 
    const Int_t NRGBs = 5; const Int_t NCont = 63;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 1.00 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.80, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable( NRGBs, stops, red, green, blue, NCont );
    gStyle->SetNumberContours( NCont );
    gPad->SetRightMargin( 0.2 );
    cvs.SetLogz(1);

    hist2D.GetXaxis()->SetTitle( "T_{DM} [MeV]" );
    hist2D.GetYaxis()->SetTitle( "T_{N} [MeV]" );

    hist2D.Draw("colz");
    cvs.SaveAs( "energy.png" );

    return;
}
