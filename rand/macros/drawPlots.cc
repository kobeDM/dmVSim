#include "inc/shinclude.h"

void drawPlots( const String& inputFile, const String& outputDir )
{
    ShUtil::ExistCreateDir( outputDir );
    SetAtlasStyle( );

    TFile file( inputFile.c_str( ) );
    TTree* pTree = dynamic_cast< TTree* >( file.Get( "tree" ) );
    if( pTree == nullptr ) return;

    double dmInjV = 0.0;
    double dmInjTheta = 0.0, dmInjPhi = 0.0;

    double nuRecTheta = 0.0, nuRecPhi = 0.0;
    double nuRecE = 0.0;

    double formFactorSq = 0.0;
    double rndm = 0.0;

    pTree->SetBranchAddress( "dmInjV",       &dmInjV );
    pTree->SetBranchAddress( "dmInjTheta",   &dmInjTheta );
    pTree->SetBranchAddress( "dmInjPhi",     &dmInjPhi );
    
    pTree->SetBranchAddress( "nuRecE",       &nuRecE );
    pTree->SetBranchAddress( "nuRecTheta",   &nuRecTheta );
    pTree->SetBranchAddress( "nuRecPhi",     &nuRecPhi );

    pTree->SetBranchAddress( "formFactorSq", &formFactorSq );
    pTree->SetBranchAddress( "rndm",         &rndm );

    TH2D* pHist2DDMDir    = new TH2D( "hist2DDMDir",    "hist2DDMDir",      100, -180, 180, 100, -90, 90 );
    TH2D* pHist2DDMDirV   = new TH2D( "hist2DDMDirV",   "hist2DDMDirV",     100, -180, 180, 100, 0, 1.0 ); // v/c
    TH2D* pHist2DDMCosGCV = new TH2D( "hist2DDMCosGCV", "hist2D2DDMCosGCV", 100, -1.0, 1.0, 100, 0, 1.0 ); // v/c

    TH2D* pHist2DNuDir    = new TH2D( "hist2DNuDir",    "hist2DNuDir",      100, -180, 180, 100, -90, 90 );
    TH2D* pHist2DNuDirE   = new TH2D( "hist2DNuDirE",   "hist2DNuDirE",     100, -180, 180, 100, 0, 0.01 ); // GeV
    TH2D* pHist2DNuCosGCE = new TH2D( "hist2DNuCosGCE", "hist2DNuCosGCE",   100, -1.0, 1.0, 100, 0.0, 10.0 ); // GeV
    // TH2D* pHist2DNuDirE = new TH2D( "hist2DNuDirE", "hist2DNuDirE", 100, -180, 180, 100, 0, 1.0 ); // GeV

    TH1D* pHistNuCosGC_E100keV = new TH1D( "histNuCosGC_E100keV", "histNuCosGC_E100keV", 100, -1.0, 1.0 );
    TH1D* pHistNuCosGC_E1MeV   = new TH1D( "histNuCosGC_E1MeV",   "histNuCosGC_E1MeV",   100, -1.0, 1.0 );
    TH1D* pHistNuCosGC_E10MeV  = new TH1D( "histNuCosGC_E10MeV",  "histNuCosGC_E10MeV",  100, -1.0, 1.0 );
    TH1D* pHistNuCosGC_E100MeV = new TH1D( "histNuCosGC_E100MeV", "histNuCosGC_E100MeV", 100, -1.0, 1.0 );
    TH1D* pHistNuCosGC_EMore   = new TH1D( "histNuCosGC_EMore",   "histNuCosGC_EMore",   100, -1.0, 1.0 );

    int totEvt = pTree->GetEntries( );
    for( int evt = 0; evt < totEvt; ++evt ) {
        ShUtil::PrintProgressBar( evt, totEvt );
        pTree->GetEntry( evt );

        double dmPhiCorr = ( dmInjPhi > TMath::Pi( ) ) ? dmInjPhi - 2.0 * TMath::Pi( ) : dmInjPhi;
        double nuPhiCorr = ( nuRecPhi > TMath::Pi( ) ) ? nuRecPhi - 2.0 * TMath::Pi( ) : nuRecPhi;

        pHist2DDMDir->Fill( dmPhiCorr * 180.0 / TMath::Pi( ), dmInjTheta * 180.0 / TMath::Pi( ) );
        pHist2DDMDirV->Fill( dmPhiCorr * 180.0 / TMath::Pi( ), dmInjV / ( TMath::C( ) * 0.001 ) );

        double dmCosGC = cos( dmInjTheta ) *cos( dmPhiCorr );
        pHist2DDMCosGCV->Fill( dmCosGC, dmInjV / ( TMath::C( ) * 0.001 ) );
        
        if( rndm <= formFactorSq ) {
            pHist2DNuDir->Fill( nuPhiCorr * 180.0 / TMath::Pi( ), nuRecTheta * 180.0 / TMath::Pi( ) );
            pHist2DNuDirE->Fill( nuPhiCorr * 180.0 / TMath::Pi( ), nuRecE );

            double nuRecCosGC = cos( nuRecTheta ) *cos( nuPhiCorr );
            // if( nuRecE > 0.01 ) pHist2DNuCosGCE->Fill( nuRecCosGC, nuRecE );
            pHist2DNuCosGCE->Fill( nuRecCosGC, nuRecE );

            if     ( nuRecE > 0.1     ) pHistNuCosGC_EMore->Fill( nuRecCosGC );
            else if( nuRecE > 0.01    ) pHistNuCosGC_E100MeV->Fill( nuRecCosGC );
            else if( nuRecE > 0.001   ) pHistNuCosGC_E10MeV->Fill( nuRecCosGC );
            else if( nuRecE > 0.0001  ) pHistNuCosGC_E1MeV->Fill( nuRecCosGC );
            else if( nuRecE > 0.00001 ) pHistNuCosGC_E100keV->Fill( nuRecCosGC );
        }
    }

    TCanvas cvs( "cvs", "cvs", 800, 600 );

    pHistNuCosGC_EMore  ->Scale( 1.0 / pHistNuCosGC_EMore  ->Integral( ) );
    pHistNuCosGC_E100MeV->Scale( 1.0 / pHistNuCosGC_E100MeV->Integral( ) );
    pHistNuCosGC_E10MeV ->Scale( 1.0 / pHistNuCosGC_E10MeV ->Integral( ) );
    pHistNuCosGC_E1MeV  ->Scale( 1.0 / pHistNuCosGC_E1MeV  ->Integral( ) );
    pHistNuCosGC_E100keV->Scale( 1.0 / pHistNuCosGC_E100keV->Integral( ) );
    
    pHistNuCosGC_EMore  ->SetLineColor( kGreen+1 );
    pHistNuCosGC_E100MeV->SetLineColor( kBlue    );
    pHistNuCosGC_E10MeV ->SetLineColor( kViolet  );
    pHistNuCosGC_E1MeV  ->SetLineColor( kRed     );
    pHistNuCosGC_E100keV->SetLineColor( kBlack   );

    pHistNuCosGC_EMore  ->GetXaxis( )->SetTitle( "cos#it{#gamma}_{GC}" );
    pHistNuCosGC_EMore  ->GetYaxis( )->SetTitle( "A.U." );
    pHistNuCosGC_EMore  ->GetYaxis( )->SetRangeUser( 0.0, 0.1 );

    pHistNuCosGC_EMore  ->Draw( "hist" );
    pHistNuCosGC_E100MeV->Draw( "histsame" );
    pHistNuCosGC_E10MeV ->Draw( "histsame" );
    pHistNuCosGC_E1MeV  ->Draw( "histsame" );
    pHistNuCosGC_E100keV->Draw( "histsame" );

    TLegend* pLeg = ShTUtil::CreateLegend( 0.2, 0.6, 0.7, 0.9 );
    pLeg->AddEntry( pHistNuCosGC_EMore,   "100 MeV < E_{rec}",          "l" );
    pLeg->AddEntry( pHistNuCosGC_E100MeV, "10 MeV < E_{rec} < 100 MeV", "l" );
    pLeg->AddEntry( pHistNuCosGC_E10MeV,  "1 MeV < E_{rec} < 10 MeV",   "l" );
    pLeg->AddEntry( pHistNuCosGC_E1MeV,   "100 keV < E_{rec} < 1 MeV",  "l" );
    pLeg->AddEntry( pHistNuCosGC_E100keV, "10 keV < E_{rec} < 100 keV", "l" );
    pLeg->Draw( );

    cvs.SaveAs( Form( "%s/cosGC.png", outputDir.c_str( ) ) );

    cvs.SetGridx( 1 ); cvs.SetGridy( 1 );

    const Int_t NRGBs = 5; const Int_t NCont = 63;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 1.00 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.80, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable( NRGBs, stops, red, green, blue, NCont );
    gStyle->SetNumberContours( NCont );

    pHist2DDMDir->GetXaxis( )->SetTitle( "#it{#phi}_{DM}" );
    pHist2DDMDir->GetYaxis( )->SetTitle( "#it{#theta}_{DM}" );
    pHist2DDMDir->GetZaxis( )->SetTitle( "A.U." );
    pHist2DDMDir->Draw( "aitoff" );
    cvs.SaveAs( Form( "%s/dmDir.png", outputDir.c_str( ) ) );

    pHist2DNuDir->GetXaxis( )->SetTitle( "#it{#phi}_{N}" );
    pHist2DNuDir->GetYaxis( )->SetTitle( "#it{#theta}_{N}" );
    pHist2DNuDir->GetZaxis( )->SetTitle( "A.U." );
    pHist2DNuDir->Draw( "aitoff" );
    cvs.SaveAs( Form( "%s/nuDir.png", outputDir.c_str( ) ) );

    gPad->SetRightMargin( 0.2 );

    pHist2DDMDirV->GetXaxis( )->SetTitle( "#it{#phi}_{DM}" );
    pHist2DDMDirV->GetYaxis( )->SetTitle( "#it{v}_{DM} / #it{c}" );
    pHist2DDMDirV->GetZaxis( )->SetTitle( "A.U." );
    pHist2DDMDirV->Draw( "colz" );
    cvs.SaveAs( Form( "%s/dmDirV.png", outputDir.c_str( ) ) );

    cvs.SetLogz( 1 );
    pHist2DNuDirE->GetXaxis( )->SetTitle( "#it{#phi}_{N}" );
    pHist2DNuDirE->GetYaxis( )->SetTitle( "#it{E}_{N} [GeV]" );
    pHist2DNuDirE->GetZaxis( )->SetTitle( "A.U." );
    pHist2DNuDirE->Draw( "colz" );
    cvs.SaveAs( Form( "%s/nuDirE.png", outputDir.c_str( ) ) );

    pHist2DDMCosGCV->GetXaxis( )->SetTitle( "cos#it{#gamma}_{GC}" );
    pHist2DDMCosGCV->GetYaxis( )->SetTitle( "#it{v}_{DM} / #it{c}" );
    pHist2DDMCosGCV->GetZaxis( )->SetTitle( "A.U." );
    pHist2DDMCosGCV->Draw( "colz" );
    cvs.SaveAs( Form( "%s/dmCosGCV.png", outputDir.c_str( ) ) );

    pHist2DNuCosGCE->GetXaxis( )->SetTitle( "cos#it{#gamma}_{GC}" );
    pHist2DNuCosGCE->GetYaxis( )->SetTitle( "#it{E}_{N} [GeV]" );
    pHist2DNuCosGCE->GetZaxis( )->SetTitle( "A.U." );
    pHist2DNuCosGCE->Draw( "colz" );
    cvs.SaveAs( Form( "%s/nuCosGCE.png", outputDir.c_str( ) ) );


    return;
}

