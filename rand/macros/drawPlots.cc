#include "inc/shinclude.h"

void drawPlots( const String& inputFile, const String& outputDir )
{
    ShUtil::ExistCreateDir( outputDir );
    SetAtlasStyle( );

    TFile file( inputFile.c_str( ) );
    TTree* pTree = dynamic_cast< TTree* >( file.Get( "tree" ) );
    if( pTree == nullptr ) return;

    double dmM = 0.0;
    double dmInjV = 0.0;
    double dmInjTheta = 0.0, dmInjPhi = 0.0;

    double nuRecTheta = 0.0, nuRecPhi = 0.0;
    double nuRecE = 0.0;

    double formFactorSq = 0.0;
    double rndm = 0.0;

    double invWeight = 0.0, totalRateSI = 0.0, totalRateSD = 0.0;

    pTree->SetBranchAddress( "dmM",          &dmM );

    pTree->SetBranchAddress( "dmInjV",       &dmInjV );
    pTree->SetBranchAddress( "dmInjTheta",   &dmInjTheta );
    pTree->SetBranchAddress( "dmInjPhi",     &dmInjPhi );
    
    pTree->SetBranchAddress( "nuRecE",       &nuRecE );
    pTree->SetBranchAddress( "nuRecTheta",   &nuRecTheta );
    pTree->SetBranchAddress( "nuRecPhi",     &nuRecPhi );

    pTree->SetBranchAddress( "formFactorSq", &formFactorSq );
    pTree->SetBranchAddress( "rndm",         &rndm );

    pTree->SetBranchAddress( "invWeight",    &invWeight );
    pTree->SetBranchAddress( "totalRateSI",  &totalRateSI );
    pTree->SetBranchAddress( "totalRateSD",  &totalRateSD );


    TH2D* pHist2DDMDir    = new TH2D( "hist2DDMDir",    "hist2DDMDir",      100, -180, 180, 100, -90, 90 );
    // TH2D* pHist2DDMDirV   = new TH2D( "hist2DDMDirV",   "hist2DDMDirV",     100, -180, 180, 100, 0, 0.1 ); // v/c
    TH2D* pHist2DDMCosGCV = new TH2D( "hist2DDMCosGCV", "hist2D2DDMCosGCV", 100, -1.0, 1.0, 100, 0, 0.1 ); // v/c

    TH2D* pHist2DNuDir    = new TH2D( "hist2DNuDir",    "hist2DNuDir",      100, -180, 180, 100, -90, 90 );
    // TH2D* pHist2DNuDirE   = new TH2D( "hist2DNuDirE",   "hist2DNuDirE",     100, -180, 180, 100, 0, 0.01 ); // GeV
    TH2D* pHist2DNuCosGCE = new TH2D( "hist2DNuCosGCE", "hist2DNuCosGCE",   100, -1.0, 1.0, 100, 0.0, 0.1 ); // MeV
    // TH2D* pHist2DNuCosGCE = new TH2D( "hist2DNuCosGCE", "hist2DNuCosGCE",   100, -1.0, 1.0, 100, 0.0, 1.0 ); // GeV
    // TH2D* pHist2DNuDirE = new TH2D( "hist2DNuDirE", "hist2DNuDirE", 100, -180, 180, 100, 0, 1.0 ); // GeV

    TH1D* pHistNuCosGC_E100keV = new TH1D( "histNuCosGC_E100keV", "histNuCosGC_E100keV", 40, -1.0, 1.0 );
    TH1D* pHistNuCosGC_E1MeV   = new TH1D( "histNuCosGC_E1MeV",   "histNuCosGC_E1MeV",   40, -1.0, 1.0 );
    TH1D* pHistNuCosGC_E10MeV  = new TH1D( "histNuCosGC_E10MeV",  "histNuCosGC_E10MeV",  40, -1.0, 1.0 );
    TH1D* pHistNuCosGC_E100MeV = new TH1D( "histNuCosGC_E100MeV", "histNuCosGC_E100MeV", 40, -1.0, 1.0 );
    TH1D* pHistNuCosGC_EMore   = new TH1D( "histNuCosGC_EMore",   "histNuCosGC_EMore",   40, -1.0, 1.0 );

    TH1D* pHistNuCosGCRateSI_E100keV = new TH1D( "histNuCosGCRateSI_E100keV", "histNuCosGCRateSI_E100keV", 40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSI_E1MeV   = new TH1D( "histNuCosGCRateSI_E1MeV",   "histNuCosGCRateSI_E1MeV",   40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSI_E10MeV  = new TH1D( "histNuCosGCRateSI_E10MeV",  "histNuCosGCRateSI_E10MeV",  40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSI_E100MeV = new TH1D( "histNuCosGCRateSI_E100MeV", "histNuCosGCRateSI_E100MeV", 40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSI_EMore   = new TH1D( "histNuCosGCRateSI_EMore",   "histNuCosGCRateSI_EMore",   40, -1.0, 1.0 );

    TH1D* pHistNuCosGCRateSD_E100keV = new TH1D( "histNuCosGCRateSD_E100keV", "histNuCosGCRateSD_E100keV", 40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSD_E1MeV   = new TH1D( "histNuCosGCRateSD_E1MeV",   "histNuCosGCRateSD_E1MeV",   40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSD_E10MeV  = new TH1D( "histNuCosGCRateSD_E10MeV",  "histNuCosGCRateSD_E10MeV",  40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSD_E100MeV = new TH1D( "histNuCosGCRateSD_E100MeV", "histNuCosGCRateSD_E100MeV", 40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSD_EMore   = new TH1D( "histNuCosGCRateSD_EMore",   "histNuCosGCRateSD_EMore",   40, -1.0, 1.0 );

    int totEvt = pTree->GetEntries( );
    for( int evt = 0; evt < totEvt; ++evt ) {
        ShUtil::PrintProgressBar( evt, totEvt );
        pTree->GetEntry( evt );

        double dmPhiCorr = ( dmInjPhi > TMath::Pi( ) ) ? dmInjPhi - 2.0 * TMath::Pi( ) : dmInjPhi;
        double nuPhiCorr = ( nuRecPhi > TMath::Pi( ) ) ? nuRecPhi - 2.0 * TMath::Pi( ) : nuRecPhi;

        pHist2DDMDir->Fill( dmPhiCorr * 180.0 / TMath::Pi( ), dmInjTheta * 180.0 / TMath::Pi( ) );
        // pHist2DDMDirV->Fill( dmPhiCorr * 180.0 / TMath::Pi( ), dmInjV / ( TMath::C( ) * 0.001 ) );

        double dmCosGC = cos( dmInjTheta ) *cos( dmPhiCorr );
        pHist2DDMCosGCV->Fill( dmCosGC, dmInjV / ( TMath::C( ) * 0.001 ) );
        
        if( rndm <= formFactorSq ) {
            pHist2DNuDir->Fill( nuPhiCorr * 180.0 / TMath::Pi( ), nuRecTheta * 180.0 / TMath::Pi( ) );
            // pHist2DNuDirE->Fill( nuPhiCorr * 180.0 / TMath::Pi( ), nuRecE );

            double nuRecCosGC = cos( nuRecTheta ) *cos( nuPhiCorr );
            // if( nuRecE > 0.01 ) pHist2DNuCosGCE->Fill( nuRecCosGC, nuRecE );
            pHist2DNuCosGCE->Fill( nuRecCosGC, nuRecE*1000.0 );

            if     ( nuRecE > 0.1     ) { 
                pHistNuCosGC_EMore->Fill( nuRecCosGC );
                pHistNuCosGCRateSI_EMore->Fill( nuRecCosGC, invWeight * totalRateSI / (double)totEvt ); 
                pHistNuCosGCRateSD_EMore->Fill( nuRecCosGC, invWeight * totalRateSD / (double)totEvt ); 
            }
            else if( nuRecE > 0.01    ) {
                pHistNuCosGC_E100MeV->Fill( nuRecCosGC );
                pHistNuCosGCRateSI_E100MeV->Fill( nuRecCosGC, invWeight * totalRateSI / (double)totEvt );
                pHistNuCosGCRateSD_E100MeV->Fill( nuRecCosGC, invWeight * totalRateSD / (double)totEvt );
            }
            else if( nuRecE > 0.001   ) {
                pHistNuCosGC_E10MeV->Fill( nuRecCosGC );
                pHistNuCosGCRateSI_E10MeV->Fill( nuRecCosGC, invWeight * totalRateSI / (double)totEvt ); 
                pHistNuCosGCRateSD_E10MeV->Fill( nuRecCosGC, invWeight * totalRateSD / (double)totEvt ); 
            }
            else if( nuRecE > 0.0001  ) { 
                pHistNuCosGC_E1MeV->Fill( nuRecCosGC );
                pHistNuCosGCRateSI_E1MeV->Fill( nuRecCosGC, invWeight * totalRateSI / (double)totEvt );
                pHistNuCosGCRateSD_E1MeV->Fill( nuRecCosGC, invWeight * totalRateSD / (double)totEvt );
            }
            else if( nuRecE > 0.00001 ) { 
                pHistNuCosGC_E100keV->Fill( nuRecCosGC );
                pHistNuCosGCRateSI_E100keV->Fill( nuRecCosGC, invWeight * totalRateSI / (double)totEvt );
                pHistNuCosGCRateSD_E100keV->Fill( nuRecCosGC, invWeight * totalRateSD / (double)totEvt );
            }
        }
    }

    TCanvas cvs( "cvs", "cvs", 800, 600 );

    if( pHistNuCosGC_EMore  ->Integral( ) > 0.0 ) pHistNuCosGC_EMore  ->Scale( 1.0 / pHistNuCosGC_EMore  ->Integral( ) );
    if( pHistNuCosGC_E100MeV->Integral( ) > 0.0 ) pHistNuCosGC_E100MeV->Scale( 1.0 / pHistNuCosGC_E100MeV->Integral( ) );
    if( pHistNuCosGC_E10MeV ->Integral( ) > 0.0 ) pHistNuCosGC_E10MeV ->Scale( 1.0 / pHistNuCosGC_E10MeV ->Integral( ) );
    if( pHistNuCosGC_E1MeV  ->Integral( ) > 0.0 ) pHistNuCosGC_E1MeV  ->Scale( 1.0 / pHistNuCosGC_E1MeV  ->Integral( ) );
    if( pHistNuCosGC_E100keV->Integral( ) > 0.0 ) pHistNuCosGC_E100keV->Scale( 1.0 / pHistNuCosGC_E100keV->Integral( ) );
    
    pHistNuCosGC_EMore  ->SetLineColor( kGreen+1 );
    pHistNuCosGC_E100MeV->SetLineColor( kBlue    );
    pHistNuCosGC_E10MeV ->SetLineColor( kViolet  );
    pHistNuCosGC_E1MeV  ->SetLineColor( kRed     );
    pHistNuCosGC_E100keV->SetLineColor( kBlack   );

    double max_histNuCosGC = 0.0;
    if( pHistNuCosGC_EMore->GetMaximum( )   > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGC_EMore->GetMaximum( );
    if( pHistNuCosGC_E100MeV->GetMaximum( ) > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGC_E100MeV->GetMaximum( );
    if( pHistNuCosGC_E10MeV->GetMaximum( )  > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGC_E10MeV->GetMaximum( );
    if( pHistNuCosGC_E1MeV->GetMaximum( )   > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGC_E1MeV->GetMaximum( );
    if( pHistNuCosGC_E100keV->GetMaximum( ) > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGC_E100keV->GetMaximum( );

    pHistNuCosGC_EMore  ->GetYaxis( )->SetRangeUser( 0.0, max_histNuCosGC * 1.2 );
    pHistNuCosGC_EMore  ->GetXaxis( )->SetTitle( "cos#it{#gamma}_{GC}" );
    pHistNuCosGC_EMore  ->GetYaxis( )->SetTitle( "A.U." );

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


    pHistNuCosGCRateSI_EMore  ->SetLineColor( kGreen+1 );
    pHistNuCosGCRateSI_E100MeV->SetLineColor( kBlue    );
    pHistNuCosGCRateSI_E10MeV ->SetLineColor( kViolet  );
    pHistNuCosGCRateSI_E1MeV  ->SetLineColor( kRed     );
    pHistNuCosGCRateSI_E100keV->SetLineColor( kBlack   );

    max_histNuCosGC = 0.0;
    if( pHistNuCosGCRateSI_EMore->GetMaximum( )   > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGCRateSI_EMore->GetMaximum( );
    if( pHistNuCosGCRateSI_E100MeV->GetMaximum( ) > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGCRateSI_E100MeV->GetMaximum( );
    if( pHistNuCosGCRateSI_E10MeV->GetMaximum( )  > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGCRateSI_E10MeV->GetMaximum( );
    if( pHistNuCosGCRateSI_E1MeV->GetMaximum( )   > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGCRateSI_E1MeV->GetMaximum( );
    if( pHistNuCosGCRateSI_E100keV->GetMaximum( ) > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGCRateSI_E100keV->GetMaximum( );

    pHistNuCosGCRateSI_EMore  ->GetYaxis( )->SetRangeUser( 0.0, max_histNuCosGC * 1.2 );
    pHistNuCosGCRateSI_EMore  ->GetXaxis( )->SetTitle( "cos#it{#gamma}_{GC}" );
    pHistNuCosGCRateSI_EMore  ->GetYaxis( )->SetTitle( "Events / sec / kg" );

    pHistNuCosGCRateSI_EMore  ->Draw( "hist" );
    pHistNuCosGCRateSI_E100MeV->Draw( "histsame" );
    pHistNuCosGCRateSI_E10MeV ->Draw( "histsame" );
    pHistNuCosGCRateSI_E1MeV  ->Draw( "histsame" );
    pHistNuCosGCRateSI_E100keV->Draw( "histsame" );
    pLeg->Draw( );
    
    ShTUtil::CreateDrawText( 0.225, 0.53, "#sigma_{#chi - p} = 10^{-30} cm^{2}, Spin Independent" );
    ShTUtil::CreateDrawText( 0.225, 0.45, Form("#it{m}_{#chi} = %1.1lf GeV", dmM) );
    cvs.SaveAs( Form( "%s/cosGCSI.png", outputDir.c_str( ) ) );


    pHistNuCosGCRateSD_EMore  ->SetLineColor( kGreen+1 );
    pHistNuCosGCRateSD_E100MeV->SetLineColor( kBlue    );
    pHistNuCosGCRateSD_E10MeV ->SetLineColor( kViolet  );
    pHistNuCosGCRateSD_E1MeV  ->SetLineColor( kRed     );
    pHistNuCosGCRateSD_E100keV->SetLineColor( kBlack   );

    max_histNuCosGC = 0.0;
    if( pHistNuCosGCRateSD_EMore->GetMaximum( )   > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGCRateSD_EMore->GetMaximum( );
    if( pHistNuCosGCRateSD_E100MeV->GetMaximum( ) > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGCRateSD_E100MeV->GetMaximum( );
    if( pHistNuCosGCRateSD_E10MeV->GetMaximum( )  > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGCRateSD_E10MeV->GetMaximum( );
    if( pHistNuCosGCRateSD_E1MeV->GetMaximum( )   > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGCRateSD_E1MeV->GetMaximum( );
    if( pHistNuCosGCRateSD_E100keV->GetMaximum( ) > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGCRateSD_E100keV->GetMaximum( );

    pHistNuCosGCRateSD_EMore  ->GetYaxis( )->SetRangeUser( 0.0, max_histNuCosGC * 1.2 );
    pHistNuCosGCRateSD_EMore  ->GetXaxis( )->SetTitle( "cos#it{#gamma}_{GC}" );
    pHistNuCosGCRateSD_EMore  ->GetYaxis( )->SetTitle( "Events / sec / kg" );

    pHistNuCosGCRateSD_EMore  ->Draw( "hist" );
    pHistNuCosGCRateSD_E100MeV->Draw( "histsame" );
    pHistNuCosGCRateSD_E10MeV ->Draw( "histsame" );
    pHistNuCosGCRateSD_E1MeV  ->Draw( "histsame" );
    pHistNuCosGCRateSD_E100keV->Draw( "histsame" );
    pLeg->Draw( );

    ShTUtil::CreateDrawText( 0.225, 0.53, "#sigma_{#chi - p} = 10^{-30} cm^{2}, Spin Dependent" );
    ShTUtil::CreateDrawText( 0.225, 0.45, Form("#it{m}_{#chi} = %1.1lf GeV", dmM) );
    cvs.SaveAs( Form( "%s/cosGCSD.png", outputDir.c_str( ) ) );


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

    // pHist2DDMDirV->GetXaxis( )->SetTitle( "#it{#phi}_{DM}" );
    // pHist2DDMDirV->GetYaxis( )->SetTitle( "#it{v}_{DM} / #it{c}" );
    // pHist2DDMDirV->GetZaxis( )->SetTitle( "A.U." );
    // pHist2DDMDirV->Draw( "colz" );
    // cvs.SaveAs( Form( "%s/dmDirV.png", outputDir.c_str( ) ) );

    cvs.SetLogz( 1 );
    // pHist2DNuDirE->GetXaxis( )->SetTitle( "#it{#phi}_{N}" );
    // pHist2DNuDirE->GetYaxis( )->SetTitle( "#it{E}_{N} [GeV]" );
    // pHist2DNuDirE->GetZaxis( )->SetTitle( "A.U." );
    // pHist2DNuDirE->Draw( "colz" );
    // cvs.SaveAs( Form( "%s/nuDirE.png", outputDir.c_str( ) ) );

    pHist2DDMCosGCV->GetXaxis( )->SetTitle( "cos#it{#gamma}_{GC}" );
    pHist2DDMCosGCV->GetYaxis( )->SetTitle( "#it{v}_{DM} / #it{c}" );
    pHist2DDMCosGCV->GetZaxis( )->SetTitle( "A.U." );
    pHist2DDMCosGCV->Draw( "colz" );
    cvs.SaveAs( Form( "%s/dmCosGCV.png", outputDir.c_str( ) ) );

    pHist2DNuCosGCE->GetXaxis( )->SetTitle( "cos#it{#gamma}_{GC}" );
    pHist2DNuCosGCE->GetYaxis( )->SetTitle( "#it{E}_{N} [MeV]" );
    pHist2DNuCosGCE->GetZaxis( )->SetTitle( "A.U." );
    pHist2DNuCosGCE->Draw( "colz" );
    cvs.SaveAs( Form( "%s/nuCosGCE.png", outputDir.c_str( ) ) );

    return;
}

