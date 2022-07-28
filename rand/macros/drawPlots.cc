#include "inc/shinclude.h"

// const double EXPOSURE = 0.15 * 1.3e+06; // 
const double EXPOSURE = 0.155 * 0.78 * 1.0 * 31536000.0; // [SF6 density (20 Torr) : kg/m3] * [F occupancy for SF6] * volume [m3] * 1 year [sec]

const double CORRECTION = 1.0; // should not considered any corrections if our calculation is correct!!!

void drawPlots( const String& inputFile, const String& outputDir, const String& plotName )
{
    ShUtil::ExistCreateDir( outputDir );
    
    String subDir = outputDir + "/";
    subDir += plotName;
    ShUtil::ExistCreateDir( subDir );

    SetAtlasStyle( );

    TFile file( inputFile.c_str( ) );
    TTree* pTree = dynamic_cast< TTree* >( file.Get( "tree" ) );
    if( pTree == nullptr ) return;

    // find text file
    String txtFileName = outputDir + "/events.txt";
    DEBUG(txtFileName);
    std::ofstream ofs;
    if( ShUtil::ExistFile( txtFileName ) == true ) {
        ofs.open( txtFileName, std::ios_base::app );
    }
    else {
        ofs.open( txtFileName, std::ios_base::out  );
        ofs << "#plotname\tnTotalSI\tnPlusSI\tnFOMSI\tnFOMSIError\tnTotalSD\tnPlusSD\tnFOMSD\tnFOMSDError" << std::endl;
    }

    if( ofs.is_open( ) == false ) return;

    double dmM = 0.0;
    double dmInjV = 0.0;
    double dmInjTheta = 0.0, dmInjPhi = 0.0;

    double nuRecTheta = 0.0, nuRecPhi = 0.0;
    double nuRecCosTheta = 0.0;
    double nuRecE = 0.0;

    double formFactorSq = 0.0;
    double rndm = 0.0;

    double invWeight = 0.0, totalRateSI = 0.0, totalRateSD = 0.0;

    int atom = 0;

    pTree->SetBranchAddress( "dmM",           &dmM );

    pTree->SetBranchAddress( "dmInjV",        &dmInjV );
    pTree->SetBranchAddress( "dmInjTheta",    &dmInjTheta );
    pTree->SetBranchAddress( "dmInjPhi",      &dmInjPhi );
    
    pTree->SetBranchAddress( "nuRecE",        &nuRecE );
    pTree->SetBranchAddress( "nuRecTheta",    &nuRecTheta );
    pTree->SetBranchAddress( "nuRecCosTheta", &nuRecCosTheta );
    pTree->SetBranchAddress( "nuRecPhi",      &nuRecPhi );

    pTree->SetBranchAddress( "formFactorSq",  &formFactorSq );
    pTree->SetBranchAddress( "rndm",          &rndm );

    pTree->SetBranchAddress( "invWeight",     &invWeight );
    pTree->SetBranchAddress( "totalRateSI",   &totalRateSI );
    pTree->SetBranchAddress( "totalRateSD",   &totalRateSD );

    pTree->SetBranchAddress( "atom",          &atom );

    TH2D* pHist2DDMDir    = new TH2D( "hist2DDMDir",    "hist2DDMDir",      180, -180.0, 180.0, 179, -89.5, 89.5 );
    TH2D* pHist2DDMCosGCV = new TH2D( "hist2DDMCosGCV", "hist2D2DDMCosGCV", 180, -1.0, 1.0, 100, 0.0, 0.1 ); // v/c

    TH2D* pHist2DNuDir    = new TH2D( "hist2DNuDir",    "hist2DNuDir",      180, -180.0, 180.0, 179, -89.5, 89.5 );
    TH2D* pHist2DNuCosGCE = new TH2D( "hist2DNuCosGCE", "hist2DNuCosGCE",   100, -1.0, 1.0, 100, 0.0, 0.1 ); // MeV

    TH1D* pHistNuCosGC_E100keV = new TH1D( "histNuCosGC_E100keV", "histNuCosGC_E100keV", 40, -1.0, 1.0 );
    TH1D* pHistNuCosGC_E1MeV   = new TH1D( "histNuCosGC_E1MeV",   "histNuCosGC_E1MeV",   40, -1.0, 1.0 );
    TH1D* pHistNuCosGC_E10MeV  = new TH1D( "histNuCosGC_E10MeV",  "histNuCosGC_E10MeV",  40, -1.0, 1.0 );
    TH1D* pHistNuCosGC_E100MeV = new TH1D( "histNuCosGC_E100MeV", "histNuCosGC_E100MeV", 40, -1.0, 1.0 );
    TH1D* pHistNuCosGC_E1GeV   = new TH1D( "histNuCosGC_E1GeV",   "histNuCosGC_E1GeV",   40, -1.0, 1.0 );
    TH1D* pHistNuCosGC_EMore   = new TH1D( "histNuCosGC_EMore",   "histNuCosGC_EMore",   40, -1.0, 1.0 );

    TH1D* pHistNuCosGCRateSI_E100keV = new TH1D( "histNuCosGCRateSI_E100keV", "histNuCosGCRateSI_E100keV", 40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSI_E600keV = new TH1D( "histNuCosGCRateSI_E600keV", "histNuCosGCRateSI_E600keV", 40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSI_E1MeV   = new TH1D( "histNuCosGCRateSI_E1MeV",   "histNuCosGCRateSI_E1MeV",   40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSI_E3MeV   = new TH1D( "histNuCosGCRateSI_E3MeV",   "histNuCosGCRateSI_E3MeV",   40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSI_E10MeV  = new TH1D( "histNuCosGCRateSI_E10MeV",  "histNuCosGCRateSI_E10MeV",  40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSI_E100MeV = new TH1D( "histNuCosGCRateSI_E100MeV", "histNuCosGCRateSI_E100MeV", 40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSI_E1GeV   = new TH1D( "histNuCosGCRateSI_E1GeV",   "histNuCosGCRateSI_E1GeV",   40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSI_EMore   = new TH1D( "histNuCosGCRateSI_EMore",   "histNuCosGCRateSI_EMore",   40, -1.0, 1.0 );

    TH1D* pHistNuCosGCRateSD_E100keV = new TH1D( "histNuCosGCRateSD_E100keV", "histNuCosGCRateSD_E100keV", 40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSD_E600keV = new TH1D( "histNuCosGCRateSD_E600keV", "histNuCosGCRateSD_E600keV", 40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSD_E1MeV   = new TH1D( "histNuCosGCRateSD_E1MeV",   "histNuCosGCRateSD_E1MeV",   40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSD_E3MeV   = new TH1D( "histNuCosGCRateSD_E3MeV",   "histNuCosGCRateSD_E3MeV",   40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSD_E10MeV  = new TH1D( "histNuCosGCRateSD_E10MeV",  "histNuCosGCRateSD_E10MeV",  40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSD_E100MeV = new TH1D( "histNuCosGCRateSD_E100MeV", "histNuCosGCRateSD_E100MeV", 40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSD_E1GeV   = new TH1D( "histNuCosGCRateSD_E1GeV",   "histNuCosGCRateSD_E1GeV",   40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSD_EMore   = new TH1D( "histNuCosGCRateSD_EMore",   "histNuCosGCRateSD_EMore",   40, -1.0, 1.0 );

    TH1D* pHistNuE_600keV    = new TH1D( "histNuE_600keV",   "histNuE_600keV",   100, 0.0, 600 );
    TH1D* pHistNuE_3000keV   = new TH1D( "histNuE_3000keV",  "histNuE_3000keV",  100, 0.0, 3000 );
    TH1D* pHistNuE_10000keV  = new TH1D( "histNuE_10000keV", "histNuE_10000keV", 100, 0.0, 10000 );

    int totEvt = pTree->GetEntries( );
    for( int evt = 0; evt < totEvt; ++evt ) {
        ShUtil::PrintProgressBar( evt, totEvt );
        pTree->GetEntry( evt );

        double dmPhiCorr = ( dmInjPhi > TMath::Pi( ) ) ? dmInjPhi - 2.0 * TMath::Pi( ) : dmInjPhi;
        double nuPhiCorr = ( nuRecPhi > TMath::Pi( ) ) ? nuRecPhi - 2.0 * TMath::Pi( ) : nuRecPhi;

        pHist2DDMDir->Fill( dmPhiCorr * 180.0 / TMath::Pi( ), dmInjTheta * 180.0 / TMath::Pi( ) );

        double dmCosGC = cos( dmInjTheta ) *cos( dmPhiCorr );
        pHist2DDMCosGCV->Fill( dmCosGC, dmInjV / ( TMath::C( ) * 0.001 ) );
        
        if( rndm <= formFactorSq ) {
            if( nuRecCosTheta < 0.001 ) nuRecCosTheta = 0.001;
            pHist2DNuDir->Fill( nuPhiCorr * 180.0 / TMath::Pi( ), nuRecTheta * 180.0 / TMath::Pi( ), 1.0 / nuRecCosTheta ); // debug

            double nuRecCosGC = cos( nuRecTheta ) *cos( nuPhiCorr );

            if     ( nuRecE > 1.0     ) {  // > 1 GeV
                pHistNuCosGC_EMore->Fill( nuRecCosGC );
                pHistNuCosGCRateSI_EMore->Fill( nuRecCosGC, invWeight * totalRateSI / (double)totEvt * CORRECTION ); 
                pHistNuCosGCRateSD_EMore->Fill( nuRecCosGC, invWeight * totalRateSD / (double)totEvt * CORRECTION ); 
            }
            else if( nuRecE > 0.1     ) { // 100 MeV < m < 1 GeV
                pHistNuCosGC_E1GeV->Fill( nuRecCosGC );
                pHistNuCosGCRateSI_E1GeV->Fill( nuRecCosGC, invWeight * totalRateSI / (double)totEvt * CORRECTION ); 
                pHistNuCosGCRateSD_E1GeV->Fill( nuRecCosGC, invWeight * totalRateSD / (double)totEvt * CORRECTION ); 
            }
            else if( nuRecE > 0.01    ) { // 10 MeV < m < 100 MeV
                pHistNuCosGC_E100MeV->Fill( nuRecCosGC );
                pHistNuCosGCRateSI_E100MeV->Fill( nuRecCosGC, invWeight * totalRateSI / (double)totEvt * CORRECTION );
                pHistNuCosGCRateSD_E100MeV->Fill( nuRecCosGC, invWeight * totalRateSD / (double)totEvt * CORRECTION );
            }
            else if( nuRecE > 0.001   ) { // 1 MeV < m < 10 MeV
                pHistNuCosGC_E10MeV->Fill( nuRecCosGC );
                pHistNuCosGCRateSI_E10MeV->Fill( nuRecCosGC, invWeight * totalRateSI / (double)totEvt * CORRECTION ); 
                pHistNuCosGCRateSD_E10MeV->Fill( nuRecCosGC, invWeight * totalRateSD / (double)totEvt * CORRECTION ); 
            }
            else if( nuRecE > 0.0001  ) { // 100 keV < m < 1 MeV
                pHistNuCosGC_E1MeV->Fill( nuRecCosGC );
                pHistNuCosGCRateSI_E1MeV->Fill( nuRecCosGC, invWeight * totalRateSI / (double)totEvt * CORRECTION );
                pHistNuCosGCRateSD_E1MeV->Fill( nuRecCosGC, invWeight * totalRateSD / (double)totEvt * CORRECTION );
            }
            else if( nuRecE > 0.00001 ) { // 10 keV < m < 100 keV
                pHistNuCosGC_E100keV->Fill( nuRecCosGC );
                pHistNuCosGCRateSI_E100keV->Fill( nuRecCosGC, invWeight * totalRateSI / (double)totEvt * CORRECTION );
                pHistNuCosGCRateSD_E100keV->Fill( nuRecCosGC, invWeight * totalRateSD / (double)totEvt * CORRECTION );
            }

            if( nuRecE > 0.0001 && nuRecE < 0.0006 ) {
                pHistNuCosGCRateSI_E600keV->Fill( nuRecCosGC, invWeight * totalRateSI / (double)totEvt * CORRECTION );
                pHistNuCosGCRateSD_E600keV->Fill( nuRecCosGC, invWeight * totalRateSD / (double)totEvt * CORRECTION );
            }
            if( nuRecE > 0.001 && nuRecE < 0.003 ) {
                pHistNuCosGCRateSI_E3MeV->Fill( nuRecCosGC, invWeight * totalRateSI / (double)totEvt * CORRECTION );
                pHistNuCosGCRateSD_E3MeV->Fill( nuRecCosGC, invWeight * totalRateSD / (double)totEvt * CORRECTION );
            }

            if( nuRecE < 0.0006 ) pHistNuE_600keV->Fill  ( nuRecE * 1000000.0, invWeight * totalRateSI / (double)totEvt * CORRECTION );
            if( nuRecE < 0.003  ) pHistNuE_3000keV->Fill ( nuRecE * 1000000.0, invWeight * totalRateSI / (double)totEvt * CORRECTION );
            if( nuRecE < 0.01   ) pHistNuE_10000keV->Fill( nuRecE * 1000000.0, invWeight * totalRateSI / (double)totEvt * CORRECTION );
        }
    }

    String atomStr = "";
    if     ( atom == 10 ) atomStr = "F";
    else if( atom == 11 ) atomStr = "Ag";
    else if( atom == 12 ) atomStr = "p";

    TCanvas cvs( "cvs", "cvs", 800, 600 );

    if( pHistNuCosGC_EMore  ->Integral( ) > 0.0 ) pHistNuCosGC_EMore  ->Scale( 1.0 / pHistNuCosGC_EMore  ->Integral( ) );
    if( pHistNuCosGC_E1GeV  ->Integral( ) > 0.0 ) pHistNuCosGC_E1GeV  ->Scale( 1.0 / pHistNuCosGC_E1GeV  ->Integral( ) );
    if( pHistNuCosGC_E100MeV->Integral( ) > 0.0 ) pHistNuCosGC_E100MeV->Scale( 1.0 / pHistNuCosGC_E100MeV->Integral( ) );
    if( pHistNuCosGC_E10MeV ->Integral( ) > 0.0 ) pHistNuCosGC_E10MeV ->Scale( 1.0 / pHistNuCosGC_E10MeV ->Integral( ) );
    if( pHistNuCosGC_E1MeV  ->Integral( ) > 0.0 ) pHistNuCosGC_E1MeV  ->Scale( 1.0 / pHistNuCosGC_E1MeV  ->Integral( ) );
    if( pHistNuCosGC_E100keV->Integral( ) > 0.0 ) pHistNuCosGC_E100keV->Scale( 1.0 / pHistNuCosGC_E100keV->Integral( ) );
    
    pHistNuCosGC_EMore  ->SetLineColor( kCyan+1  );
    pHistNuCosGC_E1GeV  ->SetLineColor( kGreen+1 );
    pHistNuCosGC_E100MeV->SetLineColor( kBlue    );
    pHistNuCosGC_E10MeV ->SetLineColor( kViolet  );
    pHistNuCosGC_E1MeV  ->SetLineColor( kRed     );
    pHistNuCosGC_E100keV->SetLineColor( kBlack   );

    double max_histNuCosGC = 0.0;
    if( pHistNuCosGC_EMore->GetMaximum( )   > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGC_EMore->GetMaximum( );
    if( pHistNuCosGC_E1GeV->GetMaximum( )   > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGC_E1GeV->GetMaximum( );
    if( pHistNuCosGC_E100MeV->GetMaximum( ) > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGC_E100MeV->GetMaximum( );
    if( pHistNuCosGC_E10MeV->GetMaximum( )  > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGC_E10MeV->GetMaximum( );
    if( pHistNuCosGC_E1MeV->GetMaximum( )   > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGC_E1MeV->GetMaximum( );
    if( pHistNuCosGC_E100keV->GetMaximum( ) > max_histNuCosGC ) max_histNuCosGC = pHistNuCosGC_E100keV->GetMaximum( );

    pHistNuCosGC_EMore  ->GetYaxis( )->SetRangeUser( 0.0, max_histNuCosGC * 1.2 );
    pHistNuCosGC_EMore  ->GetXaxis( )->SetTitle( "cos#it{#gamma}_{GC}" );
    pHistNuCosGC_EMore  ->GetYaxis( )->SetTitle( "A.U." );

    pHistNuCosGC_EMore  ->Draw( "hist" );
    pHistNuCosGC_E1GeV  ->Draw( "histsame" );
    pHistNuCosGC_E100MeV->Draw( "histsame" );
    pHistNuCosGC_E10MeV ->Draw( "histsame" );
    pHistNuCosGC_E1MeV  ->Draw( "histsame" );
    pHistNuCosGC_E100keV->Draw( "histsame" );

    TLegend* pLeg = ShTUtil::CreateLegend( 0.2, 0.6, 0.7, 0.9 );
    pLeg->AddEntry( pHistNuCosGC_EMore,   "1 GeV < E_{rec}",            "l" );
    pLeg->AddEntry( pHistNuCosGC_E1GeV,   "100 MeV < E_{rec} < 1 GeV",  "l" );
    pLeg->AddEntry( pHistNuCosGC_E100MeV, "10 MeV < E_{rec} < 100 MeV", "l" );
    pLeg->AddEntry( pHistNuCosGC_E10MeV,  "1 MeV < E_{rec} < 10 MeV",   "l" );
    pLeg->AddEntry( pHistNuCosGC_E1MeV,   "100 keV < E_{rec} < 1 MeV",  "l" );
    pLeg->AddEntry( pHistNuCosGC_E100keV, "10 keV < E_{rec} < 100 keV", "l" );
    pLeg->Draw( );

    cvs.SaveAs( Form( "%s/%s_cosGC.png", subDir.c_str( ), plotName.c_str( ) ) );

    TH1D* pHistNuCosGCRateSI_tot = new TH1D( "histNuCosGCRateSI_tot",  "histNuCosGCRateSI_tot",  40, -1.0, 1.0 );
    TH1D* pHistNuCosGCRateSD_tot = new TH1D( "histNuCosGCRateSD_tot",  "histNuCosGCRateSD_tot",  40, -1.0, 1.0 );
    TLegend* pLeg_SISD = ShTUtil::CreateLegend( 0.2, 0.65, 0.7, 0.9 );
    if( atom == 10 ) { // F
        pHistNuCosGCRateSI_E100keV->SetLineColor( kRed   );
        pHistNuCosGCRateSI_E1MeV->SetLineColor( kViolet   );
        pHistNuCosGCRateSI_E3MeV->SetLineColor( kBlue   );
        pHistNuCosGCRateSI_tot->SetLineColor( kBlack   );
        
        pHistNuCosGCRateSI_tot->Add( pHistNuCosGCRateSI_E100keV );
        pHistNuCosGCRateSI_tot->Add( pHistNuCosGCRateSI_E1MeV   );
        pHistNuCosGCRateSI_tot->Add( pHistNuCosGCRateSI_E3MeV   );

        pHistNuCosGCRateSD_E100keV->SetLineColor( kRed   );
        pHistNuCosGCRateSD_E1MeV->SetLineColor( kViolet   );
        pHistNuCosGCRateSD_E3MeV->SetLineColor( kBlue   );
        pHistNuCosGCRateSD_tot->SetLineColor( kBlack   );
        
        pHistNuCosGCRateSD_tot->Add( pHistNuCosGCRateSD_E100keV );
        pHistNuCosGCRateSD_tot->Add( pHistNuCosGCRateSD_E1MeV   );
        pHistNuCosGCRateSD_tot->Add( pHistNuCosGCRateSD_E3MeV   );

        pLeg_SISD->AddEntry( pHistNuCosGCRateSI_E3MeV,   "1 MeV < E_{rec} < 3 MeV",    "l" );
        pLeg_SISD->AddEntry( pHistNuCosGCRateSI_E1MeV,   "100 keV < E_{rec} < 1 MeV",  "l" );
        pLeg_SISD->AddEntry( pHistNuCosGCRateSI_E100keV, "10 keV < E_{rec} < 100 keV", "l" );
    }
    else if( atom == 11 ) { // Ag
        pHistNuCosGCRateSI_E100keV->SetLineColor( kRed   );
        pHistNuCosGCRateSI_E600keV->SetLineColor( kViolet   );
        pHistNuCosGCRateSI_tot->SetLineColor( kBlack   );
        
        pHistNuCosGCRateSI_tot->Add( pHistNuCosGCRateSI_E100keV );
        pHistNuCosGCRateSI_tot->Add( pHistNuCosGCRateSI_E600keV );

        pHistNuCosGCRateSD_E100keV->SetLineColor( kRed   );
        pHistNuCosGCRateSD_E600keV->SetLineColor( kViolet   );
        pHistNuCosGCRateSD_tot->SetLineColor( kBlack   );
        
        pHistNuCosGCRateSD_tot->Add( pHistNuCosGCRateSD_E100keV );
        pHistNuCosGCRateSD_tot->Add( pHistNuCosGCRateSD_E600keV );

        pLeg_SISD->AddEntry( pHistNuCosGCRateSI_E600keV, "100 keV < E_{rec} < 600 keV", "l" );
        pLeg_SISD->AddEntry( pHistNuCosGCRateSI_E100keV, "10 keV < E_{rec} < 100 keV", "l" );
    }
    if( atom == 12 ) { // H
        pHistNuCosGCRateSI_E100keV->SetLineColor( kRed   );
        pHistNuCosGCRateSI_E1MeV->SetLineColor( kViolet   );
        pHistNuCosGCRateSI_E10MeV->SetLineColor( kBlue   );
        pHistNuCosGCRateSI_tot->SetLineColor( kBlack   );
        
        pHistNuCosGCRateSI_tot->Add( pHistNuCosGCRateSI_E100keV );
        pHistNuCosGCRateSI_tot->Add( pHistNuCosGCRateSI_E1MeV   );
        pHistNuCosGCRateSI_tot->Add( pHistNuCosGCRateSI_E10MeV   );

        pHistNuCosGCRateSD_E100keV->SetLineColor( kRed   );
        pHistNuCosGCRateSD_E1MeV->SetLineColor( kViolet   );
        pHistNuCosGCRateSD_E10MeV->SetLineColor( kBlue   );
        pHistNuCosGCRateSD_tot->SetLineColor( kBlack   );
        
        pHistNuCosGCRateSD_tot->Add( pHistNuCosGCRateSD_E100keV );
        pHistNuCosGCRateSD_tot->Add( pHistNuCosGCRateSD_E1MeV   );
        pHistNuCosGCRateSD_tot->Add( pHistNuCosGCRateSD_E10MeV   );

        pLeg_SISD->AddEntry( pHistNuCosGCRateSI_E10MeV,  "1 MeV < E_{rec} < 10 MeV",   "l" );
        pLeg_SISD->AddEntry( pHistNuCosGCRateSI_E1MeV,   "100 keV < E_{rec} < 1 MeV",  "l" );
        pLeg_SISD->AddEntry( pHistNuCosGCRateSI_E100keV, "10 keV < E_{rec} < 100 keV", "l" );
    }
    pLeg_SISD->AddEntry( pHistNuCosGCRateSI_tot,     "All",            "l" );

    max_histNuCosGC = pHistNuCosGCRateSI_tot->GetMaximum( );
    pHistNuCosGCRateSI_tot->GetXaxis( )->SetTitle( "cos#it{#gamma}_{GC}" );
    pHistNuCosGCRateSI_tot->GetYaxis( )->SetTitle( "Events / sec / kg" );
    pHistNuCosGCRateSI_tot->GetYaxis( )->SetRangeUser( 0.0, max_histNuCosGC * 1.4 );

    pHistNuCosGCRateSI_tot    ->Draw( "hist" );
    if( atom == 10 ) {
        pHistNuCosGCRateSI_E3MeV  ->Draw( "histsame" );
        pHistNuCosGCRateSI_E1MeV  ->Draw( "histsame" );
        pHistNuCosGCRateSI_E100keV->Draw( "histsame" );
    }
    else if( atom == 11 ) {
        pHistNuCosGCRateSI_E600keV->Draw( "histsame" );
        pHistNuCosGCRateSI_E100keV->Draw( "histsame" );
    }
    else if( atom == 12 ) {
        pHistNuCosGCRateSI_E10MeV ->Draw( "histsame" );
        pHistNuCosGCRateSI_E1MeV  ->Draw( "histsame" );
        pHistNuCosGCRateSI_E100keV->Draw( "histsame" );
    }
    pLeg_SISD->Draw( );
    
    ShTUtil::CreateDrawText( 0.225, 0.61, "#sigma_{#chi - p} = 10^{-32} cm^{2}, Spin Independent" );
    if     ( dmM < 0.001 ) ShTUtil::CreateDrawText( 0.225, 0.53, Form("#it{m}_{#chi} = %1.1lf keV, %s recoil", dmM*1000000.0, atomStr.c_str( ) ) );
    else if( dmM < 1.0   ) ShTUtil::CreateDrawText( 0.225, 0.53, Form("#it{m}_{#chi} = %1.1lf MeV, %s recoil", dmM*1000.0, atomStr.c_str( ) ) );
    else                   ShTUtil::CreateDrawText( 0.225, 0.53, Form("#it{m}_{#chi} = %1.1lf GeV, %s recoil", dmM, atomStr.c_str( ) ) );

    cvs.SaveAs( Form( "%s/%s_cosGCSI.png", subDir.c_str( ), plotName.c_str( ) ) );

    max_histNuCosGC = pHistNuCosGCRateSD_tot->GetMaximum( );
    pHistNuCosGCRateSD_tot->GetYaxis( )->SetRangeUser( 0.0, max_histNuCosGC * 1.4 );
    pHistNuCosGCRateSD_tot->GetXaxis( )->SetTitle( "cos#it{#gamma}_{GC}" );
    pHistNuCosGCRateSD_tot->GetYaxis( )->SetTitle( "Events / sec / kg" );

    pHistNuCosGCRateSD_tot    ->Draw( "hist" );
    if( atom == 10 ) {
        pHistNuCosGCRateSD_E3MeV  ->Draw( "histsame" );
        pHistNuCosGCRateSD_E1MeV  ->Draw( "histsame" );
        pHistNuCosGCRateSD_E100keV->Draw( "histsame" );
    }
    else if( atom == 11 ) {
        pHistNuCosGCRateSD_E600keV->Draw( "histsame" );
        pHistNuCosGCRateSD_E100keV->Draw( "histsame" );
    }
    else if( atom == 12 ) {
        pHistNuCosGCRateSD_E10MeV ->Draw( "histsame" );
        pHistNuCosGCRateSD_E1MeV  ->Draw( "histsame" );
        pHistNuCosGCRateSD_E100keV->Draw( "histsame" );
    }
    pLeg_SISD->Draw( );

    ShTUtil::CreateDrawText( 0.225, 0.61, "#sigma_{#chi - p} = 10^{-32} cm^{2}, Spin Dependent" );
    if     ( dmM < 0.001 ) ShTUtil::CreateDrawText( 0.225, 0.53, Form("#it{m}_{#chi} = %1.1lf keV, %s recoil", dmM*1000000.0, atomStr.c_str( ) ) );
    else if( dmM < 1.0   ) ShTUtil::CreateDrawText( 0.225, 0.53, Form("#it{m}_{#chi} = %1.1lf MeV, %s recoil", dmM*1000.0, atomStr.c_str( ) ) );
    else                   ShTUtil::CreateDrawText( 0.225, 0.53, Form("#it{m}_{#chi} = %1.1lf GeV, %s recoil", dmM, atomStr.c_str( ) ) );
    cvs.SaveAs( Form( "%s/%s_cosGCSD.png", subDir.c_str( ), plotName.c_str( ) ) );

    if( atom == 10 ) {
        pHistNuE_3000keV->GetXaxis()->SetTitle( "Recoil energy [keV]" );
        pHistNuE_3000keV->GetYaxis()->SetTitle( "Events / sec / kg" );
        pHistNuE_3000keV->Draw( "hist" );
    }
    else if( atom == 11 ) {
        pHistNuE_600keV->GetXaxis()->SetTitle( "Recoil energy [keV]" );
        pHistNuE_600keV->GetYaxis()->SetTitle( "Events / sec / kg" );
        pHistNuE_600keV->Draw( "hist" );
    }
    else if( atom == 12 ) {
        pHistNuE_10000keV->GetXaxis()->SetTitle( "Recoil energy [keV]" );
        pHistNuE_10000keV->GetYaxis()->SetTitle( "Events / sec / kg" );
        pHistNuE_10000keV->Draw( "hist" );
    }

    cvs.SaveAs( Form( "%s/%s_nuE.png", subDir.c_str( ), plotName.c_str( ) ) );

    // 2D plots
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
    cvs.SaveAs( Form( "%s/%s_dmDir.png", subDir.c_str( ), plotName.c_str( ) ) );

    pHist2DNuDir->GetXaxis( )->SetTitle( "#it{#phi}_{N}" );
    pHist2DNuDir->GetYaxis( )->SetTitle( "#it{#theta}_{N}" );
    pHist2DNuDir->GetZaxis( )->SetTitle( "A.U." );
    pHist2DNuDir->Draw( "aitoff" );
    // pHist2DNuDir->Draw( "colz" ); // debug
    cvs.SaveAs( Form( "%s/%s_nuDir.png", subDir.c_str( ), plotName.c_str( ) ) );

    // input txt
    double nTotalSI = pHistNuCosGCRateSI_tot->Integral( 1, 40 ) * EXPOSURE;
    double nTotalSD = pHistNuCosGCRateSD_tot->Integral( 1, 40 ) * EXPOSURE;

    double nPlusSI = pHistNuCosGCRateSI_tot->Integral( 21, 40 ) * EXPOSURE;
    double nPlusSD = pHistNuCosGCRateSD_tot->Integral( 21, 40 ) * EXPOSURE;

    if( atom == 11 ) { // only in case of Ag, threshold is set to be 100 keV
        //        pHistNuCosGCRateSI_E600keV->
        nTotalSI = pHistNuCosGCRateSI_E600keV->Integral( 1, 40 ) * EXPOSURE;
        nTotalSD = pHistNuCosGCRateSD_E600keV->Integral( 1, 40 ) * EXPOSURE;
        nPlusSI  = pHistNuCosGCRateSI_E600keV->Integral( 21, 40 ) * EXPOSURE;
        nPlusSD  = pHistNuCosGCRateSD_E600keV->Integral( 21, 40 ) * EXPOSURE;
    }

    double nFOMSI = nPlusSI / nTotalSI;
    double nFOMSD = nPlusSD / nTotalSD;

    double nFOMSIError = sqrt( nTotalSI*nTotalSI*nPlusSI + nTotalSI*nPlusSI*nPlusSI) / nTotalSI / nTotalSI;
    double nFOMSDError = sqrt( nTotalSD*nTotalSD*nPlusSD + nTotalSD*nPlusSD*nPlusSD) / nTotalSD / nTotalSD;

    std::cout << "total SI: " << nTotalSI << ",\tplus SI: " << nPlusSI << ", \tFOM SI: " << nFOMSI << ", pm: " << nFOMSIError << std::endl;
    std::cout << "total SD: " << nTotalSD << ",\tplus SD: " << nPlusSD << ", \tFOM SD: " << nFOMSD << ", pm: " << nFOMSDError << std::endl;

    ofs << plotName << "\t" << nTotalSI << "\t" << nPlusSI << "\t" << nFOMSI << "\t" << nFOMSIError << "\t" << nTotalSD << "\t" << nPlusSD << "\t" << nFOMSD << "\t" << nFOMSDError << "\t" << std::endl;

    return;
}
