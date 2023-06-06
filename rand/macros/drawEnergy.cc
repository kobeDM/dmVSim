#include "inc/shinclude.h"

void drawEnergy( const String& inputWIMP, const String& inputCRDMList )
{
    SetAtlasStyle( );

    TFile eneWIMPfile( inputWIMP.c_str( ) );
    TH1F* pHistEne5   = dynamic_cast< TH1F* >( eneWIMPfile.Get( "hist_SI_5_Xe" ) );
    TH1F* pHistEne10  = dynamic_cast< TH1F* >( eneWIMPfile.Get( "hist_SI_10_Xe" ) );
    TH1F* pHistEne25  = dynamic_cast< TH1F* >( eneWIMPfile.Get( "hist_SI_25_Xe" ) );
    TH1F* pHistEne50  = dynamic_cast< TH1F* >( eneWIMPfile.Get( "hist_SI_50_Xe" ) );
    TH1F* pHistEne100 = dynamic_cast< TH1F* >( eneWIMPfile.Get( "hist_SI_100_Xe" ) );
    TH1F* pHistEne200 = dynamic_cast< TH1F* >( eneWIMPfile.Get( "hist_SI_200_Xe" ) );

    // TH1F* pHistEne5   = dynamic_cast< TH1F* >( eneWIMPfile.Get( "hist_SD_5_F" ) );
    // TH1F* pHistEne10  = dynamic_cast< TH1F* >( eneWIMPfile.Get( "hist_SD_10_F" ) );
    // TH1F* pHistEne25  = dynamic_cast< TH1F* >( eneWIMPfile.Get( "hist_SD_25_F" ) );
    // TH1F* pHistEne50  = dynamic_cast< TH1F* >( eneWIMPfile.Get( "hist_SD_50_F" ) );
    // TH1F* pHistEne100 = dynamic_cast< TH1F* >( eneWIMPfile.Get( "hist_SD_100_F" ) );
    // TH1F* pHistEne200 = dynamic_cast< TH1F* >( eneWIMPfile.Get( "hist_SD_200_F" ) );
    // TH1F* pHistEne300 = dynamic_cast< TH1F* >( eneWIMPfile.Get( "hist_SD_300_F" ) );
    // TH1F* pHistEne400 = dynamic_cast< TH1F* >( eneWIMPfile.Get( "hist_SD_400_F" ) );

    if( pHistEne5   == nullptr || pHistEne10  == nullptr || pHistEne25  == nullptr || pHistEne50  == nullptr ||
        pHistEne100 == nullptr || pHistEne200 == nullptr /*|| pHistEne300 == nullptr || pHistEne400 == nullptr*/ ) return;

    pHistEne5->Scale( 1.0 / 60.0 / 60.0 / 24.0 );
    pHistEne10->Scale( 1.0 / 60.0 / 60.0 / 24.0 );
    pHistEne25->Scale( 1.0 / 60.0 / 60.0 / 24.0 );
    pHistEne50->Scale( 1.0 / 60.0 / 60.0 / 24.0 );
    pHistEne100->Scale( 1.0 / 60.0 / 60.0 / 24.0 );
    pHistEne200->Scale( 1.0 / 60.0 / 60.0 / 24.0 );
    // pHistEne300->Scale( 1.0 / 60.0 / 60.0 / 24.0 );
    // pHistEne400->Scale( 1.0 / 60.0 / 60.0 / 24.0 );

    pHistEne5->SetLineColor( 1 );
    pHistEne10->SetLineColor( 2 );
    pHistEne25->SetLineColor( 3 );
    pHistEne50->SetLineColor( 4 );
    pHistEne100->SetLineColor( 5 );
    pHistEne200->SetLineColor( 6 );
    // pHistEne300->SetLineColor( 5 );
    // pHistEne400->SetLineColor( 6 );

    TCanvas cvsWIMP( "cvsWIMP", "cvsWIMP", 800, 600 );
    cvsWIMP.SetLogy( 1 );
    // pHistEne5->GetXaxis( )->SetRangeUser( 0.0, 200.0 );
    // pHistEne5->GetXaxis( )->SetTitle( "Recoil energy [keV]" );
    // pHistEne5->Draw( );
    // pHistEne10->Draw( "same" );

    pHistEne5->GetXaxis( )->SetRangeUser( 0.0, 100.0 );
    pHistEne5->GetXaxis( )->SetTitle( "Recoil energy [keV]" );
    pHistEne5->GetYaxis( )->SetRangeUser( 0.0002, 20.0 );
    pHistEne5->GetYaxis( )->SetTitle( "Events/kg/keV/sec" );

    // pHistEne25->GetXaxis( )->SetRangeUser( 0.0, 200.0 );
    // pHistEne25->GetXaxis( )->SetTitle( "Recoil energy [keV]" );
    // pHistEne25->GetYaxis( )->SetRangeUser( 0.00000002, 0.0002 );
    // pHistEne25->GetYaxis( )->SetTitle( "Events/kg/keV/sec" );
    pHistEne5->Draw( "hist" );

    pHistEne10->Draw( "histsame" );
    pHistEne25->Draw( "histsame" );
    pHistEne50->Draw( "histsame" );
    pHistEne100->Draw( "histsame" );
    pHistEne200->Draw( "histsame" );
    // pHistEne300->Draw( "histsame" );
    // pHistEne400->Draw( "histsame" );

    TLegend* pLeg = ShTUtil::CreateLegend( 0.6, 0.6, 0.9, 0.9 );
    // pLeg->AddEntry( pHistEne5,   "#it{m}_{#chi} = 5 GeV", "f" );
    // pLeg->AddEntry( pHistEne10,  "#it{m}_{#chi} = 10 GeV", "f" );
    pLeg->AddEntry( pHistEne5,  "#it{m}_{#chi} = 5 GeV", "f" );
    pLeg->AddEntry( pHistEne10,  "#it{m}_{#chi} = 10 GeV", "f" );
    pLeg->AddEntry( pHistEne25,  "#it{m}_{#chi} = 25 GeV", "f" );
    pLeg->AddEntry( pHistEne50,  "#it{m}_{#chi} = 50 GeV", "f" );
    pLeg->AddEntry( pHistEne100, "#it{m}_{#chi} = 100 GeV", "f" );
    pLeg->AddEntry( pHistEne200, "#it{m}_{#chi} = 200 GeV", "f" );
    // pLeg->AddEntry( pHistEne300, "#it{m}_{#chi} = 300 GeV", "f" );
    // pLeg->AddEntry( pHistEne400, "#it{m}_{#chi} = 400 GeV", "f" );
    pLeg->Draw( );

    ShTUtil::CreateDrawText( 0.2, 0.85, "Normal WIMP" );
    ShTUtil::CreateDrawText( 0.2, 0.78, "Toy MC" );
    cvsWIMP.SaveAs( "testWIMP.png" );

    std::list< String > fileList;
    if( ShUtil::GetLines( inputCRDMList, &fileList ) == false ) return;

    // std::vector< TH1D* > histArray;
    std::map< double, TH1D* > histTable;
    // histArray.reserve( fileList.size( ) );

    double recEnergy = 0.0;
    double dmM = 0.0;
    double totalRateSI = 0.0, totalRateSD = 0.0;
    for( auto fileCRDM : fileList ) {
        TFile file( fileCRDM.c_str( ) );
        TTree* pTree = dynamic_cast< TTree* >( file.Get( "tree" ) );
        if( pTree == nullptr ) continue;
        
        std::cout << "File: " << fileCRDM << std::endl;
        
        pTree->SetBranchAddress( "nuRecE", &recEnergy );
        pTree->SetBranchAddress( "dmM", &dmM );
        pTree->SetBranchAddress( "totalRateSI",   &totalRateSI );
        pTree->SetBranchAddress( "totalRateSD",   &totalRateSD );

        pTree->GetEntry( 1 );
        TH1D* pHist = new TH1D( Form("histCRDM_%lfGeV", dmM), Form("histCRDM_%lfGeV", dmM), 1000, 0.0, 1000.0 );
        // TH1D* pHist = new TH1D( Form("histCRDM_%lfGeV", dmM), Form("histCRDM_%lfGeV", dmM), 200, 0.0, 200.0 );
        pHist->SetDirectory( nullptr );
        int totEvt = pTree->GetEntries( );
        for( int evtId = 0; evtId < totEvt; ++evtId ) {
            ShUtil::PrintProgressBar( evtId, totEvt );
            pTree->GetEntry( evtId );
            pHist->Fill( recEnergy * 1000000.0, totalRateSI / (double)totEvt ); // keV
            // pHist->Fill( recEnergy * 1000000.0, totalRateSD / (double)totEvt ); // keV
        }
        histTable.insert( std::make_pair( dmM, pHist ) );
    }

    TCanvas cvsCRDM( "cvsCRDM", "cvsCRDM", 800, 600 );
    cvsCRDM.SetLogy( 1 );
    bool doFirst = false;
    TLegend* pLegCRDM = ShTUtil::CreateLegend( 0.6, 0.6, 0.9, 0.9 );
    int color = 1;
    for( auto pair : histTable ) {
        double mass = pair.first;
        TH1D* pHist = pair.second;
        if( pHist == nullptr ) continue;
        pHist->SetLineColor( color );

        if( mass < 0.001 ) {
            pLegCRDM->AddEntry( pHist, Form("#it{m}_{#chi} = %0.2lf keV", mass*1000000.0), "f" );
        }
        else if( mass < 1.0 ) {
            pLegCRDM->AddEntry( pHist, Form("#it{m}_{#chi} = %0.2lf MeV", mass*1000.0), "f" );
        }
        else {
            pLegCRDM->AddEntry( pHist, Form("#it{m}_{#chi} = %0.2lf GeV", mass), "f" );
        }

        if( doFirst == false ) {
            pHist->SetLineColor( color );
            pHist->GetXaxis( )->SetTitle( "Recoil energy [keV]" );
            pHist->GetXaxis( )->SetRangeUser( 0.0, 100.0 );
            pHist->GetYaxis( )->SetTitle( "Events/kg/keV/sec" );
            pHist->GetYaxis( )->SetRangeUser( 2e-10, 1e-3 );
            pHist->Draw("hist");
            doFirst = true;
        }
        else { 
            pHist->Draw("histsame");
        }
        ++color;
        
    }
    pLegCRDM->Draw( );

    ShTUtil::CreateDrawText( 0.2, 0.85, "CRDM" );
    ShTUtil::CreateDrawText( 0.2, 0.78, "Toy MC" );

    cvsCRDM.SaveAs( "testCRDM.png" );

    return;
}
