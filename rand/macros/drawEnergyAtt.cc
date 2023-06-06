#include "inc/shinclude.h"

String getAtom( const int& atom );

void drawEnergyAtt( const String& inputCRDMList )
{
    SetAtlasStyle( );
    std::list< String > fileList;
    if( ShUtil::GetLines( inputCRDMList, &fileList ) == false ) return;

    std::map< String, TH1D* > histTable;

    double recEnergy = 0.0;
    double dmM = 0.0;
    double totalRateSI = 0.0, totalRateSD = 0.0;
    int  atom = 0;
    for( auto fileCRDM : fileList ) {
        bool isAtt = false;
        if( fileCRDM.find( "att" ) != String::npos ) isAtt = true;

        TFile file( fileCRDM.c_str( ) );
        TTree* pTree = dynamic_cast< TTree* >( file.Get( "tree" ) );
        if( pTree == nullptr ) continue;
        
        std::cout << "File: " << fileCRDM << std::endl;
        
        pTree->SetBranchAddress( "nuRecE", &recEnergy );
        pTree->SetBranchAddress( "dmM", &dmM );
        pTree->SetBranchAddress( "atom", &atom );
        pTree->SetBranchAddress( "totalRateSI",   &totalRateSI );
        pTree->SetBranchAddress( "totalRateSD",   &totalRateSD );

        pTree->GetEntry( 1 );
        String histName = isAtt ? Form("histCRDM_%lfGeV_att", dmM) : Form("histCRDM_%lfGeV", dmM);


        TH1D* pHist = new TH1D( histName.c_str( ), histName.c_str( ), 1000, 0.0, 1000.0 );
        // TH1D* pHist = new TH1D( Form("histCRDM_%lfGeV", dmM), Form("histCRDM_%lfGeV", dmM), 200, 0.0, 200.0 );
        pHist->SetDirectory( nullptr );
        int totEvt = pTree->GetEntries( );
        for( int evtId = 0; evtId < totEvt; ++evtId ) {
            ShUtil::PrintProgressBar( evtId, totEvt );
            pTree->GetEntry( evtId );
            if( atom == 10 ) pHist->Fill( recEnergy * 1000000.0, totalRateSD / (double)totEvt ); // keV
            else             pHist->Fill( recEnergy * 1000000.0, totalRateSI / (double)totEvt ); // keV
            // pHist->Fill( recEnergy * 1000000.0, totalRateSD / (double)totEvt ); // keV
        }

        String legStr = "";
        if     ( dmM < 0.001 ) legStr = Form("#it{m}_{#chi} = %0.0lf keV, %s", dmM*1000000.0, getAtom( atom ).c_str( ) );
        else if( dmM < 1.0   ) legStr = Form("#it{m}_{#chi} = %0.0lf MeV, %s", dmM*1000.0, getAtom( atom ).c_str( ) );
        else                   legStr = Form("#it{m}_{#chi} = %0.0lf GeV, %s", dmM, getAtom( atom ).c_str( ) );
        if( isAtt == true ) legStr += "_att";

        histTable.insert( std::make_pair( legStr, pHist ) );
    }

    TCanvas cvsCRDM( "cvsCRDM", "cvsCRDM", 800, 600 );
    cvsCRDM.SetLogy( 1 );
    bool doFirst = false;
    TLegend* pLegCRDM = ShTUtil::CreateLegend( 0.6, 0.6, 0.9, 0.9 );
    int color = 1;
    for( auto pair : histTable ) {
        String legStr = pair.first;
        TH1D* pHist = pair.second;
        if( pHist == nullptr ) continue;
        pHist->SetLineColor( color );

        pLegCRDM->AddEntry( pHist, legStr.c_str( ), "f" );

        if( doFirst == false ) {
            pHist->SetLineColor( color );
            pHist->GetXaxis( )->SetTitle( "Recoil energy [keV]" );
            pHist->GetXaxis( )->SetRangeUser( 0.0, 100.0 );
            pHist->GetYaxis( )->SetTitle( "Events/kg/keV/sec" );
            pHist->GetYaxis( )->SetRangeUser( 2e-12, 1e-3 );
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
    ShTUtil::CreateDrawText( 0.2, 0.69, "NFW profile" );

    cvsCRDM.SaveAs( "eneCRDMAtt.png" );
    cvsCRDM.SaveAs( "eneCRDMAtt.pdf" );
    cvsCRDM.SaveAs( "eneCRDMAtt.eps" );
    // cvsCRDM.SaveAs( "eneCRDMAtom.png" );
    // cvsCRDM.SaveAs( "eneCRDMAtom.pdf" );
    // cvsCRDM.SaveAs( "eneCRDMAtom.eps" );

    return;
}

String getAtom( const int& atom )
{
    String retVal = "";
    if     ( atom == 0 ) retVal = "C";
    else if( atom == 1 ) retVal = "S";
    else if( atom == 2 ) retVal = "Br";
    else if( atom == 3 ) retVal = "I";
    else if( atom == 10 ) retVal = "F";
    else if( atom == 11 ) retVal = "Ag";
    else if( atom == 12 ) retVal = "p";
    else if( atom == 13 ) retVal = "Xe";
    
    return retVal;
}
