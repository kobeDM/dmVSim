#include "inc/shinclude.h"

// const double EXPOSURE_F  = 0.155 * 0.78  * 1.0 * 31536000.0; // [SF6 density (20 Torr) : kg/m3] * [F fraction in SF6] * volume [m3] * 1 year [sec]
// const double EXPOSURE_p  = 0.1   * 0.016       * 31536000.0; // [NIT density           : kg   ] * [p fraction in NIT]               * 1 year [sec]
// const double EXPOSURE_Ag = 0.1   * 0.44        * 31536000.0; // [NIT density           : kg   ] * [p fraction in NIT]               * 1 year [sec]

// assuming CYGNUS-1000 and 5 kg/year NIT
const double EXPOSURE_org  = 0.155 * 0.78  * 1.0    * 31536000.0; // [SF6 density (20 Torr) : kg/m3] * [F fraction in SF6] * volume [m3] * 1 year [sec]
const double EXPOSURE_F    = 0.155 * 0.78  * 1000.0 * 31536000.0; // [SF6 density (20 Torr) : kg/m3] * [F fraction in SF6] * 1000   [m3] * 1 year [sec]
const double EXPOSURE_p    = 5.0   * 0.016          * 31536000.0; // [NIT density           : kg   ] * [p fraction in NIT]               * 1 year [sec]
const double EXPOSURE_Ag   = 5.0   * 0.44           * 31536000.0; // [NIT density           : kg   ] * [p fraction in NIT]               * 1 year [sec]

class CRDMParam
{
public:
    
    String m_profile;
    String m_element;
    double m_mass; // [GeV]
    double m_nTotalSI;
    double m_nPlusSI;
    double m_nFOMSI;
    double m_nFOMSIError;
    double m_nTotalSD;
    double m_nPlusSD;
    double m_nFOMSD;
    double m_nFOMSDError;
};

bool getParams( const String& input, String& profile, String& element, double& mass );
bool getProfile( const String& input, String& output );
bool getElement( const String& input, String& output );
bool getMass( const String& input, double& output );

void drawSensitivity( const String& targetDir )
{
    SetAtlasStyle( );

    if( ShUtil::ExistDir( targetDir ) == false ) {
        ShUtil::Cerr( Form( "failed to find %s", targetDir.c_str( ) ) );
        abort( );
    }
    
    String filePath = targetDir + "/events.txt";
    std::ifstream ifs( filePath );
    if( ifs.is_open( ) == false ) {
        ShUtil::Cerr( Form( "failed to open %s", filePath.c_str( ) ) );
        abort( );
    }

    std::vector< CRDMParam* > paramArr;
    while( !ifs.eof( ) ) {
        String line = "";
        std::getline( ifs, line );
        if( line.length( ) <= 0 || strncmp( line.c_str( ), "#", 1 ) == 0 ) continue;
        
        StringStream ss( line );
        
        CRDMParam* pParam = new CRDMParam( );
        String plotName = "";
        ss >> plotName
           >> pParam->m_nTotalSI >> pParam->m_nPlusSI >> pParam->m_nFOMSI >> pParam->m_nFOMSIError
           >> pParam->m_nTotalSD >> pParam->m_nPlusSD >> pParam->m_nFOMSD >> pParam->m_nFOMSDError;

        if( getParams( plotName, pParam->m_profile, pParam->m_element, pParam->m_mass ) == false ) continue;
        paramArr.push_back( pParam );
    }

    /*
    // debug
    for( auto pParam : paramArr ) {
        if( pParam == nullptr ) continue;
        std::cout << pParam->m_profile << "\t"
                  << pParam->m_element << "\t"
                  << pParam->m_mass    << "\t"
                  << pParam->m_nTotalSI << "\t"
                  << pParam->m_nPlusSI << "\t"
                  << pParam->m_nTotalSD << "\t"
                  << pParam->m_nPlusSD << "\n";
    }
    */

    double mass[10]            = {};
    double massNFW[10]         = {};
    double massPIT[10]         = {};
    double massEIN[10]         = {};
    double FOM_mass_F_NFW[10]  = {};
    double FOM_mass_F_EIN[10]  = {};
    double FOM_mass_F_PIT[10]  = {};
    double FOM_mass_p_NFW[10]  = {};
    double FOM_mass_p_EIN[10]  = {};
    double FOM_mass_p_PIT[10]  = {};
    double FOM_mass_Ag_NFW[10] = {};
    double FOM_mass_Ag_EIN[10] = {};
    double FOM_mass_Ag_PIT[10] = {};

    double FOM_mass_F_NFW_err[10]  = {};
    double FOM_mass_F_EIN_err[10]  = {};
    double FOM_mass_F_PIT_err[10]  = {};
    double FOM_mass_p_NFW_err[10]  = {};
    double FOM_mass_p_EIN_err[10]  = {};
    double FOM_mass_p_PIT_err[10]  = {};
    double FOM_mass_Ag_NFW_err[10] = {};
    double FOM_mass_Ag_EIN_err[10] = {};
    double FOM_mass_Ag_PIT_err[10] = {};

    int massIdx = 0;
    // double scaling = 20.0;
    double scaling = 1.0;
    for( auto pParam : paramArr ) {
        if( pParam == nullptr ) continue;

        // calculate statistics including exposure (considering 1 year measurement using C/N-1.0 or NIT)
        double numTot  = 1.0, numPlus = 1.0;
        if( pParam->m_element == "F" ) {
            numTot = pParam->m_nTotalSD * EXPOSURE_F / EXPOSURE_org;
            numPlus = pParam->m_nPlusSD * EXPOSURE_F / EXPOSURE_org;
        }
        else if( pParam->m_element == "p" ) {
            numTot = pParam->m_nTotalSI * EXPOSURE_p / EXPOSURE_org;
            numPlus = pParam->m_nPlusSI * EXPOSURE_p / EXPOSURE_org;
        }
        else if( pParam->m_element == "Ag" ) {
            numTot = pParam->m_nTotalSI * EXPOSURE_Ag / EXPOSURE_org;
            numPlus = pParam->m_nPlusSI * EXPOSURE_Ag / EXPOSURE_org;
        }
        
        // assume 5 year data taking
        numTot  *= scaling;
        numPlus *= scaling;

        double ratio = numPlus / numTot;
        double FOM = 2.0*ratio - 1.0;
        double FOMErr = 2.0*sqrt(ratio*(1.0 - ratio)/numTot);

        if( pParam->m_mass > 0.002 ) massIdx = 0; // reset index
        mass[massIdx] = pParam->m_mass * 1000.0; // GeV -> MeV

        if( pParam->m_profile == "NFW" ) {
            if( pParam->m_element == "F" ) {
                FOM_mass_F_NFW[massIdx] = FOM;
                FOM_mass_F_NFW_err[massIdx] = FOMErr;
            }
            else if( pParam->m_element == "p" ) {
                FOM_mass_p_NFW[massIdx] = FOM;
                FOM_mass_p_NFW_err[massIdx] = FOMErr;
            }
            else if( pParam->m_element == "Ag" ) {
                FOM_mass_Ag_NFW[massIdx] = FOM;
                FOM_mass_Ag_NFW_err[massIdx] = FOMErr;
            }
            massNFW[massIdx] = mass[massIdx] * 0.8;
        }
        else if( pParam->m_profile == "IT" ) {
            if( pParam->m_element == "F" ) {
                FOM_mass_F_PIT[massIdx] = FOM;
                FOM_mass_F_PIT_err[massIdx] = FOMErr;
            }
            else if( pParam->m_element == "p" ) {
                FOM_mass_p_PIT[massIdx] = FOM;
                FOM_mass_p_PIT_err[massIdx] = FOMErr;
            }
            else if( pParam->m_element == "Ag" ) {
                FOM_mass_Ag_PIT[massIdx] = FOM;
                FOM_mass_Ag_PIT_err[massIdx] = FOMErr;
            }
            massPIT[massIdx] = mass[massIdx] * 1.0;
        }
        else if( pParam->m_profile == "EIN" ) {
            if( pParam->m_element == "F" ) {
                FOM_mass_F_EIN[massIdx] = FOM;
                FOM_mass_F_EIN_err[massIdx] = FOMErr;
            }
            else if( pParam->m_element == "p" ) {
                FOM_mass_p_EIN[massIdx] = FOM;
                FOM_mass_p_EIN_err[massIdx] = FOMErr;
            }
            else if( pParam->m_element == "Ag" ) {
                FOM_mass_Ag_EIN[massIdx] = FOM;
                FOM_mass_Ag_EIN_err[massIdx] = FOMErr;
            }
            massEIN[massIdx] = mass[massIdx] * 1.2;
        }
        
        ++massIdx;
    }

    std::cout << FOM_mass_F_NFW[2] << std::endl;
    std::cout << FOM_mass_F_NFW[4] << std::endl;

    // sensitivity plots
    double exposureF[100]  = {};
    double exposureP[100]  = {};
    double exposureAg[100] = {};
    double Z_exposure_F_NFW_mass1[100] = {};
    double Z_exposure_F_NFW_mass2[100] = {};
    double Z_exposure_F_EIN_mass1[100] = {};
    double Z_exposure_F_EIN_mass2[100] = {};
    double Z_exposure_F_PIT_mass1[100] = {};
    double Z_exposure_F_PIT_mass2[100] = {};
    double Z_exposure_p_NFW_mass1[100] = {};
    double Z_exposure_p_NFW_mass2[100] = {};
    double Z_exposure_p_EIN_mass1[100] = {};
    double Z_exposure_p_EIN_mass2[100] = {};
    double Z_exposure_p_PIT_mass1[100] = {};
    double Z_exposure_p_PIT_mass2[100] = {};
    double Z_exposure_Ag_NFW_mass1[100] = {};
    double Z_exposure_Ag_NFW_mass2[100] = {};
    double Z_exposure_Ag_EIN_mass1[100] = {};
    double Z_exposure_Ag_EIN_mass2[100] = {};
    double Z_exposure_Ag_PIT_mass1[100] = {};
    double Z_exposure_Ag_PIT_mass2[100] = {};

    double exposureF_year  = EXPOSURE_F  / 31536000.0; // kg * year
    double exposureP_year  = EXPOSURE_p  / 31536000.0; // kg * year
    double exposureAg_year = EXPOSURE_Ag / 31536000.0; // kg * year
    double den = 20.0;
    for( int i = 0; i < 100; ++i ) {
        exposureF[i]  = exposureF_year  * static_cast<double>(i) * scaling / den;
        exposureP[i]  = exposureP_year  * static_cast<double>(i) * scaling / den;
        exposureAg[i] = exposureAg_year * static_cast<double>(i) * scaling / den;
        
        // mass1 = 100 keV, mass2 = 1 keV
        Z_exposure_F_NFW_mass1[i] = (FOM_mass_F_NFW[2] / FOM_mass_F_NFW_err[2]) * sqrt( static_cast<double>(i) / den );
        Z_exposure_F_NFW_mass2[i] = (FOM_mass_F_NFW[4] / FOM_mass_F_NFW_err[4]) * sqrt( static_cast<double>(i) / den );
        Z_exposure_F_PIT_mass1[i] = (FOM_mass_F_PIT[2] / FOM_mass_F_PIT_err[2]) * sqrt( static_cast<double>(i) / den );
        Z_exposure_F_PIT_mass2[i] = (FOM_mass_F_PIT[4] / FOM_mass_F_PIT_err[4]) * sqrt( static_cast<double>(i) / den );
        Z_exposure_F_EIN_mass1[i] = (FOM_mass_F_EIN[2] / FOM_mass_F_EIN_err[2]) * sqrt( static_cast<double>(i) / den );
        Z_exposure_F_EIN_mass2[i] = (FOM_mass_F_EIN[4] / FOM_mass_F_EIN_err[4]) * sqrt( static_cast<double>(i) / den );

        Z_exposure_p_NFW_mass1[i] = (FOM_mass_p_NFW[2] / FOM_mass_p_NFW_err[2]) * sqrt( static_cast<double>(i) / den );
        Z_exposure_p_NFW_mass2[i] = (FOM_mass_p_NFW[4] / FOM_mass_p_NFW_err[4]) * sqrt( static_cast<double>(i) / den );
        Z_exposure_p_PIT_mass1[i] = (FOM_mass_p_PIT[2] / FOM_mass_p_PIT_err[2]) * sqrt( static_cast<double>(i) / den );
        Z_exposure_p_PIT_mass2[i] = (FOM_mass_p_PIT[4] / FOM_mass_p_PIT_err[4]) * sqrt( static_cast<double>(i) / den );
        Z_exposure_p_EIN_mass1[i] = (FOM_mass_p_EIN[2] / FOM_mass_p_EIN_err[2]) * sqrt( static_cast<double>(i) / den );
        Z_exposure_p_EIN_mass2[i] = (FOM_mass_p_EIN[4] / FOM_mass_p_EIN_err[4]) * sqrt( static_cast<double>(i) / den );

        Z_exposure_Ag_NFW_mass1[i] = (FOM_mass_Ag_NFW[2] / FOM_mass_Ag_NFW_err[2]) * sqrt( static_cast<double>(i) / den );
        Z_exposure_Ag_NFW_mass2[i] = (FOM_mass_Ag_NFW[4] / FOM_mass_Ag_NFW_err[4]) * sqrt( static_cast<double>(i) / den );
        Z_exposure_Ag_PIT_mass1[i] = (FOM_mass_Ag_PIT[2] / FOM_mass_Ag_PIT_err[2]) * sqrt( static_cast<double>(i) / den );
        Z_exposure_Ag_PIT_mass2[i] = (FOM_mass_Ag_PIT[4] / FOM_mass_Ag_PIT_err[4]) * sqrt( static_cast<double>(i) / den );
        Z_exposure_Ag_EIN_mass1[i] = (FOM_mass_Ag_EIN[2] / FOM_mass_Ag_EIN_err[2]) * sqrt( static_cast<double>(i) / den );
        Z_exposure_Ag_EIN_mass2[i] = (FOM_mass_Ag_EIN[4] / FOM_mass_Ag_EIN_err[4]) * sqrt( static_cast<double>(i) / den );
    }



    // fluorine
    TGraphErrors gFOM_mass_F_NFW( massIdx, massNFW, FOM_mass_F_NFW, nullptr, FOM_mass_F_NFW_err );
    TGraphErrors gFOM_mass_F_PIT( massIdx, massPIT, FOM_mass_F_PIT, nullptr, FOM_mass_F_PIT_err );
    TGraphErrors gFOM_mass_F_EIN( massIdx, massEIN, FOM_mass_F_EIN, nullptr, FOM_mass_F_EIN_err );

    gFOM_mass_F_NFW.SetMarkerStyle( 20 );
    gFOM_mass_F_PIT.SetMarkerStyle( 20 );
    gFOM_mass_F_EIN.SetMarkerStyle( 20 );
    gFOM_mass_F_NFW.SetMarkerColor( kRed );
    gFOM_mass_F_PIT.SetMarkerColor( kBlue );
    gFOM_mass_F_EIN.SetMarkerColor( kViolet );
    gFOM_mass_F_NFW.SetLineColor( kRed );
    gFOM_mass_F_PIT.SetLineColor( kBlue );
    gFOM_mass_F_EIN.SetLineColor( kViolet );

    TMultiGraph mg_mass_F;
    mg_mass_F.SetTitle(";#it{m}_{DM} [MeV];Asymmetry");
    mg_mass_F.Add(&gFOM_mass_F_NFW);
    mg_mass_F.Add(&gFOM_mass_F_PIT);
    mg_mass_F.Add(&gFOM_mass_F_EIN);

    // proton
    TGraphErrors gFOM_mass_p_NFW( massIdx, massNFW, FOM_mass_p_NFW, nullptr, FOM_mass_p_NFW_err );
    TGraphErrors gFOM_mass_p_PIT( massIdx, massPIT, FOM_mass_p_PIT, nullptr, FOM_mass_p_PIT_err );
    TGraphErrors gFOM_mass_p_EIN( massIdx, massEIN, FOM_mass_p_EIN, nullptr, FOM_mass_p_EIN_err );

    gFOM_mass_p_NFW.SetMarkerStyle( 20 );
    gFOM_mass_p_PIT.SetMarkerStyle( 20 );
    gFOM_mass_p_EIN.SetMarkerStyle( 20 );
    gFOM_mass_p_NFW.SetMarkerColor( kRed );
    gFOM_mass_p_PIT.SetMarkerColor( kBlue );
    gFOM_mass_p_EIN.SetMarkerColor( kViolet );
    gFOM_mass_p_NFW.SetLineColor( kRed );
    gFOM_mass_p_PIT.SetLineColor( kBlue );
    gFOM_mass_p_EIN.SetLineColor( kViolet );

    TMultiGraph mg_mass_p;
    mg_mass_p.SetTitle(";#it{m}_{DM} [MeV];Asymmetry");
    mg_mass_p.Add(&gFOM_mass_p_NFW);
    mg_mass_p.Add(&gFOM_mass_p_PIT);
    mg_mass_p.Add(&gFOM_mass_p_EIN);

    // Ag
    TGraphErrors gFOM_mass_Ag_NFW( massIdx, massNFW, FOM_mass_Ag_NFW, nullptr, FOM_mass_Ag_NFW_err );
    TGraphErrors gFOM_mass_Ag_PIT( massIdx, massPIT, FOM_mass_Ag_PIT, nullptr, FOM_mass_Ag_PIT_err );
    TGraphErrors gFOM_mass_Ag_EIN( massIdx, massEIN, FOM_mass_Ag_EIN, nullptr, FOM_mass_Ag_EIN_err );

    gFOM_mass_Ag_NFW.SetMarkerStyle( 20 );
    gFOM_mass_Ag_PIT.SetMarkerStyle( 20 );
    gFOM_mass_Ag_EIN.SetMarkerStyle( 20 );
    gFOM_mass_Ag_NFW.SetMarkerColor( kRed );
    gFOM_mass_Ag_PIT.SetMarkerColor( kBlue );
    gFOM_mass_Ag_EIN.SetMarkerColor( kViolet );
    gFOM_mass_Ag_NFW.SetLineColor( kRed );
    gFOM_mass_Ag_PIT.SetLineColor( kBlue );
    gFOM_mass_Ag_EIN.SetLineColor( kViolet );

    TMultiGraph mg_mass_Ag;
    mg_mass_Ag.SetTitle(";#it{m}_{DM} [MeV];Asymmetry");
    mg_mass_Ag.Add(&gFOM_mass_Ag_NFW);
    mg_mass_Ag.Add(&gFOM_mass_Ag_PIT);
    mg_mass_Ag.Add(&gFOM_mass_Ag_EIN);

    TLegend* pLeg = ShTUtil::CreateLegend( 0.2, 0.2, 0.6, 0.45 );
    pLeg->AddEntry( &gFOM_mass_F_NFW, "NFW", "lep" );
    pLeg->AddEntry( &gFOM_mass_F_PIT, "PIT", "lep" );
    pLeg->AddEntry( &gFOM_mass_F_EIN, "Einasto", "lep" );

    // draw
    TCanvas* cvsMass = new TCanvas( "cvsMass", "cvsMass", 800, 600 );
    cvsMass->SetLogx(1);
    mg_mass_F.Draw("AP");
    mg_mass_F.GetXaxis()->SetLimits( 0.00002, 50.0 );
    mg_mass_F.GetYaxis()->SetRangeUser( -5.5, 5.5 );
    mg_mass_F.Draw("AP");

    TLine l(0.00002, 0.0, 50.0, 0.0 );
    l.SetLineStyle( 2 );
    l.Draw( );
    pLeg->Draw( );
    ShTUtil::CreateDrawText( 0.2, 0.87, "F recoil, SD" );
    ShTUtil::CreateDrawText( 0.2, 0.80, "#sigma_{#chi-p} = 10^{-32} cm^{2}");
    cvsMass->SaveAs(Form("%s/asymm_mass_F.png", targetDir.c_str( )));
    cvsMass->SaveAs(Form("%s/asymm_mass_F.pdf", targetDir.c_str( )));
    cvsMass->SaveAs(Form("%s/asymm_mass_F.eps", targetDir.c_str( )));

    mg_mass_p.Draw("AP");
    mg_mass_p.GetXaxis()->SetLimits( 0.00002, 50.0 );
    mg_mass_p.GetYaxis()->SetRangeUser( -5.5, 5.5 );
    mg_mass_p.Draw("AP");
    l.Draw( );
    pLeg->Draw( );
    ShTUtil::CreateDrawText( 0.2, 0.87, "p recoil, SI" );
    ShTUtil::CreateDrawText( 0.2, 0.80, "#sigma_{#chi-p} = 10^{-32} cm^{2}");
    cvsMass->SaveAs(Form("%s/asymm_mass_p.png", targetDir.c_str( )));
    cvsMass->SaveAs(Form("%s/asymm_mass_p.pdf", targetDir.c_str( )));
    cvsMass->SaveAs(Form("%s/asymm_mass_p.eps", targetDir.c_str( )));

    mg_mass_Ag.Draw("AP");
    mg_mass_Ag.GetXaxis()->SetLimits( 0.00002, 50.0 );
    mg_mass_Ag.GetYaxis()->SetRangeUser( -5.5, 5.5 );
    mg_mass_Ag.Draw("AP");
    l.Draw( );
    pLeg->Draw( );
    ShTUtil::CreateDrawText( 0.2, 0.87, "Ag recoil, SI" );
    ShTUtil::CreateDrawText( 0.2, 0.80, "#sigma_{#chi-p} = 10^{-32} cm^{2}");
    cvsMass->SaveAs(Form("%s/asymm_mass_Ag.png", targetDir.c_str( )));
    cvsMass->SaveAs(Form("%s/asymm_mass_Ag.pdf", targetDir.c_str( )));
    cvsMass->SaveAs(Form("%s/asymm_mass_Ag.eps", targetDir.c_str( )));

    delete cvsMass;

    TGraph gZ_exposure_F_NFW_mass1( 100, exposureF, Z_exposure_F_NFW_mass1 );
    TGraph gZ_exposure_F_NFW_mass2( 100, exposureF, Z_exposure_F_NFW_mass2 );
    TGraph gZ_exposure_F_PIT_mass1( 100, exposureF, Z_exposure_F_PIT_mass1 );
    TGraph gZ_exposure_F_PIT_mass2( 100, exposureF, Z_exposure_F_PIT_mass2 );
    TGraph gZ_exposure_F_EIN_mass1( 100, exposureF, Z_exposure_F_EIN_mass1 );
    TGraph gZ_exposure_F_EIN_mass2( 100, exposureF, Z_exposure_F_EIN_mass2 );

    TGraph gZ_exposure_p_NFW_mass1( 100, exposureP, Z_exposure_p_NFW_mass1 );
    TGraph gZ_exposure_p_NFW_mass2( 100, exposureP, Z_exposure_p_NFW_mass2 );
    TGraph gZ_exposure_p_PIT_mass1( 100, exposureP, Z_exposure_p_PIT_mass1 );
    TGraph gZ_exposure_p_PIT_mass2( 100, exposureP, Z_exposure_p_PIT_mass2 );
    TGraph gZ_exposure_p_EIN_mass1( 100, exposureP, Z_exposure_p_EIN_mass1 );
    TGraph gZ_exposure_p_EIN_mass2( 100, exposureP, Z_exposure_p_EIN_mass2 );

    TGraph gZ_exposure_Ag_NFW_mass1( 100, exposureAg, Z_exposure_Ag_NFW_mass1 );
    TGraph gZ_exposure_Ag_NFW_mass2( 100, exposureAg, Z_exposure_Ag_NFW_mass2 );
    TGraph gZ_exposure_Ag_PIT_mass1( 100, exposureAg, Z_exposure_Ag_PIT_mass1 );
    TGraph gZ_exposure_Ag_PIT_mass2( 100, exposureAg, Z_exposure_Ag_PIT_mass2 );
    TGraph gZ_exposure_Ag_EIN_mass1( 100, exposureAg, Z_exposure_Ag_EIN_mass1 );
    TGraph gZ_exposure_Ag_EIN_mass2( 100, exposureAg, Z_exposure_Ag_EIN_mass2 );

    gZ_exposure_F_NFW_mass1.SetLineColor( kRed );
    gZ_exposure_F_NFW_mass2.SetLineColor( kRed );
    gZ_exposure_F_PIT_mass1.SetLineColor( kBlue );
    gZ_exposure_F_PIT_mass2.SetLineColor( kBlue );
    gZ_exposure_F_EIN_mass1.SetLineColor( kViolet );
    gZ_exposure_F_EIN_mass2.SetLineColor( kViolet );

    gZ_exposure_p_NFW_mass1.SetLineColor( kRed );   
    gZ_exposure_p_NFW_mass2.SetLineColor( kRed );   
    gZ_exposure_p_PIT_mass1.SetLineColor( kBlue );  
    gZ_exposure_p_PIT_mass2.SetLineColor( kBlue );  
    gZ_exposure_p_EIN_mass1.SetLineColor( kViolet );
    gZ_exposure_p_EIN_mass2.SetLineColor( kViolet );

    gZ_exposure_Ag_NFW_mass1.SetLineColor( kRed );   
    gZ_exposure_Ag_NFW_mass2.SetLineColor( kRed );   
    gZ_exposure_Ag_PIT_mass1.SetLineColor( kBlue );  
    gZ_exposure_Ag_PIT_mass2.SetLineColor( kBlue );  
    gZ_exposure_Ag_EIN_mass1.SetLineColor( kViolet );
    gZ_exposure_Ag_EIN_mass2.SetLineColor( kViolet );

    gZ_exposure_F_NFW_mass1.SetLineWidth(2);
    gZ_exposure_F_NFW_mass2.SetLineWidth(2);
    gZ_exposure_F_PIT_mass1.SetLineWidth(2);
    gZ_exposure_F_PIT_mass2.SetLineWidth(2);
    gZ_exposure_F_EIN_mass1.SetLineWidth(2);
    gZ_exposure_F_EIN_mass2.SetLineWidth(2);
    gZ_exposure_p_NFW_mass1.SetLineWidth(2);
    gZ_exposure_p_NFW_mass2.SetLineWidth(2);
    gZ_exposure_p_PIT_mass1.SetLineWidth(2);
    gZ_exposure_p_PIT_mass2.SetLineWidth(2);
    gZ_exposure_p_EIN_mass1.SetLineWidth(2);
    gZ_exposure_p_EIN_mass2.SetLineWidth(2);
    gZ_exposure_Ag_NFW_mass1.SetLineWidth(2);
    gZ_exposure_Ag_NFW_mass2.SetLineWidth(2);
    gZ_exposure_Ag_PIT_mass1.SetLineWidth(2);
    gZ_exposure_Ag_PIT_mass2.SetLineWidth(2);
    gZ_exposure_Ag_EIN_mass1.SetLineWidth(2);
    gZ_exposure_Ag_EIN_mass2.SetLineWidth(2);

    gZ_exposure_F_NFW_mass2.SetLineStyle( 3 );
    gZ_exposure_F_PIT_mass2.SetLineStyle( 3 );
    gZ_exposure_F_EIN_mass2.SetLineStyle( 3 );
    gZ_exposure_p_NFW_mass2.SetLineStyle( 3 );
    gZ_exposure_p_PIT_mass2.SetLineStyle( 3 );
    gZ_exposure_p_EIN_mass2.SetLineStyle( 3 );
    gZ_exposure_Ag_NFW_mass2.SetLineStyle( 3 );
    gZ_exposure_Ag_PIT_mass2.SetLineStyle( 3 );
    gZ_exposure_Ag_EIN_mass2.SetLineStyle( 3 );

    TMultiGraph mg_Z_F;
    mg_Z_F.SetTitle( ";Exposure [kg year];Sensitivity [#sigma]" );
    mg_Z_F.Add( &gZ_exposure_F_NFW_mass1 );
    mg_Z_F.Add( &gZ_exposure_F_NFW_mass2 );
    mg_Z_F.Add( &gZ_exposure_F_PIT_mass1 );
    mg_Z_F.Add( &gZ_exposure_F_PIT_mass2 );
    mg_Z_F.Add( &gZ_exposure_F_EIN_mass1 );
    mg_Z_F.Add( &gZ_exposure_F_EIN_mass2 );

    TMultiGraph mg_Z_p;
    mg_Z_p.SetTitle( ";Exposure [kg year];Sensitivity [#sigma]" );
    mg_Z_p.Add( &gZ_exposure_p_NFW_mass1 );
    mg_Z_p.Add( &gZ_exposure_p_NFW_mass2 );
    mg_Z_p.Add( &gZ_exposure_p_PIT_mass1 );
    mg_Z_p.Add( &gZ_exposure_p_PIT_mass2 );
    mg_Z_p.Add( &gZ_exposure_p_EIN_mass1 );
    mg_Z_p.Add( &gZ_exposure_p_EIN_mass2 );

    TMultiGraph mg_Z_Ag;
    mg_Z_Ag.SetTitle( ";Exposure [kg year];Sensitivity [#sigma]" );
    mg_Z_Ag.Add( &gZ_exposure_Ag_NFW_mass1 );
    mg_Z_Ag.Add( &gZ_exposure_Ag_NFW_mass2 );
    mg_Z_Ag.Add( &gZ_exposure_Ag_PIT_mass1 );
    mg_Z_Ag.Add( &gZ_exposure_Ag_PIT_mass2 );
    mg_Z_Ag.Add( &gZ_exposure_Ag_EIN_mass1 );
    mg_Z_Ag.Add( &gZ_exposure_Ag_EIN_mass2 );

    TLegend* pLegZ = ShTUtil::CreateLegend( 0.19, 0.55, 0.65, 0.84 );
    pLegZ->AddEntry( &gZ_exposure_F_NFW_mass1, "NFW   #it{m}_{DM} = 100 keV", "l" );
    pLegZ->AddEntry( &gZ_exposure_F_NFW_mass2, "NFW   #it{m}_{DM} = 1 keV", "l" );
    pLegZ->AddEntry( &gZ_exposure_F_PIT_mass1, "PIT   #it{m}_{DM} = 100 keV", "l" );
    pLegZ->AddEntry( &gZ_exposure_F_PIT_mass2, "PIT   #it{m}_{DM} = 1 keV", "l" );
    pLegZ->AddEntry( &gZ_exposure_F_EIN_mass1, "Einasto   #it{m}_{DM} = 100 keV", "l" );
    pLegZ->AddEntry( &gZ_exposure_F_EIN_mass2, "Einasto   #it{m}_{DM} = 1 keV", "l" );
    
    TCanvas* cvsZ = new TCanvas( "cvsZ", "cvsZ", 800, 600 );
    mg_Z_F.Draw( "AC" );
    mg_Z_F.GetYaxis()->SetRangeUser( 0.0, Z_exposure_F_EIN_mass2[99] * 1.5);
    mg_Z_F.Draw( "AC" );
    pLegZ->Draw( );
    ShTUtil::CreateDrawText(0.2, 0.87, "F recoil, SD" );
    ShTUtil::CreateDrawText(0.7, 0.87, "#sigma_{#chi-p} = 10^{-32} cm^{2}");
    cvsZ->SaveAs(Form("%s/z_F.png", targetDir.c_str( )));
    cvsZ->SaveAs(Form("%s/z_F.pdf", targetDir.c_str( )));
    cvsZ->SaveAs(Form("%s/z_F.eps", targetDir.c_str( )));

    mg_Z_p.Draw( "AC" );
    mg_Z_p.GetYaxis()->SetRangeUser( 0.0, Z_exposure_p_EIN_mass2[99] * 1.5);
    mg_Z_p.Draw( "AC" );
    pLegZ->Draw( );
    ShTUtil::CreateDrawText(0.2, 0.87, "p recoil, SI");
    ShTUtil::CreateDrawText(0.7, 0.87, "#sigma_{#chi-p} = 10^{-32} cm^{2}");
    cvsZ->SaveAs(Form("%s/z_p.png", targetDir.c_str( )));
    cvsZ->SaveAs(Form("%s/z_p.pdf", targetDir.c_str( )));
    cvsZ->SaveAs(Form("%s/z_p.eps", targetDir.c_str( )));

    mg_Z_Ag.Draw( "AC" );
    mg_Z_Ag.GetYaxis()->SetRangeUser( 0.0, Z_exposure_Ag_EIN_mass2[99] * 1.5);
    mg_Z_Ag.Draw( "AC" );
    pLegZ->Draw( );
    ShTUtil::CreateDrawText(0.2, 0.87, "Ag recoil, SI");
    ShTUtil::CreateDrawText(0.7, 0.87, "#sigma_{#chi-p} = 10^{-32} cm^{2}");
    cvsZ->SaveAs(Form("%s/z_Ag.png", targetDir.c_str( )));
    cvsZ->SaveAs(Form("%s/z_Ag.pdf", targetDir.c_str( )));
    cvsZ->SaveAs(Form("%s/z_Ag.eps", targetDir.c_str( )));

    delete cvsZ;

    return;
}


bool getParams( const String& input, String& profile, String& element, double& mass )
{
    bool retVal = false;

    if( input.length( ) <= 0 ) return retVal;

    std::vector< String > arr;
    if( ShUtil::Split( input, '_', &arr ) == false ) return retVal;
    
    bool profileOK = false;
    bool elementOK = false;
    bool massOK    = false;
    for( auto str : arr ) {
        if( str == "recoil" ) continue;
        
        // DEBUG(str);
        if( profileOK == false && getProfile( str, profile ) == true ) profileOK = true;
        if( elementOK == false && getElement( str, element ) == true ) elementOK = true;
        if( massOK    == false && getMass   ( str, mass    ) == true ) massOK    = true; 
    }

    if( profileOK == true && elementOK == true && massOK == true ) retVal = true;

    return retVal;
}

bool getProfile( const String& input, String& output )
{
    bool retVal = true;
    if     ( input == "NFW" ) output = "NFW";
    else if( input == "IT"  ) output = "IT";
    else if( input == "EIN" ) output = "EIN";
    else retVal = false;
    
    return retVal;
}

bool getElement( const String& input, String& output )
{
    bool retVal = true;
    if     ( input == "F"  ) output = "F";
    else if( input == "Ag" ) output = "Ag";
    else if( input == "p"  ) output = "p";
    else retVal = false;
    
    return retVal;
}


bool getMass( const String& input, double& output )
{
    bool retVal = false;
    if( input.find( "los" ) == String::npos ) return retVal;

    size_t pos = input.find( "dm" );
    String massStr = input.substr( pos+2, input.length( ) );

    retVal = true;
    if     ( massStr == "10"       ) output = 10.0;
    else if( massStr == "1"        ) output = 1.0;
    else if( massStr == "01"       ) output = 0.1;
    else if( massStr == "001"      ) output = 0.01;
    else if( massStr == "0001"     ) output = 0.001;
    else if( massStr == "00001"    ) output = 0.0001;
    else if( massStr == "000001"   ) output = 0.00001;
    else if( massStr == "0000001"  ) output = 0.000001;
    else if( massStr == "00000001" ) output = 0.0000001;
    else                             retVal = false;
    
    return retVal;
}
