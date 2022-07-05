#include "inc/shinclude.h"
#include "CRDMFunc.h"
#include "CRDMFunc.cc"

double getNFWProfile( double* x, double* par );
double getPITProfile( double* x, double* par );
double getEINProfile( double* x, double* par );

void drawProfile( )
{
    SetAtlasStyle( );
    
    TF1 funcNFW( "funcNFW", getNFWProfile, 0.0001, 100.0, 2, 1 );
    TF1 funcPIT( "funcPIT", getPITProfile, 0.0001, 100.0, 2, 1 );
    TF1 funcEIN( "funcEIN", getEINProfile, 0.0001, 100.0, 3, 1 );

    funcNFW.SetNpx(10000);
    funcPIT.SetNpx(10000);
    funcEIN.SetNpx(10000);

    funcNFW.SetLineColor( kRed );
    funcPIT.SetLineColor( kBlue );
    funcEIN.SetLineColor( kViolet );

    funcNFW.SetParameter( 0, DM_R_SCALE_NFW );
    funcNFW.SetParameter( 1, DM_RHO_SCALE_NFW );

    funcPIT.SetParameter( 0, DM_R_SCALE_PIT );
    funcPIT.SetParameter( 1, DM_RHO_SCALE_PIT );

    funcEIN.SetParameter( 0, DM_R_SCALE_EIN );
    funcEIN.SetParameter( 1, DM_RHO_SCALE_EIN );
    funcEIN.SetParameter( 2, DM_ALPHA_EIN );

    // scaling
    funcNFW.SetParameter( 1, DM_RHO_SCALE_NFW / funcNFW.Eval( SUN_DISTANCE ) );
    funcPIT.SetParameter( 1, DM_RHO_SCALE_PIT / funcPIT.Eval( SUN_DISTANCE ) );
    funcEIN.SetParameter( 1, DM_RHO_SCALE_EIN / funcEIN.Eval( SUN_DISTANCE ) );

    TCanvas cvs( "cvs", "cvs", 800, 600 );
    cvs.SetLogx(1);
    cvs.SetLogy(1);
    funcNFW.GetXaxis()->SetRangeUser(0.001, 50.0);
    funcNFW.GetYaxis()->SetRangeUser(0.01, 100000.0);
    funcNFW.GetXaxis()->SetTitle( "Distance from G.C. [kpc]" );
    funcNFW.GetYaxis()->SetTitle( "DM density [GeV / cm^{3}]" );


    // funcNFW.GetYaxis()->SetLimits(0.001, 100000.0);
    funcNFW.Draw();
    funcPIT.Draw("same");
    funcEIN.Draw("same");

    TLegend* pLeg = ShTUtil::CreateLegend( 0.55, 0.7, 0.98, 0.93 );
    pLeg->AddEntry( &funcNFW, "NFW", "l" );
    pLeg->AddEntry( &funcPIT, "Pseudo Iso-thermal", "l" );
    pLeg->AddEntry( &funcEIN, "Einasto", "l" );
    pLeg->Draw();

    cvs.SaveAs( "profile.png" );
    cvs.SaveAs( "profile.pdf" );
    cvs.SaveAs( "profile.eps" );

    return;
}


double getNFWProfile( double* x, double* par )
{
    double r = x[0];
    double dmRScale = par[0];
    double dmDScale = par[1];
    double relR     = r / dmRScale;
    return dmDScale / (relR*pow( 1 + relR, 2 )); // GeV/cm^3
}

double getPITProfile( double* x, double* par )
{
    double r = x[0];
    double dmRScale = par[0];
    double dmDScale = par[1];
    double relR     = r / dmRScale;
    return dmDScale / ( 1.0 + pow( relR, 2 ) ); // GeV/cm^3
}

double getEINProfile( double* x, double* par )
{
    double r = x[0];
    double dmRScale = par[0];
    double dmDScale = par[1];
    double einAlpha = par[2];
    double relR     = r / dmRScale;
    return dmDScale * exp( -2.0 * ( pow( relR, einAlpha ) - 1.0 ) / einAlpha ); // GeV/cm^3
}

