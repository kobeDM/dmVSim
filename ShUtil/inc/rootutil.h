//////////////////////////////////////////////////////////////////
//
// ROOT Utility
//
// 2016/4/19
// Satoshi Higashino
// satoshi.higashino@cern.ch
//
//////////////////////////////////////////////////////////////////
#ifndef SH_ROOT_UTIL_H
#define SH_ROOT_UTIL_H

// COMMON
#include "util.h"

// ROOT
#include "TString.h"
#include "TList.h"
#include "TFile.h"
#include "TKey.h"
#include "TTree.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TEfficiency.h"

#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCut.h"
#include "TArrow.h"

class ShTUtil
{
public:

    //////////////////////////////////////////////////////////////////
    //
    // Histogram
    //
    //////////////////////////////////////////////////////////////////

    // ---------------------------------------------------------------
    // Set title range
    static void SetXTitle( TH1* pHist, const String& title );
    static void SetYTitle( TH1* pHist, const String& title );
    static void SetZTitle( TH1* pHist, const String& title );

    // ---------------------------------------------------------------
    // Set range
    static void SetXRange( TH1* pHist, const Double_t& begin, const Double_t& end );
    static void SetYRange( TH1* pHist, const Double_t& begin, const Double_t& end );
    static void SetZRange( TH1* pHist, const Double_t& begin, const Double_t& end );

    // ---------------------------------------------------------------
    // Set draw style
    static void SetStyles( TH1*    pHist,
                           Color_t lineColor = 1,
                           Width_t lineWidth = 2,
                           Color_t markColor = 1,
                           Size_t  markSize  = 2,
                           Style_t markStyle = 20 );

    // ---------------------------------------------------------------
    // Set fill options
    static void SetFillOption( TH1*    pHist,
                               Color_t fillColor     = 1,
                               Style_t fillStyle     = 3001,
                               Bool_t  editLine      = false,
                               Float_t transparency  = 1.00,
                               Color_t lineColor     = 1,
                               Width_t lineWidth     = 2,
                               Style_t lineStyle     = 20 );

    // ---------------------------------------------------------------
    // Get histogram from TDirectory
    static TObject* ReadObjFromDir    ( TDirectory* pDir, const String& keyName );
    static bool     ReadObjListFromDir( TDirectory* pDir, std::list< TObject* >* pList );

    // ---------------------------------------------------------------
    // Get histogram from TDirectory
    static bool SetPointPoissonError( TGraphAsymmErrors* pGraph,
                                      const int&         index,
                                      const double&      valx,
                                      const double&      valy,
                                      const double&      weight = 1.0 );

    //////////////////////////////////////////////////////////////////
    //
    // THStack
    //
    //////////////////////////////////////////////////////////////////
    // ---------------------------------------------------------------
    // Set title range
    static void SetXTitle( THStack* pStack, const String& title );

    // ---------------------------------------------------------------
    // Set title range
    static void SetYTitle( THStack* pStack, const String& title );

    //////////////////////////////////////////////////////////////////
    //
    // Function
    //
    //////////////////////////////////////////////////////////////////

    // ---------------------------------------------------------------
    // Set draw style
    static void SetStyles( TF1*    pFunc,
                           Color_t lineColor,
                           Width_t lineWidth );

    //////////////////////////////////////////////////////////////////
    //
    // ROOT Visual Objetcs
    //
    //////////////////////////////////////////////////////////////////

    // ---------------------------------------------------------------
    // Create text in your plots
    //   note: Detail is shown in : https://root.cern.ch/doc/master/classTLatex.html
    //
    static TLatex* CreateDrawText( const double&  x,
                                   const double&  y,
                                   const String&  text,
                                   const double&  size = 0.05,
                                   const Color_t& color = 1 );

    // ---------------------------------------------------------------
    // Create legend in your plots
    //   note: Detail is shown in : https://root.cern.ch/doc/master/classTLegend.html
    //
    static TLegend* CreateLegend( const double& x1 = 0.2,
                                  const double& y1 = 0.7,
                                  const double& x2 = 0.5,
                                  const double& y2 = 0.75 );

};

#endif // SH_ROOT_UTIL_H
