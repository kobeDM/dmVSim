//////////////////////////////////////////////////////////////////
//
// ROOT Histogram Manager
//
// 2016/12/28
// Satoshi Higashino
// satoshi.higashino@cern.ch
//
//////////////////////////////////////////////////////////////////
#ifndef SH_HIST_MANAGER_H
#define SH_HIST_MANAGER_H

// COMMON
#include "inc/shinclude.h"

enum ShTSaveOption {
    pdf = 0,
    png,
    ps,
    eps,
    jpg,
    max,
};

class ShTHistManager
{

// private:
public:
    std::map< String, TH1F* >        m_histTH1F;
    std::map< String, TH2F* >        m_histTH2F;
    std::map< String, TH3F* >        m_histTH3F;
    std::map< String, TProfile* >    m_histTProfile;
    std::map< String, TF1* >         m_funcTF1;

public:

    ShTHistManager( ){ }
    ~ShTHistManager( ){ }
    
    // Create and store TH1F histogram
    bool createTH1F( const String&              name,
                     const int&                 Nbins,
                     const double&              xmin,
                     const double&              xmax,
                     const String&              title = "" );
    bool createTH1F( const String&              name,
                     const std::vector<double>& bins,
                     const String&              title = "" );
    bool createTH1F( const String&              name,
                     const int&                 nBin,
                     const double*              pBins,
                     const String&              title = "" ); // variable bins


    // Create and store TH2F histogram
    bool createTH2F( const String&              name,
                     const int&                 NbinsX,
                     const double&              xmin,
                     const double&              xmax,
                     const int&                 NBinsY,
                     const double&              ymin,
                     const double&              ymax,
                     const String&              title = "" );
    bool createTH2F( const String&              name,
                     const std::vector<double>& xbins,
                     const std::vector<double>& ybins,
                     const String&              title = "" ); 
    bool createTH2F( const String&              name,
                     const int&                 nBinX,
                     const double*              pBinsX,
                     const int&                 nBinY,
                     const double*              pBinsY,
                     const String&              title = "" ); // variable bins
    bool createTH2F( const String&              name,
                     const int&                 nBinX,
                     const double*              pBinsX,
                     const int&                 NBinsY,
                     const double&              ymin,
                     const double&              ymax,
                     const String&              title = "" ); // variable bins
    bool createTH2F( const String&              name,
                     const int&                 NBinsX,
                     const double&              xmin,
                     const double&              xmax,
                     const int&                 nBinY,
                     const double*              pBinsY,
                     const String&              title = "" ); // variable bins
    // Create and store TH3F histogram
    bool createTH3F( const String&              name,
                     const int&                 NbinsX,
                     const double&              xmin,
                     const double&              xmax,
                     const int&                 NBinsY,
                     const double&              ymin,
                     const double&              ymax,
                     const int&                 NBinsZ,
                     const double&              zmin,
                     const double&              zmax,
                     const String&              title = "" );
    bool createTH3F( const String&              name,
                     const std::vector<double>& xbins,
                     const std::vector<double>& ybins,
                     const std::vector<double>& zbins,
                     const String&              title = "" );
    bool createTH3F( const String&              name,
                     const int&                 nBinX,
                     const double*              pBinsX,
                     const int&                 nBinY,
                     const double*              pBinsY,
                     const int&                 nBinZ,
                     const double*              pBinsZ,
                     const String&              title = "" ); // variable bins

    // Create and store TProfile histogram
    bool createTProfile( const String&              name,
                         const int&                 NbinsX,
                         const double&              xmin,
                         const double&              xmax,
                         const String&              title = "" );
    bool createTProfile( const String&              name,
                         const std::vector<double>& xbins,
                         const String&              title = "" );


    // Fill existing TH1F histogram
    inline void fillTH1F( const String& name,
                          const double& x,
                          const double& w = 1.0 ) { if( hasTH1F( name ) == true ) getTH1F( name )->Fill( x, w ); }

    // Fill existing TH2F histogram
    inline void fillTH2F( const String& name,
                          const double& x,
                          const double& y,
                          const double& w = 1.0 ) { if( hasTH2F( name ) == true ) getTH2F( name )->Fill( x, y, w ); }

    // Fill existing TH3F histogram
    inline void fillTH3F( const String& name,
                          const double& x,
                          const double& y,
                          const double& z,
                          const double& w = 1.0 ) { if( hasTH3F( name ) == true ) getTH3F( name )->Fill( x, y, z, w ); }

    // Fill existing TProfile histogram
    inline void fillTProfile( const String& name,
                              const double& x,
                              const double& y,
                              const double& w = 1.0 ) { if( hasTProfile( name ) == true ) getTProfile( name )->Fill( x, y, w ); }


    // check whether a given TH1F exist in the store
    inline bool hasTH1F( const String& name ) { return m_histTH1F.count( name ) > 0; }

    // check whether a given TH2F exist in the store
    inline bool hasTH2F( const String& name ) { return m_histTH2F.count( name ) > 0; }

    // check whether a given TH3F exists in the store
    inline bool hasTH3F( const String& name ) { return m_histTH3F.count( name ) > 0; }

    // check whether a given TProfile exist in the store
    inline bool hasTProfile( const String& name ) { return m_histTProfile.count( name ) > 0; }

    // check whether a given TH1 (inclusive) exist in the store
    inline bool hasTH1( const String& name ) {
        return ( hasTH1F( name ) || hasTH2F( name ) || hasTH3F( name ) || hasTProfile( name ) ); 
    }

    // check whether a given TF1 exist in the store
    inline bool hasTF1( const String& name ) { return m_funcTF1.count( name ) > 0; }
    
    // Retrieve TH1F histogram from internal store
    TH1F* getTH1F( const String& name );

    // Retrieve TH2F histogram from internal store
    TH2F* getTH2F( const String& name );

    // Retrieve TH3F histogram from internal store
    TH3F* getTH3F( const String& name );

    // Retrieve TProfile histogram from internal store
    TProfile* getTProfile( const String& name );

    // Retrieve TF1 function from internal store
    TF1* getTF1( const String& name );

    // Add TH1F histogram to internal store
    void addTH1F( const String& name, TH1F* pHist ) { m_histTH1F.insert( std::make_pair( name, pHist ) ); }

    // Add TH2F histogram to internal store
    void addTH2F( const String& name, TH2F* pHist ) { m_histTH2F.insert( std::make_pair( name, pHist ) ); }

    // Add TH3F histogram to internal store
    void addTH3F( const String& name, TH3F* pHist ) { m_histTH3F.insert( std::make_pair( name, pHist ) ); }

    // Add TProfile histogram to internal store
    void addTProfile            ( const String& name, TProfile* pHist ) { m_histTProfile.insert( std::make_pair( name, pHist ) ); }
    
    // Add TF1 function to internal store
    void addTF1                 ( const String& name, TF1* pFit ) { m_funcTF1.insert( std::make_pair( name, pFit ) ); }

    // Retrieve List of all histograms in internal store
    bool getListOfHistograms    ( std::vector< TH1* >* pArray );

    // Retrieve iterators
    inline std::map< String, TH1F* >::iterator     beginTH1F    ( ) { return m_histTH1F.begin( ); }
    inline std::map< String, TH2F* >::iterator     beginTH2F    ( ) { return m_histTH2F.begin( ); }
    inline std::map< String, TH3F* >::iterator     beginTH3F    ( ) { return m_histTH3F.begin( ); }
    inline std::map< String, TProfile* >::iterator beginTProfile( ) { return m_histTProfile.begin( ); }
    inline std::map< String, TF1* >::iterator      beginTF1     ( ) { return m_funcTF1.begin( ); }

    inline std::map< String, TH1F* >::iterator     endTH1F      ( ) { return m_histTH1F.end( ); }
    inline std::map< String, TH2F* >::iterator     endTH2F      ( ) { return m_histTH2F.end( ); }
    inline std::map< String, TH3F* >::iterator     endTH3F      ( ) { return m_histTH3F.end( ); }
    inline std::map< String, TProfile* >::iterator endTProfile  ( ) { return m_histTProfile.end( ); }
    inline std::map< String, TF1* >::iterator      endTF1       ( ) { return m_funcTF1.end( ); }

    // Save histograms
    bool saveTH1F        ( const String& outputDir, const String& option = "" );
    bool saveTH2F        ( const String& outputDir, const String& option = "" );
    bool saveTH3F        ( const String& outputDir, const String& option = "" );
    bool saveTProfile    ( const String& outputDir, const String& option = "" );
    bool saveAll         ( const String& outputDir, const String& option = "" );
    bool saveTF1         ( const String& outputDir, const String& option = "" );

    bool saveSameTH1F    ( const String&              outputDir,
                           const String&              keyword,
                           const String&              option = "",
                           const bool&                inheritColor = true,
                           const bool&                isLogY = false,
                           const bool&                isNorm = false,
                           const bool&                poissonErr = true,
                           const bool&                drawFitFunc = false,
                           const String&              signalName = "Signal",
                           const String&              dataName = "Data" );

    bool saveStackTH1F   ( const String&              outputDir,
                           const String&              keyword,
                           std::vector< String >*     pExcludeList,
                           const bool&                drawExcludeAsSame = false,
                           const String&              option = "",
                           const bool&                inheritColor = true );
    
    // Save as XXX.root format
    bool saveToRootFile  ( TFile* pFile );

    // Set Directory
    void setDirectoryTH1F        ( TDirectory* pDir );
    void setDirectoryTH2F        ( TDirectory* pDir );
    void setDirectoryTH3F        ( TDirectory* pDir );
    void setDirectoryTProfile    ( TDirectory* pDir );
    void setDirectoryAll         ( TDirectory* pDir );

    // Normalize
    bool normalizeTH1F( const bool& fitMax = false );
    bool normalizeTH1F( const String& name,
                        const bool&   fitMax = false );
    bool normalizeTH1F( const double& normFactor );
    bool normalizeTH1F( const String& name,
                        const double& normFactor );

    bool normalizeTH2F( const double& normFactor );
    bool normalizeTH2F( const String& name,
                        const double& normFactor );

    
    double getNormFactorTH1F( const String& name,
                              const bool&   fitMax = false );

    // Ratio Plot
    TH1F* createRatioPlot( TH1F* pNum,
                           TH1F* pDen );

    bool saveRatioPlot( const String&        outputDir,
                        const String&        keyword,
                        const String&        numName,
                        const String&        denName,
                        const String&        option = "",
                        const bool&          isLogY = false );

    bool saveRatioPlot( const String&        outputDir,
                        const String&        keyword,
                        const String&        denName,
                        const String&        option = "",
                        const bool&          isLogY = false );
    
};

#endif // SH_HIST_MANAGER_H

