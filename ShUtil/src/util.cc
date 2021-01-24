//////////////////////////////////////////////////////////////////
//
// Utility
//
// 2014/1/29
// Satoshi Higashino
// satoshi.higashino@cern.ch
//
//////////////////////////////////////////////////////////////////
#include "inc/util.h"

//////////////////////////////////////////////////////////////////
//
// Converter
//
//////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------
// -> string data converter
template <typename T> String ShUtil::ToString( const T& val )
{
    String str = "";
    StringStream stream;

    stream << val;
    stream >> str;

    return str;
}

// ---------------------------------------------------------------
// -> string data converter
int ShUtil::StrToInt( const String& str )
{
    StringStream stream( str );
    int retVal = 0;
        
    // stream << str;
    stream >> retVal;

    return retVal;
}

// ---------------------------------------------------------------
// -> string data converter
double ShUtil::StrToDouble( const String& str )
{
    StringStream stream;
    double retVal = 0.0;
        
    stream << str;
    stream >> retVal;

    return retVal;
}

// ---------------------------------------------------------------
// Convert angle unit
// ---------------------------------------------------------------
// ---------------------------------------------------------------
// Radian -> Degree
double ShUtil::RadToDeg( double rad )
{
    return rad * 180.0 / PI;
}

// ---------------------------------------------------------------
// Degree -> Radian
double ShUtil::DegToRad( double deg )
{
    return deg * PI / 180.0;
}

// ---------------------------------------------------------------
// Convert timee unit
// ---------------------------------------------------------------
// ---------------------------------------------------------------
// hour -> minute
double ShUtil::HourToMin( double hour )
{
    return hour * 60.0;
}

// ---------------------------------------------------------------
// hour -> second
double ShUtil::HourToSec( double hour )
{
    return hour * 3600.0;
}

// ---------------------------------------------------------------
// minute -> hour
double ShUtil::MinToHour( double min )
{
    return min / 60.0;
}

// ---------------------------------------------------------------
// minute -> second
double ShUtil::MinToSec( double min )
{
    return min * 60.0;
}

// ---------------------------------------------------------------
// second -> hour
double ShUtil::SecToHour( double sec )
{
    return sec * 3600.0;
}

// ---------------------------------------------------------------
// second -> minute
double ShUtil::SecToMin( double sec )
{
    return sec * 60.0;
}

// ---------------------------------------------------------------
// Console Output
// ---------------------------------------------------------------
// ---------------------------------------------------------------
// only applied for string value.
// eg. MuUitl::Cout( "Example" );
// template <typename T> void ShUtil::Cout( const T& val )
// {
//     std::cout << val << std::endl;
// }

void ShUtil::Cout( const String& val )
{
    std::cout << val << std::endl;
}

// // ---------------------------------------------------------------
// // output error message
// template <typename T> void ShUtil::Cerr( const T& val )
// {
//     std::cout << "ShUtil Error: " << val << std::endl;
// }

void ShUtil::Cerr( const String& val )
{
    std::cout << "ShUtil Error: " << val << std::endl;
}

// ---------------------------------------------------------------
// output error message
// template <typename T> void ShUtil::Cwarn( const T& val )
// {
//     std::cout << "ShUtil Warning: " << val << std::endl;
// }

void ShUtil::Cwarn( const String& val )
{
    std::cout << "ShUtil Warning: " << val << std::endl;
}

// ---------------------------------------------------------------
// output error message
// template <typename T> void ShUtil::Cinfo( const T& val )
// {
//     std::cout << "ShUtil Information: " << val << std::endl;
// }

void ShUtil::Cinfo( const String& val )
{
    std::cout << "ShUtil Information: " << val << std::endl;
}

//////////////////////////////////////////////////////////////////
//
// Comparetor
//
//////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------
// Comparing "double" value
// at 10e-6 (the significance of float value)
// return :
//     1 for first > second
//    -1 for first < second
//     0 for first = second
int ShUtil::CompareDouble( const double& first, const double& second )
{
    if( fabs( first - second ) < 0.000001 ) {
        return 0;
    }
    else {
        if( first > second ) {
            return 1;
        }
        else {
            return -1;
        }
    }
}

// ---------------------------------------------------------------
// Comparing "double" value to Zero
// at 10e-6 (the significance of float value)
bool ShUtil::IsZero( const double& val )
{
    if( CompareDouble( val, 0.0 ) == 0 )
        return true;

    return false;
}

//////////////////////////////////////////////////////////////////
//
// File Checker
//
//////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------
// return:
//    true ...exist 
//    false...NOT exist 
bool ShUtil::ExistFile( const String& filePath )
{
    bool exist = false;
    FILE* fp = fopen( filePath.c_str( ), "r" );
    if( fp != nullptr ){
        exist = true;
        fclose( fp );
    }
    return exist;
}

// ---------------------------------------------------------------
// Check directry
// return:
//    true ...exist 
//    false...NOT exist 
// comment:
//    faster to use S_ISDIR( ) instead of opendir( )?
bool ShUtil::ExistDir( const String& dirPath )
{
    bool exist = false;
    DIR* dp = opendir( dirPath.c_str( ) );
    if( dp != nullptr ) {
        exist = true;
        closedir( dp );
    }
    return exist;
}

// ---------------------------------------------------------------
// Check directry: if the directory is not found, create it
void ShUtil::ExistCreateDir( const String& dirPath )
{
    if( ExistDir( dirPath ) == true ) return;

    // create directory
    if( mkdir( dirPath.c_str( ), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH ) == 0 ) {
        String output = "create directory : " + dirPath;
        Cinfo( output );
    }
    return;
}

// ---------------------------------------------------------------
// Check directry belonging to file path
// return:
//    true ...exist 
//    false...NOT exist 
bool ShUtil::ExistFilePathDir( const String& filePath )
{
    int index = filePath.rfind( "/", filePath.size( ) - 1 );
    String dirPath = filePath.substr( 0, index );

    return ShUtil::ExistDir( dirPath );
}

// ---------------------------------------------------------------
// Get file name from full path
// return: file name
String ShUtil::GetFileName( const String& path )
{
    size_t pos1 = path.rfind( '/' );
    if( pos1 != String::npos ) {
        return path.substr( pos1 + 1, path.size( ) - pos1 - 1 );
    }
        
    return path;
}

// ---------------------------------------------------------------
// Extract extension from file name
// return: file name (without extension)
String ShUtil::ExtractPathWithoutExt( const String& path )
{
    String::size_type pos;
    if( ( pos = path.find_last_of( "." ) ) == String::npos ) {
        return path;
    }

    return path.substr( 0, pos );
}

// ---------------------------------------------------------------
// get file path list from input directory
// argument:
//     [input]  inputDir      ... input directory (Don't include "/" at the end of directory name!!!)
//     [output] pFilePathList ... relative file path from input directory
//     [input]  recursive     ... if true, output file path recursively
//     [input]  fileKey       ... picking up only containing "fileKey" phrase in a file name
//     [input]  dirKey        ... picking up only containing "dirKey" phrase in a directory name
// return: 
//     true ... success
//     false... failure
bool ShUtil::GetFilePathList( const String&        inputDir,
                              std::list< String >* pFilePathList,
                              const bool&          recursive,
                              const String&        fileKey,
                              const String&        dirKey )
{
    if( pFilePathList == nullptr ) return false;
    // if( ExistDir( inputDir ) == false ) return false;

    bool retVal = true;
    DIR* dp = opendir( inputDir.c_str( ) );
    if( dp != nullptr ) {
        struct dirent* entry;
        struct stat st;
        String name = "";
            
        while( (entry = readdir( dp ) ) != nullptr ) {
            name = entry->d_name;
            String filePath = inputDir + "/";
            filePath += name;
                
            // skip "." and ".."
            if( name.length( ) <= 0 || name == "." || name == ".." ) continue;

            // skip stat file
            if( stat( filePath.c_str( ), &st ) == -1 ) continue;

            if( S_ISDIR( st.st_mode ) == true ) {
                if( recursive == true ) {
                    if( GetFilePathList( filePath, pFilePathList, recursive, fileKey, dirKey ) == false )
                        retVal = false;
                }
            }
            else {
                if( fileKey.length( ) > 0 && name.find( fileKey ) == String::npos ) continue;
                if( dirKey.length( ) > 0 && filePath.find( dirKey ) == String::npos ) continue;
                pFilePathList->push_back( filePath );
                // std::cout << "Add file: " << filePath << std::endl;
            }
        }

        closedir( dp );
    }
    else {
        Cerr( "error!!!" );
    }

    return retVal;
}


bool ShUtil::GetFilePathArr( const String&          inputDir,
                             std::vector< String >* pFilePathArr,
                             const bool&            recursive,
                             const String&          fileKey,
                             const String&          dirKey )
{
    if( pFilePathArr == nullptr ) return false;
    // if( ExistDir( inputDir ) == false ) return false;

    bool retVal = true;
    DIR* dp = opendir( inputDir.c_str( ) );
    if( dp != nullptr ) {
        struct dirent* entry;
        struct stat st;
        String name = "";
            
        while( (entry = readdir( dp ) ) != nullptr ) {
            name = entry->d_name;
            String filePath = inputDir + "/";
            filePath += name;
                
            // skip "." and ".."
            if( name.length( ) <= 0 || name == "." || name == ".." ) continue;

            // skip stat file
            if( stat( filePath.c_str( ), &st ) == -1 ) continue;

            if( S_ISDIR( st.st_mode ) == true ) {
                if( recursive == true ) {
                    if( GetFilePathArr( filePath, pFilePathArr, recursive, fileKey, dirKey ) == false )
                        retVal = false;
                }
            }
            else {
                if( fileKey.length( ) > 0 && name.find( fileKey ) == String::npos ) continue;
                if( dirKey.length( ) > 0 && filePath.find( dirKey ) == String::npos ) continue;
                pFilePathArr->push_back( filePath );
                // std::cout << "Add file: " << filePath << std::endl;
            }
        }

        closedir( dp );
    }
    else {
        Cerr( "error!!!" );
    }

    return retVal;
}


// ---------------------------------------------------------------
// get lines from a text file
// argument:
//     [input]  inputFile  ... input filename
//     [output] pList      ... line list
//     [input]  commSyntax ... comment syntax (default: #)
// return: 
//     true ... success
//     false... failure
bool ShUtil::GetLines( const String&        inputFile,
                       std::list< String >* pList,
                       const String&        commSyntax )
{
    if( pList == nullptr ) return false;

    std::ifstream ifs( inputFile );
    if( ifs.is_open( ) == false ) return false;

    while( !ifs.eof( ) ) {
        String line = "";
        std::getline( ifs, line );
        if( line.length( ) <= 0 || strncmp( line.c_str( ), commSyntax.c_str( ), 1 ) == 0 ) continue;
            
        pList->push_back( line );
    }

    if( pList->size( ) <= 0 ) Cwarn( "ShUtil::getLines( ) ... no lines are retrieved." );
        
    return true;
}



//////////////////////////////////////////////////////////////////
//
// Calculation
//
//////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------
// Get Sum of Squares
// 
double ShUtil::GetSumOfSqrt( const double& val1, const double& val2 )
{
    return sqrt( val1 * val1 + val2 * val2 );
}
double ShUtil::GetSumOfSqrt( const double& val1, const double& val2, const double& val3 )
{
    return sqrt( val1*val1 + val2*val2 + val3*val3 );
}

//////////////////////////////////////////////////////////////////
//
// String Operation
//
//////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------
// Replacement
// 
String ShUtil::Replace( const String& inputStr, const String& targetStr, const String& replacedStr )
{
    String retVal = "";
    size_t pos = inputStr.find( targetStr );
    if( pos != String::npos ) {
        retVal = inputStr;
        retVal.replace( pos, targetStr.length( ), replacedStr.c_str( ) );
    }
    else {
        retVal = inputStr;
    }

    return retVal;
}

// ---------------------------------------------------------------
// Split (by white spaces)
// 
bool ShUtil::Split( const String& inputStr, const char delim, std::vector< String >* pArr )
{
    if( pArr == nullptr ) return false;
    StringStream ss( inputStr );
    String element = "";
    
    while( std::getline( ss, element, delim ) ) pArr->push_back( element );
        
    return true;
}


//////////////////////////////////////////////////////////////////
//
// Progress Bar
//
//////////////////////////////////////////////////////////////////
void ShUtil::PrintProgressBar( const int& index, const int& total )
{
    if( index % 100 == 0 ) {
        String printBar = " [";
        double progress = static_cast< double >( index ) / static_cast< double >( total );
        for( int bar = 0; bar < 20; ++bar ) {
            double currentFraction = static_cast< double >( bar ) * 0.05;
            if( progress > currentFraction ) printBar += "/";
            else printBar += ".";
        }
        printBar += "] ";
        double percent = 100.0 * progress;
        StringStream percentSS;
        percentSS << std::setprecision( 2 ) << percent;
        String text = printBar + " ";
        text += percentSS.str( );
        std::cout << std::flush; 
        std::cout << text << "%\r" << std::flush; 
    }
    return;
}
    

