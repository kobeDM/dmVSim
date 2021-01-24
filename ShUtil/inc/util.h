//////////////////////////////////////////////////////////////////
//
// Utility
//
// 2014/1/29
// Satoshi Higashino
// satoshi.higashino@cern.ch
//
//////////////////////////////////////////////////////////////////
#ifndef SH_UTIL_H
#define SH_UTIL_H

// COMMON
#include "inc/def.h"

// File operation
#include <dirent.h>

// io manipulator
#include <iomanip>

class ShUtil
{
public:

    //////////////////////////////////////////////////////////////////
    //
    // Converter
    //
    //////////////////////////////////////////////////////////////////

    // ---------------------------------------------------------------
    // -> string data converter
    template <typename T> static String ToString( const T& val );

    // ---------------------------------------------------------------
    // -> string data converter
    static int StrToInt( const String& str );

        // ---------------------------------------------------------------
    // -> string data converter
    static double StrToDouble( const String& str );

    // ---------------------------------------------------------------
    // Convert angle unit
    // ---------------------------------------------------------------
    // ---------------------------------------------------------------
    // Radian -> Degree
    static double RadToDeg( double rad );

    // ---------------------------------------------------------------
    // Degree -> Radian
    static double DegToRad( double deg );

    // ---------------------------------------------------------------
    // Convert timee unit
    // ---------------------------------------------------------------
    // ---------------------------------------------------------------
    // hour -> minute
    static double HourToMin( double hour );

    // ---------------------------------------------------------------
    // hour -> second
    static double HourToSec( double hour );

    // ---------------------------------------------------------------
    // minute -> hour
    static double MinToHour( double min );

    // ---------------------------------------------------------------
    // minute -> second
    static double MinToSec( double min );

    // ---------------------------------------------------------------
    // second -> hour
    static double SecToHour( double sec );

    // ---------------------------------------------------------------
    // second -> minute
    static double SecToMin( double sec );

    // ---------------------------------------------------------------
    // Console Output
    // ---------------------------------------------------------------
    // ---------------------------------------------------------------
    // only applied for string value.
    // eg. MuUitl::Cout( "Example" );
    // template <typename T> static void Cout( const T& val );
    static void Cout( const String& val );

    // ---------------------------------------------------------------
    // output error message
    // template <typename T> static void Cerr( const T& val );
    static void Cerr( const String& val );

    // ---------------------------------------------------------------
    // output error message
    // template <typename T> static void Cwarn( const T& val );
    static void Cwarn( const String& val );

    // ---------------------------------------------------------------
    // output error message
    // template <typename T> static void Cinfo( const T& val );
    static void Cinfo( const String& val );

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
    static int CompareDouble( const double& first, const double& second );

    // ---------------------------------------------------------------
    // Comparing "double" value to Zero
    // at 10e-6 (the significance of float value)
    static bool IsZero( const double& val );

    //////////////////////////////////////////////////////////////////
    //
    // File Checker
    //
    //////////////////////////////////////////////////////////////////
    // ---------------------------------------------------------------
    // return:
    //    true ...exist 
    //    false...NOT exist 
    static bool ExistFile( const String& filePath );

    // ---------------------------------------------------------------
    // Check directry
    // return:
    //    true ...exist 
    //    false...NOT exist 
    // comment:
    //    faster to use S_ISDIR( ) instead of opendir( )?
    static bool ExistDir( const String& dirPath );

    // ---------------------------------------------------------------
    // Check directry: if the directory is not found, create it
    static void ExistCreateDir( const String& dirPath );

    // ---------------------------------------------------------------
    // Check directry belonging to file path
    // return:
    //    true ...exist 
    //    false...NOT exist 
    static bool ExistFilePathDir( const String& filePath );

    // ---------------------------------------------------------------
    // Get file name from full path
    // return: file name
    static String GetFileName( const String& path );

    // ---------------------------------------------------------------
    // Extract extension from file name
    // return: file name (without extension)
    static String ExtractPathWithoutExt( const String& path );

    // ---------------------------------------------------------------
    // get file path list/array from input directory
    // argument:
    //     [input]  inputDir           ... input directory (Don't include "/" at the end of directory name!!!)
    //     [output] pFilePathList(Arr) ... relative file path from input directory
    //     [input]  recursive          ... if true, output file path recursively
    //     [input]  fileKey            ... picking up only containing "fileKey" phrase in a file name
    //     [input]  dirKey             ... picking up only containing "dirKey" phrase in a directory name
    // return: 
    //     true ... success
    //     false... failure
    static bool GetFilePathList( const String&        inputDir,
                                 std::list< String >* pFilePathList,
                                 const bool&          recursive    = true,
                                 const String&        fileKey      = "",
                                 const String&        dirKey       = "" );

    static bool GetFilePathArr( const String&          inputDir,
                                std::vector< String >* pFilePathArr,
                                const bool&            recursive    = true,
                                const String&          fileKey      = "",
                                const String&          dirKey       = "" );
    
    // ---------------------------------------------------------------
    // get lines from a text file
    // argument:
    //     [input]  inputFile  ... input filename
    //     [output] pList      ... line list
    //     [input]  commSyntax ... comment syntax (default: #)
    // return: 
    //     true ... success
    //     false... failure
    static bool GetLines( const String&        inputDir,
                          std::list< String >* pList,
                          const String&        commSyntax = "#" );

    //////////////////////////////////////////////////////////////////
    //
    // Calculation
    //
    //////////////////////////////////////////////////////////////////
    // ---------------------------------------------------------------
    // Get Sum of Squares
    // 
    static double GetSumOfSqrt( const double& val1, const double& val2 );
    static double GetSumOfSqrt( const double& val1, const double& val2, const double& val3 );

    //////////////////////////////////////////////////////////////////
    //
    // String Operation
    //
    //////////////////////////////////////////////////////////////////
    // ---------------------------------------------------------------
    // Replacement
    // 
    static String Replace( const String& inputStr, const String& targetStr, const String& replacedStr );

    // ---------------------------------------------------------------
    // Split (by white spaces)
    // 
    static bool Split( const String& inputStr, const char delim, std::vector< String >* pArr );

    //////////////////////////////////////////////////////////////////
    //
    // Progress Bar
    //
    //////////////////////////////////////////////////////////////////
    static void PrintProgressBar( const int& index, const int& total );
    
};

#endif // SH_UTIL_H

