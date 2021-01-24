//////////////////////////////////////////////////////////////////
//
// Parameter definition
//
// 2013/8/27
// Satoshi Higashino
//
// Include this header at all files & use below macros.
// Please define each patameter here. 
//
//////////////////////////////////////////////////////////////////
#ifndef SH_DEF_H
#define SH_DEF_H

// ---------------------------------------------------------------
// Include standard C++ header files
// ---------------------------------------------------------------

// sys
#include <sys/types.h>
#include <sys/stat.h>

// C++ basics
#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>

#include <cmath>
#include <ctime>
#include <cstdio>

// STL
#include <set>
#include <vector>
#include <map>
#include <list>

// typedef std namespace data types
typedef std::string       String;
typedef std::stringstream StringStream;

// ---------------------------------------------------------------
// Constant parameter definition
// ---------------------------------------------------------------

// mathematical constant
const double PI                 = 3.141593;

// for debug
#define DEBUG(val) std::cout<<"Debugging : "<<val<<std::endl;

#endif // SH_DEF_H
