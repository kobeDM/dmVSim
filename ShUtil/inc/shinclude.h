#ifndef SH_INCLUDE_H
#define SH_INCLUDE_H


#include "inc/util.h"
#include "inc/rootutil.h"
#include "inc/physutil.h"

#include "inc/AtlasLabels.h"
#include "inc/AtlasUtils.h"
#include "inc/AtlasStyle.h"
#include "inc/histmanager.h"

#ifdef __CLING__
#include "src/util.cc"
#include "src/rootutil.cc"
#include "src/physutil.cc"

#include "src/AtlasLabels.C"
#include "src/AtlasUtils.C"
#include "src/AtlasStyle.C"
#include "src/histmanager.C"
#endif

#endif // SH_INCLUDE_H
