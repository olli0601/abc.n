#include "nabc_globals.h"
#include <limits>

char nabcGlobals::BUFFER[256]= "";
double nabcGlobals::NABC_DBL_EPSILON= std::numeric_limits<double>::epsilon();
double nabcGlobals::NABC_DBL_MIN= std::numeric_limits<double>::min();
//char nabcGlobals::EQ_TEST_NAMES={"mutost"};
