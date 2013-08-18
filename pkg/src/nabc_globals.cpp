#include "nabc_globals.h"
#include <limits>

char nabcGlobals::BUFFER[256]= "";
double nabcGlobals::NABC_DBL_EPSILON= std::numeric_limits<double>::epsilon();
double nabcGlobals::NABC_DBL_MIN= std::numeric_limits<double>::min();

EqTestMapValue nabcGlobals::EqTestMapEntries[] =
{
    
	EqTestMapValue( "mutost",  MUTOST_ONE_SAMPLE )
    
};

EqTestMap nabcGlobals::s_mapEqTestValue(&EqTestMapEntries[MUTOST_ONE_SAMPLE], &EqTestMapEntries[NUMBER_OF_EQ_TEST]);

//EqTestMap nabcGlobals::s_mapEqTestValue(&EqTestMapEntries[MUTOST_ONE_SAMPLE], &EqTestMapEntries[MUTOST_ONE_SAMPLE]);

