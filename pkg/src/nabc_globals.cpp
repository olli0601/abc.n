#include "nabc_globals.h"
#include <cmath>

char nabcGlobals::BUFFER[256]= "";
double nabcGlobals::NABC_DBL_EPSILON= std::numeric_limits<double>::epsilon();
double nabcGlobals::NABC_DBL_MIN= std::numeric_limits<double>::min();
double nabcGlobals::NABC_DBL_TOL= std::pow(nabcGlobals::NABC_DBL_EPSILON,0.25);

EqTestMapValue nabcGlobals::EqTestMapEntries[] =
{
    
	EqTestMapValue( "mutost",  MUTOST_ONE_SAMPLE )
    
};

EqTestMap nabcGlobals::s_mapEqTestValue(&EqTestMapEntries[MUTOST_ONE_SAMPLE], &EqTestMapEntries[NUMBER_OF_EQ_TEST]);


//KlArgMapValue nabcGlobals::KlArgMapEntries[] =
//{
//    
//	KlArgMapValue( "mutost",  MUTOST_ONE_SAMPLE )
//    
//};
//
//EqTestMap nabcGlobals::s_mapEqTestValue(&EqTestMapEntries[MUTOST_ONE_SAMPLE], &EqTestMapEntries[NUMBER_OF_EQ_TEST]);
//

//EqTestMap nabcGlobals::s_mapEqTestValue(&EqTestMapEntries[MUTOST_ONE_SAMPLE], &EqTestMapEntries[MUTOST_ONE_SAMPLE]);

