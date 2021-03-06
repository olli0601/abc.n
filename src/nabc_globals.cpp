#include "nabc_globals.h"

char nabcGlobals::BUFFER[256]= "";
double nabcGlobals::NABC_DBL_EPSILON= std::numeric_limits<double>::epsilon();
double nabcGlobals::NABC_DBL_MIN= std::numeric_limits<double>::min();
double nabcGlobals::NABC_DBL_TOL= std::pow(nabcGlobals::NABC_DBL_EPSILON,0.25);

EqTestMapValue nabcGlobals::EqTestMapEntries[] =
{
    
	EqTestMapValue( "mutost",  MUTOST_ONE_SAMPLE )
    
};

EqTestMap nabcGlobals::s_mapEqTestValue(&EqTestMapEntries[MUTOST_ONE_SAMPLE], &EqTestMapEntries[NUMBER_OF_EQ_TEST]);


CalibrationMapValue nabcGlobals::CalibrationMapEntries[] =
{
    
	CalibrationMapValue( "tau_no_max_pow",  TAU_NO_MAX_POW ),
	CalibrationMapValue( "m_and_tau_yes_max_pow",  M_AND_TAU_YES_MAX_POW )
    
    
};

CalibrationMap nabcGlobals::s_mapCalibrationValue(&CalibrationMapEntries[TAU_NO_MAX_POW], &CalibrationMapEntries[NUMBER_OF_CALIBRATION]);
