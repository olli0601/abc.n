/** \file nabc_error_handling.h
    \brief This file contains inline functions for error handling.
*/

#ifndef NABC_ERROR_HANDLING_H_
#define NABC_ERROR_HANDLING_H_

#include "nabc_globals.h"
#include <cstdarg>
#include <sstream> 
#include <stdexcept>

#ifdef __cplusplus
    extern "C" {
#endif

#define ERROR_ON(condition,message) if (condition){ error("%s:%s:%d:%s",__FILE__,__FUNCTION__,__LINE__,message);}//provided by Rinternal.h
#define WARNING_ON(condition,message) if (condition){ warning("%s:%s:%d:%s",__FILE__,__FUNCTION__,__LINE__,message);}//provided by Rinternal.h

#ifdef __cplusplus
    }
#endif

#endif /*NABC_ERROR_HANDLING_H_*/
