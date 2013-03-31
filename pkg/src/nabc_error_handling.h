/** \file nabc_error_handling.h
    \brief This file contains inline functions for error handling.
*/

#ifndef NABC_ERROR_HANDLING_H_
#define NABC_ERROR_HANDLING_H_

#include <cstdarg>
#include <sstream>
#include "nabc_globals.h"

#ifdef __cplusplus
extern "C" {
#endif

/**\brief Generate error message and exit.
 * \params calling Character array informing in which function an error occurs.
 * \params format Character array on how the error message is to be formatted.
 */
inline void postIfError(const char* format, ...)
{
	va_list args;
	va_start(args, format);
	//Rvprintf(format, args);
	vprintf(format, args);
	va_end(args);
	fflush(stdout);
 	std::exit(EXIT_FAILURE);
 	//throw std::runtime_error(nabcGlobals::BUFFER);
}

#define POST_ERROR(X,Y) postIfError(__FILE__,X, Y)/**< macro to check for error messages */
#define FAIL_ON(condition, message) if (condition){ postIfError(message, condition); }/**< macro to throw error messages */



 #ifdef __cplusplus
}
#endif

#endif /*NABC_ERROR_HANDLING_H_*/
