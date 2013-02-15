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
inline void postIfError(char* calling, char* format,...)
{
 	va_list args;
 	va_start (args, format);
  	vsprintf(nabcGlobals::BUFFER,format, args);
 	va_end(args);

 	char tmp[std::strlen(calling)+3+std::strlen(nabcGlobals::BUFFER)];
 	std::strcat(tmp,calling);
 	std::strcat(tmp,"\t");
 	std::strcat(tmp,nabcGlobals::BUFFER);
 	std::strcat(tmp,"\n");
 	printf("%s",tmp);

 	std::exit(EXIT_FAILURE);
 	//throw std::runtime_error(nabcGlobals::BUFFER);
}

#define POST_ERROR(X,Y) postIfError(__FILE__,X, Y)/**< macro to check for error messages */
#define FAIL_ON(condition, message) if (condition){ POST_ERROR(message,condition); }/**< macro to throw error messages */



 #ifdef __cplusplus
}
#endif

#endif /*NABC_ERROR_HANDLING_H_*/
