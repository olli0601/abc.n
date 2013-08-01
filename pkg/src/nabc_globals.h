/** \file nabc_globals.h
    \brief This file contains all type definitions, macros, enumerators and global constants.
*/

#ifndef NABC_GLOBALS_H_
#define NABC_GLOBALS_H_

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <stdint.h>

typedef double flt_ty; /**< type for scientific accuracy in rate constants */
#ifdef USE_UNSIGNEDLONG
	typedef unsigned long int ulint;/**< type for population counts */
#else
typedef unsigned long long ulint;/**< type for population counts */

//typedef uint64_t ulint;/**< type for population counts */
#endif
typedef unsigned int uint;

#define NEW(x) (x *)malloc(sizeof(x))									/**< macro to allocate new memory*/
#define NEW_ARY(x,y) (x *)malloc((y)*sizeof(x))					/**< macro to allocate a contiguous memory array*/
#define NEW_ZERO_ARY(x,y) (x *)calloc(y,sizeof(x))			/**< macro to allocate a contiguous memory array that is initialized to zero*/
#define RLC_ARY(x,y,z) x = (z *)realloc((x),(y)*sizeof(z))		/**< macro to re-allocate a contiguous memory array; this preserves the content while resizing*/
#define DELETE(x) free(x)														/**< macro to free memory*/
#define MOVE(ty, fr, to, n) (ty *)memmove(to, fr, (n)*sizeof(ty))	/** < macro to move an array */	//void * memmove (void *to, const void *from, size_t size)
#define COPY_INTO_SEPARATE(ty, fr, to, n)	(ty *)memcpy(to, fr, (n)*sizeof(ty))		/** < macro to copy an array into a non-overlapping array; both arrays are assumed to be at least 'n' units long */			//void * memcpy ( void * destination, const void * source, size_t num );
#define CONST(x,y) (const_cast<x>(y))									/**< macro to make cast const*/
#define CAST(x,y) (static_cast<x>(y))									/**< macro to cast between types*/
#define MIN(x,y) ((x < y) ? x : y) 										/**< macro to return the minimum of two values */
#define MAX(x,y) ((x < y) ? y : x)										/**< macro to return the maximum of two values */
#define ABSDIFF(x,y) ((y-x>0) ? y-x : x-y) 								/**< macro to return the absolute difference of two values */
#define ABS(x) ((x>0) ? x : -x) 										/**< macro to return the absolute value of a number */


/**\brief All global variables. */
struct nabcGlobals
{
	//user-defined variables that may change
	static char BUFFER[256];/**<String buffer to perform string operations*/
    static double NABC_DBL_EPSILON;
    static double NABC_DBL_MIN;
};

typedef struct{
    double nx;
    double sx;
    double ssn; /* sx/sqrt(nx) */
    double df;
    double norm;
    int give_log;
    
    double alpha;
    double tau_up;
    double sT; /* sy/sqrt(ny) */
    
} basic_arg;

typedef struct{
    
    double (*p)(double, void *);
    double (*q)(double, void *);
    void *p_arg;
    void *q_arg;
    
} kl_integrand_arg;

typedef struct{
    
} kl_arg;




#endif /*NABC_GLOBALS_H_*/
