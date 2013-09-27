/** \file nabc_Rcall.h
 \brief Main file that provides all functions that are callable from R.
 */

#ifndef NABC_RCALL_H_
#define NABC_RCALL_H_

#include "nabc_globals.h"
#include "nabc_error_handling.h"
#include "nabc_KLdiv.h"
#include "nabc_mutost.h"
#include "nabc_chisqstretch.h"


extern "C" {
    
    SEXP abcIntersectLevelSets(SEXP m1, SEXP m2, SEXP s);
}


#endif /*NABC_RCALL_H_*/
