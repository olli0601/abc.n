//
//  nabc_utils.h
//  nABC
//
//  Created by anton camacho on 06/09/13.
//
//

#ifndef __nABC__nabc_utils__
#define __nABC__nabc_utils__

#include "nabc_globals.h"
#include "nabc_error_handling.h"

void oprintf(const char * format, ...);

void oprintff(const char * format, ...);

void oprinta(double const * const start,const int &n,std::ostream& os);

void oseq_nout(const double &a, const double &b, const int &n, double * const ans);

void ovar(const int &n, double * const x, double * const fx, const double &mean, double &var);

int oIsZero(const int &n, double * const x);

void printBArg(basic_arg *arg);


#endif /* defined(__nABC__nabc_utils__) */
