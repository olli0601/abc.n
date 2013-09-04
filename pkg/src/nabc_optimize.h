//
//  nabc_optimize.h
//  nABC
//
//  Created by anton camacho on 14/06/13.
//
//

#ifndef __nABC__nabc_optimize__
#define __nABC__nabc_optimize__

#include <iostream>
#include <cmath>
#include <Rmath.h>
#include "nabc_error_handling.h"
#include "nabc_globals.h"

double Brent_fmin(double ax, double bx, double (*f)(double, void *),
                  void *info, double tol);

#endif /* defined(__nABC__nabc_optimize__) */
