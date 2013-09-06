//
//  nabc_chisqstretch.h
//  nABC
//
//  Created by anton camacho on 06/09/13.
//
//

#ifndef __nABC__nabc_chisqstretch__
#define __nABC__nabc_chisqstretch__

#include "nabc_globals.h"
#include "nabc_error_handling.h"
#include "nabc_utils.h"


void abcScaledChiSq_criticalregion(	const double scale, const double df, const double tl, const double tu, const double alpha, const double tol, const double incit,
                                                 double &maxit, double &c1, double &c2, double &mxpw, double &error	);

#endif /* defined(__nABC__nabc_chisqstretch__) */
