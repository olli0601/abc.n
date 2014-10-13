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

extern "C" {
/*\brief Estimate the rejection interval for one-sample dispersion equivalence numerically.
 * \param args SEXP structure of a vector of doubles with
 * 			nx	length of the observed data set
 ny	length of the simulated data set
 tl	lower tau
 tu	upper tau
 alpha	significance level of the test
 tol		numerical tolerance
 maxit	maximal number of iterations to be performed
 incit	increment to be used at each iteration
 *\return SEXP vector of length 5, holding estimates of C1, C2, the power at rho=1, the tolerance, and the number of iterations.
 */
SEXP abcScaledChiSq(SEXP args);
}

#endif /* defined(__nABC__nabc_chisqstretch__) */
