/** \file nabc_fun.h
    \brief Main file that provides all functions that are callable from R.
*/

#ifndef NABC_FUN_H_
#define NABC_FUN_H_

#include <R.h>
#include <Rinternals.h>


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


#endif /*NABC_FUN_H_*/