/** \file nabc_fun.h
 \brief Main file that provides all functions that are callable from R.
 */

#ifndef NABC_FUN_H_
#define NABC_FUN_H_

#include <R.h>
#include <Rinternals.h>
#include "nabc_qng.h"
#include "nabc_optimize.h"

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
    
    SEXP abcMuTOST_nsim(SEXP args);
    
    SEXP abcMuTOST_taulowup_pw(SEXP args);
    
    SEXP abcMuTOST_taulowup_var(SEXP args);
    
    SEXP abcMuTOST_pwvar(SEXP args);
    
    SEXP abcMuTOST_pow(SEXP arg_rho, SEXP arg_df, SEXP arg_tau_up, SEXP arg_sT, SEXP arg_alpha);
    
    SEXP abcMuTOST_sulkl(SEXP arg_rho, SEXP arg_nx, SEXP arg_sx, SEXP arg_norm, SEXP arg_log);
    
    SEXP abcIntersectLevelSets(SEXP m1, SEXP m2, SEXP s);
    
    SEXP abcMuTOST_sulkl_integrate_qng(SEXP arg_lower, SEXP arg_upper, SEXP arg_abs_tol, SEXP arg_rel_tol, SEXP arg_nx, SEXP arg_sx, SEXP arg_norm, SEXP arg_log);
    
    SEXP abcMuTOST_pow_integrate_qng(SEXP arg_lower, SEXP arg_upper, SEXP arg_abs_tol, SEXP arg_rel_tol, SEXP arg_df, SEXP arg_sT, SEXP arg_tau_up, SEXP arg_alpha, SEXP arg_norm, SEXP arg_log);
    
    SEXP abcMuTOST_KL(SEXP arg_nx, SEXP arg_sx, SEXP arg_ny, SEXP arg_sy, SEXP arg_mx_pw, SEXP arg_alpha, SEXP arg_calibrate_tau_up, SEXP arg_tau_up, SEXP arg_pow_scale);
    
    SEXP abcCalibrate_minimize_KL(SEXP arg_test_name, SEXP arg_calibration_name, SEXP list_KL_args, SEXP arg_max_it);
    
}


#endif /*NABC_FUN_H_*/
