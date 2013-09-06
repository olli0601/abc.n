//
//  nabc_mutost.h
//  nABC
//
//  Created by anton camacho on 06/09/13.
//
//

#ifndef __nABC__nabc_mutost__
#define __nABC__nabc_mutost__

#include "nabc_globals.h"
#include "nabc_error_handling.h"
#include "nabc_utils.h"
#include "nabc_integrate.h"
#include "nabc_KLdiv.h"


void abcMuTOST_pow(const int &nrho, double * const rho, const double &df, const double &tau_up, const double &sT, const double &alpha, double *ans);

double abcMuTOST_pow_scalar(double x, void *arg);

void abcMuTOST_sulkl(const int &nrho, double * const rho, const double &nx, const double &sx, const double &norm, const int &give_log,double *ans);

double abcMuTOST_sulkl_scalar(double x, void *arg_void);

void abcMuTOST_taulowup_pw(	const double &mxpw, const double &df, const double &sT, const double &tau_ub, const double &alpha, const double &rho_eq, const double &tol,
                                         double &maxit, double &tau_u, double &curr_mxpw, double &error);

void abcMuTOST_taulowup_var(	const double &slkl, const double &df, const double &sT, const double &tau_ub, const double &alpha, const double &rho_eq, const double &tol,
                                          double &maxit, double &tau_u, double &curr_pwv, double &error);


void abcMUTOST_pwvar(	const double &mxpw, const double &df, const double &sT, const double &tau_ub, const double &alpha, const double &rho_eq, const double &tol,
                                   double &maxit, double &curr_pwv);


void abcMuTOST_nsim(	const double &nobs, const double &slkl, const double &mxpw, const double &sSim, const double &tau_ub, const double &alpha, const double &rho_eq, const double &tol,
                                  double &maxit, double &nsim, double &tau_u, double &curr_pwv, double &curr_pw, double &error);

void abcMuTOST_KL(const double &nx, const double &sx, const double &ny, const double &sy, const double &mx_pw, const double &alpha, const int &calibrate_tau_up, double &tau_up,const double &pow_scale, double &curr_mxpw, double &KL_div);

void abcMuTOST_KL(kl_arg *arg);



#endif /* defined(__nABC__nabc_mutost__) */
