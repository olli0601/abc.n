//
//  nabc_mutost.h
//  nABC
//
//
//

#ifndef __nABC__nabc_mutost__
#define __nABC__nabc_mutost__

#include "nabc_globals.h"
#include "nabc_error_handling.h"
#include "nabc_utils.h"
#include "nabc_integrate.h"
#include "nabc_KLdiv.h"

typedef struct
{
    double nx;//pass
    double sx;
    double ny;
    double sy;

    double ssn; /* sx/sqrt(nx) */
    double sT; /* sy/sqrt(ny) */
    double df;
    double norm;
    int give_log;

    double mx_pw;
    double alpha;
    int calibrate_tau_up;
    double tau_up;

    double pow_scale;
    double curr_mxpw;
    double KL_div;
} arg_mutost;

void printBArg(arg_mutost *arg);

void abcMuTOST_power(const int &nrho, double * const rho, const double &df, const double &c_up, const double &sT, double *ans);

void abcMuTOST_pow(const int &nrho, double * const rho, const double &df, const double &tau_up, const double &sT, const double &alpha, double *ans);

double abcMuTOST_power_scalar(double x, void *arg);

double abcMuTOST_pow_scalar(double x, void *arg);

void abcMuTOST_sulkl(const int &nrho, double * const rho, const double &nx, const double &sx, const double &norm, const int &give_log,double *ans);

extern "C"
{
    SEXP abcMuTOST_pwvar(SEXP args);
    SEXP abcMuTOST_power(SEXP arg_rho, SEXP arg_df, SEXP arg_c_up, SEXP arg_sT);
    SEXP abcMuTOST_pow(SEXP arg_rho, SEXP arg_df, SEXP arg_tau_up, SEXP arg_sT, SEXP arg_alpha);
    SEXP abcMuTOST_sulkl(SEXP arg_rho, SEXP arg_nx, SEXP arg_sx, SEXP arg_norm, SEXP arg_log);

    SEXP abc_mutost_integrate_sulkl(SEXP arg_lower, SEXP arg_upper, SEXP arg_abs_tol, SEXP arg_rel_tol, SEXP arg_nx, SEXP arg_sx, SEXP arg_norm, SEXP arg_log);
    SEXP abc_mutost_integrate_power(SEXP arg_lower, SEXP arg_upper, SEXP arg_abs_tol, SEXP arg_rel_tol, SEXP arg_df, SEXP arg_sT, SEXP arg_c_up, SEXP arg_norm, SEXP arg_log);
    SEXP abc_mutost_integrate_pow(SEXP arg_lower, SEXP arg_upper, SEXP arg_abs_tol, SEXP arg_rel_tol, SEXP arg_df, SEXP arg_sT, SEXP arg_tau_up, SEXP arg_alpha, SEXP arg_norm, SEXP arg_log);
    SEXP abc_mutost_get_KL(SEXP arg_nx, SEXP arg_sx, SEXP arg_ny, SEXP arg_sy, SEXP arg_mx_pw, SEXP arg_alpha, SEXP arg_calibrate_tau_up, SEXP arg_tau_up, SEXP arg_pow_scale);
    SEXP abc_mutost_calibrate_tauup_for_mxpw(SEXP args);
    SEXP abc_mutost_calibrate_powertighter(SEXP list_KL_args, SEXP arg_max_it);
    SEXP abc_mutost_calibrate_powerbroader(SEXP list_KL_args, SEXP arg_max_it);
}


#endif /* defined(__nABC__nabc_mutost__) */
