//
//  nabc_ftest.h
//  nABC
//
//
//

#ifndef __nABC__nabc_ftest__
#define __nABC__nabc_ftest__

#include "nabc_globals.h"
#include "nabc_error_handling.h"
#include "nabc_utils.h"

typedef struct
{
	double p;
    double nx;
    double t2x;
    double ny;

    double tau;
    double norm;
    int give_log;

    double mx_pw;
    double alpha;
    double tau_up;

    double pow_scale;
    double curr_mxpw;
    double pw_error;
    double KL_div;
} arg_ftest;

void nabcFTEST_printArg(arg_ftest *arg);

extern "C"
{
    SEXP abcFTEST_pow(SEXP arg_rho, SEXP arg_tau, SEXP arg_ny, SEXP arg_p, SEXP arg_alpha);
    SEXP abcFTEST_sulkl(SEXP arg_rho, SEXP arg_tx, SEXP arg_nx, SEXP arg_p, SEXP arg_norm, SEXP arg_log);
    SEXP abcFTEST_calibrate_tau_for_mxpw(SEXP args);
    SEXP abcFTEST_calibrate_KL(SEXP list_KL_args, SEXP arg_max_it);
}


#endif /* defined(__nABC__nabc_ftest__) */
