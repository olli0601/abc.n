//
//  nabc_KLdiv.h
//  nABC
//
//  Created by anton camacho on 06/09/13.
//
//

#ifndef __nABC__nabc_KLdiv__
#define __nABC__nabc_KLdiv__

#include "nabc_globals.h"
#include "nabc_error_handling.h"
#include "nabc_optimize.h"
#include "nabc_integrate.h"
//#include "nabc_mutost.h"
//#include "nabc_chisqstretch.h"


double abcKL_integrand(double x,void *arg_void);

void abc_generic_calibrate_tauup_for_KL(void (*KL_divergence)(void*), double (*KL_optimize)(double, void*), void *KL_arg, const int &max_it);

void abc_generic_calibrate_yn_for_KL(void (*KL_divergence)(void*), double (*KL_optimize)(double, void*), void *KL_arg, const int &max_it);

#endif /* defined(__nABC__nabc_KLdiv__) */
