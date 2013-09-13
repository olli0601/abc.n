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
#include "nabc_mutost.h"
#include "nabc_chisqstretch.h"


double abcKL_integrand(double x,void *arg_void);

void abcCalibrate_tau_nomxpw_yesKL(void (*KL_divergence)(kl_arg *), kl_arg *KL_arg, const int &max_it);

void abcCalibrate_m_and_tau_yesmxpw_yesKL(void (*KL_divergence)(kl_arg *), kl_arg *KL_arg, const int &max_it);


#endif /* defined(__nABC__nabc_KLdiv__) */
