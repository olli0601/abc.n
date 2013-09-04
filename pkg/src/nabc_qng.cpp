/* integration/qng.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "nabc_qng.h"

static double
rescale_error (double err, const double result_abs, const double result_asc)
{
    err = fabs(err) ;
    
    if (result_asc != 0 && err != 0)
    {
        double scale = pow((200 * err / result_asc), 1.5) ;
        
        if (scale < 1)
        {
            err = result_asc * scale ;
        }
        else
        {
            err = result_asc ;
        }
    }
    if (result_abs > nabcGlobals::NABC_DBL_MIN / (50 * nabcGlobals::NABC_DBL_EPSILON))
    {
        double min_err = 50 * nabcGlobals::NABC_DBL_EPSILON * result_abs ;
        
        if (min_err > err)
        {
            err = min_err ;
        }
    }
    
    return err ;
}


int
nabc_integration_qng (double (*f)(double, void *), void *f_arg,
                      double a, double b,
                      double epsabs, double epsrel,
                      double * result, double * abserr, int * neval)
{
    
    //basic_arg *arg=(basic_arg *) f_arg;
    //std::cout<<"nabc_integration_qng\t"<<arg->ssn<<'\t'<<arg->norm<<'\t'<<arg->give_log<<'\t'<<a<<'\t'<<b<<'\t'<<epsabs<<'\t'<<epsrel<<std::endl;
    
    double fv1[5], fv2[5], fv3[5], fv4[5];
    double savfun[21];  /* array of function values which have been computed */
    double res10, res21, res43, res87;    /* 10, 21, 43 and 87 point results */
    double result_kronrod, err ;
    double resabs; /* approximation to the integral of abs(f) */
    double resasc; /* approximation to the integral of abs(f-i/(b-a)) */
    
    const double half_length =  0.5 * (b - a);
    const double abs_half_length = fabs (half_length);
    const double center = 0.5 * (b + a);
    const double f_center = (*f)(center, f_arg);
    
    int k ;
    
    if (epsabs <= 0 && (epsrel < 50 * nabcGlobals::NABC_DBL_EPSILON || epsrel < 0.5e-28))
    {
        * result = 0;
        * abserr = 0;
        * neval = 0;
        WARNING_ON(1,"tolerance cannot be acheived with given epsabs and epsrel");
    };
    
    /* Compute the integral using the 10- and 21-point formula. */
    
    res10 = 0;
    res21 = w21b[5] * f_center;
    resabs = w21b[5] * fabs (f_center);
    
    for (k = 0; k < 5; k++)
    {
        const double abscissa = half_length * x1[k];
        const double fval1 = (*f)(center + abscissa, f_arg);
        const double fval2 = (*f)(center - abscissa, f_arg);
        const double fval = fval1 + fval2;
        res10 += w10[k] * fval;
        res21 += w21a[k] * fval;
        resabs += w21a[k] * (fabs (fval1) + fabs (fval2));
        savfun[k] = fval;
        fv1[k] = fval1;
        fv2[k] = fval2;
    }
    
    for (k = 0; k < 5; k++)
    {
        const double abscissa = half_length * x2[k];
        const double fval1 = (*f)(center + abscissa, f_arg);
        const double fval2 = (*f)(center - abscissa, f_arg);
        const double fval = fval1 + fval2;
        res21 += w21b[k] * fval;
        resabs += w21b[k] * (fabs (fval1) + fabs (fval2));
        savfun[k + 5] = fval;
        fv3[k] = fval1;
        fv4[k] = fval2;
    }
    
    resabs *= abs_half_length ;
    
    {
        const double mean = 0.5 * res21;
        
        resasc = w21b[5] * fabs (f_center - mean);
        
        for (k = 0; k < 5; k++)
        {
            resasc +=
            (w21a[k] * (fabs (fv1[k] - mean) + fabs (fv2[k] - mean))
             + w21b[k] * (fabs (fv3[k] - mean) + fabs (fv4[k] - mean)));
        }
        resasc *= abs_half_length ;
    }
    
    result_kronrod = res21 * half_length;
    
    err = rescale_error ((res21 - res10) * half_length, resabs, resasc) ;
    
    /*   test for convergence. */
    
    if (err < epsabs || err < epsrel * fabs (result_kronrod))
    {
        * result = result_kronrod ;
        * abserr = err ;
        * neval = 21;
        //std::cout<<"res21\t"<<result_kronrod<<'\t'<<err<<std::endl;
        
        return 1;
    }
    
    /* compute the integral using the 43-point formula. */
    
    res43 = w43b[11] * f_center;
    
    for (k = 0; k < 10; k++)
    {
        res43 += savfun[k] * w43a[k];
    }
    
    for (k = 0; k < 11; k++)
    {
        const double abscissa = half_length * x3[k];
        const double fval = ((*f)(center + abscissa, f_arg)
                             + (*f)(center - abscissa, f_arg));
        res43 += fval * w43b[k];
        savfun[k + 10] = fval;
    }
    
    /*  test for convergence */
    
    result_kronrod = res43 * half_length;
    err = rescale_error ((res43 - res21) * half_length, resabs, resasc);
    
    if (err < epsabs || err < epsrel * fabs (result_kronrod))
    {
        * result = result_kronrod ;
        * abserr = err ;
        * neval = 43;
        //std::cout<<"res43\t"<<result_kronrod<<'\t'<<err<<std::endl;
        
        return 1;
    }
    
    /* compute the integral using the 87-point formula. */
    
    res87 = w87b[22] * f_center;
    
    for (k = 0; k < 21; k++)
    {
        res87 += savfun[k] * w87a[k];
    }
    
    for (k = 0; k < 22; k++)
    {
        const double abscissa = half_length * x4[k];
        res87 += w87b[k] * ((*f)(center + abscissa, f_arg)
                            + (*f)(center - abscissa, f_arg));
    }
    
    /*  test for convergence */
    
    result_kronrod = res87 * half_length ;
    
    err = rescale_error ((res87 - res43) * half_length, resabs, resasc);
    
    if (err < epsabs || err < epsrel * fabs (result_kronrod))
    {
        * result = result_kronrod ;
        * abserr = err ;
        * neval = 87;
        //std::cout<<"res87\t"<<result_kronrod<<'\t'<<err<<std::endl;
        
        return 1;
    }
    
    /* failed to converge */
    
    * result = result_kronrod ;
    * abserr = err ;
    * neval = 87;
    
    WARNING_ON(1,"failed to reach tolerance with highest-order rule");
    
    return 0;
    
}
