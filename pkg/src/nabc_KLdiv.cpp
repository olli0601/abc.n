//
//  nabc_KLdiv.cpp
//  nABC
//
//  Created by anton camacho on 06/09/13.
//
//

#include "nabc_KLdiv.h"


double abcKL_integrand(double x,void *arg_void)
{
    kl_integrand_arg *arg=(kl_integrand_arg *) arg_void;
    const double log_P=(*(arg->p))(x,arg->p_arg);
    const double log_Q=(*(arg->q))(x,arg->q_arg);
    
    return (log_P-log_Q)*exp(log_P);
    
}

void abc_generic_calibrate_tauup_for_KL(void (*KL_divergence)(void*), double (*KL_optimize)(double, void*), void *KL_arg_void, const int &max_it)
{
    double previous_KL_div, next_KL_div, tau_up_lb;
    int curr_it;
    arg_mutost *KL_arg= (arg_mutost *)  KL_arg_void;
    //do not calibrate tau_up for given max.pw
    KL_arg->calibrate_tau_up = 0;
    // current KL
    KL_arg->tau_up	= KL_arg->tau_up;
    KL_divergence( KL_arg );
    previous_KL_div = KL_arg->KL_div;
    // next KL
    KL_arg->tau_up*= 2;
    KL_divergence( KL_arg );
    next_KL_div = KL_arg->KL_div;

    curr_it = max_it;
    while((next_KL_div < previous_KL_div) && curr_it--)
    {
        previous_KL_div = next_KL_div;
        KL_arg->tau_up*= 2;
        KL_divergence( KL_arg );
        next_KL_div = KL_arg->KL_div;
    }
    ERROR_ON(!curr_it ,"could not find upper bound for tau.u within the maximum number of iterations.");
    tau_up_lb	= (curr_it == max_it) ? KL_arg->tau_up/2 : KL_arg->tau_up/4;
    //std::cout<<"tau_up_lb:"<<KL_arg->tau_up/4<<"\t tau_up_ub:"<<KL_arg->tau_up<<"\t curr_it:"<<curr_it<<std::endl;
    
    //minimize KL between KL_arg->tau_up/4 and KL_arg->tau_up
    //if curr_it==max_it then the lb is obtained by dividing ub by 2, by 4 otherwise
    KL_arg->tau_up = Brent_fmin(tau_up_lb, KL_arg->tau_up, KL_optimize, KL_arg, nabcGlobals::NABC_DBL_TOL, 0);
    
    //compute final KL_div
    KL_divergence(KL_arg);
}

void abc_generic_calibrate_yn_for_KL(void (*KL_divergence)(void*), double (*KL_optimize)(double, void*), void *KL_arg_void, const int &max_it)
{
    
    double test_KL_div,current_KL_div,ny_lb;
    int curr_it;
    arg_mutost *KL_arg= (arg_mutost *)  KL_arg_void;
   //calibrate tau_up when calibrating yn
    KL_arg->calibrate_tau_up = 1;

    // test KL
    KL_arg->ny--;
    KL_divergence( KL_arg );
    test_KL_div= KL_arg->KL_div;
    //test_KL_div = KL_divergence_switch_arg(KL_arg->ny - 1,&KL_switch_arg);
    //std::cout<<"next KL div switch:"<<next_KL_div<<"\tKL_switch_arg.KL_arg->tau_up:"<<KL_switch_arg.KL_arg->tau_up<<"\tKL_arg->tau_up:"<<KL_arg->tau_up<<std::endl;

    // current KL
    KL_arg->ny++;								// ny +1 since it has been modified by previous call to KL_divergence_switch_arg
    KL_divergence( KL_arg );
    current_KL_div = KL_arg->KL_div;
    //current_KL_div = KL_divergence_switch_arg(KL_arg->ny+1,&KL_switch_arg);
    //std::cout<<"previous KL div switch:"<<previous_KL_div<<"\tKL_switch_arg.KL_arg->tau_up:"<<KL_switch_arg.KL_arg->tau_up<<"\tKL_arg->tau_up:"<<KL_arg->tau_up<<std::endl;
    if (test_KL_div < current_KL_div) {
        ny_lb = 1;
    }
    else
    {
        curr_it 	= max_it;
        //test_KL_div = KL_divergence_switch_arg(2*KL_arg->ny,&KL_switch_arg);
        KL_arg->ny*= 2;
        KL_divergence( KL_arg );
        test_KL_div	= KL_arg->KL_div;
        while ((test_KL_div < current_KL_div) && curr_it--)
        {
            current_KL_div 	= test_KL_div;
            KL_arg->ny*= 2;
            KL_divergence( KL_arg );
            test_KL_div		= KL_arg->KL_div;
            //test_KL_div 	= KL_divergence_switch_arg(2*KL_arg->ny,&KL_switch_arg);
        }
        ERROR_ON(!curr_it ,"could not find upper bound for n.of.y within the maximum number of iterations.");
        
        //if curr_it==max_it then the lb is obtained by dividing ub by 2, by 4 otherwise
        ny_lb = (curr_it == max_it)?KL_arg->ny/2:KL_arg->ny/4;
        
    }
    //minimize KL between ny_lb and KL_arg->ny
    KL_arg->ny = Brent_fmin(ny_lb, KL_arg->ny, KL_optimize, KL_arg, nabcGlobals::NABC_DBL_TOL, 1);//
    //compute final KL_div
    KL_divergence(KL_arg);
}





