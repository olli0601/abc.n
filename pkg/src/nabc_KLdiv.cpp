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



double KL_divergence_switch_arg(double x, void* KL_switch_arg){
    
    //std::cout<<"x:"<<x<<std::endl;
    
    kl_switch_arg *arg=(kl_switch_arg *) KL_switch_arg;
    
    //switch arg
    switch (arg->KL_arg_value) {
        case TAU_UP:
            arg->KL_arg->tau_up = x;
            break;
            
        case NY:
            arg->KL_arg->ny = round(x);
            break;
            
        default:
            error("%s:%s:%d: edit switch to include KlArgValue %d",__FILE__,__FUNCTION__,__LINE__,nabcGlobals::BUFFER);
            break;
    }
    
    //compute KL divergence
    arg->KL_divergence(arg->KL_arg);
    
    return arg->KL_arg->KL_div;
    
}


void abcCalibrate_tau_nomxpw_yesKL(void (*KL_divergence)(kl_arg *), kl_arg *KL_arg, const int &max_it)
{
    
    double previous_KL_div,next_KL_div;
    int curr_it;
    //initialize kl_switch_arg
    kl_switch_arg KL_switch_arg;
    //choose parameter to calibrate
    KL_switch_arg.KL_arg_value = TAU_UP;
    KL_switch_arg.KL_arg = KL_arg;
    KL_switch_arg.KL_divergence = KL_divergence;
    
    //do not calibrate tau_up when computing KL
    KL_arg->calibrate_tau_up = 0;
    
    // current KL
    previous_KL_div = KL_divergence_switch_arg(KL_arg->tau_up,&KL_switch_arg);
    //std::cout<<"previous KL div switch:"<<previous_KL_div<<"\tKL_switch_arg.KL_arg->tau_up:"<<KL_switch_arg.KL_arg->tau_up<<"\tKL_arg->tau_up:"<<KL_arg->tau_up<<std::endl;
    
    
    // next KL
    next_KL_div = KL_divergence_switch_arg(2*KL_arg->tau_up,&KL_switch_arg);
    //std::cout<<"next KL div switch:"<<next_KL_div<<"\tKL_switch_arg.KL_arg->tau_up:"<<KL_switch_arg.KL_arg->tau_up<<"\tKL_arg->tau_up:"<<KL_arg->tau_up<<std::endl;
    
    curr_it = max_it;
    
    
    while ((next_KL_div < previous_KL_div) && curr_it--) {
        previous_KL_div = next_KL_div;
        next_KL_div = KL_divergence_switch_arg(2*KL_arg->tau_up,&KL_switch_arg);
    }
    
    ERROR_ON(!curr_it ,"could not find upper bound for tau.u within the maximum number of iterations.");
    
    //std::cout<<"tau_up_lb:"<<KL_arg->tau_up/4<<"\t tau_up_ub:"<<KL_arg->tau_up<<"\t curr_it:"<<curr_it<<std::endl;
    
    
    //minimize KL between KL_arg->tau_up/4 and KL_arg->tau_up
    //if curr_it==max_it then the lb is obtained by dividing ub by 2, by 4 otherwise
    KL_arg->tau_up = Brent_fmin((curr_it == max_it)?KL_arg->tau_up/2:KL_arg->tau_up/4, KL_arg->tau_up, &KL_divergence_switch_arg, &KL_switch_arg, nabcGlobals::NABC_DBL_TOL, 0);
    
    //compute final KL_div
    KL_divergence(KL_arg);
    
}


void abcCalibrate_m_and_tau_yesmxpw_yesKL(void (*KL_divergence)(kl_arg *), kl_arg *KL_arg, const int &max_it)
{
    
    double test_KL_div,current_KL_div,ny_lb;
    int curr_it;
    //initialize kl_switch_arg
    kl_switch_arg KL_switch_arg;
    //choose parameter to calibrate
    KL_switch_arg.KL_arg_value = NY;
    KL_switch_arg.KL_arg = KL_arg;
    KL_switch_arg.KL_divergence = KL_divergence;
    
    //calibrate tau_up when computing KL
    KL_arg->calibrate_tau_up = 1;
    
    
    // test KL
    test_KL_div = KL_divergence_switch_arg(KL_arg->ny - 1,&KL_switch_arg);
    //std::cout<<"next KL div switch:"<<next_KL_div<<"\tKL_switch_arg.KL_arg->tau_up:"<<KL_switch_arg.KL_arg->tau_up<<"\tKL_arg->tau_up:"<<KL_arg->tau_up<<std::endl;
    
    
    // current KL
    current_KL_div = KL_divergence_switch_arg(KL_arg->ny+1,&KL_switch_arg);// ny +1 since it has been modified by previous call to KL_divergence_switch_arg
    //std::cout<<"previous KL div switch:"<<previous_KL_div<<"\tKL_switch_arg.KL_arg->tau_up:"<<KL_switch_arg.KL_arg->tau_up<<"\tKL_arg->tau_up:"<<KL_arg->tau_up<<std::endl;
    
    
    if (test_KL_div < current_KL_div) {
        ny_lb = 1;
    }
    else
    {
        curr_it = max_it;
        
        test_KL_div = KL_divergence_switch_arg(2*KL_arg->ny,&KL_switch_arg);
        
        while ((test_KL_div < current_KL_div) && curr_it--) {
            current_KL_div = test_KL_div;
            test_KL_div = KL_divergence_switch_arg(2*KL_arg->ny,&KL_switch_arg);
        }
        
        ERROR_ON(!curr_it ,"could not find upper bound for n.of.y within the maximum number of iterations.");
        
        //if curr_it==max_it then the lb is obtained by dividing ub by 2, by 4 otherwise
        ny_lb = (curr_it == max_it)?KL_arg->ny/2:KL_arg->ny/4;
        
    }
    
    //minimize KL between ny_lb and KL_arg->ny
    KL_arg->ny = Brent_fmin(ny_lb, KL_arg->ny, &KL_divergence_switch_arg, &KL_switch_arg, nabcGlobals::NABC_DBL_TOL, 1);//
    
    //compute final KL_div
    KL_divergence(KL_arg);
    
    
}





