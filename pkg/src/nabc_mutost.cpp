//
//  nabc_mutost.cpp
//  nABC
//
//
//

#include "nabc_mutost.h"

void abcMuTOST_pow(const int &nrho, double * const rho, const double &df, const double &tau_up, const double &sT, const double &alpha, double *ans)
{
	int n= nrho;
	const int LOG= 0, LOWERTAIL=1;
	double *xans= ans, *xrho=rho;
	const double QU= qt(alpha, df, LOWERTAIL, LOG), TIS=tau_up/sT;
    
	for(; n--; xrho++, xans++)
	{
		*xans= pnt( TIS + QU, df, *xrho/sT, LOWERTAIL, LOG) - pnt( -TIS - QU, df, *xrho/sT, LOWERTAIL, LOG);
		*xans= *xans<0 ? 0 : *xans;
	}
}

double abcMuTOST_pow_scalar(double x, void *arg)
{
    basic_arg *a=(basic_arg *) arg;
	double ans;
	const double QU_TIS= qt(a->alpha, a->df, 1, 0)+a->tau_up/a->sT;
    
    ans= pnt( QU_TIS, a->df, x/a->sT, 1, 0) - pnt( -QU_TIS, a->df, x/a->sT, 1, 0);
    
    /* Can be 0 at some point! */
    /* Usually this happens due to numerical inacurracy in the tail. */
    /* To avoid infinity we assume that ans=nabcGlobals::NABC_DBL_MIN at these points. */
    
    ans= ans<=0 ? nabcGlobals::NABC_DBL_MIN : ans/a->norm;
    
    return a->give_log ? log(ans) : ans;
    
}


void abcMuTOST_sulkl(const int &nrho, double * const rho, const double &nx, const double &sx, const double &norm, const int &give_log,double *ans)
{
	int n= nrho;
	const double ssn= sx/sqrt(nx),df=nx-1;
	double *xans= ans, *xrho=rho;
    
	for(; n--; xrho++, xans++)
	{
		*xans= dt( *xrho/ssn, df, give_log);
        *xans= give_log ? *xans-log(ssn*norm) : *xans/ssn/norm;
	}
}

double abcMuTOST_sulkl_scalar(double x, void *arg_void)
{
    basic_arg *arg=(basic_arg *) arg_void;
    //std::cout<<"abcMuTOST_sulkl_scalar input:\nssn\t"<<arg->ssn<<"\ndf\t"<<arg->df<<"\ngive_log\t"<<arg->give_log<<"\nnorm\t"<<arg->norm<<"\nx\t"<<x<<std::endl;
    
	double ans;
    
    ans= dt( x/arg->ssn, arg->df, arg->give_log);
    ans= arg->give_log ? ans-log(arg->ssn*arg->norm) : ans/arg->ssn/arg->norm;
	
    return ans;
}

void abcMuTOST_taulowup_pw(	const double &mxpw, const double &df, const double &sT, const double &tau_ub, const double &alpha, const double &rho_eq, const double &tol,
                                         double &maxit, double &tau_u, double &curr_mxpw, double &error)
{
	const int NRHO=1;
	int n= CAST(int,maxit);
	const double DIGITS= std::ldexp(1,33);
	double tau_ubd= tau_ub/2, tau_lbd=0, *rho= NULL;
    
	rho= NEW_ARY(double,NRHO);
	*rho= rho_eq;
	for(curr_mxpw=0; 	n-- && curr_mxpw<mxpw; 		)
	{
		tau_ubd*=2;
		abcMuTOST_pow(NRHO, rho, df, tau_ubd, sT, alpha, &curr_mxpw);
		//std::cout<<"H1 "<<curr_mxpw<<'\t'<<mxpw<<'\t'<<tau_ubd<<'\t'<<n<<std::endl;
	}
	if(n<0)	std::cout<<"nabcMuTOST_taulowup_pw: "<<tau_ub<<'\t'<<tau_ubd<<'\t'<<curr_mxpw<<'\t'<<mxpw<<'\t'<<df<<'\t'<<sT<<'\t'<<maxit<<std::endl;
	FAIL_ON(n<0,"\nabcMuTOST_taulowup_pw: error at 1a %c");
	for(	error=1, n= CAST(int,maxit);
        n-- && (ABS(error)>tol) && std::floor(tau_ubd*DIGITS)!=std::floor(tau_lbd*DIGITS);
        )
	{
		tau_u= (tau_lbd+tau_ubd)/2;
		abcMuTOST_pow(NRHO, rho, df, tau_u, sT, alpha, &curr_mxpw);
		error= curr_mxpw-mxpw;
		if(error<0)
			tau_lbd= tau_u;
		else
			tau_ubd= tau_u;
		//std::cout<<"H2 "<<curr_mxpw<<'\t'<<mxpw<<'\t'<<tau_ubd<<'\t'<<tau_u<<'\t'<<n<<std::endl;
	}
	if(n<0)
		oprintff("\nabcMuTOST_tau_taulowup: reached max it %i current pw %g requ pw %g",n,curr_mxpw,mxpw); //could not find Rvprintf ??
	maxit= n+1;
	DELETE(rho);
}

void abcMuTOST_taulowup_var(	const double &slkl, const double &df, const double &sT, const double &tau_ub, const double &alpha, const double &rho_eq, const double &tol,
                                          double &maxit, double &tau_u, double &curr_pwv, double &error)
{
	const int NRHO= 1024;
	int n= CAST(int,maxit);
	const double DIGITS= std::ldexp(1,33), S2LKL= slkl*slkl, MEAN=0.;
	double tau_ubd= tau_ub/2, tau_lbd=0, *rho= NULL, *pw=NULL;
    
	rho= NEW_ARY(double,NRHO);
	pw= NEW_ARY(double,NRHO);
	for(curr_pwv=0; 	n-- && curr_pwv<S2LKL; 		)
	{
		tau_ubd*=2;
		oseq_nout(-2*tau_ubd, 2*tau_ubd, NRHO, rho);
		abcMuTOST_pow(NRHO, rho, df, tau_ubd, sT, alpha, pw);
		if(oIsZero(NRHO,pw))//must increase tau_u to get non-zero power
			curr_pwv= 0;
		else
			ovar(NRHO, rho, pw, MEAN, curr_pwv);
		//FAIL_ON(CAST(double,n+1)==maxit && curr_pwv>S2LKL,"abcMuTOST_taulowup_var: variance of power assumed to be smaller than variance of summary likelihood %c");
        //std::cout<<"H1 "<<curr_pwv<<'\t'<<S2LKL<<'\t'<<tau_ubd<<'\t'<<n<<std::endl;
		//oprinta(pw, NRHO, std::cout);
	}
	FAIL_ON(n<0,"\nabcMuTOST_taulowup_var: error at 1a %c");
	for(	error=1, n= CAST(int, maxit);
        n-- && (ABS(error)>tol) && std::floor(tau_ubd*DIGITS)!=std::floor(tau_lbd*DIGITS);
        )
	{
		tau_u= (tau_lbd+tau_ubd)/2;
		oseq_nout(-2*tau_u, 2*tau_u, NRHO, rho);
		abcMuTOST_pow(NRHO, rho, df, tau_u, sT, alpha, pw);
		if(oIsZero(NRHO,pw))//must increase tau_u to get non-zero power
			curr_pwv= 0;
		else
			ovar(NRHO, rho, pw, MEAN, curr_pwv);
		error= curr_pwv-S2LKL;
		if(error<0)
			tau_lbd= tau_u;
		else
			tau_ubd= tau_u;
        //std::cout<<"H2 "<<curr_pwv<<'\t'<<S2LKL<<'\t'<<tau_u<<'\t'<<error<<'\t'<<maxit<<'\t'<<tau_lbd<<'\t'<<tau_ubd<<'\t'<<tol<<'\t'<< (std::floor(tau_ubd*DIGITS)!=std::floor(tau_lbd*DIGITS)) <<std::endl;
	}
    //oprinta(pw, NRHO, std::cout);
	if(n<0)
		oprintff("\nabcMuTOST_tau_taulowup_var: reached max it %i current pw variance %g requ pw variance %g",n,curr_pwv,S2LKL);
	maxit= n+1;
	DELETE(rho);
	DELETE(pw);
}

void abcMUTOST_pwvar(	const double &mxpw, const double &df, const double &sT, const double &tau_ub, const double &alpha, const double &rho_eq, const double &tol,
                                   double &maxit, double &curr_pwv)
{
	const int NRHO= 1024;
	const double MEAN= 0.;
	double curr_pw= 1., tau_u=1., error=0., *rho= NULL, *pw=NULL;
    
	rho= NEW_ARY(double,NRHO);
	pw= NEW_ARY(double,NRHO);
    
	abcMuTOST_taulowup_pw(mxpw, df, sT, tau_ub, alpha, rho_eq, tol, maxit, tau_u, curr_pw, error);
	oseq_nout(-tau_u, tau_u, NRHO, rho);
	abcMuTOST_pow(NRHO, rho, df, tau_u, sT, alpha, pw);
	if(oIsZero(NRHO,pw))//must increase nsim to get non-zero power
		curr_pwv= R_PosInf;
	else
		ovar(NRHO, rho, pw, MEAN, curr_pwv);
	DELETE(rho);
	DELETE(pw);
}

void abcMuTOST_nsim(	const double &nobs, const double &slkl, const double &mxpw, const double &sSim, const double &tau_ub, const double &alpha, const double &rho_eq, const double &tol,
                                  double &maxit, double &nsim, double &tau_u, double &curr_pwv, double &curr_pw, double &error)
{
	const int NRHO= 1024;
	int n= CAST(int,maxit);
	const double S2LKL= slkl*slkl, MEAN=0.;
	double nsim_lb=nobs-1, nsim_ub=std::ceil(nobs/2) /* ceiling to make sure that nsim_ub>=nobs*/, xtau_ub= tau_ub, *rho= NULL, *pw=NULL, *cali_tau=NULL;
    
	rho= NEW_ARY(double,NRHO);
	pw= NEW_ARY(double,NRHO);
	cali_tau= NEW_ARY(double,2);
    //std::cout<<"\nnobs"<<nobs<<"\nslkl"<<slkl<<"\nmxpw"<<mxpw<<"\nsSim"<<sSim<<"\ntau_ub"<<tau_ub<<"\nalpha"<<alpha<<"\nrho_eq"<<rho_eq<<"\ntol"<<tol<<"\nmaxit"<<maxit<<std::endl;
	for(	curr_pwv=2*S2LKL;
        n-- && curr_pwv>S2LKL;
        )
	{
		nsim_ub*= 2;
		*(cali_tau+1)= maxit;
		abcMuTOST_taulowup_pw(mxpw, nsim_ub-1, sSim/std::sqrt(nsim_ub), xtau_ub, alpha, rho_eq, tol, *(cali_tau+1) /*maxit*/, tau_u /*tau_u*/, curr_pw, *cali_tau /*error*/);
		oseq_nout(-tau_u, tau_u, NRHO, rho);
		abcMuTOST_pow(NRHO, rho, nsim_ub-1, tau_u, sSim/std::sqrt(nsim_ub), alpha, pw);
		if(oIsZero(NRHO,pw))//must increase nsim to get non-zero power
			curr_pwv= 2*S2LKL;
		else
			ovar(NRHO, rho, pw, MEAN, curr_pwv);
        //std::cout<<"H1 "<<curr_pwv<<'\t'<<S2LKL<<'\t'<<tau_u<<'\t'<<nsim_ub<<std::endl;
		//FAIL_ON(CAST(double,n+1)==maxit && curr_pwv<S2LKL,"abcMuTOST_nsim: variance of power assumed to be larger than variance of summary likelihood %c");
	}
	FAIL_ON(n<0,"\nabcMuTOST_nsim: error at 1a %c");
	for(	error=1, nsim= nsim_ub, n= CAST(int,maxit);
        n-- && (ABS(error)>tol) && (nsim_lb+1)!=nsim_ub;
        )
	{
        nsim= std::floor( (nsim_lb+nsim_ub)/2 );
        *(cali_tau+1)= maxit;
        xtau_ub= tau_u;
        abcMuTOST_taulowup_pw(mxpw, nsim-1, sSim/std::sqrt(nsim), xtau_ub, alpha, rho_eq, tol, *(cali_tau+1) /*maxit*/, tau_u /*tau_u*/, curr_pw, *cali_tau /*error*/);
        oseq_nout(-tau_u, tau_u, NRHO, rho);
        abcMuTOST_pow(NRHO, rho, nsim-1, tau_u, sSim/std::sqrt(nsim), alpha, pw);
        if(oIsZero(NRHO,pw))//must increase nsim to get non-zero power
            curr_pwv= 2*S2LKL;
        else
            ovar(NRHO, rho, pw, MEAN, curr_pwv);
        
        error= curr_pwv-S2LKL;
        if(error<0)
            nsim_ub= nsim;
        else
            nsim_lb= nsim;
        //std::cout<<"H2 "<<curr_pwv<<'\t'<<S2LKL<<'\t'<<tau_u<<'\t'<<nsim<<'\t'<<curr_pw<<'\t'<<nsim_lb<<'\t'<<nsim_ub<<std::endl;
	}
	if(n<0)
		oprintff("\nabcMuTOST_nsim: reached max it %i current pw variance %g requ pw variance %g",n,curr_pwv,S2LKL);
    
	maxit= n+1;
	DELETE(rho);
	DELETE(pw);
	DELETE(cali_tau);
}

void abcMuTOST_KL(const double &nx, const double &sx, const double &ny, const double &sy, const double &mx_pw, const double &alpha, const int &calibrate_tau_up, double &tau_up,const double &pow_scale, double &curr_mxpw, double &KL_div)
{
    
    //calibrate tau_up
    if(calibrate_tau_up)
    {
        double *cali_tau= NEW_ARY(double,5);
        *(cali_tau+1)= tau_up;//tau_ub
        *(cali_tau+2)= 0;//rho_eq
        *(cali_tau+3)= 1e-5;//tol
        *(cali_tau+4)= 100;//maxit
        
        //std::cout<<"abcMuTOST_taulowup_pw at a1:\t"<<tau_up<<'\t'<<curr_mxpw<<'\t'<<*cali_tau<<std::endl;
        
        abcMuTOST_taulowup_pw(mx_pw, ny-1, sy/std::sqrt(ny), *(cali_tau+1), alpha, *(cali_tau+2), *(cali_tau+3), *(cali_tau+4), tau_up/*tau_up*/, curr_mxpw/*curr_mxpw*/, *cali_tau/*error*/);
        
        //std::cout<<"abcMuTOST_taulowup_pw at a2:\t"<<tau_up<<'\t'<<curr_mxpw<<'\t'<<*cali_tau<<std::endl;
        
        DELETE(cali_tau);
    }
    
    //compute pow_norm
    int neval;
    double lower,upper,abs_tol,rel_tol,abserr;
    upper=tau_up*pow_scale;
    lower=-upper;
    rel_tol=nabcGlobals::NABC_DBL_TOL;
    abs_tol=rel_tol;
    
    //create arg for pow
    basic_arg pow_arg;
    pow_arg.df=ny-1;
    pow_arg.sT=sy/std::sqrt(ny);
    pow_arg.tau_up=tau_up;
    pow_arg.alpha=alpha;
    pow_arg.norm=1;
    pow_arg.give_log=0;
    
    //compute actual max power if tau.u is not calibrated
    if (!calibrate_tau_up) {
        curr_mxpw = abcMuTOST_pow_scalar(0,&pow_arg);
    }
    
    //std::cout<<"abcMuTOST_taulowup_pw at b1:\t"<<lower<<'\t'<<pow_arg.norm<<'\t'<<pow_arg.give_log<<'\t'<<neval<<std::endl;
    
    //printBArg(&pow_arg);
    nabc_integration_qng(abcMuTOST_pow_scalar, &pow_arg, lower, upper, abs_tol, rel_tol, &pow_arg.norm /*updated norm*/, &abserr, &neval);
    
    //std::cout<<"abcMuTOST_taulowup_pw at b2:\t"<<lower<<'\t'<<pow_arg.norm<<'\t'<<abserr<<'\t'<<neval<<std::endl;
    
    //create arg for sulkl
    basic_arg sulkl_arg;
    sulkl_arg.nx=nx;
    sulkl_arg.df=nx-1;
    sulkl_arg.sx=sx;
    sulkl_arg.ssn=sx/std::sqrt(nx);
    sulkl_arg.norm=pt(upper/sulkl_arg.ssn,sulkl_arg.df,1,0)-pt(lower/sulkl_arg.ssn,sulkl_arg.df,1,0);
    sulkl_arg.give_log=0;
    
    //put give_log to 1
    sulkl_arg.give_log=1;
    pow_arg.give_log=1;
    
    kl_integrand_arg KL_arg;
    KL_arg.p=&abcMuTOST_sulkl_scalar;
    KL_arg.q=&abcMuTOST_pow_scalar;
    KL_arg.p_arg=&sulkl_arg;
    KL_arg.q_arg=&pow_arg;
    
    //    printf("pow arguments\n");
    //    printBArg(&pow_arg);
    //    printf("\nsulkl arguments\n");
    //    printBArg(&sulkl_arg);
    //    printf("Start KL integration\n");
    
    nabc_integration_qng(abcKL_integrand, &KL_arg, lower, upper, abs_tol, rel_tol, &KL_div, &abserr, &neval);
    
    //std::cout<<"abcMuTOST_taulowup_pw at c2:\t"<<KL_div<<'\t'<<abserr<<'\t'<<neval<<std::endl;
    
    
}


void abcMuTOST_KL(kl_arg *arg)
{
    //std::cout<<"abcMuTOST_KL 1a:\t"<<arg->nx<<'\t'<<arg->sx<<'\t'<<arg->ny<<'\t'<<arg->sy<<'\t'<<arg->mx_pw<<'\t'<<arg->alpha<<'\t'<<arg->calibrate_tau_up<<'\t'<<arg->tau_up<<'\t'<<arg->pow_scale<<'\t'<<arg->curr_mxpw<<'\t'<<arg->KL_div<<std::endl;
    
    abcMuTOST_KL(arg->nx, arg->sx, arg->ny, arg->sy, arg->mx_pw, arg->alpha, arg->calibrate_tau_up, arg->tau_up, arg->pow_scale, arg->curr_mxpw, arg->KL_div);
    
    //std::cout<<"abcMuTOST_KL 1b:\t"<<arg->nx<<'\t'<<arg->sx<<'\t'<<arg->ny<<'\t'<<arg->sy<<'\t'<<arg->mx_pw<<'\t'<<arg->alpha<<'\t'<<arg->calibrate_tau_up<<'\t'<<arg->tau_up<<'\t'<<arg->pow_scale<<'\t'<<arg->curr_mxpw<<'\t'<<arg->KL_div<<std::endl;
    
}



