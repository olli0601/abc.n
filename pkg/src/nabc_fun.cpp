#include <iostream>
#include <cmath>
#include <Rmath.h>
#include "nabc_error_handling.h"
#include "nabc_fun.h"
#include "nabc_globals.h"

static inline void oprintf(const char * format, ...)
{
	va_list args;
	va_start(args, format);
	//Rvprintf(format, args);
	vprintf(format, args);
	va_end(args);
}

static inline void oprintff(const char * format, ...)
{
	va_list args;
	va_start(args, format);
	//Rvprintf(format, args);
	vprintf(format, args);
	va_end(args);
	fflush(stdout);
}

static inline void oprinta(double const * const start,const int &n,std::ostream& os)
{
	double const * xfltty= start;
	int m= n;
    
	os<<"c(";
	if(m>1)
		for(m--;m--;)
			os<<*xfltty++<<", ";
	if(n>0)
		os<<*xfltty;
	os<<")";
}

static inline void oseq_nout(const double &a, const double &b, const int &n, double * const ans)
{
	FAIL_ON(n<2,"oseq_nout: error at 1a %c");
	double *xans=NULL, *yans=NULL;
	int xn= n-1;
	const double by= (b-a)/xn;
    
	yans= xans= ans;
	*xans++= a;
	for(; 		xn--; 		*xans++= *yans++ + by);
}

static inline void ovar(const int &n, double * const x, double * const fx, const double &mean, double &var)
{
	FAIL_ON(n<1,"ovar: error at 1a %c");
	int xn= n;
	double *xx= x, *xfx= fx, norm=0;
    
	if(mean==0)
		for(var=0, norm=0;		xn--;		xx++)
		{
			var+= *xx * *xx * *xfx;
			norm+= *xfx++;//assume that x is equidistant
		}
	else
		for(var=0, norm=0;		xn--;		xx++)
		{
			var+= (*xx-mean) * (*xx-mean) * *xfx;
			norm+= *xfx++;//assume that x is equidistant
		}
	var/= norm;
}


static inline int oIsZero(const int &n, double * const x)
{
	int xn= n, isZero=1;
	double *xx= x;
	for(; 		isZero && xn--;		xx++)
		isZero= CAST(int,*xx==0);
    //isZero= (	*xx > std::numeric_limits<double>::epsilon() || *xx < -std::numeric_limits<double>::epsilon()	)?0:1;
    //isZero= CAST(int,*xx==0);
	return isZero;
}

static inline void abcMuTOST_pow(const int &nrho, double * const rho, const double &df, const double &tau_up, const double &sT, const double &alpha, double *ans)
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


static inline void abcMuTOST_sulkl(const int &nrho, double * const rho, const double &nx, const double &sx, const double &norm, const int &give_log,double *ans)
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

double abcKL_integrand(double x,void *arg_void)
{
    kl_integrand_arg *arg=(kl_integrand_arg *) arg_void;
    const double log_P=(*(arg->p))(x,arg->p_arg);
    const double log_Q=(*(arg->q))(x,arg->q_arg);
    
    return (log_P-log_Q)*exp(log_P);
    
}




static inline void abcScaledChiSq_criticalregion(	const double scale, const double df, const double tl, const double tu, const double alpha, const double tol, const double incit,
                                                 double &maxit, double &c1, double &c2, double &mxpw, double &error	)
{
	const int LOWERTAIL= 1;
	const double tuisc= tu/scale, scitu= scale/tu, scitl= scale/tl;
	double c1l= 0, c1r= 0, level= 0;
	//std::cout<<"abcScaledChiSq input:\nscale\t"<<scale<<"\ndf\t"<<df<<"\ntl\t"<<tl<<"\ntu\t"<<tu<<"\nalpha\t"<<alpha<<"\ntol\t"<<tol<<"\nmaxit\t"<<maxit<<"\nincit\t"<<incit<<std::endl;
    
	//determine [c1l,c2] such that p(R|tl)>alpha --> then, [c1l,c1l+incit] contains the true c1.
	//start with [c1,c2] sth p(R|tl)<alpha:		c1l= std::sqrt( tl*tu )		and substract first increment
	for(level= 0, c1l= std::sqrt( tl*tu )+incit;			level<alpha;			)
	{
        c1l-=incit;
        c2=	alpha + pchisq(c1l*scitu, df, LOWERTAIL, 0);
        c2= tuisc * qchisq(c2, df, LOWERTAIL, 0);
        level= pchisq(c2*scitl, df, LOWERTAIL, 0) - pchisq(c1l*scitl, df, LOWERTAIL, 0);
        //std::cout<<c1l<<'\t'<<c2<<'\t'<<level<<'\t'<<(level-alpha)<<std::endl; std::cout<<"OK1"<<std::endl;
	}
    //initialize c1 in [c1l,c2], set to c1l for simplicity.
	for(	c1= c1l, c1r= c1l+incit;
        (ABSDIFF(level,alpha))>=tol	&&	maxit--;
        )
	{
        c1= (c1l + c1r) * 0.5;
        c2= tuisc * qchisq(		alpha + pchisq(c1*scitu, df, LOWERTAIL, 0), 			df, LOWERTAIL, 0);
        level= pchisq(c2*scitl, df, LOWERTAIL, 0) - pchisq(c1*scitl, df, LOWERTAIL, 0);
        if(level <= alpha)
            c1r= c1;	//if interval too much on the right, then update C1R
        else
            c1l= c1;	//if interval too much on the left, then update C1L
        //std::cout<<c1<<'\t'<<c2<<'\t'<<c1r<<'\t'<<c1l<<'\t'<<ABSDIFF(level,alpha)<<'\t'<<tol<<std::endl;	std::cout<<"OK2"<<std::endl;
	}
    
	//output
	//c1 updated
	//c2 updated
	mxpw= pchisq(c2*scale,df,LOWERTAIL,0) - pchisq(c1*scale,df,LOWERTAIL,0);//power at rho=1
	error= ABSDIFF(level,alpha);
	//maxit updated
}

static inline void abcMuTOST_taulowup_pw(	const double &mxpw, const double &df, const double &sT, const double &tau_ub, const double &alpha, const double &rho_eq, const double &tol,
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

static inline void abcMuTOST_taulowup_var(	const double &slkl, const double &df, const double &sT, const double &tau_ub, const double &alpha, const double &rho_eq, const double &tol,
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

static inline void abcMUTOST_pwvar(	const double &mxpw, const double &df, const double &sT, const double &tau_ub, const double &alpha, const double &rho_eq, const double &tol,
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

static inline void abcMuTOST_nsim(	const double &nobs, const double &slkl, const double &mxpw, const double &sSim, const double &tau_ub, const double &alpha, const double &rho_eq, const double &tol,
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

SEXP abcScaledChiSq(SEXP args)
{
	FAIL_ON(! Rf_isReal(args) ,"abcScaledChiSqu: error at 1a %c");
	double scale=0, df=0, tl= 1, tu= 1, alpha= 0.01, tol= 1e-10, maxit= 100, incit= 0.05;
	double *xans= NULL;
	SEXP ans;
    
	//convert SEXP into C
	FAIL_ON(length(args)!=8,"abcScaledChiSqu: error at 1b %c");
	FAIL_ON((scale= REAL(args)[0])<3,"abcScaledChiSqu: error at 1c %c");
	FAIL_ON((df= REAL(args)[1])<2,"abcScaledChiSqu: error at 1d %c");
	FAIL_ON((tl= REAL(args)[2])>1 || tl<0,"abcScaledChiSqu: error at 1e %c");
	FAIL_ON((tu= REAL(args)[3])<1,"abcScaledChiSqu: error at 1f %c");
	FAIL_ON((alpha= REAL(args)[4])>1 || alpha<0,"abcScaledChiSqu: error at 1g %c");
	tol= REAL(args)[5];
	maxit= REAL(args)[6];
	incit= REAL(args)[7];
    
	//output
	PROTECT(ans=  allocVector(REALSXP,5));
	xans= REAL(ans);
	*(xans+4)= maxit;
    
	//compute stuff
	abcScaledChiSq_criticalregion(scale, df, tl, tu, alpha, tol, incit, *(xans+4) /*maxit will be updated*/, *xans /*c1*/, *(xans+1) /*c2*/, *(xans+2) /*mxpw*/, *(xans+3) /*error*/);
    
	UNPROTECT(1);
	return ans;
}

SEXP abcMuTOST_pow(SEXP arg_rho, SEXP arg_df, SEXP arg_tau_up, SEXP arg_sT, SEXP arg_alpha)
{
	FAIL_ON(!Rf_isReal(arg_rho) ,"abcMuTOST_pow: error at 1a %c");
    
	int nrho= Rf_length(arg_rho);
	double *xrho= REAL(arg_rho), *xans=NULL;
	double df= ::Rf_asReal(arg_df), tau_up= ::Rf_asReal(arg_tau_up), sT= ::Rf_asReal(arg_sT), alpha= ::Rf_asReal(arg_alpha);
	SEXP ans;
    
	PROTECT(ans=  allocVector(REALSXP,nrho));
	xans= REAL(ans);
	//std::cout<<"abcMuTOST_pow_c\t"<<nrho<<'\t'<<*xrho<<'\t'<<df<<'\t'<<tau_up<<'\t'<<sT<<'\t'<<alpha<<std::endl;
	abcMuTOST_pow(nrho, xrho, df, tau_up, sT, alpha, xans);
    
	UNPROTECT(1);
	return ans;
}

SEXP abcMuTOST_sulkl_integrate_qng(SEXP arg_lower, SEXP arg_upper, SEXP arg_abs_tol, SEXP arg_rel_tol, SEXP arg_nx, SEXP arg_sx, SEXP arg_norm, SEXP arg_log)
{
    
    basic_arg f_arg;
    f_arg.nx=asReal(arg_nx);
    f_arg.df=f_arg.nx-1;
    f_arg.sx=asReal(arg_sx);
    f_arg.ssn=f_arg.sx/sqrt(f_arg.nx);
    f_arg.norm=asReal(arg_norm);
    f_arg.give_log=asInteger(arg_log);
    
    
	double lower= asReal(arg_lower),upper= asReal(arg_upper), abs_tol=asReal(arg_abs_tol),rel_tol=asReal(arg_rel_tol) ;
    double abserr;
    double *xans=NULL;
    int neval,res;
    
	SEXP ans;
    
	PROTECT(ans=  allocVector(REALSXP,1));
	xans= REAL(ans);
	res=nabc_integration_qng(abcMuTOST_sulkl_scalar, &f_arg, lower, upper, abs_tol, rel_tol, xans, &abserr, &neval);
	
    //std::cout<<"abcMuTOST_sulkl_integrate_qng\nnx\t"<<f_arg.nx<<"\nsx\t"<<f_arg.sx<<"\ndf\t"<<f_arg.df<<"\nssn\t"<<f_arg.ssn<<"\nnorm\t"<<f_arg.norm<<"\ngive_log\t"<<f_arg.give_log<<"\nreturn\t"<<res<<"\nans\t"<<*xans<<std::endl;
    
	UNPROTECT(1);
	return ans;
}


SEXP abcMuTOST_pow_integrate_qng(SEXP arg_lower, SEXP arg_upper, SEXP arg_abs_tol, SEXP arg_rel_tol, SEXP arg_df, SEXP arg_sT, SEXP arg_tau_up, SEXP arg_alpha, SEXP arg_norm, SEXP arg_log)
{
    
    basic_arg f_arg;
    f_arg.df=asReal(arg_df);
    f_arg.sT=asReal(arg_sT);
    f_arg.tau_up=asReal(arg_tau_up);
    f_arg.alpha=asReal(arg_alpha);
    f_arg.norm=asReal(arg_norm);
    f_arg.give_log=asInteger(arg_log);
    
    
	double lower= asReal(arg_lower),upper= asReal(arg_upper), abs_tol=asReal(arg_abs_tol),rel_tol=asReal(arg_rel_tol) ;
    double abserr;
    double *xans=NULL;
    int neval,res;
    
	SEXP ans;
    
	PROTECT(ans=  allocVector(REALSXP,1));
	xans= REAL(ans);
	res=nabc_integration_qng(abcMuTOST_pow_scalar, &f_arg, lower, upper, abs_tol, rel_tol, xans, &abserr, &neval);
	
    //std::cout<<"abcMuTOST_pow_integrate_qng\ndf\t"<<f_arg.df<<"\nsT\t"<<f_arg.sT<<"\ntau_up\t"<<f_arg.tau_up<<"\nalpha\t"<<f_arg.alpha<<"\nnorm\t"<<f_arg.norm<<"\ngive_log\t"<<f_arg.give_log<<"\nreturn\t"<<res<<"\nans\t"<<*xans<<std::endl;
    
	UNPROTECT(1);
	return ans;
}


SEXP abcMuTOST_sulkl(SEXP arg_rho, SEXP arg_nx, SEXP arg_sx, SEXP arg_norm, SEXP arg_log)
{
	int nrho= Rf_length(arg_rho),give_log=asInteger(arg_log);
	double *xrho= REAL(arg_rho), *xans=NULL;
	double nx= asReal(arg_nx), sx= asReal(arg_sx), norm= asReal(arg_norm);
	SEXP ans;
    
	PROTECT(ans=  allocVector(REALSXP,nrho));
	xans= REAL(ans);
	//std::cout<<"abcMuTOST_sulkl_c\t"<<nrho<<'\t'<<*xrho<<'\t'<<nx<<'\t'<<sx<<'\t'<<norm<<'\t'<<log<<std::endl;
	abcMuTOST_sulkl(nrho, xrho, nx, sx, norm, give_log, xans);
    
	UNPROTECT(1);
	return ans;
}


SEXP abcMuTOST_pwvar(SEXP args)
{
	FAIL_ON(!Rf_isReal(args) ,"abcMuTOST_pwvar: error at 1a %c");
	double mxpw=0, df=0, sT= 1, tu_ub= 1, alpha= 0.01, rho_eq=0, tol= 1e-10, maxit= 100;
	double *xans= NULL;
	SEXP ans;
    
	//convert SEXP into C
	FAIL_ON(length(args)!=8,					"abcMuTOST_pwvar: error at 1b %c");
	FAIL_ON((mxpw= REAL(args)[0])>1 || mxpw<=0,	"abcMuTOST_pwvar: error at 1c %c");
	FAIL_ON((df= REAL(args)[1])<2,				"abcMuTOST_pwvar: error at 1d %c");
	FAIL_ON((sT= REAL(args)[2])<=0,				"abcMuTOST_pwvar: error at 1e %c");
	FAIL_ON((tu_ub= REAL(args)[3])<0,			"abcMuTOST_pwvar: error at 1f %c");
	FAIL_ON((alpha= REAL(args)[4])>1 || alpha<0,"abcMuTOST_pwvar: error at 1g %c");
	rho_eq= REAL(args)[5];
	tol= REAL(args)[6];
	maxit= REAL(args)[7];
    
	PROTECT(ans=  allocVector(REALSXP,2));
	xans= REAL(ans);
	*(xans+1)= maxit;
	abcMUTOST_pwvar(mxpw, df, sT, tu_ub, alpha, rho_eq, tol, *(xans+1) /*maxit*/, *xans/*curr_pwv*/);
    
	UNPROTECT(1);
	return ans;
}

void printBArg(basic_arg *arg)
{
    printf("nx=%f\nsx=%f\nssn=%f\ndf=%f\nnorm=%f\ngive_log=%d\nalpha=%f\ntau_up=%f\nsT=%f\n",arg->nx,arg->sx,arg->ssn,arg->df,arg->norm,arg->give_log,arg->alpha,arg->tau_up,arg->sT);
    
}

static inline void abcMuTOST_KL(const double &nx, const double &sx, const double &ny, const double &sy, const double &mx_pw, const double &alpha, const int &calibrate_tau_up, double &tau_up,const double &pow_scale, double &curr_mxpw, double &KL_div)
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
    rel_tol=std::pow(nabcGlobals::NABC_DBL_EPSILON,0.25);
    abs_tol=rel_tol;
    
    //create arg for pow
    basic_arg pow_arg;
    pow_arg.df=ny-1;
    pow_arg.sT=sy/std::sqrt(ny);
    pow_arg.tau_up=tau_up;
    pow_arg.alpha=alpha;
    pow_arg.norm=1;
    pow_arg.give_log=0;
    
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
    KL_arg.p=abcMuTOST_sulkl_scalar;
    KL_arg.q=abcMuTOST_pow_scalar;
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


static inline void abcMuTOST_KL(kl_arg *arg)
{
    abcMuTOST_KL(arg->nx, arg->sx, arg->ny, arg->sy, arg->mx_pw, arg->alpha, arg->calibrate_tau_up, arg->tau_up, arg->pow_scale, arg->curr_mxpw, arg->KL_div);
}


SEXP abcMuTOST_KL(SEXP arg_nx, SEXP arg_sx, SEXP arg_ny, SEXP arg_sy, SEXP arg_mx_pw, SEXP arg_alpha, SEXP arg_calibrate_tau_up, SEXP arg_tau_up, SEXP arg_pow_scale)
{
    //SEXP to C
    double nx=asReal(arg_nx);
    double sx=asReal(arg_sx);
    double ny=asReal(arg_ny);
    double sy=asReal(arg_sy);
    double mx_pw=asReal(arg_mx_pw);
    double alpha=asReal(arg_alpha);
    int calibrate_tau_up=asLogical(arg_calibrate_tau_up);
    double tau_up=asReal(arg_tau_up);
    double pow_scale=asReal(arg_pow_scale);
    
    //answer
    double *xans= NULL;
	SEXP ans;
	PROTECT(ans= allocVector(REALSXP,3));
	xans= REAL(ans);
    xans[1]=tau_up;
    xans[2]=mx_pw;
    
    //Call C function
    abcMuTOST_KL(nx,sx, ny, sy,mx_pw,alpha, calibrate_tau_up, xans[1], pow_scale, xans[2],xans[0]);
    
	UNPROTECT(1);
	return ans;
    
}

SEXP getListElement(SEXP list, const char *str)
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    
    for (int i = 0; i < length(list); i++)
        if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
            elmt = VECTOR_ELT(list, i);
            break;
        }
    return elmt;
}


SEXP abcCalibrate_tau_nomxpw_yesKL(SEXP arg_test_name, SEXP list_KL_args, SEXP arg_tau_up_lb, SEXP arg_max_it)
{
    //SEXP to C
    //    double nx=asReal(arg_nx);
    //    double sx=asReal(arg_sx);
    //    double ny=asReal(arg_ny);
    //    double sy=asReal(arg_sy);
    //    double mx_pw=asReal(arg_mx_pw);
    //    double alpha=asReal(arg_alpha);
    //    int calibrate_tau_up=asLogical(arg_calibrate_tau_up);
    //    double tau_up=asReal(arg_tau_up);
    //    double pow_scale=asReal(arg_pow_scale);

    EqTestValue eq_test_val;

    //get element from list_KL_args and put them in a kl_arg structure
    //by default just add the arguments you need, un-needed arguments will just be NULL if not present in list_KL_args
    kl_arg arg;
    arg.nx = asReal(getListElement(list_KL_args,"n.of.x"));
    arg.sx = asReal(getListElement(list_KL_args,"s.of.x"));
    arg.ny = asReal(getListElement(list_KL_args,"n.of.y"));
    arg.sy = asReal(getListElement(list_KL_args,"s.of.y"));
    arg.mx_pw = asInteger(getListElement(list_KL_args,"mx.pw"));
    arg.alpha = asReal(getListElement(list_KL_args,"alpha"));
    
    strcpy(nabcGlobals::BUFFER,CHAR(STRING_ELT(arg_test_name, 0)));
    
    try{
        eq_test_val=s_mapEqTestValue.at(nabcGlobals::BUFFER);
    }
    catch (const std::out_of_range& oor) {
        error("%s is not in the lookup table of equivalence test, see documentation\n",nabcGlobals::BUFFER);
    }
    
    //given the arg_test_name choose right function and ptr to it
    switch (eq_test_val) {
        case MUTOST_ONE_SAMPLE:
            
            break;
            
        default:
            error("%s:%s:%d: edit switch to include test %s",__FILE__,__FUNCTION__,__LINE__,nabcGlobals::BUFFER);
            break;
    }
    
    
    //answer
    double *xans= NULL;
    SEXP ans;
    PROTECT(ans= allocVector(REALSXP,3));
    xans= REAL(ans);
    xans[1]=arg.nx;
    xans[2]=arg.ny;
    
    //Call C function
    //abcMuTOST_KL(nx,sx, ny, sy,mx_pw,alpha, calibrate_tau_up, xans[1], pow_scale, xans[2],xans[0]);
    
    UNPROTECT(1);
    return ans;
    
}



SEXP abcMuTOST_taulowup_pw(SEXP args)
{
	FAIL_ON(!Rf_isReal(args) ,"abcMuTOST_tau_taulowup_pw: error at 1a %c");
	double mxpw=0, df=0, sT= 1, tu_ub= 1, alpha= 0.01, rho_eq=0, tol= 1e-10, maxit= 100;
	double *xans= NULL;
	SEXP ans;
    
	//convert SEXP into C
	FAIL_ON(length(args)!=8,					"abcMuTOST_tau_taulowup_pw: error at 1b %c");
	FAIL_ON((mxpw= REAL(args)[0])>1 || mxpw<=0,	"abcMuTOST_tau_taulowup_pw: error at 1c %c");
	FAIL_ON((df= REAL(args)[1])<2,				"abcMuTOST_tau_taulowup_pw: error at 1d %c");
	FAIL_ON((sT= REAL(args)[2])<=0,				"abcMuTOST_tau_taulowup_pw: error at 1e %c");
	FAIL_ON((tu_ub= REAL(args)[3])<0,			"abcMuTOST_tau_taulowup_pw: error at 1f %c");
	FAIL_ON((alpha= REAL(args)[4])>1 || alpha<0,"abcMuTOST_tau_taulowup_pw: error at 1g %c");
	rho_eq= REAL(args)[5];
	tol= REAL(args)[6];
	maxit= REAL(args)[7];
    
	PROTECT(ans=  allocVector(REALSXP,5));
	xans= REAL(ans);
	*(xans+4)= maxit;
	abcMuTOST_taulowup_pw(mxpw, df, sT, tu_ub, alpha, rho_eq, tol, *(xans+4) /*maxit*/, *(xans+1)/*tau_u*/, *(xans+2)/*curr_mxpw*/, *(xans+3)/*error*/);
	*xans= -*(xans+1);
    
	UNPROTECT(1);
	return ans;
}

SEXP abcMuTOST_taulowup_var(SEXP args)
{
	FAIL_ON(!Rf_isReal(args) ,"abcMuTOST_tau_taulowup_var: error at 1a %c");
	double sLkl=0, df=0, sT= 1, tu_ub= 1, alpha= 0.01, rho_eq=0, tol= 1e-10, maxit= 100;
	double *xans= NULL;
	SEXP ans;
    
	//convert SEXP into C
	FAIL_ON(length(args)!=8,					"abcMuTOST_tau_taulowup_var: error at 1b %c");
	FAIL_ON((sLkl= REAL(args)[0])<=0,			"abcMuTOST_tau_taulowup_var: error at 1c %c");
	FAIL_ON((df= REAL(args)[1])<2,				"abcMuTOST_tau_taulowup_var: error at 1d %c");
	FAIL_ON((sT= REAL(args)[2])<=0,				"abcMuTOST_tau_taulowup_var: error at 1e %c");
	FAIL_ON((tu_ub= REAL(args)[3])<0,			"abcMuTOST_tau_taulowup_var: error at 1f %c");
	FAIL_ON((alpha= REAL(args)[4])>1 || alpha<0,"abcMuTOST_tau_taulowup_var: error at 1g %c");
	rho_eq= REAL(args)[5];
	tol= REAL(args)[6];
	maxit= REAL(args)[7];
    
	PROTECT(ans=  allocVector(REALSXP,5));
	xans= REAL(ans);
	*(xans+4)= maxit;
	abcMuTOST_taulowup_var(sLkl, df, sT, tu_ub, alpha, rho_eq, tol, *(xans+4) /*maxit*/, *(xans+1)/*tau_u*/, *(xans+2)/*curr_pwv*/, *(xans+3)/*error*/);
	*xans= -*(xans+1);
    
	UNPROTECT(1);
	return ans;
}

SEXP abcMuTOST_nsim(SEXP args)
{
	FAIL_ON(!Rf_isReal(args) ,"abcMuTOST_nsim: error at 1a %c");
	double nobs= 1, sLkl=0, mxpw=0, sSim=0, tu_ub= 1, alpha= 0.01, rho_eq=0, tol= 1e-10, maxit= 100;
	double *xans= NULL;
	SEXP ans;
    
	//convert SEXP into C
	FAIL_ON(length(args)!=9,					"abcMuTOST_nsim: error at 1b %c");
	FAIL_ON((nobs= REAL(args)[0])<=0,			"abcMuTOST_nsim: error at 1c %c");
	FAIL_ON((sLkl= REAL(args)[1])<=0,			"abcMuTOST_nsim: error at 1d %c");
	FAIL_ON((mxpw= REAL(args)[2])<=0 || mxpw>1,	"abcMuTOST_nsim: error at 1e %c");
	FAIL_ON((sSim= REAL(args)[3])<=0,			"abcMuTOST_nsim: error at 1f %c");
	FAIL_ON((tu_ub= REAL(args)[4])<0,			"abcMuTOST_nsim: error at 1f %c");
	FAIL_ON((alpha= REAL(args)[5])>1 || alpha<0,"abcMuTOST_nsim: error at 1g %c");
	rho_eq= REAL(args)[6];
	tol= REAL(args)[7];
	maxit= REAL(args)[8];
    
	PROTECT(ans=  allocVector(REALSXP,7));
	xans= REAL(ans);
	*(xans+6)= maxit;
	abcMuTOST_nsim(nobs, sLkl, mxpw, sSim, tu_ub, alpha, rho_eq, tol, *(xans+6) /*maxit*/, *xans /*nsim*/,*(xans+2)/*tau_u*/, *(xans+3)/*curr_pwv*/, *(xans+4)/*curr_pw*/, *(xans+5)/*error*/);
	*(xans+1)= -*(xans+2);
    
	UNPROTECT(1);
	return ans;
}

SEXP abcIntersectLevelSets(SEXP m1, SEXP m2, SEXP s)
{
	FAIL_ON(! Rf_isMatrix(m1) ,"abcIntersectLevelSets: error at 1a %c");
	FAIL_ON(! Rf_isMatrix(m2) ,"abcIntersectLevelSets: error at 1b %c");
	FAIL_ON(! Rf_isReal(s) ,"abcIntersectLevelSets: error at 1c %c");
    
	int i,j1,j2;
	int nr1= Rf_nrows(m1), nr2= Rf_nrows(m2), nc1= Rf_ncols(m1), nc2= Rf_ncols(m2), ns= Rf_length(s);
	double tmp, *xd= NULL, *ym2=NULL, *xm1= NULL, *xm2= NULL, *xs= NULL, *is=NULL;
	SEXP ans;
    
	FAIL_ON(nr1!=nr2,"abcIntersectLevelSets: error at 2a %c");
	FAIL_ON(nr1!=ns,"abcIntersectLevelSets: error at 2b %c");
    
	PROTECT(ans=  Rf_allocMatrix(REALSXP, nc1, nc2));
	for(i= nc1*nc2, xd= REAL(ans); 	i--;  *xd++= 0.);
    
	is= NEW_ARY(double,ns);
	for(i=ns, xd=is, xs=REAL(s);		i--;	*xd++= 1 / *xs++);
    
	//std::cout<<"n"<<nr1<<"col"<<nc1<<"\t"<<*xm1<<'\t'<<*(xm1+1)<<'\t'<<*(xm1+2)<<'\t'<<*(xm1+3)<<std::endl;
    
	//xd= (dist= oNFLTYM_NewFill(nc2,nc1,0.,error))->rdata;	//rdata is byrow
	for(j2= nc2, xm2= REAL(m2), xd= REAL(ans);		j2--;		xm2+=nr1)
		for(j1= nc1, xm1= REAL(m1);		j1--;		xd++)
			for(i=nr1, xs= is, ym2=xm2;		i--;	)
			{
				tmp= (*xm1++ - *ym2++) * *xs++;
				(*xd)+= tmp * tmp;
			}
	DELETE(is);
	UNPROTECT(1);
	return ans;
}
