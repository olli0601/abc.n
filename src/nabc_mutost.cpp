//
//  nabc_mutost.cpp
//  nABC
//
//
//
#include "nabc_utils.h"
#include "nabc_mutost.h"
#include "nabc_Rcall.h"

void printBArg(arg_mutost *arg)
{
	printf("nx=%f\nsx=%f\nssn=%f\ndf=%f\nnorm=%f\ngive_log=%d\nalpha=%f\ntau_up=%f\nsT=%f\n",arg->nx,arg->sx,arg->ssn,arg->df,arg->norm,arg->give_log,arg->alpha,arg->tau_up,arg->sT);

}

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

void abcMuTOST_power(const int &nrho, double * const rho, const double &df, const double &c_up, const double &sT, double *ans)
{
	int n= nrho;
	const int LOG= 0, LOWERTAIL=1;
	double *xans= ans, *xrho=rho;
	const double CIS=c_up/sT;

	for(; n--; xrho++, xans++)
	{
		*xans= pnt( CIS, df, *xrho/sT, LOWERTAIL, LOG) - pnt( -CIS, df, *xrho/sT, LOWERTAIL, LOG);
		*xans= *xans<0 ? 0 : *xans;
	}
}

double abcMuTOST_power_scalar(double x, void *arg)
{
	arg_mutost *a=(arg_mutost *) arg;
	double ans;
	const double CIS= a->tau_up/a->sT;	//we set a->tau_up to c_up

	ans= pnt( CIS, a->df, x/a->sT, 1, 0) - pnt( -CIS, a->df, x/a->sT, 1, 0);

    /* Can be 0 at some point! */
    /* Usually this happens due to numerical inacurracy in the tail. */
    /* To avoid infinity we assume that ans=nabcGlobals::NABC_DBL_MIN at these points. */

	ans= ans<=0 ? nabcGlobals::NABC_DBL_MIN : ans/a->norm;

	return a->give_log ? log(ans) : ans;
}

double abcMuTOST_pow_scalar(double x, void *arg)
{
	arg_mutost *a=(arg_mutost *) arg;
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
	const double ssn= sx/sqrt(nx);
	double *xans= ans, *xrho=rho;

	for(; n--; xrho++, xans++)
	{
		*xans= dnorm( *xrho, 0, ssn, give_log);
		*xans= give_log ? *xans-log(norm) : *xans/norm;
	}
}

static inline double abcMuTOST_sulkl_scalar(double x, void *arg_void)
{
	arg_mutost *arg=(arg_mutost *) arg_void;
    //std::cout<<"abcMuTOST_sulkl_scalar input:\nssn\t"<<arg->ssn<<"\ndf\t"<<arg->df<<"\ngive_log\t"<<arg->give_log<<"\nnorm\t"<<arg->norm<<"\nx\t"<<x<<std::endl;

	double ans;

	//ans= dt( x/arg->ssn, arg->df, arg->give_log);
	//ans= arg->give_log ? ans-log(arg->ssn*arg->norm) : ans/arg->ssn/arg->norm;
	ans= dnorm( x, 0, arg->ssn, arg->give_log);
	ans= arg->give_log ? ans-log(arg->norm) : ans/arg->norm;
	
	return ans;
}

SEXP abcMuTOST_pow(SEXP arg_rho, SEXP arg_df, SEXP arg_tau_up, SEXP arg_sT, SEXP arg_alpha)
{
	ERROR_ON(!Rf_isReal(arg_rho) ,"abcMuTOST_pow: error at 1a ");

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

SEXP abcMuTOST_power(SEXP arg_rho, SEXP arg_df, SEXP arg_c_up, SEXP arg_sT)
{
	ERROR_ON(!Rf_isReal(arg_rho) ,"abcMuTOST_power: error at 1a ");

	int nrho= Rf_length(arg_rho);
	double *xrho= REAL(arg_rho), *xans=NULL;
	double df= ::Rf_asReal(arg_df), c_up= ::Rf_asReal(arg_c_up), sT= ::Rf_asReal(arg_sT);
	SEXP ans;

	PROTECT(ans=  allocVector(REALSXP,nrho));
	xans= REAL(ans);
	//std::cout<<"abcMuTOST_pow_c\t"<<nrho<<'\t'<<*xrho<<'\t'<<df<<'\t'<<tau_up<<'\t'<<sT<<std::endl;
	abcMuTOST_power(nrho, xrho, df, c_up, sT, xans);

	UNPROTECT(1);
	return ans;
}

SEXP abc_mutost_integrate_sulkl(SEXP arg_lower, SEXP arg_upper, SEXP arg_abs_tol, SEXP arg_rel_tol, SEXP arg_nx, SEXP arg_sx, SEXP arg_norm, SEXP arg_log)
{

	arg_mutost f_arg;
	f_arg.nx		= asReal(arg_nx);
	f_arg.df		= f_arg.nx-1;
	f_arg.sx		= asReal(arg_sx);
	f_arg.ssn		= f_arg.sx/sqrt(f_arg.nx);
	f_arg.norm		= asReal(arg_norm);
	f_arg.give_log	= asInteger(arg_log);
	double lower= asReal(arg_lower),upper= asReal(arg_upper), abs_tol=asReal(arg_abs_tol),rel_tol=asReal(arg_rel_tol) ;
	double abserr;
	double *xans=NULL;
	int neval,res;
	SEXP ans;

	PROTECT(ans=  allocVector(REALSXP,1));
	xans	= REAL(ans);
	res		= nabc_integration_qng(abcMuTOST_sulkl_scalar, &f_arg, lower, upper, abs_tol, rel_tol, xans, &abserr, &neval);
    //std::cout<<"abc_mutost_integrate_sulkl\nnx\t"<<f_arg.nx<<"\nsx\t"<<f_arg.sx<<"\ndf\t"<<f_arg.df<<"\nssn\t"<<f_arg.ssn<<"\nnorm\t"<<f_arg.norm<<"\ngive_log\t"<<f_arg.give_log<<"\nreturn\t"<<res<<"\nans\t"<<*xans<<std::endl;
	UNPROTECT(1);
	return ans;
}

SEXP abc_mutost_integrate_power(SEXP arg_lower, SEXP arg_upper, SEXP arg_abs_tol, SEXP arg_rel_tol, SEXP arg_df, SEXP arg_sT, SEXP arg_c_up, SEXP arg_norm, SEXP arg_log)
{
	double lower= asReal(arg_lower), upper= asReal(arg_upper), abs_tol=asReal(arg_abs_tol), rel_tol=asReal(arg_rel_tol) ;
	double abserr;
	double *xans=NULL;
	int neval,res;

	arg_mutost f_arg;
	f_arg.df		= asReal(arg_df);
	f_arg.sT		= asReal(arg_sT);
	f_arg.tau_up	= asReal(arg_c_up);
	f_arg.norm		= asReal(arg_norm);
	f_arg.give_log	= asInteger(arg_log);

	SEXP ans;

	PROTECT(ans=  allocVector(REALSXP,1));
	xans	= REAL(ans);
	res		= nabc_integration_qng(abcMuTOST_power_scalar, &f_arg, lower, upper, abs_tol, rel_tol, xans, &abserr, &neval);
    //std::cout<<"abc_mutost_integrate_pow\ndf\t"<<f_arg.df<<"\nsT\t"<<f_arg.sT<<"\ntau_up\t"<<f_arg.tau_up<<"\nalpha\t"<<f_arg.alpha<<"\nnorm\t"<<f_arg.norm<<"\ngive_log\t"<<f_arg.give_log<<"\nreturn\t"<<res<<"\nans\t"<<*xans<<std::endl;
	UNPROTECT(1);
	return ans;
}

SEXP abc_mutost_integrate_pow(SEXP arg_lower, SEXP arg_upper, SEXP arg_abs_tol, SEXP arg_rel_tol, SEXP arg_df, SEXP arg_sT, SEXP arg_tau_up, SEXP arg_alpha, SEXP arg_norm, SEXP arg_log)
{

	arg_mutost f_arg;
	f_arg.df		= asReal(arg_df);
	f_arg.sT		= asReal(arg_sT);
	f_arg.tau_up	= asReal(arg_tau_up);
	f_arg.alpha		= asReal(arg_alpha);
	f_arg.norm		= asReal(arg_norm);
	f_arg.give_log	= asInteger(arg_log);

	double lower= asReal(arg_lower), upper= asReal(arg_upper), abs_tol=asReal(arg_abs_tol), rel_tol=asReal(arg_rel_tol) ;
	double abserr;
	double *xans=NULL;
	int neval,res;

	SEXP ans;

	PROTECT(ans=  allocVector(REALSXP,1));
	xans	= REAL(ans);
	res		= nabc_integration_qng(abcMuTOST_pow_scalar, &f_arg, lower, upper, abs_tol, rel_tol, xans, &abserr, &neval);
    //std::cout<<"abc_mutost_integrate_pow\ndf\t"<<f_arg.df<<"\nsT\t"<<f_arg.sT<<"\ntau_up\t"<<f_arg.tau_up<<"\nalpha\t"<<f_arg.alpha<<"\nnorm\t"<<f_arg.norm<<"\ngive_log\t"<<f_arg.give_log<<"\nreturn\t"<<res<<"\nans\t"<<*xans<<std::endl;
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

static inline void abc_mutost_calibrate_tauup_for_mxpw(	const double &mxpw, const double &df, const double &sT, const double &tau_ub, const double &alpha, const double &rho_eq, const double &tol,
	double &maxit, double &tau_u, double &curr_mxpw, double &pw_error)
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
	if(n<0)	
		std::cout<<"abc_mutost_calibrate_tauup_for_mxpw: "<<tau_ub<<'\t'<<tau_ubd<<'\t'<<curr_mxpw<<'\t'<<mxpw<<'\t'<<df<<'\t'<<sT<<'\t'<<maxit<<std::endl;
	
	ERROR_ON(n<0,"abc_mutost_calibrate_tauup_for_mxpw: error at 1a");
	
	for(	pw_error=1, n= CAST(int,maxit);
		n-- && (ABS(pw_error)>tol) && std::floor(tau_ubd*DIGITS)!=std::floor(tau_lbd*DIGITS);
		)
	{
		tau_u= (tau_lbd+tau_ubd)/2;
		abcMuTOST_pow(NRHO, rho, df, tau_u, sT, alpha, &curr_mxpw);
		pw_error= curr_mxpw-mxpw;
		if(pw_error<0)
			tau_lbd= tau_u;
		else
			tau_ubd= tau_u;
		//std::cout<<"H2 "<<curr_mxpw<<'\t'<<mxpw<<'\t'<<tau_ubd<<'\t'<<tau_u<<'\t'<<n<<std::endl;
	}
	if(n<0)
		oprintff("\nabc_mutost_calibrate_tauup_for_mxpw: reached max it %i current pw %g requ pw %g",n,curr_mxpw,mxpw); //could not find Rvprintf ??
	maxit= n+1;
	DELETE(rho);
}

SEXP abc_mutost_calibrate_tauup_for_mxpw(SEXP args)
{
	ERROR_ON(!Rf_isReal(args) ,"abcMuTOST_tau_taulowup_pw: error at 1a ");
	double mxpw=0, df=0, sT= 1, tu_ub= 1, alpha= 0.01, rho_eq=0, tol= 1e-10, maxit= 100;
	double *xans= NULL;
	SEXP ans;

	//convert SEXP into C
	ERROR_ON(length(args)!=8,					"abcMuTOST_tau_taulowup_pw: error at 1b ");
	ERROR_ON((mxpw= REAL(args)[0])>1 || mxpw<=0,	"abcMuTOST_tau_taulowup_pw: error at 1c ");
	ERROR_ON((df= REAL(args)[1])<2,				"abcMuTOST_tau_taulowup_pw: error at 1d ");
	ERROR_ON((sT= REAL(args)[2])<=0,				"abcMuTOST_tau_taulowup_pw: error at 1e ");
	ERROR_ON((tu_ub= REAL(args)[3])<0,			"abcMuTOST_tau_taulowup_pw: error at 1f ");
	ERROR_ON((alpha= REAL(args)[4])>1 || alpha<0,"abcMuTOST_tau_taulowup_pw: error at 1g ");
	rho_eq= REAL(args)[5];
	tol= REAL(args)[6];
	maxit= REAL(args)[7];

	PROTECT(ans=  allocVector(REALSXP,5));
	xans= REAL(ans);
	*(xans+4)= maxit;
	abc_mutost_calibrate_tauup_for_mxpw(mxpw, df, sT, tu_ub, alpha, rho_eq, tol, *(xans+4) /*maxit*/, *(xans+1)/*tau_u*/, *(xans+2)/*curr_mxpw*/, *(xans+3)/*error*/);
	*xans= -*(xans+1);

	UNPROTECT(1);
	return ans;
}

static void abc_mutost_get_KL(const double &nx, const double &sx, const double &ny, const double &sy, const double &mx_pw, const double &alpha, const int &calibrate_tau_up, double &tau_up,const double &pow_scale, double &curr_mxpw, double &KL_div)
{
	int neval;
	double lower, upper, abs_tol, rel_tol, abserr;
	double *cali_tau= NULL;
	arg_mutost pow_arg;
	arg_mutost sulkl_arg;
	kl_integrand_arg KL_arg;
    //calibrate tau_up
	if(calibrate_tau_up)
	{
		cali_tau			= NEW_ARY(double,5);
        *(cali_tau+1)		= tau_up;//tau_ub
        *(cali_tau+2)		= 0;//rho_eq
        *(cali_tau+3)		= 1e-5;//tol
        *(cali_tau+4)		= 100;//maxit
        //std::cout<<"abc_mutost_calibrate_tauup_for_mxpw at a1:\t"<<tau_up<<'\t'<<curr_mxpw<<'\t'<<*cali_tau<<std::endl;
        abc_mutost_calibrate_tauup_for_mxpw(mx_pw, ny-1, sy/std::sqrt(ny), *(cali_tau+1), alpha, *(cali_tau+2), *(cali_tau+3), *(cali_tau+4), tau_up/*tau_up*/, curr_mxpw/*curr_mxpw*/, *cali_tau/*error*/);
        //std::cout<<"abc_mutost_calibrate_tauup_for_mxpw at a2:\t"<<tau_up<<'\t'<<curr_mxpw<<'\t'<<*cali_tau<<std::endl;
        DELETE(cali_tau);
    }
    //compute pow_norm
    upper	= tau_up*pow_scale;
    lower	= -upper;
    rel_tol	= nabcGlobals::NABC_DBL_TOL;
    abs_tol	= rel_tol;
    //create arg for pow
    pow_arg.df		= ny-1;
    pow_arg.sT		= sy/std::sqrt(ny);
    pow_arg.tau_up	= tau_up;
    pow_arg.alpha	= alpha;
    pow_arg.norm	= 1;
    pow_arg.give_log= 0;
    //compute actual max power if tau.u is not calibrated
    if (!calibrate_tau_up)
    {
    	curr_mxpw 	= abcMuTOST_pow_scalar(0,&pow_arg);
    }
    //std::cout<<"abc_mutost_calibrate_tauup_for_mxpw at b1:\t"<<lower<<'\t'<<pow_arg.norm<<'\t'<<pow_arg.give_log<<'\t'<<neval<<std::endl;
    //printBArg(&pow_arg);
    nabc_integration_qng(abcMuTOST_pow_scalar, &pow_arg, lower, upper, abs_tol, rel_tol, &pow_arg.norm /*updated norm*/, &abserr, &neval);
    //std::cout<<"abc_mutost_calibrate_tauup_for_mxpw at b2:\t"<<lower<<'\t'<<pow_arg.norm<<'\t'<<abserr<<'\t'<<neval<<std::endl;
    //create arg for sulkl
    sulkl_arg.nx		= nx;
    sulkl_arg.df		= nx-1;
    sulkl_arg.sx		= sx;
    sulkl_arg.ssn		= sx/std::sqrt(nx);
    sulkl_arg.norm		= pt(upper/sulkl_arg.ssn,sulkl_arg.df,1,0)-pt(lower/sulkl_arg.ssn,sulkl_arg.df,1,0);
    sulkl_arg.give_log	= 0;
    
    //put give_log to 1
    sulkl_arg.give_log	= 1;
    pow_arg.give_log	= 1;
    KL_arg.p			= &abcMuTOST_sulkl_scalar;
    KL_arg.q			= &abcMuTOST_pow_scalar;
    KL_arg.p_arg		= &sulkl_arg;
    KL_arg.q_arg		= &pow_arg;
    //    printf("pow arguments\n");
    //    printBArg(&pow_arg);
    //    printf("\nsulkl arguments\n");
    //    printBArg(&sulkl_arg);
    //    printf("Start KL integration\n");
    nabc_integration_qng(abcKL_integrand, &KL_arg, lower, upper, abs_tol, rel_tol, &KL_div, &abserr, &neval);
    // std::cout<<"abc_mutost_get_KL at c2:\t"<<KL_div<<'\t'<<abserr<<'\t'<<neval<<std::endl;
}

static void abc_mutost_get_KL(void *arg_void)
{
	arg_mutost *arg= (arg_mutost *) arg_void;
    //std::cout<<"abc_mutost_get_KL 1a:\t"<<arg->nx<<'\t'<<arg->sx<<'\t'<<arg->ny<<'\t'<<arg->sy<<'\t'<<arg->mx_pw<<'\t'<<arg->alpha<<'\t'<<arg->calibrate_tau_up<<'\t'<<arg->tau_up<<'\t'<<arg->pow_scale<<'\t'<<arg->curr_mxpw<<'\t'<<arg->KL_div<<std::endl;
	abc_mutost_get_KL(arg->nx, arg->sx, arg->ny, arg->sy, arg->mx_pw, arg->alpha, arg->calibrate_tau_up, arg->tau_up, arg->pow_scale, arg->curr_mxpw, arg->KL_div);
    //std::cout<<"abc_mutost_get_KL 1b:\t"<<arg->nx<<'\t'<<arg->sx<<'\t'<<arg->ny<<'\t'<<arg->sy<<'\t'<<arg->mx_pw<<'\t'<<arg->alpha<<'\t'<<arg->calibrate_tau_up<<'\t'<<arg->tau_up<<'\t'<<arg->pow_scale<<'\t'<<arg->curr_mxpw<<'\t'<<arg->KL_div<<std::endl;
}

SEXP abc_mutost_get_KL(SEXP arg_nx, SEXP arg_sx, SEXP arg_ny, SEXP arg_sy, SEXP arg_mx_pw, SEXP arg_alpha, SEXP arg_calibrate_tau_up, SEXP arg_tau_up, SEXP arg_pow_scale)
{
    //SEXP to C
	double nx			= asReal(arg_nx);
	double sx			= asReal(arg_sx);
	double ny			= asReal(arg_ny);
	double sy			= asReal(arg_sy);
	double mx_pw		= asReal(arg_mx_pw);
	double alpha		= asReal(arg_alpha);
	int calibrate_tau_up= asLogical(arg_calibrate_tau_up);
	double tau_up		= asReal(arg_tau_up);
	double pow_scale	= asReal(arg_pow_scale);
    //answer
	double *xans		= NULL;
	SEXP ans;
	PROTECT(ans= allocVector(REALSXP,3));
	xans	= REAL(ans);
	xans[1]	= tau_up;
	xans[2]	= mx_pw;
    //Call C function
	abc_mutost_get_KL(nx,sx, ny, sy,mx_pw,alpha, calibrate_tau_up, xans[1], pow_scale, xans[2],xans[0]);
	UNPROTECT(1);
	return ans;
}

static double abc_mutost_Brentfun_optimize_tauup_for_KL(double x, void* KL_arg)
{
    //std::cout<<"x:"<<x<<std::endl;
	arg_mutost *arg	= (arg_mutost *) KL_arg;
	arg->tau_up	= x;
    abc_mutost_get_KL(arg);				//compute KL divergence
    return arg->KL_div;
}

SEXP abc_mutost_calibrate_powertighter(SEXP list_KL_args, SEXP arg_max_it)
{
	int max_it	= asInteger(arg_max_it);
	arg_mutost arg;
	double *xans= NULL;
	SEXP ans, ans_names;

	arg.nx 			= asReal(getListElement(list_KL_args,"n.of.x"));
	arg.sx 			= asReal(getListElement(list_KL_args,"s.of.x"));
	arg.ny 			= asReal(getListElement(list_KL_args,"n.of.y"));
	arg.sy 			= asReal(getListElement(list_KL_args,"s.of.y"));
	arg.mx_pw 		= asReal(getListElement(list_KL_args,"mx.pw"));
	arg.alpha 		= asReal(getListElement(list_KL_args,"alpha"));
	arg.tau_up 		= asReal(getListElement(list_KL_args,"tau.u"));
	arg.pow_scale 	= asReal(getListElement(list_KL_args,"pow_scale"));

	abc_generic_calibrate_tauup_for_KL(abc_mutost_get_KL, abc_mutost_Brentfun_optimize_tauup_for_KL, &arg, max_it);

	PROTECT(ans= allocVector(REALSXP,4));
	PROTECT(ans_names= allocVector(STRSXP,4));
	xans	= REAL(ans);
	xans[0]	= arg.KL_div;
	xans[1]	= round(arg.ny);
	xans[2]	= arg.tau_up;
	xans[3]	= arg.curr_mxpw;
	SET_STRING_ELT(ans_names,0,mkChar("KL_div"));
	SET_STRING_ELT(ans_names,1,mkChar("n.of.y"));
	SET_STRING_ELT(ans_names,2,mkChar("tau.u"));
	SET_STRING_ELT(ans_names,3,mkChar("pw.cmx"));
	setAttrib(ans, R_NamesSymbol, ans_names);
	UNPROTECT(2);
	return ans;
}

static double abc_mutost_Brentfun_optimize_yn_for_KL(double x, void* KL_arg)
{
    //std::cout<<"x:"<<x<<std::endl;
	arg_mutost *arg	=(arg_mutost *) KL_arg;
	arg->ny 	= round(x);
    abc_mutost_get_KL(arg);				//compute KL divergence
    return arg->KL_div;
}

SEXP abc_mutost_calibrate_powerbroader(SEXP list_KL_args, SEXP arg_max_it)
{
	int max_it	= asInteger(arg_max_it);
	double *xans= NULL;
	SEXP ans, ans_names;
	arg_mutost arg;

	arg.nx 			= asReal(getListElement(list_KL_args,"n.of.x"));
	arg.sx 			= asReal(getListElement(list_KL_args,"s.of.x"));
	arg.ny 			= asReal(getListElement(list_KL_args,"n.of.y"));
	arg.sy 			= asReal(getListElement(list_KL_args,"s.of.y"));
	arg.mx_pw 		= asReal(getListElement(list_KL_args,"mx.pw"));
	arg.alpha 		= asReal(getListElement(list_KL_args,"alpha"));
	arg.tau_up 		= asReal(getListElement(list_KL_args,"tau.u"));
	arg.pow_scale 	= asReal(getListElement(list_KL_args,"pow_scale"));

	abc_generic_calibrate_yn_for_KL(abc_mutost_get_KL, abc_mutost_Brentfun_optimize_yn_for_KL, &arg, max_it);

	PROTECT(ans= allocVector(REALSXP,4));
	PROTECT(ans_names= allocVector(STRSXP,4));
	xans	= REAL(ans);
	xans[0]	= arg.KL_div;
	xans[1]	= round(arg.ny);
	xans[2]	= arg.tau_up;
	xans[3]	= arg.curr_mxpw;
	SET_STRING_ELT(ans_names,0,mkChar("KL_div"));
	SET_STRING_ELT(ans_names,1,mkChar("n.of.y"));
	SET_STRING_ELT(ans_names,2,mkChar("tau.u"));
	SET_STRING_ELT(ans_names,3,mkChar("pw.cmx"));
	setAttrib(ans, R_NamesSymbol, ans_names);
	UNPROTECT(2);
	return ans;
}

/*
static inline void abcMUTOST_pwvar(	const double &mxpw, const double &df, const double &sT, const double &tau_ub, const double &alpha, const double &rho_eq, const double &tol,
                                   double &maxit, double &curr_pwv)
{
	const int NRHO= 1024;
	const double MEAN= 0.;
	double curr_pw= 1., tau_u=1., error=0., *rho= NULL, *pw=NULL;

	rho= NEW_ARY(double,NRHO);
	pw= NEW_ARY(double,NRHO);

	abc_mutost_calibrate_tauup_for_mxpw(mxpw, df, sT, tau_ub, alpha, rho_eq, tol, maxit, tau_u, curr_pw, error);
	oseq_nout(-tau_u, tau_u, NRHO, rho);
	abcMuTOST_pow(NRHO, rho, df, tau_u, sT, alpha, pw);
	if(oIsZero(NRHO,pw))//must increase nsim to get non-zero power
		curr_pwv= R_PosInf;
	else
		ovar(NRHO, rho, pw, MEAN, curr_pwv);
	DELETE(rho);
	DELETE(pw);
}

SEXP abcMuTOST_pwvar(SEXP args)
{
	ERROR_ON(!Rf_isReal(args) ,"abcMuTOST_pwvar: error at 1a ");
	double mxpw=0, df=0, sT= 1, tu_ub= 1, alpha= 0.01, rho_eq=0, tol= 1e-10, maxit= 100;
	double *xans= NULL;
	SEXP ans;

	//convert SEXP into C
	ERROR_ON(length(args)!=8,					"abcMuTOST_pwvar: error at 1b ");
	ERROR_ON((mxpw= REAL(args)[0])>1 || mxpw<=0,	"abcMuTOST_pwvar: error at 1c ");
	ERROR_ON((df= REAL(args)[1])<2,				"abcMuTOST_pwvar: error at 1d ");
	ERROR_ON((sT= REAL(args)[2])<=0,				"abcMuTOST_pwvar: error at 1e ");
	ERROR_ON((tu_ub= REAL(args)[3])<0,			"abcMuTOST_pwvar: error at 1f ");
	ERROR_ON((alpha= REAL(args)[4])>1 || alpha<0,"abcMuTOST_pwvar: error at 1g ");
	rho_eq= REAL(args)[5];
	tol= REAL(args)[6];
	maxit= REAL(args)[7];

	PROTECT(ans=  allocVector(REALSXP,2));
	xans= REAL(ans);
	*(xans+1)= maxit;
	abcMUTOST_pwvar(mxpw, df, sT, tu_ub, alpha, rho_eq, tol, *(xans+1), *xans);

	UNPROTECT(1);
	return ans;
}
*/

/*
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
		//ERROR_ON(CAST(double,n+1)==maxit && curr_pwv>S2LKL,"abcMuTOST_taulowup_var: variance of power assumed to be smaller than variance of summary likelihood ");
        //std::cout<<"H1 "<<curr_pwv<<'\t'<<S2LKL<<'\t'<<tau_ubd<<'\t'<<n<<std::endl;
		//oprinta(pw, NRHO, std::cout);
	}
	ERROR_ON(n<0,"abcMuTOST_taulowup_var: error at 1a ");
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
		oprintff("abcMuTOST_tau_taulowup_var: reached max it %i current pw variance %g requ pw variance %g",n,curr_pwv,S2LKL);
	maxit= n+1;
	DELETE(rho);
	DELETE(pw);
}
*/
/*
void abcMuTOST_nsim(	const double &nobs, const double &slkl, const double &mxpw, const double &sSim, const double &tau_ub, const double &alpha, const double &rho_eq, const double &tol,
                                  double &maxit, double &nsim, double &tau_u, double &curr_pwv, double &curr_pw, double &error)
{
	const int NRHO= 1024;
	int n= CAST(int,maxit);
	const double S2LKL= slkl*slkl, MEAN=0.;
	double nsim_lb=nobs-1, nsim_ub=std::ceil(nobs/2), xtau_ub= tau_ub, *rho= NULL, *pw=NULL, *cali_tau=NULL;

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
		abc_mutost_calibrate_tauup_for_mxpw(mxpw, nsim_ub-1, sSim/std::sqrt(nsim_ub), xtau_ub, alpha, rho_eq, tol, *(cali_tau+1), tau_u, curr_pw, *cali_tau);
		oseq_nout(-tau_u, tau_u, NRHO, rho);
		abcMuTOST_pow(NRHO, rho, nsim_ub-1, tau_u, sSim/std::sqrt(nsim_ub), alpha, pw);
		if(oIsZero(NRHO,pw))//must increase nsim to get non-zero power
			curr_pwv= 2*S2LKL;
		else
			ovar(NRHO, rho, pw, MEAN, curr_pwv);
        //std::cout<<"H1 "<<curr_pwv<<'\t'<<S2LKL<<'\t'<<tau_u<<'\t'<<nsim_ub<<std::endl;
	}
	ERROR_ON(n<0,"abcMuTOST_nsim: error at 1a ");
	for(	error=1, nsim= nsim_ub, n= CAST(int,maxit);
        n-- && (ABS(error)>tol) && (nsim_lb+1)!=nsim_ub;
        )
	{
        nsim= std::floor( (nsim_lb+nsim_ub)/2 );
        *(cali_tau+1)= maxit;
        xtau_ub= tau_u;
        abc_mutost_calibrate_tauup_for_mxpw(mxpw, nsim-1, sSim/std::sqrt(nsim), xtau_ub, alpha, rho_eq, tol, *(cali_tau+1), tau_u, curr_pw, *cali_tau);
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
        std::cout<<"H2 "<<curr_pwv<<'\t'<<S2LKL<<'\t'<<tau_u<<'\t'<<nsim<<'\t'<<curr_pw<<'\t'<<nsim_lb<<'\t'<<nsim_ub<<std::endl;
	}
	if(n<0)
		oprintff("abcMuTOST_nsim: reached max it %i current pw variance %g requ pw variance %g",n,curr_pwv,S2LKL);

	maxit= n+1;
	DELETE(rho);
	DELETE(pw);
	DELETE(cali_tau);
}
*/
