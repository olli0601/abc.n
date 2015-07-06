//
//  nabc_mutost.cpp
//  nABC
//
//
//
#include "nabc_utils.h"
#include "nabc_ftest.h"
#include "nabc_Rcall.h"
#include "nabc_integrate.h"
#include "nabc_KLdiv.h"

void nabcFTEST_printArg(arg_ftest *arg)
{
	printf("nx=%f\tt2x=%f\tp=%f\nny=%f\ttau=%f\n norm=%f\tgive_log=%d\tmx_pw=%f\talpha=%f\n",arg->nx,arg->t2x,arg->p,arg->ny,arg->tau,arg->norm,arg->give_log,arg->mx_pw,arg->alpha);
}

static inline void abcFTEST_pow(const int &nrho, double * const rho, const double &tau, const double &ny, const double &p, const double &alpha, double *ans)
{
	int n= nrho;
	const int LOG= 0, LOWERTAIL=1;
	double *xans= ans, *xrho=rho;
	const double CRIT= qnf(alpha, p, ny-p, ny*tau, LOWERTAIL, LOG);

	for(; n--; xrho++, xans++)
	{
		*xans= pnf(CRIT, p, ny-p, *xrho*ny, LOWERTAIL, LOG);
		*xans= *xans<0 ? 0 : *xans;
	}
}

static inline double abcFTEST_pow_scalar(double x, void *arg)
{
	arg_ftest *a=(arg_ftest *) arg;
	double ans= pnf(qnf(a->alpha, a->p, a->ny-a->p, a->ny*a->tau, 1, 0), a->p, a->ny-a->p, x*a->ny, 1, 0);
    /* Can be 0 at some point! */
    /* Usually this happens due to numerical inacurracy in the tail. */
    /* To avoid infinity we assume that ans=nabcGlobals::NABC_DBL_MIN at these points. */
	ans= ans<=0 ? nabcGlobals::NABC_DBL_MIN : ans/a->norm;
	return a->give_log ? log(ans) : ans;
}

SEXP abcFTEST_pow(SEXP arg_rho, SEXP arg_tau, SEXP arg_ny, SEXP arg_p, SEXP arg_alpha)
{
	ERROR_ON(!Rf_isReal(arg_rho) ,"abcFTEST_pow: error at 1a ");
	int nrho= Rf_length(arg_rho);
	double *xrho= REAL(arg_rho), *xans=NULL;
	double tau= ::Rf_asReal(arg_tau), ny= ::Rf_asReal(arg_ny), p= ::Rf_asReal(arg_p), alpha= ::Rf_asReal(arg_alpha);
	SEXP ans;
	PROTECT(ans=  allocVector(REALSXP,nrho));
	xans= REAL(ans);
	//std::cout<<"abcMuTOST_pow_c\t"<<nrho<<'\t'<<*xrho<<'\t'<<df<<'\t'<<tau_up<<'\t'<<sT<<'\t'<<alpha<<std::endl;
	abcFTEST_pow(nrho, xrho, tau, ny, p, alpha, xans);
	UNPROTECT(1);
	return ans;
}

static inline void abcFTEST_sulkl(const int &nrho, double * const rho, const double &t2x, const double &nx, const double &p, const double &norm, const int &give_log, double *ans)
{
	int n= nrho;
	const double XV= t2x*(nx-p)/(p*(nx-1));
	double *xans= ans, *xrho=rho;
	for(; n--; xrho++, xans++)
	{
		*xans= dnf(XV, p, nx-p, *xrho*nx, give_log);
		*xans= give_log ? *xans-log(norm) : *xans/norm;
	}
}

static inline double abcFTEST_sulkl_scalar(double x, void *arg_void)
{
	arg_ftest *a=(arg_ftest *) arg_void;
	double ans;
    //std::cout<<"abcMuTOST_sulkl_scalar input:\nssn\t"<<arg->ssn<<"\ndf\t"<<arg->df<<"\ngive_log\t"<<arg->give_log<<"\nnorm\t"<<arg->norm<<"\nx\t"<<x<<std::endl;
	ans= dnf(	a->t2x*(a->nx-a->p)/(a->p*(a->nx-1)),	a->p,	a->nx-a->p,		a->nx*x, 	0);
	ans= a->give_log ? ans-log(a->norm) : ans/a->norm;
	return ans;
}

SEXP abcFTEST_sulkl(SEXP arg_rho, SEXP arg_tx, SEXP arg_nx, SEXP arg_p, SEXP arg_norm, SEXP arg_log)
{
	int nrho= Rf_length(arg_rho), give_log=asInteger(arg_log);
	double *xrho= REAL(arg_rho), *xans=NULL;
	double nx= asReal(arg_nx), tx= asReal(arg_tx), p= asReal(arg_p), norm= asReal(arg_norm);
	SEXP ans;
	PROTECT(ans=  allocVector(REALSXP,nrho));
	xans= REAL(ans);
	//std::cout<<"abcMuTOST_sulkl_c\t"<<nrho<<'\t'<<*xrho<<'\t'<<nx<<'\t'<<sx<<'\t'<<norm<<'\t'<<log<<std::endl;
	abcFTEST_sulkl(nrho, xrho, tx, nx, p, norm, give_log, xans);
	UNPROTECT(1);
	return ans;
}

static inline double abcFTEST_criticalvalue(	const double &tau, const double &ny, const double &p, const double &alpha )
{
	return p * (ny-1) * qnf(alpha, p, ny-p, ny*tau, 1, 0) / (ny-p);
}

static inline void abcFTEST_calibrate_tau_for_mxpw(	const double &mxpw, const double &ny, const double &p, const double &tau_ub, const double &alpha, const double &tol,
	double &maxit, double &tau, double &curr_mxpw, double &pw_error)
{
	ERROR_ON(mxpw<=0.5,"abcFTEST_calibrate_tau_for_mxpw: mx.pw<0.5");
	ERROR_ON(mxpw>=0.99,"abcFTEST_calibrate_tau_for_mxpw: mx.pw>0.99");
	ERROR_ON(alpha<=0,"abcFTEST_calibrate_tau_for_mxpw: alpha<=0");
	ERROR_ON(alpha>=0.25,"abcFTEST_calibrate_tau_for_mxpw: alpha>=0.25");
	ERROR_ON(tau_ub<=0,"abcFTEST_calibrate_tau_for_mxpw: tau_ub<=0");
	ERROR_ON(maxit<=10,"abcFTEST_calibrate_tau_for_mxpw: maxit<=10");
	ERROR_ON(tol>=0.01,"abcFTEST_calibrate_tau_for_mxpw: tol>=0.01");
	ERROR_ON(ny<2,"abcFTEST_calibrate_tau_for_mxpw: ny<2");
	ERROR_ON(p<2,"abcFTEST_calibrate_tau_for_mxpw: p<2");
	ERROR_ON(p>ny,"abcFTEST_calibrate_tau_for_mxpw: p>ny");

	const double PWL= mxpw-tol, PWU= mxpw+tol, DIGITS= std::ldexp(1,30);
	double	tau_u=0, tau_l=0;
	double	tmp1=0, tmp2=0;
	int n= CAST(int, maxit);
	arg_ftest arg;
	arg.p		= p;
	arg.ny		= ny;
	arg.tau		= tau_ub;
	arg.norm	= 1;
	arg.give_log= 0;
	arg.alpha	= alpha;

	// find upper tolerance bound
	for(; n-- && abcFTEST_pow_scalar(0., &arg)<PWU ; arg.tau++);
	ERROR_ON(n<0,"abcFTEST_calibrate_tau_for_mxpw: could not find upper tau bound");
	tau_u= arg.tau;
	// find lower tolerance bound
	for(	n= CAST(int, maxit), arg.tau=(tau_l+tau_u)/2;
			n-- && abcFTEST_pow_scalar(0., &arg)>PWL ;
			arg.tau=(arg.tau+tau_l)/2	);
	ERROR_ON(n<0,"abcFTEST_calibrate_tau_for_mxpw: could not find lower tau bound");
	tau_l= arg.tau;
	ERROR_ON(tau_l>=tau_u, "abcFTEST_calibrate_tau_for_mxpw: found tau_l>=tau_u");
	// binary search to find tau
	for(	n= CAST(int,maxit), pw_error=1;
			n-- && (ABS(pw_error)>tol) && std::floor(tau_u*DIGITS)!=std::floor(tau_l*DIGITS);
			)
		{
			arg.tau		= (tau_l+tau_u)/2;
			curr_mxpw	= abcFTEST_pow_scalar(0., &arg);
			pw_error	= curr_mxpw-mxpw;
			if(pw_error<0)
				tau_l= arg.tau;
			else
				tau_u= arg.tau;
		}
	ERROR_ON(n<0,"abcFTEST_calibrate_tau_for_mxpw: could not find tau in binary search");
	maxit	= n+1;
	tau		= arg.tau;
}

SEXP abcFTEST_calibrate_tau_for_mxpw(SEXP args)
{
	double mxpw=0, ny=0, p= 0, tau_ub= 0, alpha= 0, tol= 0, maxit= 0, tau=0, curr_mxpw=0, pw_error=0, crit=0;
	double *xans= NULL;
	SEXP ans, ans_names;

	//(mx.pw, ny, p, tau.ub, alpha = 0.01, tol = 1e-5, max.it = 100, verbose = 0)
	mxpw 	= asReal(getListElement(args,"mx.pw"));
	ny 		= asReal(getListElement(args,"ny"));
	p 		= asReal(getListElement(args,"p"));
	tau_ub 	= asReal(getListElement(args,"tau.ub"));
	alpha 	= asReal(getListElement(args,"alpha"));
	tol 	= asReal(getListElement(args,"tol"));
	maxit 	= asReal(getListElement(args,"max.it"));
	abcFTEST_calibrate_tau_for_mxpw(mxpw, ny, p, tau_ub, alpha, tol, maxit, tau, curr_mxpw, pw_error);
	crit	= abcFTEST_criticalvalue(tau, ny, p, alpha );

	//convert SEXP into C
	PROTECT(ans=  allocVector(REALSXP,5));
	PROTECT(ans_names= allocVector(STRSXP,5));
	xans	= REAL(ans);
	xans[0]	= crit;
	xans[1]	= tau;
	xans[2]	= curr_mxpw;
	xans[3]	= pw_error;
	xans[4]	= maxit;
	SET_STRING_ELT(ans_names,0,mkChar("c"));
	SET_STRING_ELT(ans_names,1,mkChar("tau"));
	SET_STRING_ELT(ans_names,2,mkChar("curr.pw"));
	SET_STRING_ELT(ans_names,3,mkChar("error.pw"));
	SET_STRING_ELT(ans_names,4,mkChar("max.it"));
	setAttrib(ans, R_NamesSymbol, ans_names);
	UNPROTECT(2);
	return ans;
}

static inline void abcFTEST_get_KL(const double &nx, const double &t2x, const double &ny, const double &p, double &tau, const double &mx_pw, const double &alpha,  const double &pow_scale, double &curr_mxpw, double &pw_error, double &KL_div)
{
	int neval=0;
	double lower=0, upper=0, abs_tol=0, rel_tol=0, abserr=0, tmp1=0;
	double *cali_tau= NULL;
	arg_ftest pow_arg;
	arg_ftest sulkl_arg;
	kl_integrand_arg KL_arg;

	//calibrate tau_up
	cali_tau			= NEW_ARY(double, 3);
	*cali_tau			= 4*t2x;
    *(cali_tau+1)		= 1e-5;//tol
    *(cali_tau+2)		= 100;//maxit
    //std::cout<<"abc_mutost_calibrate_tauup_for_mxpw at a1:\t"<<tau_up<<'\t'<<curr_mxpw<<'\t'<<*cali_tau<<std::endl;
    abcFTEST_calibrate_tau_for_mxpw(mx_pw, ny, p, *cali_tau /*tau_ub*/, alpha, *(cali_tau+1) /*tol*/, *(cali_tau+2) /*maxit*/, tau, curr_mxpw, pw_error);
    //std::cout<<"abc_mutost_calibrate_tauup_for_mxpw at a2:\t"<<tau_up<<'\t'<<curr_mxpw<<'\t'<<*cali_tau<<std::endl;
    DELETE(cali_tau);

    //compute pow_norm
    upper	= tau*pow_scale;
    lower	= 0;
    rel_tol	= nabcGlobals::NABC_DBL_TOL;
    abs_tol	= rel_tol;
    //create arg for pow
    pow_arg.p			= p;
    pow_arg.ny			= ny;
    pow_arg.tau			= tau;
    pow_arg.norm		= 1;
    pow_arg.give_log	= 0;
    pow_arg.alpha		= alpha;
    //std::cout<<"abc_mutost_calibrate_tauup_for_mxpw at b1:\t"<<lower<<'\t'<<pow_arg.norm<<'\t'<<pow_arg.give_log<<'\t'<<neval<<std::endl;
    //printBArg(&pow_arg);
    nabc_integration_qng(abcFTEST_pow_scalar, &pow_arg, lower, upper, abs_tol, rel_tol, &pow_arg.norm /*updated norm*/, &abserr, &neval);
    //std::cout<<"abc_mutost_calibrate_tauup_for_mxpw at b2:\t"<<lower<<'\t'<<pow_arg.norm<<'\t'<<abserr<<'\t'<<neval<<std::endl;
    //create arg for sulkl
    sulkl_arg.p			= p;
    sulkl_arg.nx		= nx;
    sulkl_arg.t2x		= t2x;
    sulkl_arg.norm		= 1;
    sulkl_arg.give_log	= 0;
    nabc_integration_qng(abcFTEST_sulkl_scalar, &sulkl_arg, lower, upper, abs_tol, rel_tol, &sulkl_arg.norm /*updated norm*/, &abserr, &neval);

    //put give_log to 1
    sulkl_arg.give_log	= 1;
    pow_arg.give_log	= 1;
    KL_arg.p			= &abcFTEST_sulkl_scalar;
    KL_arg.q			= &abcFTEST_pow_scalar;
    KL_arg.p_arg		= &sulkl_arg;
    KL_arg.q_arg		= &pow_arg;
    if(1)
    {
    	printf("pow arguments\n");
    	nabcFTEST_printArg(&pow_arg);
    	printf("\nsulkl arguments\n");
    	nabcFTEST_printArg(&sulkl_arg);
    	printf("Start KL integration\n");
    }
    nabc_integration_qng(abcKL_integrand, &KL_arg, lower, upper, abs_tol, rel_tol, &KL_div, &abserr, &neval);
    //std::cout<<"abc_mutost_calibrate_tauup_for_mxpw at c2:\t"<<KL_div<<'\t'<<abserr<<'\t'<<neval<<std::endl;
}

static void abcFTEST_get_KL(void *arg_void)
{
	arg_ftest *arg= (arg_ftest *) arg_void;
    //std::cout<<"abc_mutost_get_KL 1a:\t"<<arg->nx<<'\t'<<arg->sx<<'\t'<<arg->ny<<'\t'<<arg->sy<<'\t'<<arg->mx_pw<<'\t'<<arg->alpha<<'\t'<<arg->calibrate_tau_up<<'\t'<<arg->tau_up<<'\t'<<arg->pow_scale<<'\t'<<arg->curr_mxpw<<'\t'<<arg->KL_div<<std::endl;
	abcFTEST_get_KL(arg->nx, arg->t2x, arg->ny, arg->p, arg->tau, arg->mx_pw, arg->alpha, arg->pow_scale, arg->curr_mxpw, arg->pw_error, arg->KL_div);
    //std::cout<<"abc_mutost_get_KL 1b:\t"<<arg->nx<<'\t'<<arg->sx<<'\t'<<arg->ny<<'\t'<<arg->sy<<'\t'<<arg->mx_pw<<'\t'<<arg->alpha<<'\t'<<arg->calibrate_tau_up<<'\t'<<arg->tau_up<<'\t'<<arg->pow_scale<<'\t'<<arg->curr_mxpw<<'\t'<<arg->KL_div<<std::endl;
}

static double abcFTEST_Brentfun_optimize_yn_for_KL(double x, void* KL_arg)
{
    //std::cout<<"x:"<<x<<std::endl;
	arg_ftest *arg	= (arg_ftest *) KL_arg;
	arg->ny 		= round(x);
	abcFTEST_get_KL(arg->nx, arg->t2x, arg->ny, arg->p, arg->tau, arg->mx_pw, arg->alpha, arg->pow_scale, arg->curr_mxpw, arg->pw_error, arg->KL_div);
    return arg->KL_div;
}

/*
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
*/
/*
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
*/

static inline void abcFTEST_calibrate_yn_for_KL(void (*KL_divergence)(void*), double (*KL_optimize)(double, void*), void *KL_arg_void, const int &max_it)
{

    double test_KL_div, current_KL_div, ny_lb, ny_ub;
    int curr_it;
    arg_ftest *KL_arg= (arg_ftest *)  KL_arg_void;

    // test KL
    KL_arg->ny++;
    KL_divergence( KL_arg );
    test_KL_div= KL_arg->KL_div;
std::cout<<"test KL at ny+1:"<<test_KL_div<<"\tKL_arg->ny:"<<KL_arg->ny<<"\tKL_arg->tau:"<<KL_arg->tau<<std::endl;

    // current KL
    KL_arg->ny--;
    KL_divergence( KL_arg );
    current_KL_div = KL_arg->KL_div;
    //current_KL_div = KL_divergence_switch_arg(KL_arg->ny+1,&KL_switch_arg);
std::cout<<"current KL:"<<current_KL_div<<"\tKL_arg->ny:"<<KL_arg->ny<<"\tKL_arg->tau:"<<KL_arg->tau<<std::endl;
	// if KL is increasing for ny+1, then find optimum value between [1,ny]
	// else determine upper bound and then find optimum value between [ny, ny_ub]
	if (test_KL_div > current_KL_div)
    {
        ny_lb = 1;
        ny_ub = KL_arg->ny;
    }
    else
    {
        curr_it 	= max_it;
        KL_arg->ny*= 2;
        KL_divergence( KL_arg );
        test_KL_div	= KL_arg->KL_div;
        while ((test_KL_div < current_KL_div) && curr_it--)
        {
            current_KL_div 	= test_KL_div;
            KL_arg->ny*= 2;
            KL_divergence( KL_arg );
            test_KL_div		= KL_arg->KL_div;
std::cout<<"upper trial KL:"<<test_KL_div<<'\t'<<current_KL_div<<"\tKL_arg->ny:"<<KL_arg->ny<<"\tKL_arg->tau:"<<KL_arg->tau<<std::endl;
			ERROR_ON(KL_arg->ny>1e7 ,"could not find upper bound for n.of.y with ny<1e7.");
            //test_KL_div 	= KL_divergence_switch_arg(2*KL_arg->ny,&KL_switch_arg);
        }
        ERROR_ON(!curr_it ,"could not find upper bound for n.of.y within the maximum number of iterations.");

        //if curr_it==max_it then the lb is obtained by dividing ub by 2, by 4 otherwise
        ny_lb = (curr_it == max_it) ? KL_arg->ny/2 : KL_arg->ny/4;
        ny_ub = KL_arg->ny;
    }
std::cout<<"ny lower upper:"<<ny_lb<<'\t'<<ny_ub<<std::endl;
    //minimize KL between ny_lb and KL_arg->ny
    KL_arg->ny = Brent_fmin(ny_lb, ny_ub, KL_optimize, KL_arg, nabcGlobals::NABC_DBL_TOL, 1);//
    //compute final KL_div
    KL_divergence(KL_arg);
}

SEXP abcFTEST_calibrate_KL(SEXP list_KL_args, SEXP arg_max_it)
{
	int max_it	= asInteger(arg_max_it);
	double *xans= NULL;
	SEXP ans, ans_names;
	arg_ftest arg;
//(t2.x, n.of.x, p, n.of.y=n.of.x, mx.pw=0.9, alpha=0.01, max.it=100, debug=FALSE, plot=FALSE)
	arg.nx 			= asReal(getListElement(list_KL_args,"n.of.x"));
	arg.t2x			= asReal(getListElement(list_KL_args,"t2.x"));
	arg.p 			= asReal(getListElement(list_KL_args,"p"));
	arg.ny 			= asReal(getListElement(list_KL_args,"n.of.y"));
	arg.mx_pw 		= asReal(getListElement(list_KL_args,"mx.pw"));
	arg.alpha 		= asReal(getListElement(list_KL_args,"alpha"));
	arg.pow_scale 	= 1.5;

	abcFTEST_calibrate_yn_for_KL(abcFTEST_get_KL, abcFTEST_Brentfun_optimize_yn_for_KL, &arg, max_it);
	//c(n.of.y=n.of.y, tau=tmp['tau'], c=tmp['c'], pw.cmx=tmp['pw.cmx'], KL_div=tmp['KL_div'])
	PROTECT(ans= allocVector(REALSXP,5));
	PROTECT(ans_names= allocVector(STRSXP,5));
	xans	= REAL(ans);
	xans[0]	= round(arg.ny);
	xans[1]	= arg.tau;
	xans[2]	= abcFTEST_criticalvalue(	arg.tau, round(arg.ny), arg.p, arg.alpha );;
	xans[3]	= arg.curr_mxpw;
	xans[4]	= arg.KL_div;
	SET_STRING_ELT(ans_names,0,mkChar("n.of.y"));
	SET_STRING_ELT(ans_names,1,mkChar("tau"));
	SET_STRING_ELT(ans_names,2,mkChar("c"));
	SET_STRING_ELT(ans_names,3,mkChar("pw.cmx"));
	SET_STRING_ELT(ans_names,4,mkChar("KL_div"));
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
