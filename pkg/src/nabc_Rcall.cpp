#include "nabc_Rcall.h"



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

SEXP abcCalibrate_minimize_KL(SEXP arg_test_name, SEXP arg_calibration_name, SEXP list_KL_args, SEXP arg_max_it)
{
    //SEXP to C
    int max_it=asInteger(arg_max_it);
    
    EqTestValue eq_test_val;
    CalibrationValue calibration_val;
    
    void (*KL_divergence)(kl_arg *);
    
    //answer
    double *xans= NULL;
    SEXP ans, ans_names;
    
    //get element from list_KL_args and put them in a kl_arg structure
    //by default just add the arguments you need, un-needed arguments will just be NULL if not present in list_KL_args
    kl_arg arg;
    arg.nx = asReal(getListElement(list_KL_args,"n.of.x"));
    arg.sx = asReal(getListElement(list_KL_args,"s.of.x"));
    arg.ny = asReal(getListElement(list_KL_args,"n.of.y"));
    arg.sy = asReal(getListElement(list_KL_args,"s.of.y"));
    arg.mx_pw = asReal(getListElement(list_KL_args,"mx.pw"));
    arg.alpha = asReal(getListElement(list_KL_args,"alpha"));
    arg.tau_up = asReal(getListElement(list_KL_args,"tau.u"));
    arg.pow_scale = asReal(getListElement(list_KL_args,"pow_scale"));
    
    strcpy(nabcGlobals::BUFFER,CHAR(STRING_ELT(arg_test_name, 0)));
    
    try{
        eq_test_val=nabcGlobals::s_mapEqTestValue.at(nabcGlobals::BUFFER);
    }
    catch (const std::out_of_range& oor) {
        error("%s is not in the lookup table of equivalence test, see documentation\n",nabcGlobals::BUFFER);
    }
    
    //given the eq_test_val choose right function and ptr to it
    switch (eq_test_val) {
        case MUTOST_ONE_SAMPLE:
            KL_divergence=&abcMuTOST_KL;
            break;
            
        default:
            error("%s:%s:%d: edit switch to include test %s",__FILE__,__FUNCTION__,__LINE__,nabcGlobals::BUFFER);
            break;
    }
    
    
        
    //Call C function according to type_of_calibration
    strcpy(nabcGlobals::BUFFER,CHAR(STRING_ELT(arg_calibration_name, 0)));
    
    try{
        calibration_val=nabcGlobals::s_mapCalibrationValue.at(nabcGlobals::BUFFER);
    }
    catch (const std::out_of_range& oor) {
        error("%s is not in the lookup table of calibration, see documentation\n",nabcGlobals::BUFFER);
    }
    
    //given the calibration_val run appropriate calibration function
    switch (calibration_val) {
        case TAU_NO_MAX_POW:
            abcCalibrate_tau_nomxpw_yesKL(KL_divergence, &arg, max_it);
            
            PROTECT(ans= allocVector(REALSXP,3));
            PROTECT(ans_names= allocVector(STRSXP,3));
            
            xans= REAL(ans);
            
            xans[0]=arg.KL_div;
            xans[1]=arg.tau_up;
            xans[2]=arg.curr_mxpw;

            SET_STRING_ELT(ans_names,0,mkChar("KL_div"));
            SET_STRING_ELT(ans_names,1,mkChar("tau.u"));
            SET_STRING_ELT(ans_names,2,mkChar("pw.cmx"));

            break;
            
        case M_AND_TAU_YES_MAX_POW:
            abcCalibrate_m_and_tau_yesmxpw_yesKL(KL_divergence, &arg, max_it);
            
            PROTECT(ans= allocVector(REALSXP,4));
            PROTECT(ans_names= allocVector(STRSXP,4));

            xans= REAL(ans);
            
            xans[0]=arg.KL_div;
            xans[1]=round(arg.ny);
            xans[2]=arg.tau_up;
            xans[3]=arg.curr_mxpw;
            
            SET_STRING_ELT(ans_names,0,mkChar("KL_div"));
            SET_STRING_ELT(ans_names,1,mkChar("n.of.y"));
            SET_STRING_ELT(ans_names,2,mkChar("tau.u"));
            SET_STRING_ELT(ans_names,3,mkChar("pw.cmx"));

            break;
            
        default:
            error("%s:%s:%d: edit switch to include calibration %s",__FILE__,__FUNCTION__,__LINE__,nabcGlobals::BUFFER);
            break;
    }

    setAttrib(ans, R_NamesSymbol, ans_names);
    
    UNPROTECT(2);
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
