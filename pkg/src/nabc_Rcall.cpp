#include "nabc_Rcall.h"

/*
SEXP abcMuTOST_taulowup_var(SEXP args)
{
	ERROR_ON(!Rf_isReal(args) ,"abcMuTOST_tau_taulowup_var: error at 1a ");
	double sLkl=0, df=0, sT= 1, tu_ub= 1, alpha= 0.01, rho_eq=0, tol= 1e-10, maxit= 100;
	double *xans= NULL;
	SEXP ans;
    
	//convert SEXP into C
	ERROR_ON(length(args)!=8,					"abcMuTOST_tau_taulowup_var: error at 1b ");
	ERROR_ON((sLkl= REAL(args)[0])<=0,			"abcMuTOST_tau_taulowup_var: error at 1c ");
	ERROR_ON((df= REAL(args)[1])<2,				"abcMuTOST_tau_taulowup_var: error at 1d ");
	ERROR_ON((sT= REAL(args)[2])<=0,				"abcMuTOST_tau_taulowup_var: error at 1e ");
	ERROR_ON((tu_ub= REAL(args)[3])<0,			"abcMuTOST_tau_taulowup_var: error at 1f ");
	ERROR_ON((alpha= REAL(args)[4])>1 || alpha<0,"abcMuTOST_tau_taulowup_var: error at 1g ");
	rho_eq= REAL(args)[5];
	tol= REAL(args)[6];
	maxit= REAL(args)[7];
    
	PROTECT(ans=  allocVector(REALSXP,5));
	xans= REAL(ans);
	*(xans+4)= maxit;
	abcMuTOST_taulowup_var(sLkl, df, sT, tu_ub, alpha, rho_eq, tol, *(xans+4), *(xans+1), *(xans+2), *(xans+3));
	*xans= -*(xans+1);
    
	UNPROTECT(1);
	return ans;
}
*/
/*
SEXP abcMuTOST_nsim(SEXP args)
{
	ERROR_ON(!Rf_isReal(args) ,"abcMuTOST_nsim: error at 1a ");
	double nobs= 1, sLkl=0, mxpw=0, sSim=0, tu_ub= 1, alpha= 0.01, rho_eq=0, tol= 1e-10, maxit= 100;
	double *xans= NULL;
	SEXP ans;
    
	//convert SEXP into C
	ERROR_ON(length(args)!=9,					"abcMuTOST_nsim: error at 1b ");
	ERROR_ON((nobs= REAL(args)[0])<=0,			"abcMuTOST_nsim: error at 1c ");
	ERROR_ON((sLkl= REAL(args)[1])<=0,			"abcMuTOST_nsim: error at 1d ");
	ERROR_ON((mxpw= REAL(args)[2])<=0 || mxpw>1,	"abcMuTOST_nsim: error at 1e ");
	ERROR_ON((sSim= REAL(args)[3])<=0,			"abcMuTOST_nsim: error at 1f ");
	ERROR_ON((tu_ub= REAL(args)[4])<0,			"abcMuTOST_nsim: error at 1f ");
	ERROR_ON((alpha= REAL(args)[5])>1 || alpha<0,"abcMuTOST_nsim: error at 1g ");
	rho_eq= REAL(args)[6];
	tol= REAL(args)[7];
	maxit= REAL(args)[8];
    
	PROTECT(ans=  allocVector(REALSXP,7));
	xans= REAL(ans);
	*(xans+6)= maxit;
	abcMuTOST_nsim(nobs, sLkl, mxpw, sSim, tu_ub, alpha, rho_eq, tol, *(xans+6) , *xans ,*(xans+2), *(xans+3), *(xans+4), *(xans+5));
	*(xans+1)= -*(xans+2);
    
	UNPROTECT(1);
	return ans;
}
*/

SEXP abcIntersectLevelSets(SEXP m1, SEXP m2, SEXP s)
{
	ERROR_ON(! Rf_isMatrix(m1) ,"abcIntersectLevelSets: error at 1a ");
	ERROR_ON(! Rf_isMatrix(m2) ,"abcIntersectLevelSets: error at 1b ");
	ERROR_ON(! Rf_isReal(s) ,"abcIntersectLevelSets: error at 1c ");
    
	int i,j1,j2;
	int nr1= Rf_nrows(m1), nr2= Rf_nrows(m2), nc1= Rf_ncols(m1), nc2= Rf_ncols(m2), ns= Rf_length(s);
	double tmp, *xd= NULL, *ym2=NULL, *xm1= NULL, *xm2= NULL, *xs= NULL, *is=NULL;
	SEXP ans;
    
	ERROR_ON(nr1!=nr2,"abcIntersectLevelSets: error at 2a ");
	ERROR_ON(nr1!=ns,"abcIntersectLevelSets: error at 2b ");
    
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
