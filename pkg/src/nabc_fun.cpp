#include <cmath>
#include <Rmath.h>
#include "nabc_error_handling.h"
#include "nabc_fun.h"

SEXP abcScaledChiSq(SEXP args)
{
	FAIL_ON(! Rf_isReal(args) ,"abcScaledChiSqu: error at 1a %c");
	const int LOWERTAIL= 1;
	int error= 0;
	double scale=0, df=0, tl= 1, tu= 1, tuisc= 1, scitu= 1, scitl= 1, alpha= 0.01, tol= 1e-10, maxit= 100, incit= 0.05;
	double c1, c2, c1l, c1r, level;
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
	tuisc= tu/scale;
	scitu= scale/tu;
	scitl= scale/tl;
//std::cout<<"abcScaledChiSq input:\ndf\t"<<df<<"\ntl\t"<<tl<<"\ntu\t"<<tu<<"\nalpha\t"<<alpha<<"\ntol\t"<<tol<<"\nmaxit\t"<<maxit<<"\nincit\t"<<incit<<std::endl;

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
	}
//std::cout<<c1<<'\t'<<c2<<'\t'<<c1r<<'\t'<<c1l<<'\t'<<ABSDIFF(level,alpha)<<'\t'<<tol<<std::endl;	std::cout<<"OK2"<<std::endl;
	//output
	PROTECT(ans=  allocVector(REALSXP,5));
	xans= REAL(ans);
	*xans++= c1;
	*xans++= c2;
	*xans++= pchisq(c2*scale,df,LOWERTAIL,0) - pchisq(c1*scale,df,LOWERTAIL,0);//power at rho=1
	*xans++= ABSDIFF(level,alpha);
	*xans=	 maxit;

	UNPROTECT(1);
	return ans;
}


