//
//  nabc_utils.cpp
//  nABC
//
//
//

#include "nabc_utils.h"


void oprintf(const char * format, ...)
{
	va_list args;
	va_start(args, format);
	//Rvprintf(format, args);
	vprintf(format, args);
	va_end(args);
}

void oprintff(const char * format, ...)
{
	va_list args;
	va_start(args, format);
	//Rvprintf(format, args);
	vprintf(format, args);
	va_end(args);
	fflush(stdout);
}

void oprinta(double const * const start,const int &n,std::ostream& os)
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

void oseq_nout(const double &a, const double &b, const int &n, double * const ans)
{
	FAIL_ON(n<2,"oseq_nout: error at 1a %c");
	double *xans=NULL, *yans=NULL;
	int xn= n-1;
	const double by= (b-a)/xn;
    
	yans= xans= ans;
	*xans++= a;
	for(; 		xn--; 		*xans++= *yans++ + by);
}

void ovar(const int &n, double * const x, double * const fx, const double &mean, double &var)
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


int oIsZero(const int &n, double * const x)
{
	int xn= n, isZero=1;
	double *xx= x;
	for(; 		isZero && xn--;		xx++)
		isZero= CAST(int,*xx==0);
    //isZero= (	*xx > std::numeric_limits<double>::epsilon() || *xx < -std::numeric_limits<double>::epsilon()	)?0:1;
    //isZero= CAST(int,*xx==0);
	return isZero;
}



void printBArg(basic_arg *arg)
{
    printf("nx=%f\nsx=%f\nssn=%f\ndf=%f\nnorm=%f\ngive_log=%d\nalpha=%f\ntau_up=%f\nsT=%f\n",arg->nx,arg->sx,arg->ssn,arg->df,arg->norm,arg->give_log,arg->alpha,arg->tau_up,arg->sT);
    
}







