#' this file contains all R functions of the abc-n package
#' @import nortest locfit
#' @useDynLib nabc

NABC.DEFAULT.ANS<- {tmp<- c(0,50, 1, NA, NA, NA, 0, 0,0,1,1,1,NA); names(tmp)<- c("lkl", "error", "pval","link.mc.obs","link.mc.sim", "rho.mc", "cil", "cir","al","ar","pfam.pval","ts.pval","mx.pow"); tmp}

#------------------------------------------------------------------------------------------------------------------------
#' Test if summary values are normally distributed
#' @param x 			summary values
#' @param norma.test 	name of function with which normality of the summary values is tested
#' @return p value of the test
#' @examples nabc.get.pfam.pval(rnorm(1e4),"shapiro.test")
nabc.get.pfam.pval<- function(x,normal.test)
{
	if(!normal.test%in%c("shapiro.test","lillie.test","cvm.test","ad.test","pearson.test","sf.test"))		
		stop("nabc.get.pfam.pval: error at 1a")
	if(normal.test%in%c("lillie.test","cvm.test","ad.test","pearson.test","sf.test"))		
		require(nortest, lib.loc=LIB.LOC)
	
	ifelse(	any(diff(x)>0), 
			ifelse(normal.test%in%c("shapiro.test","sf.test") && length(x)>5000, 
					eval(call(normal.test,x[1:5000]))$p.value, 
					eval(call(normal.test,x))$p.value), 
			0)
}
#------------------------------------------------------------------------------------------------------------------------
#' Compute power of the asymptotic equivalence test for autocorrelations at lag 1 
#' @param rho 	true difference in simulated and observed autocorrelation at lag 1
#' @param tau.u	upper tolerance of the equivalence region
#' @param alpha	level of the equivalence test
#' @param s 	standard deviation of the test statistic
#' @return power of the asymptotic test. this is approximate because the test is asymptotic
#' @examples	tau.u<- 0.09
#' 	tau.l<- -tau.u
#' 	sim.n<-	5e3
#' 	rho<- seq(tau.l,tau.u,0.001)
#' 	pw<- nabc.acf.equivalence.pow(rho, tau.u, alpha, 1/sqrt(floor(sim.n/3)-3))
nabc.acf.equivalence.pow<- function(rho, tau.u, alpha, s)
{ 
	tmp<- pnorm( (tau.u-rho)/s + qnorm(alpha) ) - pnorm( -(tau.u+rho)/s - qnorm(alpha) )
	tmp[tmp<=0]<- 0
	tmp
}
#------------------------------------------------------------------------------------------------------------------------
#' Compute the autocorrelation in a time series along with some other info
#' @param x 		time series (simply vector)
#' @param leave.out thinning, how many values in the pair sequence (x_i,x_i-1) should be left out. Defaults to zero.
#' @return	\item{cor}{autocorrelation in the thinned sequence}
#' 	\item{z}{Z-transformation of the autocorrelation (this is atanh of "cor")}
#' 	\item{n}{Number of pairs (x_i,x_i-1) after thinning}
#' @examples nabc.acf.equivalence.cor(rnorm(100,0,1), leave.out=2)
nabc.acf.equivalence.cor<- function(x, leave.out=0)
{
	tmp<-	rbind( x[-1], x[-length(x)]	)
	tmp<-	tmp[,seq.int(1,ncol(tmp),by=1+leave.out)]
	tmp2<-	cor(tmp[1,],tmp[2,])
	ans<-	c(tmp2,		.5 * log( (1+tmp2)/(1-tmp2) ), 		ncol(tmp)	)
	names(ans)<- c("cor","z","n")
	ans
}
#------------------------------------------------------------------------------------------------------------------------
#' Compute the ABC tolerances of the asymptotic equivalence test for autocorrelations at lag 1
#' @param tau.l	lower tolerance of the equivalence region
#' @param tau.u	upper tolerance of the equivalence region
#' @param n		number of pairs (x_i,x_i-1) after thinning of the time series x_1, x_2, ...
#' @param alpha	level of the equivalence test
#' @return vector of length 2, first entry is lower ABC tolerance, second entry is upper ABC tolerance
#' @examples	tau.u<- 0.09
#' 	tau.l<- -tau.u
#' 	sim.n<-	5e3
#' leave.out<- 2
#' nabc.acf.equivalence.abctol(tau.l, tau.u, floor(sim.n / (1+leave.out)), 0.01)
nabc.acf.equivalence.abctol<- function(tau.l, tau.u, n, alpha)
{		
	z.isd<- sqrt( n - 3 )	 
	c(min(tau.l*z.isd+qnorm(1-alpha),0), max(0,tau.u*z.isd+qnorm(alpha)))
}
#------------------------------------------------------------------------------------------------------------------------
#' Calibrate the equivalence region of the asymptotic equivalence test for autocorrelations at lag 1 for given maximum power
#' @param mx.pw		maximum power at the point of reference (rho.star).
#' @param tau.up.ub	guess on an upper bound on the upper tolerance of the equivalence region
#' @param n			number of pairs (x_i,x_i-1) after thinning of the time series x_1, x_2, ...
#' @param alpha		level of the equivalence test
#' @param rho.star	point of reference. Defaults to the point of equality rho.star=0.
#' @param tol		this algorithm stops when the actual maximum power is less than 'tol' from 'mx.pw'
#' @param max.it	this algorithm stops prematurely when the number of iterations to find the equivalence region exceeds 'max.it'
#' @return	vector of length 4
#' 	\item{1}{lower tolerance of the equivalence region}		
#' 	\item{2}{upper tolerance of the equivalence region}
#' 	\item{3}{actual maximum power associated with the equivalence region}
#' 	\item{4}{error ie abs(actual power - mx.pw)}
#' @examples tau.u<- 0.09
#' 	tau.l<- -tau.u
#' 	sim.n<-	5e3
#' 	leave.out<- 2
#' 	nabc.chisqstretch.tau.lowup(0.9, 2, floor(sim.n / (1+leave.out)), 0.01)
nabc.acf.equivalence.tau.lowup<- function(mx.pw, tau.up.ub, n, alpha, rho.star=0, tol= 1e-5, max.it=100)
{
	curr.mx.pw<- 0
	tau.up.ub<- tau.up.ub/2
	tmp<- max.it
	while(curr.mx.pw<mx.pw && tmp>0)
	{
		tmp<- tmp-1
		tau.up.ub<- 2*tau.up.ub
		curr.mx.pw<- nabc.acf.equivalence.pow(rho.star, tau.up.ub, alpha, 1/sqrt(n-3))		
	}
	if(tmp==0)	stop("could not find tau.up.ub")	
#print(tau.up.ub); stop()
	tau.up.lb<- 0
	error<- 1	
	while(abs(error)>tol && round(tau.up.lb,d=10)!=round(tau.up.ub,d=10) && max.it>0)
	{
		max.it<- max.it-1
		tau.up<- (tau.up.lb + tau.up.ub)/2
		curr.mx.pw<- nabc.acf.equivalence.pow(rho.star, tau.up, alpha, 1/sqrt(n-3))
		error<- curr.mx.pw - mx.pw
#print(c(curr.mx.pw, tau.up, tau.up.lb, tau.up.ub, max.it))
		if(error<0)
			tau.up.lb<- tau.up
		else
			tau.up.ub<- tau.up
#print(c(abs(error), round(tau.up.lb,d=10)!=round(tau.up.ub,d=10)) )	
	}
	if(max.it==0)	warning("reached max.it")
	c(-tau.up,tau.up,curr.mx.pw,abs(error))
}
#------------------------------------------------------------------------------------------------------------------------
#' Perform the asymptotic equivalence test for autocorrelations at lag 1
#' @param sim			simulated summary values
#' @param obs			observed summary values
#' @param args			argument that contains the equivalence region and the level of the test (see Examples). This is the preferred method for specifying arguments and overwrites the dummy default values
#' @param verbose		flag if detailed information on the computations should be printed to standard out
#' @param alpha			level of the equivalence test
#' @param leave.out		thinning, how many values in the pair sequence (x_i,x_i-1) should be left out. Defaults to zero.
#' @param normal.test	name of function with which normality of the summary values is tested
#' @return	vector containing
#' \item{error}{test statistic. Here, instead of T we return the p-value of the TOST.}
#' \item{cil}{lower ABC tolerance. Here, instead of c^- we return 0.}
#' \item{cir}{upper ABC tolerance. Here, instead of c^+ we return 'alpha'.}
#' \item{al}{free entry. Here set to c^-.}
#' \item{ar}{free entry. Here set to c^+.}
#' \item{mx.pw}{Maximum power at the point of equality}
#' \item{rho.mc}{sample estimate of 'rho'}
#' @examples leave.out<- 2
#' tau.u<- 0.09
#' alpha<- 0.01
#' n<- 5e3
#' sigma<- 1
#' a<- 0.1
#' args<- paste("acfequiv",leave.out,tau.u,alpha,sep='/')	
#' x<-	rnorm(n+1,0,sigma)
#' x<- x[-1] + x[-(n+1)]*a 
#' y<-	rnorm(n+1,0,sigma)
#' y<- y[-1] + y[-(n+1)]*a
#' nabc.acf.equivalence(y,x,args)
nabc.acf.equivalence<- function(sim, obs, args=NA, verbose= FALSE, alpha=0, leave.out=0, normal.test= "sf.test")
{
	#verbose<- 0
	if(any(is.na(sim)))	stop("nabc.acf.equivalence: error at 1a")
	if(any(is.na(obs)))	stop("nabc.acf.equivalence: error at 1b")
	if(length(obs)!=length(sim))		stop("nabc.acf.equivalence: error at 1c")
	if(length(sim)<5)	stop("nabc.acf.equivalence: error at 1d")
	if(!is.na(args))
	{
		args<- strsplit(args,'/')[[1]]
		if(length(args)==4)	
		{
			leave.out<- as.numeric( args[2] )
			tau.u<- as.numeric( args[3] )
			tau.l<- -tau.u
			alpha<- as.numeric( args[4] )
		}
		else if(length(args)==5)
		{
			leave.out<- as.numeric( args[2] )
			tau.l<- as.numeric( args[3] )
			tau.u<- as.numeric( args[4] )
			alpha<- as.numeric( args[5] )
		}
		else 
			stop("nabc.acf.equivalence: error at 1e")
		args<- args[1]
	}
	if(alpha<0 || alpha>1)		stop("nabc.acf.equivalence: error at 1f")
	if(tau.u<0 )		stop("nabc.acf.equivalence: error at 1g")
	if(tau.l>0 )		stop("nabc.acf.equivalence: error at 1h")
	if(leave.out<0)		stop("nabc.acf.equivalence: error at 1i")
	ans<- NABC.DEFAULT.ANS	
	ans["pfam.pval"]<-	nabc.get.pfam.pval(sim,normal.test) 						
	
	#compute z-scores for sim and obs
	z.sim<- nabc.acf.equivalence.cor(sim, leave.out=leave.out)
	z.obs<- nabc.acf.equivalence.cor(obs, leave.out=leave.out)
	z.isd<- sqrt(z.sim["n"]-3)
	#perform TOST -- this now only depends on z.sim z.obs z.isd   
	tmp<- c(z.obs["z"] - z.sim["z"] - tau.l, z.obs["z"] - z.sim["z"] - tau.u, z.obs["z"] - z.sim["z"]) * z.isd #compute ZL and ZU for TOST, add Z for tau=0
	#print(z.sim); print(z.obs); print(z.isd); print(tmp)	
	which.not.reject<- which( c(pnorm( tmp[1] )<1-alpha, pnorm( tmp[2] )>alpha) )
	#only if length==0, TOST will be rejected ie ABC accept
	if(length(which.not.reject)%in%c(0,2))
	{
		#decide which of the two pvalues are to be reported --> figure out which one is closer to boundary
		#both upper and lower test statistics indicate mean difference < -tau and mean difference > tau  	OR 		mean difference >= -tau and mean difference <= tau
		which.not.reject<- which.min(abs(c(1-alpha-pnorm( tmp[1] ), pnorm( tmp[2] )-alpha )))
	}
	#the pvalue of ZU is the lower tail, but for ZL it is the upper tail, so..
	ans["error"]<-			ifelse(which.not.reject==1, 		1-pnorm( tmp[which.not.reject] ),		pnorm( tmp[which.not.reject] ) )
	ans[c("cil","cir")]<-	c(0,alpha)		
	ans[c("al","ar")]<- 	c(min(tau.l*z.isd+qnorm(1-alpha),0), max(0,tau.u*z.isd+qnorm(alpha)))	#CL and CU of the rejection region of the standardized test statistic
	ans["mx.pow"]<- 		diff(pnorm( ans[c("al","ar")]))		
	ans["pval"]<-			( pnorm(tmp[3],0,1) - (1-ans["mx.pow"])/2 ) / ans["mx.pow"]				#rescaled p-value that is expected to follow U(0,1) under the point null hypothesis
	ans["lkl"]<- 			dnorm(tmp[3],0,1)
	ans["link.mc.sim"]<- 	z.sim["z"]
	ans["link.mc.obs"]<- 	z.obs["z"]
	ans["rho.mc"]<- 		z.sim["z"] - z.obs["z"]
	if(verbose)	cat(paste(paste("\n{",args,"<-list(sim.cor=",z.sim["cor"] ,", obs.cor=",z.obs["cor"] ,", ZL=",tmp[1] ,", ZU=",tmp[2] ,", alpha=",alpha,", tau.l=",tau.l,", tau.u=",tau.u,", ",sep=''),paste(names(ans), ans, collapse=', ', sep='='),")}",sep=''))

	ans
}
#------------------------------------------------------------------------------------------------------------------------
#' Compute power of the exact equivalence test for dispersion 
#' @param rho 	true ratio in simulated variance / observed variance
#' @param df	degrees of freedom
#' @param cl	lower ABC tolerance
#' @param cu	upper ABC tolerance 
#' @return power of the exact test. this is exact.
#' @examples	alpha<- 0.01						
#'	tau.up<- 1.09
#'	yn<- 5e3
#'	tau.low<- nabc.chisqstretch.tau.low(tau.up, yn-1, alpha)
#'	rej<- .Call("abcScaledChiSq",	c(yn-1,tau.low,tau.up,alpha,1e-10,100,0.05)	)
#'	rho<- seq(tau.low,tau.up,by=0.001)
#'	nabc.chisqstretch.pow(rho,yn-1,rej[1],rej[2])		
nabc.chisqstretch.pow<- function(rho, df, cl, cu)
{
	pchisq( cu/rho*df, df ) - pchisq( cl/rho*df, df)
}
#------------------------------------------------------------------------------------------------------------------------
#' Calibrate the lower tolerance interval of the equivalence region for the test of dispersion equivalence
#' @param tau.up	upper tolerance of the equivalence region
#' @param df		degrees of freedom
#' @param alpha		level of the equivalence test
#' @param rho.star	point of reference. Defaults to the point of equality rho.star=1
#' @param tol		this algorithm stops when the actual point of reference is less than 'tol' from 'rho.star'
#' @param max.it	this algorithm stops prematurely when the number of iterations to find the equivalence region exceeds 'max.it'
#' @return tau.low, lower tolerance of the equivalence region
#' @examples	tau.u<- 2.2
#'  yn<- 60
#'	tau.l<- nabc.chisqstretch.tau.low(tau.u, yn-1, 0.01)
nabc.chisqstretch.tau.low<- function(tau.up, df, alpha, rho.star=1, tol= 1e-5, max.it=100) 
{
	tau.low.lb<-	1/tau.up		#1/tau.up gives rho.max<1 so we know we only need to go up
	tau.low.ub<- 	1				#tau.low must be between [tau.low.lb, tau.low.ub]
	error<- 1
	while(abs(error)>tol && round(tau.low.lb,d=10)!=round(tau.low.ub,d=10) && max.it>0)
	{
		max.it<- max.it-1
		tau.low<- (tau.low.lb + tau.low.ub)/2
		rej<- .Call("abcScaledChiSq",	c(df,tau.low,tau.up,alpha,1e-10,100,0.05)	)
		rho<- seq(tau.low,tau.up,by=0.001)
		pw<- nabc.chisqstretch.pow(rho,df,rej[1],rej[2])
#print( c(rho[ which.max(pw) ],pw[ which.max(pw) ], tau.low.lb, tau.low.ub,round(tau.low.lb,d=10)==round(tau.low.ub,d=10) ))	
		error<- rho[ which.max(pw) ] - rho.star
#print( error )			
		if(error<0)
			tau.low.lb<- tau.low
		else
			tau.low.ub<- tau.low			
	}
	if(max.it==0)	warning("nabc.chisqstretch.tau.low: reached max.it")
	tau.low
}
#------------------------------------------------------------------------------------------------------------------------
#' Calibrate the equivalence region for the test of dispersion equivalence for given maximum power
#' @param mx.pw		maximum power at the point of reference (rho.star).
#' @param tau.up.ub	guess on an upper bound on the upper tolerance of the equivalence region
#' @param df		degrees of freedom
#' @param alpha		level of the equivalence test
#' @param rho.star	point of reference. Defaults to the point of equality rho.star=1.
#' @param tol		this algorithm stops when the actual maximum power is less than 'tol' from 'mx.pw'
#' @param max.it	this algorithm stops prematurely when the number of iterations to find the equivalence region exceeds 'max.it'
#' @return	vector of length 4
#' 	\item{1}{lower tolerance of the equivalence region}		
#' 	\item{2}{upper tolerance of the equivalence region}
#' 	\item{3}{actual maximum power associated with the equivalence region}
#' 	\item{4}{error ie abs(actual power - mx.pw)}
#' @examples yn<- 60
#' 	nabc.chisqstretch.tau.lowup(0.9, 2.5, yn-1, 0.01)
nabc.chisqstretch.tau.lowup<- function(mx.pw, tau.up.ub, df, alpha, rho.star=1, tol= 1e-5, max.it=100)
{
	curr.mx.pw	<- 0
	tau.up.ub	<- tau.up.ub/2
	tmp			<- max.it
	while(curr.mx.pw<mx.pw && tmp>0)
	{
		tmp			<- tmp-1
		tau.up.ub	<- 2*tau.up.ub
		tau.low		<- nabc.chisqstretch.tau.low(tau.up.ub, df, alpha, rho.star=rho.star, tol=tol, max.it=max.it)
		rej			<- .Call("abcScaledChiSq",	c(df,tau.low,tau.up.ub,alpha,1e-10,100,0.05)	)
		curr.mx.pw	<- nabc.chisqstretch.pow(rho.star, df, rej[1], rej[2])		
	}
	if(tmp==0)	stop("nabc.chisqstretch.tau.lowup: could not find tau.up.ub")	
	#print(tau.up.ub)
	tau.up.lb	<- 1
	error		<- 1	
	while(abs(error)>tol && round(tau.up.lb,d=10)!=round(tau.up.ub,d=10) && max.it>0)
	{
		max.it	<- max.it-1
		tau.up	<- (tau.up.lb + tau.up.ub)/2
		tau.low	<- nabc.chisqstretch.tau.low(tau.up, df, alpha, rho.star=rho.star, tol=tol, max.it=max.it)
		rej		<- .Call("abcScaledChiSq",	c(df,tau.low,tau.up,alpha,1e-10,100,0.05)	)
		curr.mx.pw<- nabc.chisqstretch.pow(rho.star, df, rej[1], rej[2])
		error	<- curr.mx.pw - mx.pw
#print(c(curr.mx.pw, tau.low, tau.up, tau.up.lb, tau.up.ub, max.it))
		if(error<0)
			tau.up.lb<- tau.up
		else
			tau.up.ub<- tau.up
#print(c(abs(error), round(tau.up.lb,d=10)!=round(tau.up.ub,d=10)) )	
	}
	if(max.it==0)	warning("nabc.chisqstretch.tau.lowup: reached max.it")
	c(tau.low,tau.up,curr.mx.pw,abs(error))
}
#------------------------------------------------------------------------------------------------------------------------
#' Calibrate the number of simulated summary values and the equivalence region for the test of dispersion equivalence
#' @param n.of.x	number of observed summary values
#' @param s.of.Sx	standard deviation in the observed summary likelihood
#' @param mx.pw		maximum power at the point of reference (rho.star).
#' @param alpha		level of the equivalence test
#' @param tau.up.ub	guess on an upper bound on the upper tolerance of the equivalence region
#' @param tol		this algorithm stops when the actual variation in the ABC approximation to the summary likelihood is less than 'tol' from 's.of.Sx*s.of.Sx'
#' @param max.it	this algorithm stops prematurely when the number of iterations to calibrate the number of simulated data points exceeds 'max.it'
#' @return	vector of length 8
#' 	\item{1}{number of simulated summary values}
#' 	\item{2}{lower tolerance of the equivalence region}		
#' 	\item{3}{upper tolerance of the equivalence region}
#' 	\item{4}{lower ABC tolerance c^-}		
#' 	\item{5}{upper ABC tolerance c^+}
#' 	\item{6}{actual variation of the power}
#' 	\item{7}{actual maximum power associated with the equivalence region}
#' 	\item{8}{error ie abs(actual variation - variation in the observed summary likelihood)}
#' @examples xn<- 60; alpha	<- 0.01; prior.u	<- 3; prior.l	<- 1/3; tau.u<- 2.5; xsig2	<- 1		
#' 	#summary likelihood of sigma2 given sample mean and sum of squares
#' 	th		<- seq(prior.l,prior.u,length.out=1e3)
#' 	shape		<- (xn-2)/2	 
#' 	scale		<- xsig2*xn*xn/(xn-1)/2
#' 	y		<- densigamma(th, shape, scale)
#' 	var.Sx	<- scale*scale/((shape-1)*(shape-1)*(shape-2))		
#' 	#abc approximation to summary likelihood
#' 	nabc.chisqstretch.n.of.y(xn, sqrt(var.Sx), 0.9, alpha, tau.u.ub=tau.u)
#' 	yn		<- tmp[1]
#' 	tau.l	<- tmp[2]
#' 	tau.u	<- tmp[3]
#' 	c.l		<- tmp[4]
#' 	c.u		<- tmp[5]
#' 	y2		<- nabc.chisqstretch.pow(th, yn-1, c.l, c.u)
#' #plot the summary likelihood and the abc approximation
#' plot(th,y/mean(y),ylim=range(c(y/mean(y),y2/mean(y2))),type='l')
#' lines(th,y2/mean(y2),col="blue")
nabc.chisqstretch.n.of.y<- function(n.of.x, s.of.Sx, mx.pw, alpha, tau.u.ub=2, tol= 1e-5, max.it=100)
{
	s2.of.Sx	<- s.of.Sx*s.of.Sx
	pw.cvar		<- 2*s2.of.Sx
	curr.mx.pw	<- 0
	yn.ub		<- round( n.of.x/2 )
	curr.it		<- max.it
	tau.u		<- tau.u.ub
	while(pw.cvar>s2.of.Sx && curr.it>0)
	{
		curr.it	<- curr.it-1
		yn.ub	<- 2*yn.ub
#print(c(mx.pw,yn.ub))	
		tmp		<- nabc.chisqstretch.tau.lowup(mx.pw, 2*tau.u, yn.ub-1, alpha )
#print(c("OK ",tmp))		
		tau.l	<- tmp[1]
		tau.u	<- tmp[2]
		pw.cmx	<- tmp[3]
		rho		<- seq(tau.l/2,2*tau.u,length.out=1e3)
		tmp		<- .Call("abcScaledChiSq",	c(yn.ub-1,tau.l,tau.u,alpha,1e-10,100,0.05)	)
		pw		<- nabc.chisqstretch.pow(rho, yn.ub-1, tmp[1], tmp[2])
		pw.cme	<- sum(rho*pw) / sum(pw)
		pw.cvar	<- sum((rho-pw.cme)*(rho-pw.cme)*pw) / sum(pw)			
#print(c(yn.ub, tau.l, tau.u, pw.cvar, s2.of.Sx, pw.cmx ))		
	}
	if(curr.it==0)	stop("nabc.chisqstretch.n.of.y: could not find upper bound for yn")	
#print(c(yn.ub, tau.u, pw.cvar, s2.of.Sx, pw.cmx ))
	yn.lb	<- n.of.x
	error	<- 1	
	while(abs(error)>tol && (yn.lb+1)!=yn.ub && max.it>0)
	{
		max.it	<- max.it-1
		yn		<- round( (yn.lb + yn.ub)/2 )
		tmp		<- nabc.chisqstretch.tau.lowup(mx.pw, 2*tau.u, yn-1, alpha )
		tau.l	<- tmp[1]
		tau.u	<- tmp[2]
		pw.cmx	<- tmp[3]
		rho		<- seq(tau.l/2,2*tau.u,length.out=1e3)
		tmp		<- .Call("abcScaledChiSq",	c(yn-1,tau.l,tau.u,alpha,1e-10,100,0.05)	)
		c.l		<- tmp[1]
		c.u		<- tmp[2]
		pw		<- nabc.chisqstretch.pow(rho, yn-1, c.l, c.u)
		pw.cme	<- sum(rho*pw) / sum(pw)
		pw.cvar	<- sum((rho-pw.cme)*(rho-pw.cme)*pw) / sum(pw)
		error	<- pw.cvar - s2.of.Sx
#print(c(pw.cvar, yn, yn.lb, yn.ub, max.it))
		if(error<0)
			yn.ub<- yn
		else
			yn.lb<- yn
#print(c(abs(error), (yn.lb+1)!=yn.ub) )	
	}
	c(yn,tau.l,tau.u,c.l,c.u,pw.cvar,pw.cmx,abs(error))
}
#------------------------------------------------------------------------------------------------------------------------
#' Perform the exact test for dispersion equivalence when the summary values are normally distributed
#' @param sim			simulated summary values
#' @param obs.mc		variance of the observed summary values
#' @param args			argument that contains the equivalence region and the level of the test (see Examples). This is the preferred method for specifying arguments and overwrites the dummy default values
#' @param verbose		flag if detailed information on the computations should be printed to standard out
#' @param tau.l			lower tolerance of the equivalence region
#' @param tau.u			upper tolerance of the equivalence region
#' @param guess.tau.l	guess on the lower tolerance of the equivalence region. Used when the tolerances are annealed and calibration is numerically unstable.
#' @param alpha			level of the equivalence test
#' @param leave.out		thinning, how many values in the pair sequence (x_i,x_i-1) should be left out. Defaults to zero.
#' @param normal.test	name of function with which normality of the summary values is tested
#' @return	vector containing
#' \item{error}{test statistic, here var(sim)/obs.mc}
#' \item{cil}{lower ABC tolerance c^-}
#' \item{cir}{upper ABC tolerance c^+}
#' \item{mx.pw}{Maximum power at the point of equality}
#' \item{rho.mc}{log(var(sim) / obs.mc)}
#' @examples alpha<- 0.01; xn<- yn<- 60; xsigma2<- 1; tau.u<- 2.2
#'	tau.l<- nabc.chisqstretch.tau.low(tau.u, yn-1, alpha)		
#'	args<- paste("chisqstretch",tau.l,tau.u,alpha,sep='/')
#'	x<- rnorm(xn,0,sd=sqrt(xsigma2))
#'	y<- rnorm(yn,0,sd=sqrt(xsigma2))
#'	nabc.chisqstretch(y, var(x), args=args, verbose= 0)
nabc.chisqstretch<- function(sim, obs.mc, args=NA, verbose= FALSE, tau.l=1, tau.u=1, guess.tau.l=0, alpha=0, normal.test= "sf.test")
{
	#verbose<- 1
	#sim<- rnorm(100, 8.1, 1.1)
	if(any(is.na(sim)))	stop("get.dist.chisqstretch: error at 1a")
	if(length(obs.mc)>1 || is.na(obs.mc))	stop("get.dist.chisqstretch: error at 1b")
	
	df.sim<- length(sim) - 1
	if(!is.na(args))
	{
		args<- strsplit(args,'/')[[1]]
		if(length(args)==3)
		{
			guess.tau.l<- as.numeric( args[2] )
			tau.u<- tau.l<- NA
			alpha<- as.numeric( args[3] )
		}		
		if(length(args)==4)	
		{
			guess.tau.l<- as.numeric( args[2] )
			tau.u<- as.numeric( args[3] )
			tau.l<- NA
			alpha<- as.numeric( args[4] )
		}
		else if(length(args)==5)
		{
			guess.tau.l<- as.numeric( args[2] )
			tau.l<- as.numeric( args[3] )
			tau.u<- as.numeric( args[4] )
			alpha<- as.numeric( args[5] )
		}
		else 
			stop("get.dist.chisqstretch: error at 1A")
		if(is.na(tau.u))
		{
			tmp<- nabc.chisqstretch.tau.lowup(0.9, 2, df.sim, alpha)
			tau.l<- tmp[1]
			tau.u<- tmp[2]						
		}
		else if(is.na(tau.l) && guess.tau.l && tau.u>200)		#if tau.u too large, then "nabc.chisqstretch.tau.low" may take a while 
			tau.l<- 1/tau.u
		else if(is.na(tau.l))
			tau.l<- nabc.chisqstretch.tau.low(tau.u, df.sim, alpha)
		#print(c(df.sim,tau.l,tau.u,alpha))			
		args<- args[1]
	}
	if(alpha<0 || alpha>1)		stop("get.dist.chisqstretch: error at 1e")
	if(tau.u<1 )		stop("get.dist.chisqstretch: error at 1f")
	if(tau.l>1 )		stop("get.dist.chisqstretch: error at 1g")
		 
	ans<- NABC.DEFAULT.ANS	
	ans["pfam.pval"]<-	nabc.get.pfam.pval(sim,normal.test) 
			
	#get confidence intervals by numerical approximation	
	tmp<- .Call("abcScaledChiSq",	c(df.sim,tau.l,tau.u,alpha,1e-10,100,0.05)	)
	if(tmp[4]>1e-10)	stop("get.dist.chisqstretch: error at 3a")
	ans[c("cil","cir","mx.pow")]<- tmp[1:3]
	
	
	ans["error"]<- 			var(sim) / obs.mc
	ans["lkl"]<- 			dchisq(ans["error"]*df.sim,df.sim)	
	ans["pval"]<- 			pchisq(ans["error"]*df.sim,df.sim)
	ans[c("al","ar")]<- 	c(0, 1 - diff( pchisq(ans[c("cil","cir")]*df.sim,df.sim) ) )
	ans["pval"]<-			( ans["pval"] - ans["ar"]/2 ) / ( 1 - ans["ar"] )
	ans["link.mc.sim"]<- 	var(sim)
	ans["link.mc.obs"]<- 	obs.mc
	ans["rho.mc"]<- log(var(sim) / obs.mc)
	if(verbose)	cat(paste(paste("\n{",args,"<-list(sim.var=",var(sim) ," , obs.var=",obs.mc," , alpha=",alpha," , tau.l=",tau.l," , tau.u=",tau.u,", log.ciu=",log(ans["cir"]),", ",sep=''),paste(names(ans), ans, collapse=', ', sep='='),")}",sep=''))
	ans
}	
#------------------------------------------------------------------------------------------------------------------------
nabc.get.locfit.links<- function(th.d, m, th.thin= 1, th.sep=100)
{
	require(locfit)
	th.names	<- colnames(m)[1:th.d]
	su.idx		<- (th.d+1):ncol(m)
	su.names	<- colnames(m)[su.idx]
	colnames(m)[1:th.d]<- paste("th",1:th.d,sep='')
	colnames(m)[su.idx]<- paste("su",seq_along(su.idx),sep='')
	links	<- lapply(seq_along(su.idx),function(k)
			{
				colnames(m)[su.idx]<- paste("su",seq_along(su.idx),sep='')
				colnames(m)[th.d+k]<- "rho"				
				lnk<-try(locfit(rho~th1:th2, data=m[seq.int(1,nrow(m),by=th.thin),],mint=20,deg=2,alpha=2,kern="epan"))
				if(inherits(lnk, "try-error"))
					lnk<- locfit(rho~th1:th2, data=m[seq.int(1,nrow(m),by=th.thin),],maxk=200,mint=20,deg=2,alpha=2,kern="epan")	
				#do locpoly regression on loc poly regression to save memory space
				theta<- expand.grid(lfmarg(lnk$box, rep(th.sep, lnk$mi["d"])))
				colnames(theta)<- paste("th",1:th.d,sep='')
				rhok<- predict(lnk, theta)
				lnk<-try(locfit(rho~th1:th2, data=cbind(theta,rho=rhok),mint=20,deg=2,alpha=2,kern="epan"))
				if(inherits(lnk, "try-error"))
					lnk<- locfit(rho~th1:th2, data=cbind(theta,rho=rhok),maxk=200,mint=20,deg=2,alpha=2,kern="epan")
				lnk	
			})
	links
}
#------------------------------------------------------------------------------------------------------------------------
nabc.get.locfit.jac<- function(th.d, m, th.thin= 1, th.sep=100)
{
	require(locfit)
	th.idx		<- 1:th.d
	th.names	<- colnames(m)[th.idx]
	su.idx		<- (th.d+1):ncol(m)
	su.names	<- colnames(m)[su.idx]
	colnames(m)[1:th.d]<- paste("th",th.idx,sep='')
	colnames(m)[su.idx]<- paste("su",seq_along(su.idx),sep='')
	jac	<- lapply(seq_along(su.idx),function(k)
			{
				jacrow<- lapply(th.idx,function(d)
						{
							colnames(m)[su.idx]<- paste("su",seq_along(su.idx),sep='')
							colnames(m)[th.d+k]<- "rho"				
							partialder<-try(locfit(rho~th1:th2, data=m[seq.int(1,nrow(m),by=th.thin),],mint=20,deg=2,alpha=2,kern="epan",deriv=paste("th",d,sep='')))
							if(inherits(partialder, "try-error"))
								partialder<- locfit(rho~th1:th2, data=m[seq.int(1,nrow(m),by=th.thin),],maxk=200,mint=20,deg=2,alpha=2,kern="epan",deriv=paste("th",d,sep=''))	
							#do locpoly regression on loc poly regression to save memory space
							theta<- expand.grid(lfmarg(partialder$box, rep(th.sep, partialder$mi["d"])))
							colnames(theta)<- paste("th",1:th.d,sep='')
							rhok<- predict(partialder, theta)
							partialder<-try(locfit(rho~th1:th2, data=cbind(theta,rho=rhok),mint=20,deg=2,alpha=2,kern="epan",deriv=paste("th",d,sep='')))
							if(inherits(partialder, "try-error"))
								partialder<- locfit(rho~th1:th2, data=cbind(theta,rho=rhok),maxk=200,mint=20,deg=2,alpha=2,kern="epan",deriv=paste("th",d,sep=''))
							partialder	
						})		
				names(jacrow)<- paste("dth",th.idx,sep='')
				jacrow
			})
	names(jac)<- paste("dL",su.idx,sep='')
	jac
}
#------------------------------------------------------------------------------------------------------------------------
checkjac<- function(a,sig2,ax,vx)
{
	if(length(sig2)>1 || length(ax)>1 || length(vx)>1)	stop("checkjac not vectorized")
	c( 2*a*sig2/vx, (1+a*a)/vx, (1-a*a)/(1+a*a+a^4), rep(0,length(a)))	
}
#------------------------------------------------------------------------------------------------------------------------
checklink<- function(a,sig2,ax,vx)
{
	c((1+a*a)*sig2/vx, project.nABC.movingavg.a2rho(a)-project.nABC.movingavg.a2rho(ax) )
}
#------------------------------------------------------------------------------------------------------------------------
nabc.estimate.jac<- function( links, th.eval, th.dh, ax, vx )
{
	th.d	<- ncol(th.eval)
	#print(th.d)
	if(length(th.dh)==1)	th.dh<- rep(th.dh,ncol(th.eval))
	if(length(th.dh)!=th.d)	stop("nabc.estimate.Jac: incorrect length of dh")	
	if(th.d!=2)	stop("nabc.estimate.Jac: incorrect th.d")		
	th.dh	<- matrix(c(th.dh[1],0, -th.dh[1],0, 0,th.dh[2], 0,-th.dh[2]), nrow=4, byrow=1)	#TODO th.d for higher th.d
#	print(th.dh); #plot(links[[2]]); print(th.eval)
	#th.eval<- th.eval[150:160,]
	jac		<- apply(th.eval,1,function(theta)
			{
#print(theta)
				if(0)			#exact to double check
				{
					diffq<- t(sapply(seq_along(links),function(k)
						{				
							tmp<- matrix(rep(theta,length(nrow(th.dh))),nrow=nrow(th.dh),ncol=th.d,byrow=1)+th.dh
#print(tmp)							
							rhok<- sapply(seq_len(nrow(tmp)),function(i)	checklink(tmp[i,1],tmp[i,2],ax, vx)		)[k,]
#print(rhok)
#stop()
							tmp2<- matrix(links[[k]]$box,ncol=th.d,byrow=1)
							tmp2<- apply(tmp,1,function(row){		all(row>=links[[k]]$box[1:th.d]	& row<=links[[k]]$box[(th.d+1):(2*th.d)]) 	})
							rhok[!tmp2]<- NA
							#print(tmp); print( tmp2 ); print(rhok);  							
							if(length(rhok)!=nrow(th.dh))	stop("nabc.estimate.Jac: unexpected behavior of predict.locfit - inappropriate links?")							
							rhok<- -apply(matrix(rhok,nrow=2), 2, diff)					#L_k(theta+h_d) - L_k(theta-h_d) for all d
							tmp<- rhok/apply(matrix(abs(th.dh[,1]+th.dh[,2]),nrow=2),2,sum)	#difference quotient for L_k and all d
							#print(tmp)
							tmp
						}))
#print(diffq); print(det(diffq))		
				}
#	
#stop()
				if(1)
				{
				diffq<- t(sapply(seq_along(links),function(k)
						{				
							tmp<- matrix(rep(theta,length(nrow(th.dh))),nrow=nrow(th.dh),ncol=th.d,byrow=1)+th.dh
							rhok<- predict(links[[k]], tmp )
#print(rhok)										
							tmp2<- matrix(links[[k]]$box,ncol=th.d,byrow=1)
							tmp2<- apply(tmp,1,function(row){		all(row>=links[[k]]$box[1:th.d]	& row<=links[[k]]$box[(th.d+1):(2*th.d)]) 	})
							rhok[!tmp2]<- NA
							#print(tmp); print( tmp2 ); print(rhok);  							
							if(length(rhok)!=nrow(th.dh))	stop("nabc.estimate.Jac: unexpected behavior of predict.locfit - inappropriate links?")							
							rhok<- -apply(matrix(rhok,nrow=2), 2, diff)					#L_k(theta+h_d) - L_k(theta-h_d) for all d
							tmp<- rhok/apply(matrix(abs(th.dh[,1]+th.dh[,2]),nrow=2),2,sum)	#difference quotient for L_k and all d
							#print(tmp)
							tmp
						}))
#print(diffq); print(det(diffq))		
				}

#stop()
				rownames(diffq)<- paste("L",seq_along(links),sep='')
				colnames(diffq)<- paste("th",1:th.d,sep='')
				#diffq[2,2]<- 0
				abs(det(diffq))								
			})
	jac
}
#------------------------------------------------------------------------------------------------------------------------
nabc.get.jacobian.2d<- function(m, th.mode, th.sep= rep(10,2), th.thin= 2000)
{
	require(locfit)
	th.d		<- 2	
	th.names	<- colnames(m)[1:th.d]
	su.idx		<- (th.d+1):ncol(m)
	su.names	<- colnames(m)[su.idx]
	colnames(m)[1:th.d]<- paste("th",1:th.d,sep='')
	colnames(m)[su.idx]<- paste("su",seq_along(su.idx),sep='')
	
	#print(colnames(m))
	th.range<- apply(m[,1:th.d],2,range)
	#print(th.range)
	th.dh	<- apply(th.range,2,diff) / th.sep
	#print(th.mode)
	th.eval	<- (th.range-rep(th.mode,each=2))/2		#+rep(th.mode,each=2)	
	th.eval	<- cbind(	th.mode+c(th.eval[1,1],0),th.mode+c(th.eval[2,1],0),
						th.mode+c(0,th.eval[1,2]),th.mode+c(0,th.eval[2,2]), th.mode		)
	th.dh	<- matrix(c(th.dh[1],0, -th.dh[1],0, 0,th.dh[2], 0,-th.dh[2]), ncol=4)				
	#print(th.dh)
	#print(th.eval)
	#remains to compute difference quotient to all rows in 'th.eval', using the dh in 'th.dh'
	#now first construct approximations to L_k
	links	<- lapply(seq_along(su.idx),function(k)
			{
				colnames(m)[su.idx]<- paste("su",seq_along(su.idx),sep='')
				colnames(m)[th.d+k]<- "rho"
				locfit(rho~th1:th2, data=m[seq.int(1,nrow(m),by=th.thin),])				
			})
	jac		<- apply(th.eval,2,function(theta)
			{
				#print(theta)
				diffq<- sapply(seq_along(links),function(k)
					{
						rhok<- apply(th.dh,2,function(h)		#the rho_k's for each of theta+h_1, theta-h_1, theta+h_2, theta-h_2, ...
							{																
								predict(links[[k]], matrix(theta+h, nrow=1, dimnames=list(c(),paste("th",1:th.d,sep=''))))								
							})					
						rhok<- -apply(matrix(rhok,nrow=2), 2, diff)					#L_k(theta+h_d) - L_k(theta-h_d) for all d
						rhok/apply(matrix(abs(th.dh[1,]+th.dh[2,]),nrow=2),2,sum)	#difference quotient for L_k and all d												
					})
				colnames(diffq)<- paste("L",seq_along(links),sep='')
				rownames(diffq)<- paste("th",1:th.d,sep='')
				abs(det(t(diffq)))				
			})
	ans		<- rbind(th.eval, jac)
	rownames(ans)<- c(th.names,"jac")
	ans
}
#------------------------------------------------------------------------------------------------------------------------
get.dist.fstretch<- function(sim, obs, args=NA, verbose= FALSE, alpha=0, tau.l=1, tau.u=1, plot=FALSE, xlab= NA, nbreaks=40, normal.test= "sf.test")
{
	#verbose<- 1
	#sim<- rnorm(100, 8.1, 1.1)
	if(any(is.na(sim)))	stop("get.dist.fstretch: error at 1a")
	if(any(is.na(obs)))	stop("get.dist.fstretch: error at 1b")	
	if(length(obs)!=length(sim))		stop("get.dist.fstretch: error at 1c")
	if(!is.na(args))
	{
		args<- strsplit(args,'/')[[1]]
		if(length(args)==3)	
		{
			tau.u<- as.numeric( args[2] )
			tau.l<- 1/tau.u
			alpha<- as.numeric( args[3] )
		}
		else if(length(args)==4)
		{
			tau.l<- as.numeric( args[2] )
			tau.u<- as.numeric( args[3] )
			alpha<- as.numeric( args[4] )
		}
		else 
			stop("get.dist.fstretch: error at 1A")
		args<- args[1]
	}
	if(alpha<0 || alpha>1)		stop("get.dist.fstretch: error at 1e")
	if(tau.u<1 )		stop("get.dist.fstretch: error at 1f")
	if(tau.l>1 )		stop("get.dist.fstretch: error at 1g")
	if(!normal.test%in%c("shapiro.test","lillie.test","cvm.test","ad.test","pearson.test","sf.test"))		stop("get.dist.fstretch: error at 1h")
	if(normal.test%in%c("lillie.test","cvm.test","ad.test","pearson.test","sf.test"))		require(nortest, lib.loc=LIB.LOC)
	df.sim<- length(sim) - 1
	df.obs<- length(obs) - 1 
	ans<- NABC.DEFAULT.ANS
	
	ans["pfam.pval"]<- ifelse(any(diff(sim)>0), ifelse(normal.test%in%c("shapiro.test","sf.test") && length(sim)>5000, eval(call(normal.test,sim[1:5000]))$p.value, eval(call(normal.test,sim))$p.value), 0)

	#get confidence intervals by numerical approximation; this is the fstretch R command by Wellek
	if(tau.l==1/tau.u && df.sim==df.obs)
	{
		fstretch.solvefor.cir<- function(cir, tau.u, n) pf( cir / tau.u,  n, n) - pf( 1 / (cir*tau.u),  n, n) - alpha
		catch<- try({
			ans["cir"]<- uniroot(	fstretch.solvefor.cir,	c(1,100), tol=.Machine$double.eps^0.5, tau.u=tau.u, n=df.sim)$root
			ans["cil"]<- 1/ans["cir"]
		})
		if(inherits(catch, "try-error"))
		{
			geterrmessage()						#the upper and lower bound of 'fstretch.solvefor.cir' may not cross 0 if tau.u is too large
			ans[c("cil","cir")]<- c(0,1e5)
		}
		ans["mx.pow"]<- pf(ans["cir"],df.obs, df.sim) - pf(ans["cil"],df.obs,df.sim)
	}
	else
	{
		tmp<- .Call("abcScaledF",	c(df.obs,df.sim,tau.l,tau.u,alpha,1e-10,100,0.05)	)
		if(tmp[4]>1e-10)	stop("get.dist.fstretch: error at 3a")
		ans[c("cil","cir","mx.pow")]<- tmp[1:3]
		
	}
	ans["error"]<- 			var(obs) / var(sim)
	ans["lkl"]<- 			df(ans["error"],df.obs,df.sim)	
	ans["pval"]<- 			pf(ans["error"],df.obs,df.sim)
	ans[c("al","ar")]<- 	c(0, 1 - diff( pf(ans[c("cil","cir")],df.obs, df.sim) ) )
	ans["pval"]<-			( ans["pval"] - ans["ar"]/2 ) / ( 1 - ans["ar"] )
	ans["link.mc.sim"]<- 	var(sim)
	ans["link.mc.obs"]<- 	var(obs)
	ans["rho.mc"]<- log(var(obs) / var(sim))
	if(verbose)	cat(paste(paste("\n{",args,"<-list(sim.var=",var(sim) ," , obs.var=",var(obs)," , alpha=",alpha," , tau.l=",tau.l," , tau.u=",tau.u,", log.ciu=",log(ans["cir"]),", ",sep=''),paste(names(ans), ans, collapse=', ', sep='='),")}",sep=''))
	ans
}
#------------------------------------------------------------------------------------------------------------------------
get.dist.mwu.equivalence<- function(sim, obs, args= NA, verbose= FALSE, alpha= 0.05, tau= 0, tau.translate= 0, plot= FALSE, xlab= NA, nbreaks= 40)
{
	verbose<- 1
	get.dist.mwu.equivalence.psd<- function(sim,obs)
	{
		if(length(sim)==1) sd(obs)
		else if(length(obs)==1)	sd(sim)
		else	sqrt(1/length(sim)+1/length(obs)) * sqrt(		(  (length(sim)-1)*var(sim)+(length(obs)-1)*var(obs)  )		/		(length(sim)+length(obs)-2)				)
	}
	get.dist.mwu.equivalence.tau<- function(x,c=8)
	{
		ans<- numeric(length(x));
		tmp<- x<c;
		ans[tmp]<- pnorm(x[tmp])-.5;
		ans[!tmp]<- (x[!tmp]-(c-4))/8;
		ans
	}
	#compute two sample t-test on either z-scores or untransformed data points
	if(any(is.na(sim)))	stop("get.dist.mwu.equivalence: error at 1a")
	if(any(is.na(obs)))	stop("get.dist.mwu.equivalence: error at 1b")
	if(length(sim)<1)		stop("get.dist.mwu.equivalence: error at 1c")
	if(length(obs)<1)		stop("get.dist.mwu.equivalence: error at 1d")
	if(length(obs)<2 && length(sim)<2)		stop("get.dist.mwu.equivalence: error at 1e")
	#if not missing, arguments 'args' always overwrite 'alpha' and 'nc'
	#expect args string of the form studenttXXX/beta/non-centrality
	if(!is.na(args))
	{
		args<- strsplit(args,'/')[[1]]
		if(length(args)!=5)	stop("get.dist.mwu.equivalence: error at 2a")
		standardize<- as.numeric( args[2] )
		tau.translate<- as.numeric( args[3] )
		tau<- as.numeric( args[4] )
		alpha<- as.numeric( args[5] )
		args<- args[1]
	}
	if(!standardize%in%c(0,1))	stop("get.dist.mwu.equivalence: error at 2b")
	if(alpha<0 || alpha>1)		stop("get.dist.mwu.equivalence: error at 2c")
	if(tau<=0 )		stop("get.dist.mwu.equivalence: error at 2d")
	if(tau.translate && standardize)	stop("get.dist.mwu.equivalence: error at 2e")
	ans<- NABC.DEFAULT.ANS
	#transform tau from  cdf tolerance of pi+ into location tolerance based on 6.10 and 6.2 in Wellek 2010
	#we want tau for |nu_k(x)-nu_k(theta)|<tau_k. so need an estimate of sigma. we pool over sim and obs, assuming equal sample sizes.
	if(!tau.translate && !standardize)
	{
		if(tau>.5)	tau<- 0.5
		if(tau<=0)					stop("get.dist.mwu.equivalence: error at 3a")
		tau.loc<- qnorm(.5+tau)*get.dist.mwu.equivalence.psd(sim,obs)
	}
	else if(tau.translate)
	{
		tau.loc<- tau
		tau<- tau.loc / get.dist.mwu.equivalence.psd(sim,obs)
		if(tau>8)		#allow for tau>0.5 -- this does not make much sense for MWU but can be useful for annealing purposes. in any case raise a warning
		{
			tau<- get.dist.mwu.equivalence.tau(tau)
			options(warn=1)
			warning(paste("get.dist.mwu.equivalence: tau is larger than 0.5:",tau))
			options(warn=2)
		}
		else
			tau<- pnorm( tau ) - 0.5
	}
	else
		tau.loc<- tau
	if(length(sim)==1)	sim<- obs-median(obs)+sim		#apply shift hypothesis
	if(length(obs)==1)	obs<- sim-median(sim)+obs		#apply shift hypothesis
#print(obs); print(sim)
	if(0)
	{
		#compute mwu test statistic
		w<- sum( sapply(obs, function(x){			length(which(x>=sim))		}) )
		#compute PIxxy
		tmp<-   obs[-1]
		tmp2<- matrix(ncol= length(obs)-1, nrow= length(obs)-1)
		tmp2<- col(tmp2) >= row(tmp2)
		tmp2<-	c(  sapply(seq_along(tmp),function(i){     tmp[tmp2[i,]]      }), recursive=1 )	#elements of upper triangular matrix with entries x[2] .. x[m],       x[3] .. x[m], 	..., x[m,m]		in this order
		tmp<- rep(obs[-length(obs)], times= seq(length(obs)-1,1,-1))										#m-1 times x[1], m-2 times x[2], ..., 1 times x[m-1]	in this order
		#compute PIxxy and append to W
		w<- c(w, sum(sapply(seq_along(tmp), function(i){		length(which(tmp[i]>=sim & tmp2[i]>=sim))		})))
		#compute PIxyy
		tmp<-   sim[-1]
		tmp2<- matrix(ncol= length(sim)-1, nrow= length(sim)-1)
		tmp2<- col(tmp2) >= row(tmp2)
		tmp2<-	c(  sapply(seq_along(tmp),function(i){     tmp[tmp2[i,]]      }), recursive=1 )	#elements of upper triangular matrix with entries x[2] .. x[m],       x[3] .. x[m], 	..., x[m,m]		in this order
		tmp<- rep(sim[-length(sim)], times= seq(length(sim)-1,1,-1))										#m-1 times x[1], m-2 times x[2], ..., 1 times x[m-1]	in this order
		#compute PIxyy and append to W
		w<- c(w, sum(sapply(seq_along(tmp), function(i){		length(which(obs>=tmp[i] & obs>=tmp2[i]))		})))
		#check with simple implementation taken from mawi.R in Wellek 2010
		#m<- length(obs); n<- length(obs); x<- obs; y<- sim; wxy<- pihxxy<- pihxyy<- 0
		#for (i in 1:m) for (j in 1:n) wxy <- wxy + trunc(0.5*(sign(x[i] - y[j]) + 1))
		#for (i in 1:m) for (j1 in 1:(n-1)) for (j2 in (j1+1):n) pihxyy <- pihxyy + trunc(0.5*(sign(x[i] - max(y[j1],y[j2])) + 1))
		#for (i1 in 1:(m-1)) for (i2 in (i1+1):m) for (j in 1:n) pihxxy <- pihxxy + trunc(0.5*(sign(min(x[i1],x[i2]) - y[j]) + 1))
		#normalize
		tmp<- length(obs)
		tmp2<- length(sim)
		w<- w * c(  1/(tmp*tmp2),    2/(tmp*(tmp-1)*tmp2),    	2/(tmp2*(tmp2-1)*tmp)   )
		#estimate sd of mwu statistic
		w<- c(w,     	sqrt( (w[1] - (tmp+tmp2-1)*w[1]*w[1] + (tmp-1)*w[2]  + (tmp2-1)*w[3])/(tmp*tmp2) )    )
	}
	else		#faster than vectorized version. w[1]= W+  w[2]= (m-1)PIxxy   w[3]= (n-1)PIxyy  w[4]= sd(W+)   w[5]=W-  w[6]=sd(W-)
		w<- .Call("abcMWUE",obs,sim)
	#can compute lkl and pval at tau=0 explicitly for small sample sizes
	if(0)
	{	#this is slow and currently not needed
		options(warn=1)
		require(coin, warn.conflicts = FALSE)
		options(warn=2)

		tmp<- data.frame( val= c(obs,sim), group= c(rep("1",length(obs)), rep("2",length(sim))) )
		tmp<- wilcox_test(val ~ group, data = tmp,    distribution = "exact" )		#this assumes unpaired data
		idx<-  which(support(tmp)<=as.numeric(statistic(tmp,"standardized")))
		ans[c("lkl","pval")]<-	c( dperm( tmp, support(tmp)[	ifelse(!length(idx), 1, idx[length(idx)])	] ), pvalue(tmp) )
	}
#print(w)
	#compute test statistic and critical bound, assuming symmetric tolerance.
	#the CI regions in (6.19) are numerically not stable (large ncp). We use (4.5) and construct a TOST.
	#this leads to a numerically stable procedure
	if(is.na(w[4]) && sd(sim)==sd(obs))							#if variance in obs & sim equal, then sigma(w-)=sigma(w+), and sigma(w-) might not be NA
		w[4]<- w[6]
	if(is.na(w[4]) || !w[4] )
	{
		tmp<- rep(NA,2)
		ans["error"]<- ifelse(tau>=0.5, 1, 2)*alpha			#always accept when tau >= 0.5
	}
	else if(standardize)
	{
		tmp<- c( (w[1] - 0.5) / w[4] + tau, (w[1] - 0.5) / w[4] - tau)		#here, tau is the same as tau/sigma(W_+) below, and hence independent of y 
		which.not.reject<- which( c(pnorm( tmp[1] )<1-alpha, pnorm( tmp[2] )>alpha) )		
		if(length(which.not.reject)%in%c(0,2))
			which.not.reject<- which.min(abs(c(1-alpha-pnorm( tmp[1] ), pnorm( tmp[2] )-alpha )))
		ans["error"]<-	ifelse(which.not.reject==1, 		1-pnorm( tmp[which.not.reject] ),		pnorm( tmp[which.not.reject] ) )
	}
	else
	{
		tmp<- c( (w[1] - 0.5 + tau) / w[4], (w[1] - 0.5 - tau) / w[4] )		#compute ZL= ( W_+ - 0.5 + tau ) / sigma(W_+) and ZU for TOST; 		
		which.not.reject<- which( c(pnorm( tmp[1] )<1-alpha, pnorm( tmp[2] )>alpha) )
		#only if length==0, TOST will be rejected ie ABC accept
		if(length(which.not.reject)%in%c(0,2))
		{
			#decide which of the two pvalues are to be reported --> figure out which one is closer to boundary
			#both upper and lower test statistics indicate mean difference < -tau and mean difference > tau  	OR 		mean difference >= -tau and mean difference <= tau
			which.not.reject<- which.min(abs(c(1-alpha-pnorm( tmp[1] ), pnorm( tmp[2] )-alpha )))
		}
		#the pvalue of ZU is the lower tail, but for ZL it is the upper tail, so..
		ans["error"]<-	ifelse(which.not.reject==1, 		1-pnorm( tmp[which.not.reject] ),		pnorm( tmp[which.not.reject] ) )
	}
	ans[c("al","ar")]<- ans[c("cil","cir")]<-	c(0,alpha)
	ans["mx.pow"]<- tau.loc			#just for now
	ans["link.mc.sim"]<- 	-w[1]
	ans["link.mc.obs"]<- 	w[1]
	ans["rho.mc"]<- 		w[1] - 0.5
	
	if(verbose)	cat(paste(paste("\n{",args,"<-list(sim.median=",median(sim) ,", obs.median=",median(obs) ,", ZL=",tmp[1] ,", ZU=",tmp[2] ,", CU=",qnorm( alpha ) ,", alpha=",alpha,", tau.cdf=",tau,", tau.loc=",tau.loc,", ",sep=''),paste(names(ans), ans, collapse=', ', sep='='),")}",sep=''))
	if(plot)
	{
		breaks<-  range(c(sim,obs))
		breaks[1]<- breaks[1]*ifelse(breaks[1]<0, 1.2, 0.8)
		breaks[2]<- breaks[2]*ifelse(breaks[2]>0, 1.2, 0.8)
		breaks<- seq(from= breaks[1], to= breaks[2], by= (breaks[2]-breaks[1])/nbreaks)
		xlim<- range(c(sim,obs))
		sim.h<- hist(sim,	breaks= breaks,	plot= F)
		obs.h<- hist(obs,	breaks= breaks,	plot= F)
		ylim<- c(0,max(c(obs.h$intensities, sim.h$intensities),na.rm=TRUE)*1.1)
		plot(1,1,xlab=xlab, xlim= xlim,type='n',ylim=ylim,ylab="probability",main="")
		plot(sim.h,freq=F,add=TRUE,col=myFadeCol("#0080FFFF",0.6),border=myFadeCol("#0080FFFF",0.6))
		plot(obs.h,freq=F,add=TRUE)
		points( sim, rep(ylim[2],length(sim)),col=myFadeCol("#0080FFFF",0.6),pch=20,cex=2.5 )
		points( obs, rep(ylim[2],length(obs)),pch=2,cex=1.5 )
	}
	ans
}
#------------------------------------------------------------------------------------------------------------------------
nabc.mutost.pow<- function(rho, df, tau.u, s.of.T, alpha, rtn.fun= FALSE)
{ 
	x<-	rho
	if(length(x)<10)
		x<- seq(-2*tau.u,2*tau.u,length.out=1e3)
	if(length(rho)<10 && any(rho>2*tau.u))	stop("unexpected rho")
	
	ncp	<- x/s.of.T	
	#may not be a stable numerical approximation
	suppressWarnings({
		tmp	<- pt( tau.u/s.of.T + qt(alpha,df) ,	df, ncp	) - pt( -tau.u/s.of.T - qt(alpha,df) , df, ncp)
	})
	tmp[tmp<=0]<- 0
	#print(c(tau.u,s.of.T,rho)); print(x); print(tmp)
	#fix numerical approximation	
	if(length(which(tmp!=0))<2)
		pw.fun<- approxfun(x=x, y=rep(0,length(x)), method="linear", 0, 0, rule=2 )
	else
		pw.fun<- approxfun(x=x[which(tmp!=0)], y=tmp[which(tmp!=0)], method="linear", 0, 0, rule=2 )
	if(rtn.fun)
		return( pw.fun )
	else
		return( pw.fun(rho) )
}
#------------------------------------------------------------------------------------------------------------------------
#determine both tau.low and tau.up sth max power at rho.star is as given
#this may not work here because the tost is not unbiased
nabc.mutost.onesample.tau.lowup<- function(mx.pw, df, s.of.T, tau.up.ub, alpha, rho.star=0, tol= 1e-5, max.it=100)
{
	curr.mx.pw	<- 0
	tau.up.ub	<- tau.up.ub/2
	tmp			<- max.it
	while(curr.mx.pw<mx.pw && tmp>0)
	{
		tmp			<- tmp-1
		tau.up.ub	<- 2*tau.up.ub
		curr.mx.pw	<- nabc.mutost.pow(rho.star, df, tau.up.ub, s.of.T, alpha)		
	}
	if(tmp==0)	stop("nabc.mutost.onesample.tau.lowup: could not find tau.up.ub")	
#print(tau.up.ub); stop()
	tau.up.lb	<- 0
	error		<- 1	
	while(abs(error)>tol && round(tau.up.lb,d=10)!=round(tau.up.ub,d=10) && max.it>0)
	{
		max.it		<- max.it-1
		tau.up		<- (tau.up.lb + tau.up.ub)/2
		curr.mx.pw	<- nabc.mutost.pow(rho.star, df, tau.up, s.of.T, alpha) 
		error		<- curr.mx.pw - mx.pw
#print(c(curr.mx.pw, tau.up, tau.up.lb, tau.up.ub, max.it))
		if(error<0)
			tau.up.lb<- tau.up
		else
			tau.up.ub<- tau.up
#print(c(abs(error), round(tau.up.lb,d=10)!=round(tau.up.ub,d=10)) )	
	}
	if(max.it==0)	warning("nabc.mutost.onesample.tau.lowup: reached max.it")
	c(-tau.up,tau.up,curr.mx.pw,abs(error))
}
#------------------------------------------------------------------------------------------------------------------------
nabc.mutost.onesample.n.of.y<- function(n.of.x, s.of.Sx, mx.pw, s.of.y, alpha, tau.u.ub=2, tol= 1e-5, max.it=100)
{
	
	s2.of.Sx	<- s.of.Sx*s.of.Sx
	pw.cvar		<- 2*s2.of.Sx
	curr.mx.pw	<- 0
	yn.ub		<- round( n.of.x/2 )
	curr.it		<- max.it
	tau.u		<- tau.u.ub
	while(pw.cvar>s2.of.Sx && curr.it>0)
	{
		curr.it	<- curr.it-1
		yn.ub	<- 2*yn.ub
#print(c(mx.pw,yn.ub))		
		tmp		<- nabc.mutost.onesample.tau.lowup(mx.pw, yn.ub-1, s.of.y/sqrt(yn.ub), 2*tau.u, alpha)
#print(c("OK ",tmp))		
		tau.u	<- tmp[2]
		pw.cmx	<- tmp[3]
		rho		<- seq(-2*tau.u,2*tau.u,length.out=1e3)
		pw.fun	<- nabc.mutost.pow(rho, yn.ub-1, tau.u, s.of.y/sqrt(yn.ub), alpha, rtn.fun=1)
		pw.cvar	<- sum(rho*rho*pw.fun(rho)) / sum(pw.fun(rho))			#mean is 0
#print(c(yn.ub, tau.u, pw.cvar, s2.of.Sx, pw.cmx ))		
	}
	if(curr.it==0)	stop("nabc.mutost.onesample.n.of.y: could not find upper bound for yn")	
#print(c(yn.ub, tau.u, pw.cvar, s2.of.Sx, pw.cmx ))
	yn.lb	<- n.of.x
	error	<- 1	
	while(abs(error)>tol && (yn.lb+1)!=yn.ub && max.it>0)
	{
		max.it	<- max.it-1
		yn		<- round( (yn.lb + yn.ub)/2 )
		tmp		<- nabc.mutost.onesample.tau.lowup(mx.pw, yn-1, s.of.y/sqrt(yn), 2*tau.u, alpha)
		tau.u	<- tmp[2]
		pw.cmx	<- tmp[3]
		rho		<- seq(-2*tau.u,2*tau.u,length.out=1e3)		
		pw.fun	<- nabc.mutost.pow(rho, yn-1, tau.u, s.of.y/sqrt(yn), alpha, rtn.fun=1)
		pw.cvar	<- sum(rho*rho*pw.fun(rho)) / sum(pw.fun(rho))		 
		error	<- pw.cvar - s2.of.Sx
#print(c(pw.cvar, yn, yn.lb, yn.ub, max.it))
		if(error<0)
			yn.ub<- yn
		else
			yn.lb<- yn
#print(c(abs(error), (yn.lb+1)!=yn.ub) )	
	}
	c(yn,-tau.u,tau.u,pw.cvar,pw.cmx,abs(error))
}
#------------------------------------------------------------------------------------------------------------------------
nabc.generic.tost<- function(tost.args, tau.l, tau.u, alpha, tost.distr="t")
{
	ans<- numeric(7)
	names(ans)<- c("error","p.error","lkl","cl","cu","ass.alpha","ass.pval")
	if(!tost.distr%in%c("t"))	stop("unexpected tost.distr")
	
	which.not.reject<- which( c(pt( tost.args[1], tost.args[4] )<1-alpha, pt( tost.args[2], tost.args[4] )>alpha))
	if(length(which.not.reject)%in%c(0,2))		#both upper and lower test statistics indicate mean difference < -tau and mean difference > tau  	OR 		mean difference >= -tau and mean difference <= tau
	{
		#figure out which one is closer to boundary, and schedule that this one is reported
		which.not.reject<- which.min(abs(c(1-alpha-pt( tost.args[1], tost.args[4] ), pt( tost.args[2], tost.args[4] )-alpha )))
	}
	
	ans["error"]	<- tost.args[3]
	ans["p.error"]	<- ifelse(	which.not.reject==1, 		
								1-pt( tost.args[which.not.reject], tost.args[4] ),		
								pt( tost.args[which.not.reject], tost.args[4] )		)
	ans["lkl"]		<- dt( tost.args[3], tost.args[4])		
	ans["cl"]		<- min(tau.l/tost.args[5]+qt( 1-alpha, tost.args[4] ),0)
	ans["cu"]		<- max(0, tau.u/tost.args[5]+qt( alpha, tost.args[4] ))
	ans["ass.alpha"]<- 1-diff(pt( ans[c("cl","cu")], tost.args[4]))											#the corresponding alpha quantile of the traditional t-test with acceptance region [c-,c+] and above s, df
	ans["ass.pval"]	<- ( pt( tost.args[3], tost.args[4] ) - ans["ass.alpha"]/2 ) / ( 1 - ans["ass.alpha"] )						#rescaled p-value that is expected to follow U(0,1) under the point null hypothesis
	
	#POWER[[length(POWER)+1]]<<- c(tau.l, tau.u, tmp[6], sum(moments[,1]), tmp[4], tmp[5] )	#PowerTOST:::.power.TOST(alpha=alpha, tau.l, tau.u, seq(tau.l, tau.u, length.out= 1e3), tmp[6], sum(moments[,1]), tmp[4], bk = 4)	
	ans	
}
#------------------------------------------------------------------------------------------------------------------------
nabc.mutost.onesample<- function(sim, obs.mean, std.sd=NA, args= NA, verbose= FALSE, alpha= 0, tau.u= 0, tau.l= -tau.u, mx.pw=0.9, annealing=1, plot= FALSE, xlab= NA, nbreaks= 40, normal.test= "sf.test")
{
	#verbose<- 1
	ans<- NABC.DEFAULT.ANS
	#compute two sample t-test on either z-scores or untransformed data points
	if(any(is.na(sim)))			stop("nabc.mutost: error at 1a")
	if(any(is.na(obs.mean)))	stop("nabc.mutost: error at 1b")
	if(!is.na(args))
	{
		args<- strsplit(args,'/')[[1]]
		if(length(args)==4)
		{
			standardize	<- as.numeric( args[2] )						
			tau.u		<- as.numeric( args[3] )
			tau.l		<- -tau.u
			alpha		<- as.numeric( args[4] )
		}
		else if(length(args)==5)
		{
			standardize	<- as.numeric( args[2] )
			if(standardize==2)
			{
				annealing<- as.numeric( args[3] )
				mx.pw	<- as.numeric( args[4] )
			}
			else
			{
				tau.l	<- as.numeric( args[3] )
				tau.u	<- as.numeric( args[4] )
			}
			alpha		<- as.numeric( args[5] )
		}
		else
			stop("nabc.mutost: error at 1c")
		args<- args[1]
	}
	if(!standardize%in%c(0,1,2))	stop("incorrect standardize")
	if(standardize==1 && is.na(std.sd)) stop("std.sd not specified but required")	
	if(alpha<0 || alpha>1)		stop("incorrect alpha")
	if(tau.u<0 || tau.l>0)		stop("incorrect tau.u or tau.l")
	if(annealing<1)				stop("incorrect annealing parameter")
#print(standardize)
#standardize<- 0
	ans["pfam.pval"]<-	nabc.get.pfam.pval(sim,normal.test)	
	if(!any(diff(sim)>0))	return(ans)
	sim.mean	<- mean(sim)
	if(standardize==1)
		sim		<- (sim-sim.mean)/sd(sim)*std.sd+sim.mean
	sim.n		<- length(sim)	
	sim.sd		<- sd(sim)	
	if(standardize==2)
	{	
		#print(sim); 
		#print(c(mx.pw,sim.n,sim.sd,alpha))
		tmp		<- nabc.mutost.onesample.tau.lowup(mx.pw, sim.n-1, sim.sd/sqrt(sim.n), 2, alpha)
		if(tmp[4]>0.09)	stop("nabc.mutost: tau.up not accurate")		
		tau.l	<- tmp[1]*annealing
		tau.u	<- tmp[2]*annealing				
		#print(c(annealing,mx.pw,tau.l,tau.u))
		#rho<- seq(tau.l,tau.u,length.out=1e3); y<- nabc.mutost.pow(rho, sim.n-1, tau.u, sim.sd/sqrt(sim.n), alpha); plot(rho,y,type='l')		
	}
	tmp			<- c(	sqrt(sim.n)*(sim.mean-obs.mean-tau.l) / sim.sd,			#[1]	T-	test statistic for -tau (lower test); estimate of the common std dev is simply the std dev in the sample whose sample size is > 1
						sqrt(sim.n)*(sim.mean-obs.mean-tau.u) / sim.sd,			#[2]	T+	test statistic for tau (upper test); estimate of the common std dev is simply the std dev in the sample whose sample size is > 1
						sqrt(sim.n)*(sim.mean-obs.mean) / sim.sd,				#[3]	T	test statistic for equality; estimate of the common std dev is simply the std dev in the sample whose sample size is > 1
						sim.n-1,												#[4] 	degrees of freedom
						sim.sd/sqrt(sim.n),										#[5]	estimate of the std dev of the test statistic is simply the std dev in the sample whose sample size is > 1 divided by that sample size
						sim.sd )												#[6]  	standard deviation of the sample
	tost.ans	<-	nabc.generic.tost(tmp, tau.l, tau.u, alpha, tost.distr="t")
	#print(tost.ans)
	ans[c("error","cil","cir")]<- c(tost.ans["p.error"], 0, alpha)
	ans[c("lkl","pval")]<-  tost.ans[c("lkl","ass.pval")]
	ans[c("al","ar")]	<- 	c(0,alpha)								
	ans["mx.pow"]		<-	nabc.mutost.pow(0, tmp[4], tau.u, tmp[5], alpha) 				
	ans["link.mc.sim"]	<- 	sim.mean
	ans["link.mc.obs"]	<- 	obs.mean
	ans["rho.mc"]		<- 	sim.mean - obs.mean
	ans
}
#------------------------------------------------------------------------------------------------------------------------
nabc.mutost<- function(sim, obs, args= NA, verbose= FALSE, alpha= 0, tau.u= 0, tau.l= -tau.u, plot= FALSE, xlab= NA, nbreaks= 40, normal.test= "sf.test")
{
	#verbose<- 1
	options(warn=1)
	require(PowerTOST)
	options(warn=2)
	ans<- NABC.DEFAULT.ANS
	#compute two sample t-test on either z-scores or untransformed data points
	if(any(is.na(sim)))	stop("nabc.mutost: error at 1a")
	if(any(is.na(obs)))	stop("nabc.mutost: error at 1b")
	#if not missing, arguments 'args' always overwrite 'alpha' and 'nc'
	#expect args string of the form studenttXXX/beta/non-centrality
	if(!is.na(args))
	{
		args<- strsplit(args,'/')[[1]]
		if(length(args)==4)
		{
			standardize<- as.numeric( args[2] )
			tau.u<- as.numeric( args[3] )
			tau.l<- -as.numeric( args[3] )
			alpha<- as.numeric( args[4] )
		}
		else if(length(args)==5)
		{
			standardize<- as.numeric( args[2] )
			tau.l<-  as.numeric( args[3] )
			tau.u<- as.numeric( args[4] )
			alpha<- as.numeric( args[5] )
		}
		else
			stop("nabc.mutost: error at 1c")
		args<- args[1]
	}
#print(standardize)
	if(!standardize%in%c(0,1))	stop("nabc.mutost: error at 1e")
	if(alpha<0 || alpha>1)		stop("nabc.mutost: error at 1f")
	if(tau.u<0 || tau.l>0)		stop("nabc.mutost: error at 1g")
#print(standardize)
#standardize<- 0
	if(length(sim)>3)		
		tmp<- sim
	else 
		tmp<- obs
	ans["pfam.pval"]<-	nabc.get.pfam.pval(tmp,normal.test)
				
	if(!any(diff(sim)>0))
	{	
		length(sim)<- 1
		if(length(obs)<2)	
			return(ans)
	}
	moments<- matrix(NA, 2, 3)	#store number of samples, mean and var in columns
	moments[, 1]<- c( length(sim), length(obs) )	
	if(standardize==1 && moments[1,1]>2 && moments[2,1]>2)
	{	
		sim<- mean(sim)
		moments[1,1]<- 1
		#moments[, 3]<- c(sd(sim), sd(obs))										#second attempt did not work - too much information lost
		#sim<- sim * moments[2,3]/moments[1,3]		
		##sim<- sim / ifelse(moments[1,1]>2,	moments[1,3], moments[2,3])		#first attempt did not work - essentially testing mu/sigma which is more difficult
		##obs<- obs / ifelse(moments[2,1]>2,	moments[2,3], moments[1,3])
	}
	if(all(moments[,1]<2))	stop("nabc.mutost: error at 3a")
	#if(nc!=0)		stop("get.dist.studentt: non-centrality not yet implemented")
	#convert tau from log(obs/sim)<tau to obs-sim<tau	
	moments[,2]<- c( mean(sim), mean(obs) )
	moments[,3]<- c( var(sim), var(obs) )
	if(moments[1,1]>2 && moments[2,1]>2 && standardize==0)						#use unequal sample size, unequal variance
	{
		tmp<- moments[,3]/moments[,1]												#temporarily store 		var(sim)/sim.n, var(obs)/obs.n
		tmp<- c(	(diff( moments[, 2] )-tau.l) / sqrt( sum(tmp) ),				#[1] t test statistic	 for -tau (lower test)
						(diff( moments[, 2] )-tau.u) / sqrt( sum(tmp) ),			#[2] t test statistic	 for tau (upper test)
						diff( moments[, 2] ) / sqrt( sum(tmp) ),					#[3] t test statistic	 for equality
						sum(tmp)^2 / (		tmp[2]^2/(moments[2,1]-1)	+ tmp[1]^2/(moments[1,1]-1)		),				#[4] degrees of freedom: Welch Satterthwaite approximation
						sqrt( sum(tmp) ),											#[5] estimate of standard deviation of the test statistic
						sqrt(sum(moments[,3]*(moments[,1]-1))/(sum(moments[,1])-2))	#[6] common std dev, assuming equal variance. only for power calculations
						)
	}
	else if(moments[1,1]>2 && moments[2,1]>2 && standardize==1)					#estimate common std dev as usual, then use unequal sample size, equal variance
	{
		tmp<- sqrt(sum(moments[,3]*(moments[,1]-1))/(sum(moments[,1])-2)*sum(1/moments[,1]))		#std dev of test statistic
		tmp<- c(	(diff( moments[, 2] )-tau.l) / tmp,								#[1] t test statistic	 for -tau (lower test)
						(diff( moments[, 2] )-tau.u) / tmp,							#[2] t test statistic	 for tau (upper test)
						diff( moments[, 2] ) / tmp,									#[3] t test statistic	 for equality
						sum(moments[,1])-2,											#[4] degrees of freedom
						tmp,														#[5] estimate of standard deviation of the test statistic
						sqrt(sum(moments[,3]*(moments[,1]-1))/(sum(moments[,1])-2))	#[6] estimate of standard deviation of the pooled sample
		)
	}
	else if(moments[1,1]==1 || moments[2,1]==1)					#set common std dev to the std dev in the sample whose sample size is > 1, and use unequal sample size, pooled variance
	{															#this is the same, if standardized or not
		tmp<- !is.na(moments[,3])
		tmp<- c(	(diff( moments[, 2] )-tau.l) / sqrt( moments[tmp,3] / moments[tmp,1] ),			#[1]	test statistic for -tau (lower test); estimate of the common std dev is simply the std dev in the sample whose sample size is > 1
						(diff( moments[, 2] )-tau.u) / sqrt( moments[tmp,3] / moments[tmp,1] ),		#[2]	test statistic for tau (upper test); estimate of the common std dev is simply the std dev in the sample whose sample size is > 1
						diff( moments[, 2] ) / sqrt( moments[tmp,3] / moments[tmp,1] ),				#[3]	test statistic for equality; estimate of the common std dev is simply the std dev in the sample whose sample size is > 1
						moments[tmp,1]-1,															#[4] 	degrees of freedom
						sqrt( moments[tmp,3]/moments[tmp,1] ),										#[5]	estimate of the std dev of the test statistic is simply the std dev in the sample whose sample size is > 1 divided by that sample size
						sqrt( moments[tmp,3] )														#[6]  	standard deviation of the sample
						)
		args<- paste(args,".equalvariance",sep='')
	}
	else	stop("nabc.mutost: error at 3b")
	#cat( paste( "\nloc.cdf= ",pnorm( tau / ( tmp[6] * sqrt(2) ) ) - 0.5 ) )

	which.not.reject<- which( c(pt( tmp[1], tmp[4] )<1-alpha, pt( tmp[2], tmp[4] )>alpha))
	if(length(which.not.reject)%in%c(0,2))		#both upper and lower test statistics indicate mean difference < -tau and mean difference > tau  	OR 		mean difference >= -tau and mean difference <= tau
	{
		#figure out which one is closer to boundary, and schedule that this one is reported
		which.not.reject<- which.min(abs(c(1-alpha-pt( tmp[1], tmp[4] ), pt( tmp[2], tmp[4] )-alpha )))
	}
	
	
	#POWER[[length(POWER)+1]]<<- c(tau.l, tau.u, tmp[6], sum(moments[,1]), tmp[4], tmp[5] )	#PowerTOST:::.power.TOST(alpha=alpha, tau.l, tau.u, seq(tau.l, tau.u, length.out= 1e3), tmp[6], sum(moments[,1]), tmp[4], bk = 4)	
	
	#print(which.not.reject)
	ans[c("lkl","pval","error")]<- c( 		dt( tmp[3], tmp[4]),
											pt( tmp[3], tmp[4] ),
											ifelse(which.not.reject==1, 		1-pt( tmp[which.not.reject], tmp[4] ),		pt( tmp[which.not.reject], tmp[4] ) )				)
	ans[c("al","ar")]	<- 	c(min(tau.l/tmp[5]+qt( 1-alpha, tmp[4] ),0), max(0, tau.u/tmp[5]+qt( alpha, tmp[4] ))		)		#'ar' is upper boundary point c+ of TOST for the traditional two-sided t-test statistic T= (bar(x)-bar(y)) / s, and likewise 'al' for c-								
	ans["mx.pow"]		<-	1-diff(pt( ans[c("al","ar")], tmp[4]))	#the corresponding alpha quantile of the traditional t-test with acceptance region [c-,c+] and above s, df						
	ans["pval"]			<-	( ans["pval"] - ans["mx.pow"]/2 ) / ( 1 - ans["mx.pow"] )							#rescaled p-value that is expected to follow U(0,1) under the point null hypothesis
	ans[c("cil","cir")]	<-	c(0, alpha)
	ans["mx.pow"]		<-	nabc.mutost.pow(0, tmp[4], tau.u, tmp[5], alpha) 				
	ans["link.mc.sim"]	<- 	moments[1,2]
	ans["link.mc.obs"]	<- 	moments[2,2]
	ans["rho.mc"]		<- 	diff(moments[,2]) / tmp[5]
	#ans["al"]<- tmp[3]

	if(verbose)	cat(paste(paste("\n{",args,"<-list(sim.mean=",moments[1, 2] ,", obs.mean=",moments[2,2] ,", sim.n=",moments[1, 1] ,", obs.n=",moments[2,1] ,", per.capita.s=",tmp[6],", df=",tmp[4],", alpha=",alpha,", tau.l=",tau.l,", tau.u=",tau.u,", tl=",tmp[1],", tu=",tmp[2],", cil=", tau.l/tmp[5]+qt( 1-alpha, tmp[4] ),", ciu=", tau.u/tmp[5]-qt( 1-alpha, tmp[4] ),", ",sep=''),paste(names(ans), ans, collapse=', ', sep='='),")}",sep=''))
	if(plot)
	{
		breaks<-  range(c(sim,obs))
		breaks[1]<- breaks[1]*ifelse(breaks[1]<0, 1.2, 0.8)
		breaks[2]<- breaks[2]*ifelse(breaks[2]>0, 1.2, 0.8)
		breaks<- seq(from= breaks[1], to= breaks[2], by= (breaks[2]-breaks[1])/nbreaks)
		xlim<- range(c(sim,obs))
		if(moments[1,1]>2 && moments[2,1]>2)
		{
			sim.h<- hist(sim,	breaks= breaks,	plot= F)
			obs.h<- hist(obs,	breaks= breaks,	plot= F)
		}
		else
		{
			if(which(!is.na(moments[,3]))==1)			#1 if sim has > 1 sample size, 2 if obs has > 1 sample size
				sim.h<-obs.h<- hist(sim,	breaks= breaks,	plot= F)
			else
				sim.h<-obs.h<- hist(obs,	breaks= breaks,	plot= F)
		}
		ylim<- c(0,max(c(obs.h$intensities, sim.h$intensities),na.rm=TRUE)*1.1)
		plot(1,1,xlab=xlab, xlim= xlim,type='n',ylim=ylim,ylab="probability",main="")
		if(moments[1,1]>2 && moments[2,1]>2)
		{
			plot(sim.h,freq=F,add=TRUE,col=myFadeCol("#0080FFFF",0.6),border=myFadeCol("#0080FFFF",0.6))
			lines(seq( xlim[1], xlim[2], by= diff(xlim)/100 ), dnorm( 	seq( xlim[1], xlim[2], by= diff(xlim)/100 ), moments[1,2],  sqrt(moments[1,3])), col=myFadeCol("#0080FFFF",0.6) )
			points( sim, rep(ylim[2],length(sim)),col=myFadeCol("#0080FFFF",0.6),pch=20,cex=2.5 )
			points( obs, rep(ylim[2],length(obs)),pch=2,cex=1.5 )
		}
		plot(obs.h,freq=F,add=TRUE)
		if(moments[1,1]>2 && moments[2,1]>2)
			lines(seq( xlim[1], xlim[2], by= diff(xlim)/100 ), dnorm( 	seq( xlim[1], xlim[2], by= diff(xlim)/100 ), moments[2,2],  sqrt(moments[2,3])) )
		else
		{
			lines(seq( xlim[1], xlim[2], by= diff(xlim)/100 ), dnorm( 	seq( xlim[1], xlim[2], by= diff(xlim)/100 ), moments[!is.na(moments[,3]),2],  tmp[5]) )
			abline(v= moments[is.na(moments[,3]),2], col=myFadeCol("#0080FFFF",0.6))
		}
	}
	ans
}
#---------------------------------------------------------------------------------------------------------------------------
nABC.getlevelset.2d<- function(df, lnk.name, theta.names, rho.eq=0, rho.eq.sep=100, rho.eq.q= 0.075, theta.sep= 100, method="quantile",plot=0, verbose=0)
{
	if(!method%in%c("quantile","fixed"))	stop("nABC.getlevelset.2d: error at 1a")
	if(method=="quantile")
	{
		rho.eq.sep.u<- length(df[,lnk.name])
		rho.eq.sep.l<- 1		
	}
	require(locfit)
	tmp<- paste("locfit(",lnk.name,'~',paste(theta.names,collapse=':',sep=''),", data=df,maxk=200)",sep='')
	lnk.locfit<- eval(parse(text=tmp))
	lnk.xrange<- lfmarg(lnk.locfit$box, rep(ifelse(lnk.locfit$mi["d"]==1, 100, theta.sep), lnk.locfit$mi["d"]))
	names(lnk.xrange)<- theta.names
	lnk.pred<- locfit:::preplot.locfit(lnk.locfit, lnk.xrange, band = "none", tr = NULL, what = "coef", get.data = 0, f3d = 0)
	lnk.len<- sapply(lnk.xrange,length)
	lnk.pred<- array(lnk.pred$fit, dim=lnk.len)
	lnk.xsep<- sapply(lnk.xrange, function(x)	diff(x[1:2])	)
	if(method=="quantile")
	{	
		while(1)
		{		
			rho.eq.sep<- mean(c(rho.eq.sep.l,rho.eq.sep.u))
			rho.eq.eps<- range(df[,lnk.name]) /  rho.eq.sep
			lnk.eqidx<- which( lnk.pred<(rho.eq+rho.eq.eps) & lnk.pred>(rho.eq-rho.eq.eps) )
			if(verbose) cat(paste("\nrho.eq.seps are",rho.eq.sep.l,rho.eq.sep.u,"\tlevel set length is",length(lnk.eqidx)))
			if(	abs(length(lnk.eqidx)/length(df[,lnk.name])-rho.eq.q)<rho.eq.q/20 ||
					rho.eq.sep.l==rho.eq.sep.u	)
				break
			if(length(lnk.eqidx)/length(df[,lnk.name])<rho.eq.q)
				rho.eq.sep.u<- rho.eq.sep
			else
				rho.eq.sep.l<- rho.eq.sep
		}
	}
	else
	{
		rho.eq.eps<- range(df[,lnk.name]) / rho.eq.sep
		lnk.eqidx<- which( lnk.pred<(rho.eq+rho.eq.eps) & lnk.pred>(rho.eq-rho.eq.eps) )
	}	
	#print(lnk.locfit)
	#print(lnk.len)
	#print(lnk.xsep)
	#print(dim(lnk.pred))
	#print(lnk.eqidx)	
	lnk.eqx<- matrix(NA,nrow=length(lnk.len),ncol=length(lnk.eqidx),dimnames=list(names(lnk.xrange),c()))
	for(i in rev(seq_along(lnk.len)))
	{
		lnk.eqx[i,]<- (lnk.eqidx-1) %/% prod(lnk.len[-length(lnk.len)]) + 1		#idx wrt ith theta
		lnk.eqidx<- (lnk.eqidx-1) %% prod(lnk.len[-length(lnk.len)]) + 1
		lnk.len<- lnk.len[-length(lnk.len)]
	}
	lnk.eqx<- sapply(seq_len(ncol(lnk.eqx)),function(j)
			{
				sapply(seq_len(nrow(lnk.eqx)),function(i)		lnk.xrange[[i]][ lnk.eqx[i,j] ]	)
			}) 	
	rownames(lnk.eqx)<- names(lnk.xrange)
	#print(lnk.eqx)
	lnk.eqx	
}
#---------------------------------------------------------------------------------------------------------------------------
nABC.getlevelsetintersection.2d<- function(lsets,theta.names,rho.eq.sep= 25,plot=0,verbose=0)
{
	theta.range<- sapply(seq_len(nrow(lsets[[1]])),function(d)
			{
				range(c(sapply(seq_along(lsets),function(k){		range(lsets[[k]][d,])	}),recursive=1))		
			})
	colnames(theta.range)<- theta.names
	theta.sep<- apply(theta.range,2,diff)/rho.eq.sep 
	
	theta.intersection<- lsets[[1]]
	lsets.intersection<- vector("list",length(lsets)-1)
	for(k in seq_along(lsets)[-1])	
	{
		if(verbose)	cat(paste("\nnumber of theta in intersection",ncol(theta.intersection),"\nnow intersect w summary",names(lsets)[k],"\n"))		
		tmp<- .Call("abcIntersectLevelSets", lsets[[k]], theta.intersection, theta.sep)				
		tmp<- which(tmp<sum(theta.sep*theta.sep))				
		theta.close<- matrix(	c( (tmp-1) %/% ncol(lsets[[k]]) + 1, (tmp-1) %% ncol(lsets[[k]]) + 1 ), nrow=2,ncol=length(tmp),byrow=1,dimnames=list(c("1","k"),c()) )		
		tmp<- unique(theta.close["1",])
		lsets.intersection[[k-1]]<- theta.intersection[,unique(theta.close["1",]),drop=0]				
		if(!length(tmp))
			break
		theta.intersection<- theta.intersection[,tmp,drop=0]
	}	
	if(plot)
	{
		plot(	theta.intersection[1,],theta.intersection[2,],
				xlim=range(theta.intersection[1,]),ylim=range(theta.intersection[2,]),xlab=rownames(theta.intersection)[1],ylab=rownames(theta.intersection)[2],
				type='p',pch=19,col=myFadeCol("black",0.3))
	}
	if(verbose) cat(paste("\nfinal number of theta in intersection",ncol(theta.intersection),"\n"))
	theta.intersection
}
#---------------------------------------------------------------------------------------------------------------------------
nABC.getlevelset.3d<- function(df, lnk.name, theta.names, rho.eq=0, rho.eq.sep=100, rho.eq.q= 0.075, theta.sep= 100, method="quantile",plot=0)
{	
	if(!method%in%c("quantile","fixed"))	stop("nABC.getlevelset.3d: error at 1a")
	if(method=="quantile")
	{
		rho.eq.sep.u<- length(df[,lnk.name])
		rho.eq.sep.l<- 1		
	}			
	tmp<- paste("locfit(",lnk.name,'~',paste(theta.names,collapse=':',sep=''),", data=df)",sep='')								
	lnk.locfit<- eval(parse(text=tmp))					
	lnk.xrange<- lfmarg(lnk.locfit$box, rep(ifelse(lnk.locfit$mi["d"]==1, 100, theta.sep), lnk.locfit$mi["d"]))
	names(lnk.xrange)<- theta.names
	lnk.pred<- locfit:::preplot.locfit(lnk.locfit, lnk.xrange, band = "none", tr = NULL, what = "coef", get.data = 0, f3d = 0)
	lnk.len<- sapply(lnk.xrange,length)
	lnk.pred<- array(lnk.pred$fit, dim=lnk.len)
	lnk.xsep<- sapply(lnk.xrange, function(x)	diff(x[1:2])	)
	
	if(method=="quantile")
	{	
		while(1)
		{		
			rho.eq.sep<- mean(c(rho.eq.sep.l,rho.eq.sep.u))
			rho.eq.eps<- range(df[,lnk.name]) /  rho.eq.sep
			lnk.eqidx<- which( lnk.pred<(rho.eq+rho.eq.eps) & lnk.pred>(rho.eq-rho.eq.eps) )
			cat(paste("\nrho.eq.seps are",rho.eq.sep.l,rho.eq.sep.u,"\tlevel set length is",length(lnk.eqidx)))
			if(	abs(length(lnk.eqidx)/length(df[,lnk.name])-rho.eq.q)<rho.eq.q/20 ||
					rho.eq.sep.l==rho.eq.sep.u	)
				break
			if(length(lnk.eqidx)/length(df[,lnk.name])<rho.eq.q)
				rho.eq.sep.u<- rho.eq.sep
			else
				rho.eq.sep.l<- rho.eq.sep
		}
	}
	else
	{
		rho.eq.eps<- range(df[,lnk.name]) / rho.eq.sep
		lnk.eqidx<- which( lnk.pred<(rho.eq+rho.eq.eps) & lnk.pred>(rho.eq-rho.eq.eps) )
	}
	lnk.eqx<- matrix(NA,nrow=length(lnk.len),ncol=length(lnk.eqidx),dimnames=list(names(lnk.xrange),c()))
	for(i in rev(seq_along(lnk.len)))
	{
		lnk.eqx[i,]<- (lnk.eqidx-1) %/% prod(lnk.len[-length(lnk.len)]) + 1		#idx wrt ith theta
		lnk.eqidx<- (lnk.eqidx-1) %% prod(lnk.len[-length(lnk.len)]) + 1
		lnk.len<- lnk.len[-length(lnk.len)]
	}
	lnk.eqx<- sapply(seq_len(ncol(lnk.eqx)),function(j)
			{
				sapply(seq_len(nrow(lnk.eqx)),function(i)		lnk.xrange[[i]][ lnk.eqx[i,j] ]	)
			}) 	
	rownames(lnk.eqx)<- names(lnk.xrange)
	if(plot)
	{
		ABC.CI.MMCMC.plot.trellis.levelset(lnk.locfit, zlab=theta.names[length(theta.names)], ysep=rho.eq.eps)		
	}
	lnk.eqx
}