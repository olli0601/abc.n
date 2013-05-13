
#------------------------------------------------------------------------------------------------------------------------
#' Compute power of the asymptotic equivalence test for autocorrelations at lag 1
#' @export 
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
#' @export
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
#' @export
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
#' @export
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
#' 	nabc.acf.equivalence.tau.lowup(0.9, 2, floor(sim.n / (1+leave.out)), 0.01)
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
#' @export
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
	ans["error"]				<- ifelse(which.not.reject==1, 		1-pnorm( tmp[which.not.reject] ),		pnorm( tmp[which.not.reject] ) )
	ans[c("cil","cir")]			<- c(0,alpha)		
	ans[c("tl","tr","nsim")]	<- c(tau.l, tau.u, z.sim["n"]) 
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