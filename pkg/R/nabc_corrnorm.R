
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
#' 	pw<- nabc.tosz.pow(rho, tau.u, alpha, 1/sqrt(floor(sim.n/3)-3))
nabc.tosz.pow<- function(rho, tau.u, alpha, sigma, norm = 1, support= c(-Inf,Inf), log=FALSE)
{ 
	stopifnot(sigma>0,norm>0,alpha>0, alpha<0.5, support[1]<=support[2])
	ans 			<- rho
	in_support 		<- (rho >= support[1] & rho <= support[2])
	ans[!in_support]<- ifelse(log,-Inf,0)
	if(any(in_support)) 
	{
		ans[in_support]		<- (pnorm( (tau.u-rho[in_support])/sigma + qnorm(alpha) ) - pnorm( -(tau.u+rho[in_support])/sigma - qnorm(alpha) ))/norm
		tmp					<- which(ans[in_support]<0)
		if(length(tmp))
			ans[tmp]		<- 0
		if(log)
			ans[in_support]	<- log(ans[in_support])
	}
	ans
}
#------------------------------------------------------------------------------------------------------------------------
nabc.tosz.pow.norm<- function(tau.u, sigma, alpha=0.01, support=c(-Inf,Inf))
{
	ans	<- integrate(nabc.tosz.pow, lower=support[1], upper=support[2], tau.u=tau.u, sigma=sigma, alpha=alpha, norm=1, support=support, log=FALSE)
	ans$value
}
#------------------------------------------------------------------------------------------------------------------------
nabc.tosz.sulkl<- function(rho, sigma, norm = 1, support= c(-Inf,Inf), log=FALSE)
{ 
	stopifnot(sigma>0,norm>0,support[1]<=support[2])	
	ans 			<- rho
	in_support 		<- (rho >= support[1] & rho <= support[2])
	ans[!in_support]<- ifelse(log,-Inf,0)
	if(any(in_support)) 
	{
		if(log)
			ans[in_support]	<- dnorm(rho[in_support], 0, sd=sigma, log=T)-log(norm)
		else
			ans[in_support] <- dnorm(rho[in_support], 0, sd=sigma)/norm		
	}	
	ans
}
#------------------------------------------------------------------------------------------------------------------------
nabc.tosz.sulkl.norm<- function(Tsd, support= c(-Inf,Inf))
{
	diff(pnorm(support,mean=0,sd=Tsd))
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
nabc.acf.equivalence.cor<- function(x, leave.out=0, len= ceiling(length(x)/(1+leave.out)) )
{
	tmp	<- rbind( x[-1], x[-length(x)]	)
	tmp	<- tmp[,seq.int(1,ncol(tmp),by=1+leave.out)]
	tmp	<- tmp[,seq_len(len)]
	if(any(is.na(tmp)))		stop("Unexpected NA in tmp")
	tmp2<- cor(tmp[1,],tmp[2,])
	ans	<- c(tmp2,		.5 * log( (1+tmp2)/(1-tmp2) ), 		ncol(tmp)	)
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
#' nabc.tosz.criticalregion(tau.l, tau.u, floor(sim.n / (1+leave.out)), 0.01)
nabc.tosz.criticalregion<- function(tau.u, sigma, alpha=0.01)
{				 
	c(cl=min(-tau.u/sigma+qnorm(1-alpha),0), cu=max(0,tau.u/sigma+qnorm(alpha)))
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
#' 	nabc.tosz.calibrate.tolerances(0.9, 2, 1/sqrt(floor(sim.n / (1+leave.out))-3), 0.01)
nabc.tosz.calibrate.tolerances<- function(mx.pw, tau.up.ub, sigma, alpha=0.01, rho.star=0, tol= 1e-5, max.it=100, pow_scale=2, verbose=0)
{
	stopifnot(mx.pw>0, mx.pw<1, sigma>0, alpha>0, alpha<0.5)
	curr.mx.pw	<- 0
	tau.up.ub	<- tau.up.ub/2
	tmp			<- max.it
	while(curr.mx.pw<mx.pw && tmp>0)
	{
		tmp			<- tmp-1
		tau.up.ub	<- 2*tau.up.ub
		rho			<- seq(-tau.up.ub*pow_scale, tau.up.ub*pow_scale, len=1024)
		pw			<- nabc.tosz.pow(rho, tau.up.ub, alpha, sigma)
		curr.mx.pw	<- max(pw)				 
		if(verbose)	cat(paste("\ntrial upper bound",tau.up.ub,"with power",curr.mx.pw,"at rho=",rho[ which.max(pw) ]))
	}
	if(tmp==0)	stop("could not find tau.up.ub")		
	if(verbose)	cat(paste("\nFound upper bound",tau.up.ub,"with power",curr.mx.pw,"at rho=",rho[ which.max(pw) ]))
	tau.up.lb	<- 0
	error		<- 1	
	while(abs(error)>tol && round(tau.up.lb,d=10)!=round(tau.up.ub,d=10) && max.it>0)
	{
		max.it		<- max.it-1
		tau.up		<- (tau.up.lb + tau.up.ub)/2
		rho			<- seq(-tau.up*pow_scale, tau.up*pow_scale, len=1024)
		pw			<- nabc.tosz.pow(rho, tau.up, alpha, sigma)
		curr.mx.pw	<- max(pw)				 
		error		<- curr.mx.pw - mx.pw
		if(verbose)	cat(paste("\ntrial tau.u",tau.up,"with power",curr.mx.pw,"at rho=",rho[ which.max(pw) ],"search interval",tau.up.lb,tau.up.ub,"error",error))
		if(error<0)
			tau.up.lb	<- tau.up
		else
			tau.up.ub	<- tau.up
#print(c(abs(error), round(tau.up.lb,d=10)!=round(tau.up.ub,d=10)) )	
	}
	if(max.it==0)	warning("reached max.it")
	ans				<- c(-tau.up, tau.up, curr.mx.pw, abs(error))	
	names(ans)		<- c("tau.l", "tau.u", "curr.mx.pw", "error")
	c(ans, nabc.tosz.criticalregion(tau.up, sigma, alpha))
}
#------------------------------------------------------------------------------------------------------------------------
#	mx.pw=0.9; alpha=0.01; pow_scale=2; debug = 0; calibrate.tau.u=T; plot = F; verbose=0
nabc.tosz.calibrate.tolerances.getkl <- function(s.of.x, s.of.T, tau.u, mx.pw=0.9, alpha=0.01, pow_scale=2, debug = 0, calibrate.tau.u=T, plot = F, verbose=0) 
{
	tau.l<- pw.cmx<- error<- cl<- cu<- NA	
	if(calibrate.tau.u) 
	{			
		g(tau.l, tau.u, pw.cmx,	error, cl, cu)	%<-% nabc.tosz.calibrate.tolerances(mx.pw, tau.u, s.of.T, alpha=alpha, pow_scale=pow_scale, verbose=verbose)	#tau.u is taken as upper bound on calibrated tau.u			
		if (abs(pw.cmx - mx.pw) > 0.09) 	stop("tau.up not accurate")			
	}
	#else tau.u is taken as final tau.u
		
	#truncate pow and compute pow_norm	
	rho			<- seq(-tau.u*pow_scale,tau.u*pow_scale, len=1024)	
	support		<- range(rho)
	pw.norm		<- nabc.tosz.pow.norm(tau.u, s.of.T, alpha=alpha, support=support)
	pw			<- nabc.tosz.pow(rho, tau.u, alpha, s.of.T, norm=pw.norm, support=support, log=FALSE)	
	lkl.norm	<- nabc.tosz.sulkl.norm(s.of.x, support=support)	
	lkl			<- nabc.tosz.sulkl(rho, s.of.x, norm=lkl.norm, support=support, log=FALSE)
	
	pw.arg		<- list(tau.u=tau.u, sigma=s.of.T, alpha=alpha, norm=pw.norm, support=support)
	lkl.arg		<- list(sigma=s.of.x, norm=lkl.norm, support=support)
	tmp 		<- integrate(nabc.kl.integrand, lower=support[1], upper=support[2], dP=nabc.tosz.sulkl, dQ=nabc.tosz.pow, P_arg=lkl.arg, Q_arg=pw.arg)	
	if (tmp$message != "OK") 	warning(tmp$message)
	KL_div		<- tmp$value
	
	if (plot) 
	{
		library(ggplot2)
		library(reshape2)
		df_lkl 				<- data.frame(x=rho, yes=lkl, no=lkl*lkl.norm )
		df_lkl$distribution <- "summary likelihood"
		df_pow 				<- data.frame(x=rho, yes=pw, no=pw*pw.norm)
		df_pow$distribution <- "ABC approximation"
		gdf 				<- rbind(df_pow, df_lkl)
		gdf					<- melt(gdf,id.vars=c("x","distribution"))
		p 					<- ggplot(data = gdf, aes(x = x, y = value, colour = distribution, linetype=variable))
		p					<- p+geom_vline(xintercept=c(-tau.u,tau.u),linetype="dotted")
		p					<- p+geom_hline(yintercept= mx.pw,linetype="dotted")
		p					<- p+ geom_line()
		p					<- p+scale_linetype("truncated and\nstandardized?")
		p					<- p+xlab(expression(rho))+ylab("")
		p 					<- p + ggtitle(paste("n.of.y=", 1/s.of.T^2+3,"\ntau.u=", tau.u,"\nKL=", KL_div))
		print(p)
	}
	pw.cmx 	<- ifelse(calibrate.tau.u, pw.cmx, 	nabc.tosz.pow(rho=1, tau.u, s.of.T, alpha, norm=1, support=support))
	cu	 	<- ifelse(calibrate.tau.u, cu, 		nabc.tosz.criticalregion(tau.u, s.of.T, alpha))
	cl	 	<- ifelse(calibrate.tau.u, cl, 		-cu)
	c(KL_div=KL_div, tau.l= -tau.u, tau.u=tau.u, pw.cmx = pw.cmx, cl=cl, cu=cu)	
}
#------------------------------------------------------------------------------------------------------------------------
#n.of.y=n.of.x; n2s= function(n){ 1/sqrt(n-3) }; mx.pw=0.9; alpha=0.01; max.it=100; pow_scale=2; debug=F; plot=F
nabc.tosz.calibrate<- function(n.of.x, n.of.y=n.of.x, n2s= function(n){ 1/sqrt(n-3) }, s2n=function(s){ (1/s)^2+3 }, mx.pw=0.9, alpha=0.01, max.it=100, pow_scale=2, debug=F, plot=F)
{	
	KL.of.yn_ub	<- KL.of.yn		<- error <- curr.mx.pw <- tau.low <- cl <- cu	<- NA
	s.of.x		<- n2s(n.of.x)
	s.of.T		<- n2s(n.of.y)
	#KL for initial n.of.y	
	KL.of.yn		<- nabc.tosz.calibrate.tolerances.getkl(s.of.x, s.of.T, 3*s.of.T, mx.pw=mx.pw, alpha=alpha, pow_scale=pow_scale)["KL_div"]	
	#KL always decreases from n.of.x. Find upper bound yn.ub such that KL first increases again.	
	curr.it 		<- max.it
	yn.ub 			<- 2 * n.of.y	
	KL.of.yn_ub		<- nabc.tosz.calibrate.tolerances.getkl(s.of.x, n2s(yn.ub), 3*n2s(yn.ub), mx.pw=mx.pw, alpha=alpha, pow_scale=pow_scale)["KL_div"]
	
	while (KL.of.yn_ub < KL.of.yn && curr.it > 0) 
	{
		#print(c(yn.ub, KL.of.yn_ub, KL.of.yn, curr.it))
		curr.it 		<- curr.it - 1
		KL.of.yn 		<- KL.of.yn_ub
		yn.ub 			<- 2 * yn.ub
		KL.of.yn_ub		<- nabc.tosz.calibrate.tolerances.getkl(s.of.x, n2s(yn.ub), 3*n2s(yn.ub), mx.pw=mx.pw, alpha=alpha, pow_scale=pow_scale)["KL_div"]				
		if(debug)	cat(paste("\ntrial upper bound m=",yn.ub,"with KL",KL.of.yn_ub))
	}			
	if (curr.it == 0) 	stop("could not find upper bound for yn")					
	if(debug)			cat(paste("\nFound upper bound m=",yn.ub,"with KL",KL.of.yn_ub))
	yn.lb				<- ifelse(curr.it==max.it, yn.ub/2, yn.ub/4)
	if(debug)	cat(paste("\nupper and lower bounds on m:",yn.lb, yn.ub))
	
	
	
	KL_args				<- list(s.of.x=s.of.x, tau.u=3*s.of.x, mx.pw=mx.pw, alpha=alpha, pow_scale=pow_scale, calibrate.tau.u=T, plot=F)	
	tmp 				<- optimize(nabc.kl.optimize, interval = n2s(c(yn.ub, yn.lb)), x_name="s.of.T", is_integer=F, KL_divergence = "nabc.tosz.calibrate.tolerances.getkl", KL_args=KL_args, verbose=debug, tol=1e-5)	
	n.of.y 				<- round(s2n(tmp$minimum))
	
	g(KL_div, tau.l, tau.u, pw.cmx, cl, cu)	%<-%	nabc.tosz.calibrate.tolerances.getkl(s.of.x, n2s(n.of.y), 3*s.of.x, mx.pw=mx.pw, alpha=alpha, pow_scale=pow_scale,calibrate.tau.u=T, plot=plot)
			
	c(n.of.y=n.of.y, tau.l=tau.l, tau.u=tau.u, pw.cmx=pw.cmx, KL_div=KL_div, cl=cl, cu=cu)		
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