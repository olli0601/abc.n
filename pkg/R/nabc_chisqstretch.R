
#' Compute the density of the (possible truncated) power of the equivalence test for population dispersion of normal summary values
#' @param rho vector of quantile
#' @param scale	scaling of T apart from rho, either n-1 for unbiased ABC or n for exact MAP
#' @param df	degrees of freedom
#' @param c.l,c.u lower and upper ABC tolerance
#' @param norm normalization constant for the truncated power function.
#' @param support vector of dimension 2. Support of the truncated power function.
#' @param log logical; if \code{TRUE}, densities d are given as log(d). 
#' @note The summary likelihood can be truncated to \code{support} and then standardized with \code{norm}.
#' For computational efficiency, both \code{norm} and \code{support} must be provided although each one can be derived (numerically) from the other.
#' @export
#' 	
dchisqstretch_pow <- function(rho, scale, df, c.l, c.u, norm=1, support=c(0,Inf), log=FALSE){
	
	ans<-rho
	in_support<- (rho >= support[1] & rho <= support[2])
	ans[!in_support] <- 0
	
	if(any(in_support)){			
		ans[in_support] <- (pchisq( c.u/rho[in_support]*scale, df ) - pchisq( c.l/rho[in_support]*scale, df))/norm
	}
	
	if(log){
		ans<-log(ans)
	}
	
	return(ans)
}
#------------------------------------------------------------------------------------------------------------------------
#' Compute power of the exact equivalence test for dispersion
#' @export 
#' @param rho 	true ratio in simulated variance / observed variance
#' @param scale	scaling of T apart from rho, either n-1 for unbiased ABC or n for exact MAP
#' @param df	degrees of freedom
#' @param cl	lower ABC tolerance
#' @param cu	upper ABC tolerance 
#' @return power of the exact test. this is exact.
#' @examples	alpha<- 0.01						
#'	tau.up<- 1.09
#'	yn<- 5e3
#'	tau.low<- nabc.chisqstretch.tau.low(tau.up, yn-1, alpha)
#'	rej<- .Call("abcScaledChiSq",	c(yn-1,yn-1,tau.low,tau.up,alpha,1e-10,100,0.05)	)
#'	rho<- seq(tau.low,tau.up,by=0.001)
#'	nabc.chisqstretch.pow(rho,yn-1,yn-1,rej[1],rej[2])		
nabc.chisqstretch.pow<- function(rho, scale, df, cl, cu)
{
	pchisq( cu/rho*scale, df ) - pchisq( cl/rho*scale, df)
}
#------------------------------------------------------------------------------------------------------------------------
#' Compute the density of the (possibly truncated) summary likelihood for population dispersion of normal summary values
#' @inheritParams nabc.mutost.kl
#' @param rho vector of quantile
#' @param norm scalar, 0<\code{norm}<=1, normalization constant for the truncated summary likelihood.
#' @param support vector of dimension 2, support of the truncated summary likelihood.
#' @param log logical; if \code{TRUE}, densities d are given as log(d). 
#' @note The summary likelihood can be truncated to \code{support} and then standardized with \code{norm}.
#' For computational efficiency, both \code{norm} and \code{support} must be provided although each one can be derived from the other. 
#' \code{support=qigamma(c(1-norm,1+norm)/2,(n.of.x-2)/2,s.of.x^2*(n.of.x-1)/2)} and \code{norm=diff(pigamma(support,(n.of.x-2)/2,s.of.x^2*(n.of.x-1)/2)}.
#' @import pscl
#' @export	
#'
nabc.chisqstretch.sulkl<- function(rho, n.of.x, s.of.x, norm = 1, support= c(0,Inf), log=FALSE) 
{
	alpha	<- (n.of.x-2)/2	 
	beta	<- s.of.x^2*(n.of.x-1)/2
	if(norm>1) 	stop("norm must be <= 1")
	ans 				<- rho
	in_support 			<- (rho >= support[1] & rho <= support[2])
	ans[!in_support]	<- 0
	if (any(in_support)) 	
		ans[in_support]	<- densigamma(rho[in_support], alpha, beta)/norm	
	if(log)
		ans				<- log(ans)
	return(ans)
}
#------------------------------------------------------------------------------------------------------------------------
#' Compute Kullback-Leibler divergence between the summary likelihood and the power function of chisqstretch 
#' @param n.of.x number of observed summary values
#' @param s.of.x standard deviation of observed summary values
#' @param n.of.y number of simulated summary values
#' @param s.of.y standard deviation of simulated summary values
#' @param mx.pw maximum power at the point of reference (rho.star=0) (only when \code{calibrate.tau.u==TRUE}).
#' @param alpha level of the equivalence test
#' @param calibrate.tau.u if \code{TRUE} the upper tolerance of the equivalence region (\code{tau.u}) is calibrated so that power at the point of reference is equal to \code{mx.pw}
#' @param tau.u	upper tolerance of the equivalence region. If \code{calibrate.tau.u==TRUE}, \code{tau.u} is just a guess on an upper bound on the upper tolerance of the equivalence region.
#' @param for.mle	calibrate so that the mode of the power is at the MLE
#' @param pow_scale scale for the support of the standardized power. The power is truncated between \code{[tau.l/pow_scale,tau.u*pow_scale]} and then standardized.
#' @param debug flag if C implementation is used
#' @param plot whether to plot the two distributions
#' @return	vector of length 6
#' 	\item{KL_div}{the Kullback Leibler divergence}	
#' 	\item{tau.l}{lower tolerance of the equivalence region}	
#' 	\item{tau.u}{upper tolerance of the equivalence region}
#' 	\item{c.l}{lower ABC tolerance}	
#' 	\item{c.u}{upper ABC tolerance}	
#' 	\item{pw.cmx}{actual maximum power associated with the equivalence region}
#' @note Whatever the value of \code{calibrate.tau.u}, the lower tolerance of the equivalence region (\code{tau.l}) is always numerically calibrated using \link{nabc.chisqstretch.tau.low}.
#' @export
#' @import ggplot2 reshape2 pscl
#' @examples
#' 
#' nabc.chisqstretch.kl(n.of.x=60,s.of.x=0.1,n.of.y=60,s.of.y=0.3, mx.pw=0.9,
#' alpha=0.01, calibrate.tau.u=T, tau.u=1, plot=T)
#'
nabc.chisqstretch.kl <- function(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = F, tau.u=1, for.mle=0, pow_scale=1.5, debug = 0, plot = F) 
{
	require(pscl)	
	scale	<- ifelse(for.mle, n.of.y, n.of.y-1)		
	df 		<- n.of.y-1
	
	if (calibrate.tau.u) 
	{
		#calibrate tau.u constrained on yn, alpha and mx.pw	
		tmp		<- nabc.chisqstretch.tau.lowup(mx.pw, tau.u, df, alpha, for.mle=for.mle )
		tau.l	<- tmp[1]
		tau.u	<- tmp[2]
		pw.cmx	<- tmp[3]
		c.l		<- tmp[5]
		c.u		<- tmp[6]	
		if (abs(pw.cmx - mx.pw) > 0.09) 	stop("tau.up not accurate")			
	}
	else
	{
		tau.l 	<- nabc.chisqstretch.tau.low(tau.u, df, alpha, for.mle= for.mle) 
		rej		<- .Call("abcScaledChiSq",	c(scale,df,tau.l,tau.u,alpha,1e-10,100,0.05))
		c.l 	<- rej[1]
		c.u 	<- rej[2]
	}
	
	#truncate pow and compute pow_norm	
	pow_support <- c(tau.l/pow_scale, tau.u*pow_scale) 
	rho 		<- seq(pow_support[1], pow_support[2], length.out = 1000)
	pow			<- nabc.chisqstretch.pow(rho, scale, df, c.l, c.u)
	pow_norm 	<- sum(pow) * diff(rho)[1]
	pow 		<- pow/pow_norm
	
	#compute the norm of lkl, given its support 
	alpha		<- (n.of.x-2)/2	 
	beta		<- s.of.x^2*(n.of.x-1)/2
	lkl_support	<- pow_support
	lkl_norm	<- diff(pigamma(lkl_support,alpha,beta))
	
	integral_range	<- pow_support			
	lkl_arg			<- list(n.of.x= n.of.x, s.of.x= s.of.x, norm = lkl_norm, support = lkl_support)
	pow_arg			<- list(scale = scale, df = df, c.l=c.l, c.u=c.u, norm=pow_norm, support=pow_support)	
	tmp 			<- integrate(nabc.kl.integrand, lower = integral_range[1], upper = integral_range[2], dP=nabc.chisqstretch.sulkl,dQ=dchisqstretch_pow,P_arg=lkl_arg,Q_arg=pow_arg)
	KL_div			<- tmp$value
	if (tmp$message != "OK") 
	{
		warning(tmp$message)
	}
	if (plot) 
	{
		require(reshape)
		require(ggplot2)
		rho_lkl 			<- seq(lkl_support[1], lkl_support[2], length.out = 1000)
		lkl					<- nabc.chisqstretch.sulkl(rho_lkl, n.of.x, s.of.x, lkl_norm, lkl_support)
		df_lkl 				<- data.frame(x = rho_lkl, no = lkl*lkl_norm ,yes = lkl)
		df_lkl$distribution <- "summary likelihood"
		df_pow 				<- data.frame(x = rho, no = pow*pow_norm ,yes = pow)
		df_pow$distribution <- "ABC power"
		df 					<- rbind(df_pow, df_lkl)
		gdf					<- melt(df,id.vars=c("x","distribution"))
		p 					<- ggplot(data = gdf, aes(x = x, y = value, colour = distribution,linetype=variable))
		p					<- p+geom_vline(xintercept=c(tau.l,tau.u),linetype="dotted")
		p					<- p+geom_hline(yintercept= mx.pw,linetype="dotted")
		p					<- p+ geom_line()
		p					<- p+scale_linetype("truncated?")
		p					<- p+xlab(expression(rho))+ylab("")
		p 					<- p + ggtitle(paste("n.of.y=", n.of.y, "\ntau.l=", tau.l,"\ntau.u=", tau.u,"\nKL=", KL_div))
		print(p)
	}
	pw.cmx 	<- ifelse(calibrate.tau.u, pw.cmx, nabc.chisqstretch.pow(rho=1, scale, df, c.l, c.u))	
	ans 	<- c(KL_div = KL_div, tau.l = tau.l, tau.u = tau.u, c.l = c.l, c.u = c.u, pw.cmx = pw.cmx)
	return(ans)
}
#------------------------------------------------------------------------------------------------------------------------
#' Calibrate the lower tolerance interval of the equivalence region for the test of dispersion equivalence
#' @export
#' @param tau.up	upper tolerance of the equivalence region
#' @param df		degrees of freedom
#' @param alpha		level of the equivalence test
#' @param rho.star	point of reference. Defaults to the point of equality rho.star=1
#' @param tol		this algorithm stops when the actual point of reference is less than 'tol' from 'rho.star'
#' @param max.it	this algorithm stops prematurely when the number of iterations to find the equivalence region exceeds 'max.it'
#' @param for.mle	calibrate so that the mode of the power is at the MLE
#' @return tau.low, lower tolerance of the equivalence region
#' @examples	tau.u<- 2.2
#'  yn<- 60
#'	tau.l<- nabc.chisqstretch.tau.low(tau.u, yn-1, 0.01)
nabc.chisqstretch.tau.low<- function(tau.up, df, alpha, rho.star=1, tol= 1e-5, max.it=100, for.mle=0) 
{
	tau.low.lb	<-	1/tau.up		#1/tau.up gives rho.max<1 so we know we only need to go up
	tau.low.ub	<- 	1				#tau.low must be between [tau.low.lb, tau.low.ub]
	error		<- 1
	scale		<- ifelse(for.mle,df+1,df)
	while(abs(error)>tol && round(tau.low.lb,d=10)!=round(tau.low.ub,d=10) && max.it>0)
	{
		max.it<- max.it-1
		tau.low<- (tau.low.lb + tau.low.ub)/2
		rej<- .Call("abcScaledChiSq",	c(scale,df,tau.low,tau.up,alpha,1e-10,100,0.05)	)
		rho<- seq(tau.low,tau.up,by=0.001)
		pw<- nabc.chisqstretch.pow(rho,scale,df,rej[1],rej[2])
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
#' @export
#' @param mx.pw		maximum power at the point of reference (rho.star).
#' @param tau.up.ub	guess on an upper bound on the upper tolerance of the equivalence region
#' @param df		degrees of freedom
#' @param alpha		level of the equivalence test
#' @param rho.star	point of reference. Defaults to the point of equality rho.star=1.
#' @param tol		this algorithm stops when the actual maximum power is less than 'tol' from 'mx.pw'
#' @param max.it	this algorithm stops prematurely when the number of iterations to find the equivalence region exceeds 'max.it'
#' @param for.mle	calibrate so that the mode of the power is at the MLE
#' @return	vector of length 6
#' 	\item{1}{lower tolerance of the equivalence region}		
#' 	\item{2}{upper tolerance of the equivalence region}
#' 	\item{3}{actual maximum power associated with the equivalence region}
#' 	\item{4}{error ie abs(actual power - mx.pw)}
#' 	\item{5}{lower point of critical region}
#' 	\item{6}{upper point of critical region}
#' @examples yn<- 60
#' 	nabc.chisqstretch.tau.lowup(0.9, 2.5, yn-1, 0.01)
nabc.chisqstretch.tau.lowup<- function(mx.pw, tau.up.ub, df, alpha, rho.star=1, tol= 1e-5, max.it=100, for.mle=0)
{
	curr.mx.pw	<- 0
	tau.up.ub	<- tau.up.ub/2
	tmp			<- max.it
	scale		<- ifelse(for.mle,df+1,df)
	while(curr.mx.pw<mx.pw && tmp>0)
	{
		tmp			<- tmp-1
		tau.up.ub	<- 2*tau.up.ub
		tau.low		<- nabc.chisqstretch.tau.low(tau.up.ub, df, alpha, rho.star=rho.star, tol=tol, max.it=max.it, for.mle=for.mle)
		rej			<- .Call("abcScaledChiSq",	c(scale,df,tau.low,tau.up.ub,alpha,1e-10,100,0.05)	)
		curr.mx.pw	<- nabc.chisqstretch.pow(rho.star, scale, df, rej[1], rej[2])		
	}
	if(tmp==0)	stop("nabc.chisqstretch.tau.lowup: could not find tau.up.ub")	
	#print(tau.up.ub)
	tau.up.lb	<- 1
	error		<- 1	
	while(abs(error)>tol && round(tau.up.lb,d=10)!=round(tau.up.ub,d=10) && max.it>0)
	{
		max.it	<- max.it-1
		tau.up	<- (tau.up.lb + tau.up.ub)/2
		tau.low	<- nabc.chisqstretch.tau.low(tau.up, df, alpha, rho.star=rho.star, tol=tol, max.it=max.it, for.mle=for.mle)
		rej		<- .Call("abcScaledChiSq",	c(scale,df,tau.low,tau.up,alpha,1e-10,100,0.05)	)
		curr.mx.pw<- nabc.chisqstretch.pow(rho.star, scale, df, rej[1], rej[2])
		error	<- curr.mx.pw - mx.pw
#print(c(curr.mx.pw, tau.low, tau.up, tau.up.lb, tau.up.ub, max.it))
		if(error<0)
			tau.up.lb<- tau.up
		else
			tau.up.ub<- tau.up
#print(c(abs(error), round(tau.up.lb,d=10)!=round(tau.up.ub,d=10)) )	
	}
	if(max.it==0)	warning("nabc.chisqstretch.tau.lowup: reached max.it")
	c(tau.low,tau.up,curr.mx.pw,abs(error), rej[1], rej[2])
}
#------------------------------------------------------------------------------------------------------------------------
#' Calibrate the number of simulated summary values and the equivalence region for the test of dispersion equivalence
#' @export
#' @param n.of.x	number of observed summary values
#' @param s.of.Sx	standard deviation in the observed summary likelihood
#' @param mx.pw		maximum power at the point of reference (rho.star).
#' @param alpha		level of the equivalence test
#' @param tau.up.ub	guess on an upper bound on the upper tolerance of the equivalence region
#' @param tol		this algorithm stops when the actual variation in the ABC approximation to the summary likelihood is less than 'tol' from 's.of.Sx*s.of.Sx'
#' @param max.it	this algorithm stops prematurely when the number of iterations to calibrate the number of simulated data points exceeds 'max.it'
#' @param for.mle	calibrate so that the mode of the power is at the MLE
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
#' 	y2		<- nabc.chisqstretch.pow(th, yn-1, yn-1, c.l, c.u)
#' #plot the summary likelihood and the abc approximation
#' plot(th,y/mean(y),ylim=range(c(y/mean(y),y2/mean(y2))),type='l')
#' lines(th,y2/mean(y2),col="blue")
nabc.chisqstretch.n.of.y<- function(n.of.x, s.of.Sx, mx.pw, alpha, tau.u.ub=2, tol= 1e-5, max.it=100, for.mle=0)
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
		tmp		<- nabc.chisqstretch.tau.lowup(mx.pw, 2*tau.u, yn.ub-1, alpha, for.mle=for.mle )
#print(c("OK ",tmp))		
		tau.l	<- tmp[1]
		tau.u	<- tmp[2]
		pw.cmx	<- tmp[3]
		c.l		<- tmp[5]
		c.u		<- tmp[6]
		rho		<- seq(tau.l/2,2*tau.u,length.out=1e3)
		scale	<- ifelse(for.mle,yn.ub,yn.ub-1)		
		pw		<- nabc.chisqstretch.pow(rho, scale, yn.ub-1, c.l, c.u)
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
		tmp		<- nabc.chisqstretch.tau.lowup(mx.pw, 2*tau.u, yn-1, alpha, for.mle=for.mle )
		tau.l	<- tmp[1]
		tau.u	<- tmp[2]
		pw.cmx	<- tmp[3]
		c.l		<- tmp[5]
		c.u		<- tmp[6]		
		rho		<- seq(tau.l/2,2*tau.u,length.out=1e3)
		scale	<- ifelse(for.mle,yn,yn-1)
		pw		<- nabc.chisqstretch.pow(rho, scale, yn-1, c.l, c.u)
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
#' @export
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
#' @param for.mle		calibrate so that the mode of the power is at the MLE
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
nabc.chisqstretch<- function(sim, obs.mc, args=NA, verbose= FALSE, tau.l=1, tau.u=1, guess.tau.l=0, alpha=0, normal.test= "sf.test", for.mle=0)
{
	#verbose<- 1
	#sim<- rnorm(100, 8.1, 1.1)
	if(any(is.na(sim)))	stop("get.dist.chisqstretch: error at 1a")
	if(length(obs.mc)>1 || is.na(obs.mc))	stop("get.dist.chisqstretch: error at 1b")
	
	df.sim<- length(sim) - 1
	if(!is.na(args))
	{
		args<- strsplit(args,'/')[[1]]
		if(length(args)==4)
		{
			for.mle		<- as.numeric( args[2] )
			guess.tau.l	<- as.numeric( args[3] )
			tau.u		<- tau.l	<- NA
			alpha		<- as.numeric( args[4] )
		}		
		if(length(args)==5)	
		{
			for.mle		<- as.numeric( args[2] )
			guess.tau.l	<- as.numeric( args[3] )
			tau.u		<- as.numeric( args[4] )
			tau.l		<- NA
			alpha		<- as.numeric( args[5] )
		}
		else if(length(args)==6)
		{
			for.mle		<- as.numeric( args[2] )
			guess.tau.l	<- as.numeric( args[3] )
			tau.l		<- as.numeric( args[4] )
			tau.u		<- as.numeric( args[5] )
			alpha		<- as.numeric( args[6] )
		}
		else 
			stop("get.dist.chisqstretch: error at 1A")
		if(is.na(tau.u))
		{
			tmp<- nabc.chisqstretch.tau.lowup(0.9, 2, df.sim, alpha, for.mle=for.mle)
			tau.l<- tmp[1]
			tau.u<- tmp[2]						
		}
		else if(is.na(tau.l) && guess.tau.l && tau.u>200)		#if tau.u too large, then "nabc.chisqstretch.tau.low" may take a while 
			tau.l<- 1/tau.u
		else if(is.na(tau.l))
			tau.l<- nabc.chisqstretch.tau.low(tau.u, df.sim, alpha, for.mle=for.mle)
		#print(c(df.sim,tau.l,tau.u,alpha))			
		args<- args[1]
	}
	if(alpha<0 || alpha>1)		stop("get.dist.chisqstretch: error at 1e")
	if(tau.u<1 )		stop("get.dist.chisqstretch: error at 1f")
	if(tau.l>1 )		stop("get.dist.chisqstretch: error at 1g")
		 
	ans					<- NABC.DEFAULT.ANS	
	ans["pfam.pval"]	<- nabc.get.pfam.pval(sim,normal.test) 
			
	#get confidence intervals by numerical approximation
	scale				<- ifelse(for.mle,df.sim+1,df.sim)
	tmp					<- .Call("abcScaledChiSq",	c(scale,df.sim,tau.l,tau.u,alpha,1e-10,100,0.05)	)
	if(tmp[4]>1e-10)	stop("get.dist.chisqstretch: error at 3a")
	ans[c("cil","cir","mx.pow")]<- tmp[1:3]
	
	
	ans["error"]		<- var(sim) / obs.mc
	ans[c("tl","tr")]	<- c(tau.l, tau.u)
	ans["nsim"]			<- length(sim)
	ans["lkl"]			<- dchisq(ans["error"]*scale,df.sim)	
	ans["pval"]			<- pchisq(ans["error"]*scale,df.sim)
	ans[c("al","ar")]	<- c(0, 1 - diff( pchisq(ans[c("cil","cir")]*scale,df.sim) ) )
	ans["pval"]			<- ( ans["pval"] - ans["ar"]/2 ) / ( 1 - ans["ar"] )
	ans["link.mc.sim"]	<- var(sim)
	ans["link.mc.obs"]	<- obs.mc
	ans["rho.mc"]		<- log(var(sim) / obs.mc )
	if(verbose)	cat(paste(paste("\n{",args,"<-list(sim.var=",var(sim) ," , obs.var=",obs.mc," , alpha=",alpha," , tau.l=",tau.l," , tau.u=",tau.u,", log.ciu=",log(ans["cir"]),", ",sep=''),paste(names(ans), ans, collapse=', ', sep='='),")}",sep=''))
	ans
}	
