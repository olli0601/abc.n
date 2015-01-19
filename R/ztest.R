#' @title \code{ztest} power function 
#' @description Compute the power of the one-sample equivalence test for population means of normal summary values with known population variance.
#' @export
#' @inheritParams ztest.calibrate
#' @param rho 		vector of quantile
#' @param sigma 	Standard deviation of the test statistic.
#' @param norm 		Normalization constant for the truncated power function.
#' @param support 	Support of the truncated power function (vector of dimension 2).
#' @param log 		If \code{TRUE}, the power function is returned on the log scale. 
#' @note The power function can be truncated to \code{support} and then standardized with \code{norm}.
#' If one of these is set, the other must be provided too.
#' @seealso \code{\link{ztest.pow.norm}} 
#' @example example/ex.ztest.pow.R
#' @references  http://arxiv.org/abs/1305.4283
ztest.pow<- function(rho, tau.u, alpha, sigma, norm = 1, support= c(-Inf,Inf), log=FALSE)
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
#' @title Area under the \code{ztest} power function
#' @export
#' @description This function computes the area under the power function \code{ztest.pow}.
#' @inheritParams ztest.pow
#' @seealso \code{\link{ztest.pow}}
ztest.pow.norm<- function(tau.u, sigma, alpha=0.01, support=c(-Inf,Inf))
{
	ans	<- integrate(ztest.pow, lower=support[1], upper=support[2], tau.u=tau.u, sigma=sigma, alpha=alpha, norm=1, support=support, log=FALSE)
	ans$value
}
#------------------------------------------------------------------------------------------------------------------------
# @title Density of the summary likelihood of the autocorrelation error for normal summary values
# @description		The density of the summary likelihood of the autocorrelation error is normal.
# @param rho 		Auxiliary error parameter (vector).
# @param sigma 	Standard deviation of the test statistic.
# @param norm 		scalar, 0<\code{norm}<=1, normalization constant for the truncated summary likelihood.
# @param support 	vector of dimension 2, support of the truncated summary likelihood.
# @param log 		logical; if \code{TRUE}, densities d are given as log(d). 
# @note The summary likelihood can be truncated to \code{support} and then standardized with \code{norm}.
# For computational efficiency, both \code{norm} and \code{support} must be provided although each one can be derived from the other. 
# @seealso \code{\link{chisqstretch.calibrate}}, \code{\link{ztest.sulkl.norm}}
# @return Summary likelihood for the error parameter \code{rho}
# @export	
ztest.sulkl<- function(rho, sigma, norm = 1, support= c(-Inf,Inf), log=FALSE)
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
# @title Area under the density of the summary likelihood of the autocorrelation error for normal summary values
# @inheritParams ztest.sulkl
# @param Tsd	standard deviation of the test statistic
# @return Integral of the summary likelihood on the support
# @seealso \code{\link{ztest.sulkl}}
ztest.sulkl.norm<- function(Tsd, support= c(-Inf,Inf))
{
	diff(pnorm(support,mean=0,sd=Tsd))
}
#------------------------------------------------------------------------------------------------------------------------
ztest.criticalregion<- function(tau.u, sigma, alpha=0.01)
{				 
	c(cl=min(-tau.u/sigma+qnorm(1-alpha),0), cu=max(0,tau.u/sigma+qnorm(alpha)))
}
#------------------------------------------------------------------------------------------------------------------------
# Calibrate the equivalence region of the asymptotic equivalence test for autocorrelations at lag 1 for given maximum power
# @examples tau.u<- 0.09
# 	tau.l<- -tau.u
# 	sim.n<-	5e3
# 	leave.out<- 2
# 	ztest.calibrate.tolerances(0.9, 2, 1/sqrt(floor(sim.n / (1+leave.out))-3), 0.01)
ztest.calibrate.tolerances<- function(mx.pw, tau.up.ub, sigma, alpha=0.01, rho.star=0, tol= 1e-5, max.it=100, pow_scale=2, verbose=0)
{
	stopifnot(mx.pw>0, mx.pw<1, sigma>0, alpha>0, alpha<0.5, tau.up.ub>0, max.it>10, tol<0.01, pow_scale>1)
	curr.mx.pw	<- 0
	tau.up.ub	<- tau.up.ub/2
	tmp			<- max.it
	while(curr.mx.pw<mx.pw && tmp>0)
	{
		tmp			<- tmp-1
		tau.up.ub	<- 2*tau.up.ub
		rho			<- seq(-tau.up.ub*pow_scale, tau.up.ub*pow_scale, len=1024)
		pw			<- ztest.pow(rho, tau.up.ub, alpha, sigma)
		curr.mx.pw	<- max(pw)				 
		if(verbose)	cat(paste("\ntrial upper bound",tau.up.ub,"with power",curr.mx.pw,"at rho=",rho[ which.max(pw) ]))
	}
	if(tmp==0)	stop("could not find tau.up.ub")		
	if(verbose)	cat(paste("\nFound upper bound",tau.up.ub,"with power",curr.mx.pw,"at rho=",rho[ which.max(pw) ]))
	tau.up.lb	<- 0
	error		<- 1	
	while(abs(error)>tol && round(tau.up.lb,digits=10)!=round(tau.up.ub,digits=10) && max.it>0)
	{
		max.it		<- max.it-1
		tau.up		<- (tau.up.lb + tau.up.ub)/2
		rho			<- seq(-tau.up*pow_scale, tau.up*pow_scale, len=1024)
		pw			<- ztest.pow(rho, tau.up, alpha, sigma)
		curr.mx.pw	<- max(pw)				 
		error		<- curr.mx.pw - mx.pw
		if(verbose)	cat(paste("\ntrial tau.u",tau.up,"with power",curr.mx.pw,"at rho=",rho[ which.max(pw) ],"search interval",tau.up.lb,tau.up.ub,"error",error))
		if(error<0)
			tau.up.lb	<- tau.up
		else
			tau.up.ub	<- tau.up
#print(c(abs(error), round(tau.up.lb,digits=10)!=round(tau.up.ub,digits=10)) )	
	}
	if(max.it==0)	warning("reached max.it")
	ans				<- c(-tau.up, tau.up, curr.mx.pw, abs(error))	
	names(ans)		<- c("tau.l", "tau.u", "curr.mx.pw", "error")
	c(ans, ztest.criticalregion(tau.up, sigma, alpha))
}
#------------------------------------------------------------------------------------------------------------------------
#	mx.pw=0.9; alpha=0.01; pow_scale=2; debug = 0; calibrate.tau.u=T; plot = F; verbose=0
# @param s.of.x	standard deviation of the test statistic for the observed summary values; used to construct the summary likelihood
# @param s.of.T	standard deviation of the test statistic for the simulted summary values; used to construct the power function
# @param calibrate.tau.u	If \code{calibrate.tau.u==TRUE} the upper tolerance of the equivalence region (\code{tau.u}) is calibrated so that power at the point of reference is equal to \code{mx.pw}
# @param plot		Logical. If \code{plot==TRUE}, the power of the calibrated test is plotted along with the summary likelihood.
# @return	vector of length 6
# 	\item{KL_div}{KL divergence between the power and the summary likelihood}
# 	\item{tau.l}{lower tolerance of the equivalence region}		
# 	\item{tau.u}{upper tolerance of the equivalence region}
# 	\item{pw.cmx}{actual maximum power associated with the equivalence region}
# 	\item{cl}{lower point of the critical region, i.e. lower standard ABC tolerance}		
# 	\item{cu}{upper point of the critical region, i.e. upper standard ABC tolerance}
# @seealso \code{\link{ztest.calibrate}}
ztest.getkl <- function(s.of.x, s.of.T, tau.u, mx.pw=0.9, alpha=0.01, pow_scale=2, calibrate.tau.u=T, plot = F, verbose=0) 
{
	tau.l<- pw.cmx<- error<- cl<- cu<- NA	
	if(calibrate.tau.u) 
	{			
		g(tau.l, tau.u, pw.cmx,	error, cl, cu)	%<-% ztest.calibrate.tolerances(mx.pw, tau.u, s.of.T, alpha=alpha, pow_scale=pow_scale, verbose=verbose)	#tau.u is taken as upper bound on calibrated tau.u			
		if (abs(pw.cmx - mx.pw) > 0.09) 	stop("tau.up not accurate")			
	}
	#else tau.u is taken as final tau.u
		
	#truncate pow and compute pow_norm	
	rho			<- seq(-tau.u*pow_scale,tau.u*pow_scale, len=1024)	
	support		<- range(rho)
	pw.norm		<- ztest.pow.norm(tau.u, s.of.T, alpha=alpha, support=support)
	pw			<- ztest.pow(rho, tau.u, alpha, s.of.T, norm=pw.norm, support=support, log=FALSE)	
	lkl.norm	<- ztest.sulkl.norm(s.of.x, support=support)	
	lkl			<- ztest.sulkl(rho, s.of.x, norm=lkl.norm, support=support, log=FALSE)
	
	pw.arg		<- list(tau.u=tau.u, sigma=s.of.T, alpha=alpha, norm=pw.norm, support=support)
	lkl.arg		<- list(sigma=s.of.x, norm=lkl.norm, support=support)
	tmp 		<- integrate(kl.integrand, lower=support[1], upper=support[2], dP=ztest.sulkl, dQ=ztest.pow, P_arg=lkl.arg, Q_arg=pw.arg)	
	if (tmp$message != "OK") 	warning(tmp$message)
	KL_div		<- tmp$value
	
	if (plot) 
	{
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
	pw.cmx 	<- ifelse(calibrate.tau.u, pw.cmx, 	ztest.pow(rho=1, tau.u, s.of.T, alpha, norm=1, support=support))
	cu	 	<- ifelse(calibrate.tau.u, cu, 		ztest.criticalregion(tau.u, s.of.T, alpha))
	cl	 	<- ifelse(calibrate.tau.u, cl, 		-cu)
	c(KL_div=KL_div, tau.l= -tau.u, tau.u=tau.u, pw.cmx = pw.cmx, cl=cl, cu=cu)	
}
#------------------------------------------------------------------------------------------------------------------------
ztest.calibrate.kl<- function(n.of.x, n.of.y=n.of.x, n2s= function(n){ 1/sqrt(n-3) }, s2n=function(s){ (1/s)^2+3 }, mx.pw=0.9, alpha=0.01, max.it=100, pow_scale=2, debug=F, plot=F)
{	
	KL.of.yn_ub	<- KL.of.yn		<- error <- curr.mx.pw <- tau.low <- cl <- cu	<- NA
	s.of.x		<- n2s(n.of.x)
	s.of.T		<- n2s(n.of.y)
	#KL for initial n.of.y	
	KL.of.yn		<- ztest.getkl(s.of.x, s.of.T, 3*s.of.T, mx.pw=mx.pw, alpha=alpha, pow_scale=pow_scale)["KL_div"]	
	#KL always decreases from n.of.x. Find upper bound yn.ub such that KL first increases again.	
	curr.it 		<- max.it
	yn.ub 			<- 2 * n.of.y	
	KL.of.yn_ub		<- ztest.getkl(s.of.x, n2s(yn.ub), 3*n2s(yn.ub), mx.pw=mx.pw, alpha=alpha, pow_scale=pow_scale)["KL_div"]
	
	while (KL.of.yn_ub < KL.of.yn && curr.it > 0) 
	{
		#print(c(yn.ub, KL.of.yn_ub, KL.of.yn, curr.it))
		curr.it 		<- curr.it - 1
		KL.of.yn 		<- KL.of.yn_ub
		yn.ub 			<- 2 * yn.ub
		KL.of.yn_ub		<- ztest.getkl(s.of.x, n2s(yn.ub), 3*n2s(yn.ub), mx.pw=mx.pw, alpha=alpha, pow_scale=pow_scale)["KL_div"]				
		if(debug)	cat(paste("\ntrial upper bound m=",yn.ub,"with KL",KL.of.yn_ub))
	}			
	if (curr.it == 0) 	stop("could not find upper bound for yn")					
	if(debug)			cat(paste("\nFound upper bound m=",yn.ub,"with KL",KL.of.yn_ub))
	yn.lb				<- ifelse(curr.it==max.it, yn.ub/2, yn.ub/4)
	if(debug)	cat(paste("\nupper and lower bounds on m:",yn.lb, yn.ub))
	
	
	
	KL_args				<- list(s.of.x=s.of.x, tau.u=3*s.of.x, mx.pw=mx.pw, alpha=alpha, pow_scale=pow_scale, calibrate.tau.u=T, plot=F)	
	tmp 				<- optimize(kl.optimize, interval = n2s(c(yn.ub, yn.lb)), x_name="s.of.T", is_integer=F, KL_divergence = "ztest.getkl", KL_args=KL_args, verbose=debug, tol=1e-5)	
	n.of.y 				<- round(s2n(tmp$minimum))
	
	g(KL_div, tau.l, tau.u, pw.cmx, cl, cu)	%<-%	ztest.getkl(s.of.x, n2s(n.of.y), 3*s.of.x, mx.pw=mx.pw, alpha=alpha, pow_scale=pow_scale,calibrate.tau.u=T, plot=plot)
	
	c(n.of.y=n.of.y, tau.l=tau.l, tau.u=tau.u, pw.cmx=pw.cmx, KL_div=KL_div, cl=cl, cu=cu)		
}
#------------------------------------------------------------------------------------------------------------------------
ztest.plot<- function(n.of.y, sigma, c.u, tau.u, alpha, pow_scale=1.5)
{
	tmp			<- data.frame(rho= seq(-pow_scale*tau.u, pow_scale*tau.u, length.out=1024))	
	tmp$power	<- ztest.pow(tmp$rho, tau.u, alpha, sigma )
	
	p	<- ggplot(tmp, aes(x=rho, y=power)) + geom_line() + labs(x=expression(rho), y='Power\n(ABC acceptance probability)') +
			scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
			geom_vline(xintercept = c(-tau.u, tau.u), linetype = "dotted") +
			#geom_vline(xintercept = c(-c.u, c.u), linetype = "dashed") +
			ggtitle(paste("n.of.y=", n.of.y, "\ntau.u=", tau.u, "\nc.u=", c.u))
	print(p)
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Calibrating \code{ztest} for ABC
#' @description Calibrate the one-sample equivalence test for population means of normal summary values with known population variance for ABC inference.
#' This is an asymptotic test that can be used for several testing problems, for example testing if autocorrelations are similar. 
#' 
#' Different types of calibrations are available, see Notes for details:
#' \enumerate{ 
#'  \item (\code{what=ALPHA}) compute the ABC false positive rate for given critical region,
#'  \item (\code{what=CR}) calibrate the critical region for given ABC false positive rate,
#'  \item (\code{what=MXPW}) calibrate the critical region and the equivalence region for given ABC false positive rate and maximum power,
#'  \item (\code{what=KL}) calibrate the critical region, the equivalence region and the number of simulated summary values for given ABC false positive rate, maximum power and sample standard deviation of the observed data.
#' }
#' The calibration KL should be used. Here, the KL calibration does not require multiple i. i. d. instances of observed summary statistics
#' at each ABC iteration because this is an asymptotic test. 
#' 
#' Depending on the type of calibration, some of the following inputs must be specified (see Examples). 
#' @export 
#' @param n.of.x 	Number of observed summary values  
#' @param n.of.y 	Number of simulated summary values
#' @param n2s 		Function to calculate the standard deviation of the test statistic from n.of.y
#' @param s2n 		Function to calculate n.of.y from the standard deviation of the test statistic
#' @param what		Character string to indicate the type of calibration to be performed 
#' @param c.u		Upper boundary point of the critical region 
#' @param tau.u		Upper boundary point of the equivalence region  
#' @param mx.pw 	Maximum power at the point of equality
#' @param alpha 	Level of the equivalence test
#' @param max.it  	Maximum number of optimization steps at each calibration hierarchy 
#' @param pow_scale Scale for the support of the standardized power. The power is truncated to \code{pow_scale*[-tau.u,tau.u]} and then standardized
#' @param tol		Required error tolerance in calibrating the actual maximum power to the requested maximum power 
#' @param plot  	Flag to plot calibrations
#' @param debug		Flag to switch off C implementation
#' @param plot_debug	Flag to plot at each calibration iteration
#' @param verbose	Flag to run in verbose mode
#' @return	vector
#' @seealso \code{\link{mutost.calibrate}}, \code{\link{vartest.calibrate}}
#' @note 
#' \enumerate{	
#'  \item (\code{what=ALPHA}) This calibration requires the inputs \code{c.u}, \code{tau.u} with \code{0<tau.u} and \code{c.u>0}. 
#' 				The output contains the corresponding ABC false positive rate \code{alpha}.
#' 				This option does not specify any of the free ABC parameters, but may be useful to determine the ABC
#' 				false positive rate for uncalibrated ABC routines.
#'  \item (\code{what=CR}) This calibration requires the input \code{tau.u}, \code{alpha} with \code{tau.u>0} and default \code{alpha=0.01}. 
#' 				The output contains the corresponding critical region \code{[c.l, c.u]}, which corresponds to the ABC tolerance region typically denoted by \code{[-epsilon, epsilon]}. 
#' 				The resulting critical region may be empty under this test, if stochasticity is too large. 
#' 				This is an intermediate calibration step and may result in unsuitable power properties (see Examples).
#'  \item (\code{what=MXPW}) This calibration requires the inputs \code{alpha}, \code{mx.pw}, with default values 0.01 and 0.9 respectively.
#' 				The output contains the corresponding critical region \code{[c.l, c.u]} (to be used in ABC, see Notes on (2)), and 
#' 				the corresponding equivalence region \code{[tau.l, tau.u]} that gives a suitable ABC accept/reject probability if the simulated summary values are close to the observed summary values.
#' 				As a check to the numerical calibrations, the actual power at the point of equality is returned (\code{pw.cmx}).
#' \item (\code{what=KL}) This is the default calibration routine for ABC inference. 
#' 				Inputs are \code{n.of.x} (always known), \code{alpha}, \code{mx.pw}, with default values 0.01 and 0.9 respectively.
#' 				The output consists of the corresponding critical region \code{[c.l, c.u]} (to be used in ABC, see Notes on (2)), the equivalence
#' 				region \code{[tau.l, tau.u]}, and the number of simulated summary values needed (\code{n.of.y}). 
#' 				As a check to the numerical calibrations, the KL divergence is returned (\code{KL}). It is desirable to compare the power to the 
#' 				summary likelihood in terms of the KL divergence, see References.
#' }
#' Note that the underlying test statistic only depends on \code{n.of.y}, which can be calibrated
#' before ABC is run. It is not necessary to re-calibrate this test during ABC runtime. 
#' 
#' The lower boundary point of the critical region \code{c.l} is always fixed to \code{-c.u}, because this choice implies that
#' the power is maximized at the point of equality \code{rho=0}. This is also commonly used in uncalibrated ABC routines.
#' @example example/ex.ztest.calibrate.R
#' @references  http://arxiv.org/abs/1305.4283
ztest.calibrate<- function(	n.of.x=NA, n.of.y=NA, n2s=NULL, s2n=NULL, what='MXPW', 
							mx.pw=0.9, alpha=0.01, c.u=NA, tau.u=NA, tol=1e-5, max.it=100, pow_scale=1.5, debug=FALSE, plot=FALSE, verbose=FALSE)
{
	stopifnot(what%in%c('ALPHA','CR','MXPW','KL'))
	if(what=='ALPHA')
	{
		stopifnot(c.u>0, tau.u>0, n.of.y>1, alpha>0, alpha<1 )
		ans			<- pnorm( c.u-tau.u/n2s(n.of.y) )
		names(ans)	<- c('alpha')
	}
	if(what=='CR')
	{
		stopifnot(tau.u>0, n.of.y>1, alpha>0, alpha<1)
		ans			<- ztest.criticalregion(tau.u, n2s(n.of.y), alpha=alpha)
		names(ans)	<- c('c.l','c.u')
		if(plot)
			ztest.plot(n.of.y, n2s(n.of.y), ans['c.u'], tau.u, alpha)
	}
	if(what=='MXPW')
	{						
		tau.u.ub	<- 3*n2s(2*n.of.y)
		stopifnot(mx.pw>0, mx.pw<1, n.of.y>1, alpha>0, alpha<0.5, tau.up.ub>0, max.it>10, tol<0.01, pow_scale>1)
		tmp			<- ztest.calibrate.tolerances(mx.pw, tau.u.ub, n2s(n.of.y), alpha=alpha, rho.star=0, tol=tol, max.it=max.it, pow_scale=pow_scale, verbose=verbose)
		ans			<- c(tmp[5], tmp[6], tmp[1], tmp[2], tmp[3], tmp[4])
		names(ans)	<- c('c.l','c.u','tau.l','tau.u','pw.cmx','pw.error')
		if(plot)
			ztest.plot(n.of.y, n2s(n.of.y), ans['c.u'], ans['tau.u'], alpha, pow_scale=pow_scale)			
	}
	if(what=='KL')
	{
		tmp	<- ztest.calibrate.kl(n.of.x, n.of.y=n.of.x, n2s=n2s, s2n=s2n, mx.pw=mx.pw, alpha=alpha, max.it=max.it, pow_scale=pow_scale, debug=debug, plot=plot)		
		ans	<- c(tmp[6], tmp[7], tmp[2], tmp[3], tmp[1], tmp[4], tmp[5])
		names(ans)	<- c('c.l','c.u','tau.l','tau.u','n.of.y','pw.cmx','KL.div')								 
	}	
	ans	
}
