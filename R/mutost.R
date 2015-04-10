#' @title \code{mutost} power function 
#' @description Compute the power of the one-sample equivalence test for population means of normal summary values
#' @export
#' @import data.table pscl reshape2 ggplot2 ash nortest
#' @inheritParams mutost.calibrate
#' @param rho 		vector of quantile
#' @param df		Degrees of freedom
#' @param s.of.T	Standard deviation of the test statistic
#' @param norm 		Normalization constant for the truncated power function.
#' @param support 	Support of the truncated power function (vector of dimension 2).
#' @param log 		If \code{TRUE}, the power function is returned on the log scale. 
#' @note To compute the power, either \code{c.u} or \code{tau.u} and \code{alpha} are required. If \code{tau.u} and \code{alpha} 
#' are given, then the value of \code{c.u} is ignored.  
#'  
#' The power function can be truncated to \code{support} and then standardized with \code{norm}.
#' If one of these is set, the other must be provided too. 
#' @example example/ex.mutost.pow.R
#' @references  http://arxiv.org/abs/1305.4283
mutost.pow <- function(rho, df, s.of.T, c.u=NA, tau.u=NA, alpha=NA, norm=1, support=c(-Inf,Inf), log=FALSE)
{	
	stopifnot( (!is.na(tau.u) & !is.na(alpha)) | (!is.na(c.u) & is.na(tau.u) & is.na(alpha)))
	if(!is.na(tau.u) & !is.na(alpha))
		c.u			<- max(0, qt(alpha, df)*s.of.T+tau.u)
	stopifnot(c.u>=0)
	ans				<- rho
	in_support		<- (rho >= support[1] & rho <= support[2])
	ans[!in_support]<- 0
	
	if(any(in_support))
	{
		suppressWarnings({ #suppress numerical inaccuracy warnings
					ans[in_support] <- .Call("abcMuTOST_power", rho[in_support], df, c.u, s.of.T)/norm
				})
	}	
	if(log)
		ans			<- log(ans)	
	ans
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Area under the \code{mutost} power function
#' @export
#' @description This function computes the area under the power function \code{mutost.pow}.
#' @inheritParams mutost.pow
#' @seealso \code{\link{mutost.pow}}
#' @references  http://arxiv.org/abs/1305.4283
mutost.pow.norm<- function(df, s.of.T, c.u, tau.u=NA, alpha=NA, norm=1, support=c(-Inf,Inf), log=FALSE)
{
	stopifnot( (!is.na(tau.u) & !is.na(alpha)) | (is.na(tau.u) & is.na(alpha)) )
	if(!is.na(tau.u) & !is.na(alpha))
		c.u			<- max(0, qt(alpha, df)*s.of.T+tau.u)
	stopifnot(c.u>0)	
	.Call("abc_mutost_integrate_power", support[1], support[2], .Machine$double.eps^0.25, .Machine$double.eps^0.25, as.double(df), s.of.T, c.u, norm, log)
}	
#------------------------------------------------------------------------------------------------------------------------
# Generic two one sided test (TOST).
# 
# Perform a generic two one sided test. This is an internal function.
# @param tost.args vector of arguments for generic TOST
# @param tau.l		lower tolerance of equivalence region
# @param tau.u		upper tolerance of equivalence region
# @param alpha		level of equivalence test
# @param tost.distr	name of distribution of tost
# @return vector of length 7
generic.tost <- function(tost.args, tau.l, tau.u, alpha, tost.distr="t")
{
	ans<- numeric(7)
	names(ans)<- c("error","p.error","lkl","cl","cu","ass.alpha","ass.pval")
	if(!tost.distr%in%c("t"))	
		stop("unexpected tost.distr")
	
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
# Density of the summary likelihood for \code{mutost}
# 
# Compute the density of the (possibly truncated) summary likelihood for population means of normal summary values
# @inheritParams mutost.calibrate
# @inheritParams mutost.pow
# @param norm 		scalar, 0<\code{norm}<=1, normalization constant for the truncated summary likelihood.
# @param support 	vector of dimension 2, support of the truncated summary likelihood. 
# @note The summary likelihood can be truncated to \code{support} and then standardized with \code{norm}.
# For computational efficiency, both \code{norm} and \code{support} must be provided although each one can be derived from the other.
# \code{support=s.of.x/sqrt(n.of.x)*qt(c(1-norm,1+norm)/2,n.of.x-1)} and \code{norm=diff(pt(support/(s.of.x/sqrt(n.of.x)),n.of.x-1))}.
#
mutost.sulkl <- function(rho, n.of.x, s.of.x, norm = 1, support= c(-Inf,Inf), log=FALSE, debug=0) 
{
	
	stopifnot(n.of.x>0,s.of.x>0,norm<=1,norm>0,support[1]<=support[2])
	
	ans 			<- rho
	in_support 		<- (rho >= support[1] & rho <= support[2])
	ans[!in_support]<- ifelse(log,-Inf,0)
	
	
	if(debug){
		#R code
		if (any(in_support)) 
		{
			if(log)
				ans[in_support] <- dnorm(rho[in_support], sd=s.of.x/sqrt(n.of.x), log=TRUE)-log(norm)
			else
				ans[in_support] <- dnorm(rho[in_support], sd=s.of.x/sqrt(n.of.x))/norm		
		}
		
	}else{
		#C code
		if (any(in_support)){
			ans[in_support]<-.Call("abcMuTOST_sulkl",rho[in_support], n.of.x, s.of.x, norm, log)
		}
		
	}
	return(ans)
}
#------------------------------------------------------------------------------------------------------------------------
# Compute the Kullback-Leibler divergence for the \code{mutost}
#
# Compute Kullback-Leibler divergence between the summary likelihood and the power function of the test of location equivalence 
# @inheritParams mutost.calibrate
# @param calibrate.tau.u if \code{TRUE} the upper tolerance of the equivalence region (\code{tau.u}) is calibrated so that power at the point of reference is equal to \code{mx.pw}
# @param legend.title character; title of the plot (only when \code{plot==TRUE}).
# @examples
# 
# 	mutost.getkl(n.of.x=60,s.of.x=0.1,n.of.y=60,s.of.y=0.3, mx.pw=0.9, alpha=0.01, calibrate.tau.u=TRUE, tau.u=1, plot=TRUE)
#
mutost.getkl <- function(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = F, tau.u = 1, pow_scale = 1.5, debug = 0, plot = F, legend.title='') 
{
	
	
	stopifnot(n.of.x > 1, s.of.x > 0, n.of.y > 1, s.of.y > 0, mx.pw > 0, mx.pw<=1, alpha > 0, alpha<=0.5, tau.u>0, pow_scale > 0)
	
	if (!debug)		#ALL IN C 
	{						
		suppressWarnings({ #suppress numerical inaccuracy warnings
					ans <- .Call("abc_mutost_get_KL", n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u, tau.u, pow_scale)
				})
		
		KL_div <- ans[1]
		tau.u <- ans[2]
		pw.cmx <- ans[3]
		
		if(plot)
		{
			ssn 			<- s.of.x/sqrt(n.of.x)
			lkl_support 	<- pow_support <- c(-tau.u, tau.u) * pow_scale
			lkl_norm 		<- diff(pt(lkl_support/ssn, n.of.x-1))
			suppressWarnings({ #suppress numerical inaccuracy warnings
						pow_norm 	<- .Call("abc_mutost_integrate_pow", pow_support[1], pow_support[2],.Machine$double.eps^0.25,.Machine$double.eps^0.25,as.double(n.of.y-1),s.of.y/sqrt(n.of.y),tau.u,alpha,1,0)
					})			
		}
	} 
	else			#ALL IN R 
	{		
		if (calibrate.tau.u) 
		{
			#calibrate tau.u constrained on yn, alpha and mx.pw	
			tmp 	<- mutost.calibrate.mxpw(mx.pw, n.of.y - 1, s.of.y/sqrt(n.of.y), tau.u, alpha)
			tau.u 	<- tmp[2]
			pw.cmx 	<- tmp[3]
			if (abs(pw.cmx - mx.pw) > 0.09) 	
				stop("tau.up not accurate")
		}
		#truncate pow and compute pow_norm
		pow_support <- c(-tau.u, tau.u) * pow_scale
		#pow_norm <- integrate(mutost.pow, lower = pow_support[1], upper = pow_support[2], df=n.of.y-1, s.of.T=s.of.y/sqrt(n.of.y), tau.u= tau.u, alpha= alpha, norm=1, support= pow_support, log=FALSE)
		suppressWarnings({ #suppress numerical inaccuracy warnings
					pow_norm <- .Call("abc_mutost_integrate_pow", pow_support[1], pow_support[2],.Machine$double.eps^0.25,.Machine$double.eps^0.25,as.double(n.of.y-1),s.of.y/sqrt(n.of.y),tau.u,alpha,1,0)
				})		
		#compute the norm of lkl, given its support 
		ssn 			<- s.of.x/sqrt(n.of.x)
		df 				<- n.of.x - 1
		lkl_support 	<- pow_support
		lkl_norm 		<- diff(pt(lkl_support/ssn, df))
		integral_range	<- pow_support
		lkl_arg 		<- list(n.of.x = n.of.x, s.of.x = s.of.x, norm = lkl_norm, support = lkl_support)
		pow_arg 		<- list(df = n.of.y - 1, s.of.T = s.of.y/sqrt(n.of.y), tau.u = tau.u, alpha = alpha, norm = pow_norm, support = pow_support)
		tmp 			<- integrate(kl.integrand, lower = integral_range[1], upper = integral_range[2], dP = mutost.sulkl, dQ = mutost.pow, P_arg = lkl_arg, Q_arg = pow_arg)
		KL_div 			<- tmp$value
		if (tmp$message != "OK")	
			warning(tmp$message)
		pw.cmx 			<- ifelse(calibrate.tau.u, pw.cmx, mutost.pow(rho = 0, n.of.y - 1, s.of.y/sqrt(n.of.y), tau.u=tau.u, alpha=alpha))
		#print(c(calibrate.tau.u, n.of.y, s.of.y, s.of.y/sqrt(n.of.y), alpha, tau.u, pw.cmx, mutost.pow(rho = 0, n.of.y - 1, s.of.y/sqrt(n.of.y), tau.u=tau.u, alpha=alpha)))		
	}
	
	if (plot) 
	{
		rho 				<- seq(lkl_support[1], lkl_support[2], length.out = 1000)
		lkl 				<- mutost.sulkl(rho, n.of.x, s.of.x, lkl_norm, lkl_support)
		df_lkl 				<- data.frame(x = rho, no = lkl * lkl_norm, yes = lkl)
		df_lkl$distribution <- "summary likelihood"
		pow					<- mutost.pow(rho, df=n.of.y-1, s.of.T=s.of.y/sqrt(n.of.y), tau.u=tau.u, alpha=alpha, norm=pow_norm, support= pow_support, log=FALSE)
		df_pow 				<- data.frame(x = rho, no = pow * pow_norm, yes = pow)
		df_pow$distribution <- "ABC power"
		df 					<- rbind(df_pow, df_lkl)
		gdf 				<- melt(df, id.vars = c("x", "distribution"))
		p 					<- ggplot(data = gdf, aes(x = x, y = value, colour = distribution, linetype = variable))
		p 					<- p + geom_vline(xintercept = c(-tau.u, tau.u), linetype = "dotted")
		p 					<- p + geom_hline(yintercept = mx.pw, linetype = "dotted")
		p 					<- p + geom_line()
		p 					<- p + scale_linetype("truncated and\nstandardized?")
		p 					<- p + xlab(expression(rho)) + ylab("")
		p 					<- p + ggtitle(paste(ifelse(!is.na(legend.title), legend.title, ''),"n.of.y=", n.of.y, "\ntau.u=", tau.u, "\nKL=", KL_div))
		print(p)
	}
	ans <- c(KL_div = KL_div, n.of.y=n.of.y, tau.u = tau.u, pw.cmx = pw.cmx)
	ans
}
#------------------------------------------------------------------------------------------------------------------------
# Plot the calibrated \code{mutost}
mutost.plot<- function(n.of.y, s.of.y, c.u, tau.u, alpha)
{
	tmp			<- data.frame(rho= seq(-1.5*tau.u, 1.5*tau.u, length.out=1024))	
	tmp$power	<- mutost.pow(tmp$rho, n.of.y-1, s.of.y/sqrt(n.of.y), c.u=c.u )
	
	p	<- ggplot(tmp, aes(x=rho, y=power)) + geom_line() + labs(x=expression(rho), y='Power\n(ABC acceptance probability)') +
			scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
			geom_vline(xintercept = c(-tau.u, tau.u), linetype = "dotted") +
			geom_vline(xintercept = c(-c.u, c.u), linetype = "dashed") +
			ggtitle(paste("n.of.y=", n.of.y, "\ntau.u=", tau.u, "\nc.u=", c.u))
	print(p)
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Calibrating \code{mutost} for ABC
#' @description Calibrate the one-sample equivalence test for population means of normal summary values for ABC inference.
#' The one-sample \code{mutost} can be used within ABC to test the null hypothesis that the underlying population mean of the simulated 
#' summary values is not similar to the sample mean of the observed summary values. It is applicable when the simulated and 
#' observed summary values follow a normal distribution, or when normality cannot be rejected. 
#' 
#' Different types of calibrations are available, see Notes for details:
#' \enumerate{ 
#'  \item (\code{what=ALPHA}) compute the ABC false positive rate for given critical region,
#'  \item (\code{what=CR}) calibrate the critical region for given ABC false positive rate,
#'  \item (\code{what=MXPW}) calibrate the critical region and the equivalence region for given ABC false positive rate and maximum power,
#'  \item (\code{what=KL}) calibrate the critical region, the equivalence region and the number of simulated summary values for given ABC false positive rate, maximum power and sample standard deviation of the observed data.
#' }
#' 
#' In the ideal case, the calibration KL is used. However, the KL calibration requires multiple i. i. d. instances of observed summary statistics
#' at each ABC iteration. If this is not available, the MXPW calibration should be used.
#' 
#' Depending on the type of calibration, some of the following inputs must be specified (see Examples).
#' @export 
#' @import data.table pscl reshape2 ggplot2 ash nortest
#' @param n.of.x 	Number of observed summary values 
#' @param s.of.x 	Standard deviation of observed summary values 
#' @param n.of.y 	Number of simulated summary values
#' @param s.of.y 	Standard deviation of simulated summary values
#' @param what		Character string to indicate the type of calibration to be performed 
#' @param c.u		Upper boundary point of the critical region 
#' @param tau.u		Upper boundary point of the equivalence region 
#' @param tau.u.ub	Guess on the upper boundary point of the equivalence region 
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
#' @seealso \code{\link{vartest.calibrate}}, \code{\link{ztest.calibrate}}, \code{\link{ratetest.calibrate}}
#' @note 
#' \enumerate{	
#'  \item (\code{what=ALPHA}) This calibration requires the inputs \code{c.u}, \code{tau.u} with \code{c.u<tau.u} and \code{c.u>0}. 
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
#' \item (\code{what=KL}) This calibration can be used when a set of observed summary values is available. It is desirable because it specifies the number of simulated summary 
#' 				values so that the power is very close to the desired summary likelihood in terms of the KL divergence. 
#' 				The inputs are \code{alpha}, \code{mx.pw}, with default values 0.01 and 0.9 respectively.
#' 				values so that the power is very close to the desired summary likelihood in terms of the KL divergence. Thus, the output  
#' 				consists of the corresponding critical region \code{[c.l, c.u]} (to be used in ABC, see Notes on (2)), the equivalence
#' 				region \code{[tau.l, tau.u]}, and the number of simulated summary values needed (\code{n.of.y}). As a check to the numerical calibrations, 
#' 				the KL divergence is returned (\code{KL}). It is desirable to compare the power to the summary likelihood in terms of the KL divergence, see References.
#' }
#' Note that the underlying test statistic depends on the sample standard deviation of the simulated summary values, \code{s.of.y}.
#' Consequently, the free ABC parameters must be re-calibrated when \code{s.of.y} changes. 
#' 
#' The lower boundary point of the critical region \code{c.l} is always fixed to \code{-c.u}, because this choice implies that
#' the power is maximized at the point of equality \code{rho=0}. This is also commonly used in uncalibrated ABC routines.
#' @example example/ex.mutost.calibrate.R
#' @references  http://arxiv.org/abs/1305.4283
mutost.calibrate<- function(  	n.of.x=NA, s.of.x=NA, n.of.y=n.of.x, s.of.y=NA, what='MXPW',
		c.u=NA, tau.u=NA, tau.u.ub=NA, mx.pw=0.9, alpha=0.01, max.it=100, pow_scale=1.5, tol=1e-5, debug=FALSE, plot=FALSE, plot_debug=FALSE, verbose=FALSE)
{	
	stopifnot(what%in%c('ALPHA','CR','MXPW','KL'))
	if(what=='ALPHA')
	{
		stopifnot(c.u>0, c.u<tau.u, tau.u>0, n.of.y>1, alpha>0, alpha<1, s.of.y>0)
		ans			<- pt( (c.u-tau.u)*sqrt(n.of.y)/s.of.y, n.of.y-1 )
		names(ans)	<- c('alpha')
	}
	if(what=='CR')
	{
		stopifnot(tau.u>0, n.of.y>1, alpha>0, alpha<1, s.of.y>0)
		ans			<- max(0, tau.u + s.of.y/sqrt(n.of.y)*qt(alpha, n.of.y-1))
		ans			<- c(-ans, ans)
		names(ans)	<- c('c.l','c.u')
		if(plot)
			mutost.plot(n.of.y, s.of.y, ans['c.u'], tau.u, alpha)
	}
	if(what=='MXPW')
	{		
		if(is.na(tau.u.ub))
			tau.u.ub<- 3*s.of.y		
		ans			<- mutost.calibrate.mxpw(mx.pw, n.of.y-1, s.of.y/sqrt(n.of.y), tau.u.ub, alpha, rho.star=0, tol=tol, max.it=max.it, debug=debug)
		names(ans)	<- c('tau.l','tau.u','pw.cmx','pw.err')
		tmp			<- max(0, ans['tau.u'] + s.of.y/sqrt(n.of.y)*qt(alpha, n.of.y-1))
		tmp			<- c(-tmp, tmp)
		names(tmp)	<- c('c.l','c.u')
		ans			<- c(tmp, ans)
		if(plot)
			mutost.plot(n.of.y, s.of.y, ans['c.u'], ans['tau.u'], alpha)
	}
	if(what=='KL')
	{
		if(is.na(tau.u.ub))
			tau.u.ub<- 3*s.of.y
		ans			<- mutost.calibrate.kl(  n.of.x, s.of.x, n.of.y, s.of.y, tau.u.ub, 
				mx.pw=mx.pw, alpha=alpha, max.it=max.it, pow_scale=pow_scale, debug=debug, plot=plot, plot_debug=plot_debug, verbose=verbose)
		tmp			<- max(0, ans['tau.u'] + s.of.y/sqrt(ans['n.of.y'])*qt(alpha, ans['n.of.y']-1))
		tmp			<- c(-tmp, tmp)
		names(tmp)	<- c('c.l','c.u')
		ans			<- c(tmp, ans)							 
	}	
	ans
}
#------------------------------------------------------------------------------------------------------------------------
# Calibrate the \code{mutost} with option \code{what=KL}
# @inheritParams mutost.calibrate
mutost.calibrate.kl<- function(  n.of.x, s.of.x, n.of.y, s.of.y, tau.u.ub, 
		mx.pw=0.9, alpha=0.01, max.it=100, pow_scale=1.5, debug=FALSE, plot=FALSE, plot_debug=FALSE, verbose=FALSE)
{	
	KL_args 	<- list(n.of.x=length(obs), s.of.x= sd(obs), n.of.y=length(sim), s.of.y=sd(sim), mx.pw=mx.pw, alpha=alpha, tau.u=tau.u.ub, pow_scale=1.5, debug=debug)
	
	stopifnot(all(c("n.of.x", "s.of.x", "n.of.y", "s.of.y", "mx.pw", "alpha", "tau.u", "pow_scale") %in% names(KL_args)), max.it > 0)
	with(KL_args, stopifnot(n.of.x > 1, s.of.x > 0, n.of.y > 1, s.of.y > 0, mx.pw > 0, mx.pw<=1, alpha > 0, alpha<=0.5, tau.u>0, pow_scale > 0))
	#print(KL_args)
	
	#see if power function "too tight" for yn=xn -> as a marker we use that KL(yn-1) < KL(yn). 
	#In this case, 	the only way to minimize KL is to use yn<xn, which is not what we like.
	#				instead, we give up on mx.pw=0.9	
	KL_divergence			<-  "mutost.getkl"		
	KL_args$calibrate.tau.u <- T
	KL_args$plot 			<- F			
	KL.of.yn 				<- do.call(KL_divergence, KL_args)["KL_div"]
	
	
	n.of.y 					<- KL_args$n.of.y	
	KL_args$n.of.y 			<- n.of.y - 1
	KL.of.yn_m1 			<- do.call(KL_divergence, KL_args)["KL_div"]		
	decrease_n.of.y			<- as.logical(KL.of.yn_m1 < KL.of.yn)
	if(verbose)	
		cat(paste("\ninitial m=",n.of.y, "KL(m)=",KL.of.yn, "KL(m-1)=",KL.of.yn_m1))
	KL_args$n.of.y 			<- n.of.y	
	
	
	if (!debug)		#all in C 
	{		
		if(!decrease_n.of.y)		#case power function not "too tight", adjust yn for fixed mx.pw
		{
			suppressWarnings({ #suppress numerical inaccuracy warnings
						ans 	<- .Call("abc_mutost_calibrate_powerbroader", KL_args, as.integer(max.it))
					})
		}
		else						#case power function "too tight", force yn=xn and give up mx.pw
		{
			KL_args$tau.u 	<- as.double(KL_args$s.of.y / 5)
			stopifnot(KL_args$tau.u>0)
			suppressWarnings({ #suppress numerical inaccuracy warnings
						ans <- .Call("abc_mutost_calibrate_powertighter", KL_args, as.integer(max.it))
					})
		}
	} 
	else			#all in R
	{
		if(!decrease_n.of.y){
			#case power function not "too tight", adjust yn for fixed mx.pw 
			#we have KL(yn-1)>KL(yn), ie KL decreases as yn-1 is incremented. Find upper bound yn.ub such that KL first increases again. 						
			curr.it 		<- max.it
			yn.ub 			<- 2 * (KL_args$n.of.y - 1)
			KL_args$n.of.y 	<- yn.ub
			KL.of.yn_ub 	<- do.call(KL_divergence, KL_args)["KL_div"]
			while (KL.of.yn_ub < KL.of.yn && curr.it > 0) 
			{
				#print(c(yn.ub, KL.of.yn_ub, KL.of.yn, curr.it))
				curr.it 		<- curr.it - 1
				KL.of.yn 		<- KL.of.yn_ub
				yn.ub 			<- 2 * yn.ub
				KL_args$n.of.y 	<- yn.ub
				KL.of.yn_ub 	<- do.call(KL_divergence, KL_args)["KL_div"]				
			}			
			if (curr.it == 0) 	
				stop("could not find upper bound for yn")					
			yn.lb 				<-ifelse(curr.it == max.it,KL_args$n.of.y/2,KL_args$n.of.y/4)		#previous 'yn.ub' could be exactly the minimum, so need to choose 'yn.ub/4' as lower bound 
			#KL is minimized in the open set (yn.lb,yn.ub)						
			if(verbose)	
				cat(paste("\nupper and lower bounds on m:",yn.lb, yn.ub))
			if (plot_debug) 
				pdf("KL_optimization.pdf", onefile = T)		
			KL_args$plot 			<- plot_debug
			KL_args["n.of.y"] 		<- NULL
			tmp 					<- optimize(kl.optimize, interval = c(yn.lb, yn.ub), x_name = "n.of.y", is_integer = T, KL_divergence = KL_divergence, KL_args=KL_args, verbose=verbose, tol = 1)
			if (plot_debug) 
				dev.off()
			KL_args$n.of.y 			<- round(tmp$minimum)
			KL_args$plot 			<- plot
			ans 					<- do.call(KL_divergence, KL_args)
		}
		else					#case power function "too tight", force yn=xn and give up mx.pw
		{
			KL_args$calibrate.tau.u <- F
			KL_args$tau.u 			<- KL_args$s.of.y / 5
			# print(KL_args$tau.u)
			KL_args$plot 			<- plot_debug
			if (plot_debug)
				pdf("KL_initial.pdf", onefile = T)
			KL.of.tau.u 			<- do.call(KL_divergence, KL_args)["KL_div"]
			if (plot_debug)
				dev.off()
			
			#to find tau.u.ub increase tau.u until the KL increases too
			tau.u.ub 				<- 2*KL_args$tau.u
			curr.it 				<- max.it			
			KL_args$tau.u 			<- tau.u.ub
			KL_args$plot 			<- F
			KL.of.tau.u.ub 			<- do.call(KL_divergence, KL_args)["KL_div"]			
			while (KL.of.tau.u.ub < KL.of.tau.u && curr.it > 0) 
			{
				#print(c(tau.u.ub, KL.of.tau.u.ub , KL.of.tau.u, curr.it))
				curr.it 			<- curr.it - 1
				KL.of.tau.u 		<- KL.of.tau.u.ub
				tau.u.ub 			<- 2 * tau.u.ub
				KL_args$tau.u 		<- tau.u.ub
				KL.of.tau.u.ub 		<- do.call(KL_divergence, KL_args)["KL_div"]
			}
			if(curr.it == 0) 	
				stop("could not find upper bound for tau.u")
			#print(c(tau.u.ub, KL.of.tau.u.ub , KL.of.tau.u, curr.it))
			tau.u.lb				<- ifelse(curr.it==max.it, tau.u.ub/2, tau.u.ub/4)
			#minimize KL_tau.u between [tau.u.lb, tau.u.ub]
			if(verbose)	
				cat(paste("\nupper and lower bounds on tau.u:",tau.u.lb, tau.u.ub))
			if(plot_debug)
				pdf("KL_optimization.pdf", onefile = T)
			KL_args$plot 			<- plot_debug
			KL_args["tau.u"]		<- NULL
			tmp 					<- optimize(kl.optimize, interval = c(tau.u.lb,tau.u.ub), x_name="tau.u" ,KL_divergence= KL_divergence, KL_args = KL_args, verbose=verbose)
			if(plot_debug)
				dev.off()
			KL_args$tau.u 	<- tmp$minimum
			KL_args$plot 	<- plot	
			ans 			<- do.call(KL_divergence, KL_args)
		}			
	}
	names(ans)	<- gsub('KL_div','KL',names(ans))
	return(ans)
}
#------------------------------------------------------------------------------------------------------------------------
# Calibrate the \code{mutost} with option \code{what=MXPW}
# @inheritParams mutost.calibrate
# @inheritParams mutost.pow
# @param rho.star	Point of equality
# @examples yn<- 60; ysigma2<- 1
#	mutost.calibrate.mxpw(0.9, yn-1, sqrt(ysigma2/yn), 2, 0.01 )
mutost.calibrate.mxpw<- function(mx.pw, df, s.of.T, tau.up.ub, alpha, rho.star=0, tol= 1e-5, max.it=100, debug=0)
{
	stopifnot(mx.pw<1, mx.pw>0, df>1, s.of.T>0, tau.up.ub>0, alpha>0, alpha<1)
	if(!debug)
	{
		suppressWarnings({	#suppress numerical inaccuracy warnings
					ans<- .Call("abc_mutost_calibrate_tauup_for_mxpw",c(mx.pw, df, s.of.T, tau.up.ub, alpha, rho.star, tol, max.it))
				})
		return(ans)
	}
	#else do R implementation
	curr.mx.pw	<- 0
	tau.up.ub	<- tau.up.ub/2
	tmp			<- max.it
	while(curr.mx.pw<mx.pw && tmp>0)
	{
		tmp			<- tmp-1
		tau.up.ub	<- 2*tau.up.ub
#print( list(rho.star, df, tau.up.ub, s.of.T, alpha) )
		curr.mx.pw	<- .Call("abcMuTOST_pow", as.double(rho.star), df, tau.up.ub, s.of.T, alpha)
		#print( curr.mx.pw )
		#curr.mx.pw	<- mutost.pow(rho.star, df, tau.up.ub, s.of.T, alpha)
		#print( curr.mx.pw )		
	}
	if(tmp==0)	stop("mutost.calibrate.mxpw: could not find tau.up.ub")	
#print(tau.up.ub); stop()
	tau.up.lb	<- 0
	error		<- 1	
	while(abs(error)>tol && round(tau.up.lb,d=10)!=round(tau.up.ub,d=10) && max.it>0)
	{
		max.it		<- max.it-1
		tau.up		<- (tau.up.lb + tau.up.ub)/2
		curr.mx.pw	<- .Call("abcMuTOST_pow", as.double(rho.star), df, tau.up, s.of.T, alpha)
		#print( curr.mx.pw )
		#curr.mx.pw	<- mutost.pow(rho.star, df, tau.up, s.of.T, alpha)
		#print( curr.mx.pw )
		error		<- curr.mx.pw - mx.pw
#print(c(curr.mx.pw, tau.up, tau.up.lb, tau.up.ub, max.it))
		if(error<0)
			tau.up.lb<- tau.up
		else
			tau.up.ub<- tau.up
#print(c(abs(error), round(tau.up.lb,d=10)!=round(tau.up.ub,d=10)) )	
	}
	if(max.it==0)	warning("mutost.calibrate.mxpw: reached max.it")
#	stop("HERE")
	c(-tau.up,tau.up,curr.mx.pw,abs(error))
}
#------------------------------------------------------------------------------------------------------------------------
# Exact TOST for location equivalence.
# 
# Perform the exact TOST for location equivalence when the summary values are normally distributed
# @param sim			simulated summary values
# @param obs			observed summary values
# @param args			argument that contains the equivalence region and the level of the test (see Examples). This is the preferred method for specifying arguments and overwrites the dummy default values
# @param verbose		flag if detailed information on the computations should be printed to standard out
# @param tau.u			upper tolerance of the equivalence region
# @param tau.l			lower tolerance of the equivalence region
# @param alpha			level of the equivalence test
# @param mx.pw			maximum power at the point of equality
# @param annealing		inflation factor of tolerances of the equivalence region
# @param sd.tolerance	tolerance on the difference between the standart deviation of the simulated data and the one used in power calibrations. If this difference if too high, the actual maximum of the power might be considerably different from mx.pw.	
# @param normal.test	name of function with which normality of the summary values is tested
# @param plot			logical; if \code{TRUE} the summary likelihood and calibrated power are ploted.
# @param legend.txt 	character; title of the plot (only when \code{plot==TRUE}).
# @return	vector of length 16 containing
# \item{lkl}{likelihood of T (test statistic for equality) under the TOST}
# \item{error}{test statistic, here p-value of TOST}
# \item{pval}{rescaled p-value of TOST that is expected to follow U(0,1) under the point null hypothesis}
# \item{link.mc.obs}{mean(obs)}
# \item{link.mc.sim}{mean(sim)}
# \item{rho.mc}{mean(sim) - mean(obs)}
# \item{cil}{lower ABC tolerance, here 0}
# \item{cir}{upper ABC tolerance, here alpha}
# \item{tl}{lower tolerance of the equivalence region}
# \item{tr}{upper tolerance of the equivalence region}
# \item{al}{lower ABC tolerance, here 0}
# \item{ar}{upper ABC tolerance, here alpha}
# \item{pfam.pval}{p-value of the normal test for the simulated summary values}
# \item{nsim}{number of simulated summary values}
# \item{mx.pow}{maximum power at the point of equality}
# \item{rho.pow}{power at the point rho.mc}
# @example example/ex.mutost.R
mutost.onesample<- function(sim, obs, args= NA, verbose= FALSE, tau.u= 0, tau.l= -tau.u, alpha= 0, mx.pw=0.9, annealing=1, sd.tolerance=0.05, normal.test= "sf.test", plot=0, legend.txt="")
{
	ans <- c(0, 50, 1, NA, NA, NA, 0, 0, 0, 0, 0, 1, 1, NA, NA, NA)
	names(ans)<- c("lkl", "error", "pval","link.mc.obs","link.mc.sim", "rho.mc", "cil", "cir","tl","tr","al","ar","pfam.pval","nsim","mx.pow","rho.pow")
	
	#compute two sample t-test on either z-scores or untransformed data points
	if(any(is.na(sim)))		
		stop("unexpected NA in sim")
	if(any(is.na(obs)))		
		stop("unexpected NA in obs")
	if(!is.na(args)){			
		# stop("args missing")
		#cat(print.v(sim, prefix="simu", print.char=0))
		#cat(print.v(obs, prefix="obs", print.char=0))
		args<- strsplit(args,'/')[[1]]
		if(length(args)==3)
		{
			annealing	<- as.numeric( args[2] )	#annealing must be at pos 2 otherwise 'ANNEAL.CI.setTau' goes wrong
			obs.n		<- length(obs)
			obs.sd		<- sd(obs)
			tau.u.ub	<- 3*obs.sd
			alpha		<- as.numeric( args[3] )			
		}
		else if(length(args)==4)
		{
			annealing	<- as.numeric( args[2] )
			obs.n		<- length(obs)
			obs.sd		<- sd(obs)
			tau.u.ub	<- as.numeric( args[3] )
			alpha		<- as.numeric( args[4] )			
		}
		else if(length(args)==6)
		{
			annealing	<- as.numeric( args[2] )		
			obs.n		<- as.numeric( args[3] )	#read when obs data just a single point
			obs.sd		<- as.numeric( args[4] )	#read when obs data just a single point
			tau.u.ub	<- as.numeric( args[5] )
			alpha		<- as.numeric( args[6] )			
		}
		else stop("args with unexpected length")
	}else{
		obs.n		<- length(obs)
		obs.sd		<- sd(obs)
		tau.u.ub	<- 3*obs.sd
	}
	
	if(alpha<0 || alpha>1)		
		stop("incorrect alpha")
	if(tau.u.ub<0 )				
		stop("incorrect args for tau.u.ub")
	if(annealing<1)				
		stop("incorrect annealing parameter")
	
	args		<- args[1]
	if(verbose)	
		cat(paste("\ninput call:",args,"annealing=",annealing,"obs.sd=",obs.sd,"tau.u.ub=",tau.u.ub,"alpha=",alpha))
	
	# mx.pw		<- 0.9
	
	ans["pfam.pval"]	<-	abccheck.normal(sim, normal.test)	
	if(!any(diff(sim)>0))	
		return(ans)		
	obs.mean	<- mean(obs)	
	if(verbose)		
		cat(paste("\nn=",obs.n,"obs.mean",obs.mean,"obs.sd=",obs.sd,"sim.sd",sd(sim)))	
	
	if(obs.n<2)	
		stop("length of observed summaries too small, or set 'obs.n' explicitly using args")		
	#trial run on full sd(sim)
	n.sd.refine	<- 8
	tmp			<- sd(sim)
	while(n.sd.refine>0)
	{
		KL_args 		<- list(n.of.x=obs.n, s.of.x= obs.sd, n.of.y=min( obs.n, length(sim) ), s.of.y=tmp, mx.pw=mx.pw, alpha=alpha, tau.u=tau.u.ub, pow_scale=1.5, debug=0)
		abc.param		<- mutost.calibrate( KL_args, 100, debug=0, plot_debug=0, plot=plot, verbose=0)
		if(abc.param["n.of.y"]>length(sim))
		{
			cat(paste("\nnot enough simulated summary values available:",abc.param["n.of.y"], "requested, and used:", length(sim)))
			abc.param	<- mutost.getkl(obs.n, obs.sd, length(sim), sd(sim), mx.pw, alpha, calibrate.tau.u=TRUE, tau.u=3*obs.sd, debug=0, plot=FALSE)
		}		
#print(abc.param); print(sd(sim[seq_len(abc.param["n.of.y"])])); print(tmp); print(sd.tolerance)
		if( abs( log( sd(sim[seq_len(abc.param["n.of.y"])]) / tmp ) ) > sd.tolerance)
		{
			tmp			<- sd(sim[seq_len(abc.param["n.of.y"])])
			n.sd.refine	<- n.sd.refine-1
		}
		else
			n.sd.refine	<-	0		
		#print(n.sd.refine)
	}
	if( abs( log( sd(sim[seq_len(abc.param["n.of.y"])]) / tmp ) ) > sd.tolerance)
		cat(paste("\nWARNING: actual sd.sim differs from the one used in power calibrations so much that pw.max might be considerably different from",mx.pw, tmp, sd(sim[seq_len(abc.param["n.of.y"])]) ))		
	else if(abs(abc.param["pw.cmx"]-mx.pw)>0.09 && abc.param["n.of.y"]>obs.n  && abc.param["n.of.y"]<length(sim))
	{
		print(abc.param)
		stop("unexpected difference in max power when m>n - tau.up not accurate")
	}
	sim.n		<- abc.param["n.of.y"]		
	sim.mean	<- mean(sim[1:sim.n])
	sim.sd		<- sd(sim[1:sim.n])
	tau.u		<- abc.param["tau.u"]*annealing
	tau.l		<- -tau.u		
	if(verbose)		
		cat(paste("\nsim.mean",sim.mean,"sd sim=",sim.sd,"sd sim long=",tmp,"\n Free ABC parameters calibrated to m=",sim.n,"tau.u=",abc.param["tau.u"],"annealed tau.u=",tau.u))	
	if(plot)
		tmp		<- mutost.getkl(obs.n, obs.sd, sim.n, sd(sim), mx.pw, alpha, calibrate.tau.u=FALSE, tau.u=abc.param["tau.u"], debug=0, plot=TRUE, legend.title=legend.txt)
	
	tmp			<- c(	sqrt(sim.n)*(sim.mean-obs.mean-tau.l) / sim.sd,			#[1]	T-	test statistic for -tau (lower test); estimate of the common std dev is simply the std dev in the sample whose sample size is > 1
			sqrt(sim.n)*(sim.mean-obs.mean-tau.u) / sim.sd,			#[2]	T+	test statistic for tau (upper test); estimate of the common std dev is simply the std dev in the sample whose sample size is > 1
			sqrt(sim.n)*(sim.mean-obs.mean) / sim.sd,				#[3]	T	test statistic for equality; estimate of the common std dev is simply the std dev in the sample whose sample size is > 1
			sim.n-1,												#[4] 	degrees of freedom
			sim.sd/sqrt(sim.n),										#[5]	estimate of the std dev of the test statistic is simply the std dev in the sample whose sample size is > 1 divided by that sample size
			sim.sd )												#[6]  	standard deviation of the sample
	tost.ans	<-	generic.tost(tmp, tau.l, tau.u, alpha, tost.distr="t")
	#print("")
	#print(tost.ans)
	
	ans[c("error","cil","cir")]	<- c(tost.ans["p.error"], 0, alpha)
	ans[c("tl","tr","nsim")]	<- c(tau.l,tau.u,sim.n)
	ans[c("lkl","pval")]<-  tost.ans[c("lkl","ass.pval")]
	ans[c("al","ar")]	<- 	c(0,alpha)								
	#ans["mx.pow"]		<-	mutost.pow(0, tmp[4], tmp[5], tau.u,  alpha)
	ans["mx.pow"]		<-	mutost.pow(0, tmp[4], tmp[5], tau.u=abc.param["tau.u"],  alpha=alpha)		#debugging	
	ans["link.mc.sim"]	<- 	sim.mean
	ans["link.mc.obs"]	<- 	obs.mean
	ans["rho.mc"]		<- 	sim.mean - obs.mean
	ans["rho.pow"]		<-	mutost.pow(ans["rho.mc"], tmp[4], tmp[5], tau.u=tau.u, alpha=alpha )
	if(verbose){	cat("\n"); print(ans)	}
	ans
}
