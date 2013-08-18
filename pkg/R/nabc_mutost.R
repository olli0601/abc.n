	
#' Compute the density of the (possible truncated) power of the equivalence test for population means of normal summary values
#' @inheritParams nabc.mutost.kl
#' @inheritParams nabc.mutost.pow
#' @param rho vector of quantile
#' @param norm normalization constant for the truncated power function.
#' @param support vector of dimension 2. Support of the truncated power function.
#' @param log logical; if \code{TRUE}, densities d are given as log(d). 
#' @note The summary likelihood can be truncated to \code{support} and then standardized with \code{norm}.
#' For computational efficiency, both \code{norm} and \code{support} must be provided although each one can be derived (numerically) from the other.
#' @export
#' 	
dMuTOST_pow <- function(rho, df, s.of.T, tau.u, alpha, norm=1, support=c(-Inf,Inf), log=FALSE){
	
	ans<-rho
	in_support<- (rho >= support[1] & rho <= support[2])
	ans[!in_support] <- 0
	
	if(any(in_support)){
		suppressWarnings({ #suppress numerical inaccuracy warnings
					ans[in_support] <- .Call("abcMuTOST_pow", rho[in_support], df, tau.u, s.of.T, alpha)/norm
				})
	}
	
	if(log){
		ans<-log(ans)
	}
	
	return(ans)
}

#------------------------------------------------------------------------------------------------------------------------
#' Perform a generic two one sided test. This is an internal function.
#' @export
#' @param tost.args vector of arguments for generic TOST
#' @param tau.l		lower tolerance of equivalence region
#' @param tau.u		upper tolerance of equivalence region
#' @param alpha		level of equivalence test
#' @param tost.distr	name of distribution of tost
#' @return vector of length 7
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
#' Compute power of the equivalence test for population means of normal summary values
#' @export 
#' @param rho 		true difference in simulated and observed population means
#' @param df		degrees of freedom of the simulated summary values
#' @param tau.u		upper tolerance of the equivalence region
#' @param s.of.T	standard deviation of the test statistic
#' @param alpha		level of the equivalence test
#' @param rtn.fun 	indicator if a function to compute the power should be returned. Defaults to 0.
#' @param force		if TRUE, enforce power computation outside of acceptance region 
#' @return approximate power of the exact test. this is approximate because the standard deviation of the normal model for the simulated summary values is not known.
#' @examples	prior.u<- 5; prior.l<- -prior.u; tau.u	<- 0.75; yn<- 60; ysigma2<- 1; alpha<- 0.01
#' rho	<- seq(prior.l,prior.u,length.out=1e3)
#' nabc.mutost.pow(rho, yn-1, tau.u, sqrt(ysigma2/yn), alpha)
nabc.mutost.pow<- function(rho, df, tau.u, s.of.T, alpha, rtn.fun= FALSE,force= FALSE)
{ 
	x<-	rho
	if(length(x)<10)
		x<- seq(-2*tau.u,2*tau.u,length.out=1e3)
	if(length(rho)<10 && any(rho>2*tau.u) & !force)	stop("unexpected rho")
	
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
#' Compute the density of the (possibly truncated) summary likelihood for population means of normal summary values
#' @inheritParams nabc.mutost.kl
#' @param rho vector of quantile
#' @param norm scalar, 0<\code{norm}<=1, normalization constant for the truncated summary likelihood.
#' @param support vector of dimension 2, support of the truncated summary likelihood.
#' @param log logical; if \code{TRUE}, densities d are given as log(d). 
#' @note The summary likelihood can be truncated to \code{support} and then standardized with \code{norm}.
#' For computational efficiency, both \code{norm} and \code{support} must be provided although each one can be derived from the other.
#' \code{support=s.of.x/sqrt(n.of.x)*qt(c(1-norm,1+norm)/2,n.of.x-1)} and \code{norm=diff(pt(support/(s.of.x/sqrt(n.of.x)),n.of.x-1))}.
#' @export	
#'
nabc.mutost.sulkl <- function(rho, n.of.x, s.of.x, norm = 1, support= c(-Inf,Inf), log=FALSE,debug=0) 
{
	
	stopifnot(n.of.x>0,s.of.x>0,norm<=1,norm>0,support[1]<=support[2])
	
	ans 			<- rho
	in_support 		<- (rho >= support[1] & rho <= support[2])
	ans[!in_support]<- ifelse(log,-Inf,0)


	if(debug){
		#R code
		ssn				<- s.of.x/sqrt(n.of.x)
		df 				<- n.of.x - 1
		if (any(in_support)) 
		{
			if(log)
				ans[in_support] <- dt(rho[in_support]/ssn, df,log=T)-log(ssn*norm)
			else
				ans[in_support] <- dt(rho[in_support]/ssn, df)/ssn/norm		
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
#' Compute Kullback-Leibler divergence between the summary likelihood and the power function of mutost 
#' @param n.of.x number of observed summary values
#' @param s.of.x standard deviation of observed summary values
#' @param n.of.y number of simulated summary values
#' @param s.of.y standard deviation of simulated summary values
#' @param mx.pw maximum power at the point of reference (rho.star=0) (only when \code{calibrate.tau.u==TRUE}).
#' @param alpha level of the equivalence test
#' @param calibrate.tau.u if \code{TRUE} the upper tolerance of the equivalence region (\code{tau.u}) is calibrated so that power at the point of reference is equal to \code{mx.pw}
#' @param tau.u	upper tolerance of the equivalence region. If \code{calibrate.tau.u==TRUE}, it's just a guess on an upper bound on the upper tolerance of the equivalence region.
#' @param pow_scale scale for the support of the standardized power. The power is truncated between \code{pow_scale*[-tau.u,tau.u]} and then standardized.
#' @param debug flag if C implementation is used
#' @param plot whether to plot the two distributions
#' @return	vector of length 3
#' 	\item{KL_div}{the Kullback Leibler divergence}		
#' 	\item{tau.u}{upper tolerance of the equivalence region}
#' 	\item{pw.cmx}{actual maximum power associated with the equivalence region}
#' @export
#' @import ggplot2 reshape2
#' @examples
#' 
#' nabc.mutost.kl(n.of.x=60,s.of.x=0.1,n.of.y=60,s.of.y=0.3, mx.pw=0.9,
#' alpha=0.01, calibrate.tau.u=T, tau.u=1, plot=T)
#'
#------------------------------------------------------------------------------------------------------------------------
nabc.mutost.kl <- function(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = F, tau.u = 1, pow_scale = 1.5, debug = 0, 
	plot = F) {

	#TODO check argument with stopifnot
	
	if (!debug) {
				
		#ALL IN C
		suppressWarnings({ #suppress numerical inaccuracy warnings
			ans <- .Call("abcMuTOST_KL", n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u, tau.u, pow_scale)
		})

		KL_div <- ans[1]
		tau.u <- ans[2]
		pw.cmx <- ans[3]
		
		if(plot){
			lkl_support <- pow_support <- c(-tau.u, tau.u) * pow_scale
			lkl_norm <- diff(pt(lkl_support/ssn, df))
			suppressWarnings({ #suppress numerical inaccuracy warnings
				pow_norm <- .Call("abcMuTOST_pow_integrate_qng", pow_support[1], pow_support[2],.Machine$double.eps^0.25,.Machine$double.eps^0.25,as.double(n.of.y-1),s.of.y/sqrt(n.of.y),tau.u,alpha,1,0)
			})			
		}

	
	} else {
		#ALL IN R
		if (calibrate.tau.u) {
			#calibrate tau.u constrained on yn, alpha and mx.pw	
			tmp <- nabc.mutost.onesample.tau.lowup.pw(mx.pw, n.of.y - 1, s.of.y/sqrt(n.of.y), tau.u, alpha)
			tau.u <- tmp[2]
			pw.cmx <- tmp[3]
			if (abs(pw.cmx - mx.pw) > 0.09) 
				stop("tau.up not accurate")
		}

		#truncate pow and compute pow_norm
		pow_support <- c(-tau.u, tau.u) * pow_scale
		#pow_norm <- integrate(dMuTOST_pow, lower = pow_support[1], upper = pow_support[2], df=n.of.y-1, s.of.T=s.of.y/sqrt(n.of.y), tau.u= tau.u, alpha= alpha, norm=1, support= pow_support, log=FALSE)
		suppressWarnings({ #suppress numerical inaccuracy warnings
				pow_norm <- .Call("abcMuTOST_pow_integrate_qng", pow_support[1], pow_support[2],.Machine$double.eps^0.25,.Machine$double.eps^0.25,as.double(n.of.y-1),s.of.y/sqrt(n.of.y),tau.u,alpha,1,0)
			})		
		#compute the norm of lkl, given its support 
		ssn <- s.of.x/sqrt(n.of.x)
		df <- n.of.x - 1
		lkl_support <- pow_support
		lkl_norm <- diff(pt(lkl_support/ssn, df))

		integral_range <- pow_support

		lkl_arg <- list(n.of.x = n.of.x, s.of.x = s.of.x, norm = lkl_norm, support = lkl_support)
		pow_arg <- list(df = n.of.y - 1, s.of.T = s.of.y/sqrt(n.of.y), tau.u = tau.u, alpha = alpha, norm = pow_norm, support = pow_support)

		tmp <- integrate(nabc.kl.integrand, lower = integral_range[1], upper = integral_range[2], dP = nabc.mutost.sulkl, dQ = dMuTOST_pow, 
			P_arg = lkl_arg, Q_arg = pow_arg)
		KL_div <- tmp$value

		if (tmp$message != "OK") {
			warning(tmp$message)
		}
		
		pw.cmx <- ifelse(calibrate.tau.u, pw.cmx, dMuTOST_pow(rho = 0, n.of.y - 1, s.of.y/sqrt(n.of.y), tau.u, alpha))

	}

	if (plot) {
		require(reshape)
		require(ggplot2)
		rho <- seq(lkl_support[1], lkl_support[2], length.out = 1000)
		lkl <- nabc.mutost.sulkl(rho, n.of.x, s.of.x, lkl_norm, lkl_support)
		df_lkl <- data.frame(x = rho, no = lkl * lkl_norm, yes = lkl)
		df_lkl$distribution <- "summary likelihood"
		pow<-dMuTOST_pow(rho, df=n.of.y-1, s.of.T=s.of.y/sqrt(n.of.y), tau.u, alpha, norm=pow_norm, support= pow_support, log=FALSE)
		df_pow <- data.frame(x = rho, no = pow * pow_norm, yes = pow)
		df_pow$distribution <- "ABC power"
		df <- rbind(df_pow, df_lkl)
		gdf <- melt(df, id.vars = c("x", "distribution"))
		p <- ggplot(data = gdf, aes(x = x, y = value, colour = distribution, linetype = variable))
		p <- p + geom_vline(xintercept = c(-tau.u, tau.u), linetype = "dotted")
		p <- p + geom_hline(yintercept = mx.pw, linetype = "dotted")
		p <- p + geom_line()
		p <- p + scale_linetype("truncated and\nstandardized?")
		p <- p + xlab(expression(rho)) + ylab("")
		p <- p + ggtitle(paste("n.of.y=", n.of.y, "\ntau.u=", tau.u, "\nKL=", KL_div))
		print(p)
	}

	ans <- c(KL_div = KL_div, tau.u = tau.u, pw.cmx = pw.cmx)
	return(ans)
}
#------------------------------------------------------------------------------------------------------------------------
#' Calibrate the equivalence region for the test of location equivalence for given maximum power
#' @export
#' @param mx.pw		maximum power at the point of reference (rho.star).
#' @param df		degrees of freedom
#' @param s.of.T	standard deviation of the test statistic
#' @param tau.up.ub	guess on an upper bound on the upper tolerance of the equivalence region
#' @param alpha		level of the equivalence test
#' @param rho.star	point of reference. Defaults to the point of equality rho.star=0.
#' @param tol		this algorithm stops when the actual maximum power is less than 'tol' from 'mx.pw'
#' @param max.it	this algorithm stops prematurely when the number of iterations to find the equivalence region exceeds 'max.it'
#' @param debug		Flag if C implementation is used.
#' @return	vector of length 4
#' 	\item{1}{lower tolerance of the equivalence region}		
#' 	\item{2}{upper tolerance of the equivalence region}
#' 	\item{3}{actual maximum power associated with the equivalence region}
#' 	\item{4}{error ie abs(actual power - mx.pw)}
#' @examples yn<- 60; ysigma2<- 1; alpha<- 0.01
#'	nabc.mutost.onesample.tau.lowup.pw(0.9, yn-1, sqrt(ysigma2/yn), 2, alpha )
nabc.mutost.onesample.tau.lowup.pw<- function(mx.pw, df, s.of.T, tau.up.ub, alpha, rho.star=0, tol= 1e-5, max.it=100, debug=0)
{
	if(!debug)
	{
		suppressWarnings({	#suppress numerical inaccuracy warnings
					ans<- .Call("abcMuTOST_taulowup_pw",c(mx.pw, df, s.of.T, tau.up.ub, alpha, rho.star, tol, max.it))
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
		#curr.mx.pw	<- nabc.mutost.pow(rho.star, df, tau.up.ub, s.of.T, alpha)
		#print( curr.mx.pw )		
	}
	if(tmp==0)	stop("nabc.mutost.onesample.tau.lowup.pw: could not find tau.up.ub")	
#print(tau.up.ub); stop()
	tau.up.lb	<- 0
	error		<- 1	
	while(abs(error)>tol && round(tau.up.lb,d=10)!=round(tau.up.ub,d=10) && max.it>0)
	{
		max.it		<- max.it-1
		tau.up		<- (tau.up.lb + tau.up.ub)/2
		curr.mx.pw	<- .Call("abcMuTOST_pow", as.double(rho.star), df, tau.up, s.of.T, alpha)
		#print( curr.mx.pw )
		#curr.mx.pw	<- nabc.mutost.pow(rho.star, df, tau.up, s.of.T, alpha)
		#print( curr.mx.pw )
		error		<- curr.mx.pw - mx.pw
#print(c(curr.mx.pw, tau.up, tau.up.lb, tau.up.ub, max.it))
		if(error<0)
			tau.up.lb<- tau.up
		else
			tau.up.ub<- tau.up
#print(c(abs(error), round(tau.up.lb,d=10)!=round(tau.up.ub,d=10)) )	
	}
	if(max.it==0)	warning("nabc.mutost.onesample.tau.lowup.pw: reached max.it")
#	stop("HERE")
	c(-tau.up,tau.up,curr.mx.pw,abs(error))
}
#------------------------------------------------------------------------------------------------------------------------
#' Calibrate the number of simulated summary values and the equivalence region for the test of location equivalence
#' @export
#' @param n.of.x	 number of observed summary values
#' @param s.of.Sx	standard deviation in the observed summary likelihood
#' @param mx.pw		maximum power at the point of reference (rho.star).
#' @param s.of.y 	standard deviation in the simulated summary values
#' @param alpha		level of the equivalence test
#' @param tau.up.ub	guess on an upper bound on the upper tolerance of the equivalence region
#' @param tol		this algorithm stops when the actual variation in the ABC approximation to the summary likelihood is less than 'tol' from 's.of.Sx*s.of.Sx'
#' @param max.it 	this algorithm stops prematurely when the number of iterations to calibrate the number of simulated data points exceeds 'max.it'
#' @param debug		Flag if C implementation is used.
#' @return	vector of length 6
#' 	\item{1}{number of simulated summary values}
#' 	\item{2}{lower tolerance of the equivalence region}		
#' 	\item{3}{upper tolerance of the equivalence region}
#' 	\item{4}{actual variation of the power}
#' 	\item{5}{actual maximum power associated with the equivalence region}
#' 	\item{6}{error ie abs(actual variation - variation in the observed summary likelihood)}
#' @examples prior.u<- 2; prior.l<- -prior.u; tau.u<- 0.75; xn<- yn<- 60; xmu<- 0.5; xsigma2<- ysigma2<- 2; alpha<- 0.01
#'	rho<- seq(prior.l,prior.u,length.out=1e3)
#' 	#summary likelihood 		
#'	y<-	dnorm(rho,0,sqrt(xsigma2/xn))
#'	y<- y / diff(pnorm(c(prior.l,prior.u),0,sqrt(xsigma2/xn)))
#' 	#abc approximation to summary likelihood based on equivalence test 
#'	tmp	<- nabc.mutost.onesample.n.of.y(xn, sqrt(xsigma2/xn), 0.9, sqrt(ysigma2), alpha, tau.u.ub=2*tau.u )
#'	yn	<- tmp[1]
#'	tau.u	<- tmp[3]						
#'	y2<- nabc.mutost.pow(rho, yn-1, tau.u, sqrt(ysigma2/yn), alpha)
#'	rho2<- rho[which(y2!=0)]
#'	y2<- y2[which(y2!=0)]
#'	y2<- y2/sum(diff(rho2)*y2[-1])	
#'	#plot summary likelihood and abc approximation thereof
#'	plot(1,1,type='n',xlim=range(rho),ylim=range(c(y,y2)),xlab=expression(rho))
#'	lines(rho,y,col="red")
#'	lines(rho2,y2,col="blue")			
#'	abline(v=0,col="red")			
nabc.mutost.onesample.n.of.y<- function(n.of.x, s.of.Sx, mx.pw, s.of.y, alpha, tau.u.ub=2, tol= 1e-5, max.it=100, debug=0)
{
	if(!debug)
	{
		suppressWarnings({	#suppress numerical inaccuracy warnings
					ans<- .Call("abcMuTOST_nsim",c(n.of.x, s.of.Sx, mx.pw, s.of.y, tau.u.ub, alpha, 0, tol, max.it))
				})
		return(ans)
	}
	
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
		tmp		<- nabc.mutost.onesample.tau.lowup.pw(mx.pw, yn.ub-1, s.of.y/sqrt(yn.ub), 2*tau.u, alpha, debug=0)
#print(c("OK ",tmp))		
		tau.u	<- tmp[2]
		pw.cmx	<- tmp[3]
		rho		<- seq(-2*tau.u,2*tau.u,length.out=1e3)
		#pw.fun	<- nabc.mutost.pow(rho, yn.ub-1, tau.u, s.of.y/sqrt(yn.ub), alpha, rtn.fun=1)
		suppressWarnings({	#suppress numerical inaccuracy warnings
					pw		<- .Call("abcMuTOST_pow", rho, yn.ub-1, tau.u, s.of.y/sqrt(yn.ub), alpha)
				})	
		#pw.cvar	<- sum(rho*rho*pw.fun(rho)) / sum(pw.fun(rho))			#mean is 0
		pw.cvar	<- sum(rho*rho*pw) / sum(pw)			#mean is 0	
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
		tmp		<- nabc.mutost.onesample.tau.lowup.pw(mx.pw, yn-1, s.of.y/sqrt(yn), 2*tau.u, alpha, debug=0)
		tau.u	<- tmp[2]
		pw.cmx	<- tmp[3]
		rho		<- seq(-2*tau.u,2*tau.u,length.out=1e3)
		suppressWarnings({
					pw		<- .Call("abcMuTOST_pow", rho, yn-1, tau.u, s.of.y/sqrt(yn), alpha)
				})
		#pw.fun	<- nabc.mutost.pow(rho, yn-1, tau.u, s.of.y/sqrt(yn), alpha, rtn.fun=1)
		#pw.cvar	<- sum(rho*rho*pw.fun(rho)) / sum(pw.fun(rho))
		pw.cvar	<- sum(rho*rho*pw) / sum(pw)
		error	<- pw.cvar - s2.of.Sx
#print(c(pw.cvar, yn, yn.lb, yn.ub, max.it))
		if(error<0)
			yn.ub<- yn
		else
			yn.lb<- yn
#print(c(abs(error), (yn.lb+1)!=yn.ub) )	
	}
	c(yn,-tau.u,tau.u,pw.cvar,pw.cmx,abs(error), max.it)
}
#------------------------------------------------------------------------------------------------------------------------
#' Calibrate the equivalence region for the test of location equivalence for given variance of the summary likelihood
#' @export
#' @param s.of.Sx	standard deviation of the summary likelihood
#' @param df		degrees of freedom
#' @param s.of.T	standard deviation of the test statistic
#' @param tau.up.ub	guess on an upper bound on the upper tolerance of the equivalence region
#' @param alpha		level of the equivalence test
#' @param rho.star	point of reference. Defaults to the point of equality rho.star=0.
#' @param tol		this algorithm stops when the actual maximum power is less than 'tol' from 'mx.pw'
#' @param max.it	this algorithm stops prematurely when the number of iterations to find the equivalence region exceeds 'max.it'
#' @param debug		Flag if C implementation is used.
#' @return	vector of length 4
#' 	\item{1}{lower tolerance of the equivalence region}		
#' 	\item{2}{upper tolerance of the equivalence region}
#' 	\item{3}{actual variance associated with the power}
#' 	\item{4}{error ie abs(actual var(power) - var(summary likelihood))}
#' @examples yn<- 60; ysigma2<- 1; alpha<- 0.01
#'	nabc.mutost.onesample.tau.lowup.var(0.002, yn-1, sqrt(ysigma2/yn), 2, alpha )
nabc.mutost.onesample.tau.lowup.var<- function(s.of.Sx, df, s.of.T, tau.up.ub, alpha, rho.star=0, tol= 1e-5, max.it=100, debug=0)
{
	if(!debug)
	{
		suppressWarnings({	#suppress numerical inaccuracy warnings
					ans<- .Call("abcMuTOST_taulowup_var",c(s.of.Sx, df, s.of.T, tau.up.ub, alpha, rho.star, tol, max.it))
				})
		return(ans)
	}
	s2.of.Sx	<- s.of.Sx*s.of.Sx
	pw.cvar		<- 0	
	tmp			<- max.it		
	tau.up.ub	<- tau.up.ub/2
	tmp			<- max.it
	while(pw.cvar<s2.of.Sx && tmp>0)
	{
		tmp			<- tmp-1
		tau.up.ub	<- 2*tau.up.ub		#increase variance until larger than summary likelihood
		rho			<- seq(-2*tau.up.ub,2*tau.up.ub,length.out=1024)
		suppressWarnings({	#suppress numerical inaccuracy warnings
					pw			<- .Call("abcMuTOST_pow", rho, df, tau.up.ub, s.of.T, alpha)
				})
		pw.cvar		<- ifelse(	sum(pw)<EPS,	0, 	sum(rho*rho*pw) / sum(pw)	)	#mean is 0				
	}
	if(tmp==0)	stop("nabc.mutost.onesample.tau.lowup: could not find tau.up.ub")	
	
#print(c(yn.ub, tau.u, pw.cvar, s2.of.Sx, pw.cmx ))	
	tau.up.lb	<- 0
	error		<- 1	
	while(abs(error)>tol && round(tau.up.lb,d=10)!=round(tau.up.ub,d=10) && max.it>0)
	{
		max.it		<- max.it-1
		tau.u		<- (tau.up.lb + tau.up.ub)/2
		rho			<- seq(-2*tau.u,2*tau.u,length.out=1024)
		suppressWarnings({	#suppress numerical inaccuracy warnings
					pw			<- .Call("abcMuTOST_pow", rho, df, tau.u, s.of.T, alpha)
					#print(pw)
				})
		pw.cvar		<- ifelse(	sum(pw)<EPS,	0,	sum(rho*rho*pw) / sum(pw)		)		#mean is 0
		error		<- pw.cvar - s2.of.Sx
#print(c("H2",pw.cvar,s2.of.Sx,tau.u, error, max.it, tau.up.lb, tau.up.ub,round(tau.up.lb,d=10)!=round(tau.up.ub,d=10) ))
		if(error<0)
			tau.up.lb<- tau.u
		else
			tau.up.ub<- tau.u		
#print(c(abs(error), (yn.lb+1)!=yn.ub) )	
	}
	c(-tau.u,tau.u,pw.cvar,abs(error),max.it)
}
#------------------------------------------------------------------------------------------------------------------------
#' Perform the exact TOST for location equivalence when the summary values are normally distributed
#' @export
#' @param sim			simulated summary values
#' @param obs		observed summary values
#' @param args			argument that contains the equivalence region and the level of the test (see Examples). This is the preferred method for specifying arguments and overwrites the dummy default values
#' @param verbose		flag if detailed information on the computations should be printed to standard out
#' @param s.of.x		standard deviation of the observed summary values
#' @param tau.u			upper tolerance of the equivalence region
#' @param tau.l			lower tolerance of the equivalence region
#' @param alpha			level of the equivalence test
#' @param mx.pw			maximum power at the point of equality
#' @param annealing		inflation factor of tolerances of the equivalence region
#' @param normal.test	name of function with which normality of the summary values is tested
#' @return	vector containing
#' \item{error}{test statistic, here p-value of TOST}
#' \item{cil}{lower ABC tolerance, here 0}
#' \item{cir}{upper ABC tolerance, here alpha}
#' \item{mx.pw}{Maximum power at the point of equality}
#' \item{rho.mc}{mean(sim) - obs.mean}
#' @examples tau.u<- 0.5; tau.l<- -tau.u; alpha<- 0.01; xn<- yn<- 60; xmu<- ymu<- 0.5; xsigma2<- ysigma2<- 2
#'	args<- paste("mutost",1,tau.u,alpha,sep='/')
#'	x<- rnorm(xn,xmu,sd=sqrt(xsigma2))
#'	y<- rnorm(yn,ymu,sd=sqrt(ysigma2))
#'	nabc.mutost.onesample(y, x, args= args, verbose= 0)
nabc.mutost.onesample<- function(sim, obs, obs.n=NA, obs.sd=NA, args= NA, verbose= FALSE, tau.u= 0, tau.l= -tau.u, alpha= 0, mx.pw=0.9, annealing=1, normal.test= "sf.test", plot=0, legend.txt="")
{
	verbose<- 1
	ans<- NABC.DEFAULT.ANS
	#compute two sample t-test on either z-scores or untransformed data points
	if(any(is.na(sim)))			stop("nabc.mutost: error at 1a")
	if(any(is.na(obs)))	stop("nabc.mutost: error at 1b")
	if(!is.na(args))
	{
		args<- strsplit(args,'/')[[1]]
		if(length(args)==4)
		{
			standardize	<- as.numeric( args[2] )
			if(standardize!=1)	stop("standardize must be 1")
			tau.u		<- as.numeric( args[3] )
			tau.l		<- -tau.u
			alpha		<- as.numeric( args[4] )
		}
		else if(length(args)==6)
		{
			standardize	<- as.numeric( args[2] )		
			if(!standardize%in%c(2,3,4,5))	stop("standardize must be 2-5 if 6 args specified")			
			annealing	<- as.numeric( args[3] )
			mx.pw		<- ifelse(standardize!=5,	as.numeric( args[4] ), 	0.9)
			obs.sd		<- ifelse(standardize!=5,	obs.sd,					as.numeric( args[4] ))
			tau.u.ub	<- as.numeric( args[5] )
			alpha		<- as.numeric( args[6] )
		}
		else
			stop("nabc.mutost: error at 1c")
		args<- args[1]
	}
	if(!standardize%in%c(0,1,2,3,4,5))	stop("incorrect standardize")		
	if(standardize==1 && is.na(sd(obs)) && is.na(obs.sd))	stop("cannot use standardize==1 when sd(obs) is undefined and obs.sd missing")
	if(standardize==5 && is.na(obs.sd))	stop("cannot use standardize==5 when obs.sd missing")
	if(alpha<0 || alpha>1)		stop("incorrect alpha")
	if(tau.u<0 || tau.l>0)		stop("incorrect tau.u or tau.l")
	if(annealing<1)				stop("incorrect annealing parameter")
#print(standardize)
#standardize<- 0
	ans["pfam.pval"]<-	nabc.get.pfam.pval(sim,normal.test)	
	if(!any(diff(sim)>0))	return(ans)
	obs.n		<- ifelse(!is.na(obs.n),obs.n,length(obs))
	obs.mean	<- mean(obs)
	if(obs.n<2)	stop("length of observed summaries too small, or set 'obs.n' explicitly")		
	sim.n		<- length(sim)	
	
	if(standardize==1)
	{
		if(sim.n>obs.n)
			sim.n	<- obs.n
		obs.sd	<- ifelse(!is.na(sd(obs)),sd(obs),obs.sd)			 	
		sim.mean<- mean(sim[seq.int(1,sim.n)])
		sim.sd	<- sd(sim[seq.int(1,sim.n)])
		sim		<- (sim[seq.int(1,sim.n)]-sim.mean)/sim.sd*obs.sd+sim.mean					
	}
	else if(standardize==2)
	{	
		#print(sim); 
		#print(c(mx.pw,sim.n,sim.sd,alpha))
		if(sim.n>obs.n)
			sim.n	<- obs.n		
		sim.mean<- mean(sim[seq.int(1,sim.n)])
		sim.sd	<- sd(sim[seq.int(1,sim.n)])					
		tmp		<- nabc.mutost.onesample.tau.lowup.pw(mx.pw, sim.n-1, sim.sd/sqrt(sim.n), 2*tau.u.ub, alpha)
		if(tmp[4]>0.09)	stop("tau.up not accurate")		
		tau.l	<- tmp[1]*annealing
		tau.u	<- tmp[2]*annealing	
		if(verbose) 
			cat(paste("\nstd is 2 and sim.n is",sim.n,"annealing is",annealing,mx.pw,tau.l,tau.u))		
		#print(c(annealing,mx.pw,tau.l,tau.u))
		#rho<- seq(tau.l,tau.u,length.out=1e3); y<- nabc.mutost.pow(rho, sim.n-1, tau.u, sim.sd/sqrt(sim.n), alpha); plot(rho,y,type='l')		
	}
	else if(standardize==4)
	{	
		#print(sim); 
		#print(c(mx.pw,sim.n,sim.sd,alpha))
		if(sim.n>obs.n)
			sim.n	<- obs.n		
		sim.mean<- mean(sim[seq.int(1,sim.n)])
		sim.sd	<- sd(sim[seq.int(1,sim.n)])					
		tmp		<- nabc.mutost.onesample.tau.lowup.pw(mx.pw, sim.n-1, sim.sd/sqrt(sim.n), 2*tau.u.ub, alpha)
		if(tmp[4]>0.09)	stop("tau.up not accurate")		
		tau.l	<- max(tmp[1],-tau.u.ub)*annealing
		tau.u	<- min(tmp[2],tau.u.ub)*annealing	
		if(verbose) 
			cat(paste("\nstd is 4 and sim.n is",sim.n,"annealing is",annealing,mx.pw,tau.l,tau.u))		
		#print(c(annealing,mx.pw,tau.l,tau.u))
		#rho<- seq(tau.l,tau.u,length.out=1e3); y<- nabc.mutost.pow(rho, sim.n-1, tau.u, sim.sd/sqrt(sim.n), alpha); plot(rho,y,type='l')		
	}
	else if(standardize==3)
	{
#cat(print.v(sim,print.char=0)); cat(print.v(obs,print.char=0))
		obs.sd		<- ifelse(obs.n>length(obs),sd(sim[1:obs.n]),sd(obs))	
		#instead of the true std dev of the summary likelihood, decided to use the std dev of the summary likelihood truncated to reasonable mass
		tmp			<- seq( -2*obs.sd/sqrt(obs.n),2*obs.sd/sqrt(obs.n), len=1024 )
		su.lkl		<- dt(tmp/obs.sd*sqrt(obs.n), obs.n-1)
		s.of.lkl	<- sqrt( sum( tmp*tmp*su.lkl ) / sum(su.lkl) )			
		#s.of.lkl	<- obs.sd * sqrt( (obs.n-1)/(obs.n-3)/obs.n	)			#assuming empirical Bayes prior on sig2 with df0=n-1, S^2_0=S^2(x) / (n-1)
		sim.sd		<- sd(sim)
		suppressWarnings({
					s.of.pw	<- sqrt( .Call("abcMuTOST_pwvar",c(mx.pw, obs.n, sim.sd/sqrt(obs.n), tau.u.ub, alpha, 0, tol= s.of.lkl*s.of.lkl*1e-5, 100))[1] )	#base case: sim.n=obs.n -- for simplicity use sim.sd
				})
		if(s.of.pw>=s.of.lkl)	#adjust sim.n
		{
#print(c(obs.n,s.of.lkl, mx.pw, sim.sd, alpha, tau.u.ub))
			tmp		<- nabc.mutost.onesample.n.of.y(obs.n, s.of.lkl, mx.pw, sim.sd, alpha, tau.u.ub=2*tau.u.ub, tol= s.of.lkl*s.of.lkl*1e-5)		#for simplicity keep sim.sd fixed even if we use shorter 'sim' overall
			if(abs(tmp[5]-mx.pw)>0.09)	stop("tau.up not accurate")
			sim.n	<- tmp[1]
			options(warn=1)
			if(sim.n>length(sim))
			{
				warning(paste("not enough simulated summary values",sim.n,length(sim)))
				sim.n<- length(sim)
			}
			if(sim.n<obs.n)
			{
				print(c(obs.n,s.of.lkl, mx.pw, sim.sd, alpha, tau.u.ub))
				print.v(sim,print.char=0)
				print.v(obs,print.char=0)
				stop()
			}
			options(warn=2)
			sim.mean<- mean(sim[1:sim.n])
			tau.l	<- tmp[2]*annealing
			tau.u	<- tmp[3]*annealing
			if(verbose)
				cat(paste("\nstd is 3, larger sim.n, and sim.n obs.n is",sim.n,obs.n,"sd sim/obs",sim.sd,obs.sd,"sd pw/lkl",s.of.pw,s.of.lkl))			
		}
		else					#adjust tau.u so that the variance of the summary likelihood is matched even if that means the max pw is > 0.9
		{
			if(sim.n>obs.n)
				sim.n	<- obs.n
			sim.mean<- mean(sim[seq.int(1,sim.n)])
			sim.sd	<- sd(sim[seq.int(1,sim.n)])
			tmp		<- nabc.mutost.onesample.tau.lowup.var(s.of.lkl, sim.n-1, sim.sd/sqrt(sim.n), 2*tau.u.ub, alpha, 0, tol= s.of.lkl*s.of.lkl*1e-8, debug=0)
			print(tmp)
			if(tmp[4]>0.09)	stop("tau.up not accurate")		
			tau.l	<- tmp[1]*annealing
			tau.u	<- tmp[2]*annealing
			if(verbose)
				cat(paste("\nstd is 3, broader mx.pw, and sim.n obs.n is",sim.n,obs.n,"sd sim/obs",sim.sd,obs.sd,"sd pw/pwmatched/lkl",s.of.pw,sqrt(tmp[3]),s.of.lkl))
		}		
	}
	else if(standardize==5)		#obs.sd is specified
	{
#cat(print.v(sim,print.char=0)); cat(print.v(obs,print.char=0))
		tmp			<- seq( -2*obs.sd/sqrt(obs.n),2*obs.sd/sqrt(obs.n), len=1024 )
		su.lkl		<- dt(tmp/obs.sd*sqrt(obs.n), obs.n-1)
		s.of.lkl	<- sqrt( sum( tmp*tmp*su.lkl ) / sum(su.lkl) )					
		#s.of.lkl	<- obs.sd * sqrt( (obs.n-1)/(obs.n-3)/obs.n	)			#assuming empirical Bayes prior on sig2 with df0=n-1, S^2_0=S^2(x) / (n-1)
		sim.sd		<- sd(sim)
		suppressWarnings({
					s.of.pw	<- sqrt( .Call("abcMuTOST_pwvar",c(mx.pw, obs.n, sim.sd/sqrt(obs.n), tau.u.ub, alpha, 0, tol= s.of.lkl*s.of.lkl*1e-5, 100))[1] )	#base case: sim.n=obs.n -- for simplicity use sim.sd
				})
		if(s.of.pw>=s.of.lkl)	#adjust sim.n
		{
#print(c(obs.n,s.of.lkl, mx.pw, sim.sd, alpha, tau.u.ub))
			tmp		<- nabc.mutost.onesample.n.of.y(obs.n, s.of.lkl, mx.pw, sim.sd, alpha, tau.u.ub=2*tau.u.ub, tol= s.of.lkl*s.of.lkl*1e-5)		#for simplicity keep sim.sd fixed even if we use shorter 'sim' overall
			if(abs(tmp[5]-mx.pw)>0.09)	stop("tau.up not accurate")
			sim.n	<- tmp[1]
			options(warn=1)
			if(sim.n>length(sim))
			{
				warning(paste("not enough simulated summary values",sim.n,length(sim)))
				sim.n<- length(sim)
			}
			if(sim.n<obs.n)
			{
				print(c(obs.n,s.of.lkl, mx.pw, sim.sd, alpha, tau.u.ub))
				print.v(sim,print.char=0)
				print.v(obs,print.char=0)
				stop()
			}
			options(warn=2)
			sim.mean<- mean(sim[1:sim.n])
			tau.l	<- tmp[2]*annealing
			tau.u	<- tmp[3]*annealing
			if(verbose)
				cat(paste("\nstd is 5, adj sim.n, and sim.n obs.n is",sim.n,obs.n,"sd sim/obs",sim.sd,obs.sd,"sd pw/lkl",s.of.pw,s.of.lkl))			
		}
		else					#adjust tau.u so that the variance of the summary likelihood is matched even if that means the max pw is > 0.9
		{
			if(sim.n>obs.n)
				sim.n	<- obs.n
			sim.mean<- mean(sim[seq.int(1,sim.n)])
			sim.sd	<- sd(sim[seq.int(1,sim.n)])
			tmp		<- nabc.mutost.onesample.tau.lowup.var(s.of.lkl, sim.n-1, sim.sd/sqrt(sim.n), 2*tau.u.ub, alpha, 0, tol= s.of.lkl*s.of.lkl*1e-5)
			if(tmp[4]>0.09)	stop("tau.up not accurate")		
			tau.l	<- tmp[1]*annealing
			tau.u	<- tmp[2]*annealing
			if(verbose)
				cat(paste("\nstd is 5, broader mx.pw, and sim.n obs.n is",sim.n,obs.n,"sd sim/obs",sim.sd,obs.sd,"sd pw/lkl",sqrt(tmp[3]),s.of.lkl))
		}		
	}
	else if(standardize==6)
	{
		#variables I need:
		n.of.x = length(obs) #QOLLI or obs.n? why obs.n can also be specified in the argument? isn't it always length(obs).
		s.of.x = sd(obs)
		n.of.y = sim.n
		s.of.y = sd(sim)
		
		
		#we start from the case n.of.y=n.of.x (but - for simplicity? - we use sd(sim))	
		#check if KL decreases when the number of simulations increases by one
		#this should be equivalent to check that s.of.pow > s.of.lkl
		KL.sim.n <- nabc.mutost.kl(n.of.x, s.of.x, n.of.x, s.of.y, mx.pw, alpha, calibrate.tau.u = T, tau.u = tau.u.ub)
		KL.sim.n_p1 <- nabc.mutost.kl(n.of.x, s.of.x, n.of.x + 1, s.of.y, mx.pw, alpha, calibrate.tau.u = T, tau.u = tau.u.ub)
		
		increase_n.of.y <- as.logical(KL.sim.n_p1["KL_div"] < KL.sim.n["KL_div"])
		
		if (increase_n.of.y) {
			#increase sim.n so that the KL.div between the summary likelihood and the power is minimised
			tmp <- nabc.calibrate.m.and.tau.yesmxpw.yesKL("nabc.mutost.kl", args = list(n.of.x = n.of.x, s.of.x = s.of.x, n.of.y = n.of.y, s.of.y = s.of.y, 
							mx.pw = mx.pw, alpha = alpha, calibrate.tau.u = T, tau.u = tau.u.ub))
			
			if (abs(tmp["pw.cmx"] - mx.pw) > 0.09) 
				stop("tau.up not accurate")
			sim.n <- tmp["n.of.y"]
			
			options(warn = 1)
			if (sim.n > length(sim)) {
				warning(paste("not enough simulated summary values", sim.n, length(sim)))
				sim.n <- length(sim)
			}
			if (sim.n < obs.n) {
				#QOLLI not sure what this if() does..
				print(c(obs.n, s.of.x, mx.pw, s.of.y, alpha, tau.u.ub))
				print.v(sim, print.char = 0)
				print.v(obs, print.char = 0)
				stop()
			}
			options(warn = 2)
			
			sim.mean <- mean(sim[1:sim.n])
			#QOLLI shall we also recompute sim.sd with sim.n?
			#sim.sd	<- sd(sim[1:sim.n])
			
			tau.u <- tmp["tau.u"] * annealing
			tau.l <- -tau.u
			
		} else {
			#keep n.of.y=n.of.x
			#increase tau.u so that the KL.div between the summary likelihood and the power is minimized even if that means the max pw is > mx.pw
			#use half of tau.u calibrated for n.of.y = n.of.x as lower bound (just to be sure in case sim.sd changes a bit when sim.n <- n.of.y)
			if (sim.n > n.of.x) {
				sim.n <- n.of.x
			}
			tau.u.lb <- KL.sim.n["tau.u"]
			
			sim.mean <- mean(sim[1:sim.n])
			sim.sd <- sd(sim[1:sim.n])
			
			tmp <- nabc.calibrate.tau.nomxpw.yesKL(test_name="mutost", KL_args=list(n.of.x= n.of.x, s.of.x= s.of.x, n.of.y=n.of.y, s.of.y=s.of.y, mx.pw=mx.pw, alpha=alpha), tau.u.lb= tau.u.lb)
			
			#we don't adjust mx.pw
			#if (abs(tmp["pw.cmx"] - mx.pw) > 0.09) 
			#	stop("tau.up not accurate")
			
			tau.u <- tmp["tau.u"] * annealing
			tau.l <- -tau.u
			
		}
	}
	
	tmp			<- c(	sqrt(sim.n)*(sim.mean-obs.mean-tau.l) / sim.sd,			#[1]	T-	test statistic for -tau (lower test); estimate of the common std dev is simply the std dev in the sample whose sample size is > 1
			sqrt(sim.n)*(sim.mean-obs.mean-tau.u) / sim.sd,			#[2]	T+	test statistic for tau (upper test); estimate of the common std dev is simply the std dev in the sample whose sample size is > 1
			sqrt(sim.n)*(sim.mean-obs.mean) / sim.sd,				#[3]	T	test statistic for equality; estimate of the common std dev is simply the std dev in the sample whose sample size is > 1
			sim.n-1,												#[4] 	degrees of freedom
			sim.sd/sqrt(sim.n),										#[5]	estimate of the std dev of the test statistic is simply the std dev in the sample whose sample size is > 1 divided by that sample size
			sim.sd )												#[6]  	standard deviation of the sample
	tost.ans	<-	nabc.generic.tost(tmp, tau.l, tau.u, alpha, tost.distr="t")
	#print("")
	#print(tost.ans)
	ans[c("error","cil","cir")]	<- c(tost.ans["p.error"], 0, alpha)
	ans[c("tl","tr","nsim")]	<- c(tau.l,tau.u,sim.n)
	ans[c("lkl","pval")]<-  tost.ans[c("lkl","ass.pval")]
	ans[c("al","ar")]	<- 	c(0,alpha)								
	ans["mx.pow"]		<-	nabc.mutost.pow(0, tmp[4], tau.u, tmp[5], alpha) 				
	ans["link.mc.sim"]	<- 	sim.mean
	ans["link.mc.obs"]	<- 	obs.mean
	ans["rho.mc"]		<- 	sim.mean - obs.mean
	ans["rho.pow"]		<-	nabc.mutost.pow(ans["rho.mc"], tmp[4], tau.u, tmp[5], alpha,force=T) 				
	
	if(plot)
	{
		rho		<- seq(-4*tau.u,4*tau.u,length.out=1e3)
		pw		<- nabc.mutost.pow(rho, sim.n-1, tau.u, tmp[5], alpha)			
		pw		<- pw/(sum(pw)*diff(rho)[1])
		ylim	<- range(pw)
		su.lkl	<- NA
		if(standardize%in%c(3,5))
		{
			su.lkl		<- dt(rho/obs.sd*sqrt(obs.n), obs.n-1)
			su.lkl		<- su.lkl / (sum(su.lkl)*diff(rho)[1])
			#print( sum(pw)*diff(rho)[1] )
			#print( sqrt( sum( rho^2*pw ) / sum(pw) ) )
			#print( sqrt( sum( rho^2*su.lkl ) / sum(su.lkl) )	) 
			ylim		<- range(c(ylim,su.lkl))
		}	
		plot(1,1,type='n',bty='n',xlim=range(rho),ylim=ylim,ylab="power density",xlab=expression(rho))			
		lines(rho,pw,lty=2)
		if(!any(is.na(su.lkl)))
			lines(rho,su.lkl,lty=1)				
		legend("topright",bty='n',legend=legend.txt)
		legend("topleft",bty='n',legend=c("power","sulkl"),lty=c(2,1))
	}
	#print(ans)
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