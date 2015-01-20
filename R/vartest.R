#' @title \code{vartest} power function
#' @description Computes the power of the one-sample chi sqare test to
#' test the equivalence of population variances of normal summary values.
#' @param rho 		Vector of error quantiles
#' @param scale		Scaling parameter; typically set to the number of observed summary values
#' @param df		Degrees of freedom; typically set to \code{m-1} where \code{m} is the number of simulated summary values 
#' @param c.l		Lower boundary point of the critical region (equivalent to the lower ABC tolerance \code{epsilon^-})
#' @param c.u 		Upper boundary point of the critical region (equivalent to the upper ABC tolerance \code{epsilon^+})
#' @param trafo		Parameter transformation to translate the power
#' @param norm 		Normalization constant for the truncated power function
#' @param support 	Support of the truncated power function
#' @param log 		If \code{TRUE}, the power function is returned on the log scale. 
#' @note The power function can be truncated to \code{support} and then standardized with \code{norm}.
#' If one of these is set, the other must be provided too.
#' @seealso \code{\link{vartest.pow.norm}}
#' @example example/ex.chisqstretch.pow.R
#' @export
#' @import data.table pscl reshape2 ggplot2 ash nortest
vartest.pow <- function(rho, scale, df, c.l, c.u, norm=1, trafo=1, support=c(0,Inf), log=FALSE)
{	
	ans					<- rep(0,length(rho))
	in_support			<- (rho >= support[1] & rho <= support[2])	
	if(any(in_support))			
		ans[in_support] <- (pchisq( c.u/rho[in_support]*trafo*scale, df ) - pchisq( c.l/rho[in_support]*trafo*scale, df))/norm
	if(log)
		ans<-log(ans)
	ans
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Area under the \code{vartest} power function
#' @export
#' @description This function computes the area under the power function \code{vartest.pow}.
#' @inheritParams vartest.pow
#' @seealso \code{\link{vartest.pow}}
vartest.pow.norm<- function(scale, df, c.l, c.u, trafo= 1, support=c(0,Inf))
{
	ans <- integrate(vartest.pow, lower=support[1], upper=support[2], scale=scale, df=df, c.l=c.l, c.u=c.u, norm=1, trafo=trafo, support=support, log=FALSE)
	ans$value
}	
#------------------------------------------------------------------------------------------------------------------------
# @title Density of the summary likelihood of the variance for normal summary values
# @description		The density of the summary likelihood of the variance for normal summary values is scaled inverse chi square.
# @param rho 		Auxiliary error parameter
# @param n.of.x	Number of observed summary values
# @param s.of.x	Standard deviation of the summary values
# @param trafo		Parameter transformation to compute the summary likelihood on the \code{rho} error space. 		
# @param norm 		scalar, 0<\code{norm}<=1, normalization constant for the truncated summary likelihood.
# @param support 	vector of dimension 2, support of the truncated summary likelihood.
# @param log 		logical; if \code{TRUE}, densities d are given as log(d). 
# @note The summary likelihood can be truncated to \code{support} and then standardized with \code{norm}.
# For computational efficiency, both \code{norm} and \code{support} must be provided although each one can be derived from the other. 
# \code{support=qigamma(c(1-norm,1+norm)/2,(n.of.x-2)/2,s.of.x^2*(n.of.x-1)/2)} and \code{norm=diff(pigamma(support,(n.of.x-2)/2,s.of.x^2*(n.of.x-1)/2)}.
# @seealso \code{\link{vartest.calibrate}}, \code{\link{vartest.getkl}} for an example.
vartest.sulkl<- function(rho, n.of.x, s.of.x, trafo=(n.of.x-1)/n.of.x*s.of.x*s.of.x, norm = 1, support= c(0,Inf), log=FALSE) 
{
	alpha	<- (n.of.x-2)/2	 
	beta	<- s.of.x^2*(n.of.x-1)/2
	ans 				<- rho
	in_support 			<- (rho >= support[1] & rho <= support[2])
	ans[!in_support]	<- 0
	#to match to sulkl, need to apply transformation nu=trafo_fun(rho)= rho*trafo where trafo= mle.nux= (n.of.x-1)/n.of.x * s.of.x
	#dput(alpha); dput(beta); dput(rho[in_support] * trafo)
	if (any(in_support)) 			
		ans[in_support]	<- densigamma(rho[in_support] * trafo, alpha, beta)/norm
	if(log)
		ans				<- log(ans)
	return(ans)
}
#------------------------------------------------------------------------------------------------------------------------
# @title Area under the summary likelihood of the variance for normal summary values
# @description This function computes the area under the summary likelihood \code{vartest.sulkl}.
# @inheritParams vartest.sulkl
# @seealso \code{\link{vartest.calibrate}}
vartest.su.lkl.norm	<- function(n.of.x, s.of.x, trafo=1, support=c(0,Inf))
{
	ans	<- integrate(vartest.sulkl, lower=support[1], upper=support[2],  n.of.x=n.of.x, s.of.x=s.of.x, norm=1, trafo=trafo , support=support, log=FALSE)	
	ans$value
}
#------------------------------------------------------------------------------------------------------------------------
# @title KL divergence between the summary likelihood and the power function of \code{vartest}
# @description Compute the Kullback-Leibler divergence between the summary likelihood and the power function of the \code{vartest} equivalence test.
# The KL divergence is required to calibrate the number of simulated data points of the test.
# @param mx.pw 			Power at the point of reference rho.star=1 (only used when \code{calibrate.tau.u==TRUE}).
# @param alpha 			Level of the equivalence test
# @param calibrate.tau.u	If \code{TRUE} the upper tolerance of the equivalence region (\code{tau.u}) is calibrated so that power at the point of reference is equal to \code{mx.pw}
# @param tau.u				Upper tolerance of the equivalence region. If \code{calibrate.tau.u==TRUE}, \code{tau.u} is just a guess on an upper bound on the upper tolerance of the equivalence region to speed up calibration.
# @param pow_scale 	 	Used to set the support of the pdf associated to the power function. The power is truncated between \code{[tau.l/pow_scale,tau.u*pow_scale]} and then standardized.
# @param plot 				Logical. If \code{plot==TRUE}, the power of the calibrated test is plotted along with the summary likelihood.
# @return	vector of length 6
# 	\item{KL_div}{the Kullback Leibler divergence}	
# 	\item{tau.l}{lower tolerance of the equivalence region}	
# 	\item{tau.u}{upper tolerance of the equivalence region}
# 	\item{c.l}{lower point of the critical region, i.e. lower standard ABC tolerance}	
# 	\item{c.u}{upper point of the critical region, i.e. upper standard ABC tolerance}	
# 	\item{pw.cmx}{actual maximum power at the point of equality}
# @note Whatever the value of \code{calibrate.tau.u}, the lower tolerance of the equivalence region (\code{tau.l}) is always numerically calibrated so that the mode of the power function is at the point of equality rho.star.
vartest.getkl <- function(n.of.x, s.of.x, scale, df, tau.u, mx.pw=0.9, alpha=0.01, pow_scale=1.5, calibrate.tau.u=T, plot = F) 
{
	tau.l<- pw.cmx<- error<- c.l<- c.u<- NA	
	if(calibrate.tau.u)	#calibrate tau.u constrained on yn, alpha and mx.pw 
	{			
		g(tau.l, tau.u, pw.cmx,	error, c.l, c.u)	%<-%	vartest.calibrate.tauup( mx.pw, tau.u, scale, df, alpha )						#tau.u is taken as upper bound on calibrated tau.u
		if (abs(pw.cmx - mx.pw) > 0.09) 	stop("tau.up not accurate")			
	}
	else
	{
		g(tau.l, c.l, c.u, error)	%<-%	vartest.calibrate.taulow(tau.u, scale, df, alpha )	#tau.u is taken as final tau.u
	}
	
	#truncate pow and compute pow_norm	
	pow_support 	<- c(tau.l/pow_scale, tau.u*pow_scale) 	
	pow_norm 		<- vartest.pow.norm(scale, df, c.l, c.u, trafo=1, support=pow_support)
	#compute the norm of lkl, given its support 
	lkl_support	<- pow_support	
	#print(c(n.of.x, s.of.x, (n.of.x-1)/n.of.x*s.of.x*s.of.x)); print(lkl_support)
	lkl_norm	<- vartest.su.lkl.norm(n.of.x, s.of.x, trafo=(n.of.x-1)/n.of.x*s.of.x*s.of.x, support=lkl_support)
	integral_range	<- pow_support			
	lkl_arg			<- list(n.of.x= n.of.x, s.of.x= s.of.x, trafo= (n.of.x-1)/n.of.x*s.of.x*s.of.x, norm = lkl_norm, support = lkl_support)
	pow_arg			<- list(scale = scale, df = df, c.l=c.l, c.u=c.u, norm=pow_norm, support=pow_support, trafo= 1)	
	tmp 			<- integrate(kl.integrand, lower = integral_range[1], upper = integral_range[2], dP=vartest.sulkl, dQ=vartest.pow, P_arg=lkl_arg, Q_arg=pow_arg)
	KL_div			<- tmp$value
	if (tmp$message != "OK") 
	{
		warning(tmp$message)
	}
	if (plot) 
	{
		rho_lkl 			<- seq(lkl_support[1], lkl_support[2], length.out = 1000)
		lkl					<- vartest.sulkl(rho_lkl, n.of.x, s.of.x, trafo= (n.of.x-1)/n.of.x*s.of.x*s.of.x, norm=lkl_norm, support=lkl_support)
		df_lkl 				<- data.frame(x = rho_lkl, yes = lkl, no = lkl*lkl_norm )
		df_lkl$distribution <- "summary likelihood"
		rho_pow	 			<- seq(pow_support[1], pow_support[2], length.out = 1000)
		pow					<- vartest.pow(rho_pow, scale, df, c.l, c.u, trafo= 1, norm=pow_norm)		
		df_pow 				<- data.frame(x = rho_pow, yes = pow, no = pow*pow_norm)
		df_pow$distribution <- "ABC approximation"
		gdf 				<- rbind(df_pow, df_lkl)
		gdf					<- melt(gdf,id.vars=c("x","distribution"))
		p 					<- ggplot(data = gdf, aes(x = x, y = value, colour = distribution,linetype=variable))
		p					<- p+geom_vline(xintercept=c(tau.l,tau.u),linetype="dotted")
		p					<- p+geom_hline(yintercept= mx.pw,linetype="dotted")
		p					<- p+ geom_line()
		p					<- p+scale_linetype("truncated and\nstandardized?")
		p					<- p+xlab(expression(rho))+ylab("")
		p 					<- p + ggtitle(paste("n.of.y=", df+1, "\ntau.l=", tau.l,"\ntau.u=", tau.u,"\nKL=", KL_div))
		print(p)
	}
	pw.cmx 	<- ifelse(calibrate.tau.u, pw.cmx, vartest.pow(rho=1, scale, df, c.l, c.u))	
	c(KL_div = KL_div, tau.l = tau.l, tau.u = tau.u, c.l = c.l, c.u = c.u, pw.cmx = pw.cmx)	
}
#------------------------------------------------------------------------------------------------------------------------
# @title Calibrate the lower tolerance interval of the equivalence region for \code{vartest}
# @description This function calibrates the lower tolerance interval of the equivalence region for the \code{vartest} equivalence test
# so that the mode of the power function is at \code{rho.star=1}.
# @return vector of length 4	
# 	\item{tau.low}{lower tolerance of the equivalence region}	
# 	\item{cl}{lower point of the critical region, i.e. lower standard ABC tolerance}	
# 	\item{cu}{upper point of the critical region, i.e. upper standard ABC tolerance}	
# 	\item{error}{actual error between the mode of the power function and rho.star}
# @seealso \code{\link{vartest.calibrate}}
# @example example/ex.vartest.calibrate.taulow.R 
vartest.calibrate.taulow<- function(tau.up, scale, df, alpha=0.01, rho.star=1, tol= 1e-5, max.it=100, pow_scale=1.5, verbose=0) 
{	
	rho			<- seq(1/(tau.up*pow_scale), tau.up*pow_scale, len=1024)
	tau.low.lb	<- 2/tau.up				
	tmp			<- max.it
	c.rho.max	<- Inf
	while(c.rho.max>rho.star && tmp>0)
	{
		tmp							<- tmp-1
		tau.low.lb					<- tau.low.lb/2
		rej							<- .Call("abcScaledChiSq",	c(scale,df,tau.low.lb,tau.up,alpha,1e-10,max.it,0.05)	)
		if(rej[4]>tol)	stop("compute tau.low.lb: rejection region does not have level alpha within tolerance")
		pw							<- vartest.pow(rho,scale,df,rej[1],rej[2])
		c.rho.max					<- rho[ which.max(pw) ]
		if(verbose)	cat(paste("\ntrial lower bound",tau.low.lb,"with current rho.max",c.rho.max,"critical region",rej[1],rej[2],"error in level is",rej[4]))
	}
	tau.low.ub	<- ifelse(tmp+1<max.it,2*tau.low.lb,tau.up)
	if(verbose)	cat(paste("\nregion for tau.low is",tau.low.lb,tau.low.ub))	
	error		<- 1
	while(abs(error)>tol && round(tau.low.lb,digits=10)!=round(tau.low.ub,digits=10) && max.it>0)
	{
		max.it	<- max.it-1
		tau.low	<- (tau.low.lb + tau.low.ub)/2
#print(c(tau.low, tau.up))
		rej		<- .Call("abcScaledChiSq",	c(scale,df,tau.low,tau.up,alpha,1e-10,max.it,0.05)	)
		if(rej[4]>tol)	stop("compute tau.low: rejection region does not have level alpha within tolerance")
		pw		<- vartest.pow(rho,scale,df,rej[1],rej[2])
#print( c(rho[ which.max(pw) ],pw[ which.max(pw) ], tau.low.lb, tau.low.ub,round(tau.low.lb,digits=10)==round(tau.low.ub,digits=10) ))	
		error	<- rho[ which.max(pw) ] - rho.star
		if(verbose)	cat(paste("\ntrial tau.l=",tau.low,"pw.max is",max(pw),"at",rho[ which.max(pw) ], "it",max.it))
#print( error )			
		if(error<0)
			tau.low.lb<- tau.low
		else
			tau.low.ub<- tau.low			
	}
	if(max.it==0)	warning("vartest.calibrate.taulow: reached max.it")
	c(tau.low=tau.low, cl=rej[1], cu=rej[2], error=error)
}
#------------------------------------------------------------------------------------------------------------------------
# @title Calibrate the upper tolerance interval of the equivalence region for \code{vartest}
# @description This function calibrates the upper tolerance interval of the equivalence region for the \code{vartest} equivalence test
# so that the mode of the power function is at \code{rho.star=1} and so that the power function at \code{rho.star} euqals \code{mx.pw}. 
# This involves recursive recursive calls to re-calibrate the lower tolerance region.
# @return	vector of length 6
# 	\item{tau.low}{lower tolerance of the equivalence region}		
# 	\item{tau.up}{upper tolerance of the equivalence region}
# 	\item{curr.mx.pw}{actual maximum power associated with the equivalence region}
# 	\item{error}{actual error between the power at rho.star and mx.pw}
# 	\item{cl}{lower point of the critical region, i.e. lower standard ABC tolerance}	
# 	\item{cu}{upper point of the critical region, i.e. upper standard ABC tolerance}	
vartest.calibrate.tauup<- function(mx.pw, tau.up.ub, scale, df, alpha=0.01, rho.star=1, tol= 1e-5, max.it=100, pow.scale=1.5, verbose=0)
{
	tau.low		<- cl <- cu	<- NA
	error		<- curr.mx.pw	<- 0
	tau.up.ub	<- tau.up.ub/2	
	tmp			<- max.it		
	while(curr.mx.pw<mx.pw && tmp>0)
	{
		tmp							<- tmp-1
		tau.up.ub					<- 2*tau.up.ub
		g(tau.low, cl, cu, error)	%<-%	vartest.calibrate.taulow(tau.up.ub, scale, df, alpha, rho.star=rho.star, tol=tol, max.it=max.it)
		rho							<- seq(tau.low/pow.scale, tau.up.ub*pow.scale, len=1024)
		pw							<- vartest.pow(rho, scale, df, cl, cu)
		curr.mx.pw					<- max(pw)		
		if(verbose)	cat(paste("\ntrial upper bound",tau.up.ub,"with power",curr.mx.pw,"at rho=",rho[ which.max(pw) ]))
	}
	if(tmp==0)	stop("could not find tau.up.ub")
	if(verbose)	cat(paste("\nFound upper bound",tau.up.ub,"with power",curr.mx.pw,"at rho=",rho[ which.max(pw) ]))
	tau.up.lb	<- 1
	error		<- 1	
	
	while(abs(error)>tol && round(tau.up.lb,digits=10)!=round(tau.up.ub,digits=10) && max.it>0)
	{
		max.it						<- max.it-1
		tau.up						<- (tau.up.lb + tau.up.ub)/2
		g(tau.low, cl, cu, error)	%<-%	vartest.calibrate.taulow(tau.up, scale, df, alpha, rho.star=rho.star, tol=tol, max.it=max.it)
		rho							<- seq(tau.low/pow.scale, tau.up*pow.scale, len=1024)		
		pw							<- vartest.pow(rho, scale, df, cl, cu)
		curr.mx.pw					<- max(pw)		
		error						<- curr.mx.pw - mx.pw
		if(verbose)	cat(paste("\ntrial tau.u",tau.up,"with power",curr.mx.pw,"at rho=",rho[ which.max(pw) ],"search interval",tau.up.lb,tau.up.ub,"error",error))
		if(error<0)
			tau.up.lb<- tau.up
		else
			tau.up.ub<- tau.up
#print(c(abs(error), round(tau.up.lb,digits=10)!=round(tau.up.ub,digits=10)) )	
	}
	if(max.it==0)	warning("vartest.calibrate.tauup: reached max.it")
	c(tau.low=tau.low, tau.up=tau.up, curr.mx.pw=curr.mx.pw,	error=abs(error), cl=cl, cu=cu)
}
#------------------------------------------------------------------------------------------------------------------------
# @title Calibrate the power function of \code{vartest}
# @description This function calibrates the power function of \code{vartest} so that its mode coincides 
# with the mode of the summary likelihood and so that its KL divergence to the summary likelihood is 
# minimized. The function minimizes the KL divergences and includes recursive calls to re-calibrate the
# upper and lower tolerance regions for every new proposed number of simulated summary values.  
# @example example/ex.vartest.calibrate.R
# @example example/ex.vartest.abcreject.R
vartest.calibrate.kl<- function(n.of.x, s.of.x, scale=n.of.x, n.of.y=n.of.x, mx.pw=0.9, alpha=0.01, max.it=100, debug=F, plot=F)
{	
	KL.of.yn_ub<- KL.of.yn<- error <- curr.mx.pw <- tau.low <- cl <- cu	<- NA		
	#KL for initial n.of.y
	KL.of.yn		<- vartest.getkl(n.of.x, s.of.x, scale, n.of.y-1, 3*s.of.x, mx.pw=mx.pw, alpha=alpha, pow_scale=1.5, calibrate.tau.u=T, plot=F)["KL_div"]	
	#KL always decreases from n.of.x. Find upper bound yn.ub such that KL first increases again.	
	curr.it 		<- max.it
	yn.ub 			<- 2 * n.of.y		
	KL.of.yn_ub		<- vartest.getkl(n.of.x, s.of.x, scale, yn.ub-1, 3*s.of.x, mx.pw=mx.pw, alpha=alpha, pow_scale=1.5, calibrate.tau.u=T, plot=F)["KL_div"]		
	while (KL.of.yn_ub < KL.of.yn && curr.it > 0) 
	{
		#print(c(yn.ub, KL.of.yn_ub, KL.of.yn, curr.it))
		curr.it 		<- curr.it - 1
		KL.of.yn 		<- KL.of.yn_ub
		yn.ub 			<- 2 * yn.ub
		KL.of.yn_ub		<- vartest.getkl(n.of.x, s.of.x, scale, yn.ub-1, 3*s.of.x, mx.pw=mx.pw, alpha=alpha, pow_scale=1.5, calibrate.tau.u=T, plot=F)["KL_div"]
		if(debug)	cat(paste("\ntrial upper bound m=",yn.ub,"with KL",KL.of.yn_ub))
	}			
	if (curr.it == 0) 	stop("could not find upper bound for yn")					
	if(debug)	cat(paste("\nFound upper bound m=",yn.ub,"with KL",KL.of.yn_ub))
	yn.lb	<- ifelse(curr.it==max.it, yn.ub/2, yn.ub/4)
	if(debug)	cat(paste("\nupper and lower bounds on m:",yn.lb, yn.ub))
	
	KL_args					<- list(n.of.x=n.of.x, s.of.x=s.of.x, scale=scale, tau.u=3*s.of.x, mx.pw=mx.pw, alpha=alpha, calibrate.tau.u=T, plot=F)	
	tmp 					<- optimize(kl.optimize, interval = c(yn.lb-1, yn.ub-1), x_name = "df", is_integer = T, KL_divergence = "vartest.getkl", KL_args = KL_args, verbose = debug, tol = 1)
	
	n.of.y 										<- round(tmp$minimum)+1
	g(KL_div, tau.l, tau.u, c.l, c.u, pw.cmx)	%<-%	vartest.getkl(n.of.x, s.of.x, scale, n.of.y-1, 3*s.of.x, mx.pw=mx.pw, alpha=alpha, pow_scale=1.5, calibrate.tau.u=T, plot=plot)
	c(n.of.y=n.of.y, tau.l=tau.l, tau.u=tau.u, cl=c.l, cu=c.u, pw.cmx=pw.cmx, KL_div=KL_div)		
}
#------------------------------------------------------------------------------------------------------------------------
vartest.plot<- function(scale, df, c.l, c.u, tau.l, tau.u, pow_scale=1.5)
{
	pow_support <- c(tau.l/pow_scale, tau.u*pow_scale) 	
	pow_norm 	<- vartest.pow.norm(scale, df, c.l, c.u, trafo=1, support=pow_support)	
	tmp			<- data.frame(rho=seq(pow_support[1], pow_support[2], length.out = 1024))	
	tmp$power	<- vartest.pow(tmp$rho, scale, df, c.l, c.u, trafo= 1, norm=pow_norm)*pow_norm	
	
	p	<- ggplot(tmp, aes(x=rho, y=power)) + geom_line() + labs(x=expression(rho), y='Power\n(ABC acceptance probability)') +
			scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
			scale_x_continuous(limits=c(0,pow_support[2])) +
			geom_vline(xintercept = c(tau.l, tau.u), linetype = "dotted") +
			geom_vline(xintercept = c(c.l, c.u), linetype = "dashed") +
			ggtitle(paste("n.of.y=", n.of.y, "\ntau.l=", round(tau.l,d=5), " tau.u=", round(tau.u,d=5), "\nc.l=", round(c.l,d=5), " c.u=", round(c.u,d=5)))
	print(p)
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Calibrating \code{vartest} for ABC
#' @description Calibrate the one-sample equivalence test for population variances of normal summary values for ABC inference.
#' The one-sample \code{vartest} can be used to test the null hypothesis that the underlying population variance of the simulated summary values 
#' is not similar to the variance of the observed summary values. It is applicable when the simulated and observed summary values follow a normal 
#' distribution, or when normality cannot be rejected.
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
#' @param scale 	Scale parameter of the test statistic, usually \code{n.of.x}
#' @param what		Character string to indicate the type of calibration to be performed
#' @param c.l		Lower boundary point of the critical region 
#' @param c.u		Upper boundary point of the critical region
#' @param tau.l		Lower boundary point of the equivalence region 
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
#' @seealso \code{\link{mutost.calibrate}}, \code{\link{ztest.calibrate}}, \code{\link{ratetest.calibrate}}
#' @note 
#' \enumerate{	
#'  \item (\code{what=ALPHA}) This calibration requires the inputs \code{c.l}, \code{c.u}, \code{tau.l}, \code{tau.u} with \code{c.l>tau.l}, \code{c.u<tau.u}, \code{tau.u>1}, \code{tau.l<1}. 
#' 				The output contains the corresponding ABC false positive rate \code{alpha}.
#' 				This option does not specify any of the free ABC parameters, but may be useful to determine the ABC
#' 				false positive rate for uncalibrated ABC routines.
#'  \item (\code{what=CR}) This calibration requires the inputs \code{tau.l}, \code{tau.u}, \code{alpha} with \code{tau.l<1}, \code{tau.u>1} and default \code{alpha=0.01}. 
#' 				The output contains the corresponding critical region \code{[c.l, c.u]}, which corresponds to the ABC tolerance region typically denoted by \code{[-epsilon, epsilon]}. 
#' 				This is an intermediate calibration step and may result in unsuitable power properties (see Examples).				
#'  \item (\code{what=MXPW}) This calibration requires the inputs \code{alpha}, \code{mx.pw}, with default values 0.01 and 0.9 respectively.
#' 				The output contains the corresponding critical region \code{[c.l, c.u]} (to be used in ABC, see Notes on (2)), and 
#' 				the corresponding equivalence region \code{[tau.l, tau.u]} that gives a suitable ABC accept/reject probability if the simulated summary values are close to the observed summary values.
#' 				As a check to the numerical calibrations, the actual power at the point of equality is returned (\code{pw.cmx}). 
#' \item (\code{what=KL}) This calibration can be used when a set of observed summary values is available. It is desirable because it specifies the number of simulated summary 
#' 				values so that the power is very close to the desired summary likelihood in terms of the KL divergence. 
#' 				The inputs are \code{alpha}, \code{mx.pw}, with default values 0.01 and 0.9 respectively.
#' 				The output consists of the corresponding critical region \code{[c.l, c.u]} (to be used in ABC, see Notes on (2)), the equivalence
#' 				region \code{[tau.l, tau.u]}, and the number of simulated summary values needed (\code{n.of.y}). As a check to the numerical calibrations, 
#' 				the KL divergence is returned (\code{KL}). It is desirable to compare the power to the summary likelihood in terms of the KL divergence, see References.
#' }
#' Note that the underlying test statistic only depends on \code{n.of.x}, \code{n.of.y}, \code{s.of.x}, which are all 
#' known before ABC is run. Consequently, the free ABC parameters are calibrated once, before ABC is started. 
#' 
#' The lower boundary point of the critical region \code{c.l} is calibrated numerically, so that
#' the power is maximized at the point of equality \code{rho=1}. The calibrated \code{c.l} does not equal 1/\code{c.u}. 
#' @example example/ex.chisqstretch.calibrate.R
#' @references  http://arxiv.org/abs/1305.4283
vartest.calibrate<- function(n.of.x=NA, s.of.x=NA, n.of.y=NA, what='MXPW', scale=n.of.x, mx.pw=0.9, alpha=0.01, tau.l=NA, tau.u=NA, tau.u.ub=NA, c.l=NA, c.u=NA, max.it=100, tol= 1e-5, pow_scale=1.5, debug=FALSE, plot=FALSE, verbose=FALSE)
{
	stopifnot(what%in%c('ALPHA','CR','MXPW_AT_EQU','MXPW','KL'))
	if(what=='ALPHA')
	{		
		stopifnot(scale>0, c.u<tau.u, tau.u>1, tau.l<1, c.l>tau.l, n.of.y>2, alpha>0, alpha<1)
		ans	<- pchisq(scale*c.u/tau.u, n.of.y-1)-pchisq(scale*c.l/tau.u, n.of.y-1)
		names(ans)	<- 'alpha'
		if(plot)
			vartest.plot(scale, n.of.y-1, c.l, c.u, tau.l, tau.u, pow_scale=pow_scale)
	}
	if(what=='CR')
	{	
		stopifnot(scale>0, tau.u>1, tau.l<1, n.of.y>2, alpha>0, alpha<1, pow_scale>1)
		tmp	<- .Call("abcScaledChiSq",	c(scale, n.of.y-1, tau.l, tau.u, alpha, 1e-10, max.it, 0.05)	)
		ans	<- c(tmp[1], tmp[2], tau.l, tau.u, tmp[3], tmp[4], tmp[5])
		names(ans)<- c('c.l','c.u','tau.l','tau.u','mx.cpw','alpha.err','n.it')
		if(plot)
			vartest.plot(scale, n.of.y-1, ans['c.l'], ans['c.u'], tau.l, tau.u, pow_scale=pow_scale)		
	}
	if(what=='MXPW_AT_EQU')
	{
		stopifnot(scale>0, tau.u>1, n.of.y>2, alpha>0, alpha<1, pow_scale>1, max.it>10, tol<0.2)
		tmp	<- vartest.calibrate.taulow(tau.u, scale, n.of.y-1, alpha=alpha, rho.star=1, tol=tol, max.it=max.it, pow_scale=pow_scale, verbose=verbose)
		ans	<- c(tmp[2],tmp[3],tmp[1],tau.u,tmp[4])
		names(ans)<- c('c.l','c.u','tau.l','tau.u','eq.error')
		if(plot)
			vartest.plot(scale, n.of.y-1, ans['c.l'], ans['c.u'], ans['tau.l'], tau.u, pow_scale=pow_scale)		
	}	
	if(what=='MXPW')
	{
		stopifnot(scale>0, tau.u.ub>1, n.of.y>2, alpha>0, alpha<1, pow_scale>1, max.it>10, tol<0.2)
		tmp <- vartest.calibrate.tauup(mx.pw, tau.u.ub, scale, n.of.y-1, alpha=alpha, rho.star=1, tol=tol, max.it=max.it, pow.scale=pow_scale, verbose=verbose)
		ans	<- c(tmp[5], tmp[6], tmp[1], tmp[2], tmp[3], tmp[4])
		names(ans)<- c('c.l','c.u','tau.l','tau.u','pw.cmx','pw.error')
		if(plot)
			vartest.plot(scale, n.of.y-1, ans['c.l'], ans['c.u'], ans['tau.l'], ans['tau.u'], pow_scale=pow_scale)
	}
	if(what=='KL')
	{
		tmp	<- vartest.calibrate.kl(n.of.x, s.of.x, scale=scale, n.of.y=n.of.x, mx.pw=mx.pw, alpha=alpha, max.it=max.it, debug=debug, plot=plot)
		ans	<- c(tmp[4], tmp[5], tmp[2], tmp[3], tmp[1], tmp[6], tmp[7])
		names(ans)	<- c('c.l','c.u','tau.l','tau.u','n.of.y','pw.cmx','KL.div')
	}
	ans
}