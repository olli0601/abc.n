
#' @title Power of the variance test for normal summary values
#' @description This function computes the power function of the two-sided, one-sample scaled chi sqare test for 
#' testing the equivalence of variances for normal summary values.
#' @param rho 		Auxiliary error parameter (vector).
#' @param scale		scaling parameter; typically \code{n}, the number of observed summary values
#' @param df		degrees of freedom; typically \code{m-1} where \code{m} is the number of simulated summary values 
#' @param c.l		lower point of the critical region of the test (equivalent to the lower standard ABC tolerance \code{epsilon^-})
#' @param c.u 		upper point of the critical region of the test (equivalent to the upper standard ABC tolerance \code{epsilon^+})
#' @param norm 		normalization constant for the truncated power function.
#' @param support 	vector of dimension 2. Support of the truncated power function.
#' @param log l		ogical; if \code{TRUE}, densities d are given as log(d). 
#' @note The summary likelihood can be truncated to \code{support} and then standardized with \code{norm}.
#' For computational efficiency, both \code{norm} and \code{support} must be provided although each one can be derived (numerically) from the other.
#' @seealso \code{\link{chisqstretch.pow.norm}}, \code{\link{chisqstretch.calibrate}}
#' @example example/ex.chisqstretch.pow.R
#' @export
#' 	
chisqstretch.pow <- function(rho, scale, df, c.l, c.u, norm=1, trafo=1, support=c(0,Inf), log=FALSE){
	
	ans					<- rep(0,length(rho))
	in_support			<- (rho >= support[1] & rho <= support[2])	
	if(any(in_support))			
		ans[in_support] <- (pchisq( c.u/rho[in_support]*trafo*scale, df ) - pchisq( c.l/rho[in_support]*trafo*scale, df))/norm
	if(log)
		ans<-log(ans)
	ans
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Area under the power function of the variance test for normal summary values
#' @description This function computes the area under the power function \code{chisqstretch.pow.norm}.
#' @inheritParams chisqstretch.pow
#' @seealso \code{\link{chisqstretch.pow}}, \code{\link{chisqstretch.calibrate}}
chisqstretch.pow.norm<- function(scale, df, c.l, c.u, trafo= 1, support=c(0,Inf))
{
	ans <- integrate(chisqstretch.pow, lower=support[1], upper=support[2], scale=scale, df=df, c.l=c.l, c.u=c.u, norm=1, trafo=trafo, support=support, log=FALSE)
	ans$value
}	
#------------------------------------------------------------------------------------------------------------------------
#' @title Density of the summary likelihood of the variance for normal summary values
#' @description		The density of the summary likelihood of the variance for normal summary values is scaled inverse chi square.
#' @param rho 		Auxiliary error parameter
#' @param n.of.x	Number of observed summary values
#' @param s.of.x	Standard deviation of the summary values
#' @param trafo		Parameter transformation to compute the summary likelihood on the \code{rho} error space. 		
#' @param norm 		scalar, 0<\code{norm}<=1, normalization constant for the truncated summary likelihood.
#' @param support 	vector of dimension 2, support of the truncated summary likelihood.
#' @param log 		logical; if \code{TRUE}, densities d are given as log(d). 
#' @note The summary likelihood can be truncated to \code{support} and then standardized with \code{norm}.
#' For computational efficiency, both \code{norm} and \code{support} must be provided although each one can be derived from the other. 
#' \code{support=qigamma(c(1-norm,1+norm)/2,(n.of.x-2)/2,s.of.x^2*(n.of.x-1)/2)} and \code{norm=diff(pigamma(support,(n.of.x-2)/2,s.of.x^2*(n.of.x-1)/2)}.
#' @seealso \code{\link{chisqstretch.calibrate}}, \code{\link{chisqstretch.calibrate.tolerances.getkl}} for an example.
#' @import pscl
#' @export	
#'
chisqstretch.sulkl<- function(rho, n.of.x, s.of.x, trafo=(n.of.x-1)/n.of.x*s.of.x*s.of.x, norm = 1, support= c(0,Inf), log=FALSE) 
{
	require(pscl)
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
#' @title Area under the summary likelihood of the variance for normal summary values
#' @description This function computes the area under the summary likelihood \code{chisqstretch.sulkl}.
#' @inheritParams chisqstretch.sulkl
#' @seealso \code{\link{chisqstretch.calibrate}}
chisqstretch.su.lkl.norm	<- function(n.of.x, s.of.x, trafo=1, support=c(0,Inf))
{
	ans	<- integrate(chisqstretch.sulkl, lower=support[1], upper=support[2],  n.of.x=n.of.x, s.of.x=s.of.x, norm=1, trafo=trafo , support=support, log=FALSE)	
	ans$value
}
#------------------------------------------------------------------------------------------------------------------------
#' @title KL divergence between the summary likelihood and the power function of \code{chisqstretch}
#' @description Compute the Kullback-Leibler divergence between the summary likelihood and the power function of the \code{chisqstretch} equivalence test.
#' The KL divergence is required to calibrate the number of simulated data points of the test.
#' @inheritParams chisqstretch.sulkl 
#' @inheritParams chisqstretch.pow
#' @param mx.pw 			Power at the point of reference rho.star=1 (only used when \code{calibrate.tau.u==TRUE}).
#' @param alpha 			Level of the equivalence test
#' @param calibrate.tau.u	If \code{TRUE} the upper tolerance of the equivalence region (\code{tau.u}) is calibrated so that power at the point of reference is equal to \code{mx.pw}
#' @param tau.u				Upper tolerance of the equivalence region. If \code{calibrate.tau.u==TRUE}, \code{tau.u} is just a guess on an upper bound on the upper tolerance of the equivalence region to speed up calibration.
#' @param pow_scale 	 	Used to set the support of the pdf associated to the power function. The power is truncated between \code{[tau.l/pow_scale,tau.u*pow_scale]} and then standardized.
#' @param plot 				Logical. If \code{plot==TRUE}, the power of the calibrated test is plotted along with the summary likelihood.
#' @return	vector of length 6
#' 	\item{KL_div}{the Kullback Leibler divergence}	
#' 	\item{tau.l}{lower tolerance of the equivalence region}	
#' 	\item{tau.u}{upper tolerance of the equivalence region}
#' 	\item{c.l}{lower point of the critical region, i.e. lower standard ABC tolerance}	
#' 	\item{c.u}{upper point of the critical region, i.e. upper standard ABC tolerance}	
#' 	\item{pw.cmx}{actual maximum power at the point of equality}
#' @note Whatever the value of \code{calibrate.tau.u}, the lower tolerance of the equivalence region (\code{tau.l}) is always numerically calibrated so that the mode of the power function is at the point of equality rho.star.
#' @export
#' @import ggplot2 reshape2 pscl
#' @example example/ex.chisqstretch.calibrate.tolerances.getkl.R
#' 
chisqstretch.calibrate.tolerances.getkl <- function(n.of.x, s.of.x, scale, df, tau.u, mx.pw=0.9, alpha=0.01, pow_scale=1.5, calibrate.tau.u=T, plot = F) 
{
	tau.l<- pw.cmx<- error<- c.l<- c.u<- NA	
	if(calibrate.tau.u)	#calibrate tau.u constrained on yn, alpha and mx.pw 
	{			
		g(tau.l, tau.u, pw.cmx,	error, c.l, c.u)	%<-%	chisqstretch.calibrate.tauup( mx.pw, tau.u, scale, df, alpha )						#tau.u is taken as upper bound on calibrated tau.u
		if (abs(pw.cmx - mx.pw) > 0.09) 	stop("tau.up not accurate")			
	}
	else
	{
		g(tau.l, c.l, c.u, error)	%<-%	chisqstretch.calibrate.taulow(tau.u, scale, df, alpha )	#tau.u is taken as final tau.u
	}
	
	#truncate pow and compute pow_norm	
	pow_support <- c(tau.l/pow_scale, tau.u*pow_scale) 	
	pow_norm 	<- chisqstretch.pow.norm(scale, df, c.l, c.u, trafo=1, support=pow_support)
	#compute the norm of lkl, given its support 
	lkl_support	<- pow_support	
	#print(c(n.of.x, s.of.x, (n.of.x-1)/n.of.x*s.of.x*s.of.x)); print(lkl_support)
	lkl_norm	<- chisqstretch.su.lkl.norm(n.of.x, s.of.x, trafo=(n.of.x-1)/n.of.x*s.of.x*s.of.x, support=lkl_support)
	integral_range	<- pow_support			
	lkl_arg			<- list(n.of.x= n.of.x, s.of.x= s.of.x, trafo= (n.of.x-1)/n.of.x*s.of.x*s.of.x, norm = lkl_norm, support = lkl_support)
	pow_arg			<- list(scale = scale, df = df, c.l=c.l, c.u=c.u, norm=pow_norm, support=pow_support, trafo= 1)	
	tmp 			<- integrate(kl.integrand, lower = integral_range[1], upper = integral_range[2], dP=chisqstretch.sulkl, dQ=chisqstretch.pow, P_arg=lkl_arg, Q_arg=pow_arg)
	KL_div			<- tmp$value
	if (tmp$message != "OK") 
	{
		warning(tmp$message)
	}
	if (plot) 
	{
		library(ggplot2)
		library(reshape2)
		rho_lkl 			<- seq(lkl_support[1], lkl_support[2], length.out = 1000)
		lkl					<- chisqstretch.sulkl(rho_lkl, n.of.x, s.of.x, trafo= (n.of.x-1)/n.of.x*s.of.x*s.of.x, norm=lkl_norm, support=lkl_support)
		df_lkl 				<- data.frame(x = rho_lkl, yes = lkl, no = lkl*lkl_norm )
		df_lkl$distribution <- "summary likelihood"
		rho_pow	 			<- seq(pow_support[1], pow_support[2], length.out = 1000)
		pow					<- chisqstretch.pow(rho_pow, scale, df, c.l, c.u, trafo= 1, norm=pow_norm)		
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
	pw.cmx 	<- ifelse(calibrate.tau.u, pw.cmx, chisqstretch.pow(rho=1, scale, df, c.l, c.u))	
	c(KL_div = KL_div, tau.l = tau.l, tau.u = tau.u, c.l = c.l, c.u = c.u, pw.cmx = pw.cmx)	
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Calibrate the lower tolerance interval of the equivalence region for \code{chisqstretch}
#' @description This function calibrates the lower tolerance interval of the equivalence region for the \code{chisqstretch} equivalence test
#' so that the mode of the power function is at \code{rho.star=1}.
#' @export
#' @inheritParams 	chisqstretch.pow
#' @param tau.up	Upper tolerance of the equivalence region
#' @param rho.star	point of reference. Defaults to the point of equality \code{rho.star=1}.
#' @param tol		this algorithm stops when the actual point of reference is less than 'tol' from 'rho.star'
#' @param max.it	this algorithm stops prematurely when the number of iterations to find the equivalence region exceeds 'max.it'
#' @param pow_scale	Used to set the support of the power function. The power is truncated between \code{[tau.l/pow_scale,tau.u*pow_scale]} and then standardized.
#' @param verbose	Logical. If \code{verbose==TRUE}, details of the calibration are printed to the console. 
#' @return vector of length 4	
#' 	\item{tau.low}{lower tolerance of the equivalence region}	
#' 	\item{cl}{lower point of the critical region, i.e. lower standard ABC tolerance}	
#' 	\item{cu}{upper point of the critical region, i.e. upper standard ABC tolerance}	
#' 	\item{error}{actual error between the mode of the power function and rho.star}
#' @seealso \code{\link{chisqstretch.calibrate}}
#' @example example/ex.chisqstretch.calibrate.taulow.R
#' 
chisqstretch.calibrate.taulow<- function(tau.up, scale, df, alpha=0.01, rho.star=1, tol= 1e-5, max.it=100, pow_scale=1.5, verbose=0) 
{	
	rho			<- seq(1/(tau.up*pow_scale), tau.up*pow_scale, len=1024)
	tau.low.lb	<- 2/tau.up				
	tmp			<- max.it
	c.rho.max	<- Inf
	while(c.rho.max>rho.star && tmp>0)
	{
		tmp							<- tmp-1
		tau.low.lb					<- tau.low.lb/2
		rej							<- .Call("abcScaledChiSq",	c(scale,df,tau.low.lb,tau.up,alpha,1e-10,100,0.05)	)
		if(rej[4]>tol)	stop("compute tau.low.lb: rejection region does not have level alpha within tolerance")
		pw							<- chisqstretch.pow(rho,scale,df,rej[1],rej[2])
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
		rej		<- .Call("abcScaledChiSq",	c(scale,df,tau.low,tau.up,alpha,1e-10,100,0.05)	)
		if(rej[4]>tol)	stop("compute tau.low: rejection region does not have level alpha within tolerance")
		pw		<- chisqstretch.pow(rho,scale,df,rej[1],rej[2])
#print( c(rho[ which.max(pw) ],pw[ which.max(pw) ], tau.low.lb, tau.low.ub,round(tau.low.lb,digits=10)==round(tau.low.ub,digits=10) ))	
		error	<- rho[ which.max(pw) ] - rho.star
		if(verbose)	cat(paste("\ntrial tau.l=",tau.low,"pw.max is",max(pw),"at",rho[ which.max(pw) ], "it",max.it))
#print( error )			
		if(error<0)
			tau.low.lb<- tau.low
		else
			tau.low.ub<- tau.low			
	}
	if(max.it==0)	warning("chisqstretch.calibrate.taulow: reached max.it")
	c(tau.low=tau.low, cl=rej[1], cu=rej[2], error=error)
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Calibrate the upper tolerance interval of the equivalence region for \code{chisqstretch}
#' @description This function calibrates the upper tolerance interval of the equivalence region for the \code{chisqstretch} equivalence test
#' so that the mode of the power function is at \code{rho.star=1} and so that the power function at \code{rho.star} euqals \code{mx.pw}. 
#' This involves recursive recursive calls to re-calibrate the lower tolerance region.
#' @export
#' @inheritParams 	chisqstretch.pow
#' @param mx.pw		maximum power at the point of reference (rho.star).
#' @param tau.up.ub	guess on an upper bound on the upper tolerance of the equivalence region
#' @param rho.star	point of reference. Defaults to the point of equality rho.star=1.
#' @param tol		this algorithm stops when the actual maximum power is less than 'tol' from 'mx.pw'
#' @param max.it	this algorithm stops prematurely when the number of iterations to find the equivalence region exceeds 'max.it'
#' @param pow_scale	Used to set the support of the power function. The power is truncated between \code{[tau.l/pow_scale,tau.u*pow_scale]} and then standardized.
#' @return	vector of length 6
#' 	\item{tau.low}{lower tolerance of the equivalence region}		
#' 	\item{tau.up}{upper tolerance of the equivalence region}
#' 	\item{curr.mx.pw}{actual maximum power associated with the equivalence region}
#' 	\item{error}{actual error between the power at rho.star and mx.pw}
#' 	\item{cl}{lower point of the critical region, i.e. lower standard ABC tolerance}	
#' 	\item{cu}{upper point of the critical region, i.e. upper standard ABC tolerance}	
#' @example example/ex.chisqstretch.calibrate.tauup.R
#' @seealso \code{\link{chisqstretch.calibrate}}
#'
chisqstretch.calibrate.tauup<- function(mx.pw, tau.up.ub, scale, df, alpha=0.01, rho.star=1, tol= 1e-5, max.it=100, pow.scale=1.5, verbose=0)
{
	tau.low		<- cl <- cu	<- NA
	error		<- curr.mx.pw	<- 0
	tau.up.ub	<- tau.up.ub/2	
	tmp			<- max.it		
	while(curr.mx.pw<mx.pw && tmp>0)
	{
		tmp							<- tmp-1
		tau.up.ub					<- 2*tau.up.ub
		g(tau.low, cl, cu, error)	%<-%	chisqstretch.calibrate.taulow(tau.up.ub, scale, df, alpha, rho.star=rho.star, tol=tol, max.it=max.it)
		rho							<- seq(tau.low*pow.scale, tau.up.ub*pow.scale, len=1024)
		pw							<- chisqstretch.pow(rho, scale, df, cl, cu)
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
		g(tau.low, cl, cu, error)	%<-%	chisqstretch.calibrate.taulow(tau.up, scale, df, alpha, rho.star=rho.star, tol=tol, max.it=max.it)
		rho							<- seq(tau.low*pow.scale, tau.up*pow.scale, len=1024)		
		pw							<- chisqstretch.pow(rho, scale, df, cl, cu)
		curr.mx.pw					<- max(pw)		
		error						<- curr.mx.pw - mx.pw
		if(verbose)	cat(paste("\ntrial tau.u",tau.up,"with power",curr.mx.pw,"at rho=",rho[ which.max(pw) ],"search interval",tau.up.lb,tau.up.ub,"error",error))
		if(error<0)
			tau.up.lb<- tau.up
		else
			tau.up.ub<- tau.up
#print(c(abs(error), round(tau.up.lb,digits=10)!=round(tau.up.ub,digits=10)) )	
	}
	if(max.it==0)	warning("chisqstretch.calibrate.tauup: reached max.it")
	c(tau.low=tau.low, tau.up=tau.up, curr.mx.pw=curr.mx.pw,	error=abs(error), cl=cl, cu=cu)
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Calibrate the power function of \code{chisqstretch}
#' @description This function calibrates the power function of \code{chisqstretch} so that its mode coincides 
#' with the mode of the summary likelihood and so that its KL divergence to the summary likelihood is 
#' minimized. The function minimizes the KL divergences and includes recursive calls to re-calibrate the
#' upper and lower tolerance regions for every new proposed number of simulated summary values.  
#' @export
#' @inheritParams 	chisqstretch.pow
#' @inheritParams 	chisqstretch.sulkl
#' @param plot		Logical. If \code{plot==TRUE}, the calibrated power function is plotted along with the summary likelihood.
#' @param debug		Logical. If \code{debug==TRUE}, detailed optimization output for the number of simulated summary values is printed to the console.
#' @param mx.pw		maximum power at the point of equality \code{rho.star}.
#' @param max.it	this algorithm stops prematurely when the number of iterations to calibrate the number of simulated data points exceeds 'max.it'
#' @return	vector of length 7
#' 	\item{n.of.y}{number of simulated summary values}
#' 	\item{tau.l}{lower tolerance of the equivalence region}		
#' 	\item{tau.u}{upper tolerance of the equivalence region}
#' 	\item{c.l}{lower point of the critical region, i.e. lower standard ABC tolerance}	
#' 	\item{c.u}{upper point of the critical region, i.e. upper standard ABC tolerance}	
#' 	\item{pw.cmax}{actual power at rho.star}
#' 	\item{KL_div}{actual KL divergence between the ABC approximation and the summary likelihood}
#' @example example/ex.chisqstretch.calibrate.R
#' @example example/ex.chisqstretch.abcreject.R
#' @seealso \code{\link{chisqstretch.calibrate.tauup}}, \code{\link{chisqstretch.calibrate.taulow}}, \code{\link{chisqstretch.pow}}
chisqstretch.calibrate<- function(n.of.x, s.of.x, scale=n.of.x, n.of.y=n.of.x, mx.pw=0.9, alpha=0.01, max.it=100, debug=F, plot=F)
{	
	KL.of.yn_ub<- KL.of.yn<- error <- curr.mx.pw <- tau.low <- cl <- cu	<- NA		
	#KL for initial n.of.y
	KL.of.yn		<- chisqstretch.calibrate.tolerances.getkl(n.of.x, s.of.x, scale, n.of.y-1, 3*s.of.x, mx.pw=mx.pw, alpha=alpha, pow_scale=1.5, calibrate.tau.u=T, plot=F)["KL_div"]	
	#KL always decreases from n.of.x. Find upper bound yn.ub such that KL first increases again.	
	curr.it 		<- max.it
	yn.ub 			<- 2 * n.of.y		
	KL.of.yn_ub		<- chisqstretch.calibrate.tolerances.getkl(n.of.x, s.of.x, scale, yn.ub-1, 3*s.of.x, mx.pw=mx.pw, alpha=alpha, pow_scale=1.5, calibrate.tau.u=T, plot=F)["KL_div"]		
	while (KL.of.yn_ub < KL.of.yn && curr.it > 0) 
	{
		#print(c(yn.ub, KL.of.yn_ub, KL.of.yn, curr.it))
		curr.it 		<- curr.it - 1
		KL.of.yn 		<- KL.of.yn_ub
		yn.ub 			<- 2 * yn.ub
		KL.of.yn_ub		<- chisqstretch.calibrate.tolerances.getkl(n.of.x, s.of.x, scale, yn.ub-1, 3*s.of.x, mx.pw=mx.pw, alpha=alpha, pow_scale=1.5, calibrate.tau.u=T, plot=F)["KL_div"]
		if(debug)	cat(paste("\ntrial upper bound m=",yn.ub,"with KL",KL.of.yn_ub))
	}			
	if (curr.it == 0) 	stop("could not find upper bound for yn")					
	if(debug)	cat(paste("\nFound upper bound m=",yn.ub,"with KL",KL.of.yn_ub))
	yn.lb	<- ifelse(curr.it==max.it, yn.ub/2, yn.ub/4)
	if(debug)	cat(paste("\nupper and lower bounds on m:",yn.lb, yn.ub))
	
	KL_args					<- list(n.of.x=n.of.x, s.of.x=s.of.x, scale=scale, tau.u=3*s.of.x, mx.pw=mx.pw, alpha=alpha, calibrate.tau.u=T, plot=F)	
	tmp 					<- optimize(kl.optimize, interval = c(yn.lb-1, yn.ub-1), x_name = "df", is_integer = T, KL_divergence = "chisqstretch.calibrate.tolerances.getkl", KL_args = KL_args, verbose = debug, tol = 1)
	
	n.of.y 										<- round(tmp$minimum)+1
	g(KL_div, tau.l, tau.u, c.l, c.u, pw.cmx)	%<-%	chisqstretch.calibrate.tolerances.getkl(n.of.x, s.of.x, scale, n.of.y-1, 3*s.of.x, mx.pw=mx.pw, alpha=alpha, pow_scale=1.5, calibrate.tau.u=T, plot=plot)
	c(n.of.y=n.of.y, tau.l=tau.l, tau.u=tau.u, cl=c.l, cu=c.u, pw.cmx=pw.cmx, KL_div=KL_div)		
}
