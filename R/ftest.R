#' @title \code{ftest} power function 
#' @description Compute the power of the one-sample multivariate equivalence test for population means of multivariate normal summary values with unknown population variance.
#' @export
#' @import data.table ggplot2
#' @inheritParams ftest.calibrate
#' @param rho 		Vector of quantiles
#' @param support 	Support of the truncated power function (vector of dimension 2).
#' @param log 		If \code{TRUE}, the power function is returned on the log scale. 
#' @note The power function can be truncated to \code{support}.
#' @example example/ex.ftest.pow.R
#' @references  http://arxiv.org/abs/1305.4283
ftest.pow <- function(rho, tau, n.of.y, p, alpha=0.01, support=c(0, Inf), log=FALSE, norm=1)
{ 
	stopifnot(alpha > 0, alpha < 0.5, support[1] <= support[2], tau > 0, support[1] >= 0)
	stopifnot(n.of.y %% 1 == 0, p %% 1 == 0, n.of.y > 1, p > 0, p < n.of.y)
	ans 			<- rho
	in_support 		<- (rho >= support[1] & rho <= support[2])
	ans[!in_support]<- ifelse(log, -Inf, 0)
	if(any(in_support))
	{
		crit.val            <- qf(alpha, p, n.of.y - p, n.of.y * tau)
		ans[in_support]		<- pf(crit.val, p, n.of.y - p, n.of.y * rho[in_support])/norm
		tmp					<- which(ans[in_support] < 0)
		if(length(tmp))
			ans[tmp]		<- 0
		if(log)
			ans[in_support]	<- log(ans[in_support])
	}
	ans
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Area under the \code{ftest} power function
#' @export
#' @description This function computes the area under the power function \code{ftest.pow}.
#' @inheritParams ftest.pow
#' @seealso \code{\link{ftest.pow}}
ftest.pow.norm<- function(tau, n.of.y, p, support=c(0,Inf))
{
	ans <- integrate(ftest.pow, lower=support[1], upper=support[2], tau=tau, n.of.y=n.of.y, p=p, norm=1, support=support, log=FALSE)
	ans$value
}	
#------------------------------------------------------------------------------------------------------------------------
# density of Hotelling T2 (partial likelihood) on observed data
ftest.sulkl	<- function(rho, t2.x, n.of.x, p, norm = 1, support= c(0,Inf), log=FALSE) 
{
	stopifnot(p>1, n.of.x>1)
	ans 				<- rho
	in_support 			<- (rho >= support[1] & rho <= support[2])
	ans[!in_support]	<- 0
	if(any(in_support)) 			
		ans[in_support]	<- df( t2.x * (n.of.x-p) / (p*(n.of.x-1)), df1=p, df2=n.of.x-p, ncp= n.of.x*rho[in_support]) / norm
	if(log)
		ans				<- log(ans)
	return(ans)
}
#------------------------------------------------------------------------------------------------------------------------
# normalising constant of density of Hotelling T2 (partial likelihood) on observed data
ftest.sulkl.norm	<- function(t2.x, n.of.x, p, support=c(0,Inf))
{	
	ans	<- integrate(ftest.sulkl, lower=support[1], upper=support[2],  t2.x=t2.x, n.of.x=n.of.x, p=p, norm=1, support=support, log=FALSE)	
	ans$value
}
#------------------------------------------------------------------------------------------------------------------------
# calculate the critical value for a given equivalence region
ftest.criticalvalue <- function(tau, ny, p, alpha = 0.01)
{				 
	c <- p * (ny - 1) * qf(alpha, p, ny - p, ny * tau) / (ny - p)
	names(c) <- "c"
	c
}
#------------------------------------------------------------------------------------------------------------------------
# Calibrate the equivalence region for given maximum power
ftest.calibrate.tau <- function(mx.pw, ny, p, tau.ub, alpha = 0.01, tol = 1e-5, max.it = 100, verbose = 0, use.R = 0)
{
	stopifnot(mx.pw > 0, mx.pw < 1, alpha > 0, alpha < 0.5, tau.ub > 0, max.it > 10, tol < 0.01)
	stopifnot(ny %% 1 == 0, p %% 1 == 0, ny > 1, p > 0, p < ny)
	stopifnot((mx.pw - tol) > 0, (mx.pw + tol) < 1)
	if(!use.R)
	{
		cali.args	<- list(mx.pw=mx.pw, ny=ny, p=p, tau.ub=tau.ub, alpha=alpha, tol=tol, max.it=max.it)	
		ans			<- .Call("abcFTEST_calibrate_tau_for_mxpw", cali.args)
	}
	if(use.R)
	{
		pw.lb <- mx.pw - tol
		pw.ub <- mx.pw + tol
		
		tau.lb <- 0
		tmp	<- max.it
		while(ftest.pow(0, tau.ub, ny, p, alpha) < pw.ub & tmp > 0)
		{
			tau.ub <- tau.ub + 1
			tmp <- tmp - 1
		}
		if(tmp == 0) stop("could not set viable upper bound for tau")	

			tau <- (tau.ub + tau.lb) / 2
		tmp	<- max.it
		while(ftest.pow(0, tau, ny, p, alpha) > pw.lb & tmp > 0)
		{
			tau <- (tau + tau.lb) / 2
			tmp <- tmp - 1
		}
		if(tmp == 0) stop("could not set viable lower bound for tau")
			else tau.lb <- tau
		
		if(ftest.pow(0, tau.lb, ny, p, alpha) > ftest.pow(0, tau.ub, ny, p, alpha)) stop("non-viable bounds generated")

			tmp	<- max.it
		while(tmp > 0 & round(tau.lb, digits = 10) != round(tau.ub, digits = 10))
		{
			tau <- seq(tau.lb, tau.ub, len = 1024)
			curr.pw <- sapply(tau, function(x) ftest.pow(0, x, ny, p, alpha))
			in_tol <- (curr.pw >= pw.lb & curr.pw <= pw.ub)
			if(all(!in_tol))
			{
				in_tol <- (curr.pw < pw.lb)
				if(!all(!in_tol)) tau.lb <- tau[max(which(in_tol))]
				in_tol <- (curr.pw > pw.ub)
				if(!all(!in_tol)) tau.ub <- tau[min(which(in_tol))]
			}
			else
			{   
				tau.ub <- min(tau[in_tol])
				min_in_tol <- which(1:length(curr.pw) < min(which(in_tol)))
				if(length(min_in_tol) > 0) tau.lb <- tau[max(min_in_tol)]
			}
			if(verbose)	cat(paste("\ntrial upper bound", tau.ub, "with power", ftest.pow(0, tau.ub, ny, p, alpha), " lower bound", tau.lb, "with power", ftest.pow(0, tau.lb, alpha, ny, p)))
				tmp <- tmp - 1
		}
		if(tmp == 0)	stop("could not find tau")
			tau <- round(tau.lb, digits = 10)
		curr.pw <- ftest.pow(0, tau, ny, p, alpha)
		if(verbose)	cat(paste("\nFound tau", tau, "with power", curr.pw, "at rho=0\n"))
			ans				<- c(ftest.criticalvalue(tau=tau, ny=ny, p=p, alpha=alpha), tau, curr.pw, abs(curr.pw - mx.pw))	
		names(ans)		<- c("c", "tau", "curr.pw", "error.pw")	
	}	
	ans
}
#------------------------------------------------------------------------------------------------------------------------
#	compute KL divergence between power of F test and summary likelihood
ftest.getkl <- function(n.of.x, t2.x, n.of.y, p, tau.u, mx.pw=0.9, alpha=0.01, pow_scale=1.5, plot = F) 
{	
	#	calibrate tau	
	tmp				<- ftest.calibrate.tau( mx.pw, n.of.y, p, tau.u, alpha )						#tau.u is taken as upper bound on calibrated tau.u
	crit			<- tmp['c']
	tau				<- tmp['tau']
	curr.pw			<- tmp['curr.pw']	
	if (abs(curr.pw - mx.pw) > 0.09) 	
		stop("tau.up not accurate")			
	#	truncate pow and compute the normalizing constant pow_norm	
	pow_support 	<- c(0, tau*pow_scale) 	
	pow_norm 		<- ftest.pow.norm(tau, n.of.y, p, support=pow_support)
	#	truncate lkl and compute the normalizing constant lkl_norm 
	lkl_support		<- pow_support	
	#print(c(n.of.x)); print(lkl_support)
	lkl_norm		<- ftest.sulkl.norm(t2.x, n.of.x, p, support=lkl_support)
	integral_range	<- pow_support			
	lkl_arg			<- list(n.of.x=n.of.x, t2.x=t2.x, p=p, norm=lkl_norm, support=lkl_support)
	pow_arg			<- list(tau=tau, n.of.y=n.of.y, p=p, alpha=alpha, norm=pow_norm, support=pow_support)	
	tmp 			<- integrate(kl.integrand, lower=integral_range[1], upper=integral_range[2], dP=ftest.sulkl, dQ=ftest.pow, P_arg=lkl_arg, Q_arg=pow_arg)
	KL_div			<- tmp$value
	if (tmp$message != "OK") 
	{
		warning(tmp$message)
	}
	if (plot) 
	{
		rho_lkl 			<- seq(lkl_support[1], lkl_support[2], length.out = 1000)
		lkl					<- ftest.sulkl(rho_lkl, t2.x, n.of.x, p, norm=lkl_norm, support=lkl_support)
		df_lkl 				<- data.table(X = rho_lkl, density = lkl, power = lkl*lkl_norm, TYPE="summary likelihood" )
		rho_pow	 			<- seq(pow_support[1], pow_support[2], length.out = 1000)
		pow					<- ftest.pow(rho_pow, tau, n.of.y, p, alpha=alpha, support=pow_support, norm=pow_norm)		
		df_pow 				<- data.table(X = rho_pow, density = pow, power = pow*pow_norm, TYPE="ABC approximation")
		gdf 				<- rbind(df_pow, df_lkl)
		gdf					<- melt(gdf, id.vars=c("X","TYPE"))
		gdf					<- subset(gdf, !(TYPE=='summary likelihood' & variable=='power'))
		set(gdf, NULL, 'TYPE', gdf[, factor(TYPE, levels=c('summary likelihood',"ABC approximation"), labels=c('summary likelihood',"ABC approximation"))])
		pp	<- ggplot(gdf, aes(x=X, y=value, group=TYPE, colour=TYPE)) +
		geom_ribbon(data=subset(gdf, TYPE=='summary likelihood'), aes(ymax=value, ymin=0), fill='grey70', guide=FALSE) +				
		geom_vline(xintercept= tau,linetype="dotted") + 
		geom_hline(yintercept= mx.pw,linetype="dotted") + 
		geom_line() +
		scale_y_continuous() +				 
		scale_colour_manual(values=c('black','grey70')) + 
		labs(x= expression(rho), y="", linetype="Normalized", colour='') +
		facet_wrap(~variable, scales='free') +
		theme_bw() + theme(legend.position='bottom') #+ guides(colour=guide_legend(ncol=2))
		#p 					<- p + ggtitle(paste("n.of.y=", df+1, "\ntau.l=", tau.l,"\ntau.u=", tau.u,"\nKL=", KL_div))
		print(pp)
	}
	tmp			<- c(KL_div, tau, crit, curr.pw,abs(curr.pw - mx.pw))
	names(tmp)	<- c('KL_div', 'tau', 'c', 'pw.cmx', 'err.pw')
	tmp
}
#------------------------------------------------------------------------------------------------------------------------
#	calibrate Ftest to match summary likelihood
ftest.calibrate.kl<- function(t2.x, n.of.x, p, n.of.y=n.of.x, mx.pw=0.9, alpha=0.01, max.it=100, use.R= FALSE, debug=FALSE, plot=FALSE)
{				
	if(!use.R)
	{
		KL.args		<- list(n.of.x=n.of.x, t2.x=t2.x, p=p, n.of.y=n.of.y, mx.pw=mx.pw, alpha=alpha)
		ans			<- .Call("abcFTEST_calibrate_KL", KL.args, as.integer(max.it))
	}
	if(use.R)
	{
		#	KL for initial n.of.y
		KL.of.yn		<- ftest.getkl(n.of.x, t2.x, n.of.y, p, 4*t2.x, mx.pw=mx.pw, alpha=alpha, pow_scale=1.5, plot=F)["KL_div"]		
		#	KL always decreases from n.of.x. Find upper bound yn.ub such that KL first increases again.	
		curr.it 		<- max.it
		yn.ub 			<- 2 * n.of.y		
		KL.of.yn_ub		<- ftest.getkl(n.of.x, t2.x, yn.ub, p, 4*t2.x, mx.pw=mx.pw, alpha=alpha, pow_scale=1.5, plot=F)["KL_div"]
		if(debug)
			cat(paste('\nKL for n.of.y',KL.of.yn,'KL for 2 x n.of.y',KL.of.yn_ub))
		while (KL.of.yn_ub < KL.of.yn && curr.it > 0) 
		{
			#print(c(yn.ub, KL.of.yn_ub, KL.of.yn, curr.it))
			curr.it 		<- curr.it - 1
			KL.of.yn 		<- KL.of.yn_ub
			yn.ub 			<- 2 * yn.ub
			KL.of.yn_ub		<- ftest.getkl(n.of.x, t2.x, yn.ub, p, 4*t2.x,  mx.pw=mx.pw, alpha=alpha, pow_scale=1.5, plot=F)["KL_div"]
			if(debug)	cat(paste("\ntrial upper bound n.of.y=",yn.ub,"with KL",KL.of.yn_ub))
		}			
	if (curr.it == 0) 	stop("could not find upper bound for n.of.y")					
		if(debug)	
			cat(paste("\nFound upper bound n.of.y=",yn.ub,"with KL",KL.of.yn_ub))
		yn.lb			<- ifelse(curr.it==max.it, yn.ub/2, yn.ub/4)
		if(debug)	
			cat(paste("\nupper and lower bounds on n.of.y:",yn.lb, yn.ub))
		
		KL_args				<- list(	n.of.x=n.of.x, t2.x=t2.x, p=p, tau.u=4*t2.x, mx.pw=mx.pw, alpha=alpha, plot=F)	
		tmp 				<- optimize(kl.optimize, interval=c(yn.lb-1, yn.ub-1), x_name="n.of.y", is_integer=TRUE, KL_divergence="ftest.getkl", KL_args=KL_args, verbose=debug, tol=1)	
		n.of.y 				<- round(tmp$minimum)+1		
		tmp					<- ftest.getkl(n.of.x, t2.x, n.of.y, p, 4*t2.x, mx.pw=mx.pw, alpha=alpha, pow_scale=1.5, plot=plot)
		ans					<- c(n.of.y=n.of.y, tau=tmp['tau'], c=tmp['c'], pw.cmx=tmp['pw.cmx'], KL_div=tmp['KL_div'])	
	}			
	ans
}
#------------------------------------------------------------------------------------------------------------------------
ftest.plot <- function(n.of.y, p, tau, alpha, pow_scale = 1.5)
{
	tmp			<- data.frame(rho = seq(0, pow_scale * tau, length.out = 1024))	
	tmp$power	<- ftest.pow(tmp$rho, tau, n.of.y, p, alpha)
	
	p	<- ggplot(tmp, aes(x = rho, y = power)) + geom_line() + labs(x = expression(rho), y = 'Power\n(ABC acceptance probability)') +
	scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
	geom_vline(xintercept = tau, linetype = "dotted") +
	ggtitle(substitute(paste("n = ", n.of.y, ", p = ", p, ", ", " tau = ", tau, sep = ""), list(n.of.y = n.of.y, p = p, tau = tau)))
	print(p)
}

#------------------------------------------------------------------------------------------------------------------------
#' @title Calibrating \code{ftest} for ABC
#' @description Calibrate the one-sample equivalence test for population means of multivariate normal summary values with unknown population covariance matrix for ABC inference.
#' 
#' Different types of calibrations are available, see Notes for details:
#' \enumerate{ 
#'  \item (\code{what = ALPHA}) compute the ABC false positive rate for given critical region,
#'  \item (\code{what = CR}) calibrate the critical region for given ABC false positive rate,
#'  \item (\code{what = MXPW}) calibrate the critical region and the equivalence region for given ABC false positive rate and maximum power.
#' }
#' 
#' Depending on the type of calibration, some of the following inputs must be specified (see Examples). 
#' @export 
#' @import data.table ggplot2
#' @param n 		Number of replicate simulations
#' @param p 		Number of variables
#' @param what		Character string to indicate the type of calibration to be performed 
#' @param c		Upper boundary point of the critical region 
#' @param tau		Upper boundary point of the equivalence region  
#' @param mx.pw 	Maximum power at the point of equality
#' @param alpha 	Level of the equivalence test
#' @param max.it  	Maximum number of optimization steps at each calibration hierarchy
#' @param tol		Required error tolerance in calibrating the actual maximum power to the requested maximum power 
#' @param plot  	Flag to plot calibrations
#' @param pow_scale Scale for the support of the power function (used for plotting)
#' @param verbose	Flag to run in verbose mode
#' @return	vector
#' @seealso \code{\link{mutost.calibrate}}, \code{\link{vartest.calibrate}}, \code{\link{ratetest.calibrate}}, \code{\link{ztest.calibrate}}
#' @note 
#' \enumerate{	
#'  \item (\code{what = ALPHA}) This calibration requires the inputs \code{c}, \code{tau} with \code{0 < tau} and \code{c > 0}. 
#' 				The output contains the corresponding ABC false positive rate \code{alpha}.
#' 				This option does not specify any of the free ABC parameters, but may be useful to determine the ABC
#' 				false positive rate for uncalibrated ABC routines.
#'  \item (\code{what = CR}) This calibration requires the input \code{tau}, \code{alpha} with \code{tau > 0} and default \code{alpha = 0.01}. 
#' 				The output contains the corresponding critical value \code{c}, which defines critical region \code{[0, c]}, which corresponds to the ABC tolerance region typically denoted by \code{[0, tau]}.
#'  \item (\code{what=MXPW}) This calibration requires the inputs \code{alpha}, \code{mx.pw}, with default values 0.01 and 0.9 respectively.
#' 				The output contains the corresponding critical value \code{c}, corresponding to critical region \code{[0, c]} (to be used in ABC). It also contains 
#' 				the corresponding squared distance \code{tau}, defining equivalence region \code{[0, tau]} that gives a suitable ABC accept/reject probability if the simulated summary values are close to the observed summary values.
#' 				As a check to the numerical calibrations, the actual power at the point of equality is returned (\code{pw.cmx}).
#' }
#' @example example/ex.ftest.calibrate.R
#' @references  http://arxiv.org/abs/1305.4283
ftest.calibrate<- function(n.of.x=NA, t2.x=NA, n.of.y = NA, p = NA, what = 'MXPW', 
	mx.pw = 0.9, alpha = 0.01, c = NA, tau = NA, tol = 1e-5, 
	max.it = 100, pow_scale = 1.5, use.R= FALSE, debug=FALSE, plot=FALSE, verbose=FALSE)
{
	stopifnot(what %in% c('ALPHA', 'CR', 'MXPW', 'KL'))
	if(what == 'ALPHA')
	{
		stopifnot(c > 0, tau > 0)
		stopifnot(n.of.y %% 1 == 0, p %% 1 == 0, n.of.y > 1, p > 0, p < n.of.y)
		ans         <- pf((n.of.y - p) * c / (p * (n.of.y - 1)), p, n.of.y - p, n.of.y * tau)
		names(ans)	<- c('alpha')
		if(plot)
			ftest.plot(n.of.y = n.of.y, p = p, tau = tau, alpha = ans['alpha'], pow_scale = pow_scale)
	}
	if(what == 'CR')
	{
		stopifnot(tau > 0, alpha > 0, alpha < 1)
		stopifnot(n.of.y %% 1 == 0, p %% 1 == 0, n.of.y > 1, p > 0, p < n.of.y)
		ans			<- ftest.criticalvalue(tau = tau, n.of.y = n.of.y, p = p, alpha = alpha)
		names(ans)	<- c('c')
		if(plot)
			ftest.plot(n.of.y = n.of.y, p = p, tau = tau, alpha = alpha, pow_scale = pow_scale)
	}
	if(what == 'MXPW')
	{
		stopifnot(mx.pw > 0, mx.pw < 1, alpha > 0, alpha < 0.5, max.it > 10, tol < 0.01, pow_scale > 1)
		stopifnot(n.of.y %% 1 == 0, p %% 1 == 0, n.of.y > 1, p > 0, p < n.of.y)
		tmp			<- ftest.calibrate.tau(mx.pw, ny = n.of.y, p = p, tau.ub = ifelse(is.na(tau), 1, tau), alpha = alpha, tol = tol, max.it = max.it, use.R = use.R, verbose = verbose)
		ans			<- c(tmp['c'], tmp['tau'], tmp['curr.pw'], tmp['error.pw'])
		names(ans)	<- c('c','tau', 'pw.cmx', 'pw.error')
		if(plot)
			ftest.plot(n.of.y=n.of.y, p = p, tau = ans['tau'], alpha = alpha, pow_scale = pow_scale)		
	}	
	if(what == 'KL')
	{
		stopifnot(!is.na(t2.x), n.of.x>2, p>1, p<n.of.x, alpha>0, alpha<1, max.it > 10, tol < 0.01, pow_scale > 1)
		stopifnot(n.of.x %% 1 == 0, p %% 1 == 0)		
		tmp			<- ftest.calibrate.kl(t2.x, n.of.x, p, n.of.y=n.of.x, mx.pw=mx.pw, alpha=alpha, max.it=max.it, use.R=use.R, debug=debug, plot=plot)
		ans			<- c(tmp[3], tmp[2], tmp[1], tmp[4], tmp[5])
		names(ans)	<- c('c','tau','n.of.y','pw.cmx','KL.div')
	}
	ans	
}