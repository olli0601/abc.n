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
ftest.pow <- function(rho2, tau2, alpha = 0.01, n, p, support = c(0, Inf), log = FALSE)
{ 
	stopifnot(alpha > 0, alpha < 0.5, support[1] <= support[2], tau2 > 0, support[1] >= 0)
	stopifnot(n %% 1 == 0, p %% 1 == 0, n > 1, p > 0, p < n)
	ans 			<- rho2
	in_support 		<- (rho2 >= support[1] & rho2 <= support[2])
	ans[!in_support]<- ifelse(log, -Inf, 0)
	if(any(in_support))
	{
		crit.val            <- qf(alpha, p, n - p, n * tau2)
		ans[in_support]		<- pf(crit.val, p, n - p, n * rho2[in_support])
		tmp					<- which(ans[in_support] < 0)
		if(length(tmp))
			ans[tmp]		<- 0
		if(log)
			ans[in_support]	<- log(ans[in_support])
	}
	ans
}
#------------------------------------------------------------------------------------------------------------------------
# calculate the critical value for a given equivalence region
ftest.criticalvalue <- function(tau2, n, p, alpha = 0.01)
{				 
	c <- p * (n - 1) * qf(alpha, p, n - p, n * tau2) / (n - p)
	names(c) <- "c"
	c
}
#------------------------------------------------------------------------------------------------------------------------
# Calibrate the equivalence region for given maximum power
ftest.calibrate.tolerances <- function(mx.pw, n, p, tau2.ub, alpha = 0.01, tol = 1e-5, max.it = 100, verbose = 0)
{
	stopifnot(mx.pw > 0, mx.pw < 1, alpha > 0, alpha < 0.5, tau2.ub > 0, max.it > 10, tol < 0.01)
	stopifnot(n %% 1 == 0, p %% 1 == 0, n > 1, p > 0, p < n)
	stopifnot((mx.pw - tol) > 0, (mx.pw + tol) < 1)
	
	pw.lb <- mx.pw - tol
	pw.ub <- mx.pw + tol
	
	tau2.lb <- 0
	tmp	<- max.it
	while(ftest.pow(0, tau2.ub, alpha, n, p) < pw.ub & tmp > 0)
	{
		tau2.ub <- tau2.ub + 1
		tmp <- tmp - 1
	}
	if(tmp == 0) stop("could not set viable upper bound for tau2")	
	
	tau2 <- (tau2.ub + tau2.lb) / 2
	tmp	<- max.it
	while(ftest.pow(0, tau2, alpha, n, p) > pw.lb & tmp > 0)
	{
		tau2 <- (tau2 + tau2.lb) / 2
		tmp <- tmp - 1
	}
	if(tmp == 0) stop("could not set viable lower bound for tau2")
	else tau2.lb <- tau2
	
	if(ftest.pow(0, tau2.lb, alpha, n, p) > ftest.pow(0, tau2.ub, alpha, n, p)) stop("non-viable bounds generated")
	
	tmp	<- max.it
	while(tmp > 0 & round(tau2.lb, digits = 10) != round(tau2.ub, digits = 10))
	{
		tau2 <- seq(tau2.lb, tau2.ub, len = 1024)
		curr.pw <- sapply(tau2, function(x) ftest.pow(0, x, alpha, n, p))
		in_tol <- (curr.pw >= pw.lb & curr.pw <= pw.ub)
		if(all(!in_tol))
		{
			in_tol <- (curr.pw < pw.lb)
			if(!all(!in_tol)) tau2.lb <- tau2[max(which(in_tol))]
			in_tol <- (curr.pw > pw.ub)
			if(!all(!in_tol)) tau2.ub <- tau2[min(which(in_tol))]
		}
		else
		{   
			tau2.ub <- min(tau2[in_tol])
			min_in_tol <- which(1:length(curr.pw) < min(which(in_tol)))
			if(length(min_in_tol) > 0) tau2.lb <- tau2[max(min_in_tol)]
		}
		if(verbose)	cat(paste("\ntrial upper bound", tau2.ub, "with power", ftest.pow(0, tau2.ub, alpha, n, p), " lower bound", tau2.lb, "with power", ftest.pow(0, tau2.lb, alpha, n, p)))
		tmp <- tmp - 1
	}
	if(tmp == 0)	stop("could not find tau2")
	tau2 <- round(tau2.lb, digits = 10)
	curr.pw <- ftest.pow(0, tau2, alpha, n, p)
	if(verbose)	cat(paste("\nFound tau2", tau2, "with power", curr.pw, "at rho=0\n"))
	ans				<- c(tau2, curr.pw, abs(curr.pw - mx.pw))	
	names(ans)		<- c("tau2", "curr.pw", "error")
	c(ans, ftest.criticalvalue(tau2 = tau2, n = n, p = p, alpha = alpha))
}

#------------------------------------------------------------------------------------------------------------------------
ftest.plot <- function(n, p, tau2, alpha, pow_scale = 1.5)
{
	tmp			<- data.frame(rho2 = seq(0, pow_scale * tau2, length.out = 1024))	
	tmp$power	<- ftest.pow(tmp$rho2, tau2, alpha, n, p)
	
	p	<- ggplot(tmp, aes(x = rho2, y = power)) + geom_line() + labs(x = expression(rho^2), y = 'Power\n(ABC acceptance probability)') +
			scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
			geom_vline(xintercept = tau2, linetype = "dotted") +
			ggtitle(substitute(paste("n = ", n, ", p = ", p, ", ", tau^2, " = ", tau2, sep = ""), list(n = n, p = p, tau2 = tau2)))
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
#' @param tau2 		Squared distance that defines the equivalence region
#' @param n 		Number of replicate simulations
#' @param p 		Number of variables
#' @param what		Character string to indicate the type of calibration to be performed 
#' @param c		Upper boundary point of the critical region 
#' @param tau2		Upper boundary point of the equivalence region  
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
#'  \item (\code{what = ALPHA}) This calibration requires the inputs \code{c}, \code{tau2} with \code{0 < tau2} and \code{c > 0}. 
#' 				The output contains the corresponding ABC false positive rate \code{alpha}.
#' 				This option does not specify any of the free ABC parameters, but may be useful to determine the ABC
#' 				false positive rate for uncalibrated ABC routines.
#'  \item (\code{what = CR}) This calibration requires the input \code{tau2}, \code{alpha} with \code{tau2 > 0} and default \code{alpha = 0.01}. 
#' 				The output contains the corresponding critical value \code{c}, which defines critical region \code{[0, c]}, which corresponds to the ABC tolerance region typically denoted by \code{[0, tau2]}.
#'  \item (\code{what=MXPW}) This calibration requires the inputs \code{alpha}, \code{mx.pw}, with default values 0.01 and 0.9 respectively.
#' 				The output contains the corresponding critical value \code{c}, corresponding to critical region \code{[0, c]} (to be used in ABC). It also contains 
#' 				the corresponding squared distance \code{tau2}, defining equivalence region \code{[0, tau2]} that gives a suitable ABC accept/reject probability if the simulated summary values are close to the observed summary values.
#' 				As a check to the numerical calibrations, the actual power at the point of equality is returned (\code{pw.cmx}).
#' }
#' @example example/ex.ftest.calibrate.R
#' @references  http://arxiv.org/abs/1305.4283
ftest.calibrate<- function(n = NA, p = NA, what = 'MXPW', 
		mx.pw = 0.9, alpha = 0.01, c = NA, tau2 = NA, tol = 1e-5, 
		max.it = 100, pow_scale = 1.5, plot = FALSE, verbose = FALSE)
{
	stopifnot(what %in% c('ALPHA', 'CR', 'MXPW'))
	if(what == 'ALPHA')
	{
		stopifnot(c > 0, tau2 > 0)
		stopifnot(n %% 1 == 0, p %% 1 == 0, n > 1, p > 0, p < n)
		ans         <- pf((n - p) * c / (p * (n - 1)), p, n - p, n * tau2)
		names(ans)	<- c('alpha')
		if(plot)
			ftest.plot(n = n, p = p, tau2 = tau2, alpha = ans['alpha'], pow_scale = pow_scale)
	}
	if(what == 'CR')
	{
		stopifnot(tau2 > 0, alpha > 0, alpha < 1)
		stopifnot(n %% 1 == 0, p %% 1 == 0, n > 1, p > 0, p < n)
		ans			<- ftest.criticalvalue(tau2 = tau2, n = n, p = p, alpha = alpha)
		names(ans)	<- c('c')
		if(plot)
			ftest.plot(n = n, p = p, tau2 = tau2, alpha = alpha, pow_scale = pow_scale)
	}
	if(what == 'MXPW')
	{
		stopifnot(mx.pw > 0, mx.pw < 1, alpha > 0, alpha < 0.5, max.it > 10, tol < 0.01, pow_scale > 1)
		stopifnot(n %% 1 == 0, p %% 1 == 0, n > 1, p > 0, p < n)
		tmp			<- ftest.calibrate.tolerances(mx.pw, n = n, p = p, alpha = alpha, tau2.ub = ifelse(is.na(tau2), 1, tau2), tol = tol, max.it = max.it, verbose = verbose)
		ans			<- c(tmp[4], tmp[1], tmp[2], tmp[3])
		names(ans)	<- c('c','tau2', 'pw.cmx', 'pw.error')
		if(plot)
			ftest.plot(n = n, p = p, tau2 = ans['tau2'], alpha = alpha, pow_scale = pow_scale)		
	}	
	ans	
}