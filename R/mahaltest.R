mahaltest.plot <- function(p, df, c.l, c.u, tau.l, tau.u, pow_scale = 1.5)
{
	pow_support <- c(tau.l / pow_scale, tau.u * pow_scale) 	
	pow_norm 	<- vartest.pow.norm(1, df, c.l, c.u, trafo = 1, support = pow_support)	
	tmp			<- data.frame(rho=seq(pow_support[1], pow_support[2], length.out = 1024))	
	tmp$power	<- vartest.pow(tmp$rho, 1, df, c.l, c.u, trafo = 1, norm = pow_norm) * pow_norm	
	
	p	<- ggplot(tmp, aes(x = rho, y = power)) + geom_line() + labs(x = expression(rho), y = 'Power\n(ABC acceptance probability)') +
			scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
			scale_x_continuous(limits = c(0, pow_support[2])) +
			geom_vline(xintercept = c(tau.l, tau.u), linetype = "dotted") +
			geom_vline(xintercept = c(c.l, c.u), linetype = "dashed") +
			ggtitle(paste0("p=", p, ", n.of.y=", df + p, "\ntau.l=", round(tau.l, d = 5), " tau.u=", round(tau.u, d = 5), "\nc.l=", round(c.l, d = 5), " c.u=", round(c.u, d = 5)))
	print(p)
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Calibrating \code{mahaltest} for ABC
#' @description Calibrate the one-sample equivalence test for Mahalanobis distance of 
#' multivariate normal summary values for ABC inference.
#' The one-sample \code{mahaltest} can be used to test the null hypothesis that the underlying population means of the simulated summary values 
#' are similar to the observed summary values, given an unknown covariance matrix. It is applicable when the simulated and observed summary values follow a 
#' multivariate normal distribution, or when normality cannot be rejected.
#' 
#' Different types of calibrations are available, see Notes for details:
#' \enumerate{ 
#'  \item (\code{what=ALPHA}) compute the ABC false positive rate for given critical region,
#'  \item (\code{what=CR}) calibrate the critical region for given ABC false positive rate,
#'  \item (\code{what=MXPW}) calibrate the critical region and the equivalence region for given ABC false positive rate and maximum power.
#' }
#' 
#' Depending on the type of calibration, some of the following inputs must be specified (see Examples).
#' @export 
#' @import data.table pscl reshape2 ggplot2 ash nortest
#' @param n.of.x 	Number of observations
#' @param p      	Dimensionality of the underlying MVN distribution
#' @param n.of.y 	Number of simulations
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
#' @param pow_scale Scale over which to calibrate power function and summary likelihood
#' @param plot  	Flag to plot calibrations
#' @param debug		Flag to switch off C implementation
#' @param plot_debug	Flag to plot at each calibration iteration
#' @param verbose	Flag to run in verbose mode
#' @return	vector
#' @seealso \code{\link{mutost.calibrate}}, \code{\link{ztest.calibrate}}, \code{\link{ratetest.calibrate}}, \code{\link{vartest.calibrate}}
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
#'              Sometimes the default initial value for \code{n.of.y} returns an error. In this case the argument \code{n.of.y} can be used to set an initial value (start as close to \code{p + 2} as possible).
#' 				The output consists of the corresponding critical region \code{[c.l, c.u]} (to be used in ABC, see Notes on (2)), the equivalence
#' 				region \code{[tau.l, tau.u]}, and the number of simulated summary values needed (\code{n.of.y}). As a check to the numerical calibrations, 
#' 				the KL divergence is returned (\code{KL}). It is desirable to compare the power to the summary likelihood in terms of the KL divergence, see References.
#' }
#' Note that the underlying test statistic only depends on \code{n.of.x} and \code{p}, which are 
#' known before ABC is run. Consequently, the free ABC parameters are calibrated once, before ABC is started. 
#' 
#' The lower boundary point of the critical region \code{c.l} is calibrated numerically, so that
#' the power is maximized at the point of equality \code{rho=p}. The calibrated \code{c.l} does not equal 1/\code{c.u}.
#' @references  http://arxiv.org/abs/1305.4283
mahaltest.calibrate <- function(n.of.x = NA, p = NA, n.of.y = NA, what = 'MXPW', mx.pw = 0.9, alpha = 0.01, tau.l = NA, tau.u = NA, tau.u.ub = NA, c.l = NA, c.u = NA, max.it = 100, tol = 1e-5, pow_scale = 1.5, debug = FALSE, plot = FALSE, verbose = FALSE)
{
	stopifnot(what %in% c('ALPHA', 'CR', 'MXPW', 'KL'))
	if(what == 'ALPHA')
	{
	    stop("Currently not implemented for what = 'ALPHA'")
#		stopifnot(c.u < tau.u, tau.u > 1, tau.l < 1, c.l > tau.l, alpha > 0, alpha < 1, p > 1, n.of.x > 0)
#		ans	<- pchisq(n.of.x * c.u / tau.u, p) - pchisq(n.of.x * c.l / tau.u, p)
#		names(ans)	<- 'alpha'
#		if(plot)
#			mahaltest.plot(n.of.x, p, c.l, c.u, tau.l, tau.u, pow_scale = pow_scale)
	}
	if(what == 'CR')
	{	
	    stop("Currently not implemented for what = 'CR'")
		stopifnot(scale>0, tau.u>1, tau.l<1, n.of.y>2, alpha>0, alpha<1, pow_scale>1)
#		tmp	<- .Call("abcScaledChiSq",	c(scale, n.of.y-df, tau.l, tau.u, alpha, 1e-10, max.it, 0.05)	)
#		ans	<- c(tmp[1], tmp[2], tau.l, tau.u, tmp[3], tmp[4], tmp[5])
#		names(ans)<- c('c.l','c.u','tau.l','tau.u','mx.cpw','alpha.err','n.it')
#		if(plot)
#			mahaltest.plot(scale, n.of.y-df, ans['c.l'], ans['c.u'], tau.l, tau.u, pow_scale=pow_scale)		
	}
	if(what == 'MXPW')
	{
    	print("NEED TO AMEND THIS FOR F-LIKELIHOOD")
		stopifnot(n.of.y %% 1 == 0, p %% 1 == 0, (n.of.y - p) > 0, tau.u.ub > (p - 2), p > 1, alpha > 0, alpha < 0.5, pow_scale > 1, max.it > 10, tol < 0.2, mx.pw > 0, mx.pw < 1)
		if(p > 2)
		{
		    tmp <- vartest.calibrate.tauup(mx.pw, tau.u.ub, 1, n.of.y - p, alpha = alpha, rho.star = p - 2, tol = tol, max.it = max.it, pow.scale = pow_scale, verbose = verbose)
    		ans	<- c(tmp[5], tmp[6], tmp[1], tmp[2], tmp[3], tmp[4])
    		names(ans) <- c('c.l', 'c.u', 'tau.l', 'tau.u', 'pw.cmx', 'pw.error')
		    if(plot) mahaltest.plot(p, n.of.y - p, ans['c.l'], ans['c.u'], ans['tau.l'], ans['tau.u'], pow_scale = pow_scale)
        }
		else 
		{
		    tmp			<- ftest.calibrate.tau(mx.pw, ny = n.of.y, p = p, tau.ub = tau.u.ub, alpha = alpha, tol = tol, max.it = max.it, use.R = F, verbose = verbose)
    		ans			<- c(tmp['c'], tmp['tau'], tmp['curr.pw'], tmp['error.pw'])
    		names(ans)	<- c('c','tau', 'pw.cmx', 'pw.error')
		    if(plot)
			    ftest.plot(n.of.y=n.of.y, p = p, tau = ans['tau'], alpha = alpha, pow_scale = pow_scale)
		}
	}
	if(what == 'KL')
	{
		stopifnot(p > 1, alpha > 0, alpha < 1, pow_scale > 1, max.it > 10, tol < 0.2)
		ans	<- mahaltest.calibrate.kl(n.of.x, p, n.of.y = n.of.y, mx.pw = mx.pw, alpha = alpha, max.it = max.it, debug = debug, plot = plot, pow_scale = pow_scale, tol = tol)
	}
	ans
}
#------------------------------------------------------------------------------------------------------------------------
# @title Calibrate the power function of \code{mahaltest}
# @description This function calibrates the power function of \code{mahaltest} so that its mode coincides 
# with the mode of the summary likelihood and so that its KL divergence to the summary likelihood is 
# minimized. The function minimizes the KL divergences and includes recursive calls to re-calibrate the
# upper and lower tolerance regions for every new proposed number of simulated summary values.
mahaltest.calibrate.kl <- function(n.of.x, p, n.of.y = p + 2, mx.pw = 0.9, alpha = 0.01, max.it = 100, pow_scale = 1.5, debug = F, plot = F, tol = 1e-5)
{	
    if(!is.na(n.of.x) & n.of.x <= p) {
#        print("Since n < p, must calibrate according to chi-squared summary likelihood")
#        n.of.x <- NA
        stop("n must be greater than p")
    } else {
        if(is.na(n.of.x)) print("Calibrating according to chi-squared summary likelihood")
        else print("Calibrating according to an F-distributed summary likelihood")
    }
    if(p == 2) print("Calibrating according to an F-test statistic")
    else print("Calibrating according to a chi-squared statistic")
    
	KL.of.yn_ub <- KL.of.yn <- error <- curr.mx.pw <- tau.low <- cl <- cu <- NA	
	if(is.na(n.of.y)) n.of.y <- p + 2 #ifelse(is.na(n.of.x), p + 2, n.of.x)
#	if(n.of.y < n.of.x) stop("n.of.y < n.of.x")
	#set rho.star based on degrees-of-freedom and form of summary likelihood
	if(p > 2) {
	    rho.star <- ifelse(is.na(n.of.x), p - 2, (p - 2) * (n.of.x - p) / (p * (n.of.x - p + 2)))
	} else rho.star <- 0
	#KL for initial n.of.y
	tau.u <- 6 * (rho.star + 3)
	
	KL.of.yn		<- mahaltest.getkl(n.of.x, p, n.of.y, rho.star, tau.u, mx.pw = mx.pw, alpha = alpha, pow_scale = pow_scale, calibrate.tau.u = T, plot = F, max.it = max.it, tol = tol)["KL_div"]
	KL.of.yn_ub		<- mahaltest.getkl(n.of.x, p, n.of.y + 1, rho.star, tau.u, mx.pw = mx.pw, alpha = alpha, pow_scale = pow_scale, calibrate.tau.u = T, plot = F, max.it = max.it, tol = tol)["KL_div"]
	print("Do we need this additional check for decreasing KL?")	
	if(KL.of.yn_ub < KL.of.yn)
	{
	    #KL always decreases from n.of.x. Find upper bound yn.ub such that KL first increases again.
	    curr.it 		<- max.it
	    yn.ub 			<- 2 * n.of.y		
	    KL.of.yn_ub		<- mahaltest.getkl(n.of.x, p, yn.ub, rho.star, tau.u, mx.pw = mx.pw, alpha = alpha, pow_scale = pow_scale, calibrate.tau.u = T, plot = F, max.it = max.it, tol = tol)["KL_div"]		
	    while (KL.of.yn_ub < KL.of.yn && curr.it > 0) 
	    {
		    curr.it 		<- curr.it - 1
		    KL.of.yn 		<- KL.of.yn_ub
		    yn.ub 			<- 2 * yn.ub
		    KL.of.yn_ub		<- mahaltest.getkl(n.of.x, p, yn.ub, rho.star, tau.u, mx.pw = mx.pw, alpha = alpha, pow_scale = pow_scale, calibrate.tau.u = T, plot = F, max.it = max.it, tol = tol)["KL_div"]
		    if(debug)	cat(paste("\ntrial upper bound m=", yn.ub, "with KL", KL.of.yn_ub))
	    }			
	    if (curr.it == 0) 	stop("could not find upper bound for yn")					
	    if(debug)	cat(paste("\nFound upper bound m=", yn.ub, "with KL", KL.of.yn_ub))
	    yn.lb	<- ifelse(curr.it == max.it, yn.ub / 2, yn.ub / 4)
	    if(debug)	cat(paste("\nupper and lower bounds on m:", yn.lb, yn.ub))
	
	    KL_args					<- list(n.of.x = n.of.x, p = p, rho.star = rho.star, tau.u = tau.u, mx.pw = mx.pw, alpha = alpha, calibrate.tau.u = T, plot = F, max.it = max.it, tol = tol)	
	    tmp 					<- optimize(kl.optimize, interval = c(yn.lb - 1, yn.ub - 1), x_name = "n.of.y", is_integer = T, KL_divergence = "mahaltest.getkl", KL_args = KL_args, verbose = debug, tol = 1)
	
	    n.of.y 										<- round(tmp$minimum) + 1
	} else print("KL doesn't decrease")
	if(p == 2)
	{
	    g(KL_div, tau, c1, pw.cmx, err.pw)	%<-%	mahaltest.getkl(n.of.x, p, n.of.y, rho.star, tau.u, mx.pw = mx.pw, alpha = alpha, pow_scale = pow_scale, calibrate.tau.u = T, plot = plot, max.it = max.it, tol = tol)
	    ans <- c(c = c1, tau = tau, n.of.y = n.of.y, pw.cmx = pw.cmx, KL_div = KL_div)
	}
	else
	{
	    g(KL_div, tau.l, tau.u, c.l, c.u, pw.cmx)	%<-%	mahaltest.getkl(n.of.x, p, n.of.y, rho.star, tau.u, mx.pw = mx.pw, alpha = alpha, pow_scale = pow_scale, calibrate.tau.u = T, plot = plot, max.it = max.it, tol = tol)
	    ans <- c(cl = c.l, cu = c.u, tau.l = tau.l, tau.u = tau.u, n.of.y = n.of.y, pw.cmx = pw.cmx, KL_div = KL_div)
	}
	ans
}
#------------------------------------------------------------------------------------------------------------------------
# @title KL divergence between the summary likelihood and the power function of \code{mahaltest}
# @description Compute the Kullback-Leibler divergence between the summary likelihood and the power function of the \code{mahaltest} equivalence test.
# The KL divergence is required to calibrate the number of simulated data points of the test.
# @param n.of.x         The number of observations.
# @param p              The dimensionality of the underlying MVN distribution.
# @param n.of.y         The number of simulations.
# @param rho.star       The location of the mode of the underlying summary likelihood.
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
mahaltest.getkl <- function(n.of.x, p, n.of.y, rho.star, tau.u, mx.pw = 0.9, alpha = 0.01, pow_scale = 1.5, calibrate.tau.u = T, plot = F, max.it = 100, tol = 1e-5) 
{
    #if rho.star = 0 then use F-statistic, else use chi-squared statistic
    #if n.of.x = NA then use chi-squared summary likelihood, else use F-dist summary likelihood
    if(rho.star == 0) ans <- mahaltest.getkl.F(n.of.x = n.of.x, n.of.y = n.of.y, p = p, tau.u = tau.u, mx.pw = mx.pw, alpha = alpha, pow_scale = pow_scale, plot = plot)
    else ans <- mahaltest.getkl.chi(n.of.x = n.of.x, p = p, n.of.y = n.of.y, rho.star = rho.star, tau.u = tau.u, mx.pw = mx.pw, alpha = alpha, pow_scale = pow_scale, calibrate.tau.u = calibrate.tau.u, plot = plot, max.it = max.it, tol = tol)
    ans
} 
	
#calibrate according to KL divergence using a chi-squared test statistic
mahaltest.getkl.chi <- function(n.of.x, p, n.of.y, rho.star, tau.u, mx.pw = 0.9, alpha = 0.01, pow_scale = 1.5, calibrate.tau.u = T, plot = F, max.it = 100, tol = 1e-5) 
{
	tau.l <- pw.cmx <- error <- c.l <- c.u <- NA
	#set parameters for calibration
	df <- n.of.y - p
	scale <- 1
	stopifnot(df > 0)
	if(calibrate.tau.u)	#calibrate tau.u constrained on yn, alpha and mx.pw 
	{			
		g(tau.l, tau.u, pw.cmx,	error, c.l, c.u)	%<-%	vartest.calibrate.tauup( mx.pw, tau.u, scale, df, alpha, rho.star = rho.star, max.it = max.it, tol = tol )						#tau.u is taken as upper bound on calibrated tau.u
#		print(c(tau.l, tau.u, pw.cmx, mx.pw))
		if (abs(pw.cmx - mx.pw) > 0.09) 	stop("tau.up not accurate")			
	}
	else
	{
		g(tau.l, c.l, c.u, error)	%<-%	vartest.calibrate.taulow(tau.u, scale, df, alpha, rho.star = rho.star, max.it = max.it, tol = tol )	#tau.u is taken as final tau.u
	}
	
	#truncate pow and compute pow_norm	
	pow_support 	<- c(tau.l / pow_scale, tau.u * pow_scale) 	
	pow_norm 		<- vartest.pow.norm(scale, df, c.l, c.u, trafo = 1, support = pow_support)
	#compute the norm of lkl, given its support 
	lkl_support		<- pow_support
	lkl_norm		<- mahaltest.sulkl.norm(n.of.x, p, support = lkl_support)
	integral_range	<- pow_support			
	lkl_arg			<- list(n.of.x = n.of.x, p = p, norm = lkl_norm, support = lkl_support)
	pow_arg			<- list(scale = scale, df = df, c.l = c.l, c.u = c.u, norm = pow_norm, support = pow_support, trafo = 1)	
	tmp 			<- integrate(kl.integrand, lower = integral_range[1], upper = integral_range[2], dP = mahaltest.sulkl, dQ = vartest.pow, P_arg = lkl_arg, Q_arg = pow_arg)
	KL_div			<- tmp$value
	if (tmp$message != "OK") 
	{
		warning(tmp$message)
	}
	if (plot) 
	{
		rho_lkl 			<- seq(lkl_support[1], lkl_support[2], length.out = 1000)
		lkl					<- mahaltest.sulkl(rho_lkl, n.of.x, p, norm = lkl_norm, support = lkl_support)
		rho_lkl             <- c(min(rho_lkl), rho_lkl, max(rho_lkl))
		lkl                 <- c(0, lkl, 0)
		lkl_norm            <- c(0, lkl_norm, 0)
		df_lkl 				<- data.frame(x = rho_lkl, yes = lkl, no = lkl * lkl_norm)
		df_lkl$distribution <- "summary likelihood"
		rho_pow	 			<- seq(pow_support[1], pow_support[2], length.out = 1000)
		pow					<- vartest.pow(rho_pow, scale, df, c.l, c.u, trafo = 1, norm = pow_norm)		
#		rho_pow             <- c(min(rho_pow), rho_pow, max(rho_pow))
#		pow                 <- c(0, pow, 0)
#		pow_norm            <- c(0, pow_norm, 0)
		df_pow 				<- data.frame(x = rho_pow, yes = pow, no = pow * pow_norm)
		df_pow$distribution <- "ABC approximation"
		gdf 				<- rbind(df_pow, df_lkl)
		gdf					<- melt(gdf, id.vars = c("x", "distribution"))
		p 					<- ggplot(data = gdf, aes(x = x, y = value, colour = distribution, linetype = variable)) +
								geom_polygon(data = subset(gdf, distribution == 'summary likelihood'), fill = 'grey70') +
								geom_vline(xintercept = rho.star, colour = 'black', linetype = "dotted") +
								geom_vline(xintercept = c(tau.l, tau.u), linetype = "dotted") + 
								geom_hline(yintercept = mx.pw, linetype = "dotted") + geom_line() +
								scale_y_continuous(lim = c(-0.02, 2.3), expand = c(0, 0)) +
								#scale_colour_brewer(palette='Set1', guide=FALSE) + 
								scale_linetype_manual(values = c('solid', 'longdash'), guide = FALSE) + 
								scale_colour_manual(values = c('black', 'transparent'), guide = FALSE) + 
								labs(x = expression(rho), y = "", linetype = "Normalized", colour='Distribution') +
								theme_bw() + theme(legend.position = 'bottom') #+ guides(colour=guide_legend(ncol=2))
		#p 					<- p + ggtitle(paste("n.of.y=", df+1, "\ntau.l=", tau.l,"\ntau.u=", tau.u,"\nKL=", KL_div))
		print(p)
	}
	pw.cmx 	<- ifelse(calibrate.tau.u, pw.cmx, vartest.pow(rho = rho.star, scale, df, c.l, c.u))	
	c(KL_div = KL_div, tau.l = tau.l, tau.u = tau.u, c.l = c.l, c.u = c.u, pw.cmx = pw.cmx)	
}

#calibrate according to KL divergence using an F-test statistic	
mahaltest.getkl.F <- function(n.of.x, n.of.y, p, tau.u, mx.pw = 0.9, alpha = 0.01, pow_scale = 1.5, plot = F) 
{	
	#	calibrate tau	
	tmp				<- ftest.calibrate.tau( mx.pw, n.of.y, p, tau.u, alpha )						#tau.u is taken as upper bound on calibrated tau.u
	crit			<- tmp['c']
	tau				<- tmp['tau']
	curr.pw			<- tmp['curr.pw']	
	if (abs(curr.pw - mx.pw) > 0.09) 	
		stop("tau.up not accurate")			
	#	truncate pow and compute the normalizing constant pow_norm	
	pow_support 	<- c(0, tau * pow_scale) 	
	pow_norm 		<- ftest.pow.norm(tau, n.of.y, p, support = pow_support)
	#	truncate lkl and compute the normalizing constant lkl_norm 
	lkl_support		<- pow_support
	lkl_norm		<- mahaltest.sulkl.norm(n.of.x, p, support = lkl_support)
	integral_range	<- pow_support			
	lkl_arg			<- list(n.of.x = n.of.x, p = p, norm = lkl_norm, support = lkl_support)
	pow_arg			<- list(tau = tau, n.of.y = n.of.y, p = p, alpha = alpha, norm = pow_norm, support=pow_support)	
	tmp 			<- integrate(kl.integrand, lower = integral_range[1], upper = integral_range[2], dP = mahaltest.sulkl, dQ = ftest.pow, P_arg = lkl_arg, Q_arg = pow_arg)
	KL_div			<- tmp$value
	if (tmp$message != "OK") 
	{
		warning(tmp$message)
	}
	if (plot) 
	{
		rho_lkl 			<- seq(lkl_support[1], lkl_support[2], length.out = 1000)
		lkl					<- mahaltest.sulkl(rho_lkl, n.of.x, p, norm = lkl_norm, support = lkl_support)
		rho_lkl             <- c(min(rho_lkl), rho_lkl, max(rho_lkl))
		lkl                 <- c(0, lkl, 0)
		lkl_norm            <- c(0, lkl_norm, 0)
		df_lkl 				<- data.frame(x = rho_lkl, yes = lkl, no = lkl * lkl_norm)
		df_lkl$distribution <- "summary likelihood"
		rho_pow	 			<- seq(pow_support[1], pow_support[2], length.out = 1000)
		pow					<- ftest.pow(rho_pow, tau, n.of.y, p, alpha = alpha, support = pow_support, norm = pow_norm)
		df_pow 				<- data.frame(x = rho_pow, yes = pow, no = pow * pow_norm)
		df_pow$distribution <- "ABC approximation"
		gdf 				<- rbind(df_pow, df_lkl)
		gdf					<- melt(gdf, id.vars = c("x", "distribution"))
		p 					<- ggplot(data = gdf, aes(x = x, y = value, colour = distribution, linetype = variable)) +
								geom_polygon(data = subset(gdf, distribution == 'summary likelihood'), fill = 'grey70') +
								geom_vline(xintercept = tau, linetype = "dotted") + 
								geom_hline(yintercept = mx.pw, linetype = "dotted") + geom_line() +
								scale_y_continuous(lim = c(-0.02, 2.3), expand = c(0, 0)) +
								#scale_colour_brewer(palette='Set1', guide=FALSE) + 
								scale_linetype_manual(values = c('solid', 'longdash'), guide = FALSE) + 
								scale_colour_manual(values = c('black', 'transparent'), guide = FALSE) + 
								labs(x = expression(rho), y = "", linetype = "Normalized", colour='Distribution') +
								theme_bw() + theme(legend.position = 'bottom') #+ guides(colour=guide_legend(ncol=2))
		#p 					<- p + ggtitle(paste("n.of.y=", df+1, "\ntau.l=", tau.l,"\ntau.u=", tau.u,"\nKL=", KL_div))
		print(p)
	}
	tmp			<- c(KL_div, tau, crit, curr.pw,abs(curr.pw - mx.pw))
	names(tmp)	<- c('KL_div', 'tau', 'c', 'pw.cmx', 'err.pw')
	tmp
}
#------------------------------------------------------------------------------------------------------------------------
# @title Density of the summary likelihood of the Mahalanobis distance for multivariate normal summary values
# @description		The density of the scaled Mahalanobis distance is chi-squared on p degrees-of-freedom.
# @param rho 		Auxiliary error parameter
# @param n.of.x     Number of observations (if missing then use chi-squared summary likelihood, else use an F-distributed summary likelihood
# @param p  		Dimensionality of the underlying MVN distribution
# @param norm 		scalar, 0<\code{norm}<=1, normalization constant for the truncated summary likelihood.
# @param support 	vector of dimension 2, support of the truncated summary likelihood.
# @param log 		logical; if \code{TRUE}, densities d are given as log(d). 
# @note The summary likelihood can be truncated to \code{support} and then standardized with \code{norm}.
# For computational efficiency, both \code{norm} and \code{support} must be provided although each one can be derived from the other. 
# @seealso \code{\link{mahaltest.calibrate}}, \code{\link{mahaltest.getkl}} for an example.
mahaltest.sulkl <- function(rho, n.of.x, p, norm = 1, support= c(0, Inf), log = FALSE) 
{
	ans 				<- rho
	in_support 			<- (rho >= support[1] & rho <= support[2])
	ans[!in_support]	<- 0
    if (any(in_support)) 
	{			
	    if(is.na(n.of.x)) ans[in_support]	<- dchisq(rho[in_support], df = p) / norm
		else ans[in_support]	<- df(rho[in_support], p, n.of.x - p) / norm
    }	    
		    
	if(log)
		ans				<- log(ans)
	return(ans)
}
#------------------------------------------------------------------------------------------------------------------------
# @title Area under the summary likelihood of the Mahalanobis distance for multivariate normal summary values
# @description This function computes the area under the summary likelihood \code{mahaltest.sulkl}.
# @inheritParams mahaltest.sulkl
# @seealso \code{\link{mahaltest.calibrate}}
mahaltest.sulkl.norm	<- function(n.of.x, p, support = c(0, Inf))
{
    ans	<- integrate(mahaltest.sulkl, lower = support[1], upper = support[2], n.of.x = n.of.x, p = p, norm = 1, support = support, log = FALSE)	
	ans$value
}
