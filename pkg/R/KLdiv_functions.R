

#' Integrand of the Kullback-Leibler divergence D(P||Q) for generic distributions
#' @param x value at which integrand is evaluated. Can be a vector.
#' @param dP,dQ name of the functions that compute the density of P and Q
#' @param P_arg,Q_arg list of arguments for \code{dP} and \code{dQ}
#' @export
#'
nabc.kl.integrand<-function(x,dP,dQ,P_arg,Q_arg)
{
	#this function must accept x as a vector	
	log_P_x	<- do.call(dP,c(list(x,log=T),P_arg))				
	log_Q_x	<- do.call(dQ,c(list(x,log=T),Q_arg))
		
	if(any(tmp<-(log_Q_x == -Inf)))
	{
		warning("Q(x) is 0 at some point! Usually this happens due to numerical inacurracy in the tail of Q. For simplicity we assume that log(Q(x))=-100 at these points.")
		log_Q_x[tmp]<- -100
	}
	ans		<- (log_P_x-log_Q_x)*exp(log_P_x)
	if(any(!is.finite(ans)))
	{
		print(data.frame(x,P_x=exp(log_P_x),log_P_x,log_Q_x,ans))
		stop()
	}			
	return(ans)
}

#' A wrapper to minimize \code{KL_divergence} over the parameter \code{x_name} using the function \link{optimize}
#' @param x_value value tested.
#' @param x_name name of the parameter over which minimization is performed.
#' @param is_integer if \code{TRUE}, \code{x_value} is rounded to the closest integer. See \link{round}.
#' @param KL_divergence character, name of the function that computes the KL divergence.
#' @param args additional arguments to be passed to \code{KL_divergence}.
#' @param verbose logical, if \code{TRUE}, print warnings.
#' @return the minimized Kullback-Leibler divergence (scalar).
#' @export
nabc.kl.optimize<- function(x_value, x_name, is_integer=FALSE, KL_divergence, args, verbose=FALSE) 
{
	if (verbose)	print(x_value)
	if(is_integer)
		x_value 	<- round(x_value)	
	

	args 			<- c(list(x_name = x_value),args)
	names(args)[1] 	<- x_name
	tmp				<- do.call(KL_divergence,args)	
	if(0 || !verbose)								#TODO anton check
		options(warn=0)
	
	return(tmp["KL_div"])
}

#------------------------------------------------------------------------------------------------------------------------
#' Calibrate the equivalence region for a specified test of equivalence by minimising the Kullback-Leibler divergence between the power function and the summary likelihood.
#' @param KL_divergence character, name of the function to compute the KL divergence for the test of equivalence. This function must have at least two arguments: \code{tau.u}, the upper tolerance of the equivalence region, and \code{plot}, which flag some plotting options. 
#' @inheritParams nabc.calibrate.m.and.tau.yesmxpw.yesKL
#' @param tau.u.lb	guess on an lower bound on the upper tolerance of the equivalence region
#' @return	vector of variable length. The first value is \code{n.of.y}: the number of simulated summaries. The rest is returned from \code{KL_divergence}.
#' @note This procedure is not relevant for the test of dispersion equivalence for normally distributed variables (chisqstretch).
#' @export
#' @import stats
#' @example inst/examples/ex.nabc.adjust.tau.lowup.KL.R
#'
nabc.calibrate.tau.nomxpw.yesKL <- function(KL_divergence, args, tau.u.lb=1, max.it=100,debug = 0, plot = F) {
	
	if(!tau.u.lb){stop("tau.u.lb must be >0")}
	
	args$calibrate.tau.u <- F
	args$tau.u <- tau.u.lb
	args$plot <- debug
	if (debug) {
		cairo_pdf("KL_initial.pdf", onefile = T)
	}
	KL.of.tau.u <- do.call(KL_divergence, args)["KL_div"]
	if (debug) {
		dev.off()
	}
	
	#to find tau.u.ub increase tau.u until the KL increases too
	tau.u.ub < 2*tau.u.lb
	curr.it <- max.it
	
	args$tau.u <- tau.u.ub
	args$plot <- F
	KL.of.tau.u.ub <- do.call(KL_divergence, args)["KL_div"]
	
	while (KL.of.tau.u.ub < KL.of.tau.u && curr.it > 0) {
		curr.it <- curr.it - 1
		KL.of.tau.u <- KL.of.tau.u.ub
		tau.u.ub <- 2 * tau.u.ub
		args$tau.u <- tau.u.ub
		KL.of.tau.u.ub <- do.call(KL_divergence, args)["KL_div"]
	}
	
	if (curr.it == 0) 
		stop("nabc.adjust.tau.lowup.KL: could not find upper bound for tau.u")
	
	#minimize KL_tau.u between [tau.u.lb - tau.u.ub]	
	if (debug) {
		cairo_pdf("KL_optimization.pdf", onefile = T)
	}
	args$plot <- debug
	args <- args[-which(names(args) == "tau.u")]
	tmp <- optimize(nabc.kl.optimize, interval = c(tau.u.lb,tau.u.ub), x_name="tau.u" ,KL_divergence= KL_divergence, args = args, verbose = debug)
	if (debug) {
		dev.off()
	}
	tau.u <- tmp$minimum
	args$tau.u <- tau.u
	args$plot <- plot	
	tmp <- do.call(KL_divergence, args)
	
	return(c(n.of.y = n.of.y, tmp))
}

#------------------------------------------------------------------------------------------------------------------------
#' Calibrate the number of simulated summary values and the equivalence region for a specified test of equivalence by minimising the Kullback-Leibler divergence from the standardized power function to the summary likelihood.
#' @param KL_divergence character, name of the function to compute the KL divergence for the test of equivalence. This function must have at least two arguments: \code{n.of.y}, the number of simulated summaries, and \code{plot}, which flag some plotting options. 
#' @param args list of argument to be passed to \code{KL_divergence}. 
#' @param max.it this algorithm stops prematurely when the number of iterations to calibrate the number of simulated data points exceeds 'max.it'
#' @param debug flag if C implementation is used.
#' @param plot if \code{TRUE}, plot the result of minimization. Only if implemented in \code{KL_divergence}.
#' @return	vector of variable length. The first value is \code{n.of.y}: the adjusted number of simulated summaries. The rest is returned from \code{KL_divergence}.
#' @export
#' @import stats
#' @example inst/examples/ex.nabc.adjust.n.of.y.KL.R
#' 
nabc.calibrate.m.and.tau.yesmxpw.yesKL <- function(KL_divergence, args, max.it = 100, debug = 0, plot = F) 
{	
	args$plot 			<- F
	args$calibrate.tau.u<- T
	
	if (debug) 
	{
		cairo_pdf("KL_initial.pdf", onefile = T)
	}	
	KL.of.yn 			<- do.call(KL_divergence, args)["KL_div"]
	if (debug) 
	{
		dev.off()
	}
	n.of.y				<- args$n.of.y 
	args$n.of.y 		<- n.of.y - 1
	args$plot 			<- F
	KL.of.yn_m1 		<- do.call(KL_divergence, args)["KL_div"]
	
	decrease_n.of.y 	<- as.logical(KL.of.yn_m1 < KL.of.yn)
	
	if (decrease_n.of.y) 
	{
		#optimize between 1 and n.of.y
		yn.lb 			<- 1
		yn.ub 			<- n.of.y
	} 
	else 
	{
		#find upper bound for optimize
		yn.lb 			<- n.of.y
		curr.it 		<- max.it
		yn.ub 			<- 2 * n.of.y
		
		args$n.of.y 	<- yn.ub
		KL.of.yn_ub 	<- do.call(KL_divergence, args)["KL_div"]
		
		while (KL.of.yn_ub < KL.of.yn && curr.it > 0) 
		{
			curr.it 	<- curr.it - 1
			KL.of.yn 	<- KL.of.yn_ub
			yn.ub 		<- 2 * yn.ub
			args$n.of.y <- yn.ub
			KL.of.yn_ub <- do.call(KL_divergence, args)["KL_div"]
		}		
		if (curr.it == 0) 
			stop("nabc.calibrate.m.and.tau.yesmxpw.yesKL: could not find upper bound for yn")
	}
	if (debug) 
	{
		cairo_pdf("KL_optimization.pdf", onefile = T)
	}
	args$plot 	<- debug
	args 		<- args[-which(names(args) == "n.of.y")]
	tmp 		<- optimize(nabc.kl.optimize, interval = c(yn.lb, yn.ub), x_name = "n.of.y", is_integer = T, KL_divergence = KL_divergence, args = args, verbose = debug, tol = 1)	
	if (debug) 
	{
		dev.off()
	}
	yn 			<- round(tmp$minimum)
	
	args$n.of.y <- yn
	args$plot 	<- plot
	tmp 		<- do.call(KL_divergence, args)
	
	return(c(n.of.y = yn, tmp))
}