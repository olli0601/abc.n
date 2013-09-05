

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
		log_Q_x[tmp]<- log(.Machine$double.xmin)
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
#' @param KL_args additional arguments to be passed to \code{KL_divergence}.
#' @param verbose logical, if \code{TRUE}, print warnings.
#' @return the minimized Kullback-Leibler divergence (scalar).
#' @export
nabc.kl.optimize<- function(x_value, x_name, is_integer=FALSE, KL_divergence, KL_args, verbose=FALSE) 
{
	if (verbose)	print(x_value)
	if(is_integer)
		x_value 	<- round(x_value)	
	

	KL_args 			<- c(list(x_name = x_value),KL_args)
	names(KL_args)[1] 	<- x_name
	tmp				<- do.call(KL_divergence,KL_args)	
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
nabc.calibrate.tau.nomxpw.yesKL <- function(test_name, KL_args, tau.u.lb=1, max.it=100,debug = 0, plot = F, plot_debug = FALSE) {
	
	stopifnot(tau.u.lb>0,max.it>0)

	if(test_name=="mutost"){
		
		#check KL_args
		stopifnot(all(c("n.of.x","s.of.x","n.of.y","s.of.y","mx.pw","alpha","pow_scale")%in%names(KL_args)))
		with(KL_args,stopifnot(n.of.x>0,s.of.x>0, n.of.y>0, s.of.y>0,mx.pw>0, alpha>0, pow_scale>0))
		
	}else{
		stop(test_name,"is currently not supported.")
	}	
	
	if(!debug){
		#all in C
		KL_args$tau.u <- as.double(tau.u.lb)	
		suppressWarnings({ #suppress numerical inaccuracy warnings
			tmp <- .Call("abcCalibrate_minimize_KL", test_name, "tau_no_max_pow", KL_args, as.integer(max.it))
		})
		
		names(tmp)<-c("KL_div","tau.u","pw.cmx")
		
	}else{
	
	#all in R
	KL_divergence=switch(test_name,mutost="nabc.mutost.kl")
	
	KL_args$calibrate.tau.u <- F
	KL_args$tau.u <- tau.u.lb
	KL_args$plot <- plot_debug
	if (plot_debug) {
		cairo_pdf("KL_initial.pdf", onefile = T)
	}
	
	
	KL.of.tau.u <- do.call(KL_divergence, KL_args)["KL_div"]
	
	
	if (plot_debug) {
		dev.off()
	}
	
	#to find tau.u.ub increase tau.u until the KL increases too
	tau.u.ub <- 2*tau.u.lb
	curr.it <- max.it
	
	KL_args$tau.u <- tau.u.ub
	KL_args$plot <- F
	KL.of.tau.u.ub <- do.call(KL_divergence, KL_args)["KL_div"]
		
	while (KL.of.tau.u.ub < KL.of.tau.u && curr.it > 0) {
		curr.it <- curr.it - 1
		KL.of.tau.u <- KL.of.tau.u.ub
		tau.u.ub <- 2 * tau.u.ub
		KL_args$tau.u <- tau.u.ub
		KL.of.tau.u.ub <- do.call(KL_divergence, KL_args)["KL_div"]
	}
	
	
		
	if (curr.it == 0) 
		stop("could not find upper bound for tau.u")
	
	if(curr.it < max.it){
		tau.u.lb <- tau.u.ub/4
	}
	
	#minimize KL_tau.u between [tau.u.lb - tau.u.ub]	
	if (plot_debug) {
		cairo_pdf("KL_optimization.pdf", onefile = T)
	}
	KL_args$plot <- plot_debug
	KL_args["tau.u"]<-NULL
	tmp <- optimize(nabc.kl.optimize, interval = c(tau.u.lb,tau.u.ub), x_name="tau.u" ,KL_divergence= KL_divergence, KL_args = KL_args, verbose = debug)
	if (plot_debug) {
		dev.off()
	}
	tau.u <- tmp$minimum
	KL_args$tau.u <- tau.u
	KL_args$plot <- plot	
	tmp <- do.call(KL_divergence, KL_args)
	
	}
	
	return(c(n.of.y = n.of.y, tmp))
	
}

#------------------------------------------------------------------------------------------------------------------------
#' Calibrate the number of simulated summary values and the equivalence region for a specified test of equivalence by minimising the Kullback-Leibler divergence from the standardized power function to the summary likelihood.
#' @param KL_divergence character, name of the function to compute the KL divergence for the test of equivalence. This function must have at least two arguments: \code{n.of.y}, the number of simulated summaries, and \code{plot}, which flag some plotting options. 
#' @param KL_args list of argument to be passed to \code{KL_divergence}. 
#' @param max.it this algorithm stops prematurely when the number of iterations to calibrate the number of simulated data points exceeds 'max.it'
#' @param debug flag if C implementation is used.
#' @param plot if \code{TRUE}, plot the result of minimization. Only if implemented in \code{KL_divergence}.
#' @return	vector of variable length. The first value is \code{n.of.y}: the adjusted number of simulated summaries. The rest is returned from \code{KL_divergence}.
#' @export
#' @import stats
#' @example inst/examples/ex.nabc.adjust.n.of.y.KL.R
#' 
nabc.calibrate.m.and.tau.yesmxpw.yesKL <- function(test_name, KL_args, max.it = 100, debug = 0, plot = F) {
	
	stopifnot(max.it > 0)

	if (test_name == "mutost") {

		#check KL_args
		stopifnot(all(c("n.of.x", "s.of.x", "n.of.y", "s.of.y", "mx.pw", "alpha", "tau.u", "pow_scale") %in% names(KL_args)))
		with(KL_args, stopifnot(n.of.x > 0, s.of.x > 0, n.of.y > 0, s.of.y > 0, mx.pw > 0, alpha > 0, tau.u>0, pow_scale > 0))

	} else {
		stop(test_name, "is currently not supported.")
	}

	if (!debug) {
		#all in C
		suppressWarnings({ #suppress numerical inaccuracy warnings
			tmp <- .Call("abcCalibrate_tau_nomxpw_yesKL", test_name, KL_args, as.integer(max.it))
		})

		names(tmp) <- c("KL_div", "tau.u", "pw.cmx")

	} else {

		#all in R
		KL_divergence = switch(test_name, mutost = "nabc.mutost.kl")

		KL_args$plot <- F
		KL_args$calibrate.tau.u <- T

		if (debug) {
			cairo_pdf("KL_initial.pdf", onefile = T)
		}
		KL.of.yn <- do.call(KL_divergence, KL_args)["KL_div"]
		if (debug) {
			dev.off()
		}
		n.of.y <- KL_args$n.of.y
		KL_args$n.of.y <- n.of.y - 1
		KL_args$plot <- F
		KL.of.yn_m1 <- do.call(KL_divergence, KL_args)["KL_div"]

		decrease_n.of.y <- as.logical(KL.of.yn_m1 < KL.of.yn)

		if (decrease_n.of.y) {
			#optimize between 1 and n.of.y
			yn.lb <- 1
			yn.ub <- n.of.y
		} else {
			#find upper bound for optimize
			yn.lb <- n.of.y
			curr.it <- max.it
			yn.ub <- 2 * n.of.y

			KL_args$n.of.y <- yn.ub
			KL.of.yn_ub <- do.call(KL_divergence, KL_args)["KL_div"]

			while (KL.of.yn_ub < KL.of.yn && curr.it > 0) {
				curr.it <- curr.it - 1
				KL.of.yn <- KL.of.yn_ub
				yn.ub <- 2 * yn.ub
				KL_args$n.of.y <- yn.ub
				KL.of.yn_ub <- do.call(KL_divergence, KL_args)["KL_div"]
			}
			if (curr.it == 0) 
				stop("nabc.calibrate.m.and.tau.yesmxpw.yesKL: could not find upper bound for yn")
		}
		if (debug) {
			cairo_pdf("KL_optimization.pdf", onefile = T)
		}
		KL_args$plot <- debug
		KL_args["n.of.y"] <- NULL
		tmp <- optimize(nabc.kl.optimize, interval = c(yn.lb, yn.ub), x_name = "n.of.y", is_integer = T, KL_divergence = KL_divergence, 
			KL_args = KL_args, verbose = debug, tol = 1)
		if (debug) {
			dev.off()
		}
		yn <- round(tmp$minimum)

		KL_args$n.of.y <- yn
		KL_args$plot <- plot
		tmp <- do.call(KL_divergence, KL_args)
	}

	return(c(n.of.y = yn, tmp))
}