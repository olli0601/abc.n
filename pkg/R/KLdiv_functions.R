
#' Integrand of the Kullback-Leibler divergence for generic distribution
#' @param x value at which integrand is evaluated. Can be a vector.
#' @param P_dist name of the function that compute the density of P
#' @param Q_dist name of the function that compute the density of Q
#' @param P_arg list of arguments for \code{P_dist}
#' @param Q_arg list of arguments for \code{Q_dist}
integrand_KL_divergence_1D<-function(x,P_dist,Q_dist,P_arg,Q_arg){
	
	P_x<-do.call(P_dist,c(list(x=x),P_arg))					
	log_P_x<-do.call(P_dist,c(list(x=x,log=T),P_arg))				
	log_Q_x<-do.call(Q_dist,c(list(x=x,log=T),Q_arg))
	
	return((log_P_x-log_Q_x)*P_x)
}

#' Integrand of the Kullback-Leibler divergence for mutost
#' @param rho value at which integrand is evaluated. Can be a vector.
#' @inheritParams KL_divergence_1D_mutost
integrand_KL_divergence_1D_mutost <- function(rho, n.of.x, s.of.x, n.of.y, s.of.y, tau.u, alpha, lkl_norm = 1, pow_norm = 1, lkl_support, pow_support) {

	ssn<-s.of.x/sqrt(n.of.x)
	P_rho <- dt(rho/ssn, df = n.of.x - 1)/ssn/lkl_norm
	P_rho[rho<lkl_support[1] | rho>lkl_support[2]]<-0
	
	suppressWarnings({ #suppress numerical inaccuracy warnings
		Q_rho <- .Call("abcMuTOST_pow", rho, n.of.y - 1, tau.u, s.of.y/sqrt(n.of.y), alpha)/pow_norm
	})
	Q_rho[rho<pow_support[1] | rho>pow_support[2]]<-0

	ans <- log(P_rho/Q_rho) * P_rho

	if (any(!is.finite(ans))) {
		#print(data.frame(rho=rho,P=P_rho,Q=Q_rho,ans=ans))
		#warning("infinite KL integrand replaced by 0 (usually this happens when Q_rho=0, i.e. in the tail)")
		ans[!is.finite(ans)] <- 0
	}

	if (any(is.na(ans))) {
		#warning("NA KL integrand replaced by 0")
		ans[is.na(ans)] <- 0
	}

	return(ans)
}

#' Compute the Kullback-Leibler divergence D(P||Q) from the 1D-distribution P and Q, when the density is available.
#' @inheritParams integrand_KL_divergence_1D
#' @param lower lower bound for integration. Can be infinite.
#' @param upper upper bound for integration. Can be infinite.
#' @export 
#' @import stats
#' @examples
#' KL_divergence_1D(P_dist="dexp",Q_dist="dexp",P_arg=list(rate=0.2),Q_arg=list(rate=0.4),lower=0,upper=Inf)
#' #theoretical divergence = log(0.2/0.4)+(0.4-0.2)-1 = 1-log(2) = 0.307
KL_divergence_1D<-function(P_dist,Q_dist,P_arg,Q_arg,lower,upper){

	KL_d<-integrate(integrand_KL_divergence_1D,lower,upper,P_dist,Q_dist,P_arg,Q_arg)
	
	if(KL_d$message!="OK") warning(KL_d$message)

	return(KL_d$value)
}

#' Compute the Kullback-Leibler divergence D(P||Q) where P is the summary likelihood and Q is the mutost power.
#' @inheritParams KL_divergence_mutost
#' @param tau.u	upper tolerance of the equivalence region
#' @param lkl_norm constant for normalization for the truncated summary likelihood (default to 1)
#' @param pow_norm constant for normalization for the truncated power (default to 1)
#' @param lower lower bound for integration. Can be infinite.
#' @param upper upper bound for integration. Can be infinite.
#' @export 
#' @import stats
#'
KL_divergence_1D_mutost <- function(n.of.x, s.of.x, n.of.y, s.of.y, tau.u, alpha, lkl_norm = 1, pow_norm = 1, lkl_support, pow_support, lower, upper) {

	KL_d <- integrate(integrand_KL_divergence_1D_mutost, lower, upper, n.of.x, s.of.x, n.of.y, s.of.y, tau.u, alpha, lkl_norm, pow_norm, lkl_support, pow_support)

	if (KL_d$message != "OK") {
		warning(KL_d$message)
	}

	return(KL_d$value)
}


#' Compute Kullback-Leibler divergence between the summary likelihood and the power function of mutost after calibration of tau.u 
#' @param n.of.x number of observed summary values
#' @param s.of.x standard deviation of observed summary values
#' @param n.of.y number of simulated summary values
#' @param s.of.y standard deviation of simulated summary values
#' @param mx.pw maximum power at the point of reference (rho.star=0) (only when \code{calibrate.tau.u==TRUE}).
#' @param calibrate.tau.u if \code{TRUE} the upper tolerance of the equivalence region (\code{tau.u}) is calibrated so that power at the point of reference is equal to \code{mx.pw}
#' @param tau.u	upper tolerance of the equivalence region. If \code{calibrate.tau.u==TRUE}, it's just a guess on an upper bound on the upper tolerance of the equivalence region.
#' @param alpha level of the equivalence test
#' @param debug flag if C implementation is used
#' @param plot whether to plot the two distributions
#' @return	vector of length 3
#' 	\item{KL_div}{the Kullback Leibler divergence}		
#' 	\item{tau.u}{upper tolerance of the equivalence region}
#' 	\item{pw.cmx}{actual maximum power associated with the equivalence region}
#' @note If \code{calibrate.tau.u==TRUE}
#' @export
#' @examples
#' KL_divergence_mutost(n.of.x=60,s.of.x=0.1,n.of.y=60,s.of.y=0.3, mx.pw=0.9,
#' tau.u.ub=1, alpha=0.01,plot=T)
#'
KL_divergence_mutost <- function(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, calibrate.tau.u=F,tau.u, alpha, lkl_norm, lkl_support, debug = 0, plot = F) {

	if(calibrate.tau.u){
		#calibrate tau.u constrained on yn, alpha and mx.pw	
		tmp <- nabc.mutost.onesample.tau.lowup.pw(mx.pw, n.of.y - 1, s.of.y/sqrt(n.of.y), tau.u, alpha, debug)
		tau.u <- tmp[2]
		pw.cmx <- tmp[3]
		if(abs(pw.cmx-mx.pw)>0.09)	stop("tau.up not accurate")
	}

	#truncate pow to [-2*tau.u,2*tau.u] and compute pow_norm
	pow_support<- c(-2*tau.u,2*tau.u)
	rho <- seq(pow_support[1], pow_support[2], length.out = 1000)
	suppressWarnings({ #suppress numerical inaccuracy warnings
		pow <- .Call("abcMuTOST_pow", rho, n.of.y - 1, tau.u, s.of.y/sqrt(n.of.y), alpha)
	})
	pow_norm <- sum(pow) * diff(rho)[1]
	pow<-pow/pow_norm
	
	integral_range<-range(c(lkl_support,pow_support))	
	
	KL_div <- KL_divergence_1D_mutost(n.of.x, s.of.x, n.of.y, s.of.y, tau.u, alpha, lkl_norm = lkl_norm, pow_norm = pow_norm, lkl_support= lkl_support, pow_support= pow_support, lower = integral_range[1], 
		upper = integral_range[2])

	if (plot) {
		ssn<-s.of.x/sqrt(n.of.x)
		rho_lkl<- seq(lkl_support[1], lkl_support[2], length.out = 1000)
		lkl <- dt(rho_lkl/ssn, n.of.x - 1)/ssn/lkl_norm
		plot(rho_lkl, lkl, t = "l")
		lines(rho, pow, col = "red")
		title(main = paste("n.of.y=", n.of.y,"\ntau.u=",tau.u ,"\nKL=", KL_div))
	}

	pw.cmx <- ifelse(calibrate.tau.u,pw.cmx,.Call("abcMuTOST_pow", rho=0, n.of.y - 1, tau.u, s.of.y/sqrt(n.of.y), alpha))
	
	ans <- c(KL_div = KL_div, tau.u = tau.u, pw.cmx= pw.cmx)
	return(ans)
}

#' A wrapper to minimize \code{KL_divergence_mutost} over n.of.y using the function \code{optimize}
#' @inheritParams KL_divergence_mutost
#' @export
KL_divergence_mutost_optimize_n.of.y <- function(n.of.y, n.of.x, s.of.x, s.of.y, mx.pw, calibrate.tau.u,tau.u , alpha, lkl_norm, lkl_support, debug = 0, plot = F) {
	
	if(plot){print(n.of.y)}
	
	n.of.y <- round(n.of.y)
	tmp <- KL_divergence_mutost(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, calibrate.tau.u,tau.u , alpha,lkl_norm, lkl_support, debug, plot)
	return(tmp["KL_div"])
}

#' A wrapper to minimize \code{KL_divergence_mutost} over tau.u using the function \code{optimize}
#' @inheritParams KL_divergence_mutost
#' @export
KL_divergence_mutost_optimize_tau.u <- function(tau.u , n.of.x, s.of.x, n.of.y,s.of.y, mx.pw, calibrate.tau.u=F, alpha, lkl_norm, lkl_support, debug = 0, plot = F) {
	
	if(plot){print(tau.u)}
	tmp <- KL_divergence_mutost(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, calibrate.tau.u=F,tau.u , alpha,lkl_norm, lkl_support, debug, plot)
	return(tmp["KL_div"])
}






