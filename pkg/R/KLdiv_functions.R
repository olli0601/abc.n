
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
#' @inheritParams nabc.mutost.onesample.n.of.y
#' @param n.of.y number of simulated summary values
#' @param tau.u upper tolerance of the equivalence region
#' @param lkl_norm constant for normalization for the truncated summary likelihood (default to 1)
#' @param pow_norm constant for normalization for the truncated power (default to 1)
integrand_KL_divergence_mutost<-function(rho,n.of.x,n.of.y,tau.u,s.of.y,alpha,lkl_norm=1,pow_norm=1){
	
	
	P_x<-dt(rho,df=n.of.x-1)/lkl_norm
	suppressWarnings({	#suppress numerical inaccuracy warnings
	Q_x<-.Call("abcMuTOST_pow", rho, n.of.y-1, tau.u, s.of.y/sqrt(n.of.y), alpha)/pow_norm
	})
	
	ans<-log(P_x/Q_x)*P_x
	ans[is.na(ans)]<-0
	
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

	KL_d<-integrate(integrand_KL_divergence_mutost,lower,upper,P_dist,Q_dist,P_arg,Q_arg)
	
	if(KL_d$message!="OK") warning(KL_d$message)

	return(KL_d$value)
}

#' Compute the Kullback-Leibler divergence D(P||Q) where P is the summary likelihood and Q is the mutost power.
#' @inheritParams integrand_KL_divergence_mutost
#' @param lower lower bound for integration. Can be infinite.
#' @param upper upper bound for integration. Can be infinite.
#' @export 
#' @import stats
#'
KL_divergence_1D_mutost<-function(n.of.x,n.of.y,tau.u,s.of.y,alpha,lkl_norm=1,pow_norm=1,lower,upper){

	KL_d<-integrate(integrand_KL_divergence_mutost,lower,upper,n.of.x,n.of.y,tau.u,s.of.y,alpha,lkl_norm,pow_norm)
	
	if(KL_d$message!="OK") warning(KL_d$message)

	return(KL_d$value)
}


#' Compute Kullback-Leibler divergence between the summary likelihood and the power function of mutost after calibration of tau.u 
#' @inheritParams KL_divergence_1D_mutost
#' @param tau.u.ub	guess on an upper bound on the upper tolerance of the equivalence region
#' @param debug Flag if C implementation is used.
#' @export
#' @examples
#' KL_divergence_mutost_n(n.of.x=20,n.of.y=30, mx.pw=0.9, s.of.y=5, tau.u.ub=6, alpha=0.01, debug = 0)
#'
KL_divergence_mutost_n <- function(n.of.x,n.of.y, mx.pw, s.of.y, tau.u.ub, alpha, debug = 0) {

	#calibrate tau.u constrained on yn, alpha and mx.pw	
	tmp <- nabc.mutost.onesample.tau.lowup.pw(mx.pw, n.of.y - 1, s.of.y/sqrt(n.of.y), tau.u.ub, alpha, debug)
	tau.u <- tmp[2]
	pw.cmx <- tmp[3]

	#truncate lkl and pow to -2*tau.u 2*tau.u
	lkl_norm<-pt(2*tau.u,df= n.of.x-1)-pt(-2*tau.u,df= n.of.x-1)
	rho		<- seq(-2*tau.u,2*tau.u,length.out=1e3)
	suppressWarnings({	#suppress numerical inaccuracy warnings
		pw	<- .Call("abcMuTOST_pow", rho, n.of.y-1, tau.u, s.of.y/sqrt(n.of.y), alpha)
	})	
	pow_norm<-sum(diff(rho)*pw[-length(pw)])
	
	KL_div<-KL_divergence_1D_mutost(n.of.x,n.of.y,tau.u,s.of.y,alpha,lkl_norm= lkl_norm,pow_norm= pow_norm,lower=-2*tau.u,upper=2*tau.u)

	return(KL_div)
}








