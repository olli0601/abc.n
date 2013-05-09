
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
integrand_KL_divergence_1D_mutost<-function(rho,n.of.x,s.of.x,n.of.y,s.of.y,tau.u,alpha,lkl_norm=1,pow_norm=1){
	
	
	P_rho<-dt(rho*sqrt(n.of.x)/s.of.x,df=n.of.x-1)/lkl_norm

	suppressWarnings({	#suppress numerical inaccuracy warnings
	Q_rho<-.Call("abcMuTOST_pow", rho, n.of.y-1, tau.u, s.of.y/sqrt(n.of.y), alpha)/pow_norm
	})
	
	ans<-log(P_rho/Q_rho)*P_rho
	
	if(any(!is.finite(ans))){
		#print(data.frame(rho=rho,P=P_rho,Q=Q_rho,ans=ans))
		warning("infinite KL integrand replaced by 0 (usually this happens when Q_rho=0, i.e. in the tail)")
		ans[!is.finite(ans)]<-0
	}
	
	if(any(is.na(ans))){
		warning("NA KL integrand replaced by 0")
		ans[is.na(ans)]<-0
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

	KL_d<-integrate(integrand_KL_divergence_1D_mutost,lower,upper,P_dist,Q_dist,P_arg,Q_arg)
	
	if(KL_d$message!="OK") warning(KL_d$message)

	return(KL_d$value)
}

#' Compute the Kullback-Leibler divergence D(P||Q) where P is the summary likelihood and Q is the mutost power.
#' @inheritParams KL_divergence_mutost_tau.u
#' @param tau.u	upper tolerance of the equivalence region
#' @param lkl_norm constant for normalization for the truncated summary likelihood (default to 1)
#' @param pow_norm constant for normalization for the truncated power (default to 1)
#' @param lower lower bound for integration. Can be infinite.
#' @param upper upper bound for integration. Can be infinite.
#' @export 
#' @import stats
#'
KL_divergence_1D_mutost<-function(n.of.x,s.of.x,n.of.y,s.of.y,tau.u,alpha,lkl_norm=1,pow_norm=1,lower,upper){

	KL_d<-integrate(integrand_KL_divergence_mutost,lower,upper,n.of.x,s.of.x,n.of.y,s.of.y,tau.u,alpha,lkl_norm,pow_norm)
	
	if(KL_d$message!="OK") warning(KL_d$message)

	return(KL_d$value)
}


#' Compute Kullback-Leibler divergence between the summary likelihood and the power function of mutost after calibration of tau.u 
#' @param n.of.x number of observed summary values
#' @param s.of.x standard deviation of observed summary values
#' @param n.of.y number of simulated summary values
#' @param s.of.y standard deviation of simulated summary values
#' @param mx.pw maximum power at the point of reference (rho.star=0)
#' @param tau.u.ub	guess on an upper bound on the upper tolerance of the equivalence region
#' @param alpha level of the equivalence test
#' @param debug flag if C implementation is used
#' @param plot whether to plot the two distributions
#' @return	vector of length 3
#' 	\item{KL_div}{the Kullback Leibler divergence}		
#' 	\item{tau.u}{upper tolerance of the equivalence region}
#' 	\item{pw.cmx}{actual maximum power associated with the equivalence region}
#' @export
#' @examples
#' KL_divergence_mutost_n(n.of.x=60,s.of.x=0.1,n.of.y=60,s.of.y=0.3, mx.pw=0.9, tau.u.ub=1, alpha=0.01,plot=T)
#'
KL_divergence_mutost_tau.u <- function(n.of.x,s.of.x,n.of.y,s.of.y, mx.pw, tau.u.ub, alpha, debug = 0, plot = F) {

	#calibrate tau.u constrained on yn, alpha and mx.pw	
	tmp <- nabc.mutost.onesample.tau.lowup.pw(mx.pw, n.of.y - 1, s.of.y/sqrt(n.of.y), tau.u.ub, alpha, debug)
	tau.u <- tmp[2]
	pw.cmx <- tmp[3]

	#truncate lkl and pow to -2*tau.u 2*tau.u
	rho		<- seq(-2*tau.u,2*tau.u,length.out=1e3)
	su.lkl<-dt(rho/s.of.x*sqrt(n.of.x),n.of.x-1)
	lkl_norm<- sum(su.lkl)*diff(rho)[1]
	suppressWarnings({	#suppress numerical inaccuracy warnings
		pow	<- .Call("abcMuTOST_pow", rho, n.of.y-1, tau.u, s.of.y/sqrt(n.of.y), alpha)
	})	
	pow_norm<-sum(pow)*diff(rho)[1]

	KL_div<-KL_divergence_1D_mutost(n.of.x,s.of.x,n.of.y,s.of.y,tau.u,alpha,lkl_norm= lkl_norm,pow_norm= pow_norm,lower=-2*tau.u,upper=2*tau.u)
	
	if(plot){
	plot(rho,su.lkl/lkl_norm,t='l')
	lines(rho,pow/pow_norm,col='red')
	title(main=paste("n.of.y=",n.of.y,"KL=",KL_div))
	}
	
	ans<-c(KL_div=KL_div,tau.u=tau.u,pw.cmx=pw.cmx)	
	return(ans)
}

#' A wrapper to minimize \code{KL_divergence_mutost_n} over n.of.y using the function \code{optimize}
#' @inheritParams KL_divergence_mutost_tau.u
#' @export
KL_optimize<-function(n.of.y,n.of.x,s.of.x,s.of.y, mx.pw, tau.u.ub, alpha, debug = 0,plot=F){
	
	n.of.y<-round(n.of.y)
	tmp<-KL_divergence_mutost_tau.u(n.of.x,s.of.x,n.of.y,s.of.y, mx.pw, tau.u.ub, alpha, debug,plot)
	return(tmp["KL_div"])
}








