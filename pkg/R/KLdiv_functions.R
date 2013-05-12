
#' Compute the density of the (possible truncated) power of the equivalence test for population dispersion of normal summary values
#' @inheritParams KL_divergence_mutost
#' @inheritParams nabc.chisqstretch.pow
#' @param rho vector of quantile
#' @param norm normalization constant for the truncated power function.
#' @param support vector of dimension 2. Support of the truncated power function.
#' @param log logical; if \code{TRUE}, densities d are given as log(d). 
#' @note The summary likelihood can be truncated to \code{support} and then standardized with \code{norm}.
#' For computational efficiency, both \code{norm} and \code{support} must be provided although each one can be derived (numerically) from the other.
#' @export
#' 	
dchisqstretch_pow <- function(rho, scale, df, cl, cu, norm=1, support=c(0,Inf), log=FALSE){
	
	ans<-rho
	in_support<- (rho >= support[1] & rho <= support[2])
	ans[!in_support] <- 0
	
	if(any(in_support)){			
		ans[in_support] <- (pchisq( cu/rho[in_support]*scale, df ) - pchisq( cl/rho[in_support]*scale, df))/norm
	}
	
	if(log){
		ans<-log(ans)
	}

	return(ans)
}

#' Compute the density of the (possibly truncated) summary likelihood for population dispersion of normal summary values
#' @inheritParams KL_divergence_mutost
#' @param rho vector of quantile
#' @param norm scalar, 0<\code{norm}<=1, normalization constant for the truncated summary likelihood.
#' @param support vector of dimension 2, support of the truncated summary likelihood.
#' @param log logical; if \code{TRUE}, densities d are given as log(d). 
#' @note The summary likelihood can be truncated to \code{support} and then standardized with \code{norm}.
#' For computational efficiency, both \code{norm} and \code{support} must be provided although each one can be derived from the other. 
#' \code{support=qigamma(c(1-norm,1+norm)/2,(n.of.x-2)/2,s.of.x^2*(n.of.x-1)/2)} and \code{norm=diff(pigamma(support,(n.of.x-2)/2,s.of.x^2*(n.of.x-1)/2)}.
#' @export	
#'
dchisqstretch_lkl <- function(rho, n.of.x, s.of.x, norm = 1, support= c(0,Inf), log=FALSE) {

	alpha<- (n.of.x-2)/2	 
	beta<- s.of.x^2*(n.of.x-1)/2

	if(norm>1){
		stop("norm must be <= 1")
	}

	ans <- rho
	in_support <- (rho >= support[1] & rho <= support[2])
	ans[!in_support] <- 0

	if (any(in_support)) {
		ans[in_support] <- densigamma(rho[in_support], alpha, beta)/norm
	}

	if(log){
		ans<-log(ans)
	}
	
	return(ans)
}

#' Compute the density of the (possible truncated) power of the equivalence test for population means of normal summary values
#' @inheritParams KL_divergence_mutost
#' @inheritParams nabc.mutost.pow
#' @param rho vector of quantile
#' @param norm normalization constant for the truncated power function.
#' @param support vector of dimension 2. Support of the truncated power function.
#' @param log logical; if \code{TRUE}, densities d are given as log(d). 
#' @note The summary likelihood can be truncated to \code{support} and then standardized with \code{norm}.
#' For computational efficiency, both \code{norm} and \code{support} must be provided although each one can be derived (numerically) from the other.
#' @export
#' 	
dMuTOST_pow <- function(rho, df, s.of.T, tau.u, alpha, norm=1, support=c(-Inf,Inf), log=FALSE){
	
	ans<-rho
	in_support<- (rho >= support[1] & rho <= support[2])
	ans[!in_support] <- 0
	
	if(any(in_support)){
		suppressWarnings({ #suppress numerical inaccuracy warnings
			ans[in_support] <- .Call("abcMuTOST_pow", rho[in_support], df, tau.u, s.of.T, alpha)/norm
		})
	}
	
	if(log){
		ans<-log(ans)
	}

	return(ans)
}

#' Compute the density of the (possibly truncated) summary likelihood for population means of normal summary values
#' @inheritParams KL_divergence_mutost
#' @param rho vector of quantile
#' @param norm scalar, 0<\code{norm}<=1, normalization constant for the truncated summary likelihood.
#' @param support vector of dimension 2, support of the truncated summary likelihood.
#' @param log logical; if \code{TRUE}, densities d are given as log(d). 
#' @note The summary likelihood can be truncated to \code{support} and then standardized with \code{norm}.
#' For computational efficiency, both \code{norm} and \code{support} must be provided although each one can be derived from the other.
#' \code{support=s.of.x/sqrt(n.of.x)*qt(c(1-norm,1+norm)/2,n.of.x-1)} and \code{norm=diff(pt(support/(s.of.x/sqrt(n.of.x)),n.of.x-1))}.
#' @export	
#'
dMuTOST_lkl <- function(rho, n.of.x, s.of.x, norm = 1, support= c(-Inf,Inf), log=FALSE) {

	ssn <- s.of.x/sqrt(n.of.x)
	df <- n.of.x - 1

	if(norm>1){
		stop("norm must be <= 1")
	}

	ans <- rho
	in_support <- (rho >= support[1] & rho <= support[2])
	ans[!in_support] <- ifelse(log,-Inf,0)

	if (any(in_support)) {
		if(log){
			ans[in_support] <- dt(rho[in_support]/ssn, df,log=T)-log(ssn*norm)
		}else{
			ans[in_support] <- dt(rho[in_support]/ssn, df)/ssn/norm
		}
	}

	return(ans)
}


#' Compute Kullback-Leibler divergence between the summary likelihood and the power function of chisqstretch 
#' @param n.of.x number of observed summary values
#' @param s.of.x standard deviation of observed summary values
#' @param n.of.y number of simulated summary values
#' @param s.of.y standard deviation of simulated summary values
#' @param mx.pw maximum power at the point of reference (rho.star=0) (only when \code{calibrate.tau.u==TRUE}).
#' @param alpha level of the equivalence test
#' @param calibrate.tau.u if \code{TRUE} the upper tolerance of the equivalence region (\code{tau.u}) is calibrated so that power at the point of reference is equal to \code{mx.pw}
#' @param tau.u	upper tolerance of the equivalence region. If \code{calibrate.tau.u==TRUE}, it's just a guess on an upper bound on the upper tolerance of the equivalence region.
#' @param for.mle	calibrate so that the mode of the power is at the MLE
#' @param pow_scale scale for the support of the standardized power. The power is truncated between \code{[tau.l/pow_scale,tau.u*pow_scale]} and then standardized.
#' @param debug flag if C implementation is used
#' @param plot whether to plot the two distributions
#' @return	vector of length 3
#' 	\item{KL_div}{the Kullback Leibler divergence}		
#' 	\item{tau.u}{upper tolerance of the equivalence region}
#' 	\item{pw.cmx}{actual maximum power associated with the equivalence region}
#' @note Whatever the value of \code{calibrate.tau.u}, the lower tolerance of the equivalence region (\code{tau.l}) is always numerically calibrated using \link{nabc.chisqstretch.tau.low}.
#' @export
#' @import ggplot2 reshape2
#' @examples
#' 
#' KL_divergence_mutost(n.of.x=60,s.of.x=0.1,n.of.y=60,s.of.y=0.3, mx.pw=0.9,
#' alpha=0.01, calibrate.tau.u=T, tau.u=1, plot=T)
#'
KL_divergence_chisqstretch <- function(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = F, tau.u, for.mle=0, pow_scale=2, debug = 0, plot = F) {

	if (calibrate.tau.u) {
		#calibrate tau.u constrained on yn, alpha and mx.pw	
		tmp		<- nabc.chisqstretch.tau.lowup(mx.pw, tau.u, n.of.y-1, alpha, for.mle=for.mle )
		tau.l	<- tmp[1]
		tau.u	<- tmp[2]
		pw.cmx	<- tmp[3]
		c.l		<- tmp[5]
		c.u		<- tmp[6]	
		if (abs(pw.cmx - mx.pw) > 0.09) 
			stop("tau.up not accurate")
			
	}else{
		#find tau.l, c.l and c.u 
	}

	#truncate pow and compute pow_norm	
	pow_support <- c(tau.l/pow_scale, tau.u*pow_scale) 
	rho <- seq(pow_support[1], pow_support[2], length.out = 1000)
	scale	<- ifelse(for.mle, n.of.y, n.of.y-1)		
	df <- n.of.y-1
	pow		<- nabc.chisqstretch.pow(rho, scale, df, c.l, c.u)
	pow_norm <- sum(pow) * diff(rho)[1]
	pow <- pow/pow_norm

	#compute the norm of lkl, given its support 
	alpha<- (n.of.x-2)/2	 
	beta<- s.of.x^2*(n.of.x-1)/2
	lkl_support<-pow_support
	lkl_norm<-diff(pigamma(support,alpha,beta))

	integral_range<-pow_support	
		
	lkl_arg<-list(n.of.x= n.of.x, s.of.x= s.of.x, norm = lkl_norm, support = lkl_support)
	pow_arg<-list(scale = scale, df = df, cl=cl, cu=cu, norm=pow_norm, support=pow_support)

	tmp <- integrate(integrand_KL_divergence_1D, lower = integral_range[1], upper = integral_range[2], dP=dchisqstretch_lkl,dQ=dchisqstretch_pow,P_arg=lkl_arg,Q_arg=pow_arg)

	if (tmp$message != "OK") {
		warning(tmp$message)
	}
	
	KL_div<-tmp$value
	
	if (plot) {
		rho_lkl <- seq(lkl_support[1], lkl_support[2], length.out = 1000)
		lkl<-dchisqstretch_lkl(rho_lkl, n.of.x, s.of.x, lkl_norm, lkl_support)
		df_lkl <- data.frame(x = rho_lkl, no = lkl*lkl_norm ,yes = lkl)
		df_lkl$distribution <- "summary likelihood"
		df_pow <- data.frame(x = rho, no = pow*pow_norm ,yes = pow)
		df_pow$distribution <- "ABC power"
		df <- rbind(df_pow, df_lkl)
		gdf<-melt(df,id.vars=c("x","distribution"))
		p <- ggplot(data = gdf, aes(x = x, y = value, colour = distribution,linetype=variable))
		p<-p+geom_vline(xintercept=c(-tau.u,tau.u),linetype="dotted")
		p<-p+geom_hline(yintercept= mx.pw,linetype="dotted")
		p<-p+ geom_line()
		p<-p+scale_linetype("truncated?")
		p<-p+xlab(expression(rho))+ylab("")
		p <- p + ggtitle(paste("n.of.y=", n.of.y, "\ntau.u=", tau.u, "\ntau.l=", tau.l,"\nKL=", KL_div))
		print(p)
	}

	pw.cmx <- ifelse(calibrate.tau.u, pw.cmx, nabc.chisqstretch.pow(rho=1, scale, df, c.l, c.u))
	
	ans <- c(KL_div = KL_div, tau.l = tau.l, tau.u = tau.u, c.l = c.l, c.u = c.u, pw.cmx = pw.cmx)
	return(ans)
}

#' Compute Kullback-Leibler divergence between the summary likelihood and the power function of mutost 
#' @param n.of.x number of observed summary values
#' @param s.of.x standard deviation of observed summary values
#' @param n.of.y number of simulated summary values
#' @param s.of.y standard deviation of simulated summary values
#' @param mx.pw maximum power at the point of reference (rho.star=0) (only when \code{calibrate.tau.u==TRUE}).
#' @param alpha level of the equivalence test
#' @param calibrate.tau.u if \code{TRUE} the upper tolerance of the equivalence region (\code{tau.u}) is calibrated so that power at the point of reference is equal to \code{mx.pw}
#' @param tau.u	upper tolerance of the equivalence region. If \code{calibrate.tau.u==TRUE}, it's just a guess on an upper bound on the upper tolerance of the equivalence region.
#' @param pow_scale scale for the support of the standardized power. The power is truncated between \code{pow_scale*[-tau.u,tau.u]} and then standardized.
#' @param debug flag if C implementation is used
#' @param plot whether to plot the two distributions
#' @return	vector of length 3
#' 	\item{KL_div}{the Kullback Leibler divergence}		
#' 	\item{tau.u}{upper tolerance of the equivalence region}
#' 	\item{pw.cmx}{actual maximum power associated with the equivalence region}
#' @export
#' @import ggplot2 reshape2
#' @examples
#' 
#' KL_divergence_mutost(n.of.x=60,s.of.x=0.1,n.of.y=60,s.of.y=0.3, mx.pw=0.9,
#' alpha=0.01, calibrate.tau.u=T, tau.u=1, plot=T)
#'
KL_divergence_mutost <- function(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = F, tau.u, pow_scale=1.5, debug = 0, plot = F) {

	if (calibrate.tau.u) {
		#calibrate tau.u constrained on yn, alpha and mx.pw	
		tmp <- nabc.mutost.onesample.tau.lowup.pw(mx.pw, n.of.y - 1, s.of.y/sqrt(n.of.y), tau.u, alpha)
		tau.u <- tmp[2]
		pw.cmx <- tmp[3]
		if (abs(pw.cmx - mx.pw) > 0.09) 
			stop("tau.up not accurate")
	}

	#truncate pow and compute pow_norm
	pow_support <- c(-tau.u, tau.u)* pow_scale
	rho <- seq(pow_support[1], pow_support[2], length.out = 1000)
	pow <- dMuTOST_pow(rho, n.of.y-1, s.of.y/sqrt(n.of.y), tau.u, alpha)
	pow_norm <- sum(pow) * diff(rho)[1]
	pow <- pow/pow_norm

	#compute the norm of lkl, given its support 
	ssn <- s.of.x/sqrt(n.of.x)
	df <- n.of.x - 1
	lkl_support<-pow_support
	lkl_norm<-diff(pt(lkl_support/ssn,df))

	integral_range<-pow_support	
	
	lkl_arg<-list(n.of.x= n.of.x, s.of.x= s.of.x, norm = lkl_norm, support = lkl_support)
	pow_arg<-list(df=n.of.y-1, s.of.T=s.of.y/sqrt(n.of.y), tau.u= tau.u, alpha= alpha, norm=pow_norm, support=pow_support)

	tmp <- integrate(integrand_KL_divergence_1D, lower = integral_range[1], upper = integral_range[2], dP=dMuTOST_lkl,dQ=dMuTOST_pow,P_arg=lkl_arg,Q_arg=pow_arg)

	if (tmp$message != "OK") {
		warning(tmp$message)
	}
	
	KL_div<-tmp$value
	
	if (plot) {
		rho_lkl <- seq(lkl_support[1], lkl_support[2], length.out = 1000)
		lkl<-dMuTOST_lkl(rho_lkl, n.of.x, s.of.x, lkl_norm, lkl_support)
		df_lkl <- data.frame(x = rho_lkl, no = lkl*lkl_norm ,yes = lkl)
		df_lkl$distribution <- "summary likelihood"
		df_pow <- data.frame(x = rho, no = pow*pow_norm ,yes = pow)
		df_pow$distribution <- "ABC power"
		df <- rbind(df_pow, df_lkl)
		gdf<-melt(df,id.vars=c("x","distribution"))
		p <- ggplot(data = gdf, aes(x = x, y = value, colour = distribution,linetype=variable))
		p<-p+geom_vline(xintercept=c(-tau.u,tau.u),linetype="dotted")
		p<-p+geom_hline(yintercept= mx.pw,linetype="dotted")
		p<-p+ geom_line()
		p<-p+scale_linetype("truncated?")
		p<-p+xlab(expression(rho))+ylab("")
		p <- p + ggtitle(paste("n.of.y=", n.of.y, "\ntau.u=", tau.u, "\nKL=", KL_div))
		print(p)
	}

	pw.cmx <- ifelse(calibrate.tau.u, pw.cmx, dMuTOST_pow(rho=0, n.of.y-1, s.of.y/sqrt(n.of.y), tau.u, alpha))
	
	ans <- c(KL_div = KL_div, tau.u = tau.u, pw.cmx = pw.cmx)
	return(ans)
}

#' Integrand of the Kullback-Leibler divergence D(P||Q) for generic distributions
#' @param x value at which integrand is evaluated. Can be a vector.
#' @param dP,dQ name of the functions that compute the density of P and Q
#' @param P_arg,Q_arg list of arguments for \code{dP} and \code{dQ}
#' @export
#'
integrand_KL_divergence_1D<-function(x,dP,dQ,P_arg,Q_arg){
	#this function must accept x as a vector
	
	log_P_x<-do.call(dP,c(list(x,log=T),P_arg))				
	log_Q_x<-do.call(dQ,c(list(x,log=T),Q_arg))
		
	if(any(tmp<-(log_Q_x == -Inf))){
		warning("Q(x) is 0 at some point! Usually this happens due to numerical inacurracy in the tail of Q. For simplicity we assume that log(Q(x))=-100 at these points.")
		log_Q_x[tmp]<- -100
	}
	
	ans<-(log_P_x-log_Q_x)*exp(log_P_x)

	if(any(!is.finite(ans))){
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
KL_divergence_optimize_1D <- function(x_value, x_name, is_integer=FALSE, KL_divergence, args, verbose=FALSE) {

	if (verbose) {
		print(x_value)
	}

	if(is_integer){
		x_value <- round(x_value)	
	}

	args <- c(list(x_name = x_value),args)
	names(args)[1] <- x_name
	
	if(!verbose){
	}
	
	tmp<-do.call(KL_divergence,args)
	
	if(!verbose){
		options(warn=0)
	}
	
	return(tmp["KL_div"])
}
