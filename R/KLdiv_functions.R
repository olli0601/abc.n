

#' @title Integrand of the Kullback-Leibler divergence D(P||Q) for generic distributions
#' @export
#' @param x value at which integrand is evaluated. Can be a vector.
#' @param dP name of the functions that compute the density of P
#' @param dQ name of the functions that compute the density of Q
#' @param P_arg list of arguments for \code{dP}
#' @param Q_arg list of arguments for \code{dQ}
#' @return KL divergence
#'
kl.integrand<-function(x, dP, dQ, P_arg, Q_arg)
{
	#this function must accept x as a vector	
	log_P_x	<- do.call(dP,c(list(x,log=T),P_arg))				
	log_Q_x	<- do.call(dQ,c(list(x,log=T),Q_arg))
		
	if(any(tmp<-(log_Q_x == -Inf)))
	{
		cat("\nQ(x) is 0 at some point! Usually this happens due to numerical inacurracy in the tail of Q. For simplicity we assume that log(Q(x))=-100 at these points.")
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

#' @title Compute the KL divergence for two dimensional densities P and Q
#' @param df1	data.table with two columes th1 and th2, random draws from P
#' @param df2	data.table with two columes th1 and th2, random draws from Q
#' @param nbin	number of bins for both densities
#' @return KL divergence
#' @export
kl.2D <- function(df1, df2, nbin=100)
{
	df1.lim		<- df1[, lapply(.SD, range)]
	df2.lim		<- df2[, lapply(.SD, range)]
	df.lim		<- as.matrix( rbind(df1.lim, df2.lim)[, lapply(.SD, range)] )
	df1.bins	<- bin2(as.matrix(df1), ab=t(df.lim), nbin=rep(nbin,2))
	df1.bins$nc	<- df1.bins$nc / sum( df1.bins$nc )
	df2.bins	<- bin2(as.matrix(df2), ab=t(df.lim), nbin=rep(nbin,2))
	df2.bins$nc	<- df2.bins$nc / sum( df2.bins$nc )
	
	kl			<- data.table(one= as.numeric(log( df1.bins$nc / df2.bins$nc )), two= as.numeric( log( df2.bins$nc / df1.bins$nc )) )
	set(kl, which(is.nan(kl[,one])), 'one', 0. )
	set(kl, which(is.nan(kl[,two])), 'two', 0. )
	set(kl, which(is.infinite(kl[,one])), 'one', NA )
	set(kl, which(is.infinite(kl[,two])), 'two', NA )
	set(kl, NULL, 'one', kl[, one] * as.numeric(df1.bins$nc) )
	set(kl, NULL, 'two', kl[, two] * as.numeric(df2.bins$nc) )
	set(kl, which(as.numeric(df1.bins$nc)==0), 'one', 0. )
	set(kl, which(as.numeric(df2.bins$nc)==0), 'two', 0. )
	
	ans			<- list(one= sum(kl[,one], na.rm=1), two= sum(kl[,two], na.rm=1), zero.denom.one= sum(is.na(kl[,one])), zero.denom.two= sum(is.na(kl[,two]))	)
	ans
}

#' A wrapper to minimize \code{KL_divergence} over the parameter \code{x_name} using the function \link{optimize}
#' @export
#' @param x_value value tested.
#' @param x_name name of the parameter over which minimization is performed.
#' @param is_integer if \code{TRUE}, \code{x_value} is rounded to the closest integer. See \link{round}.
#' @param KL_divergence character, name of the function that computes the KL divergence.
#' @param KL_args additional arguments to be passed to \code{KL_divergence}.
#' @param verbose logical, if \code{TRUE}, print warnings.
#' @return the minimized Kullback-Leibler divergence (scalar).
kl.optimize<- function(x_value, x_name, is_integer=FALSE, KL_divergence, KL_args, verbose=FALSE) 
{
	if (verbose)		cat(paste("\nOptimizing x=",x_value))
	if(is_integer)
		x_value 		<- round(x_value)	
	KL_args 			<- c(list(x_name = x_value),KL_args)
	names(KL_args)[1] 	<- x_name
	tmp					<- do.call(KL_divergence,KL_args)	
	if(0 || !verbose)								#TODO anton check
		options(warn=0)	
	return(tmp["KL_div"])
}