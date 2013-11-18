

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