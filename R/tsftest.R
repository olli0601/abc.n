#' @title \code{tsftest} power function 
#' @description Compute the power of the two-sided one-sample multivariate equivalence test for population means of multivariate normal summary values with unknown population variance.
#' @export
#' @import data.table ggplot2
#' @param rho 		Vector of quantiles
#' @param cl		Lower boundary point of the critical region
#' @param cu		Upper boundary point of the critical region
#' @param n.of.y 	Number of replicate simulations
#' @param p 		Number of variables
#' @param support 	Support of the truncated power function (vector of dimension 2).
#' @param log 		If \code{TRUE}, the power function is returned on the log scale. 
#' @note The power function can be truncated to \code{support}.
#' @example example/ex.tsftest.pow.R
#' @references  http://arxiv.org/abs/1305.4283
#------------------------------------------------------------------------------------------------------------------------
tsftest.pow <- function(rho, cl, cu, n.of.y, p, support=c(0, Inf), log=FALSE, norm=1)
{ 
	stopifnot(cl>0, cu>0, cl<cu, support[1] <= support[2], support[1] >= 0)
	stopifnot(n.of.y %% 1 == 0, p %% 1 == 0, n.of.y > 1, p > 0, p < n.of.y)
	ans 			<- rho
	in_support 		<- (rho >= support[1] & rho <= support[2])
	ans[!in_support]<- ifelse(log, -Inf, 0)
	if(any(in_support))
	{
		ans[in_support]		<- pf( (n.of.y-p)/(p*(n.of.y-1))*cu, p, n.of.y - p, n.of.y * rho[in_support] ) - pf( (n.of.y-p)/(p*(n.of.y-1))*cl, p, n.of.y - p, n.of.y * rho[in_support] )
		ans[in_support]		<- ans[in_support]/norm		
		tmp					<- which(ans[in_support] < 0)
		if(length(tmp))
			ans[tmp]		<- 0
		if(log)
			ans[in_support]	<- log(ans[in_support])
	}
	ans
}
