#=============================================================================================================================

ratetest.rejectint <- function(alpha=0.01, tau.l, tau.u, m, tol=1.0e-10, itmax=100){
  error <- -alpha  # error term to record the difference
  c.l <- sqrt(tau.l * tau.u)   #start value
  while (error < 0)
    {  c.l <- c.l - 0.05
  h <- alpha + pgamma(m * c.l/tau.u, m, 1)
  c.u <- qgamma(h,m,1)*tau.u/m
  error <- pgamma(m*c.u/tau.l, m, 1) - pgamma(m * c.l/tau.l, m, 1) - alpha
}
c1L <- c.l
c1R <- c.l + 0.05
iter <- 0
while (abs(error) >= tol && iter <= itmax)
  {  iter <- iter + 1
c.l <- (c1L + c1R) / 2
h <- alpha + pgamma(m*c.l/tau.u, m, 1)
c.u <- qgamma(h, m, 1) * tau.u /m
error <- pgamma(m*c.u/tau.l, m, 1) - pgamma(m*c.l/tau.l, m, 1) - alpha
if (error <= 0) 
 c1R <- c.l   else
c1L <- c.l
}

ans <- c(c.l, c.u)
names(ans)  <- c("c.l","c.u")

return(ans)
}

#========================================================================================================================
#' @title \code{ratetest} power function
#' @description Computes the power function of the two-sided, one-sample gamma test for testing if simulated and observed 
#' summary values occur at similar rates. This test is applicable when the observed and simulated summary values are Exponentially
#' distributed, or when Exponentiality cannot be rejected. The testing problem is described in terms of the scale parameter of the Exponential distribution, the reciprocal of the
#' rate parameter. 
#' @param rho 		Vector of error quantiles
#' @param c.l		Lower boundary point of the critical region (equivalent to the lower ABC tolerance \code{epsilon^-})
#' @param c.u 		Upper boundary point of the critical region (equivalent to the upper ABC tolerance \code{epsilon^+})
#' @param m       	Number of simulated values
#' @param trafo		Parameter transformation to translate the power
#' @param norm 		Normalization constant for the truncated power function
#' @param support 	Support of the truncated power function
#' @param log 		If \code{TRUE}, the power function is returned on the log scale. 
#' @note The power function can be truncated to \code{support} and then standardized with \code{norm}.
#' If one of these is set, the other must be provided too.
#' @seealso \code{\link{ratetest.calibrate}}, \code{\link{ratetest.pow.norm}}
#' @example example/ex.ratetest.pow.R
#' @export
#' @import data.table pscl reshape2 ggplot2 ash nortest
ratetest.pow <- function(rho, c.l, c.u, m,  norm=1, trafo=1, support=c(0, Inf), log=FALSE){
  ans <- rep(0, length(rho))
  in_support <- (rho >= support[1] & rho <= support[2])  
  if(any(in_support)){
    ans[in_support] <- ( pgamma(m*c.u/rho[in_support] *trafo, m, 1) - pgamma(m*c.l/rho[in_support] *trafo, m, 1) )/norm 
  }
  if(log){
    ans<-log(ans)
  }
  return(ans)
}


#========================================================================================================================
#' @title Area under the \code{ratetest} power function
#' @export
#' @description This function computes the area under the power function \code{ratetest.pow}.
#' @inheritParams ratetest.pow 
#' @seealso \code{\link{ratetest.pow}}, \code{\link{ratetest.calibrate}}
ratetest.pow.norm<- function(c.l, c.u, m, trafo= 1, support=c(0,Inf)){

  ans <- integrate(ratetest.pow, lower=support[1], upper=support[2], c.l=c.l, c.u=c.u, m=m, norm=1, trafo=trafo, support=support, log=FALSE)
  ans$value
  
}	

#==========================================================================================================================

ratetest.sulkl<- function(rho, n.of.x, mean.x, trafo = mean.x, norm = 1, support= c(0,Inf), log=FALSE) {

  ans 				<- rho
  in_support 			<- (rho >= support[1] & rho <= support[2])
  
  # parameters of inv-gamma
  x <- rho[in_support]
  shape <- n.of.x - 1
  scale <- n.of.x

  log_dens <- shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - (scale/x)

  if (!log) {
    ans[!in_support] <- 0
    ans[in_support] <-  exp(log_dens)
  } else {
    ans[!in_support] <- -Inf
    ans[in_support] <- log_dens
  }

  return(ans)
  
}


#=========================================================================================================================

ratetest.sulkl.norm  <- function(n.of.x, mean.x, trafo = 1, support=c(0,Inf)){

  ans	<- integrate(ratetest.sulkl, lower=support[1], upper=support[2], n.of.x=n.of.x, mean.x = mean.x, trafo=trafo, norm=1, support=support, log=FALSE)	
  ans$value
  
}


#===========================================================================================================================================

ratetest.calibrate.taulow<- function(tau.up, n.of.y, alpha=0.01, rho.star=1, tol= 1e-5, max.it=100, pow_scale=1.5, verbose=FALSE){

  rho			<- seq(1/(tau.up*pow_scale), tau.up*pow_scale, len=1024) 
  tau.low.lb	<- 2/tau.up  #randomly chosen a smaller lower boundary				
  tmp			<- max.it  # control max number of iterations
  c.rho.max	<- Inf
  
  while(c.rho.max > rho.star && tmp>0){
    tmp					<- tmp-1
    tau.low.lb	<- tau.low.lb/2
    rej    <-   ratetest.rejectint(alpha=alpha, tau.l=tau.low.lb,  tau.u = tau.up, m=n.of.y , tol=1.0e-10, itmax=100)
    #if(rej[4]>tol)	stop("compute tau.low.lb: rejection region does not have level alpha within tolerance")
    pw   <- ratetest.pow(rho, rej[1], rej[2], m = n.of.y,  norm=1, trafo=1, support=c(0, Inf), log=FALSE)
    c.rho.max	<- rho[ which.max(pw) ]
    if(verbose) {cat(paste("\ntrial lower bound",tau.low.lb,"with current rho.max",c.rho.max,"critical region",rej[1],rej[2] )) } #,"error in level is",rej[4]))
}

tau.low.ub	<- ifelse(tmp+1<max.it,  2*tau.low.lb,  tau.up)
if(verbose)	cat(paste("\nregion for tau.low is",tau.low.lb,tau.low.ub))	
  error		<- 1

while(abs(error)>tol && round(tau.low.lb, digits=10)!=round(tau.low.ub, digits=10) &&  max.it>0){
  max.it	<- max.it-1
  tau.low	<- (tau.low.lb + tau.low.ub)/2
  rej    <-   ratetest.rejectint(alpha=alpha, tau.l=tau.low,  tau.u = tau.up, m=n.of.y , tol=1.0e-10, itmax=100)
    #if(rej[4]>tol)	stop("compute tau.low: rejection region does not have level alpha within tolerance")
  pw  <-  ratetest.pow(rho, rej[1], rej[2], m=n.of.y,  norm=1, trafo=1, support=c(0, Inf), log=FALSE)
    #print( c(rho[ which.max(pw) ],pw[ which.max(pw) ], tau.low.lb, tau.low.ub,round(tau.low.lb,digits=10)==round(tau.low.ub,digits=10) ))	
  error	<- rho[ which.max(pw) ] - rho.star
  if(verbose) {cat(paste("\ntrial tau.l=",tau.low,"pw.max is",max(pw),"at",rho[ which.max(pw) ], "it",max.it)) }	
  if(error<0)
    tau.low.lb<- tau.low
  else
    tau.low.ub<- tau.low			
}

if(max.it==0)	warning("vartest.calibrate.taulow: reached maximum iteration")

  c(tau.low=tau.low, cl=rej[1], cu=rej[2], error=error)
}




#==========================================================================================================================

ratetest.calibrate.tauup <- function(mx.pw, tau.up.ub, n.of.y, alpha=0.01, rho.star=1, tol= 1e-5, max.it=100, pow_scale=1.5, verbose=FALSE){

  tau.low		<- cl <- cu	<- NA
  error		<- curr.mx.pw	<- 0
  tau.up.ub	<- tau.up.ub/2
  tmp		<- max.it		
  
  while(curr.mx.pw < mx.pw && tmp>0){
    tmp		<- tmp-1
    tau.up.ub 	<- 2*tau.up.ub
    
    calcu1 <- ratetest.calibrate.taulow(tau.up.ub, n.of.y, alpha=alpha, rho.star=rho.star, tol= tol, max.it=max.it)
    tau.low <- calcu1[1]
    cl <- calcu1[2]
    cu <- calcu1[3]
    error <- calcu1[4]

    rho		<- seq(tau.low/pow_scale, tau.up.ub*pow_scale, len=1024)   
    
    pw							<- ratetest.pow(rho, cl, cu, m=n.of.y)
    curr.mx.pw					<- max(pw)		
    if(verbose){cat(paste("\ntrial upper bound",tau.up.ub,"with power",curr.mx.pw,"at rho=",rho[ which.max(pw) ]))}
  }
  
  if(tmp==0){stop("could not find tau.up.ub")}
  if(verbose){cat(paste("\nFound upper bound",tau.up.ub,"with power",curr.mx.pw,"at rho=",rho[ which.max(pw) ]))}
  
  tau.up.lb	<- 1  # randomly chosen the lower boundary for tau.up
  error		<- 1	
  
  while(abs(error)>tol && round(tau.up.lb,digits=10)!=round(tau.up.ub,digits=10) && max.it>0){
    max.it 	<- max.it-1
    tau.up	<- (tau.up.lb + tau.up.ub)/2
    
    calcu2 <-	ratetest.calibrate.taulow(tau.up, n.of.y, alpha=alpha, rho.star=rho.star, tol= tol, max.it=max.it)
    tau.low <- calcu2[1]
    cl <- calcu2[2]
    cu <- calcu2[3]
    error <- calcu2[4]
    
    rho		<- seq(tau.low/pow_scale, tau.up*pow_scale, len=1024)		
    pw		<- ratetest.pow(rho, cl, cu, m=n.of.y)
    curr.mx.pw  	<- max(pw)		
    error 	<- curr.mx.pw - mx.pw  
    if(verbose){cat(paste("\ntrial tau.u",tau.up,"with power",curr.mx.pw,"at rho=",rho[ which.max(pw) ],"search interval",tau.up.lb,tau.up.ub,"error",error))}
    if(error<0)
      tau.up.lb<- tau.up
    else
      tau.up.ub<- tau.up
    #print(c(abs(error), round(tau.up.lb,digits=10)!=round(tau.up.ub,digits=10)) )	
  }
  if(max.it==0)	warning("vartest.calibrate.tauup: reached max.it")

    c(tau.low=tau.low, tau.up=tau.up, curr.mx.pw=curr.mx.pw,	error=abs(error), cl=cl, cu=cu)
}



#==========================================================================================================================

ratetest.getkl <- function(mean.x, n.of.x, n.of.y, tau.u, mx.pw=0.9, alpha=0.01, pow_scale=1.5, calibrate.tau.u=TRUE, plot = FALSE, legend.title=''){

  tau.l<- pw.cmx<- error<- c.l<- c.u<- NA	
  
  if(calibrate.tau.u) {	
    #calibrate tau.u constrained on yn, alpha and mx.pw 	
    res.taup <-	ratetest.calibrate.tauup(mx.pw = mx.pw, tau.u, n.of.y, alpha=alpha)	#tau.u is taken as upper bound on calibrated tau.u
    tau.l <- res.taup[1]
    tau.u <- res.taup[2]
    pw.cmx <- res.taup[3]
    error <- res.taup[4]
    c.l <- res.taup[5]
    c.u <- res.taup[6]

    if (abs(pw.cmx - mx.pw) > 0.09) {
      stop("tau.up not accurate")         
    }

  } else {
    res.taulow <-	ratetest.calibrate.taulow(tau.u, n.of.y, alpha=alpha )	#tau.u is taken as final tau.u
    tau.l <- res.taulow[1]
    c.l <- res.taulow[2]
    c.u <- res.taulow[3]
    error <- res.taulow[4]  
  }

  #truncate pow and compute pow_norm	
  pow_support <- c(tau.l/pow_scale, tau.u*pow_scale)  # this is to decide the truncated part of power function 	
  pow_norm 	<- ratetest.pow.norm(c.l, c.u, m = n.of.y, trafo=1, support=pow_support)

  #compute the norm of lkl (likelihood), given its support 
  lkl_support	<- pow_support	
  #print(c(n.of.x, s.of.x, (n.of.x-1)/n.of.x*s.of.x*s.of.x)); print(lkl_support)
  lkl_norm	<- ratetest.sulkl.norm(n.of.x, mean.x, trafo=mean.x, support=lkl_support)
  integral_range	<- pow_support	

  lkl_arg			<- list(n.of.x=n.of.x, mean.x = mean.x, trafo=mean.x, norm = lkl_norm, support= lkl_support)
  pow_arg			<- list(c.l=c.l, c.u=c.u, m=n.of.y, norm=pow_norm, trafo=1, support=pow_support)	
  tmp 			<- integrate( kl.integrand, lower = integral_range[1], upper = integral_range[2], dP=ratetest.sulkl, dQ=ratetest.pow, P_arg=lkl_arg, Q_arg=pow_arg)
  KL_div			<- tmp$value

  if (tmp$message != "OK"){
    warning(tmp$message)
  }

  if (plot) {

    rho         <- seq(lkl_support[1], lkl_support[2], length.out = 1000)

    lkl         <- ratetest.sulkl(rho =  rho, n.of.x = n.of.x, mean.x = mean.x, trafo=mean.x, norm = lkl_norm, support= lkl_support, log=FALSE)    
    df_lkl        <- data.frame(x = rho, no = lkl * lkl_norm, yes = lkl, distribution = "summary likelihood")

    pow <- ratetest.pow(rho = rho, c.l = c.l, c.u = c.u, m = n.of.y,  norm = pow_norm, trafo=1, support=pow_support, log=FALSE)
    df_pow        <- data.frame(x = rho, no = pow * pow_norm, yes = pow, distribution = "ABC power")
    
    df          <- rbind(df_pow, df_lkl)
    gdf         <- melt(df, id.vars = c("x", "distribution"))

    p           <- ggplot(data = gdf, aes(x = x, y = value, colour = distribution, linetype = variable))
    p           <- p + geom_vline(xintercept = c(tau.l, tau.u), linetype = "dotted")
    p           <- p + geom_hline(yintercept = mx.pw, linetype = "dotted")
    p           <- p + geom_line()
    p           <- p + scale_linetype("truncated and\nstandardized?")
    p           <- p + xlab(expression(rho)) + ylab("")
    p           <- p + ggtitle(paste(ifelse(!is.na(legend.title), legend.title, ''),"n.of.y=", n.of.y, "\ntau.u=", tau.u, "\nKL=", KL_div))
    print(p)

  }
  
  pw.cmx 	<- ifelse(calibrate.tau.u, pw.cmx, ratetest.pow(rho=1, c.l, c.u, m=n.of.y))	
  
  return(c(KL_div = KL_div, tau.l = tau.l, tau.u = tau.u, c.l = c.l, c.u = c.u, pw.cmx = pw.cmx))

}



#=========================================================================================================================================

ratetest.calibrate.kl <- function(n.of.x, mean.x, n.of.y, mx.pw=0.9, alpha=0.05, max.it=100, debug =FALSE, plot=FALSE ){

  KL.of.yn_ub <- KL.of.yn <- error <- curr.mx.pw <- tau.low <- cl <- cu	<- NA		
  
  #KL for initial n.of.y
  KL.of.yn <- ratetest.getkl(mean.x, n.of.x, n.of.y, tau.u = 3*mean.x, mx.pw=mx.pw, alpha=alpha, pow_scale=1.5, calibrate.tau.u=TRUE, plot = FALSE)["KL_div"] 
  
  #KL always decreases from n.of.x. Find upper bound yn.ub such that KL first increases again.	
  curr.it 		<- max.it
  yn.ub 			<- 2 * n.of.y	  # start value 	
  KL.of.yn_ub <- ratetest.getkl(mean.x, n.of.x, n.of.y=yn.ub, tau.u = 3*mean.x, mx.pw = mx.pw, alpha = alpha, pow_scale=1.5, calibrate.tau.u=TRUE, plot=FALSE)["KL_div"]  	
  

  while (KL.of.yn_ub < KL.of.yn && curr.it > 0){

    #print(c(yn.ub, KL.of.yn_ub, KL.of.yn, curr.it))
    curr.it 		<- curr.it - 1
    KL.of.yn 		<- KL.of.yn_ub
    yn.ub 			<- 2 * yn.ub
    KL.of.yn_ub		<- ratetest.getkl(mean.x, n.of.x, n.of.y=yn.ub, tau.u = 3*mean.x, mx.pw = mx.pw, alpha = alpha, pow_scale=1.5, calibrate.tau.u=TRUE, plot=FALSE)["KL_div"]
    
    if(debug){cat(paste("\ntrial upper bound m=",yn.ub,"with KL",KL.of.yn_ub))}
  }			
  
  if (curr.it == 0) {stop("could not find upper bound for yn")}		
  
  if(debug){cat(paste("\nFound upper bound m=",yn.ub,"with KL",KL.of.yn_ub))}
  
  yn.lb	<- ifelse(curr.it==max.it, yn.ub/2, yn.ub/4)
  
  if(debug){cat(paste("\nupper and lower bounds on m:",yn.lb, yn.ub))}
  
  
  KL_args <- list(mean.x=mean.x, n.of.x=n.of.x, tau.u = 3*mean.x, mx.pw = mx.pw, alpha = alpha, calibrate.tau.u=TRUE, plot=FALSE)
  
  tmp 		<- optimize(kl.optimize, interval = c(yn.lb, yn.ub), x_name = "n.of.y", is_integer = TRUE, KL_divergence = "ratetest.getkl", KL_args = KL_args, verbose = debug, tol = 1)
  
  n.of.y 	<- round(tmp$minimum)
  
  klres  <-	ratetest.getkl(mean.x, n.of.x=n.of.x, n.of.y = n.of.y, tau.u=3*mean.x, mx.pw=mx.pw, alpha=alpha, pow_scale=1.5, calibrate.tau.u=TRUE, plot=plot)
  KL_div <- klres[1]
  tau.l <- klres[2]
  tau.u <- klres[3]
  c.l <- klres[4]
  c.u <- klres[5]
  pw.cmx <- klres[6]
  
  c(n.of.y=n.of.y, tau.l=tau.l, tau.u=tau.u, cl=c.l, cu=c.u, pw.cmx=pw.cmx, KL_div=KL_div)		
}


#=============================================================================================================

ratetest.plot<- function(n.of.y, c.l, c.u, tau.l, tau.u, pow_scale=1.5){

  pow_support <- c(tau.l/pow_scale, tau.u*pow_scale) 	
  pow_norm 	<- ratetest.pow.norm(c.l=c.l, c.u=c.u, m=n.of.y, trafo=1, support=pow_support)	
  
  tmp			<- data.frame(rho=seq(pow_support[1], pow_support[2], length.out = 1024))	
  tmp$power	<- ratetest.pow(tmp$rho, c.l=c.l, c.u=c.u, m=n.of.y, norm=pow_norm, trafo= 1)*pow_norm	
  
  p	<- ggplot(tmp, aes(x=rho, y=power)) + geom_line() + labs(x=expression(rho), y='Power\n(ABC acceptance probability)') +
  scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
  scale_x_continuous(limits=c(0,pow_support[2])) +
  geom_vline(xintercept = c(tau.l, tau.u), linetype = "dotted") +
  geom_vline(xintercept = c(c.l, c.u), linetype = "dashed") +
  ggtitle(paste("n.of.y=", n.of.y, "\ntau.l=", round(tau.l,d=5), " tau.u=", round(tau.u,d=5), "\nc.l=", round(c.l,d=5), " c.u=", round(c.u,d=5)))
  print(p)
}




#==========================================================================================================
#' @title Calibrating \code{ratetest} for ABC
#' @description Calibrate the one-sample equivalence test for testing if simulated and observed 
#' summary values occur at similar rates. This test is applicable when the observed and simulated summary values are Exponentially
#' distributed, or when Exponentiality cannot be rejected. The testing problem is described in terms of the scale parameter of the Exponential distribution, the reciprocal of the
#' rate parameter.
#'  
#' Different types of calibrations are available, see Notes for details:
#' \enumerate{ 
#'  \item (\code{what=ALPHA}) compute the ABC false positive rate for given critical region,
#'  \item (\code{what=CR}) calibrate the critical region for given ABC false positive rate,
#'  \item (\code{what=MXPW}) calibrate the critical region and the equivalence region for given ABC false positive rate and maximum power,
#'  \item (\code{what=KL}) calibrate the critical region, the equivalence region and the number of simulated summary values for given ABC false positive rate, maximum power and sample standard deviation of the observed data.
#' }
#' 
#' The calibration KL should be used. Typically, the KL calibration requires multiple i. i. d. instances of observed summary statistics
#' at each ABC iteration. Here, the maximum likelihood estimate of the scale parameter is the sample mean, and can be computed from a single 
#' summary value. However, based on a single summary value, it is not possible to determine if the Exponential distribution is appropriate.
#' It is recommended to use this test with multiple summary values, and to evaluate on the fly if the Exponential distribution is appropriate. 
#' 
#' Depending on the type of calibration, some of the following inputs must be specified (see Examples).
#' @export 
#' @import data.table pscl reshape2 ggplot2 ash nortest
#' @param n.of.x   	Number of observed summary values 
#' @param n.of.y 	Number of simulated summary values
#' @param mean.x  	Mean of observed summary values 
#' @param what		Character string to indicate the type of calibration to be performed
#' @param c.l		Lower boundary point of the critical region 
#' @param c.u		Upper boundary point of the critical region
#' @param tau.l		Lower boundary point of the equivalence region 
#' @param tau.u		Upper boundary point of the equivalence region 
#' @param tau.u.ub	Guess on the upper boundary point of the equivalence region 
#' @param mx.pw 	Maximum power at the point of equality
#' @param alpha 	Level of the equivalence test
#' @param max.it  	Maximum number of optimization steps at each calibration hierarchy 
#' @param pow_scale   Scale for the support of the standardized power. The power is truncated to \code{[tau.l/pow_scale,tau.u*pow_scale]} and then standardized
#' @param tol		Required error tolerance in calibrating the actual maximum power to the requested maximum power 
#' @param plot  	Flag to plot calibrations
#' @param debug		Flag to switch off C implementation
#' @param plot_debug	Flag to plot at each calibration iteration
#' @param verbose   	Flag to run in verbose mode
#' @return	vector
#' @seealso \code{\link{mutost.calibrate}}, \code{\link{ztest.calibrate}}, \code{\link{vartest.calibrate}}
#' @note 
#' \enumerate{	
#'  \item (\code{what=ALPHA}) This calibration requires the inputs \code{n.of.y}, \code{c.l}, \code{c.u}, \code{tau.l}, \code{tau.u} with \code{c.l>tau.l}, \code{c.u<tau.u}, \code{tau.u>1}, \code{tau.l<1}. 
#' 				The output contains the corresponding ABC false positive rate \code{alpha}.
#' 				This option does not specify any of the free ABC parameters, but may be useful to determine the ABC
#' 				false positive rate for uncalibrated ABC routines.
#'  \item (\code{what=CR}) This calibration requires the inputs \code{tau.l}, \code{tau.u}, \code{alpha} with \code{tau.l<1}, \code{tau.u>1} and default \code{alpha=0.01}. 
#' 				The output contains the corresponding critical region \code{[c.l, c.u]}, which corresponds to the ABC tolerance region typically denoted by \code{[-epsilon, epsilon]}. 
#' 				This is an intermediate calibration step and may result in unsuitable power properties (see Examples).				
#'  \item (\code{what=MXPW}) This calibration requires the inputs \code{alpha} and \code{mx.pw}, with default values 0.01 and 0.9 respectively. 
#' 				The output contains the corresponding critical region \code{[c.l, c.u]} (to be used in ABC, see Notes on (2)), and 
#' 				the corresponding equivalence region \code{[tau.l, tau.u]} that gives a suitable ABC accept/reject probability if the simulated 
#' 				summary values are close to the observed summary values.
#' 				As a check to the numerical calibrations, the actual power at the point of equality is returned (\code{pw.cmx}). 
#' \item (\code{what=KL}) This calibration can be used when a set of observed summary values is available. It is desirable because it specifies the number of simulated summary 
#' 				values so that the power is very close to the desired summary likelihood in terms of the KL divergence. 
#' 				The inputs are \code{alpha} (default is 0.01), \code{mx.pw} (default is 0.9), \code{n.of.x}, \code{mean.x}, \code{n.of.y}. 
#' 				The output consists of the corresponding critical region \code{[c.l, c.u]} (to be used in ABC, see Notes on (2)), the equivalence
#' 				region \code{[tau.l, tau.u]}, and the number of simulated summary values needed (\code{n.of.y}). As a check to the numerical calibrations, 
#' 				the KL divergence is returned (\code{KL}). It is desirable to compare the power to the summary likelihood in terms of the KL divergence, see References.
#' }
#' Note that the underlying test statistic only depends on \code{n.of.x}, \code{n.of.y}, \code{mean.x}, which are all 
#' known before ABC is run. Consequently, the free ABC parameters are calibrated once, before ABC is started. 
#' 
#' The lower boundary point of the critical region \code{c.l} is calibrated numerically, so that
#' the power is maximized at the point of equality \code{rho=1}. The calibrated \code{c.l} does not equal 1/\code{c.u}. 
#' @example example/ex.ratetest.calibrate.R
#' @references  http://arxiv.org/abs/1305.4283
ratetest.calibrate<- function(n.of.x=NA, mean.x=NA, n.of.y=NA, what='MXPW', mx.pw=0.9, alpha=0.01, tau.l=NA, tau.u=NA, tau.u.ub=NA, c.l=NA, c.u=NA, max.it=100, tol= 1e-5, pow_scale=1.5, debug=FALSE, plot=FALSE, verbose=FALSE)
{

  stopifnot(what%in%c('ALPHA','CR','MXPW_AT_EQU','MXPW','KL'))
  
  if(what=='ALPHA'){

    stopifnot(c.u<=tau.u, tau.u>1, tau.l<1, c.l>=tau.l, n.of.y>2, alpha>0, alpha<1)
    
    ans	<- pgamma(n.of.y*c.u/tau.l, n.of.y, 1)-pgamma(n.of.y*c.l/tau.l, n.of.y, 1)
    names(ans)	<- 'alpha'
    if(plot)
      ratetest.plot(n.of.y, c.l, c.u, tau.l, tau.u, pow_scale=pow_scale)

  }
  
  if(what=='CR'){

    stopifnot(tau.u>1, tau.l<1, n.of.y>2, alpha>0, alpha<1, pow_scale>1)
    
    tmp <- ratetest.rejectint(alpha, tau.l, tau.u, n.of.y, tol=1e-10, itmax=max.it)
    ans <- c(tmp[1], tmp[2])
    names(ans)<- c('c.l','c.u')
    
    if(plot)
      ratetest.plot(n.of.y, ans['c.l'], ans['c.u'], tau.l, tau.u, pow_scale=pow_scale)		

  }
  
  if(what=='MXPW_AT_EQU'){   

    stopifnot(tau.u>1, n.of.y>2, alpha>0, alpha<1, pow_scale>1, max.it>10, tol<0.2)
    
    tmp	<- ratetest.calibrate.taulow(tau.u, n.of.y, alpha=alpha, rho.star=1, tol=tol, max.it=max.it, pow_scale=pow_scale, verbose=verbose)
    
    ans	<- c(tmp[2],tmp[3],tmp[1],tau.u,tmp[4])
    names(ans)<- c('c.l','c.u','tau.l','tau.u','eq.error')
    if(plot)
      ratetest.plot(n.of.y, ans['c.l'], ans['c.u'], ans['tau.l'], tau.u, pow_scale=pow_scale)		
    
  }	
  
  if(what=='MXPW'){

    stopifnot(tau.u.ub>1, n.of.y>2, alpha>0, alpha<1, pow_scale>1, max.it>10, tol<0.2)
    
    tmp <- ratetest.calibrate.tauup(mx.pw, tau.u.ub, n.of.y, alpha=alpha, rho.star=1, tol=tol, max.it=max.it, pow_scale=pow_scale, verbose=verbose)
    ans	<- c(tmp[5], tmp[6], tmp[1], tmp[2], tmp[3], tmp[4])
    names(ans)<- c('c.l','c.u','tau.l','tau.u','pw.cmx','pw.error')
    
    if(plot)
      ratetest.plot(n.of.y, ans['c.l'], ans['c.u'], ans['tau.l'], ans['tau.u'], pow_scale=pow_scale)
  }
  
  if(what=='KL'){

    tmp	<- ratetest.calibrate.kl(n.of.x, mean.x, n.of.y=n.of.x, mx.pw=mx.pw, alpha=alpha, max.it=max.it, debug=debug, plot=plot)
    ans	<- c(tmp[4], tmp[5], tmp[2], tmp[3], tmp[1], tmp[6], tmp[7])
    names(ans)	<- c('c.l','c.u','tau.l','tau.u','n.of.y','pw.cmx','KL.div')
    
    # if(plot)
      # ratetest.plot(ans['n.of.y'], ans['c.l'], ans['c.u'], ans['tau.l'], ans['tau.u'], pow_scale=pow_scale)

  }
  
  ans
}


