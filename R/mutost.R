#' Density of the power for muTOST
#' 
#' Compute the density of the (possible truncated) power of the equivalence test for population means of normal summary values
#' @inheritParams mutost.calibrate.tolerances.getkl
#' @param rho vector of quantile
#' @param df		degrees of freedom
#' @param s.of.T	standard deviation of the test statistic
#' @param tau.u		upper tolerance of equivalence region
#' @param alpha		level of equivalence test
#' @param norm normalization constant for the truncated power function.
#' @param support vector of dimension 2. Support of the truncated power function.
#' @param log logical; if \code{TRUE}, densities d are given as log(d). 
#' @note The summary likelihood can be truncated to \code{support} and then standardized with \code{norm}.
#' For computational efficiency, both \code{norm} and \code{support} must be provided although each one can be derived (numerically) from the other.
#' @export
#' @example example/ex.mutost.pow.R
#' 	
mutost.pow <- function(rho, df, s.of.T, tau.u, alpha, norm=1, support=c(-Inf,Inf), log=FALSE)
{	
	ans				<- rho
	in_support		<- (rho >= support[1] & rho <= support[2])
	ans[!in_support]<- 0
	
	if(any(in_support))
	{
		suppressWarnings({ #suppress numerical inaccuracy warnings
			ans[in_support] <- .Call("abcMuTOST_pow", rho[in_support], df, tau.u, s.of.T, alpha)/norm
		})
	}	
	if(log)
		ans			<- log(ans)	
	ans
}

#------------------------------------------------------------------------------------------------------------------------
#' Generic two one sided test (TOST).
#' 
#' Perform a generic two one sided test. This is an internal function.
#' @param tost.args vector of arguments for generic TOST
#' @param tau.l		lower tolerance of equivalence region
#' @param tau.u		upper tolerance of equivalence region
#' @param alpha		level of equivalence test
#' @param tost.distr	name of distribution of tost
#' @return vector of length 7
generic.tost <- function(tost.args, tau.l, tau.u, alpha, tost.distr="t")
{
	ans<- numeric(7)
	names(ans)<- c("error","p.error","lkl","cl","cu","ass.alpha","ass.pval")
	if(!tost.distr%in%c("t"))	
		stop("unexpected tost.distr")

	which.not.reject<- which( c(pt( tost.args[1], tost.args[4] )<1-alpha, pt( tost.args[2], tost.args[4] )>alpha))
	if(length(which.not.reject)%in%c(0,2))		#both upper and lower test statistics indicate mean difference < -tau and mean difference > tau  	OR 		mean difference >= -tau and mean difference <= tau
	{
		#figure out which one is closer to boundary, and schedule that this one is reported
		which.not.reject<- which.min(abs(c(1-alpha-pt( tost.args[1], tost.args[4] ), pt( tost.args[2], tost.args[4] )-alpha )))
	}
	
	ans["error"]	<- tost.args[3]
	ans["p.error"]	<- ifelse(	which.not.reject==1, 		
		1-pt( tost.args[which.not.reject], tost.args[4] ),		
		pt( tost.args[which.not.reject], tost.args[4] )		)
	ans["lkl"]		<- dt( tost.args[3], tost.args[4])		
	ans["cl"]		<- min(tau.l/tost.args[5]+qt( 1-alpha, tost.args[4] ),0)
	ans["cu"]		<- max(0, tau.u/tost.args[5]+qt( alpha, tost.args[4] ))
	ans["ass.alpha"]<- 1-diff(pt( ans[c("cl","cu")], tost.args[4]))											#the corresponding alpha quantile of the traditional t-test with acceptance region [c-,c+] and above s, df
	ans["ass.pval"]	<- ( pt( tost.args[3], tost.args[4] ) - ans["ass.alpha"]/2 ) / ( 1 - ans["ass.alpha"] )						#rescaled p-value that is expected to follow U(0,1) under the point null hypothesis
	
	#POWER[[length(POWER)+1]]<<- c(tau.l, tau.u, tmp[6], sum(moments[,1]), tmp[4], tmp[5] )	#PowerTOST:::.power.TOST(alpha=alpha, tau.l, tau.u, seq(tau.l, tau.u, length.out= 1e3), tmp[6], sum(moments[,1]), tmp[4], bk = 4)	
	ans	
}

#' Density of the summary likelihood for muTOST
#' 
#' Compute the density of the (possibly truncated) summary likelihood for population means of normal summary values
#' @inheritParams mutost.calibrate.tolerances.getkl
#' @param rho vector of quantile
#' @param norm scalar, 0<\code{norm}<=1, normalization constant for the truncated summary likelihood.
#' @param support vector of dimension 2, support of the truncated summary likelihood.
#' @param log logical; if \code{TRUE}, densities d are given as log(d). 
#' @note The summary likelihood can be truncated to \code{support} and then standardized with \code{norm}.
#' For computational efficiency, both \code{norm} and \code{support} must be provided although each one can be derived from the other.
#' \code{support=s.of.x/sqrt(n.of.x)*qt(c(1-norm,1+norm)/2,n.of.x-1)} and \code{norm=diff(pt(support/(s.of.x/sqrt(n.of.x)),n.of.x-1))}.
#' @export	
#'
mutost.sulkl <- function(rho, n.of.x, s.of.x, norm = 1, support= c(-Inf,Inf), log=FALSE, debug=0) 
{
	
	stopifnot(n.of.x>0,s.of.x>0,norm<=1,norm>0,support[1]<=support[2])
	
	ans 			<- rho
	in_support 		<- (rho >= support[1] & rho <= support[2])
	ans[!in_support]<- ifelse(log,-Inf,0)


	if(debug){
		#R code
		ssn				<- s.of.x/sqrt(n.of.x)
		df 				<- n.of.x - 1
		if (any(in_support)) 
		{
			if(log)
				ans[in_support] <- dt(rho[in_support]/ssn, df,log=TRUE)-log(ssn*norm)
			else
				ans[in_support] <- dt(rho[in_support]/ssn, df)/ssn/norm		
		}
		
	}else{
		#C code
		if (any(in_support)){
			ans[in_support]<-.Call("abcMuTOST_sulkl",rho[in_support], n.of.x, s.of.x, norm, log)
		}
		
	}
	return(ans)
}
#------------------------------------------------------------------------------------------------------------------------
#' Kullback-Leibler divergence for the muTOST
#'
#' Compute Kullback-Leibler divergence between the summary likelihood and the power function of the test of location equivalence 
#' @param n.of.x number of observed summary values
#' @param s.of.x standard deviation of observed summary values
#' @param n.of.y number of simulated summary values
#' @param s.of.y standard deviation of simulated summary values
#' @param mx.pw maximum power at the point of reference (rho.star=0) (only when \code{calibrate.tau.u==TRUE}).
#' @param alpha level of the equivalence test
#' @param calibrate.tau.u if \code{TRUE} the upper tolerance of the equivalence region (\code{tau.u}) is calibrated so that power at the point of reference is equal to \code{mx.pw}
#' @param tau.u	upper tolerance of the equivalence region. If \code{calibrate.tau.u==TRUE}, it's just a guess on an upper bound on the upper tolerance of the equivalence region.
#' @param pow_scale scale for the support of the standardized power. The power is truncated to \code{pow_scale*[-tau.u,tau.u]} and then standardized.
#' @param debug flag if C implementation is used
#' @param plot  logical; if \code{TRUE} the summary likelihood and calibrated power are ploted.
#' @param legend.title character; title of the plot (only when \code{plot==TRUE}).
#' @return	vector of length 4
#' 	\item{KL_div}{the Kullback Leibler divergence}		
#' 	\item{n.of.y}{number of simulated summary values}		
#' 	\item{tau.u}{upper tolerance of the equivalence region}
#' 	\item{pw.cmx}{actual maximum power associated with the equivalence region}
#' @export
#' @import ggplot2 reshape2
#' @examples
#' 
#' mutost.calibrate.tolerances.getkl(n.of.x=60,s.of.x=0.1,n.of.y=60,s.of.y=0.3, mx.pw=0.9, alpha=0.01, calibrate.tau.u=TRUE, tau.u=1, plot=TRUE)
#'
mutost.calibrate.tolerances.getkl <- function(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = F, tau.u = 1, pow_scale = 1.5, debug = 0, plot = F, legend.title='') 
{
	

	stopifnot(n.of.x > 1, s.of.x > 0, n.of.y > 1, s.of.y > 0, mx.pw > 0, mx.pw<=1, alpha > 0, alpha<=0.5, tau.u>0, pow_scale > 0)

	if (!debug)		#ALL IN C 
	{						
		suppressWarnings({ #suppress numerical inaccuracy warnings
			ans <- .Call("abc_mutost_get_KL", n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u, tau.u, pow_scale)
		})

		KL_div <- ans[1]
		tau.u <- ans[2]
		pw.cmx <- ans[3]
		
		if(plot)
		{
			ssn 			<- s.of.x/sqrt(n.of.x)
			lkl_support 	<- pow_support <- c(-tau.u, tau.u) * pow_scale
			lkl_norm 		<- diff(pt(lkl_support/ssn, n.of.x-1))
			suppressWarnings({ #suppress numerical inaccuracy warnings
				pow_norm 	<- .Call("abc_mutost_integrate_pow", pow_support[1], pow_support[2],.Machine$double.eps^0.25,.Machine$double.eps^0.25,as.double(n.of.y-1),s.of.y/sqrt(n.of.y),tau.u,alpha,1,0)
			})			
		}
	} 
	else			#ALL IN R 
	{		
		if (calibrate.tau.u) 
		{
			#calibrate tau.u constrained on yn, alpha and mx.pw	
			tmp 	<- mutost.calibrate.tolerances(mx.pw, n.of.y - 1, s.of.y/sqrt(n.of.y), tau.u, alpha)
			tau.u 	<- tmp[2]
			pw.cmx 	<- tmp[3]
			if (abs(pw.cmx - mx.pw) > 0.09) 	
				stop("tau.up not accurate")
		}
		#truncate pow and compute pow_norm
		pow_support <- c(-tau.u, tau.u) * pow_scale
		#pow_norm <- integrate(mutost.pow, lower = pow_support[1], upper = pow_support[2], df=n.of.y-1, s.of.T=s.of.y/sqrt(n.of.y), tau.u= tau.u, alpha= alpha, norm=1, support= pow_support, log=FALSE)
		suppressWarnings({ #suppress numerical inaccuracy warnings
			pow_norm <- .Call("abc_mutost_integrate_pow", pow_support[1], pow_support[2],.Machine$double.eps^0.25,.Machine$double.eps^0.25,as.double(n.of.y-1),s.of.y/sqrt(n.of.y),tau.u,alpha,1,0)
		})		
		#compute the norm of lkl, given its support 
		ssn 			<- s.of.x/sqrt(n.of.x)
		df 				<- n.of.x - 1
		lkl_support 	<- pow_support
		lkl_norm 		<- diff(pt(lkl_support/ssn, df))
		integral_range	<- pow_support
		lkl_arg 		<- list(n.of.x = n.of.x, s.of.x = s.of.x, norm = lkl_norm, support = lkl_support)
		pow_arg 		<- list(df = n.of.y - 1, s.of.T = s.of.y/sqrt(n.of.y), tau.u = tau.u, alpha = alpha, norm = pow_norm, support = pow_support)
		tmp 			<- integrate(kl.integrand, lower = integral_range[1], upper = integral_range[2], dP = mutost.sulkl, dQ = mutost.pow, P_arg = lkl_arg, Q_arg = pow_arg)
		KL_div 			<- tmp$value
		if (tmp$message != "OK")	
			warning(tmp$message)
		pw.cmx 			<- ifelse(calibrate.tau.u, pw.cmx, mutost.pow(rho = 0, n.of.y - 1, s.of.y/sqrt(n.of.y), tau.u, alpha))
		#print(c(calibrate.tau.u, n.of.y, s.of.y, s.of.y/sqrt(n.of.y), alpha, tau.u, pw.cmx, mutost.pow(rho = 0, n.of.y - 1, s.of.y/sqrt(n.of.y), tau.u, alpha)))		
	}

	if (plot) 
	{
		rho 				<- seq(lkl_support[1], lkl_support[2], length.out = 1000)
		lkl 				<- mutost.sulkl(rho, n.of.x, s.of.x, lkl_norm, lkl_support)
		df_lkl 				<- data.frame(x = rho, no = lkl * lkl_norm, yes = lkl)
		df_lkl$distribution <- "summary likelihood"
		pow					<- mutost.pow(rho, df=n.of.y-1, s.of.T=s.of.y/sqrt(n.of.y), tau.u, alpha, norm=pow_norm, support= pow_support, log=FALSE)
		df_pow 				<- data.frame(x = rho, no = pow * pow_norm, yes = pow)
		df_pow$distribution <- "ABC power"
		df 					<- rbind(df_pow, df_lkl)
		gdf 				<- melt(df, id.vars = c("x", "distribution"))
		p 					<- ggplot(data = gdf, aes(x = x, y = value, colour = distribution, linetype = variable))
		p 					<- p + geom_vline(xintercept = c(-tau.u, tau.u), linetype = "dotted")
		p 					<- p + geom_hline(yintercept = mx.pw, linetype = "dotted")
		p 					<- p + geom_line()
		p 					<- p + scale_linetype("truncated and\nstandardized?")
		p 					<- p + xlab(expression(rho)) + ylab("")
		p 					<- p + ggtitle(paste(ifelse(!is.na(legend.title), legend.title, ''),"n.of.y=", n.of.y, "\ntau.u=", tau.u, "\nKL=", KL_div))
		print(p)
	}
	ans <- c(KL_div = KL_div, n.of.y=n.of.y, tau.u = tau.u, pw.cmx = pw.cmx)
	ans
}
#------------------------------------------------------------------------------------------------------------------------
#'Calibration for muTOST
#'
#'Calibrate the equivalence region for the test of location equivalence for given maximum power
#' @param KL_args list of arguments needed to run \code{\link{mutost.calibrate.tolerances.getkl}} (see note below)
#' @param  max.it  this algorithm stops prematurely when the number of iterations to find the equivalence region exceeds 'max.it'
#' @param  debug  flag if C implementation is used
#' @param  plot  logical; if \code{TRUE} the summary likelihood and calibrated power are ploted (only when \code{debug==TRUE}).
#' @param  plot_debug  logical; if \code{TRUE} the summary likelihood and power are ploted at each calibration step (only when \code{debug==TRUE}).
#' @param  verbose flag if detailed information on the computations should be printed to standard out
#' @export
#' @note \code{KL_args} is a named list with the following elements: \code{n.of.x}, \code{s.of.x}, \code{n.of.y}, \code{s.of.y}, \code{mx.pw}, \code{alpha}, \code{tau.u} and \code{pow_scale}.
#' @seealso \code{\link{mutost.calibrate.tolerances.getkl}}
#' @return	vector of length 4
#' 	\item{KL_div}{the Kullback Leibler divergence}		
#' 	\item{n.of.y}{number of simulated summary values}		
#' 	\item{tau.u}{upper tolerance of the equivalence region}
#' 	\item{pw.cmx}{actual maximum power associated with the equivalence region}
#' @example example/ex.mutost.calibrate.R
mutost.calibrate<- function(KL_args, max.it = 100, debug = 0, plot = FALSE, plot_debug = FALSE, verbose=FALSE) 
{	

	stopifnot(max.it > 0)
	stopifnot(all(c("n.of.x", "s.of.x", "n.of.y", "s.of.y", "mx.pw", "alpha", "tau.u", "pow_scale") %in% names(KL_args)))
	with(KL_args, stopifnot(n.of.x > 1, s.of.x > 0, n.of.y > 1, s.of.y > 0, mx.pw > 0, mx.pw<=1, alpha > 0, alpha<=0.5, tau.u>0, pow_scale > 0))
	#print(KL_args)
	
	#see if power function "too tight" for yn=xn -> as a marker we use that KL(yn-1) < KL(yn). 
	#In this case, 	the only way to minimize KL is to use yn<xn, which is not what we like.
	#				instead, we give up on mx.pw=0.9	
	KL_divergence			<-  "mutost.calibrate.tolerances.getkl"		
	KL_args$calibrate.tau.u <- T
	KL_args$plot 			<- F			
	KL.of.yn 				<- do.call(KL_divergence, KL_args)["KL_div"]

	
	n.of.y 					<- KL_args$n.of.y	
	KL_args$n.of.y 			<- n.of.y - 1
	KL.of.yn_m1 			<- do.call(KL_divergence, KL_args)["KL_div"]		
	decrease_n.of.y			<- as.logical(KL.of.yn_m1 < KL.of.yn)
	if(verbose)	
		cat(paste("\ninitial m=",n.of.y, "KL(m)=",KL.of.yn, "KL(m-1)=",KL.of.yn_m1))
	KL_args$n.of.y 			<- n.of.y	
	
	
	if (!debug)		#all in C 
	{		
		if(!decrease_n.of.y)		#case power function not "too tight", adjust yn for fixed mx.pw
		{
			suppressWarnings({ #suppress numerical inaccuracy warnings
				ans 	<- .Call("abc_mutost_calibrate_powerbroader", KL_args, as.integer(max.it))
			})
		}
		else						#case power function "too tight", force yn=xn and give up mx.pw
		{
			KL_args$tau.u 	<- as.double(KL_args$s.of.y / 5)
			stopifnot(KL_args$tau.u>0)
			suppressWarnings({ #suppress numerical inaccuracy warnings
				ans <- .Call("abc_mutost_calibrate_powertighter", KL_args, as.integer(max.it))
			})
		}
	} 
	else			#all in R
	{
		if(!decrease_n.of.y){
			#case power function not "too tight", adjust yn for fixed mx.pw 
			#we have KL(yn-1)>KL(yn), ie KL decreases as yn-1 is incremented. Find upper bound yn.ub such that KL first increases again. 						
			curr.it 		<- max.it
			yn.ub 			<- 2 * (KL_args$n.of.y - 1)
			KL_args$n.of.y 	<- yn.ub
			KL.of.yn_ub 	<- do.call(KL_divergence, KL_args)["KL_div"]
			while (KL.of.yn_ub < KL.of.yn && curr.it > 0) 
			{
					#print(c(yn.ub, KL.of.yn_ub, KL.of.yn, curr.it))
				curr.it 		<- curr.it - 1
				KL.of.yn 		<- KL.of.yn_ub
				yn.ub 			<- 2 * yn.ub
				KL_args$n.of.y 	<- yn.ub
				KL.of.yn_ub 	<- do.call(KL_divergence, KL_args)["KL_div"]				
			}			
			if (curr.it == 0) 	
				stop("could not find upper bound for yn")					
			yn.lb 				<-ifelse(curr.it == max.it,KL_args$n.of.y/2,KL_args$n.of.y/4)		#previous 'yn.ub' could be exactly the minimum, so need to choose 'yn.ub/4' as lower bound 
			#KL is minimized in the open set (yn.lb,yn.ub)						
			if(verbose)	
				cat(paste("\nupper and lower bounds on m:",yn.lb, yn.ub))
			if (plot_debug) 
				pdf("KL_optimization.pdf", onefile = T)		
			KL_args$plot 			<- plot_debug
			KL_args["n.of.y"] 		<- NULL
			tmp 					<- optimize(kl.optimize, interval = c(yn.lb, yn.ub), x_name = "n.of.y", is_integer = T, KL_divergence = KL_divergence, KL_args = KL_args, verbose = debug, tol = 1)
			if (plot_debug) 
				dev.off()
			KL_args$n.of.y 			<- round(tmp$minimum)
			KL_args$plot 			<- plot
			ans 					<- do.call(KL_divergence, KL_args)
		}
		else					#case power function "too tight", force yn=xn and give up mx.pw
		{
			KL_args$calibrate.tau.u <- F
			KL_args$tau.u 			<- KL_args$s.of.y / 5
			# print(KL_args$tau.u)
			KL_args$plot 			<- plot_debug
			if (plot_debug)
				pdf("KL_initial.pdf", onefile = T)
			KL.of.tau.u 			<- do.call(KL_divergence, KL_args)["KL_div"]
			if (plot_debug)
				dev.off()

			#to find tau.u.ub increase tau.u until the KL increases too
			tau.u.ub 				<- 2*KL_args$tau.u
			curr.it 				<- max.it			
			KL_args$tau.u 			<- tau.u.ub
			KL_args$plot 			<- F
			KL.of.tau.u.ub 			<- do.call(KL_divergence, KL_args)["KL_div"]			
			while (KL.of.tau.u.ub < KL.of.tau.u && curr.it > 0) 
			{
				#print(c(tau.u.ub, KL.of.tau.u.ub , KL.of.tau.u, curr.it))
				curr.it 			<- curr.it - 1
				KL.of.tau.u 		<- KL.of.tau.u.ub
				tau.u.ub 			<- 2 * tau.u.ub
				KL_args$tau.u 		<- tau.u.ub
				KL.of.tau.u.ub 		<- do.call(KL_divergence, KL_args)["KL_div"]
			}
			if(curr.it == 0) 	
				stop("could not find upper bound for tau.u")
			#print(c(tau.u.ub, KL.of.tau.u.ub , KL.of.tau.u, curr.it))
			tau.u.lb				<- ifelse(curr.it==max.it, tau.u.ub/2, tau.u.ub/4)
			#minimize KL_tau.u between [tau.u.lb, tau.u.ub]
			if(verbose)	
				cat(paste("\nupper and lower bounds on tau.u:",tau.u.lb, tau.u.ub))
			if(plot_debug)
				pdf("KL_optimization.pdf", onefile = T)
			KL_args$plot 			<- plot_debug
			KL_args["tau.u"]		<- NULL
			tmp 					<- optimize(kl.optimize, interval = c(tau.u.lb,tau.u.ub), x_name="tau.u" ,KL_divergence= KL_divergence, KL_args = KL_args, verbose = debug)
			if(plot_debug)
				dev.off()
			KL_args$tau.u 	<- tmp$minimum
			KL_args$plot 	<- plot	
			ans 			<- do.call(KL_divergence, KL_args)
		}			
	}
	return(ans)
}
#------------------------------------------------------------------------------------------------------------------------
#' Calibrate the equivalence region for muTOST
#' 
#' Calibrate the equivalence region for the test of location equivalence for given maximum power
#' @export
#' @param mx.pw		maximum power at the point of reference (rho.star).
#' @param df		degrees of freedom
#' @param s.of.T	standard deviation of the test statistic
#' @param tau.up.ub	guess on an upper bound on the upper tolerance of the equivalence region
#' @param alpha		level of the equivalence test
#' @param rho.star	point of reference. Defaults to the point of equality rho.star=0.
#' @param tol		this algorithm stops when the actual maximum power is less than 'tol' from 'mx.pw'
#' @param max.it	this algorithm stops prematurely when the number of iterations to find the equivalence region exceeds 'max.it'
#' @param debug		Flag if C implementation is used.
#' @return	vector of length 4
#' 	\item{1}{lower tolerance of the equivalence region}		
#' 	\item{2}{upper tolerance of the equivalence region}
#' 	\item{3}{actual maximum power associated with the equivalence region}
#' 	\item{4}{error ie abs(actual power - mx.pw)}
#' @examples yn<- 60; ysigma2<- 1; alpha<- 0.01
#'	mutost.calibrate.tolerances(0.9, yn-1, sqrt(ysigma2/yn), 2, alpha )
mutost.calibrate.tolerances<- function(mx.pw, df, s.of.T, tau.up.ub, alpha, rho.star=0, tol= 1e-5, max.it=100, debug=0)
{
	if(!debug)
	{
		suppressWarnings({	#suppress numerical inaccuracy warnings
			ans<- .Call("abc_mutost_calibrate_tauup_for_mxpw",c(mx.pw, df, s.of.T, tau.up.ub, alpha, rho.star, tol, max.it))
		})
		return(ans)
	}
	#else do R implementation
	curr.mx.pw	<- 0
	tau.up.ub	<- tau.up.ub/2
	tmp			<- max.it
	while(curr.mx.pw<mx.pw && tmp>0)
	{
		tmp			<- tmp-1
		tau.up.ub	<- 2*tau.up.ub
#print( list(rho.star, df, tau.up.ub, s.of.T, alpha) )
		curr.mx.pw	<- .Call("abcMuTOST_pow", as.double(rho.star), df, tau.up.ub, s.of.T, alpha)
		#print( curr.mx.pw )
		#curr.mx.pw	<- mutost.pow(rho.star, df, tau.up.ub, s.of.T, alpha)
		#print( curr.mx.pw )		
	}
	if(tmp==0)	stop("mutost.calibrate.tolerances: could not find tau.up.ub")	
#print(tau.up.ub); stop()
		tau.up.lb	<- 0
	error		<- 1	
	while(abs(error)>tol && round(tau.up.lb,d=10)!=round(tau.up.ub,d=10) && max.it>0)
	{
		max.it		<- max.it-1
		tau.up		<- (tau.up.lb + tau.up.ub)/2
		curr.mx.pw	<- .Call("abcMuTOST_pow", as.double(rho.star), df, tau.up, s.of.T, alpha)
		#print( curr.mx.pw )
		#curr.mx.pw	<- mutost.pow(rho.star, df, tau.up, s.of.T, alpha)
		#print( curr.mx.pw )
		error		<- curr.mx.pw - mx.pw
#print(c(curr.mx.pw, tau.up, tau.up.lb, tau.up.ub, max.it))
		if(error<0)
			tau.up.lb<- tau.up
		else
			tau.up.ub<- tau.up
#print(c(abs(error), round(tau.up.lb,d=10)!=round(tau.up.ub,d=10)) )	
	}
	if(max.it==0)	warning("mutost.calibrate.tolerances: reached max.it")
#	stop("HERE")
		c(-tau.up,tau.up,curr.mx.pw,abs(error))
}
#------------------------------------------------------------------------------------------------------------------------
#' Exact TOST for location equivalence.
#' 
#' Perform the exact TOST for location equivalence when the summary values are normally distributed
#' @export
#' @param sim			simulated summary values
#' @param obs			observed summary values
#' @param args			argument that contains the equivalence region and the level of the test (see Examples). This is the preferred method for specifying arguments and overwrites the dummy default values
#' @param verbose		flag if detailed information on the computations should be printed to standard out
#' @param tau.u			upper tolerance of the equivalence region
#' @param tau.l			lower tolerance of the equivalence region
#' @param alpha			level of the equivalence test
#' @param mx.pw			maximum power at the point of equality
#' @param annealing		inflation factor of tolerances of the equivalence region
#' @param sd.tolerance	tolerance on the difference between the standart deviation of the simulated data and the one used in power calibrations. If this difference if too high, the actual maximum of the power might be considerably different from mx.pw.	
#' @param normal.test	name of function with which normality of the summary values is tested
#' @param plot			logical; if \code{TRUE} the summary likelihood and calibrated power are ploted.
#' @param legend.txt 	character; title of the plot (only when \code{plot==TRUE}).
#' @return	vector of length 16 containing
#' \item{lkl}{likelihood of T (test statistic for equality) under the TOST}
#' \item{error}{test statistic, here p-value of TOST}
#' \item{pval}{rescaled p-value of TOST that is expected to follow U(0,1) under the point null hypothesis}
#' \item{link.mc.obs}{mean(obs)}
#' \item{link.mc.sim}{mean(sim)}
#' \item{rho.mc}{mean(sim) - mean(obs)}
#' \item{cil}{lower ABC tolerance, here 0}
#' \item{cir}{upper ABC tolerance, here alpha}
#' \item{tl}{lower tolerance of the equivalence region}
#' \item{tr}{upper tolerance of the equivalence region}
#' \item{al}{lower ABC tolerance, here 0}
#' \item{ar}{upper ABC tolerance, here alpha}
#' \item{pfam.pval}{p-value of the normal test for the simulated summary values}
#' \item{nsim}{number of simulated summary values}
#' \item{mx.pow}{maximum power at the point of equality}
#' \item{rho.pow}{power at the point rho.mc}
#' @example example/ex.mutost.R
mutost.onesample<- function(sim, obs, args= NA, verbose= FALSE, tau.u= 0, tau.l= -tau.u, alpha= 0, mx.pw=0.9, annealing=1, sd.tolerance=0.05, normal.test= "sf.test", plot=0, legend.txt="")
{
	ans <- c(0, 50, 1, NA, NA, NA, 0, 0, 0, 0, 0, 1, 1, NA, NA, NA)
	names(ans)<- c("lkl", "error", "pval","link.mc.obs","link.mc.sim", "rho.mc", "cil", "cir","tl","tr","al","ar","pfam.pval","nsim","mx.pow","rho.pow")

	#compute two sample t-test on either z-scores or untransformed data points
	if(any(is.na(sim)))		
		stop("unexpected NA in sim")
	if(any(is.na(obs)))		
		stop("unexpected NA in obs")
	if(!is.na(args)){			
		# stop("args missing")
		#cat(print.v(sim, prefix="simu", print.char=0))
		#cat(print.v(obs, prefix="obs", print.char=0))
		args<- strsplit(args,'/')[[1]]
		if(length(args)==3)
		{
			annealing	<- as.numeric( args[2] )	#annealing must be at pos 2 otherwise 'ANNEAL.CI.setTau' goes wrong
			obs.n		<- length(obs)
			obs.sd		<- sd(obs)
			tau.u.ub	<- 3*obs.sd
			alpha		<- as.numeric( args[3] )			
		}
		else if(length(args)==4)
		{
			annealing	<- as.numeric( args[2] )
			obs.n		<- length(obs)
			obs.sd		<- sd(obs)
			tau.u.ub	<- as.numeric( args[3] )
			alpha		<- as.numeric( args[4] )			
		}
		else if(length(args)==6)
		{
			annealing	<- as.numeric( args[2] )		
			obs.n		<- as.numeric( args[3] )	#read when obs data just a single point
			obs.sd		<- as.numeric( args[4] )	#read when obs data just a single point
			tau.u.ub	<- as.numeric( args[5] )
			alpha		<- as.numeric( args[6] )			
		}
		else stop("args with unexpected length")
	}else{
		obs.n		<- length(obs)
		obs.sd		<- sd(obs)
		tau.u.ub	<- 3*obs.sd
	}

	if(alpha<0 || alpha>1)		
		stop("incorrect alpha")
	if(tau.u.ub<0 )				
		stop("incorrect args for tau.u.ub")
	if(annealing<1)				
		stop("incorrect annealing parameter")

	args		<- args[1]
	if(verbose)	
		cat(paste("\ninput call:",args,"annealing=",annealing,"obs.sd=",obs.sd,"tau.u.ub=",tau.u.ub,"alpha=",alpha))
	
	# mx.pw		<- 0.9
	
	ans["pfam.pval"]	<-	abccheck.normal(sim, normal.test)	
	if(!any(diff(sim)>0))	
		return(ans)		
	obs.mean	<- mean(obs)	
	if(verbose)		
		cat(paste("\nn=",obs.n,"obs.mean",obs.mean,"obs.sd=",obs.sd,"sim.sd",sd(sim)))	

	if(obs.n<2)	
		stop("length of observed summaries too small, or set 'obs.n' explicitly using args")		
	#trial run on full sd(sim)
	n.sd.refine	<- 8
	tmp			<- sd(sim)
	while(n.sd.refine>0)
	{
		KL_args 		<- list(n.of.x=obs.n, s.of.x= obs.sd, n.of.y=min( obs.n, length(sim) ), s.of.y=tmp, mx.pw=mx.pw, alpha=alpha, tau.u=tau.u.ub, pow_scale=1.5, debug=0)
		abc.param		<- mutost.calibrate( KL_args, 100, debug=0, plot_debug=0, plot=plot, verbose=0)
		if(abc.param["n.of.y"]>length(sim))
		{
			cat(paste("\nnot enough simulated summary values available:",abc.param["n.of.y"], "requested, and used:", length(sim)))
			abc.param	<- mutost.calibrate.tolerances.getkl(obs.n, obs.sd, length(sim), sd(sim), mx.pw, alpha, calibrate.tau.u=TRUE, tau.u=3*obs.sd, debug=0, plot=FALSE)
		}		
#print(abc.param); print(sd(sim[seq_len(abc.param["n.of.y"])])); print(tmp); print(sd.tolerance)
		if( abs( log( sd(sim[seq_len(abc.param["n.of.y"])]) / tmp ) ) > sd.tolerance)
		{
			tmp			<- sd(sim[seq_len(abc.param["n.of.y"])])
			n.sd.refine	<- n.sd.refine-1
		}
		else
			n.sd.refine	<-	0		
		#print(n.sd.refine)
	}
	if( abs( log( sd(sim[seq_len(abc.param["n.of.y"])]) / tmp ) ) > sd.tolerance)
		cat(paste("\nWARNING: actual sd.sim differs from the one used in power calibrations so much that pw.max might be considerably different from",mx.pw, tmp, sd(sim[seq_len(abc.param["n.of.y"])]) ))		
	else if(abs(abc.param["pw.cmx"]-mx.pw)>0.09 && abc.param["n.of.y"]>obs.n  && abc.param["n.of.y"]<length(sim))
	{
		print(abc.param)
		stop("unexpected difference in max power when m>n - tau.up not accurate")
	}
	sim.n		<- abc.param["n.of.y"]		
	sim.mean	<- mean(sim[1:sim.n])
	sim.sd		<- sd(sim[1:sim.n])
	tau.u		<- abc.param["tau.u"]*annealing
	tau.l		<- -tau.u		
	if(verbose)		
		cat(paste("\nsim.mean",sim.mean,"sd sim=",sim.sd,"sd sim long=",tmp,"\n Free ABC parameters calibrated to m=",sim.n,"tau.u=",abc.param["tau.u"],"annealed tau.u=",tau.u))	
	if(plot)
		tmp		<- mutost.calibrate.tolerances.getkl(obs.n, obs.sd, sim.n, sd(sim), mx.pw, alpha, calibrate.tau.u=FALSE, tau.u=abc.param["tau.u"], debug=0, plot=TRUE, legend.title=legend.txt)

	tmp			<- c(	sqrt(sim.n)*(sim.mean-obs.mean-tau.l) / sim.sd,			#[1]	T-	test statistic for -tau (lower test); estimate of the common std dev is simply the std dev in the sample whose sample size is > 1
		sqrt(sim.n)*(sim.mean-obs.mean-tau.u) / sim.sd,			#[2]	T+	test statistic for tau (upper test); estimate of the common std dev is simply the std dev in the sample whose sample size is > 1
		sqrt(sim.n)*(sim.mean-obs.mean) / sim.sd,				#[3]	T	test statistic for equality; estimate of the common std dev is simply the std dev in the sample whose sample size is > 1
		sim.n-1,												#[4] 	degrees of freedom
		sim.sd/sqrt(sim.n),										#[5]	estimate of the std dev of the test statistic is simply the std dev in the sample whose sample size is > 1 divided by that sample size
		sim.sd )												#[6]  	standard deviation of the sample
	tost.ans	<-	generic.tost(tmp, tau.l, tau.u, alpha, tost.distr="t")
	#print("")
	#print(tost.ans)

	ans[c("error","cil","cir")]	<- c(tost.ans["p.error"], 0, alpha)
	ans[c("tl","tr","nsim")]	<- c(tau.l,tau.u,sim.n)
	ans[c("lkl","pval")]<-  tost.ans[c("lkl","ass.pval")]
	ans[c("al","ar")]	<- 	c(0,alpha)								
	#ans["mx.pow"]		<-	mutost.pow(0, tmp[4], tmp[5], tau.u,  alpha)
	ans["mx.pow"]		<-	mutost.pow(0, tmp[4], tmp[5], abc.param["tau.u"],  alpha)		#debugging	
	ans["link.mc.sim"]	<- 	sim.mean
	ans["link.mc.obs"]	<- 	obs.mean
	ans["rho.mc"]		<- 	sim.mean - obs.mean
	ans["rho.pow"]		<-	mutost.pow(ans["rho.mc"], tmp[4], tmp[5], tau.u, alpha )
	if(verbose){	cat("\n"); print(ans)	}
	ans
}