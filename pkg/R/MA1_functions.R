
#------------------------------------------------------------------------------------------------------------------------
#' @title Transform \code{sig2} to the variance error of the MA(1) process
#' @param s2	\code{sig2} parameter
#' @param a		\code{a} parameter
#' @param vx	MLE of the variance of the observed summary values
#' @return Variance error of the MA(1) process
#' @export
ma.sig22rho <- function(s2, a, vx = 1) 
{
	(1 + a * a) * s2/vx
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Transform the variance error of the MA(1) process to \code{sig2}
#' @param x		the variance error of the MA(1) process 
#' @param a		\code{a} parameter
#' @param vx	MLE of the variance of the observed summary values
#' @return \code{sig2} parameter
#' @export
ma.rho2sig2 <- function(x, a, vx = 1) 
{
	x * vx/(1 + a * a)
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Transform \code{a} to the autocorrelation of the MA(1) process
#' @param a	\code{a} parameter
#' @return Autocorrelation parameter of the MA(1) process
#' @export
ma.a2nu<- function(a) 
{
	a/(1+a*a)
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Transform the autocorrelation of the MA(1) process into the \code{a} parameter
#' @param x	Autocorrelation parameter of the MA(1) process
#' @return \code{a} parameter
#' @export
ma.nu2a <- function(x) 
{
	x[x > 0.5] <- 0.5
	x[x < -0.5] <- -0.5
	(1 - sqrt(1 - 4 * x * x))/(2 * x)
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Transform \code{a} to the autocorrelation error of the MA(1) process
#' @param a		\code{a} parameter
#' @param ax	MLE of the autocorrelation of the observed summary values
#' @return Autocorrelation error parameter of the MA(1) process
#' @export
ma.a2rho <- function(a, ax = 0) 
{
	atanh(a/(1 + a * a)) - atanh(ax)
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Transform the autocorrelation error of the MA(1) process to the \code{a} parameter
#' @param x		Autocorrelation error parameter
#' @param ax	MLE of the autocorrelation of the observed summary values
#' @return \code{a} parameter
#' @export
ma.rho2a <- function(x, ax = 0) 
{
	x <- tanh(x + atanh(ax))
	x[x > 0.5] <- 0.5
	x[x < -0.5] <- -0.5
	(1 - sqrt(1 - 4 * x * x))/(2 * x)
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Compute the autocorrelation in a time series along with some other info
#' @export
#' @param x 		numerical vector of time series values
#' @param leave.out thinning, how many values in the pair sequence (x_i,x_i-1) should be left out. Defaults to zero.
#' @param len		number of (x_i,x_i-1) pairs 
#' @return vector of length 3	
#' 	\item{cor}{autocorrelation in the thinned sequence}
#' 	\item{z}{Z-transformation of the autocorrelation (this is atanh of "cor")}
#' 	\item{n}{Number of pairs (x_i,x_i-1) after thinning}
#' @examples ma.cor(rnorm(100,0,1), leave.out=2)
ma.cor<- function(x, leave.out=0, len= ceiling(length(x)/(1+leave.out)) )
{
	tmp	<- rbind( x[-1], x[-length(x)]	)
	tmp	<- tmp[,seq.int(1,ncol(tmp),by=1+leave.out)]
	if(leave.out>0){
		tmp	<- tmp[,seq_len(min(ncol(tmp),len))]		
	}
	if(any(is.na(tmp)))		stop("Unexpected NA in tmp")
	tmp2<- cor(tmp[1,],tmp[2,])
	ans	<- c(tmp2,		.5 * log( (1+tmp2)/(1-tmp2) ), 		ncol(tmp)	)
	names(ans)<- c("cor","z","n")
	ans
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Compute the 2d mode of a 2D density 
#' @description This function computes the 2d mode of a 2D density based on average shifted histograms
#' @param	x			first dimension of random variable 
#' @param	y			second dimension of random variable
#' @param	xlim		range of first dimension of random variable
#' @param	ylim		range of second dimension of random variable
#' @param	xlab		name of first dimension of random variable
#' @param	ylab		name of second dimension of random variable
#' @param	n.hists		number of average shifted histograms
#' @param	nbin		number of bins for histogram
#' @param	nlevels		number of contour levels
#' @param	width.infl	bandwidth inflation of kernel density estimate
#' @param	gridsize	gridsize of kernel density estimate
#' @param	method		either 'kde' or 'ash'
#' @param	plot		Logical
#' @param	contour.col	color of contours
#' @param	cols		color gradient of density estimate	
#' @param	...			further options
#' @import ash
#' @import KernSmooth
#' @import fields
#' @export 
ma.get.2D.mode<- function(x,y,xlim=NA,ylim=NA,xlab='x',ylab='y',n.hists=5,nbin=2, nlevels=5, width.infl=0.25, gridsize=c(100,100), method="kde", plot=0, contour.col="black", cols= head( rev(gray(seq(0,.95,len=trunc(50*1.4)))), 50), ...)
{
	if(any(is.na(xlim)))	xlim<- range(x)*1.05
	if(any(is.na(ylim)))	ylim<- range(y)*1.05
	if(method=="kde")
	{
		bins<- bin2(cbind(x, y), ab=rbind(xlim,ylim),nbin=nbin*c(nclass.Sturges(x),nclass.Sturges(y)))
		f <- ash2(bins,rep(n.hists,2))
		mxidx<- c( (which.max(f$z)-1)%%nrow(f$z)+1, (which.max(f$z)-1)%/%ncol(f$z)+1 ) #row, col
		mx<- c(		mean( f$x[ c(mxidx[1],ifelse(mxidx[1]<length(f$x),mxidx[1]+1,mxidx[1])) ] ),
				mean( f$y[ c(mxidx[2],ifelse(mxidx[2]<length(f$y),mxidx[2]+1,mxidx[2])) ] )		)
		if(plot==1)
		{
			image(f$x,f$y,f$z, col=cols, ...)				
			contour(f$x,f$y,f$z,add=TRUE, nlevels= nlevels, col=contour.col)
			points(mx, col="red", pch=19)
		}
	}
	else if(method=="ash")
	{
		x.bw<- width.infl*diff(summary(x)[c(2,5)])
		y.bw<- width.infl*diff(summary(y)[c(2,5)])
		f <- bkde2D(cbind(x, y), range.x=list(xlim,ylim), bandwidth=c(x.bw,y.bw), gridsize=gridsize)
		mxidx<- c( (which.max(f$fhat)-1)%%nrow(f$fhat)+1, (which.max(f$fhat)-1)%/%ncol(f$fhat)+1 ) #row, col	
		#print(mxidx)
		#print(f$x2[ c(mxidx[2],ifelse(mxidx[2]<length(f$x2),mxidx[2]+1,mxidx[2])) ] )
		mx<- c(		mean( f$x1[ c(mxidx[1],ifelse(mxidx[1]<length(f$x1),mxidx[1]+1,mxidx[1])) ] ),
				mean( f$x2[ c(mxidx[2],ifelse(mxidx[2]<length(f$x2),mxidx[2]+1,mxidx[2])) ] )		)
		if(plot)
		{
			image(f$x1,f$x2,f$fhat, col=cols,xlab=xlab,ylab=ylab )
			contour(f$x1, f$x2, f$fhat, nlevels= nlevels, add=1, col=contour.col, ...)			
		}
	}	
	mx
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Compute the 2d mode of a 2D density 
#' @description This function computes the 2d mode of a 2D density based on average shifted histograms
#' @import KernSmooth
#' @import fields
#' @inheritParams ma.get.2D.mode
#' @export 
ma.add.contour<- function(x,y,xlim=NA,ylim=NA, nlevels=5, width.infl=0.25, gridsize=c(100,100), contour.col="black", ...)
{
	if(any(is.na(xlim)))	xlim<- range(x)*1.05
	if(any(is.na(ylim)))	ylim<- range(y)*1.05	
	x.bw<- width.infl*diff(summary(x)[c(2,5)])
	y.bw<- width.infl*diff(summary(y)[c(2,5)])
	f <- bkde2D(cbind(x, y), range.x=list(xlim,ylim), bandwidth=c(x.bw,y.bw), gridsize=gridsize)	
	contour(f$x1, f$x2, f$fhat, nlevels= nlevels, add=1, col=contour.col, ...)						
}
#------------------------------------------------------------------------------------------------------------------------
#' @title Generate an MA(1) pseudo data set 
#' @description This function generates an MA(1) pseudo data set subject to constraints on the sample autocorrelation, the sample variance and the MLE of the time series.
#' All these statistics of the time series are within a certain tolerance value of the corresponding population level statistics.
#' @param n				length of time series	
#' @param mu			mean of time series
#' @param a				\code{a} parameter
#' @param sd			\code{sigma2} parameter
#' @param leave.out.a	thinning of summary values for the autocorrelation test
#' @param leave.out.s2	thinning of summary values for the variance test
#' @param verbose		print verbose output to the console
#' @param tol			tolerance with which expected MAP is matched in pseudo data
#' @param return_eps_0	logical. If \code{return_eps_0==TRUE}, u0 used to create the time series is returned
#' @return MA(1) time series
#' @export 
ma.get.pseudo.data<- function(n, mu, a, sd, leave.out.a=2, leave.out.s2=1, verbose=0, tol=2e-3, return_eps_0=FALSE )
{
	verbose				<- 1
	m					<- ceiling( max(1, n / 5000) )
	m					<- 1e4 / m
	#true values to match in observed data
	rho.a.0				<- ma.a2rho(a)
	rho.sig2.0			<- sd*sd *( 1 + a*a )
	ans					<- NA
	index.leave.out.s2	<- seq.int(1,n,by=1+leave.out.s2)
	while(any(is.na(ans)))
	{
		x				<- rnorm(m*n+1,mu,sd=sd)
		u0	 			<- x[seq(1,m*n,by=n)]
		x				<- x[-1] + x[-(m*n+1)]*a
		x				<- matrix(x,ncol=m)	
		#adjust variance of time series
		x_std_cte 		<- sapply(seq_len(ncol(x)),function(i)		sqrt(	((1+a*a)*sd*sd)/(var(x[index.leave.out.s2,i])*(n-1)/n)	)	)
		x				<- sapply(seq_len(ncol(x)),function(i)		x[,i]*x_std_cte[i]							)			#set desired variance		
		#compute sample autocorrelation with leave.out
		if(leave.out.a)				
			rho.a.x.lo		<- apply(x,2,function(col)	ma.cor(col, leave.out=leave.out.a)["z"] )
		else
			rho.a.x.lo		<- rep(rho.a.0, ncol(x))
		#compute sample autocorrelation without leave.out 
		rho.a.x.nlo			<- apply(x,2,function(col)	ma.cor(col, leave.out=0)["z"] )
		#compute sample variance without leave.out
		if(leave.out.s2)	
			rho.sig2.x.nlo	<- apply(x,2,function(col) var(col) )
		else
			rho.sig2.x.nlo	<- rep(rho.sig2.0, ncol(x))
		
		error		<- sapply(seq_len(ncol(x)),function(i)
				{ 
					tmp<- arima(x[,i], order=c(0,0,1), include.mean=0, method="CSS-ML")
					c( tmp[["coef"]][1], tmp[["sigma2"]] )
				})
		error		<- t( abs(error-c(a,sd*sd)) )
		colnames(error)	<- c("error.arima.a","error.arima.sd")
		error		<- cbind( data.table( error ), data.table(	error.abc.a.lo=abs(rho.a.x.lo - rho.a.0 ), 
						error.abc.a.nlo=abs(rho.a.x.nlo - rho.a.0 ),
						error.abc.sig2.nlo=abs(rho.sig2.x.nlo - rho.sig2.0),
						u0= u0/x_std_cte), dummy=seq_along(u0) )
		error		<- subset(error, error.arima.a<tol & error.arima.sd<tol &  error.abc.a.lo<tol  &  error.abc.a.nlo<tol & error.abc.sig2.nlo<tol)
		if(nrow(error))
		{
			print(error)
			ans		<- error[1,dummy]
			ans		<- x[,ans]
		}
		else if(verbose)	cat(paste("\nerror above",tol))						
	}	
	
	if(verbose && leave.out.a)	cat(paste("\n atanh(cor) of leave-out time series is",ma.cor(ans, leave.out=leave.out.a)["z"],"and should be",rho.a.0))	
	if(verbose && leave.out.s2)	cat(paste("\n var*(n-1)/n of leave-out time series is",var(ans[index.leave.out.s2])*(n-1)/n,"and should be",rho.sig2.0))
	tmp				<- arima(ans, order=c(0,0,1), include.mean=0, method="CSS-ML")
	if(verbose)	cat(paste("\narima MLE of time series is a=",tmp[["coef"]][1],"sig2=",tmp[["sigma2"]]))
	
	if(return_eps_0)
	{
		unthinned <- list(MLE=data.frame(a=tmp[["coef"]][1],sig2=tmp[["sigma2"]]),s_stat=data.frame(variance= var(ans),autocorr= ma.cor(ans)[["cor"]]))
		thinned <- list(MLE=data.frame(a=NA,sig2=NA),s_stat=data.frame(variance= var(ans[index.leave.out.s2]),autocorr= ma.cor(ans,leave.out=leave.out.a)[["cor"]]))
		ans_error <- as.data.frame(error[1,c(error.arima.a, error.arima.sd, error.abc.a.lo,error.abc.a.nlo,error.abc.sig2.nlo)])
		ans <- list(param=data.frame(a=a,sig2=sd*sd),eps_0=error[1,u0],n=n,thin=data.frame(variance= leave.out.s2,autocorr= leave.out.a), unthinned=unthinned, thinned=thinned, precision=list(tol=tol,error=ans_error),x=ans)
	}
	ans
}
#------------------------------------------------------------------------------------------------------------------------
#' Perform the asymptotic equivalence test for autocorrelations at lag 1
#' @param sim			simulated summary values
#' @param obs			observed summary values
#' @param args			argument that contains the equivalence region and the level of the test (see Examples). This is the preferred method for specifying arguments and overwrites the dummy default values
#' @param verbose		flag if detailed information on the computations should be printed to standard out
#' @param alpha			level of the equivalence test
#' @param leave.out		thinning, how many values in the pair sequence (x_i,x_i-1) should be left out. Defaults to zero.
#' @param normal.test	name of function with which normality of the summary values is tested
#' @return	vector containing
#' \item{error}{test statistic. Here, instead of T we return the p-value of the TOST.}
#' \item{cil}{lower ABC tolerance. Here, instead of c^- we return 0.}
#' \item{cir}{upper ABC tolerance. Here, instead of c^+ we return 'alpha'.}
#' \item{al}{free entry. Here set to c^-.}
#' \item{ar}{free entry. Here set to c^+.}
#' \item{mx.pw}{Maximum power at the point of equality}
#' \item{rho.mc}{sample estimate of 'rho'}
#' @examples leave.out<- 2
#' tau.u<- 0.09
#' alpha<- 0.01
#' n<- 5e3
#' sigma<- 1
#' a<- 0.1
#' args<- paste("acfequiv",leave.out,tau.u,alpha,sep='/')	
#' x<-	rnorm(n+1,0,sigma)
#' x<- x[-1] + x[-(n+1)]*a 
#' y<-	rnorm(n+1,0,sigma)
#' y<- y[-1] + y[-(n+1)]*a
#' ma.equivalence(y,x,args)
ma.equivalence<- function(sim, obs, args=NA, verbose= FALSE, alpha=0, leave.out=0, normal.test= "sf.test")
{
	#verbose<- 0
	if(any(is.na(sim)))	stop("ma.equivalence: error at 1a")
	if(any(is.na(obs)))	stop("ma.equivalence: error at 1b")
	if(length(obs)!=length(sim))		stop("ma.equivalence: error at 1c")
	if(length(sim)<5)	stop("ma.equivalence: error at 1d")
	if(!is.na(args))
	{
		args<- strsplit(args,'/')[[1]]
		if(length(args)==4)	
		{
			leave.out<- as.numeric( args[2] )
			tau.u<- as.numeric( args[3] )
			tau.l<- -tau.u
			alpha<- as.numeric( args[4] )
		}
		else if(length(args)==5)
		{
			leave.out<- as.numeric( args[2] )
			tau.l<- as.numeric( args[3] )
			tau.u<- as.numeric( args[4] )
			alpha<- as.numeric( args[5] )
		}
		else 
			stop("ma.equivalence: error at 1e")
		args<- args[1]
	}
	if(alpha<0 || alpha>1)		stop("ma.equivalence: error at 1f")
	if(tau.u<0 )		stop("ma.equivalence: error at 1g")
	if(tau.l>0 )		stop("ma.equivalence: error at 1h")
	if(leave.out<0)		stop("ma.equivalence: error at 1i")
	ans<- NABC.DEFAULT.ANS	
	ans["pfam.pval"]<-	abccheck.normal(sim,normal.test) 						
	
	#compute z-scores for sim and obs
	z.sim<- ma.cor(sim, leave.out=leave.out)
	z.obs<- ma.cor(obs, leave.out=leave.out)
	z.isd<- sqrt(z.sim["n"]-3)
	#perform TOST -- this now only depends on z.sim z.obs z.isd   
	tmp<- c(z.obs["z"] - z.sim["z"] - tau.l, z.obs["z"] - z.sim["z"] - tau.u, z.obs["z"] - z.sim["z"]) * z.isd #compute ZL and ZU for TOST, add Z for tau=0
	#print(z.sim); print(z.obs); print(z.isd); print(tmp)	
	which.not.reject<- which( c(pnorm( tmp[1] )<1-alpha, pnorm( tmp[2] )>alpha) )
	#only if length==0, TOST will be rejected ie ABC accept
	if(length(which.not.reject)%in%c(0,2))
	{
		#decide which of the two pvalues are to be reported --> figure out which one is closer to boundary
		#both upper and lower test statistics indicate mean difference < -tau and mean difference > tau  	OR 		mean difference >= -tau and mean difference <= tau
		which.not.reject<- which.min(abs(c(1-alpha-pnorm( tmp[1] ), pnorm( tmp[2] )-alpha )))
	}
	#the pvalue of ZU is the lower tail, but for ZL it is the upper tail, so..
	ans["error"]				<- ifelse(which.not.reject==1, 		1-pnorm( tmp[which.not.reject] ),		pnorm( tmp[which.not.reject] ) )
	ans[c("cil","cir")]			<- c(0,alpha)		
	ans[c("tl","tr","nsim")]	<- c(tau.l, tau.u, z.sim["n"]) 
	ans[c("al","ar")]<- 	c(min(tau.l*z.isd+qnorm(1-alpha),0), max(0,tau.u*z.isd+qnorm(alpha)))	#CL and CU of the rejection region of the standardized test statistic
	ans["mx.pow"]<- 		diff(pnorm( ans[c("al","ar")]))		
	ans["pval"]<-			( pnorm(tmp[3],0,1) - (1-ans["mx.pow"])/2 ) / ans["mx.pow"]				#rescaled p-value that is expected to follow U(0,1) under the point null hypothesis
	ans["lkl"]<- 			dnorm(tmp[3],0,1)
	ans["link.mc.sim"]<- 	z.sim["z"]
	ans["link.mc.obs"]<- 	z.obs["z"]
	ans["rho.mc"]<- 		z.sim["z"] - z.obs["z"]
	if(verbose)	cat(paste(paste("\n{",args,"<-list(sim.cor=",z.sim["cor"] ,", obs.cor=",z.obs["cor"] ,", ZL=",tmp[1] ,", ZU=",tmp[2] ,", alpha=",alpha,", tau.l=",tau.l,", tau.u=",tau.u,", ",sep=''),paste(names(ans), ans, collapse=', ', sep='='),")}",sep=''))
	
	ans
}