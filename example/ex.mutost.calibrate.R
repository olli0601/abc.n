yn			<- 60
ymean 		<- 1.3
ysigma 		<- 0.8
sim 		<- rnorm(yn, ymean, ysigma)
sim 		<- (sim - mean(sim))/sd(sim) * ysigma + ymean

# Example 1A: calibrate critical region and power of ABC accept/reject step (default)
# this requires to specify alpha, mx.pw (calibration parameters), 
# n.of.y, s.of.y (summary parameters), 
# tau.u.ub (for numerical optimization)
# note: this is the default calibration because it specifies all ABC parameters
# for sensible calibration parameters, and only requires minimal information on 
# the observed summary values. If a set of observed summary values is available,
# use the calibrations in Example 3.

tau.u.ub 	<- 0.5
ans 		<- mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), 
	tau.u.ub=tau.u.ub, what='MXPW', mw.pw=0.9, alpha=0.01, plot=TRUE)
						
# Example 1B: calibrate critical region and power of ABC accept/reject step
# note: critical region depends on ysigma

ysigma 		<- 1.3
tau.u.ub 	<- 0.5
sim 		<- rnorm(yn, ymean, ysigma)
sim 		<- (sim - mean(sim))/sd(sim) * ysigma + ymean
ans 		<- mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), 
	tau.u.ub=tau.u.ub, what='MXPW', mx.pw=0.9, alpha=0.01, plot=TRUE)
						
# Example 2A: calculate ABC false positive rate for given ABC tolerance
# this requires to specify c.u, tau.u (ad-hoc ABC parameters), n.of.y, s.of.y (summary parameters)
# note: can be useful to compute the ABC false positive rate for uncalibrated ABC routines

cus			<- c(0.1, 0.2, 0.3, 0.4, 0.5)
tau.u		<- 0.6817	
ans 		<- sapply(cus, function(c.u)	mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), c.u=c.u, tau.u=tau.u, what='ALPHA')	) 

# Example 2B: calibrate critical region for given ABC false positive rate and equivalence region
# this requires to specify alpha (calibration parameters), tau.u (ad-hoc ABC parameter), n.of.y, s.of.y (summary parameters)
# note: this is just an intermediate calibration and may result in unsuitable power properties

ans 		<- mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), tau.u=1.2,
	what='CR', alpha=0.01, plot=TRUE)

# Example 3A: calibrate critical region, power of ABC accept/reject step, and #simulated data points
# this requires to specify alpha, mx.pw (calibration parameters), 
# n.of.x, s.of.x, s.of.y (summary parameters), 
# tau.u.ub (for numerical optimization)
# note: the advantage here is that the KL divergence is also minimised, but we also
# need to have a set of observed summary values
#
# case sd(sim)>sd(obs): the power function is not "too tight" for yn=xn.

xn 			<- 60
xmean 		<- 1
xsigma 		<- 1
obs 		<- rnorm(xn, xmean, xsigma)
obs 		<- (obs - mean(obs))/sd(obs) * xsigma + xmean
yn			<- NA	#this is now a calibration output
ans 		<- mutost.calibrate(n.of.x=length(obs), s.of.x= sd(obs), n.of.y=length(sim), s.of.y=sd(sim), tau.u.ub=5,
	what='KL', mx.pw=0.9, alpha=0.01, plot=TRUE, debug=TRUE)

# case sd(sim)<sd(obs): the power function is "too tight" for yn=xn.

ysigma 		<- 0.7
sim 		<- rnorm(yn, ymean, ysigma)
sim 		<- (sim - mean(sim))/sd(sim) * ysigma + ymean
ans 		<- mutost.calibrate(n.of.x=length(obs), s.of.x= sd(obs), n.of.y=length(sim), s.of.y=sd(sim), 
	tau.u.ub=3, what='KL', mx.pw=0.9, alpha=0.01, plot=TRUE, debug=TRUE)

\dontrun{
abc.presim.uprior.mu<- function(abc.nit, xn, xmean, xsigma, prior.l, prior.u, ysigma, yn=NA )		
{		
	ans			<- vector("list",5)
	names(ans)	<- c("x","xn","xmean","xsigma","sim")
	obs 		<- rnorm(xn, xmean, xsigma)
	obs 		<- (obs - mean(obs))/sd(obs) * xsigma + xmean
	ans[["x"]]			<- obs
	ans[["xmean"]]		<- xmean
	ans[["xsigma"]]		<- xsigma
	
	ans[["sim"]]		<- sapply(1:abc.nit, function(i)
			{					
				ymu		<- runif(1, prior.l, prior.u)
				y		<- rnorm(yn, ymu, sd=ysigma)
				tmp		<- c(yn, ymu, ysigma, mean(y), sd(y) )									
				tmp					
			})								
	rownames(ans[["sim"]])	<- c('m','ymu','ysigma','ysmean','yssd')
	ans
}


abc.presim	<- abc.presim.uprior.mu( 	abc.nit=1e7, xn=60, xmean=1.34, xsigma=1.4, 
		prior.l=1.34-5, prior.u=1.34+5, ysigma=1.4, yn=60 )
abc.df		<- as.data.table(t(abc.presim$sim))								
abc.df[, it:=seq_len(nrow(abc.df))]								
tmp			<- abc.df[,  as.list( mutost.calibrate(	n.of.y=m, s.of.y=yssd, tau.u.ub=1, what='MXPW', mx.pw=0.9, alpha=0.01)[1:4] ), by='it']
abc.df		<- merge(abc.df, tmp, by='it')								
abc.df[, T:= abc.df[, ysmean-1.34]]
abc.df[, ABC_C90:= abc.df[, c.l<=T & T<=c.u]]

ggplot( subset(abc.df, ABC_C90), aes(x=ymu-abc.presim$xmean)) + geom_histogram()
ggplot( subset(abc.df, ABC_C90), aes(x=ymu-abc.presim$xmean)) + geom_histogram(aes(y= ..density..))
}


